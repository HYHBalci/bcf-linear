library(Rcpp)
sourceCpp("src/horseshoe_samplers.cpp")
set.seed(123)
library(MASS) 
loglikeTauInt <- function(
    tau_int,
    # Baseline interaction info
    beta_int_base,
    int_pairs_base,        # Typically a data.frame or matrix with columns (Var1, Var2)
    # Main-effect shrinkage scales (tau) used in var = sigma^2 * tau_int * tau[iVar] * tau[jVar]
    tau_main,
    sigma,
    include_treatment_int = FALSE,
    beta_int_trt = numeric(0),
    tau_trt      = numeric(0),   # Not used directly here, but included for consistency
    int_pairs_trt = data.frame()
) {
  # 1) Check Uniform(0.01, 1.0) boundary
  if (tau_int < 0.01 || tau_int > 1.0) {
    return(-Inf)
  }
  
  logp    <- 0.0
  log2pi  <- log(2.0 * pi)
  
  # 2) Baseline interaction terms
  # var_ij = tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
  # log N(0, var_ij) = -0.5 * [ log(2π var_ij) + beta^2 / var_ij ]
  for (k in seq_len(nrow(int_pairs_base))) {
    iVar <- int_pairs_base[k, 1]  # or int_pairs_base$Var1[k], if it's a data.frame with named columns
    jVar <- int_pairs_base[k, 2]  # or int_pairs_base$Var2[k]
    
    var_ij <- tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
    beta2  <- beta_int_base[k]^2
    
    logp <- logp -
      0.5 * (log2pi + log(var_ij)) -
      0.5 * (beta2 / var_ij)
  }
  
  # 3) Treatment interaction terms (if included)
  if (include_treatment_int && length(beta_int_trt) > 0) {
    for (k in seq_len(nrow(int_pairs_trt))) {
      iVar <- int_pairs_trt[k, 1]
      jVar <- int_pairs_trt[k, 2]
      
      var_ij <- tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
      beta2  <- beta_int_trt[k]^2
      
      logp <- logp -
        0.5 * (log2pi + log(var_ij)) -
        0.5 * (beta2 / var_ij)
    }
  }
  
  return(logp)
}

sample_linear_part <- function(y, z, X, intTreat = TRUE, iter = 1000, burnin = 500) {
  n <- length(y)
  p <- ncol(X)
  
  # Generate interaction terms
  int_pairs <- expand.grid(1:p, 1:p)
  int_pairs <- int_pairs[int_pairs$Var1 <= int_pairs$Var2, ]  # Avoid duplicates
  p_int <- nrow(int_pairs)
  
  X_int <- do.call(
    cbind,
    lapply(seq_len(nrow(int_pairs)), function(k) {
      X[, int_pairs$Var1[k]] * X[, int_pairs$Var2[k]]
    })
  )
  
  # Initialize parameters
  alpha <- 0.0
  sigma <- 1
  beta <- rep(0.0, p)
  tau <- rep(1.0, p)
  beta_int <- rep(0.0, p_int)
  tau_int <- 0.5
  
  # Storage for posterior samples
  alpha_samples <- numeric(iter)
  beta_samples <- matrix(0, iter, p)
  tau_samples <- matrix(0, iter, p)
  beta_int_samples <- matrix(0, iter, p_int)
  tau_int_samples <- numeric(iter)
  sigma_samples <- numeric(iter)
  
  for (i in 1:(burnin+iter)) {
    # Sample alpha (intercept)
    r_alpha <- y - (z * alpha + X %*% beta + X_int %*% beta_int)
    alpha <- sample_alpha(n, r_alpha, z, sigma, alpha_prior_sd = 10.0)
    
    # Sample main effects (beta)
    for (j in 1:p) {
      r_beta <- y - (z * alpha + X[, -j] %*% beta[-j] + X_int %*% beta_int)
      beta[j] <- sample_beta_j(n, r_beta, z, X[, j], tau[j], sigma)
      tau[j] <- sample_tau_j_slice(tau[j], beta[j], sigma)
    }
    
    # Sample interaction effects (if enabled)
    if (intTreat) {
      for (k in 1:p_int) {
        iVar <- int_pairs$Var1[k]
        jVar <- int_pairs$Var2[k]
        r_beta_int <- y - (z * alpha + X %*% beta + X_int[, -k] %*% beta_int[-k])
        beta_int[k] <- sample_beta_j(n, r_beta_int, z, X[, iVar] * X[, jVar], sqrt(tau_int * tau[iVar] * tau[jVar]), sigma)
      }
      
      # Sample tau_int (global shrinkage for interactions)
      currentTauInt <- tau_int
      proposedTauInt <- runif(1, 0.01, 1.0)
      
      logPosteriorCurrent <- loglikeTauInt(currentTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      logPosteriorProposed <- loglikeTauInt(proposedTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      print(logPosteriorProposed)
      logAccept <- logPosteriorProposed - logPosteriorCurrent
      if (log(runif(1)) < logAccept) tau_int <- proposedTauInt
    }
    
    # Sample sigma
    resid <- y - (z * alpha + X %*% beta + X_int %*% beta_int)
    sigma <- sqrt(sample_sigma2_ig(n, resid,1.0, 0.001))
    
    # Store samples
    if (i > burnin) {
      idx <- i - burnin
      alpha_samples[idx] <- alpha
      beta_samples[idx, ] <- beta
      tau_samples[idx, ] <- tau
      beta_int_samples[idx, ] <- beta_int
      tau_int_samples[idx] <- tau_int
      sigma_samples[idx] <- sigma
    }
  }
  
  return(list(alpha = alpha_samples, beta = beta_samples, tau = tau_samples,
              beta_int = beta_int_samples, tau_int = tau_int_samples, sigma = sigma_samples))
}


set.seed(123)

n <- 500    # Number of observations
p <- 5      # Number of predictors

# 1) Generate design matrix X ~ N(0,1)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)

# 2) Generate *only* strictly upper (or lower) triangular interactions
int_pairs <- expand.grid(1:p, 1:p)
int_pairs <- int_pairs[int_pairs$Var1 < int_pairs$Var2, ]
# Now we have 10 pairs for p=5
X_int <- do.call(
  cbind,
  lapply(seq_len(nrow(int_pairs)), function(k) {
    X[, int_pairs$Var1[k]] * X[, int_pairs$Var2[k]]
  })
)

# 3) True interaction coefficients (10 total)
beta_int_true <- c(0.5, -0.7, 0.3, 0.0, 0.0, -0.2, 0.0, 0.4, -0.5, 0.6)

# 4) True main-effect coefficients
beta_true <- c(1.5, -2.0, 0.8, 0.0, 1.2)  # length = p = 5

# 5) Treatment assignment
z <- rbinom(n, 1, 0.5)

# 6) Intercept, noise
alpha_true <- 2.0
sigma_true <- 1.0

# 7) Generate y
y <- alpha_true +
  z * 1.0 +                       # optional treatment effect
  X %*% beta_true +               # main effects
  X_int %*% beta_int_true +       # interactions
  rnorm(n, 0, sigma_true)

# 9) Run the sampler
result <- sample_linear_part(y, z, X, intTreat = TRUE, iter = 5000, burnin = 2000)

hist(result$tau_int)

### check with coda 

library(coda)

chain_1 <- as.mcmc(result$beta)
effectiveSize(chain_1)