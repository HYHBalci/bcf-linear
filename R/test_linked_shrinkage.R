library(Rcpp)
sourceCpp("src/horseshoe_samplers.cpp")
set.seed(123)
library(MASS) 

rinvgamma <- function(shape, scale){
  if (shape <= 0.0 || scale <= 0.0) {
    stop("Shape and scale must be positive.");
  }
  
  g = rgamma(1, shape, scale)
  out = 1.0 / g
  return(out)
}

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

sample_linear_part <- function(y, z, X, intTreat = TRUE, iter = 1000, burnin = 500, horseshoe = FALSE) {
  n <- length(y)
  p <- ncol(X)
  
  # Generate interaction terms 
  int_pairs <- expand.grid(1:p, 1:p)
  int_pairs <- int_pairs[int_pairs$Var1 < int_pairs$Var2, ]  # Avoid duplicates
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
  tau_glob <- 1.0
  beta <- rep(0.0, p)
  if(horseshoe){
    tau <- rep(1.0, p + p_int)
  } else{ 
    tau <- rep(1.0, p)
  }
  beta_int <- rep(0.0, p_int)
  tau_int <- 0.5
  
  nu <- 1 / (tau^2)
  xi <- 0.01
  
  # Storage for posterior samples
  alpha_samples <- numeric(iter)
  xi_samples <- numeric(iter)
  tau_glob_samples <- numeric(iter)
  beta_samples <- matrix(0, iter, p)
  if(!horseshoe){
    tau_samples <- matrix(0, iter, p)
  } else {
    tau_samples <- matrix(0, iter, p + p_int)
  }
  nu_samples <- matrix(0, iter, p + p_int)
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
      beta[j] <- sample_beta_j(n, r_beta, z, X[, j], 1 / (tau[j]*tau_glob), sigma)
      if(horseshoe){
        tau[j] <- sqrt(rinvgamma(1, (1 / nu[j]) + (beta[j] * beta[j]) / (2*tau_glob*tau_glob*sigma*sigma)));
      } else{
        tau[j] <- sample_tau_j_slice(tau[j], beta[j], j, beta_int = beta_int, tau = tau, tau_int = tau_int, sigma = sigma)
      }
    }
    
    # Sample interaction effects (if enabled)
    if (intTreat) {
      for (k in 1:p_int) {
        iVar <- int_pairs$Var1[k]
        jVar <- int_pairs$Var2[k]
        r_beta_int <- y - (z * alpha + X %*% beta + X_int[, -k] %*% beta_int[-k])
        if(horseshoe){
          beta_int[k] <- sample_beta_j(n, r_beta_int, z, X[, iVar] * X[, jVar], tau[p+k]*tau_glob, sigma)
        } else {
          beta_int[k] <- sample_beta_j(n, r_beta_int, z, X[, iVar] * X[, jVar], sqrt(tau_int * tau[iVar] * tau[jVar]), sigma)
        }
        if(horseshoe){
          for(o in 1:(p + p_int)){
            nu[o] <- rinvgamma(1, 1 + 1 / (tau[o]*tau[o]))
          } 
          tau[p+k] <- sqrt(rinvgamma(1, (1 / nu[p + k]) + (beta_int[k] * beta_int[k]) / (2*tau_glob*tau_glob*sigma*sigma)));
        }
      }
      
      if(!horseshoe) {
      # Sample tau_int (global shrinkage for interactions)
      currentTauInt <- tau_int
      proposedTauInt <- runif(1, 0.01, 1.0)
      
      logPosteriorCurrent <- loglikeTauInt(currentTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      logPosteriorProposed <- loglikeTauInt(proposedTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      logAccept <- logPosteriorProposed - logPosteriorCurrent
      if (log(runif(1)) < logAccept) tau_int <- proposedTauInt
      } else {
        sum = 0
        for(j in 1:(p + p_int)){
          if(j <= p){
            sum = sum +  beta[j]^2 / tau[j]^2
          } else{
            sum = sum + beta_int[j-p]^2 / tau[j]^2
          }
        }
        xi <- rinvgamma(1, 1 + 1 / (tau_glob*tau_glob))
        tau_glob <- sqrt(rinvgamma((p + p_int +1)/2, 1/xi +  (1/(2*sigma*sigma))*sum))
        
      }
    }
    
    # Sample sigma
    resid <- y - (z * alpha + X %*% beta + X_int %*% beta_int)
    
    if(horseshoe){
      # Lambda matrix: Λ = diag(lambda1^2, ..., lambdap^2)
      lambda_sq <- tau^2  # assume 'lambda' is a vector of lambda_j
      Lambda_inv <- diag(1 / lambda_sq)
      beta_tot <- c(beta, beta_int)
      # Compute scale parameter
      scale <- as.numeric(0.5 * crossprod(resid) + 0.5 * t(beta_tot) %*% Lambda_inv %*% beta_tot)
      
      # Shape parameter
      shape <- (n + p + p_int) / 2
      
      # Sample from inverse-gamma
      sigma <- sqrt(rinvgamma(shape, scale))
    } else {
      sigma <- sqrt(sample_sigma2_ig(n, resid,1.0, 0.001))
    }
    
    # Store samples
    if (i > burnin) {
      idx <- i - burnin
      alpha_samples[idx] <- alpha
      beta_samples[idx, ] <- beta
      tau_glob_samples[idx] <- tau_glob
      nu_samples[idx, ] <- nu
      xi_samples[idx] <- xi 
      print(tau)
      tau_samples[idx, ] <- tau
      beta_int_samples[idx, ] <- beta_int
      tau_int_samples[idx] <- tau_int
      sigma_samples[idx] <- sigma
    }
  }
  
  return(list(alpha = alpha_samples, beta = beta_samples, tau = tau_samples,
              beta_int = beta_int_samples, tau_int = tau_int_samples, sigma = sigma_samples, xi = xi_samples, tau_glob = tau_glob_samples, nu = nu_samples))
}

