library(Rcpp)
Rcpp::sourceCpp("src/horseshoe_samplers.cpp")

# =======================================================================
# MCMC routine for a univariate regression with Horseshoe prior on the slope
# =======================================================================

run_mcmc_univariate_horseshoe <- function(
    y,          # outcome vector
    x,          # predictor vector (same length as y)
    n_iter = 2000, 
    burn   = 500,
    thin = 10,
    # hyperparameters
    # Priors for intercept
    alpha_prior_sd = 10.0,
    # Inverse-Gamma prior for sigma^2
    sigma2_prior_shape = 1.0,
    sigma2_prior_rate  = 0.001,
    # Slice sampler settings for tau
    step_out = 0.5,
    max_steps = 50
) {
  # 0) Basic checks
  N <- length(y)

  # 1) Storage
  # We'll store draws for (beta0, beta, tau, sigma)
  out_alpha <- numeric(n_iter)
  out_beta  <- matrix(0, nrow = n_iter, ncol = p)
  out_tau   <- matrix(0, nrow = n_iter, ncol = p)
  out_lambda <- matrix(1.0, nrow = n_iter, ncol = p)
  out_sigma <- numeric(n_iter)
  
  # 2) Initialize parameters
  alpha <- 0.0
  beta  <- rep(0.0, p) # Local shrinkage for each beta
  tau   <- rep(1.0, p)           # Global shrinkage
  sigma <- 1.0           # Noise standard deviation
  
  # 3) MCMC Loop
  for (iter in seq_len(n_iter)) {
    # -------------------------------------------
    # (A) Update alpha (intercept)
    # Residual excluding alpha
    r_alpha <- y - (X %*% beta)
    
    alpha <- sample_alpha(
      N          = N,
      r_alpha    = r_alpha,
      z_ = rep(1.0, N),
      sigma      = sigma,
      alpha_prior_sd = alpha_prior_sd
    )
    
    # -------------------------------------------
    # (B) Update each beta_j
    for (j in seq_len(p)) {
      # Partial residual excluding beta_j contribution
      r_beta <- y - alpha - X[, -j, drop=FALSE] %*% beta[-j]
      
      # Sample beta_j using Rcpp function
      beta[j] <- sample_beta_j(
        N        = N,
        r_beta   = r_beta,
        z        = rep(1.0, N),  # Set all z to 1 for continuous predictors
        w_j      = X[, j],
        tau_j    = tau[j],  # Local-global shrinkage
        sigma    = sigma
      )
      
      # -------------------------------------------
      
    # -------------------------------------------
    # (D) Update global scale tau
    tau[j] <- sample_tau_j_slice(
      tau_old = tau[j],
      beta_j  = beta[j],
      sigma   = sigma,
      step_out = step_out,
      max_steps = max_steps
    )
    }
    # -------------------------------------------
    # (E) Update sigma^2 (error variance)
    resid <- y - (alpha + X %*% beta)
    sigma2 <- sample_sigma2_ig(
      N           = N,
      resid       = resid,
      shape_prior = sigma2_prior_shape,
      rate_prior  = sigma2_prior_rate
    )
    sigma <- sqrt(sigma2)
    
    # -------------------------------------------
    # Store results
    out_alpha[iter] <- alpha
    out_beta[iter, ] <- beta
    out_tau[iter, ]  <- tau
    out_sigma[iter] <- sigma}
  
  # 4) Return posterior draws after burn-in and thinning
  keep_indices <- seq(from = burn + 1, to = n_iter, by = thin)
  
  return (list(
    alpha  = out_alpha[keep_indices],
    beta   = out_beta[keep_indices, , drop = FALSE],
    lambda = out_lambda[keep_indices, , drop = FALSE],
    tau    = out_tau[keep_indices],
    sigma  = out_sigma[keep_indices]
  ))
}

# =======================================================================
# Example usage:
# =======================================================================
# Set up the simulation
set.seed(123)
N <- 200

true_beta0 <- 1.0
true_beta  <- c(2.0, 3.0, 0.0, -1, 0, 0, 4, 0)
p <- length(true_beta)

# Generate predictor matrix X
X <- matrix(rnorm(N * p, mean = 0, sd = 1), nrow = N, ncol = p)

# Generate the response variable y
true_sigma <- 1.0
y <- true_beta0 + X %*% true_beta + rnorm(N, 0, true_sigma)

# Run MCMC
res <- run_mcmc_univariate_horseshoe(y, X, n_iter = 50000, burn = 1000)

# Display results
cat("Posterior mean of alpha:", mean(res$alpha), "\n")
cat("Posterior mean of beta:\n")
print(colMeans(res$beta))
cat("Posterior mean of tau:", mean(res$tau), "\n")
cat("Posterior mean of sigma:", mean(res$sigma), "\n")

# Visualize results
hist(res$beta[,4], main = "Posterior distribution of sigma", xlab = "sigma")