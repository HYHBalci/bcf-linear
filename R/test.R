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
    # hyperparameters
    sigma_prior_nu   = 3.0,  # e.g. ~ InvGamma(nu/2, nu*lambda/2)
    sigma_prior_lambda = 1.0,
    beta0_prior_sd   = 10.0  # for intercept ~ N(0, beta0_prior_sd^2)
) {
  # 0) Basic checks
  N <- length(y)
  if(length(x) != N) stop("x and y must have the same length")
  
  # 1) Storage
  # We'll store draws for (beta0, beta, tau, sigma)
  out_beta0 <- numeric(n_iter)
  out_beta  <- numeric(n_iter)
  out_tau   <- numeric(n_iter)
  out_sigma <- numeric(n_iter)
  
  # 2) Initialize parameters
  beta0 <- 0.0
  beta  <- 0.0
  tau   <- 1.0
  sigma <- 1.0
  
  # 3) The main MCMC loop
  for(iter in seq_len(n_iter)) {
    
    # --- (A) Update beta0 (intercept) ----------------------------------
    # Model: y_i = beta0 + beta*x_i + e_i,  e_i ~ N(0, sigma^2).
    # If we consider partial residuals wrt beta*x_i, i.e.:
    #    r_i = y_i - beta*x_i,
    # We have r_i ~ Normal(beta0, sigma^2).
    # With prior beta0 ~ N(0, beta0_prior_sd^2),
    # the posterior is also normal (conjugate).
    #
    # posterior variance = 1 / (n/sigma^2 + 1/beta0_prior_var)
    # posterior mean    = var * ( (sum(r_i)/sigma^2) )
    r_i <- y - beta*x  # partial residual for intercept
    prior_var0 <- beta0_prior_sd^2
    prec_data  <- N / (sigma^2)
    prec_prior <- 1.0 / prior_var0
    postVar0   <- 1.0 / (prec_data + prec_prior)
    postMean0  <- postVar0 * (sum(r_i) / (sigma^2))
    beta0 <- rnorm(1, mean = postMean0, sd = sqrt(postVar0))
    
    # --- (B) Update beta (slope) via sample_beta_j ---------------------
    # We'll treat beta as if we have "z[i] = 1" for all i (like the "treatment" is always 1),
    # and w_j[i] = (x[i]), then partial residual is: r_beta[i] = y[i] - beta0.
    # So r_beta[i] ~ z[i]*w_j[i]*beta + e_i with prior beta ~ Normal(0, sigma^2 tau^2).
    r_beta <- y - beta0  # partial residual ignoring beta's contribution
    # We'll create dummy "z" (all 1) and "w_j" (the x vector)
    z_vec  <- rep(1.0, N)
    w_j    <- x
    # call the Rcpp function sample_beta_j
    beta_new <- sample_beta_j(
      N         = N, 
      r_beta    = r_beta, 
      z         = z_vec, 
      w_j       = w_j,
      tau_j     = tau, 
      sigma     = sigma
    )
    beta <- beta_new
    
    # --- (C) Update tau (local scale) via slice sampler ---------------
    # We call sample_tau_j_slice(tau_old, beta, sigma)
    tau_new <- sample_tau_j_slice(
      tau_old = tau, 
      beta_j  = beta, 
      sigma   = sigma
    )
    tau <- tau_new
    
    # --- (D) Update sigma^2 via Inverse-Gamma (conjugate) -------------
    # residual e_i = y[i] - [beta0 + beta*x[i]]
    e <- y - (beta0 + beta*x)
    rss <- sum(e^2)
    # prior: sigma^2 ~ InvGamma(nu/2, nu*lambda/2)
    # posterior shape = (nu + N)/2
    # posterior rate  = (nu*lambda + rss)/2
    shape_post <- 0.5 * (sigma_prior_nu + N)
    rate_post  <- 0.5 * (sigma_prior_nu*sigma_prior_lambda + rss)
    # sample: sigma^2 from IG(shape, rate)
    # in R: 1 / rgamma(1, shape, rate=...)
    sigma2_new <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    sigma <- sqrt(sigma2_new)
    
    # 4) Store draws
    out_beta0[iter] <- beta0
    out_beta[iter]  <- beta
    out_tau[iter]   <- tau
    out_sigma[iter] <- sigma
  }
  
  # 5) Return posterior draws (minus burn-in)
  keep <- seq(from = burn+1, to = n_iter)
  res <- list(
    beta0 = out_beta0[keep],
    beta  = out_beta[keep],
    tau   = out_tau[keep],
    sigma = out_sigma[keep]
  )
  return(res)
}

# =======================================================================
# Example usage:
# =======================================================================
set.seed(123)
# Generate synthetic data: y = 1 + 2*x + e
N <- 10000
x <- rnorm(N, 0, 1)
true_beta0 <- 1.0
true_beta  <- 2.0
true_sigma <- 1.0
y <- true_beta0 + true_beta*x + rnorm(N, 0, true_sigma)

# Run MCMC
res <- run_mcmc_univariate_horseshoe(y, x, n_iter=3000, burn=1000)

cat("Posterior means:\n")
cat("beta0 =", mean(res$beta0), "\n")
cat("beta  =", mean(res$beta),  "\n")
cat("tau   =", mean(res$tau),   "\n")
cat("sigma =", mean(res$sigma), "\n")
