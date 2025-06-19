# --- Simulation Data Generator for Simple Logistic Model ---
simulate_simple_logistic_data <- function(n = 500, p = 2, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate covariates X (e.g., p continuous covariates)
  X_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X_mat) <- paste0("X", 1:p)
  
  # Generate binary treatment Z (0/1) - simple RCT for this example
  Z_vec <- rbinom(n, 1, 0.5)
  
  # Define true parameters
  true_alpha <- 0.2
  true_beta <- rnorm(p, mean = 0, sd = 0.75) # True coefficients for X main effects
  true_aleph <- 0.5                          # True coefficient for Z main effect (additive to X effects)
  true_gamma <- rnorm(p, mean = 0, sd = 0.5) # True coefficients for Z*X interaction terms
  
  # Construct linear predictor
  # eta = alpha + X*beta + Z*aleph + Z*(X*gamma)
  eta <- true_alpha + 
    as.vector(X_mat %*% true_beta) + 
    Z_vec * true_aleph + 
    as.vector((X_mat * Z_vec) %*% true_gamma) # Element-wise product for Z*(X*gamma) then sum
  
  # Calculate probabilities
  prob_y <- 1 / (1 + exp(-eta))
  
  # Generate binary outcome Y
  y_vec <- rbinom(n, 1, prob_y)
  
  # Return data and true parameters
  list(
    y_vec = y_vec,
    X_mat = X_mat,
    Z_vec = Z_vec,
    true_params = list(
      alpha = true_alpha,
      beta = true_beta,
      aleph = true_aleph,
      gamma = true_gamma
    ),
    true_eta = eta,
    true_prob_y = prob_y
  )
}

#-------------------------------------------------------------
# --- Simplified Logistic Gibbs Sampler (Pure R Version) ---
simple_logistic_gibbs_R <- function(y_vec, X_mat, Z_vec,
                                    n_iter, burn_in,
                                    prior_mean = 0, # Common prior mean for all coeffs
                                    prior_sd = 10.0,  # Common prior SD for all coeffs
                                    init_alpha = 0.0,
                                    init_aleph = 0.0,
                                    seed = 456) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!requireNamespace("BayesLogit", quietly = TRUE)) {
    stop("Package 'BayesLogit' needed for rpg function. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for mvrnorm. Please install it.", call. = FALSE)
  }
  
  N <- nrow(X_mat)
  p <- ncol(X_mat)
  
  if (length(y_vec) != N) stop("Length of y_vec must match number of rows in X_mat.")
  if (length(Z_vec) != N) stop("Length of Z_vec must match number of rows in X_mat.")
  
  # --- Construct the full design matrix for X_full * all_coeffs ---
  # Columns: Intercept, X1..Xp (for beta), Z (for aleph), Z*X1..Z*Xp (for gamma)
  X_for_beta <- X_mat
  X_for_gamma <- X_mat * Z_vec # Element-wise multiplication, Z_vec will be recycled row-wise
  
  X_full <- cbind(
    Intercept = rep(1, N),   # For alpha
    X_for_beta,              # For beta_vec
    Z_Intercept = Z_vec,     # For aleph
    X_for_gamma              # For gamma_vec
  )
  
  # Check for perfect multicollinearity early (simplistic check)
  if (qr(X_full)$rank < ncol(X_full)) {
    warning("Design matrix X_full might be rank deficient. Consider checking for multicollinearity.")
  }
  
  num_coeffs_total <- 1 + p + 1 + p # alpha, beta_vec, aleph, gamma_vec
  
  # --- Initialize Parameters ---
  alpha <- init_alpha
  beta_vec <- rep(0, p)
  aleph <- init_aleph
  gamma_vec <- rep(0, p)
  
  all_coeffs <- c(alpha, beta_vec, aleph, gamma_vec)
  
  # --- Priors ---
  # Assuming independent Normal priors N(prior_mean, prior_sd^2) for all coefficients
  prior_precision_scalar <- 1 / (prior_sd^2)
  P_inv_diag_mat <- diag(prior_precision_scalar, num_coeffs_total)
  # If specific prior means are needed for P_inv_beta0:
  prior_mean_vector <- rep(prior_mean, num_coeffs_total)
  P_inv_beta0 <- prior_precision_scalar * prior_mean_vector
  
  
  # --- MCMC Storage ---
  num_samples_to_store <- n_iter - burn_in
  if (num_samples_to_store <= 0) {
    stop("n_iter must be greater than burn_in.")
  }
  
  alpha_samples <- numeric(num_samples_to_store)
  beta_samples <- matrix(0, nrow = num_samples_to_store, ncol = p)
  aleph_samples <- numeric(num_samples_to_store)
  gamma_samples <- matrix(0, nrow = num_samples_to_store, ncol = p)
  
  current_sample_idx <- 0
  
  # --- MCMC Loop ---
  for (iter in 1:n_iter) {
    if (iter %% 100 == 0) message(paste("R Sampler - Iteration:", iter, "/", n_iter))
    
    # 1. Update Linear Predictor eta_vec based on current coefficients
    eta_vec <- as.vector(X_full %*% all_coeffs)
    
    # 2. Sample Polya-Gamma latent variables omega_vec
    # Ensure eta_vec has no NAs or Infs which can cause rpg to fail
    if(any(is.na(eta_vec)) || any(is.infinite(eta_vec))) {
      warning(paste("NA or Inf in eta_vec at iteration", iter, "- using previous omega_vec or stopping"))
      # Potentially skip this iteration or use previous omega_vec if robust handling is needed
      # For simplicity, we might error out or use previous omega if it exists from iter > 1
      if(iter == 1) stop("NA/Inf in eta_vec at first iteration.")
    } else {
      omega_vec <- BayesLogit::rpg(num = N, h = 1, z = eta_vec) # abs() not needed for z here
    }
    omega_vec[omega_vec < 1e-9] <- 1e-9 # Safeguard for numerical stability
    
    # 3. Construct transformed response Y_star_vec for Gaussian likelihood
    kappa_vec <- y_vec - 0.5
    Y_star_vec <- kappa_vec / omega_vec
    
    # 4. Sample all regression coefficients (alpha, beta, aleph, gamma) jointly
    Omega_diag_mat <- diag(omega_vec) # Observation-specific precisions
    
    # Posterior precision: X_full' * Omega * X_full + Prior_Precision
    posterior_precision_mat <- t(X_full) %*% Omega_diag_mat %*% X_full + P_inv_diag_mat
    
    # Posterior mean component from data: X_full' * Omega * Y_star
    # Add prior mean component: P_inv_diag_mat * prior_mean_vector (or just P_inv_beta0)
    rhs_for_mean <- t(X_full) %*% Omega_diag_mat %*% Y_star_vec + P_inv_beta0
    
    # Attempt Cholesky decomposition for stability
    chol_attempt <- try(chol(posterior_precision_mat), silent = TRUE)
    if (inherits(chol_attempt, "try-error")) {
      jitter_val <- 1e-6 * mean(diag(posterior_precision_mat)) # Relative jitter
      posterior_precision_mat_jittered <- posterior_precision_mat + diag(jitter_val, num_coeffs_total)
      chol_attempt <- try(chol(posterior_precision_mat_jittered))
      if (inherits(chol_attempt, "try-error")) {
        warning(paste("Cholesky decomposition failed at iteration", iter, "even with jitter. Using previous coefficients."))
        # Skip updating coefficients for this iteration if Cholesky still fails
        if (iter == 1) stop("Cholesky failed at first iteration with jitter.")
        # (all_coeffs remains from previous iteration)
      } else {
        posterior_covariance_mat <- chol2inv(chol_attempt)
        posterior_mean_vec <- posterior_covariance_mat %*% rhs_for_mean
        all_coeffs <- MASS::mvrnorm(n = 1, mu = as.vector(posterior_mean_vec), Sigma = posterior_covariance_mat)
      }
    } else {
      posterior_covariance_mat <- chol2inv(chol_attempt)
      posterior_mean_vec <- posterior_covariance_mat %*% rhs_for_mean
      all_coeffs <- MASS::mvrnorm(n = 1, mu = as.vector(posterior_mean_vec), Sigma = posterior_covariance_mat)
    }
    
    # Unpack coefficients
    alpha     <- all_coeffs[1]
    beta_vec  <- all_coeffs[2:(1+p)]
    aleph     <- all_coeffs[(2+p)]
    gamma_vec <- all_coeffs[(3+p):(2+p+p)]
    
    # Store Samples after burn-in
    if (iter > burn_in) {
      current_sample_idx <- iter - burn_in
      alpha_samples[current_sample_idx] <- alpha
      beta_samples[current_sample_idx, ] <- beta_vec
      aleph_samples[current_sample_idx] <- aleph
      gamma_samples[current_sample_idx, ] <- gamma_vec
    }
  } # End MCMC loop
  
  # Return named list of samples
  list(
    alpha = alpha_samples,
    beta = beta_samples,
    aleph = aleph_samples,
    gamma = gamma_samples
  )
}
#--------------------------------------------------------------------
# --- Example: Simulate, Fit, and Analyze ---

# 1. Simulate Data
sim_p <- 3 # Number of covariates for X
sim_data_list <- simulate_simple_logistic_data(n = 3000, p = sim_p, seed = 2024)
y_sim <- sim_data_list$y_vec
X_sim <- sim_data_list$X_mat
Z_sim <- -0.5 + sim_data_list$Z_vec

true_params <- sim_data_list$true_params

cat("\n--- True Parameters from Simulation ---\n")
print(true_params)

# 2. Fit the Simplified Gibbs Sampler
# Note: R Gibbs samplers are slow. For real use, more iterations are needed.
# For this example, using fewer iterations for speed.
n_iterations_example <- 10000 # Total iterations
n_burnin_example <- 1000    # Burn-in

cat(sprintf("\n--- Running Simplified Gibbs Sampler (N=%d, p=%d, Iter=%d, Burn=%d) ---\n", 
            nrow(X_sim), ncol(X_sim), n_iterations_example, n_burnin_example))

# You can adjust prior_sd. A larger SD means a less informative prior.
# A very small SD would shrink estimates strongly towards prior_mean (0 by default).
mcmc_samples <- simple_logistic_gibbs_R(
  y_vec = y_sim,
  X_mat = X_sim,
  Z_vec = Z_sim,
  n_iter = n_iterations_example,
  burn_in = n_burnin_example,
  prior_sd = 5.0, # Adjusted prior SD
  seed = 101
)

cat("\n--- MCMC Sampling Complete ---\n")

# 3. Analyze Results
cat("\n--- Posterior Means vs True Parameters ---\n")
cat("Alpha (Intercept): True =", true_params$alpha, 
    ", Estimated =", mean(mcmc_samples$alpha), "\n")

cat("\nBeta (Main X effects):\n")
beta_comparison <- data.frame(
  True = true_params$beta,
  Estimated_Mean = colMeans(mcmc_samples$beta)
)
rownames(beta_comparison) <- paste0("Beta_X", 1:sim_p)
print(beta_comparison)

cat("\nAleph (Main Z effect / Treatment Intercept):\n True =", true_params$aleph, 
    ", Estimated =", mean(mcmc_samples$aleph), "\n")

cat("\nGamma (Z*X interaction effects):\n")
gamma_comparison <- data.frame(
  True = true_params$gamma,
  Estimated_Mean = colMeans(mcmc_samples$gamma)
)
rownames(gamma_comparison) <- paste0("Gamma_Z_X", 1:sim_p)
print(gamma_comparison)

# --- Plotting some trace plots (optional, good for checking convergence) ---
# Requires more iterations for meaningful plots usually

if (length(mcmc_samples$alpha) > 1) {
  par(mfrow = c(2, 3)) # Arrange plots
  plot(mcmc_samples$alpha, type = "l", main = "Trace: Alpha (Intercept)", ylab = "alpha")
  abline(h = true_params$alpha, col = "red", lwd = 2)
  
  if (sim_p > 0) {
    plot(mcmc_samples$beta[,1], type = "l", main = "Trace: Beta_X1", ylab = "beta1")
    abline(h = true_params$beta[1], col = "red", lwd = 2)
  }
  if (sim_p > 0) {
    plot(mcmc_samples$beta[,2], type = "l", main = "Trace: Beta_X2", ylab = "beta2")
    abline(h = true_params$beta[2], col = "red", lwd = 2)
  }
  
  plot(mcmc_samples$aleph, type = "l", main = "Trace: Aleph (Z Intercept)", ylab = "aleph")
  abline(h = true_params$aleph, col = "red", lwd = 2)
  
  if (sim_p > 0) {
    plot(mcmc_samples$gamma[,1], type = "l", main = "Trace: Gamma_Z_X1", ylab = "gamma1")
    abline(h = true_params$gamma[1], col = "red", lwd = 2)
  }
  if (sim_p > 0) {
    plot(mcmc_samples$gamma[,2], type = "l", main = "Trace: Gamma_Z_X2", ylab = "gamma2")
    abline(h = true_params$gamma[2], col = "red", lwd = 2)
  }
  par(mfrow = c(1, 1)) # Reset plot layout
}

cat("\n--- Analysis Example Complete ---\n")
cat("For better results and convergence checks, increase n_iter and burn_in significantly.\n")
