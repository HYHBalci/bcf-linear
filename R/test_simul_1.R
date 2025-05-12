source('R/simul_1.R')
library(stochtree)
# Progress bar: range from 0 to n_simul
pb <- txtProgressBar(min = 0, max = n_simul, style = 3)

# A wrapper function that updates the progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Define model specs
n_values <- c(250, 500)
heterogeneity <- c(TRUE, FALSE)
linearity <- c(TRUE, FALSE)


compute_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper)
  interval_length <- mean(ci_upper - ci_lower)
  return(c(rmse, coverage, interval_length))
}

interaction_pairs <- function(num_covariates) {
  combn(1:num_covariates, 2)
}

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = T, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = i, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = T
)
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 7, RCT = T, test = F)
weights <- rep(1,n_obser)

nbcf_fit <- bcf_linear(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 5000, 
  general_params = general_params_default
)
# Compute means of x for correct alpha rescalin
X <-  as.matrix(sapply(data[, c(1:6)], as.numeric))
means_x <- colMeans(X)

# Compute means of interaction terms jointly
interaction_means <- apply(interaction_pairs(ncol(X)), 2, function(idx) {
  mean(X[, idx[1]] * X[, idx[2]])  # Joint mean of interaction term X_j * X_k
})

# Rescale coefficients
alpha_samples <- as.vector(t(nbcf_fit$alpha))
beta_samples  <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta[chain, , ]))
beta_int_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta_int[chain, , ]))

# Properly rescale slopes
beta_samples  <- beta_samples  * sd(data$y)
beta_int_samples <- beta_int_samples * sd(data$y)

# Properly rescale intercept with joint interaction means
alpha_samples <- alpha_samples * sd(data$y)
X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
# Vectorized approach
ipairs <- interaction_pairs(6)
tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                        nrow = n_obser, byrow = FALSE) +
  X %*% t(beta_samples)

for (idx in 1:ncol(ipairs)) {
  j <- ipairs[1, idx]
  k <- ipairs[2, idx]
  tau_posterior <- tau_posterior +
    (X[, j] * X[, k]) %*% t(beta_int_samples[, idx])
}

mean(tau_posterior)
median(tau_posterior)