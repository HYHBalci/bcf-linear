# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
source('R/simul_1.R')
# library(stochtree) # Uncomment if standardize_X_by_index or other functions require it

# --- Source necessary files ---
# Ensure this path points to the file containing generate_data_2 and standardize_X_by_index
# source("R/simul_1.R")

# --- Helper Functions (copied from original script) ---
compute_mode <- function(x) {
  # Using mean as a proxy for mode per original code.
  return(mean(x))
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(data.frame(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

interaction_pairs <- function(num_covariates, boolean_vector) {
  interaction_list <- list()
  if (num_covariates > 1) {
    for (j in 1:(num_covariates - 1)) {
      for (k in (j + 1):num_covariates) {
        if (boolean_vector[j] || boolean_vector[k])
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
      }
    }
  }
  if (length(interaction_list) == 0) return(matrix(nrow = 2, ncol = 0))
  return(do.call(cbind, interaction_list))
}

# 2. CONFIGURATION FOR SINGLE RUN ANALYSIS
# --------------------------------------------------------------------------



# 3. CORE EVALUATION LOGIC
# --------------------------------------------------------------------------
scenario_n <- 500
data <- generate_data_2(scenario_n, is_te_hetero = F, is_mu_nonlinear = T, seed = 31, RCT = F, z_diff = 0.5, BCF = F,  sigma_sq =1)

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = NULL, 
  sigma2_global_shape = 0, sigma2_global_scale = 0, 
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = F, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = -1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = TRUE, 
  probit_outcome_model = FALSE
)
nbcf_fit <- bcf(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 0, 
  num_mcmc = 500,
  general_params = general_params_default
)
true_cate <- data$tau
true_ate <- mean(true_cate)
tau_posterior <- nbcf_fit$tau_hat_train

# --- Compute Point Estimates and Intervals ---
cat("Calculating point estimates and credible intervals...\n")
# CATE estimates
tau_estimate <- apply(tau_posterior, 1, compute_mode)
ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)

# ATE estimates
ate_draws <- colMeans(tau_posterior)
est_ate <- mean(ate_draws)
ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)

# --- Calculate Final Metrics ---
ate_metrics <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
cate_metrics <- compute_metrics(true_cate, tau_estimate, ci_tau_lower, ci_tau_upper)

# --- Display Results ---
cat("\n--- Evaluation Results ---\n")
cat("ATE Metrics:\n")
print(ate_metrics)
cat("\nCATE Metrics:\n")
print(cate_metrics)



# 4. VISUALIZATION FOR SINGLE RUN
# --------------------------------------------------------------------------
# Create a data frame for plotting true vs. estimated CATE values
plot_data <- data.frame(
  true_cate = true_cate,
  estimated_cate = tau_estimate,
  ci_lower = ci_tau_lower,
  ci_upper = ci_tau_upper
)

# Scatter plot of True CATE vs. Estimated CATE
cate_plot <- ggplot(plot_data, aes(x = true_cate, y = estimated_cate)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  labs(
    title = "True vs. Estimated CATE",
    subtitle = paste("File:"),
    x = "True Treatment Effect",
    y = "Estimated Treatment Effect (Posterior Mean)"
  ) +
  theme_minimal()

print(cate_plot)