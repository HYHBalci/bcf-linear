# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
source('R/simul_1.R')

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
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, sample_sigma2_global = TRUE, sigma2_global_init = NULL, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001, variable_weights = NULL, 
  propensity_covariate = "none", adaptive_coding = FALSE, rfx_prior_var = NULL, random_seed = -1, 
  keep_burnin = FALSE, keep_gfr = FALSE, keep_every = 1, num_chains = 1, verbose = TRUE, 
  probit_outcome_model = FALSE, standardize_cov = FALSE, interaction_rule = "continuous"
)
scenario_n <- 750
data <- generate_data_2(scenario_n, is_te_hetero = F, is_mu_nonlinear = T, seed = 31, RCT = T, z_diff = 0.5, BCF = F,  sigma_sq = 1)
nbcf_fit <- bcf_restricted_score_test(X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
                y_train = as.numeric(data$y),
                Z_train = as.numeric(data$z), 
                propensity_train = rep(0, length(data$z)),
                num_gfr = 25, 
                num_burnin = 500, 
                num_mcmc = 2000,
                general_params = general_params_default)
