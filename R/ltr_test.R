# ==============================================================================
# TWO-STAGE SEMI-PARAMETRIC LTR: SIMULATION PIPELINE
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
library(DescTools)

# Source your internal files
source('R/simul_1.R')

num_chains <- 1
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = T, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = F, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE,   #30
  keep_every = 1, num_chains = num_chains, verbose = T, 
  sample_global_prior = "half-cauchy", unlink = T, 
  propensity_seperate = "none", gibbs = T, step_out = 0.5, 
  max_steps = 150, save_output = F, probit_outcome_model = F, 
  interaction_rule = "continuous_or_binary", standardize_cov = T, 
  simple_prior = F, save_partial_residual = F, regularize_ATE = F,
  sigma_residual = 0, hn_scale = 0, use_ncp = F, n_tijn = 1
)

# ------------------------------------------------------------------------------
# 1. SIMULATE THE DATA
# ------------------------------------------------------------------------------
cat("1. Simulating Data...\n")
scenario_n <- 500
data <- generate_data_2(scenario_n, is_te_hetero = T, is_mu_nonlinear = T, 
                        seed = 12, RCT = F, z_diff = 0.5, BCF = F,  sigma_sq = 1)

# Check signal-to-noise ratio of the outcome space
cat("Signal-to-noise check (sd(y_hat)/sd(y)):", sd(data$y_hat) / sd(data$y), "\n")

# Prepare matrices based on your exact column specifications
X_mat <- as.matrix(sapply(data[, c(1:6)], as.numeric))
Y_vec <- as.numeric(data$y)
Z_vec <- as.numeric(data$z)
pi_vec <- as.numeric(data$pi_x)

# ------------------------------------------------------------------------------
# 2. FIT THE MAGNITUDE MODEL (STAGE 1)
# ------------------------------------------------------------------------------
cat("2. Fitting BCF Linear Probit (Stage 1)...\n")
nbcf_fit <- bcf_linear_probit(
  X_train = X_mat,
  y_train = Y_vec,
  Z_train = Z_vec,
  propensity_train = pi_vec,
  num_gfr = 31, 
  num_burnin = 0, 
  num_mcmc = 500, 
  general_params = general_params_default
)

# ------------------------------------------------------------------------------
# 3. EXTRACT SCORE RESIDUALS (THE GROUND TRUTH PROXY)
# ------------------------------------------------------------------------------
cat("3. Extracting Score Residuals...\n")
# Extract the posterior mean of the prognostic function mu(x)
mu_hat <- colMeans(nbcf_fit$mu_hat_train) 

# Calculate the partial residual for the treatment function (R**)
R_star_star <- Y_vec - mu_hat

# Isolate the continuous treatment signal (tau_tilde)
# Since you used z_diff = 0.5 in generate_data_2, Z should be centered at -0.5 / 0.5.
# If Z_vec is still 0/1, we explicitly convert it to symmetric coding here.
Z_sym <- ifelse(Z_vec > 0, 0.5, -0.5) 
tau_tilde <- R_star_star / Z_sym

# ------------------------------------------------------------------------------
# 4. PASS TO LTR_SCORECARD (STAGE 2)
# ------------------------------------------------------------------------------
cat("4. Fitting Learning-to-Rank Scorecard (Stage 2)...\n")
# We pass your X_mat and tau_tilde to the updated LTR function
if (!exists("ltr_method")) {
  ltr_method <- "PG"
}
ltr_fit <- ltr_scorecard(
  tau_tilde = tau_tilde, 
  X = X_mat, 
  method = ltr_method,                                  # Configurable Pairwise Mean Squared Error or PG
  interaction_rule = general_params_default$interaction_rule, # Passed directly from your params
  unlink = general_params_default$unlink,
  mini_batch_size = 90000,                          # Active set to prevent O(N^2)
  n_iter = 2000, 
  burn_in = 2000, 
  seed = 12
)

# ------------------------------------------------------------------------------
# 5. EVALUATE RANKINGS (CORRECTED FOR INTERACTIONS)
# ------------------------------------------------------------------------------
cat("5. Evaluating Prescriptive Utility...\n")

# Reconstruct the interaction design matrix
interaction_rule <- general_params_default$interaction_rule
P_main <- ncol(X_mat)

are_continuous <- rep(0L, P_main)
for (j in 1:P_main) {
  unique_vals <- length(unique(X_mat[,j]))
  if (interaction_rule == "continuous") {
    if (unique_vals > 2) are_continuous[j] <- 1L
  } else if (interaction_rule == "continuous_or_binary") {
    if (unique_vals >= 2) are_continuous[j] <- 1L
  } else if (interaction_rule == "all") {
    are_continuous[j] <- 1L
  }
}

X_int <- matrix(nrow = nrow(X_mat), ncol = 0)

# Loop through identically to ltr_scorecard.R logic to rebuild columns
for (i in 1:(P_main - 1)) {
  for (j in (i + 1):P_main) {
    if (are_continuous[i] == 1 || are_continuous[j] == 1) {
      # Element-wise multiplication for the interaction term
      X_int <- cbind(X_int, X_mat[, i] * X_mat[, j])
    }
  }
}

X_full <- cbind(X_mat, X_int)

# Extract posterior means for ALL effects (main + interactions)
beta_hat <- colMeans(ltr_fit$beta)

# Compute the final predicted score
if (ncol(X_full) == length(beta_hat)) {
  predicted_rank_score <- X_full %*% beta_hat
  cat(sprintf("Successfully integrated %d main effects and %d interaction terms into the final scorecard.\n", P_main, ncol(X_int)))
} else {
  warning("Dimension mismatch between reconstructed design matrix and beta_hat. Using only main effects as fallback.")
  beta_main_hat <- beta_hat[1:P_main]
  predicted_rank_score <- X_mat %*% beta_main_hat
}

# ------------------------------------------------------------------------------
# CALCULATE FINAL METRICS
# ------------------------------------------------------------------------------
if ("tau" %in% colnames(data)) {
  true_tau <- data$tau
  
  # Metric A: Spearman Rank Correlation (Overall sorting accuracy)
  spearman_corr <- cor(predicted_rank_score, true_tau, method = "spearman")
  cat("Spearman Rank Correlation (Predicted vs True):", round(spearman_corr, 3), "\n")
  # Metric B: Top 10% Precision (Clinical Targeting Accuracy)
  q90_pred <- quantile(predicted_rank_score, 0.90)
  q90_true <- quantile(true_tau, 0.90)
  
  targeted_by_model <- which(predicted_rank_score >= q90_pred)
  true_top_responders <- which(true_tau >= q90_true)
  
  precision_top_10 <- length(intersect(targeted_by_model, true_top_responders)) / length(targeted_by_model)
  cat("Precision in Top 10% Targeting:", round(precision_top_10 * 100, 1), "%\n")
} else {
  cat("True CATE ('tau_vec') not found in data frame. Skipping benchmark evaluation.\n")
}