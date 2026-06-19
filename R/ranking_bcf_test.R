# ==============================================================================
# RANKING BCF: SIMULATION PIPELINE
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
library(DescTools)

# Source your internal files
source('R/simul_1.R')

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

# Convert treatment from 0/1 to symmetric coding if required, but ranking_bcf expects 0/1 generally 
# or continuous Z. Z_vec here is binary.

# ------------------------------------------------------------------------------
# 2. FIT THE RANKING BCF MODEL
# ------------------------------------------------------------------------------
cat("2. Fitting Ranking BCF...\n")

# Choose kl_divergence_method: "listnet" or "distributional"
# We will use "listnet" for this test.
fit_ranking <- ranking_bcf(
  X_train = X_mat,
  y_train = Y_vec,
  Z_train = Z_vec,
  propensity_train = pi_vec,
  kl_divergence_method = "listnet",
  num_burnin = 100, 
  num_mcmc = 200, 
  prognostic_forest_params = list(num_trees = 50),
  treatment_effect_forest_params = list(num_trees = 50)
)

# ------------------------------------------------------------------------------
# 3. EVALUATE RANKINGS
# ------------------------------------------------------------------------------
cat("3. Evaluating Prescriptive Utility...\n")

# Extract posterior mean of the treatment effect (tau)
tau_hat_matrix <- fit_ranking$tau_hat_train
tau_hat_mean <- colMeans(tau_hat_matrix)

predicted_rank_score <- tau_hat_mean

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
  
  if(length(targeted_by_model) > 0) {
    precision_top_10 <- length(intersect(targeted_by_model, true_top_responders)) / length(targeted_by_model)
    cat("Precision in Top 10% Targeting:", round(precision_top_10 * 100, 1), "%\n")
  } else {
    cat("Precision in Top 10% Targeting: N/A (no targets identified)\n")
  }
} else {
  cat("True CATE ('tau_vec') not found in data frame. Skipping benchmark evaluation.\n")
}
