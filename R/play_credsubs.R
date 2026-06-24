# --------------------------------------------------------------------------
# Play script: Combining bcf_linear_probit with the credsubs package
# This script generates data, fits the semi-parametric BCF model,
# reconstructs the posterior CATE distribution, and applies the 
# Bayesian Credible Subgroups (Schnell et al.) methodology.
# --------------------------------------------------------------------------

# 1. SETUP
library(dplyr)
library(ggplot2)
library(stochtree) # Assuming bcf_linear_probit is available here

# Install credsubs if you haven't already:
# install.packages("credsubs")
library(credsubs)

source('R/simul_1.R')

# --------------------------------------------------------------------------
# 2. GENERATE DATA
# --------------------------------------------------------------------------
cat("Generating data...\n")
scenario_n <- 500 # Using a smaller sample for a quick play script
data <- generate_data_2(n = scenario_n, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, 
                        seed = 1234, RCT = FALSE, z_diff = 0.5, BCF = FALSE, sigma_sq = 1)

X_train <- as.matrix(sapply(data[, c(1:6)], as.numeric))
y_train <- as.numeric(data$y)
Z_train <- as.numeric(data$z)
pi_hat <- as.numeric(data$pi_x)

# --------------------------------------------------------------------------
# 3. FIT THE BCF-LINEAR-PROBIT MODEL
# --------------------------------------------------------------------------
cat("Fitting the bcf_linear_probit model...\n")
num_chains <- 1
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE,
  keep_every = 1, num_chains = num_chains, verbose = FALSE, 
  sample_global_prior = "half-cauchy", unlink = TRUE, propensity_seperate = "none", 
  gibbs = TRUE, step_out = 0.5, max_steps = 150, save_output = FALSE, 
  probit_outcome_model = FALSE, interaction_rule = "continuous_or_binary",
  standardize_cov = FALSE, simple_prior = TRUE, save_partial_residual = FALSE, 
  regularize_ATE = FALSE, sigma_residual = 0, hn_scale = 0, use_ncf = FALSE, n_tijn = 1
)

nbcf_fit <- bcf_linear_probit(
  X_train = X_train,
  y_train = y_train,
  Z_train = Z_train,
  propensity_train = pi_hat,
  num_gfr = 10, 
  num_burnin = 500, 
  num_mcmc = 1000, 
  general_params = general_params_default
)

# --------------------------------------------------------------------------
# 4. RECONSTRUCT POSTERIOR SAMPLES OF CATE (tau)
# --------------------------------------------------------------------------
cat("Reconstructing posterior samples of CATE...\n")
X_info <- standardize_X_by_index(X_initial = X_train, 
                                 process_data = general_params_default$standardize_cov, 
                                 interaction_rule = general_params_default$interaction_rule, 
                                 cat_coding_method = "difference")

alpha_samples <- as.vector(t(nbcf_fit$alpha))
beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))

sd_y <- sd(data$y)
alpha_samples <- alpha_samples * sd_y
beta_samples <- beta_samples * sd_y
if (nrow(beta_samples) > ncol(beta_samples)) {
  beta_samples <- t(beta_samples)
}

# Assume p_beta aligns directly with X_info$X_final for this script (from one_simul.R logic)
X_target <- cbind(1, X_info$X_final) 

tau_posterior <- matrix(rep(alpha_samples, each = scenario_n), nrow = scenario_n, byrow = FALSE)
# Add main effects (Note: excluding intercept from X_target if alpha is added separately)
tau_posterior <- tau_posterior + X_info$X_final %*% beta_samples

# Add interactions
if (!is.null(nbcf_fit$interaction_pairs) && ncol(nbcf_fit$interaction_pairs) > 0) {
  beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))
  beta_int_samples <- beta_int_samples * sd_y
  if (nrow(beta_int_samples) > ncol(beta_int_samples)) {
    beta_int_samples <- t(beta_int_samples)
  }
  
  int_pairs <- nbcf_fit$interaction_pairs
  num_interactions <- ncol(int_pairs)
  X_int <- matrix(0, nrow = nrow(X_info$X_final), ncol = num_interactions)
  for (k in 1:num_interactions) {
    X_int[, k] <- X_info$X_final[, int_pairs[1, k]] * X_info$X_final[, int_pairs[2, k]]
  }
  tau_posterior <- tau_posterior + (X_int %*% beta_int_samples)
}

# tau_posterior is now [n_obs x n_mcmc]. 
# credsubs expects [n_mcmc x n_obs] when providing the evaluated surface directly.
tau_mcmc <- t(tau_posterior) 

# --------------------------------------------------------------------------
# 5. APPLY CREDIBLE SUBGROUPS METHODOLOGY (credsubs package)
# --------------------------------------------------------------------------
cat("Applying Bayesian Credible Subgroups (credsubs)...\n")

# We want to find the subgroups where the treatment effect is strictly greater than a threshold.
# Let's say we want to identify patients where the CATE > 1
threshold_val <- 1.0

# --------------------------------------------------------------------------
# A. Pure Bayes (PB) Method via Surface Evaluation
# --------------------------------------------------------------------------
# Providing the evaluated surface directly implies PB is the most straightforward application.
cred_pb <- credsubs(params = tau_mcmc, threshold = threshold_val, cred.level = 0.80)

# --------------------------------------------------------------------------
# B. HPD and RCS Methods via Parametric Input
# --------------------------------------------------------------------------
cat("Applying HPD and RCS methods...\n")
# To use HPD and RCS, credsubs requires the raw posterior parameter draws and the design matrix.

# Construct the full parameter matrix [n_mcmc x (1 + p_main + p_int)]
if (exists("beta_int_samples")) {
  params_matrix <- cbind(matrix(alpha_samples, ncol = 1), t(beta_samples), t(beta_int_samples))
  design_matrix <- cbind(1, X_info$X_final, X_int)
} else {
  params_matrix <- cbind(matrix(alpha_samples, ncol = 1), t(beta_samples))
  design_matrix <- cbind(1, X_info$X_final)
}

# HPD guarantees coverage over the entire infinite covariate space
cred_hpd <- credsubs(params = params_matrix, design = design_matrix, threshold = threshold_val, method = "HPD", cred.level = 0.80)

# RCS restricts the coverage guarantee to the specific points provided in the design matrix
cred_rcs <- credsubs(params = params_matrix, design = design_matrix, threshold = threshold_val, method = "RCS", cred.level = 0.80)

# --------------------------------------------------------------------------
# 6. CLASSIFY AND COMPARE RESULTS
# --------------------------------------------------------------------------
# Add the boolean classification vectors to the data
data$in_D_pb <- cred_pb$exclusive
data$in_D_hpd <- cred_hpd$exclusive
data$in_D_rcs <- cred_rcs$exclusive

data <- data %>%
  mutate(
    class_pb = case_when(in_D_pb ~ "1. Benefit (D)", !cred_pb$inclusive ~ "3. No Benefit (S_c)", TRUE ~ "2. Uncertain"),
    class_hpd = case_when(in_D_hpd ~ "1. Benefit (D)", !cred_hpd$inclusive ~ "3. No Benefit (S_c)", TRUE ~ "2. Uncertain"),
    class_rcs = case_when(in_D_rcs ~ "1. Benefit (D)", !cred_rcs$inclusive ~ "3. No Benefit (S_c)", TRUE ~ "2. Uncertain")
  )

# Print Summary
cat("\n--- Credible Subgroups (Benefit 'D') Summary (Threshold =", threshold_val, ") ---\n")
cat("Pure Bayes (PB)   patients in D:", sum(data$in_D_pb), "\n")
cat("Restricted (RCS)  patients in D:", sum(data$in_D_rcs), "\n")
cat("Highest PD (HPD)  patients in D:", sum(data$in_D_hpd), "\n")

# --------------------------------------------------------------------------
# 7. VISUALIZATION (Using Pure Bayes Classification as Example)
# --------------------------------------------------------------------------
# Plot True CATE vs. Estimated CATE, colored by Credible Subgroup classification
data$est_cate <- rowMeans(tau_posterior)

p1 <- ggplot(data, aes(x = tau, y = est_cate, color = class_pb)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = threshold_val, linetype = "dotted", color = "red", size=1) +
  geom_vline(xintercept = threshold_val, linetype = "dotted", color = "red", size=1) +
  scale_color_manual(values = c("1. Benefit (D)" = "darkgreen", 
                                "2. Uncertain" = "gold3", 
                                "3. No Benefit (S_c)" = "darkred")) +
  labs(
    title = "CATE Estimates by PB Credible Subgroup",
    subtitle = paste("Threshold (delta) =", threshold_val, "| Credible Level (1-alpha) = 80%"),
    x = "True CATE",
    y = "Estimated Posterior Mean CATE",
    color = "Subgroup"
  ) +
  theme_minimal()

# Save the plot
ggsave("R/credible_subgroups_plot.png", p1, width = 8, height = 6)
cat("\nPlot saved to R/credible_subgroups_plot.png\n")
