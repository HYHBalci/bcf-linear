# test_subgroup_sechidis.R
# This script evaluates the subgroup testing framework of Kostas Sechidis
# against the Bayesian Causal Forest (BCF) model using tolerance/credible intervals.

# 1. SETUP
# --------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)

# Load helper functions
source('R/simul_1.R')
source('R/forGiorgio.R')

# Helper to construct interaction pairs exactly as in one_simul.R
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

# 2. CONFIGURATION
# --------------------------------------------------------------------------
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
  standardize_cov = TRUE, simple_prior = FALSE, save_partial_residual = FALSE, 
  regularize_ATE = FALSE, sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1
)

# Simulation parameters
B <- 5 # Number of simulation iterations
scenario_n <- 500
tau_target_percentile <- 0.50 # Select tau such that roughly 50% of people have true CATE >= tau

# Store results
results_df <- data.frame(
  Iteration = integer(),
  Tau_Threshold = numeric(),
  True_Subgroup_Size = integer(),
  
  Giorgio_Est_Size = integer(),
  Giorgio_Type_I_Error = logical(),
  Giorgio_FSR = numeric(),
  Giorgio_Power = numeric(),
  
  Sechidis_Est_Size = integer(),
  Sechidis_Type_I_Error = logical(),
  Sechidis_FSR = numeric(),
  Sechidis_Power = numeric()
)

# 3. SIMULATION LOOP
# --------------------------------------------------------------------------
cat("Starting simulation loop...\n")

for (b in 1:B) {
  cat(sprintf("\n--- Iteration %d / %d ---\n", b, B))
  
  # A. Generate Data
  # We use generate_data_2 with the same settings as one_simul.R
  # using different seeds for each iteration
  data <- generate_data_2(scenario_n, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, 
                          seed = 12 + b, RCT = FALSE, z_diff = 0.5, BCF = FALSE, sigma_sq = 1)
  
  true_cate <- data$tau
  # Dynamically pick tau so approx 50% of the population forms the true subgroup S_tau
  tau <- quantile(true_cate, probs = tau_target_percentile)
  S_tau <- which(true_cate >= tau)
  
  # Ensure the target subgroup is non-empty
  if (length(S_tau) == 0) {
    cat("  Warning: True subgroup S_tau is empty. Skipping.\n")
    next
  }
  
  cat(sprintf("  True subgroup size: %d / %d (%.1f%%) for threshold tau = %.3f\n", 
              length(S_tau), scenario_n, length(S_tau)/scenario_n * 100, tau))
  
  # B. Fit BCF Model
  X_train <- as.matrix(sapply(data[, c(1:6)], as.numeric))
  y_train <- as.numeric(data$y)
  Z_train <- as.numeric(data$z)
  propensity_train <- as.numeric(data$pi_x)
  
  nbcf_fit <- bcf_linear_probit(
    X_train = X_train,
    y_train = y_train,
    Z_train = Z_train,
    propensity_train = propensity_train,
    num_gfr = 31, 
    num_burnin = 1000, 
    num_mcmc = 2000, 
    general_params = general_params_default
  )
  
  # C. Reconstruct Posterior CATE Distribution
  # The exact method used in one_simul.R
  X_info <- standardize_X_by_index(X_initial = X_train, process_data = FALSE, 
                                   interaction_rule = "continuous", cat_coding_method = "difference")
  boolean_vector <- as.logical(as.numeric(X_info$X_final_var_info$is_continuous)) + 
                    as.logical(as.numeric(X_info$X_final_var_info$is_binary))
  ipairs <- interaction_pairs(ncol(X_train), boolean_vector)
  
  alpha_samples <- as.vector(t(nbcf_fit$alpha))
  beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
  beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))
  
  sd_y <- sd(y_train)
  alpha_samples <- alpha_samples * sd_y
  beta_samples <- beta_samples * sd_y
  beta_int_samples <- beta_int_samples * sd_y
  
  tau_posterior <- matrix(rep(alpha_samples, each = scenario_n),
                          nrow = scenario_n, byrow = FALSE)
  tau_posterior <- tau_posterior + X_train %*% t(beta_samples)
  
  if (ncol(ipairs) > 0) {
    for (idx in 1:ncol(ipairs)) {
      j <- ipairs[1, idx]
      k <- ipairs[2, idx]
      interaction_term_values <- X_train[, j] * X_train[, k]
      tau_posterior <- tau_posterior + interaction_term_values %*% t(beta_int_samples[, idx, drop = FALSE])
    }
  }
  
  # tau_posterior is now (scenario_n x num_mcmc)
  
  # D. Subgroup Selection via Tolerance Intervals (Giorgio's Method)
  cat("  Computing simultaneous tolerance intervals (Giorgio's Method)...\n")
  
  # findconfglob searches for the marginal confidence level that yields global (simultaneous) 
  # coverage at a certain tolerance. Setting tol = 0.05 targets 95% global confidence.
  confmarg <- findconfglob(tau_posterior, tol = 0.05)
  cat(sprintf("    Found marginal confidence level: %.4f\n", confmarg))
  
  # Extract the intervals at this marginal confidence level
  intervals <- intfrompost1(tau_posterior, conf = confmarg)[[1]]
  lower_bounds_giorgio <- intervals[, 1]
  
  # The estimated subgroup L_alpha comprises individuals whose lower bound >= tau
  L_alpha_giorgio <- which(lower_bounds_giorgio >= tau)
  
  # E. Subgroup Selection via GLM Simultaneous Confidence Bands (Sechidis/Wan et al. Method)
  cat("  Computing simultaneous confidence bands (Sechidis GLM Method)...\n")
  # Convert posterior CATEs into point estimates (pseudo-observations)
  Y_star <- rowMeans(tau_posterior)
  df_train <- as.data.frame(X_train)
  lm_fit <- lm(Y_star ~ ., data = df_train)
  
  # Get predictions and standard errors for the fitted model
  preds <- predict(lm_fit, se.fit = TRUE)
  mu_hat <- preds$fit
  se_hat <- preds$se.fit
  
  # Simulate critical value c1 over the covariate region (using observed covariates)
  # We use the mvrnorm to draw from the asymptotic distribution of beta_hat
  Sigma_hat <- vcov(lm_fit)
  p_vars <- length(coef(lm_fit))
  Z_sim <- MASS::mvrnorm(1000, mu = rep(0, p_vars), Sigma = Sigma_hat)
  
  X_design <- model.matrix(lm_fit)
  
  # Avoid division by zero if standard error is numerically zero
  safe_se_hat <- ifelse(se_hat < 1e-8, 1e-8, se_hat)
  sim_values <- X_design %*% t(Z_sim)
  sim_max <- apply(sim_values, 2, function(col) max(col / safe_se_hat))
  
  alpha <- 0.05
  c1_hat <- quantile(sim_max, probs = 1 - alpha)
  
  # Calculate lower bounds for Sechidis method
  lower_bounds_sechidis <- mu_hat - c1_hat * se_hat
  L_alpha_sechidis <- which(lower_bounds_sechidis >= tau)
  
  cat(sprintf("  Estimated subgroup sizes - Giorgio: %d, Sechidis: %d (out of %d)\n", 
              length(L_alpha_giorgio), length(L_alpha_sechidis), scenario_n))
  
  # F. Calculate Performance Metrics
  
  # Giorgio's Method Metrics
  type_I_error_giorgio <- any(!(L_alpha_giorgio %in% S_tau))
  fsr_giorgio <- if (length(L_alpha_giorgio) > 0) mean(!(L_alpha_giorgio %in% S_tau)) else 0
  power_giorgio <- sum(L_alpha_giorgio %in% S_tau) / length(S_tau)
  
  # Sechidis GLM Method Metrics
  type_I_error_sechidis <- any(!(L_alpha_sechidis %in% S_tau))
  fsr_sechidis <- if (length(L_alpha_sechidis) > 0) mean(!(L_alpha_sechidis %in% S_tau)) else 0
  power_sechidis <- sum(L_alpha_sechidis %in% S_tau) / length(S_tau)
  
  cat(sprintf("    Giorgio Method  -> Type I Error: %s | FSR: %.3f | Power: %.3f\n", type_I_error_giorgio, fsr_giorgio, power_giorgio))
  cat(sprintf("    Sechidis Method -> Type I Error: %s | FSR: %.3f | Power: %.3f\n", type_I_error_sechidis, fsr_sechidis, power_sechidis))
  
  # Store metrics
  results_df <- rbind(results_df, data.frame(
    Iteration = b,
    Tau_Threshold = tau,
    True_Subgroup_Size = length(S_tau),
    
    Giorgio_Est_Size = length(L_alpha_giorgio),
    Giorgio_Type_I_Error = type_I_error_giorgio,
    Giorgio_FSR = fsr_giorgio,
    Giorgio_Power = power_giorgio,
    
    Sechidis_Est_Size = length(L_alpha_sechidis),
    Sechidis_Type_I_Error = type_I_error_sechidis,
    Sechidis_FSR = fsr_sechidis,
    Sechidis_Power = power_sechidis
  ))
}

# 4. SUMMARY OF SIMULATION
# --------------------------------------------------------------------------
cat("\n======================================================\n")
cat("SIMULATION SUMMARY OVER", B, "ITERATIONS\n")
cat("======================================================\n")
summary_stats <- data.frame(
  Metric = c("Mean Type I Error Rate", "Mean False Selection Rate (FSR)", "Mean Power"),
  Giorgio_Method = c(
    mean(results_df$Giorgio_Type_I_Error),
    mean(results_df$Giorgio_FSR),
    mean(results_df$Giorgio_Power)
  ),
  Sechidis_Method = c(
    mean(results_df$Sechidis_Type_I_Error),
    mean(results_df$Sechidis_FSR),
    mean(results_df$Sechidis_Power)
  )
)
print(summary_stats)
print(results_df)

# ==============================================================================
# 5. BONUS: APPLIED SETTING ON ACTG175
# ==============================================================================
cat("\n======================================================\n")
cat("APPLIED SETTING: ACTG175 SUBGROUP SELECTION\n")
cat("======================================================\n")
library(speff2trial)

# 1. Load Data
data(ACTG175, package = "speff2trial")

actg_sub <- ACTG175 %>% 
  filter(arms %in% c(1, 3)) %>%
  mutate(A = ifelse(arms == 1, 1, 0)) %>%
  filter(!is.na(cd420))

Y_actg <- actg_sub$cd420
Z_actg <- actg_sub$A

X_actg_raw <- actg_sub %>% 
  dplyr::select(age, wtkg, karnof, cd40, cd80, gender, homo, race, symptom, drugs, hemo, z30)

# Preprocessing from ACTG_complete.R
results_actg <- standardize_X_by_index(
  X_initial = X_actg_raw, 
  process_data = TRUE, 
  interaction_rule = "continuous_or_binary", 
  cat_coding_method = "difference"
)
X_train_actg <- results_actg$X_final
n_actg <- nrow(X_train_actg)

# 2. Fit Model
cat("Fitting BCF model on ACTG175...\n")
fit_actg <- bcf_linear_probit(
  X_train = X_train_actg, y_train = Y_actg, Z_train = Z_actg - 0.5,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
  general_params = general_params_default # Borrowing from the above
)

# 3. Extract CATEs via the patched predict function (assuming predict_linear_bcf_patched exists)
# Because predict_linear_bcf_patched requires more setup, we can also extract from the posterior directly
# using the same manual method as the simulation.
cat("Extracting posterior CATEs...\n")

alpha_actg <- as.vector(t(fit_actg$alpha))
beta_actg <- do.call(rbind, lapply(1:num_chains, function(chain) fit_actg$Beta[chain, , ]))

# Note: In ACTG_complete.R, interactions might be enabled. 
# We'll stick to main effects for simplicity in this applied demonstration, 
# or use the prediction function if we source ACTG_complete.R. Let's do the manual extraction:
sd_y_actg <- sd(Y_actg)
alpha_actg <- alpha_actg * sd_y_actg
beta_actg <- beta_actg * sd_y_actg

tau_post_actg <- matrix(rep(alpha_actg, each = n_actg), nrow = n_actg, byrow = FALSE)
tau_post_actg <- tau_post_actg + X_train_actg %*% t(beta_actg)

# Interaction effects
if (!is.null(fit_actg$interaction_pairs) && ncol(fit_actg$interaction_pairs) > 0) {
  beta_int_actg <- do.call(rbind, lapply(1:num_chains, function(chain) fit_actg$Beta_int[chain, , ]))
  beta_int_actg <- beta_int_actg * sd_y_actg
  ipairs_actg <- fit_actg$interaction_pairs
  for (idx in 1:ncol(ipairs_actg)) {
    j <- ipairs_actg[1, idx]
    k <- ipairs_actg[2, idx]
    int_terms <- X_train_actg[, j] * X_train_actg[, k]
    tau_post_actg <- tau_post_actg + int_terms %*% t(beta_int_actg[, idx, drop = FALSE])
  }
}

# 4. Find subgroups
cat("Computing simultaneous intervals (Giorgio's Method)...\n")
confmarg_actg <- findconfglob(tau_post_actg, tol = 0.05)
int_actg <- intfrompost1(tau_post_actg, conf = confmarg_actg)[[1]]
lb_actg_giorgio <- int_actg[, 1]

cat("Computing simultaneous confidence bands (Sechidis GLM Method)...\n")
Y_star_actg <- rowMeans(tau_post_actg)
df_train_actg <- as.data.frame(X_train_actg)
lm_fit_actg <- lm(Y_star_actg ~ ., data = df_train_actg)

preds_actg <- predict(lm_fit_actg, se.fit = TRUE)
mu_hat_actg <- preds_actg$fit
se_hat_actg <- preds_actg$se.fit

Sigma_hat_actg <- vcov(lm_fit_actg)
p_vars_actg <- length(coef(lm_fit_actg))
Z_sim_actg <- MASS::mvrnorm(1000, mu = rep(0, p_vars_actg), Sigma = Sigma_hat_actg)

X_design_actg <- model.matrix(lm_fit_actg)
safe_se_hat_actg <- ifelse(se_hat_actg < 1e-8, 1e-8, se_hat_actg)
sim_values_actg <- X_design_actg %*% t(Z_sim_actg)
sim_max_actg <- apply(sim_values_actg, 2, function(col) max(col / safe_se_hat_actg))

c1_hat_actg <- quantile(sim_max_actg, probs = 1 - 0.05)
lb_actg_sechidis <- mu_hat_actg - c1_hat_actg * se_hat_actg

# We test multiple clinically relevant thresholds
tau_values_to_test <- c(0, 10, 20, 30)

cat("\nSubgroup Selection Results (ACTG175):\n")
for (t_val in tau_values_to_test) {
  L_giorgio <- which(lb_actg_giorgio >= t_val)
  L_sechidis <- which(lb_actg_sechidis >= t_val)
  cat(sprintf("  Threshold tau = %2d :\n", t_val))
  cat(sprintf("    Giorgio Method : Selected %4d / %4d patients (%.1f%%)\n", 
              length(L_giorgio), n_actg, length(L_giorgio)/n_actg * 100))
  cat(sprintf("    Sechidis Method: Selected %4d / %4d patients (%.1f%%)\n", 
              length(L_sechidis), n_actg, length(L_sechidis)/n_actg * 100))
}

cat("\nScript completed successfully.\n")
