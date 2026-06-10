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

# Simulation parameters (Expanded Grid for Rigorous Testing)
B <- 30 # Number of simulation draws/iterations per scenario
tau_target_percentile <- 0.50 # Select tau such that roughly 50% of people have true CATE >= tau

scenario_grid <- expand.grid(
  sample_size = c(250, 500, 1000),
  sigma_sq = c(1.0, 3.0, 5.0),
  is_te_hetero = c(TRUE, FALSE),
  is_mu_nonlinear = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

# Master dataframe for all results
master_results_df <- data.frame()

# 3. SIMULATION LOOP
# --------------------------------------------------------------------------
cat("Starting simulation loop over", nrow(scenario_grid), "scenarios...\n")

for (s in 1:nrow(scenario_grid)) {
  scenario_n <- scenario_grid$sample_size[s]
  s_sigma_sq <- scenario_grid$sigma_sq[s]
  s_te_hetero <- scenario_grid$is_te_hetero[s]
  s_mu_nonlin <- scenario_grid$is_mu_nonlinear[s]
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("SCENARIO %d / %d: n = %d, sigma_sq = %.1f, TE_Hetero = %s, Mu_Nonlin = %s\n", 
              s, nrow(scenario_grid), scenario_n, s_sigma_sq, s_te_hetero, s_mu_nonlin))
  cat(sprintf("======================================================\n"))
  
  for (b in 1:B) {
    cat(sprintf("\n--- Iteration %d / %d ---\n", b, B))
    
    # A. Generate Data
    data <- generate_data_2(scenario_n, is_te_hetero = s_te_hetero, is_mu_nonlinear = s_mu_nonlin, 
                            seed = s * 1000 + b, RCT = FALSE, z_diff = 0.5, BCF = FALSE, sigma_sq = s_sigma_sq)
    
    true_cate <- data$tau
    tau <- quantile(true_cate, probs = tau_target_percentile)
    S_tau <- which(true_cate >= tau)
    
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
      num_burnin = 2000, 
      num_mcmc = 5000, 
      general_params = general_params_default
    )
    
    # C. Reconstruct Posterior CATE Distribution
    X_info <- standardize_X_by_index(X_initial = X_train, process_data = general_params_default$standardize_cov, 
                                     interaction_rule = general_params_default$interaction_rule, cat_coding_method = "difference")
  
    alpha_samples <- as.vector(t(nbcf_fit$alpha))
    beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
    
    sd_y <- sd(y_train)
    alpha_samples <- alpha_samples * sd_y
    beta_samples <- beta_samples * sd_y
    
    if (nrow(beta_samples) > ncol(beta_samples)) {
      beta_samples <- t(beta_samples)
    }
    
    tau_posterior <- matrix(rep(alpha_samples, each = scenario_n),
                            nrow = scenario_n, byrow = FALSE)
    
    X_target <- X_train
    p_beta <- nrow(beta_samples)
    p_train <- ncol(X_train)
    p_final <- ncol(X_info$X_final)
    
    if (p_beta == p_final) {
      X_target <- X_info$X_final
    } else if (p_beta == p_train + 1) {
      if (general_params_default$propensity_seperate == "none" && general_params_default$propensity_covariate != "none") {
        X_target <- cbind(X_train, propensity_train)
      } else {
        X_target <- cbind(1, X_train)
      }
    } else if (p_beta == p_train + 2) {
      X_target <- cbind(1, X_train, propensity_train)
    } else if (p_beta == p_final + 1) {
      X_target <- cbind(1, X_info$X_final)
    } else if (p_beta == p_final + 2) {
      X_target <- cbind(1, X_info$X_final, propensity_train)
    } else if (p_beta != p_train) {
      stop(paste("Cannot align X_train (ncol =", p_train, ") with beta_samples (nrow =", p_beta, ")"))
    }
    
    tau_posterior <- tau_posterior + X_target %*% beta_samples
    
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
    
    # D. Subgroup Selection via Tolerance Intervals (Giorgio's Method)
    cat("  Computing simultaneous tolerance intervals (Giorgio's Method)...\n")
    confmarg <- findconfglob(tau_posterior, tol = 0.05)
    cat(sprintf("    Found marginal confidence level: %.4f\n", confmarg))
    
    intervals <- intfrompost1(tau_posterior, conf = confmarg)[[1]]
    lower_bounds_giorgio <- intervals[, 1]
    L_alpha_giorgio <- which(lower_bounds_giorgio >= tau)
    
    # E. Subgroup Selection via Simultaneous Confidence Bands (Sechidis Method on BCF)
    cat("  Computing simultaneous confidence bands (Sechidis Method on BCF)...\n")
    mu_hat_bcf <- rowMeans(tau_posterior)
    se_hat_bcf <- apply(tau_posterior, 1, sd)
    safe_se_hat_bcf <- ifelse(se_hat_bcf < 1e-8, 1e-8, se_hat_bcf)
    
    Z_post <- (mu_hat_bcf - tau_posterior) / safe_se_hat_bcf
    M_b <- apply(Z_post, 2, max)
    
    alpha <- 0.05
    c1_bcf <- quantile(M_b, probs = 1 - alpha)
    
    lower_bounds_sechidis <- mu_hat_bcf - c1_bcf * safe_se_hat_bcf
    L_alpha_sechidis <- which(lower_bounds_sechidis >= tau)
    
    cat(sprintf("  Estimated subgroup sizes - Giorgio: %d, Sechidis: %d (out of %d)\n", 
                length(L_alpha_giorgio), length(L_alpha_sechidis), scenario_n))
    
    # F. Calculate Performance Metrics
    type_I_error_giorgio <- as.numeric(any(!(L_alpha_giorgio %in% S_tau)))
    fsr_giorgio <- if (length(L_alpha_giorgio) > 0) mean(!(L_alpha_giorgio %in% S_tau)) else 0
    power_giorgio <- sum(L_alpha_giorgio %in% S_tau) / length(S_tau)
    
    type_I_error_sechidis <- as.numeric(any(!(L_alpha_sechidis %in% S_tau)))
    fsr_sechidis <- if (length(L_alpha_sechidis) > 0) mean(!(L_alpha_sechidis %in% S_tau)) else 0
    power_sechidis <- sum(L_alpha_sechidis %in% S_tau) / length(S_tau)
    
    cat(sprintf("    Giorgio Method  -> Type I Error: %s | FSR: %.3f | Power: %.3f\n", as.logical(type_I_error_giorgio), fsr_giorgio, power_giorgio))
    cat(sprintf("    Sechidis Method -> Type I Error: %s | FSR: %.3f | Power: %.3f\n", as.logical(type_I_error_sechidis), fsr_sechidis, power_sechidis))
    
    # Store metrics
    master_results_df <- rbind(master_results_df, data.frame(
      Scenario_ID = s,
      Sample_Size = scenario_n,
      Sigma_Sq = s_sigma_sq,
      TE_Hetero = s_te_hetero,
      Mu_Nonlin = s_mu_nonlin,
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
}

# 4. SUMMARY OF SIMULATION AND EXPORT
# --------------------------------------------------------------------------
cat("\n======================================================\n")
cat("SIMULATION SUMMARY COMPLETED\n")
cat("======================================================\n")

# Aggregate results
summary_stats <- master_results_df %>%
  group_by(Scenario_ID, Sample_Size, Sigma_Sq, TE_Hetero, Mu_Nonlin) %>%
  summarize(
    Iter_Count = n(),
    Giorgio_Mean_TypeI = mean(Giorgio_Type_I_Error),
    Giorgio_Mean_FSR = mean(Giorgio_FSR),
    Giorgio_Mean_Power = mean(Giorgio_Power),
    Sechidis_Mean_TypeI = mean(Sechidis_Type_I_Error),
    Sechidis_Mean_FSR = mean(Sechidis_FSR),
    Sechidis_Mean_Power = mean(Sechidis_Power),
    .groups = "drop"
  )

print(summary_stats)

# Save to local directory
output_dir <- "run_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, "sechidis_simulation_results.csv")
write.csv(master_results_df, output_file, row.names = FALSE)
write.csv(summary_stats, file.path(output_dir, "sechidis_simulation_summary.csv"), row.names = FALSE)

cat(sprintf("\nResults successfully saved to %s\n", output_dir))

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
  num_gfr = 50, num_burnin = 2000, num_mcmc = 5000,
  general_params = general_params_default # Borrowing from the above
)

# 3. Extract CATEs via the patched predict function (assuming predict_linear_bcf_patched exists)
# Because predict_linear_bcf_patched requires more setup, we can also extract from the posterior directly
# using the same manual method as the simulation.
cat("Extracting posterior CATEs...\n")

alpha_actg <- as.vector(t(fit_actg$alpha))
beta_actg <- do.call(rbind, lapply(1:num_chains, function(chain) fit_actg$Beta[chain, , ]))

sd_y_actg <- sd(Y_actg)
alpha_actg <- alpha_actg * sd_y_actg
beta_actg <- beta_actg * sd_y_actg

  if (nrow(beta_actg) > ncol(beta_actg)) {
    beta_actg <- t(beta_actg)
  }
  
  tau_post_actg <- matrix(rep(alpha_actg, each = n_actg), nrow = n_actg, byrow = FALSE)
  
  X_target_actg <- X_train_actg
  p_beta_actg <- nrow(beta_actg)
  p_train_actg <- ncol(X_train_actg)
  
  if (p_beta_actg == p_train_actg + 1) {
    if (general_params_default$propensity_seperate == "none" && general_params_default$propensity_covariate != "none") {
      X_target_actg <- cbind(X_train_actg, propensity_actg)
    } else {
      X_target_actg <- cbind(1, X_train_actg)
    }
  } else if (p_beta_actg == p_train_actg + 2) {
    X_target_actg <- cbind(1, X_train_actg, propensity_actg)
  } else if (p_beta_actg != p_train_actg) {
    stop(paste("Cannot align X_train_actg (ncol =", p_train_actg, ") with beta_actg (nrow =", p_beta_actg, ")"))
  }
  
  tau_post_actg <- tau_post_actg + X_target_actg %*% beta_actg

# Interaction effects
if (!is.null(fit_actg$interaction_pairs) && ncol(fit_actg$interaction_pairs) > 0) {
  beta_int_actg <- do.call(rbind, lapply(1:num_chains, function(chain) fit_actg$Beta_int[chain, , ]))
  beta_int_actg <- beta_int_actg * sd_y_actg
  if (nrow(beta_int_actg) > ncol(beta_int_actg)) {
    beta_int_actg <- t(beta_int_actg)
  }
  ipairs_actg <- fit_actg$interaction_pairs
  for (idx in 1:ncol(ipairs_actg)) {
    j <- ipairs_actg[1, idx]
    k <- ipairs_actg[2, idx]
    int_terms <- X_train_actg[, j] * X_train_actg[, k]
    beta_int_actg_idx <- if(is.matrix(beta_int_actg)) beta_int_actg[idx, , drop = FALSE] else matrix(beta_int_actg[idx], nrow=1)
    tau_post_actg <- tau_post_actg + matrix(int_terms, ncol=1) %*% beta_int_actg_idx
  }
}

# 4. Find subgroups
cat("Computing simultaneous intervals (Giorgio's Method)...\n")
confmarg_actg <- findconfglob(tau_post_actg, tol = 0.05)
int_actg <- intfrompost1(tau_post_actg, conf = confmarg_actg)[[1]]
lb_actg_giorgio <- int_actg[, 1]

cat("Computing simultaneous confidence bands (Sechidis Method on BCF)...\n")
mu_hat_actg <- rowMeans(tau_post_actg)
se_hat_actg <- apply(tau_post_actg, 1, sd)
safe_se_hat_actg <- ifelse(se_hat_actg < 1e-8, 1e-8, se_hat_actg)

Z_post_actg <- (mu_hat_actg - tau_post_actg) / safe_se_hat_actg
M_b_actg <- apply(Z_post_actg, 2, max)

alpha_actg <- 0.05
c1_actg <- quantile(M_b_actg, probs = 1 - alpha_actg)
lb_actg_sechidis <- mu_hat_actg - c1_actg * safe_se_hat_actg

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
