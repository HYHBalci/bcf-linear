# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(stochtree)

# Ensure the file containing your data generation function is available
source("R/simul_1.R")

# --- Helper Functions ---
compute_mode <- function(x) {
  return(mean(x))
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

interaction_pairs <- function(num_covariates, boolean_vector) {
  interaction_list <- list()
  if (num_covariates > 1) {
    for (j in 1:(num_covariates - 1)) {
      for (k in (j + 1):num_covariates) {
        # This rule includes an interaction if at least one of the pair is continuous or binary
        if (boolean_vector[j] || boolean_vector[k]) {
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
        }
      }
    }
  }
  if (length(interaction_list) == 0) return(matrix(nrow = 2, ncol = 0))
  return(do.call(cbind, interaction_list))
}


# 2. SIMULATION SPECIFICATIONS
# --------------------------------------------------------------------------
n_simul <- 50
n_values <- c(250, 500)
heterogeneity_opts <- c(TRUE, FALSE)
linearity_opts <- c(TRUE, FALSE)
num_chains <- 2

# --- Pre-calculate EXPECTED dimensions ---
# This creates a single source of truth for the expected number of parameters.
# Any loaded file that does not conform to this will be skipped.
X_template <- as.matrix(generate_data_2(500, T, T, seed = 1)[, 1:6])
res <- standardize_X_by_index(X_initial = X_template, process_data = F, interaction_rule = "continuous_or_binary", cat_coding_method = "difference")
boolean <- as.logical(as.numeric(res$X_final_var_info$is_binary) + as.numeric(res$X_final_var_info$is_continuous))

expected_num_betas <- ncol(X_template)
ipairs <- interaction_pairs(ncol(X_template), boolean)
expected_num_beta_ints <- ncol(ipairs)


# 3. INITIALIZE RESULTS STORAGE
# --------------------------------------------------------------------------
# Summary table
results <- expand.grid(
  n = n_values,
  heterogeneity = heterogeneity_opts,
  linearity = linearity_opts
) %>%
  mutate(
    rmse_ate = NA, cover_ate = NA, len_ate = NA,
    rmse_cate = NA, cover_cate = NA, len_cate = NA,
    alpha_coverage = NA, param_tpr = NA, param_fpr = NA,
    param_fdr = NA, param_coverage = NA
  )

# List for detailed results for plotting
all_scenario_results_list <- list()


# 4. SIMULATION LOOPS
# --------------------------------------------------------------------------
for (n_obser in n_values) {
  for (het in heterogeneity_opts) {
    for (lin in linearity_opts) {
      
      result_matrix <- matrix(NA, nrow = n_simul, ncol = 11)
      colnames(result_matrix) <- c(
        "ATE_RMSE", "ATE_Coverage", "ATE_Length",
        "CATE_RMSE", "CATE_Coverage", "CATE_Length",
        "Alpha_Coverage", "Param_TPR", "Param_FPR",
        "Param_FDR", "Param_Coverage"
      )
      
      cat(paste("\nRunning Scenario: n =", n_obser, "| Heterogeneity =", het, "| Linearity =", lin, "\n"))
      
      for (i in 1:n_simul) {
        file_name <- sprintf(
          "D:/new ok 2025 horseshoe/Block_horse_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
          ifelse(het, "T", "F"), ifelse(lin, "T", "F"), n_obser, i
        )
        
        
          if (!file.exists(file_name)) {
            cat("File not found, skipping:", file_name, "\n")
            next
          } else {
            load(file_name) # loads object: nbcf_fit
          }
          
          # --- !! KEY ROBUSTNESS CHECK !! ---
          # Verify that the dimensions of the loaded object match expectations
          if (dim(nbcf_fit$Beta)[3] != expected_num_betas || dim(nbcf_fit$Beta_int)[3] != expected_num_beta_ints) {
            cat("!! WARNING: Dimension mismatch in", basename(file_name), "- SKIPPING.\n")
            cat("  Expected:", expected_num_betas, "main,", expected_num_beta_ints, "interaction effects.\n")
            cat("  Found:   ", dim(nbcf_fit$Beta)[3], "main,", dim(nbcf_fit$Beta_int)[3], "interaction effects.\n")
            next # Skip this iteration
          }
          
          # --- ATE/CATE Processing ---
          data <- generate_data_2(n_obser, het, lin, seed = i, RCT = FALSE, z_diff = 0.5)
          X <- as.matrix(data[, 1:6])
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          
          alpha_samples <- as.vector(t(nbcf_fit$alpha))
          beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
          beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))
          
          sd_y <- sd(data$y)
          alpha_samples <- alpha_samples * sd_y
          beta_samples <- beta_samples * sd_y
          beta_int_samples <- beta_int_samples * sd_y
          
          tau_posterior <- matrix(rep(alpha_samples, each = n_obser), nrow = n_obser, byrow = FALSE) + X %*% t(beta_samples)
          
          if (ncol(ipairs) > 0) {
            for (idx in 1:ncol(ipairs)) {
              j <- ipairs[1, idx]
              k <- ipairs[2, idx]
              tau_posterior <- tau_posterior + (X[, j] * X[, k]) %*% t(beta_int_samples[, idx, drop = FALSE])
            }
          }
          
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          ci_tau <- apply(tau_posterior, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
          
          ate_draws <- colMeans(tau_posterior)
          est_ate <- mean(ate_draws)
          ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
          
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau[1, ], ci_tau[2, ])
          
          # --- Parameter-Level Analysis ---
          if (het) { # Heterogeneous case: tau = 1 + 2*x2*x4
            true_alpha <- 1
            true_beta <- rep(0, expected_num_betas)
            true_beta_int <- rep(0, expected_num_beta_ints)
            interaction_idx <- which(ipairs[1, ] == 2 & ipairs[2, ] == 4)
            if (length(interaction_idx) > 0) true_beta_int[interaction_idx] <- 2
          } else { # Homogeneous case: tau = 3
            true_alpha <- 3
            true_beta <- rep(0, expected_num_betas)
            true_beta_int <- rep(0, expected_num_beta_ints)
          }
          
          ci_alpha <- quantile(alpha_samples, probs = c(0.025, 0.975))
          ci_beta <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))
          ci_beta_int <- apply(beta_int_samples, 2, quantile, probs = c(0.025, 0.975))
          
          alpha_coverage <- as.numeric(true_alpha >= ci_alpha[1] & true_alpha <= ci_alpha[2])
          all_params_true <- c(true_beta, true_beta_int)
          all_params_ci_lower <- c(ci_beta[1, ], ci_beta_int[1, ])
          all_params_ci_upper <- c(ci_beta[2, ], ci_beta_int[2, ])
          
          param_coverage <- mean(all_params_true >= all_params_ci_lower & all_params_true <= all_params_ci_upper, na.rm = TRUE)
          
          is_truly_nonzero <- all_params_true != 0
          is_selected <- !(all_params_ci_lower < 0 & all_params_ci_upper > 0)
          
          TP <- sum(is_selected & is_truly_nonzero, na.rm = TRUE)
          FP <- sum(is_selected & !is_truly_nonzero, na.rm = TRUE)
          P <- sum(is_truly_nonzero, na.rm = TRUE)
          N <- sum(!is_truly_nonzero, na.rm = TRUE)
          
          TPR <- if (P > 0) TP / P else 0
          FPR <- if (N > 0) FP / N else 0
          FDR <- if ((TP + FP) > 0) FP / (TP + FP) else 0
          
          param_metrics <- c(alpha_coverage, TPR, FPR, FDR, param_coverage)
          result_matrix[i, ] <- c(ate_vec, cate_vec, param_metrics)
          
        
      }
      
      # --- Update Summary and Detailed Results ---
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      if(length(idx) == 1) {
        results[idx, 4:14] <- colMeans(result_matrix, na.rm = TRUE)
      }
      
      scenario_df <- as.data.frame(result_matrix) %>%
        mutate(
          n = n_obser,
          heterogeneity = ifelse(het, "Heterogeneous", "Homogeneous"),
          linearity = ifelse(lin, "Linear", "Nonlinear")
        )
      all_scenario_results_list[[length(all_scenario_results_list) + 1]] <- scenario_df
    }
  }
}

# 5. FINAL OUTPUTS
# --------------------------------------------------------------------------
save(results, file = 'results_linked_with_params.RData')
cat("\n\n--- Summary Results Table ---\n")
print(results)

# Combine all detailed results into a single data frame
all_results_df <- bind_rows(all_scenario_results_list)

# --- Boxplots for ATE/CATE Metrics ---
if (nrow(all_results_df) > 0) {
  cat("\n\n--- Generating Boxplots for ATE/CATE Metrics ---\n")
  results_long_effects <- all_results_df %>%
    na.omit() %>%
    pivot_longer(
      cols = c(ends_with("RMSE"), ends_with("Length")),
      names_to = c("effect_type", "metric"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    mutate(metric = factor(metric, levels = c("RMSE", "Length")))
  
  metrics_plot <- ggplot(results_long_effects, aes(x = heterogeneity, y = value, fill = linearity)) +
    geom_boxplot(alpha = 0.8) +
    facet_grid(metric ~ n, scales = "free_y") +
    scale_fill_brewer(palette = "Pastel1") +
    labs(
      title = "Distribution of ATE/CATE Performance Metrics",
      x = "Treatment Effect Type", y = "Metric Value", fill = "Functional Form"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  print(metrics_plot)
  
  # --- Boxplots for Parameter-Level Metrics ---
  cat("\n\n--- Generating Boxplots for Parameter-Level Metrics ---\n")
  results_long_params <- all_results_df %>%
    na.omit() %>%
    pivot_longer(
      cols = c("Param_TPR", "Param_FPR", "Param_FDR", "Param_Coverage", "Alpha_Coverage"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = factor(metric, levels = c("Param_TPR", "Param_FPR", "Param_FDR", "Param_Coverage", "Alpha_Coverage")))
  
  param_plot <- ggplot(results_long_params, aes(x = as.factor(n), y = value, fill = linearity)) +
    geom_boxplot(alpha = 0.8) +
    facet_grid(metric ~ heterogeneity, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Distribution of Parameter-Level Performance Metrics",
      x = "Sample Size (n)", y = "Metric Value", fill = "Functional Form"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    geom_hline(data = . %>% filter(metric %in% c("Param_Coverage", "Alpha_Coverage")),
               aes(yintercept = 0.95), linetype = "dashed", color = "red")
  print(param_plot)
} else {
  cat("\n\nNo valid simulation results to plot.\n")
}