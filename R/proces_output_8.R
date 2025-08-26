# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(stochtree)
source("R/simul_1.R")

# --- Helper Functions (as provided in your script) ---
compute_mode <- function(x) {
  if (any(is.na(x))) return(NA)
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower = NULL, ci_upper = NULL) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  if (is.null(ci_lower) || is.null(ci_upper)) {
    return(c(rmse = rmse))
  }
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

# --- Simulation specifications ---
n_simul <- 50
n_values <- c(250, 500)
heterogeneity_opts <- c(TRUE, FALSE)
linearity_opts <- c(TRUE, FALSE)

# INITIALIZE the summary table (with new columns for predictive performance)
results <- expand.grid(n = n_values,
                       heterogeneity = heterogeneity_opts,
                       non_linearity = !linearity_opts) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA,
         # New columns for predictive RMSE
         rmse_y_obs = NA, rmse_y0 = NA, rmse_y1 = NA)

# INITIALIZE a list to store the new detailed data for plotting
all_scenario_results_list <- list()


# 2. SIMULATION LOOPS
# --------------------------------------------------------------------------
# Loop over all model settings
for (n_obser in n_values) {
  for (het in heterogeneity_opts) {
    for (lin in linearity_opts) {
      
      # This matrix will temporarily hold results for the current scenario (now with 9 columns)
      result_matrix <- matrix(NA, nrow = n_simul, ncol = 9)
      colnames(result_matrix) <- c("ATE_RMSE", "ATE_Coverage", "ATE_Length",
                                   "CATE_RMSE", "CATE_Coverage", "CATE_Length",
                                   "Y_obs_RMSE", "Y0_RMSE", "Y1_RMSE")
      
      cat(paste("\nRunning Scenario: n =", n_obser, "| Heterogeneity =", het, "| Non-Linearity =", !lin, "\n"))
      
      for (i in 1:n_simul) {
        file_name <- sprintf(
          "D:/BCF_STANDARD_2/BCF_STANDARD_2_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
          ifelse(het, "T", "F"), ifelse(lin, "T", "F"), n_obser, i
        )
        
        tryCatch({
          if (!file.exists(file_name)) {
            stop(paste("File not found:", file_name))
          } else {
            load(file_name) # loads object: nbcf_fit
          }
          
          # --- Generate the true data for this simulation run ---
          data <- generate_data_2(n_obser, het, lin, seed = i, RCT = FALSE)
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          
          # --- ATE and CATE Metrics (Original Logic) ---
          tau_posterior <- nbcf_fit$tau_hat_train
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
          ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
          
          ate_draws <- colMeans(tau_posterior)
          est_ate <- mean(ate_draws)
          ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
          
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper)
          
          # --- Predictive Performance Metrics (New Logic) ---
          mu_posterior <- nbcf_fit$mu_hat_train
          
          # Construct posteriors for potential outcomes Y(0) and Y(1)
          y0_posterior <- mu_posterior
          y1_posterior <- mu_posterior + tau_posterior
          
          # Get point estimates (posterior mean)
          est_y0 <- rowMeans(y0_posterior)
          est_y1 <- rowMeans(y1_posterior)
          
          # Construct point estimate for the observed outcome Y
          est_y_obs <- est_y0 * (1 - data$z) + est_y1 * data$z
          
          # Get true values
          true_y0 <- data$mu
          true_y1 <- data$mu + data$tau
          
          # Compute predictive RMSE
          rmse_y_obs <- sqrt(mean((data$y - est_y_obs)^2))
          rmse_y0 <- sqrt(mean((true_y0 - est_y0)^2))
          rmse_y1 <- sqrt(mean((true_y1 - est_y1)^2))
          
          # Store all results for this simulation run
          result_matrix[i, ] <- c(ate_vec, cate_vec, rmse_y_obs, rmse_y0, rmse_y1)
          
        }, error = function(e) {
          # Error message
          cat("Error processing file:", file_name, "-", conditionMessage(e), "\n")
        })
      }
      
      # --- UPDATE THE SUMMARY TABLE ---
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$non_linearity == !lin)
      results[idx, 4:6] <- colMeans(result_matrix[, 1:3, drop = FALSE], na.rm = TRUE) # ATE
      results[idx, 7:9] <- colMeans(result_matrix[, 4:6, drop = FALSE], na.rm = TRUE) # CATE
      results[idx, 10:12] <- colMeans(result_matrix[, 7:9, drop = FALSE], na.rm = TRUE) # Predictive
      
      # --- STORE DETAILED RESULTS FOR PLOTTING ---
      scenario_df <- as.data.frame(result_matrix) %>%
        mutate(
          n = n_obser,
          heterogeneity = ifelse(het, "Heterogeneous", "Homogeneous"),
          linearity = ifelse(lin, "Linear", "Non-Linear")
        )
      all_scenario_results_list[[length(all_scenario_results_list) + 1]] <- scenario_df
    }
  }
}

# 3. FINAL OUTPUTS
# --------------------------------------------------------------------------
save(results, file = 'results_BCF_standard_2_with_preds.RData')

# --- 3A. Print the updated summary results table ---
cat("\n\n--- Summary Results Table ---\n")
print(results, digits = 3)


# --- 3B. Create and print the boxplots from detailed results ---
cat("\n\n--- Generating Boxplots ---\n")

# Combine the list of data frames into one
all_results_df <- bind_rows(all_scenario_results_list)

# Pivot data to a long format for easier plotting
results_long <- all_results_df %>%
  na.omit() %>%
  pivot_longer(
    cols = ends_with("RMSE"), # Focus on RMSE for clarity
    names_to = "estimand",
    names_pattern = "(.*)_RMSE",
    values_to = "rmse_value"
  ) %>%
  mutate(
    estimand = factor(estimand, levels = c("ATE", "CATE", "Y_obs", "Y0", "Y1"),
                      labels = c("ATE", "CATE (τ)", "Observed (Y)", "Y(0)", "Y(1)"))
  )

# Create the plot object
metrics_plot <- ggplot(results_long, aes(x = heterogeneity, y = rmse_value, fill = linearity)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  facet_grid(estimand ~ n, scales = "free_y") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribution of RMSE Across Simulation Scenarios",
    x = "Treatment Effect Type",
    y = "Root Mean Squared Error (RMSE)",
    fill = "Functional Form"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Explicitly print the plot
print(metrics_plot)