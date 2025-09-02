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

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper, file_name) {
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
        if (boolean_vector[j] || boolean_vector[k])
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
      }
    }
  }
  if (length(interaction_list) == 0) return(matrix(nrow=2, ncol=0))
  return(do.call(cbind, interaction_list))
}

# --- Simulation specifications ---
n_simul <- 50
n_values <- c(250, 500)
heterogeneity_opts <- c(TRUE, FALSE)
linearity_opts <- c(TRUE,FALSE)

data <- generate_data_2(500, T, T, seed = 1848, RCT = FALSE)

X <- as.matrix(data[, 1:6])
res <- standardize_X_by_index(X_initial = X, process_data = F, interaction_rule = "continuous_or_binary", cat_coding_method = "difference")

boolean <- as.logical(as.numeric(res$X_final_var_info$is_binary) + as.numeric(res$X_final_var_info$is_continuous))


# INITIALIZE the summary table (as in your original code)
results <- expand.grid(n = n_values,
                       heterogeneity = heterogeneity_opts,
                       linearity = linearity_opts) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

# INITIALIZE a list to store the new detailed data for plotting
all_scenario_results_list <- list()


# 2. SIMULATION LOOPS
# --------------------------------------------------------------------------
# Loop over all model settings
for (n_obser in n_values) {
  for (het in heterogeneity_opts) {
    for (lin in linearity_opts) {
      
      # This matrix will temporarily hold results for the current scenario
      result_matrix <- matrix(NA, nrow = n_simul, ncol = 6)
      colnames(result_matrix) <- c("ATE_RMSE", "ATE_Coverage", "ATE_Length",
                                   "CATE_RMSE", "CATE_Coverage", "CATE_Length")
      
      cat(paste("\nRunning Scenario: n =", n_obser, "| Heterogeneity =", het, "| Linearity =", lin, "\n"))
      
      for (i in 1:n_simul) {
        file_name <- sprintf(
          "D:/block_horseshoe/Block_horse_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
          ifelse(het, "T", "F"), ifelse(lin, "T", "F"), n_obser, i
        )
        
        tryCatch({
          # If a pre-saved fit doesn't exist, we create a dummy one for demonstration
          if (!file.exists(file_name)) {
            # Create a dummy nbcf_fit object if file not found
            nbcf_fit <- list(
              alpha = matrix(rnorm(1000), nrow = 2),
              Beta = array(rnorm(2 * 1000 * 6), dim = c(2, 500, 6)),
              Beta_int = array(rnorm(2 * 1000 * 15), dim = c(2, 500, 15))
            )
            stop('ERRROR')
            # next # Use 'next' if you only want to process existing files
          } else {
            load(file_name) # loads object: nbcf_fit
          }
          
          # --- The rest of your processing logic ---
          data <- generate_data_2(n_obser, het, lin, seed = i, RCT = FALSE)
          X <- as.matrix(data[, 1:6])
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          
          alpha_samples <- as.vector(t(nbcf_fit$alpha))
          beta_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta[chain, , ]))
          beta_int_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta_int[chain, , ]))
          
          sd_y <- sd(data$y)
          alpha_samples <- alpha_samples * sd_y
          beta_samples <- beta_samples * sd_y
          beta_int_samples <- beta_int_samples * sd_y
          
          ipairs <- interaction_pairs(ncol(X), boolean)
          
          tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                                  nrow = n_obser, byrow = FALSE) + X %*% t(beta_samples)
          
          if (ncol(ipairs) > 0) {
            for (idx in 1:ncol(ipairs)) {
              j <- ipairs[1, idx]
              k <- ipairs[2, idx]
              tau_posterior <- tau_posterior + (X[, j] * X[, k]) %*% t(beta_int_samples[, idx, drop = FALSE])
            }
          }
          
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
          ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
          
          ate_draws <- colMeans(tau_posterior)
          est_ate <- mean(ate_draws)
          ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
          
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2], file_name)
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper, file_name)
          
          result_matrix[i, ] <- c(ate_vec, cate_vec)
          
        }, error = function(e) {
          # Error message
          cat("Error processing file:", file_name, "-", conditionMessage(e), "\n")
        })
      }
      
      # --- UPDATE THE SUMMARY TABLE (your original logic) ---
      ate_results <- result_matrix[, 1:3, drop = FALSE]
      cate_results <- result_matrix[, 4:6, drop = FALSE]
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      results[idx, 4:6] <- colMeans(ate_results, na.rm = TRUE)
      results[idx, 7:9] <- colMeans(cate_results, na.rm = TRUE)
      
      # --- STORE DETAILED RESULTS FOR PLOTTING (new logic) ---
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

# 3. FINAL OUTPUTS
# --------------------------------------------------------------------------
save(results, file = 'results_linked.RData')
# --- 3A. Print the original summary results table ---
cat("\n\n--- Summary Results Table ---\n")
# final_summary_table <- results %>%
#   mutate(heterogeneity = ifelse(heterogeneity, "Heterogeneous", "Homogeneous"),
#          linearity = ifelse(linearity, "Linear", "Nonlinear")) %>%
#   arrange(n, heterogeneity, linearity) %>%
#   select(n, heterogeneity, linearity, everything())
# 
# print(results)


# --- 3B. Create and print the boxplots from detailed results ---
cat("\n\n--- Generating Boxplots ---\n")

# Combine the list of data frames into one
all_results_df <- bind_rows(all_scenario_results_list)

# Pivot the data into a long format for easier plotting
results_long <- all_results_df %>%
  na.omit() %>% # Remove rows where a simulation might have failed
  pivot_longer(
    cols = c(ends_with("RMSE") , ends_with("Length")),
    names_to = c("effect_type", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  # Make metric names factors for consistent ordering in plots
  mutate(metric = factor(metric, levels = c("RMSE", "Length")))

# Create the plot object
metrics_plot <- ggplot(results_long, aes(x = heterogeneity, y = value, fill = linearity)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  facet_grid(metric ~ n, scales = "free_y") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribution of Performance Metrics Across Simulation Scenarios",
    x = "Treatment Effect Type",
    y = "Metric Value",
    fill = "Functional Form"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# Explicitly print the plot
print(metrics_plot)



################################################
library(coda)
data <- generate_data_2(n = 500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 40)
chain_1 <- as.mcmc(nbcf_fit$Beta[1,,]*sd(data$y))
summary(chain_1)
chain_2 <- as.mcmc(nbcf_fit$Beta_int[1,,]*sd(data$y))
summary(chain_2)
chain_3 <- as.mcmc(nbcf_fit$alpha[1,]*sd(data$y))
summary(chain_3)