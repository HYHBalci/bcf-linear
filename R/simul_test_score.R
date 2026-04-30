# 1. SETUP & IMPORTS
# --------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree) # Assumes bcf_restricted_score_test is available
library(MASS)

# Load your Data Generating Process (DGP) functions
source('R/simul_1.R')

# Set working directory if needed (defaults to current R session directory)
# setwd("path/to/your/save/folder")

# 2. SIMULATION CONFIGURATION
# --------------------------------------------------------------------------
num_sims <- 50          # Number of datasets to simulate per configuration (increase for final paper)

# Define the grid of all scenarios we want to test
sim_grid <- expand.grid(
  n = c(250, 500, 750, 1500),
  is_mu_nonlinear = c(FALSE, TRUE),
  is_te_hetero = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

# Add readable labels for the plots and tables
sim_grid <- sim_grid %>%
  mutate(
    mu_label = ifelse(is_mu_nonlinear, "Non-Linear DGP", "Linear DGP"),
    te_label = ifelse(is_te_hetero, "Heterogeneous (Alternative)", "Homogeneous (Null)"),
    config_name = paste(te_label, "|", mu_label, "| n =", n)
  )

# Initialize list to store results
results_list <- list()

# Define general BCF parameters (verbose = FALSE to prevent console flooding)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, sample_sigma2_global = TRUE, 
  sigma2_global_init = NULL, sigma2_global_shape = 1, sigma2_global_scale = 0.001, 
  variable_weights = NULL, propensity_covariate = "none", adaptive_coding = FALSE, 
  rfx_prior_var = NULL, random_seed = -1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = FALSE, probit_outcome_model = FALSE, 
  standardize_cov = FALSE, interaction_rule = "continuous"
)

# 3. CORE SIMULATION LOOP
# --------------------------------------------------------------------------
total_configs <- nrow(sim_grid)

for (s in 1:total_configs) {
  current_config <- sim_grid[s, ]
  
  cat(sprintf("\n=== Starting Config %d of %d: %s ===\n", s, total_configs, current_config$config_name))
  
  for (i in 1:num_sims) {
    cat(sprintf("  Simulating dataset %d of %d...\r", i, num_sims))
    
    # Set a unique seed for each simulation to ensure proper variance across runs
    sim_seed <- 10000 * s + i 
    
    # Generate Data (RCT setup, Z is centered using z_diff = 0.5)
    data <- generate_data_2(
      n = current_config$n, 
      is_te_hetero = current_config$is_te_hetero, 
      is_mu_nonlinear = current_config$is_mu_nonlinear, 
      seed = sim_seed, 
      RCT = TRUE, 
      z_diff = 0.5, 
      BCF = FALSE,  
      sigma_sq = 1
    )
    
    # Extract covariate matrix (x1, x2, x3, x4, x5_1, x5_2)
    X_mat <- as.matrix(sapply(data[, c(1:6)], as.numeric))
    
    # Fit the Restricted Score Test Model
    nbcf_fit <- tryCatch({
      bcf_restricted_score_test(
        X_train = X_mat,
        y_train = as.numeric(data$y),
        Z_train = as.numeric(data$z), 
        propensity_train = rep(0, nrow(data)), 
        num_gfr = 25, 
        num_burnin = 500, 
        num_mcmc = 2000,
        general_params = general_params_default
      )
    }, error = function(e) {
      cat("\n  [Error in sim", i, ":", e$message, "]\n")
      return(NULL)
    })
    
    if (is.null(nbcf_fit)) next # Skip to next iteration if model failed
    
    # Extract Posterior P-values
    quad_pvals <- as.vector(nbcf_fit$pval_quad)
    max_pvals <- as.vector(nbcf_fit$pval_max)
    
    # Calculate Posterior Median P-values
    median_quad <- median(quad_pvals, na.rm = TRUE)
    median_max <- median(max_pvals, na.rm = TRUE)
    
    # Append to results
    results_list[[length(results_list) + 1]] <- data.frame(
      n = current_config$n,
      mu_type = current_config$mu_label,
      te_type = current_config$te_label,
      sim_id = i,
      median_pval_quad = median_quad,
      median_pval_max = median_max,
      reject_quad = ifelse(median_quad < 0.05, 1, 0),
      reject_max = ifelse(median_max < 0.05, 1, 0)
    )
  }
}
cat("\n\n=== SIMULATION COMPLETE ===\n")

# 4. AGGREGATE AND SUMMARIZE RESULTS
# --------------------------------------------------------------------------
final_results <- bind_rows(results_list)

# Calculate summary statistics grouped by Sample Size, DGP type, and Null/Alt
summary_stats <- final_results %>%
  group_by(te_type, mu_type, n) %>%
  summarise(
    simulations_run = n(),
    avg_median_quad = mean(median_pval_quad, na.rm = TRUE),
    avg_median_max = mean(median_pval_max, na.rm = TRUE),
    rejection_rate_quad = mean(reject_quad, na.rm = TRUE),
    rejection_rate_max = mean(reject_max, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(te_type, mu_type, n)

print(summary_stats)


# 5. VISUALIZATION
# --------------------------------------------------------------------------

# Plot 1: Rejection Rates (Power & Type I Error) across Sample Sizes
power_plot <- summary_stats %>%
  ggplot(aes(x = n, y = rejection_rate_quad, color = mu_type, group = mu_type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ te_type, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Score Test Rejection Rates (Quadratic Statistic)",
    subtitle = "Left: Type I Error (False Positives) | Right: Statistical Power (True Positives)",
    x = "Sample Size (n)",
    y = "Rejection Rate (p < 0.05)",
    color = "DGP Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot 2: P-Value Distributions (Faceted by n and DGP)
p_val_plot <- final_results %>%
  ggplot(aes(x = median_pval_quad, fill = te_type)) +
  geom_histogram(alpha = 0.7, bins = 20, color = "black", position = "identity") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.8) +
  facet_grid(mu_type ~ n) +
  labs(
    title = "Distribution of Median Quadratic P-Values across Configurations",
    subtitle = "Columns: Sample Size | Rows: Prognostic Function Type",
    x = "Posterior Median P-Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Homogeneous (Null)" = "grey50", "Heterogeneous (Alternative)" = "blue")) +
  theme(legend.position = "bottom", legend.title = element_blank())


# 6. SAVE ARTIFACTS LOCALLY
# --------------------------------------------------------------------------
cat("\nSaving all artifacts to the local working directory...\n")

# Save tables as .RData for easy reloading later
save(final_results, summary_stats, file = "score_test_simulation_results.RData")

# Save visualizations as high-quality PNGs
ggsave("score_test_power_typeI.png", plot = power_plot, width = 10, height = 6, dpi = 300)
ggsave("score_test_pval_distributions.png", plot = p_val_plot, width = 12, height = 8, dpi = 300)

cat("Success! Saved:\n")
cat(" - score_test_simulation_results.RData\n")
cat(" - score_test_power_typeI.png\n")
cat(" - score_test_pval_distributions.png\n")