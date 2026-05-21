# 1. SETUP & IMPORTS
# --------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree) # Assumes updated bcf_restricted_score_test is available
library(MASS)
library(future)
library(furrr)

# Load your Data Generating Process (DGP) functions
source('R/simul_1.R')

# Set working directory if needed
# setwd("path/to/your/save/folder")

# 2. SIMULATION CONFIGURATION & PARALLEL SETUP
# --------------------------------------------------------------------------
num_sims <- 100          # Number of datasets to simulate per configuration

# Set up parallel backend (leaves 1 core free to keep your system responsive)
plan(multisession, workers = availableCores() - 1)

# Define the base grid
sim_grid <- expand.grid(
  n = c(250, 500, 750, 1500),
  is_mu_nonlinear = c(FALSE, TRUE),
  is_te_hetero = c(FALSE, TRUE),
  stringsAsFactors = FALSE
) %>%
  mutate(
    mu_label = ifelse(is_mu_nonlinear, "Non-Linear DGP", "Linear DGP"),
    te_label = ifelse(is_te_hetero, "Heterogeneous (Alternative)", "Homogeneous (Null)"),
    config_name = paste(te_label, "|", mu_label, "| n =", n),
    config_id = row_number()
  )

# Flatten the nested loops into a single task grid (one row = one simulation)
full_task_grid <- sim_grid %>%
  crossing(sim_id = 1:num_sims) %>%
  mutate(sim_seed = 10000 * config_id + sim_id)

# Define general BCF parameters 
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, sample_sigma2_global = TRUE, 
  sigma2_global_init = NULL, sigma2_global_shape = 1, sigma2_global_scale = 0.001, 
  variable_weights = NULL, propensity_covariate = "none", adaptive_coding = FALSE, 
  rfx_prior_var = NULL, random_seed = -1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = FALSE, probit_outcome_model = FALSE, 
  standardize_cov = FALSE, interaction_rule = "continuous"
)

# 3. DEFINE THE WORKER FUNCTION
# --------------------------------------------------------------------------
run_single_sim <- function(n, is_te_hetero, is_mu_nonlinear, mu_label, te_label, sim_id, sim_seed, config_name) {
  
  # Generate Data
  data <- generate_data_2(
    n = n, 
    is_te_hetero = is_te_hetero, 
    is_mu_nonlinear = is_mu_nonlinear, 
    seed = sim_seed, 
    RCT = TRUE, 
    z_diff = 0.5, 
    BCF = FALSE,  
    sigma_sq = 1
  )
  
  X_mat <- as.matrix(sapply(data[, c(1:6)], as.numeric))
  
  # Pass the explicit seed to the BCF params to ensure exact reproducibility 
  params_for_run <- general_params_default
  params_for_run$random_seed <- sim_seed
  
  # Fit Model
  nbcf_fit <- tryCatch({
    bcf_restricted_score_test(
      X_train = X_mat,
      y_train = as.numeric(data$y),
      Z_train = as.numeric(data$z), 
      propensity_train = rep(0, nrow(data)), 
      num_gfr = 25, 
      num_burnin = 500, 
      num_mcmc = 2000,
      use_rao_blackwell = TRUE,
      general_params = params_for_run
    )
  }, error = function(e) {
    return(NULL)
  })
  
  # Handle failures gracefully
  if (is.null(nbcf_fit)) {
    return(data.frame(
      n = n, mu_type = mu_label, te_type = te_label, sim_id = sim_id,
      pval_quad = NA, pval_max = NA, reject_quad = NA, reject_max = NA
    ))
  }
  
  # Extract P-values based on model output structure
  if (!is.null(nbcf_fit$rb_pval_quad)) {
    val_quad <- median(nbcf_fit$rb_pval_quad, na.rm = TRUE)
    val_max  <- median(nbcf_fit$rb_pval_max, na.rm = TRUE)
  } else {
    val_quad <- median(as.vector(nbcf_fit$pval_quad_samples), na.rm = TRUE)
    val_max  <- median(as.vector(nbcf_fit$pval_max_samples), na.rm = TRUE)
  }
  
  # Return single row dataframe
  data.frame(
    n = n,
    mu_type = mu_label,
    te_type = te_label,
    sim_id = sim_id,
    pval_quad = val_quad,
    pval_max = val_max,
    reject_quad = ifelse(val_quad < 0.05, 1, 0),
    reject_max = ifelse(val_max < 0.05, 1, 0)
  )
}

# 4. EXECUTE IN PARALLEL
# --------------------------------------------------------------------------
cat("\n=== STARTING PARALLEL SIMULATION ===\n")
cat(sprintf("Running %d total tasks across %d workers...\n", nrow(full_task_grid), nbrOfWorkers()))

# Define exports and packages for parallel workers
worker_options <- furrr_options(
  globals = c("run_single_sim", "generate_data_2", "general_params_default"),
  packages = c("dplyr", "MASS", "stochtree"),
  seed = TRUE
)

# Map the function over the grid rows with a progress bar
final_results <- future_pmap_dfr(
  .l = list(
    n = full_task_grid$n,
    is_te_hetero = full_task_grid$is_te_hetero,
    is_mu_nonlinear = full_task_grid$is_mu_nonlinear,
    mu_label = full_task_grid$mu_label,
    te_label = full_task_grid$te_label,
    sim_id = full_task_grid$sim_id,
    sim_seed = full_task_grid$sim_seed,
    config_name = full_task_grid$config_name
  ),
  .f = run_single_sim,
  .options = worker_options, 
  .progress = TRUE
)

# Return to sequential processing when done to free up RAM
plan(sequential) 
cat("\n=== SIMULATION COMPLETE ===\n\n")

# 5. AGGREGATE AND SUMMARIZE RESULTS
# --------------------------------------------------------------------------
# Calculate summary statistics grouped by Sample Size, DGP type, and Null/Alt
summary_stats <- final_results %>%
  group_by(te_type, mu_type, n) %>%
  summarise(
    simulations_run = sum(!is.na(pval_quad)),
    failed_models = sum(is.na(pval_quad)),
    avg_pval_quad = mean(pval_quad, na.rm = TRUE),
    avg_pval_max = mean(pval_max, na.rm = TRUE),
    rejection_rate_quad = mean(reject_quad, na.rm = TRUE),
    rejection_rate_max = mean(reject_max, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(te_type, mu_type, n)

print(summary_stats)

# 6. SAVE ARTIFACTS LOCALLY
# --------------------------------------------------------------------------
cat("\nSaving artifacts to the local working directory...\n")

# Save tables as .RData
save(final_results, summary_stats, file = "score_test_simulation_results.RData")
cat(" - Saved: score_test_simulation_results.RData\n")

# Optional: Generate and save visualizations
power_plot <- summary_stats %>%
  ggplot(aes(x = n, y = rejection_rate_quad, color = mu_type, group = mu_type)) +
  geom_line(linewidth = 1) +
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

ggsave("score_test_power_typeI.png", plot = power_plot, width = 10, height = 6, dpi = 300)
cat(" - Saved: score_test_power_typeI.png\n")