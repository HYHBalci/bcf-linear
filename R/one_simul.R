# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
source('R/simul_1.R')
# library(stochtree) # Uncomment if standardize_X_by_index or other functions require it

# --- Source necessary files ---
# Ensure this path points to the file containing generate_data_2 and standardize_X_by_index
# source("R/simul_1.R")

# --- Helper Functions (copied from original script) ---
compute_mode <- function(x) {
  # Using mean as a proxy for mode per original code.
  return(median(x))
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(data.frame(rmse = rmse, coverage = coverage, interval_length = interval_length))
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
  if (length(interaction_list) == 0) return(matrix(nrow = 2, ncol = 0))
  return(do.call(cbind, interaction_list))
}

# 2. CONFIGURATION FOR SINGLE RUN ANALYSIS
# --------------------------------------------------------------------------



# 3. CORE EVALUATION LOGIC
# --------------------------------------------------------------------------
num_chains <- 1
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = T, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE,   #30
  keep_every = 1, num_chains = num_chains, verbose = T, 
  sample_global_prior = "half-cauchy", unlink = T, propensity_seperate = "none", gibbs = F, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = F, interaction_rule = "continuous_or_binary",standardize_cov = F, simple_prior = T, save_partial_residual = F, regularize_ATE = F,
  sigma_residual = 0, hn_scale = 0, use_ncf = F
)
prognostic_forest_params_new <- list(
  num_trees = 50, 
  beta = 3.0, 
  max_depth = 5, 
  alpha = 0.95, 
  min_samples_leaf = 5,
  sample_sigma2_leaf = TRUE, 
  sigma2_leaf_init = NULL, 
  sigma2_leaf_shape = 3, 
  sigma2_leaf_scale = NULL, 
  keep_vars = NULL, 
  drop_vars = NULL
)
#'all'
#
scenario_n <- 2000
data <- generate_data_2(scenario_n, is_te_hetero = T, is_mu_nonlinear = T, seed = 1848, RCT = F, z_diff = 0.5, BCF = F,  sigma_sq =1)


nbcf_fit <- bcf_linear_probit(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z),
  propensity_train = as.numeric(data$pi_x),
  num_gfr = 25, 
  num_burnin = 1500, 
  num_mcmc = 1000, 
  general_params = general_params_default,
  prognostic_forest_params = prognostic_forest_params_new
)

save(nbcf_fit, file = 'nbcf_fit_one_BCF_linear.RData')
X <- as.matrix(data[, 1:6])
true_cate <- data$tau
true_ate <- mean(true_cate)
  
# --- Prepare Interaction Features ---
# Calculate boolean vector based on the characteristics of the loaded data X
# This determines which covariates are included in interaction calculations.
X_info <- standardize_X_by_index(X_initial = X, process_data = F, interaction_rule = "continuous", cat_coding_method = "difference")
boolean_vector <- as.logical(as.numeric(X_info$X_final_var_info$is_continuous)) + as.logical(as.numeric(X_info$X_final_var_info$is_binary))
ipairs <- interaction_pairs(ncol(X), boolean_vector)
  
 # --- Extract and Process Posterior Samples ---
cat("Extracting and processing posterior samples...\n")
  # Combine chains from model fit object
alpha_samples <- as.vector(t(nbcf_fit$alpha))
beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))
  
  # Rescale samples if necessary (matching original script logic)
sd_y <- sd(data$y)
alpha_samples <- alpha_samples * sd_y
beta_samples <- beta_samples * sd_y
beta_int_samples <- beta_int_samples * sd_y
  
  # --- Calculate Posterior CATE Distribution ---
  cat("Reconstructing CATE posterior distribution...\n")
  # Start with intercept (alpha) for all observations
  tau_posterior <- matrix(rep(alpha_samples, each = scenario_n),
                          nrow = scenario_n, byrow = FALSE)
  
  # Add main effects: X * beta
  tau_posterior <- tau_posterior + X %*% t(beta_samples)
  
  # Add interaction effects: sum(X_j * X_k * beta_jk)
  if (ncol(ipairs) > 0) {
    for (idx in 1:ncol(ipairs)) {
      j <- ipairs[1, idx]
      k <- ipairs[2, idx]
      interaction_term_values <- X[, j] * X[, k]
      tau_posterior <- tau_posterior + interaction_term_values %*% t(beta_int_samples[, idx, drop = FALSE])
    
    }
  }
  
  # --- Compute Point Estimates and Intervals ---
  cat("Calculating point estimates and credible intervals...\n")
  # CATE estimates
  tau_estimate <- apply(tau_posterior, 1, compute_mode)
  ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # ATE estimates
  ate_draws <- colMeans(tau_posterior)
  est_ate <- mean(ate_draws)
  ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # --- Calculate Final Metrics ---
  ate_metrics <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
  cate_metrics <- compute_metrics(true_cate, tau_estimate, ci_tau_lower, ci_tau_upper)
  
  # --- Display Results ---
  cat("\n--- Evaluation Results ---\n")
  cat("ATE Metrics:\n")
  print(ate_metrics)
  cat("\nCATE Metrics:\n")
  print(cate_metrics)



# 4. VISUALIZATION FOR SINGLE RUN
# --------------------------------------------------------------------------
# Create a data frame for plotting true vs. estimated CATE values
plot_data <- data.frame(
  true_cate = true_cate,
  estimated_cate = tau_estimate,
  ci_lower = ci_tau_lower,
  ci_upper = ci_tau_upper
)

  
# --- Create the text label for the plot ---
rmse_label <- paste(
    paste("ATE RMSE:", round(ate_metrics$rmse, 3)),
    paste("CATE RMSE:", round(cate_metrics$rmse, 3)),
    sep = "\n" # Put them on separate lines
  )
  
# Scatter plot of True CATE vs. Estimated CATE with 95% credible intervals
cate_plot <- ggplot(plot_data, aes(x = true_cate, y = estimated_cate)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0, alpha = 0.3, color = "gray50") +
    geom_point(alpha = 0.6, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
    # --- ADDED: Annotate plot with RMSE values ---
    annotate(
      geom = "text",
      x = -Inf,        # Position horizontally at the far left
      y = Inf,         # Position vertically at the very top
      label = rmse_label,
      hjust = -0.1,    # Adjust to push text right from the edge
      vjust = 1.5,     # Adjust to push text down from the edge
      size = 4,
      color = "black",
      fontface = "bold"
    ) +
    labs(
      title = "True vs. Estimated CATE with 95% Credible Intervals for BCF-linear. local shrinkage.",
      subtitle = "Points represent posterior means, error bars represent 95% credible intervals. n = 500",
      x = "True Treatment Effect",
      y = "Estimated Treatment Effect"
    ) +
    theme_minimal() +
    coord_fixed() # Ensures the 1:1 line is at a 45-degree angle
  
print(cate_plot)

cat("\nGenerating plot to investigate errors for mu(X)...\n")

# --- Extract true values and posterior draws for mu ---
true_mu <- data$mu
mu_posterior <- nbcf_fit$mu_hat_train

# --- Compute Point Estimates and Intervals for mu ---
cat("Calculating point estimates and credible intervals for mu...\n")
mu_estimate <- apply(mu_posterior, 1, compute_mode)
ci_mu_lower <- apply(mu_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_mu_upper <- apply(mu_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)

# --- Calculate Metrics for mu ---
mu_metrics <- compute_metrics(true_mu, mu_estimate, ci_mu_lower, ci_mu_upper)
cat("\nPrognostic Effect (mu) Metrics:\n")
print(mu_metrics)

# --- Prepare data for plotting ---
mu_plot_data <- data.frame(
  true_mu = true_mu,
  estimated_mu = mu_estimate,
  ci_lower = ci_mu_lower,
  ci_upper = ci_mu_upper
)

# --- Create the text label for the plot ---
mu_rmse_label <- paste("Mu RMSE:", round(mu_metrics$rmse, 3))

# --- Create the scatter plot ---
mu_error_plot <- ggplot(mu_plot_data, aes(x = true_mu, y = estimated_mu)) +
  # Add credible intervals as error bars
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0, alpha = 0.3, color = "gray50") +
  # Add point estimates
  geom_point(alpha = 0.6, color = "darkgreen") +
  # Add the line of perfect correspondence (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # Add a linear fit to the points
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  # Add the RMSE metric as an annotation
  annotate(
    geom = "text",
    x = -Inf, 
    y = Inf,  
    label = mu_rmse_label,
    hjust = -0.1, 
    vjust = 1.5,  
    size = 4,
    color = "black",
    fontface = "bold"
  ) +
  labs(
    title = "True vs. Estimated Prognostic Effect (μ) with 95% Credible Intervals",
    subtitle = "Comparing model estimates of μ(X) against true values. n = 500",
    x = "True Prognostic Effect (μ)",
    y = "Estimated Prognostic Effect (μ)"
  ) +
  theme_minimal() +
  coord_fixed()

print(mu_error_plot)

shapley <- compute_shapley_all(X = as.matrix(sapply(data[, c(1:6)], as.numeric)), beta_post = nbcf_fit$Beta[1,,], nbcf_fit$Beta_int[1,,], boolean_interaction =c(T,T,T,T,F,F))


########################################### PREDICT MU ##################################
# Preprocess X for the forest
X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
Z <- as.numeric(data$z)
object <- nbcf_fit
model_params <- object$model_params
y_std <- model_params$outcome_scale
y_bar <- model_params$outcome_mean
num_samples <- model_params$num_samples
n_obs <- nrow(X)

propensity <- as.numeric(data$pi_x)

# --- 1. Prepare Data & Predict for FOREST Model (mu_hat) ---


train_set_metadata_forest <- object$train_set_metadata
X_forest <- preprocessPredictionData(X, train_set_metadata_forest)

# Add propensities if needed by the *forest*
X_forest_combined <- X_forest
if (model_params$propensity_covariate != "none") {
  X_forest_combined <- cbind(X_forest, propensity)
}

# Create prediction dataset for forest
# Z is just a placeholder, mu_hat prediction doesn't use it.
forest_dataset_pred <- createForestDataset(X_forest_combined, Z) 

# Get mu_hat predictions (n_obs x num_samples)
mu_hat <- object$forests_mu$predict(forest_dataset_pred) * y_std + y_bar