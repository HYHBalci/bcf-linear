# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
source('R/simul_1.R')

# --- Helper Functions (copied from original script) ---
compute_mode <- function(x) {
  # Using mean as a proxy for mode per original code.
  return(mean(x))
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
scenario_n <- 500
data <- generate_data_2(scenario_n, is_te_hetero = T, is_mu_nonlinear = T, seed = 25, RCT = F, z_diff = 0.5, BCF = T,  sigma_sq =1)

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = NULL, 
  sigma2_global_shape = 0, sigma2_global_scale = 0, 
  variable_weights = NULL, propensity_covariate = "none", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = -1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = TRUE, 
  probit_outcome_model = FALSE
)

nbcf_fit <- bcf(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z) - as.numeric(data$pi_x), 
  propensity_train = as.numeric(data$pi_x),
  num_gfr = 25, 
  num_burnin = 0, 
  num_mcmc = 500,
  general_params = general_params_default
)
true_cate <- data$tau
true_ate <- mean(true_cate)
tau_posterior <- nbcf_fit$tau_hat_train
save(nbcf_fit, file = "nbcf_fit_one_BCF.RData")
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
    title = "True vs. Estimated CATE with 95% Credible Intervals for BCF-standard.",
    subtitle = "Points represent posterior means, error bars represent 95% credible intervals. n = 500.",
    x = "True Treatment Effect",
    y = "Estimated Treatment Effect"
  ) +
  theme_minimal() +
  coord_fixed() # Ensures the 1:1 line is at a 45-degree angle

print(cate_plot)

#####################################################################################

# 5. FRIEDMAN PARTIAL DEPENDENCE PLOT (with 95% Credible Intervals)
# --------------------------------------------------------------------------
cat("\nGenerating Friedman partial dependence plot with credible intervals...\n")

# --- Define the grid ---
grid_size <- 25
x1_grid <- seq(min(data[, 1]), max(data[, 1]), length.out = grid_size)
x4_grid <- c(-1, 1) # x4 is binary

# --- Create matrices to store the results ---
# We now need to store the point estimate (mean) and the interval bounds
pd_mean <- matrix(NA, nrow = grid_size, ncol = length(x4_grid))
pd_lower <- matrix(NA, nrow = grid_size, ncol = length(x4_grid))
pd_upper <- matrix(NA, nrow = grid_size, ncol = length(x4_grid))

# --- Loop through each point in the grid ---
for (i in 1:grid_size) {
  for (j in 1:length(x4_grid)) {
    # 1. Create a temporary copy of the original data
    temp_data <- data
    
    # 2. Set x1 and x4 to the current grid values
    temp_data[, 2] <- x1_grid[i]
    temp_data[, 4] <- x4_grid[j]
    
    # Create the new covariate matrix
    X_new <- as.matrix(sapply(temp_data[, c(1:6)], as.numeric))
    
    # 3. Predict potential outcomes for BOTH treated and control groups
    preds_treated <- predict(
      nbcf_fit, 
      X = X_new,
      Z = rep(1, nrow(X_new)),
      propensity = as.numeric(temp_data$pi_x)
    )
    preds_control <- predict(
      nbcf_fit,
      X = X_new,
      Z = rep(0, nrow(X_new)),
      propensity = as.numeric(temp_data$pi_x)
    )
    
    # 4. Calculate the posterior distribution of the CATE
    # This is a matrix of [num_observations x num_mcmc_samples]
    cate_posterior_draws <- preds_treated$y_hat - preds_control$y_hat
    
    # 5. To get the PDP uncertainty, we average over observations for each MCMC draw
    # This gives the posterior distribution of the average CATE for this grid point.
    pdp_posterior_dist <- colMeans(cate_posterior_draws)
    
    # 6. Store the mean (for the line) and quantiles (for the ribbon)
    pd_mean[i, j] <- mean(pdp_posterior_dist)
    ci <- quantile(pdp_posterior_dist, probs = c(0.025, 0.975))
    pd_lower[i, j] <- ci[1]
    pd_upper[i, j] <- ci[2]
  }
}

# --- Prepare the results for plotting ---
pd_df <- expand.grid(x1 = x1_grid, x4 = x4_grid)
pd_df$cate <- as.vector(pd_mean)
pd_df$ci_lower <- as.vector(pd_lower)
pd_df$ci_upper <- as.vector(pd_upper)

# --- Create the line plot with a confidence ribbon ---
pdp_plot <- ggplot(pd_df, aes(x = x1, y = cate, color = as.factor(x4), group = as.factor(x4))) +
  # Add the ribbon layer for the credible interval
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = as.factor(x4)), alpha = 0.2, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "Value of x4", aesthetics = c("color", "fill")) +
  labs(
    title = "Partial Dependence Plot with 95% Credible Intervals",
    subtitle = "Shows the marginal effect of x2 on CATE at each level of binary x4.",
    x = "Value of x2 (Continuous)",
    y = "Average CATE (Partial Dependence)"
  ) +
  theme_minimal()

print(pdp_plot)

# 6. FRIEDMAN PARTIAL DEPENDENCE PLOT FOR MU (with 95% Credible Intervals)
# --------------------------------------------------------------------------
cat("\nGenerating Friedman partial dependence plot for mu(x) with credible intervals...\n")

# --- Define the grid for x3 ---
grid_size <- 25
# Assuming the 3rd column of your original 'data' dataframe is x3
x3_grid <- seq(min(data[, 3]), max(data[, 3]), length.out = grid_size)

# --- Create vectors to store the results ---
# We are only varying one predictor, so vectors are sufficient.
pd_mean_mu <- numeric(grid_size)
pd_lower_mu <- numeric(grid_size)
pd_upper_mu <- numeric(grid_size)

# --- Loop through each point in the grid ---
for (i in 1:grid_size) {
  # 1. Create a temporary copy of the original data
  temp_data <- data
  
  # 2. Set x3 to the current grid value for all observations
  temp_data[, 3] <- x3_grid[i]
  
  # Create the new covariate matrix
  X_new <- as.matrix(sapply(temp_data[, c(1:6)], as.numeric))
  
  # 3. Predict outcomes. For mu, we only need to predict once as it's the prognostic
  #    effect, which is independent of the treatment assignment Z.
  preds <- predict(
    nbcf_fit, 
    X = X_new,
    # Z is a required argument but does not affect mu_hat
    Z = rep(0, nrow(X_new)), 
    propensity = as.numeric(temp_data$pi_x)
  )
  
  # 4. Extract the posterior distribution of the prognostic effect (mu)
  # This is a matrix of [num_observations x num_mcmc_samples]
  mu_posterior_draws <- preds$mu_hat
  
  # 5. To get the PDP uncertainty, we average over observations for each MCMC draw
  # This gives the posterior distribution of the average mu for this grid point.
  pdp_mu_posterior_dist <- colMeans(mu_posterior_draws)
  
  # 6. Store the mean (for the line) and quantiles (for the ribbon)
  pd_mean_mu[i] <- mean(pdp_mu_posterior_dist)
  ci_mu <- quantile(pdp_mu_posterior_dist, probs = c(0.025, 0.975))
  pd_lower_mu[i] <- ci_mu[1]
  pd_upper_mu[i] <- ci_mu[2]
}

# --- Prepare the results for plotting ---
pd_mu_df <- data.frame(
  x3 = x3_grid,
  mu = pd_mean_mu,
  ci_lower = pd_lower_mu,
  ci_upper = pd_upper_mu
)

# --- Create the line plot with a confidence ribbon ---
pdp_mu_plot <- ggplot(pd_mu_df, aes(x = x3, y = mu)) +
  # Add the ribbon layer for the credible interval
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "skyblue", alpha = 0.3, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.1, color = "navy") +
  geom_point(size = 2, color = "navy") +
  labs(
    title = "Partial Dependence Plot for Prognostic Effect (μ)",
    subtitle = "Shows the marginal effect of x3 on the expected outcome.",
    x = "Value of x3 (Continuous)",
    y = "Average Prognostic Effect (Partial Dependence on μ)"
  ) +
  theme_minimal()

print(pdp_mu_plot)

# 7. ADD TRUE MU FUNCTION TO THE PLOT
# --------------------------------------------------------------------------
cat("\nAdding the true mu function to the PDP plot...\n")

# --- Calculate the true mu values based on the provided formula ---
# The formula is: mu = -6 + g_x5 + 6*abs(x3 - 1)
# To plot this against x3, we need to integrate out the effect of x5.
# We do this by taking the average of the g_x5 term over all observations.

# Assuming g_x5 = sin(pi * x5) and x5 is the 5th column of the data
# If your g(x) function is different, replace it here.
avg_g_x5 <- mean(sin(pi * data[, 5]))

# Now, calculate the true mu for each point on our x3 grid
true_mu_values <- -6 + avg_g_x5 + 6 * abs(pd_mu_df$x3 - 1)

# Add the true values to our plotting dataframe
pd_mu_df$true_mu <- true_mu_values


# --- Create the plot with both estimated PDP and true function ---
pdp_mu_plot_with_true <- ggplot(pd_mu_df, aes(x = x3)) +
  # Add the ribbon layer for the credible interval of the estimate
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "skyblue", alpha = 0.3) +
  
  # Add the line for the estimated PDP
  geom_line(aes(y = mu, color = "Estimated PDP"), linewidth = 1.1) +
  geom_point(aes(y = mu, color = "Estimated PDP"), size = 2) +
  
  # Add the line for the true mu function
  geom_line(aes(y = true_mu, color = "True Function"), linewidth = 1.2, linetype = "dashed") +
  
  # Manually define the colors and legend
  scale_color_manual(
    name = "Function", 
    values = c("Estimated PDP" = "navy", "True Function" = "red")
  ) +
  
  labs(
    title = "Partial Dependence Plot vs. True Prognostic Function (μ)",
    subtitle = "Comparing the NBCF estimate with the true marginal effect of x3.",
    x = "Value of x3 (Continuous)",
    y = "Average Prognostic Effect (μ)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(pdp_mu_plot_with_true)

# 8. VISUALIZATION OF MU(X) ERRORS
# --------------------------------------------------------------------------
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
    subtitle = "Comparing model estimates of μ(X) against true values. n = 500.",
    x = "True Prognostic Effect (μ)",
    y = "Estimated Prognostic Effect (μ)"
  ) +
  theme_minimal() +
  coord_fixed()

print(mu_error_plot)
