# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(Rcpp)
# --- Source necessary files ---
# Make sure these paths are correct for your project structure
source('R/simul_1.R') 
source('R/simul_1.R')
source('R/test_logit_2.R')
source('R/old_linear_linear.R')

# --- Helper Functions ---
compute_mode <- function(x) {
  # Using mean as a proxy for mode, consistent with the original script's approach.
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

# 2. CONFIGURATION AND DATA GENERATION
# --------------------------------------------------------------------------
cat("Generating simulated data...\n")
set.seed(32) # for reproducibility
scenario_n <- 1000
data <- generate_data_2(
  scenario_n, 
  is_te_hetero = TRUE, 
  is_mu_nonlinear = F, # Set to FALSE for a linear prognostic surface
  seed = 131, 
  RCT = FALSE, 
  z_diff = 0.5
)

X <- as.matrix(data[, 1:6])
true_cate <- data$tau
true_mu <- data$mu
true_ate <- mean(true_cate)


# 3. FIT THE BAYESIAN LINEAR MODEL
# --------------------------------------------------------------------------
cat("Fitting the Bayesian linear model with horseshoe priors...\n")
fit_grouped_hs <- fit_grouped_horseshoes_R(
  y_vec = as.numeric(data$y),
  X_mat = X,
  Z_vec = as.numeric(data$z),
  propensity_train = as.numeric(data$pi_x),
  family = "gaussian",
  n_iter = 4000,
  burn_in = 1000,
  num_chains = 2,
  propensity_as_covariate = TRUE,
  method_tau_prognostic = "halfCauchy", tau_prognostic_init = 1,
  method_tau_treatment = "halfCauchy", tau_treatment_init = 1,
  method_tau_overall = "fixed", tau_overall_init = 1,
  alpha_global_prior_sd = 5.0,
  aleph_prior_sd = 5.0,
  thin = 1,
  seed = 231,
  verbose = TRUE,
  ping = 1000,
  standardize_cov = FALSE,
  interaction_rule = "continuous_or_binary",
  cat_coding_method = "difference"
)


# 4. PROCESS POSTERIOR SAMPLES
# --------------------------------------------------------------------------
cat("Extracting and processing posterior samples...\n")
# --- Extract posterior samples for CATE (tau) components ---
aleph_samples <- fit_grouped_hs$aleph       # Intercept for tau(x)
gamma_samples <- fit_grouped_hs$gamma       # Main effects for tau(x)
gamma_int_samples <- fit_grouped_hs$gamma_int # Interaction effects for tau(x)

# --- Extract posterior samples for Prognostic (mu) components ---
alpha_samples <- fit_grouped_hs$alpha       # Intercept for mu(x)
beta_samples <- fit_grouped_hs$beta        # Main effects for mu(x)
beta_int_samples <- fit_grouped_hs$beta_int  # Interaction effects for mu(x)

# --- Prepare Interaction Features ---
X_info <- standardize_X_by_index(X_initial = X, process_data = F, interaction_rule = "continuous", cat_coding_method = "difference")
boolean_vector <- as.logical(as.numeric(X_info$X_final_var_info$is_continuous)) + as.logical(as.numeric(X_info$X_final_var_info$is_binary))
ipairs <- interaction_pairs(ncol(X), boolean_vector)


# --- Reconstruct CATE (tau) Posterior Distribution ---
cat("Reconstructing CATE posterior distribution...\n")
tau_posterior <- matrix(rep(aleph_samples, each = scenario_n), nrow = scenario_n, byrow = FALSE)
tau_posterior <- tau_posterior + X %*% t(gamma_samples)
if (ncol(ipairs) > 0) {
  for (idx in 1:ncol(ipairs)) {
    j <- ipairs[1, idx]
    k <- ipairs[2, idx]
    interaction_term_values <- X[, j] * X[, k]
    tau_posterior <- tau_posterior + interaction_term_values %*% t(gamma_int_samples[, idx, drop = FALSE])
  }
}

# --- Reconstruct Prognostic (mu) Posterior Distribution ---
cat("Reconstructing prognostic (mu) posterior distribution...\n")
mu_posterior <- matrix(rep(alpha_samples, each = scenario_n), nrow = scenario_n, byrow = FALSE)
mu_posterior <- mu_posterior + X %*% t(beta_samples)
if (ncol(ipairs) > 0) {
  for (idx in 1:ncol(ipairs)) {
    j <- ipairs[1, idx]
    k <- ipairs[2, idx]
    interaction_term_values <- X[, j] * X[, k]
    mu_posterior <- mu_posterior + interaction_term_values %*% t(beta_int_samples[, idx, drop = FALSE])
  }
}


# 5. COMPUTE ESTIMATES, INTERVALS, AND METRICS
# --------------------------------------------------------------------------
cat("Calculating point estimates, credible intervals, and metrics...\n")
# --- CATE estimates ---
tau_estimate <- apply(tau_posterior, 1, compute_mode)
ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
cate_metrics <- compute_metrics(true_cate, tau_estimate, ci_tau_lower, ci_tau_upper)

# --- ATE estimates ---
ate_draws <- colMeans(tau_posterior)
est_ate <- mean(ate_draws)
ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
ate_metrics <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])

# --- Prognostic (mu) estimates ---
mu_estimate <- apply(mu_posterior, 1, compute_mode)
ci_mu_lower <- apply(mu_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_mu_upper <- apply(mu_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
mu_metrics <- compute_metrics(true_mu, mu_estimate, ci_mu_lower, ci_mu_upper)

# --- Display Results ---
cat("\n--- Evaluation Results ---\n")
cat("ATE Metrics:\n")
print(ate_metrics)
cat("\nCATE Metrics:\n")
print(cate_metrics)
cat("\nPrognostic Effect (mu) Metrics:\n")
print(mu_metrics)


# 6. VISUALIZATION
# --------------------------------------------------------------------------
# --- Plot for CATE ---
cat("\nGenerating plot for True vs. Estimated CATE...\n")
plot_data_cate <- data.frame(
  true_cate = true_cate,
  estimated_cate = tau_estimate,
  ci_lower = ci_tau_lower,
  ci_upper = ci_tau_upper
)

rmse_label_cate <- paste(
  paste("ATE RMSE:", round(ate_metrics$rmse, 3)),
  paste("CATE RMSE:", round(cate_metrics$rmse, 3)),
  sep = "\n"
)

cate_plot <- ggplot(plot_data_cate, aes(x = true_cate, y = estimated_cate)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0, alpha = 0.3, color = "gray50") +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  annotate(
    geom = "text", x = -Inf, y = Inf, label = rmse_label_cate,
    hjust = -0.1, vjust = 1.5, size = 4, color = "black", fontface = "bold"
  ) +
  labs(
    title = "True vs. Estimated CATE (Bayesian Linear Model)",
    subtitle = paste("Points are posterior means, error bars are 95% credible intervals. n =", scenario_n),
    x = "True Treatment Effect (CATE)",
    y = "Estimated Treatment Effect (CATE)"
  ) +
  theme_minimal() +
  coord_fixed()

print(cate_plot)

# --- Plot for Prognostic Effect (mu) ---
cat("\nGenerating plot for True vs. Estimated Prognostic Effect (mu)...\n")
plot_data_mu <- data.frame(
  true_mu = true_mu,
  estimated_mu = mu_estimate,
  ci_lower = ci_mu_lower,
  ci_upper = ci_mu_upper
)

rmse_label_mu <- paste("Mu RMSE:", round(mu_metrics$rmse, 3))

mu_plot <- ggplot(plot_data_mu, aes(x = true_mu, y = estimated_mu)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.0, alpha = 0.3, color = "gray50") +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  annotate(
    geom = "text", x = -Inf, y = Inf, label = rmse_label_mu,
    hjust = -0.1, vjust = 1.5, size = 4, color = "black", fontface = "bold"
  ) +
  labs(
    title = "True vs. Estimated Prognostic Effect (μ) (Bayesian Linear Model)",
    subtitle = paste("Comparing model estimates of μ(X) against true values. n =", scenario_n),
    x = "True Prognostic Effect (μ)",
    y = "Estimated Prognostic Effect (μ)"
  ) +
  theme_minimal() +
  coord_fixed()

print(mu_plot)

shapley <- compute_shapley_all(X = X, beta_post = fit_grouped_hs$gamma , fit_grouped_hs$gamma_int, boolean_interaction =c(T,T,T,T,F,F))
print(plot_shapley_summary(shapley))