library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")
data <- generate_data_2(500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 18, RCT = T, z_diff = F)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = TRUE, global_shrinkage = T, unlink = T, propensity_seperate = F
)

nbcf_fit <- bcf_linear( X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
                        y_train = as.numeric(data$y),
                        Z_train = as.numeric(data$z), num_gfr = 25, num_burnin =  1000, num_mcmc = 4000, general_params = general_params_default )

library(coda)

library(coda)
for(i in 1:15) {
chain_1 <- as.mcmc(nbcf_fit$Beta_int[1,,i] * sd(data$y))

# Plot density
dens <- density(chain_1)

plot(dens,
     main = expression(paste("Posterior of ", beta[i])),
     xlab = expression(alpha[i]),
     col = "blue",
     lwd = 2)

# Add vertical line at the true value (1)
abline(v = 0, col = "red", lwd = 2, lty = 2)

# Optional: Add a legend
legend("topright",
       legend = c("Posterior", "True Value = 0"),
       col = c("blue", "red"),
       lty = c(1, 2),
       lwd = 2)
}
library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = T, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1884, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = T, 
  global_shrinkage = F, unlink = F, propensity_seperate = F, gibbs = F, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = F, interaction_rule = "all", standardize_cov = F, return_shapley = T
)
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 40, RCT = FALSE, z_diff = T)

nbcf_fit <- bcf_linear_probit(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 2000, 
  general_params = general_params_default
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ---- Run for many observations
compute_shapley_all <- function(X, beta_post, beta_int_post, indices = NULL) {
  if (is.null(indices)) indices <- sample(1:nrow(X), 100)
  all_shapleys <- purrr::map_dfr(indices, ~shapley(.x, X, beta_post, beta_int_post))
  return(all_shapleys)
}
save(nbcf_fit, file = "place_holder.RData")
# ---- Example usage
# Assume: X is your centered covariate matrix
# beta_post: chains x samples x p
# beta_int_post: chains x samples x p_int
# Example: run for 100 observations

shapley_all_df <- compute_shapley_all(as.matrix(sapply(data[, c(1:6)], as.numeric)), nbcf_fit$Beta, nbcf_fit$Beta_int)
# ---- Plotting summary
plot_shapley_summary <- function(shapley_all_df) {
  shapley_summary <- shapley_all_df %>%
    group_by(feature, obs) %>%
    summarise(
      mean = mean(shapley),
      .groups = "drop"
    )
  
  ggplot(shapley_summary, aes(x = mean, y = feature)) +
    geom_jitter(width = 0, height = 0.2, alpha = 0.6) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red") +
    labs(
      title = "Shapley Value Distribution Across Observations",
      x = "Shapley Value (ϕ_ij)",
      y = "Feature"
    ) +
    theme_minimal()
}
library(stochtree)
plot_shapley_vs_covariate <- function(X, beta_post, beta_int_post, indices = NULL, feature = NULL) {
  if (is.null(indices)) indices <- sample(1:nrow(X), 100)
  
  shapley_list <- purrr::map_dfr(indices, function(i) {
    res <- tryCatch(shapley(i, X, beta_post, beta_int_post), error = function(e) NULL)
    if (is.null(res) || is.null(res$summary)) return(NULL)
    
    tibble::tibble(
      obs = i,
      feature = colnames(X),
      x_value = as.numeric(X[i, ]),
      shapley = as.numeric(res$summary[, "mean"])
    )
  })
  
  if (nrow(shapley_list) == 0) stop("No valid Shapley values could be computed.")
  
  if (!is.null(feature)) {
    if (is.numeric(feature)) feature <- colnames(X)[feature]
    shapley_list <- dplyr::filter(shapley_list, feature == !!feature)
  }
  
  p <- ggplot(shapley_list, aes(x = x_value, y = shapley)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "red", linetype = "dashed") +
    labs(
      title = "Shapley Value vs. Covariate",
      x = "Covariate Value",
      y = "Shapley Value (ϕ_ij)"
    ) +
    theme_minimal()
  
  if (is.null(feature)) {
    p <- p + facet_wrap(~feature, scales = "free_x")
  }
  
  return(p)
}


######################################################################## 
#              PROBIT ##################################################
########################################################################

library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1234, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = T, 
  global_shrinkage = T, unlink = T, propensity_seperate = F, gibbs = T, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = F, interaction_rule = "continuous_or_binary",standardize_cov = F
)
#'all'
#
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 10, RCT = FALSE, z_diff = T)

nbcf_fit <- bcf_linear_probit(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 4000, 
  general_params = general_params_default
)


###################### BCF ##########################################################################

 
###########################################################################

library(coda)
chain_beta <- as.mcmc(nbcf_fit$Beta_int[1,,]*sd(data$y))
traceplot(chain_beta)
summary(chain_beta)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")

library(stochtree)
data <- simulate_logit_direct(n = 1000, mu_is_nonlinear = FALSE, f_is_heterogeneous =  T)
fit <- linked_shrinkage_logistic_gibbs_R(y_vec = as.numeric(data$Y), X_mat = as.matrix(sapply(data[, c(1:6)], as.numeric)), Z_vec = as.numeric(data$T), n_iter = 5000, burn_in = 1000)


######################## LINEAR-LINEAR #################################
library(stochtree)
source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")
data <- generate_data_2(500, is_te_hetero = F, is_mu_nonlinear = F, seed = 40, RCT = FALSE, z_diff = F, contrast_binary = T)
X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
y <- as.numeric(data$y)
z <- as.numeric(data$z)
tau_model <- list(gibbs = TRUE, global_shrink = TRUE, unlink = TRUE, p_int = 14)
mu_model <- list(gibbs = TRUE, global_shrink = TRUE, unlink = TRUE, p_int = 14)
res <- horseshoe_sampler(X = X, y = y, z = z, num_iterations = 2000, burn_in = 0, tau_model = tau_model, mu_model = mu_model, propensity_separate = T, standardize_cov = F, interaction_rule = 'continuous_or_binary')
