source('R/simul_1.R')
library(stochtree)

n_obser <- 500
het <- TRUE
lin <- FALSE
i <- 1

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper, file_name) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}


general_params_default <- list(
          cutpoint_grid_size = 100, standardize = TRUE, 
          sample_sigma2_global = TRUE, sigma2_global_init = 1, 
          variable_weights = NULL, propensity_covariate = "mu", 
          adaptive_coding = TRUE, control_coding_init = -0.5, 
          treated_coding_init = 0.5, rfx_prior_var = NULL, 
          keep_burnin = FALSE, keep_gfr = TRUE, 
          keep_every = 1, num_chains = 1, verbose = FALSE, 
          probit_outcome_model = FALSE
        )
        
data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = FALSE, z_diff = FALSE, tian = FALSE)
true_cate <- data$tau
true_ate <- mean(true_cate)
        
nbcf_fit <- bcf(
          X_train = data[, c(1:6)],
          y_train = as.numeric(data$y),
          Z_train = as.numeric(data$z), 
          num_gfr = 40, 
          num_burnin = 0, 
          num_mcmc = 500, 
          general_params = general_params_default
        )
  
tau_posterior <- nbcf_fit$tau_hat_train
tau_mean <- apply(tau_posterior, 1, mean)
ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
        
        
compute_metrics(true_cate, tau_mean, ci_tau_lower, ci_tau_upper, file_name)
        
plot(true_cate,tau_mean)
abline(0,1)
        