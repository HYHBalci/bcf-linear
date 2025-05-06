library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")
data <- generate_data_2(500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 18, RCT = F, z_diff = T)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = TRUE, global_shrinkage = F, unlink = F, propensity_seperate = F
)

nbcf_fit <- bcf_linear( X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
                        y_train = as.numeric(data$y),
                        Z_train = as.numeric(data$z), num_gfr = 25, num_burnin =  1000, num_mcmc = 4000, general_params = general_params_default )

load()
library(coda)
# chain_1 <- as.mcmc(nbcf_fit$alpha[1,]*sd(data$y))
# traceplot(chain_1)
# effectiveSize(chain_1)
# effectiveSize(chain_1)
# hist(nbcf_fit$Beta_int[1,,9]*sd(data$y))
# hist(nbcf_fit$sigma2_samples)


# Assuming this is your MCMC sample:
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

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = T, 
  global_shrinkage = FALSE, unlink = FALSE, propensity_seperate = F, step_out = 0.5, max_steps = 50
)
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 1848, RCT = FALSE, z_diff = T)
prognostic_forest_params_default <- list(
  num_trees = 250, alpha = 0.95, beta = 2, 
  min_samples_leaf = 5, max_depth = 10, 
  sample_sigma2_leaf = TRUE, sigma2_leaf_init = NULL, 
  sigma2_leaf_shape = 10, sigma2_leaf_scale = 0.1 * var(data$y) / 250, 
  keep_vars = NULL, drop_vars = NULL
)
nbcf_fit <- bcf_linear(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 7000, 
  general_params = general_params_default,
  prognostic_forest_params = prognostic_forest_params_default
)

set.seed(1848)
mu_hat_train <- nbcf_fit$mu_hat_train
mu_hat_real <- data$y_hat

library(coda)
chain_1 <- as.mcmc(nbcf_fit$Beta_int[2,,])
effectiveSize(chain_1)
traceplot(chain_1)
for(i in 1:25){
hist(mu_hat_train[i,])
abline(v = mu_hat_real[i],          # x-position (a single number or vector)
       col = "red",          # line colour
       lwd = 2,              # line width
       lty = 2) 
}