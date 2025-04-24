library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")
data <- generate_data_2(500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 18, RCT = F)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = TRUE, global_shrinkage = T, unlink = F, propensity_seperate = T
)

nbcf_fit <- bcf_linear( X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
                        y_train = as.numeric(data$y),
                        Z_train = as.numeric(data$z), num_gfr = 25, num_burnin =  1000, num_mcmc = 4000, general_params = general_params_default )

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