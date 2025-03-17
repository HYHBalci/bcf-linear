library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")
data <- generate_data_2(500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 18, RCT = T)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = TRUE
)

nbcf_fit <- bcf_linear( X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
                        y_train = as.numeric(data$y),
                        Z_train = as.numeric(data$z), num_gfr = 25, num_burnin =  1000, num_mcmc = 10000, general_params = general_params_default )

library(coda)
chain_1 <- as.mcmc(nbcf_fit$Beta_int)
effectiveSize(chain_1)
traceplot(chain_1)
hist(nbcf_fit$Beta[1, ,7]*sd(data$y))

beta_2_mcmc_list <- mcmc.list(
  mcmc(nbcf_fit$Beta[1, , 1]),  
  mcmc(nbcf_fit$Beta[2, , 1]) 
)
effectiveSize(nbcf_fit$Beta[1, , ])
traceplot(mcmc(nbcf_fit$Beta_int[2, , 1]*sd(data$y)))
gelman.diag(beta_2_mcmc_list)
# Compute convergence diagnostics