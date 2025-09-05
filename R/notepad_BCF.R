source('R/simul_1.R')
source('R/test_logit_2.R')
source('R/old_linear_linear.R')
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 40, RCT = FALSE, z_diff = 0, contrast_binary = T)
X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
y <- as.numeric(data$y)
z <- as.numeric(data$z)

start.time <- Sys.time()
test_linked <- linear_linear_function(
    y_vec = y, X_mat = X, Z_vec = z,
    family = "gaussian",
    n_iter = 4000, burn_in = 1000,
    prognostic_shrinkage = "horseshoe",
    treatment_shrinkage = "horseshoe",
    standardize_cov = F,
    interaction_rule = "continuous_or_binary",
    num_chains = 1, propensity_as_covariate = TRUE,
    thin = 1, seed = 123, verbose = TRUE, ping = 31, regularize_ate = F
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
summary(test_linked$sigma_sq)
###############################################################
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 1848, RCT = FALSE, z_diff = 0, tian = F)
fit_grouped_hs <- fit_grouped_horseshoes_R(
  y_vec = as.numeric(data$y),
  X_mat = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  Z_vec = as.numeric(data$z),
  family = "gaussian",
  n_iter = 4000, 
  burn_in = 1000,
  num_chains = 2,
  propensity_as_covariate = T,
  method_tau_prognostic = "halfCauchy", tau_prognostic_init = 0.1,
  method_tau_treatment = "halfCauchy", tau_treatment_init = 0.1,
  method_tau_overall = "fixed", tau_overall_init = 1,
  alpha_global_prior_sd = 5.0,
  aleph_prior_sd = 5.0,
  thin = 1,
  seed = i,
  verbose = F,
  ping = 1,
  standardize_cov = F,
  interaction_rule ="continuous_or_binary",
  cat_coding_method = "difference"
)
summary(fit_grouped_hs$gamma)



###############################################################

library(stochtree)
library(coda)
source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")

general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = FALSE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 30, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = T, 
  global_shrinkage = T, unlink = T, propensity_seperate = "none", gibbs = T, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = F, interaction_rule = "continuous_or_binary",standardize_cov = F, simple_prior = F, save_partial_residual = T
)
#'all'
#
data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = F, seed = 45, RCT = T, z_diff = 0, BCF = F,  sigma_sq =1)
data$z <- data$z -0.5
nbcf_fit <- bcf_linear_probit(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 1500, 
  general_params = general_params_default
)

residuals <- rowMeans(nbcf_fit$partial_residuals[,,])

propensity_train <- rowMeans(nbcf_fit$bart_propensity_model$y_hat_train)
mean(residuals)
plot(residuals, propensity_train)

# hist(nbcf_fit$gamma[1,]*sd(data$y))
# load("D:/block_horseshoe/Block_horse_fit_heter_T_linear_F_n_500_sim_2.Rdata")
# data <- generate_data_2(500, is_te_hetero =  T, is_mu_nonlinear = F, seed = 3, RCT = FALSE)
chain_1 <- as.mcmc(nbcf_fit$Beta_int[1,,]*sd(data$y))
chain_2 <- as.mcmc(nbcf_fit$Beta[1,,]*sd(data$y))
chain_3 <- as.mcmc(nbcf_fit$alpha[1,]*sd(data$y))
summary(chain_1)
summary(chain_2)
summary(as.mcmc(nbcf_fit$alpha[1,]*sd(data$y)))
traceplot(chain_1)
traceplot(chain_3)
#################################################################################################################
#FUCK AROUND WITH BCF!!!#########################################################################################
#################################################################################################################

library(stochtree)

source("C:/Users/P094412/Documents/bcf-linear/R/simul_1.R")

data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = T, seed = 1848, RCT = FALSE, z_diff = 0, contrast_binary = T, BCF = T, sigma_sq = 2)
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = NULL, 
  sigma2_global_shape = 0, sigma2_global_scale = 0, 
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = TRUE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = -1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = TRUE, 
  probit_outcome_model = FALSE
)
nbcf_fit <- bcf(
  X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
  y_train = as.numeric(data$y),
  Z_train = as.numeric(data$z), 
  num_gfr = 25, 
  num_burnin = 0, 
  num_mcmc = 500,
  general_params = general_params_default
)

hist(nbcf_fit$sigma2_global_samples)