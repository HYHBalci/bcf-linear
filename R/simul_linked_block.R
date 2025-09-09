source('R/simul_1.R')
library(stochtree)

n_simul <- 50
heter <- c(TRUE,FALSE)
linear <-c(TRUE,FALSE)
n  <- c(250,500)

for(het in heter){
  for(lin in linear){
    for(n_obser in n){
        for(i in 1:n_simul){
          set.seed(i)
          general_params_default <- list(
            cutpoint_grid_size = 100, standardize = TRUE, 
            sample_sigma2_global = T, sigma2_global_init = 1, 
            sigma2_global_shape = 1, sigma2_global_scale = 0.001,
            variable_weights = NULL, propensity_covariate = "mu", 
            adaptive_coding = FALSE, control_coding_init = -0.5, 
            treated_coding_init = 0.5, rfx_prior_var = NULL, 
            random_seed = i, keep_burnin = FALSE, keep_gfr = FALSE, 
            keep_every = 1, num_chains = num_chains, verbose = T, 
            global_shrinkage = T, unlink = T, propensity_seperate = "none", gibbs = T, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = F, interaction_rule = "continuous_or_binary",standardize_cov = F, simple_prior = F, save_partial_residual = F, regularize_ATE = F
          )
          data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = F, z_diff = 0.5, BCF = F,  sigma_sq =1)
          nbcf_fit <- bcf_linear_probit(
            X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
            y_train = as.numeric(data$y),
            Z_train = as.numeric(data$z), 
            propensity_train = as.numeric(data$pi_x),
            num_gfr = 25, 
            num_burnin = 1000, 
            num_mcmc = 2000, 
            general_params = general_params_default
          )
          # Generate a dynamic filename based on model settings
          filename <- sprintf("E:/block_linked/Block_link_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
                              ifelse(het, "T", "F"), 
                              ifelse(lin, "T", "F"), 
                              n_obser, 
                              i)
          print(filename)
          save(nbcf_fit, file = filename)}
      }
    }
  }
