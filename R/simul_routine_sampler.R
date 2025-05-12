source('R/simul_1.R')
library(stochtree)

n_simul <- 50
heter <- c(TRUE, FALSE)
linear <-c(TRUE, FALSE)
n  <- c(250, 500)
index <- 0
step_outs <- c(0.1,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
for(step_out in step_outs){
        set.seed(1848)
        index <- index + 1
        general_params_default <- list(
          cutpoint_grid_size = 100, standardize = TRUE, 
          sample_sigma2_global = TRUE, sigma2_global_init = 1, 
          sigma2_global_shape = 1, sigma2_global_scale = 0.001,
          variable_weights = NULL, propensity_covariate = "mu", 
          adaptive_coding = TRUE, control_coding_init = -0.5, 
          treated_coding_init = 0.5, rfx_prior_var = NULL, 
          random_seed = i, keep_burnin = FALSE, keep_gfr = FALSE, 
          keep_every = 1, num_chains = 2, verbose = F, global_shrinkage = F, unlink = F, 
          propensity_seperate = F, step_out = step_out, max_steps = 50
        )
        data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = F, seed = i, RCT = T)
        
        nbcf_fit <- bcf_linear(
          X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
          y_train = as.numeric(data$y),
          Z_train = as.numeric(data$z), 
          num_gfr = 25, 
          num_burnin = 1000, 
          num_mcmc = 7500, 
          general_params = general_params_default
        )
        
        # Generate a dynamic filename based on model settings
        filename <- sprintf("BCF_fit_step_sim_%d.Rdata", 
                            index)
        print(filename)
        # Save the model with dynamic naming
        save(nbcf_fit, file = filename)
        flush.console()
      }
