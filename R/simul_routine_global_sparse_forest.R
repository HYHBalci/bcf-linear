source('R/simul_1.R')
library(stochtree)

n_simul <- 50
heter <- c(TRUE)
linear <- c(TRUE,FALSE)
n <- c(500)

for (het in heter) {
  for (lin in linear) {
    for (n_obser in n) {
      for (i in 1:n_simul) {
        set.seed(i)
        general_params_default <- list(
          cutpoint_grid_size = 100, standardize = TRUE, 
          sample_sigma2_global = TRUE, sigma2_global_init = 1, 
          sigma2_global_shape = 1, sigma2_global_scale = 0.001,
          variable_weights = NULL, propensity_covariate = "mu", 
          adaptive_coding = TRUE, control_coding_init = -0.5, 
          treated_coding_init = 0.5, rfx_prior_var = NULL, 
          random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
          keep_every = 1, num_chains = 2, verbose = F, 
          global_shrinkage = FALSE, unlink = FALSE, propensity_seperate = F, step_out = 0.5, max_steps = 50
        )
        data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = FALSE, z_diff = T)
        
        nbcf_fit <- bcf_linear(
          X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
          y_train = as.numeric(data$y),
          Z_train = as.numeric(data$z), 
          num_gfr = 25, 
          num_burnin = 1000, 
          num_mcmc = 5000, 
          general_params = general_params_default        )
        
        filename <- sprintf("test_heter_%s_linear_%s_n_%d_sim_%d.RData", 
                            ifelse(het, "T", "F"), 
                            ifelse(lin, "T", "F"), 
                            n_obser, 
                            i)
        # Save using qs
        save(nbcf_fit, file = filename)
        print(filename)
      }
    }
  }
}