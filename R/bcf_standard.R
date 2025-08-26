source('R/simul_1.R')
library(stochtree)

n_simul <- 50
heter <- c(FALSE)
linear <-c(TRUE,FALSE)
n  <- c(250, 500)

for(het in heter){
  for(lin in linear){
    for(n_obser in n){
      for(i in 1:n_simul){
        set.seed(i)
        general_params_default <- list(
          cutpoint_grid_size = 100, standardize = TRUE, 
          sample_sigma2_global = TRUE, sigma2_global_init = 1, 
          sigma2_global_shape = 1, sigma2_global_scale = 0.001, 
          variable_weights = NULL, propensity_covariate = "mu", 
          adaptive_coding = F, control_coding_init = 0, 
          treated_coding_init = 1, rfx_prior_var = NULL, 
          random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
          keep_every = 1, num_chains = 2, verbose = FALSE, 
          probit_outcome_model = FALSE
        )
        data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = FALSE, z_diff = F, tian = F)
        nbcf_fit <- bcf(
          X_train = as.matrix(sapply(data[, c(1:6)], as.numeric)),
          y_train = as.numeric(data$y),
          Z_train = as.numeric(data$z), 
          num_gfr = 25, 
          num_burnin = 1000, 
          num_mcmc = 4000, 
          general_params = general_params_default
        )
        # Generate a dynamic filename based on model settings
        filename <- sprintf("D:/BCF_STANDARD_2/BCF_STANDARD_2_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
                            ifelse(het, "T", "F"), 
                            ifelse(lin, "T", "F"), 
                            n_obser, 
                            i)
        print(filename)
        save(nbcf_fit, file = filename)}
    }
  }
}
