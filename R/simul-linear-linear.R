source('R/simul_1.R')
source('R/test_logit_2.R')
library(stochtree)

n_simul <- 50
heter <- c(TRUE,FALSE)
linear <-c(TRUE,FALSE)
n  <- c(250, 500)

for(het in heter){
  for(lin in linear){
    for(n_obser in n){
      for(i in 1:n_simul){
        set.seed(i)
        data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = FALSE, z_diff = F, tian = F)
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
          ping = 10000
        )
        
        # Generate a dynamic filename based on model settings
        filename <- sprintf("D:/linear-linear/linear_linear_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
                            ifelse(het, "T", "F"), 
                            ifelse(lin, "T", "F"), 
                            n_obser, 
                            i)
        print(filename)
        save(fit_grouped_hs, file = filename)}
    }
  }
}
