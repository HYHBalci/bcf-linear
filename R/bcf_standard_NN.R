# Load necessary libraries
source('R/simul_1.R')
library(stochtree)
library(nnet) 

# --- Simulation Parameters ---
n_simul <- 50
heter <- c(TRUE, FALSE)
linear <- c(TRUE, FALSE)
n <- c(250, 500)

# --- Main Simulation Loop ---
for(het in heter){
  for(lin in linear){
    for(n_obser in n){
      for(i in 1:n_simul){
        set.seed(i)
        
        # Define general parameters for the BCF model
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
        
        # Generate simulation data
        data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = FALSE, z_diff = T, tian = T)
        
        # Prepare data for modeling
        X_train_matrix <- as.matrix(sapply(data[, c(1:6)], as.numeric))
        y_train_vector <- as.numeric(data$y)
        Z_train_vector <- as.numeric(data$z)
        
        # --- Propensity Score Estimation (MODIFIED TO OVERFIT) ---

        ps_model_data <- data.frame(z = Z_train_vector, X_train_matrix)
        
        ps_model <- nnet(
          z ~ ., 
          data = ps_model_data,
          size = 10,     # Increased hidden units to make model complex
          maxit = 200,  # Increased iterations
          decay = 0,     # No regularization
          linout = FALSE,# Use logistic activation for binary outcome
          trace = FALSE  # Suppress fitting output
        )
        
        # Predict the propensity scores (probabilities of treatment)
        propensity_scores <- as.numeric(predict(ps_model, type = "raw"))
        
        # --- BCF Model Fitting ---
        nbcf_fit <- bcf(
          X_train = X_train_matrix,
          y_train = y_train_vector,
          Z_train = Z_train_vector,
          propensity_train = propensity_scores, 
          num_gfr = 25,
          num_burnin = 1000,
          num_mcmc = 4000,
          general_params = general_params_default
        )
        
        # --- Save Results ---
        # Generate a dynamic filename based on model settings
        filename <- sprintf("D:/BCF_STANDARD_2_NN/BCF_STANDARD_prop_NN_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
                            ifelse(het, "T", "F"),
                            ifelse(lin, "T", "F"),
                            n_obser,
                            i)
        
        message(paste("Saving results to:", filename))
        save(nbcf_fit, file = filename)
      }
    }
  }
}