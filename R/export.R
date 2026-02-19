generate_data_medical <- function(n = 250,
                                  is_te_hetero = TRUE,
                                  is_mu_nonlinear = TRUE, 
                                  seed = 1848, 
                                  RCT = FALSE, 
                                  # New Argument to toggle specific scenarios
                                  scenario = "default", # Options: "default", "complex_interaction", "medical_saturation"
                                  z_diff = FALSE, 
                                  contrast_binary = TRUE, 
                                  BCF = FALSE, 
                                  sigma_sq = 1) {
  set.seed(seed)
  
  # -- 1. Generate covariates --
  # x1: Continuous (e.g., Biomarker / Age / BMI)
  x1 <- rnorm(n, mean=0, sd=1)
  # x2: Continuous (e.g., Blood Pressure)
  x2 <- rnorm(n, mean=0, sd=1)
  # x3: Continuous (e.g., Lab value)
  x3 <- rnorm(n, mean=0, sd=1)
  # x4: Binary (e.g., Sex or Comorbidity)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
  # x5: Categorical (e.g., Hospital Site or Region)
  x5_raw <- sample(1:3, size=n, replace=TRUE, prob=c(1/3,1/3,1/3))
  
  # Define function for categorical effect
  g_func <- function(x5) {
    out <- rep(NA, length(x5))
    out[x5 == 1] <-  2
    out[x5 == 2] <- -1
    out[x5 == 3] <- -4
    return(out)
  }
  
  g_x5 <- g_func(x5_raw)
  
  # -- 2. Prognostic function mu(x) --
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  # -- 3. Treatment effect tau(x) --
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    
    if(scenario == "complex_interaction"){ 
      # Previous "test" logic: High complexity polynomial interaction
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
      
    } else if(scenario == "medical_saturation"){
      # NEW SCENARIO: Non-linear Sigmoid Saturation
      # x1 acts as a biomarker. 
      # Patients with low x1 get a small baseline benefit (0.5).
      # As x1 increases, benefit grows rapidly but caps at 4.5 (saturation).
      tau_vec <- 0.5 + 4 / (1 + exp(-2.5 * x1))
      
    } else {
      # Default linear interaction
      tau_vec <- 1 + 2*x2*x4
    }
  }
  
  # -- 4. Compute standard deviation 's' of mu --
  s <- sd(mu)
  
  # -- 5. Propensity function --
  u_i <- runif(n, 0, 1)
  Phi <- function(z) pnorm(z, mean=0, sd=1)
  
  if (RCT) {
    pi_x <- rep(0.5, n)
  } else {
    pi_x <- 0.8 * Phi((3*mu)/s - 0.5*x1) + 0.05 + (u_i / 10)
  }
  
  pi_x <- pmin(pmax(pi_x, 0), 1)
  z <- rbinom(n, size=1, prob=pi_x)
  
  # -- 6. Treatment assignment & Outcome --
  eps <- rnorm(n, 0, sqrt(sigma_sq))
  
  if(BCF){
    z_binary <- z
  }
  
  # Handle Z difference/centering if requested
  if(is.logical(z_diff)){
    delta <- if(z_diff) 0.5 else 0
  } else {
    delta <- z_diff
  }
  
  if(z_diff != FALSE){
    z <- z - delta
  }
  
  # Calculate Outcome
  y <- mu + z*tau_vec + eps
  y_hat <- mu + z*tau_vec
  
  # -- 7. Formatting Output --
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1] 
  colnames(x5_dev) <- c("x5_1", "x5_2")
  
  if(BCF){
    z <- z_binary
  }
  
  df <- data.frame(
    x1 = x1, 
    x2 = x2, 
    x3 = x3,
    x4 = x4,
    x5_1 = x5_dev[, 1], 
    x5_2 = x5_dev[, 2], 
    z  = z,
    y  = y,
    mu = mu,
    pi_x = pi_x, 
    tau = tau_vec,
    y_hat = y_hat
  )
  
  return(df)
}

library(stochtree)

# -- Simulation Settings --
n_simul <- 50
heter   <- c(TRUE, FALSE)
linear  <- c(TRUE, FALSE)
n_vals  <- c(250, 500, 750, 1000, 1500, 3000)
scenarios <- c("default", "medical_saturation") # Add scenarios here
num_chains <- 1 # Ensure this is defined for the params list

# -- Main Loop --
for (scen in scenarios) {
  for (het in heter) {
    
    # Optimization: If TE is Homogeneous (heter=FALSE), the scenario logic 
    # (default vs medical) is irrelevant as tau is fixed at 3. 
    # Skip non-default scenarios to avoid duplicate simulations.
    if (het == FALSE && scen != "default") next
    
    for (lin in linear) {
      for (n_obser in n_vals) {
        for (i in 1:n_simul) {
          
          set.seed(i)
          
          # Define General Params
          general_params_default <- list(
            cutpoint_grid_size = 100, standardize = TRUE, 
            sample_sigma2_global = T, sigma2_global_init = 1, 
            sigma2_global_shape = 1, sigma2_global_scale = 0.001,
            variable_weights = NULL, propensity_covariate = "mu", 
            adaptive_coding = FALSE, control_coding_init = -0.5, 
            treated_coding_init = 0.5, rfx_prior_var = NULL, 
            random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE,   #30
            keep_every = 1, num_chains = num_chains, verbose = T, 
            sample_global_prior = "half-cauchy", unlink = T, propensity_seperate = "none", gibbs = F, step_out = 0.5, max_steps = 150, save_output = F, probit_outcome_model = F, interaction_rule = "continuous_or_binary",standardize_cov = F, simple_prior = T, save_partial_residual = F, regularize_ATE = F,
            sigma_residual = 0, hn_scale = 0, use_ncp = F, n_tijn = 1
          )
          # -- GENERATE DATA --
          # Using the new medical function with scenario argument
          data <- generate_data_medical(
            n = n_obser, 
            is_te_hetero = het, 
            is_mu_nonlinear = lin, 
            seed = i, 
            RCT = FALSE, 
            scenario = scen,  # New Argument
            z_diff = 0.5, 
            BCF = FALSE,  
            sigma_sq = 1
          )
          
          # -- FIT MODEL --
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
          
          # -- SAVE OUTPUT --
          # Filename now includes the scenario to differentiate files
          filename <- sprintf(
            "E:/block_linked/Block_link_fit_scen_%s_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
            scen,
            ifelse(het, "T", "F"), 
            ifelse(lin, "T", "F"), 
            n_obser, 
            i
          )
          
          print(filename)
          save(nbcf_fit, file = filename)
          
        } 
      } # end n
    } # end linear
  } # end heter
} # end scenarios