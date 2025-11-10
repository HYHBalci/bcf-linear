
generate_data <- function(n = 250,
                          is_te_hetero = FALSE,  # toggles homogeneous vs heterogeneous
                          is_mu_nonlinear = TRUE, seed = 1848,RCT = F) {
  set.seed(seed)
  # -- 1. Generate covariates --
  x1 <- rnorm(n, mean=0, sd=1)
  x2 <- rnorm(n, mean=0, sd=1)
  x3 <- rnorm(n, mean=0, sd=1)
  
  # x4 is binary
  x4 <- rbinom(n, size=1, prob=0.5)
  
  # x5 is unordered categorical with 3 levels: 1,2,3 (equal prob)
  # We'll store as integer 1,2,3 but treat it as factor if needed
  x5_raw <- sample(1:3, size=n, replace=TRUE, prob=c(1/3,1/3,1/3))
  
  # Helper function g(x5)
  g_func <- function(x5) {
    # x5 takes values in {1,2,3}
    # g(1)=2, g(2)=-1, g(3)=-4
    out <- rep(NA, length(x5))
    out[x5 == 1] <-  2
    out[x5 == 2] <- -1
    out[x5 == 3] <- -4
    return(out)
  }
  
  # -- 2. Prognostic function mu(x) (linear vs. nonlinear) --
  g_x5  <- g_func(x5_raw)
  if(!is_mu_nonlinear) {
    # linear: mu(x) = 1 + g(x5) + x1*x3
    mu <- 1 + g_x5 + x1*x3
  } else {
    # nonlinear: mu(x) = -6 + g(x5) + 6|x3 - 1|
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  # -- 3. Treatment effect tau(x) (homogeneous vs heterogeneous) --
  if(!is_te_hetero) {
    # homogeneous: tau(x) = 3
    tau_vec <- rep(3, n)
  } else {
    # heterogeneous: tau(x) = 1 + 2*x2*x5
    # x5 is numeric here, but recall it's actually categorical
    tau_vec <- 1 +  2*x2*x1
  }
  
  # -- 4. Compute standard deviation 's' of mu over the sample
  s <- sd(mu)
  
  # -- 5. Propensity function pi(x) = 0.8 Phi(3 mu(x)/s - 0.5 x1) + 0.05 + u_i/10
  #    where u_i ~ Uniform(0,1)
  u_i <- runif(n, 0, 1)
  # standard Normal CDF:
  Phi <- function(z) pnorm(z, mean=0, sd=1)
  
  if(RCT){
    pi_x <- rep(0.5,n)
  } else{
    pi_x <- 0.8 * Phi( (3*mu)/s - 0.5*x1 ) + 0.05 + (u_i / 10)
  }

  # clamp or ensure pi_x in [0,1], just in case
  pi_x <- pmin(pmax(pi_x, 0), 1)
  
  # -- 6. Treatment assignment z ~ Bernoulli(pi_x)
  z <- rbinom(n, size=1, prob=pi_x)
  
  # -- 7. Outcome Y = mu + z*tau + eps, eps ~ N(0,1)
  eps <- rnorm(n, 0, 1)
  y <- mu + z*tau_vec + eps
  y_hat <- mu + z*tau_vec
  # Convert x5 into a factor
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))  
  
  # Create dummy variables
  x5_dummy <- model.matrix(~ x5_factor - 1)  # Removes intercept to get separate dummies
  
  # Combine with the original data
  df <- data.frame(
    x1 = x1, 
    x2 = x2, 
    x3 = x3,
    x4 = x4,
    x5_1 = x5_dummy[, 1],  # Dummy variable for level 1
    x5_2 = x5_dummy[, 2],  # Dummy variable for level 2
    x5_3 = x5_dummy[, 3],  # Dummy variable for level 3
    z  = z,
    y  = y,
    mu = mu,
    pi_x = pi_x, 
    tau = tau_vec
  )
}

generate_data_2 <- function(n = 250,
                            is_te_hetero = FALSE,  
                            is_mu_nonlinear = TRUE, seed = 1848, RCT = FALSE, test = FALSE, z_diff = F, tian = F, contrast_binary = T, BCF = F, sigma_sq = 1) {
  set.seed(seed)
  
  # -- 1. Generate covariates --
  x1 <- rnorm(n, mean=0, sd=1)
  x2 <- rnorm(n, mean=0, sd=1)
  x3 <- rnorm(n, mean=0, sd=1)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
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
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
    # mu <- ifelse(x1 > 0.5, 10*(x3 < 0) + 2*(x3 > 0), 5*(x2 < 0) + 2*(x2 > 0))
  }
  
  # -- 3. Treatment effect tau(x) --
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    if(test){
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
    } else {
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
  # -- 6. Treatment assignment --
  eps <- rnorm(n, 0, sqrt(sigma_sq))
  if(BCF){
    z_binary <- z
  }
  if(z_diff == T){
    delta <- 0.5
  } else if(z_diff == F){
    delta <- 0 
  } else {
    delta <- z_diff
  }
  
  if(z_diff){
    z <- z - delta
  }
 
  y <- mu + z*tau_vec + eps
  y_hat <- mu + z*tau_vec
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1]  # Remove intercept (which is all ones)
  colnames(x5_dev) <- c("x5_1", "x5_2")
  if(BCF){
    z <- z_binary
  }
  df <- data.frame(
    x1 = x1, 
    x2 = x2, 
    x3 = x3,
    x4 = x4,
    x5_1 = x5_dev[, 1],  # Deviance coding for level 1
    x5_2 = x5_dev[, 2],  # Deviance coding for level 2
    z  = z,
    y  = y,
    mu = mu,
    pi_x = pi_x, 
    tau = tau_vec,
    y_hat = y_hat
  )
  
  return(df)
}

generate_data_3 <- function(n = 250,
                            is_te_hetero = FALSE,  # toggles homogeneous vs heterogeneous
                            is_mu_nonlinear = TRUE, seed = 1848, RCT = FALSE, test = FALSE, z_diff = F, tian = F, contrast_binary = T) {
  set.seed(seed)
  
  # -- 1. Generate covariates --
  x1 <- rnorm(n, mean=0, sd=1)
  x2 <- rnorm(n, mean=0, sd=1)
  x3 <- rnorm(n, mean=0, sd=1)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
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
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  # -- 3. Treatment effect tau(x) --
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    if(test){
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
    } else {
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
  # -- 6. Treatment assignment --
  eps <- rnorm(n, 0, 1)
  if(tian){
    y <- mu -(1-z)*tau_vec + z*tau_vec + eps
    y_hat <- mu -(1-z)*tau_vec + z*tau_vec
  } else {
    y <- mu + z*tau_vec + eps
    y_hat <- mu + z*tau_vec
  }
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1]  # Remove intercept (which is all ones)
  
  colnames(x5_dev) <- c("x5_1", "x5_2")
  if(z_diff){
    z <- z - 0.5 
  }
  df <- data.frame(
    x1 = x1, 
    x2 = x2, 
    x3 = x3,
    x4 = x4,
    x5_1 = x5_dev[, 1],  # Deviance coding for level 1
    x5_2 = x5_dev[, 2],  # Deviance coding for level 2
    z  = z,
    y  = y,
    mu = mu,
    pi_x = pi_x, 
    tau = tau_vec,
    y_hat = y_hat
  )
  
  return(df)
}

print_interaction_pairs <- function(p) {
  print("Printing all unique (j, k) pairs where j < k")
  interaction_pairs <- combn(1:p, 2)
  
  # Convert to a data frame for readability
  pairs_df <- data.frame(
    Index = 1:ncol(interaction_pairs),
    Covariate_1 = interaction_pairs[1, ],
    Covariate_2 = interaction_pairs[2, ]
  )
  # Return data frame for further use if needed
  return(pairs_df)
}

simulate_probit_direct <- function(n = 250,
                                   mu_is_nonlinear = TRUE,
                                   f_is_heterogeneous = TRUE,
                                   confound_treatment = FALSE,
                                   seed = 1848) {
  set.seed(seed)
  
  # --- 1. Generate Covariates (fixed set: x1-x3 continuous, x4 binary, x5 categorical) ---
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 0, sd = 1)
  x3 <- rnorm(n, mean = 0, sd = 1)
  x4 <- rbinom(n, size = 1, prob = 0.5) # Binary 0/1
  x5_raw <- sample(1:3, size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3)) # Categorical 1,2,3
  
  # For potential linear terms: sum-coded dummies for x5
  # Reference level will be the last one (level 3)
  x5_factor <- factor(x5_raw, levels = 1:3)
  contrasts(x5_factor) <- contr.sum(3)
  x5_dummies_mat <- model.matrix(~ x5_factor, data = data.frame(x5_factor))[, -1, drop = FALSE] # Drop intercept
  colnames(x5_dummies_mat) <- c("x5_s1", "x5_s2") # Sum-coded: level1 vs mean, level2 vs mean
  x5_s1 <- x5_dummies_mat[, "x5_s1"]
  x5_s2 <- x5_dummies_mat[, "x5_s2"]
  
  # --- 2. Define True Prognostic Function mu(X) ---
  mu_X <- rep(0, n)
  beta_mu_linear_true <- NULL # Store if mu is linear
  X_names_for_mu_linear <- NULL
  
  # Define g_func for categorical effect in non-linear mu (as in generate_data_2)
  g_func_x5 <- function(val) {
    out <- numeric(length(val))
    out[val == 1] <- 1.5
    out[val == 2] <- -0.5
    out[val == 3] <- -1.0
    return(out)
  }
  
  if (mu_is_nonlinear) {
    mu_X <- 0.2 + g_func_x5(x5_raw) + sin(pi * x1 * 0.8) + 1.2 * (x2 * 0.7 - 0.5)^2 - 0.6 * x3 + 0.4 * x4
  } else {
    # Linear mu(X)
    intercept_mu <- 0.1
    beta_mu1 <- 0.7
    beta_mu2 <- -0.5
    beta_mu3 <- 0.45 # for x3
    beta_mu4 <- 0.6 # for x4 (0/1)
    beta_mu_x5_s1 <- 0.5 # for x5_s1 (sum-coded dummy)
    beta_mu_x5_s2 <- -0.243 # for x5_s2 (sum-coded dummy)
    
    mu_X <- intercept_mu + beta_mu1*x1 + beta_mu2*x2 + beta_mu3*x3 + beta_mu4*x4 + 
      beta_mu_x5_s1*x5_s1 + beta_mu_x5_s2*x5_s2
    
    beta_mu_linear_true <- c(intercept_mu, beta_mu1, beta_mu2, beta_mu3, beta_mu4, beta_mu_x5_s1, beta_mu_x5_s2)
    X_names_for_mu_linear <- c("intercept_mu", "x1", "x2", "x3", "x4", "x5_s1", "x5_s2")
  }
  mu_X <- mu_X - mean(mu_X) # Center mu_X
  
  # --- 3. Define True Linear Treatment Effect Modifier f(X) ---
  f_intercept <- 2 # Base treatment effect on the latent scale
  beta_f_true <- c(f_intercept)
  X_names_for_f <- c("intercept_f")
  
  f_X_hetero_terms <- rep(0, n)
  
  if (f_is_heterogeneous) {
    beta_f_x1 <- 5
    beta_f_x4 <- 5
    beta_f_x3 <- 2 
    # Example: f(X) = f_intercept + beta_f_x1*x1 + beta_f_x4*x4
    f_X_hetero_terms <- beta_f_x1 * x1 + beta_f_x4 * x4  + beta_f_x3 * x3 
    
    beta_f_true <- c(beta_f_true, beta_f_x1, beta_f_x4, beta_f_x3)
    X_names_for_f <- c(X_names_for_f, "x1", "x4")
  }
  f_X <- f_intercept + f_X_hetero_terms
  
  # --- 4. Generate Binary Treatment T ---
  propensity_T <- rep(0.5, n) # Default for RCT
  
  if (confound_treatment) {
    # Propensity model inspired by generate_data_2: Phi((3*mu)/s - 0.5*x1) + noise
    # Let's use a simpler linear probit model for propensity for clarity
    eta_T_intercept <- -0.2
    gamma_T_mu <- 0.3
    gamma_T_x1 <- 0.4
    gamma_T_x3 <- -0.25
    
    mu_sd <- sd(mu_X)
    mu_X_scaled_for_prop <- if(mu_sd > 1e-9) (mu_X / mu_sd) else mu_X
    
    eta_T <- eta_T_intercept + gamma_T_mu * mu_X_scaled_for_prop + gamma_T_x1 * x1 + gamma_T_x3 * x3
    propensity_T <- pnorm(eta_T)
    propensity_T <- pmin(pmax(propensity_T, 0.01), 0.99) # Bound propensities
  }
  T_treatment <- rbinom(n, 1, propensity_T)
  
  # --- 5. Calculate Linear Predictor (eta_Y) for Y ---
  linear_predictor_Y <- mu_X + T_treatment * f_X
  
  # --- 6. Calculate Probabilities P(Y=1) ---
  prob_Y_true <- pnorm(linear_predictor_Y)
  
  # --- 7. Generate Binary Outcomes Y ---
  Y_outcome <- rbinom(n, 1, prob_Y_true)
  
  # --- 8. Assemble and Return Data ---
  sim_data <- data.frame(
    Y = Y_outcome,
    T = T_treatment,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4, # Binary 0/1
    x5_s1 = x5_s1, # Sum-coded dummy for categorical
    x5_s2 = x5_s2, # Sum-coded dummy for categorical
    x5_raw = x5_factor, # Original factor for reference or direct use by some models
    mu_X_true = mu_X,
    f_X_true = f_X, 
    propensity_T_true = propensity_T,
    linear_predictor_Y_true = linear_predictor_Y,
    prob_Y_true = prob_Y_true
  )
  
  attr(sim_data, "beta_f_true") <- beta_f_true
  attr(sim_data, "X_names_for_f") <- X_names_for_f
  if (!mu_is_nonlinear) {
    attr(sim_data, "beta_mu_linear_true") <- beta_mu_linear_true
    attr(sim_data, "X_names_for_mu_linear") <- X_names_for_mu_linear
  }
  
  return(sim_data)
}

simulate_logit_direct <- function(n = 250,
                                  mu_is_nonlinear = FALSE,
                                  f_is_heterogeneous = TRUE,
                                  confound_treatment = FALSE,
                                  seed = 1848) {
  set.seed(seed)
  
  # --- 1. Generate Covariates (fixed set: x1-x3 continuous, x4 binary, x5 categorical) ---
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 0, sd = 1)
  x3 <- rnorm(n, mean = 0, sd = 1)
  x4 <- rbinom(n, size = 1, prob = 0.5) # Binary 0/1
  x5_raw <- sample(1:3, size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3)) # Categorical 1,2,3
  
  # For potential linear terms: sum-coded dummies for x5
  x5_factor <- factor(x5_raw, levels = 1:3)
  contrasts(x5_factor) <- contr.sum(3)
  x5_dummies_mat <- model.matrix(~ x5_factor, data = data.frame(x5_factor))[, -1, drop = FALSE] # Drop intercept
  colnames(x5_dummies_mat) <- c("x5_s1", "x5_s2") 
  x5_s1 <- x5_dummies_mat[, "x5_s1"]
  x5_s2 <- x5_dummies_mat[, "x5_s2"]
  
  # --- 2. Define True Prognostic Function mu(X) ---
  mu_X <- rep(0, n)
  beta_mu_linear_true <- NULL 
  X_names_for_mu_linear <- NULL
  
  g_func_x5 <- function(val) {
    out <- numeric(length(val))
    out[val == 1] <- 1.5
    out[val == 2] <- -0.5
    out[val == 3] <- -1.0
    return(out)
  }
  
  if (mu_is_nonlinear) {
    mu_X <- 0.2 + g_func_x5(x5_raw) + sin(pi * x1 * 0.8) + 1.2 * (x2 * 0.7 - 0.5)^2 - 0.6 * x3 + 0.4 * x4
  } else {
    intercept_mu <- 2
    beta_mu1 <- 2
    beta_mu2 <- -3
    beta_mu3 <- 4.5 
    beta_mu4 <- 2.6
    beta_mu_x5_s1 <- 0.5 
    beta_mu_x5_s2 <- -0.2 
    
    mu_X <- intercept_mu + beta_mu1*x1 + beta_mu2*x2  + beta_mu4*x4 + 
      beta_mu_x5_s1*x5_s1 + beta_mu_x5_s2*x5_s2 + 0.25 *x1*x2 + 0.5*x1*x4 + 0.33*x1*x5_s2
    
    beta_mu_linear_true <- c(intercept_mu, beta_mu1, beta_mu2, beta_mu3, beta_mu4, beta_mu_x5_s1, beta_mu_x5_s2)
    X_names_for_mu_linear <- c("intercept_mu", "x1", "x2", "x3", "x4", "x5_s1", "x5_s2")
  }
  mu_X <- mu_X - mean(mu_X) # Center mu_X
  
  # --- 3. Define True Linear Treatment Effect Modifier f(X) ---
  f_intercept <- 2 # Base treatment effect on the latent/linear predictor scale
  beta_f_true <- c(f_intercept)
  X_names_for_f <- c("intercept_f")
  
  f_X_hetero_terms <- rep(0, n)
  
  if (f_is_heterogeneous) {
    beta_f_x1 <- 4
    beta_f_x4 <- 3.5
    
    f_X_hetero_terms <- beta_f_x1 * x1 + beta_f_x4 * x4 
    
    beta_f_true <- c(beta_f_true, beta_f_x1, beta_f_x4)
    X_names_for_f <- c(X_names_for_f, "x1", "x4")
  }
  f_X <- f_intercept + f_X_hetero_terms
  
  # --- 4. Generate Binary Treatment T ---
  propensity_T <- rep(0.5, n) 
  
  if (confound_treatment) {
    eta_T_intercept <- -0.2
    gamma_T_mu <- 0.3
    gamma_T_x1 <- 0.4
    gamma_T_x3 <- -0.25
    
    mu_sd <- sd(mu_X)
    mu_X_scaled_for_prop <- if(mu_sd > 1e-9) (mu_X / mu_sd) else mu_X
    
    eta_T <- eta_T_intercept + gamma_T_mu * mu_X_scaled_for_prop + gamma_T_x1 * x1 + gamma_T_x3 * x3
    # Using plogis for propensity score model for consistency, though pnorm was used before
    propensity_T <- plogis(eta_T) 
    propensity_T <- pmin(pmax(propensity_T, 0.01), 0.99) # Bound propensities
  }
  T_treatment <- rbinom(n, 1, propensity_T)
  
  # --- 5. Calculate Linear Predictor (eta_Y) for Y ---
  linear_predictor_Y <- mu_X + T_treatment * f_X
  prob_Y_true <- plogis(linear_predictor_Y)
  
  # --- 7. Generate Binary Outcomes Y ---
  Y_outcome <- rbinom(n, 1, prob_Y_true)
  
  # --- 8. Assemble and Return Data ---
  sim_data <- data.frame(
    Y = Y_outcome,
    T = T_treatment,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5_s1 = x5_s1,
    x5_s2 = x5_s2,
    x5_raw = x5_factor, 
    mu_X_true = mu_X,
    f_X_true = f_X, 
    propensity_T_true = propensity_T,
    linear_predictor_Y_true = linear_predictor_Y,
    prob_Y_true = prob_Y_true
  )
  
  attr(sim_data, "beta_f_true") <- beta_f_true
  attr(sim_data, "X_names_for_f") <- X_names_for_f
  if (!mu_is_nonlinear) {
    attr(sim_data, "beta_mu_linear_true") <- beta_mu_linear_true
    attr(sim_data, "X_names_for_mu_linear") <- X_names_for_mu_linear
  }
  
  return(sim_data)
}