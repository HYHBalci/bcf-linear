
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