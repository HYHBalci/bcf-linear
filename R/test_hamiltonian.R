library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/linked_shrinkage_ll_lg.cpp")

set.seed(42)

# Number of samples and predictors
n <- 1000  # Sample size
p <- 3   # Number of main effects

# Generate design matrix X (without intercept)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# True coefficients
alpha_true <- 1.0                  
beta_true <- c(2.0, -1.5, 0.5)      
beta_inter_true <- c(1.2, -0.8, 0.3)  # Interactions: X1X2, X1X3, X2X3

# Generate interaction terms
interaction_terms <- cbind(X[,1] * X[,2], X[,1] * X[,3], X[,2] * X[,3])

# True sigma
sigma_true <- 1.0  

# Generate response variable
y <- alpha_true + X %*% beta_true + interaction_terms %*% beta_inter_true + rnorm(n, sd = sigma_true)

# Combine main effects and interaction terms
X_full <- cbind(X, interaction_terms)

# Get parameter dimensions
pC2 <- p * (p - 1) / 2

leapfrog <- function(param, momentum, step_size, num_steps, X, y) {
  # First half-step for momentum
  momentum <- momentum + 0.5 * step_size * grad_log_posterior_linked_shrinkage(param, X, y)
  
  # Full Leapfrog steps
  for (i in 1:num_steps) {
    # Full position step
    param <- param + step_size * momentum
    if (i != num_steps) {
      # Full momentum step (except for last)
      momentum <- momentum + step_size * grad_log_posterior_linked_shrinkage(param, X, y)
    }
  }
  
  # Final half-step for momentum
  momentum <- momentum + 0.5 * step_size * grad_log_posterior_linked_shrinkage(param, X, y)
  
  return(list(param = param, momentum = momentum))
}

hmc_sampler <- function(init_param,n_warmup, n_samples, step_size, num_steps, X, y) {
  samples <- matrix(NA, nrow = n_samples, ncol = length(init_param))
  param <- init_param
  print(param)
  accept <- 0 
  for (i in 1:n_samples) {
    # Sample random momentum
    momentum <- rnorm(length(init_param), mean = 0, sd = 1)
    
    # Compute initial Hamiltonian
    H_old <- -log_posterior_linked_shrinkage(param, X, y) + sum(momentum^2) / 2
    
    # Perform Leapfrog integration
    leap <- leapfrog(param, momentum, step_size, num_steps, X, y)
    param_new <- leap$param
    momentum_new <- leap$momentum
    
    # Compute new Hamiltonian
    H_new <- -log_posterior_linked_shrinkage(param_new, X, y) + sum(momentum_new^2) / 2
    
    # Metropolis accept/reject
    print(H_old)
    print(H_new)
    if(is.nan(H_new)){
      print(param_new)
    }
    if (runif(1) < exp(H_old - H_new)) {
      param <- param_new  # Accept move
      accept <- accept + 1
    }
    
    samples[i, ] <- param
  }
  
  return(list(samples = samples, accept_rate = round(accept/n_samples, 2)))
}

# Initial parameter values (starting at zero)
init_param <- c(
  0,                            # α (Intercept)
  rnorm(p, 0, 0.1),             # β_j (Main Effects)
  rnorm(pC2, 0, 0.1),           # β_{jk} (Interaction Effects)
  rep(log(1), p),             # log(τ_j), ensuring τ_j > 0
  log(0.5),                   # logit(τ_int), ensuring τ_int ∈ (0.01,1)
  log(1)                        # log(σ²), ensuring σ² > 0
)

# Run HMC
out <- hmc_sampler(init_param, n_samples = 5000, step_size = 0.005, num_steps = 10, X = X, y = y)
samples <- out$samples
print(out$accept_rate)
# Posterior mean estimates
apply(samples, 2, mean)

par(mfrow = c(3, 3))
for (i in 1:9) {
  hist(samples[, i], breaks = 30, main = paste("Posterior: Param", i))
}

library(coda)
sample_mcmc <- as.mcmc(samples[-(1:500), ,drop=FALSE])
traceplot(sample_mcmc)
effectiveSize(sample_mcmc)