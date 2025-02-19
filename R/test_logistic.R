# Load necessary libraries
library(Rcpp)
library(RcppArmadillo)

# Source the compiled C++ function (if needed)
Rcpp::sourceCpp("src/bayesian_logistic_horseshoe.cpp", rebuild = TRUE)
# Generate synthetic data
set.seed(42)
n <- 3000  # Number of observations
p <- 5    # Number of predictors

# Simulate X matrix with random normal values
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Simulate treatment assignment (binary 0 or 1)
z <- rbinom(n, 1, 0.5)

# Generate true coefficients
true_alpha <- 0.5
true_beta_base <- runif(p, -1, 1)  # Main effect coefficients
true_beta_trt <- runif(p, -1, 1)   # Treatment effect coefficients

# Generate linear predictor
lin_pred <- X %*% true_beta_base + z * (true_alpha + X %*% true_beta_trt)

# Convert to probabilities using logistic function
prob <- 1 / (1 + exp(-lin_pred))

# Generate binary response variable y
y <- rbinom(n, 1, prob)

# Run MCMC function
mcmc_results <- run_mcmc_logistic_treat_inter(
  y = y,
  X = X,
  z = z,
  n_iter = 6000,
  burn = 1000,
  thin = 10
)

# Check output
str(mcmc_results)  # Structure of the returned list

# Plot trace plots to check convergence
par(mfrow = c(2, 2))  # Arrange plots in a grid
plot(mcmc_results$alpha, type = "l", main = "Alpha Trace Plot", ylab = "Alpha")
plot(mcmc_results$sigma, type = "l", main = "Sigma Trace Plot", ylab = "Sigma")
plot(mcmc_results$beta_base[,1], type = "l", main = "Beta_base[1] Trace", ylab = "Beta_base[1]")
plot(mcmc_results$beta_trt[,1], type = "l", main = "Beta_trt[1] Trace", ylab = "Beta_trt[1]")

# Summary statistics
summary(mcmc_results$alpha)
summary(mcmc_results$sigma)
summary(mcmc_results$beta_base)
summary(mcmc_results$beta_trt)
