library(Rcpp)
library(ggplot2)
source('R/simul_1.R')
source('R/test_linked_shrinkage.R')
sourceCpp("src/horseshoe_samplers.cpp")

set.seed(123)

# Parameters
n <- 500     # number of samples
p <- 5        # number of predictors
alpha_true <- 2.0  # intercept

# Main effect coefficients (true)
beta_true <- c(1.0, -0.5, 0.8, -1.2, 0.3)

# Generate predictors
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Generate all pairwise interaction terms (i < j)
int_pairs <- expand.grid(1:p, 1:p)
int_pairs <- int_pairs[int_pairs$Var1 < int_pairs$Var2, ]
p_int <- nrow(int_pairs)

# Create interaction matrix
X_int <- do.call(
  cbind,
  lapply(1:p_int, function(k) {
    X[, int_pairs$Var1[k]] * X[, int_pairs$Var2[k]]
  })
)

# Interaction coefficients (some non-zero, others zero for sparsity)
beta_int_true <- rep(0, p_int)
nonzero_idx <- sample(1:p_int, 5)  # randomly select 5 non-zero interactions
beta_int_true[nonzero_idx] <- runif(5, -0.7, 0.7)

# Response generation
y <- alpha_true +
  X %*% beta_true +
  X_int %*% beta_int_true +
  rnorm(n, sd = 1.0)

posterior <- sample_linear_part(y, rep(1,n), X, intTreat = TRUE, iter = 5000, burnin = 1000, horseshoe =  TRUE)
