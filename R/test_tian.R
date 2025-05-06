source('R/test_linked_shrinkage.R')
source('R/simul_1.R')

#### TEST (OLD)


set.seed(123)

n <- 500    # Number of observations
p <- 5      # Number of predictors

# 1) Generate design matrix X ~ N(0,1)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)

# 2) Generate *only* strictly upper (or lower) triangular interactions
int_pairs <- expand.grid(1:p, 1:p)
int_pairs <- int_pairs[int_pairs$Var1 < int_pairs$Var2, ]
# Now we have 10 pairs for p=5
X_int <- do.call(
  cbind,
  lapply(seq_len(nrow(int_pairs)), function(k) {
    X[, int_pairs$Var1[k]] * X[, int_pairs$Var2[k]]
  })
)

# 3) True interaction coefficients (10 total)
beta_int_true <- c(0.5, -0.7, 0.3, 0.0, 0.0, -0.2, 0.0, 0.4, -0.5, 0.6)

# 4) True main-effect coefficients
beta_true <- c(1.5, -2.0, 0.8, 0.0, 1.2)  # length = p = 5

# 5) Treatment assignment
z <- rbinom(n, 1, 0.5)

# 6) Intercept, noise
alpha_true <- 2.0
sigma_true <- 1.0

# 7) Generate y
y <- alpha_true +
  z * 1.0 +                       # optional treatment effect
  X %*% beta_true +               # main effects
  X_int %*% beta_int_true +       # interactions
  rnorm(n, 0, sigma_true)

# 9) Run the sampler
result <- sample_linear_part(y, z, X, intTreat = TRUE, iter = 5000, burnin = 2000)

hist(result$tau_int)

### check with coda 

library(coda)

chain_1 <- as.mcmc(result$beta)
effectiveSize(chain_1)


#### TIAN SIMULATION
n_obv <- 500
data <- generate_data_2(n = n_obv, is_te_hetero = T, is_mu_nonlinear = F, RCT = T, z_diff = T, tian = T)

X_tian <-as.matrix(sapply(data[, c(1:6)], as.numeric))
result <- sample_linear_part(data$y, data$z, X_tian, intTreat = TRUE, iter = 5000, burnin = 2000)
result$beta
summary(result$beta)