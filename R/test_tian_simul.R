# Load required package
library(ggplot2)

# --- Data Generation Function ---
generate_data_2 <- function(n = 250,
                            is_te_hetero = FALSE,
                            is_mu_nonlinear = TRUE,
                            seed = 1848,
                            RCT = FALSE,
                            test = FALSE,
                            z_diff = FALSE,
                            tian = FALSE) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rbinom(n, 1, 0.5)
  x5_raw <- sample(1:3, n, replace = TRUE)
  
  g_func <- function(x5) {
    out <- rep(NA, length(x5))
    out[x5 == 1] <- 2
    out[x5 == 2] <- -1
    out[x5 == 3] <- -4
    return(out)
  }
  
  g_x5 <- g_func(x5_raw)
  
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1 * x3
  } else {
    mu <- -6 + g_x5 + 6 * abs(x3 - 1)
  }
  
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    if (test) {
      tau_vec <- 1 + 4 * x1 + 3 * x2 + 2 * x2 * x1
    } else {
      tau_vec <- 1 + 2 * x2 * x4
    }
  }
  
  s <- sd(mu)
  u_i <- runif(n)
  Phi <- function(z) pnorm(z)
  
  if (RCT) {
    pi_x <- rep(0.5, n)
  } else {
    pi_x <- 0.8 * Phi((3 * mu) / s - 0.5 * x1) + 0.05 + (u_i / 10)
  }
  
  pi_x <- pmin(pmax(pi_x, 0), 1)
  z <- rbinom(n, 1, pi_x)
  eps <- rnorm(n)
  
  if (tian) {
    y <- mu - (1 - z) * tau_vec + z * tau_vec + eps
    y_hat <- mu - (1 - z) * tau_vec + z * tau_vec
  } else {
    y <- mu + z * tau_vec + eps
    y_hat <- mu + z * tau_vec
  }
  
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor)[, -1]
  colnames(x5_dev) <- c("x5_1", "x5_2")
  
  if (z_diff) {
    z <- z - 0.5
  }
  
  df <- data.frame(x1, x2, x3, x4, x5_1 = x5_dev[, 1], x5_2 = x5_dev[, 2],
                   z, y, mu, pi_x, tau = tau_vec, y_hat)
  return(df)
}

# --- Call DGP ---
library(Rcpp)
library(ggplot2)
source('R/simul_1.R')
source('R/test_linked_shrinkage.R')
sourceCpp("src/horseshoe_samplers.cpp")

set.seed(42)
n <- 500
data <- generate_data_2(n = n, is_te_hetero = TRUE, is_mu_nonlinear = TRUE,
                        RCT = TRUE, z_diff = TRUE, tian = TRUE)

X <- as.matrix(data[, c("x1", "x2", "x3", "x4", "x5_1", "x5_2")])
y <- 
z <- data$z

# --- Run Tian Model ---
posterior <- sample_linear_part(data$y, rep(1,n), X, intTreat = TRUE, iter = 3000, burnin = 1000, horseshoe =  TRUE)

# --- Construct Interaction Terms ---
p <- ncol(X)
interaction_pairs <- expand.grid(1:p, 1:p)
interaction_pairs <- interaction_pairs[interaction_pairs$Var1 <= interaction_pairs$Var2, ]

X_int <- do.call(
  cbind,
  lapply(seq_len(nrow(interaction_pairs)), function(k) {
    X[, interaction_pairs$Var1[k]] * X[, interaction_pairs$Var2[k]]
  })
)

# --- Posterior CATE Estimation ---
n_samples <- nrow(posterior$beta)
n_obs <- nrow(X)
posterior_tau_hat <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (s in 1:n_samples) {
  alpha_s <- posterior$alpha[s]
  beta_s <- posterior$beta[s, ]
  beta_int_s <- posterior$beta_int[s, ]
  posterior_tau_hat[s, ] <- alpha_s + X %*% beta_s + X_int %*% beta_int_s
}

tau_hat_mean <- colMeans(posterior_tau_hat)

# --- Plot Estimated vs. True CATE ---
df_plot <- data.frame(tau_true = data$tau, tau_hat = tau_hat_mean)

ggplot(df_plot, aes(x = tau_true, y = tau_hat)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    title = "Estimated CATE vs True CATE (Tian Adjusted Covariates)",
    x = "True Treatment Effect",
    y = "Estimated Posterior Mean (CATE)"
  ) +
  theme_minimal()

# --- Diagnostics (optional) ---
cat("Correlation between true and estimated tau:", cor(df_plot$tau_true, df_plot$tau_hat), "\n")
cat("RMSE:", sqrt(mean((df_plot$tau_true - df_plot$tau_hat)^2)), "\n")
