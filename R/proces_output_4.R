
library(doSNOW)
library(foreach)
library(dplyr)
library(tidyr)
library(MASS)
source("R/simul_1.R")  # your generate_data_2 function

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

n_simul <- 50  # for example

# Progress bar: range from 0 to n_simul
pb <- txtProgressBar(min = 0, max = n_simul, style = 3)

# A wrapper function that updates the progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Define model specs
n_values <- c(500)
heterogeneity <- c(TRUE)
linearity <- c(TRUE)

results <- expand.grid(n = n_values,
                       heterogeneity = heterogeneity,
                       linearity = linearity) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

compute_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper)
  interval_length <- mean(ci_upper - ci_lower)
  return(c(rmse, coverage, interval_length))
}

interaction_pairs <- function(num_covariates) {
  combn(1:num_covariates, 2)
}

# Loop over model settings
for (n_obser in n_values) {
  for (het in heterogeneity) {
    for (lin in linearity) {
      # We'll store the 6 metrics per simulation in a matrix
      # [n_simul x 6] => (ate_rmse, ate_cover, ate_len, cate_rmse, cate_cover, cate_len)
      result_matrix <-
        foreach(i = 1:n_simul,
                .combine = rbind,
                .packages = c("dplyr","tidyr","MASS"),
                .options.snow = opts) %dopar% {
                  
                  # Update progress bar
                  # (the foreach machinery calls progress automatically,
                  #  so we don't need to call it explicitly here)
                  
                  file_name <- sprintf(
                    "test_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
                    ifelse(het, "T", "F"),
                    ifelse(lin, "T", "F"),
                    n_obser, i
                  )
                  
                  # If file doesn't exist, return NA row
                  if (!file.exists(file_name)) {
                    return(rep(NA, 6))
                  }
                  
                  load(file_name)  # Load BCF fit
                  set.seed(i)
                  
                  data <- generate_data_2(
                    n_obser, 
                    is_te_hetero = het, 
                    is_mu_nonlinear = lin,
                    seed = i, 
                    RCT = F
                  )
                  
                  X <- as.matrix(data[, 1:6])
                  z <- data$z
                  true_cate <- data$tau
                  
                  X_treated <- X[z == 1, , drop = FALSE]
                  true_cate_treated <- true_cate[z == 1]
                  true_ate <- mean(true_cate)
                  
                  # Compute means of x for correct alpha rescaling
                  means_x <- colMeans(X)
                  
                  # Rescale coefficients
                  alpha_samples <- as.vector(t(nbcf_fit$alpha))
                  beta_samples  <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta[chain, , ]))
                  beta_int_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta_int[chain, , ]))
                  
                  # Properly rescale
                  beta_samples  <- beta_samples  * sd(data$y)
                  beta_int_samples <- beta_int_samples * sd(data$y)
                  
                  alpha_samples <- alpha_samples * sd(data$y)
                  num_samples    <- length(alpha_samples)
                  num_covariates <- ncol(X)
                  ipairs <- interaction_pairs(num_covariates)
                  
                  # Vectorized approach
                  tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                                          nrow = n_obser, byrow = FALSE) +
                    X %*% t(beta_samples)
                  
                  for (idx in 1:ncol(ipairs)) {
                    j <- ipairs[1, idx]
                    k <- ipairs[2, idx]
                    tau_posterior <- tau_posterior +
                      (X[, j] * X[, k]) %*% t(beta_int_samples[, idx])
                  }
                  
                  # Posterior summary
                  tau_mode <- apply(tau_posterior, 1, compute_mode)
                  
                  ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025)
                  ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975)
                  
                  # Posterior draws of the ATE
                  ate_draws <- colMeans(tau_posterior)  # colMeans, not rowMeans
                  
                  # Point estimate and 95% CI for the ATE
                  est_ate <- mean(ate_draws)
                  ci_ate  <- quantile(ate_draws, probs = c(0.025, 0.975))
                  
                  ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
                  cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper)
                  
                  # Return a 6-element vector
                  c(ate_vec, cate_vec)
                }
      
      # result_matrix is [n_simul x 6], each row => one iteration
      ate_results  <- result_matrix[, 1:3, drop = FALSE]
      cate_results <- result_matrix[, 4:6, drop = FALSE]
      
      # Store final means
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      results[idx, 4:6] <- colMeans(ate_results, na.rm = TRUE)
      results[idx, 7:9] <- colMeans(cate_results, na.rm = TRUE)
    }
  }
}

# Stop cluster and close progress bar
close(pb)
stopCluster(cl)

# Print final table
results %>%
  mutate(heterogeneity = ifelse(heterogeneity, "Heterogeneous", "Homogeneous"),
         linearity = ifelse(linearity, "Linear", "Nonlinear")) %>%
  arrange(n, heterogeneity, linearity) %>%
  select(n, heterogeneity, linearity, everything()) %>%
  print(n = Inf)
results
`

`