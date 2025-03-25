library(doParallel)   # For parallel backend
library(foreach)      # For parallel foreach loops
library(dplyr)
library(tidyr)
library(MASS)

source('R/simul_1.R')  # Load `generate_data_2`


# Define model specifications
n_simul <- 50  # Number of simulations
heterogeneity <- c(TRUE, FALSE)  # Homogeneous vs. Heterogeneous effects
linearity <- c(TRUE, FALSE)  # Linear vs. Nonlinear models
n_values <- c(250, 500)  # Sample sizes

# Initialize storage for results
results <- expand.grid(n = n_values, heterogeneity = heterogeneity, linearity = linearity) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

# Function to compute mode from posterior samples
compute_mode <- function(x) {
  d <- density(x)  # Kernel density estimation
  d$x[which.max(d$y)]
}

# Function to compute RMSE, Coverage, and Interval Length
compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper)
  interval_length <- mean(ci_upper - ci_lower)
  return(c(rmse, coverage, interval_length))
}

# Precompute interaction pairs (no squared terms)
interaction_pairs <- function(num_covariates) {
  combn(1:num_covariates, 2)  # Generates all (j,k) pairs
}

# Loop through each model specification
for (n_obser in n_values) {
  for (het in heterogeneity) {
    for (lin in linearity) {
      
      # Initialize storage for simulation results
      ate_results <- matrix(NA, n_simul, 3)
      cate_results <- matrix(NA, n_simul, 3)
      
      for (i in 1:n_simul) {
        file_name <- sprintf("D:/simulations 2/BCF_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
                             ifelse(het, "T", "F"), 
                             ifelse(lin, "T", "F"), 
                             n_obser, i)
        print(file_name)
        
        if (file.exists(file_name)) {
          load(file_name)  # Load saved BCF fit
          
          set.seed(i)  # Ensure consistency across simulations
          data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = TRUE)
          
          X <- as.matrix(data[, c(1:6)])
          z <- data$z
          true_cate <- data$tau  # True CATE from data generation
          
          X_treated <- X[z == 1, , drop = FALSE]  
          true_cate_treated <- true_cate[z == 1]
          true_ate <- mean(true_cate_treated)
          
          # ---- Reconstruct Posterior Tau (CATE) ----
          alpha_samples <- as.vector(t(nbcf_fit$alpha))  # Merge chains
          beta_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta[chain, , ]))
          beta_int_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta_int[chain, , ]))
          
          alpha_samples <- alpha_samples*sd(data$y)
          beta_samples <- beta_samples*sd(data$y)
          beta_int_samples <- beta_int_samples*sd(data$y)
          
          num_samples <- length(alpha_samples)
          num_covariates <- ncol(X)
          interaction_idx <- interaction_pairs(num_covariates)  # Precompute (j,k) pairs
          
          # ---- **Vectorized Matrix Computation for Tau** ----
          tau_posterior <- matrix(NA, nrow = n_obser, ncol = num_samples)
          
          # Compute **all tau samples** in a single matrix operation
          tau_posterior <- matrix(rep(alpha_samples, each = n_obser), nrow = n_obser, byrow = FALSE) +
            X %*% t(beta_samples)  # Apply main effect contributions
          
          # Compute **interaction effects all at once**
          for (idx in 1:ncol(interaction_idx)) {
            j <- interaction_idx[1, idx]
            k <- interaction_idx[2, idx]
            tau_posterior <- tau_posterior + (X[, j] * X[, k]) %*% t(beta_int_samples[, idx])
          }
          
          # ---- Compute Posterior Summary ----
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          
          # Compute **95% credible intervals**
          ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025)
          ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975)
          
          # Compute **ATE from mean of posterior tau**
          est_ate <- mean(tau_mode)
          ci_ate <- quantile(rowMeans(tau_posterior), probs = c(0.025, 0.975))
          
          # Compute **ATE metrics**
          ate_results[i, ] <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          
          # Compute **CATE metrics**
          cate_results[i, ] <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper)
        } else {
          message("File not found: ", file_name)
        }
      }
      
      # ---- Store Results ----
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      results[idx, 4:6] <- colMeans(ate_results, na.rm = TRUE)
      results[idx, 7:9] <- colMeans(cate_results, na.rm = TRUE)
      
    }
  }
}

# ---- Print Final Table ----
results %>%
  mutate(heterogeneity = ifelse(heterogeneity, "Heterogeneous", "Homogeneous"),
         linearity = ifelse(linearity, "Linear", "Nonlinear")) %>%
  arrange(n, heterogeneity, linearity) %>%
  select(n, heterogeneity, linearity, everything()) %>%
  print(n = Inf)
