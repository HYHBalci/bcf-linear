library(dplyr)
library(tidyr)

source("R/simul_1.R")  # for generate_data_2()

n_simul <- 50

# Define design grid
n_values <- c(250, 500)
heterogeneity <- c(TRUE, FALSE)
linearity <- c(TRUE, FALSE)

# Storage for results
results <- expand.grid(n = n_values,
                       heterogeneity = heterogeneity,
                       linearity = linearity) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

# Helper functions
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

# Main evaluation loop
for (n_obser in n_values) {
  for (het in heterogeneity) {
    for (lin in linearity) {
      
      result_matrix <- matrix(NA, nrow = n_simul, ncol = 6)
      
      for (i in 1:n_simul) {
        cat(sprintf("Running sim %d for n = %d, heter = %s, linear = %s\n", i, n_obser, het, lin))
        
        file_name <- sprintf("D:/tian/tian_heter_%s_linear_%s_n_%d_sim_%d.RData",
                             ifelse(het, "T", "F"),
                             ifelse(lin, "T", "F"),
                             n_obser, i)
        
        if (!file.exists(file_name)) next
        
        tryCatch({
          load(file_name)
          if (!exists("posterior")) stop("Posterior object missing in file")
          
          set.seed(i)
          
          data <- generate_data_2(n_obser,
                                  is_te_hetero = het,
                                  is_mu_nonlinear = lin,
                                  seed = i,
                                  RCT = TRUE,
                                  z_diff = TRUE,
                                  tian = TRUE)
          
          X <- as.matrix(data[, c("x1", "x2", "x3", "x4", "x5_1", "x5_2")])
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          n_obs <- nrow(X)
          
          interaction_pairs <- expand.grid(Var1 = 1:ncol(X), Var2 = 1:ncol(X)) %>%
            filter(Var1 <= Var2)
          X_int <- do.call(cbind, lapply(seq_len(nrow(interaction_pairs)), function(k) {
            X[, interaction_pairs$Var1[k]] * X[, interaction_pairs$Var2[k]]
          }))
          
          n_samples <- nrow(posterior$beta)
          tau_post <- matrix(NA, nrow = n_samples, ncol = n_obs)
          
          for (s in 1:n_samples) {
            
            alpha_s <- posterior$alpha[s]
            beta_s <- as.numeric(posterior$beta[s, ])
            beta_int_s <- as.numeric(posterior$beta_int[s, ])
            if (ncol(X) != length(beta_s)) {
              stop(sprintf("Mismatch: ncol(X) = %d, length(beta_s) = %d", ncol(X), length(beta_s)))
            }
            if (ncol(X_int) != length(beta_int_s)) {
              stop(sprintf("Mismatch: ncol(X_int) = %d, length(beta_int_s) = %d", ncol(X_int), length(beta_int_s)))
            }
            tau_post[s, ] <- alpha_s + X %*% beta_s + X_int %*% beta_int_s
          }
          
          tau_mode <- apply(tau_post, 2, compute_mode)
          ci_lower <- apply(tau_post, 2, quantile, probs = 0.025)
          ci_upper <- apply(tau_post, 2, quantile, probs = 0.975)
          
          ate_samples <- rowMeans(tau_post)
          est_ate <- mean(ate_samples)
          ci_ate <- quantile(ate_samples, probs = c(0.025, 0.975))
          
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_lower, ci_upper)
          
          result_matrix[i, ] <- c(ate_vec, cate_vec)
          
        }, error = function(e) {
          cat(sprintf("Simulation %d failed: %s\n", i, e$message))
        })
      }
      
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      ate_results <- result_matrix[, 1:3, drop = FALSE]
      cate_results <- result_matrix[, 4:6, drop = FALSE]
      results[idx, 4:6] <- colMeans(ate_results, na.rm = TRUE)
      results[idx, 7:9] <- colMeans(cate_results, na.rm = TRUE)
    }
  }
}

# Display results
results %>%
  mutate(heterogeneity = ifelse(heterogeneity, "Heterogeneous", "Homogeneous"),
         linearity = ifelse(linearity, "Linear", "Nonlinear")) %>%
  arrange(n, heterogeneity, linearity) %>%
  select(n, heterogeneity, linearity, everything()) %>%
  print(n = Inf)
