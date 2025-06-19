library(dplyr)
library(tidyr)
library(MASS)

# Source data generation
source("R/simul_1.R")  # your generate_data_2 function

n_simul <- 50  # number of simulations

# Define model specs
n_values <- c(250, 500)
heterogeneity <- c(TRUE, FALSE)
linearity <- c(TRUE, FALSE)

results <- expand.grid(n = n_values,
                       heterogeneity = heterogeneity,
                       linearity = linearity) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

compute_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper, file_name) {
  rmse <- sqrt(mean((true_values - estimates)^2))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper)
  interval_length <- mean(ci_upper - ci_lower)
  if (rmse > 10) {
    msg <- paste(Sys.time(), "| High RMSE (>10) in", file_name, "\n")
    cat(msg, file = "simulation_errors.txt", append = TRUE)
  }
  return(c(rmse, coverage, interval_length))
}

interaction_pairs <- function(num_covariates) {
  combn(1:num_covariates, 2)
}

# Loop over model settings
for (n_obser in n_values) {
  for (het in heterogeneity) {
    for (lin in linearity) {
      
      result_matrix <- matrix(NA, nrow = n_simul, ncol = 6)
      
      for (i in 1:n_simul) {
        file_name <- sprintf(
          "D:/simulhorseshoe/HORSESHOE_fit_heter_%s_linear_%s_n_%d_sim_%d.Rdata",
          ifelse(het, "T", "F"),
          ifelse(lin, "T", "F"),
          n_obser, i
        )
        
        tryCatch({
          
          if (!file.exists(file_name)) next
          load(file_name)  # loads object: nbcf_fit
          print(file_name)
          set.seed(i)
          data <- generate_data_2(n_obser, het, lin, seed = i, RCT = FALSE)
          
          X <- as.matrix(data[, 1:6])
          z <- data$z
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          
          # Posterior samples
          alpha_samples <- as.vector(t(nbcf_fit$alpha))
          beta_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta[chain, , ]))
          beta_int_samples <- do.call(rbind, lapply(1:2, function(chain) nbcf_fit$Beta_int[chain, , ]))
          
          # Rescale
          sd_y <- sd(data$y)
          alpha_samples <- alpha_samples * sd_y
          beta_samples <- beta_samples * sd_y
          beta_int_samples <- beta_int_samples * sd_y
          
          num_samples <- length(alpha_samples)
          num_covariates <- ncol(X)
          ipairs <- interaction_pairs(num_covariates)
          
          # Posterior samples of tau(x)
          tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                                  nrow = n_obser, byrow = FALSE) +
            X %*% t(beta_samples)
          
          for (idx in 1:ncol(ipairs)) {
            j <- ipairs[1, idx]
            k <- ipairs[2, idx]
            tau_posterior <- tau_posterior +
              (X[, j] * X[, k]) %*% t(beta_int_samples[, idx])
          }
          
          # Posterior summaries
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025)
          ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975)
          
          ate_draws <- colMeans(tau_posterior)
          est_ate <- mean(ate_draws)
          ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975))
          
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2], file_name)
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper, file_name)
          
          result_matrix[i, ] <- c(ate_vec, cate_vec)
          
        }, error = function(e) {
          msg <- paste(Sys.time(), "| ERROR in simulation", file_name, ":", conditionMessage(e), "\n")
          cat(msg, file = "simulation_errors.txt", append = TRUE)
        })
      }
      
      ate_results <- result_matrix[, 1:3, drop = FALSE]
      cate_results <- result_matrix[, 4:6, drop = FALSE]
      
      idx <- which(results$n == n_obser & results$heterogeneity == het & results$linearity == lin)
      results[idx, 4:6] <- colMeans(ate_results, na.rm = TRUE)
      results[idx, 7:9] <- colMeans(cate_results, na.rm = TRUE)
    }
  }
}

# Final output
results %>%
  mutate(heterogeneity = ifelse(heterogeneity, "Heterogeneous", "Homogeneous"),
         linearity = ifelse(linearity, "Linear", "Nonlinear")) %>%
  arrange(n, heterogeneity, linearity) %>%
  select(n, heterogeneity, linearity, everything()) %>%
  print(n = Inf)

results
