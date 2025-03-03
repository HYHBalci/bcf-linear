# Load required libraries
library(dplyr)
library(tidyr)
source('R/simul_1.R')  # Load `generate_data` function

# Define model specifications
n_simul <- 50  # Number of simulations
heterogeneity <- c(TRUE, FALSE)  # Homogeneous vs. Heterogeneous effects
linearity <- c(TRUE, FALSE)  # Linear vs. Nonlinear models
n_values <- c(250, 500)  # Sample sizes

# Initialize storage for results
results <- expand.grid(n = n_values, heterogeneity = heterogeneity, linearity = linearity) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

# Function to compute mode
compute_mode <- function(x) {
  uniq_x <- unique(na.omit(x))
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

# Function to compute RMSE, Coverage, and Interval Length
compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper)
  interval_length <- mean(ci_upper - ci_lower)
  return(c(rmse, coverage, interval_length))
}

# Loop through each model specification
for (n_obser in n_values) {
  for (het in heterogeneity) {
    for (lin in linearity) {
      
      # Initialize storage for simulation results
      ate_results <- matrix(NA, n_simul, 3)
      cate_results <- matrix(NA, n_simul, 3)
      
      for (i in 1:n_simul) {
        file_name <- paste0("simulations output/bcf_out_het_", het, "_lin_", lin, "_n_", n_obser, "_sim_", i, ".RData")
        
        # Check if file exists
        if (file.exists(file_name)) {
          load(file_name)  # Load the saved BCF output
          
          set.seed(i)  # Ensure consistency across simulations
          data <- generate_data(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i)
          
          X <- as.matrix(data[, c(1:7)])
          z <- data$z
          true_cate <- data$tau  # Use true CATE from generate_data function
          
          X_treated <- X[z == 1, , drop = FALSE]  
          true_cate_treated <- true_cate[z == 1]
          true_ate <- mean(true_cate_treated)
          
          # ---- 2️⃣ Extract Posterior Estimates from BCF ----
          if (!is.null(bcf_out$alpha)) {
            muy <- bcf_out$muy
            sdy <- bcf_out$sdy
            bcf_out$alpha <- sdy * bcf_out$alpha + muy
            
            est_ate <- compute_mode(bcf_out$alpha)  # Mode estimate for ATE
            ci_ate <- quantile(bcf_out$alpha, probs = c(0.025, 0.975))  # Bayesian Credible Interval
            
            # Compute ATE metrics
            ate_results[i, ] <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          }
          
          if (!is.null(bcf_out$beta) && !is.null(bcf_out$beta_int)) {
            # ---- 3️⃣ Extract Coefficients Using Mode ----
            beta_estimates <- apply(bcf_out$beta, 2, compute_mode)  # Mode of posterior estimates for main effects
            
            # Construct interaction matrix dynamically
            p <- ncol(X)
            beta_int_matrix <- matrix(0, nrow = p, ncol = p)
            
            interaction_pairs <- t(combn(1:p, 2))  # Get all possible interaction pairs
            
            # Assign interaction effects from `beta_int`
            beta_int_values <- apply(bcf_out$beta_int, 2, compute_mode)
            
            # Ensure `beta_int_values` aligns correctly
            for (idx in 1:nrow(interaction_pairs)) {
              i <- interaction_pairs[idx, 1]
              j <- interaction_pairs[idx, 2]
              if (idx <= length(beta_int_values)) {  # Ensure we don't go out of bounds
                beta_int_matrix[i, j] <- beta_int_values[idx]
                beta_int_matrix[j, i] <- beta_int_values[idx]  # Symmetric interactions
              }
            }
            
            # Compute Estimated CATE using interaction effects properly
            est_cate <- compute_mode(bcf_out$alpha) + X_treated %*% beta_estimates + rowSums((X_treated %*% beta_int_matrix) * X_treated)
            
            # Compute Bayesian Credible Intervals dynamically
            ci_cate_lower <- quantile(bcf_out$alpha, 0.025)
            ci_cate_upper <- quantile(bcf_out$alpha, 0.975)
            
            # Compute CATE metrics for treated observations
            cate_results[i, ] <- compute_metrics(true_cate_treated, est_cate, ci_cate_lower, ci_cate_upper)
          }
        } else {
          message("File not found: ", file_name)
        }
      }
      
      # ---- 4️⃣ Store Results ----
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

