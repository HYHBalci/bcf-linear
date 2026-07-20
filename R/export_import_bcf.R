# 1. SETUP
# --------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1A. Data Generation Function (Matches Fitting Script Exactly) ---
generate_data_medical <- function(n = 250,
                                  is_te_hetero = TRUE,
                                  is_mu_nonlinear = TRUE, 
                                  seed = 1848, 
                                  RCT = FALSE, 
                                  scenario = "default", 
                                  z_diff = FALSE, 
                                  contrast_binary = TRUE, 
                                  BCF = FALSE, 
                                  sigma_sq = 1) {
  set.seed(seed)
  
  x1 <- rnorm(n, mean=0, sd=1)
  x2 <- rnorm(n, mean=0, sd=1)
  x3 <- rnorm(n, mean=0, sd=1)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
  x5_raw <- sample(1:3, size=n, replace=TRUE, prob=c(1/3,1/3,1/3))
  
  g_func <- function(x5) {
    out <- rep(NA, length(x5))
    out[x5 == 1] <-  2
    out[x5 == 2] <- -1
    out[x5 == 3] <- -4
    return(out)
  }
  g_x5 <- g_func(x5_raw)
  
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    if(scenario == "complex_interaction"){ 
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
      
    } else if(scenario == "medical_saturation"){
      tau_vec <- 0.5 + 4 / (1 + exp(-2.5 * x1))
      
    } else if(scenario == "mild_saturation"){
      tau_vec <- 0.5 + 4 / (1 + exp(-0.6 * x1))
      
    } else {
      tau_vec <- 1 + 2*x2*x4
    }
  }
  
  s <- sd(mu)
  
  u_i <- runif(n, 0, 1)
  Phi <- function(z) pnorm(z, mean=0, sd=1)
  
  if (RCT) {
    pi_x <- rep(0.5, n)
  } else {
    pi_x <- 0.8 * Phi((3*mu)/s - 0.5*x1) + 0.05 + (u_i / 10)
  }
  
  pi_x <- pmin(pmax(pi_x, 0), 1)
  z <- rbinom(n, size=1, prob=pi_x)
  
  eps <- rnorm(n, 0, sqrt(sigma_sq))
  
  if(BCF){
    z_binary <- z
  }
  
  if(is.logical(z_diff)){
    delta <- if(z_diff) 0.5 else 0
  } else {
    delta <- z_diff
  }
  
  if(z_diff != FALSE){
    z <- z - delta
  }
  
  y <- mu + z*tau_vec + eps
  y_hat <- mu + z*tau_vec
  
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1] 
  colnames(x5_dev) <- c("x5_1", "x5_2")
  
  if(BCF){
    z <- z_binary
  }
  
  df <- data.frame(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    x5_1 = x5_dev[, 1], x5_2 = x5_dev[, 2], 
    z  = z, y  = y, mu = mu, pi_x = pi_x, 
    tau = tau_vec, y_hat = y_hat
  )
  
  return(df)
}

# --- 1B. Helper Functions ---
compute_mode <- function(x) {
  if (any(is.na(x))) return(NA)
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

compute_threshold_metrics <- function(true_cate, ci_lower, threshold) {
  true_positive_mask <- true_cate >= threshold
  model_confident_mask <- ci_lower >= threshold
  
  denom_sens <- sum(true_positive_mask, na.rm = TRUE)
  sensitivity <- ifelse(denom_sens > 0, 
                        sum(true_positive_mask & model_confident_mask, na.rm = TRUE) / denom_sens, 
                        NA)
  
  denom_prec <- sum(model_confident_mask, na.rm = TRUE)
  precision <- ifelse(denom_prec > 0, 
                      sum(true_positive_mask & model_confident_mask, na.rm = TRUE) / denom_prec, 
                      NA) 
  
  return(c(sensitivity = sensitivity, precision = precision))
}


# 2. EVALUATION SETUP
# --------------------------------------------------------------------------
n_simul <- 50
n_values <- c(250, 500, 750, 1000, 1500, 3000) 
heterogeneity_opts <- c(TRUE, FALSE)
linearity_opts <- c(TRUE, FALSE)
scenarios <- c("mild_saturation", "medical_saturation", "default") 

# FIXED INITIALIZATION: Properly keeps homogeneous logic intact
results <- expand.grid(
  n = n_values,
  heterogeneity = heterogeneity_opts,
  linearity = linearity_opts,
  scenario = scenarios
) %>%
  # Keep all heterogeneous runs, BUT for homogeneous, only keep the "default" scenario
  filter(heterogeneity == TRUE | (heterogeneity == FALSE & scenario == "default")) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA,
         sens_top10 = NA, prec_top10 = NA,
         sens_top25 = NA, prec_top25 = NA,
         sens_top50 = NA, prec_top50 = NA)


# 3. EVALUATION LOOP
# --------------------------------------------------------------------------
for (scen in scenarios) {
  for (het in heterogeneity_opts) {
    
    # Skip redundant homogeneous scenarios
    if (het == FALSE && scen != "default") next
    
    for (lin in linearity_opts) {
      for (n_obser in n_values) {
        
        result_matrix <- matrix(NA, nrow = n_simul, ncol = 12)
        
        cat(sprintf("Evaluating: Scen=%s | n=%d | Het=%s | Lin=%s\n", 
                    scen, n_obser, het, lin))
        
        for (i in 1:n_simul) {
          
          # Match exactly the saving format from your fitting script
          file_name <- sprintf(
            "E:/tout bcf/Block_link_fit_scen_%s_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
            scen, ifelse(het, "T", "F"), ifelse(lin, "T", "F"), n_obser, i
          )
          
          if (!file.exists(file_name)) {
            next 
          } else {
            load(file_name) # Loads nbcf_fit
          }
          
          # Generate exact matched TRUE data
          data <- generate_data_medical(
            n = n_obser, is_te_hetero = het, is_mu_nonlinear = lin, 
            seed = i, RCT = FALSE, scenario = scen, z_diff = 0.5, BCF = TRUE, sigma_sq = 1
          )
          
          true_cate <- data$tau
          true_ate <- mean(true_cate)
          
          # Extract BCF Posteriors
          tau_posterior <- nbcf_fit$tau_hat_train
          
          tau_mode <- apply(tau_posterior, 1, compute_mode)
          ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
          ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
          
          ate_draws <- colMeans(tau_posterior)
          est_ate <- mean(ate_draws)
          ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
          
          # Calculate Core Metrics
          ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2])
          cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper)
          
          # Calculate Ranking Metrics (Only valid for Heterogeneous models)
          if(het == TRUE) {
            thresh_10 <- quantile(true_cate, probs = 0.90, na.rm = TRUE)
            thresh_25 <- quantile(true_cate, probs = 0.75, na.rm = TRUE)
            thresh_50 <- quantile(true_cate, probs = 0.50, na.rm = TRUE)
            
            m10 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_10)
            m25 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_25)
            m50 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_50)
          } else {
            m10 <- c(sensitivity = NA, precision = NA)
            m25 <- c(sensitivity = NA, precision = NA)
            m50 <- c(sensitivity = NA, precision = NA)
          }
          
          # Store in matrix
          result_matrix[i, ] <- c(ate_vec, cate_vec, 
                                  m10["sensitivity"], m10["precision"],
                                  m25["sensitivity"], m25["precision"],
                                  m50["sensitivity"], m50["precision"])
        }
        
        # Aggregate means and update the final results table
        idx <- which(results$n == n_obser & 
                       results$heterogeneity == het & 
                       results$linearity == lin & 
                       results$scenario == scen)
        
        if(length(idx) > 0) {
          results[idx, 5:16] <- colMeans(result_matrix, na.rm = TRUE)
        }
      }
    }
  }
}

# 4. FINAL OUTPUTS
# --------------------------------------------------------------------------
save(results, file = 'E:/tout bcf/results_BCF_Final_Evaluation.RData')
cat("\n\n--- Evaluation Complete. Results Table Saved ---\n")
print(head(results))