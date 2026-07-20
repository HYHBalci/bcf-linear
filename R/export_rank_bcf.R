# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(stochtree)
library(DescTools)
library(grid)
library(gridExtra)

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
  
  # -- 1. Generate covariates --
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
  
  # -- 2. Prognostic function mu(x) --
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  # -- 3. Treatment effect tau(x) --
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
  
  # -- 4. Compute standard deviation 's' of mu --
  s <- sd(mu)
  
  # -- 5. Propensity function --
  u_i <- runif(n, 0, 1)
  Phi <- function(z) pnorm(z, mean=0, sd=1)
  
  if (RCT) {
    pi_x <- rep(0.5, n)
  } else {
    pi_x <- 0.8 * Phi((3*mu)/s - 0.5*x1) + 0.05 + (u_i / 10)
  }
  
  pi_x <- pmin(pmax(pi_x, 0), 1)
  z <- rbinom(n, size=1, prob=pi_x)
  
  # -- 6. Treatment assignment & Outcome --
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
  
  # -- 7. Formatting Output --
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1] 
  colnames(x5_dev) <- c("x5_1", "x5_2")
  
  if(BCF){
    z <- z_binary
  }
  
  df <- data.frame(
    x1 = x1, 
    x2 = x2, 
    x3 = x3,
    x4 = x4,
    x5_1 = x5_dev[, 1], 
    x5_2 = x5_dev[, 2], 
    z  = z,
    y  = y,
    mu = mu,
    pi_x = pi_x, 
    tau = tau_vec,
    y_hat = y_hat
  )
  
  return(df)
}

# --- Helper Functions  ---
compute_mode <- function(x) {
  if (any(is.na(x))) return(NA)
  d <- density(x)
  d$x[which.max(d$y)]
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper, file_name) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

# Threshold Confidence Metrics
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

# --- Simulation specifications ---
n_simul <- 50
n_values <- c(250, 500, 750, 1000, 1500, 3000) 
heterogeneity_opts <- c(TRUE,FALSE)
linearity_opts <- c(TRUE, FALSE)
scenarios <- c("medical_saturation", "mild_saturation") 

# INITIALIZE the summary table (Expanded with Threshold Columns)
results <- expand.grid(
  n = n_values,
  heterogeneity = heterogeneity_opts,
  linearity = linearity_opts,
  scenario = scenarios
) %>%
  filter(!(heterogeneity == FALSE & scenario != "default")) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA,
         sens_top10 = NA, prec_top10 = NA,
         sens_top25 = NA, prec_top25 = NA,
         sens_top50 = NA, prec_top50 = NA)

# INITIALIZE a list to store the detailed data for plotting
all_scenario_results_list <- list()

# 2. SIMULATION LOOPS
# --------------------------------------------------------------------------
for (scen in scenarios) {
  for (het in heterogeneity_opts) {
    for (lin in linearity_opts) {
      for (n_obser in n_values) {
        
        # 12 columns to include threshold metrics
        result_matrix <- matrix(NA, nrow = n_simul, ncol = 12)
        colnames(result_matrix) <- c("ATE_RMSE", "ATE_Coverage", "ATE_Length",
                                     "CATE_RMSE", "CATE_Coverage", "CATE_Length",
                                     "Sens_Top10", "Prec_Top10", 
                                     "Sens_Top25", "Prec_Top25", 
                                     "Sens_Top50", "Prec_Top50")
        
        cat(sprintf("\nProcessing: Scen=%s | n=%d | Het=%s | Lin=%s\n", 
                    scen, n_obser, het, lin))
        
        for (i in 1:n_simul) {
          # Make sure path matches where your BCF standard outputs are stored
          file_name <- sprintf(
            "E:/tout bcf/Block_link_fit_scen_%s_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
            scen, ifelse(het, "T", "F"), ifelse(lin, "T", "F"), n_obser, i
          )
          
          tryCatch({
            if (!file.exists(file_name)) {
              cat("  [Missing] File not found:", file_name, "\n")
              next 
            } else {
              load(file_name) 
            }
            
            # Generate True Data with BCF logic enabled
            data <- generate_data_medical(
              n = n_obser, is_te_hetero = het, is_mu_nonlinear = lin, 
              seed = i, RCT = FALSE, scenario = scen, z_diff = 0.5, BCF = TRUE, sigma_sq = 1
            )
            
            true_cate <- data$tau
            true_ate <- mean(true_cate)
            
            # --- Extract Standard BCF Posteriors ---
            tau_posterior <- nbcf_fit$tau_hat_train
            
            # Calculate modes and bounds
            tau_mode <- apply(tau_posterior, 1, compute_mode)
            ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
            ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
            
            ate_draws <- colMeans(tau_posterior)
            est_ate <- mean(ate_draws)
            ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
            
            # Core Metrics
            ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2], file_name)
            cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper, file_name)
            
            # --- Threshold Ranking Logic ---
            thresh_10 <- quantile(true_cate, probs = 0.90, na.rm = TRUE)
            thresh_25 <- quantile(true_cate, probs = 0.75, na.rm = TRUE)
            thresh_50 <- quantile(true_cate, probs = 0.50, na.rm = TRUE)
            
            m10 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_10)
            m25 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_25)
            m50 <- compute_threshold_metrics(true_cate, ci_tau_lower, thresh_50)
            
            result_matrix[i, ] <- c(ate_vec, cate_vec, 
                                    m10["sensitivity"], m10["precision"],
                                    m25["sensitivity"], m25["precision"],
                                    m50["sensitivity"], m50["precision"])
            
          }, error = function(e) {
            cat("Error processing file:", file_name, "-", conditionMessage(e), "\n")
          })
        }
        
        # --- UPDATE THE SUMMARY TABLE ---
        idx <- which(results$n == n_obser & 
                       results$heterogeneity == het & 
                       results$linearity == lin & 
                       results$scenario == scen)
        
        if(length(idx) > 0) {
          results[idx, 5:16] <- colMeans(result_matrix, na.rm = TRUE)
        }
        
        # --- STORE DETAILED RESULTS FOR PLOTTING ---
        scenario_df <- as.data.frame(result_matrix) %>%
          mutate(
            scenario = scen,
            n = factor(n_obser), 
            heterogeneity = ifelse(het, "Heterogeneous", "Homogeneous"),
            linearity = ifelse(lin, "Linear", "Nonlinear")
          )
        all_scenario_results_list[[length(all_scenario_results_list) + 1]] <- scenario_df
      }
    }
  }
}

# 3. FINAL OUTPUTS
# --------------------------------------------------------------------------
save(results, file = 'E:/tout bcf/results_linked_analysis_ranking_BCF.RData')
cat("\n\n--- Summary Results Table Saved ---\n")

cat("\n\n--- Generating Boxplots ---\n")

all_results_df <- bind_rows(all_scenario_results_list)

# Subset just the standard metrics for the initial plot (filters out threshold metrics cleanly)
results_long <- all_results_df %>%
  na.omit() %>%
  pivot_longer(
    cols = c(ends_with("RMSE") , ends_with("Length")),
    names_to = c("effect_type", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("RMSE", "Length")))

metrics_plot <- ggplot(results_long, aes(x = n, y = value, fill = linearity)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  facet_grid(metric ~ scenario + heterogeneity, scales = "free_y") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribution of Performance Metrics Across Medical Scenarios (Standard BCF)",
    x = "Sample Size (n)",
    y = "Metric Value",
    fill = "Functional Form (Mu)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

print(metrics_plot)