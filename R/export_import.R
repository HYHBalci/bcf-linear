# 1. SETUP
# --------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(stochtree)
library(DescTools)

generate_data_medical <- function(n = 250,
                                  is_te_hetero = TRUE,
                                  is_mu_nonlinear = TRUE, 
                                  seed = 1848, 
                                  RCT = FALSE, 
                                  # New Argument to toggle specific scenarios
                                  scenario = "default", # Options: "default", "complex_interaction", "medical_saturation"
                                  z_diff = FALSE, 
                                  contrast_binary = TRUE, 
                                  BCF = FALSE, 
                                  sigma_sq = 1) {
  set.seed(seed)
  
  # -- 1. Generate covariates --
  # x1: Continuous (e.g., Biomarker / Age / BMI)
  x1 <- rnorm(n, mean=0, sd=1)
  # x2: Continuous (e.g., Blood Pressure)
  x2 <- rnorm(n, mean=0, sd=1)
  # x3: Continuous (e.g., Lab value)
  x3 <- rnorm(n, mean=0, sd=1)
  # x4: Binary (e.g., Sex or Comorbidity)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
  # x5: Categorical (e.g., Hospital Site or Region)
  x5_raw <- sample(1:3, size=n, replace=TRUE, prob=c(1/3,1/3,1/3))
  
  # Define function for categorical effect
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
      # Previous "test" logic: High complexity polynomial interaction
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
      
    } else if(scenario == "medical_saturation"){
      # NEW SCENARIO: Non-linear Sigmoid Saturation
      # x1 acts as a biomarker. 
      # Patients with low x1 get a small baseline benefit (0.5).
      # As x1 increases, benefit grows rapidly but caps at 4.5 (saturation).
      tau_vec <- 0.5 + 4 / (1 + exp(-2.5 * x1))
      
    } else {
      # Default linear interaction
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
  
  # Handle Z difference/centering if requested
  if(is.logical(z_diff)){
    delta <- if(z_diff) 0.5 else 0
  } else {
    delta <- z_diff
  }
  
  if(z_diff != FALSE){
    z <- z - delta
  }
  
  # Calculate Outcome
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
  d <- density(x)
  mode_result <- optimise(
    f = function(val) {
      -approxfun(d$x, d$y)(val)
    },
    interval = range(x)
  )
  return(mode_result$minimum)
}

compute_metrics <- function(true_values, estimates, ci_lower, ci_upper, file_name) {
  rmse <- sqrt(mean((true_values - estimates)^2, na.rm = TRUE))
  coverage <- mean(true_values >= ci_lower & true_values <= ci_upper, na.rm = TRUE)
  interval_length <- mean(ci_upper - ci_lower, na.rm = TRUE)
  return(c(rmse = rmse, coverage = coverage, interval_length = interval_length))
}

interaction_pairs <- function(num_covariates, boolean_vector) {
  interaction_list <- list()
  if (num_covariates > 1) {
    for (j in 1:(num_covariates - 1)) {
      for (k in (j + 1):num_covariates) {
        if (boolean_vector[j] || boolean_vector[k])
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
      }
    }
  }
  if (length(interaction_list) == 0) return(matrix(nrow=2, ncol=0))
  return(do.call(cbind, interaction_list))
}

# --- Simulation specifications ---
n_simul <- 50
n_values <- c(250, 500, 750, 1000, 1500, 3000) # Updated to match your new grid
heterogeneity_opts <- c(TRUE, FALSE)
linearity_opts <- c(TRUE, FALSE)
scenarios <- c("default", "medical_saturation") # New scenario options
num_chains <- 1 

# --- Determine Interaction Booleans ---
# Generate a dummy dataset to extract covariate properties
dummy_data <- generate_data_medical(n = 500, seed = 1848, scenario = "default")
X_dummy <- as.matrix(dummy_data[, 1:6])
res <- standardize_X_by_index(
  X_initial = X_dummy, 
  process_data = FALSE, 
  interaction_rule = "continuous_or_binary", 
  cat_coding_method = "difference"
)

boolean <- as.logical(as.numeric(res$X_final_var_info$is_binary) + 
                        as.numeric(res$X_final_var_info$is_continuous))

# INITIALIZE the summary table
# We use filter to drop rows that your simulation skipped (Homogeneous + medical_saturation)
results <- expand.grid(
  n = n_values,
  heterogeneity = heterogeneity_opts,
  linearity = linearity_opts,
  scenario = scenarios
) %>%
  filter(!(heterogeneity == FALSE & scenario != "default")) %>%
  mutate(rmse_ate = NA, cover_ate = NA, len_ate = NA,
         rmse_cate = NA, cover_cate = NA, len_cate = NA)

# INITIALIZE a list to store the detailed data for plotting
all_scenario_results_list <- list()


# 2. SIMULATION LOOPS
# --------------------------------------------------------------------------
for (scen in scenarios) {
  for (het in heterogeneity_opts) {
    
    # Optimization match: Skip if TE is Homogeneous and scenario is not default
    if (het == FALSE && scen != "default") next
    
    for (lin in linearity_opts) {
      for (n_obser in n_values) {
        
        result_matrix <- matrix(NA, nrow = n_simul, ncol = 6)
        colnames(result_matrix) <- c("ATE_RMSE", "ATE_Coverage", "ATE_Length",
                                     "CATE_RMSE", "CATE_Coverage", "CATE_Length")
        
        cat(sprintf("\nProcessing: Scen=%s | n=%d | Het=%s | Lin=%s\n", 
                    scen, n_obser, het, lin))
        
        for (i in 1:n_simul) {
          # Updated file path and name to match your new simulation output
          file_name <- sprintf(
            "E:/tout!/Block_link_fit_scen_%s_heter_%s_linear_%s_n_%d_sim_%d.Rdata", 
            scen,
            ifelse(het, "T", "F"), 
            ifelse(lin, "T", "F"), 
            n_obser, 
            i
          )
          
          tryCatch({
            if (!file.exists(file_name)) {
              cat("  [Missing] File not found:", file_name, "\n")
              next 
            } else {
              load(file_name) # loads object: nbcf_fit
            }
            
            # Use the NEW data generation function
            data <- generate_data_medical(
              n = n_obser, 
              is_te_hetero = het, 
              is_mu_nonlinear = lin, 
              seed = i, 
              RCT = FALSE, 
              scenario = scen, 
              z_diff = 0.5, 
              BCF = FALSE,  
              sigma_sq = 1
            )
            
            X <- as.matrix(data[, 1:6])
            true_cate <- data$tau
            true_ate <- mean(true_cate)
            
            alpha_samples <- as.vector(t(nbcf_fit$alpha))
            beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
            beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))
            
            sd_y <- sd(data$y)
            alpha_samples <- alpha_samples * sd_y
            beta_samples <- beta_samples * sd_y
            beta_int_samples <- beta_int_samples * sd_y
            
            ipairs <- interaction_pairs(ncol(X), boolean)
            
            tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                                    nrow = n_obser, byrow = FALSE) + X %*% t(beta_samples)
            
            if (ncol(ipairs) > 0) {
              for (idx in 1:ncol(ipairs)) {
                j <- ipairs[1, idx]
                k <- ipairs[2, idx]
                tau_posterior <- tau_posterior + (X[, j] * X[, k]) %*% t(beta_int_samples[, idx, drop = FALSE])
              }
            }
            
            tau_mode <- apply(tau_posterior, 1, compute_mode)
            ci_tau_lower <- apply(tau_posterior, 1, quantile, probs = 0.025, na.rm = TRUE)
            ci_tau_upper <- apply(tau_posterior, 1, quantile, probs = 0.975, na.rm = TRUE)
            
            ate_draws <- colMeans(tau_posterior)
            est_ate <- mean(ate_draws)
            ci_ate <- quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE)
            
            ate_vec <- compute_metrics(true_ate, est_ate, ci_ate[1], ci_ate[2], file_name)
            cate_vec <- compute_metrics(true_cate, tau_mode, ci_tau_lower, ci_tau_upper, file_name)
            
            result_matrix[i, ] <- c(ate_vec, cate_vec)
            
          }, error = function(e) {
            cat("Error processing file:", file_name, "-", conditionMessage(e), "\n")
          })
        }
        
        # --- UPDATE THE SUMMARY TABLE ---
        ate_results <- result_matrix[, 1:3, drop = FALSE]
        cate_results <- result_matrix[, 4:6, drop = FALSE]
        
        idx <- which(results$n == n_obser & 
                       results$heterogeneity == het & 
                       results$linearity == lin & 
                       results$scenario == scen)
        
        if(length(idx) > 0) {
          results[idx, 5:7] <- colMeans(ate_results, na.rm = TRUE)
          results[idx, 8:10] <- colMeans(cate_results, na.rm = TRUE)
        }
        
        # --- STORE DETAILED RESULTS FOR PLOTTING ---
        scenario_df <- as.data.frame(result_matrix) %>%
          mutate(
            scenario = scen,
            n = factor(n_obser), # Factorized n for better plotting
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
# Updated save path to match your E: drive setup
save(results, file = 'E:/tout!/results_linked_analysis.RData')
cat("\n\n--- Summary Results Table Saved ---\n")


# --- 3B. Create and print the boxplots from detailed results ---
cat("\n\n--- Generating Boxplots ---\n")

all_results_df <- bind_rows(all_scenario_results_list)

results_long <- all_results_df %>%
  na.omit() %>%
  pivot_longer(
    cols = c(ends_with("RMSE") , ends_with("Length")),
    names_to = c("effect_type", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("RMSE", "Length")))

# Updated Plot: Facet by metric AND scenario to handle the new complexity
metrics_plot <- ggplot(results_long, aes(x = n, y = value, fill = linearity)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  # We facet by metric (rows) and scenario + heterogeneity (columns)
  facet_grid(metric ~ scenario + heterogeneity, scales = "free_y") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribution of Performance Metrics Across Medical Scenarios",
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