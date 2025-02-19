# Load required libraries
library(dplyr)
library(tidyr)

# Define model specifications
n_simul <- 50
heter <- c(TRUE, FALSE)
linear <- c(TRUE, FALSE)
n_values <- c(250, 500)

# Initialize empty lists to store results for Beta, Beta_int, and Alpha
summary_list_beta <- list()
summary_list_beta_int <- list()
summary_list_alpha <- list()

# Loop through each model specification
for (het in heter) {
  for (lin in linear) {
    for (n_obser in n_values) {
      # Initialize storage for posterior summaries
      beta_summaries <- list()
      beta_int_summaries <- list()
      alpha_summaries <- list()
      
      for (i in 1:n_simul) {
        # Construct file name
        file_name <- paste0("simulations output/bcf_out_het_", het, "_lin_", lin, "_n_", n_obser, "_sim_", i, ".RData")
        
        # Check if file exists
        if (file.exists(file_name)) {
          load(file_name)  # Load the saved BCF output
          
          # ---- Extract and Summarize Beta ----
          if (!is.null(bcf_out$beta)) {
            posterior_samples <- bcf_out$beta  # Extract beta posterior samples
            
            num_params <- ncol(posterior_samples)  # Ensure the correct orientation
            param_names <- paste0("beta_", 1:num_params)
            
            beta_summary <- apply(posterior_samples, 1, function(param_draws) {
              data.frame(
                mean_value = mean(param_draws),
                sd_value = sd(param_draws),
                lower_ci = quantile(param_draws, 0.025),
                upper_ci = quantile(param_draws, 0.975)
              )
            }) %>%
              bind_rows(.id = "param_index") %>%
              mutate(param_name = param_names[as.numeric(param_index)])  # Assign beta names
            
            beta_summary$sim_id <- i
            beta_summaries[[i]] <- beta_summary
          }
          
          # ---- Extract and Summarize Beta_int ----
          if (!is.null(bcf_out$beta_int)) {
            posterior_samples <- bcf_out$beta_int  # Extract beta_int posterior samples
            
            num_params <- ncol(posterior_samples)
            param_names <- paste0("beta_int_", 1:num_params)
            
            beta_int_summary <- apply(posterior_samples, 1, function(param_draws) {
              data.frame(
                mean_value = mean(param_draws),
                sd_value = sd(param_draws),
                lower_ci = quantile(param_draws, 0.025),
                upper_ci = quantile(param_draws, 0.975)
              )
            }) %>%
              bind_rows(.id = "param_index") %>%
              mutate(param_name = param_names[as.numeric(param_index)])  # Assign beta_int names
            
            beta_int_summary$sim_id <- i
            beta_int_summaries[[i]] <- beta_int_summary
          }
          
          # ---- Extract and Summarize Alpha ----
          if (!is.null(bcf_out$alpha) && length(bcf_out$alpha) > 1) {  # If alpha is a vector
            muy <- bcf_out$muy
            sdy <- bcf_out$sdy
            bcf_out$alpha <- sdy*bcf_out$alpha + muy
            alpha_summary <- data.frame(
              param_name = "alpha",
              mean_value = mean(bcf_out$alpha),
              sd_value = sd(bcf_out$alpha),
              lower_ci = quantile(bcf_out$alpha, 0.025),
              upper_ci = quantile(bcf_out$alpha, 0.975),
              sim_id = i
            )
            alpha_summaries[[i]] <- alpha_summary
          }
        } else {
          message("File not found: ", file_name)
        }
      }
      
      # ---- Combine results across simulations for Beta ----
      if (length(beta_summaries) > 0) {
        summary_df_beta <- bind_rows(beta_summaries) %>%
          group_by(param_name) %>%
          summarise(
            mean_beta_avg = mean(mean_value),
            mean_beta_sd = sd(mean_value),
            overall_sd_beta = mean(sd_value),
            lower_ci_avg = mean(lower_ci),
            upper_ci_avg = mean(upper_ci),
            .groups = "drop"
          ) %>%
          mutate(heterogeneity = het, linearity = lin, sample_size = n_obser)
        
        summary_list_beta[[paste0("het_", het, "_lin_", lin, "_n_", n_obser)]] <- summary_df_beta
      }
      
      # ---- Combine results across simulations for Beta_int ----
      if (length(beta_int_summaries) > 0) {
        summary_df_beta_int <- bind_rows(beta_int_summaries) %>%
          group_by(param_name) %>%
          summarise(
            mean_beta_int_avg = mean(mean_value),
            mean_beta_int_sd = sd(mean_value),
            overall_sd_beta_int = mean(sd_value),
            lower_ci_avg = mean(lower_ci),
            upper_ci_avg = mean(upper_ci),
            .groups = "drop"
          ) %>%
          mutate(heterogeneity = het, linearity = lin, sample_size = n_obser)
        
        summary_list_beta_int[[paste0("het_", het, "_lin_", lin, "_n_", n_obser)]] <- summary_df_beta_int
      }
      
      # ---- Combine results across simulations for Alpha ----
      if (length(alpha_summaries) > 0) {
        summary_df_alpha <- bind_rows(alpha_summaries) %>%
          summarise(
            mean_alpha_avg = mean(mean_value),
            mean_alpha_sd = sd(mean_value),
            overall_sd_alpha = mean(sd_value),
            lower_ci_avg = mean(lower_ci),
            upper_ci_avg = mean(upper_ci),
            .groups = "drop"
          ) %>%
          mutate(param_name = "alpha", heterogeneity = het, linearity = lin, sample_size = n_obser)
        
        summary_list_alpha[[paste0("het_", het, "_lin_", lin, "_n_", n_obser)]] <- summary_df_alpha
      }
    }
  }
}

# ---- Combine results into final data frames ----
final_summary_beta <- bind_rows(summary_list_beta)
final_summary_beta_int <- bind_rows(summary_list_beta_int)
final_summary_alpha <- bind_rows(summary_list_alpha)
# ---- Save summaries to CSV files ----

write.csv(final_summary_beta, "posterior_summary_beta.csv", row.names = FALSE)
write.csv(final_summary_beta_int, "posterior_summary_beta_int.csv", row.names = FALSE)
write.csv(final_summary_alpha, "posterior_summary_alpha.csv", row.names = FALSE)

# ---- Print summary tables ----
print(final_summary_beta, n = 64)
print(final_summary_beta_int, Width = 9, n = 250)
print(final_summary_alpha)
