# Load necessary libraries
library(ggplot2)

# --- Function to calculate risk shift from a model fit ---
# This encapsulates the logic you provided to avoid repeating code.
calculate_risk_shift <- function(fit_object, X_matrix, is_continuous_map) {
  
  # Extract posterior samples from the fit object
  num_chains <- dim(fit_object$Beta)[1]
  alpha_samples <- as.vector(t(fit_object$alpha))
  beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) fit_object$Beta[chain, , ]))
  beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) fit_object$Beta_int[chain, , ]))
  
  # Recreate the interaction pairs based on the provided X matrix
  p_mod <- ncol(X_matrix)
  int_pairs <- list()
  col_idx_counter <- 1
  if (p_mod > 1) {
    for (i in 1:(p_mod - 1)) {
      for (j in (i + 1):p_mod) {
        if (is_continuous_map[i] || is_continuous_map[j]) {
          int_pairs[[col_idx_counter]] <- c(i, j)
          col_idx_counter <- col_idx_counter + 1
        }
      }
    }
  }
  
  n_obser <- nrow(X_matrix)
  
  # --- Calculate posterior of tau(x) on the probit scale ---
  # Start with the intercept (alpha) and main effects (beta)
  tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                          nrow = n_obser, byrow = FALSE) +
    X_matrix %*% t(beta_samples)
  
  # Add the interaction effects
  if (length(int_pairs) > 0 && ncol(beta_int_samples) > 0) {
    # Correctly convert list of pairs to a 2-row matrix for easy iteration
    ipairs_mat <- do.call(cbind, int_pairs)
    
    for (idx in 1:ncol(ipairs_mat)) {
      j <- ipairs_mat[1, idx]
      k <- ipairs_mat[2, idx]
      
      # Ensure the interaction term is a matrix for matrix multiplication
      interaction_term_matrix <- as.matrix(X_matrix[, j] * X_matrix[, k])
      
      # Add the effect for this interaction across all MCMC samples
      tau_posterior <- tau_posterior +
        interaction_term_matrix %*% t(as.matrix(beta_int_samples[, idx]))
    }
  }
  
  # --- Calculate posterior of the risk shift ---
  mu_posterior_samples <- fit_object$mu_hat_train
  tau_posterior_samples <- tau_posterior
  
  # Convert from probit scale to probability scale
  prob_control_posterior <- pnorm(mu_posterior_samples)
  prob_treated_posterior <- pnorm(mu_posterior_samples + tau_posterior_samples)
  
  # Calculate the risk difference
  risk_shift_posterior <- prob_treated_posterior - prob_control_posterior
  
  # Return the posterior mean for each individual
  return(rowMeans(risk_shift_posterior))
}


# --- Main Analysis Script ---

# Load the first model fit (assuming it's in the workspace or loaded from .RData)
load('japanese_fit_linked_glob.RData') 
# Rename the object to avoid overwriting
nbcf_fit_linked <- nbcf_fit 

# Load the second model fit
load('japanese_fit.RData')
# This object is now named nbcf_fit

# --- Prepare the necessary X matrix and metadata ---
# This part needs to be consistent with how the models were originally fit.
# Assuming X_initial was the original data before processing.
# We need to run the standardization/processing step to get the final X matrix.
# This assumes a function `standardize_X_by_index` is available.
# If not, we'll define a placeholder.


data <- read.csv("C:/Users/P094412/Documents/PhD project/Data/data_Japan_clean.csv")
data_ihdp <- read.csv("C:/Users/P094412/Documents/PhD project/Data/ihdp_data.csv")
library(stochtree)
library(dplyr)

Y <- as.vector(data$ponv)
Z <- as.vector(data$dexamethasone_use)
X_initial <- data[, !(names(data) %in% c("ponv", "dexamethasone_use"))]

output <- standardize_X_by_index(X_initial, interaction_rule = "continuous", cat_coding_method = "difference")
X <- output$X_final
is_continuous_map_final_X <- output$X_final_var_info$is_continuous


# --- Calculate Mean Risk Shift for Both Models ---
cat("Calculating risk shift for 'linked shrinkage' model...\n")
risk_shift_linked <- calculate_risk_shift(nbcf_fit_linked, X, is_continuous_map_final_X)

cat("Calculating risk shift for standard BCF model...\n")
risk_shift_normal <- calculate_risk_shift(nbcf_fit, X, is_continuous_map_final_X)

# --- Create Comparison Dataframe and Plot ---
comparison_df <- data.frame(
  Risk_Shift_Standard_BCF = risk_shift_normal,
  Risk_Shift_Linked_Shrinkage = risk_shift_linked
)

cat("\n--- Comparison Plot ---\n")
# Create the scatter plot using ggplot2
comparison_plot <- ggplot(comparison_df, aes(x = Risk_Shift_Standard_BCF, y = Risk_Shift_Linked_Shrinkage)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Comparison of Individual Treatment Effect Estimates",
    subtitle = "Horseshoe vs. Linked Shrinkage Model",
    x = "Posterior Mean Risk Shift (Horseshoe)",
    y = "Posterior Mean Risk Shift (Linked Shrinkage)"
  ) +
  theme_minimal(base_size = 14) +
  coord_fixed(ratio = 1, xlim = range(comparison_df), ylim = range(comparison_df)) # Ensure 1:1 aspect ratio

# Print the plot
print(comparison_plot)

# Print a summary of the correlation
correlation <- cor(comparison_df$Risk_Shift_Standard_BCF, comparison_df$Risk_Shift_Linked_Shrinkage)
cat(sprintf("\nCorrelation between the two sets of estimates: %.4f\n", correlation))

