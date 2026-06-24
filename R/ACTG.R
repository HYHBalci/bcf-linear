# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(speff2trial) # Contains the ACTG175 dataset
library(stochtree)   # BCF and Semi-Parametric BCF methods
source('R/shapley_aux.R')
# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PATCHED PREDICTION
# ==============================================================================

# Custom preprocessing function
standardize_X_by_index_new <- function(X_initial,
                                       process_data = TRUE,
                                       interaction_rule = c("continuous", "continuous_or_binary", "all"),
                                       cat_coding_method = c("sum", "difference")) {
  
  interaction_rule <- match.arg(interaction_rule)
  cat_coding_method <- match.arg(cat_coding_method)
  if (!is.data.frame(X_initial) && !is.matrix(X_initial)) stop("X_initial must be a data.frame or a matrix.")
  
  X_df <- as.data.frame(X_initial)
  colnames(X_df) <- trimws(colnames(X_df))
  original_colnames <- colnames(X_df)
  
  if (process_data) {
    all_cols_info <- lapply(1:ncol(X_df), function(j_idx) {
      col_data <- X_df[[j_idx]]
      valid_data <- col_data[!is.na(col_data)]
      n_unique <- length(unique(valid_data))
      
      is_continuous <- FALSE; is_binary <- FALSE; is_categorical <- FALSE
      processed_col_data <- NULL; should_include <- FALSE
      
      if (n_unique > 1) {
        if (is.logical(col_data)) col_data <- as.character(col_data)
        if (n_unique == 2) {
          factor_col <- factor(col_data)
          processed_col_data <- ifelse(factor_col == levels(factor_col)[1], -1, 1)
          is_binary <- TRUE; should_include <- TRUE
        } else if (is.character(col_data) || is.factor(col_data) || (is.numeric(col_data) && n_unique < 5)) {
          col_as_factor <- factor(col_data)
          contrasts(col_as_factor) <- if (cat_coding_method == "difference") MASS::contr.sdif(nlevels(col_as_factor)) else contr.sum(nlevels(col_as_factor))
          processed_col_data <- col_as_factor
          is_categorical <- TRUE; should_include <- TRUE
        } else if (is.numeric(col_data)) {
          processed_col_data <- col_data
          is_continuous <- TRUE; should_include <- TRUE
        }
      }
      list(processed_data = processed_col_data, is_continuous = is_continuous, 
           is_binary = is_binary, is_categorical = is_categorical, 
           original_idx = j_idx, col_name = original_colnames[j_idx], should_include = should_include)
    })
    
    valid_cols <- Filter(function(x) x$should_include, all_cols_info)
    X_intermediate_list <- lapply(valid_cols, `[[`, "processed_data")
    names(X_intermediate_list) <- sapply(valid_cols, `[[`, "col_name")
    X_intermediate_df <- as.data.frame(X_intermediate_list, check.names = FALSE)
    
    intermediate_map <- data.frame(
      col_name_intermediate = names(X_intermediate_list),
      is_continuous = sapply(valid_cols, `[[`, "is_continuous"),
      is_binary = sapply(valid_cols, `[[`, "is_binary"),
      stringsAsFactors = FALSE
    )
    
    cont_cols <- intermediate_map$col_name_intermediate[intermediate_map$is_continuous]
    X_intermediate_df[cont_cols] <- lapply(X_intermediate_df[cont_cols], function(x) as.vector(scale(x)))
    
    X_final_with_int <- model.matrix(~ ., data = X_intermediate_df)
    X_final <- X_final_with_int[, -1, drop = FALSE]
    attr_assign <- attr(X_final_with_int, "assign")[-1]
    matched_rows <- match(colnames(X_intermediate_df)[attr_assign], intermediate_map$col_name_intermediate)
    
    X_final_var_info <- data.frame(
      is_continuous = intermediate_map$is_continuous[matched_rows],
      is_binary = intermediate_map$is_binary[matched_rows],
      col_name_final = colnames(X_final),
      stringsAsFactors = FALSE
    )
    non_continous_idx_cpp <- which(!X_final_var_info$is_continuous) - 1
  } else {
    X_final <- as.matrix(X_df)
    X_final_var_info <- data.frame(is_continuous = TRUE, is_binary = FALSE, col_name_final = colnames(X_final))
    non_continous_idx_cpp <- numeric(0)
  }
  
  p_int_calculated <- 0
  n_final_cols <- ncol(X_final)
  if (n_final_cols > 1) {
    is_candidate <- switch(interaction_rule,
                           "all" = rep(TRUE, n_final_cols),
                           "continuous_or_binary" = (X_final_var_info$is_continuous | X_final_var_info$is_binary),
                           "continuous" = X_final_var_info$is_continuous)
    is_candidate[is.na(is_candidate)] <- FALSE
    for (i in 1:(n_final_cols - 1)) {
      for (j in (i + 1):n_final_cols) {
        if (is_candidate[i] || is_candidate[j]) p_int_calculated <- p_int_calculated + 1
      }
    }
  }
  return(list(X_final = X_final, p_int = p_int_calculated, 
              non_continous_idx_cpp = non_continous_idx_cpp, X_final_var_info = X_final_var_info))
}

# Patched Predict Function to Extract MCMC Draws from bcf_linear_probit
predict_linear_bcf_patched <- function(object, X, Z, propensity = NULL, rfx_group_ids = NULL, rfx_basis = NULL, ...) {
  if ((!is.data.frame(X)) && (!is.matrix(X))) stop("X must be a matrix or dataframe")
  
  train_set_metadata <- object$train_set_metadata
  X_processed <- preprocessPredictionData(X, train_set_metadata)
  X_linear <- X # Linear component uses raw model matrix
  
  if (!is.null(Z) && is.null(dim(Z))) Z <- as.matrix(as.numeric(Z))
  if (!is.null(propensity) && is.null(dim(propensity))) propensity <- as.matrix(propensity)
  
  if ((object$model_params$propensity_covariate != "none") || (object$model_params$propensity_seperate != "none")) {
    if (is.null(propensity)) {
      if (object$model_params$internal_propensity_model) {
        propensity <- rowMeans(predict(object$bart_propensity_model, X_processed)$y_hat)
        if (is.null(dim(propensity))) propensity <- as.matrix(propensity)
      } else {
        stop("Propensity scores must be provided for this model.")
      }
    }
  }
  
  X_forest <- X_linear
  if (object$model_params$propensity_covariate != "none") X_forest <- cbind(X_linear, propensity)
  forest_dataset_pred <- createForestDataset(X_forest, Z)  
  
  prop_sep <- object$model_params$propensity_seperate
  if (!is.null(prop_sep) && prop_sep == "tau") X_linear <- cbind(X_linear, propensity)
  
  num_chains <- dim(object$Beta)[1]
  num_mcmc <- dim(object$Beta)[2]
  total_samples <- num_chains * num_mcmc
  
  alpha_samples <- as.vector(t(object$alpha))
  beta_samples <- matrix(aperm(object$Beta, c(2, 1, 3)), nrow = total_samples)
  
  has_interactions <- !is.null(object$Beta_int) && (dim(object$Beta_int)[3] > 0)
  if (has_interactions) beta_int_samples <- matrix(aperm(object$Beta_int, c(2, 1, 3)), nrow = total_samples)
  
  # Linear combination reconstructs the CATE parameter
  linear_pred <- matrix(rep(alpha_samples, each=nrow(X)), nrow=nrow(X)) + 
    as.matrix(X_linear) %*% t(beta_samples)
  
  if (has_interactions) {
    int_pairs <- object$interaction_pairs
    num_interactions <- ncol(int_pairs)
    X_int <- matrix(0, nrow = nrow(X_linear), ncol = num_interactions)
    for (k in 1:num_interactions) {
      X_int[, k] <- X_linear[, int_pairs[1, k]] * X_linear[, int_pairs[2, k]]
    }
    linear_pred <- linear_pred + (X_int %*% t(beta_int_samples))
  }
  
  y_std <- object$model_params$outcome_scale
  y_bar <- object$model_params$outcome_mean
  mu_hat <- object$forests_mu$predict(forest_dataset_pred) * y_std + y_bar
  
  # Rescale to original units
  tau_hat_total <- if (object$model_params$standardize) linear_pred * y_std else linear_pred
  
  result <- list(
    mu_hat = mu_hat,
    tau_hat = tau_hat_total, 
    y_hat = mu_hat + (tau_hat_total * as.vector(Z))
  )
  return(result)
}


# ==============================================================================
# 3. PREPARE ACTG175 DATA (Section 4 Setup)
# ==============================================================================
data(ACTG175, package = "speff2trial")

# Subset to: ddI (arm=3, mapped to Control=0) vs ZDV+ddI (arm=1, mapped to Treated=1)
actg_sub <- ACTG175 %>% 
  filter(arms %in% c(1, 3)) %>%
  mutate(A = ifelse(arms == 1, 1, 0)) %>%
  filter(!is.na(cd420))

Y_actg <- actg_sub$cd420
Z_actg <- actg_sub$A

X_actg_raw <- actg_sub %>% 
  dplyr::select(age, wtkg, karnof, cd40, cd80, gender, homo, race, symptom, drugs, hemo, z30)

results_actg <- standardize_X_by_index_new(
  X_initial = X_actg_raw, 
  process_data = TRUE, 
  interaction_rule = "continuous_or_binary", 
  cat_coding_method = "sum"
)
X_train_mat <- results_actg$X_final


# ==============================================================================
# 4. FIT MODELS
# ==============================================================================
cat("\nFitting Models... (This may take a minute)\n")

# A. Semi-Parametric Probit BCF Parameters
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = FALSE, 
  sample_global_prior = "half-cauchy", unlink = TRUE, 
  propensity_seperate = "none", gibbs = TRUE, step_out = 0.5, 
  max_steps = 150, save_output = FALSE, probit_outcome_model = FALSE, 
  interaction_rule = "continuous_or_binary", standardize_cov = FALSE, 
  simple_prior = FALSE, save_partial_residual = FALSE, regularize_ATE = FALSE,
  sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1
)

# Fit Semi-Parametric Model
fit_nbcf <- bcf_linear_probit(
  X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
  general_params = general_params_default
)

# Fit Standard BCF Model
fit_bcf <- bcf(
  X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg, 
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000
)


# ==============================================================================
# 5. INFERENCE (CATEs and ATE)
# ==============================================================================
# Extract MCMC draws. For nbcf, we use the patched predict function passing the train data
tau_draws_bcf <- fit_bcf$tau_hat_train
tau_draws_nbcf <- predict_linear_bcf_patched(fit_nbcf, X = X_train_mat, Z = as.matrix(Z_actg))$tau_hat

# Point Estimates (Mean over MCMC draws per patient)
cate_hat_bcf <- rowMeans(tau_draws_bcf)
cate_hat_nbcf <- rowMeans(tau_draws_nbcf)

# ATE Distributions (Average across patients for each MCMC draw)
ate_draws_bcf <- colMeans(tau_draws_bcf)
ate_draws_nbcf <- colMeans(tau_draws_nbcf)

# Calculate Mean and 95% Credible Intervals for ATE
ate_res_bcf <- c(Mean = mean(ate_draws_bcf), quantile(ate_draws_bcf, probs = c(0.025, 0.975)))
ate_res_nbcf <- c(Mean = mean(ate_draws_nbcf), quantile(ate_draws_nbcf, probs = c(0.025, 0.975)))

cat("\n--- IN-SAMPLE INFERENCE RESULTS ---\n")
cat(sprintf("BCF ATE: %.2f [95%% CI: %.2f, %.2f]\n", ate_res_bcf[1], ate_res_bcf[2], ate_res_bcf[3]))
cat(sprintf("Semi-Parametric ATE: %.2f [95%% CI: %.2f, %.2f]\n", ate_res_nbcf[1], ate_res_nbcf[2], ate_res_nbcf[3]))
cat(sprintf("CATE Correlation (r): %.3f\n", cor(cate_hat_bcf, cate_hat_nbcf)))


# ==============================================================================
# 6. GGPLOT VISUALIZATIONS
# ==============================================================================
plot_data <- data.frame(
  CD4_Baseline = actg_sub$cd40,
  CATE_BCF = cate_hat_bcf,
  CATE_NBCF = cate_hat_nbcf
)

# A. CATE Correlation Plot (with Total Least Squares line)
pca_fit <- prcomp(~ CATE_BCF + CATE_NBCF, data = plot_data)
tls_slope <- pca_fit$rotation["CATE_NBCF", "PC1"] / pca_fit$rotation["CATE_BCF", "PC1"]
tls_intercept <- pca_fit$center["CATE_NBCF"] - tls_slope * pca_fit$center["CATE_BCF"]

p_cor <- ggplot(plot_data, aes(x = CATE_BCF, y = CATE_NBCF)) +
  geom_point(alpha = 0.5, color = "darkblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_abline(slope = tls_slope, intercept = tls_intercept, color = "orange", size = 1.2) +
  labs(title = paste0("CATE Correlation (In-Sample, r = ", round(cor(cate_hat_bcf, cate_hat_nbcf), 3), ")"),
       subtitle = "Orange: Total Least Squares | Red Dashed: 1:1 Line",
       x = "Estimated CATE (Standard BCF)", 
       y = "Estimated CATE (Semi-Parametric BCF)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# B. Density Distributions of Estimated CATEs
density_data <- data.frame(
  CATE = c(plot_data$CATE_BCF, plot_data$CATE_NBCF),
  Model = rep(c("Standard BCF", "Semi-Parametric BCF"), each = nrow(plot_data))
)

p_dist <- ggplot(density_data, aes(x = CATE, fill = Model, color = Model)) +
  geom_density(alpha = 0.4, size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  scale_fill_manual(values = c("tomato", "#008080")) +
  scale_color_manual(values = c("tomato", "#008080")) +
  labs(title = "Distribution of Estimated CATEs",
       x = "Treatment Effect (Difference in CD4 at 20 Weeks)",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

# C. Treatment Effect Heterogeneity by Baseline CD4
p_het <- ggplot(plot_data) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_BCF, color = "Standard BCF"), method = "loess", se = FALSE, size = 1.2) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BCF"), method = "loess", se = FALSE, size = 1.2) +
  geom_point(aes(x = CD4_Baseline, y = CATE_BCF, color = "Standard BCF"), alpha = 0.2) +
  geom_point(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BCF"), alpha = 0.2) +
  scale_color_manual(values = c("Standard BCF" = "tomato", "Semi-Parametric BCF" = "#008080")) +
  labs(title = "Treatment Effect Heterogeneity",
       subtitle = "Estimated CATE vs Baseline CD4 T-cell Count",
       x = "Baseline CD4 Count",
       y = "Estimated CATE",
       color = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

# Print the plots
print(p_cor)
print(p_dist)
print(p_het)

# ==============================================================================
# 7. SHAPLEY ANALYSIS & VISUALIZATION (SEMI-PARAMETRIC PROBIT BCF)
# ==============================================================================

cat("\n--- STARTING SHAPLEY ANALYSIS ---\n")

# Extract the design matrix, posterior betas, and interaction pairs from the fitted model
X_actg <- results_actg$X_final
beta_post <- fit_nbcf$Beta
beta_int_post <- fit_nbcf$Beta_int

# Determine Interaction Pairs 
# Check if the model explicitly tracked pairs; if not, generate them based on your metadata.
if (!is.null(fit_nbcf$interaction_pairs)) {
  ipairs <- fit_nbcf$interaction_pairs
} else {
  # If the model didn't store them, reconstruct using the same rule applied during data prep
  var_info <- results_actg$X_final_var_info
  is_candidate <- var_info$is_continuous | var_info$is_binary # Assuming "continuous_or_binary"
  ipairs <- interaction_pairs(ncol(X_actg), is_candidate)
}

# 1. Compute Shapley values for the entire training cohort
cat("Computing global Shapley values for all individuals...\n")
shapley_results <- compute_shapley_all(
  X = X_actg, 
  beta_post = beta_post, 
  beta_int_post = beta_int_post, 
  ipairs = ipairs
)

# 2. Generate SHAP Summary Plot
cat("Generating SHAP Summary plot...\n")
p_shap_summary <- plot_shapley_contribution(shapley_results)
print(p_shap_summary)

# 3. Generate Feature Importance Breakdown (Main vs. Interaction)
cat("Generating Feature Importance plot...\n")
p_shap_importance <- plot_shapley_importance_breakdown(shapley_results)
print(p_shap_importance)

# 4. Generate Tolerance/Credible Intervals for the top feature
# First, identify the feature with the highest mean absolute Shapley value
top_feature <- shapley_results %>%
  group_by(feature) %>%
  summarize(mean_abs_importance = mean(abs(shapley_total))) %>%
  arrange(desc(mean_abs_importance)) %>%
  slice(1) %>%
  pull(feature)

cat(sprintf("Generating Tolerance Intervals for the top feature: %s...\n", top_feature))

# We'll calculate the tolerance interval for a subset of individuals (e.g., first 500) 
# to keep the computation reasonable for the plot.
plot_indices <- 1:min(500, nrow(X_actg)) 

# Extract the full MCMC posterior matrix for the top feature
post_matrix_top_feat <- get_shapley_posterior_single(
  feature_name = top_feature, 
  indices = plot_indices, 
  X = X_actg, 
  beta_post = beta_post, 
  beta_int_post = beta_int_post, 
  ipairs = ipairs
)

# Generate the Tolerance plot
p_tolerance <- plot_tolerance_intervals(
  post_matrix = post_matrix_top_feat, 
  feature_values = X_actg[plot_indices, top_feature], 
  feature_name = top_feature
)
print(p_tolerance)

cat("\n--- SHAPLEY ANALYSIS COMPLETE ---\n")

# ==============================================================================
# 8. PROFILING THE NEGATIVE CATE TAIL (SEMI-PARAMETRIC BCF)
# ==============================================================================
library(tidyr)
library(ggplot2)
cat("\n--- PROFILING NEGATIVE RESPONDERS ---\n")

# Combine raw demographics with the estimated CATEs
patient_profiles <- actg_sub %>%
  dplyr::select(cd40, cd80, age, wtkg, karnof, gender, homo, race, symptom, drugs, hemo, z30) %>%
  mutate(
    CATE_Standard = cate_hat_bcf,
    CATE_SemiParametric = cate_hat_nbcf,
    # Flag patients who have a predicted negative effect (harm/no benefit)
    Response_Group = ifelse(CATE_SemiParametric < 0, "Negative Effect (CATE < 0)", "Positive Effect (CATE >= 0)")
  )

# 1. See how many patients fall into this tail
cat("Count of patients in each response group:\n")
print(table(patient_profiles$Response_Group))

# 2. Summarize demographic differences between the groups
# (Note: ACTG175 binary coding -> gender: 1=male, homo: 1=yes, race: 1=non-white, z30: 1=prior ZDV)
group_summary <- patient_profiles %>%
  group_by(Response_Group) %>%
  summarise(
    Total_Patients = n(),
    Mean_Baseline_CD4 = mean(cd40, na.rm = TRUE),
    Mean_Baseline_CD8 = mean(cd80, na.rm = TRUE),
    Mean_Age = mean(age, na.rm = TRUE),
    Mean_Weight = mean(wtkg, na.rm = TRUE),
    Pct_Homosexual_Activity = mean(homo == 1, na.rm = TRUE) * 100,
    Pct_Male = mean(gender == 1, na.rm = TRUE) * 100,
    Pct_NonWhite = mean(race == 1, na.rm = TRUE) * 100,
    Pct_Prior_ZDV = mean(z30 == 1, na.rm = TRUE) * 100,
    Pct_Symptomatic = mean(symptom == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  # Pivot for easier side-by-side reading in the console
  pivot_longer(cols = -Response_Group, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Response_Group, values_from = Value)

cat("\nDemographic Comparison (Negative vs Positive Responders):\n")
print(group_summary)

# ==============================================================================
# OPTIONAL: QUICK VISUALIZATION OF THE DIFFERENCES
# ==============================================================================
# Let's plot Baseline CD4 and Weight, since they were top main/interaction drivers

p_box_cd4 <- ggplot(patient_profiles, aes(x = Response_Group, y = cd40, fill = Response_Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("tomato", "#008080")) +
  labs(title = "Baseline CD4 by Response Group", y = "Baseline CD4 Count") +
  theme_minimal() + theme(legend.position = "none")

p_box_wt <- ggplot(patient_profiles, aes(x = Response_Group, y = wtkg, fill = Response_Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("tomato", "#008080")) +
  labs(title = "Patient Weight by Response Group", y = "Weight (kg)") +
  theme_minimal() + theme(legend.position = "none")

grid.arrange(p_box_cd4, p_box_wt, ncol = 2)

# ==============================================================================
# 9. QINI CURVE EVALUATION
# ==============================================================================

cat("\n--- GENERATING QINI CURVES ---\n")

# 1. Define the Qini function (from your uploaded oos_MCMC.R script)
get_qini <- function(df, model_col, model_name) {
  df %>%
    arrange(desc(.data[[model_col]])) %>%
    mutate(
      cum_y_1 = cumsum(Y * Z),
      cum_y_0 = cumsum(Y * (1 - Z)),
      uplift = cum_y_1 - (cum_y_0 * (sum(Z) / sum(1 - Z))),
      frac = row_number() / n(),
      Model = model_name
    ) %>%
    dplyr::select(frac, uplift, Model)
}

# 2. Prepare the evaluation dataframe
# Note: get_qini expects the outcome column to be exactly 'Y' and treatment exactly 'Z'
qini_eval_df <- data.frame(
  Y = Y_actg,               # Outcome: CD4 at 20 weeks
  Z = Z_actg,               # Treatment indicator (1 = ZDV+ddI, 0 = ddI)
  CATE_BCF = cate_hat_bcf,  # Estimates from Standard BCF
  CATE_NBCF = cate_hat_nbcf # Estimates from Semi-Parametric BCF
)

# 3. Calculate the Qini data for both models
qini_bcf  <- get_qini(qini_eval_df, "CATE_BCF", "Standard BCF")
qini_nbcf <- get_qini(qini_eval_df, "CATE_NBCF", "Semi-Parametric BCF")

# Combine for plotting
qini_plot_data <- bind_rows(qini_bcf, qini_nbcf)

# 4. Generate the Plot
p_qini <- ggplot(qini_plot_data, aes(x = frac, y = uplift, color = Model)) +
  geom_line(linewidth = 1.2) +
  # Add a random assignment baseline (straight line from 0 to the final uplift)
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = tail(qini_bcf$uplift, 1)), 
    color = "black", linetype = "dashed", linewidth = 0.8
  ) +
  labs(
    title = "Qini Curve: Standard BCF vs Semi-Parametric BCF",
    subtitle = "Cumulative Gain by Prioritizing Patients with Highest Predicted CATEs",
    x = "Fraction of Population Treated",
    y = "Cumulative Gain (Uplift in CD4 Count)"
  ) +
  scale_color_manual(values = c("Standard BCF" = "tomato", "Semi-Parametric BCF" = "#008080")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Display the plot
print(p_qini)

# 5. Optional: Calculate Area Under the Qini Curve (AUQC)
get_auqc <- function(qini_df) {
  x <- qini_df$frac
  y <- qini_df$uplift
  # Trapezoidal rule for Area Under the Curve
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

auqc_bcf <- get_auqc(qini_bcf)
auqc_nbcf <- get_auqc(qini_nbcf)

cat(sprintf("AUQC - Standard BCF: %.2f\n", auqc_bcf))
cat(sprintf("AUQC - Semi-Parametric BCF: %.2f\n", auqc_nbcf))
# ==============================================================================
# 10. SMOOTHED QINI CURVE (FOR PRESENTATION)
# ==============================================================================

cat("\n--- GENERATING SMOOTHED QINI CURVES ---\n")

# Extract the final uplift value to use in the annotation
final_uplift <- tail(qini_bcf$uplift, 1)

p_qini_smooth <- ggplot(qini_plot_data, aes(x = frac, y = uplift, color = Model)) +
  
  # LOESS smoother 
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
  
  # CLEAN FIX: Use annotate() to draw a single baseline without inheriting the main dataset
  annotate(
    "segment", 
    x = 0, y = 0, 
    xend = 1, yend = final_uplift, 
    color = "black", linetype = "dashed", linewidth = 0.8
  ) +
  
  labs(
    title = "Smoothed Qini Curve: Standard BCF vs Semi-Parametric BCF",
    subtitle = "Smoothed cumulative gain by prioritizing highest predicted CATEs",
    x = "Fraction of Population Treated",
    y = "Cumulative Gain (Uplift in CD4 Count)"
  ) +
  scale_color_manual(values = c("Standard BCF" = "tomato", "Semi-Parametric BCF" = "#008080")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Display the smoothed plot
print(p_qini_smooth)