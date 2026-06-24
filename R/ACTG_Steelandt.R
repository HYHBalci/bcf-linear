# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(speff2trial) # Contains the ACTG175 dataset
library(stochtree)   # Semi-Parametric BCF method
library(ranger)      # Random Forest for the Vansteelandt DR-Learner
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

predict_linear_bcf_patched <- function(object, X, Z, propensity = NULL, rfx_group_ids = NULL, rfx_basis = NULL, ...) {
  if ((!is.data.frame(X)) && (!is.matrix(X))) stop("X must be a matrix or dataframe")
  
  train_set_metadata <- object$train_set_metadata
  X_processed <- preprocessPredictionData(X, train_set_metadata)
  X_linear <- X 
  
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
# 4. FIT MODELS (Semi-Parametric vs Vansteelandt DR-Learner)
# ==============================================================================
cat("\nFitting Models... (This may take a minute)\n")

# ---------------------------------------------------------
# MODEL A: Semi-Parametric Probit BCF 
# ---------------------------------------------------------
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

fit_nbcf <- bcf_linear_probit(
  X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
  general_params = general_params_default
)

# ---------------------------------------------------------
# MODEL B: Vansteelandt DR-Learner (Algorithm SS-B)
# ---------------------------------------------------------
# Note: ACTG175 is randomized, so true propensity is known (pi = 522/1083 ~ 0.48)
pi_hat <- mean(Z_actg) 
df_dr <- data.frame(Y = Y_actg, A = Z_actg, X_train_mat)

# We use 5-fold cross-fitting to generate the pseudo-outcomes as defined in Section 4
set.seed(123)
K_folds <- 5
fold_ids <- sample(rep(1:K_folds, length.out = nrow(df_dr)))
pseudo_Y <- numeric(nrow(df_dr))

for (k in 1:K_folds) {
  train_idx <- which(fold_ids != k)
  val_idx <- which(fold_ids == k)
  
  # Learn outcome nuisance model E[Y|A,X]
  rf_mu <- ranger(Y ~ ., data = df_dr[train_idx, ])
  
  df_val_1 <- df_dr[val_idx, ]; df_val_1$A <- 1
  df_val_0 <- df_dr[val_idx, ]; df_val_0$A <- 0
  
  mu_1_hat <- predict(rf_mu, data = df_val_1)$predictions
  mu_0_hat <- predict(rf_mu, data = df_val_0)$predictions
  mu_A_hat <- ifelse(df_dr$A[val_idx] == 1, mu_1_hat, mu_0_hat)
  
  # Calculate DR Pseudo-Outcome (Equation 1 in Vansteelandt paper)
  pseudo_Y[val_idx] <- (df_dr$A[val_idx] - pi_hat) / (pi_hat * (1 - pi_hat)) * (df_dr$Y[val_idx] - mu_A_hat) + (mu_1_hat - mu_0_hat)
}

# Final DR-Learner Step: Regress pseudo-outcome on covariates to get CATE
df_pseudo <- data.frame(pseudo_Y = pseudo_Y, X_train_mat)
fit_dr <- ranger(pseudo_Y ~ ., data = df_pseudo)


# ==============================================================================
# 5. INFERENCE (CATEs and ATE)
# ==============================================================================

# Extract CATEs for Semi-Parametric
tau_draws_nbcf <- predict_linear_bcf_patched(fit_nbcf, X = X_train_mat, Z = as.matrix(Z_actg))$tau_hat
cate_hat_nbcf <- rowMeans(tau_draws_nbcf)
ate_draws_nbcf <- colMeans(tau_draws_nbcf)
ate_res_nbcf <- c(Mean = mean(ate_draws_nbcf), quantile(ate_draws_nbcf, probs = c(0.025, 0.975)))

# Extract CATEs & ATE for DR-Learner
cate_hat_dr <- fit_dr$predictions
ate_res_dr <- mean(pseudo_Y) # AIPW estimator of ATE

cat("\n--- IN-SAMPLE INFERENCE RESULTS ---\n")
cat(sprintf("Semi-Parametric ATE: %.2f [95%% CI: %.2f, %.2f]\n", ate_res_nbcf[1], ate_res_nbcf[2], ate_res_nbcf[3]))
cat(sprintf("Vansteelandt DR-Learner ATE (AIPW): %.2f\n", ate_res_dr))
cat(sprintf("CATE Correlation (r): %.3f\n", cor(cate_hat_nbcf, cate_hat_dr)))


# ==============================================================================
# 6. GGPLOT VISUALIZATIONS
# ==============================================================================
plot_data <- data.frame(
  CD4_Baseline = actg_sub$cd40,
  CATE_NBCF = cate_hat_nbcf,
  CATE_DR = cate_hat_dr
)

# A. CATE Correlation Plot (with Total Least Squares line)
pca_fit <- prcomp(~ CATE_NBCF + CATE_DR, data = plot_data)
tls_slope <- pca_fit$rotation["CATE_DR", "PC1"] / pca_fit$rotation["CATE_NBCF", "PC1"]
tls_intercept <- pca_fit$center["CATE_DR"] - tls_slope * pca_fit$center["CATE_NBCF"]

p_cor <- ggplot(plot_data, aes(x = CATE_NBCF, y = CATE_DR)) +
  geom_point(alpha = 0.5, color = "darkblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_abline(slope = tls_slope, intercept = tls_intercept, color = "orange", size = 1.2) +
  labs(title = paste0("CATE Correlation (r = ", round(cor(cate_hat_nbcf, cate_hat_dr), 3), ")"),
       subtitle = "Orange: Total Least Squares | Red Dashed: 1:1 Line",
       x = "Estimated CATE (Semi-Parametric BCF)", 
       y = "Estimated CATE (Vansteelandt DR-Learner)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# B. Density Distributions of Estimated CATEs
density_data <- data.frame(
  CATE = c(plot_data$CATE_NBCF, plot_data$CATE_DR),
  Model = rep(c("Semi-Parametric BCF", "Vansteelandt DR-Learner"), each = nrow(plot_data))
)

p_dist <- ggplot(density_data, aes(x = CATE, fill = Model, color = Model)) +
  geom_density(alpha = 0.4, size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  scale_fill_manual(values = c("#008080", "purple")) +
  scale_color_manual(values = c("#008080", "purple")) +
  labs(title = "Distribution of Estimated CATEs",
       x = "Treatment Effect (Difference in CD4 at 20 Weeks)",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

# C. Treatment Effect Heterogeneity by Baseline CD4
p_het <- ggplot(plot_data) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BCF"), method = "loess", se = FALSE, size = 1.2) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_DR, color = "Vansteelandt DR-Learner"), method = "loess", se = FALSE, size = 1.2) +
  geom_point(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BCF"), alpha = 0.2) +
  geom_point(aes(x = CD4_Baseline, y = CATE_DR, color = "Vansteelandt DR-Learner"), alpha = 0.2) +
  scale_color_manual(values = c("Semi-Parametric BCF" = "#008080", "Vansteelandt DR-Learner" = "purple")) +
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
# 7. SHAPLEY ANALYSIS vs. VANSTEELANDT TE-VIMs
# ==============================================================================
cat("\n--- STARTING IMPORTANCE ANALYSIS ---\n")

# A. Compute Shapley for Semi-Parametric BCF
X_actg <- results_actg$X_final
if (!is.null(fit_nbcf$interaction_pairs)) {
  ipairs <- fit_nbcf$interaction_pairs
} else {
  var_info <- results_actg$X_final_var_info
  is_candidate <- var_info$is_continuous | var_info$is_binary
  ipairs <- interaction_pairs(ncol(X_actg), is_candidate)
}

shapley_results <- compute_shapley_all(
  X = X_actg, beta_post = fit_nbcf$Beta, 
  beta_int_post = fit_nbcf$Beta_int, ipairs = ipairs
)
p_shap_importance <- plot_shapley_importance_breakdown(shapley_results) + 
  labs(title = "Semi-Parametric BCF: Shapley Importance")
print(p_shap_importance)

# B. Compute Vansteelandt Leave-One-Out (LOO) TE-VIMs for the DR-Learner
# Represents the increase in MSE when a variable is removed from the CATE conditioning set
mse_full <- mean((pseudo_Y - cate_hat_dr)^2)
te_vims <- numeric(ncol(X_train_mat))
names(te_vims) <- colnames(X_train_mat)

for(j in 1:ncol(X_train_mat)) {
  col_to_drop <- names(te_vims)[j]
  df_loo <- df_pseudo[, -which(colnames(df_pseudo) == col_to_drop)]
  rf_loo <- ranger(pseudo_Y ~ ., data = df_loo)
  mse_loo <- mean((pseudo_Y - rf_loo$predictions)^2)
  te_vims[j] <- mse_loo - mse_full
}

te_vim_df <- data.frame(Feature = names(te_vims), TE_VIM = te_vims) %>%
  arrange(desc(TE_VIM))

p_tevims <- ggplot(te_vim_df, aes(x = TE_VIM, y = reorder(Feature, TE_VIM))) +
  geom_col(fill = "purple", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Vansteelandt DR-Learner: LOO TE-VIMs", 
       subtitle = "Increase in Mean Squared Error (MSE) when removed",
       x = "TE-VIM (MSE Difference)", y = "Feature") +
  theme_minimal(base_size = 14)
print(p_tevims)


# ==============================================================================
# 8. PROFILING THE NEGATIVE CATE TAIL 
# ==============================================================================
cat("\n--- PROFILING NEGATIVE RESPONDERS ---\n")

patient_profiles <- actg_sub %>%
  dplyr::select(cd40, cd80, age, wtkg, karnof, gender, homo, race, symptom, drugs, hemo, z30) %>%
  mutate(
    CATE_NBCF = cate_hat_nbcf,
    CATE_DR = cate_hat_dr,
    Harm_NBCF = ifelse(CATE_NBCF < 0, 1, 0),
    Harm_DR = ifelse(CATE_DR < 0, 1, 0)
  )

cat("Cross-Tabulation of Predicted Harm (CATE < 0) Agreement:\n")
print(table("Semi-Parametric BCF Harm" = patient_profiles$Harm_NBCF, 
            "DR-Learner Harm" = patient_profiles$Harm_DR))


# ==============================================================================
# 9 & 10. QINI CURVE EVALUATION (RAW & SMOOTHED)
# ==============================================================================
cat("\n--- GENERATING QINI CURVES ---\n")

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

qini_eval_df <- data.frame(
  Y = Y_actg, Z = Z_actg, 
  CATE_NBCF = cate_hat_nbcf, CATE_DR = cate_hat_dr
)

qini_nbcf <- get_qini(qini_eval_df, "CATE_NBCF", "Semi-Parametric BCF")
qini_dr <- get_qini(qini_eval_df, "CATE_DR", "Vansteelandt DR-Learner")
qini_plot_data <- bind_rows(qini_nbcf, qini_dr)

# 10. Smoothed Qini Curve
final_uplift <- tail(qini_nbcf$uplift, 1)

p_qini_smooth <- ggplot(qini_plot_data, aes(x = frac, y = uplift, color = Model)) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = final_uplift, 
           color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Smoothed Qini Curve: Semi-Parametric vs Vansteelandt DR",
    subtitle = "Cumulative gain by prioritizing highest predicted CATEs",
    x = "Fraction of Population Treated",
    y = "Cumulative Gain (Uplift in CD4 Count)"
  ) +
  scale_color_manual(values = c("Semi-Parametric BCF" = "#008080", "Vansteelandt DR-Learner" = "purple")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p_qini_smooth)

# AUQC
get_auqc <- function(qini_df) {
  x <- qini_df$frac; y <- qini_df$uplift
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
cat(sprintf("AUQC - Semi-Parametric BCF: %.2f\n", get_auqc(qini_nbcf)))
cat(sprintf("AUQC - Vansteelandt DR-Learner: %.2f\n", get_auqc(qini_dr)))

# ==============================================================================
# 11. SAVE PLOTS TO LOCAL DIRECTORY
# ==============================================================================
cat("\n--- SAVING PLOTS TO DIRECTORY ---\n")

# Define the target directory path (using forward slashes for R compatibility)
plot_dir <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/plots/ACTG"

# Check if the directory exists; if not, create it recursively
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat(sprintf("Created new directory at: %s\n", plot_dir))
}

# Define standard dimensions for PhD/Publication quality
plot_width  <- 8
plot_height <- 6
plot_dpi    <- 300

# Save each plot 
ggsave(filename = file.path(plot_dir, "01_CATE_Correlation.png"), 
       plot = p_cor, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(plot_dir, "02_CATE_Density_Distribution.png"), 
       plot = p_dist, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(plot_dir, "03_Treatment_Heterogeneity.png"), 
       plot = p_het, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(plot_dir, "04_SemiParametric_Shapley_Importance.png"), 
       plot = p_shap_importance, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(plot_dir, "05_Vansteelandt_LOO_TE_VIMs.png"), 
       plot = p_tevims, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(plot_dir, "06_Qini_Curve_Comparison.png"), 
       plot = p_qini_smooth, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

cat(sprintf("All 6 plots have been successfully saved to: %s\n", plot_dir))
cat("\n================ ANALYSIS COMPLETE ================\n")