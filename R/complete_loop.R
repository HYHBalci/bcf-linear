# ==============================================================================
# 1. LOAD LIBRARIES & EXTERNAL SCRIPTS
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(stochtree)   
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(ranger)
library(coda)
library(rpart)
library(rpart.plot)
library(BayesLogit)
library(MASS)

# Source auxiliary scripts (Ensure these are in your working directory)
source('R/shapley_aux.R')
source('R/old_linear_linear.R') # Loads fit_grouped_horseshoes_R and its helpers

# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PREDICTION
# ==============================================================================
standardize_X_by_index_new <- function(X_initial, process_data = TRUE, interaction_rule = c("continuous", "continuous_or_binary", "all"), cat_coding_method = c("sum", "difference")) {
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
      list(processed_data = processed_col_data, is_continuous = is_continuous, is_binary = is_binary, is_categorical = is_categorical, original_idx = j_idx, col_name = original_colnames[j_idx], should_include = should_include)
    })
    
    valid_cols <- Filter(function(x) x$should_include, all_cols_info)
    X_intermediate_list <- lapply(valid_cols, `[[`, "processed_data")
    names(X_intermediate_list) <- sapply(valid_cols, `[[`, "col_name")
    X_intermediate_df <- as.data.frame(X_intermediate_list, check.names = FALSE)
    
    intermediate_map <- data.frame(col_name_intermediate = names(X_intermediate_list), is_continuous = sapply(valid_cols, `[[`, "is_continuous"), is_binary = sapply(valid_cols, `[[`, "is_binary"), stringsAsFactors = FALSE)
    
    cont_cols <- intermediate_map$col_name_intermediate[intermediate_map$is_continuous]
    X_intermediate_df[cont_cols] <- lapply(X_intermediate_df[cont_cols], function(x) as.vector(scale(x)))
    X_final_with_int <- model.matrix(~ ., data = X_intermediate_df)
    X_final <- X_final_with_int[, -1, drop = FALSE]
    attr_assign <- attr(X_final_with_int, "assign")[-1]
    matched_rows <- match(colnames(X_intermediate_df)[attr_assign], intermediate_map$col_name_intermediate)
    
    X_final_var_info <- data.frame(is_continuous = intermediate_map$is_continuous[matched_rows], is_binary = intermediate_map$is_binary[matched_rows], col_name_final = colnames(X_final), stringsAsFactors = FALSE)
    non_continous_idx_cpp <- which(!X_final_var_info$is_continuous) - 1
  } else {
    X_final <- as.matrix(X_df)
    X_final_var_info <- data.frame(is_continuous = TRUE, is_binary = FALSE, col_name_final = colnames(X_final))
    non_continous_idx_cpp <- numeric(0)
  }
  
  p_int_calculated <- 0; n_final_cols <- ncol(X_final)
  if (n_final_cols > 1) {
    is_candidate <- switch(interaction_rule, "all" = rep(TRUE, n_final_cols), "continuous_or_binary" = (X_final_var_info$is_continuous | X_final_var_info$is_binary), "continuous" = X_final_var_info$is_continuous)
    is_candidate[is.na(is_candidate)] <- FALSE
    for (i in 1:(n_final_cols - 1)) { for (j in (i + 1):n_final_cols) { if (is_candidate[i] || is_candidate[j]) p_int_calculated <- p_int_calculated + 1 } }
  }
  return(list(X_final = X_final, p_int = p_int_calculated, non_continous_idx_cpp = non_continous_idx_cpp, X_final_var_info = X_final_var_info))
}

predict_linear_bcf_patched <- function(object, X, Z, propensity = NULL) {
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
      } else stop("Propensity scores must be provided for this model.")
    }
  }
  
  X_forest <- X_linear
  if (object$model_params$propensity_covariate != "none") X_forest <- cbind(X_linear, propensity)
  forest_dataset_pred <- createForestDataset(X_forest, Z)  
  prop_sep <- object$model_params$propensity_seperate
  if (!is.null(prop_sep) && prop_sep == "tau") X_linear <- cbind(X_linear, propensity)
  
  num_chains <- dim(object$Beta)[1]; num_mcmc <- dim(object$Beta)[2]; total_samples <- num_chains * num_mcmc
  alpha_samples <- as.vector(t(object$alpha))
  beta_samples <- matrix(aperm(object$Beta, c(2, 1, 3)), nrow = total_samples)
  
  has_interactions <- !is.null(object$Beta_int) && (dim(object$Beta_int)[3] > 0)
  if (has_interactions) beta_int_samples <- matrix(aperm(object$Beta_int, c(2, 1, 3)), nrow = total_samples)
  
  linear_pred <- matrix(rep(alpha_samples, each=nrow(X)), nrow=nrow(X)) + as.matrix(X_linear) %*% t(beta_samples)
  
  if (has_interactions) {
    int_pairs <- object$interaction_pairs; num_interactions <- ncol(int_pairs)
    X_int <- matrix(0, nrow = nrow(X_linear), ncol = num_interactions)
    for (k in 1:num_interactions) { X_int[, k] <- X_linear[, int_pairs[1, k]] * X_linear[, int_pairs[2, k]] }
    linear_pred <- linear_pred + (X_int %*% t(beta_int_samples))
  }
  
  y_std <- object$model_params$outcome_scale; y_bar <- object$model_params$outcome_mean
  mu_hat <- object$forests_mu$predict(forest_dataset_pred) * y_std + y_bar
  tau_hat_total <- if (object$model_params$standardize) linear_pred * y_std else linear_pred
  
  return(list(mu_hat = mu_hat, tau_hat = tau_hat_total, y_hat = mu_hat + (tau_hat_total * as.vector(Z))))
}

# Functions for Tolerance Bands
postsum <- function(posters, conf = 0.05) { qs <- quantile(posters, p = c(conf/2, 1 - conf/2)); return(c(qs[1], mean(posters), qs[2])) }
intfrompost1 <- function(post, conf = 0.05){ t(apply(post, 1, postsum, conf = conf)) }
conftol <- function(margconf = 0.05, post, tol = 0.05){
  int <- intfrompost1(post, conf = margconf)
  difl <- ((post - int[,1]) >= 0) * ((post - int[,3]) <= 0)
  return(mean(colMeans(difl) >= 1 - tol))
} 
findconfglob <- function(post, tol = 0.05){
  grid <- seq(0.0005, 0.05, by = 0.0005)
  conftols <- sapply(grid, conftol, post = post, tol = tol)
  wm <- which(conftols - 0.95 < 0)[1] - 1
  if(is.na(wm) || wm < 1) return(0.0005) 
  return(grid[wm])
}

build_tolerance_plot <- function(tau_draws, alpha_draws, title_suffix) {
  het_draws <- t(t(tau_draws) - alpha_draws)
  confglob <- findconfglob(het_draws, tol = 0.02)
  cred_int <- intfrompost1(het_draws, conf = 0.05)
  tol_int  <- intfrompost1(het_draws, conf = confglob)
  
  het_df <- data.frame(Mean = cred_int[, 2], Cred_Lower = cred_int[, 1], Cred_Upper = cred_int[, 3], Tol_Lower = tol_int[, 1], Tol_Upper = tol_int[, 3]) %>%
    mutate(Significance = case_when(Tol_Lower > 0 ~ "Significant Positive", Tol_Upper < 0 ~ "Significant Negative", TRUE ~ "Not Significant")) %>%
    arrange(Mean) %>% mutate(Ordered_ID = row_number())
  
  pct_sig <- sum(het_df$Significance != "Not Significant") / nrow(het_df) * 100
  
  ggplot(het_df, aes(x = Ordered_ID)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_ribbon(aes(ymin = Tol_Lower, ymax = Tol_Upper, fill = "Tolerance Interval"), alpha = 0.6) +
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper, fill = "95% Credible Interval"), alpha = 0.85) +
    geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
    geom_point(data = filter(het_df, Significance != "Not Significant"), aes(y = Mean, color = Significance), size = 1.8, stroke = 0.5) +
    scale_fill_manual(values = c("Tolerance Interval" = "#D0D1E6", "95% Credible Interval" = "#74A9CF")) +
    scale_color_manual(values = c("Significant Positive" = "#005a5a", "Significant Negative" = "#cc0000")) +
    labs(title = paste("Individual Treatment Effect Heterogeneity", title_suffix),
         subtitle = paste0("Deep Blue: 95% CI | Light Blue: Tolerance Band (tol=0.02) | ", round(pct_sig, 1), "% significant"),
         x = "Patients", y = "\u0394 CATE") +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# ==============================================================================
# 3. MASTER ANALYSIS FUNCTION (RUNS A SINGLE DATASET COMPLETELY)
# ==============================================================================
analyze_competition_dataset <- function(file_path, base_out_dir) {
  
  # 1. Dataset ID and Directory Setup
  dataset_id <- trimws(tools::file_path_sans_ext(basename(file_path)))
  cat(sprintf("\n\n=======================================================\n"))
  cat(sprintf("▶ PROCESSING: %s\n", dataset_id))
  cat(sprintf("=======================================================\n"))
  
  plot_dir <- file.path(base_out_dir, dataset_id, "plots")
  model_dir <- file.path(base_out_dir, dataset_id, "models")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
  
  # 2. Load Data and Detect Outcome Type
  comp_data <- read.csv(file_path)
  Y_comp <- comp_data$Y
  Z_comp <- as.numeric(comp_data$W)
  
  # DYNAMIC COVARIATE SELECTION (Finds all columns like x1, x2 ... x50)
  x_cols <- grep("^x[0-9]+$", names(comp_data), value = TRUE)
  X_comp_raw <- comp_data[, x_cols, drop = FALSE]
  
  cat(sprintf("Detected %d covariates.\n", length(x_cols)))
  
  # DYNAMIC OUTCOME DETECTION
  is_binary <- length(unique(na.omit(Y_comp))) == 2
  cat(sprintf("Detected Outcome Type: %s\n", ifelse(is_binary, "BINARY", "CONTINUOUS")))
  
  results_comp <- standardize_X_by_index_new(X_initial = X_comp_raw, process_data = TRUE, interaction_rule = "continuous_or_binary", cat_coding_method = "sum")
  X_train_mat <- results_comp$X_final
  pi_hat <- mean(Z_comp) 
  
  # 3. Fit Models
  general_params_default <- list(
    cutpoint_grid_size = 100, standardize = !is_binary, sample_sigma2_global = !is_binary, sigma2_global_init = 1, 
    sigma2_global_shape = 1, sigma2_global_scale = 0.001, variable_weights = NULL, propensity_covariate = "mu", 
    adaptive_coding = FALSE, control_coding_init = -0.5, treated_coding_init = 0.5, rfx_prior_var = NULL, 
    random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE, keep_every = 1, num_chains = 1, verbose = FALSE, 
    sample_global_prior = "half-cauchy", unlink = TRUE, propensity_seperate = "none", gibbs = TRUE, step_out = 0.5, 
    max_steps = 150, save_output = FALSE, interaction_rule = "continuous_or_binary", 
    standardize_cov = FALSE, simple_prior = FALSE, save_partial_residual = FALSE, regularize_ATE = FALSE,
    sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1,
    probit_outcome_model = is_binary # <-- Set Dynamically
  )
  
  # A. BCF (HORSESHOE)
  cat("  -> Fitting BCF (Horseshoe)...\n")
  fit_nbcf <- bcf_linear_probit(X_train = X_train_mat, y_train = Y_comp, Z_train = Z_comp, num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = general_params_default)
  saveRDS(fit_nbcf, file.path(model_dir, "fit_nbcf.rds"))
  
  # B. BCF (NO SHRINKAGE)
  cat("  -> Fitting BCF (No Shrinkage)...\n")
  general_params_no_shrinkage <- general_params_default
  general_params_no_shrinkage$sample_global_prior <- "none" 
  fit_nbcf_no_shrink <- bcf_linear_probit(X_train = X_train_mat, y_train = Y_comp, Z_train = Z_comp, num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = general_params_no_shrinkage)
  saveRDS(fit_nbcf_no_shrink, file.path(model_dir, "fit_nbcf_no_shrink.rds"))
  
  # C. BCF (HALF-NORMAL)
  cat("  -> Fitting BCF (Half-Normal)...\n")
  general_params_hn <- general_params_default
  general_params_hn$sample_global_prior <- "half-normal"
  general_params_hn$hn_scale <- 0.5 
  fit_nbcf_hn <- bcf_linear_probit(X_train = X_train_mat, y_train = Y_comp, Z_train = Z_comp, num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = general_params_hn)
  saveRDS(fit_nbcf_hn, file.path(model_dir, "fit_nbcf_hn.rds"))
  
  # D. DR-Learner (SuperLearner)
  cat("  -> Fitting DR-Learner (SuperLearner)...\n")
  sl_library <- c("SL.glm", "SL.glmnet", "SL.gam", "SL.xgboost", "SL.ranger")
  sl_family <- if (is_binary) binomial() else gaussian() # <-- Set Dynamically
  
  df_dr <- data.frame(Y = Y_comp, A = Z_comp, X_train_mat)
  set.seed(123)
  K_folds <- 10
  fold_ids <- sample(rep(1:K_folds, length.out = nrow(df_dr)))
  pseudo_Y <- numeric(nrow(df_dr))
  
  for (k in 1:K_folds) {
    train_idx <- which(fold_ids != k); val_idx <- which(fold_ids == k)
    sl_mu <- SuperLearner(Y = df_dr$Y[train_idx], X = df_dr[train_idx, setdiff(names(df_dr), "Y")], SL.library = sl_library, family = sl_family, cvControl = list(V = 10))
    X_val_1 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_1$A <- 1
    X_val_0 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_0$A <- 0
    mu_1_hat <- predict(sl_mu, newdata = X_val_1)$pred
    mu_0_hat <- predict(sl_mu, newdata = X_val_0)$pred
    pseudo_Y[val_idx] <- (df_dr$A[val_idx] - pi_hat) / (pi_hat * (1 - pi_hat)) * (df_dr$Y[val_idx] - ifelse(df_dr$A[val_idx] == 1, mu_1_hat, mu_0_hat)) + (mu_1_hat - mu_0_hat)
  }
  fit_dr <- SuperLearner(Y = pseudo_Y, X = as.data.frame(X_train_mat), SL.library = sl_library, family = gaussian(), cvControl = list(V = 10))
  saveRDS(fit_dr, file.path(model_dir, "fit_dr.rds"))
  cate_hat_dr <- as.vector(fit_dr$SL.predict)
  
  # --- NEW: Extract and Save SuperLearner Weights ---
  cat("  -> Saving SuperLearner Weights...\n")
  sl_weights <- data.frame(
    Learner = names(fit_dr$coef),
    Weight = as.numeric(fit_dr$coef)
  )
  write.csv(sl_weights, file.path(model_dir, "SuperLearner_Weights.csv"), row.names = FALSE)
  
  # E. Grouped Horseshoe
  cat("  -> Fitting Custom Grouped Horseshoe...\n")
  ghs_family <- if (is_binary) "binomial" else "gaussian" # <-- Set Dynamically
  fit_ghs <- fit_grouped_horseshoes_R(
    y_vec = Y_comp, X_mat = X_train_mat, Z_vec = Z_comp, 
    family = ghs_family, n_iter = 3000, burn_in = 1000, num_chains = 1,
    propensity_train = rep(pi_hat, length(Z_comp)), propensity_as_covariate = T
  )
  saveRDS(fit_ghs, file.path(model_dir, "fit_ghs.rds"))
  
  # 4. Extract CATEs
  cate_hat_nbcf <- rowMeans(predict_linear_bcf_patched(fit_nbcf, X = X_train_mat, Z = as.matrix(Z_comp))$tau_hat)
  cate_hat_nbcf_no_shrink <- rowMeans(predict_linear_bcf_patched(fit_nbcf_no_shrink, X = X_train_mat, Z = as.matrix(Z_comp))$tau_hat)
  cate_hat_nbcf_hn <- rowMeans(predict_linear_bcf_patched(fit_nbcf_hn, X = X_train_mat, Z = as.matrix(Z_comp))$tau_hat)
  
  # Grouped Horseshoe CATE Extraction
  boolean_continuous <- as.vector(results_comp$X_final_var_info$is_continuous | results_comp$X_final_var_info$is_binary)
  main_interaction_indices <- create_interaction_pairs_R(ncol(X_train_mat), boolean_continuous)
  X_int_mat <- NULL
  if (length(main_interaction_indices) > 0) {
    X_int_mat <- matrix(0, nrow = nrow(X_train_mat), ncol = length(main_interaction_indices))
    for (k in seq_along(main_interaction_indices)) { X_int_mat[, k] <- X_train_mat[, main_interaction_indices[[k]][1]] * X_train_mat[, main_interaction_indices[[k]][2]] }
  }
  cate_hat_ghs <- mean(fit_ghs$aleph) + as.vector(X_train_mat %*% colMeans(fit_ghs$gamma))
  if (!is.null(X_int_mat) && ncol(fit_ghs$gamma_int) > 0) cate_hat_ghs <- cate_hat_ghs + as.vector(X_int_mat %*% colMeans(fit_ghs$gamma_int))
  
  # ---------------------------------------------------------
  # 5. Determine Most Important Feature
  # ---------------------------------------------------------
  cat("\n  -> Computing Shapley & Identifying Top Feature for Plotting...\n")
  var_info <- results_comp$X_final_var_info
  ipairs <- interaction_pairs(ncol(X_train_mat), var_info$is_continuous | var_info$is_binary)
  shap_bcf <- compute_shapley_all(X_train_mat, fit_nbcf$Beta, fit_nbcf$Beta_int, ipairs)
  
  # Extract numeric columns only to calculate row sums
  shap_numeric <- shap_bcf[, sapply(shap_bcf, is.numeric), drop = FALSE]
  top_row_idx <- which.max(rowSums(abs(shap_numeric)))
  top_feature_name <- as.character(shap_bcf$feature[top_row_idx])
  
  # Safely match to the design matrix
  if (top_feature_name %in% colnames(X_train_mat)) {
    top_feature_vector <- X_train_mat[, top_feature_name]
  } else {
    top_feature_name <- colnames(X_train_mat)[1] # Fallback
    top_feature_vector <- X_train_mat[, 1]
  }
  
  p_shap_horseshoe <- plot_shapley_importance_breakdown(shap_bcf) + labs(title = "BCF (Horseshoe): Shapley Importance")
  p_shap_hn <- plot_shapley_importance_breakdown(compute_shapley_all(X_train_mat, fit_nbcf_hn$Beta, fit_nbcf_hn$Beta_int, ipairs)) + labs(title = "BCF (Half-Normal): Shapley Importance")
  p_shap_no_shrink <- plot_shapley_importance_breakdown(compute_shapley_all(X_train_mat, fit_nbcf_no_shrink$Beta, fit_nbcf_no_shrink$Beta_int, ipairs)) + labs(title = "BCF (No Shrinkage): Shapley Importance")
  
  # 6. Generate Plots
  cat("  -> Generating Plots...\n")
  color_palette <- c("BCF (Horseshoe)" = "#008080", "BCF (Half-Normal)" = "dodgerblue", "BCF (No Shrinkage)" = "darkorange", "Vansteelandt DR" = "purple", "Grouped Horseshoe Linear" = "#E6194B")
  
  density_data <- data.frame(
    CATE = c(cate_hat_nbcf, cate_hat_nbcf_hn, cate_hat_nbcf_no_shrink, cate_hat_dr, cate_hat_ghs),
    Model = rep(names(color_palette), each = nrow(X_train_mat))
  )
  
  p_dist <- ggplot(density_data, aes(x = CATE, fill = Model, color = Model)) +
    geom_density(alpha = 0.25, linewidth = 1) + geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
    scale_fill_manual(values = color_palette) + scale_color_manual(values = color_palette) +
    labs(title = paste("Distribution of CATEs -", dataset_id), x = "Treatment Effect", y = "Density") + theme_minimal(base_size = 14)
  
  plot_data <- data.frame(Top_Feature = top_feature_vector, CATE_NBCF = cate_hat_nbcf, CATE_NBCF_NoShrink = cate_hat_nbcf_no_shrink, CATE_NBCF_HN = cate_hat_nbcf_hn, CATE_DR = cate_hat_dr, CATE_GHS = cate_hat_ghs)
  
  p_het <- ggplot(plot_data) +
    geom_smooth(aes(x = Top_Feature, y = CATE_NBCF, color = "BCF (Horseshoe)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = Top_Feature, y = CATE_NBCF_HN, color = "BCF (Half-Normal)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = Top_Feature, y = CATE_NBCF_NoShrink, color = "BCF (No Shrinkage)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = Top_Feature, y = CATE_DR, color = "Vansteelandt DR"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = Top_Feature, y = CATE_GHS, color = "Grouped Horseshoe Linear"), method = "loess", se = FALSE, linewidth = 1.2) +
    scale_color_manual(values = color_palette) +
    labs(title = paste("Treatment Heterogeneity vs", top_feature_name), x = top_feature_name, y = "Estimated CATE") + theme_minimal(base_size = 14)
  
  # Qini and Uplift Curves
  get_eval_curves <- function(df, model_col, model_name) {
    df %>% arrange(desc(.data[[model_col]])) %>% mutate(
      cum_Z = cumsum(Z), cum_C = cumsum(1 - Z), cum_Y_Z = cumsum(Y * Z), cum_Y_C = cumsum(Y * (1 - Z)),
      qini = cum_Y_Z - cum_Y_C * (sum(Z) / sum(1 - Z)),
      cum_ATE = ifelse(cum_Z > 0 & cum_C > 0, (cum_Y_Z / cum_Z) - (cum_Y_C / cum_C), 0),
      uplift = cum_ATE * (cum_Z + cum_C), frac = row_number() / n(), Model = model_name
    ) %>% dplyr::select(frac, qini, uplift, Model)
  }
  
  qini_eval_df <- data.frame(Y = Y_comp, Z = Z_comp, CATE_NBCF = cate_hat_nbcf, CATE_NBCF_HN = cate_hat_nbcf_hn, CATE_NBCF_NoShrink = cate_hat_nbcf_no_shrink, CATE_DR = cate_hat_dr, CATE_GHS = cate_hat_ghs)
  
  eval_data <- bind_rows(
    get_eval_curves(qini_eval_df, "CATE_NBCF", "BCF (Horseshoe)"), get_eval_curves(qini_eval_df, "CATE_NBCF_HN", "BCF (Half-Normal)"),
    get_eval_curves(qini_eval_df, "CATE_NBCF_NoShrink", "BCF (No Shrinkage)"), get_eval_curves(qini_eval_df, "CATE_DR", "Vansteelandt DR"), get_eval_curves(qini_eval_df, "CATE_GHS", "Grouped Horseshoe Linear")
  )
  
  p_qini <- ggplot(eval_data, aes(x = frac, y = qini, color = Model)) + geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) + scale_color_manual(values = color_palette) + labs(title = paste("Qini Curve -", dataset_id), x = "Fraction Treated", y = "Cumulative Qini") + theme_minimal()
  p_uplift <- ggplot(eval_data, aes(x = frac, y = uplift, color = Model)) + geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) + scale_color_manual(values = color_palette) + labs(title = paste("Robust Uplift -", dataset_id), x = "Fraction Treated", y = "True Cumulative Uplift") + theme_minimal()
  
  # Tolerance Bands
  p_het_tol_horseshoe <- build_tolerance_plot(predict_linear_bcf_patched(fit_nbcf, X_train_mat, Z_comp)$tau_hat, as.vector(t(fit_nbcf$alpha)), "(Horseshoe)")
  p_het_tol_hn <- build_tolerance_plot(predict_linear_bcf_patched(fit_nbcf_hn, X_train_mat, Z_comp)$tau_hat, as.vector(t(fit_nbcf_hn$alpha)), "(Half-Normal)")
  p_het_tol_noshrink <- build_tolerance_plot(predict_linear_bcf_patched(fit_nbcf_no_shrink, X_train_mat, Z_comp)$tau_hat, as.vector(t(fit_nbcf_no_shrink$alpha)), "(No Shrinkage)")
  
  # Subgroup Extraction (Tree) - DYNAMIC COVARIATES FIX
  extract_tree <- function(cate_preds, m_name, f_name) {
    tr_data <- data.frame(CATE = cate_preds)
    tr_data <- cbind(tr_data, comp_data[, x_cols, drop = FALSE]) # <-- Uses dynamic x_cols!
    
    tr <- rpart(CATE ~ ., data = tr_data, method = "anova", control = rpart.control(cp = 0.01, maxdepth = 3))
    png(filename = file.path(plot_dir, paste0("08_Subgroup_Tree_", f_name, ".png")), width = 2400, height = 1800, res = 300)
    rpart.plot(tr, type = 4, extra = 1, roundint = FALSE, main = paste("Subgroup Rule:", m_name), box.palette = c("tomato", "white", "#008080"), shadow.col = "gray", nn = TRUE)
    dev.off()
  }
  
  extract_tree(cate_hat_nbcf, "BCF (Horseshoe)", "Horseshoe")
  extract_tree(cate_hat_dr, "Vansteelandt DR", "DR_Learner")
  
  # 7. Save Everything
  cat("  -> Saving Plots to Folder...\n")
  ggsave(file.path(plot_dir, "01_Density_All.png"), plot = p_dist, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "02_Heterogeneity.png"), plot = p_het, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "03a_Shapley_Horseshoe.png"), plot = p_shap_horseshoe, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "03b_Shapley_HalfNormal.png"), plot = p_shap_hn, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "03c_Shapley_NoShrink.png"), plot = p_shap_no_shrink, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "04_Qini.png"), plot = p_qini, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "05_Uplift.png"), plot = p_uplift, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "06a_Tol_Horseshoe.png"), plot = p_het_tol_horseshoe, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "06b_Tol_HalfNormal.png"), plot = p_het_tol_hn, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, "06c_Tol_NoShrink.png"), plot = p_het_tol_noshrink, width = 10, height = 7, dpi = 300, bg = "white")
  
  cat(sprintf("✔ SUCCESS: All outputs saved for %s\n", dataset_id))
}

# ==============================================================================
# 4. EXECUTION LOOP: PROCESS ALL DATASETS IN THE FOLDER
# ==============================================================================

# Define Paths
data_directory <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/Data/Datasets/Blinded_data_sets/"
base_output_directory <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/Results/"

# Find all simulated CSV datasets (modify pattern if your files are named differently)
data_files <- list.files(path = data_directory, pattern = "dataset30_.*\\.csv", full.names = TRUE)

cat(sprintf("\nFOUND %d DATASETS TO PROCESS.\n", length(data_files)))

# Loop through each file and run the analysis
for (file_path in data_files) {
  tryCatch({
    analyze_competition_dataset(file_path, base_output_directory)
  }, error = function(e) {
    cat(sprintf("\n[ERROR] Failed processing %s:\n%s\n", basename(file_path), e$message))
  })
}

cat("\n================ ALL DATASETS COMPLETED ================\n")