# ==============================================================================
# AUTO-INJECTED ETA TIMER (adapted for CV)
# ==============================================================================
cv_seeds <- c(123, 456, 789)
K_folds <- 10
total_models_to_run <- length(cv_seeds) * K_folds * 4
models_run <- 0

bcf_linear_probit_eta <- function(...) {
  t0 <- Sys.time()
  res <- bcf_linear_probit(...)
  t1 <- Sys.time()
  dur <- as.numeric(difftime(t1, t0, units="mins"))
  cat(sprintf("\n[TIMER] Model finished in %.2f mins.\n\n", dur))
  return(res)
}

# ==============================================================================
# 1. LOAD LIBRARIES & PARALLEL SETUP
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(speff2trial) # Contains the ACTG175 dataset
library(stochtree)   # Semi-Parametric BCF and Standard BCF method
# DR-Learner Ensemble Libraries:
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(ranger)
library(coda)
library(foreach)
library(doParallel)
source('R/shapley_aux.R')

# Setup parallel cluster
num_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(num_cores, outfile = "")
registerDoParallel(cl)
cat(sprintf("\n[INFO] Registered Parallel Cluster with %d cores\n", num_cores))

# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PATCHED PREDICTIONS
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

# --- 1. PREDICT FUNCTION FOR SEMI-PARAMETRIC BART ---
predict_linear_bcf_patched <- function(object, X, Z, propensity = NULL, rfx_group_ids = NULL, rfx_basis = NULL, ...) {
  if ((!is.data.frame(X)) && (!is.matrix(X))) stop("X must be a matrix or dataframe")
  
  train_set_metadata <- object$train_set_metadata
  X_processed <- preprocessPredictionData(X, train_set_metadata)
  X_linear <- X_processed 
  
  if (!is.null(Z) && is.null(dim(Z))) Z <- as.matrix(as.numeric(Z))
  if (!is.null(propensity) && is.null(dim(propensity))) propensity <- as.matrix(propensity)
  
  prop_cov <- object$model_params$propensity_covariate
  prop_sep <- object$model_params$propensity_seperate
  has_prop_cov <- !is.null(prop_cov) && !is.na(prop_cov) && prop_cov != "none"
  has_prop_sep <- !is.null(prop_sep) && !is.na(prop_sep) && prop_sep != "none"
  
  if (has_prop_cov || has_prop_sep) {
    if (is.null(propensity)) {
      if (isTRUE(object$model_params$internal_propensity_model)) {
        propensity <- rowMeans(predict(object$bart_propensity_model, X_processed)$y_hat)
        if (is.null(dim(propensity))) propensity <- as.matrix(propensity)
      } else { stop("Propensity scores must be provided for this model.") }
    }
  }
  
  X_forest <- X_linear
  if (has_prop_cov) X_forest <- cbind(X_linear, propensity)
  forest_dataset_pred <- createForestDataset(X_forest, Z)  
  
  if (has_prop_sep && prop_sep == "tau") X_linear <- cbind(X_linear, propensity)
  
  num_chains <- dim(object$Beta)[1]
  num_mcmc <- dim(object$Beta)[2]
  total_samples <- num_chains * num_mcmc
  
  alpha_samples <- as.vector(t(object$alpha))
  beta_samples <- matrix(aperm(object$Beta, c(2, 1, 3)), nrow = total_samples)
  
  has_interactions <- !is.null(object$Beta_int) && (dim(object$Beta_int)[3] > 0)
  if (has_interactions) beta_int_samples <- matrix(aperm(object$Beta_int, c(2, 1, 3)), nrow = total_samples)
  
  p_beta <- ncol(beta_samples)
  p_linear <- ncol(X_linear)
  X_target <- X_linear
  
  if (p_beta == p_linear + 1) {
    X_target <- cbind(1, X_linear)
  } else if (p_beta != p_linear) {
    stop(paste("Cannot align X_linear (ncol =", p_linear, ") with beta_samples (ncol =", p_beta, ")"))
  }
  
  linear_pred <- matrix(rep(alpha_samples, each=nrow(X)), nrow=nrow(X)) + as.matrix(X_target) %*% t(beta_samples)
  
  if (has_interactions) {
    int_pairs <- object$interaction_pairs
    num_interactions <- ncol(int_pairs)
    X_int <- matrix(0, nrow = nrow(X_linear), ncol = num_interactions)
    for (k in 1:num_interactions) { X_int[, k] <- X_linear[, int_pairs[1, k]] * X_linear[, int_pairs[2, k]] }
    linear_pred <- linear_pred + (X_int %*% t(beta_int_samples))
  }
  
  y_std <- object$model_params$outcome_scale
  y_bar <- object$model_params$outcome_mean
  mu_hat <- object$forests_mu$predict(forest_dataset_pred) * y_std + y_bar
  is_std <- isTRUE(object$model_params$standardize)
  tau_hat_total <- if (is_std) linear_pred * y_std else linear_pred
  
  result <- list(mu_hat = mu_hat, tau_hat = tau_hat_total, y_hat = mu_hat + (tau_hat_total * as.vector(Z)))
  return(result)
}

# --- 2. PREDICT FUNCTION FOR STANDARD BCF (From Vignette) ---
predict.bcfmodel <- function(object, X, Z, propensity = NULL, rfx_group_ids = NULL, rfx_basis = NULL, ...){
  # Preprocess covariates
  if ((!is.data.frame(X)) && (!is.matrix(X))) { stop("X must be a matrix or dataframe") }
  train_set_metadata <- object$train_set_metadata
  X <- preprocessPredictionData(X, train_set_metadata)
  
  # Convert all input data to matrices if not already converted
  if ((is.null(dim(Z))) && (!is.null(Z))) { Z <- as.matrix(as.numeric(Z)) }
  if ((is.null(dim(propensity))) && (!is.null(propensity))) { propensity <- as.matrix(propensity) }
  if ((is.null(dim(rfx_basis))) && (!is.null(rfx_basis))) { rfx_basis <- as.matrix(rfx_basis) }
  
  # Safe checks for missing propensities
  prop_cov <- object$model_params$propensity_covariate
  has_prop_cov <- !is.null(prop_cov) && !is.na(prop_cov) && prop_cov != "none"
  
  if (has_prop_cov && is.null(propensity)) {
    if (!isTRUE(object$model_params$internal_propensity_model)) {
      stop("propensity must be provided for this model")
    }
    propensity <- rowMeans(predict(object$bart_propensity_model, X)$y_hat)
  }
  if (nrow(X) != nrow(Z)) { stop("X and Z must have the same number of rows") }
  if (object$model_params$num_covariates != ncol(X)) { stop("X and must have the same number of columns as the covariates used to train the model") }
  if (isTRUE(object$model_params$has_rfx) && (is.null(rfx_group_ids))) { stop("Random effect group labels (rfx_group_ids) must be provided for this model") }
  if (isTRUE(object$model_params$has_rfx_basis) && (is.null(rfx_basis))) { stop("Random effects basis (rfx_basis) must be provided for this model") }
  if (isTRUE(object$model_params$num_rfx_basis > 0) && (ncol(rfx_basis) != object$model_params$num_rfx_basis)) { stop("Random effects basis has a different dimension than the basis used to train this model") }
  
  has_rfx <- FALSE
  if (!is.null(rfx_group_ids)) {
    rfx_unique_group_ids <- object$rfx_unique_group_ids
    group_ids_factor <- factor(rfx_group_ids, levels = rfx_unique_group_ids)
    if (sum(is.na(group_ids_factor)) > 0) { stop("All random effect group labels provided in rfx_group_ids must be present in rfx_group_ids_train") }
    rfx_group_ids <- as.integer(group_ids_factor)
    has_rfx <- TRUE
  }
  
  if (isTRUE(object$model_params$has_rfx) && (is.null(rfx_basis))) {
    rfx_basis <- matrix(rep(1, nrow(X)), ncol = 1)
  }
  
  # Add propensities to covariate set if necessary (safely defining X_combined for both cases)
  if (has_prop_cov) {
    X_combined <- cbind(X, propensity)
  } else {
    X_combined <- X
  }
  
  # Create prediction datasets
  forest_dataset_pred <- createForestDataset(X_combined, Z)
  
  # Compute forest predictions
  num_samples <- object$model_params$num_samples
  y_std <- object$model_params$outcome_scale
  y_bar <- object$model_params$outcome_mean
  initial_sigma2 <- object$model_params$initial_sigma2
  mu_hat <- object$forests_mu$predict(forest_dataset_pred)*y_std + y_bar
  if (isTRUE(object$model_params$adaptive_coding)) {
    tau_hat_raw <- object$forests_tau$predict_raw(forest_dataset_pred)
    tau_hat <- t(t(tau_hat_raw) * (object$b_1_samples - object$b_0_samples))*y_std
  } else {
    tau_hat <- object$forests_tau$predict_raw(forest_dataset_pred)*y_std
  }
  if (isTRUE(object$model_params$include_variance_forest)) {
    s_x_raw <- object$forests_variance$predict(forest_dataset_pred)
  }
  
  if (isTRUE(object$model_params$has_rfx)) {
    rfx_predictions <- object$rfx_samples$predict(rfx_group_ids, rfx_basis)*y_std
  }
  
  y_hat <- mu_hat + tau_hat * as.numeric(Z)
  if (isTRUE(object$model_params$has_rfx)) y_hat <- y_hat + rfx_predictions
  
  if (isTRUE(object$model_params$include_variance_forest)) {
    if (isTRUE(object$model_params$sample_sigma2_global)) {
      sigma2_global_samples <- object$sigma2_global_samples
      variance_forest_predictions <- sapply(1:num_samples, function(i) s_x_raw[,i]*sigma2_global_samples[i])
    } else {
      variance_forest_predictions <- s_x_raw*initial_sigma2*y_std*y_std
    }
  }
  
  result <- list("mu_hat" = mu_hat, "tau_hat" = tau_hat, "y_hat" = y_hat)
  if (isTRUE(object$model_params$has_rfx)) { result[["rfx_predictions"]] = rfx_predictions }
  if (isTRUE(object$model_params$include_variance_forest)) { result[["variance_forest_predictions"]] = variance_forest_predictions }
  return(result)
}

# ==============================================================================
# 3. PREPARE ACTG175 DATA
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

# Process all data together to ensure factors and scaling remain consistent 
results_actg <- standardize_X_by_index_new(
  X_initial = X_actg_raw, 
  process_data = TRUE, 
  interaction_rule = "continuous_or_binary", 
  cat_coding_method = "sum"
)
X_all_mat <- results_actg$X_final

# ==============================================================================
# 3b. 10-FOLD CV SETUP & STORAGE (Multiple Seeds)
# ==============================================================================
for (current_seed in cv_seeds) {
  cat(sprintf("\n====================================================================\n"))
  cat(sprintf("   STARTING NEW 10-FOLD CV ITERATION (SEED: %d)\n", current_seed))
  cat(sprintf("====================================================================\n"))
  
  set.seed(current_seed)
  n_total <- nrow(X_all_mat)
  fold_ids_outer <- sample(rep(1:K_folds, length.out = n_total))

  oos_cate_nbcf <- numeric(n_total)
  oos_cate_nbcf_noshrink <- numeric(n_total)
  oos_cate_nbcf_ols <- numeric(n_total)
  oos_cate_nbcf_robust <- numeric(n_total)
  oos_cate_bcf <- numeric(n_total)
  oos_cate_dr <- numeric(n_total)
  
  oos_tau_draws_nbcf <- list()
  oos_tau_draws_nbcf_noshrink <- list()
  oos_tau_draws_nbcf_ols <- list()
  oos_tau_draws_nbcf_robust <- list()
  oos_tau_draws_bcf <- list()
  
  fit_nbcf_fold1 <- NULL
  fit_dr_fold1 <- NULL
  pseudo_Y_fold1 <- NULL
  X_train_mat_fold1 <- NULL
  X_train_df_fold1 <- NULL

  general_params_nbcf <- list(
    cutpoint_grid_size = 100, standardize = TRUE, 
    sample_sigma2_global = TRUE, sigma2_global_init = 1, 
    sigma2_global_shape = 1, sigma2_global_scale = 0.001,
    variable_weights = NULL, propensity_covariate = "none", 
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

  general_params_nbcf_noshrink <- general_params_nbcf
  general_params_nbcf_noshrink$sample_global_prior <- "none"
  
  general_params_nbcf_ols <- general_params_nbcf
  general_params_nbcf_ols$sample_global_prior <- "OLS"

  general_params_nbcf_robust <- general_params_nbcf
  general_params_nbcf_robust$sample_global_prior <- "half-cauchy"
  general_params_nbcf_robust$robust <- TRUE
  general_params_nbcf_robust$robust_nu <- 3

  general_params_bcf <- list(
    cutpoint_grid_size = 100, standardize = TRUE, 
    sample_sigma2_global = TRUE, sigma2_global_init = 1, 
    sigma2_global_shape = 1, sigma2_global_scale = 0.001, 
    variable_weights = NULL, propensity_covariate = "none", 
    adaptive_coding = TRUE, control_coding_init = 0, 
    treated_coding_init = 1, rfx_prior_var = NULL, 
    random_seed = 1848, keep_burnin = FALSE, keep_gfr = FALSE, 
    keep_every = 1, num_chains = 2, verbose = FALSE, 
    probit_outcome_model = FALSE
  )

  sl_library <- c("SL.glm", "SL.glmnet", "SL.gam", "SL.xgboost", "SL.ranger")

  # ==============================================================================
  # 4 & 5. FIT MODELS & INFERENCE IN 10-FOLD CV
  # ==============================================================================
  library(doParallel)
  registerDoParallel(cores = parallel::detectCores() - 1)
  
  cat("\nStarting 10-Fold CV (90/10 Splits)... (Grab a coffee, this takes time)\n")
  cv_start_time <- Sys.time()
  
  # Parallel execution of K_folds
  fold_results <- foreach(k = 1:K_folds, 
                          .packages = c("dplyr", "tidyr", "SuperLearner", "glmnet", "gam", "xgboost", "ranger"),
                          .export = c("bcf_linear_probit_eta", "predict_linear_bcf_patched", "predict.bcfmodel", 
                                      "general_params_nbcf", "general_params_nbcf_noshrink", "general_params_nbcf_ols", 
                                      "general_params_nbcf_robust", "general_params_bcf", "sl_library", 
                                      "X_all_mat", "Y_actg", "Z_actg", "fold_ids_outer")) %dopar% {
    
    cat(sprintf("\n[Worker] Starting Fold %d\n", k))
    test_idx <- which(fold_ids_outer == k)
    train_idx <- which(fold_ids_outer != k)
    
    X_train_mat <- X_all_mat[train_idx, , drop = FALSE]
    Y_train <- Y_actg[train_idx]
    Z_train <- Z_actg[train_idx]
    
    X_test_mat <- X_all_mat[test_idx, , drop = FALSE]
    Z_test <- Z_actg[test_idx]
    
    # --- Semi-Parametric BART ---
    cat(sprintf("  Fold %d: Fitting Semi-Parametric BART...\n", k))
    fit_nbcf <- bcf_linear_probit_eta(
      X_train = X_train_mat, y_train = Y_train, Z_train = Z_train - 0.5,
      num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
      general_params = general_params_nbcf
    )
    tau_draws_nbcf_test <- predict_linear_bcf_patched(fit_nbcf, X = X_test_mat, Z = Z_test)$tau_hat
    y_std_nbcf <- fit_nbcf$model_params$outcome_scale
    alpha_draws_nbcf <- as.vector(t(fit_nbcf$alpha)) * y_std_nbcf
    het_draws_nbcf_fold <- t(t(tau_draws_nbcf_test) - alpha_draws_nbcf)
    
    # --- Semi-Parametric BART (No Shrinkage) ---
    cat(sprintf("  Fold %d: Fitting Semi-Parametric BART (No Shrinkage)...\n", k))
    fit_nbcf_noshrink <- bcf_linear_probit_eta(
      X_train = X_train_mat, y_train = Y_train, Z_train = Z_train - 0.5,
      num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
      general_params = general_params_nbcf_noshrink
    )
    tau_draws_nbcf_noshrink_test <- predict_linear_bcf_patched(fit_nbcf_noshrink, X = X_test_mat, Z = Z_test)$tau_hat
    y_std_nbcf_noshrink <- fit_nbcf_noshrink$model_params$outcome_scale
    alpha_draws_nbcf_noshrink <- as.vector(t(fit_nbcf_noshrink$alpha)) * y_std_nbcf_noshrink
    het_draws_nbcf_noshrink_fold <- t(t(tau_draws_nbcf_noshrink_test) - alpha_draws_nbcf_noshrink)
  
    # --- Semi-Parametric BART (OLS) ---
    cat(sprintf("  Fold %d: Fitting Semi-Parametric BART (OLS)...\n", k))
    fit_nbcf_ols <- bcf_linear_probit_eta(
      X_train = X_train_mat, y_train = Y_train, Z_train = Z_train - 0.5,
      num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
      general_params = general_params_nbcf_ols
    )
    tau_draws_nbcf_ols_test <- predict_linear_bcf_patched(fit_nbcf_ols, X = X_test_mat, Z = Z_test)$tau_hat
    y_std_nbcf_ols <- fit_nbcf_ols$model_params$outcome_scale
    alpha_draws_nbcf_ols <- as.vector(t(fit_nbcf_ols$alpha)) * y_std_nbcf_ols
    het_draws_nbcf_ols_fold <- t(t(tau_draws_nbcf_ols_test) - alpha_draws_nbcf_ols)

    # --- Semi-Parametric BART (Robust) ---
    cat(sprintf("  Fold %d: Fitting Semi-Parametric BART (Robust)...\n", k))
    fit_nbcf_robust <- bcf_linear_probit_eta(
      X_train = X_train_mat, y_train = Y_train, Z_train = Z_train - 0.5,
      num_gfr = 50, num_burnin = 1000, num_mcmc = 3000,
      general_params = general_params_nbcf_robust
    )
    tau_draws_nbcf_robust_test <- predict_linear_bcf_patched(fit_nbcf_robust, X = X_test_mat, Z = Z_test)$tau_hat
    y_std_nbcf_robust <- fit_nbcf_robust$model_params$outcome_scale
    alpha_draws_nbcf_robust <- as.vector(t(fit_nbcf_robust$alpha)) * y_std_nbcf_robust
    het_draws_nbcf_robust_fold <- t(t(tau_draws_nbcf_robust_test) - alpha_draws_nbcf_robust)
    
    # --- Standard BCF ---
    cat(sprintf("  Fold %d: Fitting Standard BCF...\n", k))
    fit_bcf <- bcf(
      X_train = X_train_mat, y_train = Y_train, Z_train = Z_train, 
      num_gfr = 25, num_burnin = 1000, num_mcmc = 4000, 
      general_params = general_params_bcf
    )
    tau_draws_bcf_test <- predict.bcfmodel(fit_bcf, X = X_test_mat, Z = Z_test)$tau_hat
    tau_draws_bcf_train <- predict.bcfmodel(fit_bcf, X = X_train_mat, Z = Z_train)$tau_hat
    ate_draws_bcf_train <- colMeans(tau_draws_bcf_train)
    het_draws_bcf_fold <- t(t(tau_draws_bcf_test) - ate_draws_bcf_train)
    
    # --- DR-Learner via SuperLearner ---
    cat(sprintf("  Fold %d: Fitting DR-Learner via SuperLearner...\n", k))
    df_dr <- data.frame(Y = Y_train, A = Z_train, X_train_mat)
    pi_hat <- mean(Z_train) 
    fold_ids_inner <- sample(rep(1:10, length.out = nrow(df_dr)))
    pseudo_Y <- numeric(nrow(df_dr))
    
    for (inner_k in 1:10) {
      train_idx_cv <- which(fold_ids_inner != inner_k)
      val_idx_cv <- which(fold_ids_inner == inner_k)
      
      X_train_k <- df_dr[train_idx_cv, setdiff(names(df_dr), "Y")]
      Y_train_k <- df_dr$Y[train_idx_cv]
      
      sl_mu <- SuperLearner(
        Y = Y_train_k, X = X_train_k, 
        SL.library = sl_library, family = gaussian(),
        cvControl = list(V = 10)
      )
      
      X_val_1 <- df_dr[val_idx_cv, setdiff(names(df_dr), "Y")]; X_val_1$A <- 1
      X_val_0 <- df_dr[val_idx_cv, setdiff(names(df_dr), "Y")]; X_val_0$A <- 0
      
      mu_1_hat <- predict(sl_mu, newdata = X_val_1)$pred
      mu_0_hat <- predict(sl_mu, newdata = X_val_0)$pred
      mu_A_hat <- ifelse(df_dr$A[val_idx_cv] == 1, mu_1_hat, mu_0_hat)
      
      pseudo_Y[val_idx_cv] <- (df_dr$A[val_idx_cv] - pi_hat) / (pi_hat * (1 - pi_hat)) * (df_dr$Y[val_idx_cv] - mu_A_hat) + (mu_1_hat - mu_0_hat)
    }
    
    X_train_df <- as.data.frame(X_train_mat)
    fit_dr <- SuperLearner(
      Y = pseudo_Y, X = X_train_df, 
      SL.library = sl_library, family = gaussian(),
      cvControl = list(V = 10)
    )
    
    X_test_df <- as.data.frame(X_test_mat)
    cate_dr <- as.vector(predict(fit_dr, newdata = X_test_df)$pred)
    
    cat(sprintf("\n[Worker] Completed Fold %d\n", k))
    
    list(
      test_idx = test_idx,
      cate_nbcf = rowMeans(tau_draws_nbcf_test),
      tau_draws_nbcf = het_draws_nbcf_fold,
      cate_nbcf_noshrink = rowMeans(tau_draws_nbcf_noshrink_test),
      tau_draws_nbcf_noshrink = het_draws_nbcf_noshrink_fold,
      cate_nbcf_ols = rowMeans(tau_draws_nbcf_ols_test),
      tau_draws_nbcf_ols = het_draws_nbcf_ols_fold,
      cate_nbcf_robust = rowMeans(tau_draws_nbcf_robust_test),
      tau_draws_nbcf_robust = het_draws_nbcf_robust_fold,
      cate_bcf = rowMeans(tau_draws_bcf_test),
      tau_draws_bcf = het_draws_bcf_fold,
      cate_dr = cate_dr,
      # For fold 1
      fit_nbcf = if (k == 1) fit_nbcf else NULL,
      fit_dr = if (k == 1) fit_dr else NULL,
      pseudo_Y = if (k == 1) pseudo_Y else NULL,
      X_train_mat = if (k == 1) X_train_mat else NULL,
      X_train_df = if (k == 1) X_train_df else NULL
    )
  }

  # Aggregate results
  for (k in 1:K_folds) {
    res <- fold_results[[k]]
    test_idx <- res$test_idx
    
    oos_cate_nbcf[test_idx] <- res$cate_nbcf
    oos_tau_draws_nbcf[[k]] <- list(idx = test_idx, het_draws = res$tau_draws_nbcf)
    
    oos_cate_nbcf_noshrink[test_idx] <- res$cate_nbcf_noshrink
    oos_tau_draws_nbcf_noshrink[[k]] <- list(idx = test_idx, het_draws = res$tau_draws_nbcf_noshrink)
    
    oos_cate_nbcf_ols[test_idx] <- res$cate_nbcf_ols
    oos_tau_draws_nbcf_ols[[k]] <- list(idx = test_idx, het_draws = res$tau_draws_nbcf_ols)
    
    oos_cate_nbcf_robust[test_idx] <- res$cate_nbcf_robust
    oos_tau_draws_nbcf_robust[[k]] <- list(idx = test_idx, het_draws = res$tau_draws_nbcf_robust)
    
    oos_cate_bcf[test_idx] <- res$cate_bcf
    oos_tau_draws_bcf[[k]] <- list(idx = test_idx, het_draws = res$tau_draws_bcf)
    
    oos_cate_dr[test_idx] <- res$cate_dr
    
    if (k == 1) {
      fit_nbcf_fold1 <- res$fit_nbcf
      fit_dr_fold1 <- res$fit_dr
      pseudo_Y_fold1 <- res$pseudo_Y
      X_train_mat_fold1 <- res$X_train_mat
      X_train_df_fold1 <- res$X_train_df
    }
  }

  cat("\n--- POOLED OOS INFERENCE RESULTS ---\n")
  cat(sprintf("Semi-Parametric BART Pooled OOS ATE: %.2f\n", mean(oos_cate_nbcf)))
  cat(sprintf("Semi-Parametric BART (No Shrinkage) Pooled OOS ATE: %.2f\n", mean(oos_cate_nbcf_noshrink)))
  cat(sprintf("Semi-Parametric BART (OLS) Pooled OOS ATE: %.2f\n", mean(oos_cate_nbcf_ols)))
  cat(sprintf("Semi-Parametric BART (Robust) Pooled OOS ATE: %.2f\n", mean(oos_cate_nbcf_robust)))
  cat(sprintf("Standard BCF Pooled OOS ATE: %.2f\n", mean(oos_cate_bcf)))
  cat(sprintf("DR-Learner Pooled OOS ATE: %.2f\n", mean(oos_cate_dr)))

  # Update variable names for downstream plotting
  cate_hat_nbcf_test <- oos_cate_nbcf
  cate_hat_nbcf_noshrink_test <- oos_cate_nbcf_noshrink
  cate_hat_nbcf_ols_test <- oos_cate_nbcf_ols
  cate_hat_nbcf_robust_test <- oos_cate_nbcf_robust
  cate_hat_bcf_test <- oos_cate_bcf
  cate_hat_dr_test <- oos_cate_dr

  cd40_test <- actg_sub$cd40
  Y_test <- Y_actg
  Z_test <- Z_actg


  # ==============================================================================
  # 6. GGPLOT VISUALIZATIONS (Global & Semi-Parametric) (Pooled OOS)
  # ==============================================================================
  colors_global <- c("Semi-Parametric BART" = "#008080", "Standard BCF" = "#d95f02", "DR-Learner" = "purple")
  colors_semi <- c("Semi-Parametric BART" = "#008080", "Semi-Parametric BART (No Shrinkage)" = "#e7298a", "Semi-Parametric BART (OLS)" = "#1b9e77", "Semi-Parametric BART (Robust)" = "#d95f02")

  # A1. Density Distributions (Global)
  density_data_global <- data.frame(
    CATE = c(cate_hat_nbcf_test, cate_hat_bcf_test, cate_hat_dr_test),
    Model = rep(c("Semi-Parametric BART", "Standard BCF", "DR-Learner"), each = length(cate_hat_nbcf_test))
  )

  p_dist_global <- ggplot(density_data_global, aes(x = CATE, fill = Model, color = Model)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
    scale_fill_manual(values = colors_global) +
    scale_color_manual(values = colors_global) +
    labs(title = "Distribution of Estimated CATEs (Global OOS)",
         x = "Treatment Effect (Difference in CD4 at 20 Weeks)", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # A2. Density Distributions (Semi-Parametric)
  density_data_semi <- data.frame(
    CATE = c(cate_hat_nbcf_test, cate_hat_nbcf_noshrink_test, cate_hat_nbcf_ols_test, cate_hat_nbcf_robust_test),
    Model = rep(c("Semi-Parametric BART", "Semi-Parametric BART (No Shrinkage)", "Semi-Parametric BART (OLS)", "Semi-Parametric BART (Robust)"), each = length(cate_hat_nbcf_test))
  )

  p_dist_semi <- ggplot(density_data_semi, aes(x = CATE, fill = Model, color = Model)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
    scale_fill_manual(values = colors_semi) +
    scale_color_manual(values = colors_semi) +
    labs(title = "Distribution of Estimated CATEs (Semi-Parametric OOS)",
         x = "Treatment Effect (Difference in CD4 at 20 Weeks)", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # B1. Treatment Effect Heterogeneity (Global)
  plot_data_het_global <- data.frame(
    CD4_Baseline = cd40_test,
    CATE_NBCF = cate_hat_nbcf_test,
    CATE_BCF = cate_hat_bcf_test,
    CATE_DR = cate_hat_dr_test
  )

  p_het_global <- ggplot(plot_data_het_global) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_BCF, color = "Standard BCF"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_DR, color = "DR-Learner"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_point(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), alpha = 0.15) +
    geom_point(aes(x = CD4_Baseline, y = CATE_BCF, color = "Standard BCF"), alpha = 0.15) +
    geom_point(aes(x = CD4_Baseline, y = CATE_DR, color = "DR-Learner"), alpha = 0.15) +
    scale_color_manual(values = colors_global) +
    labs(title = "Treatment Effect Heterogeneity (Global OOS)",
         subtitle = "Estimated CATE vs Baseline CD4 T-cell Count",
         x = "Baseline CD4 Count", y = "Estimated CATE", color = "Model") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # B2. Treatment Effect Heterogeneity (Semi-Parametric)
  plot_data_het_semi <- data.frame(
    CD4_Baseline = cd40_test,
    CATE_NBCF = cate_hat_nbcf_test,
    CATE_NBCF_NOSHRINK = cate_hat_nbcf_noshrink_test,
    CATE_NBCF_OLS = cate_hat_nbcf_ols_test,
    CATE_NBCF_ROBUST = cate_hat_nbcf_robust_test
  )

  p_het_semi <- ggplot(plot_data_het_semi) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF_NOSHRINK, color = "Semi-Parametric BART (No Shrinkage)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF_OLS, color = "Semi-Parametric BART (OLS)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF_ROBUST, color = "Semi-Parametric BART (Robust)"), method = "loess", se = FALSE, linewidth = 1.2) +
    geom_point(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), alpha = 0.15) +
    geom_point(aes(x = CD4_Baseline, y = CATE_NBCF_NOSHRINK, color = "Semi-Parametric BART (No Shrinkage)"), alpha = 0.15) +
    geom_point(aes(x = CD4_Baseline, y = CATE_NBCF_OLS, color = "Semi-Parametric BART (OLS)"), alpha = 0.15) +
    geom_point(aes(x = CD4_Baseline, y = CATE_NBCF_ROBUST, color = "Semi-Parametric BART (Robust)"), alpha = 0.15) +
    scale_color_manual(values = colors_semi) +
    labs(title = "Treatment Effect Heterogeneity (Semi-Parametric OOS)",
         subtitle = "Estimated CATE vs Baseline CD4 T-cell Count",
         x = "Baseline CD4 Count", y = "Estimated CATE", color = "Model") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # C1. CATE Correlation Heatmap (Global)
  cate_df_global <- data.frame(
    "Semi-Parametric BART" = cate_hat_nbcf_test,
    "Standard BCF" = cate_hat_bcf_test,
    "DR-Learner" = cate_hat_dr_test,
    check.names = FALSE
  )

  cor_matrix_global <- cor(cate_df_global)
  cor_data_global <- as.data.frame(as.table(cor_matrix_global))
  names(cor_data_global) <- c("Model1", "Model2", "Correlation")

  p_cor_heat_global <- ggplot(cor_data_global, aes(x = Model1, y = Model2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", Correlation)), 
              color = ifelse(cor_data_global$Correlation > 0.8, "white", "black"), 
              size = 6, fontface = "bold") +
    scale_fill_gradient2(low = "#cc0000", mid = "white", high = "#008080", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    labs(title = "CATE Correlation Heatmap (Global OOS)", subtitle = "Pairwise agreement between model treatment effect estimates") +
    theme_minimal(base_size = 14) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")

  # C2. CATE Correlation Heatmap (Semi-Parametric)
  cate_df_semi <- data.frame(
    "Semi-Parametric BART" = cate_hat_nbcf_test,
    "Semi-Parametric BART (No Shrinkage)" = cate_hat_nbcf_noshrink_test,
    "Semi-Parametric BART (OLS)" = cate_hat_nbcf_ols_test,
    "Semi-Parametric BART (Robust)" = cate_hat_nbcf_robust_test,
    check.names = FALSE
  )

  cor_matrix_semi <- cor(cate_df_semi)
  cor_data_semi <- as.data.frame(as.table(cor_matrix_semi))
  names(cor_data_semi) <- c("Model1", "Model2", "Correlation")

  p_cor_heat_semi <- ggplot(cor_data_semi, aes(x = Model1, y = Model2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", Correlation)), 
              color = ifelse(cor_data_semi$Correlation > 0.8, "white", "black"), 
              size = 6, fontface = "bold") +
    scale_fill_gradient2(low = "#cc0000", mid = "white", high = "#008080", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    labs(title = "CATE Correlation Heatmap (Semi-Parametric OOS)", subtitle = "Pairwise agreement between model treatment effect estimates") +
    theme_minimal(base_size = 14) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")

  # ==============================================================================
  # 7. IMPORTANCE ANALYSIS (Evaluated on Fold 1 Parameters)
  # ==============================================================================
  cat("\n--- IMPORTANCE ANALYSIS (Using Fold 1 Trained Parameters) ---\n")

  # A. Compute Shapley for Semi-Parametric BART
  if (!is.null(fit_nbcf_fold1$interaction_pairs)) {
    ipairs <- fit_nbcf_fold1$interaction_pairs
  } else {
    var_info <- results_actg$X_final_var_info
    is_candidate <- var_info$is_continuous | var_info$is_binary
    ipairs <- interaction_pairs(ncol(X_train_mat_fold1), is_candidate)
  }

  shapley_results <- compute_shapley_all(X = X_train_mat_fold1, beta_post = fit_nbcf_fold1$Beta, beta_int_post = fit_nbcf_fold1$Beta_int, ipairs = ipairs)
  p_shap_importance <- plot_shapley_importance_breakdown(shapley_results) + labs(title = "Semi-Parametric BART: Shapley Importance (Fold 1)")
  print(p_shap_importance)


  # ==============================================================================
  # 8 & 9 & 10. ROBUST UPLIFT CURVE EVALUATION (POOLED OOS)
  # ==============================================================================
  cat("\n--- GENERATING ROBUST UPLIFT CURVES ON POOLED OOS SET ---\n")

  get_eval_curves <- function(df, model_col, model_name) {
    df %>%
      arrange(desc(.data[[model_col]])) %>%
      mutate(
        cum_Z = cumsum(Z), cum_C = cumsum(1 - Z),
        cum_Y_Z = cumsum(Y * Z), cum_Y_C = cumsum(Y * (1 - Z)),
        cum_ATE = ifelse(cum_Z > 0 & cum_C > 0, (cum_Y_Z / cum_Z) - (cum_Y_C / cum_C), 0),
        uplift = cum_ATE * (cum_Z + cum_C), frac = row_number() / n(), Model = model_name
      ) %>%
      dplyr::select(frac, uplift, Model)
  }

  # Evaluate ONLY on Holdout Data
  eval_df_test <- data.frame(Y = Y_test, Z = Z_test, 
                             CATE_NBCF = cate_hat_nbcf_test, CATE_NBCF_NOSHRINK = cate_hat_nbcf_noshrink_test, CATE_NBCF_OLS = cate_hat_nbcf_ols_test, CATE_NBCF_ROBUST = cate_hat_nbcf_robust_test,
                             CATE_BCF = cate_hat_bcf_test, CATE_DR = cate_hat_dr_test)

  uplift_nbcf <- get_eval_curves(eval_df_test, "CATE_NBCF", "Semi-Parametric BART")
  uplift_nbcf_noshrink <- get_eval_curves(eval_df_test, "CATE_NBCF_NOSHRINK", "Semi-Parametric BART (No Shrinkage)")
  uplift_nbcf_ols <- get_eval_curves(eval_df_test, "CATE_NBCF_OLS", "Semi-Parametric BART (OLS)")
  uplift_nbcf_robust <- get_eval_curves(eval_df_test, "CATE_NBCF_ROBUST", "Semi-Parametric BART (Robust)")
  uplift_bcf <- get_eval_curves(eval_df_test, "CATE_BCF", "Standard BCF")
  uplift_dr <- get_eval_curves(eval_df_test, "CATE_DR", "DR-Learner")

  calc_net_auuc <- function(df, metric) { 
    raw_auc <- sum(diff(df$frac) * (head(df[[metric]], -1) + tail(df[[metric]], -1)) / 2)
    baseline_auc <- 0.5 * 1 * tail(df[[metric]], 1)
    return(raw_auc - baseline_auc)
  }

  auuc_nbcf <- calc_net_auuc(uplift_nbcf, "uplift")
  auuc_nbcf_noshrink <- calc_net_auuc(uplift_nbcf_noshrink, "uplift")
  auuc_nbcf_ols <- calc_net_auuc(uplift_nbcf_ols, "uplift")
  auuc_nbcf_robust <- calc_net_auuc(uplift_nbcf_robust, "uplift")
  auuc_bcf <- calc_net_auuc(uplift_bcf, "uplift")
  auuc_dr <- calc_net_auuc(uplift_dr, "uplift")

  final_uplift <- tail(uplift_nbcf$uplift, 1)

  # Global Uplift
  uplift_plot_data_global <- bind_rows(uplift_nbcf, uplift_bcf, uplift_dr)
  auuc_label_global <- sprintf("Net AUUC (Model - Random):\nSemi-Parametric BART = %.1f\nStandard BCF = %.1f\nDR-Learner = %.1f", auuc_nbcf, auuc_bcf, auuc_dr)

  p_uplift_smooth_global <- ggplot(uplift_plot_data_global, aes(x = frac, y = uplift, color = Model)) +
    geom_line(linewidth = 1.2) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = final_uplift, color = "black", linetype = "dashed", linewidth = 0.8) +
    annotate("label", x = 0.95, y = min(uplift_plot_data_global$uplift, na.rm = TRUE), 
             label = auuc_label_global, hjust = 1, vjust = 0, fontface = "bold", size = 4.5, color = "black", fill = "white", alpha = 0.85, label.size = NA) +
    labs(title = "Pooled Out-of-Sample Robust Uplift Curves (Global)", subtitle = "Evaluated on the 10-Fold Cross-Validated predictions",
         x = "Fraction of Population Treated", y = "Cumulative True Uplift") +
    scale_color_manual(values = colors_global) + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # Semi-Parametric Uplift
  uplift_plot_data_semi <- bind_rows(uplift_nbcf, uplift_nbcf_noshrink, uplift_nbcf_ols, uplift_nbcf_robust)
  auuc_label_semi <- sprintf("Net AUUC (Model - Random):\nSemi-Parametric = %.1f\nNo Shrinkage = %.1f\nOLS = %.1f\nRobust = %.1f", auuc_nbcf, auuc_nbcf_noshrink, auuc_nbcf_ols, auuc_nbcf_robust)

  p_uplift_smooth_semi <- ggplot(uplift_plot_data_semi, aes(x = frac, y = uplift, color = Model)) +
    geom_line(linewidth = 1.2) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = final_uplift, color = "black", linetype = "dashed", linewidth = 0.8) +
    annotate("label", x = 0.95, y = min(uplift_plot_data_semi$uplift, na.rm = TRUE), 
             label = auuc_label_semi, hjust = 1, vjust = 0, fontface = "bold", size = 4.5, color = "black", fill = "white", alpha = 0.85, label.size = NA) +
    labs(title = "Pooled Out-of-Sample Robust Uplift Curves (Semi-Parametric)", subtitle = "Evaluated on the 10-Fold Cross-Validated predictions",
         x = "Fraction of Population Treated", y = "Cumulative True Uplift") +
    scale_color_manual(values = colors_semi) + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")


  # ==============================================================================
  # 11 & 12. UNCERTAINTY QUANTIFICATION: TOLERANCE BANDS (POOLED OOS)
  # ==============================================================================
  cat("\n--- QUANTIFYING SIGNIFICANT HETEROGENEITY (POOLED OOS PATIENTS) ---\n")

  # Optimized function to compute credible intervals and means for each patient
  intfrompost1 <- function(post, conf = 0.05) {
    p_low <- conf / 2
    p_high <- 1 - (conf / 2)
    # Vectorized quantile per row
    qs <- t(apply(post, 1, quantile, probs = c(p_low, p_high), names = FALSE))
    mns <- rowMeans(post)
    return(cbind(qs[, 1], mns, qs[, 2]))
  }

  # Optimized function to find the global tolerance band with minimal memory footprint
  findconfglob <- function(post, tol = 0.05) {
    grid <- seq(0.0005, 0.05, by = 0.0005)
    
    # Pre-sort each row once to avoid repeated quantile calls!
    post_sorted <- t(apply(post, 1, sort, method = "quick"))
    S <- ncol(post_sorted)
    
    conftols <- sapply(grid, function(margconf) {
      p_low <- margconf / 2
      p_high <- 1 - (margconf / 2)
      
      idx_low <- max(1, round(S * p_low))
      idx_high <- min(S, round(S * p_high))
      
      int_lower <- post_sorted[, idx_low]
      int_upper <- post_sorted[, idx_high]
      
      difl <- (post >= int_lower) & (post <= int_upper)
      cs <- colMeans(difl)
      return(mean(cs >= 1 - tol))
    })
    
    wm <- which(conftols - 0.95 < 0)[1] - 1
    if (is.na(wm) || wm < 1) return(0.0005) 
    return(grid[wm])
  }

  build_heterogeneity_plot <- function(het_draws, model_name) {
    confglob <- findconfglob(het_draws, tol = 0.05)
    cred_int <- intfrompost1(het_draws, conf = 0.05)
    tol_int  <- intfrompost1(het_draws, conf = confglob)
    
    het_df <- data.frame(
      Mean = cred_int[, 2], Cred_Lower = cred_int[, 1], Cred_Upper = cred_int[, 3],
      Tol_Lower = tol_int[, 1], Tol_Upper = tol_int[, 3]
    ) %>%
      mutate(Significance = case_when(Tol_Lower > 0 ~ "Significant Positive (Benefit > Baseline)", Tol_Upper < 0 ~ "Significant Negative (Benefit < Baseline)", TRUE ~ "Not Significant")) %>% 
      arrange(Mean) %>% mutate(Ordered_ID = row_number())
    
    pct_sig <- sum(het_df$Significance != "Not Significant") / nrow(het_df) * 100
    
    ggplot(het_df, aes(x = Ordered_ID)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
      geom_ribbon(aes(ymin = Tol_Lower, ymax = Tol_Upper, fill = "Tolerance Interval"), alpha = 0.6) +
      geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper, fill = "95% Credible Interval"), alpha = 0.85) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      geom_point(data = filter(het_df, Significance != "Not Significant"), aes(y = Mean, color = Significance), size = 1.8, stroke = 0.5) +
      scale_fill_manual(values = c("Tolerance Interval" = "#D0D1E6", "95% Credible Interval" = "#74A9CF")) +
      scale_color_manual(values = c("Significant Positive (Benefit > Baseline)" = "#005a5a", "Significant Negative (Benefit < Baseline)" = "#cc0000")) +
      labs(title = paste0("Out-of-Sample Heterogeneity: ", model_name), subtitle = paste0(round(pct_sig, 1), "% of pooled OOS patients show significant deviation from baseline."), x = "Pooled OOS Patients (Ordered by Effect Size)", y = "Heterogeneous Effect (\u0394 CATE)") +
      theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  }

  # A. Semi-Parametric BART (Heterogeneity = CATE_test - alpha_train)
  het_draws_nbcf_test <- matrix(0, nrow = n_total, ncol = ncol(oos_tau_draws_nbcf[[1]]$het_draws))
  for (k in 1:K_folds) {
    het_draws_nbcf_test[oos_tau_draws_nbcf[[k]]$idx, ] <- oos_tau_draws_nbcf[[k]]$het_draws
  }
  p_het_nbcf <- build_heterogeneity_plot(het_draws_nbcf_test, "Semi-Parametric BART (CATE - \u03b1)")

  # B. Standard BCF (Heterogeneity = CATE_test - ATE_train)
  het_draws_bcf_test <- matrix(0, nrow = n_total, ncol = ncol(oos_tau_draws_bcf[[1]]$het_draws))
  for (k in 1:K_folds) {
    het_draws_bcf_test[oos_tau_draws_bcf[[k]]$idx, ] <- oos_tau_draws_bcf[[k]]$het_draws
  }
  p_het_bcf <- build_heterogeneity_plot(het_draws_bcf_test, "Standard BCF (CATE - ATE)")

  # C. Semi-Parametric BART (No Shrinkage)
  het_draws_nbcf_noshrink_test <- matrix(0, nrow = n_total, ncol = ncol(oos_tau_draws_nbcf_noshrink[[1]]$het_draws))
  for (k in 1:K_folds) {
    het_draws_nbcf_noshrink_test[oos_tau_draws_nbcf_noshrink[[k]]$idx, ] <- oos_tau_draws_nbcf_noshrink[[k]]$het_draws
  }
  p_het_nbcf_noshrink <- build_heterogeneity_plot(het_draws_nbcf_noshrink_test, "Semi-Parametric BART (No Shrinkage) (CATE - α)")

  # D. Semi-Parametric BART (OLS)
  het_draws_nbcf_ols_test <- matrix(0, nrow = n_total, ncol = ncol(oos_tau_draws_nbcf_ols[[1]]$het_draws))
  for (k in 1:K_folds) {
    het_draws_nbcf_ols_test[oos_tau_draws_nbcf_ols[[k]]$idx, ] <- oos_tau_draws_nbcf_ols[[k]]$het_draws
  }
  p_het_nbcf_ols <- build_heterogeneity_plot(het_draws_nbcf_ols_test, "Semi-Parametric BART (OLS) (CATE - α)")

  # E. Semi-Parametric BART (Robust)
  het_draws_nbcf_robust_test <- matrix(0, nrow = n_total, ncol = ncol(oos_tau_draws_nbcf_robust[[1]]$het_draws))
  for (k in 1:K_folds) {
    het_draws_nbcf_robust_test[oos_tau_draws_nbcf_robust[[k]]$idx, ] <- oos_tau_draws_nbcf_robust[[k]]$het_draws
  }
  p_het_nbcf_robust <- build_heterogeneity_plot(het_draws_nbcf_robust_test, "Semi-Parametric BART (Robust) (CATE - α)")


  # ==============================================================================
  # 13. SHAPLEY DISSECTION: MAIN VS INTERACTION EFFECTS (FEATURE-LEVEL) (Evaluated on Fold 1)
  # ==============================================================================
  cat("\n--- RUNNING SHAPLEY DISSECTION ANALYSIS FOR 'cd40' (Using Fold 1 Trained Parameters) ---\n")

  # A. Prepare Data & Identify Target Feature
  feature_to_dissect <- "cd40" # Focus on Baseline CD4 count

  # Extract dissected posterior
  dissected_post <- get_shapley_posterior_dissected(
    feature_name = feature_to_dissect, 
    indices = 1:nrow(X_train_mat_fold1), 
    X = X_train_mat_fold1, 
    beta_post = fit_nbcf_fold1$Beta, 
    beta_int_post = fit_nbcf_fold1$Beta_int, 
    ipairs = ipairs
  )

  # B. Calculate Medians
  plot_df_dissect <- data.frame(
    Feature_Value = X_train_mat_fold1[, feature_to_dissect],
    Main_Med = apply(dissected_post$main, 1, median),
    Int_Med = apply(dissected_post$interaction, 1, median)
  ) %>%
    mutate(
      Total_Shapley = Main_Med + Int_Med,
      Int_Ratio = abs(Int_Med) / (abs(Main_Med) + abs(Int_Med) + 1e-8)
    )

  cat("  Generating Synergy vs. Antagonism Scatter Plot...\n")
  limit_val <- max(abs(c(plot_df_dissect$Main_Med, plot_df_dissect$Int_Med)))

  p_synergy <- ggplot(plot_df_dissect, aes(x = Main_Med, y = Int_Med, color = Feature_Value)) +
    geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dotted") +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_viridis_c(name = "Standardized\nBaseline CD4") +
    annotate("text", x = limit_val*0.8, y = limit_val*0.8, label = "Synergy\n(Amplifies)", fontface = "italic") +
    annotate("text", x = -limit_val*0.8, y = limit_val*0.8, label = "Antagonism\n(Suppresses)", fontface = "italic") +
    coord_fixed(xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) +
    labs(
      title = paste("Synergy vs Antagonism:", feature_to_dissect),
      subtitle = "Dashed line (y = -x) indicates perfect cancellation of main effect by interactions",
      x = "Main Effect Contribution",
      y = "Total Interaction Contribution"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  cat("  Generating Fractional Contribution Profile...\n")
  p_fractional <- ggplot(plot_df_dissect, aes(x = Feature_Value, y = Int_Ratio)) +
    geom_point(alpha = 0.3, color = "#2c7bb6") +
    geom_smooth(method = "loess", span = 0.3, color = "#d7191c", se = FALSE, linewidth = 1.2) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      title = paste("Interaction Dominance Profile:", feature_to_dissect),
      subtitle = "Proportion of total Shapley value driven by interaction terms",
      x = paste("Standardized", feature_to_dissect),
      y = "Interaction Ratio (|Int| / (|Main| + |Int|))"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

  cat("  Generating Dual-Density Uncertainty Plots...\n")
  quant_vals <- quantile(plot_df_dissect$Feature_Value, probs = c(0.1, 0.5, 0.9))
  idx_10 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[1]))
  idx_50 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[2]))
  idx_90 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[3]))

  post_df <- data.frame(
    Draw = 1:ncol(dissected_post$main),
    P_10_Main = dissected_post$main[idx_10, ],
    P_10_Int  = dissected_post$interaction[idx_10, ],
    P_50_Main = dissected_post$main[idx_50, ],
    P_50_Int  = dissected_post$interaction[idx_50, ],
    P_90_Main = dissected_post$main[idx_90, ],
    P_90_Int  = dissected_post$interaction[idx_90, ]
  ) %>%
    tidyr::pivot_longer(
      cols = -Draw,
      names_to = c("Patient", "Effect"),
      names_pattern = "P_(.*)_(.*)",
      values_to = "Value"
    ) %>%
    mutate(
      Patient_Label = case_when(
        Patient == "10" ~ sprintf("Patient A (10th Pctl CD4: %.2f)", plot_df_dissect$Feature_Value[idx_10]),
        Patient == "50" ~ sprintf("Patient B (Median CD4: %.2f)", plot_df_dissect$Feature_Value[idx_50]),
        Patient == "90" ~ sprintf("Patient C (90th Pctl CD4: %.2f)", plot_df_dissect$Feature_Value[idx_90])
      ),
      Effect = factor(Effect, levels = c("Main", "Int"), labels = c("Main Effect", "Interaction Effect"))
    )

  p_uncertainty <- ggplot(post_df, aes(x = Value, fill = Effect, color = Effect)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ Patient_Label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Main Effect" = "#4d4d4d", "Interaction Effect" = "#e08214")) +
    scale_color_manual(values = c("Main Effect" = "#1a1a1a", "Interaction Effect" = "#b35806")) +
    labs(
      title = paste("Posterior Uncertainty: Main vs. Interaction for", feature_to_dissect),
      subtitle = "Comparing density spreads across MCMC draws for prototypical patients",
      x = "Shapley Value Contribution",
      y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )

  # ==============================================================================
  # 14. SAVE PLOTS TO LOCAL DIRECTORY
  # ==============================================================================
  plot_dir <- sprintf("plots/ACTG_new_cross_seed_%d", current_seed)
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  plot_w <- 8; plot_h <- 6; plot_dpi <- 300

  # Global Plots
  ggsave(file.path(plot_dir, "01a_Global_CATE_Density_Distribution.png"), p_dist_global, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "02a_Global_Treatment_Heterogeneity.png"), p_het_global, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "03a_Global_CATE_Correlation_Heatmap.png"), p_cor_heat_global, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "04a_Global_Robust_Uplift.png"), p_uplift_smooth_global, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")

  # Semi-Parametric Plots
  ggsave(file.path(plot_dir, "01b_SemiParametric_CATE_Density_Distribution.png"), p_dist_semi, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "02b_SemiParametric_Treatment_Heterogeneity.png"), p_het_semi, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "03b_SemiParametric_CATE_Correlation_Heatmap.png"), p_cor_heat_semi, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "04b_SemiParametric_Robust_Uplift.png"), p_uplift_smooth_semi, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")

  # Other
  ggsave(file.path(plot_dir, "05_SemiParametric_Shapley.png"), p_shap_importance, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "06a_Heterogeneity_Bands_SemiParametric.png"), p_het_nbcf, width = 10, height = 7, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "06b_Heterogeneity_Bands_SemiParametric_NoShrinkage.png"), p_het_nbcf_noshrink, width = 10, height = 7, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "06c_Heterogeneity_Bands_SemiParametric_OLS.png"), p_het_nbcf_ols, width = 10, height = 7, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "06d_Heterogeneity_Bands_SemiParametric_Robust.png"), p_het_nbcf_robust, width = 10, height = 7, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "06e_Heterogeneity_Bands_Standard_BCF.png"), p_het_bcf, width = 10, height = 7, dpi = plot_dpi, bg = "white")
  
  # Shapley Dissection
  ggsave(file.path(plot_dir, "07a_Shapley_Synergy_Antagonism_New.png"), p_synergy, width = 8, height = 8, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "07b_Shapley_Fractional_Contribution_New.png"), p_fractional, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
  ggsave(file.path(plot_dir, "07c_Shapley_DualDensity_Uncertainty_New.png"), p_uncertainty, width = 8, height = 10, dpi = plot_dpi, bg = "white")

  cat(sprintf("All out-of-sample plots successfully saved to: %s\n", plot_dir))

} # End of multiple seeds loop

# Stop parallel cluster
stopCluster(cl)
cat("\n[INFO] Parallel Cluster Stopped. All tasks completed!\n")
