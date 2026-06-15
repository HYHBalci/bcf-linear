# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(stochtree)   # Semi-Parametric BCF method
# Vansteelandt DR-Learner Ensemble Libraries:
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(ranger)
library(coda)
library(rpart)
library(rpart.plot)

if (!require(openxlsx)) install.packages("openxlsx", repos="http://cran.us.r-project.org")
library(openxlsx)
if (!require(progress)) install.packages("progress", repos="http://cran.us.r-project.org")
library(progress)
if (!require(zip)) install.packages("zip", repos="http://cran.us.r-project.org")
library(zip)

source('R/shapley_aux.R')
source('R/old_linear_linear.R')  # Load custom linear model
source('R/utility_functions.R')

# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PATCHED PREDICTION
# ==============================================================================

standardize_X_by_index_new <- function(X_initial, process_data = TRUE, interaction_rule = c("continuous", "continuous_or_binary", "all", "none"), cat_coding_method = c("sum", "difference")) {
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
    
    X_final_var_info <- data.frame(is_continuous = intermediate_map$is_continuous[matched_rows], is_binary = intermediate_map$is_binary[matched_rows], col_name_final = colnames(X_final), stringsAsFactors = FALSE)
    non_continous_idx_cpp <- which(!X_final_var_info$is_continuous) - 1
  } else {
    X_final <- as.matrix(X_df)
    X_final_var_info <- data.frame(is_continuous = TRUE, is_binary = FALSE, col_name_final = colnames(X_final))
    non_continous_idx_cpp <- numeric(0)
  }
  
  p_int_calculated <- 0; n_final_cols <- ncol(X_final)
  if (n_final_cols > 1 && interaction_rule != "none") {
    is_candidate <- switch(interaction_rule, "all" = rep(TRUE, n_final_cols), "continuous_or_binary" = (X_final_var_info$is_continuous | X_final_var_info$is_binary), "continuous" = X_final_var_info$is_continuous)
    is_candidate[is.na(is_candidate)] <- FALSE
    for (i in 1:(n_final_cols - 1)) { for (j in (i + 1):n_final_cols) { if (is_candidate[i] || is_candidate[j]) p_int_calculated <- p_int_calculated + 1 } }
  }
  return(list(X_final = X_final, p_int = p_int_calculated, non_continous_idx_cpp = non_continous_idx_cpp, X_final_var_info = X_final_var_info))
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

# Variable Classification Function (Task 2)
classify_variables <- function(X_comp, fit_nbcf, cate_hat_dr, y_train, z_train, df_dr) {
  var_names <- colnames(X_comp)
  n_vars <- length(var_names)
  
  # BCF Importances via Inclusion Proportions
  beta_samples <- matrix(aperm(fit_nbcf$Beta, c(2, 1, 3)), nrow = dim(fit_nbcf$Beta)[1]*dim(fit_nbcf$Beta)[2])
  beta_inclusion <- colMeans(beta_samples != 0)
  
  if (!is.null(fit_nbcf$Beta_int) && dim(fit_nbcf$Beta_int)[3] > 0) {
    beta_int_samples <- matrix(aperm(fit_nbcf$Beta_int, c(2, 1, 3)), nrow = dim(fit_nbcf$Beta_int)[1]*dim(fit_nbcf$Beta_int)[2])
    int_pairs <- fit_nbcf$interaction_pairs
    var_predictive_score <- numeric(n_vars)
    for(v in 1:n_vars) {
      idx <- which(int_pairs[1,] == v | int_pairs[2,] == v)
      if(length(idx) > 0) {
        var_predictive_score[v] <- mean(rowSums(beta_int_samples[, idx, drop=FALSE] != 0) > 0)
      }
    }
  } else {
    var_predictive_score <- numeric(n_vars)
  }
  
  # DR-Learner Importances
  rf_prog <- ranger(Y ~ ., data = df_dr[, c("Y", var_names)], importance = "permutation")
  rf_prog_imp <- rf_prog$variable.importance
  
  df_pred <- data.frame(CATE = cate_hat_dr, df_dr[, var_names, drop=FALSE])
  rf_pred <- ranger(CATE ~ ., data = df_pred, importance = "permutation")
  rf_pred_imp <- rf_pred$variable.importance
  
  # Normalize scores
  prog_score <- (beta_inclusion / max(beta_inclusion + 1e-6) + rf_prog_imp / max(rf_prog_imp + 1e-6)) / 2
  pred_score <- (var_predictive_score / max(var_predictive_score + 1e-6) + rf_pred_imp / max(rf_pred_imp + 1e-6)) / 2
  
  # Thresholds (e.g. top 25% or > 0.1)
  thresh_prog <- max(0.1, quantile(prog_score, 0.75))
  thresh_pred <- max(0.1, quantile(pred_score, 0.75))
  
  is_prog <- prog_score >= thresh_prog
  is_pred <- pred_score >= thresh_pred
  
  classifications <- rep("neither", n_vars)
  classifications[is_prog & is_pred] <- "both"
  classifications[is_prog & !is_pred] <- "prognostic"
  classifications[!is_prog & is_pred] <- "predictive"
  
  res <- data.frame(Variable = var_names, Classification = classifications)
  return(res)
}

# Global Score Test for Treatment Effect Heterogeneity
# Based on Sechidis et al. "Comparing methods to assess treatment effect heterogeneity in general parametric regression models"
# Uses the score residual of the centered treatment indicator and coin::independence_test
sechidis_score_test <- function(Y, X, W, is_binary) {
  if (!requireNamespace("coin", quietly = TRUE)) {
    install.packages("coin", repos="http://cran.us.r-project.org")
  }
  
  # 1. Center the treatment indicator (Section 3.4)
  W_centered <- W - 0.5
  
  # 2. Fit base model with covariates as prognostic effects
  df <- data.frame(Y = Y, W_centered = W_centered, X)
  fam <- if(is_binary) binomial() else gaussian()
  mod <- glm(Y ~ ., data = df, family = fam)
  
  # 3. Calculate score residuals (Equation 3 / generalized version)
  # For canonical links, the score derivative w.r.t treatment parameter is (Y - \hat{\mu}) * W_centered
  mu_hat <- predict(mod, type = "response")
  s_i <- (Y - mu_hat) * W_centered
  
  # 4. Global independence permutation test
  test_data <- data.frame(s_i = s_i, X)
  form <- as.formula(paste("s_i ~", paste(colnames(X), collapse = " + ")))
  
  # Using the quadratic test statistic as evaluated in the paper, with asymptotic distribution
  test_res <- coin::independence_test(form, data = test_data, teststat = "quadratic")
  p_val <- coin::pvalue(test_res)
  
  return(p_val)
}

# Helper functions for AUUC computation
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
calc_auc <- function(df, metric) { sum(diff(df$frac) * (head(df[[metric]], -1) + tail(df[[metric]], -1)) / 2) }

# ==============================================================================
# 3. DIRECTORY SETUP & BATCH LOOP
# ==============================================================================

input_dir <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/Data/Datasets/Blinded_data_sets"
output_dir <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/Output"
plot_base_dir <- file.path(output_dir, "plots")
sub_datasets_dir <- file.path(output_dir, "submission_datasets")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_base_dir)) dir.create(plot_base_dir, recursive = TRUE)
if (!dir.exists(sub_datasets_dir)) dir.create(sub_datasets_dir, recursive = TRUE)

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# If no files found, use the sample directory from the original script logic
if (length(csv_files) == 0) {
  input_dir_alternative <- dirname("C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/Data/Datasets/Blinded_data_sets/dataset30_ 2318 .csv")
  csv_files <- list.files(input_dir_alternative, pattern = "\\.csv$", full.names = TRUE)
}

results_list <- list()
var_class_list <- list()

cat("\n--- STARTING BATCH PROCESSING ---\n")
pb <- progress_bar$new(format = "  Processing [:bar] :percent ETA: :eta | :current/:total datasets", total = length(csv_files), clear = FALSE, width = 80)

for (file_path in csv_files) {
  file_name <- basename(file_path)
  dataset_id <- tools::file_path_sans_ext(file_name)
  plot_dir <- file.path(plot_base_dir, dataset_id)
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  tryCatch({
    comp_data <- read.csv(file_path)
    
    Y_comp <- comp_data$Y
    Z_comp <- as.numeric(comp_data$W)
    X_comp_raw <- comp_data %>% dplyr::select(matches("^[xX][0-9]+$"))
    
    is_binary <- length(unique(Y_comp[!is.na(Y_comp)])) <= 2
    sl_family <- if (is_binary) binomial() else gaussian()
    
    current_interaction_rule <- if (ncol(X_comp_raw) > 15) "none" else "continuous_or_binary"
    
    results_comp <- standardize_X_by_index_new(X_initial = X_comp_raw, process_data = TRUE, interaction_rule = current_interaction_rule, cat_coding_method = "sum")
    X_train_mat <- results_comp$X_final
    
    general_params_default <- list(
      cutpoint_grid_size = 100, standardize = TRUE, sample_sigma2_global = !is_binary, sigma2_global_init = 1, 
      sigma2_global_shape = 1, sigma2_global_scale = 0.001, variable_weights = NULL, propensity_covariate = "mu", 
      adaptive_coding = FALSE, control_coding_init = -0.5, treated_coding_init = 0.5, rfx_prior_var = NULL, 
      random_seed = 1, keep_burnin = FALSE, keep_gfr = FALSE, keep_every = 1, num_chains = 1, verbose = FALSE, 
      sample_global_prior = "half-cauchy", unlink = TRUE, propensity_seperate = "none", gibbs = TRUE, step_out = 0.5, 
      max_steps = 150, save_output = FALSE, probit_outcome_model = is_binary, interaction_rule = current_interaction_rule, 
      standardize_cov = FALSE, simple_prior = FALSE, save_partial_residual = FALSE, regularize_ATE = FALSE,
      sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1
    )
    
    # Create an inner progress bar for the dataset
    pb_inner <- progress_bar$new(
      format = "    Dataset steps [:bar] :percent ETA: :eta",
      total = 5, clear = TRUE, width = 60
    )
    
    # BCF Models (Probit for binary, Linear for continuous)
    pb_inner$message(sprintf("[%s] Fitting BCF Horseshoe model...", Sys.time()))
    pb_inner$tick(0)
    fit_nbcf <- bcf_linear_probit(X_train = X_train_mat, y_train = Y_comp, Z_train = Z_comp, num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = general_params_default)
    
    # Standard BCF
    pb_inner$message(sprintf("[%s] Fitting Standard BCF model...", Sys.time()))
    pb_inner$tick(1)
    general_params_standard_bcf <- general_params_default
    general_params_standard_bcf$propensity_covariate <- "none"
    general_params_standard_bcf$adaptive_coding <- TRUE
    
    fit_standard_bcf <- bcf(
      X_train = as.matrix(X_train_mat),
      y_train = Y_comp,
      Z_train = Z_comp,
      propensity_train = NULL,
      num_gfr = 50, 
      num_burnin = 1000, 
      num_mcmc = 3000,
      general_params = general_params_standard_bcf
    )
    
    # DR-Learner
    pb_inner$message(sprintf("[%s] Fitting DR-Learner...", Sys.time()))
    pb_inner$tick(1)
    sl_library <- c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.ranger")
    if (ncol(X_comp_raw) <= 15) {
      sl_library <- c(sl_library, "SL.gam")
    }
    pi_hat <- mean(Z_comp) 
    df_dr <- data.frame(Y = Y_comp, A = Z_comp, X_train_mat)
    
    K_folds <- 10
    fold_ids <- sample(rep(1:K_folds, length.out = nrow(df_dr)))
    pseudo_Y <- numeric(nrow(df_dr))
    
    dr_mu_coefs <- matrix(NA, nrow=K_folds, ncol=length(sl_library))
    colnames(dr_mu_coefs) <- sl_library
    
    for (k in 1:K_folds) {
      train_idx <- which(fold_ids != k); val_idx <- which(fold_ids == k)
      X_train_k <- df_dr[train_idx, setdiff(names(df_dr), "Y")]; Y_train_k <- df_dr$Y[train_idx]
      
      sl_mu <- SuperLearner(Y = Y_train_k, X = X_train_k, SL.library = sl_library, family = sl_family, cvControl = list(V = 10))
      dr_mu_coefs[k, ] <- sl_mu$coef
      
      X_val_1 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_1$A <- 1
      X_val_0 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_0$A <- 0
      
      mu_1_hat <- predict(sl_mu, newdata = X_val_1)$pred
      mu_0_hat <- predict(sl_mu, newdata = X_val_0)$pred
      mu_A_hat <- ifelse(df_dr$A[val_idx] == 1, mu_1_hat, mu_0_hat)
      pseudo_Y[val_idx] <- (df_dr$A[val_idx] - pi_hat) / (pi_hat * (1 - pi_hat)) * (df_dr$Y[val_idx] - mu_A_hat) + (mu_1_hat - mu_0_hat)
    }
    X_train_df <- as.data.frame(X_train_mat)
    fit_dr <- SuperLearner(Y = pseudo_Y, X = X_train_df, SL.library = sl_library, family = gaussian(), cvControl = list(V = 10))
    cate_hat_dr <- as.vector(fit_dr$SL.predict)
    dr_tau_weights <- fit_dr$coef
    dr_mu_weights <- colMeans(dr_mu_coefs)
    
    # Inference
    tau_draws_nbcf <- predict_linear_bcf_patched(fit_nbcf, X = X_train_mat, Z = as.matrix(Z_comp))$tau_hat
    cate_hat_nbcf <- rowMeans(tau_draws_nbcf)
    ate_draws_nbcf <- colMeans(tau_draws_nbcf)
    ate_res_nbcf <- mean(ate_draws_nbcf)
    
    tau_draws_standard_bcf <- fit_standard_bcf$tau_hat_train
    cate_hat_standard_bcf <- rowMeans(tau_draws_standard_bcf)
    ate_draws_standard_bcf <- colMeans(tau_draws_standard_bcf)
    ate_res_standard_bcf <- mean(ate_draws_standard_bcf)
    
    # Fit custom Old Linear model (Source script assumed loaded)
    pb_inner$message(sprintf("[%s] Fitting Old Linear Model...", Sys.time()))
    pb_inner$tick(1)
    n_iter <- 2000; burn_in <- 1000
    fit_old <- fit_grouped_horseshoes_R(y_vec = Y_comp, X_mat = X_train_mat, Z_vec = Z_comp, family = ifelse(is_binary, "binomial", "gaussian"), n_iter = n_iter, burn_in = burn_in, num_chains = 1, propensity_as_covariate = FALSE, standardize_cov = FALSE, interaction_rule = current_interaction_rule, cat_coding_method = "difference", sample_sigma = TRUE, sigma_init = 1)
    
    aleph_hat   <- mean(fit_old$aleph)
    gamma_hat   <- colMeans(fit_old$gamma)
    gamma_int_hat <- if (!is.null(fit_old$gamma_int)) colMeans(fit_old$gamma_int) else rep(0, ncol(X_train_mat))
    handled_old <- standardize_X_by_index_new(X_train_mat, process_data = FALSE, interaction_rule = current_interaction_rule, cat_coding_method = "difference")
    X_std   <- handled_old$X_final
    X_int_std <- if (!is.null(handled_old$X_int_final)) handled_old$X_int_final else matrix(0, nrow = nrow(X_std), ncol = 0)
    
    # Calculate CATE point estimate
    if (is.null(fit_old$gamma_int)) {
      treat_part_old <- aleph_hat + X_std %*% gamma_hat
    } else {
      treat_part_old <- aleph_hat + X_std %*% gamma_hat + X_int_std %*% gamma_int_hat
    }
    cate_hat_old <- as.numeric(treat_part_old)
    
    # Calculate posterior draws of CATE
    if (is.null(fit_old$gamma_int)) {
      tau_draws_old <- sweep(X_std %*% t(fit_old$gamma), 2, fit_old$aleph, "+")
    } else {
      tau_draws_old <- sweep(X_std %*% t(fit_old$gamma) + X_int_std %*% t(fit_old$gamma_int), 2, fit_old$aleph, "+")
    }
    ate_draws_old <- colMeans(tau_draws_old)
    ate_res_old <- mean(cate_hat_old)
    
    ate_res_dr <- mean(pseudo_Y)
    
    # Save Posteriors
    posteriors_list <- list(
      beta = fit_nbcf$Beta,
      beta_int = fit_nbcf$Beta_int,
      alpha = fit_nbcf$alpha,
      gamma = fit_old$gamma,
      gamma_int = fit_old$gamma_int,
      aleph = fit_old$aleph,
      tau_standard_bcf = tau_draws_standard_bcf
    )
    saveRDS(posteriors_list, file = file.path(plot_dir, paste0("posteriors_", dataset_id, ".rds")))
    
    # Task 2: Variable Classification
    var_classes <- classify_variables(X_train_mat, fit_nbcf, cate_hat_dr, Y_comp, Z_comp, df_dr)
    var_classes$Dataset <- dataset_id
    var_class_list[[dataset_id]] <- var_classes
    predictive_vars <- var_classes$Variable[var_classes$Classification %in% c("predictive", "both")]
    
    # Task 3 & 4: Subgroup Extraction
    get_subgroup <- function(cate_preds, model_name) {
      tree_data <- data.frame(CATE = cate_preds, X_comp_raw)
      subgroup_tree <- rpart(CATE ~ ., data = tree_data, method = "anova", control = rpart.control(cp = 0.01, maxdepth = 3))
      png(filename = file.path(plot_dir, paste0("08_Subgroup_Tree_", gsub("[ ()]", "_", model_name), ".png")), width = 2400, height = 1800, res = 300)
      rpart.plot(subgroup_tree, type = 4, extra = 1, roundint = FALSE, main = paste("Subgroup Rule:", model_name), box.palette = c("tomato", "white", "#008080"), shadow.col = "gray", nn = TRUE)
      dev.off()
      tree_data$Leaf_Node <- as.factor(subgroup_tree$where)
      best_node <- (tree_data %>% group_by(Leaf_Node) %>% summarize(m = mean(CATE)) %>% arrange(desc(m)))$Leaf_Node[1]
      return(ifelse(tree_data$Leaf_Node == best_node, 1, 0))
    }
    
    sub_nbcf <- get_subgroup(cate_hat_nbcf, "BCF (Horseshoe)")
    sub_standard_bcf <- get_subgroup(cate_hat_standard_bcf, "Standard BCF")
    sub_dr <- get_subgroup(cate_hat_dr, "Vansteelandt DR-Learner")
    sub_old <- get_subgroup(cate_hat_old, "Old Linear Model")
    
    # Majority Votes
    vote_bcf <- sub_nbcf + sub_standard_bcf + sub_old
    vote_dr <- sub_dr
    
    final_S_bcf <- ifelse(vote_bcf >= 2, 1, 0)
    final_S_dr <- vote_dr
    final_S <- final_S_bcf # Main export uses BCF ensemble
    
    # Export datasets
    out_df <- comp_data
    out_df$S_BCF <- final_S_bcf
    out_df$S_DR <- final_S_dr
    out_df$S <- final_S_bcf # Backward compatibility
    
    if ("ptid" %in% names(out_df)) out_df <- out_df[order(out_df$ptid), ]
    write.csv(out_df, file.path(sub_datasets_dir, file_name), row.names = FALSE)
    
    subgroup_prop <- mean(final_S)
    subgroup_ate <- mean(cate_hat_nbcf[final_S == 1])
    
    # Visualizations
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.3f", Correlation)), color = ifelse(cor_data$Correlation > 0.8, "white", "black"), size = 5, fontface = "bold") +
      scale_fill_gradient2(low = "#cc0000", mid = "white", high = "#008080", midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation") +
      labs(title = paste("CATE Correlation Heatmap:", dataset_id)) + theme_minimal(base_size=14) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1, face="bold"), axis.text.y=element_text(face="bold"))
    
    ggsave(file.path(plot_dir, paste0("04_CATE_Correlation_Heatmap_", dataset_id, ".png")), plot = p_cor_heat, width = 8, height = 6, dpi = 150, bg = "white")
    
    # 05_Robust AUUC Uplift Curves
    eval_df <- data.frame(
      Y = Y_comp,
      Z = Z_comp,
      Y = Y_comp, Z = Z_comp,
      CATE_NBCF = cate_hat_nbcf,
      CATE_Standard_BCF = cate_hat_standard_bcf,
      CATE_DR = cate_hat_dr,
      CATE_Old = cate_hat_old
    )
    
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
    
    uplift_nbcf <- get_eval_curves(eval_df, "CATE_NBCF", "BCF (Horseshoe)")
    uplift_standard_bcf <- get_eval_curves(eval_df, "CATE_Standard_BCF", "Standard BCF")
    uplift_dr <- get_eval_curves(eval_df, "CATE_DR", "Vansteelandt DR")
    uplift_old <- get_eval_curves(eval_df, "CATE_Old", "Old Linear")
    
    uplift_plot_data <- bind_rows(uplift_nbcf, uplift_standard_bcf, uplift_dr, uplift_old)
    
    auuc_nbcf <- calc_auc(uplift_nbcf, "uplift")
    auuc_standard_bcf <- calc_auc(uplift_standard_bcf, "uplift")
    auuc_dr <- calc_auc(uplift_dr, "uplift")
    auuc_old <- calc_auc(uplift_old, "uplift")
    
    auuc_label <- sprintf("AUUC:\nBCF (Horseshoe) = %.1f\nStandard BCF = %.1f\nVansteelandt DR = %.1f\nOld Linear = %.1f", auuc_nbcf, auuc_standard_bcf, auuc_dr, auuc_old)
    final_uplift <- tail(uplift_nbcf$uplift, 1)
    
    p_uplift_smooth <- ggplot(uplift_plot_data, aes(x = frac, y = uplift, color = Model)) +
      geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
      annotate("segment", x = 0, y = 0, xend = 1, yend = final_uplift, color = "black", linetype = "dashed", linewidth = 0.8) +
      annotate("label", x = 0.95, y = min(uplift_plot_data$uplift, na.rm = TRUE), label = auuc_label, hjust = 1, vjust = 0, fontface = "bold", size = 3.5, color = "black", fill = "white", alpha = 0.85, label.size = NA) +
      scale_color_manual(values = color_palette) +
      labs(title = paste("Robust AUUC Uplift:", dataset_id), x = "Fraction of Population Treated", y = "Cumulative True Uplift") +
      theme_minimal() + theme(legend.position = "bottom")
      
    ggsave(file.path(plot_dir, paste0("05_AUUC_Uplift_", dataset_id, ".png")), plot = p_uplift_smooth, width = 8, height = 6, dpi = 150, bg = "white")
    
    # Global Score Test for TEH (already defined earlier)
    score_p_val <- sechidis_score_test(Y_comp, X_train_mat, Z_comp, is_binary)
    score_test_significant <- !is.na(score_p_val) && (score_p_val < 0.05)

    # ---------- Interval Analyses & Permutation Tests ----------
    pb_inner$message(sprintf("[%s] Running Interval Tests and AUUC Permutations...", Sys.time()))
    pb_inner$tick(1)
    het_draws_nbcf <- t(t(tau_draws_nbcf) - ate_draws_nbcf)
    het_draws_standard_bcf <- t(t(tau_draws_standard_bcf) - ate_draws_standard_bcf)
    het_draws_old <- t(t(tau_draws_old) - ate_draws_old)
    
    sig_nbcf <- check_significant_individuals(het_draws_nbcf)
    sig_standard_bcf <- check_significant_individuals(het_draws_standard_bcf)
    sig_old <- check_significant_individuals(het_draws_old)
    
    perm_bcf_shoe <- auuc_permutation_test(eval_df, "CATE_NBCF", B=1000)
    perm_standard_bcf <- auuc_permutation_test(eval_df, "CATE_Standard_BCF", B=1000)
    perm_dr <- auuc_permutation_test(eval_df, "CATE_DR", B=1000)
    perm_old <- auuc_permutation_test(eval_df, "CATE_Old", B=1000)
    
    # Export summary text file
    summary_file <- file.path(plot_dir, paste0("results_summary_", dataset_id, ".txt"))
    sink(summary_file)
    cat("==================================================================\n")
    cat(sprintf("DATASET: %s\n", dataset_id))
    if (current_interaction_rule == "none") {
      cat("NOTE: Covariate space > 15. Using ONLY MAIN EFFECTS for tau(x) (interaction_rule = 'none').\n")
    } else {
      cat(sprintf("Interaction Rule: %s\n", current_interaction_rule))
    }
    cat("==================================================================\n\n")

    cat("--- TOP COVARIATES ---\n")
    cat("Predictive/Both:\n")
    print(var_classes$Variable[var_classes$Classification %in% c("predictive", "both")])
    cat("Prognostic:\n")
    print(var_classes$Variable[var_classes$Classification %in% c("prognostic")])
    cat("\n")

    cat("--- DR-LEARNER SUPERLEARNER WEIGHTS ---\n")
    cat("Tau(x) Model Weights (Final):\n")
    print(dr_tau_weights)
    cat("\nMu(x) Model Weights (Averaged across 10 folds):\n")
    print(dr_mu_weights)
    cat("\n")

    cat("--- AUUC SCORES & PERMUTATION TESTS (B=1000) ---\n")
    cat(sprintf("BCF Horseshoe   : AUUC = %.1f, p-val = %.4f\n", auuc_nbcf, perm_bcf_shoe$p_val))
    cat(sprintf("Standard BCF    : AUUC = %.1f, p-val = %.4f\n", auuc_standard_bcf, perm_standard_bcf$p_val))
    cat(sprintf("Old Linear      : AUUC = %.1f, p-val = %.4f\n", auuc_old, perm_old$p_val))
    cat(sprintf("DR-Learner      : AUUC = %.1f, p-val = %.4f\n", auuc_dr, perm_dr$p_val))
    cat("\n")

    cat("--- SECHIDIS GLOBAL SCORE TEST ---\n")
    cat(sprintf("p-value: %.4f\n\n", score_p_val))

    cat("--- SIGNIFICANT INDIVIDUALS (Differing from ATE) ---\n")
    cat(sprintf("Total Individuals: %d\n\n", sig_nbcf$n_total))
    
    print_sig <- function(name, sig_list) {
      cat(sprintf("%s:\n", name))
      cat(sprintf("  Tolerance Interval: %d (%.1f%%)\n", sig_list$tolerance_sig, 100 * sig_list$tolerance_sig / sig_list$n_total))
      cat(sprintf("  Wan Interval      : %d (%.1f%%)\n", sig_list$wan_sig, 100 * sig_list$wan_sig / sig_list$n_total))
      cat(sprintf("  Sechidis Interval : %d (%.1f%%)\n\n", sig_list$sechidis_sig, 100 * sig_list$sechidis_sig / sig_list$n_total))
    }
    
    print_sig("BCF Horseshoe", sig_nbcf)
    print_sig("Standard BCF", sig_standard_bcf)
    print_sig("Old Linear Model", sig_old)
    
    sink()
    
    # Compile Results
    pb_inner$tick(1)
    
    heterogeneity_present <- (sum(final_S) > 0) && (sum(final_S) < length(final_S))
    final_heterogeneity_decision <- heterogeneity_present | score_test_significant
    
    results_list[[dataset_id]] <- data.frame(
      Dataset = dataset_id,
      EndpointType = ifelse(is_binary, "Binary", "Continuous"),
      Subgroup_Heterogeneity = heterogeneity_present,
      ScoreTest_PValue = score_p_val,
      ScoreTest_Significant = score_test_significant,
      Final_Heterogeneity_Decision = final_heterogeneity_decision,
      ATE_Overall = ate_res_nbcf,
      Subgroup_Proportion = subgroup_prop,
      Subgroup_ATE = subgroup_ate,
      Predictive_Variables = paste(predictive_vars, collapse = ", ")
    )
    
  }, error = function(e) {
    cat(sprintf("\nError processing %s: %s\n", file_name, e$message))
    traceback(x = e)
  })
  
  pb$tick()
}

# ==============================================================================
# 5. EXPORT RESULTS & ZIP
# ==============================================================================

if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  var_class_df <- bind_rows(var_class_list)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Dataset_Level_Results")
  writeData(wb, "Dataset_Level_Results", results_df)
  
  addWorksheet(wb, "Variable_Classifications")
  writeData(wb, "Variable_Classifications", var_class_df)
  
  saveWorkbook(wb, file.path(output_dir, "Submission_Results.xlsx"), overwrite = TRUE)
}

# Zip datasets
if (dir.exists(sub_datasets_dir)) {
  zip::zipr(zipfile = file.path(output_dir, "submission_datasets.zip"), files = list.files(sub_datasets_dir, full.names=TRUE))
  cat(sprintf("\nZipped datasets successfully to: %s\n", file.path(output_dir, "submission_datasets.zip")))
}

cat("\n================ BATCH ANALYSIS COMPLETE ================\n")