# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(speff2trial) # Contains the ACTG175 dataset
library(stochtree)   # Semi-Parametric BCF and BCF method
# DR-Learner Ensemble Libraries:
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(ranger)
library(coda)
source('R/shapley_aux.R', local = TRUE)

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

# --- 2. PREDICT FUNCTION FOR BCF (From Vignette) ---
predict.bcfmodel <- function(object, X, Z, propensity = NULL, rfx_group_ids = NULL, rfx_basis = NULL, ...){
  # Preprocess covariates
  print('1')
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
  print('flag1')
  result <- list("mu_hat" = mu_hat, "tau_hat" = tau_hat, "y_hat" = y_hat)
  print('flag2')
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

results_actg <- standardize_X_by_index_new(
  X_initial = X_actg_raw, 
  process_data = TRUE, 
  interaction_rule = "continuous_or_binary", 
  cat_coding_method = "sum"
)
X_train_mat <- results_actg$X_final

# ==============================================================================
# 4. FIT MODELS (Semi-Parametric vs BCF vs DR-Learner)
# ==============================================================================
cat("\nFitting Models... (Grab a coffee, this takes time)\n")

# ---------------------------------------------------------
# MODEL A: Semi-Parametric BART 
# ---------------------------------------------------------
cat("  Fitting Semi-Parametric BART...\n")
general_params_nbcf <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "none", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 123, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = FALSE, 
  sample_global_prior = "half-cauchy", unlink = TRUE, 
  propensity_seperate = "none", gibbs = TRUE, step_out = 0.5, 
  max_steps = 150, save_output = FALSE, probit_outcome_model = FALSE, 
  interaction_rule = "continuous_or_binary", standardize_cov = FALSE, 
  simple_prior = FALSE, save_partial_residual = FALSE, regularize_ATE = FALSE,
  sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1
)

fit_nbcf <- bcf_linear_probit(
  X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5,
  num_gfr = 50, num_burnin = 2000, num_mcmc = 4000,
  general_params = general_params_nbcf
)

# ---------------------------------------------------------
# MODEL B: BCF (stochtree)
# ---------------------------------------------------------
cat("  Fitting BCF...\n")
general_params_bcf <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001, 
  variable_weights = NULL, propensity_covariate = "none", 
  adaptive_coding = TRUE, control_coding_init = 0, 
  treated_coding_init = 1, rfx_prior_var = NULL, 
  random_seed = 123, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 2, verbose = FALSE, 
  probit_outcome_model = FALSE
)

fit_bcf <- bcf(
  X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg, 
  num_gfr = 25, num_burnin = 2000, num_mcmc = 4000, 
  general_params = general_params_bcf
)

# ---------------------------------------------------------
# MODEL C: DR-Learner via SuperLearner
# ---------------------------------------------------------
cat("  Fitting DR-Learner via SuperLearner...\n")
sl_library <- c("SL.glm", "SL.glmnet", "SL.gam", "SL.xgboost", "SL.ranger")

pi_hat <- mean(Z_actg) 
df_dr <- data.frame(Y = Y_actg, A = Z_actg, X_train_mat)

set.seed(123)
K_folds <- 10
fold_ids <- sample(rep(1:K_folds, length.out = nrow(df_dr)))
pseudo_Y <- numeric(nrow(df_dr))

for (k in 1:K_folds) {
  cat(sprintf("    Processing Cross-Fit Fold %d of %d...\n", k, K_folds))
  train_idx <- which(fold_ids != k)
  val_idx <- which(fold_ids == k)
  
  X_train_k <- df_dr[train_idx, setdiff(names(df_dr), "Y")]
  Y_train_k <- df_dr$Y[train_idx]
  
  sl_mu <- SuperLearner(
    Y = Y_train_k, X = X_train_k, 
    SL.library = sl_library, family = gaussian(),
    cvControl = list(V = 10)
  )
  
  X_val_1 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_1$A <- 1
  X_val_0 <- df_dr[val_idx, setdiff(names(df_dr), "Y")]; X_val_0$A <- 0
  
  mu_1_hat <- predict(sl_mu, newdata = X_val_1)$pred
  mu_0_hat <- predict(sl_mu, newdata = X_val_0)$pred
  mu_A_hat <- ifelse(df_dr$A[val_idx] == 1, mu_1_hat, mu_0_hat)
  
  pseudo_Y[val_idx] <- (df_dr$A[val_idx] - pi_hat) / (pi_hat * (1 - pi_hat)) * (df_dr$Y[val_idx] - mu_A_hat) + (mu_1_hat - mu_0_hat)
}

cat("    Fitting final CATE model on pseudo-outcomes...\n")
X_train_df <- as.data.frame(X_train_mat)
fit_dr <- SuperLearner(
  Y = pseudo_Y, X = X_train_df, 
  SL.library = sl_library, family = gaussian(),
  cvControl = list(V = 10)
)
cate_hat_dr <- as.vector(fit_dr$SL.predict)
cat("\n--- DR-LEARNER SUPERLEARNER ENSEMBLE WEIGHTS ---\n")
print(fit_dr$coef)
# ==============================================================================
# 5. INFERENCE (CATEs and ATE)
# ==============================================================================

# Extract CATEs for Semi-Parametric BART
tau_draws_nbcf <- predict_linear_bcf_patched(fit_nbcf, X = X_train_mat, Z = Z_actg)$tau_hat
cate_hat_nbcf <- rowMeans(tau_draws_nbcf)
ate_draws_nbcf <- colMeans(tau_draws_nbcf)
ate_res_nbcf <- c(Mean = mean(ate_draws_nbcf), quantile(ate_draws_nbcf, probs = c(0.025, 0.975)))

# Extract CATEs for BCF explicitly passing Z
tau_draws_bcf <- predict.bcfmodel(fit_bcf, X = X_train_mat, Z = Z_actg)$tau_hat
cate_hat_bcf <- rowMeans(tau_draws_bcf)
ate_draws_bcf <- colMeans(tau_draws_bcf)
ate_res_bcf <- c(Mean = mean(ate_draws_bcf), quantile(ate_draws_bcf, probs = c(0.025, 0.975)))

# Extract ATE for DR-Learner
ate_res_dr <- mean(pseudo_Y) 

cat("\n--- IN-SAMPLE INFERENCE RESULTS ---\n")
cat(sprintf("Semi-Parametric BART ATE: %.2f [95%% CI: %.2f, %.2f]\n", ate_res_nbcf[1], ate_res_nbcf[2], ate_res_nbcf[3]))
cat(sprintf("BCF ATE: %.2f [95%% CI: %.2f, %.2f]\n", ate_res_bcf[1], ate_res_bcf[2], ate_res_bcf[3]))
cat(sprintf("DR-Learner ATE (AIPW): %.2f\n", ate_res_dr))

# ==============================================================================
# 6. GGPLOT VISUALIZATIONS (3-Way Compare)
# ==============================================================================
# A. Density Distributions of Estimated CATEs
density_data <- data.frame(
  CATE = c(cate_hat_nbcf, cate_hat_bcf, cate_hat_dr),
  Model = rep(c("Semi-Parametric BART", "BCF", "DR-Learner"), each = length(cate_hat_nbcf))
)

colors_3way <- c("Semi-Parametric BART" = "#008080", "BCF" = "#d95f02", "DR-Learner" = "purple")

p_dist <- ggplot(density_data, aes(x = CATE, fill = Model, color = Model)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
  scale_fill_manual(values = colors_3way) +
  scale_color_manual(values = colors_3way) +
  labs(title = "Distribution of Estimated CATEs",
       x = "Treatment Effect (Difference in CD4 at 20 Weeks)", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# B. Treatment Effect Heterogeneity by Baseline CD4
plot_data_het <- data.frame(
  CD4_Baseline = actg_sub$cd40,
  CATE_NBCF = cate_hat_nbcf,
  CATE_BCF = cate_hat_bcf,
  CATE_DR = cate_hat_dr
)

p_het <- ggplot(plot_data_het) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), method = "loess", se = FALSE, linewidth = 1.2) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_BCF, color = "BCF"), method = "loess", se = FALSE, linewidth = 1.2) +
  geom_smooth(aes(x = CD4_Baseline, y = CATE_DR, color = "DR-Learner"), method = "loess", se = FALSE, linewidth = 1.2) +
  geom_point(aes(x = CD4_Baseline, y = CATE_NBCF, color = "Semi-Parametric BART"), alpha = 0.15) +
  geom_point(aes(x = CD4_Baseline, y = CATE_BCF, color = "BCF"), alpha = 0.15) +
  geom_point(aes(x = CD4_Baseline, y = CATE_DR, color = "DR-Learner"), alpha = 0.15) +
  scale_color_manual(values = colors_3way) +
  labs(title = "Treatment Effect Heterogeneity",
       subtitle = "Estimated CATE vs Baseline CD4 T-cell Count",
       x = "Baseline CD4 Count", y = "Estimated CATE", color = "Model") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p_dist)
print(p_het)

# C. CATE Correlation Heatmap
cate_df <- data.frame(
  "Semi-Parametric BART" = cate_hat_nbcf,
  "BCF" = cate_hat_bcf,
  "DR-Learner" = cate_hat_dr,
  check.names = FALSE
)

cor_matrix <- cor(cate_df)
cor_data <- as.data.frame(as.table(cor_matrix))
names(cor_data) <- c("Model1", "Model2", "Correlation")

p_cor_heat <- ggplot(cor_data, aes(x = Model1, y = Model2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), 
            color = ifelse(cor_data$Correlation > 0.8, "white", "black"), 
            size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#cc0000", mid = "white", high = "#008080", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  labs(title = "CATE Correlation Heatmap", subtitle = "Pairwise agreement between model treatment effect estimates") +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")

print(p_cor_heat)

# ==============================================================================
# 7. IMPORTANCE ANALYSIS
# ==============================================================================
cat("\n--- STARTING IMPORTANCE ANALYSIS ---\n")

# A. Compute Shapley for Semi-Parametric BART
X_actg <- results_actg$X_final
if (!is.null(fit_nbcf$interaction_pairs)) {
  ipairs <- fit_nbcf$interaction_pairs
} else {
  var_info <- results_actg$X_final_var_info
  is_candidate <- var_info$is_continuous | var_info$is_binary
  ipairs <- interaction_pairs(ncol(X_actg), is_candidate)
}

shapley_results <- compute_shapley_all(X = X_actg, beta_post = fit_nbcf$Beta, beta_int_post = fit_nbcf$Beta_int, ipairs = ipairs)
p_shap_importance <- plot_shapley_importance_breakdown(shapley_results) + labs(title = "Semi-Parametric BART: Shapley Importance")
print(p_shap_importance)

# B. Compute Shapley for BCF using shapr (assuming independence)
cat("\nCalculating BCF Shapley values using 'shapr' package...\n")
library(shapr)

# Prediction wrapper for shapr
predict_model.bcf <- function(x, newdata) {
  # newdata is expected to be a data.frame matching X_train_mat
  # Since BCF CATE (tau_hat) shouldn't strictly depend on Z for its own internal structure,
  # but predict.bcfmodel requires Z, we provide a dummy Z.
  newdata_df <- as.data.frame(newdata)
  dummy_z <- rep(1, nrow(newdata_df))
  res <- predict.bcfmodel(x, X = newdata_df, Z = dummy_z)
  rowMeans(res$tau_hat)
}

# Compute empirical Shapley values, using the mean CATE as baseline
explanation_bcf <- explain(
  model = fit_bcf,
  x_train = as.data.frame(X_train_mat),
  x_explain = as.data.frame(X_train_mat),
  approach = "empirical", 
  phi0 = mean(cate_hat_bcf),
  predict_model = predict_model.bcf
)

cat("\n--- BCF SHAPLEY VALUES (MEAN ABSOLUTE) ---\n")
# Extract Shapley values (excluding the "none" baseline column)
bcf_shap_vals <- as.matrix(explanation_bcf$dt[, -1])
shap_abs_mean <- colMeans(abs(bcf_shap_vals))
print(sort(shap_abs_mean, decreasing = TRUE))

# # B. Compute Vansteelandt Leave-One-Out (LOO) TE-VIMs for the DR-Learner
# cat("\nCalculating Vansteelandt LOO TE-VIMs using SuperLearner...\n")
# mse_full <- mean((pseudo_Y - cate_hat_dr)^2)
# te_vims <- numeric(ncol(X_train_mat))
# names(te_vims) <- colnames(X_train_mat)
# 
# for(j in 1:ncol(X_train_mat)) {
#   col_to_drop <- names(te_vims)[j]
#   cat(sprintf("  Evaluating removal of %s (%d of %d)...\n", col_to_drop, j, ncol(X_train_mat)))
#   X_loo <- X_train_df[, -which(colnames(X_train_df) == col_to_drop), drop = FALSE]
#   sl_loo <- SuperLearner(Y = pseudo_Y, X = X_loo, SL.library = sl_library, family = gaussian(), cvControl = list(V = 10))
#   mse_loo <- mean((pseudo_Y - sl_loo$SL.predict)^2)
#   te_vims[j] <- mse_loo - mse_full
# }
# 
# te_vim_df <- data.frame(Feature = names(te_vims), TE_VIM = te_vims) %>% arrange(desc(TE_VIM))
# 
# p_tevims <- ggplot(te_vim_df, aes(x = TE_VIM, y = reorder(Feature, TE_VIM))) +
#   geom_col(fill = "purple", alpha = 0.8) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   labs(title = "DR-Learner: LOO TE-VIMs", subtitle = "Increase in Mean Squared Error (MSE) when removed", x = "TE-VIM (MSE Difference)", y = "Feature") +
#   theme_minimal(base_size = 14)
# print(p_tevims)

# ==============================================================================
# 8 & 9 & 10. ROBUST UPLIFT CURVE EVALUATION (AUUC)
# ==============================================================================
cat("\n--- GENERATING ROBUST UPLIFT CURVES (AVOIDING SAMPLING BIAS) ---\n")

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

eval_df <- data.frame(Y = Y_actg, Z = Z_actg, CATE_NBCF = cate_hat_nbcf, CATE_BCF = cate_hat_bcf, CATE_DR = cate_hat_dr)

uplift_nbcf <- get_eval_curves(eval_df, "CATE_NBCF", "Semi-Parametric BART")
uplift_bcf <- get_eval_curves(eval_df, "CATE_BCF", "BCF")
uplift_dr <- get_eval_curves(eval_df, "CATE_DR", "DR-Learner")
uplift_plot_data <- bind_rows(uplift_nbcf, uplift_bcf, uplift_dr)

calc_auc <- function(df, metric) { sum(diff(df$frac) * (head(df[[metric]], -1) + tail(df[[metric]], -1)) / 2) }

auuc_nbcf <- calc_auc(uplift_nbcf, "uplift")
auuc_bcf <- calc_auc(uplift_bcf, "uplift")
auuc_dr <- calc_auc(uplift_dr, "uplift")

auuc_label <- sprintf("AUUC:\nSemi-Parametric BART = %.1f\nBCF = %.1f\nDR-Learner = %.1f", auuc_nbcf, auuc_bcf, auuc_dr)

final_uplift <- tail(uplift_nbcf$uplift, 1)
p_uplift_smooth <- ggplot(uplift_plot_data, aes(x = frac, y = uplift, color = Model)) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = final_uplift, 
           color = "black", linetype = "dashed", linewidth = 0.8) +
  
  # --- UPDATED: Moved to Bottom-Right ---
  # x = 0.95 places it near the right edge
  # y = min(uplift) places it at the very bottom of the data bounds
  # hjust = 1 (right-aligned), vjust = 0 (bottom-aligned) prevents cutoff
  annotate("label", x = 0.95, y = min(uplift_plot_data$uplift, na.rm = TRUE), 
           label = auuc_label, hjust = 1, vjust = 0, 
           fontface = "bold", size = 4.5, color = "black", 
           fill = "white", alpha = 0.85, label.size = NA) +
  
  labs(title = "Robust Uplift",
       x = "Fraction of Population Treated", y = "Cumulative True Uplift (Adjusted CD4 Gain)") +
  scale_color_manual(values = colors_3way) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p_uplift_smooth)
# ==============================================================================
# 11 & 12. UNCERTAINTY QUANTIFICATION: TOLERANCE BANDS
# ==============================================================================
cat("\n--- QUANTIFYING SIGNIFICANT HETEROGENEITY (TOLERANCE BANDS) ---\n")

postsum <- function(posters, conf = 0.05) {
  qs <- quantile(posters, p = c(conf/2, 1 - conf/2))
  mn <- mean(posters)
  return(c(qs[1], mn, qs[2]))
}

intfrompost1 <- function(post, conf = 0.05){ t(apply(post, 1, postsum, conf = conf)) }

conftol <- function(margconf = 0.05, post, tol = 0.05){
  int <- intfrompost1(post, conf = margconf)
  difl <- ((post - int[,1]) >= 0) * ((post - int[,3]) <= 0)
  cs <- colMeans(difl)
  return(mean(cs >= 1 - tol))
} 

findconfglob <- function(post, tol = 0.05){
  grid <- seq(0.0005, 0.05, by = 0.0005)
  conftols <- sapply(grid, conftol, post = post, tol = tol)
  wm <- which(conftols - 0.95 < 0)[1] - 1
  if(is.na(wm) || wm < 1) return(0.0005) 
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
    scale_fill_manual(values = c("Tolerance Interval" = "#F0E442", "95% Credible Interval" = "#0072B2")) +
    scale_color_manual(values = c("Significant Positive (Benefit > Baseline)" = "#005a5a", "Significant Negative (Benefit < Baseline)" = "#cc0000")) +
    labs(title = paste0("Individual Treatment Effect Heterogeneity: ", model_name), subtitle = paste0(round(pct_sig, 1), "% of patients show significant deviation from baseline."), x = "Patients (Ordered by Effect Size)", y = "Heterogeneous Effect (\u0394 CATE)") +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
}

# A. Semi-Parametric BART (Heterogeneity = CATE - alpha)
alpha_draws_nbcf <- as.vector(t(fit_nbcf$alpha))
het_draws_nbcf <- t(t(tau_draws_nbcf) - alpha_draws_nbcf)
p_het_nbcf <- build_heterogeneity_plot(het_draws_nbcf, "Semi-Parametric BART (CATE - \u03b1)")

# B. BCF (Heterogeneity = CATE - ATE)
het_draws_bcf <- t(t(tau_draws_bcf) - ate_draws_bcf)
p_het_bcf <- build_heterogeneity_plot(het_draws_bcf, "BCF (CATE - ATE)")

print(p_het_nbcf)
print(p_het_bcf)

# ==============================================================================
# 13. SAVE PLOTS TO LOCAL DIRECTORY
# ==============================================================================
plot_dir_default <- "C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/PhD project/plots/ACTG"
if (!exists("plot_dir")) {
  plot_dir <- plot_dir_default
}
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

plot_w <- 8; plot_h <- 6; plot_dpi <- 300

ggsave(file.path(plot_dir, "01_CATE_Density_Distribution.png"), p_dist, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "02a_Treatment_Heterogeneity.png"), p_het, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "02b_CATE_Correlation_Heatmap.png"), p_cor_heat, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "03_SemiParametric_Shapley.png"), p_shap_importance, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "05_Robust_Uplift_Curve_Comparison.png"), p_uplift_smooth, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "06a_Heterogeneity_Bands_SemiParametric.png"), p_het_nbcf, width = 10, height = 7, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "06b_Heterogeneity_Bands_Standard_BCF.png"), p_het_bcf, width = 10, height = 7, dpi = plot_dpi, bg = "white")

cat(sprintf("All plots successfully saved to: %s\n", plot_dir)) 

# ==============================================================================
# 14. SHAPLEY DISSECTION: MAIN VS INTERACTION EFFECTS (FEATURE-LEVEL)
# ==============================================================================
cat("\n--- RUNNING SHAPLEY DISSECTION ANALYSIS FOR 'cd40' ---\n")

# A. Prepare Data & Identify Target Feature
X_actg <- results_actg$X_final
feature_to_dissect <- "cd40" # Focus on Baseline CD4 count

# Determine interaction pairs
if (!is.null(fit_nbcf$interaction_pairs)) {
  ipairs <- fit_nbcf$interaction_pairs
} else {
  var_info <- results_actg$X_final_var_info
  is_candidate <- var_info$is_continuous | var_info$is_binary
  ipairs <- interaction_pairs(ncol(X_actg), is_candidate)
}

# Extract dissected posterior
dissected_post <- get_shapley_posterior_dissected(
  feature_name = feature_to_dissect, 
  indices = 1:nrow(X_actg), 
  X = X_actg, 
  beta_post = fit_nbcf$Beta, 
  beta_int_post = fit_nbcf$Beta_int, 
  ipairs = ipairs
)

# B. Calculate Medians
plot_df_dissect <- data.frame(
  Feature_Value = X_actg[, feature_to_dissect],
  Main_Med = apply(dissected_post$main, 1, median),
  Int_Med = apply(dissected_post$interaction, 1, median)
) %>%
  mutate(
    Total_Shapley = Main_Med + Int_Med,
    # Calculate Interaction Ratio: |Int| / (|Main| + |Int|)
    Int_Ratio = abs(Int_Med) / (abs(Main_Med) + abs(Int_Med) + 1e-8) # 1e-8 to prevent div by zero
  )

# ------------------------------------------------------------------------------
# PLOT 1: Synergy vs. Antagonism Scatter Plot
# ------------------------------------------------------------------------------
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

print(p_synergy)

# ------------------------------------------------------------------------------
# PLOT 2: Fractional Contribution Profile
# ------------------------------------------------------------------------------
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

print(p_fractional)

# ------------------------------------------------------------------------------
# PLOT 3: Dual-Density Ridge Plots for Prototypical Individuals
# ------------------------------------------------------------------------------
cat("  Generating Dual-Density Uncertainty Plots...\n")
# Find individuals at the 10th, 50th, and 90th percentiles of cd40
quant_vals <- quantile(plot_df_dissect$Feature_Value, probs = c(0.1, 0.5, 0.9))
idx_10 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[1]))
idx_50 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[2]))
idx_90 <- which.min(abs(plot_df_dissect$Feature_Value - quant_vals[3]))

# Extract full posteriors for these 3 patients
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
    subtitle = "Comparing density spreads across 5,000 MCMC draws for prototypical patients",
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

print(p_uncertainty)

# ------------------------------------------------------------------------------
# 15. SAVE DISSECTION PLOTS
# ------------------------------------------------------------------------------
cat("  Saving dissection plots to local directory...\n")

ggsave(file.path(plot_dir, "07a_Shapley_Synergy_Antagonism.png"), p_synergy, width = 8, height = 8, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "07b_Shapley_Fractional_Contribution.png"), p_fractional, width = plot_w, height = plot_h, dpi = plot_dpi, bg = "white")
ggsave(file.path(plot_dir, "07c_Shapley_DualDensity_Uncertainty.png"), p_uncertainty, width = 8, height = 10, dpi = plot_dpi, bg = "white")

cat("\n==============================================================================\n")
cat("ALL ANALYSES COMPLETE.\n")
cat("==============================================================================\n")