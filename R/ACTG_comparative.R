
# ==========================================
# AUTO-INJECTED ETA TIMER
# ==========================================
total_models_to_run <- 5
models_run <- 0

bcf_linear_probit_eta <- function(...) {
  t0 <- Sys.time()
  res <- bcf_linear_probit(...)
  t1 <- Sys.time()
  dur <- as.numeric(difftime(t1, t0, units="mins"))
  models_run <<- models_run + 1
  models_left <- total_models_to_run - models_run
  eta_mins <- dur * models_left
  cat(sprintf("\n[TIMER] Model finished in %.2f mins. %d models left. ETA: %.2f mins\n\n", dur, models_left, eta_mins))
  return(res)
}
# ==========================================

# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(speff2trial) # Contains the ACTG175 dataset
library(stochtree)   # Semi-Parametric BCF method
source('R/shapley_aux.R', local = TRUE)

dir.create("plots", showWarnings = FALSE)

# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PATCHED PREDICTIONS
# ==============================================================================
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

# Uplift functions
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
# 3. PREPARE ACTG175 DATA
# ==============================================================================
cat("\n--- PREPARING ACTG175 DATA ---\n")
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
# 4. FIT MODELS
# ==============================================================================
cat("\nFitting Models on ACTG175 Data...\n")

# Workaround for the scoping bug in stochtree's bcf_linear_probit.R for robust models
probit_outcome_model <<- FALSE

general_params_base <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = TRUE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "none", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 123, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = FALSE, 
  unlink = TRUE, propensity_seperate = "none", gibbs = TRUE, 
  step_out = 0.5, max_steps = 150, save_output = FALSE, 
  probit_outcome_model = FALSE, interaction_rule = "continuous_or_binary", 
  standardize_cov = FALSE, simple_prior = FALSE, 
  save_partial_residual = FALSE, regularize_ATE = FALSE,
  sigma_residual = 0, hn_scale = 0, use_ncp = FALSE, n_tijn = 1
)

# 1. Standard (Half-Cauchy)
cat("  Fitting Standard (Half-Cauchy)...\n")
params_std <- general_params_base
params_std$sample_global_prior <- "half-cauchy"
fit_std <- bcf_linear_probit_eta(X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5, num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = params_std)

# 2. None
cat("  Fitting No Shrinkage (none)...\n")
params_none <- general_params_base
params_none$sample_global_prior <- "none"
fit_none <- bcf_linear_probit_eta(X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5, num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = params_none)

# 3. OLS
cat("  Fitting OLS Shrinkage...\n")
params_ols <- general_params_base
params_ols$sample_global_prior <- "OLS"
fit_ols <- bcf_linear_probit_eta(X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5, num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = params_ols)

# 4. Robust t-distribution
cat("  Fitting Robust t-distribution...\n")
params_tdist <- general_params_base
params_tdist$sample_global_prior <- "half-cauchy"
params_tdist$robust <- TRUE
params_tdist$robust_nu <- 3
fit_tdist <- bcf_linear_probit_eta(X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5, num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = params_tdist)

# 5. Hybrid
cat("  Fitting Hybrid Shrinkage...\n")
params_hybrid <- general_params_base
params_hybrid$sample_global_prior <- "hybrid"
fit_hybrid <- bcf_linear_probit_eta(X_train = X_train_mat, y_train = Y_actg, Z_train = Z_actg - 0.5, num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = params_hybrid)


# ==============================================================================
# 5. INFERENCE (CATEs and ATE)
# ==============================================================================
tau_std <- predict_linear_bcf_patched(fit_std, X = X_train_mat, Z = Z_actg)$tau_hat
cate_std <- rowMeans(tau_std)
ate_std <- mean(colMeans(tau_std))

tau_none <- predict_linear_bcf_patched(fit_none, X = X_train_mat, Z = Z_actg)$tau_hat
cate_none <- rowMeans(tau_none)
ate_none <- mean(colMeans(tau_none))

tau_ols <- predict_linear_bcf_patched(fit_ols, X = X_train_mat, Z = Z_actg)$tau_hat
cate_ols <- rowMeans(tau_ols)
ate_ols <- mean(colMeans(tau_ols))

tau_tdist <- predict_linear_bcf_patched(fit_tdist, X = X_train_mat, Z = Z_actg)$tau_hat
cate_tdist <- rowMeans(tau_tdist)
ate_tdist <- mean(colMeans(tau_tdist))

tau_hybrid <- predict_linear_bcf_patched(fit_hybrid, X = X_train_mat, Z = Z_actg)$tau_hat
cate_hybrid <- rowMeans(tau_hybrid)
ate_hybrid <- mean(colMeans(tau_hybrid))

cat("\n--- IN-SAMPLE INFERENCE RESULTS ---\n")
cat(sprintf("Standard (Half-Cauchy) ATE: %.2f\n", ate_std))
cat(sprintf("No Shrinkage ATE: %.2f\n", ate_none))
cat(sprintf("OLS ATE: %.2f\n", ate_ols))
cat(sprintf("Robust t-distribution ATE: %.2f\n", ate_tdist))
cat(sprintf("Hybrid ATE: %.2f\n", ate_hybrid))

ate_results <- data.frame(
  Method = c("Standard", "None", "OLS", "Robust t-dist", "Hybrid"),
  ATE = c(ate_std, ate_none, ate_ols, ate_tdist, ate_hybrid)
)
write.csv(ate_results, "plots/ACTG_comparative_ate.csv", row.names = FALSE)

# ==============================================================================
# 6. VISUALIZATION & UPLIFT
# ==============================================================================
colors_methods <- c("Standard" = "#1b9e77", "None" = "#d95f02", "OLS" = "#7570b3", "Robust t-dist" = "#e7298a", "Hybrid" = "#66a61e")

# Density Plot
density_data <- data.frame(
  CATE = c(cate_std, cate_none, cate_ols, cate_tdist, cate_hybrid),
  Model = rep(c("Standard", "None", "OLS", "Robust t-dist", "Hybrid"), each = length(cate_std))
)
p_dist <- ggplot(density_data, aes(x = CATE, fill = Model, color = Model)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
  scale_fill_manual(values = colors_methods) +
  scale_color_manual(values = colors_methods) +
  labs(title = "Distribution of Estimated CATEs", x = "Treatment Effect (Difference in CD4 at 20 Weeks)", y = "Density") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("plots/ACTG_comparative_density.png", p_dist, width = 8, height = 6)

# Uplift
eval_df <- data.frame(Y = Y_actg, Z = Z_actg, 
                      CATE_STD = cate_std, CATE_NONE = cate_none, CATE_OLS = cate_ols,
                      CATE_TDIST = cate_tdist, CATE_HYBRID = cate_hybrid)

uplift_std <- get_eval_curves(eval_df, "CATE_STD", "Standard")
uplift_none <- get_eval_curves(eval_df, "CATE_NONE", "None")
uplift_ols <- get_eval_curves(eval_df, "CATE_OLS", "OLS")
uplift_tdist <- get_eval_curves(eval_df, "CATE_TDIST", "Robust t-dist")
uplift_hybrid <- get_eval_curves(eval_df, "CATE_HYBRID", "Hybrid")

auuc_vals <- c(
  Standard = calc_auc(uplift_std, "uplift"),
  None = calc_auc(uplift_none, "uplift"),
  OLS = calc_auc(uplift_ols, "uplift"),
  Robust_tdist = calc_auc(uplift_tdist, "uplift"),
  Hybrid = calc_auc(uplift_hybrid, "uplift")
)

uplift_plot_data <- bind_rows(uplift_std, uplift_none, uplift_ols, uplift_tdist, uplift_hybrid)
auuc_label <- sprintf("AUUC:\nStandard = %.1f\nNone = %.1f\nOLS = %.1f\nRobust t-dist = %.1f\nHybrid = %.1f", 
                      auuc_vals["Standard"], auuc_vals["None"], auuc_vals["OLS"], auuc_vals["Robust_tdist"], auuc_vals["Hybrid"])

p_uplift <- ggplot(uplift_plot_data, aes(x = frac, y = uplift, color = Model)) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = tail(uplift_std$uplift, 1), color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 0.95, y = min(uplift_plot_data$uplift, na.rm = TRUE), label = auuc_label, hjust = 1, vjust = 0, fontface = "bold", size = 4, alpha = 0.85) +
  labs(title = "Robust Uplift Curves (ACTG175)", x = "Fraction of Population Treated", y = "Cumulative True Uplift (Adjusted CD4 Gain)") +
  scale_color_manual(values = colors_methods) +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("plots/ACTG_comparative_uplift.png", p_uplift, width = 8, height = 6)

# Correlation Heatmap
cate_df_all <- data.frame(
  "Standard" = cate_std,
  "None" = cate_none,
  "OLS" = cate_ols,
  "Robust t-dist" = cate_tdist,
  "Hybrid" = cate_hybrid,
  check.names = FALSE
)

cor_matrix <- cor(cate_df_all)
cor_data <- as.data.frame(as.table(cor_matrix))
names(cor_data) <- c("Model1", "Model2", "Correlation")

p_cor_heat <- ggplot(cor_data, aes(x = Model1, y = Model2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), color = ifelse(cor_data$Correlation > 0.8, "white", "black"), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#cc0000", mid = "white", high = "#008080", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Pearson\nCorrelation") +
  labs(title = "CATE Correlation Heatmap", subtitle = "Pairwise agreement between model priors on ACTG175") +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.text.x = element_text(angle = 15, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")

ggsave("plots/ACTG_comparative_correlation.png", p_cor_heat, width = 8, height = 6)

cat("\nResults saved to 'plots/' directory.\n")
cat("\nDone!\n")
