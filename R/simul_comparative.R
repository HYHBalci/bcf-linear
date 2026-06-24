# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(stochtree)   # Semi-Parametric BCF method
source('R/shapley_aux.R', local = TRUE)

# ==============================================================================
# 2. HELPER FUNCTIONS: PREPROCESSING & PATCHED PREDICTIONS
# ==============================================================================
# From ACTG_new.R
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

# ==============================================================================
# 3. DATA SIMULATION (From simul_1.R)
# ==============================================================================
generate_data_medical <- function(n = 250,
                                  is_te_hetero = TRUE,
                                  is_mu_nonlinear = TRUE, 
                                  seed = 1848, 
                                  RCT = FALSE, 
                                  scenario = "default", 
                                  z_diff = FALSE, 
                                  contrast_binary = TRUE, 
                                  BCF = FALSE, 
                                  sigma_sq = 1) {
  set.seed(seed)
  
  x1 <- rnorm(n, mean=0, sd=1)
  x2 <- rnorm(n, mean=0, sd=1)
  x3 <- rnorm(n, mean=0, sd=1)
  x4 <- rbinom(n, size=1, prob=0.5)
  if(contrast_binary){
    x4 <- 2 * x4 - 1
  }
  x5_raw <- sample(1:3, size=n, replace=TRUE, prob=c(1/3,1/3,1/3))
  
  g_func <- function(x5) {
    out <- rep(NA, length(x5))
    out[x5 == 1] <-  2
    out[x5 == 2] <- -1
    out[x5 == 3] <- -4
    return(out)
  }
  
  g_x5 <- g_func(x5_raw)
  
  if (!is_mu_nonlinear) {
    mu <- 1 + g_x5 + x1*x3
  } else {
    mu <- -6 + g_x5 + 6*abs(x3 - 1)
  }
  
  if (!is_te_hetero) {
    tau_vec <- rep(3, n)
  } else {
    if(scenario == "complex_interaction"){ 
      tau_vec <- 1 + 4*x1 + 3*x2 + 2*x2*x1
    } else if(scenario == "medical_saturation"){
      tau_vec <- 0.5 + 4 / (1 + exp(-2.5 * x1))
    } else {
      tau_vec <- 1 + 2*x2*x4
    }
  }
  
  s <- sd(mu)
  u_i <- runif(n, 0, 1)
  Phi <- function(z) pnorm(z, mean=0, sd=1)
  
  if (RCT) {
    pi_x <- rep(0.5, n)
  } else {
    pi_x <- 0.8 * Phi((3*mu)/s - 0.5*x1) + 0.05 + (u_i / 10)
  }
  
  pi_x <- pmin(pmax(pi_x, 0), 1)
  z <- rbinom(n, size=1, prob=pi_x)
  eps <- rnorm(n, 0, sqrt(sigma_sq))
  
  if(BCF){
    z_binary <- z
  }
  
  if(is.logical(z_diff)){
    delta <- if(z_diff) 0.5 else 0
  } else {
    delta <- z_diff
  }
  
  if(z_diff != FALSE){
    z <- z - delta
  }
  
  y <- mu + z*tau_vec + eps
  y_hat <- mu + z*tau_vec
  
  x5_factor <- factor(x5_raw, levels = c(1, 2, 3))
  contrasts(x5_factor) <- contr.sum(3)
  x5_dev <- model.matrix(~ x5_factor, data = data.frame(x5_factor))
  x5_dev <- x5_dev[, -1] 
  colnames(x5_dev) <- c("x5_1", "x5_2")
  
  if(BCF){
    z <- z_binary
  }
  
  df <- data.frame(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    x5_1 = x5_dev[, 1], x5_2 = x5_dev[, 2], 
    z = z, y = y, mu = mu, pi_x = pi_x, 
    tau = tau_vec, y_hat = y_hat
  )
  
  return(df)
}

cat("\n--- GENERATING SIMULATED DATA ---\n")
n_sim <- 500
df_sim <- generate_data_medical(n = n_sim, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, RCT = FALSE, z_diff = 0.5)

X_sim_mat <- df_sim %>% dplyr::select(x1, x2, x3, x4, x5_1, x5_2) %>% as.matrix()
Y_sim <- df_sim$y
Z_sim <- df_sim$z
tau_true <- df_sim$tau

# ==============================================================================
# 4. FIT MODELS (Standard, None, OLS, T-distribution)
# ==============================================================================
cat("\nFitting Models on Simulated Data...\n")

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
fit_std <- bcf_linear_probit(
  X_train = X_sim_mat, y_train = Y_sim, Z_train = Z_sim,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = params_std
)

# 2. None
cat("  Fitting No Shrinkage (none)...\n")
params_none <- general_params_base
params_none$sample_global_prior <- "none"
fit_none <- bcf_linear_probit(
  X_train = X_sim_mat, y_train = Y_sim, Z_train = Z_sim,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = params_none
)

# 3. OLS
cat("  Fitting OLS Shrinkage...\n")
params_ols <- general_params_base
params_ols$sample_global_prior <- "OLS"
fit_ols <- bcf_linear_probit(
  X_train = X_sim_mat, y_train = Y_sim, Z_train = Z_sim,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = params_ols
)

# 4. Robust t-distribution
cat("  Fitting Robust t-distribution...\n")
params_tdist <- general_params_base
params_tdist$sample_global_prior <- "half-cauchy"
params_tdist$robust <- TRUE
params_tdist$robust_nu <- 3
fit_tdist <- bcf_linear_probit(
  X_train = X_sim_mat, y_train = Y_sim, Z_train = Z_sim,
  num_gfr = 50, num_burnin = 1000, num_mcmc = 3000, general_params = params_tdist
)

# ==============================================================================
# 5. COMPARE METHODS (RMSE for CATEs)
# ==============================================================================
cat("\n--- COMPARISON OF METHODS (CATE RMSE) ---\n")

cate_std <- rowMeans(predict_linear_bcf_patched(fit_std, X = X_sim_mat, Z = Z_sim)$tau_hat)
cate_none <- rowMeans(predict_linear_bcf_patched(fit_none, X = X_sim_mat, Z = Z_sim)$tau_hat)
cate_ols <- rowMeans(predict_linear_bcf_patched(fit_ols, X = X_sim_mat, Z = Z_sim)$tau_hat)
cate_tdist <- rowMeans(predict_linear_bcf_patched(fit_tdist, X = X_sim_mat, Z = Z_sim)$tau_hat)

rmse <- function(true, est) sqrt(mean((true - est)^2))

results_rmse <- data.frame(
  Method = c("Standard (Half-Cauchy)", "None", "OLS", "Robust t-distribution"),
  RMSE = c(rmse(tau_true, cate_std), rmse(tau_true, cate_none), rmse(tau_true, cate_ols), rmse(tau_true, cate_tdist))
)
print(results_rmse)
write.csv(results_rmse, "simul_comparative_rmse.csv", row.names = FALSE)

plot_data <- data.frame(
  True_CATE = rep(tau_true, 4),
  Estimated_CATE = c(cate_std, cate_none, cate_ols, cate_tdist),
  Method = rep(c("Standard", "None", "OLS", "t-dist"), each = n_sim)
)

p_comp <- ggplot(plot_data, aes(x = True_CATE, y = Estimated_CATE, color = Method)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(title = "True vs Estimated CATE (Simulated Data)", x = "True CATE", y = "Estimated CATE")

print(p_comp)
ggsave("simul_comparative_plot.png", p_comp, width = 8, height = 6)

cat("\nSaving model fits to RDS files...\n")
saveRDS(fit_std, "simul_fit_std.rds")
saveRDS(fit_none, "simul_fit_none.rds")
saveRDS(fit_ols, "simul_fit_ols.rds")
saveRDS(fit_tdist, "simul_fit_tdist.rds")

cat("\nResults saved to 'simul_comparative_rmse.csv' and 'simul_comparative_plot.png'\n")
cat("Model fits saved to 'simul_fit_*.rds' files.\n")
cat("\nDone!\n")
