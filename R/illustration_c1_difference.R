# illustration_c1_difference.R
# This script illustrates the difference between methods for computing 
# simultaneous bounds/intervals for subgroup selection using Semi-Parametric BART and DR-Learners.

library(dplyr)
library(tidyr)
library(ggplot2)
library(stochtree)
library(MASS)
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(ranger)

# Load helper functions
source('R/simul_1.R')
source('R/forGiorgio.R')

# ==============================================================================
# HELPER FUNCTIONS (REFACTORING)
# ==============================================================================

#' Patched Semi-Parametric BART Predict Function (from ACTG_new.R)
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

#' Construct pseudo-outcomes using 5-fold cross-fitted SuperLearner
#' @param X Covariate matrix
#' @param Y Outcome vector
#' @param A Treatment vector (must be 0 or 1)
#' @param pi_hat Propensity scalar or vector
#' @param K_folds Number of cross-fitting folds
fit_dr_pseudo_outcomes <- function(X, Y, A, pi_hat = 0.5, K_folds = 5, seed = 123) {
  set.seed(seed)
  df_dr <- data.frame(Y = Y, A = A, X)
  fold_ids <- sample(rep(1:K_folds, length.out = nrow(df_dr)))
  pseudo_Y <- numeric(nrow(df_dr))
  
  sl_library <- c("SL.glm", "SL.glmnet", "SL.gam", "SL.xgboost", "SL.ranger")
  
  for (k in 1:K_folds) {
    train_idx <- which(fold_ids != k); val_idx <- which(fold_ids == k)
    
    X_train_1 <- df_dr[train_idx[df_dr$A[train_idx] == 1], setdiff(names(df_dr), c("Y", "A")), drop=FALSE]
    Y_train_1 <- df_dr$Y[train_idx[df_dr$A[train_idx] == 1]]
    X_train_0 <- df_dr[train_idx[df_dr$A[train_idx] == 0], setdiff(names(df_dr), c("Y", "A")), drop=FALSE]
    Y_train_0 <- df_dr$Y[train_idx[df_dr$A[train_idx] == 0]]
    
    sl_mu_1 <- SuperLearner(Y = Y_train_1, X = X_train_1, SL.library = sl_library, family = gaussian(), cvControl = list(V = 2))
    sl_mu_0 <- SuperLearner(Y = Y_train_0, X = X_train_0, SL.library = sl_library, family = gaussian(), cvControl = list(V = 2))
    
    X_val <- df_dr[val_idx, setdiff(names(df_dr), c("Y", "A")), drop=FALSE]
    mu_1_hat <- predict(sl_mu_1, newdata = X_val)$pred
    mu_0_hat <- predict(sl_mu_0, newdata = X_val)$pred
    mu_A_hat <- ifelse(df_dr$A[val_idx] == 1, mu_1_hat, mu_0_hat)
    
    pi_val <- if (length(pi_hat) > 1) pi_hat[val_idx] else pi_hat
    pseudo_Y[val_idx] <- ((df_dr$A[val_idx] - pi_val) / (pi_val * (1 - pi_val))) * (df_dr$Y[val_idx] - mu_A_hat) + (mu_1_hat - mu_0_hat)
  }
  return(pseudo_Y)
}

#' Simulate continuous maximal deviation process for DR-Learner Frequentist Bounds
compute_frequentist_c1 <- function(dr_lm, valid_coefs, X_grid_raw, num_sims = 2000) {
  dr_beta <- coef(dr_lm)[valid_coefs]
  X_mat_train <- model.matrix(dr_lm)[, valid_coefs]
  XTX_inv <- solve(t(X_mat_train) %*% X_mat_train)
  nu <- df.residual(dr_lm)
  
  sim_Z <- mvrnorm(n = num_sims, mu = rep(0, length(dr_beta)), Sigma = XTX_inv)
  sim_chisq <- rchisq(n = num_sims, df = nu)
  
  dr_grid_df <- as.data.frame(X_grid_raw)
  colnames(dr_grid_df) <- colnames(dr_lm$model)[-1]
  dr_grid_df$pseudo_Y <- 0
  X_mat_grid_dr <- model.matrix(pseudo_Y ~ ., data = dr_grid_df)[, valid_coefs]
  se_grid_dr <- sqrt(rowSums((X_mat_grid_dr %*% XTX_inv) * X_mat_grid_dr))
  
  max_devs_dr <- numeric(num_sims)
  for (s in 1:num_sims) {
    devs <- (X_mat_grid_dr %*% sim_Z[s, ]) / (sqrt(sim_chisq[s] / nu) * se_grid_dr)
    max_devs_dr[s] <- max(devs)
  }
  return(quantile(max_devs_dr, 0.95))
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
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

scenarios <- list(
  list(name="Homogeneous_tau", is_te_hetero=FALSE, test=FALSE),
  list(name="Heterogeneous_Simple_tau", is_te_hetero=TRUE, test=FALSE),
  list(name="Heterogeneous_Complex_tau", is_te_hetero=TRUE, test=TRUE)
)

sample_sizes <- c(150, 300, 600) # Small, Medium, Large

results_list <- list()
bounds_plot_data <- list()
metrics_list <- list()
ranking_list <- list()

# ==============================================================================
# MAIN EXECUTION LOOP
# ==============================================================================
for (n in sample_sizes) {
  for (scen in scenarios) {
    cat(sprintf("\n======================================================\n"))
    cat(sprintf("--- Scenario: %s | N = %d ---\n", scen$name, n))
    cat(sprintf("======================================================\n"))
    
    # 1. Data Generation
    data <- generate_data_2(n = n, is_te_hetero = scen$is_te_hetero, test = scen$test, seed = 123, RCT = TRUE, z_diff = TRUE)
    X_train <- as.matrix(sapply(data[, c(1:6)], as.numeric))
    y_train <- as.numeric(data$y)
    Z_train <- as.numeric(data$z) # Because z_diff=TRUE, Z_train is already centered (-0.5, 0.5)
    propensity_train <- as.numeric(data$pi_x)
    
    # Generate dynamic continuous grid correctly bounded by actual generated data
    grid_size <- 3000
    grid_x1 <- runif(grid_size, min(X_train[,1]), max(X_train[,1]))
    grid_x2 <- runif(grid_size, min(X_train[,2]), max(X_train[,2]))
    grid_x3 <- runif(grid_size, min(X_train[,3]), max(X_train[,3]))
    grid_x4 <- sample(unique(X_train[,4]), grid_size, replace=TRUE)
    idx <- sample(1:nrow(X_train), grid_size, replace=TRUE)
    X_grid_raw <- cbind(grid_x1, grid_x2, grid_x3, grid_x4, X_train[idx, 5], X_train[idx, 6])
    colnames(X_grid_raw) <- colnames(X_train)
    propensity_grid <- rep(0.5, grid_size)
    
    # 2. Fit Models
    cat("  Fitting Semi-Parametric BART model...\n")
    nbcf_fit <- bcf_linear_probit(X_train = X_train, y_train = y_train, Z_train = Z_train, propensity_train = propensity_train,
                                  num_gfr = 50, num_burnin = 2000, num_mcmc = 4000, general_params = general_params_nbcf)
    
    cat("  Fitting DR-Learner (SuperLearner)...\n")
    # DR Learner expects A in {0, 1}, so we pass Z_train + 0.5
    pseudo_Y <- fit_dr_pseudo_outcomes(X_train, y_train, Z_train + 0.5, pi_hat = 0.5, K_folds = 5, seed = 123)
    dr_data <- as.data.frame(X_train); dr_data$pseudo_Y <- pseudo_Y
    dr_lm <- lm(pseudo_Y ~ ., data = dr_data)
    valid_coefs <- !is.na(coef(dr_lm))
    
    # 3. Extract Posteriors
    cat("  Extracting Semi-Parametric BART Posteriors...\n")
    tau_post_obs  <- predict_linear_bcf_patched(nbcf_fit, X = X_train, Z = rep(0, nrow(X_train)))$tau_hat
    tau_post_grid <- predict_linear_bcf_patched(nbcf_fit, X = X_grid_raw, Z = rep(0, nrow(X_grid_raw)))$tau_hat
    
    # 4. Calculate Simultaneous Critical Constants
    cat("  Calculating Critical Constants...\n")
    mu_hat_obs <- rowMeans(tau_post_obs); se_hat_obs <- apply(tau_post_obs, 1, sd); se_hat_obs <- ifelse(se_hat_obs < 1e-8, 1e-8, se_hat_obs)
    c1_bayesian <- quantile(apply((mu_hat_obs - tau_post_obs) / se_hat_obs, 2, max), 0.95)
    
    mu_hat_grid <- rowMeans(tau_post_grid); se_hat_grid <- apply(tau_post_grid, 1, sd); se_hat_grid <- ifelse(se_hat_grid < 1e-8, 1e-8, se_hat_grid)
    c1_inferential <- quantile(apply((mu_hat_grid - tau_post_grid) / se_hat_grid, 2, max), 0.95)
    
    confmarg_obs  <- findconfglob(tau_post_obs, tol = 0.05)
    confmarg_grid <- findconfglob(tau_post_grid[sample(1:nrow(tau_post_grid), 1000), ], tol = 0.05)
    
    c1_dr_continuous <- compute_frequentist_c1(dr_lm, valid_coefs, X_grid_raw, num_sims = 2000)
    
    # 5. Evaluate Metrics on Observed Data
    cat("  Evaluating Metrics...\n")
    true_tau <- data$tau
    lb_bayesian_obs <- mu_hat_obs - c1_bayesian * se_hat_obs
    lb_inferential_obs <- mu_hat_obs - c1_inferential * se_hat_obs
    lb_giorgio_obs <- intfrompost1(tau_post_obs, conf = confmarg_obs)[[1]][, 1]
    lb_giorgio_cont_obs <- intfrompost1(tau_post_obs, conf = confmarg_grid)[[1]][, 1]
    
    dr_obs_df <- as.data.frame(X_train); dr_obs_df$pseudo_Y <- 0
    X_mat_obs_dr <- model.matrix(pseudo_Y ~ ., data = dr_obs_df)[, valid_coefs]
    mu_hat_dr_obs <- as.vector(X_mat_obs_dr %*% coef(dr_lm)[valid_coefs])
    dr_XTX_inv <- solve(t(model.matrix(dr_lm)[, valid_coefs]) %*% model.matrix(dr_lm)[, valid_coefs])
    se_obs_dr <- summary(dr_lm)$sigma * sqrt(rowSums((X_mat_obs_dr %*% dr_XTX_inv) * X_mat_obs_dr))
    lb_dr_obs <- mu_hat_dr_obs - c1_dr_continuous * se_obs_dr
    
    bounds_list <- list("Sechidis (Observed)" = lb_bayesian_obs, "Sechidis (Continuous)" = lb_inferential_obs,
                        "Giorgio (Observed)" = lb_giorgio_obs, "Giorgio (Continuous)" = lb_giorgio_cont_obs,
                        "DR-Learner GLM (Continuous)" = lb_dr_obs)
    estimates_list <- list("Semi-Parametric BART" = mu_hat_obs, "DR-Learner GLM" = mu_hat_dr_obs)
    
    if (scen$is_te_hetero) {
      thresholds <- c("Top 50%" = quantile(true_tau, 0.50), "Top 25%" = quantile(true_tau, 0.75), "Top 10%" = quantile(true_tau, 0.90))
      for (thresh_name in names(thresholds)) {
        tau_t <- thresholds[[thresh_name]]
        S_true <- which(true_tau >= tau_t)
        for (m_name in names(bounds_list)) {
          S_hat <- which(bounds_list[[m_name]] >= tau_t)
          detected_size <- length(S_hat)
          metrics_list[[length(metrics_list) + 1]] <- data.frame(
            Scenario = scen$name, N = n, Threshold = thresh_name, Method = m_name, Detected = detected_size, 
            Precision = if (detected_size > 0) sum(S_hat %in% S_true) / detected_size else NA,
            Recall = if (detected_size > 0) sum(S_hat %in% S_true) / length(S_true) else 0
          )
        }
      }
      for (est_name in names(estimates_list)) {
        ranking_list[[length(ranking_list) + 1]] <- data.frame(Scenario = scen$name, N = n, Method = est_name, Spearman_Cor = cor(true_tau, estimates_list[[est_name]], method = "spearman"))
      }
    }
    
    results_list[[paste(scen$name, n, sep="_")]] <- data.frame(
      Scenario = scen$name, N = n, Metric = c("c1 (Semi-Parametric BART)", "c1 (Semi-Parametric BART)", "confmarg (Tolerance)", "confmarg (Tolerance)", "c1 (DR-Learner GLM)"),
      Method = c("Observed Data", "Continuous Space", "Observed Data", "Continuous Space", "Continuous Space"),
      Value = c(c1_bayesian, c1_inferential, confmarg_obs, confmarg_grid, c1_dr_continuous)
    )
    
    # 6. Extract Plot Grid Data
    if (scen$name == "Heterogeneous_Complex_tau") {
      plot_x2 <- seq(min(X_train[,2]), max(X_train[,2]), length.out = 100)
      X_plot_raw <- cbind(rep(0, 100), plot_x2, rep(0, 100), rep(1, 100), rep(1, 100), rep(0, 100))
      colnames(X_plot_raw) <- colnames(X_train)
      
      tau_post_plot <- predict_linear_bcf_patched(nbcf_fit, X = X_plot_raw, Z = rep(0, nrow(X_plot_raw)))$tau_hat
      tau_hat_plot <- rowMeans(tau_post_plot); se_plot <- apply(tau_post_plot, 1, sd)
      
      dr_plot_df <- as.data.frame(X_plot_raw); colnames(dr_plot_df) <- colnames(dr_lm$model)[-1]; dr_plot_df$pseudo_Y <- 0
      X_mat_plot_dr <- model.matrix(pseudo_Y ~ ., data = dr_plot_df)[, valid_coefs]
      tau_hat_plot_dr <- as.vector(X_mat_plot_dr %*% coef(dr_lm)[valid_coefs])
      se_plot_dr <- summary(dr_lm)$sigma * sqrt(rowSums((X_mat_plot_dr %*% dr_XTX_inv) * X_mat_plot_dr))
      
      base_df <- data.frame(x2 = plot_x2, tau_hat = tau_hat_plot, tau_hat_dr = tau_hat_plot_dr, True_tau = 1 + 2 * plot_x2 * 1, N = n)
      bounds_plot_data[[paste(scen$name, n, sep="_")]] <- bind_rows(
        base_df %>% mutate(Method = "Sechidis (Observed)", LB = tau_hat_plot - c1_bayesian * se_plot, UB = tau_hat_plot + c1_bayesian * se_plot),
        base_df %>% mutate(Method = "Sechidis (Continuous)", LB = tau_hat_plot - c1_inferential * se_plot, UB = tau_hat_plot + c1_inferential * se_plot),
        base_df %>% mutate(Method = "Giorgio (Observed)", LB = intfrompost1(tau_post_plot, conf = confmarg_obs)[[1]][, 1], UB = intfrompost1(tau_post_plot, conf = confmarg_obs)[[1]][, 2]),
        base_df %>% mutate(Method = "Giorgio (Continuous)", LB = intfrompost1(tau_post_plot, conf = confmarg_grid)[[1]][, 1], UB = intfrompost1(tau_post_plot, conf = confmarg_grid)[[1]][, 2]),
        base_df %>% mutate(Method = "DR-Learner GLM (Continuous)", LB = tau_hat_plot_dr - c1_dr_continuous * se_plot_dr, UB = tau_hat_plot_dr + c1_dr_continuous * se_plot_dr)
      )
    }
  }
}

# ==============================================================================
# VISUALIZATIONS & OUTPUTS
# ==============================================================================
results_df <- bind_rows(results_list)
df_c1 <- results_df %>% filter(grepl("c1", Metric))
df_tol <- results_df %>% filter(Metric == "confmarg (Tolerance)")

p1 <- ggplot(df_c1, aes(x = factor(N), y = Value, fill = paste(Metric, Method, sep="\n"))) +
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width=0.7) +
  facet_wrap(~Scenario) + theme_minimal() +
  labs(title = "Comparison of Critical Constant (c1) by Sample Size", y = "Critical Constant c1", x = "Sample Size (N)", fill="Model/Space") +
  scale_fill_brewer(palette = "Set2")
ggsave("R/illustration_c1_plot.png", plot = p1, width = 12, height = 6, bg="white")

p2 <- ggplot(df_tol, aes(x = factor(N), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width=0.7) +
  facet_wrap(~Scenario) + theme_minimal() +
  labs(title = "Comparison of Tolerance Interval Margin (confmarg) by Sample Size", y = "confmarg", x = "Sample Size (N)", fill="Space") +
  scale_fill_manual(values = c("Observed Data" = "skyblue", "Continuous Space" = "coral"))
ggsave("R/illustration_tol_plot.png", plot = p2, width = 12, height = 6, bg="white")

plot_df <- bind_rows(bounds_plot_data)
p3 <- ggplot(plot_df, aes(x = x2)) +
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = Method), alpha = 0.4) +
  geom_line(aes(y = True_tau, color = "True tau"), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = tau_hat, color = "Estimated tau (Semi-Parametric BART)"), linewidth = 1) +
  geom_line(aes(y = tau_hat_dr, color = "Estimated tau (DR GLM)"), linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") + 
  facet_grid(Method ~ N, labeller = label_both) + theme_minimal() +
  labs(title = "Simultaneous Confidence Intervals for Subgroup Selection (Complex Scenario, x4 = 1)", y = "Treatment Effect (CATE)", x = "x2") +
  scale_fill_manual(values = c("Sechidis (Observed)"="coral", "Sechidis (Continuous)"="red", "Giorgio (Observed)"="lightgreen", "Giorgio (Continuous)"="darkgreen", "DR-Learner GLM (Continuous)"="purple")) +
  scale_color_manual(values = c("True tau"="black", "Estimated tau (Semi-Parametric BART)"="blue", "Estimated tau (DR GLM)"="darkblue"))
ggsave("R/illustration_bounds_plot.png", plot = p3, width = 16, height = 12, bg="white")

metrics_df <- bind_rows(metrics_list)
if (nrow(metrics_df) > 0) {
  metrics_df$Threshold <- factor(metrics_df$Threshold, levels = c("Top 50%", "Top 25%", "Top 10%"))
  p4 <- ggplot(metrics_df, aes(x = Recall, y = Precision, color = Method)) +
    geom_line(aes(group = Method), alpha = 0.5, linetype = "dashed") +
    geom_point(aes(shape = factor(N)), size = 4, alpha = 0.8) +
    facet_grid(Scenario ~ Threshold) + theme_minimal() +
    labs(title = "Precision vs. Recall across Thresholds", x = "Recall (Power)", y = "Precision (1 - False Discovery Rate)") +
    scale_color_manual(values = c("Sechidis (Observed)"="coral", "Sechidis (Continuous)"="red", "Giorgio (Observed)"="lightgreen", "Giorgio (Continuous)"="darkgreen", "DR-Learner GLM (Continuous)"="purple"))
  ggsave("R/illustration_metrics_plot.png", plot = p4, width = 14, height = 8, bg="white")
}

ranking_df <- bind_rows(ranking_list)
save(results_df, bounds_plot_data, metrics_df, ranking_df, file = "R/illustration_results.RData")
cat("\nResults successfully saved to R/illustration_results.RData\n")
