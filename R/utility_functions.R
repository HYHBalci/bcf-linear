# ==============================================================================
# utility_functions.R
# Shared helper functions for the extended competition pipeline
# ==============================================================================
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. TOLERANCE BANDS / SIGNIFICANCE HELPERS (from ACTG_complete.R)
# ------------------------------------------------------------------------------
# 1. TOLERANCE BANDS / SIGNIFICANCE HELPERS (from ACTG_complete.R)
# ------------------------------------------------------------------------------
safe_row_quantiles <- function(mat, probs) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  t(apply(mat, 1, quantile, probs = probs, na.rm = TRUE))
}

conftol <- function(margconf = 0.05, post, tol = 0.05) {
  lower <- safe_row_quantiles(post, margconf / 2)
  upper <- safe_row_quantiles(post, 1 - margconf / 2)
  difl <- (post >= as.vector(lower)) & (post <= as.vector(upper))
  cs <- colMeans(difl, na.rm = TRUE)
  return(mean(cs >= 1 - tol, na.rm = TRUE))
}

findconfglob <- function(post, tol = 0.05) {
  grid <- seq(0.0005, 0.05, by = 0.0005)
  conftols <- sapply(grid, conftol, post = post, tol = tol)
  wm <- which(conftols - 0.95 < 0)[1] - 1
  if (is.na(wm) || wm < 1) return(0.0005) 
  return(grid[wm])
}

# ------------------------------------------------------------------------------
# 2. INTERVAL CALCULATION & SIGNIFICANCE CHECK
# ------------------------------------------------------------------------------
check_significant_individuals <- function(het_draws, level = 0.95, tol = 0.05) {
  alpha_val <- 1 - level
  n_draws <- ncol(het_draws)
  
  # 1. Tolerance Interval (Globally adjusted)
  confglob <- findconfglob(het_draws, tol = tol)
  tol_lower <- safe_row_quantiles(het_draws, confglob / 2)
  tol_upper <- safe_row_quantiles(het_draws, 1 - confglob / 2)
  
  sig_tol <- sum(tol_lower > 0 | tol_upper < 0, na.rm = TRUE)
  
  # 2. Wan & Sechidis Intervals
  mu <- rowMeans(het_draws)
  sd_vals <- apply(het_draws, 1, sd)
  
  t_val <- qt(1 - alpha_val / 2, df = n_draws - 1)
  wan_se <- sd_vals / sqrt(n_draws)
  wan_lower <- mu - t_val * wan_se
  wan_upper <- mu + t_val * wan_se
  sig_wan <- sum(wan_lower > 0 | wan_upper < 0)
  
  z_val <- qnorm(1 - alpha_val / 2)
  sechidis_se <- sqrt(sd_vals^2 + mu^2 / n_draws)
  sech_lower <- mu - z_val * sechidis_se
  sech_upper <- mu + z_val * sechidis_se
  sig_sechidis <- sum(sech_lower > 0 | sech_upper < 0)
  
  return(list(
    tolerance_sig = sig_tol,
    wan_sig = sig_wan,
    sechidis_sig = sig_sechidis,
    n_total = nrow(het_draws)
  ))
}

# ------------------------------------------------------------------------------
# 3. AUUC PERMUTATION TEST
# ------------------------------------------------------------------------------
calc_auc_helper <- function(df, metric) {
  sum(diff(df$frac) * (head(df[[metric]], -1) + tail(df[[metric]], -1)) / 2)
}

auuc_permutation_test <- function(eval_df, model_col, B = 1000) {
  actual_uplift <- get_eval_curves_no_model(eval_df, model_col)
  actual_auuc <- calc_auc_helper(actual_uplift, "uplift")
  
  perm_auucs <- numeric(B)
  eval_perm <- eval_df
  for (b in 1:B) {
    eval_perm[[model_col]] <- sample(eval_df[[model_col]])
    perm_uplift <- get_eval_curves_no_model(eval_perm, model_col)
    perm_auucs[b] <- calc_auc_helper(perm_uplift, "uplift")
  }
  
  p_val <- mean(perm_auucs >= actual_auuc)
  return(list(actual_auuc = actual_auuc, p_val = p_val))
}

get_eval_curves_no_model <- function(df, model_col) {
  df %>%
    arrange(desc(.data[[model_col]])) %>%
    mutate(
      cum_Z = cumsum(Z), cum_C = cumsum(1 - Z),
      cum_Y_Z = cumsum(Y * Z), cum_Y_C = cumsum(Y * (1 - Z)),
      cum_ATE = ifelse(cum_Z > 0 & cum_C > 0, (cum_Y_Z / cum_Z) - (cum_Y_C / cum_C), 0),
      uplift = cum_ATE * (cum_Z + cum_C), frac = row_number() / n()
    ) %>%
    dplyr::select(frac, uplift)
}
