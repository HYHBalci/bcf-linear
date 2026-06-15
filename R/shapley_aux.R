# ==============================================================================
# MASTER SHAPLEY & INTERACTION ANALYSIS SCRIPT
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. REQUIRED LIBRARIES
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(forcats)
library(scales)

# ------------------------------------------------------------------------------
# 2. HELPER & MATHEMATICAL FUNCTIONS
# ------------------------------------------------------------------------------

#' Generate Interaction Pairs Based on a Boolean Vector
#' (Fallback function if the model does not provide them natively)
interaction_pairs <- function(num_covariates, boolean_vector) {
  interaction_list <- list()
  if (num_covariates > 1) {
    for (j in 1:(num_covariates - 1)) {
      for (k in (j + 1):num_covariates) {
        # || ensures continuous x categorical interactions are kept
        if (boolean_vector[j] || boolean_vector[k]) {
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
        }
      }
    }
  }
  if (length(interaction_list) == 0) return(matrix(numeric(0), nrow = 2, ncol = 0))
  return(do.call(cbind, interaction_list))
}

# --- Giorgio's Tolerance Interval Functions ---
postsum <- function(posters, conf = 0.05) {
  qs <- quantile(posters, p = c(conf/2, 1 - conf/2))
  mn <- mean(posters)
  return(c(qs[1], mn, qs[2]))
}

intfrompost1 <- function(post, conf = 0.05) {
  shtotpost <- t(apply(post, 1, postsum, conf = conf))
  return(list(tot = shtotpost))
}

conftol <- function(margconf = 0.05, post, tol = 0.05) {
  int <- intfrompost1(post, conf = margconf)[[1]]
  difl <- ((post - int[,1]) >= 0) * ((post - int[,3]) <= 0)
  cs <- colMeans(difl)
  return(mean(cs >= 1 - tol))
} 

findconfglob <- function(post, tol = 0.05) {
  grid <- seq(0.0005, 0.05, by = 0.0005)
  conftols <- sapply(grid, conftol, post = post, tol = tol)
  wm <- which(conftols - 0.95 < 0)[1] - 1
  return(grid[wm])
}

# ------------------------------------------------------------------------------
# 3. CORE SHAPLEY ENGINES
# ------------------------------------------------------------------------------

#' Highly Optimized & Vectorized Shapley Calculator for a Single Observation
shapley_lean <- function(index, x_star, cov_X, beta_samples, beta_int_samples, ipairs, p_int, num_covariates, col_names) {
  num_samples <- nrow(beta_samples)
  shapley_main_vals <- matrix(0, nrow = num_samples, ncol = num_covariates)
  shapley_int_vals <- matrix(0, nrow = num_samples, ncol = num_covariates)
  
  for (j in 1:num_covariates) {
    main_term <- beta_samples[, j] * x_star[j]
    interaction_sum <- numeric(num_samples)
    
    if (p_int > 0) {
      involved <- which(ipairs[1, ] == j | ipairs[2, ] == j)
      if (length(involved) > 0) {
        k_vec <- ifelse(ipairs[1, involved] == j, ipairs[2, involved], ipairs[1, involved])
        delta_vec <- (x_star[j] * x_star[k_vec]) - cov_X[j, k_vec]
        # Vectorized Matrix Multiplication to bypass loop
        interaction_sum <- as.numeric(beta_int_samples[, involved, drop = FALSE] %*% delta_vec)
      }
    }
    
    shapley_main_vals[, j] <- main_term
    shapley_int_vals[, j] <- 0.5 * interaction_sum
  }
  
  # Calculate medians directly to drop heavy 5000-draw matrices from RAM
  shapley_main_med <- apply(shapley_main_vals, 2, median)
  shapley_int_med  <- apply(shapley_int_vals, 2, median)
  shapley_total_med <- shapley_main_med + shapley_int_med
  
  tibble::tibble(
    obs = index,
    feature = col_names,
    feature_value = x_star,
    shapley_total = shapley_total_med,
    shapley_main = shapley_main_med,
    shapley_int = shapley_int_med
  )
}

#' Optimized Global Shapley Calculator
compute_shapley_all <- function(X, beta_post, beta_int_post, ipairs, indices = NULL) {
  if (is.null(indices)) indices <- 1:nrow(X)
  num_covariates <- ncol(X)
  col_names <- make.names(colnames(X), unique = TRUE)
  
  cat("Pre-computing covariance and flattening arrays...\n")
  cov_X_precomputed <- cov(X)
  if (nrow(ipairs) != 2 && ncol(ipairs) == 2) ipairs <- t(ipairs)
  p_int <- ncol(ipairs)
  
  is_3D <- length(dim(beta_post)) == 3
  if (is_3D) {
    beta_samples <- do.call(rbind, lapply(1:dim(beta_post)[1], function(chain) beta_post[chain, , ]))
    beta_int_samples <- if (p_int > 0) do.call(rbind, lapply(1:dim(beta_int_post)[1], function(chain) beta_int_post[chain, , ])) else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  } else {
    beta_samples <- beta_post
    beta_int_samples <- if (p_int > 0) beta_int_post else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  }
  
  if (p_int > 0 && p_int != ncol(beta_int_samples)) {
    stop(sprintf("Dimension Mismatch! interaction_pairs has %d pairs, but the posterior has %d columns.", p_int, ncol(beta_int_samples)))
  }
  
  cat("Calculating Shapley medians for", length(indices), "observations...\n")
  all_shapleys <- purrr::map_dfr(indices, function(idx) {
    shapley_lean(
      index = idx, x_star = X[idx, ], cov_X = cov_X_precomputed, 
      beta_samples = beta_samples, beta_int_samples = beta_int_samples, 
      ipairs = ipairs, p_int = p_int, num_covariates = num_covariates, col_names = col_names
    )
  }, .progress = TRUE)
  
  return(all_shapleys)
}

#' Extract Full 5000-Draw Shapley Posterior for a Single Feature (For Tolerance Intervals)
get_shapley_posterior_single <- function(feature_name, indices, X, beta_post, beta_int_post, ipairs) {
  j <- match(feature_name, colnames(X))
  if (is.na(j)) stop(paste("Feature", feature_name, "not found in X"))
  
  cov_X_precomputed <- cov(X)
  if (nrow(ipairs) != 2 && ncol(ipairs) == 2) ipairs <- t(ipairs)
  p_int <- ncol(ipairs)
  
  is_3D <- length(dim(beta_post)) == 3
  if (is_3D) {
    beta_samples <- do.call(rbind, lapply(1:dim(beta_post)[1], function(chain) beta_post[chain, , ]))
    beta_int_samples <- if (p_int > 0) do.call(rbind, lapply(1:dim(beta_int_post)[1], function(chain) beta_int_post[chain, , ])) else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  } else {
    beta_samples <- beta_post
    beta_int_samples <- if (p_int > 0) beta_int_post else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  }
  
  num_samples <- nrow(beta_samples)
  post_matrix <- matrix(0, nrow = length(indices), ncol = num_samples)
  
  involved <- which(ipairs[1, ] == j | ipairs[2, ] == j)
  k_vec <- if(length(involved) > 0) ifelse(ipairs[1, involved] == j, ipairs[2, involved], ipairs[1, involved]) else integer(0)
  
  cat("Extracting full posterior for feature:", feature_name, "...\n")
  for (idx_pos in seq_along(indices)) {
    x_star <- X[indices[idx_pos], ]
    main_term <- beta_samples[, j] * x_star[j]
    if (length(involved) > 0) {
      delta_vec <- (x_star[j] * x_star[k_vec]) - cov_X_precomputed[j, k_vec]
      interaction_sum <- as.numeric(beta_int_samples[, involved, drop = FALSE] %*% delta_vec)
    } else {
      interaction_sum <- 0
    }
    post_matrix[idx_pos, ] <- main_term + 0.5 * interaction_sum
  }
  return(post_matrix)
}

# ------------------------------------------------------------------------------
# 4. PLOTTING SUITE
# ------------------------------------------------------------------------------

#' Plot 1: SHAP Summary (Directional Contribution & Distribution)
plot_shapley_contribution <- function(shapley_data) {
  plot_data <- shapley_data %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(scaled_value = scales::rescale(feature_value)) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = shapley_total, y = forcats::fct_reorder(feature, abs(shapley_total), .fun = mean))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
    ggplot2::geom_jitter(ggplot2::aes(color = scaled_value), height = 0.15, alpha = 0.7, size = 1.5) +
    ggplot2::scale_color_gradient2(
      low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c", midpoint = 0.5,
      name = "Feature Value\n(Low to High)", breaks = c(0, 1), labels = c("Low", "High")
    ) +
    ggplot2::labs(title = "SHAP Summary: Feature Contribution & Distribution", x = "Total Shapley Value", y = "Feature") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot 2: Global Importance Breakdown (Main vs Interaction)
plot_shapley_importance_breakdown <- function(shapley_data) {
  importance_data <- shapley_data %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(Main = mean(abs(shapley_main)), Interaction = mean(abs(shapley_int)), Total_Height = Main + Interaction, .groups = "drop") %>%
    tidyr::pivot_longer(cols = c(Main, Interaction), names_to = "Effect_Type", values_to = "Importance")
  
  ggplot2::ggplot(importance_data, ggplot2::aes(x = Importance, y = forcats::fct_reorder(feature, Total_Height), fill = Effect_Type)) +
    ggplot2::geom_col(position = "stack", alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("Main" = "#4d4d4d", "Interaction" = "#e0e0e0"), name = "Effect Type") +
    ggplot2::labs(title = "Global Feature Importance Breakdown", subtitle = "Mean absolute Shapley values (I_j)", x = "Mean Absolute Contribution", y = "Feature") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom")
}

#' Plot 3: Credible vs. Tolerance Intervals for a Single Feature
plot_tolerance_intervals <- function(post_matrix, feature_values, feature_name) {
  cat("Calculating optimal global confidence level...\n")
  confglob <- findconfglob(post_matrix, tol = 0.02)
  cred_int <- intfrompost1(post_matrix, conf = 0.05)[[1]]
  tol_int  <- intfrompost1(post_matrix, conf = confglob)[[1]]
  
  plot_df <- data.frame(cov_val = feature_values, mean = cred_int[, 2], cred_low = cred_int[, 1], cred_high = cred_int[, 3], tol_low = tol_int[, 1], tol_high = tol_int[, 3]) %>%
    dplyr::arrange(cov_val) %>%
    dplyr::mutate(id_ordered = row_number())
  
  ggplot(plot_df, aes(x = id_ordered)) +
    geom_linerange(aes(ymin = tol_low, ymax = tol_high, color = "Tolerance interval"), linewidth = 0.6, linetype = "dashed", alpha = 0.8) +
    geom_linerange(aes(ymin = cred_low, ymax = cred_high, color = "Credible interval"), linewidth = 1.2, alpha = 0.9) +
    geom_point(aes(y = mean), size = 1.2, shape = 15, color = "black") +
    scale_color_manual(values = c("Credible interval" = "#2b8cbe", "Tolerance interval" = "#de2d26")) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    labs(title = paste("Shapley", feature_name), x = "Individuals (ordered by feature value)", y = "Shapley Value") +
    theme_minimal(base_size = 14) +
    theme(legend.position = c(0.15, 0.9), legend.title = element_blank(), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
}

#' Extract Full 5000-Draw Shapley Posterior Dissected into Main and Interaction Effects
get_shapley_posterior_dissected <- function(feature_name, indices, X, beta_post, beta_int_post, ipairs) {
  j <- match(feature_name, colnames(X))
  if (is.na(j)) stop(paste("Feature", feature_name, "not found in X"))
  
  cov_X_precomputed <- cov(X)
  if (nrow(ipairs) != 2 && ncol(ipairs) == 2) ipairs <- t(ipairs)
  p_int <- ncol(ipairs)
  
  is_3D <- length(dim(beta_post)) == 3
  if (is_3D) {
    beta_samples <- do.call(rbind, lapply(1:dim(beta_post)[1], function(chain) beta_post[chain, , ]))
    beta_int_samples <- if (p_int > 0) do.call(rbind, lapply(1:dim(beta_int_post)[1], function(chain) beta_int_post[chain, , ])) else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  } else {
    beta_samples <- beta_post
    beta_int_samples <- if (p_int > 0) beta_int_post else matrix(0, nrow = nrow(beta_samples), ncol = 0)
  }
  
  num_samples <- nrow(beta_samples)
  
  # Initialize separate matrices for main and interaction terms
  main_matrix <- matrix(0, nrow = length(indices), ncol = num_samples)
  int_matrix <- matrix(0, nrow = length(indices), ncol = num_samples)
  
  involved <- which(ipairs[1, ] == j | ipairs[2, ] == j)
  k_vec <- if(length(involved) > 0) ifelse(ipairs[1, involved] == j, ipairs[2, involved], ipairs[1, involved]) else integer(0)
  
  cat("Extracting dissected posterior for feature:", feature_name, "...\n")
  for (idx_pos in seq_along(indices)) {
    x_star <- X[indices[idx_pos], ]
    
    # 1. Main Effect
    main_term <- beta_samples[, j] * x_star[j]
    
    # 2. Interaction Effect
    if (length(involved) > 0) {
      delta_vec <- (x_star[j] * x_star[k_vec]) - cov_X_precomputed[j, k_vec]
      interaction_sum <- as.numeric(beta_int_samples[, involved, drop = FALSE] %*% delta_vec)
    } else {
      interaction_sum <- rep(0, num_samples)
    }
    
    main_matrix[idx_pos, ] <- main_term
    int_matrix[idx_pos, ] <- 0.5 * interaction_sum
  }
  
  return(list(
    main = main_matrix, 
    interaction = int_matrix
  ))
}