library(BayesLogit)
library(MASS)
library(stochtree)

create_interaction_pairs_R <- function(p_main_local) {
  pairs_local <- list()
  if (p_main_local < 2) return(list())
  idx <- 1
  for (j_idx in 1:(p_main_local - 1)) {
    for (k_idx in (j_idx + 1):p_main_local) {
      pairs_local[[idx]] <- c(j_idx, k_idx)
      idx <- idx + 1
    }
  }
  return(pairs_local)
}

rinvgamma_R <- function(n, shape, scale) {
  if (shape <= 0 || scale <= 0) {
    if(shape <= 0 && scale <= 0) return(rep(NaN, n))
    if(shape <= 0) return(rep(NaN, n))
    if(scale <= 0) return(rep(Inf, n)) 
  }
  return(1 / stats::rgamma(n, shape = shape, rate = scale))
}

safe_var_R <- function(x) {
  return(pmax(x, 1e-9))
}


#' Fit Custom Logistic or Linear Model with Grouped Horseshoe Priors
#'
#' Model: link(E[y]) = alpha_global + X*beta + XX*beta_int +
#'                      Z*(aleph + X*gamma + XX*gamma_int)
#' - Link can be "logit" (for binomial family) or "identity" (for gaussian family).
#' - alpha_global and aleph get Normal priors.
#' - (beta, beta_int) get Horseshoe prior with global tau_prognostic.
#' - (gamma, gamma_int) get Horseshoe prior with global tau_treatment.
#' - An optional tau_overall can shrink both groups further.
#'
#' @param y_vec Response vector. Binary (0s and 1s) for binomial, continuous for gaussian.
#' @param X_mat Covariate matrix.
#' @param Z_vec Binary treatment indicator vector (0s and 1s).
#' @param family Model family, either "binomial" (default) or "gaussian".
#' @param n_iter Total MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param num_chains (integer) number of chains to be initialized for sampling.  
#' @param method_tau_prognostic Method for global shrinkage for prognostic terms.
#' @param tau_prognostic_init Initial/fixed value for prognostic global shrinkage.
#' @param method_tau_treatment Method for global shrinkage for treatment-related terms.
#' @param tau_treatment_init Initial/fixed value for treatment-related global shrinkage.
#' @param method_tau_overall Method for overall global shrinkage (shrinks both groups).
#' @param tau_overall_init Initial/fixed value for overall global shrinkage.
#' @param alpha_global_prior_mean Mean for Normal prior on alpha_global. Default 0.
#' @param alpha_global_prior_sd SD for Normal prior on alpha_global. Default 10.
#' @param aleph_prior_mean Mean for Normal prior on aleph. Default 0.
#' @param aleph_prior_sd SD for Normal prior on aleph. Default 10.
#' @param thin Thinning parameter.
#' @param seed Random seed.
#' @param verbose Logical, print progress.
#' @param ping Frequency of progress messages.
#'
#' @return A list of MCMC samples.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma rexp runif pgamma qgamma rnorm var
fit_grouped_horseshoes_R <- function(
    y_vec, X_mat, Z_vec,
    family = c("binomial", "gaussian"),
    n_iter, burn_in, num_chains = 1, propensity_as_covariate = T,
    method_tau_prognostic = c("halfCauchy", "truncatedCauchy", "fixed"), tau_prognostic_init = 1.0,
    method_tau_treatment = c("halfCauchy", "truncatedCauchy", "fixed"), tau_treatment_init = 1.0,
    method_tau_overall = c("halfCauchy", "truncatedCauchy", "fixed"), tau_overall_init = 1.0,
    alpha_global_prior_mean = 0.0, alpha_global_prior_sd = 10.0,
    aleph_prior_mean = 0.0, aleph_prior_sd = 10.0,
    thin = 1,
    seed = 123,
    verbose = FALSE,
    ping = 1000
) {
  if (!is.null(seed)) {
    set.seed(seed)
    chain_seeds <- sample.int(1e6, num_chains)
  } else {
    chain_seeds <- rep(NA, num_chains)
  }
  
  family <- match.arg(family)
  method_tau_prognostic <- match.arg(method_tau_prognostic)
  method_tau_treatment <- match.arg(method_tau_treatment)
  method_tau_overall <- match.arg(method_tau_overall)
  
  # --- Input Validation ---
  if (length(unique(na.omit(y_vec))) <= 1) {
    stop("Outcome variable `y_vec` is constant or has no variation. Model cannot be fit.")
  }
  
  N <- nrow(X_mat)
  p_main_X_orig <- ncol(X_mat)
  X_mat_orig <- X_mat # Keep a copy of the original X_mat
  
  # --- Propensity Score Estimation (if requested) ---
  propensity_scores <- NULL
  if (propensity_as_covariate) {
    if (verbose) message("Estimating propensity scores to use as a covariate...")
    ps_num_burnin <- 10
    ps_num_total <- 50
    bart_model_propensity <- bart(X_train = X_mat_orig, y_train = as.numeric(Z_vec), X_test = NULL, 
                                  num_gfr = ps_num_total, num_burnin = 0, num_mcmc = 0)
    propensity_scores <- rowMeans(bart_model_propensity$y_hat_train[,(ps_num_burnin+1):ps_num_total])
    if (verbose) message("Propensity scores calculated.")
  }
  
  # --- 1. Construct Design Matrices (once, outside the chain loop) ---
  # Prognostic Part (beta, beta_int, and potentially propensity score)
  X_prog_list <- list()
  if (p_main_X_orig > 0) {
    if(is.null(colnames(X_mat_orig))) colnames(X_mat_orig) <- paste0("X", 1:p_main_X_orig)
    X_prog_list$beta <- X_mat_orig
  }
  if (!is.null(propensity_scores)) {
    X_prog_list$beta_propensity <- propensity_scores
  }
  
  main_interaction_indices_X <- create_interaction_pairs_R(p_main_X_orig)
  p_interaction_X <- length(main_interaction_indices_X)
  if (p_interaction_X > 0) {
    X_beta_int_mat <- matrix(0, nrow = N, ncol = p_interaction_X)
    colnames_beta_int <- character(p_interaction_X)
    for (k in 1:p_interaction_X) {
      idx_pair <- main_interaction_indices_X[[k]]
      X_beta_int_mat[, k] <- X_mat_orig[, idx_pair[1]] * X_mat_orig[, idx_pair[2]]
      colnames_beta_int[k] <- paste0("beta_int_", paste(colnames(X_mat_orig)[idx_pair], collapse="_"))
    }
    colnames(X_beta_int_mat) <- colnames_beta_int
    X_prog_list$beta_int <- X_beta_int_mat
  }
  X_prognostic <- if(length(X_prog_list) > 0) do.call(cbind, X_prog_list) else matrix(0, nrow=N, ncol=0)
  p_prognostic <- ncol(X_prognostic)
  
  # Treatment-related Part (gamma, gamma_int) - uses ORIGINAL X_mat only
  X_treat_list <- list()
  if (p_main_X_orig > 0) {
    X_gamma_mat <- X_mat_orig * Z_vec 
    colnames(X_gamma_mat) <- paste0("Z_", colnames(X_mat_orig))
    X_treat_list$gamma <- X_gamma_mat
  }
  if (p_interaction_X > 0) {
    # Re-use the prognostic interaction matrix base for consistency
    X_gamma_int_mat <- X_prog_list$beta_int * Z_vec 
    colnames(X_gamma_int_mat) <- paste0("Z_", colnames(X_prog_list$beta_int))
    X_treat_list$gamma_int <- X_gamma_int_mat
  }
  X_treatment_related <- if(length(X_treat_list) > 0) do.call(cbind, X_treat_list) else matrix(0, nrow=N, ncol=0)
  p_treatment_related <- ncol(X_treatment_related)
  
  X_hs <- cbind(X_prognostic, X_treatment_related)
  p_hs_total <- ncol(X_hs)
  
  # --- List to store results from all chains ---
  all_chains_results <- list()
  
  # --- Outer loop for running multiple chains ---
  for (chain_num in 1:num_chains) {
    if(verbose && num_chains > 1) message(sprintf("\n--- Starting Chain %d / %d ---\n", chain_num, num_chains))
    
    if (!is.na(chain_seeds[chain_num])) set.seed(chain_seeds[chain_num])
    
    # --- Initialize Parameters (inside chain loop for independence) ---
    alpha_global <- 0.0
    aleph <- 0.0
    Beta_hs <- rep(0, p_hs_total)    
    lambda_hs <- rep(1, p_hs_total)  
    current_tau_prognostic <- tau_prognostic_init
    current_tau_treatment <- tau_treatment_init
    current_tau_overall <- tau_overall_init
    
    if (family == "gaussian") {
      sigma_sq <- var(y_vec, na.rm = TRUE)
      if (is.na(sigma_sq) || sigma_sq == 0) sigma_sq <- 1.0 # Fallback
    } else { 
      y_adj <- y_vec - 0.5
    }
    
    # --- MCMC Storage for the current chain ---
    eff_samp_count <- floor((n_iter - burn_in) / thin)
    alpha_global_samples <- numeric(eff_samp_count)
    aleph_samples <- numeric(eff_samp_count)
    Beta_hs_samples <- if(p_hs_total > 0) matrix(0, nrow = p_hs_total, ncol = eff_samp_count) else matrix(0,0,0)
    lambda_hs_samples <- if(p_hs_total > 0) matrix(0, nrow = p_hs_total, ncol = eff_samp_count) else matrix(0,0,0)
    tau_prognostic_samples <- numeric(eff_samp_count)
    tau_treatment_samples <- numeric(eff_samp_count)
    tau_overall_samples <- numeric(eff_samp_count)
    if (family == "gaussian") sigma_sq_samples <- numeric(eff_samp_count)
    store_idx <- 0
    
    # --- MCMC Loop for the current chain ---
    if(verbose) message(sprintf("MCMC (%s model) is running...", family))
    for (iter_idx in 1:n_iter) {
      if (verbose && iter_idx %% ping == 0) message("Chain ", chain_num, " - Iteration: ", iter_idx, "/", n_iter)
      
      # (The entire MCMC loop logic from the previous version goes here)
      # This logic remains correct after augmenting X_mat
      # Start of inner MCMC logic
      eta_hs_part <- if(p_hs_total > 0) as.vector(X_hs %*% Beta_hs) else rep(0, N)
      eta <- alpha_global + Z_vec * aleph + eta_hs_part
      if (family == "binomial") {
        omega <- BayesLogit::rpg(num = N, h = 1, z = eta); omega[omega < 1e-9] <- 1e-9; z_aug <- y_adj / omega; w_vec <- omega; scale_factor <- 1.0
      } else { z_aug <- y_vec; w_vec <- rep(1, N); scale_factor <- sigma_sq }
      resid_for_alpha <- z_aug - (Z_vec * aleph + eta_hs_part); prior_prec_alpha <- 1 / safe_var_R(alpha_global_prior_sd^2); post_prec_alpha <- sum(w_vec) / scale_factor + prior_prec_alpha; post_mean_alpha <- (sum(w_vec * resid_for_alpha) / scale_factor + alpha_global_prior_mean * prior_prec_alpha) / post_prec_alpha; alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
      resid_for_aleph <- z_aug - (alpha_global + eta_hs_part); prior_prec_aleph <- 1 / safe_var_R(aleph_prior_sd^2); post_prec_aleph <- sum(w_vec * Z_vec^2) / scale_factor + prior_prec_aleph; post_mean_aleph <- (sum(w_vec * Z_vec * resid_for_aleph) / scale_factor + aleph_prior_mean * prior_prec_aleph) / post_prec_aleph; aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
      if (p_hs_total > 0) {
        z_for_hs <- z_aug - alpha_global - Z_vec * aleph; idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0); idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_hs_total else integer(0)
        effective_global_taus_for_beta <- rep(1, p_hs_total); if(length(idx_prog) > 0) effective_global_taus_for_beta[idx_prog] <- current_tau_prognostic * current_tau_overall; if(length(idx_treat) > 0) effective_global_taus_for_beta[idx_treat] <- current_tau_treatment * current_tau_overall
        temp_lambda_scale_term <- (Beta_hs^2) / (2 * scale_factor * safe_var_R(effective_global_taus_for_beta^2)); nu_lambda_inv <- stats::rexp(p_hs_total, rate = 1 + (1 / safe_var_R(lambda_hs^2))); lambda_hs_sq_inv <- stats::rexp(p_hs_total, rate = nu_lambda_inv + temp_lambda_scale_term); lambda_hs <- 1 / sqrt(safe_var_R(lambda_hs_sq_inv))
        if (length(idx_prog) > 0 && method_tau_prognostic != "fixed") { xi_inv_prog <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_prognostic^2)); rate_prog <- xi_inv_prog + sum(Beta_hs[idx_prog]^2 / (2 * scale_factor * safe_var_R(lambda_hs[idx_prog]^2 * current_tau_overall^2))); current_tau_prognostic <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog)) }
        if (length(idx_treat) > 0 && method_tau_treatment != "fixed") { xi_inv_treat <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_treatment^2)); rate_treat <- xi_inv_treat + sum(Beta_hs[idx_treat]^2 / (2 * scale_factor * safe_var_R(lambda_hs[idx_treat]^2 * current_tau_overall^2))); current_tau_treatment <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat)) }
        tau_vector_components <- rep(1, p_hs_total); if(length(idx_prog) > 0) tau_vector_components[idx_prog] <- current_tau_prognostic; if(length(idx_treat) > 0) tau_vector_components[idx_treat] <- current_tau_treatment
        if (method_tau_overall != "fixed") { xi_inv_overall <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_overall^2)); rate_overall <- xi_inv_overall + sum(Beta_hs^2 / (2 * scale_factor * safe_var_R(lambda_hs^2 * tau_vector_components^2))); current_tau_overall <- sqrt(rinvgamma_R(1, shape = (p_hs_total + 1)/2, scale = rate_overall)) }
        current_tau_prognostic <- max(current_tau_prognostic, 1e-6); current_tau_treatment <- max(current_tau_treatment, 1e-6); current_tau_overall <- max(current_tau_overall, 1e-6); lambda_hs[lambda_hs < 1e-6] <- 1e-6
        prior_var_beta_hs <- safe_var_R(lambda_hs^2 * effective_global_taus_for_beta^2); prior_prec_beta_hs_diag <- diag(1 / prior_var_beta_hs, p_hs_total, p_hs_total)
        XtWX_hs <- t(X_hs) %*% (w_vec * X_hs); post_prec_beta_hs_unscaled <- XtWX_hs + prior_prec_beta_hs_diag
        chol_attempt <- try(chol(post_prec_beta_hs_unscaled), silent = TRUE)
        if (inherits(chol_attempt, "try-error")) { warning(paste("Cholesky for Beta_hs failed at iter", iter_idx, "- skipping Beta_hs update."))
        } else { post_cov_beta_hs_unscaled_inv <- chol2inv(chol_attempt); rhs_mean <- t(X_hs) %*% (w_vec * z_for_hs); post_mean_beta_hs <- post_cov_beta_hs_unscaled_inv %*% rhs_mean; Beta_hs <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean_beta_hs), Sigma = scale_factor * post_cov_beta_hs_unscaled_inv) }
      }
      if (family == "gaussian") { eta_hs_part_new <- if(p_hs_total > 0) as.vector(X_hs %*% Beta_hs) else rep(0, N); eta_new <- alpha_global + Z_vec * aleph + eta_hs_part_new; residuals_for_sigma <- y_vec - eta_new; shape_sigma_post <- 0.001 + N / 2; rate_sigma_post <- 0.001 + 0.5 * sum(residuals_for_sigma^2); sigma_sq <- rinvgamma_R(1, shape = shape_sigma_post, scale = rate_sigma_post); sigma_sq <- max(sigma_sq, 1e-9) }
      # End of inner MCMC logic
      
      # Store samples for the current chain
      if (iter_idx > burn_in && (iter_idx - burn_in) %% thin == 0) {
        store_idx <- store_idx + 1
        alpha_global_samples[store_idx] <- alpha_global
        aleph_samples[store_idx] <- aleph
        if(p_hs_total > 0) {
          Beta_hs_samples[, store_idx] <- Beta_hs
          lambda_hs_samples[, store_idx] <- lambda_hs
        }
        tau_prognostic_samples[store_idx] <- current_tau_prognostic
        tau_treatment_samples[store_idx] <- current_tau_treatment
        tau_overall_samples[store_idx] <- current_tau_overall
        if (family == "gaussian") sigma_sq_samples[store_idx] <- sigma_sq
      }
    } # End MCMC loop for the current chain
    
    # Store the results of this chain in the master list
    chain_output <- list(
      alpha = alpha_global_samples,
      aleph = aleph_samples,
      Beta_hs = Beta_hs_samples,
      lambda_hs = lambda_hs_samples,
      tau_prognostic = tau_prognostic_samples,
      tau_treatment = tau_treatment_samples,
      tau_overall = tau_overall_samples
    )
    if (family == "gaussian") {
      chain_output$sigma_sq <- sigma_sq_samples
    }
    all_chains_results[[chain_num]] <- chain_output
    
  } # End outer loop for chains
  
  # --- Combine Samples from All Chains ---
  if(verbose) message("Combining samples from all chains...")
  
  combined_alpha_samples <- do.call(c, lapply(all_chains_results, `[[`, "alpha"))
  combined_aleph_samples <- do.call(c, lapply(all_chains_results, `[[`, "aleph"))
  combined_tau_prognostic <- do.call(c, lapply(all_chains_results, `[[`, "tau_prognostic"))
  combined_tau_treatment <- do.call(c, lapply(all_chains_results, `[[`, "tau_treatment"))
  combined_tau_overall <- do.call(c, lapply(all_chains_results, `[[`, "tau_overall"))
  
  if(p_hs_total > 0) {
    combined_Beta_hs_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "Beta_hs"))
    combined_lambda_hs_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_hs"))
  } else {
    combined_Beta_hs_samples <- matrix(0,0,0)
    combined_lambda_hs_samples <- matrix(0,0,0)
  }
  
  if (family == "gaussian") {
    combined_sigma_sq_samples <- do.call(c, lapply(all_chains_results, `[[`, "sigma_sq"))
  }
  
  # --- Unpack Combined Samples into Named Components ---
  # Define indices based on the structure of X_hs
  num_beta_main <- if(!is.null(X_prog_list$beta)) ncol(X_prog_list$beta) else 0
  num_beta_prop <- if(!is.null(X_prog_list$beta_propensity)) 1 else 0 # It's one column if it exists
  num_beta_int  <- if(!is.null(X_prog_list$beta_int)) ncol(X_prog_list$beta_int) else 0
  
  idx_current <- 0
  idx_beta <- if(num_beta_main > 0) (idx_current + 1):(idx_current + num_beta_main) else integer(0); idx_current <- idx_current + num_beta_main
  idx_beta_prop <- if(num_beta_prop > 0) (idx_current + 1):(idx_current + num_beta_prop) else integer(0); idx_current <- idx_current + num_beta_prop
  idx_beta_int <- if(num_beta_int > 0) (idx_current + 1):(idx_current + num_beta_int) else integer(0); idx_current <- idx_current + num_beta_int
  
  num_gamma_main <- if(!is.null(X_treat_list$gamma)) ncol(X_treat_list$gamma) else 0
  num_gamma_int  <- if(!is.null(X_treat_list$gamma_int)) ncol(X_treat_list$gamma_int) else 0
  
  idx_gamma <- if(num_gamma_main > 0) (idx_current + 1):(idx_current + num_gamma_main) else integer(0); idx_current <- idx_current + num_gamma_main
  idx_gamma_int <- if(num_gamma_int > 0) (idx_current + 1):(idx_current + num_gamma_int) else integer(0)
  
  # Unpack Beta_hs_samples
  beta_out <- if(length(idx_beta) > 0) combined_Beta_hs_samples[idx_beta, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  beta_prop_out <- if(length(idx_beta_prop) > 0) combined_Beta_hs_samples[idx_beta_prop, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  beta_int_out <- if(length(idx_beta_int) > 0) combined_Beta_hs_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_out <- if(length(idx_gamma) > 0) combined_Beta_hs_samples[idx_gamma, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_int_out <- if(length(idx_gamma_int) > 0) combined_Beta_hs_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  
  # Unpack lambda_hs_samples using the same indices
  lambda_beta_out <- if(length(idx_beta) > 0) combined_lambda_hs_samples[idx_beta, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  lambda_beta_prop_out <- if(length(idx_beta_prop) > 0) combined_lambda_hs_samples[idx_beta_prop, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  lambda_beta_int_out <- if(length(idx_beta_int) > 0) combined_lambda_hs_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  lambda_gamma_out <- if(length(idx_gamma) > 0) combined_lambda_hs_samples[idx_gamma, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  lambda_gamma_int_out <- if(length(idx_gamma_int) > 0) combined_lambda_hs_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  
  # Set names for output matrices
  if(nrow(beta_out) > 0) rownames(beta_out) <- colnames(X_prog_list$beta)
  if(nrow(beta_prop_out) > 0) rownames(beta_prop_out) <- "propensity_score"
  if(nrow(beta_int_out) > 0) rownames(beta_int_out) <- colnames(X_prog_list$beta_int)
  if(nrow(gamma_out) > 0) rownames(gamma_out) <- colnames(X_treat_list$gamma)
  if(nrow(gamma_int_out) > 0) rownames(gamma_int_out) <- colnames(X_treat_list$gamma_int)
  
  output <- list(
    alpha = combined_alpha_samples,
    beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
    beta_propensity = if(nrow(beta_prop_out) > 0) t(beta_prop_out) else matrix(NA, eff_samp_count * num_chains, 0),
    beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    aleph = combined_aleph_samples,
    gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
    gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    
    # Adding the unpacked lambda samples to the output list
    lambda_beta = if(nrow(lambda_beta_out) > 0) t(lambda_beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
    lambda_beta_propensity = if(nrow(lambda_beta_prop_out) > 0) t(lambda_beta_prop_out) else matrix(NA, eff_samp_count * num_chains, 0),
    lambda_beta_int = if(nrow(lambda_beta_int_out) > 0) t(lambda_beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    lambda_gamma = if(nrow(lambda_gamma_out) > 0) t(lambda_gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
    lambda_gamma_int = if(nrow(lambda_gamma_int_out) > 0) t(lambda_gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    
    tau_prognostic = combined_tau_prognostic,
    tau_treatment = combined_tau_treatment,
    tau_overall = combined_tau_overall
  )
  
  if (family == "gaussian") {
    output$sigma_sq = combined_sigma_sq_samples
  }
  
  return(output)
}

######################## GAUSSIAN FAMILY. ########################

# --- Example Usage (using the simple simulation) ---
source('R/simul_1.R')

data <- generate_data_2(500, is_te_hetero = F, is_mu_nonlinear = F, seed = 40, RCT = FALSE, z_diff = F, contrast_binary = T)
X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
y <- as.numeric(data$y)
z <- as.numeric(data$z)

cat("\n--- Running Grouped Horseshoe Logistic Sampler ---\n")
start.time <- Sys.time()
fit_grouped_hs <- fit_grouped_horseshoes_R(
  y_vec = y,
  X_mat = X,
  Z_vec = z,
  family = "gaussian",
  n_iter = 4000, # Increase for better results
  burn_in = 1000,
  num_chains = 3,
  propensity_as_covariate = T,
  method_tau_prognostic = "halfCauchy", tau_prognostic_init = 0.1,
  method_tau_treatment = "halfCauchy", tau_treatment_init = 0.1,
  method_tau_overall = "fixed", tau_overall_init = 1,
  alpha_global_prior_sd = 5.0,
  aleph_prior_sd = 5.0,
  thin = 1,
  seed = 103,
  verbose = TRUE,
  ping = 500
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
cat("\n--- Grouped HS: Posterior Means vs True Parameters ---\n")
cat("Alpha_global (Overall Intercept): True =", true_params_check$alpha,
    ", Estimated =", mean(fit_grouped_hs$alpha), "\n")

if(ncol(fit_grouped_hs$beta) > 0) {
  beta_comp_ghs <- data.frame(True = true_params_check$beta, Estimated = colMeans(fit_grouped_hs$beta))
  rownames(beta_comp_ghs) <- paste0("Beta_", colnames(X_sim))
  print(beta_comp_ghs)
}

# For beta_interaction, true is 0 in simple DGP
if (ncol(fit_grouped_hs$beta_interaction) > 0) {
  cat("\nBeta_interaction (Prognostic XX - True should be ~0 for simple DGP):\n")
  print(colMeans(fit_grouped_hs$beta_interaction))
}

cat("\nAleph (Treatment Intercept Modifier): True =", true_params_check$aleph,
    ", Estimated =", mean(fit_grouped_hs$aleph), "\n")

if(ncol(fit_grouped_hs$gamma) > 0) {
  gamma_comp_ghs <- data.frame(True = true_params_check$gamma, Estimated = colMeans(fit_grouped_hs$gamma))
  rownames(gamma_comp_ghs) <- paste0("Gamma_Z_", colnames(X_sim))
  print(gamma_comp_ghs)
}

# For gamma_int, true is 0 in simple DGP
if (ncol(fit_grouped_hs$gamma_int) > 0) {
  cat("\nGamma_int (Z*XX - True should be ~0 for simple DGP")
  }