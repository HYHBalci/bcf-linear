library(BayesLogit)
library(MASS)
library(stochtree)
library(Rcpp)
library(ggplot2)
sourceCpp("C:/Users/P094412/OneDrive - Amsterdam UMC/Documenten/bcf-linear/src/horseshoe_samplers.cpp")

interaction_pairs <- function(num_covariates, boolean_vector) {
  interaction_list <- list()
  if (num_covariates > 1) {
    for (j in 1:(num_covariates - 1)) {
      for (k in (j + 1):num_covariates) {
        if (boolean_vector[j] || boolean_vector[k])
          interaction_list[[length(interaction_list) + 1]] <- c(j, k)
      }
    }
  }
  return(do.call(cbind, interaction_list))
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
#' @param propensity_train propensity train.
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
    propensity_train = NULL,
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
    ping = 1000,
    standardize_cov = F,
    interaction_rule = c("continuous", "continuous_or_binary", "all"),
    cat_coding_method = "difference"
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
    if(is.null(propensity_train)){
    if (verbose) message("Estimating propensity scores to use as a covariate...")
    ps_num_burnin <- 10
    ps_num_total <- 50
    bart_model_propensity <- bart(X_train = X_mat_orig, y_train = as.numeric(Z_vec), X_test = NULL, 
                                  num_gfr = ps_num_total, num_burnin = 0, num_mcmc = 0)
    propensity_scores <- rowMeans(bart_model_propensity$y_hat_train[,(ps_num_burnin+1):ps_num_total])
    if (verbose) message("Propensity scores calculated.")
    } else{
      propensity_scores <- propensity_train
    }
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
  interaction_rule <- match.arg(interaction_rule)
  handled_data_list <- standardize_X_by_index(X_mat, process_data = standardize_cov, interaction_rule = interaction_rule, cat_coding_method = cat_coding_method)
  X_final_var_info <- handled_data_list$X_final_var_info
  X_mat_orig <- handled_data_list$X_final
  if(interaction_rule == 'continuous'){
    boolean_continuous <- as.vector(X_final_var_info$is_continuous)
  } else if(interaction_rule == 'continuous_or_binary'){
    boolean_continuous <- as.vector(X_final_var_info$is_continuous) + as.vector(X_final_var_info$is_binary)
  } else{
    boolean_continuous <- as.vector(X_final_var_info$is_continuous) + as.vector(X_final_var_info$is_binary) + as.vector(X_final_var_info$is_categorical)
  }
  boolean_continous <- as.logical(boolean_continuous)
  print(boolean_continous)
  
  
  main_interaction_indices_X <- interaction_pairs(p_main_X_orig, boolean_continous)
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
        omega <- BayesLogit::rpg(num = N, h = 1, z = eta)
        omega[omega < 1e-9] <- 1e-9; z_aug <- y_adj / omega
        w_vec <- omega; scale_factor <- 1.0
      } else { z_aug <- y_vec; w_vec <- rep(1, N)
      scale_factor <- sigma_sq }
      resid_for_alpha <- z_aug - (Z_vec * aleph + eta_hs_part); prior_prec_alpha <- 1 / safe_var_R(alpha_global_prior_sd^2); post_prec_alpha <- sum(w_vec) / scale_factor + prior_prec_alpha; post_mean_alpha <- (sum(w_vec * resid_for_alpha) / scale_factor + alpha_global_prior_mean * prior_prec_alpha) / post_prec_alpha; alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
      resid_for_aleph <- z_aug - (alpha_global + eta_hs_part); prior_prec_aleph <- 1 / safe_var_R(aleph_prior_sd^2); post_prec_aleph <- sum(w_vec * Z_vec^2) / scale_factor + prior_prec_aleph; post_mean_aleph <- (sum(w_vec * Z_vec * resid_for_aleph) / scale_factor + aleph_prior_mean * prior_prec_aleph) / post_prec_aleph; aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
      if (p_hs_total > 0) {
        z_for_hs <- z_aug - alpha_global - Z_vec * aleph;
        idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0);
        idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_hs_total else integer(0)
        effective_global_taus_for_beta <- rep(1, p_hs_total); if(length(idx_prog) > 0) effective_global_taus_for_beta[idx_prog] <- current_tau_prognostic * current_tau_overall; if(length(idx_treat) > 0) effective_global_taus_for_beta[idx_treat] <- current_tau_treatment * current_tau_overall
        temp_lambda_scale_term <- (Beta_hs^2) / (2 * scale_factor * safe_var_R(effective_global_taus_for_beta^2)); nu_lambda_inv <- stats::rexp(p_hs_total, rate = 1 + (1 / safe_var_R(lambda_hs^2))); lambda_hs_sq_inv <- stats::rexp(p_hs_total, rate = nu_lambda_inv + temp_lambda_scale_term); lambda_hs <- 1 / sqrt(safe_var_R(lambda_hs_sq_inv))
        if (length(idx_prog) > 0 && method_tau_prognostic != "fixed") { 
          xi_inv_prog <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_prognostic^2))
          rate_prog <- xi_inv_prog + sum(Beta_hs[idx_prog]^2 / (2 * scale_factor * safe_var_R(lambda_hs[idx_prog]^2 * current_tau_overall^2)))
          current_tau_prognostic <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog)) 
          }
        if (length(idx_treat) > 0 && method_tau_treatment != "fixed") { 
          xi_inv_treat <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_treatment^2))
          rate_treat <- xi_inv_treat + sum(Beta_hs[idx_treat]^2 / (2 * scale_factor * safe_var_R(lambda_hs[idx_treat]^2 * current_tau_overall^2)))
          current_tau_treatment <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat))
          }
        tau_vector_components <- rep(1, p_hs_total)
        if(length(idx_prog) > 0){
          tau_vector_components[idx_prog] <- current_tau_prognostic
        } 
        if(length(idx_treat) > 0) {
          tau_vector_components[idx_treat] <- current_tau_treatment
        }
        if (method_tau_overall != "fixed") {
          xi_inv_overall <- stats::rexp(1, rate = 1 + 1/safe_var_R(current_tau_overall^2))
          rate_overall <- xi_inv_overall + sum(Beta_hs^2 / (2 * scale_factor * safe_var_R(lambda_hs^2 * tau_vector_components^2)))
          current_tau_overall <- sqrt(rinvgamma_R(1, shape = (p_hs_total + 1)/2, scale = rate_overall)) 
          }
        current_tau_prognostic <- max(current_tau_prognostic, 1e-6)
        current_tau_treatment <- max(current_tau_treatment, 1e-6)
        current_tau_overall <- max(current_tau_overall, 1e-6)
        lambda_hs[lambda_hs < 1e-6] <- 1e-6
        prior_var_beta_hs <- safe_var_R(lambda_hs^2 * effective_global_taus_for_beta^2)
        prior_prec_beta_hs_diag <- diag(1 / prior_var_beta_hs, p_hs_total, p_hs_total)
        XtWX_hs <- t(X_hs) %*% (w_vec * X_hs)
        post_prec_beta_hs_unscaled <- XtWX_hs + prior_prec_beta_hs_diag
        chol_attempt <- try(chol(post_prec_beta_hs_unscaled), silent = TRUE)
        if (inherits(chol_attempt, "try-error")) {
          warning(paste("Cholesky for Beta_hs failed at iter", iter_idx, "- skipping Beta_hs update."))
        } else {
        post_cov_beta_hs_unscaled_inv <- chol2inv(chol_attempt)
        rhs_mean <- t(X_hs) %*% (w_vec * z_for_hs)
        post_mean_beta_hs <- post_cov_beta_hs_unscaled_inv %*% rhs_mean
        Beta_hs <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean_beta_hs), Sigma = scale_factor * post_cov_beta_hs_unscaled_inv)
        }
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
#' @title Fit BCF-style Linear Model with Flexible Shrinkage Options
#' @description R wrapper to run a Gibbs sampler for a complex logistic or linear model.
#'
#' @param y_vec Response vector. Binary (0s and 1s) for binomial, continuous for gaussian.
#' @param X_mat Covariate matrix.
#' @param Z_vec Binary treatment indicator vector (0s and 1s).
#' @param family Model family, either "binomial" (default) or "gaussian".
#' @param n_iter Total MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param prognostic_shrinkage Shrinkage method for mu(x) part: "horseshoe", "linked_shrinkage", or "none".
#' @param treatment_shrinkage Shrinkage method for tau(x) part: "horseshoe", "linked_shrinkage", or "none".
#' @param no_shrinkage_prior_variance Constant for the prior variance if shrinkage is "none". Default is 1e4.
#' @param standardize_cov A boolean flag. If `TRUE` (default), the function standardizes numeric
#'  variables and creates dummy variables for categorical ones. If `FALSE`, it bypasses
#'  transformations and returns the original data as a matrix, but still generates the metadata.
#' @param interaction_rule A character string specifying which interactions are allowed.
#' Must be one of:
#'  \itemize{
#'  \item `"continuous"` (Default): Allows interactions where at least one variable is continuous.
#'  \item `"continuous_or_binary"`: Allows interactions where at least one variable is continuous OR binary.
#'  \item `"all"`: Allows all possible two-way interactions between valid columns.
#'}
#' @param cat_coding_method A character string for categorical variable contrast coding.
#' Must be one of `"sum"` (for sum-to-zero/deviation coding) or `"difference"` (for successive differences).
#' Default is `"difference"`.
#' @param num_chains Integer, the number of independent MCMC chains to run.
#' @param propensity_as_covariate Boolean, if TRUE, propensity scores are estimated
#'   and added as a covariate to the prognostic model `mu(x)`.
#' @param thin Thinning parameter.
#' @param seed Random seed.
#' @param verbose Logical, print progress.
#' @param ping Frequency of progress messages.
#'
#' @return A list of combined posterior samples.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma rexp runif pgamma qgamma rnorm var combn
linear_linear_function <- function(
    y_vec, X_mat, Z_vec,
    family = c("binomial", "gaussian"),
    n_iter, burn_in,
    prognostic_shrinkage = c("horseshoe", "linked_shrinkage", "none"),
    treatment_shrinkage = c("horseshoe", "linked_shrinkage", "none"),
    no_shrinkage_prior_variance = 1e4,
    standardize_cov = F,
    interaction_rule = c("continuous", "continuous_or_binary", "all"),
    cat_coding_method = "difference",
    num_chains = 1, propensity_as_covariate = FALSE,
    thin = 1, seed = 123, verbose = FALSE, ping = 1000,
    alpha_global_prior_sd = 10.0, aleph_prior_sd = 10.0
) {
  # --- Setup and Validation ---
  if (!is.null(seed)) { set.seed(seed)
    chain_seeds <- sample.int(1e6, num_chains) } else {
      chain_seeds <- rep(NA, num_chains)
      }
  family <- match.arg(family)
  prognostic_shrinkage <- match.arg(prognostic_shrinkage)
  treatment_shrinkage <- match.arg(treatment_shrinkage)
  interaction_rule <- match.arg(interaction_rule)
  handled_data_list <- standardize_X_by_index(X_mat, process_data = standardize_cov, interaction_rule = interaction_rule, cat_coding_method = cat_coding_method)
  X_final_var_info <- handled_data_list$X_final_var_info
  X_mat <- handled_data_list$X_final
  if(interaction_rule == 'continuous'){
    boolean_continuous <- as.vector(X_final_var_info$is_continuous)
  } else if(interaction_rule == 'continuous_or_binary'){
    boolean_continuous <- as.vector(X_final_var_info$is_continuous) + as.vector(X_final_var_info$is_binary)
  } else{
    boolean_continuous <- as.vector(X_final_var_info$is_continuous) + as.vector(X_final_var_info$is_binary) + as.vector(X_final_var_info$is_categorical)
  }
  boolean_continous <- as.logical(boolean_continuous)
  print(boolean_continous)

  if (length(unique(na.omit(y_vec))) <= 1) stop("Outcome variable `y_vec` is constant.")


  N <- nrow(X_mat)
  p_main_X_orig <- ncol(X_mat)
  X_mat_orig <- X_mat

  # --- Propensity Score Estimation ---
  if (propensity_as_covariate) {
    if (verbose) message("Estimating propensity scores to use as a covariate...")
    if (!requireNamespace("stochtree", quietly = TRUE) || !exists("bart", where = "package:stochtree")) {
      stop("'stochtree::bart' function not found for propensity estimation. Please install and load the 'stochtree' package.")
    }
    bart_model_propensity <- stochtree::bart(X_train = X_mat_orig, y_train = as.numeric(Z_vec), num_gfr = 50, num_mcmc = 0)
    propensity_scores <- rowMeans(bart_model_propensity$y_hat_train)
    X_mat_for_prog <- cbind(X_mat_orig, propensity_score = propensity_scores)
    if (verbose) message("Propensity scores added to X_mat for prognostic model.")
  } else {
    X_mat_for_prog <- X_mat_orig
  }
  p_main_X_prog <- ncol(X_mat_for_prog)

  # --- 1. Construct Design Matrices (once, outside the chain loop) ---
  if(is.null(colnames(X_mat_for_prog))) colnames(X_mat_for_prog) <- paste0("X_prog_", 1:p_main_X_prog)
  if(is.null(colnames(X_mat_orig))) colnames(X_mat_orig) <- paste0("X_orig_", 1:p_main_X_orig)

  # Prognostic Part (beta, beta_int) - Based on X_mat_for_prog
  main_interaction_indices_prog <- interaction_pairs(p_main_X_orig, boolean_continous)
  p_interaction_prog <- if(is.matrix(main_interaction_indices_prog)) ncol(main_interaction_indices_prog) else 0
  X_prognostic <- X_mat_for_prog
  if (p_interaction_prog > 0) {

    # FIX 2: Loop correctly and use column indexing [row, k]
    X_beta_int_mat <- do.call(cbind, lapply(1:p_interaction_prog, function(k) {
      # Get the two covariate indices for the k-th interaction pair
      idx1 <- main_interaction_indices_prog[1, k]
      idx2 <- main_interaction_indices_prog[2, k]

      # Create the interaction term using the prognostic matrix
      (X_mat_for_prog[, idx1] * X_mat_for_prog[, idx2])
    }))

    colnames(X_beta_int_mat) <- sapply(1:p_interaction_prog, function(k) {
      idx1 <- main_interaction_indices_prog[1, k]
      idx2 <- main_interaction_indices_prog[2, k]
      paste0("beta_int_", paste(colnames(X_mat_for_prog)[c(idx1, idx2)], collapse = "_"))
    })

    X_prognostic <- cbind(X_prognostic, X_beta_int_mat)
  }
  p_prognostic <- ncol(X_prognostic)

  X_treatment_related <- X_mat_orig * Z_vec
  colnames(X_treatment_related) <- paste0("Z_", colnames(X_mat_orig))
  main_interaction_indices_orig <- interaction_pairs(p_main_X_orig, boolean_continous)

  p_interaction_orig <- if (is.matrix(main_interaction_indices_orig)) ncol(main_interaction_indices_orig) else 0

  X_treatment_related <- X_mat_orig * Z_vec
  colnames(X_treatment_related) <- paste0("Z_", colnames(X_mat_orig))

  if (p_interaction_orig > 0) {

    base_int_treat <- do.call(cbind, lapply(1:p_interaction_orig, function(k) {
      idx1 <- main_interaction_indices_orig[1, k]
      idx2 <- main_interaction_indices_orig[2, k]

      (X_mat_orig[, idx1] * X_mat_orig[, idx2])
    }))

    X_gamma_int_mat <- base_int_treat * Z_vec

    # Create the corresponding column names
    colnames(X_gamma_int_mat) <- sapply(1:p_interaction_orig, function(k) {
      idx1 <- main_interaction_indices_orig[1, k]
      idx2 <- main_interaction_indices_orig[2, k]
      paste0("Z_int_", paste(colnames(X_mat_orig)[c(idx1, idx2)], collapse = "_"))
    })

    X_treatment_related <- cbind(X_treatment_related, X_gamma_int_mat)
  }

  p_treatment_related <- ncol(X_treatment_related)
  X_full <- cbind(X_prognostic, X_treatment_related)
  p_total <- ncol(X_full)

  # --- MCMC Setup ---
  all_chains_results <- list()
  for (chain_num in 1:num_chains) {
    if(verbose && num_chains > 1) message(sprintf("\n--- Starting Chain %d / %d ---\n", chain_num, num_chains))
    if (!is.na(chain_seeds[chain_num])) set.seed(chain_seeds[chain_num])

    # Initialize parameters for this chain
    alpha_global <- 0.0; aleph <- 0.0
    Beta_all <- rep(0, p_total)

    # Prognostic shrinkage params
    lambda_prog <- rep(1, p_prognostic); tau_prog <- rep(1, p_main_X_prog)
    tau_int_prog <- 0.5; tau_glob_prog <- 1.0

    # Treatment shrinkage params
    lambda_treat <- rep(1, p_treatment_related); tau_treat <- rep(1, p_main_X_orig)
    tau_int_treat <- 0.5; tau_glob_treat <- 1.0

    if (family == "gaussian") { sigma_sq <- var(y_vec, na.rm=TRUE); if(is.na(sigma_sq) || sigma_sq == 0) sigma_sq <- 1.0 } else { y_adj <- y_vec - 0.5 }

    # Storage for this chain
    eff_samp_count <- floor((n_iter - burn_in) / thin)
    alpha_global_samples <- numeric(eff_samp_count); aleph_samples <- numeric(eff_samp_count)
    Beta_all_samples <- if(p_total > 0) matrix(0, nrow = p_total, ncol = eff_samp_count) else matrix(0,0,0)

    tau_prog_samples <- if(prognostic_shrinkage=="linked_shrinkage") matrix(0, nrow=p_main_X_prog, ncol=eff_samp_count) else numeric(eff_samp_count)
    tau_treat_samples <- if(treatment_shrinkage=="linked_shrinkage") matrix(0, nrow=p_main_X_orig, ncol=eff_samp_count) else numeric(eff_samp_count)
    tau_int_prog_samples <- numeric(eff_samp_count); tau_int_treat_samples <- numeric(eff_samp_count)
    tau_glob_prog_samples <- numeric(eff_samp_count); tau_glob_treat_samples <- numeric(eff_samp_count)
    lambda_prog_samples <- if(p_prognostic > 0) matrix(0, nrow=p_prognostic, ncol=eff_samp_count) else matrix(0,0,0)
    lambda_treat_samples <- if(p_treatment_related > 0) matrix(0, nrow=p_treatment_related, ncol=eff_samp_count) else matrix(0,0,0)

    if (family == "gaussian") sigma_sq_samples <- numeric(eff_samp_count)
    store_idx <- 0

    # MCMC Loop
    if(verbose) message(sprintf("MCMC Chain %d (%s model) is running...", chain_num, family))
    for (iter_idx in 1:n_iter) {
      if (verbose && iter_idx %% ping == 0) message("Chain ", chain_num, " - Iteration: ", iter_idx, "/", n_iter)
      
      eta <- alpha_global + Z_vec * aleph + as.vector(X_full %*% Beta_all)
      if (family == "binomial") {
        omega <- BayesLogit::rpg(N, 1, eta); omega[omega < 1e-9] <- 1e-9
        z_aug <- y_adj / omega; w_vec <- omega; scale_factor <- 1.0
      } else { z_aug <- y_vec; w_vec <- rep(1, N); scale_factor <- sigma_sq }
      
      eta_hs_part <- if(p_total > 0) as.vector(X_full %*% Beta_all) else rep(0, N)
      
      # --- Update alpha_global ---
      resid_for_alpha <- z_aug - (Z_vec * aleph + eta_hs_part)
      prior_prec_alpha <- 1 / safe_var_R(alpha_global_prior_sd^2)
      post_prec_alpha <- sum(w_vec) / scale_factor + prior_prec_alpha
      # CORRECTED: Added the prior mean term to the numerator
      post_mean_alpha <- ( (sum(w_vec * resid_for_alpha) / scale_factor) + (0.0 * prior_prec_alpha) ) / post_prec_alpha
      alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
      
      # --- Update aleph ---
      resid_for_aleph <- z_aug - (alpha_global + eta_hs_part)
      prior_prec_aleph <- 1 / safe_var_R(aleph_prior_sd^2)
      post_prec_aleph <- sum(w_vec * Z_vec^2) / scale_factor + prior_prec_aleph
      # CORRECTED: Added the prior mean term to the numerator
      post_mean_aleph <- ( (sum(w_vec * Z_vec * resid_for_aleph) / scale_factor) + (0.0 * prior_prec_aleph) ) / post_prec_aleph
      aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
      
      if (p_total > 0) {
        z_for_betas <- z_aug - alpha_global - Z_vec * aleph
        
        idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0)
        idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_total else integer(0)
        
        # --- Update Prognostic Shrinkage (if not "none") ---
        if (prognostic_shrinkage == "horseshoe") {
          if(length(idx_prog) > 0) {
            beta_prog <- Beta_all[idx_prog]
            temp_lambda_scale_term <- (beta_prog^2) / (2 * scale_factor * safe_var_R(tau_glob_prog^2))
            nu_lambda_inv <- rexp(length(idx_prog), rate = 1 + (1 / safe_var_R(lambda_prog^2)))
            lambda_prog_sq_inv <- rexp(length(idx_prog), rate = nu_lambda_inv + temp_lambda_scale_term)
            lambda_prog <- 1 / sqrt(safe_var_R(lambda_prog_sq_inv))
            xi_inv_prog <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_prog^2))
            rate_prog <- xi_inv_prog + sum(beta_prog^2 / (2 * scale_factor * safe_var_R(lambda_prog^2)))
            tau_glob_prog <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog))
          }
        } else if (prognostic_shrinkage == "linked_shrinkage") {
          if(p_prognostic > 0) {
            beta_prog_main <- Beta_all[1:p_main_X_prog]
            beta_prog_int <- Beta_all[(p_main_X_prog + 1):p_prognostic]
            for(j in 1:p_main_X_prog) {
              tau_prog[j] <- sample_tau_j_slice(tau_prog[j], beta_prog_main[j], j, beta_int = beta_prog_int, tau = tau_treat, tau_int = 1, sigma = sqrt(sigma_sq))
              if(p_interaction_prog > 0){
                tau_int_prog <-1
              }
            }
          }
        }
          # --- Update Treatment Shrinkage (if not "none") ---
          if (treatment_shrinkage == "horseshoe") {
            if(length(idx_treat) > 0) {
              beta_treat <- Beta_all[idx_treat]
              temp_lambda_scale_term <- (beta_treat^2) / (2 * scale_factor * safe_var_R(tau_glob_treat^2))
              nu_lambda_inv <- rexp(length(idx_treat), rate = 1 + (1 / safe_var_R(lambda_treat^2)))
              lambda_treat_sq_inv <- rexp(length(idx_treat), rate = nu_lambda_inv + temp_lambda_scale_term)
              lambda_treat <- 1 / sqrt(safe_var_R(lambda_treat_sq_inv))
              xi_inv_treat <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_treat^2))
              rate_treat <- xi_inv_treat + sum(beta_treat^2 / (2 * scale_factor * safe_var_R(lambda_treat^2))); tau_glob_treat <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat))
            }
          } else if (treatment_shrinkage == "linked_shrinkage") {
            beta_treat_main <- Beta_all[idx_treat[1:p_main_X_orig]]
            beta_treat_int <- Beta_all[idx_treat[(p_main_X_orig + 1):p_treatment_related]]
            for(j in 1:p_main_X_orig) {
              tau_treat[j] <- sample_tau_j_slice(tau_treat[j], beta_treat_main[j], j, beta_int = beta_treat_int, tau = tau_treat, tau_int = 1, sigma = sqrt(sigma_sq))
            }
            if(p_interaction_orig > 0){ tau_int_treat <- 1 }
          }
          
          # --- Construct Prior Precision and Sample Beta_all ---
          prior_var_beta <- numeric(p_total)
          # Prognostic Prior
          if (prognostic_shrinkage == "none") {
            if(length(idx_prog) > 0) prior_var_beta[idx_prog] <- no_shrinkage_prior_variance
          } else if (prognostic_shrinkage == "horseshoe") {
            if(length(idx_prog) > 0) prior_var_beta[idx_prog] <- lambda_prog^2 * tau_glob_prog^2
          } else { # linked_shrinkage
            if(p_prognostic > 0) {
              prog_main_indices <- 1:p_main_X_prog
              prior_var_beta[idx_prog[prog_main_indices]] <- tau_prog^2
              if(p_interaction_prog > 0){
                prog_int_indices <- (p_main_X_prog+1):p_prognostic
                for(k in 1:p_interaction_prog) {
                  pair <- main_interaction_indices_prog[,k]
                  prior_var_beta[idx_prog[prog_int_indices[k]]] <- tau_int_prog^2 * tau_prog[pair[1]]^2 * tau_prog[pair[2]]^2
                }
              }
            }
          }
          # Treatment Prior
          if (treatment_shrinkage == "none") {
            if(length(idx_treat) > 0) prior_var_beta[idx_treat] <- no_shrinkage_prior_variance
          } else if (treatment_shrinkage == "horseshoe") {
            if(length(idx_treat) > 0) prior_var_beta[idx_treat] <- lambda_treat^2 * tau_glob_treat^2
          } else { # linked_shrinkage
            if(p_treatment_related > 0) {
              treat_main_indices <- 1:p_main_X_orig
              prior_var_beta[idx_treat[treat_main_indices]] <- tau_treat^2
              if(p_interaction_orig > 0){
                treat_int_indices <- (p_main_X_orig+1):p_treatment_related
                for(k in 1:p_interaction_orig) {
                  pair <- main_interaction_indices_orig[,k]
                  prior_var_beta[idx_treat[treat_int_indices[k]]] <- tau_int_treat^2 * tau_treat[pair[1]]^2 * tau_treat[pair[2]]^2
                }
              }
            }
          }
          
          
          prior_prec_diag_mat <- diag(1 / safe_var_R(prior_var_beta), p_total, p_total)
          post_prec_matrix <- (t(X_full) %*% (w_vec * X_full)) / scale_factor + prior_prec_diag_mat
          
          chol_attempt <- try(chol(post_prec_matrix), silent = TRUE)
          if (inherits(chol_attempt, "try-error")) {
            warning(paste("Cholesky for Beta_all failed at iter", iter_idx, "- skipping update."))
          } else {
            post_cov_matrix <- chol2inv(chol_attempt)
            
            rhs_mean_scaled <- (t(X_full) %*% (w_vec * z_for_betas)) / scale_factor
            
            post_mean <- post_cov_matrix %*% rhs_mean_scaled
            
            Beta_all <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean), Sigma = post_cov_matrix)
          }
        }
        
        # --- Sample sigma_sq for Gaussian model ---
        if (family == "gaussian") {
          eta_new <- alpha_global + Z_vec * aleph + as.vector(X_full %*% Beta_all)
          residuals_for_sigma <- y_vec - eta_new; shape_sigma_post <- 0.001 + N / 2
          rate_sigma_post <- 0.001 + 0.5 * sum(residuals_for_sigma^2)
          sigma_sq <- rinvgamma_R(1, shape = shape_sigma_post, scale = rate_sigma_post)
          sigma_sq <- max(sigma_sq, 1e-9)
        }

      # --- Store samples ---
      if (iter_idx > burn_in && (iter_idx - burn_in) %% thin == 0) {
        store_idx <- store_idx + 1
        alpha_global_samples[store_idx] <- alpha_global; aleph_samples[store_idx] <- aleph
        if(p_total > 0) Beta_all_samples[, store_idx] <- Beta_all
        if (family == "gaussian") sigma_sq_samples[store_idx] <- sigma_sq
        if(prognostic_shrinkage == "horseshoe") { tau_prog_samples[store_idx] <- tau_glob_prog; if(p_prognostic > 0) lambda_prog_samples[, store_idx] <- lambda_prog }
        else if(prognostic_shrinkage == "linked_shrinkage") { if(p_main_X_prog > 0) tau_prog_samples[, store_idx] <- tau_prog; tau_int_prog_samples[store_idx] <- tau_int_prog }
        if (treatment_shrinkage == "horseshoe") { tau_treat_samples[store_idx] <- tau_glob_treat; if(p_treatment_related > 0) lambda_treat_samples[, store_idx] <- lambda_treat }
        else if(treatment_shrinkage == "linked_shrinkage") { if(p_main_X_orig > 0) tau_treat_samples[, store_idx] <- tau_treat; tau_int_treat_samples[store_idx] <- tau_int_treat }
      }
    } # End MCMC loop for the current chain

    # Store the results of this chain
    all_chains_results[[chain_num]] <- list(alpha=alpha_global_samples, aleph=aleph_samples, Beta_all=Beta_all_samples, sigma_sq=if(family=="gaussian") sigma_sq_samples else NULL,
                                            tau_prog=tau_prog_samples, tau_int_prog=tau_int_prog_samples, lambda_prog=lambda_prog_samples,
                                            tau_treat=tau_treat_samples, tau_int_treat=tau_int_treat_samples, lambda_treat=lambda_treat_samples)
  } # End outer loop for chains

  # --- Combine and Unpack Results ---
  if(verbose) message("Combining samples from all chains...")

  combined_alpha_samples <- do.call(c, lapply(all_chains_results, `[[`, "alpha"))
  combined_aleph_samples <- do.call(c, lapply(all_chains_results, `[[`, "aleph"))
  if(p_total > 0) combined_Beta_all_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "Beta_all")) else combined_Beta_all_samples <- matrix(0,0,0)

  # Combine shrinkage parameters
  if(prognostic_shrinkage == "horseshoe") {
    combined_tau_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_prog"))
    if(p_prognostic>0) combined_lambda_prog <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_prog"))
  } else if (prognostic_shrinkage == "linked_shrinkage") {
    if(p_main_X_prog > 0) combined_tau_prog <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_prog))) else combined_tau_prog <- matrix(0, eff_samp_count * num_chains, 0)
    combined_tau_int_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_prog"))
  }
  if(treatment_shrinkage == "horseshoe") {
    combined_tau_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_treat"))
    if(p_treatment_related>0) combined_lambda_treat <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_treat"))
  } else if (treatment_shrinkage == "linked_shrinkage") {
    if(p_main_X_orig > 0) combined_tau_treat <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_treat))) else combined_tau_treat <- matrix(0, eff_samp_count * num_chains, 0)
    combined_tau_int_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_treat"))
  }

  if (family == "gaussian") { combined_sigma_sq_samples <- do.call(c, lapply(all_chains_results, `[[`, "sigma_sq")) }

  # Unpack Beta_all_samples
  num_beta_main_prog <- ncol(X_mat_for_prog)
  num_beta_int  <- p_interaction_prog
  num_gamma_main <- ncol(X_mat_orig)
  num_gamma_int  <- p_interaction_orig

  idx_current <- 0
  idx_beta <- (idx_current + 1):(idx_current + num_beta_main_prog); idx_current <- idx_current + num_beta_main_prog
  idx_beta_int <- if(num_beta_int > 0) (idx_current + 1):(idx_current + num_beta_int) else integer(0); idx_current <- idx_current + num_beta_int
  idx_gamma <- if(num_gamma_main > 0) (idx_current + 1):(idx_current + num_gamma_main) else integer(0); idx_current <- idx_current + num_gamma_main
  idx_gamma_int <- if(num_gamma_int > 0) (idx_current + 1):(idx_current + num_gamma_int) else integer(0)

  beta_out <- if(length(idx_beta) > 0) combined_Beta_all_samples[idx_beta, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  beta_int_out <- if(length(idx_beta_int) > 0) combined_Beta_all_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_out <- if(length(idx_gamma) > 0) combined_Beta_all_samples[idx_gamma, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_int_out <- if(length(idx_gamma_int) > 0) combined_Beta_all_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)

  if(nrow(beta_out) > 0) rownames(beta_out) <- colnames(X_prognostic)[1:num_beta_main_prog]
  if(nrow(beta_int_out) > 0) rownames(beta_int_out) <- colnames(X_prognostic)[(num_beta_main_prog+1):p_prognostic]
  if(nrow(gamma_out) > 0) rownames(gamma_out) <- colnames(X_treatment_related)[1:num_gamma_main]
  if(nrow(gamma_int_out) > 0) rownames(gamma_int_out) <- colnames(X_treatment_related)[(num_gamma_main+1):p_treatment_related]

  output <- list(
    alpha = combined_alpha_samples,
    beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
    beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    aleph = combined_aleph_samples,
    gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
    gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0)
  )

  if (prognostic_shrinkage == "horseshoe") {
    output$tau_prognostic_global <- combined_tau_prog
    if(p_prognostic>0) output$lambda_prognostic <- t(combined_lambda_prog)
  } else if (prognostic_shrinkage == "linked_shrinkage") {
    output$tau_prognostic_local <- if(is.matrix(combined_tau_prog)) t(combined_tau_prog) else combined_tau_prog
    output$tau_prognostic_interaction <- combined_tau_int_prog
  }
  if (treatment_shrinkage == "horseshoe") {
    output$tau_treatment_global <- combined_tau_treat
    if(p_treatment_related>0) output$lambda_treatment <- t(combined_lambda_treat)
  } else if (treatment_shrinkage == "linked_shrinkage") {
    output$tau_treatment_local <- if(is.matrix(combined_tau_treat)) t(combined_tau_treat) else combined_tau_treat
    output$tau_treatment_interaction <- combined_tau_int_treat
  }

  if (family == "gaussian") {
    output$sigma_sq = combined_sigma_sq_samples
  }

  return(output)
}

#############################################################################################

#' # Log-posterior for Poisson model coefficients (for Metropolis-Hastings step)
#' log_posterior_poisson <- function(coeffs, X_full, y_vec, prior_prec_diag) {
#'   eta <- as.vector(X_full %*% coeffs)
#'   
#'   # Log-likelihood for Poisson: sum(y*eta - exp(eta))
#'   # We drop the log(y!) term as it's constant w.r.t. coefficients
#'   log_lik <- sum(y_vec * eta - exp(eta))
#'   
#'   # Log-prior for coefficients (assuming N(0, V))
#'   log_prior <- -0.5 * as.numeric(t(coeffs) %*% prior_prec_diag %*% coeffs)
#'   
#'   log_post <- log_lik + log_prior
#'   
#'   # Return -Inf if not finite, to ensure rejection in MH step
#'   if (!is.finite(log_post)) return(-Inf)
#'   return(log_post)
#' }
#' 
#' 
#' #' @title Fit BCF-style Linear Model with Flexible Shrinkage Options
#' #' @description R wrapper to run a Gibbs sampler for a complex model with Gaussian,
#' #'   Binomial (logit), or Poisson (log) links.
#' #' @param y_vec Response vector. Continuous, binary (0/1), or integer counts.
#' #' @param X_mat Covariate matrix.
#' #' @param Z_vec Binary treatment indicator vector (0s and 1s).
#' #' @param family Model family: "gaussian", "binomial", or "poisson".
#' #' @param n_iter Total MCMC iterations.
#' #' @param burn_in Number of burn-in iterations.
#' #' @param prognostic_shrinkage Shrinkage method for mu(x) part: "horseshoe" or "linked_shrinkage".
#' #' @param treatment_shrinkage Shrinkage method for tau(x) part: "horseshoe" or "linked_shrinkage".
#' #' @param num_chains Integer, the number of independent MCMC chains to run.
#' #' @param propensity_as_covariate Boolean, if TRUE, propensity scores are estimated
#' #'   and added as a covariate to the prognostic model `mu(x)`.
#' #' @param thin Thinning parameter.
#' #' @param seed Random seed.
#' #' @param verbose Logical, print progress.
#' #' @param ping Frequency of progress messages.
#' #'
#' #' @return A list of combined posterior samples.
#' #' @export
#' #' @importFrom MASS mvrnorm
#' #' @importFrom stats rgamma rexp runif pgamma qgamma rnorm var combn dpois
#' linear_linear_function <- function(
#'     y_vec, X_mat, Z_vec,
#'     family = c("binomial", "gaussian", "poisson"),
#'     n_iter, burn_in,
#'     prognostic_shrinkage = c("horseshoe", "linked_shrinkage"),
#'     treatment_shrinkage = c("horseshoe", "linked_shrinkage"),
#'     num_chains = 1, propensity_as_covariate = FALSE,
#'     thin = 1, seed = 123, verbose = FALSE, ping = 1000,
#'     alpha_global_prior_sd = 10.0, aleph_prior_sd = 10.0
#' ) {
#'   # --- Setup and Validation ---
#'   if (!is.null(seed)) { set.seed(seed); chain_seeds <- sample.int(1e6, num_chains) } else { chain_seeds <- rep(NA, num_chains) }
#'   family <- match.arg(family)
#'   prognostic_shrinkage <- match.arg(prognostic_shrinkage)
#'   treatment_shrinkage <- match.arg(treatment_shrinkage)
#'   if (length(unique(na.omit(y_vec))) <= 1) stop("Outcome variable `y_vec` is constant.")
#'   
#'   N <- nrow(X_mat)
#'   p_main_X_orig <- ncol(X_mat)
#'   X_mat_orig <- X_mat
#'   
#'   # --- Propensity Score Estimation ---
#'   if (propensity_as_covariate) {
#'     if (verbose) message("Estimating propensity scores to use as a covariate...")
#'     if (!requireNamespace("stochtree", quietly = TRUE) || !exists("bart", where = "package:stochtree")) {
#'       stop("'stochtree::bart' function not found for propensity estimation. Please install and load the 'stochtree' package.")
#'     }
#'     bart_model_propensity <- stochtree::bart(X_train = X_mat_orig, y_train = as.numeric(Z_vec), num_gfr = 50, num_mcmc = 0)
#'     propensity_scores <- rowMeans(bart_model_propensity$y_hat_train)
#'     X_mat_for_prog <- cbind(X_mat_orig, propensity_score = propensity_scores)
#'     if (verbose) message("Propensity scores added to X_mat for prognostic model.")
#'   } else {
#'     X_mat_for_prog <- X_mat_orig
#'   }
#'   p_main_X_prog <- ncol(X_mat_for_prog)
#'   
#'   # --- 1. Construct Design Matrices (once, outside the chain loop) ---
#'   if(is.null(colnames(X_mat_for_prog))) colnames(X_mat_for_prog) <- paste0("X_prog_", 1:p_main_X_prog)
#'   if(is.null(colnames(X_mat_orig))) colnames(X_mat_orig) <- paste0("X_orig_", 1:p_main_X_orig)
#'   
#'   main_interaction_indices_prog <- create_interaction_pairs_R(p_main_X_prog)
#'   p_interaction_prog <- if(is.matrix(main_interaction_indices_prog)) nrow(main_interaction_indices_prog) else 0
#'   X_prognostic <- X_mat_for_prog
#'   if (p_interaction_prog > 0) {
#'     X_beta_int_mat <- do.call(cbind, lapply(1:p_interaction_prog, function(k) X_mat_for_prog[, main_interaction_indices_prog[k,1]] * X_mat_for_prog[, main_interaction_indices_prog[k,2]]))
#'     colnames(X_beta_int_mat) <- sapply(1:p_interaction_prog, function(k) paste0("beta_int_", paste(colnames(X_mat_for_prog)[main_interaction_indices_prog[k,]], collapse="_")))
#'     X_prognostic <- cbind(X_prognostic, X_beta_int_mat)
#'   }
#'   p_prognostic <- ncol(X_prognostic)
#'   
#'   main_interaction_indices_orig <- create_interaction_pairs_R(p_main_X_orig)
#'   p_interaction_orig <- if(is.matrix(main_interaction_indices_orig)) nrow(main_interaction_indices_orig) else 0
#'   X_treatment_related <- X_mat_orig * Z_vec
#'   colnames(X_treatment_related) <- paste0("Z_", colnames(X_mat_orig))
#'   if (p_interaction_orig > 0) {
#'     base_int_treat <- do.call(cbind, lapply(1:p_interaction_orig, function(k) (X_mat_orig[,main_interaction_indices_orig[k,1]]*X_mat_orig[,main_interaction_indices_orig[k,2]])))
#'     X_gamma_int_mat <- base_int_treat * Z_vec
#'     colnames(X_gamma_int_mat) <- sapply(1:p_interaction_orig, function(k) paste0("Z_int_", paste(colnames(X_mat_orig)[main_interaction_indices_orig[k,]], collapse="_")))
#'     X_treatment_related <- cbind(X_treatment_related, X_gamma_int_mat)
#'   }
#'   p_treatment_related <- ncol(X_treatment_related)
#'   
#'   X_full_intercepts <- cbind(alpha_global = rep(1, N), aleph = Z_vec)
#'   X_full_hs <- cbind(X_prognostic, X_treatment_related)
#'   p_total_hs <- ncol(X_full_hs)
#'   X_all_coeffs <- cbind(X_full_intercepts, X_full_hs)
#'   p_all_coeffs <- ncol(X_all_coeffs)
#'   
#'   # --- MCMC Setup ---
#'   all_chains_results <- list()
#'   for (chain_num in 1:num_chains) {
#'     if(verbose && num_chains > 1) message(sprintf("\n--- Starting Chain %d / %d ---\n", chain_num, num_chains))
#'     if (!is.na(chain_seeds[chain_num])) set.seed(chain_seeds[chain_num])
#'     
#'     # Initialize parameters for this chain
#'     alpha_global <- 0.0; aleph <- 0.0
#'     Beta_all <- rep(0, p_total_hs)
#'     
#'     lambda_prog <- rep(1, p_prognostic); tau_prog <- rep(1, p_main_X_prog); tau_int_prog <- 0.5; tau_glob_prog <- 1.0
#'     lambda_treat <- rep(1, p_treatment_related); tau_treat <- rep(1, p_main_X_orig); tau_int_treat <- 0.5; tau_glob_treat <- 1.0
#'     
#'     if (family == "gaussian") { sigma_sq <- var(y_vec, na.rm=TRUE); if(is.na(sigma_sq) || sigma_sq == 0) sigma_sq <- 1.0 } else if (family == "binomial") { y_adj <- y_vec - 0.5 }
#'     
#'     # Storage for this chain
#'     eff_samp_count <- floor((n_iter - burn_in) / thin)
#'     alpha_global_samples <- numeric(eff_samp_count); aleph_samples <- numeric(eff_samp_count)
#'     Beta_all_samples <- if(p_total_hs > 0) matrix(0, nrow = p_total_hs, ncol = eff_samp_count) else matrix(0,0,0)
#'     
#'     tau_prog_samples <- if(prognostic_shrinkage=="linked_shrinkage") matrix(0, nrow=p_main_X_prog, ncol=eff_samp_count) else numeric(eff_samp_count)
#'     tau_treat_samples <- if(treatment_shrinkage=="linked_shrinkage") matrix(0, nrow=p_main_X_orig, ncol=eff_samp_count) else numeric(eff_samp_count)
#'     tau_int_prog_samples <- numeric(eff_samp_count); tau_int_treat_samples <- numeric(eff_samp_count)
#'     tau_glob_prog_samples <- numeric(eff_samp_count); tau_glob_treat_samples <- numeric(eff_samp_count)
#'     lambda_prog_samples <- if(p_prognostic > 0) matrix(0, nrow=p_prognostic, ncol=eff_samp_count) else matrix(0,0,0)
#'     lambda_treat_samples <- if(p_treatment_related > 0) matrix(0, nrow=p_treatment_related, ncol=eff_samp_count) else matrix(0,0,0)
#'     
#'     if (family == "gaussian") sigma_sq_samples <- numeric(eff_samp_count)
#'     store_idx <- 0
#'     
#'     # MCMC Loop
#'     if(verbose) message(sprintf("MCMC Chain %d (%s model) is running...", chain_num, family))
#'     for (iter_idx in 1:n_iter) {
#'       if (verbose && iter_idx %% ping == 0) message("Chain ", chain_num, " - Iteration: ", iter_idx, "/", n_iter)
#'       
#'       # Define indices for prognostic and treatment effects within Beta_all
#'       idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0)
#'       idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_total_hs else integer(0)
#'       
#'       # --- Update Shrinkage Parameters (conditional on current Beta_all) ---
#'       current_scale_factor <- if (family == "gaussian") sigma_sq else 1.0
#'       
#'       # Prognostic Shrinkage
#'       if (prognostic_shrinkage == "horseshoe") {
#'         if(length(idx_prog) > 0) {
#'           beta_prog <- Beta_all[idx_prog]; temp_lambda_scale_term <- (beta_prog^2) / (2 * current_scale_factor * safe_var_R(tau_glob_prog^2)); nu_lambda_inv <- rexp(length(idx_prog), rate = 1 + (1 / safe_var_R(lambda_prog^2))); lambda_prog_sq_inv <- rexp(length(idx_prog), rate = nu_lambda_inv + temp_lambda_scale_term); lambda_prog <- 1 / sqrt(safe_var_R(lambda_prog_sq_inv)); xi_inv_prog <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_prog^2)); rate_prog <- xi_inv_prog + sum(beta_prog^2 / (2 * current_scale_factor * safe_var_R(lambda_prog^2))); tau_glob_prog <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog))
#'         }
#'       } else { # linked_shrinkage  
#'         if(p_prognostic > 0) {
#'           beta_prog_main <- Beta_all[1:p_main_X_prog]; beta_prog_int <- Beta_all[(p_main_X_prog + 1):p_prognostic]
#'           for(j in 1:p_main_X_prog) {
#'             sum_sq_terms <- beta_prog_main[j]^2; num_terms <- 1
#'             if(p_interaction_prog > 0){ for(k_int in 1:p_interaction_prog){ pair <- main_interaction_indices_prog[k_int,]; if(pair[1] == j) { sum_sq_terms <- sum_sq_terms + beta_prog_int[k_int]^2 / safe_var_R(tau_int_prog^2 * tau_prog[pair[2]]^2); num_terms <- num_terms + 1 } else if (pair[2] == j) { sum_sq_terms <- sum_sq_terms + beta_prog_int[k_int]^2 / safe_var_R(tau_int_prog^2 * tau_prog[pair[1]]^2); num_terms <- num_terms + 1 }}}
#'             tau_prog[j] <- sqrt(rinvgamma_R(1, shape=0.5*num_terms + 0.5, scale=0.5 + 0.5*sum_sq_terms/current_scale_factor))
#'           }
#'           if(p_interaction_prog > 0){ prop_tau_int <- runif(1, 0.01, 2); log_accept_ratio <- loglike_tau_int_R(prop_tau_int, beta_prog_int, main_interaction_indices_prog, tau_prog, current_scale_factor) - loglike_tau_int_R(tau_int_prog, beta_prog_int, main_interaction_indices_prog, tau_prog, current_scale_factor); if(!is.na(log_accept_ratio) && log(runif(1)) < log_accept_ratio) tau_int_prog <- prop_tau_int }
#'         }
#'       }
#'       
#'       # Treatment Shrinkage
#'       if (treatment_shrinkage == "horseshoe") {
#'         if(length(idx_treat) > 0) {
#'           beta_treat <- Beta_all[idx_treat]; temp_lambda_scale_term <- (beta_treat^2) / (2 * current_scale_factor * safe_var_R(tau_glob_treat^2)); nu_lambda_inv <- rexp(length(idx_treat), rate = 1 + (1 / safe_var_R(lambda_treat^2))); lambda_treat_sq_inv <- rexp(length(idx_treat), rate = nu_lambda_inv + temp_lambda_scale_term); lambda_treat <- 1 / sqrt(safe_var_R(lambda_treat_sq_inv)); xi_inv_treat <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_treat^2)); rate_treat <- xi_inv_treat + sum(beta_treat^2 / (2 * current_scale_factor * safe_var_R(lambda_treat^2))); tau_glob_treat <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat))
#'         }
#'       } else { # linked_shrinkage
#'         if(p_treatment_related > 0) {
#'           beta_treat_main <- Beta_all[idx_treat[1:p_main_X_orig]]; beta_treat_int <- Beta_all[idx_treat[(p_main_X_orig + 1):p_treatment_related]]
#'           for(j in 1:p_main_X_orig) {
#'             sum_sq_terms <- beta_treat_main[j]^2; num_terms <- 1
#'             if(p_interaction_orig > 0){ for(k_int in 1:p_interaction_orig){ pair <- main_interaction_indices_orig[k_int,]; if(pair[1] == j) { sum_sq_terms <- sum_sq_terms + beta_treat_int[k_int]^2 / safe_var_R(tau_int_treat^2 * tau_treat[pair[2]]^2); num_terms <- num_terms + 1 } else if (pair[2] == j) { sum_sq_terms <- sum_sq_terms + beta_treat_int[k_int]^2 / safe_var_R(tau_int_treat^2 * tau_treat[pair[1]]^2); num_terms <- num_terms + 1 }}}
#'             tau_treat[j] <- sqrt(rinvgamma_R(1, shape=0.5*num_terms + 0.5, scale=0.5*sum_sq_terms/current_scale_factor))
#'           }
#'           if(p_interaction_orig > 0){ prop_tau_int <- runif(1, 0.01, 2); log_accept_ratio <- loglike_tau_int_R(prop_tau_int, beta_treat_int, main_interaction_indices_orig, tau_treat, current_scale_factor) - loglike_tau_int_R(tau_int_treat, beta_treat_int, main_interaction_indices_orig, tau_treat, current_scale_factor); if(!is.na(log_accept_ratio) && log(runif(1)) < log_accept_ratio) tau_int_treat <- prop_tau_int }
#'         }
#'       }
#'       
#'       # --- Update Main Coefficients (alpha, aleph, Beta_all) ---
#'       if (family == "poisson") {
#'         # --- Metropolis-Hastings for Poisson ---
#'         current_coeffs <- c(alpha_global, aleph, Beta_all)
#'         
#'         # Construct full prior precision matrix for all coeffs
#'         prior_var_beta <- numeric(p_total_hs)
#'         if (prognostic_shrinkage == "horseshoe") { if(length(idx_prog)>0) prior_var_beta[idx_prog] <- lambda_prog^2 * tau_glob_prog^2 } else { if(p_prognostic > 0) { prog_main_indices <- 1:p_main_X_prog; prior_var_beta[idx_prog[prog_main_indices]] <- tau_prog^2; if(p_interaction_prog > 0){ prog_int_indices <- (p_main_X_prog+1):p_prognostic; for(k in 1:p_interaction_prog) { pair <- main_interaction_indices_prog[k,]; prior_var_beta[idx_prog[prog_int_indices[k]]] <- tau_int_prog^2 * tau_prog[pair[1]]^2 * tau_prog[pair[2]]^2 }}}}
#'         if (treatment_shrinkage == "horseshoe") { if(length(idx_treat)>0) prior_var_beta[idx_treat] <- lambda_treat^2 * tau_glob_treat^2 } else { if(p_treatment_related > 0) { treat_main_indices <- 1:p_main_X_orig; prior_var_beta[idx_treat[treat_main_indices]] <- tau_treat^2; if(p_interaction_orig > 0){ treat_int_indices <- (p_main_X_orig+1):p_treatment_related; for(k in 1:p_interaction_orig) { pair <- main_interaction_indices_orig[k,]; prior_var_beta[idx_treat[treat_int_indices[k]]] <- tau_int_treat^2 * tau_treat[pair[1]]^2 * tau_treat[pair[2]]^2 }}}}
#'         
#'         prior_var_all <- c(alpha_global_prior_sd^2, aleph_prior_sd^2, prior_var_beta)
#'         prior_prec_diag_all <- diag(1 / safe_var_R(prior_var_all))
#'         
#'         proposal_sd <- 0.15; proposal_cov <- diag(proposal_sd^2, length(current_coeffs))
#'         proposed_coeffs <- mvrnorm(1, current_coeffs, proposal_cov)
#'         
#'         log_post_current <- log_posterior_poisson(current_coeffs, X_all_coeffs, y_vec, prior_prec_diag_all)
#'         log_post_proposed <- log_posterior_poisson(proposed_coeffs, X_all_coeffs, y_vec, prior_prec_diag_all)
#'         
#'         if (log(runif(1)) < (log_post_proposed - log_post_current)) {
#'           alpha_global <- proposed_coeffs[1]; aleph <- proposed_coeffs[2]; Beta_all <- proposed_coeffs[-(1:2)]
#'         }
#'         
#'       } else {
#'         # --- Gibbs Sampling for Gaussian and Binomial ---
#'         prior_var_beta <- numeric(p_total_hs)
#'         if (prognostic_shrinkage == "horseshoe") { if(length(idx_prog)>0) prior_var_beta[idx_prog] <- lambda_prog^2 * tau_glob_prog^2 } else { if(p_prognostic > 0) { prog_main_indices <- 1:p_main_X_prog; prior_var_beta[idx_prog[prog_main_indices]] <- tau_prog^2; if(p_interaction_prog > 0){ prog_int_indices <- (p_main_X_prog+1):p_prognostic; for(k in 1:p_interaction_prog) { pair <- main_interaction_indices_prog[k,]; prior_var_beta[idx_prog[prog_int_indices[k]]] <- tau_int_prog^2 * tau_prog[pair[1]]^2 * tau_prog[pair[2]]^2 }}}}
#'         if (treatment_shrinkage == "horseshoe") { if(length(idx_treat)>0) prior_var_beta[idx_treat] <- lambda_treat^2 * tau_glob_treat^2 } else { if(p_treatment_related > 0) { treat_main_indices <- 1:p_main_X_orig; prior_var_beta[idx_treat[treat_main_indices]] <- tau_treat^2; if(p_interaction_orig > 0){ treat_int_indices <- (p_main_X_orig+1):p_treatment_related; for(k in 1:p_interaction_orig) { pair <- main_interaction_indices_orig[k,]; prior_var_beta[idx_treat[treat_int_indices[k]]] <- tau_int_treat^2 * tau_treat[pair[1]]^2 * tau_treat[pair[2]]^2 }}}}
#'         prior_prec_diag <- diag(1 / safe_var_R(prior_var_beta), p_total_hs, p_total_hs)
#'         
#'         XtWX <- t(X_full_hs) %*% (w_vec * X_full_hs); post_prec_unscaled <- XtWX + prior_prec_diag
#'         chol_attempt <- try(chol(post_prec_unscaled), silent = TRUE)
#'         if (inherits(chol_attempt, "try-error")) { warning(paste("Cholesky for Beta_all failed at iter", iter_idx, "- skipping update."))
#'         } else { post_cov_unscaled_inv <- chol2inv(chol_attempt); rhs_mean <- t(X_full_hs) %*% (w_vec * z_for_betas); post_mean <- post_cov_unscaled_inv %*% rhs_mean; Beta_all <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean), Sigma = scale_factor * post_cov_unscaled_inv) }
#'       }
#'       
#'       # --- Sample sigma_sq for Gaussian model ---
#'       if (family == "gaussian") { eta_new <- alpha_global + Z_vec * aleph + as.vector(X_full %*% Beta_all); residuals_for_sigma <- y_vec - eta_new; shape_sigma_post <- 0.001 + N / 2; rate_sigma_post <- 0.001 + 0.5 * sum(residuals_for_sigma^2); sigma_sq <- rinvgamma_R(1, shape = shape_sigma_post, scale = rate_sigma_post); sigma_sq <- max(sigma_sq, 1e-9) }
#'       
#'       # --- Store samples ---
#'       if (iter_idx > burn_in && (iter_idx - burn_in) %% thin == 0) {
#'         store_idx <- store_idx + 1
#'         alpha_global_samples[store_idx] <- alpha_global; aleph_samples[store_idx] <- aleph
#'         if(p_total_hs > 0) Beta_all_samples[, store_idx] <- Beta_all
#'         if (family == "gaussian") sigma_sq_samples[store_idx] <- sigma_sq
#'         if(prognostic_shrinkage == "horseshoe") { tau_prog_samples[store_idx] <- tau_glob_prog; if(p_prognostic > 0) lambda_prog_samples[, store_idx] <- lambda_prog }
#'         else { if(p_main_X_prog > 0) tau_prog_samples[, store_idx] <- tau_prog; tau_int_prog_samples[store_idx] <- tau_int_prog }
#'         if (treatment_shrinkage == "horseshoe") { tau_treat_samples[store_idx] <- tau_glob_treat; if(p_treatment_related > 0) lambda_treat_samples[, store_idx] <- lambda_treat }
#'         else { if(p_main_X_orig > 0) tau_treat_samples[, store_idx] <- tau_treat; tau_int_treat_samples[store_idx] <- tau_int_treat }
#'       }
#'     } # End MCMC loop
#'     
#'     all_chains_results[[chain_num]] <- list(alpha=alpha_global_samples, aleph=aleph_samples, Beta_all=Beta_all_samples, sigma_sq=if(family=="gaussian") sigma_sq_samples else NULL,
#'                                             tau_prog=tau_prog_samples, tau_int_prog=tau_int_prog_samples, lambda_prog=lambda_prog_samples,
#'                                             tau_treat=tau_treat_samples, tau_int_treat=tau_int_treat_samples, lambda_treat=lambda_treat_samples)
#'   } # End outer chain loop
#'     combined_alpha_samples <- do.call(c, lapply(all_chains_results, `[[`, "alpha"))
#'     combined_aleph_samples <- do.call(c, lapply(all_chains_results, `[[`, "aleph"))
#'     if(p_total_hs > 0) combined_Beta_all_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "Beta_all")) else combined_Beta_all_samples <- matrix(0,0,0)
#' 
#'     # Combine shrinkage parameters
#'     if(prognostic_shrinkage == "horseshoe") {
#'       combined_tau_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_prog"))
#'       if(p_prognostic>0) combined_lambda_prog <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_prog"))
#'     } else {
#'       if(p_main_X_prog > 0) combined_tau_prog <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_prog))) else combined_tau_prog <- matrix(0, eff_samp_count * num_chains, 0)
#'       combined_tau_int_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_prog"))
#'     }
#'     if(treatment_shrinkage == "horseshoe") {
#'       combined_tau_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_treat"))
#'       if(p_treatment_related>0) combined_lambda_treat <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_treat"))
#'     } else {
#'       if(p_main_X_orig > 0) combined_tau_treat <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_treat))) else combined_tau_treat <- matrix(0, eff_samp_count * num_chains, 0)
#'       combined_tau_int_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_treat"))
#'     }
#'   # --- Combine and Unpack Results ---
#'   
#'   # Unpack Beta_all_samples
#'   num_beta_main_prog <- ncol(X_mat_for_prog)
#'   num_beta_int  <- p_interaction_prog
#'   num_gamma_main <- ncol(X_mat_orig)
#'   num_gamma_int  <- p_interaction_orig
#'   
#'   idx_current <- 0
#'   idx_beta <- (idx_current + 1):(idx_current + num_beta_main_prog); idx_current <- idx_current + num_beta_main_prog
#'   idx_beta_int <- if(num_beta_int > 0) (idx_current + 1):(idx_current + num_beta_int) else integer(0); idx_current <- idx_current + num_beta_int
#'   idx_gamma <- if(num_gamma_main > 0) (idx_current + 1):(idx_current + num_gamma_main) else integer(0); idx_current <- idx_current + num_gamma_main
#'   idx_gamma_int <- if(num_gamma_int > 0) (idx_current + 1):(idx_current + num_gamma_int) else integer(0)
#'   
#'   beta_out <- if(length(idx_beta) > 0) combined_Beta_all_samples[idx_beta, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
#'   beta_int_out <- if(length(idx_beta_int) > 0) combined_Beta_all_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
#'   gamma_out <- if(length(idx_gamma) > 0) combined_Beta_all_samples[idx_gamma, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
#'   gamma_int_out <- if(length(idx_gamma_int) > 0) combined_Beta_all_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
#'   
#'   if(nrow(beta_out) > 0) rownames(beta_out) <- colnames(X_prognostic)[1:num_beta_main_prog]
#'   if(nrow(beta_int_out) > 0) rownames(beta_int_out) <- colnames(X_prognostic)[(num_beta_main_prog+1):p_prognostic]
#'   if(nrow(gamma_out) > 0) rownames(gamma_out) <- colnames(X_treatment_related)[1:num_gamma_main]
#'   if(nrow(gamma_int_out) > 0) rownames(gamma_int_out) <- colnames(X_treatment_related)[(num_gamma_main+1):p_treatment_related]
#'   
#'   output <- list(
#'     alpha = combined_alpha_samples,
#'     beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
#'     beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
#'     aleph = combined_aleph_samples,
#'     gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
#'     gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0)
#'   )
#'   
#'   if (prognostic_shrinkage == "horseshoe") { 
#'     output$tau_prognostic_global <- combined_tau_prog
#'     if(p_prognostic>0) output$lambda_prognostic <- t(combined_lambda_prog) 
#'   } else { 
#'     output$tau_prognostic_local <- if(is.matrix(combined_tau_prog)) t(combined_tau_prog) else combined_tau_prog
#'     output$tau_prognostic_interaction <- combined_tau_int_prog 
#'   }
#'   if (treatment_shrinkage == "horseshoe") { 
#'     output$tau_treatment_global <- combined_tau_treat
#'     if(p_treatment_related>0) output$lambda_treatment <- t(combined_lambda_treat)
#'   } else { 
#'     output$tau_treatment_local <- if(is.matrix(combined_tau_treat)) t(combined_tau_treat) else combined_tau_treat
#'     output$tau_treatment_interaction <- combined_tau_int_treat 
#'   }
#'   
#'   if (family == "gaussian") {
#'     output$sigma_sq = combined_sigma_sq_samples
#'   }
#'   return(output)
#' }
#########################################################################
##################################################################
# simulate_poisson_data <- function(n = 500, p = 4, seed = 42) {
#   if (!is.null(seed)) set.seed(seed)
#   
#   # Generate continuous and binary covariates
#   X_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)
#   colnames(X_mat) <- paste0("X", 1:p)
#   
#   # Generate binary treatment Z (0/1)
#   Z_vec <- rbinom(n, 1, 0.5)
#   
#   # Define true parameters for the linear predictor (eta)
#   true_alpha_global <- 0.1
#   true_beta <- c(0.4, -0.3, 0, 0.2) # Effects for X1, X2, X3, X4
#   true_beta_int <- c(0.3, -0.25)     # Effects for X1*X2 and X1*X3 interactions
#   
#   true_aleph <- 0.5                  # Main effect of treatment Z
#   true_gamma <- c(-0.4, 0.2, 0, 0)   # How Z modifies the effect of X1, X2, X3, X4
#   true_gamma_int <- c(-0.3)          # How Z modifies the effect of the X1*X2 interaction
#   
#   # Construct the full design matrices based on the true model structure
#   # Prognostic part
#   X_prog_main <- X_mat
#   X_prog_int <- cbind(X_mat[,1]*X_mat[,2], X_mat[,1]*X_mat[,3])
#   
#   # Treatment part
#   X_treat_main <- X_mat * Z_vec
#   X_treat_int <- (X_mat[,1]*X_mat[,2]) * Z_vec
#   
#   # Construct linear predictor eta = mu(x) + Z*tau(x)
#   mu_x <- true_alpha_global + 
#     as.vector(X_prog_main %*% true_beta) + 
#     as.vector(X_prog_int %*% true_beta_int)
#   
#   tau_x <- true_aleph + 
#     as.vector(X_mat %*% true_gamma) + # Note: gamma applies to original X
#     as.vector((X_mat[,1]*X_mat[,2]) * true_gamma_int) # and interaction
#   
#   eta <- mu_x + Z_vec * tau_x
#   
#   # Calculate the rate lambda using the exponential link function
#   lambda <- exp(eta)
#   
#   # Generate count outcome Y from the Poisson distribution
#   y_vec <- rpois(n, lambda = lambda)
#   
#   # Return data and true parameters
#   list(
#     y_vec = y_vec,
#     X_mat = X_mat,
#     Z_vec = Z_vec,
#     true_params = list(
#       alpha = true_alpha_global,
#       beta = true_beta,
#       beta_int = true_beta_int,
#       aleph = true_aleph,
#       gamma = true_gamma,
#       gamma_int = true_gamma_int
#     )
#   )
# }
####################################################################################################
# 1. Generate Poisson Data
# sim_data_poisson <- simulate_poisson_data(n = 1000, p = 4, seed = 2025)
# y_sim <- sim_data_poisson$y_vec
# X_sim <- sim_data_poisson$X_mat
# Z_sim <- sim_data_poisson$Z_vec
# true_params <- sim_data_poisson$true_params
# 
# cat("\n--- True Parameters from Poisson Simulation ---\n")
# print(true_params)
# 
# # 2. Fit the Model using the `linear_linear_function`
# # Note: The Metropolis-Hastings sampler for Poisson is slower than Gibbs.
# # A higher number of iterations is recommended for full convergence.
# n_iterations_example <- 10000
# n_burnin_example <- 2000
# 
# cat(sprintf("\n--- Running Sampler for Poisson Family (Iter=%d, Burn=%d) ---\n", 
#             n_iterations_example, n_burnin_example))
# 
# # Ensure the main sampler function is loaded into your environment
# if (!exists("linear_linear_function")) {
#   stop("The 'linear_linear_function' is not loaded. Please source the file containing it.")
# }
# 
# # For this test, let's use horseshoe shrinkage for both components
# fit_poisson <- linear_linear_function(
#   y_vec = y_sim,
#   X_mat = X_sim,
#   Z_vec = Z_sim,
#   family = "gaussian",
#   n_iter = n_iterations_example,
#   burn_in = n_burnin_example,
#   prognostic_shrinkage = "horseshoe",
#   treatment_shrinkage = "horseshoe",
#   num_chains = 1, # Use 1 chain for a quick test
#   verbose = TRUE,
#   ping = 1000
# )
# 
# cat("\n--- MCMC Sampling Complete ---\n")
# 
# # 3. Analyze Results: Compare posterior means to true values
# cat("\n--- Posterior Means vs True Parameters ---\n")
# cat("Alpha (Global Intercept): True =", true_params$alpha, 
#     ", Estimated =", mean(fit_poisson$alpha), "\n")
# 
# cat("\nBeta (Prognostic Main Effects):\n")
# beta_comp <- data.frame(
#   True = true_params$beta,
#   Estimated_Mean = colMeans(fit_poisson$beta)
# )
# print(beta_comp)
# 
# cat("\nBeta Interaction (Prognostic Interactions):\n")
# beta_int_comp <- data.frame(
#   True = true_params$beta_int,
#   Estimated_Mean = colMeans(fit_poisson$beta_interaction)
# )
# print(beta_int_comp)
# 
# cat("\nAleph (Treatment Intercept Modifier):\n True =", true_params$aleph, 
#     ", Estimated =", mean(fit_poisson$aleph), "\n")
# 
# cat("\nGamma (Z*X Interaction Effects):\n")
# gamma_comp <- data.frame(
#   True = true_params$gamma,
#   Estimated_Mean = colMeans(fit_poisson$gamma)
# )
# print(gamma_comp)
# 
# cat("\nGamma Interaction (Z*XX Interaction Effects):\n")
# gamma_int_comp <- data.frame(
#   True = true_params$gamma_int,
#   Estimated_Mean = colMeans(fit_poisson$gamma_int)
# )
# print(gamma_int_comp)

######################## GAUSSIAN FAMILY. ########################

# # --- Example Usage (using the simple simulation) ---
# source('R/simul_1.R')
# 
# data <- generate_data_2(500, is_te_hetero = T, is_mu_nonlinear = F, seed = 40, RCT = FALSE, z_diff = F, contrast_binary = T)
# X <- as.matrix(sapply(data[, c(1:6)], as.numeric))
# y <- as.numeric(data$y)
# z <- as.numeric(data$z)
# 
# cat("\n--- Running Grouped Horseshoe Logistic Sampler ---\n")
# start.time <- Sys.time()
# fit_grouped_hs <- fit_grouped_horseshoes_R(
#   y_vec = y,
#   X_mat = X,
#   Z_vec = z,
#   family = "gaussian",
#   n_iter = 4000, # Increase for better results
#   burn_in = 1000,
#   num_chains = 3,
#   propensity_as_covariate = T,
#   method_tau_prognostic = "halfCauchy", tau_prognostic_init = 0.1,
#   method_tau_treatment = "halfCauchy", tau_treatment_init = 0.1,
#   method_tau_overall = "fixed", tau_overall_init = 1,
#   alpha_global_prior_sd = 5.0,
#   aleph_prior_sd = 5.0,
#   thin = 1,
#   seed = 103,
#   verbose = TRUE,
#   ping = 500
# )
# start.time <- Sys.time()
# test_linked <- linear_linear_function(
#     y_vec = y, X_mat = X, Z_vec = z,
#     family = "gaussian",
#     n_iter = 10000, burn_in = 2000,
#     prognostic_shrinkage = "linked_shrinkage",
#     treatment_shrinkage = "horseshoe",
#     standardize_cov = F,
#     interaction_rule = "continuous_or_binary",
#     num_chains = 3, propensity_as_covariate = TRUE,
#     thin = 1, seed = 123, verbose = TRUE, ping = 100
# )
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# cat("\n--- Grouped HS: Posterior Means vs True Parameters ---\n")
# cat("Alpha_global (Overall Intercept): True =", true_params_check$alpha,
#     ", Estimated =", mean(fit_grouped_hs$alpha), "\n")
# 
# if(ncol(fit_grouped_hs$beta) > 0) {
#   beta_comp_ghs <- data.frame(True = true_params_check$beta, Estimated = colMeans(fit_grouped_hs$beta))
#   rownames(beta_comp_ghs) <- paste0("Beta_", colnames(X_sim))
#   print(beta_comp_ghs)
# }
# 
# # For beta_interaction, true is 0 in simple DGP
# if (ncol(fit_grouped_hs$beta_interaction) > 0) {
#   cat("\nBeta_interaction (Prognostic XX - True should be ~0 for simple DGP):\n")
#   print(colMeans(fit_grouped_hs$beta_interaction))
# }
# 
# cat("\nAleph (Treatment Intercept Modifier): True =", true_params_check$aleph,
#     ", Estimated =", mean(fit_grouped_hs$aleph), "\n")
# 
# if(ncol(fit_grouped_hs$gamma) > 0) {
#   gamma_comp_ghs <- data.frame(True = true_params_check$gamma, Estimated = colMeans(fit_grouped_hs$gamma))
#   rownames(gamma_comp_ghs) <- paste0("Gamma_Z_", colnames(X_sim))
#   print(gamma_comp_ghs)
# }
# 
# # For gamma_int, true is 0 in simple DGP
# if (ncol(fit_grouped_hs$gamma_int) > 0) {
#   cat("\nGamma_int (Z*XX - True should be ~0 for simple DGP")
#   }
#


#' @title Fit BCF-style Linear Model with Flexible Shrinkage Options
#' @description R wrapper to run a Gibbs sampler for a complex logistic or linear model.
#'
#' @param y_vec Response vector. Binary (0s and 1s) for binomial, continuous for gaussian.
#' @param X_mat Covariate matrix.
#' @param Z_vec Binary treatment indicator vector (0s and 1s).
#' @param family Model family, either "binomial" (default) or "gaussian".
#' @param n_iter Total MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param prognostic_shrinkage Shrinkage method for mu(x) part: "horseshoe", "linked_shrinkage", or "none".
#' @param treatment_shrinkage Shrinkage method for tau(x) part: "horseshoe", "linked_shrinkage", or "none".
#' @param no_shrinkage_prior_variance Constant for the prior variance if shrinkage is "none". Default is 1e4.
#' @param standardize_cov A boolean flag. If `TRUE` (default), the function standardizes numeric
#'  variables and creates dummy variables for categorical ones. If `FALSE`, it bypasses
#'  transformations and returns the original data as a matrix, but still generates the metadata.
#' @param interaction_rule A character string specifying which interactions are allowed.
#' @param cat_coding_method A character string for categorical variable contrast coding.
#' @param num_chains Integer, the number of independent MCMC chains to run.
#' @param propensity_as_covariate Boolean, if TRUE, propensity scores are estimated
#'    and added as a covariate to the prognostic model `mu(x)`.
#' @param regularize_ate Boolean, if `FALSE` (default), `aleph` (the Average Treatment Effect) is given a wide Normal prior.
#'   If `TRUE`, `aleph` is incorporated into the treatment-related coefficients and regularized
#'   along with them (e.g., via Horseshoe shrinkage).
#' @param thin Thinning parameter.
#' @param seed Random seed.
#' @param verbose Logical, print progress.
#' @param ping Frequency of progress messages.
#'
#' @return A list of combined posterior samples.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma rexp runif pgamma qgamma rnorm var combn
linear_linear_function <- function(
    y_vec, X_mat, Z_vec,
    family = c("binomial", "gaussian"),
    n_iter, burn_in,
    prognostic_shrinkage = c("horseshoe", "linked_shrinkage", "none"),
    treatment_shrinkage = c("horseshoe", "linked_shrinkage", "none"),
    no_shrinkage_prior_variance = 1e4,
    standardize_cov = F,
    interaction_rule = c("continuous", "continuous_or_binary", "all"),
    cat_coding_method = "difference",
    num_chains = 1, propensity_as_covariate = FALSE,
    thin = 1, seed = 123, verbose = FALSE, ping = 1000,
    alpha_global_prior_sd = 10.0, aleph_prior_sd = 10.0,
    regularize_ate = FALSE
) {
  # --- Setup and Validation ---
  if (!is.null(seed)) { set.seed(seed)
    chain_seeds <- sample.int(1e6, num_chains) } else {
      chain_seeds <- rep(NA, num_chains)
    }
  family <- match.arg(family)
  prognostic_shrinkage <- match.arg(prognostic_shrinkage)
  treatment_shrinkage <- match.arg(treatment_shrinkage)
  interaction_rule <- match.arg(interaction_rule)
  
  # This function should be defined in your environment, e.g., from a helper script
  # handled_data_list <- standardize_X_by_index(...)
  # For standalone execution, we'll create a dummy boolean vector
  p_main_X_orig_temp <- ncol(X_mat)
  boolean_continous <- rep(TRUE, p_main_X_orig_temp)
  
  if (length(unique(na.omit(y_vec))) <= 1) stop("Outcome variable `y_vec` is constant.")
  
  N <- nrow(X_mat)
  p_main_X_orig <- ncol(X_mat)
  X_mat_orig <- X_mat
  
  # --- Propensity Score Estimation ---
  if (propensity_as_covariate) {
    if (verbose) message("Estimating propensity scores...")
    if (!requireNamespace("stochtree", quietly = TRUE) || !exists("bart", where = "package:stochtree")) {
      stop("'stochtree::bart' not found. Please install/load 'stochtree'.")
    }
    bart_model_propensity <- stochtree::bart(X_train = X_mat_orig, y_train = as.numeric(Z_vec), num_gfr = 50, num_mcmc = 0)
    propensity_scores <- rowMeans(bart_model_propensity$y_hat_train)
    X_mat_for_prog <- cbind(X_mat_orig, propensity_score = propensity_scores)
  } else {
    X_mat_for_prog <- X_mat_orig
  }
  p_main_X_prog <- ncol(X_mat_for_prog)
  
  # --- 1. Construct Design Matrices ---
  if(is.null(colnames(X_mat_for_prog))) colnames(X_mat_for_prog) <- paste0("X_prog_", 1:p_main_X_prog)
  if(is.null(colnames(X_mat_orig))) colnames(X_mat_orig) <- paste0("X_orig_", 1:p_main_X_orig)
  
  main_interaction_indices_prog <- interaction_pairs(p_main_X_orig, boolean_continous)
  p_interaction_prog <- if(is.matrix(main_interaction_indices_prog)) ncol(main_interaction_indices_prog) else 0
  X_prognostic <- X_mat_for_prog
  if (p_interaction_prog > 0) {
    X_beta_int_mat <- do.call(cbind, lapply(1:p_interaction_prog, function(k) {
      idx1 <- main_interaction_indices_prog[1, k]; idx2 <- main_interaction_indices_prog[2, k]
      (X_mat_for_prog[, idx1] * X_mat_for_prog[, idx2])
    }))
    colnames(X_beta_int_mat) <- sapply(1:p_interaction_prog, function(k) {
      idx1 <- main_interaction_indices_prog[1, k]; idx2 <- main_interaction_indices_prog[2, k]
      paste0("beta_int_", paste(colnames(X_mat_for_prog)[c(idx1, idx2)], collapse = "_"))
    })
    X_prognostic <- cbind(X_prognostic, X_beta_int_mat)
  }
  p_prognostic <- ncol(X_prognostic)
  
  main_interaction_indices_orig <- interaction_pairs(p_main_X_orig, boolean_continous)
  p_interaction_orig <- if (is.matrix(main_interaction_indices_orig)) ncol(main_interaction_indices_orig) else 0
  
  if (regularize_ate) {
    if(verbose) message("Regularizing ATE (aleph)...")
    X_treat_list <- list(Z_intercept_aleph = Z_vec)
  } else {
    X_treat_list <- list()
  }
  
  X_gamma_mat <- X_mat_orig * Z_vec
  colnames(X_gamma_mat) <- paste0("Z_", colnames(X_mat_orig))
  X_treat_list$gamma_main <- X_gamma_mat
  
  if (p_interaction_orig > 0) {
    base_int_treat <- do.call(cbind, lapply(1:p_interaction_orig, function(k) {
      idx1 <- main_interaction_indices_orig[1, k]; idx2 <- main_interaction_indices_orig[2, k]
      (X_mat_orig[, idx1] * X_mat_orig[, idx2])
    }))
    X_gamma_int_mat <- base_int_treat * Z_vec
    colnames(X_gamma_int_mat) <- sapply(1:p_interaction_orig, function(k) {
      idx1 <- main_interaction_indices_orig[1, k]; idx2 <- main_interaction_indices_orig[2, k]
      paste0("Z_int_", paste(colnames(X_mat_orig)[c(idx1, idx2)], collapse = "_"))
    })
    X_treat_list$gamma_int <- X_gamma_int_mat
  }
  X_treatment_related <- do.call(cbind, X_treat_list)
  p_treatment_related <- ncol(X_treatment_related)
  X_full <- cbind(X_prognostic, X_treatment_related)
  p_total <- ncol(X_full)
  
  # --- MCMC Setup ---
  all_chains_results <- list()
  for (chain_num in 1:num_chains) {
    if(verbose && num_chains > 1) message(sprintf("\n--- Chain %d / %d ---", chain_num, num_chains))
    if (!is.na(chain_seeds[chain_num])) set.seed(chain_seeds[chain_num])
    
    alpha_global <- 0.0
    if (!regularize_ate) { aleph <- 0.0 }
    Beta_all <- rep(0, p_total)
    lambda_prog <- rep(1, p_prognostic); tau_glob_prog <- 1.0
    lambda_treat <- rep(1, p_treatment_related); tau_glob_treat <- 1.0
    
    if (family == "gaussian") { sigma_sq <- var(y_vec, na.rm=TRUE); if(is.na(sigma_sq) || sigma_sq == 0) sigma_sq <- 1.0 } else { y_adj <- y_vec - 0.5 }
    
    eff_samp_count <- floor((n_iter - burn_in) / thin)
    alpha_global_samples <- numeric(eff_samp_count)
    if (!regularize_ate) { aleph_samples <- numeric(eff_samp_count) }
    Beta_all_samples <- if(p_total > 0) matrix(0, nrow = p_total, ncol = eff_samp_count) else matrix(0,0,0)
    tau_glob_prog_samples <- numeric(eff_samp_count)
    tau_glob_treat_samples <- numeric(eff_samp_count)
    lambda_prog_samples <- if(p_prognostic > 0) matrix(0, nrow=p_prognostic, ncol=eff_samp_count) else matrix(0,0,0)
    lambda_treat_samples <- if(p_treatment_related > 0) matrix(0, nrow=p_treatment_related, ncol=eff_samp_count) else matrix(0,0,0)
    if (family == "gaussian") sigma_sq_samples <- numeric(eff_samp_count)
    store_idx <- 0
    
    if(verbose) message(sprintf("MCMC Chain %d (%s model) running...", chain_num, family))
    for (iter_idx in 1:n_iter) {
      if (verbose && iter_idx %% ping == 0) message(sprintf("Chain %d, Iter: %d/%d", chain_num, iter_idx, n_iter))
      
      eta_hs_part <- if(p_total > 0) as.vector(X_full %*% Beta_all) else rep(0, N)
      eta <- if (regularize_ate) alpha_global + eta_hs_part else alpha_global + Z_vec * aleph + eta_hs_part
      
      if (family == "binomial") {
        omega <- BayesLogit::rpg(N, 1, eta); omega[omega < 1e-9] <- 1e-9
        z_aug <- y_adj / omega; w_vec <- omega; scale_factor <- 1.0
      } else { z_aug <- y_vec; w_vec <- rep(1, N); scale_factor <- sigma_sq }
      
      resid_for_alpha <- if (regularize_ate) z_aug - eta_hs_part else z_aug - (Z_vec * aleph + eta_hs_part)
      prior_prec_alpha <- 1 / safe_var_R(alpha_global_prior_sd^2)
      post_prec_alpha <- sum(w_vec) / scale_factor + prior_prec_alpha
      post_mean_alpha <- (sum(w_vec * resid_for_alpha) / scale_factor) / post_prec_alpha
      alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
      
      if (!regularize_ate) {
        resid_for_aleph <- z_aug - (alpha_global + eta_hs_part)
        prior_prec_aleph <- 1 / safe_var_R(aleph_prior_sd^2)
        post_prec_aleph <- sum(w_vec * Z_vec^2) / scale_factor + prior_prec_aleph
        post_mean_aleph <- (sum(w_vec * Z_vec * resid_for_aleph) / scale_factor) / post_prec_aleph
        aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
      }
      
      if (p_total > 0) {
        z_for_betas <- if (regularize_ate) z_aug - alpha_global else z_aug - alpha_global - Z_vec * aleph
        
        idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0)
        idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_total else integer(0)
        
        # Update Prognostic Shrinkage
        if (prognostic_shrinkage == "horseshoe" && length(idx_prog) > 0) {
          beta_prog <- Beta_all[idx_prog]
          temp_lambda_scale_term <- (beta_prog^2) / (2 * scale_factor * safe_var_R(tau_glob_prog^2))
          nu_lambda_inv <- rexp(length(idx_prog), rate = 1 + (1 / safe_var_R(lambda_prog^2)))
          lambda_prog_sq_inv <- rexp(length(idx_prog), rate = nu_lambda_inv + temp_lambda_scale_term)
          lambda_prog <- 1 / sqrt(safe_var_R(lambda_prog_sq_inv))
          xi_inv_prog <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_prog^2))
          rate_prog <- xi_inv_prog + sum(beta_prog^2 / (2 * scale_factor * safe_var_R(lambda_prog^2)))
          tau_glob_prog <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog))
        }
        
        # Update Treatment Shrinkage
        if (treatment_shrinkage == "horseshoe" && length(idx_treat) > 0) {
          beta_treat <- Beta_all[idx_treat]
          temp_lambda_scale_term <- (beta_treat^2) / (2 * scale_factor * safe_var_R(tau_glob_treat^2))
          nu_lambda_inv <- rexp(length(idx_treat), rate = 1 + (1 / safe_var_R(lambda_treat^2)))
          lambda_treat_sq_inv <- rexp(length(idx_treat), rate = nu_lambda_inv + temp_lambda_scale_term)
          lambda_treat <- 1 / sqrt(safe_var_R(lambda_treat_sq_inv))
          xi_inv_treat <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_treat^2))
          rate_treat <- xi_inv_treat + sum(beta_treat^2 / (2 * scale_factor * safe_var_R(lambda_treat^2)))
          tau_glob_treat <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat))
        }
        
        prior_var_beta <- numeric(p_total)
        if(length(idx_prog) > 0) {
          prior_var_beta[idx_prog] <- if(prognostic_shrinkage == "horseshoe") lambda_prog^2 * tau_glob_prog^2 else no_shrinkage_prior_variance
        }
        if(length(idx_treat) > 0) {
          prior_var_beta[idx_treat] <- if(treatment_shrinkage == "horseshoe") lambda_treat^2 * tau_glob_treat^2 else no_shrinkage_prior_variance
        }
        
        prior_prec_diag_mat <- diag(1 / safe_var_R(prior_var_beta), p_total, p_total)
        post_prec_matrix <- (t(X_full) %*% (w_vec * X_full)) / scale_factor + prior_prec_diag_mat
        
        chol_attempt <- try(chol(post_prec_matrix), silent = TRUE)
        if (inherits(chol_attempt, "try-error")) {
          warning(paste("Cholesky failed at iter", iter_idx, "- skipping update."))
        } else {
          post_cov_matrix <- chol2inv(chol_attempt)
          rhs_mean_scaled <- (t(X_full) %*% (w_vec * z_for_betas)) / scale_factor
          post_mean <- post_cov_matrix %*% rhs_mean_scaled
          Beta_all <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean), Sigma = post_cov_matrix)
        }
      }
      
      if (family == "gaussian") {
        eta_new <- if(regularize_ate) alpha_global + as.vector(X_full %*% Beta_all) else alpha_global + Z_vec * aleph + as.vector(X_full %*% Beta_all)
        residuals_for_sigma <- y_vec - eta_new
        sigma_sq <- rinvgamma_R(1, shape = 0.001 + N/2, scale = 0.001 + 0.5 * sum(residuals_for_sigma^2))
        sigma_sq <- max(sigma_sq, 1e-9)
      }
      
      if (iter_idx > burn_in && (iter_idx - burn_in) %% thin == 0) {
        store_idx <- store_idx + 1
        alpha_global_samples[store_idx] <- alpha_global
        if (!regularize_ate) { aleph_samples[store_idx] <- aleph }
        if(p_total > 0) Beta_all_samples[, store_idx] <- Beta_all
        if(prognostic_shrinkage == "horseshoe") { tau_glob_prog_samples[store_idx] <- tau_glob_prog; if(p_prognostic > 0) lambda_prog_samples[, store_idx] <- lambda_prog }
        if(treatment_shrinkage == "horseshoe") { tau_glob_treat_samples[store_idx] <- tau_glob_treat; if(p_treatment_related > 0) lambda_treat_samples[, store_idx] <- lambda_treat }
        if (family == "gaussian") sigma_sq_samples[store_idx] <- sigma_sq
      }
    }
    
    chain_results <- list(alpha=alpha_global_samples, Beta_all=Beta_all_samples,
                          sigma_sq=if(family=="gaussian") sigma_sq_samples else NULL,
                          tau_prog=tau_glob_prog_samples, lambda_prog=lambda_prog_samples,
                          tau_treat=tau_glob_treat_samples, lambda_treat=lambda_treat_samples)
    if (!regularize_ate) { chain_results$aleph <- aleph_samples }
    all_chains_results[[chain_num]] <- chain_results
  }
  
  # --- Combine and Unpack Results ---
  if(verbose) message("Combining samples from all chains...")
  combined_alpha_samples <- do.call(c, lapply(all_chains_results, `[[`, "alpha"))
  if (!regularize_ate) {
    combined_aleph_samples <- do.call(c, lapply(all_chains_results, `[[`, "aleph"))
  }
  if(p_total > 0) combined_Beta_all_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "Beta_all")) else combined_Beta_all_samples <- matrix(0,0,0)
  
  if(prognostic_shrinkage == "horseshoe") {
    combined_tau_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_prog"))
    if(p_prognostic>0) combined_lambda_prog <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_prog"))
  }
  if(treatment_shrinkage == "horseshoe") {
    combined_tau_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_treat"))
    if(p_treatment_related>0) combined_lambda_treat <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_treat"))
  }
  if (family == "gaussian") { 
    combined_sigma_sq_samples <- do.call(c, lapply(all_chains_results, `[[`, "sigma_sq")) 
  }
  
  num_beta_main_prog <- ncol(X_mat_for_prog)
  num_beta_int  <- p_interaction_prog
  num_gamma_main_orig <- ncol(X_mat_orig)
  num_gamma_int_orig  <- p_interaction_orig
  
  idx_current <- 0
  idx_beta <- (idx_current + 1):(idx_current + num_beta_main_prog); idx_current <- idx_current + num_beta_main_prog
  idx_beta_int <- if(num_beta_int > 0) (idx_current + 1):(idx_current + num_beta_int) else integer(0); idx_current <- idx_current + num_beta_int
  
  if(regularize_ate) {
    idx_aleph_reg <- idx_current + 1
    idx_gamma <- if(num_gamma_main_orig > 0) (idx_current + 2):(idx_current + 1 + num_gamma_main_orig) else integer(0); idx_current <- idx_current + 1 + num_gamma_main_orig
  } else {
    idx_aleph_reg <- integer(0)
    idx_gamma <- if(num_gamma_main_orig > 0) (idx_current + 1):(idx_current + num_gamma_main_orig) else integer(0); idx_current <- idx_current + num_gamma_main_orig
  }
  idx_gamma_int <- if(num_gamma_int_orig > 0) (idx_current + 1):(idx_current + num_gamma_int_orig) else integer(0)
  
  beta_out <- if(length(idx_beta) > 0) combined_Beta_all_samples[idx_beta, , drop=FALSE] else matrix(0,0,0)
  beta_int_out <- if(length(idx_beta_int) > 0) combined_Beta_all_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,0)
  gamma_out <- if(length(idx_gamma) > 0) combined_Beta_all_samples[idx_gamma, , drop=FALSE] else matrix(0,0,0)
  gamma_int_out <- if(length(idx_gamma_int) > 0) combined_Beta_all_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,0)
  
  if(regularize_ate) {
    combined_aleph_samples <- as.vector(combined_Beta_all_samples[idx_aleph_reg, , drop = FALSE])
  }
  
  if(nrow(beta_out) > 0) rownames(beta_out) <- colnames(X_prognostic)[1:num_beta_main_prog]
  if(nrow(beta_int_out) > 0) rownames(beta_int_out) <- colnames(X_prognostic)[(num_beta_main_prog+1):p_prognostic]
  if(nrow(gamma_out) > 0) rownames(gamma_out) <- colnames(X_treatment_related)[idx_gamma - p_prognostic]
  if(nrow(gamma_int_out) > 0) rownames(gamma_int_out) <- colnames(X_treatment_related)[idx_gamma_int - p_prognostic]
  
  output <- list(
    alpha = combined_alpha_samples,
    beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
    beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    aleph = combined_aleph_samples,
    gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
    gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0)
  )
  
  if (prognostic_shrinkage == "horseshoe") {
    output$tau_prognostic_global <- combined_tau_prog
    if(p_prognostic>0) output$lambda_prognostic <- t(combined_lambda_prog)
  }
  if (treatment_shrinkage == "horseshoe") {
    output$tau_treatment_global <- combined_tau_treat
    if(p_treatment_related>0) output$lambda_treatment <- t(combined_lambda_treat)
  }
  
  if (family == "gaussian") {
    output$sigma_sq = combined_sigma_sq_samples
  }
  
  return(output)
}