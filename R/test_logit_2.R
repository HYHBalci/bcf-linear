# Helper functions (assuming create_interaction_pairs_R, rinvgamma_R, safe_var_R are available)
# If not, here they are for completeness:
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
  return(max(x, 1e-9))
}


#' Fit Custom Logistic Model with Grouped Horseshoe Priors
#'
#' Model: logit(p) = alpha_global + X*beta + XX*beta_int + 
#'                   Z*(aleph + X*gamma + XX*gamma_int)
#' - alpha_global and aleph get Normal priors.
#' - (beta, beta_int) get Horseshoe prior with global tau_prognostic.
#' - (gamma, gamma_int) get Horseshoe prior with global tau_treatment.
#' - An optional tau_overall can shrink both groups further.
#'
#' @param y_vec Binary response vector (0s and 1s).
#' @param X_mat Covariate matrix.
#' @param Z_vec Binary treatment indicator vector (0s and 1s).
#' @param n_iter Total MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param method_tau_prognostic Method for global shrinkage for prognostic terms.
#' @param tau_prognostic_init Initial/fixed value for prognostic global shrinkage.
#' @param method_tau_treatment Method for global shrinkage for treatment-related terms.
#' @param tau_treatment_init Initial/fixed value for treatment-related global shrinkage.
#' @param method_tau_overall Method for overall global shrinkage (shrinks both groups).
#' @param tau_overall_init Initial/fixed value for overall global shrinkage.
#' @param alpha_global_prior_mean Mean for Normal prior on alpha_global. Default 0.
#' @param alpha_global_prior_sd SD for Normal prior on alpha_global. Default 10.
#' @param aleph_prior_mean Mean for Normal prior on aleph. Default 0.
#' @param aleph_prior_sd SD for Normal prior on aleph. Default 10. (Diffuse prior)
#' @param thin Thinning parameter.
#' @param seed Random seed.
#' @param verbose Logical, print progress.
#' @param ping Frequency of progress messages.
#'
#' @return A list of MCMC samples.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma rexp runif pgamma qgamma rnorm
# Assuming BayesLogit::rpg or pgdraw::pgdraw is available for Polya-Gamma sampling.

fit_grouped_horseshoes_logistic_R <- function(
    y_vec, X_mat, Z_vec,
    n_iter, burn_in,
    method_tau_prognostic = "halfCauchy", tau_prognostic_init = 1.0,
    method_tau_treatment = "halfCauchy", tau_treatment_init = 1.0,
    method_tau_overall = "halfCauchy", tau_overall_init = 1.0, # Set to "fixed" and tau_overall_init=1 if no overall tau desired
    alpha_global_prior_mean = 0.0, alpha_global_prior_sd = 10.0,
    aleph_prior_mean = 0.0, aleph_prior_sd = 10.0, 
    thin = 1,
    seed = 123,
    verbose = FALSE,
    ping = 1000
) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("BayesLogit", quietly = TRUE)) {
    stop("Package 'BayesLogit' needed for Polya-Gamma sampling. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' needed for mvrnorm. Please install it.", call. = FALSE)
  }
  
  N <- nrow(X_mat)
  p_main_X <- ncol(X_mat)
  
  # --- 1. Construct Design Matrices for Coefficient Groups ---
  # Group 0a: Global Intercept (alpha_global) - handled separately
  X_alpha_global <- matrix(1, nrow = N, ncol = 1)
  
  # Group 0b: Treatment Intercept Modifier (aleph) - handled separately
  X_aleph <- matrix(Z_vec, nrow = N, ncol = 1)
  
  # Group 1: Prognostic effects (beta, beta_int) - for Horseshoe with tau_prognostic
  X_prog_list <- list()
  if (p_main_X > 0) {
    X_prog_list$beta <- X_mat
  }
  main_interaction_indices_X <- create_interaction_pairs_R(p_main_X)
  p_interaction_X <- length(main_interaction_indices_X)
  if (p_interaction_X > 0) {
    X_beta_int_mat <- matrix(0, nrow = N, ncol = p_interaction_X)
    colnames_beta_int <- character(p_interaction_X)
    for (k in 1:p_interaction_X) {
      idx_pair <- main_interaction_indices_X[[k]]
      X_beta_int_mat[, k] <- X_mat[, idx_pair[1]] * X_mat[, idx_pair[2]]
      colnames_beta_int[k] <- paste0(colnames(X_mat)[idx_pair[1]], "_x_", colnames(X_mat)[idx_pair[2]])
    }
    colnames(X_beta_int_mat) <- colnames_beta_int
    X_prog_list$beta_int <- X_beta_int_mat
  }
  X_prognostic <- if(length(X_prog_list) > 0) do.call(cbind, X_prog_list) else matrix(0, nrow=N, ncol=0)
  p_prognostic <- ncol(X_prognostic)
  
  # Group 2: Treatment-related effects (gamma, gamma_int) - for Horseshoe with tau_treatment
  X_treat_list <- list()
  if (p_main_X > 0) {
    X_gamma_mat <- X_mat * Z_vec 
    colnames(X_gamma_mat) <- paste0("Z_", colnames(X_mat))
    X_treat_list$gamma <- X_gamma_mat
  }
  if (p_interaction_X > 0) {
    X_gamma_int_mat <- X_beta_int_mat * Z_vec 
    colnames(X_gamma_int_mat) <- paste0("Z_", colnames(X_beta_int_mat))
    X_treat_list$gamma_int <- X_gamma_int_mat
  }
  X_treatment_related <- if(length(X_treat_list) > 0) do.call(cbind, X_treat_list) else matrix(0, nrow=N, ncol=0)
  p_treatment_related <- ncol(X_treatment_related)
  
  # Combine X_prognostic and X_treatment_related for block update of Horseshoe coeffs
  X_hs <- cbind(X_prognostic, X_treatment_related)
  p_hs_total <- ncol(X_hs) # Total coeffs under Horseshoe
  
  # --- Initialize Parameters ---
  alpha_global <- 0.0
  aleph <- 0.0
  
  Beta_hs <- rep(0, p_hs_total)    
  lambda_hs <- rep(1, p_hs_total)  
  
  current_tau_prognostic <- tau_prognostic_init
  current_tau_treatment <- tau_treatment_init
  current_tau_overall <- tau_overall_init
  
  # Auxiliary variables for global taus (if halfCauchy)
  xi_inv_prog <- 1.0 
  xi_inv_treat <- 1.0
  xi_inv_overall <- 1.0
  
  Sigma2 <- 1.0 # Fixed for logistic regression with Polya-Gamma
  y_adj <- y_vec - 0.5 
  
  # --- MCMC Storage ---
  eff_samp_count <- (n_iter - burn_in) / thin
  alpha_global_samples <- numeric(eff_samp_count)
  aleph_samples <- numeric(eff_samp_count)
  Beta_hs_samples <- if(p_hs_total > 0) matrix(0, nrow = p_hs_total, ncol = eff_samp_count) else matrix(0,0,0)
  lambda_hs_samples <- if(p_hs_total > 0) matrix(0, nrow = p_hs_total, ncol = eff_samp_count) else matrix(0,0,0)
  tau_prognostic_samples <- numeric(eff_samp_count)
  tau_treatment_samples <- numeric(eff_samp_count)
  tau_overall_samples <- numeric(eff_samp_count)
  
  store_idx <- 0
  
  # --- MCMC Loop ---
  if(verbose) message("MCMC (Grouped Horseshoe Logistic) is running...")
  for (iter_idx in 1:n_iter) {
    if (verbose && iter_idx %% ping == 0) {
      message("Iteration: ", iter_idx, "/", n_iter)
    }
    
    # 1. Current linear predictor for Polya-Gamma
    eta_hs_part <- if(p_hs_total > 0) as.vector(X_hs %*% Beta_hs) else rep(0, N)
    eta <- alpha_global + Z_vec * aleph + eta_hs_part
    
    # 2. Sample Polya-Gamma latent variables omega
    omega <- BayesLogit::rpg(num = N, h = 1, z = eta) 
    omega[omega < 1e-9] <- 1e-9 
    z_pg <- (y_adj / omega) # Augmented response: (y_i - 0.5) / omega_i
    
    # 3. Sample global intercept alpha_global
    resid_for_alpha <- z_pg - (Z_vec * aleph + eta_hs_part)
    # Posterior for alpha_global ~ N(mean, var)
    # Data precision sum(omega_i * X_alpha_global_i^2) = sum(omega_i * 1^2) = sum(omega_i)
    # Data mean term sum(omega_i * X_alpha_global_i * resid_for_alpha_i) = sum(omega_i * resid_for_alpha_i)
    prior_prec_alpha <- 1 / safe_var_R(alpha_global_prior_sd^2)
    post_prec_alpha <- sum(omega) + prior_prec_alpha
    post_mean_alpha <- (sum(omega * resid_for_alpha) + alpha_global_prior_mean * prior_prec_alpha) / post_prec_alpha
    alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
    
    # 4. Sample treatment intercept aleph
    resid_for_aleph <- z_pg - (alpha_global + eta_hs_part)
    # Predictor for aleph is Z_vec
    # Data precision sum(omega_i * Z_vec_i^2)
    # Data mean term sum(omega_i * Z_vec_i * resid_for_aleph_i)
    prior_prec_aleph <- 1 / safe_var_R(aleph_prior_sd^2)
    post_prec_aleph <- sum(omega * Z_vec^2) + prior_prec_aleph # Z_vec is 0/1, so Z_vec^2 = Z_vec
    post_mean_aleph <- (sum(omega * Z_vec * resid_for_aleph) + aleph_prior_mean * prior_prec_aleph) / post_prec_aleph
    aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
    
    if (p_hs_total > 0) {
      # 5. Update Horseshoe coefficients Beta_hs and their shrinkage parameters
      # Effective response for Beta_hs part: z_pg - alpha_global - Z_vec * aleph
      z_for_hs <- z_pg - alpha_global - Z_vec * aleph
      
      # Groupings for tau updates
      idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0)
      idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):(p_prognostic + p_treatment_related) else integer(0)
      
      # Update local lambdas (lambda_hs)
      effective_global_taus_for_beta <- rep(0, p_hs_total)
      if(length(idx_prog) > 0) effective_global_taus_for_beta[idx_prog] <- current_tau_prognostic * current_tau_overall
      if(length(idx_treat) > 0) effective_global_taus_for_beta[idx_treat] <- current_tau_treatment * current_tau_overall
      
      temp_lambda_denom_sq <- safe_var_R(effective_global_taus_for_beta^2) # Vectorized
      temp_lambda_scale_term <- (Beta_hs^2) / (2 * temp_lambda_denom_sq) 
      
      nu_lambda_inv <- stats::rexp(p_hs_total, rate = 1 + (1 / safe_var_R(lambda_hs^2)))
      lambda_hs_sq_inv <- stats::rexp(p_hs_total, rate = nu_lambda_inv + temp_lambda_scale_term)
      lambda_hs <- 1 / sqrt(safe_var_R(lambda_hs_sq_inv)) # Get lambda from 1/lambda^2
      lambda_hs[lambda_hs < 1e-6] <- 1e-6
      
      
      # Update current_tau_prognostic
      if (length(idx_prog) > 0) {
        if (method_tau_prognostic == "halfCauchy") {
          xi_inv_prog <- stats::rexp(1, rate = 1 + 1 / safe_var_R(current_tau_prognostic^2))
          rate_prog <- xi_inv_prog + sum(Beta_hs[idx_prog]^2 / (2 * safe_var_R(lambda_hs[idx_prog]^2 * current_tau_overall^2)))
          current_tau_prognostic <- sqrt(rinvgamma_R(1, shape = (length(idx_prog) + 1)/2, scale = rate_prog)) # rate becomes scale for rinvgamma
        } else if (method_tau_prognostic == "truncatedCauchy") {
          # (Implementation from hs, adapted for current_tau_overall)
          tempt_prog <- sum(Beta_hs[idx_prog]^2 / (2 * safe_var_R(lambda_hs[idx_prog]^2 * current_tau_overall^2)))
          et_prog <- 1 / safe_var_R(current_tau_prognostic^2)
          utau_prog <- stats::runif(1, 0, 1 / (1 + et_prog))
          ubt_1_prog <- 1 
          ubt_2_prog <- min((1 - utau_prog)/utau_prog, length(idx_prog)^2) 
          Fubt_1_prog <- stats::pgamma(ubt_1_prog, (length(idx_prog) + 1)/2, rate = tempt_prog)
          Fubt_2_prog <- stats::pgamma(ubt_2_prog, (length(idx_prog) + 1)/2, rate = tempt_prog)
          ut_prog <- stats::runif(1, min(Fubt_1_prog, Fubt_2_prog), max(Fubt_1_prog, Fubt_2_prog))
          if(Fubt_1_prog >= Fubt_2_prog && tempt_prog > 0) { current_tau_prognostic <- sqrt(1/ubt_2_prog)
          } else if (tempt_prog > 1e-9) { # ensure tempt_prog is positive
            et_prog_samp <- stats::qgamma(ut_prog, (length(idx_prog) + 1)/2, rate = tempt_prog)
            current_tau_prognostic <- 1 / sqrt(safe_var_R(et_prog_samp))
          } # else keep old value if tempt is too small
        } # else fixed
        current_tau_prognostic <- max(current_tau_prognostic, 1e-6)
      }
      
      # Update current_tau_treatment
      if (length(idx_treat) > 0) {
        if (method_tau_treatment == "halfCauchy") {
          xi_inv_treat <- stats::rexp(1, rate = 1 + 1 / safe_var_R(current_tau_treatment^2))
          rate_treat <- xi_inv_treat + sum(Beta_hs[idx_treat]^2 / (2 * safe_var_R(lambda_hs[idx_treat]^2 * current_tau_overall^2)))
          current_tau_treatment <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat)) # rate becomes scale
        } else if (method_tau_treatment == "truncatedCauchy") {
          tempt_treat <- sum(Beta_hs[idx_treat]^2 / (2 * safe_var_R(lambda_hs[idx_treat]^2 * current_tau_overall^2)))
          et_treat <- 1 / safe_var_R(current_tau_treatment^2)
          utau_treat <- stats::runif(1, 0, 1 / (1 + et_treat))
          ubt_1_treat <- 1
          ubt_2_treat <- min((1 - utau_treat)/utau_treat, length(idx_treat)^2)
          Fubt_1_treat <- stats::pgamma(ubt_1_treat, (length(idx_treat) + 1)/2, rate = tempt_treat)
          Fubt_2_treat <- stats::pgamma(ubt_2_treat, (length(idx_treat) + 1)/2, rate = tempt_treat)
          ut_treat <- stats::runif(1, min(Fubt_1_treat, Fubt_2_treat), max(Fubt_1_treat, Fubt_2_treat))
          if(Fubt_1_treat >= Fubt_2_treat && tempt_treat > 0) { current_tau_treatment <- sqrt(1/ubt_2_treat)
          } else if (tempt_treat > 1e-9) {
            et_treat_samp <- stats::qgamma(ut_treat, (length(idx_treat) + 1)/2, rate = tempt_treat)
            current_tau_treatment <- 1 / sqrt(safe_var_R(et_treat_samp))
          }
        } # else fixed
        current_tau_treatment <- max(current_tau_treatment, 1e-6)
      }
      
      # Update current_tau_overall
      # Effective component taus (tau_prognostic or tau_treatment for each Beta_hs)
      tau_vector_components_for_overall_tau <- rep(0, p_hs_total)
      if(length(idx_prog) > 0) tau_vector_components_for_overall_tau[idx_prog] <- current_tau_prognostic
      if(length(idx_treat) > 0) tau_vector_components_for_overall_tau[idx_treat] <- current_tau_treatment
      
      if (method_tau_overall == "halfCauchy") {
        xi_inv_overall <- stats::rexp(1, rate = 1 + 1 / safe_var_R(current_tau_overall^2))
        rate_overall <- xi_inv_overall + sum(Beta_hs^2 / (2 * safe_var_R(lambda_hs^2 * tau_vector_components_for_overall_tau^2)))
        current_tau_overall <- sqrt(rinvgamma_R(1, shape = (p_hs_total + 1)/2, scale = rate_overall)) # rate becomes scale
      } else if (method_tau_overall == "truncatedCauchy") {
        tempt_overall <- sum(Beta_hs^2 / (2 * safe_var_R(lambda_hs^2 * tau_vector_components_for_overall_tau^2)))
        et_overall <- 1 / safe_var_R(current_tau_overall^2)
        utau_overall <- stats::runif(1, 0, 1 / (1 + et_overall))
        ubt_1_overall <- 1
        ubt_2_overall <- min((1 - utau_overall)/utau_overall, p_hs_total^2)
        Fubt_1_overall <- stats::pgamma(ubt_1_overall, (p_hs_total + 1)/2, rate = tempt_overall)
        Fubt_2_overall <- stats::pgamma(ubt_2_overall, (p_hs_total + 1)/2, rate = tempt_overall)
        ut_overall <- stats::runif(1, min(Fubt_1_overall,Fubt_2_overall), max(Fubt_1_overall,Fubt_2_overall))
        if(Fubt_1_overall >= Fubt_2_overall && tempt_overall >0) { current_tau_overall <- sqrt(1/ubt_2_overall)
        } else if (tempt_overall > 1e-9){
          et_overall_samp <- stats::qgamma(ut_overall, (p_hs_total + 1)/2, rate = tempt_overall)
          current_tau_overall <- 1 / sqrt(safe_var_R(et_overall_samp))
        }
      } # else fixed
      current_tau_overall <- max(current_tau_overall, 1e-6)
      
      # Update Beta_hs (all Horseshoe coefficients)
      prior_var_beta_hs <- safe_var_R(lambda_hs^2 * effective_global_taus_for_beta^2) # Sigma2=1
      prior_prec_beta_hs_diag <- diag(1 / prior_var_beta_hs, p_hs_total, p_hs_total)
      
      XtOmegaX_hs <- t(X_hs) %*% diag(omega) %*% X_hs
      post_prec_beta_hs <- XtOmegaX_hs + prior_prec_beta_hs_diag
      
      # Prior mean for Beta_hs is 0
      rhs_mean_beta_hs <- t(X_hs) %*% diag(omega) %*% z_for_hs 
      
      chol_beta_hs <- try(chol(post_prec_beta_hs), silent = TRUE)
      if (inherits(chol_beta_hs, "try-error")) {
        jitter_val_beta <- 1e-6 * mean(diag(post_prec_beta_hs))
        post_prec_beta_hs_jit <- post_prec_beta_hs + diag(jitter_val_beta, p_hs_total)
        chol_beta_hs <- try(chol(post_prec_beta_hs_jit))
        if (inherits(chol_beta_hs, "try-error")){
          warning(paste("Cholesky for Beta_hs failed at iter", iter_idx, "- skipping Beta_hs update this iter."))
          # Beta_hs remains from previous iteration
        } else {
          post_cov_beta_hs <- chol2inv(chol_beta_hs)
          post_mean_beta_hs <- post_cov_beta_hs %*% rhs_mean_beta_hs
          Beta_hs <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean_beta_hs), Sigma = post_cov_beta_hs)
        }
      } else {
        post_cov_beta_hs <- chol2inv(chol_beta_hs)
        post_mean_beta_hs <- post_cov_beta_hs %*% rhs_mean_beta_hs
        Beta_hs <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean_beta_hs), Sigma = post_cov_beta_hs)
      }
    } # end if (p_hs_total > 0)
    
    # Store samples
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
    }
  } # End MCMC loop
  
  # --- Unpack Beta_hs_samples into named components ---
  beta_out       <- if(p_prognostic > 0 && "beta" %in% names(X_prog_list)) Beta_hs_samples[colnames(X_prog_list$beta), , drop=FALSE] else matrix(0,0,eff_samp_count)
  beta_int_out   <- if(p_prognostic > 0 && "beta_int" %in% names(X_prog_list)) Beta_hs_samples[colnames(X_prog_list$beta_int), , drop=FALSE] else matrix(0,0,eff_samp_count)
  gamma_out      <- if(p_treatment_related > 0 && "gamma" %in% names(X_treat_list)) Beta_hs_samples[colnames(X_treat_list$gamma), , drop=FALSE] else matrix(0,0,eff_samp_count)
  gamma_int_out  <- if(p_treatment_related > 0 && "gamma_int" %in% names(X_treat_list)) Beta_hs_samples[colnames(X_treat_list$gamma_int), , drop=FALSE] else matrix(0,0,eff_samp_count)
  
  # For lambda_hs, need to map back to original structure too
  lambda_beta_out      <- if(p_prognostic > 0 && "beta" %in% names(X_prog_list)) lambda_hs_samples[colnames(X_prog_list$beta), , drop=FALSE] else matrix(0,0,eff_samp_count)
  lambda_beta_int_out  <- if(p_prognostic > 0 && "beta_int" %in% names(X_prog_list)) lambda_hs_samples[colnames(X_prog_list$beta_int), , drop=FALSE] else matrix(0,0,eff_samp_count)
  lambda_gamma_out     <- if(p_treatment_related > 0 && "gamma" %in% names(X_treat_list)) lambda_hs_samples[colnames(X_treat_list$gamma), , drop=FALSE] else matrix(0,0,eff_samp_count)
  lambda_gamma_int_out <- if(p_treatment_related > 0 && "gamma_int" %in% names(X_treat_list)) lambda_hs_samples[colnames(X_treat_list$gamma_int), , drop=FALSE] else matrix(0,0,eff_samp_count)
  
  return(list(
    alpha = alpha_global_samples,
    beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count, 0),
    beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count, 0),
    aleph = aleph_samples,
    gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count, 0),
    gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count, 0),
    
    lambda_beta = if(nrow(lambda_beta_out) > 0) t(lambda_beta_out) else matrix(NA, eff_samp_count, 0),
    lambda_beta_int = if(nrow(lambda_beta_int_out) > 0) t(lambda_beta_int_out) else matrix(NA, eff_samp_count, 0),
    lambda_gamma = if(nrow(lambda_gamma_out) > 0) t(lambda_gamma_out) else matrix(NA, eff_samp_count, 0),
    lambda_gamma_int = if(nrow(lambda_gamma_int_out) > 0) t(lambda_gamma_int_out) else matrix(NA, eff_samp_count, 0),
    
    tau_prognostic = tau_prognostic_samples,
    tau_treatment = tau_treatment_samples,
    tau_overall = tau_overall_samples
  ))
}

# --- Example Usage (using the simple simulation) ---
source('R/simul_1.R')
sim_data_check <- simulate_logit_direct(n=1500)
y_sim <- sim_data_check$Y
X_sim <- sim_data_check$X_mat
Z_sim <- sim_data_check$T
true_params_check <- sim_data_check$true_params

cat("\n--- Running Grouped Horseshoe Logistic Sampler ---\n")
fit_grouped_hs <- fit_grouped_horseshoes_logistic_R(
  y_vec = y_sim,
  X_mat = X_sim,
  Z_vec = Z_sim,
  n_iter = 5000, # Increase for better results
  burn_in = 2000,
  # p_lin_group1 will be calculated internally, or you can set it:
  # num_prog_main = ncol(X_sim), num_treat_main = ncol(X_sim)
  # p_lin_group1 = num_prog_main, # This is for the prognostic group
  method_tau_prognostic = "halfCauchy", tau_prognostic_init = 0.1,
  method_tau_treatment = "halfCauchy", tau_treatment_init = 0.1,
  method_tau_overall = "Fixed", tau_overall_init = 1,
  alpha_global_prior_sd = 5.0,
  aleph_prior_sd = 5.0,
  thin = 1,
  seed = 103,
  verbose = TRUE,
  ping = 500
)

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