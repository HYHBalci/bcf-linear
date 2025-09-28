library(Rcpp)
sourceCpp("src/horseshoe_samplers.cpp")
set.seed(123)
library(MASS) 

rinvgamma <- function(shape, scale){
  if (shape <= 0.0 || scale <= 0.0) {
    stop("Shape and scale must be positive.");
  }
  
  g = rgamma(1, shape, scale)
  out = 1.0 / g
  return(out)
}

loglikeTauInt <- function(
    tau_int,
    beta_int_base,
    int_pairs_base,        
    tau_main,
    sigma,
    include_treatment_int = FALSE,
    beta_int_trt = numeric(0),
    tau_trt      = numeric(0),   # Not used directly here, but included for consistency
    int_pairs_trt = data.frame()
) {
  # 1) Check Uniform(0.01, 1.0) boundary
  if (tau_int < 0.01 || tau_int > 1.0) {
    return(-Inf)
  }
  
  logp    <- 0.0
  log2pi  <- log(2.0 * pi)
  
  # 2) Baseline interaction terms
  # var_ij = tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
  # log N(0, var_ij) = -0.5 * [ log(2π var_ij) + beta^2 / var_ij ]
  for (k in seq_len(nrow(int_pairs_base))) {
    iVar <- int_pairs_base[k, 1]  # or int_pairs_base$Var1[k], if it's a data.frame with named columns
    jVar <- int_pairs_base[k, 2]  # or int_pairs_base$Var2[k]
    
    var_ij <- tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
    beta2  <- beta_int_base[k]^2
    
    logp <- logp -
      0.5 * (log2pi + log(var_ij)) -
      0.5 * (beta2 / var_ij)
  }
  
  # 3) Treatment interaction terms (if included)
  if (include_treatment_int && length(beta_int_trt) > 0) {
    for (k in seq_len(nrow(int_pairs_trt))) {
      iVar <- int_pairs_trt[k, 1]
      jVar <- int_pairs_trt[k, 2]
      
      var_ij <- tau_int * tau_main[iVar] * tau_main[jVar] * (sigma^2)
      beta2  <- beta_int_trt[k]^2
      
      logp <- logp -
        0.5 * (log2pi + log(var_ij)) -
        0.5 * (beta2 / var_ij)
    }
  }
  
  return(logp)
}

sample_linear_part <- function(y, z, X, intTreat = TRUE, iter = 1000, burnin = 500, horseshoe = FALSE) {
  n <- length(y)
  p <- ncol(X)
  
  # Generate interaction terms 
  int_pairs <- expand.grid(1:p, 1:p)
  int_pairs <- int_pairs[int_pairs$Var1 < int_pairs$Var2, ]  # Avoid duplicates
  p_int <- nrow(int_pairs)
  
  X_int <- do.call(
    cbind,
    lapply(seq_len(nrow(int_pairs)), function(k) {
      X[, int_pairs$Var1[k]] * X[, int_pairs$Var2[k]]
    })
  )
  
  # Initialize parameters
  alpha <- 0.0
  sigma <- 1
  tau_glob <- 1.0
  beta <- rep(0.0, p)
  if(horseshoe){
    tau <- rep(1.0, p + p_int)
  } else{ 
    tau <- rep(1.0, p)
  }
  beta_int <- rep(0.0, p_int)
  tau_int <- 0.5
  
  nu <- 1 / (tau^2)
  xi <- 0.01
  
  # Storage for posterior samples
  alpha_samples <- numeric(iter)
  xi_samples <- numeric(iter)
  tau_glob_samples <- numeric(iter)
  beta_samples <- matrix(0, iter, p)
  if(!horseshoe){
    tau_samples <- matrix(0, iter, p)
  } else {
    tau_samples <- matrix(0, iter, p + p_int)
  }
  nu_samples <- matrix(0, iter, p + p_int)
  beta_int_samples <- matrix(0, iter, p_int)
  tau_int_samples <- numeric(iter)
  sigma_samples <- numeric(iter)
  
  for (i in 1:(burnin+iter)) {
    # Sample alpha (intercept)
    r_alpha <- y - (z * alpha + X %*% beta + X_int %*% beta_int)
    alpha <- sample_alpha(n, r_alpha, z, sigma, alpha_prior_sd = 10.0)
    
    # Sample main effects (beta)
    for (j in 1:p) {
      r_beta <- y - (z * alpha + X[, -j] %*% beta[-j] + X_int %*% beta_int)
      beta[j] <- sample_beta_j(n, r_beta, z, X[, j], 1 / (tau[j]*tau_glob), sigma)
      if(horseshoe){
        tau[j] <- sqrt(rinvgamma(1, (1 / nu[j]) + (beta[j] * beta[j]) / (2*tau_glob*tau_glob*sigma*sigma)));
      } else{
        tau[j] <- sample_tau_j_slice(tau[j], beta[j], j, beta_int = beta_int, tau = tau, tau_int = tau_int, sigma = sigma)
      }
    }
    
    # Sample interaction effects (if enabled)
    if (intTreat) {
      for (k in 1:p_int) {
        iVar <- int_pairs$Var1[k]
        jVar <- int_pairs$Var2[k]
        r_beta_int <- y - (z * alpha + X %*% beta + X_int[, -k] %*% beta_int[-k])
        if(horseshoe){
          beta_int[k] <- sample_beta_j(n, r_beta_int, z, X[, iVar] * X[, jVar], tau[p+k]*tau_glob, sigma)
        } else {
          beta_int[k] <- sample_beta_j(n, r_beta_int, z, X[, iVar] * X[, jVar], sqrt(tau_int * tau[iVar] * tau[jVar]), sigma)
        }
        if(horseshoe){
          for(o in 1:(p + p_int)){
            nu[o] <- rinvgamma(1, 1 + 1 / (tau[o]*tau[o]))
          } 
          tau[p+k] <- sqrt(rinvgamma(1, (1 / nu[p + k]) + (beta_int[k] * beta_int[k]) / (2*tau_glob*tau_glob*sigma*sigma)));
        }
      }
      
      if(!horseshoe) {
      # Sample tau_int (global shrinkage for interactions)
      currentTauInt <- tau_int
      proposedTauInt <- runif(1, 0.01, 1.0)
      
      logPosteriorCurrent <- loglikeTauInt(currentTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      logPosteriorProposed <- loglikeTauInt(proposedTauInt, beta_int, int_pairs, tau, sigma, FALSE, numeric(0), numeric(0), list())
      logAccept <- logPosteriorProposed - logPosteriorCurrent
      if (log(runif(1)) < logAccept) tau_int <- proposedTauInt
      } else {
        sum = 0
        for(j in 1:(p + p_int)){
          if(j <= p){
            sum = sum +  beta[j]^2 / tau[j]^2
          } else{
            sum = sum + beta_int[j-p]^2 / tau[j]^2
          }
        }
        xi <- rinvgamma(1, 1 + 1 / (tau_glob*tau_glob))
        tau_glob <- sqrt(rinvgamma((p + p_int +1)/2, 1/xi +  (1/(2*sigma*sigma))*sum))
        
      }
    }
    
    # Sample sigma
    resid <- y - (z * alpha + X %*% beta + X_int %*% beta_int)
    
    if(horseshoe){
      lambda_sq <- tau^2  
      Lambda_inv <- diag(1 / lambda_sq)
      beta_tot <- c(beta, beta_int)
      scale <- as.numeric(0.5 * crossprod(resid) + 0.5 * t(beta_tot) %*% Lambda_inv %*% beta_tot)
      
      # Shape parameter
      shape <- (n + p + p_int) / 2
      
      # Sample from inverse-gamma
      sigma <- sqrt(rinvgamma(shape, scale))
    } else {
      sigma <- sqrt(sample_sigma2_ig(n, resid,1.0, 0.001))
    }
    
    # Store samples
    if (i > burnin) {
      idx <- i - burnin
      alpha_samples[idx] <- alpha
      beta_samples[idx, ] <- beta
      tau_glob_samples[idx] <- tau_glob
      nu_samples[idx, ] <- nu
      xi_samples[idx] <- xi 
      print(tau)
      tau_samples[idx, ] <- tau
      beta_int_samples[idx, ] <- beta_int
      tau_int_samples[idx] <- tau_int
      sigma_samples[idx] <- sigma
    }
  }
  
  return(list(alpha = alpha_samples, beta = beta_samples, tau = tau_samples,
              beta_int = beta_int_samples, tau_int = tau_int_samples, sigma = sigma_samples, xi = xi_samples, tau_glob = tau_glob_samples, nu = nu_samples))
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
#' @param prognostic_shrinkage Shrinkage method for mu(x) part: "horseshoe" or "linked_shrinkage".
#' @param treatment_shrinkage Shrinkage method for tau(x) part: "horseshoe" or "linked_shrinkage".
#' @param num_chains Integer, the number of independent MCMC chains to run.
#' @param propensity_as_covariate Boolean, if TRUE, propensity scores are estimated and added as a covariate to the prognostic model `mu(x)`.
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
    prognostic_shrinkage = c("horseshoe", "linked_shrinkage"),
    treatment_shrinkage = c("horseshoe", "linked_shrinkage"),
    num_chains = 1, propensity_as_covariate = FALSE,
    thin = 1, seed = 123, verbose = FALSE, ping = 1000
) {
  # --- Setup and Validation ---
  if (!is.null(seed)) { set.seed(seed); chain_seeds <- sample.int(1e6, num_chains) } else { chain_seeds <- rep(NA, num_chains) }
  family <- match.arg(family)
  prognostic_shrinkage <- match.arg(prognostic_shrinkage)
  treatment_shrinkage <- match.arg(treatment_shrinkage)
  if (length(unique(na.omit(y_vec))) <= 1) stop("Outcome variable `y_vec` is constant.")
  
  N <- nrow(X_mat)
  p_main_X_orig <- ncol(X_mat)
  X_mat_orig <- X_mat
  
  # --- Propensity Score Estimation ---
  if (propensity_as_covariate) {
    if (verbose) message("Estimating propensity scores to use as a covariate...")
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
  
  main_interaction_indices_prog <- create_interaction_pairs_R(p_main_X_prog)
  p_interaction_prog <- if(is.matrix(main_interaction_indices_prog)) nrow(main_interaction_indices_prog) else 0
  X_prognostic <- X_mat_for_prog
  if (p_interaction_prog > 0) {
    X_beta_int_mat <- do.call(cbind, lapply(1:p_interaction_prog, function(k) X_mat_for_prog[, main_interaction_indices_prog[k,1]] * X_mat_for_prog[, main_interaction_indices_prog[k,2]]))
    colnames(X_beta_int_mat) <- sapply(1:p_interaction_prog, function(k) paste0("beta_int_", paste(colnames(X_mat_for_prog)[main_interaction_indices_prog[k,]], collapse="_")))
    X_prognostic <- cbind(X_prognostic, X_beta_int_mat)
  }
  p_prognostic <- ncol(X_prognostic)
  
  main_interaction_indices_orig <- create_interaction_pairs_R(p_main_X_orig)
  p_interaction_orig <- if(is.matrix(main_interaction_indices_orig)) nrow(main_interaction_indices_orig) else 0
  X_treatment_related <- X_mat_orig * Z_vec
  colnames(X_treatment_related) <- paste0("Z_", colnames(X_mat_orig))
  if (p_interaction_orig > 0) {
    base_int_treat <- do.call(cbind, lapply(1:p_interaction_orig, function(k) (X_mat_orig[,main_interaction_indices_orig[k,1]]*X_mat_orig[,main_interaction_indices_orig[k,2]])))
    X_gamma_int_mat <- base_int_treat * Z_vec
    colnames(X_gamma_int_mat) <- sapply(1:p_interaction_orig, function(k) paste0("Z_int_", paste(colnames(X_mat_orig)[main_interaction_indices_orig[k,]], collapse="_")))
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
    alpha_global <- 0.0
    aleph <- 0.0
    Beta_all <- rep(0, p_total)
    
    # Prognostic shrinkage params
    lambda_prog <- rep(1, p_prognostic)
    tau_prog <- rep(1, p_main_X_prog)
    tau_int_prog <- 0.5
    tau_glob_prog <- 1.0
    
    # Treatment shrinkage params
    lambda_treat <- rep(1, p_treatment_related)
    tau_treat <- rep(1, p_main_X_orig)
    tau_int_treat <- 0.5
    tau_glob_treat <- 1.0
    
    if (family == "gaussian") { 
      sigma_sq <- var(y_vec) 
      } 
    else {
      y_adj <- y_vec - 0.5
      }
    
    # Storage for this chain
    eff_samp_count <- floor((n_iter - burn_in) / thin)
    alpha_global_samples <- numeric(eff_samp_count); aleph_samples <- numeric(eff_samp_count)
    Beta_all_samples <- if(p_total > 0) matrix(0, nrow = p_total, ncol = eff_samp_count) else matrix(0,0,0)
    
    # Storage for shrinkage parameters
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
      resid_for_alpha <- z_aug - (Z_vec * aleph + eta_hs_part)
      post_prec_alpha <- sum(w_vec) / scale_factor + (1 / safe_var_R(alpha_global_prior_sd^2))
      post_mean_alpha <- (sum(w_vec * resid_for_alpha) / scale_factor) / post_prec_alpha
      alpha_global <- rnorm(1, mean = post_mean_alpha, sd = sqrt(1 / post_prec_alpha))
      
      resid_for_aleph <- z_aug - (alpha_global + eta_hs_part)
      post_prec_aleph <- sum(w_vec * Z_vec^2) / scale_factor + (1 / safe_var_R(aleph_prior_sd^2))
      post_mean_aleph <- (sum(w_vec * Z_vec * resid_for_aleph) / scale_factor) / post_prec_aleph
      aleph <- rnorm(1, mean = post_mean_aleph, sd = sqrt(1 / post_prec_aleph))
      
      if (p_total > 0) {
        z_for_betas <- z_aug - alpha_global - Z_vec * aleph
        
        idx_prog <- if(p_prognostic > 0) 1:p_prognostic else integer(0)
        idx_treat <- if(p_treatment_related > 0) (p_prognostic + 1):p_total else integer(0)
        
        # --- Update Prognostic Shrinkage ---
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
        } else { 
          if(p_prognostic > 0) {
            beta_prog_main <- Beta_all[1:p_main_X_prog]
            beta_prog_int <- Beta_all[(p_main_X_prog + 1):p_prognostic]
            for(j in 1:p_main_X_prog) {
              tau_prog[j] <- sample_tau_j_slice(tau_prog[j], beta_prog_main[j], j, beta_int = beta_prog_int, tau = tau_prog, tau_int = 1, sigma = sqrt(sigma_sq))
            }
              prop_tau_int <- 1 
          }
        }
        
        # --- Update Treatment Shrinkage ---
        if (treatment_shrinkage == "horseshoe") {
          if(length(idx_treat) > 0) {
            beta_treat <- Beta_all[idx_treat]; temp_lambda_scale_term <- (beta_treat^2) / (2 * scale_factor * safe_var_R(tau_glob_treat^2))
            nu_lambda_inv <- rexp(length(idx_treat), rate = 1 + (1 / safe_var_R(lambda_treat^2)))
            lambda_treat_sq_inv <- rexp(length(idx_treat), rate = nu_lambda_inv + temp_lambda_scale_term)
            lambda_treat <- 1 / sqrt(safe_var_R(lambda_treat_sq_inv))
            xi_inv_treat <- rexp(1, rate = 1 + 1/safe_var_R(tau_glob_treat^2))
            rate_treat <- xi_inv_treat + sum(beta_treat^2 / (2 * scale_factor * safe_var_R(lambda_treat^2)))
            tau_glob_treat <- sqrt(rinvgamma_R(1, shape = (length(idx_treat) + 1)/2, scale = rate_treat))
          }
        } else { # linked_shrinkage
          if(p_treatment_related > 0) {
            beta_treat_main <- Beta_all[idx_treat[1:p_main_X_orig]]
            beta_treat_int <- Beta_all[idx_treat[(p_main_X_orig + 1):p_treatment_related]]
            for(j in 1:p_main_X_orig) {
              tau_treat[j] <- sample_tau_j_slice(tau_treat[j], beta_treat_main[j], j, beta_int = beta_treat_int, tau = tau_prog, tau_int = 1, sigma = sqrt(sigma_sq))
            }
            if(p_interaction_orig > 0){
            tau_int_treat <-1
            }
          }
        }
        
        # --- Construct Prior Precision and Sample Beta_all ---
        prior_var_beta <- numeric(p_total)
        if (prognostic_shrinkage == "horseshoe") { if(length(idx_prog) > 0) prior_var_beta[idx_prog] <- lambda_prog^2 * tau_glob_prog^2
        } else { if(p_prognostic > 0) { prog_main_indices <- 1:p_main_X_prog
        prior_var_beta[idx_prog[prog_main_indices]] <- tau_prog^2
        if(p_interaction_prog > 0){ prog_int_indices <- (p_main_X_prog+1)
        p_prognostic
        for(k in 1:p_interaction_prog) 
        { pair <- main_interaction_indices_prog[k,]
          prior_var_beta[idx_prog[prog_int_indices[k]]] <- tau_int_prog^2 * tau_prog[pair[1]]^2 * tau_prog[pair[2]]^2 }}}
          }
        if (treatment_shrinkage == "horseshoe") { 
          if(length(idx_treat) > 0) prior_var_beta[idx_treat] <- lambda_treat^2 * tau_glob_treat^2
        } else { if(p_treatment_related > 0) {
          treat_main_indices <- 1:p_main_X_orig; prior_var_beta[idx_treat[treat_main_indices]] <- tau_treat^2
          if(p_interaction_orig > 0){ 
            treat_int_indices <- (p_main_X_orig+1):p_treatment_related
            for(k in 1:p_interaction_orig) { pair <- main_interaction_indices_orig[k,]
            prior_var_beta[idx_treat[treat_int_indices[k]]] <- tau_int_treat^2 * tau_treat[pair[1]]^2 * tau_treat[pair[2]]^2 
            }}}}
        
        prior_prec_diag <- diag(1 / safe_var_R(prior_var_beta), p_total, p_total)
        XtWX <- t(X_full) %*% (w_vec * X_full)
        post_prec_unscaled <- XtWX + prior_prec_diag
        
        chol_attempt <- try(chol(post_prec_unscaled), silent = TRUE)
        if (inherits(chol_attempt, "try-error")) { 
          warning(paste("Cholesky for Beta_all failed at iter", iter_idx, "- skipping update."))
        } else { 
          post_cov_unscaled_inv <- chol2inv(chol_attempt)
          rhs_mean <- t(X_full) %*% (w_vec * z_for_betas)
          post_mean <- post_cov_unscaled_inv %*% rhs_mean
          Beta_all <- MASS::mvrnorm(n = 1, mu = as.vector(post_mean), Sigma = scale_factor * post_cov_unscaled_inv) }
      }
      
      # --- Sample sigma_sq for Gaussian model ---
      if (family == "gaussian") { 
        eta_new <- alpha_global + Z_vec * aleph + as.vector(X_full %*% Beta_all)
        residuals_for_sigma <- y_vec - eta_new
        shape_sigma_post <- 0.001 + N / 2
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
        else { if(p_main_X_prog > 0) tau_prog_samples[, store_idx] <- tau_prog; tau_int_prog_samples[store_idx] <- tau_int_prog }
        if (treatment_shrinkage == "horseshoe") { tau_treat_samples[store_idx] <- tau_glob_treat; if(p_treatment_related > 0) lambda_treat_samples[, store_idx] <- lambda_treat }
        else { if(p_main_X_orig > 0) tau_treat_samples[, store_idx] <- tau_treat; tau_int_treat_samples[store_idx] <- tau_int_treat }
      }
    } 
    
    # Store the results of this chain
    all_chains_results[[chain_num]] <- list(alpha=alpha_global_samples, aleph=aleph_samples, Beta_all=Beta_all_samples, sigma_sq=if(family=="gaussian") sigma_sq_samples else NULL,
                                            tau_prog=tau_prog_samples, tau_int_prog=tau_int_prog_samples, lambda_prog=lambda_prog_samples,
                                            tau_treat=tau_treat_samples, tau_int_treat=tau_int_treat_samples, lambda_treat=lambda_treat_samples)
  } 
  
  # --- Combine and Unpack Results ---
  if(verbose) message("Combining samples from all chains...")
  
  combined_alpha_samples <- do.call(c, lapply(all_chains_results, `[[`, "alpha"))
  combined_aleph_samples <- do.call(c, lapply(all_chains_results, `[[`, "aleph"))
  if(p_total > 0) combined_Beta_all_samples <- do.call(cbind, lapply(all_chains_results, `[[`, "Beta_all")) else combined_Beta_all_samples <- matrix(0,0,0)
  
  # Combine shrinkage parameters
  if(prognostic_shrinkage == "horseshoe") { combined_tau_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_prog")); if(p_prognostic>0) combined_lambda_prog <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_prog")) }
  else { combined_tau_prog <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_prog))); combined_tau_int_prog <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_prog")) }
  if(treatment_shrinkage == "horseshoe") { combined_tau_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_treat")); if(p_treatment_related>0) combined_lambda_treat <- do.call(cbind, lapply(all_chains_results, `[[`, "lambda_treat")) }
  else { combined_tau_treat <- do.call(rbind, lapply(all_chains_results, function(x) t(x$tau_treat))); combined_tau_int_treat <- do.call(c, lapply(all_chains_results, `[[`, "tau_int_treat")) }
  
  if (family == "gaussian") { combined_sigma_sq_samples <- do.call(c, lapply(all_chains_results, `[[`, "sigma_sq")) }
  
  # Unpack Beta_all_samples
  idx_current <- 0
  idx_beta <- if(p_prognostic > 0 && "beta" %in% names(X_prog_list)) (idx_current + 1):(idx_current + p_main_X_prog) else integer(0); idx_current <- idx_current + p_main_X_prog
  idx_beta_int <- if(p_prognostic > 0 && "beta_int" %in% names(X_prog_list)) (idx_current + 1):(idx_current + p_interaction_prog) else integer(0); idx_current <- idx_current + p_interaction_prog
  idx_gamma <- if(p_treatment_related > 0 && "gamma" %in% names(X_treat_list)) (idx_current + 1):(idx_current + p_main_X_orig) else integer(0); idx_current <- idx_current + p_main_X_orig
  idx_gamma_int <- if(p_treatment_related > 0 && "gamma_int" %in% names(X_treat_list)) (idx_current + 1):(idx_current + p_interaction_orig) else integer(0)
  
  beta_out <- if(length(idx_beta) > 0) combined_Beta_all_samples[idx_beta, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  beta_int_out <- if(length(idx_beta_int) > 0) combined_Beta_all_samples[idx_beta_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_out <- if(length(idx_gamma) > 0) combined_Beta_all_samples[idx_gamma, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  gamma_int_out <- if(length(idx_gamma_int) > 0) combined_Beta_all_samples[idx_gamma_int, , drop=FALSE] else matrix(0,0,eff_samp_count * num_chains)
  
  if(nrow(beta_out) > 0) rownames(beta_out) <- colnames(X_prog_list$beta)
  if(nrow(beta_int_out) > 0) rownames(beta_int_out) <- colnames(X_prog_list$beta_int)
  if(nrow(gamma_out) > 0) rownames(gamma_out) <- colnames(X_treat_list$gamma)
  if(nrow(gamma_int_out) > 0) rownames(gamma_int_out) <- colnames(X_treat_list$gamma_int)
  
  output <- list(
    alpha = combined_alpha_samples,
    beta = if(nrow(beta_out) > 0) t(beta_out) else matrix(NA, eff_samp_count * num_chains, 0),
    beta_interaction = if(nrow(beta_int_out) > 0) t(beta_int_out) else matrix(NA, eff_samp_count * num_chains, 0),
    aleph = combined_aleph_samples,
    gamma = if(nrow(gamma_out) > 0) t(gamma_out) else matrix(NA, eff_samp_count * num_chains, 0),
    gamma_int = if(nrow(gamma_int_out) > 0) t(gamma_int_out) else matrix(NA, eff_samp_count * num_chains, 0)
  )
  
  if (prognostic_shrinkage == "horseshoe") { output$tau_prognostic_global <- combined_tau_prog; if(p_prognostic>0) output$lambda_prognostic <- t(combined_lambda_prog) }
  else { output$tau_prognostic_local <- t(combined_tau_prog); output$tau_prognostic_interaction <- combined_tau_int_prog }
  if (treatment_shrinkage == "horseshoe") { output$tau_treatment_global <- combined_tau_treat; if(p_treatment_related>0) output$lambda_treatment <- t(combined_lambda_treat) }
  else { output$tau_treatment_local <- t(combined_tau_treat); output$tau_treatment_interaction <- combined_tau_int_treat }
  
  if (family == "gaussian") {
    output$sigma_sq = combined_sigma_sq_samples
  }
  
  return(output)
}

