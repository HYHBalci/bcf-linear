data <- read.csv("C:/Users/P094412/Documents/PhD project/Data/data_Japan_clean.csv")

library(stochtree)
library(dplyr)

Y <- as.vector(data$ponv)
Z <- as.vector(data$dexamethasone_use)
X_initial <- data[, !(names(data) %in% c("ponv", "dexamethasone_use"))]
general_params_default <- list(
  cutpoint_grid_size = 100, standardize = TRUE, 
  sample_sigma2_global = FALSE, sigma2_global_init = 1, 
  sigma2_global_shape = 1, sigma2_global_scale = 0.001,
  variable_weights = NULL, propensity_covariate = "mu", 
  adaptive_coding = FALSE, control_coding_init = -0.5, 
  treated_coding_init = 0.5, rfx_prior_var = NULL, 
  random_seed = 1234, keep_burnin = FALSE, keep_gfr = FALSE, 
  keep_every = 1, num_chains = 1, verbose = T, 
  global_shrinkage = T, unlink = T, propensity_seperate = F, gibbs =T, step_out = 0.5, max_steps = 50, save_output = F, probit_outcome_model = T
)

nbcf_fit <- bcf_linear_probit(
  X_train = X_initial,
  y_train = Y,
  Z_train = Z, 
  num_gfr = 25, 
  num_burnin = 1000, 
  num_mcmc = 6000, 
  general_params = general_params_default
)
library(coda)
chain_beta <- nbcf_fit$alpha[1,]
summary(chain_beta)

output <- standardize_X_by_index(X_initial)
X <- output$X_final
is_continuous_map_final_X <- output$X_final_var_info$is_continuous
load('my_multiple_objects_2.RData')

num_chains <- dim(nbcf_fit$Beta)[1]
alpha_samples <- as.vector(t(nbcf_fit$alpha))
beta_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta[chain, , ]))
beta_int_samples <- do.call(rbind, lapply(1:num_chains, function(chain) nbcf_fit$Beta_int[chain, , ]))

p_mod <- dim(nbcf_fit$Beta)[3]
  
interaction_cols_list <- list()
int_pairs <- list()
col_idx_counter <- 1


for (i in 1:(p_mod - 1)) {
    for (j in (i + 1):p_mod) {
      if (is_continuous_map_final_X[i] || is_continuous_map_final_X[j]) {
        interaction_cols_list[[col_idx_counter]] <- X[, i] * X[, j]
        int_pairs[[col_idx_counter]] <- c(i, j) # Store original X_final indices for naming
        col_idx_counter <- col_idx_counter + 1
      }
    }
}
ipairs <- as.matrix(int_pairs)
n_obser <- nrow(X)
# Posterior samples of tau(x)
tau_posterior <- matrix(rep(alpha_samples, each = n_obser),
                        nrow = n_obser, byrow = FALSE) +
  X %*% t(beta_samples)

for (idx in 1:ncol(ipairs)) {
  j <- ipairs[[idx]][1]
  k <- ipairs[[idx]][2]
  tau_posterior <- tau_posterior +
    (X[, j] * X[, k]) %*% t(beta_int_samples[, idx])
}

y_hat_posterior <- tau_posterior + nbcf_fit$mu_hat_train
hist(y_hat_posterior)

mu_posterior_samples <- nbcf_fit$mu_hat_train

tau_posterior_samples <- tau_posterior

# --- Step 1: Calculate Posterior Probabilities for Control and Treatment ---

prob_control_posterior <- pnorm(mu_posterior_samples)

prob_treated_posterior <- pnorm(mu_posterior_samples + tau_posterior_samples)

# --- Step 2: Calculate the Posterior of the Individual Shift in Risk ---

# The risk difference is simply the difference between the two probabilities
risk_shift_posterior <- prob_treated_posterior - prob_control_posterior

# --- Step 3: Summarize and Interpret the Results ---

# Calculate the posterior mean risk shift for each individual (a point estimate)
individual_mean_risk_shift <- rowMeans(risk_shift_posterior)

# Calculate a 95% credible interval for the risk shift for each individual
individual_ci_risk_shift <- t(apply(risk_shift_posterior, 1, quantile, probs = c(0.025, 0.975)))
colnames(individual_ci_risk_shift) <- c("CI_Lower_2.5%", "CI_Upper_97.5%")

# Combine into a results data.frame
results_df <- data.frame(
  SubjectID = 1:nrow(X), # Assuming X is in the same order
  Mean_Risk_Shift = individual_mean_risk_shift,
  CI_Lower = individual_ci_risk_shift[,1],
  CI_Upper = individual_ci_risk_shift[,2]
)

# --- View and Use the Results ---
cat("--- Individual Shift in Risk (Posterior Summary) ---\n")
print(head(results_df))

# How to interpret for a specific individual (e.g., the first one):
first_person_effect <- results_df[1, ]
cat("\n--- Example Interpretation for First Individual ---\n")
print(first_person_effect)

if (first_person_effect$CI_Lower > 0) {
  cat("Result: There is strong evidence (>97.5% probability) that the treatment INCREASES the risk for this individual.\n")
} else if (first_person_effect$CI_Upper < 0) {
  cat("Result: There is strong evidence (>97.5% probability) that the treatment DECREASES the risk for this individual.\n")
} else {
  cat("Result: The 95% credible interval includes zero, so there is no strong evidence of a treatment effect for this individual.\n")
}

