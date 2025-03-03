#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

/////////////////////////////////////////////////////////////////////////////////////
// 1) Helper: compute_eta
//
//   param layout (size = 2p + pC2 + 2):
//     param[0]              = alpha
//     param[1..p]           = beta_j           (p main effects)
//     param[p+1..p + pC2]   = beta_{jk}, j<k   (pC2 interaction effects)
//     param[p + pC2 + 1..p + pC2 + p] = tau_j  (p half-Cauchy scales for main effects)
//     param[p + pC2 + p]    = tau_int
//     param[p + pC2 + p + 1]= sigma^2
//
// This function only returns the linear predictor, i.e. alpha + sum_j beta_j x_ij
// + sum_{j<k} beta_{jk} x_ij x_ik
/////////////////////////////////////////////////////////////////////////////////////
static arma::vec compute_eta(const arma::vec &param,
                             const arma::mat &X)
{
  int n = X.n_rows;
  int p = X.n_cols;
   
  double alpha = param[0];
  // main effects
  arma::vec beta = param.subvec(1, p);  // length p (1..p inclusive)
   
  // interactions 
  int pC2 = p*(p-1)/2;
  int start_inter = p + 1; // param index of first beta_{jk}
  int end_inter   = start_inter + pC2 - 1;
  arma::vec beta_int = param.subvec(start_inter, end_inter);
   
  // Build each eta_i
  arma::vec eta(n, fill::zeros);
  eta.fill(alpha);
   
  // main effects
  eta += X * beta;
   
  // interactions
  // We'll loop over j<k in row-major style
  int idx = 0;
  for(int j=0; j<p; j++){
    for(int k=j+1; k<p; k++){
      double b_jk = beta_int[idx];
      for(int i=0; i<n; i++){
        eta[i] += b_jk * X(i,j)*X(i,k);
      } 
      idx++;
    } 
  }
  return eta;
} 

inline double transform_tau_int(double w) {
    return 0.01 + (0.99 * exp(w)) / (1 + exp(w));
}
/////////////////////////////////////////////////////////////////////////////////////
// 2) log_posterior_linked_shrinkage
//    * includes likelihood + all priors as in the discussion:
//
// likelihood: y_i ~ Normal(eta_i, sigma^2)
//   => sum_i [-1/2 log(2pi) - 1/2 log(sigma^2) - (y_i - eta_i)^2/(2 sigma^2)]
//
// prior on alpha ~ N(0, 10^2)
// prior on beta_j ~ N(0, sigma^2 tau_j^2)
// prior on beta_{jk} ~ N(0, sigma^2 tau_j tau_k tau_int)
// prior on tau_j ~ half-Cauchy(0,1) => log p(tau_j) = log(2/pi) - log(1+tau_j^2)
// prior on tau_int ~ Uniform(0.01,1) => log p(tau_int) = - log(0.99)  if in (0.01,1)
// prior on sigma^2 ~ IG( shape=1, rate=0.001 )
//
// param layout repeated:
//   param[0] = alpha
//   param[1..p] = beta_j
//   param[p+1..p+pC2] = beta_{jk}
//   param[p+pC2+1.. p+pC2+p] = tau_j
//   param[p+pC2+p] = tau_int
//   param[p+pC2+p+1] = sigma^2
/////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double log_posterior_linked_shrinkage(const arma::vec &param,
                                      const arma::mat &X,
                                      const arma::vec &y)
{

  int n = X.n_rows;
  int p = X.n_cols;
  int pC2 = p * (p - 1) / 2;
  
  // Extract parameters
  double alpha = param[0];
  arma::vec beta_main = param.subvec(1, p);
  
  int start_inter = p + 1;
  int end_inter = start_inter + pC2 - 1;
  arma::vec beta_int = param.subvec(start_inter, end_inter);
  
  int start_tau = end_inter + 1;
  int end_tau = start_tau + (p - 1);
  arma::vec log_tau_vec = param.subvec(start_tau, end_tau);
  
  double w_tau_int = param[end_tau + 1];
  double log_sigma = param[end_tau + 2];
  
  arma::vec tau_vec = exp(log_tau_vec);  // **Ensures τ_j > 0**
  double tau_int = transform_tau_int(w_tau_int);  // **Ensures τ_int ∈ (0.01, 1)**
  double sigma2 = exp(log_sigma);  // **Ensures σ² > 0**
  
  if (sigma2 < 1e-10) sigma2 = 1e-10;  // Avoid division by zero
  tau_vec = clamp(tau_vec, 1e-10, 1e10);  // Avoid underflow/overflow
  
  // Compute log-likelihood
  arma::vec eta = compute_eta(param, X);
  arma::vec resid = y - eta;
  double ssr = arma::dot(resid, resid);
  
  double ll = -0.5 * n * std::log(2.0 * M_PI) - 0.5 * n *log_sigma - 0.5 * ssr / sigma2;
  
  // Priors
  double logprior_alpha = -0.5 * std::log(2.0 * M_PI * 100.0) - (alpha * alpha) / (2.0 * 100.0);
  double logprior_main = 0.0;
  for (int j = 0; j < p; j++) {
    double bj = beta_main[j];
    double tau_j = tau_vec[j];
    logprior_main += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma - 0.5 * std::log(tau_j * tau_j)
      - (bj * bj) / (2.0 * sigma2 * tau_j * tau_j);
  }
  
  double logprior_inter = 0.0;
  int idx = 0;
  for (int j = 0; j < p; j++) {
    double tau_j = tau_vec[j];
    for (int k = j + 1; k < p; k++) {
      double tau_k = tau_vec[k];
      double b_jk = beta_int[idx];
      double denom = sigma2 * tau_j * tau_k * tau_int;
      logprior_inter += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma - 0.5 * std::log(tau_j)
        - 0.5 * std::log(tau_k) - 0.5 * std::log(tau_int) - (b_jk * b_jk) / (2.0 * denom);
      idx++;
    }
  }
  
  double logprior_tau = 0.0;
  for (int j = 0; j < p; j++) {
    double tau_j = tau_vec[j];
    logprior_tau += std::log(2.0 / M_PI) - std::log(1.0 + tau_j * tau_j);
  }
  
  double logprior_tauint = -std::log(0.99);
  double logprior_sigma = std::log(0.001) - std::lgamma(1.0) - 2.0 * log_sigma - 0.001 / sigma2;
  
  return ll + logprior_alpha + logprior_main + logprior_inter + logprior_tau + logprior_tauint + logprior_sigma;
}

/////////////////////////////////////////////////////////////////////////////////////
// 3) grad_log_posterior_linked_shrinkage
//    gradient of the above log-posterior w.r.t. param
//    param layout is the same, we add partial derivatives from
//      * likelihood
//      * prior on alpha
//      * prior on main betas
//      * prior on interaction betas
//      * prior on tau_j (half-Cauchy)
//      * prior on tau_int (uniform => zero derivative inside domain)
//      * prior on sigma^2 (inverse gamma)
/////////////////////////////////////////////////////////////////////////////////////

// Helper: partial derivative of the linear predictor wrt alpha/beta_j/beta_{jk}
// for the *likelihood* portion. We'll do so by reusing compute_eta + residuals.
static arma::vec compute_residual(const arma::vec &param,
                                  const arma::mat &X,
                                  const arma::vec &y) {
  arma::vec eta = compute_eta(param, X);
  return (y - eta); // length n
}

// [[Rcpp::export]]
arma::vec grad_log_posterior_linked_shrinkage(const arma::vec &param,
                                               const arma::mat &X,
                                               const arma::vec &y)
 {
   int n = X.n_rows;
   int p = X.n_cols;
   int pC2 = p * (p - 1) / 2;
   
   // Extract parameters
   double alpha_val = param[0];
   arma::vec beta_main = param.subvec(1, p);
   
   int start_inter = p + 1;
   int end_inter = start_inter + pC2 - 1;
   arma::vec beta_int = param.subvec(start_inter, end_inter);
   
   int start_tau = end_inter + 1;
   int end_tau = start_tau + (p - 1);
   arma::vec log_tau_vec = param.subvec(start_tau, end_tau);
   
   double w_tau_int = param[end_tau + 1];
   double log_sigma = param[end_tau + 2];
   
   arma::vec tau_vec = exp(log_tau_vec);  // **Ensures τ_j > 0**
   double tau_int = 0.01 + (0.99 * exp(w_tau_int)) / (1 + exp(w_tau_int));  // **Ensures τ_int ∈ (0.01, 1)**
   double sigma2 = exp(log_sigma);  // **Ensures σ² > 0**
   
   // Prevent extreme values
   if (sigma2 < 1e-10) sigma2 = 1e-10;
   tau_vec = clamp(tau_vec, 1e-10, 1e10);
   
   // Compute residuals
   arma::vec resid = compute_residual(param, X, y);
   
   double ssr = arma::dot(resid, resid);
   
   // Initialize gradient vector
   arma::vec grad(param.size(), fill::zeros);
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  A) GRADIENT w.r.t. α (Intercept)
   grad[0] = arma::accu(resid) / sigma2;
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  B) GRADIENT w.r.t. β_j (Main effects)
   for (int j = 0; j < p; j++) {
     double bj = beta_main[j];
     double tj = tau_vec[j];
     
     double grad_beta = arma::dot(resid, X.col(j)) / sigma2;
     
     grad_beta -= bj / (sigma2 * tj * tj);
     
     grad[1 + j] = grad_beta;
   }
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  C) GRADIENT w.r.t. β_{jk} (Interaction effects)
   int idx = 0;
   for (int j = 0; j < p; j++) {
     double tj = tau_vec[j];
     for (int k = j + 1; k < p; k++) {
       double tk = tau_vec[k];
       double b_jk = beta_int[idx];
       
       double grad_beta_int = arma::dot(resid, X.col(j) % X.col(k)) / sigma2;
       
       double denom = sigma2 * tj * tk * tau_int;
       grad_beta_int -= b_jk / denom;
       
       grad[start_inter + idx] = grad_beta_int;
       idx++;
     }
   }
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  D) GRADIENT w.r.t. log(σ²)
   double d_log_sigma = -0.5 * n + ssr / (2.0 * sigma2);
   
   idx = 0;
   for (int j = 0; j < p; j++) {
     double bj = beta_main[j];
     double tj = tau_vec[j];
     
     // Contribution from prior on β_j ~ N(0, σ² τ_j²)
     d_log_sigma += -0.5 + (bj * bj) / (2.0 * sigma2 * tj * tj);
     for (int i = j +1; i < p; i++) {
       double b_jk = beta_int[idx];
       double tk = tau_vec[i];
       d_log_sigma += -0.5 + (b_jk * b_jk) / (2.0 * sigma2 * tj * tk * tau_int);
       idx ++;
     }
     
   }
   
   d_log_sigma += -2.0 + 0.001 / sigma2;
   
   grad[end_tau + 2] = d_log_sigma;
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // E) GRADIENT w.r.t. log(\u03c4_j)
   idx = 0;
   for (int j = 0; j < p; j++) {
     double tj = tau_vec[j];  // Extract \u03c4_j
     double bj = beta_main[j];  // Extract \u03b2_j
     
     // 1. Gradient from prior on \u03b2_j ~ N(0, \u03c3\u00b2 \u03c4_j\u00b2)
     double grad_tau_j = -1.0 / tj + (bj * bj) / (sigma2 * tj * tj * tj);
     
     // 2. Contribution from interactions involving \u03c4_j
     int idx = 0;  // Ensure correct indexing in beta_int
     for (int k = 0; k < p; k++) {
       if (j < k) {  // Only consider \u03b2_{j,k} when j < k
         double tk = tau_vec[k];
         double b_jk = beta_int[idx];
         grad_tau_j += -0.5 / tj + (b_jk * b_jk) / (2.0 * sigma2 * tj * tj * tk * tau_int);
         idx++;  // Move to the next interaction term
       } else if (j > k) {  // When j > k, find the corresponding interaction term
         int interaction_idx = (k * (2 * p - k - 1)) / 2 + (j - k - 1);
         double tk = tau_vec[k];
         double b_jk = beta_int[interaction_idx];
         grad_tau_j += -0.5 / tj + (b_jk * b_jk) / (2.0 * sigma2 * tj * tj * tk * tau_int);
       }
     }
     
     // 3. Add half-Cauchy prior on \u03c4_j
     grad_tau_j += -2.0 * tj / (1.0 + tj * tj);  // d/d\u03c4_j log p(\u03c4_j)
     
     // 4. Chain rule for log(\u03c4_j)
     grad[start_tau + j] = grad_tau_j * tj;
   }
   
   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  F) GRADIENT w.r.t. logit(τ_int)
   double d_tau_int = 0.0;
   int idx2 = 0;
   for (int j = 0; j < p; j++) {
     double tj = tau_vec[j];
     for (int k = j + 1; k < p; k++) {
       double tk = tau_vec[k];
       double b_jk = beta_int[idx2];
       d_tau_int += -0.5 / tau_int
       + (b_jk * b_jk) / (2.0 * sigma2 * tau_vec[j] * tau_vec[k] * (tau_int * tau_int));
       idx2++;
     }
   }
   
   double tau_int_transformed = 0.99 * std::exp(w_tau_int) 
     / std::pow(1.0 + std::exp(w_tau_int), 2.0);
     grad[end_tau + 1] = d_tau_int * tau_int_transformed;
   
   return grad;
 }

// [[Rcpp::export]]
Rcpp::List leapfrogCpp(
    arma::vec param,
    arma::vec momentum,
    double step_size,
    int num_steps,
    const arma::mat &X,
    const arma::vec &y
) {
  // First half-step for momentum
  momentum = momentum + 0.5 * step_size * grad_log_posterior_linked_shrinkage(param, X, y);
  
  // Full Leapfrog steps
  for (int i = 0; i < num_steps; i++) {
    // Full position step
    param = param + step_size * momentum;
    // Full momentum step (except for last)
    if (i != num_steps - 1) {
      momentum = momentum + step_size * grad_log_posterior_linked_shrinkage(param, X, y);
    }
  }
  
  // Final half-step for momentum
  momentum = momentum + 0.5 * step_size * grad_log_posterior_linked_shrinkage(param, X, y);
  
  return Rcpp::List::create(
    Named("param") = param,
    Named("momentum") = momentum
  );
}


