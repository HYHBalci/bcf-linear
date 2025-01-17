#include <Rcpp.h>
#include <cmath>          // for std::log, std::sqrt, etc.
#include <limits>         // for std::numeric_limits<double>::infinity()

// We need M_PI for pi. If not defined on your system, define it manually:
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Rcpp;

/*
 -------------------------------------------------------------------
 (1) sample_beta_j
 
 Given:
 - N:         number of data points
 - r_beta:    partial residual vector of length N
 r_beta[i] = y[i] - m_baseline[i] - Σ_{k≠j} [ z[i]*w_{i,k}*β_k ]
 - z:         treatment indicators, 0 or 1, length N
 - w_j:       the j-th column of the moderator design (length N)
 - tau_j:     local scale for β_j
 - sigma:     residual SD
 Returns:
 - A *new sample* for β_j from its conditional Normal posterior:
 β_j | (rest) ~ Normal(postMean, postVar)
 */

// [[Rcpp::export]]
double sample_beta_j(
    int N,
    NumericVector r_beta,
    NumericVector z,
    NumericVector w_j,
    double tau_j,
    double sigma
) {
  double sum_XZ2 = 0.0;  // sum of (z[i]*w_j[i])^2
  double sum_XZr = 0.0;  // sum of (z[i]*w_j[i]) * r_beta[i]
  
  for(int i = 0; i < N; i++){
    double xz = z[i] * w_j[i];
    sum_XZ2 += xz * xz;
    sum_XZr += xz * r_beta[i];
  }
  
  // prior variance for β_j is sigma^2 * tau_j^2
  double priorVar  = (sigma * sigma) * (tau_j * tau_j);
  // data precision from likelihood
  double prec_data  = sum_XZ2 / (sigma * sigma);
  // prior precision
  double prec_prior = 1.0 / priorVar;
  
  double postVar  = 1.0 / (prec_data + prec_prior);
  double postMean = postVar * (sum_XZr / (sigma * sigma));
  
  // Draw from Normal(postMean, sqrt(postVar)) using R::rnorm
  double draw = R::rnorm(postMean, std::sqrt(postVar));
  return draw;
  }

/*
 -------------------------------------------------------------------
 (2) sample_tau_j_slice 
 
 Uses a slice sampler to update tau_j > 0 from the posterior:
 p(tau_j | beta_j, sigma) ∝  [Half-Cauchy(0,1) on tau_j]
 × N(beta_j | 0, sigma^2 * tau_j^2)
 
 We define an internal function log_p_tau(...) to compute log posterior.
 */


// A helper for the log posterior of tau_j
static double log_p_tau(double tau, double beta_j, double sigma)
{ 
  // Half-Cauchy(0,1): p(tau) = 2/pi * 1/(1 + tau^2), tau>0
  // => log p(tau) = log(2/pi) - log(1 + tau^2)
  // Normal(0, sigma^2 * tau^2): p(beta_j|tau)
  // => log p(beta_j|tau) = -0.5 * [ log(2π) + log(sigma^2 τ^2) + (beta_j^2 / (sigma^2 τ^2)) ]
   
  if(tau <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  } 
  
  double log_prior = std::log(2.0 / M_PI) - std::log(1.0 + tau * tau);
   
  // log-likelihood of beta_j
  double log_lik = -0.5 * (
    std::log(2.0 * M_PI)
    + 2.0 * std::log(sigma)
    + 2.0 * std::log(tau)
    + (beta_j * beta_j) / ( (sigma * sigma) * (tau * tau) )
  ); 
  
  return log_prior + (-0.5) * log_lik * -2.0; 
  // ^ that line is suspect. Actually, let's do it carefully:
  // The expression inside is:
  //    -0.5*( log(2π) + 2log(sigma) + 2log(tau) + (beta_j^2/(σ^2 τ^2)) )
  // We'll rewrite it properly:
   
  log_lik = -0.5 * (
    std::log(2.0 * M_PI)
    + 2.0 * std::log(sigma)
    + 2.0 * std::log(tau)
    + ( (beta_j * beta_j) / ( (sigma * sigma) * (tau * tau) ) )
  );
   
  return log_prior + log_lik;
} 

// [[Rcpp::export]]
double sample_tau_j_slice(
    double tau_old,
    double beta_j,
    double sigma,
    double step_out = 0.5,
    int max_steps   = 50
) {
  // 1) Evaluate log posterior at old tau
  double logP_old = log_p_tau(tau_old, beta_j, sigma);
  // 2) Draw a vertical level
  double u = R::runif(0.0, 1.0);
  double y_slice = logP_old + std::log(u);

  // 3) Create an interval [L, R] containing tau_old
  double L = std::max(1e-12, tau_old - step_out);
  double R = tau_old + step_out;

  // Step out to the left
  while( (L > 1e-12) && (log_p_tau(L, beta_j, sigma) > y_slice) ) {
    L = std::max(1e-12, L - step_out);
  } 
  // Step out to the right
  while( log_p_tau(R, beta_j, sigma) > y_slice ) {
    R += step_out;
  }

  // 4) Shrink the bracket until we find a sample
  for(int rep = 0; rep < max_steps; rep++) {
    double prop = R::runif(L, R);
    double lp   = log_p_tau(prop, beta_j, sigma);

    if(lp > y_slice) {
      // Accept
      return prop;
    } else { 
      // Shrink bracket
      if(prop < tau_old) {
        L = prop;
       } else {
        R = prop;
      }
    }
  }
  
  // If we never find a point in 'max_steps' tries, return old value
  return tau_old;
}

// [[Rcpp::export]]
double sample_alpha(
    int N,
    NumericVector r_alpha,     // partial residual
    double sigma,              // current residual sd
    double alpha_prior_sd = 10.0
) {
  // prior alpha ~ N(0, alpha_prior_sd^2)
  double prior_var   = alpha_prior_sd * alpha_prior_sd;
  
  // data precision = N / sigma^2
  double prec_data   = (double)N / (sigma*sigma);
  // prior precision
  double prec_prior  = 1.0 / prior_var;
   
  double prec_post   = prec_data + prec_prior;
  double var_post    = 1.0 / prec_post;
   
  // sum of partial residual
  double sum_r = 0.0;
  for(int i=0; i<N; i++){
    sum_r += r_alpha[i];
  } 
  double mean_post = var_post * ( sum_r / (sigma*sigma) );
   
  // draw from Normal(mean_post, sqrt(var_post))
  double alpha_draw = R::rnorm(mean_post, std::sqrt(var_post));
  return alpha_draw;
} 

// [[Rcpp::export]]
double sample_sigma2_ig(
    int N,
    NumericVector resid,   // vector of residuals e_i
    double shape_prior = 1.0, 
    double rate_prior  = 0.001
) {
  // compute residual sum of squares
  double rss = 0.0;
  for(int i=0; i<N; i++){
    double e = resid[i];
    rss += e*e;
  }
   
  // posterior shape
  double shape_post = shape_prior + 0.5 * (double)N;
  // posterior rate
  double rate_post  = rate_prior + 0.5 * rss;

  double scale_post = 1.0 / rate_post;
  // gamma draw
  double gamma_draw = R::rgamma(shape_post, scale_post);

  
  double sigma2_draw = 1.0 / gamma_draw;
  return sigma2_draw;
}
