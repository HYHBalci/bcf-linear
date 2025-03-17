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

double logPosteriorTauJ(
    double tau_j,                // Proposed tau_j > 0
    double beta_j, //proposal for beta_{j}
    int index,  // index of {j}
    const std::vector<double> & beta_int, // Interaction terms: beta_{j,k} for all k in 'otherIdx'
    const std::vector<double> & tau, // The vector of all tau
    double tau_int,   
    double sigma,
    const std::vector<std::pair<int,int>> &int_pairs_trt = std::vector<std::pair<int,int>>(), //vector initializing the pairs made {j,k}
    bool interaction = false
) {
  // 1) check positivity
  if (tau_j <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  } 
  
  // 2) half-Cauchy(0,1) log prior: log(2/pi) - log(1 + tau_j^2)
  double logPrior = std::log(2.0 / M_PI) - std::log(1.0 + tau_j * tau_j);
  
  // 3) main effect log-likelihood: Normal(0, sigma^2 * tau_j^2)
  double log2pi  = std::log(2.0 * M_PI);
  double var_main= (sigma * sigma) * (tau_j * tau_j);
  double logLikMain = -0.5 * ( log2pi + std::log(var_main) ) 
    - 0.5 * ( (beta_j * beta_j) / var_main ); 
  
  // 4) interaction log-likelihood
  //    For each k in otherIdx, beta_{j,k} ~ Normal(0, sigma^2 * tau_j * tau_k * tau_int)
  double logLikInter = 0.0;
  if(interaction){
    for (size_t m=0; m < int_pairs_trt.size(); m++) {
      int iVar = int_pairs_trt[m].first;
      int jVar = int_pairs_trt[m].second;
      int target_idx = -1848; // index to fetch the other tau[]
      bool sample_flag = false;
      if (iVar == index){
        target_idx = jVar;
        sample_flag = true;
      } else if (jVar == index){ 
        target_idx =iVar;
        sample_flag = true;
      }
      if(sample_flag){
      double beta_jk = beta_int[m];
      double var_jk;
      if(target_idx == index){
        var_jk  = (sigma * sigma) * tau_j * tau_j * tau_int;
      } else {
        var_jk  = (sigma * sigma) * tau_j * tau[target_idx] * tau_int;
      }
      double b2      = beta_jk * beta_jk;
      
      double ll = -0.5 * (log2pi + std::log(var_jk)) 
        -0.5 * (b2 / var_jk); 
      logLikInter += ll;
      }}}  
  
  return logPrior + logLikMain + logLikInter;
}  

// [[Rcpp::export]]
double sample_tau_j_slice(
    double tau_old,
    double beta_j,            // main effect beta_j
    int index,                // index of {j}!
    const std::vector<double> & beta_int,  // all interaction betas
    const std::vector<double> & tau,       // all tau
    double tau_int,
    double sigma,
    bool interaction = true,  // or default false
    double step_out = 0.5,
    int max_steps = 50
)
{ 
  int p_int = beta_int.size();
  int p_mod = tau.size();
  std::vector<std::pair<int,int>> int_pairs_trt;
  int_pairs_trt.reserve(p_int);
  for(int ii = 0; ii < p_mod; ii++){
    for(int jj = ii; jj < p_mod; jj++){
      int_pairs_trt.push_back(std::make_pair(ii, jj));
    }
  }
  // 1) Evaluate log posterior at old tau
  double logP_old = logPosteriorTauJ(
    tau_old,
    beta_j,
    index,
    beta_int,
    tau,
    tau_int,
    sigma,
    int_pairs_trt,
    interaction
  ); 
  
  // 2) Draw a vertical level
  double u = R::runif(0.0, 1.0);
  double y_slice = logP_old + std::log(u);
   
  // 3) Create an interval [L, R] containing tau_old
  double L = std::max(1e-12, tau_old - step_out);
  double R = tau_old + step_out;
   
  // Step out left
  int step_count = 0;
  while ( (L > 1e-12)
            && (logPosteriorTauJ(L, beta_j, index, beta_int, tau, tau_int, sigma, int_pairs_trt, interaction) > y_slice)
            && (step_count < max_steps) )
  { 
    L = std::max(1e-12, L - step_out);
    step_count++;
  }
   
  // Step out right
  step_count = 0;
  while ( (logPosteriorTauJ(R, beta_j, index, beta_int, tau, tau_int, sigma, int_pairs_trt, interaction) > y_slice)
            && (step_count < max_steps) )
  { 
    R += step_out;
    step_count++;
  }
   
  // 4) Shrink bracket until we find a valid sample
  for (int rep = 0; rep < max_steps; rep++) {
    double prop = std::max(1e-8, R::runif(L, R)); // propose in [L,R]
    double lp   = logPosteriorTauJ(
      prop,
      beta_j,
      index,
      beta_int,
      tau,
      tau_int,
      sigma,
      int_pairs_trt,
      interaction
    ); 
    
    if (lp > y_slice) {
      // Accept
      return prop;
    } else { 
      // Shrink bracket
      if (prop < tau_old) {
        L = prop; 
      } else {
        R = prop;
      }
    } 
  }
  
  // If we never find a point in 'max_steps' tries, return a slightly perturbed tau_old
  return tau_old * (1.0 + 0.01 * R::runif(-1.0, 1.0));
} 

// [[Rcpp::export]]
double sample_alpha(
    int N,
    NumericVector r_alpha, // partial residual
    NumericVector z_,      // treatment indicator
    double sigma,
    double alpha_prior_sd = 10.0
) {
  // prior alpha ~ N(0, alpha_prior_sd^2)
  double prior_var   = alpha_prior_sd * alpha_prior_sd;
  double prec_prior  = 1.0 / prior_var;
  
  // gather sum of z_i and sum of z_i * r_alpha[i]
  double sum_z = 0.0;
  double sum_rz = 0.0;
  for(int i = 0; i < N; i++){
    sum_z  += z_[i];
    sum_rz += z_[i] * r_alpha[i];
  }
   
  // data precision = sum(z_i^2) / sigma^2 = (# treated) / sigma^2 if z is 0/1
  double prec_data = sum_z / (sigma * sigma);
  double prec_post = prec_data + prec_prior;
  double var_post  = 1.0 / prec_post;
  // posterior mean
  double mean_post = var_post * ( sum_rz / (sigma * sigma) );
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
  if (N == 0) return NA_REAL;  
  
  // Compute residual sum of squares using vectorized Rcpp
  double rss = sum(resid * resid);
  
  // Compute posterior parameters
  double shape_post = shape_prior + 0.5 * N;
  double rate_post  = rate_prior + 0.5 * rss;
  
  double gamma_draw = R::rgamma(shape_post, 1.0 / rate_post);
  
  // Return inverse of gamma draw (Inverse-Gamma sample)
  return 1.0 / gamma_draw;
}

// loglikeTauInt returns the PRIOR contribution (log p(beta_int | tau_int)).
// If tau_int is outside (0.01, 1), return -∞ (invalid).
// beta_int:   current vector of interaction coefficients
// tau:        vector of main-effect scale parameters tau[i]
// sigma:      current sigma (stdev of the error term)
// int_pairs:  vector of (iVar, jVar) pairs for interactions

double loglikeTauInt(
    double tau_int,
    // Baseline interaction info
    const std::vector<double> &beta_int_base,
    const std::vector<std::pair<int,int>> &int_pairs_base,
    // Main-effect shrinkage scales (tau) used in var = sigma^2 * tau_int * tau[iVar] * tau[jVar]
    const std::vector<double> &tau_main,
    double sigma,
    
    bool include_treatment_int = false,
    const std::vector<double> &beta_int_trt = std::vector<double>(),
    const std::vector<double> &tau_trt = std::vector<double>(),
    const std::vector<std::pair<int,int>> &int_pairs_trt = std::vector<std::pair<int,int>>()
)  
{
  // 1) Check Uniform(0.01,1.0) boundary (or whatever bounds you like)
  if(tau_int < 0.01 || tau_int > 1.0){
    return -std::numeric_limits<double>::infinity();
  } 
  
  double logp = 0.0;
  const double log2pi = std::log(2.0 * M_PI);
   
  // 2) Baseline interaction terms
  //    var_ij = tau_int * tau_main[iVar] * tau_main[jVar] * (sigma*sigma)
  //    log N(0, var_ij)
  for(size_t k = 0; k < int_pairs_base.size(); k++){
    int iVar = int_pairs_base[k].first;
    int jVar = int_pairs_base[k].second;
    double var_ij = tau_int * tau_main[iVar] * tau_main[jVar] * (sigma * sigma);
    double beta2  = beta_int_base[k]*beta_int_base[k];

    // log of N(0, var_ij) = -0.5 * [ log(2π var_ij) + beta^2 / var_ij ]
    logp += -0.5 * (log2pi + std::log(var_ij)) 
      - 0.5 * (beta2 / var_ij);
  } 
  
  // 3) Treatment interaction terms (if included)
  if(include_treatment_int && !beta_int_trt.empty()) {
    for(size_t k = 0; k < int_pairs_trt.size(); k++){
      int iVar = int_pairs_trt[k].first;
      int jVar = int_pairs_trt[k].second;

      double var_ij = tau_int * tau_main[iVar] * tau_main[jVar] * (sigma * sigma);
      double beta2  = beta_int_trt[k]*beta_int_trt[k];
    
      logp += -0.5 * (log2pi + std::log(var_ij))
        - 0.5 * (beta2 / var_ij);
    } 
  }
  return logp;
}
