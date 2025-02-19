#define RCPP_NO_VISIBILITY 1
#undef attribute_visible
#define attribute_visible /* nothing */

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "horseshoe_samplers.h"  // If needed

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]

// ------------------------------------------------------------------------
// HELPER FUNCTIONS
// ------------------------------------------------------------------------

// Logistic function
inline double logistic(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// Log-likelihood for a single data point (y in {0,1})
inline double loglike_single(double y_i, double linpred_i) {
  double p = logistic(linpred_i);
  if(y_i > 0.5) {
    return std::log(p + 1e-12); 
  } else { 
    return std::log(1.0 - p + 1e-12);
  }
} 

// Sample sigma^2 from Inverse-Gamma
// shape_prior = a, rate_prior = b
// post shape = a + n/2,  post rate = b + sum(...) / 2
double sample_sigma2_ig(const arma::vec &resid, double shape_prior, double rate_prior) {
  // For logistic, 'resid' won't be the usual Gaussian residual, 
  // but let's keep the structure for demonstration.
  int n = resid.n_elem;
  double SSR = dot(resid, resid); // sum of squares (?)
  // Posterior
  double shape_post = shape_prior + 0.5*n;
  double rate_post  = rate_prior  + 0.5*SSR;
  return 1.0 / R::rgamma(shape_post, 1.0 / rate_post);
} 

// ------------------------------------------------------------------------
// Random-Walk Metropolis for parameters
// We include the prior in the acceptance ratio
// prior: param ~ Normal(0, priorVar)
// ------------------------------------------------------------------------
double update_param_metropolis(
    double old_val,
    const arma::vec &allfit_excl,  // linear predictor excluding old contribution
    const arma::vec &y,
    const arma::vec &x_contrib,    // how param influences the predictor = z[i]*X[i,j] or z[i]*X[i,iVar]*X[i,jVar]
    double proposal_sd,
    double priorVar                // = sigma^2 * tau_j^2 or sigma^2 * (tau_iVar*tau_jVar*tau_int)
) {
  double param_prop = R::rnorm(old_val, proposal_sd);
   
  // log prior ratio
  // old ~ N(0, priorVar), new ~ N(0, priorVar)
  double lp_old = -0.5 * (old_val * old_val) / priorVar;
  double lp_new = -0.5 * (param_prop * param_prop) / priorVar;
  double lp_prior_ratio = lp_new - lp_old;
   
  // log likelihood ratio
  double loglik_old = 0.0;
  double loglik_new = 0.0;
  int n = y.n_elem;
  for(int i=0; i<n; i++){
    // old linear predictor
    double lp_o = allfit_excl[i] + old_val * x_contrib[i];
    // new linear predictor
    double lp_n = allfit_excl[i] + param_prop * x_contrib[i];
     
    loglik_old += loglike_single(y[i], lp_o);
    loglik_new += loglike_single(y[i], lp_n);
  } 
  double lp_lik_ratio = loglik_new - loglik_old;
   
  double logAccept = lp_prior_ratio + lp_lik_ratio;
  if(std::log(R::runif(0.0,1.0)) < logAccept) {
    return param_prop; // accept
  } else { 
    return old_val;    // reject
  }
} 

// ------------------------------------------------------------------------
// MAIN FUNCTION
// ------------------------------------------------------------------------

// [[Rcpp::export]]
List run_mcmc_logistic_treat_inter(
    const arma::vec &y,          
    const arma::mat &X,         
    const arma::vec &z,          
     
    int n_iter = 2000,
    int burn   = 500,
    int thin   = 10,
     
    double sigma2_prior_a = 1.0,
    double sigma2_prior_b = 1e-3,
    double prop_sd_alpha  = 0.02,
    double prop_sd_bmain  = 0.02,
    double prop_sd_bint   = 0.02,
     
    double alpha_init     = 0.0,
    double sigma_init     = 1.0,
    bool interaction = true)
  {
  // Dimensions
  int n = X.n_rows;
  int p_mod = X.n_cols;      // number of main effects
  int p_int = (p_mod * (p_mod + 1)) / 2;
  int n_kept = (n_iter - burn)/thin;
  NumericMatrix beta_intOut(n_kept, p_int);
  std::vector<std::pair<int,int>> int_pairs;
  int_pairs.reserve(p_int);
  for(int ii = 0; ii < p_mod; ii++){
    for(int jj = ii; jj < p_mod; jj++){
      int_pairs.push_back(std::make_pair(ii, jj));
    }
  }
  // STORAGE for MCMC
  arma::vec   alphaOut(n_kept, fill::zeros);
  arma::mat   betaBaseOut(n_kept, p_mod, fill::zeros);
  arma::mat   betaTrtOut(n_kept, p_mod, fill::zeros);
  arma::mat   betaBaseIntOut(n_kept, p_int, fill::zeros);
  arma::mat   betaTrtIntOut(n_kept, p_int, fill::zeros);
  arma::vec   sigmaOut(n_kept, fill::zeros);
  arma::mat   tau_baseOut(n_kept, p_mod, fill::ones);
  arma::mat   tau_trtOut(n_kept, p_mod, fill::ones);
  arma::vec   tauIntOut(n_kept, fill::ones);
   
  // PARAMETERS
  double alpha = alpha_init;                // intercept for treatment portion
  arma::vec beta_base = arma::zeros(p_mod); // baseline main effects
  arma::vec beta_trt  = arma::zeros(p_mod); // treatment main effects
  arma::vec beta_base_int = arma::zeros(p_int); // baseline interactions
  arma::vec beta_trt_int  = arma::zeros(p_int); // treatment interactions
   
  // local shrinkage for each main effect dimension
  arma::vec tau_base = arma::ones(p_mod);
  arma::vec tau_trt = arma::ones(p_mod);
  // e.g. each dimension i has its own tau[i]
  // global scale for interactions
  double tau_int = 1.0; 
   
  double sigma_lin = sigma_init;     // "sigma" for prior scale
  double sigma2 = sigma_lin * sigma_lin;
   
  
  // ----------------------------------------------------------------------
  // CONSTRUCT LINEAR PREDICTOR
  //   LP[i] = baseline(i) + z[i]*[ alpha + trt(i) ]
  // baseline(i) = sum_j beta_base[j]*X(i,j) + sum_k beta_base_int[k]*[X_int1(i,k)*X_int2(i,k)]
  // trt(i)      = alpha + sum_j beta_trt[j]*X(i,j) + sum_k beta_trt_int[k]*(...)
  // We'll store total in allfit[i].
  // ----------------------------------------------------------------------
  arma::vec allfit(n, fill::zeros);
  for(int i=0; i<n; i++){
    // baseline portion
    double basePart = 0.0;
    for(int j=0; j<p_mod; j++){
      basePart += beta_base[j]*X(i,j);
    } 
    for(int k=0; k<p_int; k++){
      int iVar = int_pairs[k].first;
      int jVar = int_pairs[k].second;
      
      basePart += beta_base_int[k]*( X(i,iVar)*X(i,jVar) );
    } 
    
    // treatment portion
    double trtPart = alpha;
    for(int j=0; j<p_mod; j++){
      trtPart += beta_trt[j]*X(i,j);
    } 
    for(int k=0; k<p_int; k++){
      int iVar = int_pairs[k].first;
      int jVar = int_pairs[k].second;
      trtPart += beta_trt_int[k]*(X(i,iVar)*X(i,jVar));
    }
     
    allfit[i] = basePart + z[i]*trtPart;
  }
   
  // ----------------------------------------------------------------------
  // MCMC LOOP
  // ----------------------------------------------------------------------
  for(int iter=0; iter < n_iter; iter++){
    // ===============================
    // 1) Update baseline main effects (beta_base[j])
    // ===============================
    for(int j=0; j<p_mod; j++){
      // Remove old contribution from allfit
      for(int i=0; i<n; i++){
        allfit[i] -= beta_base[j]*X(i,j);
      } 
      // prior variance = sigma^2 * tau[j]^2
      double priorVar = sigma2 * (tau_base[j]*tau_base[j]);
      // Call random-walk Metropolis:
      // x_contrib[i] = X(i,j), no z here since it's baseline
      arma::vec x_contrib = X.col(j);
      // Build "allfit_excl" so that we can easily compute log-lik
      //   allfit_excl[i] = allfit[i] (which no longer includes old beta_base[j])
      // Then new param is added inside acceptance check
      arma::vec allfit_excl = allfit; 
      double old_val = beta_base[j];
      double new_val = update_param_metropolis(
        old_val, allfit_excl, y, x_contrib, prop_sd_bmain, priorVar
      ); 
      beta_base[j] = new_val;
      
      // Re-add new contribution
      for(int i=0; i<n; i++){
        allfit[i] += beta_base[j]*X(i,j);
      }
      
      // Optionally sample local shrinkage for baseline main effect j
      tau_base[j] = sample_tau_j_slice(tau_base[j], beta_base[j], sigma_lin);
    }
    
    // ===============================
    // 2) Update baseline interactions (beta_base_int[k])
    // ===============================
    for(int k=0; k<p_int; k++){
      int iVar = int_pairs[k].first;
      int jVar = int_pairs[k].second;
      // Remove old contribution
      for(int i=0; i<n; i++){
        double x_ij = X(i,iVar)*X(i,jVar);
        allfit[i] -= beta_base_int[k]*x_ij;
      }
      // priorVar = sigma^2 * tau[iVar]*tau[jVar]*tau_int (shared for baseline)
      // For demonstration, assume iVar = k, jVar = k is some pairing logic.
      // Typically you'd pass in the actual iVar, jVar for each k to get tau[iVar], tau[jVar].
      // We'll just do a placeholder: let's say we pick tau_int alone, or you handle iVar/jVar outside
      double priorVar = sigma2 * tau_base[iVar] * tau_base[jVar] *tau_int; 
      
      arma::vec x_contrib(n);
      for(int i=0; i<n; i++){
        x_contrib[i] = X(i,iVar)*X(i,jVar);
      }
      arma::vec allfit_excl = allfit;
      
      double old_val = beta_base_int[k];
      double new_val = update_param_metropolis(
        old_val, allfit_excl, y, x_contrib, prop_sd_bint, priorVar
      );
      beta_base_int[k] = new_val;
      
      // Re-add
      for(int i=0; i<n; i++){
        double x_ij = X(i,iVar)*X(i,jVar);
        allfit[i] += beta_base_int[k]*x_ij;
      }
    }
    
    // ===============================
    // 3) Update alpha, treatment main effects, treatment interactions
    //    We'll do these similarly, but multiplied by z[i].
    // ===============================
    
    // 3A) alpha
    {
      // Remove old alpha
      for(int i=0; i<n; i++){
        allfit[i] -= z[i]*alpha;
      }
      double priorVar = sigma2 * 1.0; // if alpha also ~ N(0, sigma^2)
      arma::vec x_contrib = z;       // so new param = alpha * z[i]
      arma::vec allfit_excl = allfit;
      double new_val = update_param_metropolis(
        alpha, allfit_excl, y, x_contrib, prop_sd_alpha, priorVar
      );
      alpha = new_val;
      // Re-add
      for(int i=0; i<n; i++){
        allfit[i] += z[i]*alpha;
      }
    }
    
    // 3B) treatment main effects (beta_trt[j])
    for(int j=0; j<p_mod; j++){
      // Remove old contribution
      for(int i=0; i<n; i++){
        allfit[i] -= z[i]*beta_trt[j]*X(i,j);
      }
      double priorVar = sigma2 * (tau_trt[j]*tau_trt[j]);
      // x_contrib[i] = z[i]*X(i,j)
      arma::vec x_contrib(n);
      for(int i=0; i<n; i++){
        x_contrib[i] = z[i]*X(i,j);
      }
      arma::vec allfit_excl = allfit;
      
      double old_val = beta_trt[j];
      double new_val = update_param_metropolis(
        old_val, allfit_excl, y, x_contrib, prop_sd_bmain, priorVar
      );
      beta_trt[j] = new_val;
      
      // Re-add
      for(int i=0; i<n; i++){
        allfit[i] += z[i]*beta_trt[j]*X(i,j);
      }
      // local shrinkage update if desired
      tau_trt[j] = sample_tau_j_slice(tau_trt[j], beta_base[j], sigma_lin);;
    }
    
    // 3C) treatment interactions (beta_trt_int[k])
    for(int k=0; k<p_int; k++){
      // Remove old
      int iVar = int_pairs[k].first;
      int jVar = int_pairs[k].second;
      
      for(int i=0; i<n; i++){
        double x_ij = X(i,iVar)*X(i,jVar);
        allfit[i] -= z[i]*beta_trt_int[k]*x_ij;
      }
      // priorVar = sigma^2 * tau[iVar]*tau[jVar]*tau_int
      double priorVar = sigma2 * tau_trt[iVar]*tau_trt[jVar]*tau_int;
      // x_contrib[i] = z[i]*X_int1(i,k)*X_int2(i,k)
      arma::vec x_contrib(n);
      for(int i=0; i<n; i++){
        x_contrib[i] = z[i]*(X(i,iVar)*X(i,jVar));
      }
      arma::vec allfit_excl = allfit;
      
      double old_val = beta_trt_int[k];
      double new_val = update_param_metropolis(
        old_val, allfit_excl, y, x_contrib, prop_sd_bint, priorVar
      );
      beta_trt_int[k] = new_val;
      
      // Re-add
      for(int i=0; i<n; i++){
        double x_ij = X(i,iVar)*X(i,jVar);
        allfit[i] += z[i]*beta_trt_int[k]*x_ij;
      }
    }
    
    // ===============================
    // 4) Update global interaction scale tau_int (placeholder)
    // ===============================
    // Example: Metropolis step for tau_int
    double currentTauInt = tau_int;
    double proposedTauInt = R::runif(0.01, 1.0); // sample from Uniform(0.01,1)
    
    double logPosteriorCurrent = loglikeTauInt(
      currentTauInt,
      arma::conv_to<std::vector<double>>::from(beta_base_int),  // Convert arma::vec to std::vector<double>
      int_pairs,
      arma::conv_to<std::vector<double>>::from(tau_base),
      sigma_lin,
      true,  // include_treatment_int
      arma::conv_to<std::vector<double>>::from(beta_trt_int),  // Convert arma::vec to std::vector<double>
      arma::conv_to<std::vector<double>>::from(tau_trt),
      int_pairs
    );
    
    double logPosteriorProposed = loglikeTauInt(
      proposedTauInt,
      arma::conv_to<std::vector<double>>::from(beta_base_int),
      int_pairs,
      arma::conv_to<std::vector<double>>::from(tau_base),
      sigma_lin,
      true,  // include_treatment_int
      arma::conv_to<std::vector<double>>::from(beta_trt_int),
      arma::conv_to<std::vector<double>>::from(tau_trt),
      int_pairs);
    
    double logAccept = logPosteriorProposed - logPosteriorCurrent;
    
    if(R::runif(0.0, 1.0) < std::exp(logAccept)) {
      tau_int = proposedTauInt;  
    } else {
      tau_int = currentTauInt;
    }
    
    // ===============================
    // 5) Update sigma^2 from IG prior (Optional in logistic)
    // ===============================
    {
      // We can define a "residual" for demonstration, 
      // though logistic doesn't have a direct Gaussian residual
      // We'll do sum of squares of [y - logistic(allfit)] for the code skeleton.
      arma::vec resid(n);
      for(int i=0; i<n; i++){
        double p_i = logistic(allfit[i]);
        resid[i]   = (y[i] - p_i); // "residual"
      }
      double new_sigma2 = sample_sigma2_ig(resid, sigma2_prior_a, sigma2_prior_b);
      sigma2 = new_sigma2;
      sigma_lin = std::sqrt(sigma2);
    }
    
    // ===============================
    // 6) Store MCMC if in keep range
    // ===============================
    if(iter >= burn && (iter - burn) % thin == 0){
      int idx = (iter - burn)/thin;
      alphaOut[idx]      = alpha;
      betaBaseOut.row(idx) = beta_base.t();
      betaTrtOut.row(idx)  = beta_trt.t();
      betaBaseIntOut.row(idx) = beta_base_int.t();
      betaTrtIntOut.row(idx)  = beta_trt_int.t();
      sigmaOut[idx]      = sigma_lin;
      tau_trtOut.row(idx)    = tau_trt.t();
      tau_baseOut.row(idx)    = tau_base.t();
      tauIntOut[idx]     = tau_int;
    }
  } // end MCMC
  
  return List::create(
    Named("alpha")         = alphaOut,
    Named("beta_base")     = betaBaseOut,
    Named("beta_trt")      = betaTrtOut,
    Named("beta_base_int") = betaBaseIntOut,
    Named("beta_trt_int")  = betaTrtIntOut,
    Named("sigma")         = sigmaOut,
    Named("tau_trt")           = tau_trtOut,
    Named("tau_base")           = tau_baseOut,
    Named("tau_int")       = tauIntOut
  );
}
