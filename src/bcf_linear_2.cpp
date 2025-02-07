#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include "horseshoe_samplers.h"

using namespace Rcpp;
// Rstudios check's suggest not ignoring these
// #pragma GCC diagnostic ignored "-Wunused-parameter"
// #pragma GCC diagnostic ignored "-Wcomment"
// #pragma GCC diagnostic ignored "-Wformat"
// #pragma GCC diagnostic ignored "-Wsign-compare"

// y = m(x) + b(x)z + e, e~N(0, sigma^2_y

//x_con is the design matrix for m. It should have n = rows
//x_mod is the design matrix for b. It should have n = rows
//data should come in sorted with all trt first, then control cases

//Tu ne cede malis, sed contra audentior ito!
//The boolean Tian indicates whether we follow the notation set up by (Tian et al., 2016), this is a test setting to see if the linear function holds. 
// [[Rcpp::export]]
List bcflineartwo(NumericVector y_, NumericVector z_, NumericVector w_,
                               NumericVector x_con_, NumericVector x_mod_,
                               List x_con_info_list, List x_mod_info_list,
                               arma::mat random_des, //needs to come in with n rows no matter what(?)
                               arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
                               double random_var_df,
                               int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
                               int ntree_mod, int ntree_con,
                               double lambda, double nu, //prior pars for sigma^2_y
                               double con_sd, // Var(m(x)) = con_sd^2 marginally a priori (approx)
                               double mod_sd, // Var(b(x)) = mod_sd^2 marginally a priori (approx)
                               double con_alpha, double con_beta,
                               double mod_alpha, double mod_beta,
                               CharacterVector treef_con_name_, CharacterVector treef_mod_name_,
                               int status_interval=100,
                               bool RJ= false, bool use_mscale=true, bool use_bscale=true, bool b_half_normal=true,
                               double trt_init = 0.0, bool verbose_sigma=false, 
                               bool no_output=false, bool tian = false)
{

  bool randeff = false;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  } 
  
  if(randeff) Rcout << "Using random effects." << std::endl;
  
  std::ofstream treef_con;
  // std::ofstream treef_mod;
   
  std::string treef_con_name = as<std::string>(treef_con_name_);
  // std::string treef_mod_name = as<std::string>(treef_mod_name_);
   
  if((not treef_con_name.empty()) && (not no_output)){
    Rcout << "Saving Trees to"  << std::endl;
    Rcout << treef_con_name  << std::endl;
    treef_con.open(treef_con_name.c_str());
  } else {    
    Rcout << "Not Saving Trees to file"  << std::endl;
  } 
  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws
   
  //double lambda = 1.0; //this one really needs to be set
  //double nu = 3.0;
  //double kfac=2.0; //original is 2.0
  
  Logger logger = Logger();
  char logBuff[100];
  bool log_level = false;
  logger.setLevel(log_level);
  logger.log("============================================================");
  logger.log(" Starting up BCF: ");
  logger.log("============================================================");
  if (log_level){
    logger.getVectorHead(y_, logBuff);
    Rcout << "y: " <<  logBuff << "\n";
    logger.getVectorHead(z_, logBuff);
    Rcout << "z: " <<  logBuff << "\n";
    logger.getVectorHead(w_, logBuff);
    Rcout << "w: " <<  logBuff << "\n";
  }
   
  
  logger.log("BCF is Weighted");
  // Rprintf(logBuff, "Updating Moderate Tree: %d of %d");
  // logger.log(logBuff);
  logger.log("");
   
  /*****************************************************************************
   * Read, format y
   *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
  double allys_y2 = 0;
  
  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys_y2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;
  
  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys_y2-n*ybar*ybar)/(n-1)); //sample standard deviation
  /*****************************************************************************
   * Read, format  weights
  *****************************************************************************/
   double* w = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
   
   
   for(int j=0; j<n; j++) {
     w[j] = w_[j];
   }
   
   
   /*****************************************************************************
    * Read, format X_con
    *****************************************************************************/
    //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
    std::vector<double> x_con;
   for(NumericVector::iterator it=x_con_.begin(); it!= x_con_.end(); ++it) {
     x_con.push_back(*it);
   }
   size_t p_con = x_con.size()/n;
   
   Rcout << "Using " << p_con << " control variables." << std::endl;
   
   //x cutpoints
   xinfo xi_con;
   
   xi_con.resize(p_con);
   for(int i=0; i<p_con; ++i) {
     NumericVector tmp = x_con_info_list[i];
     std::vector<double> tmp2;
     for(size_t j=0; j<tmp.size(); ++j) {
       tmp2.push_back(tmp[j]);
     }
     xi_con[i] = tmp2;
   } 
   
   /*****************************************************************************
    Read, format X_mod
    *****************************************************************************/
   int ntrt = 0;
   for(size_t i=0; i<n; ++i) {
     if(tian) {
       z_[i] = 1;
       ntrt = n;
     }
   }
     std::vector<double> x_mod;
     for(NumericVector::iterator it=x_mod_.begin(); it!= x_mod_.end(); ++it) {
       x_mod.push_back(*it);
     }
     size_t p_mod = x_mod.size()/n;
     
     Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;
     

     /*****************************************************************************
       Setup the model
       *****************************************************************************/
      //--------------------------------------------------
      // trees
      std::vector<tree> t_mod(ntree_mod);
      for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init/(double)ntree_mod);
      
      std::vector<tree> t_con(ntree_con);
      for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar/(double)ntree_con);
      
      //--------------------------------------------------
      //prior parameters
      // PX scale parameter for b:
      double bscale_prec = 2;
      double bscale0 = -0.5;
      double bscale1 = 0.5;
      
      double mscale_prec = 1.0;
      double mscale = 1.0;
      double delta_con = 1.0;
      double delta_mod = 1.0;
      
      double alpha_prior_sd  = 10.0; //Alpha prior variance for linear part 
      double sigma2_prior_a = 1.0;
      double sigma2_prior_b = 0.001;
      
      pinfo pi_mod;
      pi_mod.pbd = 1.0; //prob of birth/death move
      pi_mod.pb = .5; //prob of birth given  birth/death
      
      pi_mod.alpha = mod_alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
      pi_mod.beta  = mod_beta;  //2 for bart means it is harder to build big trees.
      pi_mod.tau   = mod_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod)); //sigma_mu, variance on leaf parameters
      pi_mod.sigma = shat; //resid variance is \sigma^2_y/bscale^2 in the backfitting update
      
      pinfo pi_con;
      pi_con.pbd = 1.0; //prob of birth/death move
      pi_con.pb = .5; //prob of birth given  birth/death
      
      pi_con.alpha = con_alpha;
      pi_con.beta  = con_beta;
      pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con)); //sigma_mu, variance on leaf parameters
      
      pi_con.sigma = shat/fabs(mscale); //resid variance in backfitting is \sigma^2_y/mscale^2
      
      double sigma = shat;
      
      // @Peter This is where dinfo is initialized
      
      //--------------------------------------------------
      //dinfo for control function m(x)
      //  Rcout << "ybar " << ybar << endl;
      double* allfit_con = new double[n]; //sum of fit of all trees
      for(size_t i=0;i<n;i++)
        if(tian){
          allfit_con[i] = 0;
        } else {
          allfit_con[i] = ybar;
        }
        double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
        dinfo di_con;
        di_con.n=n;
        di_con.p = p_con;
        di_con.x = &x_con[0];
        di_con.y = r_con; //the y for each draw will be the residual
        
        //--------------------------------------------------
        //dinfo for trt effect function b(x)
        double* allfit_mod = new double[n]; //sum of fit of all trees
        for(size_t i=0;i<n;i++) 
          if(tian){
            allfit_mod[i] = 0;
          } else {
            allfit_mod[i] = trt_init;
          }
          double* r_mod = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
          dinfo di_mod;
          di_mod.n=n;
          di_mod.p=p_mod;
          di_mod.x = &x_mod[0];
          di_mod.y = r_mod; //the y for each draw will be the residual
          
          //--------------------------------------------------
          //setup for random effects
          size_t random_dim = random_des.n_cols;
          int nr=1;
          if(randeff) nr = n;
          
          arma::vec r(nr); //working residuals
          arma::vec Wtr(random_dim); // W'r
          
          arma::mat WtW = random_des.t()*random_des; //W'W
          arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));
          
          // PX parameters
          arma::vec eta(random_var_ix.n_cols); //random_var_ix is num random effects by num variance components
          eta.fill(1.0);
          
          for(size_t k=0; k<nr; ++k) {
            r(k) = y[k] - allfit_con[k] - allfit_mod[k];
          }
          
          Wtr = random_des.t()*r;
          arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
          arma::vec allfit_random = random_des*gamma;
          if(!randeff) allfit_random.fill(0);
          
          //--------------------------------------------------
          //storage for the fits
          double* allfit = new double[n]; //yhat
          for(size_t i=0;i<n;i++) {
            allfit[i] = allfit_mod[i] + allfit_con[i];
            if(randeff) allfit[i] += allfit_random[i];
          }
          double* ftemp  = new double[n]; //fit of current tree
          
          NumericVector sigma_post(nd);
          NumericVector msd_post(nd);
          NumericVector bsd_post(nd);
          NumericVector b0_post(nd);
          NumericVector b1_post(nd);
          NumericMatrix m_post(nd,n);
          NumericMatrix yhat_post(nd,n);
          NumericMatrix b_post(nd,n);
          arma::mat gamma_post(nd,gamma.n_elem);
          arma::mat random_var_post(nd,random_var.n_elem);
          
          // Linear Part Posterior vectors and Matrix
          
          NumericVector alphaOut(nd), sigmaOut(nd);
          NumericMatrix betaOut(nd, p_mod);
          NumericMatrix tauOut(nd, p_mod);
          
          double alpha = 0.0;
          double sigma_lin = 1;
          std::vector<double> beta(p_mod, 0.0);
          std::vector<double> tau (p_mod, 1.0);
          
          NumericVector r_alpha(n), r_beta(n), resid(n); 
          NumericVector z(n, 1.0);
          
          //  NumericMatrix spred2(nd,dip.n);
          
          
          // The default output precision is of C++ is 5 or 6 dp, depending on compiler.
          // I don't have much justification for 32, but it seems like a sensible number
          int save_tree_precision = 32;
          
          //save stuff to tree file
          if(not treef_con_name.empty()){
            treef_con << std::setprecision(save_tree_precision) << xi_con << endl; //cutpoints
            treef_con << ntree_con << endl;  //number of trees
            treef_con << di_con.p << endl;  //dimension of x's
            treef_con << nd << endl;
            
            // treef_mod << std::setprecision(save_tree_precision) << xi_mod << endl; //cutpoints
            // treef_mod << ntree_mod << endl;  //number of trees
            // treef_mod << di_mod.p << endl;  //dimension of x's
            // treef_mod << nd << endl;
          }
        
  
  // MCMC Loop
  for(size_t iIter=0; iIter<(nd*thin+burn); iIter++) {
    if(verbose_sigma){
      if(iIter%status_interval==0) {
        Rcout << "iteration: " << iIter << " sigma/SD(y): "<< sigma_lin << endl;
      }
    }
    // -------------------------------------------
    // (A) Update alpha (intercept)
    // -------------------------------------------
    
    // Remove current alpha contribution
    for(int i=0; i<n; i++){
      allfit[i] -= alpha;
      allfit_mod[i] -= alpha;
    }
     
    // Compute residual r_alpha (excluding alpha)
    for(int i=0; i<n; i++){ 
      r_alpha[i] = y[i] - allfit[i];
    } 
    
    // Sample new alpha
    alpha = sample_alpha(n, r_alpha, z_, sigma_lin, alpha_prior_sd);
     
    // Add back newly sampled alpha
    for(int i=0; i<n; i++){ 
      allfit[i] += alpha;
      allfit_mod[i] += alpha;
    } 
    
    // -------------------------------------------
    // (B) Update each beta_j
    // -------------------------------------------
    for(int j=0; j<p_mod; j++) {
      // Remove old beta_j contribution from allfit
      double old_beta_j = beta[j];
      for(int i=0; i<n; i++){
        double old_contrib = old_beta_j * di_mod.x[i * p_mod + j];
        allfit[i] -= old_contrib;
        allfit_mod[i] -= old_contrib;
      }
      
      // Compute residual r_beta (excluding old beta_j)
      for(int i=0; i<n; i++){
        r_beta[i] = y[i] - allfit[i];
      } 
      
      // Prepare w_j for sampling
      NumericVector w_j(n);
      for(int i=0; i<n; i++){
        w_j[i] = di_mod.x[i * p_mod + j];
      }
      
      // Sample new beta_j
      if(tian) {
        beta[j] = sample_beta_j(n, r_beta, z_, w_j, 1, sigma_lin); // Normal prior (tau[j] = 1)
        tau[j] = sample_tau_j_slice(tau[j], beta[j], sigma_lin); 
      } else {
        beta[j] = sample_beta_j(n, r_beta, z_, w_j, tau[j], sigma_lin);
        tau[j] = sample_tau_j_slice(tau[j], beta[j], sigma_lin);
      }
      
      // Add new beta_j contribution back into allfit
      for(int i=0; i<n; i++){
        double contrib_diff = beta[j] * di_mod.x[i * p_mod + j];
        allfit[i] += contrib_diff; 
        allfit_mod[i] += contrib_diff;
      }
    }
     
    // -------------------------------------------
    // (C) Update sigma (residual variance)
    // -------------------------------------------
    for(int i=0; i<n; i++){
      resid[i] = y[i] - allfit[i];
    } 
    double sigma2 = sample_sigma2_ig(n, resid, sigma2_prior_a, sigma2_prior_b);
    sigma_lin = std::sqrt(sigma2);
    sigma = sigma_lin;
     
    // -------------------------------------------
    // (D) Save results for posterior samples
    // -------------------------------------------
    if((iIter >= burn) & (iIter % thin == 0)) {
      int it = (iIter - burn) / thin;
      alphaOut[it] = alpha;
      sigmaOut[it] = sigma;
      for(int j=0; j<p_mod; j++){
        betaOut(it, j) = beta[j];
        tauOut(it, j)  = tau[j];          
      }
    }
    
    if(iIter == 38) {  // Save only the first iteration for debugging
      std::ofstream debug_file("debug_iteration_38.txt");
      
      if (debug_file.is_open()) {
        debug_file << "================ DEBUG: Iteration 38 ================\n\n";
        
        // Save Scalars
        debug_file << "alpha = " << alpha << "\n";
        debug_file << "sigma = " << sigma << "\n";
        debug_file << "sigma_lin = " << sigma_lin << "\n";
        debug_file << "alpha_prior_sd = " << alpha_prior_sd << "\n";
        debug_file << "sigma2_prior_a = " << sigma2_prior_a << "\n";
        debug_file << "sigma2_prior_b = " << sigma2_prior_b << "\n";
        debug_file << "burn = " << burn << ", nd = " << nd << ", thin = " << thin << "\n";
        debug_file << "n = " << n << ", p_mod = " << p_mod << "\n";
        
        debug_file << "\n---------------- Arrays and Vectors ----------------\n";
        
        // Save Vectors
        debug_file << "\ny = [ ";
        for(int i=0; i<n; i++) debug_file << y[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nz_ = [ ";
        for(int i=0; i<n; i++) debug_file << z_[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nw_ = [ ";
        for(int i=0; i<n; i++) debug_file << w_[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nallfit = [ ";
        for(int i=0; i<n; i++) debug_file << allfit[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nallfit_mod = [ ";
        for(int i=0; i<n; i++) debug_file << allfit_mod[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nresid = [ ";
        for(int i=0; i<n; i++) debug_file << resid[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nr_alpha = [ ";
        for(int i=0; i<n; i++) debug_file << r_alpha[i] << " ";
        debug_file << "]\n";
        
        debug_file << "\nr_beta = [ ";
        for(int i=0; i<n; i++) debug_file << r_beta[i] << " ";
        debug_file << "]\n";
        
        // Save Beta and Tau
        debug_file << "\nbeta = [ ";
        for(int j=0; j<p_mod; j++) debug_file << beta[j] << " ";
        debug_file << "]\n";
        
        debug_file << "\ntau = [ ";
        for(int j=0; j<p_mod; j++) debug_file << tau[j] << " ";
        debug_file << "]\n";
        
        // Save X_mod Matrix
        debug_file << "\nX_mod (first 10 rows) = \n";
        for(int i=0; i<std::min<size_t>(10, n); i++) {  // Print only first 10 rows for readability
          for(int j=0; j<p_mod; j++) {
            debug_file << di_mod.x[i * di_mod.p + j] << " ";
          } 
          debug_file << "\n";
        }
        
        debug_file << "\n====================================================\n";
        debug_file.close();
    
  } 
    }
  }
  
   return(List::create(_["yhat_post"] = yhat_post, _["m_post"] = m_post, _["b_post"] = b_post,
                       _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post, _["b0"] = b0_post, _["b1"] = b1_post,
                       _["gamma"] = gamma_post, _["random_var_post"] = random_var_post, _["Beta"] = betaOut, _["tau"] = tauOut, _["alpha"] = alphaOut 
   ));
   
}