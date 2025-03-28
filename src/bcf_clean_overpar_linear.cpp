#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <utility>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include "horseshoe_samplers.h"
#include "linked_shrinkage_ll_lg.h"

using namespace Rcpp;
using namespace arma;
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
List bcfoverparRcppCleanLinear(NumericVector y_, NumericVector z_, NumericVector w_,
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
                         bool no_output=false, bool intTreat = true, bool hamiltonian = false,
                         double step_size = 0.005, int num_of_steps = 10, bool sparse = false)
{
  R_FlushConsole(); // Flush the console so to not overload it with all the messages printed. 
  bool randeff = false;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  } 
  
  if(randeff) Rcout << "Using random effects." << std::endl;
  
  Rcout << 'Alo?' << std::endl;
  
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
    /* Read, format  weights
     *****************************************************************************/
    double* w = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
    
    
    for(int j=0; j<n; j++) {
      w[j] = w_[j];
    }
    
    
    /*****************************************************************************
     /* Read, format X_con
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
      if(z_[i]>0) ntrt += 1;
    }
    std::vector<double> x_mod;
    for(NumericVector::iterator it=x_mod_.begin(); it!= x_mod_.end(); ++it) {
      x_mod.push_back(*it);
    }
    size_t p_mod = x_mod.size()/n;
    Rcout << "x_mod_ size: " << x_mod_.size() << ", expected: " << (n * p_mod) << std::endl;
    for (size_t i = 0; i < std::min<size_t>(5, n); i++) {
      Rcout << "x_mod_[" << i << "] = " << x_mod_[i] << std::endl;
    }
    
    // arma::mat X(x_mod_.begin(), n, p_mod, false, true);
    // 
    // Rcout << "X dimensions: " << X.n_rows << " x " << X.n_cols << std::endl;
    // for (size_t i = 0; i < n; i++) {
    //   for (size_t j = 0; j < p_mod; j++) {
    //     Rcout << "X(" << i << "," << j << ") = " << X(i, j) << "  ";
    //   } 
    //   Rcout << std::endl;
    // }
    
    Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;
    
    //x cutpoints
    // xinfo xi_mod;
    // 
    // xi_mod.resize(p_mod);
    // for(int i=0; i<p_mod; ++i) {
    //   NumericVector tmp = x_mod_info_list[i];
    //   std::vector<double> tmp2;
    //   for(size_t j=0; j<tmp.size(); ++j) {
    //     tmp2.push_back(tmp[j]);
    //   }
    //   xi_mod[i] = tmp2;
    // }
    
    //  Rcout <<"\nburn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
    //  Rcout <<"\nlambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;
    
    /*****************************************************************************
    /* Setup the model
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
    double delta_con;
    double delta_mod;
    if(sparse == false){
      delta_con = 1.0;
      delta_mod = 1.0;
    } else {
      delta_con = 2.0;
      delta_mod = 2.0;
    }
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
    for(size_t i=0;i<n;i++) allfit_con[i] = ybar;
    double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
    dinfo di_con;
    di_con.n=n;
    di_con.p = p_con;
    di_con.x = &x_con[0];
    di_con.y = r_con; //the y for each draw will be the residual
    
    //--------------------------------------------------
    //dinfo for trt effect function b(x)
    double* allfit_mod = new double[n]; //sum of fit of all trees
    for(size_t i=0;i<n;i++) allfit_mod[i] = trt_init;
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
    
    // Initialize for the interaction terms
    NumericVector tau_int_post(nd);
    
    int p_int;
    if(!hamiltonian){
      p_int = (p_mod * (p_mod + 1)) / 2;
    } else {
      p_int = (p_mod * (p_mod - 1)) / 2;
    }
    NumericMatrix beta_intOut(nd, p_int);
    std::vector<std::pair<int,int>> int_pairs;
    int_pairs.reserve(p_int);
    for(int ii = 0; ii < p_mod; ii++){
      for(int jj = ii; jj < p_mod; jj++){
        int_pairs.push_back(std::make_pair(ii, jj));
        }
    }
    // acceptance counter if there is hamiltonian MCMC.
    
    int acceptance = 0;     // Initial guess for the stepsize
    double target_accept = 0.65; // Target acceptance rate
    double H_bar = 0.0;
    double avg_log_step = -9;   // Running average of log(step_size)
    double mu = std::log(10.0 * step_size); // Bias term
    double ksi = 0.05;         // Controls stability of adaptation
    double t0 = 10.0;            // Initial slow adaptation
    double kappa = 0.75;         // Controls decay of adaptation
    
    
    //initialize parameters for HMCMC
    int total_size = 1 + p_mod + p_int + p_mod + 1 + 1;  // Total number of parameters
    arma::vec init_param(total_size, arma::fill::zeros);
    init_param[0] = 1;
    init_param.subvec(1, p_mod) = arma::vec(p_mod, arma::fill::zeros);
    init_param.subvec(p_mod + 1, p_mod + p_int) = arma::vec(p_int, arma::fill::zeros);
    init_param.subvec(p_mod + p_int + 1, p_mod + p_int + p_mod) = arma::vec(p_mod, arma::fill::ones);
    init_param[p_mod + p_int + p_mod + 1] = 0.5;
    init_param[p_mod + p_int + p_mod + 2] = 1;
    arma::vec param = init_param;
    arma::mat samples(nd, init_param.size(), arma::fill::none);
    
    double alpha = 0.0;
    double sigma_lin = 1;
    std::vector<double> beta(p_mod, 0.0);
    std::vector<double> tau (p_mod, 1.0);
    
    NumericVector r_alpha(n), r_beta(n), resid(n); 
    NumericVector z(n, 1.0);
    
    //for interaction terms
    std::vector<double> beta_int(p_int, 0.0);
    double tau_int = 0.5;
    
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
    
    //*****************************************************************************
    /* MCMC
     * note: the allfit objects are all carrying the appropriate scales
     */
    //*****************************************************************************
    Rcout << "\n============================================================\nBeginning MCMC:\n============================================================\n";
    time_t tp;
    int time1 = time(&tp);
    
    size_t save_ctr = 0;
    bool verbose_itr = true;
    
    
    double* weight      = new double[n];
    double* weight_het  = new double[n];
    
    logger.setLevel(0);
    
    bool printTrees = false;
    
    for(size_t iIter=0;iIter<(nd*thin+burn);iIter++) {
      // verbose_itr = iIter>=burn;
      verbose_itr = false;
      
      if(verbose_sigma){
        if(iIter%status_interval==0) {
          Rcout << "iteration: " << iIter << " sigma/SD(y): "<< sigma << endl;
        }
      }
      
      logger.setLevel(verbose_itr);
      
      logger.log("==============================================");
      Rprintf(logBuff, "MCMC iteration: %d of %d Start", iIter + 1, nd*thin+burn);
      logger.log(logBuff);
      Rprintf(logBuff, "sigma %f, mscale %f, bscale0 %f, bscale1 %f",sigma, mscale, bscale0, bscale1);
      logger.log(logBuff);
      logger.log("==============================================");
      if (verbose_itr){
        logger.getVectorHead(y, logBuff);
        Rcout << "           y: " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit, logBuff);
        Rcout << "Current Fit : " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit_con, logBuff);
        Rcout << "allfit_con  : " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit_mod, logBuff);
        Rcout << "allfit_mod  : " <<  logBuff << "\n";
      }
      
      for (int k=0; k<n; ++k){
        weight[k] = w[k]*mscale*mscale/(sigma * sigma); // for non-het case, weights need to be divided by sigma square to make it similar to phi
      }
      
      for(size_t k=0; k<ntrt; ++k) {
        weight_het[k] = w[k]*bscale1*bscale1/(sigma*sigma);
      }
      for(size_t k=ntrt; k<n; ++k) {
        weight_het[k] = w[k]*bscale0*bscale0/(sigma*sigma);
      }
      
      logger.log("=====================================");
      logger.log("- Tree Processing");
      logger.log("=====================================");
      
      //draw trees for m(x)
      for(size_t iTreeCon=0;iTreeCon<ntree_con;iTreeCon++) {
        
        logger.log("==================================");
        Rprintf(logBuff, "Updating Control Tree: %d of %d",iTreeCon + 1 , ntree_con);
        logger.log(logBuff);
        logger.log("==================================");
        logger.startContext();
        
        logger.log("Attempting to Print Tree Pre Update \n");
        if(verbose_itr && printTrees){
          t_con[iTreeCon].pr(xi_con);
          Rcout << "\n\n";
        }
        
        fit(t_con[iTreeCon], // tree& t
            xi_con, // xinfo& xi
            di_con, // dinfo& di
            ftemp); // std::vector<double>& fv
        
        
        logger.log("Attempting to Print Tree Post first call to fit \n");
        if(verbose_itr && printTrees){
          t_con[iTreeCon].pr(xi_con);
          Rcout << "\n\n";
        }
        
        for(size_t k=0;k<n;k++) {
          if(ftemp[k] != ftemp[k]) {
            Rcout << "control tree " << iTreeCon <<" obs "<< k<<" "<< endl;
            Rcout << t_con[iTreeCon] << endl;
            stop("nan in ftemp");
          }
          
          allfit[k]     = allfit[k]     -mscale*ftemp[k];
          allfit_con[k] = allfit_con[k] -mscale*ftemp[k];
          
          r_con[k] = (y[k]-allfit[k])/mscale;
          
          if(r_con[k] != r_con[k]) {
            Rcout << (y[k]-allfit[k]) << endl;
            Rcout << mscale << endl;
            Rcout << r_con[k] << endl;
            stop("NaN in resid");
          }
        }
        
        
        
        if(verbose_itr && printTrees){
          logger.getVectorHead(weight, logBuff);
          Rcout << "\n weight: " <<  logBuff << "\n\n";
        }
        logger.log("Starting Birth / Death Processing");
        logger.startContext();
        bd(t_con[iTreeCon], // tree& x
           xi_con, // xinfo& xi
           di_con, // dinfo& di
           weight, // phi
           pi_con, // pinfo& pi
           gen,
           logger); // RNG& gen
        logger.stopContext();
        
        logger.log("Attempting to Print Tree Post db \n");
        if(verbose_itr && printTrees){
          t_con[iTreeCon].pr(xi_con);
          Rcout << "\n";
        }
        
        if (verbose_itr && printTrees){
          logger.log("Printing Current Status of Fit");
          
          logger.getVectorHead(z_, logBuff);
          // logger.log(logBuff);
          Rcout << "\n          z : " <<  logBuff << "\n";
          
          logger.getVectorHead(y, logBuff);
          Rcout << "          y : " <<  logBuff << "\n";
          
          logger.getVectorHead(allfit, logBuff);
          Rcout << "Fit - Tree  : " <<  logBuff << "\n";
          
          logger.getVectorHead(r_con, logBuff);
          Rcout << "     r_con  : " <<  logBuff << "\n\n";
          
          Rcout <<" MScale: " << mscale << "\n";
          
          Rcout <<" bscale0 : " << bscale0 << "\n";
          
          Rcout <<" bscale1 : " << bscale1 << "\n\n";
          
        }
        logger.log("Starting To Draw Mu");
        logger.startContext();
        
        drmu(t_con[iTreeCon],  // tree& x
             xi_con, // xinfo& xi
             di_con, // dinfo& di
             pi_con, // pinfo& pi,
             weight,
             gen); // RNG& gen
        
        logger.stopContext();
        
        logger.log("Attempting to Print Tree Post drmu \n");
        if(verbose_itr  && printTrees){
          t_con[iTreeCon].pr(xi_con);
          Rcout << "\n";
        }
        
        fit(t_con[iTreeCon],
            xi_con,
            di_con,
            ftemp);
        
        for(size_t k=0;k<n;k++) {
          allfit[k] += mscale*ftemp[k];
          allfit_con[k] += mscale*ftemp[k];
        }
        
        logger.log("Attempting to Print tree Post second call to fit \n");
        
        if(verbose_itr && printTrees){
          t_con[iTreeCon].pr(xi_con);
          Rcout << "\n";
          
        }
        logger.stopContext();
      }
      
      // Treatment Effect; Model this in a linear way drawing inspiration from 'Think Before you shrink'
      
      if(!hamiltonian){
        for(int i=0; i<n; i++){
          // double sumBX = 0.0;
          // for(int j=0; j<p_mod; j++){
          //   sumBX += beta[j]*di_mod.x[i*di_mod.p + j];
          // }
          // double bscale = (i < ntrt) ? bscale1 : bscale0;
          allfit[i]     = allfit[i]     -z_[i]*alpha;
          allfit_mod[i] = allfit_mod[i] -z_[i]*alpha;
          
          r_alpha[i] = y[i] - allfit[i];
        }
        
        alpha = sample_alpha(n, r_alpha, z_, sigma, alpha_prior_sd);
        
        for(int i = 0; i < n; i++){
          allfit[i]     = allfit[i]     +z_[i]*alpha;
          allfit_mod[i] = allfit_mod[i] +z_[i]*alpha;
        }
        
        for(int j=0; j<p_mod; j++){
          // partial residual wrt beta_j
          for(int i=0; i<n; i++){
            allfit[i]     = allfit[i]     -z_[i]*beta[j]*di_mod.x[i*di_mod.p + j]; //subtract the contributions. 
            allfit_mod[i] = allfit_mod[i] -z_[i]*beta[j]*di_mod.x[i*di_mod.p + j];
            r_beta[i] = y[i] - allfit[i];
          }
          NumericVector w_j(n);
          for(int i=0; i<n; i++){
            w_j[i] = di_mod.x[i*di_mod.p + j];
          }
          beta[j] = sample_beta_j(n, r_beta, z_, w_j, tau[j], sigma);
          if(intTreat){
            tau[j] = sample_tau_j_slice(tau[j], beta[j], j, beta_int, tau, tau_int, sigma, intTreat = true);
          } else {
            tau[j] = sample_tau_j_slice(tau[j], beta[j], j, beta_int, tau, tau_int, sigma, intTreat = true);
          }
          for(int i=0; i<n; i++){
            double bscale = (i < ntrt) ? bscale1 : bscale0;
            allfit[i]     = allfit[i]     +z_[i]*beta[j]*di_mod.x[i*di_mod.p + j]; //re-add the contributions. 
            allfit_mod[i] = allfit_mod[i] +z_[i]*beta[j]*di_mod.x[i*di_mod.p + j];
          }
        }
        if(intTreat){
          for(int k = 0; k < p_int; k++) {
            // (iVar, jVar) for the pair
            int iVar = int_pairs[k].first;
            int jVar = int_pairs[k].second;
            
            // 2) Remove old contribution of beta_int[k] from allfit[]
            for(int i = 0; i < n; i++){
              double x_ij = di_mod.x[i*di_mod.p + iVar] * di_mod.x[i*di_mod.p + jVar];
              allfit[i]   -= z_[i] * beta_int[k] * x_ij;
              // partial residual for sampling
              r_beta[i]    = y[i] - allfit[i];
            }
            
            // 3) Sample the new beta_int[k] using a function that accounts for:
            //    - Prior: N(0, tau_int * tau[iVar] * tau[jVar] * sigma^2)
            //    - Data likelihood
            //    (Adapt the same approach you used for sample_beta_j(...))
            NumericVector w_j(n);
            for(int i=0; i<n; i++){
              w_j[i] = di_mod.x[i*di_mod.p + iVar] * di_mod.x[i*di_mod.p + jVar];
            }
            beta_int[k] = sample_beta_j(n, r_beta, z_, w_j, std::sqrt(tau_int*tau[iVar]*tau[jVar]), sigma);
            
            // 4) Add the *new* contribution back
            for(int i = 0; i < n; i++){
              double x_ij = di_mod.x[i*di_mod.p + iVar] * di_mod.x[i*di_mod.p + jVar];
              allfit[i]   += z_[i] * beta_int[k] * x_ij;
            }
          }
          double currentTauInt = tau_int;
          double proposedTauInt = R::runif(0.01, 1.0); // sample from Uniform(0.01,1)
          
          // Evaluate log posterior ratio = log p(data|tau_int_prop) + log p(tau_int_prop)
          //                              - (log p(data|tau_int_current) + log p(tau_int_current))
          // log p(tau_int) = 0 if uniform(0.01,1), ignoring the boundary checks.
          
          double logPosteriorCurrent = loglikeTauInt(
            currentTauInt,
            beta_int,
            int_pairs,
            tau,
            sigma,
            false,  
            std::vector<double>(),   
            std::vector<double>(),   
            std::vector<std::pair<int,int>>() 
          );
          
          double logPosteriorProposed = loglikeTauInt(
            proposedTauInt,
            beta_int,
            int_pairs,
            tau,
            sigma,
            false,
            std::vector<double>(),
            std::vector<double>(),
            std::vector<std::pair<int,int>>()
          );
          
          double logAccept = logPosteriorProposed - logPosteriorCurrent;
          
          if(R::runif(0.0, 1.0) < std::exp(logAccept)) {
            tau_int = proposedTauInt;  
          } else {
            tau_int = currentTauInt;
          }
          
        }
        
        for(int i=0; i<n; i++){
          resid[i] = y[i] - allfit[i];}
        
        double sigma2 = sample_sigma2_ig(n, resid, sigma2_prior_a, sigma2_prior_b);
        sigma_lin = std::sqrt(sigma2); 
        
        // Save results of our iteration  
        if(((iIter>=burn) & (iIter % thin==0)) )  {
          int it = (iIter - burn) / thin;
          alphaOut[it] = alpha;
          sigmaOut[it] = sigma;
          for(int j=0; j<p_mod; j++){
            betaOut(it, j) = beta[j];
            tauOut(it, j)  = tau[j];  
          } 
          if(intTreat){
            tau_int_post[it] = tau_int;
            for(int j=0; j<p_int; j++){
              beta_intOut(it, j) = beta_int[j];
            } 
          }
        } 
      } else {
        for(int i=0; i<n; i++){
          allfit[i]     = allfit[i] - allfit_mod[i];
          resid[i] = y[i] - allfit[i];
        }
        // if(iIter < burn){
        //   // Sample random momentum
        //   int i = iIter + 1;
        //   double log_step = std::log(step_size);
        //   arma::vec momentum = arma::randn<arma::vec>(init_param.n_elem);
        // 
        //   // Compute initial Hamiltonian
        //   Rcout << "Calling LL function..."<< std::endl;
        // 
        //   double H_old = -log_posterior_linked_shrinkage(param, arma::mat(x_mod_.begin(), n, p_mod, false, true), resid) + 0.5 * arma::dot(momentum, momentum);
        //   Rcout << "H_OLD "<< H_old << std::endl;
        //   // Perform Leapfrog step
        //   Rcout << "Calling Leapfrog function..."<< std::endl;
        //   Rcout << exp(avg_log_step)<< std::endl;
        //   Rcpp::List leap = leapfrogCpp(param, momentum, std::exp(avg_log_step), 10, arma::mat(x_mod_.begin(), n, p_mod, false, true), resid);
        //   arma::vec param_new = leap["param"];
        //   // Rcout << "Param New:  "<< param_new << std::endl;
        //   arma::vec momentum_new = leap["momentum"];
        // 
        //   // Compute new Hamiltonian
        //   Rcout << "Log posterior called again"<< std::endl;
        //   double H_new = -log_posterior_linked_shrinkage(param_new, arma::mat(x_mod_.begin(), n, p_mod, false, true), resid) + 0.5 * arma::dot(momentum_new, momentum_new);
        //   Rcout << "H_new "<< H_new<< std::endl;
        //   // Compute acceptance probability
        //   double alpha_cal = std::min(1.0, std::exp(H_old - H_new));
        //   Rcout << "alpha_cal "<< alpha_cal << std::endl;
        //   // Adapt step size
        //   H_bar = (1 - 1.0 / (i + t0)) * H_bar + (1.0 / (i + t0)) * (target_accept - alpha_cal);
        //   Rcout << "H_bar"<< alpha_cal << std::endl;
        //   log_step = mu - (std::sqrt(i) / ksi) * H_bar;
        //   log_step = std::max(std::log(1e-4), std::min(log_step, std::log(0.01)));
        //   Rcout << "log_step"<< log_step << std::endl;
        //   double eta = std::pow(i, -kappa);
        //   Rcout << "eta"<< eta << std::endl;
        //   avg_log_step = eta * log_step + (1 - eta) * avg_log_step;
        //   avg_log_step = std::max(std::log(1e-4), std::min(avg_log_step, std::log(1.0)));
        //   for(int i=0; i<n; i++){
        //     allfit[i]     = allfit[i] + allfit_mod[i];
        //     
        //   }
        // } else { 
        arma::vec momentum = arma::randn<arma::vec>(init_param.size());
        
        // 2) Compute initial Hamiltonian

        double H_old = -log_posterior_linked_shrinkage(param, arma::mat(x_mod_.begin(), n, p_mod, false, true), resid)+ 0.5 * arma::dot(momentum, momentum);
        // 3) Perform Leapfrog integration
        Rcpp::List leap = leapfrogCpp(param, momentum, step_size, num_of_steps ,arma::mat(x_mod_.begin(), n, p_mod, false, true), resid);
      
        arma::vec param_new = leap["param"];
        arma::vec momentum_new = leap["momentum"];

        // 4) Compute new Hamiltonian
        double H_new = -log_posterior_linked_shrinkage(param_new,arma::mat(x_mod_.begin(), n, p_mod, false, true), resid)
          + 0.5 * arma::dot(momentum_new, momentum_new);
        Rcout << "H_new"<< H_new << std::endl;
        // 5) Metropolis accept/reject
        double log_accept_ratio = H_old - H_new;
        double u = R::runif(0.0, 1.0);
        if (std::log(u) < log_accept_ratio) {
          // Accept move
          param = param_new;
          acceptance++;
        } 
        if(((iIter>=burn) & (iIter % thin==0)) )  {
          int it = (iIter - burn) / thin;
          samples.row(it) = param.t();
        }
        
        for(int i=0; i<n; i++){
          // 1) Start with intercept
          allfit_mod[i] = param(0);
          
          // 2) Add main effects
          for(int j = 0; j < p_mod; j++){
            double x_ij = x_mod[i * p_mod + j];  
            allfit_mod[i] += param(1 + j) * x_ij;
          }
          int idx = 0;
          for(int k = 0; k < (int)int_pairs.size(); k++){
            int colA = int_pairs[k].first;
            int colB = int_pairs[k].second;
            
            double x_ij = x_mod[i * p_mod + colA] * x_mod[i * p_mod + colB];
            allfit_mod[i] += param(p_mod + 1 + idx) * x_ij;
            idx++;
          }
          
          allfit[i]     = allfit[i] + allfit_mod[i];
        }
      }
      
      
      logger.setLevel(verbose_itr);
      
      logger.log("=====================================");
      logger.log("- MCMC iteration Cleanup");
      logger.log("=====================================");
      // double ww0 = 0.0, ww1 = 0.;
      // double rw0 = 0.0, rw1 = 0.;
      // double bscale_fc_var = 1/(ww1 + bscale_prec);
      // double bscale1_old = bscale1; 
      // double bscale0_old = bscale0;
      
      if(use_bscale) {
        double ww0 = 0.0, ww1 = 0.;
        double rw0 = 0.0, rw1 = 0.;
        double s2 = sigma*sigma;
        for(size_t k=0; k<n; ++k) {
          double bscale = (k<ntrt) ? bscale1 : bscale0;
          double scale_factor = (w[k]*allfit_mod[k]*allfit_mod[k])/(s2*bscale*bscale);
          
          if(scale_factor!=scale_factor) {
            Rcout << " scale_factor " << scale_factor << endl;
            stop("");
          }
          
          double randeff_contrib = randeff ? allfit_random[k] : 0.0;
          
          double r = (y[k] - allfit_con[k] - randeff_contrib)*bscale/allfit_mod[k];
          
          if(r!=r) {
            Rcout << "bscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
            stop("");
          }
          if(k<ntrt) {
            ww1 += scale_factor;
            rw1 += r*scale_factor;
          } else {
            ww0 += scale_factor;
            rw0 += r*scale_factor;
          }
        }
        logger.log("Drawing bscale 1");
        logger.startContext();
        double bscale1_old = bscale1;
        double bscale_fc_var = 1/(ww1 + bscale_prec);
        bscale1 = bscale_fc_var*rw1 + gen.normal(0., 1.)*sqrt(bscale_fc_var);
        if(verbose_itr){
          
          Rcout << "Original bscale1 : " << bscale1_old << "\n";
          Rcout << "bscale_prec : " << bscale_prec << ", ww1 : " << ww1 << ", rw1 : " << rw1 << "\n";
          Rcout << "New  bscale1 : " << bscale1 << "\n\n";
        }
        logger.stopContext();
        
        //Error probably here. 
        logger.log("Drawing bscale 0");
        logger.startContext();
        double bscale0_old = bscale0;
        bscale_fc_var = 1/(ww0 + bscale_prec);
        bscale0 = bscale_fc_var*rw0 + gen.normal(0., 1.)*sqrt(bscale_fc_var);
        if(verbose_itr){
          Rcout << "Original bscale0 : " << bscale0_old << "\n";
          Rcout << "bscale_prec : " << bscale_prec << ", ww0 : " << ww0 << ", rw0 : " << rw0 << "\n";
          Rcout << "New  bscale0 : " << bscale0 << "\n\n";
        }
        logger.stopContext();
        if(Rcpp::NumericVector::is_na(bscale1) || Rcpp::NumericVector::is_na(bscale0)) {
          Rcout << "Original bscale1 : " << bscale1_old << "\n";
          Rcout << "bscale_prec : " << bscale_prec << ", ww1 : " << ww1 << ", rw1 : " << rw1 << "\n";
          Rcout << "New  bscale1 : " << bscale1 << "\n\n";
          Rcout << "Original bscale0 : " << bscale0_old << "\n";
          Rcout << "bscale_prec : " << bscale_prec << ", ww1 : " << ww1 << ", rw1 : " << rw1 << "\n";
          Rcout << "New  bscale0 : " << bscale0 << "\n\n";
        }
        for(size_t k=0; k<ntrt; ++k) {
          allfit_mod[k] = allfit_mod[k]*bscale1/bscale1_old;
        }
        for(size_t k=ntrt; k<n; ++k) {
          allfit_mod[k] = allfit_mod[k]*bscale0/bscale0_old;
        }
        
        // if(!b_half_normal) {
        //   double ssq = 0.0;
        //   tree::npv bnv;
        //   typedef tree::npv::size_type bvsz;
        //   double endnode_count = 0.0;
        //   
        //   for(size_t iTreeMod=0;iTreeMod<ntree_mod;iTreeMod++) {
        //     bnv.clear();
        //     t_mod[iTreeMod].getbots(bnv);
        //     bvsz nb = bnv.size();
        //     for(bvsz ii = 0; ii<nb; ++ii) {
        //       double mm = bnv[ii]->getm(); //node parameter
        //       ssq += mm*mm/(pi_mod.tau*pi_mod.tau);
        //       endnode_count += 1.0;
        //     }
        //   }
        //   delta_mod = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
        // }
        if(verbose_itr){
          Rcout << "Original pi_mod.tau : " <<  pi_mod.tau << "\n";
        }
        
        pi_mod.tau   = mod_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod));
        
        if(verbose_itr){
          Rcout << "New pi_mod.tau : " <<  pi_mod.tau << "\n\n";
        }
        
      } else {
        bscale0 = -0.5;
        bscale1 =  0.5;
      }
      pi_mod.sigma = sigma;
      
      
      if(use_mscale) {
        double ww = 0.;
        double rw = 0.;
        double s2 = sigma*sigma;
        for(size_t k=0; k<n; ++k) {
          double scale_factor = (w[k]*allfit_con[k]*allfit_con[k])/(s2*mscale*mscale);
          if(scale_factor!=scale_factor) {
            Rcout << " scale_factor " << scale_factor << endl;
            stop("");
          }
          
          double randeff_contrib = randeff ? allfit_random[k] : 0.0;
          
          double r = (y[k] - allfit_mod[k]- randeff_contrib)*mscale/allfit_con[k];
          if(r!=r) {
            Rcout << "allfit_mod " << allfit_mod[k] << endl;
            Rcout << "mscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
            std::ofstream debug_file("debug_output.txt");
            if (debug_file.is_open()) {
              debug_file << "Error at index k = " << k << "\n";
              debug_file << "allfit_mod[" << k << "] = " << allfit_mod[k] << "\n";
              debug_file << "mscale = " << mscale << "\n";
              debug_file << "bscale0 = " << bscale0 << "\n";
              // debug_file << "ww1" << ww1 << "\n";
              // debug_file << "ww0" << ww0 << "\n";
              // debug_file << "bscale_fc_var" << bscale_fc_var << "\n";
              // debug_file << "bscale_prec" << bscale_prec << "\n";
              // debug_file << "bscale1_old" << bscale1_old << "\n";
              // debug_file << "bscale0_old" << bscale0_old << "\n";
              
              
              debug_file << "bscale1 = " << bscale1 << "\n";
              debug_file << "allfit_con[" << k << "] = " << allfit_con[k] << "\n";
              debug_file << "y[" << k << "] = " << y[k] << "\n";
              debug_file << "alpha" << alpha << "\n";
              debug_file << "beta = [ ";
              for(size_t i = 0; i < beta.size(); i++) {
                debug_file << beta[i] << " ";
              }
              debug_file << "]\n";
              
              debug_file.close();
            } else {
              Rcout << "Error opening debug_output.txt file!" << endl;
            }
            stop("");
          }
          ww += scale_factor;
          rw += r*scale_factor;
        }
        
        logger.log("Drawing mscale");
        
        
        double mscale_old = mscale;
        double mscale_fc_var = 1/(ww + mscale_prec);
        mscale = mscale_fc_var*rw + gen.normal(0., 1.)*sqrt(mscale_fc_var);
        if(verbose_itr){
          Rcout << "Original mscale : " << mscale_old << "\n";
          Rcout << "mscale_prec : " << mscale_prec << ", ww : " << ww << ", rw : " << rw << "\n";
          Rcout << "New  mscale : " << mscale << "\n\n";
        }
        
        
        //Rcout<< mscale_fc_var << " " << rw <<" " << mscale << endl;
        
        for(size_t k=0; k<n; ++k) {
          allfit_con[k] = allfit_con[k]*mscale/mscale_old;
        }
        
        // update delta_con
        
        double ssq = 0.0;
        tree::npv bnv;
        typedef tree::npv::size_type bvsz;
        double endnode_count = 0.0;
        
        for(size_t iTreeCon=0;iTreeCon<ntree_con;iTreeCon++) {
          bnv.clear();
          t_con[iTreeCon].getbots(bnv);
          bvsz nb = bnv.size();
          for(bvsz ii = 0; ii<nb; ++ii) {
            double mm = bnv[ii]->getm(); //node parameter
            ssq += mm*mm/(pi_con.tau*pi_con.tau);
            endnode_count += 1.0;
          }
        }
        
        delta_con = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
        if(verbose_itr){
          logger.log("Updating pi_con.tau");
          Rcout << "Original pi_con.tau : " <<  pi_con.tau << "\n";
        }
        
        pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con));
        
        if(verbose_itr){
          Rcout << "New pi_con.tau : " <<  pi_con.tau << "\n\n";
        }
        
        
      } else {
        mscale = 1.0;
      }
      pi_con.sigma = sigma/fabs(mscale); //should be sigma/abs(mscale) for backfitting
      
      //sync allfits after scale updates, if necessary. Could do smarter backfitting updates inline
      if(use_mscale || use_bscale) {
        logger.log("Sync allfits after scale updates");
        
        for(size_t k=0; k<n; ++k) {
          double randeff_contrib = randeff ? allfit_random[k] : 0.0;
          allfit[k] = allfit_con[k] + allfit_mod[k] + randeff_contrib;
        }
      }
      
      if(randeff) {
        Rcout << "==================================\n";
        Rcout << "- Random Effects \n";
        Rcout << "==================================\n";
        
        //update random effects
        for(size_t k=0; k<n; ++k) {
          r(k) = y[k] - allfit_con[k] - allfit_mod[k];
          allfit[k] -= allfit_random[k];
        }
        
        Wtr = random_des.t()*r;
        
        arma::mat adj = diagmat(random_var_ix*eta);
        //    Rcout << adj << endl << endl;
        arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
        arma::vec m = adj*Wtr/(sigma*sigma);
        //Rcout << m << Phi << endl << Sigma_inv_random;
        gamma = rmvnorm_post(m, Phi);
        
        //Rcout << "updated gamma";
        
        // Update px parameters eta
        
        arma::mat adj2 = diagmat(gamma)*random_var_ix;
        arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
        arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
        //Rcout << m << Phi << endl << Sigma_inv_random;
        eta = rmvnorm_post(m2, Phi2);
        
        //Rcout << "updated eta";
        
        // Update variance parameters
        
        arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
        //Rcout << "A";
        arma::rowvec counts = sum(random_var_ix, 0);
        //Rcout << "B";
        for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
          random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df + ssqs(ii));
        }
        //Rcout << "updated vars" << endl;
        Sigma_inv_random = diagmat(1/(random_var_ix*random_var));
        
        allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;
        
        //Rcout << "recom allfit vars" << endl;
        
        for(size_t k=0; k<n; ++k) {
          allfit[k] = allfit_con[k] + allfit_mod[k] + allfit_random(k); //+= allfit_random[k];
        }
      }
      
      // ---------------------------------------------------------
      logger.log("Draw Sigma");
      // ---------------------------------------------------------
      double rss = 0.0;
      double restemp = 0.0;
      for(size_t k=0;k<n;k++) {
        restemp = y[k]-allfit[k];
        rss += w[k]*restemp*restemp;
      }
      if(!hamiltonian){
        sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
      } else {
        sigma = sqrt(std::exp(param[p_mod + p_int + p_mod + 2]));
      }
      pi_con.sigma = sigma/fabs(mscale);
      pi_mod.sigma = sigma; // Is this another copy paste Error?
      
      // sigma = sqrt(std::exp(init_param[p_mod + p_int + p_mod + 2]));
      
      if( ((iIter>=burn) & (iIter % thin==0)) )  {
        if(not treef_con_name.empty()){
          for(size_t j=0;j<ntree_con;j++) treef_con << std::setprecision(save_tree_precision) << t_con[j] << endl; // save trees
          // for(size_t j=0;j<ntree_mod;j++) treef_mod << std::setprecision(save_tree_precision) << t_mod[j] << endl; // save trees
        }
        
        msd_post(save_ctr) = mscale;
        bsd_post(save_ctr) = bscale1-bscale0;
        b0_post(save_ctr)  = bscale0;
        b1_post(save_ctr)  = bscale1;
        
        
        gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
        random_var_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();
        
        sigma_post(save_ctr) = sigma;
        for(size_t k=0;k<n;k++) {
          m_post(save_ctr, k) = allfit_con[k];
          yhat_post(save_ctr, k) = allfit[k];
        }
        for(size_t k=0;k<n;k++) {
          double bscale = (k<ntrt) ? bscale1 : bscale0;
          b_post(save_ctr, k) = (bscale1-bscale0)*allfit_mod[k]/bscale;
        }
        //}
        save_ctr += 1;
      }
      logger.log("==============================================");
      Rprintf(logBuff, "MCMC iteration: %d of %d End", iIter + 1, nd*thin+burn);
      logger.log(logBuff);
      Rprintf(logBuff, "sigma %f, mscale %f, bscale0 %f, bscale1 %f",sigma, mscale, bscale0, bscale1);
      logger.log(logBuff);
      logger.log("==============================================");
      if (verbose_itr){
        logger.getVectorHead(y, logBuff);
        Rcout << "           y: " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit, logBuff);
        Rcout << "Current Fit : " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit_con, logBuff);
        Rcout << "allfit_con  : " <<  logBuff << "\n";
        
        logger.getVectorHead(allfit_mod, logBuff);
        Rcout << "allfit_mod  : " <<  logBuff << "\n";
      }
      
    } // end MCMC Loop
    
    int time2 = time(&tp);
    Rcout << "\n============================================================\n MCMC Complete \n============================================================\n";
    
    Rcout << "time for loop: " << time2 - time1 << endl;
    
    // t_mod.clear(); t_con.clear();
    delete[] allfit;
    delete[] allfit_mod;
    delete[] allfit_con;
    delete[] r_mod;
    
    delete[] r_con;
    delete[] ftemp;
    
    R_FlushConsole(); // Flush the console so to not overload it with all the messages printed.
    
    if(not treef_con_name.empty()){
      treef_con.close();
      // treef_mod.close();
    }
    
    if(hamiltonian){
      alphaOut = samples.col(0);
      betaOut = wrap(samples.cols(1,p_mod));
      beta_intOut = wrap(samples.cols(p_mod + 1, p_mod + p_int));
      tauOut = wrap(samples.cols(p_mod + p_int + 1, p_mod + p_int + p_mod));
      tau_int_post = samples.col(p_mod + p_int + p_mod + 1);
      sigma_post = samples.col(p_mod + p_int + p_mod + 2);
      
    }
    
    return(List::create(_["yhat_post"] = yhat_post, _["m_post"] = m_post, _["b_post"] = b_post,
                        _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post, _["b0"] = b0_post, _["b1"] = b1_post,
                        _["gamma"] = gamma_post, _["random_var_post"] = random_var_post, _["Beta"] = betaOut, _["tau"] = tauOut, _["alpha"] = alphaOut, _["tau_int"] = tau_int_post, _["beta_int"] = beta_intOut, _["acceptance_rate"] = acceptance/nd
    ));
}