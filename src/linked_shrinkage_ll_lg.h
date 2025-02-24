#ifndef LINKED_SHRINKAGE_LL_LG_H
#define LINKED_SHRINKAGE_LL_LG_H

#include <RcppArmadillo.h>

// Ensure M_PI is defined if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Namespace declaration to avoid repetition
using namespace Rcpp;
using namespace arma;

static arma::vec compute_eta(const arma::vec &param,
                             const arma::mat &X);

double log_posterior_linked_shrinkage(const arma::vec &param,
                                      const arma::mat &X,
                                      const arma::vec &y);

static arma::vec compute_residual(const arma::vec &param,
                                  const arma::mat &X,
                                  const arma::vec &y);

arma::vec grad_log_posterior_linked_shrinkage(const arma::vec &param,
                                              const arma::mat &X,
                                              const arma::vec &y);
Rcpp::List leapfrogCpp(
    arma::vec param,
    arma::vec momentum,
    double step_size,
    int num_steps,
    const arma::mat &X,
    const arma::vec &y
);


#endif