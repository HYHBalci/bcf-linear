
#define RCPP_NO_VISIBILITY 1
#undef attribute_visible
#define attribute_visible

#include <Rcpp.h>

// [[Rcpp::export]]
int test_visibility() {
  return 123;
}