#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// #ifndef Pi
// #define Pi 3.141592653589793238462643
// #endif
//
// double normalCFD(double value)
// {
//   return 0.5 * erfc(-value * M_SQRT1_2);
// }


// [[Rcpp::export]]
NumericVector rtnorm_(int n, double mean, double sd, double lower, double upper) {
  if (lower >= upper) Rcpp::stop("Upper <= Lower !");

  NumericVector ans(n);

  int success = 0;
  double trial = 0.0;

  while (success < n) {
    trial = R::rnorm(mean, sd);

    if (trial > lower && trial < upper) {
      success++;
      ans[success-1] = trial;
    }
  }

  return(ans);
}