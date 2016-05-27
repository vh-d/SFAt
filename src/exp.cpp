#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
using namespace Rcpp;

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

double normalCFD(double value)
{
  return 0.5 * erfc(-value * M_SQRT1_2);
}

// Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double ll_cs_exp(SEXP& params, SEXP& y, SEXP& X, int ineff, bool deb) {

  NumericVector pars(params);
  NumericVector yy(y);
  NumericMatrix XX(X);

  int nbetas = XX.ncol();

  double lnsigma2_u = pars[nbetas];
  double lnsigma2_v = pars[nbetas + 1];

  double sigma2_u = std::exp(lnsigma2_u);
  double sigma2_v = std::exp(lnsigma2_v);

  double sigma_u = std::sqrt(sigma2_u);
  double sigma_v = std::sqrt(sigma2_v);

  const int N = yy.size();

  // for (int k = 0; k < nbetas; k++) {
  //   Rcout << "Coeff. " << k << " :" << pars[k] << std::endl;
  // }

  // Rcout << "Sigma_u: " << sigma_u << std::endl;
  // Rcout << "Sigma_v: " << sigma_v << std::endl;

  double ll = -N * std::log(sigma_u) + 0.5 * N * (sigma2_v / sigma2_u);

  // Rcout << "Log-lik: " << ll << std::endl;

  NumericVector epsilon(N);

  int i,j = 0;

  double tmp = 0.0;
  for (i = 0;i < N; ++i) {

    tmp = 0.0;

    for (j = 0;j < nbetas; ++j){
      tmp += XX(i, j)*pars[j];
    }

    epsilon[i] = -ineff*(yy[i] - tmp);

    // Rcout << epsilon[i] << std::endl;

    ll += std::log(normalCFD(-(epsilon[i] + (sigma2_v / sigma_u)) / sigma_v)) + (epsilon[i] / sigma_u);
  }

  return -ll;
}