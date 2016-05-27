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
double ll_cs_exp_cpp(const SEXP params, const SEXP y, const SEXP X, const int ineff, const bool deb) {

  const NumericVector pars(params);
  const NumericVector yy(y);
  const NumericMatrix XX(X);

  int nbetas = XX.ncol();

  const double lnsigma2_u = pars[nbetas];
  const double lnsigma2_v = pars[nbetas + 1];

  const double sigma2_u = exp(lnsigma2_u);
  const double sigma2_v = exp(lnsigma2_v);

  const double sigma_u = sqrt(sigma2_u);
  const double sigma_v = sqrt(sigma2_v);

  const int N_obs = yy.size();

  if (deb){
    for (int k = 0; k < nbetas; k++) {
      Rcout << "Coeff. " << k << " :" << pars[k] << std::endl;
    }

    Rcout << "Sigma_u: " << sigma_u << std::endl;
    Rcout << "Sigma_v: " << sigma_v << std::endl;

    Rcout << "N: " << N_obs << std::endl;
  }

  double ll = -1 * N_obs * log(sigma_u) + 0.5 * N_obs * (sigma2_v / sigma2_u);

  // Rcout << "Log-lik: " << ll << std::endl;

  NumericVector epsilon(N_obs);

  int i,j = 0;
  for (i = 0;i < N_obs; i++) {

    double tmp = 0.0;
    double expr = 0.0;

    for (j = 0;j < nbetas; j++){
      tmp += XX(i, j)*pars[j];
    }

    epsilon[i] = -ineff*(yy[i] - tmp);

    expr = R::pnorm(-(epsilon[i] + (sigma2_v / sigma_u)) / sigma_v, 0.0, 1.0, 1, 0);

    if (expr == 0.0) {
      if (deb) Rcout << "Inf" << std::endl;
      return(R_PosInf);
    } else {
      ll += log(expr) + (epsilon[i] / sigma_u);
    }

  }

  if (deb) {
    Rcout << std::setprecision(20) << -ll << std::endl;
  }

  return -ll;
}
