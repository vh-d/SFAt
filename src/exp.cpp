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
double ll_cs_exp_cpp(const NumericVector& pars,
                     const IntegerVector& indeces,
                     const NumericVector& y,
                     const NumericMatrix& X,
                     const SEXP& CV_u,
                     const SEXP& CV_v,
                     const SEXP& CM,
                     const int ineff, const bool deb) {

  const NumericVector f_coefs(pars.begin(),
                              pars.begin() + indeces[0]);
  const NumericVector cv_u_coefs(pars.begin() + indeces[0],
                                 pars.begin() + indeces[1]);
  const NumericVector cv_v_coefs(pars.begin() + indeces[1],
                                 pars.end());

  const int nbetas = X.ncol();
  const int N_obs = y.size();
  NumericVector epsilon(N_obs);
  int i, j = 0;

  double ll = 0.0;
  double tmp = 0.0;

  if (Rf_isMatrix(CV_u) & Rf_isMatrix(CV_v)) {
    if (deb) {
      Rcout << "-------- Matrices ---------" << std::endl;
    }

    NumericMatrix CVu(CV_u);
    NumericMatrix CVv(CV_v);

    NumericVector sigma2_u = exp(CVu * cv_u_coefs);
    NumericVector sigma2_v = exp(CVv * cv_v_coefs);

    for (i = 0; i < N_obs; i++) {
      sigma2_u[i] = exp(sum(CVu(i, _) * cv_u_coefs));
      sigma2_v[i] = exp(sum(CVv(i, _) * cv_v_coefs));
    }

    const NumericVector sigma_u = sqrt(sigma2_u);
    const NumericVector sigma_v = sqrt(sigma2_v);

    if (deb){
      Rcout << "Coeffs: " << f_coefs << std::endl;
      Rcout << "Sigma_u coeffs: " << cv_u_coefs << std::endl;
      Rcout << "Sigma_v coeffs: " << cv_v_coefs << std::endl;

      Rcout << "Sigma_u: " << sigma_u << std::endl;
      Rcout << "Sigma_v: " << sigma_v << std::endl;
      Rcout << "N: " << N_obs << std::endl;
    }

    for (i = 0;i < N_obs; i++) {
      tmp = sum(X(i,_) * f_coefs);
      epsilon[i] = -ineff*(y[i] - tmp);

      ll +=
        -log(sigma_u[i]) +
        0.5 * sigma2_v[i] / sigma2_u[i] +
        log(R::pnorm(-(epsilon[i] + (sigma2_v[i] / sigma_u[i])) / sigma_v[i], 0.0, 1.0, 1, 0)) +
        (epsilon[i] / sigma_u[i]);
    }

  } else {
    if (deb) {
      Rcout << "-------- Scalars ------" << std::endl;
    }

    double CVu = as<double>(CV_u);
    double CVv = as<double>(CV_v);
    const double sigma2_u = exp(CVu * pars[indeces[0]]);
    const double sigma2_v = exp(CVv * pars[indeces[1]]);
    const double sigma_u = sqrt(sigma2_u);
    const double sigma_v = sqrt(sigma2_v);

    if (deb){
      Rcout << "Coeffs: " << f_coefs << std::endl;
      Rcout << "Sigma_u coeffs: " << cv_u_coefs << std::endl;
      Rcout << "Sigma_v coeffs: " << cv_v_coefs << std::endl;

      Rcout << "Sigma_u: " << sigma_u << std::endl;
      Rcout << "Sigma_v: " << sigma_v << std::endl;
      Rcout << "N: " << N_obs << std::endl;
    }

    for (i = 0;i < N_obs; i++) {
      tmp = sum(X(i,_) * f_coefs);

      epsilon[i] = -ineff*(y[i] - tmp);

      ll += -1*log(sigma_u)
        + 0.5 * sigma2_v / sigma2_u
        + log(R::pnorm(-(epsilon[i] + (sigma2_v / sigma_u)) / sigma_v, 0.0, 1.0, 1, 0))
        + (epsilon[i] / sigma_u);
    }
  }

  if (deb) {
    Rcout << std::setprecision(20) << -ll << std::endl;
  }

  if (!R_finite(ll)) return(-1e100);

  return ll;
}

