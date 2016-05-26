#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ll_cs_exp(NumericVector& params, NumericVector& y, NumericMatrix& X, int ineff, bool deb) {

  int nbetas = X.ncol();

  if (params.size() != nbetas + 2) {
    stop("Incorrect nuber of parameters. ",
         nbetas + 2, " needed, but ", params.size(), " supplied.");
  }

  double lnsigma2_u = params[nbetas + 1];
  double lnsigma2_v = params[nbetas + 2];

  double sigma2_u = exp(lnsigma2_u);
  double sigma2_v = exp(lnsigma2_v);

  double sigma_u = sqrt(sigma2_u);
  double sigma_v = sqrt(sigma2_v);

  // if (deb) cat("Total of ", length(params), " parameters: \n",
  //     "Betas: ", paste(params[Range(0, nbetas)]), "\n",
  //     "Sigma_u: ", sigma_u,
  //     "Sigma_v: ", sigma_v, "\n")

  NumericVector epsilon = -ineff * (y - (X*trans(params[Range(0, nbetas)])));

  const int N = y.size();

  double ll =
    - N * log(sigma_u) +
    0.5 * N * (sigma2_v / sigma2_u) +
    sum(log(pnorm(-(epsilon + (sigma2_v / sigma_u)) / sigma_v))) +
    sum(epsilon) / sigma_u;

  return(-ll);
}
