# NORMAL / EXPONENTIAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_exp_check <- function(params,
                             indeces,
                             y, X,
                             CM = NULL,
                             CV_u,
                             CV_v) {

  # extract parameters from parameter vector
  n_f_coeff <- ncol(X)
  n_cv_u_coeff <- if (is.matrix(CV_u)) ncol(CV_u) else length(CV_u)  # number of coeffs for conditional variance of inefficiency term model
  n_cv_v_coeff <- if (is.matrix(CV_v)) ncol(CV_v) else length(CV_v) # number of coeffs for conditional variance of the symmetric error model
  # no CM in exponential model

  if (length(params) != n_f_coeff + n_cv_u_coeff + n_cv_v_coeff) {
    stop("Incorrect nuber of parameters. ",
         n_f_coeff, "+", n_cv_u_coeff, "+", n_cv_v_coeff, " needed, but ", length(params), " supplied:", paste(names(params), collapse = " "))
  } else {
    return(TRUE)
  }
}

t_par_cs_exp <- function(pars){
  pars <- sqrt(exp(pars))
  names(pars) <- c("sigma_u", "sigma_v")
  return(pars)
}


# likelihood --------------------------------------------------------------


# likelihood function normal/exponential distributional assumption
# params: beta, log(sigma_u^2), log(sigma_v^2)

ll_cs_exp_contrib <- function(params,
                              indeces,
                              y, X,
                              CM = NULL,
                              CV_u,
                              CV_v,
                              ineff,
                              deb) {

  if (deb) {
    cat("Parameters: ", params, "\n")
  }

  f_coeff    <- params[              1 :indeces[1]]
  cv_u_coeff <- params[(indeces[1] + 1):indeces[2]]
  cv_v_coeff <- params[(indeces[2] + 1):indeces[3]]

  sigma_u <- as.vector(exp(CV_u %*% cv_u_coeff))
  sigma_v <- as.vector(exp(CV_v %*% cv_v_coeff))

  sigma2_u <- sigma_u^2
  sigma2_v <- sigma_v^2

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(f_coeff), "\n",
               "Sigma_u: ", cv_u_coeff, "\n",
               "Sigma_v: ", cv_v_coeff, "\n")

  epsilon <- as.vector(-ineff * (y - X %*% f_coeff))

  N <- length(y)

  lli <- -log(sigma_u) + 0.5 * sigma2_v / sigma2_u +
    log(pnorm(-(epsilon + (sigma2_v / sigma_u)) / sigma_v)) +
    epsilon / sigma_u

  return(lli)
}


# main loklikelihood function
ll_cs_exp <- function(params,
                      indeces,
                      y, X,
                      CM = NULL,
                      CV_u,
                      CV_v,
                      ineff,
                      minmax = -1,
                      deb = FALSE) {

  # compute loglik contributions
  lli <-
    ll_cs_exp_contrib(params,
                      indeces,
                      y, X,
                      CM = NULL,
                      CV_u,
                      CV_v,
                      ineff,
                      deb)

  if (deb) cat("Misbehaving loglikelihood contributions: ", "\n",
               which(!is.finite(lli)), "\n",
               lli[!is.finite(lli)],  "\n")

  if (deb) cat("Suspicious loglikelihood contributions: ", "\n",
               which(lli < quantile(lli, 0.05)), "\n",
               lli[lli < quantile(lli, 0.05)],  "\n")

  # sum to total loglikelihood
  ll <- sum(lli)
  if (deb) cat("Loglikelihood: ", ll,  "\n")

  return(minmax*ll) # return -loglikelihood for optimization of minimum
}
if (require(compiler)) ll_cs_exp <- cmpfun(ll_cs_exp)




# gradient functions -----------------------------------------

g_cs_exp_fd <- function(params,
                        indeces,
                        y, X,
                        CV_u,
                        CV_v,
                        CM,
                        ineff,
                        minmax,
                        deb) {
  n <- length(params)
  hh <- matrix(0, n, n)
  diag(hh) <- .Machine$double.eps^(1/3)

  sapply(1:n, function(i) {
    (  ll_cs_exp(params + hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb) -
       ll_cs_exp(params - hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb)) / (2 * .Machine$double.eps^(1/3))})
}
if (require(compiler)) g_cs_exp_fd <- cmpfun(g_cs_exp_fd)

g_cs_exp_analytic <- function(params,
                              indeces,
                              y, X,
                              CM = NULL,
                              CV_u,
                              CV_v,
                              ineff,
                              minmax,
                              deb = F) {

  if (deb) {
    cat("Parameters: ", params)
  }

  f_coeff    <- params[              1 :indeces[1]]
  cv_u_coeff <- params[(indeces[1] + 1):indeces[2]]
  cv_v_coeff <- params[(indeces[2] + 1):indeces[3]]

  sigma_u <- as.vector(exp(CV_u %*% cv_u_coeff))
  sigma_v <- as.vector(exp(CV_v %*% cv_v_coeff))

  sigma2_u <- sigma_u^2
  sigma2_v <- sigma_v^2

  eps <- as.vector(-ineff * (y - X %*% f_coeff))

  N <- length(y)

  exp1 <- (-eps - sigma2_v/sigma_u)/sigma_v
  pdfnorm <- dnorm(exp1)
  cdfnorm <- pnorm(exp1)

  g_f_coeff <- -ineff*(pdfnorm/cdfnorm/sigma_v - 1/sigma_u) %*% X
  g_cv_u_coeff <- (-1 - sigma2_v/sigma2_u + (pdfnorm/cdfnorm) * sigma_v/sigma_u - eps/sigma_u) %*% CV_u
  g_cv_v_coeff <- (sigma2_v/sigma2_u + pdfnorm/cdfnorm * (eps/sigma_v - sigma_v/sigma_u)) %*% CV_v

  return(minmax*c(g_f_coeff, g_cv_u_coeff, g_cv_v_coeff))
}



# efficiency --------------------------------------------------------------


u_cs_exp <- function(object, estimator) {

  # extract sigmas from the model object
  sigma_u <- as.vector(exp(object$data$CV_u %*% object$coeff_cv_u))
  sigma_v <- as.vector(exp(object$data$CV_v %*% object$coeff_cv_v))

  # derive the rest of the parameters
  sigma2_u <- sigma_u^2
  sigma2_v <- sigma_v^2

  sigma2 <- sigma2_u + sigma2_v

  mu_ast <- object$ineff * object$residuals - sigma2_v / sigma_u
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_v*dnorm(mu_ast/sigma_v)/pnorm(mu_ast/sigma_v),
              ... = stop(paste0("Unknown type ", estimator," of conditional mean estimator.")))

  return(u)
}
