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



# likelihood function normal/exponential distributional assumption
# params: beta, log(sigma_u^2), log(sigma_v^2)
ll_cs_exp <- function(params,
                      indeces,
                      y, X,
                      CM = NULL,
                      CV_u,
                      CV_v,
                      ineff,
                      minmax,
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

  ll <- sum(lli)

  if (!is.finite(ll) && minmax == -1) {
    return(sum(!is.finite(lli))*1e100)
  }

  return(minmax*ll)
}
if (require(compiler)) ll_cs_exp <- cmpfun(ll_cs_exp)

# gradient function
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


u_cs_exp <- function(object, estimator) {

  # extract sigmas from the model object
  sigma2_u <- as.vector(exp(object$data$CV_u %*% object$coeff_cv_u))
  sigma2_v <- as.vector(exp(object$data$CV_v %*% object$coeff_cv_v))

  # derive the rest of the parameters
  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

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
