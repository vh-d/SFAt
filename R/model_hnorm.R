# NORMAL / H-NORMAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_hnorm_check <- function(params,
                               indeces,
                               y, X,
                               CM = NULL,
                               CV_u,
                               CV_v) {

  # extract parameters from parameter vector
  n_f_coeff <- ncol(X)
  n_cv_u_coeff <- if (is.matrix(CV_u)) ncol(CV_u) else length(CV_u)  # number of coeffs for conditional variance of inefficiency term model
  n_cv_v_coeff <- if (is.matrix(CV_v)) ncol(CV_v) else length(CV_v) # number of coeffs for conditional variance of the symmetric error model
  # no CM in half-normal model

  if (length(params) != n_f_coeff + n_cv_u_coeff + n_cv_v_coeff) {
    stop("Incorrect nuber of parameters. ",
         n_f_coeff, "+", n_cv_u_coeff, "+", n_cv_v_coeff, " needed, but ", length(params), " supplied:", paste(names(params), collapse = " "))
  } else {
    return(TRUE)
  }
}

t_par_cs_hnorm <- function(pars){
  pars <- sqrt(exp(pars))
  names(pars) <- c("sigma_u", "sigma_v")
  return(pars)
}

# likelihood function normal/half-normal distributional assumption
# params: beta, log(sigma_u^2), log(sigma_v^2)

ll_cs_hnorm <- function(params,
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

  sigma2 <- sigma2_u + sigma2_v
  sigma <- sqrt(sigma2)

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(f_coeff), "\n",
               "Sigma_u: ", sigma_u,
               "Sigma_v: ", sigma_v, "\n")

  epsilon <- as.vector(-ineff * (y - X %*% f_coeff))

  N <- length(y)

  lli <- 0.5*log(2/pi) -log(sigma) + log(pnorm(-(epsilon * sigma_u / sigma_v) / sigma)) - (epsilon^2) / (2*sigma2)

  ll <- sum(lli)
  if (deb) cat("Loglikelihood: ", ll,  "\n")
  if (deb) print(summary(epsilon))
  if (deb) print(summary(lli))

  # if (!is.finite(ll) && minmax == -1) {
  #   # return(sum(!is.finite(lli))*1e100)
  #   return(1e150)
  # }

  return(minmax*ll)
}

if (require(compiler)) ll_cs_hnorm <- cmpfun(ll_cs_hnorm)

# gradient function
g_cs_hnorm_analytic <- function(params,
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

  sigma2 <- sigma2_u + sigma2_v
  sigma <- sqrt(sigma2)

  lambda <- sigma_u/sigma_v

  eps <- as.vector(-ineff * (y - X %*% f_coeff))

  N <- length(y)

  sus <- sigma2_u/sigma2
  svs <- sigma2_v/sigma2
  epsigma <- eps/sigma
  cdfnorm <- pnorm(lambda*epsigma)
  lpdfnorm <- lambda*dnorm(epsigma*lambda)

  g_f_coeff <- -ineff*(eps / sigma2 + lpdfnorm/(sigma*(1-cdfnorm))) %*% X
  g_cv_u_coeff <- (-sus + sus*(epsigma)^2 - ((lpdfnorm)/(1-cdfnorm))*(epsigma)*(1 - sus)) %*% CV_u
  g_cv_v_coeff <- (-svs + svs*(epsigma)^2 + ((lpdfnorm)/(1-cdfnorm))*(epsigma)*(1 + svs)) %*% CV_v

  return(minmax*c(g_f_coeff, g_cv_u_coeff, g_cv_v_coeff))
}

if (require(compiler)) g_cs_hnorm_analytic <- cmpfun(g_cs_hnorm_analytic)

# gradient function
g_cs_hnorm_fd <- function(params,
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
  # heps <- .Machine$double.eps^(1/3)
  heps <- 1e-12
  diag(hh) <- heps
  heps2 <- 2 * heps

  sapply(1:n, function(i) {
    (  ll_cs_exp(params + hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb) -
         ll_cs_exp(params - hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb)) / heps2})
}

if (require(compiler)) g_cs_hnorm_fd <- cmpfun(g_cs_hnorm_fd)

# Jondrow et al. (1982) as in Parmeter-Kumbhakar
u_cs_hnorm <- function(object, estimator) {

   # extract sigmas from the model object
  sigma2_u <- as.vector(exp(object$data$CV_u %*% object$coeff_cv_u))
  sigma2_v <- as.vector(exp(object$data$CV_v %*% object$coeff_cv_v))

  # derive the rest of the parameters
  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_u + sigma2_v

  mu_ast <- object$ineff * object$residuals * sigma2_u / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)),
              ... = stop(paste0("Unknown type ", estimator," of conditional mean estimator.")))

  return(u)
}