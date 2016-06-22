# NORMAL / T-NORMAL MODEL- CROSS SECTION DATA ----------------------------

par_cs_tnorm <- NULL

par_cs_tnorm_check <- function(params,
                               indeces,
                               y, X,
                               CM = NULL,
                               CV_u,
                               CV_v) {

  # extract parameters from parameter vector
  n_f_coeff <- ncol(X) # number of coeffs for frontier model
  n_cv_u_coeff <- if (is.matrix(CV_u)) ncol(CV_u) else length(CV_u)  # number of coeffs for conditional variance of inefficiency term model
  n_cv_v_coeff <- if (is.matrix(CV_v)) ncol(CV_v) else length(CV_v) # number of coeffs for conditional variance of the symmetric error model
  n_cm_coeff <- if (is.matrix(CM)) ncol(CM) else length(CM) # number of coeffs for conditional mean model

  if (length(params) != n_f_coeff + n_cm_coeff + n_cv_u_coeff + n_cv_v_coeff) {
    stop("Incorrect nuber of parameters. ",
         n_f_coeff, "+", n_cm_coeff, "+", n_cv_u_coeff, "+", n_cv_v_coeff, " needed, but ", length(params), " supplied.")
  }

  return(TRUE)
}

t_par_cs_tnorm <- function(pars){
  # pars <- sqrt(exp(pars))
  # names(pars) <- c("sigma_v")
  # return(pars)
  return(NULL)
}

# Advanced version of (Battese and Coelli, 1995) and (Huang and Liu, 1994) models
# heterogeneity in efficiency term: endogeneous location parameter mu
# implemented as Hadri et al. 2003
# parameters: f_coeff, cm_coeff, cv_u_coeff, cv_v_coeff
ll_cs_tnorm <- function(params,
                        indeces,
                        y, X,
                        CV_u,
                        CV_v,
                        CM,
                        ineff,
                        deb) {

  if (deb) {
    cat("Parameters: ", params)
  }

  f_coeff    <- params[             1:indeces[1]]
  cv_u_coeff <- params[(indeces[1] + 1):indeces[2]]
  cv_v_coeff <- params[(indeces[2] + 1):indeces[3]]
  cm_coeff   <- params[(indeces[3] + 1):indeces[4]]

  sigma2_u <- as.vector(exp(CV_u %*% cv_u_coeff))
  sigma2_v <- as.vector(exp(CV_v %*% cv_v_coeff))

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_v + sigma2_u # variance of composite error term
  sigma <- sqrt(sigma2)
  sigma_ast <- sigma_u * sigma_v / sigma

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "frontier coeffs: ", paste(f_coeff), "\n",
               "cm coeffs: ", paste(cm_coeff), "\n",
               "cv_u coeffs: ", paste(cv_u_coeff), "\n",
               "cv_v coeffs: ", paste(cv_v_coeff), "\n",
               "sigma_u: ", sigma_u, "\n",
               "sigma_v: ", sigma_v, "\n")

  N <- length(y)

  Zdelta <- as.vector(CM %*% cm_coeff) # fitted means of inefficiency term

  eps <- -ineff*as.vector(y - (X %*% f_coeff)) # composite error terms

  mu_ast <- (sigma2_v*Zdelta - sigma2_u*eps) / sigma2

  if (deb) {
    cat(length(Zdelta) == N, "\n",
        length(eps) == N, "\n",
        length(mu_ast) == N, "\n",
        length(sigma_ast) == N, "\n",
        length(sigma2) == N, "\n",
        length(sigma_u) == N, "\n")
  }

  result <- sum(-0.5*log(2*pi) - log(sigma) - 0.5 * ((eps + Zdelta)^2)/sigma2 +
                  log(pnorm(mu_ast / sigma_ast)) - log(pnorm(Zdelta / sigma_u)))

  if (deb) cat("Loglikelihood: ", result,  "\n")

  if (!is.finite(result)) return(-1e100)
  return(result)
}

u_cs_tnorm <- function(object, estimator) {
  # extract sigmas from the model object
  sigma2_u <- as.vector(exp(object$data$CV_u %*% object$coeff_cv_u))
  sigma2_v <- as.vector(exp(object$data$CV_v %*% object$coeff_cv_v))

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)
  sigma2 <- sigma2_u + sigma2_v

  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  mu <- as.vector(object$data$CM %*% object$coeff_cm)
  mu_ast <- (object$ineff * object$residuals * sigma2_u + sigma2_v * mu) / sigma2

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)),
              ... = stop(paste0("Unknown type ", estimator," of conditional mean estimator.")))

  return(u)
}