# NORMAL / T-NORMAL MODEL- CROSS SECTION DATA ----------------------------

par_cs_tnorm_bc95_cv <- c(lnsigma2_v = 0.5)

t_par_cs_tnorm_bc95_cv <- function(pars){
  pars <- sqrt(exp(pars))
  names(pars) <- c("sigma_v")
  return(pars)
}

# Advanced version of (Battese and Coelli, 1995) and (Huang and Liu, 1994) models
# heterogeneity in efficiency term: endogeneous location parameter mu
# implemented as Hadri et al. 2003
# parameters: f_coeff, cm_coeff, cv_coeff, sigma2_v
ll_cs_tnorm_bc95_cv <- function(params, y, X, CM, CV, ineff, deb) {

  # extract parameters from parameter vector
  n_f_coeff <- ncol(X) # number of coeffs for frontier model
  n_cm_coeff <- ncol(CM) # number of coeffs for conditional mean model
  n_cv_coeff <- ncol(CV) # number of coeffs for conditional variance model

  if (length(params) != n_f_coeff + n_cm_coeff + n_cv_coeff + 1) {
    stop("Incorrect nuber of parameters. ",
         n_f_coeff, "+", n_cm_coeff, "+", n_cv_coeff, "+ 1 needed, but ", length(params), " supplied.")
  }

  f_coeff <- params[1:n_f_coeff]
  cm_coeff <- params[(n_f_coeff + 1):(n_f_coeff + n_cm_coeff)]
  cv_coeff <- params[(n_f_coeff + n_cm_coeff + 1):(n_f_coeff + n_cm_coeff + n_cv_coeff)]

  lnsigma2_v <- params[n_f_coeff + n_cm_coeff + n_cv_coeff + 1] # variance of symmetric error term

  sigma2_u <- as.vector(exp(CV %*% cv_coeff))
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_v + sigma2_u # variance of composite error term
  sigma <- sqrt(sigma2)
  sigma_ast <- sigma_u * sigma_v / sigma

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "frontier coeffs: ", paste(f_coeff), "\n",
               "cm coeffs: ", paste(cm_coeff), "\n",
               "cv coeffs: ", paste(cv_coeff), "\n",
               "Sigma2_v: ", sigma2_v, "\n")

  N <- length(y)

  Zdelta <- as.vector(CM %*% cm_coeff) # fitted means of inefficiency term

  eps <- -ineff*as.vector(y - (X %*% f_coeff)) # composite error terms

  mu_ast <- (sigma2_v*Zdelta - sigma2_u*eps) / sigma2

  temp_expr1 <- -0.5 * sum(log(2*pi) + log(sigma2))
  temp_expr2 <- -0.5 * sum(((eps + Zdelta)^2)/sigma2)
  temp_expr3 <- sum(log(pnorm(mu_ast / sigma_ast)) - log(pnorm(Zdelta / sigma_u)))

  if (deb) cat("Log-likelihood: ",
               temp_expr1, " + ", temp_expr2, " + ",
               temp_expr3, " + ",
               # temp_expr4,
               "\n")

  # result <- temp_expr1 + temp_expr2 + temp_expr3 + temp_expr4
  result <- temp_expr1 + temp_expr2 + temp_expr3

  if (deb) cat("Loglikelihood: ", result,  "\n")

  return(-result)
}

u_cs_tnorm_bc95_cv <- function(object, estimator) {
  # extract sigmas from model object
  lnsigma2_u <- object$parameters[length(object$coeff) + length(object$cm_coeff) + 1]
  lnsigma2_v <- object$parameters[length(object$coeff) + length(object$cm_coeff) + 2]

  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)
  sigma2 <- sigma2_u + sigma2_v

  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  mu <- as.vector(object$data$CM %*% object$cm_coeff)
  mu_ast <- (object$ineff * object$residuals * sigma2_u + sigma2_v * mu) / sigma2

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))
  if (is.null(u)) stop(paste0("Unknown type ", estimator," of conditional mean estimator."))

  return(u)
}