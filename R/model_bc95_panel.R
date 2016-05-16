# NORMAL / T-NORMAL MODEL- CROSS SECTION DATA ----------------------------

par_panel_tnorm_bc95 <- c(lnsigma2_u = 0, lnsigma2_v = 0)

# Advanced version of (Battese and Coelli, 1995) and (Huang and Liu, 1994) models
# heterogeneity in efficiency term: endogeneous location parameter mu
# implemented as Hadri et al. 2003
# parameters: beta_coef, delta_coef, sigma2_v, sigma2_u
# K - matrix of panel data indeces, k - size of cross-section dimension
ll_panel_tnorm_bc95 <- function(params, y, X, CM, K, k, ineff, deb) {

  # extract parameters from parameter vector
  nbetas <- ncol(X) # number of beta coeffs
  ndeltas <- ncol(CM) # number of delta coeffs

  if (length(params) != nbetas + ndeltas + 2) {
    stop("Incorrect nuber of parameters. ",
         nbetas, "+", ndeltas, "+ 2 needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]
  delta_coef <- params[(nbetas + 1):(nbetas + ndeltas)]

  lnsigma2_u <- params[(nbetas + ndeltas + 1)] # variance of inefficiency term
  lnsigma2_v <- params[(nbetas + ndeltas + 2)] # variance of symmetric error term

  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_v + sigma2_u # variance of composite error term
  sigma <- sqrt(sigma2)
  sigma_ast <- sigma_u * sigma_v / sigma

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Deltas: ", paste(delta_coef), "\n",
               "Sigma2_u: ", sigma2_u,
               "Sigma2_v: ", sigma2_v, "\n")

  N <- length(y)

  Zdelta <- as.vector(CM %*% delta_coef) # fitted means of inefficiency term

  eps <- as.vector(y - (X %*% beta_coef)) # composite error terms

  mu_ast <- (sigma2_v*Zdelta - sigma2_u*eps) / sigma2

  temp_expr1 <- -0.5 * N * (log(2*pi) + log(sigma2))
  temp_expr2 <- -0.5 * sum((-ineff*eps + Zdelta)^2)/sigma2
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

u_panel_tnorm_bc95 <- function(object, estimator) {
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