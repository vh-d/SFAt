# NORMAL / H-NORMAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_hnorm <- c(lnsigma2_u = 0, lnsigma2_v = 0)

# likelihood function normal/half-normal distributional assumption
# params: beta, log(sigma_u^2), log(sigma_v^2)

ll_cs_hnorm <- function(params, y, X, ineff, deb) {

  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]


  lnsigma2_u <- params[nbetas + 1]
  lnsigma2_v <- params[nbetas + 2]

  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(p_beta), "\n",

               "Sigma_u^2: ", sigma2_u,
               "Sigma_v^2: ", sigma2_v, "\n")


  epsilon <- -ineff*as.vector(y - X %*% p_beta)
  sigma <- sqrt(sigma2_u + sigma2_v)

  N <- length(y)

  ll <-
    -(N*log(pi/2)) + # the const term
    -N*log(sigma) +
    +sum(log(pnorm(-(epsilon * sqrt(sigma2_u) / sqrt(sigma2_v)) / sigma))) - 0.5 / (sigma2_u + sigma2_v) * sum(epsilon^2)

  if (deb) {
    if (is.nan(ll)) print(paste(sigma))
    cat("Log-ll: ", ll, "\n")
  }

  return(-ll)
}

# Jondrow et al. (1982) as in Parmeter-Kumbhakar
u_cs_hnorm <- function(object, estimator) {
  # extract sigmas from the model object

  lnsigma2_u <- object$parameters[length(object$coeff) + 1]
  lnsigma2_v <- object$parameters[length(object$coeff) + 2]

  # derive the rest of the parameters
  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_u + sigma2_v

  mu_ast <- object$ineff * object$residuals * sigma2_u / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))

  if (is.null(u)) stop(paste0("Unknown type ", estimator," of conditional mean estimator."))

  return(u)
}