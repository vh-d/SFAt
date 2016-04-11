# MAIN FUNCTION -----------------------------------------------------------

#' sfa.fit
#'
#' Fits stochastic frontier analysis (SFA) model
#'
#' @param y dependent (production/cost) variable.
#' @param X variables of the production/cost function.
#' @param Z exogneous determinants of the mean inefficiency location.
#' @param intercept TRUE/FALSE if the intercept should be included in the main formula.
#' @param intercept_Z TRUE/FALSE if the intercept should be included in the inefficiency location formula.
#' @param structure "cs" for cross-section or "panel" for panel data model.
#' @param dist distribution of inefficiency term ("hnorm", "exp", "tnorm").
#' @param spec specifies what model of endogeneous inefficiency term should be used (currently only bc95 for cross-section implemented).
#' @param start_val starting value of model parameters to be passed to optimization routine.
#' @param ineff -1 (or 1) for production (or cost) function.
#' @param opt_method optimization method.
#' @param deb debug mode (TRUE/FALSE).
#' @param control_opt list of options for optimization routine.
#'
#' @return list
#' @export
sfa.fit <- function(y, X,
                    Z = NULL,
                    intercept = TRUE,
                    intercept_Z = TRUE,
                    structure = "cs",
                    dist = "tnorm",
                    spec = "bc95",
                    start_val = NULL,
                    ineff = -1,
                    opt_method = "BFGS",
                    deb = F, # TRUE for debug reports
                    control_opt = NULL) {

  n_betas <- ncol(X)
  n_deltas <- if (is.null(Z)) 0 else ncol(Z)
  x_names <- colnames(X)

  if (intercept) {
    n_betas <- n_betas + 1
    X <- cbind(1, X)
    x_names <- c("intercept", x_names)
  }

  if (deb) {
    cat(ifelse(intercept == T, "X", "no X"),  "intercept, ",
        n_betas, " X coefficient parameters, ",
        n_deltas, " Z coefficient parameters,", "\n")
  }


  # STARTING VALUES

  if (is.null(start_val)) {
    if (is.null(spec)) {

      # fit OLS for starting values of betas and sigma
      lmfit <- lm(y ~ X - 1)
      if (deb) print(summary(lmfit))

      start_val <- switch (dist,
                           tnorm = c(lmfit$coefficients,
                                     mu = 0,
                                     sigma = 1,
                                     lambda = 1),
                           hnorm = c(lmfit$coefficients,
                                     var_u = 1,
                                     var_v = 1),
                           exp = c(lmfit$coefficients,
                                   var_u = 1,
                                   var_v = 1))

      names(start_val)[1:length(x_names)] <- x_names
    } else {

      # fit OLS for starting values of beta parameters
      lmfit <- lm(y ~ X + I(ineff*Z) - 1)
      if (deb) print(summary(lmfit))

      start_val <-
        c(lmfit$coefficients[1:n_betas],
          if (intercept_Z) 0 else NULL,
          lmfit$coefficients[(n_betas + 1) : (n_betas + n_deltas)],
          2, 2)

      if (intercept_Z) {
        # n_deltas <- n_deltas + 1
        Z <- cbind(mu_0 = 1, Z)
      }

      names(start_val) <- c(x_names,
                            if (is.null(colnames(Z))) paste0("delta_", 1:(n_deltas + intercept_Z)) else colnames(Z),
                                "sigma2_u", "sigma2_v")
    }
  }

  if (deb) print(start_val)

  ll_fn_call <- parse(text = paste0("ll", "_",
                                    structure, "_",
                                    dist,
                                    if (is.null(spec)) NULL else paste0("_", spec)))

  # lowerb <- c(rep(-Inf, n_betas + n_deltas), 0.000001, 0.000001)

  if (is.null(Z)) {
   est <- optim(start_val,
               fn = eval(ll_fn_call),
               method = opt_method,
               control = control_opt,
               hessian = T,
               y = y,
               X = X,
               ineff = ineff,
               deb = deb)
  } else {
     est <- optim(start_val,
               fn = eval(ll_fn_call),
               method = opt_method,
               control = control_opt,
               hessian = T,
               y = y,
               X = X,
               Z = Z,
               ineff = ineff,
               deb = deb)
  }

  result <- list(coefficients = est$par[1:n_betas],
              coefficients_Z = est$par[(n_betas + 1) : (n_betas + n_deltas + intercept_Z)],
              residuals = as.vector(y - X %*% est$par[1:n_betas]),
              parameters = est$par,
              N = length(y),
              ineff = ineff,
              ineff_name = if (ineff == -1) "production" else "cost",
              data = list(y = y, X = X, Z = Z),
              call = list(ineff = ineff,
                          intercept = intercept,
                          intercept_Z = intercept_Z,
                          dist = dist,
                          spec = spec,
                          structure = structure),
              loglik = -est$val,
              hessian = est$hessian,
              lmfit = lmfit)

  class(result) <- c("SFA")

  return(result)
}


# FORMULA FUNCTION --------------------------------------------------------

#' stochastic frontier analysis
# this formula interface is not ready yet
#' @export
SFA <- function(formula,
                data = NULL,
                intercept = TRUE,
                intercept_Z = TRUE,
                dist = "hnormal",
                start_val = NULL,
                par_mu = NULL,
                ineff = -1,
                opt_method = "BFGS", ...){

  formula_ext <- Formula(formula)
  formula_length <- length(formula_ext)


  y <- as.vector(
    model.frame(
      formula(formula_ext,
              lhs = 1, rhs = 0)))

  X <- as.matrix(
    model.frame(
      formula(formula_ext,
              lhs = 0, rhs = 1)))

  # exdogenous variables
  if (formula_length[2] > 1) {
    Z <- as.matrix(
      model.frame(
        formula(formula_ext,
                lhs = 0, rhs = 2)))
  } else {
    Z = NULL
  }

  sfa.fit(y = y,
          X = X,
          Z = Z,
          intercept = intercept,
          intercept_Z = intercept_Z,
          dist = dist,
          model = model,
          start_val = start_val,
          form = form,
          method = method,
          ...)
}


# SUMMARY FUNCTION --------------------------------------------------------
#' test statistics for SFAplus model
#' @export
summary.SFA <- function(object) {

  coef_sd <- sqrt(diag(solve(object$hessian)))
  coef_tstats <- object$parameters/coef_sd
  coef_pvalues <- 2 * pt(q = abs(coef_tstats),
                         df = object$N - length(object$parameters),
                         lower.tail = FALSE)
  coef_table <- round(x = cbind(object$parameters,
                                coef_sd,
                                coef_tstats,
                                coef_pvalues),
                      digits = 3)
  colnames(coef_table) <- c("Estimate", "Std. Error","t-stat", "Pr(>|t|)")
  row.names(coef_table) <- names(object$parameters)

  print(coef_table[1:length(object$coefficients),])

  # LR test
  lrtest(object)
}

# LR test function --------------------------------------------------------
#' LR test for SFA class
#' @export
lrtest.SFA <- function(object) {

  LR_test_stat <- 2*(object$loglik - logLik(object$lmfit))
  LR_chisq_df <- length(object$parameters) - attributes(logLik(object$lmfit))$df
  if (LR_chisq_df > 1) {
    LR_pvalue <-
      0.25*pchisq(LR_test_stat, LR_chisq_df-2, lower.tail = FALSE) +
      0.5*pchisq(LR_test_stat, LR_chisq_df-1, lower.tail = FALSE) +
      0.25*pchisq(LR_test_stat, LR_chisq_df, lower.tail = FALSE)
  } else {
    LR_pvalue <-
      0.5*pchisq(LR_test_stat, LR_chisq_df-1, lower.tail = FALSE) +
      0.5*pchisq(LR_test_stat, LR_chisq_df, lower.tail = FALSE)
  }

  cat("\n")
  cat(paste0("LR Chisq: ", round(LR_test_stat, 3)), "\n")
  cat(paste0("Chisq Df: ", round(LR_chisq_df, 3)), "\n")
  cat(paste0("Pr(>Chisq): ", round(LR_pvalue, 3)))
}
