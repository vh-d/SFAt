# MAIN FUNCTION -----------------------------------------------------------

#' sfa.fit
#'
#' Fits stochastic frontier analysis (SFA) model
#'
#' @param y dependent (production/cost) variable.
#' @param X variables of the production/cost function.
#' @param CM data for exogneous determinants of the conditional mean inefficiency.
#' @param CV data for exogneous determinants of the conditional variance of inefficiency term.
#' @param K matrix of panel data indeces.
#' @param ineff -1 (or 1) for production (or cost) function.
#' @param intercept TRUE/FALSE if the intercept term should be added to the main formula.
#' @param intercept_CM TRUE/FALSE if the intercept should be added to the conditional mean inefficiency formula.
#' @param intercept_CV TRUE/FALSE if the intercept should be added to the conditional inefficiency variance formula. (not functional yet)
#' @param structure "cs" for cross-section or "panel" for panel data model. (not functional yet)
#' @param dist distribution of inefficiency term ("hnorm", "exp", "tnorm").
#' @param spec specifies what model of endogeneous inefficiency term should be used (currently only bc95 for cross-section implemented).
#' @param sv_f starting values for frontier model parameters.
#' @param sv_cm starting values for conditional mean model parameters. (not functional yet)
#' @param sv_cv starting values for conditional variance model parameters. (not functional yet)
#' @param ll allows custom log-likelihood function that will be MINIMIZED.
#' @param opt_method optimization method.
#' @param opt_control list of options for optimization routine.
#' @param deb debug mode (TRUE/FALSE).
#' @details Experiment with different optimization algorithms:
#' @details SANN - slow but better results, maxit = 1e4, tmax = 15, temp = 1
#' @details L-BFGS-B - very fast, but can crash on infinite values
#' @return list object of class SFA.
#' @export
sfa.fit <- function(y,
                    X,
                    CM = NULL,
                    CV = NULL,
                    K = NULL,
                    ineff = -1,
                    structure = "cs",
                    dist = "tnorm",
                    spec = NULL,
                    intercept = TRUE,
                    intercept_CM = TRUE,
                    intercept_CV = TRUE,
                    sv_f = NULL,
                    sv_cm = NULL,
                    sv_cv = NULL,
                    ll = NULL,
                    opt_method = "SANN",
                    opt_control = NULL,
                    deb = F, # TRUE for debug reports
                    debll = F
                    ) {

  # ---- INIT ----

  # frontier model
  fcoeff_num <- ncol(X) # number of coefficients
  fcoeff_names <- colnames(X)

  if (intercept) {
    fcoeff_num <- fcoeff_num + 1
    X <- cbind(1, X)
    fcoeff_names <- c("intercept", fcoeff_names)
  }

  # conditional inefficiency mean
  if (is.null(CM)) {
    cmcoeff_num <- 0
  } else {
    cmcoeff_num <- ncol(CM) # number of coefficients
    cmcoeff_names <- colnames(CM) # coefficients names
  }

  # conditional inefficiency variance
  if (is.null(CV)) {
    cvarcoeff_num <- 0
  } else {
    cvarcoeff_num <- ncol(CV) # number of coefficients
    cvarcoeff_names <- colnames(CV) # coefficients names
  }

  # panel data dimensions
  ispanel <- !is.null(K)
  if (ispanel) {
    stopifnot(is.matrix(K) | is.data.frame(K), dim(K) == 2, dim(K)[2] == length(y), is.integer(K[, 2]))
    K1 <- levels(as.factor(K[, 1]))
    k = nlevels(K1)
    K2 <- K[, 2]
  }

  if (deb) {
    cat(ifelse(intercept == T, "X", "no X"),  "intercept, ",
        fcoeff_num, " X coefficient parameters, ",
        cmcoeff_num, " CM coefficient parameters,", "\n")
  }

  # ---- STARTING VALUES ----

  if (is.null(sv_f)) {
    if (is.null(spec)) {

      # fit OLS for starting values of betas and sigma
      lmfit <- lm(y ~ X - 1)
      if (deb) print(summary(lmfit))

      sv_f <- switch (dist,
                      tnorm = c(lmfit$coefficients,
                                mu = 0,
                                lnsigma2_u = 0.5,
                                lnsigma2_v = 0.5),

                      hnorm = c(lmfit$coefficients,
                                lnsigma2_u = 0.5,
                                lnsigma2_v = 0.5),

                      exp = c(lmfit$coefficients,
                              lnsigma2_u = 0.5,
                              lnsigma2_v = 0.5))

      names(sv_f)[1:length(fcoeff_names)] <- fcoeff_names
    } else {

      # fit OLS for starting values of beta parameters
      lmfit <- lm(y ~ X + I(ineff*CM) - 1)
      if (deb) print(summary(lmfit))

      sv_f <-
        c(lmfit$coefficients[1:fcoeff_num],
          if (intercept_CM) 0 else NULL,
          lmfit$coefficients[(fcoeff_num + 1) : (fcoeff_num + cmcoeff_num)],
          2, 2)

      if (intercept_CM) {
        # cmcoeff_num <- cmcoeff_num + 1
        CM <- cbind(mu_0 = 1, CM)
      }

      names(sv_f) <- c(fcoeff_names,
                       if (is.null(colnames(CM))) paste0("delta_", 1:(cmcoeff_num + intercept_CM)) else colnames(CM),
                       "lnsigma2_u",
                       "lnsigma2_v")
    }
  }

  if (deb) print(sv_f)

  if (is.null(ll)) {
    ll_fn_call <- paste0("ll", "_",
                         structure, "_",
                         dist,
                         if (is.null(spec)) NULL else paste0("_", spec))
  } else { # to-do: check existence of ll function given by user
    ll_fn_call <- ll
  }

  if (deb) {
    print(head(X))
    print(head(CM))
    print(ll_fn_call)
  }

  # lowerb <- c(rep(-Inf, fcoeff_num + cmcoeff_num), 0.000001, 0.000001)

  if (is.null(CM)) {
    est <- optim(sv_f,
                 fn = eval(parse(text = ll_fn_call)),
                 method = opt_method,
                 control = opt_control,
                 hessian = T,
                 y = y,
                 X = X,
                 ineff = ineff,
                 deb = debll)
  } else {
    est <- optim(sv_f,
                 fn = eval(parse(text = ll_fn_call)),
                 method = opt_method,
                 control = opt_control,
                 hessian = T,
                 y = y,
                 X = X,
                 CM = CM,
                 ineff = ineff,
                 debll = deb)
  }

  result <- list(coeff = est$par[1:fcoeff_num],
                 cm_coeff = est$par[(fcoeff_num + 1) : (fcoeff_num + cmcoeff_num + intercept_CM)],
                 residuals = as.vector(y - X %*% est$par[1:fcoeff_num]),
                 parameters = est$par,
                 N = length(y),
                 ineff = ineff,
                 ineff_name = if (ineff == -1) "production" else "cost",
                 data = list(y = y, X = X, CM = CM),
                 call = list(ineff = ineff,
                             intercept = intercept,
                             intercept_CM = intercept_CM,
                             dist = dist,
                             spec = spec,
                             structure = structure),
                 loglik = -est$val,
                 hessian = -est$hessian,
                 lmfit = lmfit)

  class(result) <- c("SFA")

  return(result)
}


# FORMULA FUNCTION --------------------------------------------------------

# stochastic frontier analysis
# this formula interface is not ready yet
# do not export
SFA <- function(formula,
                data = NULL,
                intercept = TRUE,
                intercept_CM = TRUE,
                dist = "hnormal",
                sv_f = NULL,
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

  # exogenous variables
  if (formula_length[2] > 1) {
    CM <- as.matrix(
      model.frame(
        formula(formula_ext,
                lhs = 0, rhs = 2)))
  } else {
    CM = NULL
  }

  sfa.fit(y = y,
          X = X,
          CM = CM,
          intercept = intercept,
          intercept_CM = intercept_CM,
          dist = dist,
          model = model,
          sv_f = sv_f,
          form = form,
          method = method,
          ...)
}


# SUMMARY FUNCTION --------------------------------------------------------
#' test statistics for SFA model
#' @param object object of class SFA
#' @export
summary.SFA <- function(object) {

  coef_sd <- sqrt(-diag(solve(object$hessian)))
  coef_tstats <- object$parameters/coef_sd
  coef_pvalues <- 2 * pt(q = abs(coef_tstats),
                         df = object$N - length(object$parameters),
                         lower.tail = FALSE)

  coef_conf_low <- object$parameters - qnorm(0.975)*coef_sd
  coef_conf_high <- object$parameters + qnorm(0.975)*coef_sd

  coef_table <- cbind(object$parameters,
                      coef_sd,
                      coef_tstats,
                      coef_pvalues,
                      coef_conf_low,
                      coef_conf_high)

  colnames(coef_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "95% (low)", "95% (high)")
  row.names(coef_table) <- names(object$parameters)

  fcoef_table <- coef_table[1:length(object$coeff), , drop = F]
  mccoef_table <- coef_table[length(object$coeff) + 1 : length(object$cm_coeff), , drop = F]
  sigmas_table <- coef_table[-(1:(length(object$coeff)+length(object$cm_coeff))), , drop = F]
  sigmas_t_table <- sigmas_table
  sigmas_t_table[] <- NA
  sigmas_t_table[, 1] <- sqrt(exp(sigmas_table[, 1]))

  ans <- list(call = object$call,
              N = object$N,
              loglik = object$loglik,
              parameters = object$parameters,
              coefficients_frontier = fcoef_table,
              coefficients_mc = mccoef_table,
              sigmas = sigmas_table,
              sigmas_t = sigmas_t_table)

  class(ans) <- "summary.SFA"

  return(ans)
}

#' Print summary for SFA objects.
#' @param object object of class SFA
#' @details Prints table of summary statistics for SFA model fitted by \code{SFA.fit} funtion.
#' @export
print.summary.SFA <- function(object) {

  col_names <- colnames(object$coefficients_frontier)

  separator_head <- t(gsub(".", "=", col_names))
  separator_mid <- t(gsub(".", "-", col_names))

  names(separator_head) = paste0(rep("=", max(nchar(col_names))), collapse = "")
  names(separator_mid) = paste0(rep("-", max(nchar(col_names))), collapse = "")

  outtable <- rbind(separator_head,
                    round(object$coefficients_frontier, 3),
                    separator_mid,
                    round(object$coefficients_mc, 3),
                    separator_mid,
                    round(object$sigmas, 3),
                    separator_mid,
                    round(object$sigmas_t, 3),
                    separator_mid)

  cat("Stochastic frontier model",
      "=========================",
      sep = "\n")

  cat(if (object$call$structure == "panel") "Panel" else "Cross-section", "data.\n")
  cat("Total observations:", object$N, if (object$call$structure == "panel") "in cross-section." else NULL,
      "\n")

  cat("Log-likelihood:", object$loglik, "\n")

  cat("\n", sep = "")
  print(outtable, quote = F)

  cat("\n", sep = "")
}
