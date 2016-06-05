# MAIN FUNCTION -----------------------------------------------------------

#' sfa.fit
#'
#' Fits stochastic frontier analysis (SFA) model
#'
#' @param y dependent (production/cost) variable.
#' @param X variables of the production/cost function.
#' @param CM data for conditional mean model of the inefficiency (asymmetric error) term.
#' @param CV_u data for conditional variance model of the inefficiency (asymmetric error) term.
#' @param CV_v data for conditional variance model of the symmetric error term.
#' @param K matrix of panel data indeces.
#' @param ineff -1 (or 1) for production (or cost) function, where inefficiency decreases (or increases) the total output (or costs).
#' @param intercept TRUE/FALSE if the intercept term should be added to the main formula.
#' @param intercept_cm TRUE/FALSE if the intercept should be added to the conditional mean equation for the asymmetric term
#' @param intercept_cv_u TRUE/FALSE if the intercept should be added to the conditional inefficiency variance formula. (not functional yet)
#' @param intercept_cv_v TRUE/FALSE if the intercept should be added to the conditional inefficiency variance formula. (not functional yet)
#' @param structure "cs" for cross-section or "panel" for panel data model. (not functional yet)
#' @param dist distribution of inefficiency term ("hnorm", "exp", "tnorm").
#' @param spec specifies what model of endogeneous inefficiency term should be used (currently only bc95 for cross-section implemented).
#' @param sv_f starting values for frontier model parameters.
#' @param sv_cm starting values for conditional mean model parameters.
#' @param sv_cv_u starting values for conditional variance of the inefficiency term model parameters. (experimental)
#' @param sv_cv_v starting values for conditional variance of the symmetric term model parameters. (not functional yet)
#' @param ll allows custom log-likelihood function that will be MINIMIZED.
#' @param opt_method optimization method.
#' @param opt_control list of options for optimization routine.
#' @param deb debug mode (TRUE/FALSE).
#' @details
#' Notice that the choice of optimization method may have significant impact on the results. It is recommanded to experiment with different optimization algorithms. Recommended are:
#' \itemize{
#' \item SANN -- In general, this is the most robust method. It can be slow with larger datasets or more complex models but the results tend to be better if parameters \code{maxit, tmax, temp} are set correctly (maxit > 1e4+, tmax = 15, temp = 1).
#' \item L-BFGS-B -- Fastest, but can crashes on complex models (infinite log-likelihood values etc...). If starting values are set well, it leads to the same results as SANN but much faster.
#' \item BFGS -- Fast, but can crashes on complex models (infinite log-likelihood values etc...).
#' }
#' See help for \code{optim()} function.
#'
#' @return Returns object of the class SFA which is a list object consisting:
#' \itemize{
#' \item coeff -- coefficients for stochastic frontier model
#' \item coeff_cm -- coefficients for conditional mean of the inefficiency term model
#' \item coeff_cv_u -- coefficients for conditional variance of the inefficiency term model (heteroskedasticity in the inefficiency)
#' \item coeff_cv_v -- coefficients for conditional variance of the symmetric error term model (heteroskedasticity in the frontier model error)
#' \item residuals -- total residuals (= both u + v terms)
#' \item parameters -- vector of all parameters returend from miximization of log-likelihood
#' \item N -- total number of observations
#' \item ineff -- -1 (1 resp.) for production (cost resp.) function
#' \item ineff_name -- either "production" or "cost" string
#' \item data -- list of all data used for estimation (including unit vectors as intercepts if appropriate)
#' \item call is list of \itemize{
#'    \item intercept
#'    \item intercept_cm
#'    \item intercept_cv_u
#'    \item intercept_cv_v
#'    \item dist
#'    \item spec
#'    \item structure
#' }
#' \item loglik -- log-likehood
#' \item hessian -- hessian matrix
#' \item lmfit -- fitted linear model
#' }
#' @examples
#' See vignettes.
#'
#' @export
sfa.fit <- function(y,
                    X,
                    CM = NULL,
                    CV_u = NULL,
                    CV_v = NULL,
                    K = NULL,
                    ineff = -1L,
                    structure = "cs",
                    dist = c("tnorm", "hnorm", "exp"),
                    spec = NULL,
                    intercept = TRUE,
                    intercept_cm = TRUE,
                    intercept_cv_u = TRUE,
                    intercept_cv_v = TRUE,
                    sv_f = NULL, # to-do list of starting value vectors
                    sv_cm = NULL,
                    sv_cv_u = NULL,
                    sv_cv_v = NULL,
                    ll = NULL,
                    opt_method = "SANN",
                    opt_control = NULL,
                    deb = F, # TRUE for debug reports
                    debll = F
) {

  # ------ VALIDATE ARGUMENTS --------

  dist <- match.arg(dist)

  if (!(
    is.integer(ineff)
    & (length(ineff) == 1)
    & (ineff %in% c(-1L, 1L)))) {
    stop("ineff must be either -1 (for production function) or 1 (for cost function). Overriding ineff = -1.")
  }


  # ---- INIT ----

  # frontier model
  fcoeff_num <- ncol(X) # number of coefficients
  fcoeff_names <- colnames(X)

  if (intercept) {
    fcoeff_num <- fcoeff_num + 1
    X <- cbind(1, X)
    fcoeff_names <- c("intercept", fcoeff_names)
  }

  # inefficiency - (un)conditional mean for t-norm distribution
  if (is.null(CM)) {
    if (!intercept_cm) {
      warning("Either intercept_cm has to be set TRUE or CM has to be provided. Overriding intercept_cm = TRUE ...")
      intercept_cm = T
    }
    CM <- c(mu = 1.0)
    cm_model <- F
    coeff_cm_num <- 1
  } else {
    if (dist != "tnorm") stop("Conditional mean of inefficiency term model only possible for normal/t-normal model. ")
    if (intercept_cm) {
      CM <- cbind(mu_0 = 1, CM)
    }
    cm_model <- T
    coeff_cm_num <- ncol(CM) # number of coefficients
  }
  coeff_cm_names <- colnames(CM) # coefficients names

  # inefficiency - (un)conditional variance (heteroskedasticity in inefficiency)
  if (is.null(CV_u)) {
    if (!intercept_cv_u) {
      warning("Either intercept_cv has to be set TRUE or CV_u has to be provided. Overriding intercept_cv = TRUE ...")
      intercept_cv_u = T
    }
    CV_u <- c("lnsigma2_u" = 1)
    cv_u_model <- F
    coeff_cv_u_num <- 1
    coeff_cv_u_names <- "lnsigma2_u"
  } else {
    if (intercept_cv_u) {
      CV_u <- cbind(sigma_u_0 = 1, CV_u)
    }
    cv_u_model <- T
    coeff_cv_u_num <- ncol(CV_u) # number of coefficients
    coeff_cv_u_names <- colnames(CV_u) # coefficients names
  }

  #  symmetric error term - (un)conditional variance (heteroskedasticity in frontier model)
  if (is.null(CV_v)) {
    if (!intercept_cv_v) {
      warning("Either intercept_cv has to be set TRUE or CV_v has to be provided. Overriding intercept_cv = TRUE ...")
      intercept_cv_v = T
    }
    CV_v <- c("lnsigma2_v" = 1)
    cv_v_model <- F
    coeff_cv_v_num <- 1
    coeff_cv_v_names <- "lnsigma2_v"
  } else {
    if (intercept_cv_v) {
      CV_v <- cbind(sigma_v_0 = 1, CV_v)
    }
    cv_v_model <- T
    coeff_cv_v_num <- ncol(CV_v) # number of coefficients
    coeff_cv_v_names <- colnames(CV_v) # coefficients names
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
    cat(ifelse(intercept == T,
               "X",
               "no X"), " intercept, ", "\n",
        fcoeff_num,     " X coefficient parameters, ", "\n",
        coeff_cm_num,   " CM coefficient parameters,", "\n",
        coeff_cv_u_num, " CV_u coefficient parameters,", "\n",
        coeff_cv_v_num, " CV_v coefficient parameters,", "\n")
  }


  # ---- MODEL SPECIFICATION ----

  model_spec <- paste0(structure, "_",
                       dist,
                       if (is.null(spec)) NULL else paste0("_", spec))

  model_parameters <- eval(parse(text = paste0("par_", model_spec)))


  # ---- STARTING VALUES ----

  # fit OLS for starting values of frontier model coefficients
  lmfit <- lm(y ~ X - 1) # intercept is part of X already
  if (deb) print(summary(lmfit))

  if (is.null(sv_f)) {
      sv_f <- lmfit$coefficients
      names(sv_f)[1:length(fcoeff_names)] <- fcoeff_names
  } # to-do: else check length

  if (dist == "tnorm") {
    if (is.null(sv_cm)) {
      sv_cm <- rep(0.0, coeff_cm_num)
      names(sv_cm) <- coeff_cm_names
    } # to-do: else check length
  } else sv_cm = NULL

  if (is.null(sv_cv_u)) {
    sv_cv_u <- rep(0.0, coeff_cv_u_num)
    names(sv_cv_u) <- coeff_cv_u_names
  } # to-do: else check length

  if (is.null(sv_cv_v)) {
    sv_cv_v <- rep(0.0, coeff_cv_v_num)
    names(sv_cv_v) <- coeff_cv_v_names
  } # to-do: else check length

  # concatenate all coefficients and model parameters into a single vector for the optimization routine
  sv <- c(sv_f,
          sv_cv_u,
          sv_cv_v,
          sv_cm,
          model_parameters)

  if (deb) cat("Starting values: ",
               sv)

  if (is.null(ll)) {
    ll_fn_call <- paste0("ll", "_", model_spec)
  } else { # to-do: check existence of ll function given by user
    ll_fn_call <- ll
  }

  if (deb) {
    print(head(X))
    print(head(CM))
    print(head(CV_u))
    print(head(CV_v))
    print(ll_fn_call)
  }

  indeces <- cumsum(c(fcoeff_num,
                      coeff_cv_u_num,
                      coeff_cv_v_num,
                      coeff_cm_num))

  # validate parameter vector
  parcheck <- do.call(what = paste0("par_", model_spec, "_check"),
                      args = list(
                        params = sv,
                        indeces = indeces,
                        y = y,
                        X = X,
                        CM = CM,
                        CV_u = CV_u,
                        CV_v = CV_v))
  if (deb & parcheck) cat("Parameters check OK.", "\n")


  # ------- MLE ----------

  est <- optim(sv,
               fn = eval(parse(text = ll_fn_call)),
               method = opt_method,
               control = opt_control,
               hessian = T,
               indeces = indeces,
               y = y,
               X = X,
               CM = CM,
               CV_u = CV_u,
               CV_v = CV_v,
               ineff = ineff,
               deb = debll)


  # ---------- RETURN ------------

  coeff_frontier <- est$par[1 : (indeces[1])]
  coeff_cv_u     <- est$par[(1 + indeces[1]) : (indeces[2])]
  coeff_cv_v     <- est$par[(1 + indeces[2]) : (indeces[3])]
  coeff_cm <-  if (dist == "tnorm")
                    est$par[(1 + indeces[3]) : (indeces[4])]

  result <- list(coeff_frontier = coeff_frontier,
                 cm_model       = cm_model,
                 coeff_cm       = coeff_cm,
                 cv_u_model     = cv_u_model,
                 coeff_cv_u     = coeff_cv_u,
                 cv_v_model     = cv_v_model,
                 coeff_cv_v     = coeff_cv_v,
                 indeces        = indeces,
                 residuals      = as.vector(y - X %*% est$par[1:fcoeff_num]),
                 parameters     = est$par,
                 N              = length(y),
                 ineff          = ineff,
                 ineff_name     = if (ineff == -1) "production" else "cost",
                 data           = list(y = y,
                                       X = X,
                                       CM = CM,
                                       CV_u = CV_u,
                                       CV_v = CV_v),
                 call           = list(ineff = ineff,
                                       intercept = intercept,
                                       intercept_cm = intercept_cm,
                                       dist = dist,
                                       spec = spec,
                                       model_spec = model_spec,
                                       structure = structure),
                 loglik         = -est$val,
                 hessian        = -est$hessian,
                 lmfit          = lmfit)

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
                intercept_cm = TRUE,
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
          intercept_cm = intercept_cm,
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

  # indeces <- cumsum(c(length(object$coeff_frontier),
  #                     length(object$coeff_cv_u),
  #                     length(object$coeff_cv_v),
  #                     length(object$coeff_cm)))

  # coeff_f_table <- coef_table[1:length(object$coeff), , drop = F]
  # coeff_cm_table <- if (!is.null(object$coeff_cm)) coef_table[length(object$coeff) + 1 : length(object$coeff_cm), , drop = F] else NULL
  # sigmas_table <- coef_table[-(1:(length(object$coeff)+length(object$coeff_cm))), , drop = F]
  # sigmas_t_table <- as.matrix(do.call(paste0("t_par_", object$call$model_spec),
  #                                     args = list(pars = as.vector(sigmas_table[, 1]))))
  # colnames(sigmas_t_table)[1] <- colnames(sigmas_table)[1]
  indeces <- object$indeces
  coeff_f_table    <- coef_table[              1  : indeces[1], , drop = F]
  coeff_cv_u_table <- coef_table[(indeces[1] + 1) : indeces[2], , drop = F]
  coeff_cv_v_table <- coef_table[(indeces[2] + 1) : indeces[3], , drop = F]
  coeff_cm_table   <- if (!is.null(object$coeff_cm))
                      coef_table[(indeces[3] + 1) : indeces[4], , drop = F] else NULL
  # sigmas_t_table <- as.matrix(do.call(paste0("t_par_", object$call$model_spec),
  #                                     args = list(pars = as.vector(sigmas_table[, 1]))))
  # colnames(sigmas_t_table)[1] <- colnames(sigmas_table)[1]

  ans <- list(call = object$call,
              N = object$N,
              loglik = object$loglik,
              parameters = object$parameters,
              coefficients_frontier = coeff_f_table,
              coefficients_cm = coeff_cm_table,
              coefficients_cv_u = coeff_cv_u_table,
              coefficients_cv_v = coeff_cv_v_table)

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

  maxrowname <- max(c(nchar(rownames(object$coefficients_frontier)),
                      nchar(rownames(object$coefficients_mc)),
                      nchar(rownames(object$sigmas)),
                      nchar(rownames(object$sigmas_t))))

  rownames(separator_head) <- paste0(rep("=", maxrowname), collapse = "")
  rownames(separator_mid) <-  paste0(rep("-", maxrowname), collapse = "")

  outtable <- rbind(separator_head,
                    round(object$coefficients_frontier, 3),
                    separator_mid,
                    if (!is.null(object$coefficients_mc)) round(object$coefficients_mc, 3) else NULL,
                    if (!is.null(object$coefficients_mc)) separator_mid else NULL,
                    round(object$coefficients_cv_u, 3),
                    separator_mid,
                    round(object$coefficients_cv_v, 3),
                    separator_mid
                    # ,
                    # cbind(round(object$sigmas_t, 3), matrix("", nrow = nrow(object$sigmas_t), ncol = ncol(separator_mid)-1)),
                    # separator_mid
                    )

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
