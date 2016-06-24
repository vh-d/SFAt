# MAIN FUNCTION -----------------------------------------------------------

#' Estimate stochastic frontier analysis (SFA) models.
#'
#' @param y dependent (production/cost) variable.
#' @param X variables of the production/cost function.
#' @param CM data for conditional mean model of the inefficiency (asymmetric error) term.
#' @param CV_u data for conditional variance model of the inefficiency (asymmetric error) term.
#' @param CV_v data for conditional variance model of the symmetric error term.
#' @param ineff -1 (or 1) for production (or cost) function, where inefficiency decreases (or increases) the total output (or costs).
#' @param intercept list of logical values:
#' \describe{
#' \item{f}{TRUE if the intercept term should be added to the main formula.}
#' \item{cm}{TRUE if the intercept should be added to the conditional mean equation for the asymmetric term}
#' \item{cv_u}{TRUE if the intercept should be added to the conditional inefficiency variance formula.}
#' \item{cv_v}{TRUE if the intercept should be added to the conditional inefficiency variance formula.}
#' }
#' @param dist distribution of inefficiency term ("hnorm", "exp", "tnorm").
#' @param spec specifies what model of endogeneous inefficiency term should be used (currently only bc95 for cross-section implemented).
#' @param sv list. starting values for:
#' \describe{
#' \item{f}{frontier model coefficients}
#' \item{cm}{starting values for conditional mean model parameters.}
#' \item{cv_u}{starting values for conditional variance of the inefficiency term model parameters.}
#' \item{cv_v}{starting values for conditional variance of the symmetric term model parameters.}
#' }
#' @param ll allows custom log-likelihood function that will be MINIMIZED.
#' @param opt_method optimization method.
#' @param opt_control list of options for optimization routine.
#' @param deb debug mode (TRUE/FALSE).
#' @param debll debug mode of log likelihood functions (TRUE/FALSE).
#' @details
#' \code{sfa.fit()} is the main workhorse function that actually estimate the SFA model. The \code{SFA.formula()} and \code{SFA.list()} methods are provided for more convenient user interface.
#'
#' For cross-section data model, the following distributions are currently supported:
#' \itemize{
#' \item normal/half-normal model
#' \item normal/truncated-normal model
#' \item normal/exponential model
#' }
#' @section Estimation:
#' Within all these models heteroskedasticity in both symmetric and asymmetric error terms can be explicitly modeled. It can be done by providing matrices of explanatory variables (\code{CV_v} for the symmetric error and \code{CV_u} for the inefficiency term).
#' Conditional mean of the inefficiency term can be modeled only within the normal/t-normal model.
#' Models are estimated via maximum likelihood estimators following established literature on the topic.
#'
#' @section Starting values:
#' Starting values are by default coefficients of a linear (OLS) model estimated during within the \code{sfa.fit()} function. Or they can be supplied by user as a list of vectors.
#'
#' @section Optimization:
#' Optimization of log-likelihood functions is currently done by R's default \code{optim()} function. Notice that the choice of optimization method may have significant impact on the results and it is highly recommanded to experiment with different optimization algorithms. Recommended are:
#'
#' \describe{
#' \item{SANN}{In general, this is the most robust method. It can be slow with larger datasets or more complex models but the results tend to be better if parameters \code{maxit, tmax, temp} are set correctly (maxit > 1e4+, tmax = 15, temp = 1).}
#' \item{L-BFGS-B}{Fastest, but can crashes on complex models (infinite log-likelihood values etc...). If starting values are set well, it leads to the same results as SANN but much faster.}
#' \item{BFGS}{Fast, but can crashes on complex models (infinite log-likelihood values etc...).}
#' }
#' See help for \code{optim()} function.
#'
#' @return
#' Returns object of the class SFA which is a list object consisting:
#' \describe{
#' \item{coeff}{coefficients for stochastic frontier model}
#' \item{coeff_cm}{coefficients for conditional mean of the inefficiency term model}
#' \item{coeff_cv_u}{coefficients for conditional variance of the inefficiency term model (heteroskedasticity in the inefficiency)}
#' \item{coeff_cv_v}{coefficients for conditional variance of the symmetric error term model (heteroskedasticity in the frontier model error)}
#' \item{residuals}{total residuals (= both u + v terms)}
#' \item{parameters}{vector of all parameters returend from miximization of log-likelihood}
#' \item{N}{total number of observations}
#' \item{ineff}{-1 (1 resp.) for production (cost resp.) function}
#' \item{ineff_name}{either "production" or "cost" string}
#' \item{data}{list of all data used for estimation (including unit vectors as intercepts if appropriate)}
#' \item{call}{is list of \itemize{
#'    \item intercept
#'    \item dist
#'    \item spec
#'    \item structure
#'    \item sv
#' }}
#' \item{loglik}{Total log-likehood.}
#' \item{hessian}{A hessian matrix as returned by optim()}
#' \item{lmfit}{lm object result of fitted linear model.}
#' }
#' @examples
#' See vignettes.
#'
#' @rdname SFA
#' @export
sfa.fit <- function(y,
                    X,
                    CM = NULL,
                    CV_u = NULL,
                    CV_v = NULL,
                    # K = NULL,
                    ineff = -1L,
                    # structure = "cs",
                    dist = c("tnorm", "hnorm", "exp"),
                    spec = NULL,
                    intercept = list(f = TRUE,
                                     cm = TRUE,
                                     cv_u = TRUE,
                                     cv_v = TRUE),
                    sv = list(f = NULL,
                              cm = NULL,
                              cv_u = NULL,
                              cv_v = NULL),
                    ll = NULL,
                    opt_method = "SANN",
                    opt_control = NULL,
                    deb = F, # TRUE for debug reports
                    debll = F
) {

  structure <- "cs"
  # ------ VALIDATE ARGUMENTS --------

  dist <- match.arg(dist)

  if (!(
    is.integer(ineff)
    & (length(ineff) == 1)
    & (ineff %in% c(-1L, 1L)))) {
    stop("ineff must be either -1 (for production function) or 1 (for cost function). Overriding ineff = -1.")
  }

  # ---- MISSING DATA ----

  if (deb) cat("y:", length(y), "X:", dim(X), "CM:", dim(CM), "CV_u", dim(CV_u), "CV_v:", dim(CV_v), "\n")
  y[!is.finite(y)] <- NA
  X[!is.finite(X)] <- NA
  if (!is.null(CV_u)) CV_u[!is.finite(CV_u)] <- NA
  if (!is.null(CV_v)) CV_v[!is.finite(CV_v)] <- NA
  if (!is.null(CM))   CM[!is.finite(CM)] <- NA

  cc <- complete.cases(y, X, if (!is.null(CM)) CM, if (!is.null(CV_u)) CV_u, if (!is.null(CV_v)) CV_v)
  ccn <- sum(cc)

  if (length(y) != ccn) {
    warning("Dropping ", length(y) - ccn, " observartions due to missingness.")
    oy    <- y
    oX    <- X
    oCV_u <- CV_u
    oCV_v <- CV_v
    oCM   <- CM

    y    <- subset(y, cc)
    X    <- subset(X, cc)
    CV_u <- subset(CV_u, cc)
    CV_v <- subset(CV_v, cc)
    CM   <- subset(CM, cc)

    missingness <- T
  } else missingness <- F

  # ---- INIT ----

  # frontier model

  if (intercept$f) {
    X <- cbind("(Intercept)" = 1, X)
    # coeff_f_names <- c("(Intercept)", coeff_f_names)
  }
  coeff_f_names <- colnames(X)
  coeff_f_n     <- ncol(X) # number of coefficients in the frontier model equation

  # inefficiency - (un)conditional mean for t-norm distribution
  if (is.null(CM)) {
    if (!intercept$cm) {
      warning("Either intercept$cm has to be set TRUE or CM has to be provided. Overriding intercept$cm = TRUE ...")
      intercept$cm = T
    }
    CM <- c(mu = 1.0)
    cm_model <- F
    coeff_cm_num <- 1
  } else {
    if (dist != "tnorm") warning("Conditional mean of inefficiency term model only possible for normal/t-normal model. ")
    if (intercept$cm) {
      CM <- cbind("(Intercept)" = 1.0, CM)
    }
    cm_model <- T
    coeff_cm_num <- ncol(CM) # number of coefficients
  }
  coeff_cm_names <- colnames(CM) # coefficients names

  # inefficiency - (un)conditional variance (heteroskedasticity in inefficiency)
  if (is.null(CV_u)) {
    if (!intercept$cv_u) {
      warning("Either intercept$cv has to be set TRUE or CV_u has to be provided. Overriding intercept$cv = TRUE ...")
      intercept$cv_u = T
    }
    CV_u <- c("(Intercept)" = 1.0)
    cv_u_model <- F
    coeff_cv_u_num <- 1
    coeff_cv_u_names <- "(Intercept)"
  } else {
    if (intercept$cv_u) {
      CV_u <- cbind("(Intercept)" = 1.0, CV_u)
    }
    cv_u_model <- T
    coeff_cv_u_num <- ncol(CV_u) # number of coefficients
    coeff_cv_u_names <- colnames(CV_u) # coefficients names
  }

  #  symmetric error term - (un)conditional variance (heteroskedasticity in frontier model)
  if (is.null(CV_v)) {
    if (!intercept$cv_v) {
      warning("Either intercept$cv has to be set TRUE or CV_v has to be provided. Overriding intercept$cv = TRUE ...")
      intercept$cv_v = T
    }
    CV_v <- c("(Intercept)" = 1.0)
    cv_v_model <- F
    coeff_cv_v_num <- 1
    coeff_cv_v_names <- "(Intercept)"
  } else {
    if (intercept$cv_v) {
      CV_v <- cbind(intercept = 1.0, CV_v)
    }
    cv_v_model <- T
    coeff_cv_v_num <- ncol(CV_v) # number of coefficients
    coeff_cv_v_names <- colnames(CV_v) # coefficients names
  }

  # # panel data dimensions
  # ispanel <- !is.null(K)
  # if (ispanel) {
  #   stopifnot(is.matrix(K) | is.data.frame(K), dim(K) == 2, dim(K)[2] == length(y), is.integer(K[, 2]))
  #   K1 <- levels(as.factor(K[, 1]))
  #   k = nlevels(K1)
  #   K2 <- K[, 2]
  # }

  if (deb) {
    cat(ifelse(intercept$f == T,
               "X",
               "no X"), " intercept, ", "\n",
        coeff_f_n,      " X coefficient parameters, ", "\n",
        coeff_cm_num,   " CM coefficient parameters,", "\n",
        coeff_cv_u_num, " CV_u coefficient parameters,", "\n",
        coeff_cv_v_num, " CV_v coefficient parameters,", "\n")
  }

  if (deb) cat(str(CV_v))
  # ---- MODEL SPECIFICATION ----

  model_spec <- paste0(structure, "_",
                       dist,
                       if (is.null(spec)) NULL else paste0("_", spec))

  model_parameters <- eval(parse(text = paste0("par_", model_spec)))


  # ---- STARTING VALUES ----

  # fit OLS for starting values of frontier model coefficients
  if (deb) cat("Fitting OLS for starting values...\n")
  if (deb) cat("y length:", length(y), "X dim:", dim(X), "\n")
  if (deb) cat(str(X), "\n", str(names(X)), "\n", sum(is.na(X)), "\n")
  lmfit <- lm(y ~ X - 1) # intercept is part of X already
  if (deb) print(summary(lmfit))

  if (is.null(sv$f)) {
    sv$f <- lmfit$coefficients
    names(sv$f)[1:length(coeff_f_names)] <- coeff_f_names
  } # to-do: else check length

  if (dist == "tnorm") {
    if (is.null(sv$cm)) {
      sv$cm <- rep(0.0, coeff_cm_num)
      names(sv$cm) <- coeff_cm_names
    } # to-do: else check length
  } else sv$cm = NULL

  if (is.null(sv$cv_u)) {
    sv$cv_u <- rep(5.0, coeff_cv_u_num)
    names(sv$cv_u) <- coeff_cv_u_names
  } # to-do: else check length

  if (is.null(sv$cv_v)) {
    sv$cv_v <- rep(5.0, coeff_cv_v_num)
    names(sv$cv_v) <- coeff_cv_v_names
  } # to-do: else check length

  # concatenate all coefficients and model parameters into a single vector for the optimization routine
  svv <- c(sv$f,
           sv$cv_u,
           sv$cv_v,
           sv$cm,
           model_parameters)

  if (deb) cat("Starting values: ", svv)

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

  indeces <- cumsum(c(coeff_f_n,
                      coeff_cv_u_num,
                      coeff_cv_v_num,
                      coeff_cm_num))

  # validate parameter vector
  parcheck <- do.call(what = paste0("par_", model_spec, "_check"),
                      args = list(
                        params = svv,
                        indeces = indeces,
                        y = y,
                        X = X,
                        CV_u = CV_u,
                        CV_v = CV_v,
                        CM = CM))
  if (deb & parcheck) cat("Parameters check OK.", "\n")


  # ------- MLE ----------
  require(maxLik)

  est <- maxLik(start = svv,
               logLik = eval(parse(text = ll_fn_call)),
               method = opt_method,
               control = opt_control,
               indeces = indeces,
               y = y,
               X = X,
               CM = CM,
               CV_u = CV_u,
               CV_v = CV_v,
               ineff = ineff,
               deb = debll)


  if (deb) {
    cat(est$par, "\n")
    cat(coeff_f_names, coeff_cv_u_names, coeff_cv_v_names, if (!is.null(CM)) coeff_cm_names else NULL, "\n")
  }

  names(est$estimate) <- c(coeff_f_names, coeff_cv_u_names, coeff_cv_v_names, if (dist == "tnorm") coeff_cm_names else NULL)


  # ---------- RETURN ------------

  if (missingness) {
    y <- oy
    X <- oX
    CV_u <- oCV_u
    CV_v <- oCV_v
    CM <- oCM
  }


  coeff_frontier <- est$estimate[1 : (indeces[1])]
  coeff_cv_u     <- est$estimate[(1 + indeces[1]) : (indeces[2])]
  coeff_cv_v     <- est$estimate[(1 + indeces[2]) : (indeces[3])]
  coeff_cm <-  if (dist == "tnorm")
    est$estimate[(1 + indeces[3]) : (indeces[4])]

  result <- list(coeff_frontier = coeff_frontier,
                 cm_model       = cm_model,
                 coeff_cm       = coeff_cm,
                 cv_u_model     = cv_u_model,
                 coeff_cv_u     = coeff_cv_u,
                 cv_v_model     = cv_v_model,
                 coeff_cv_v     = coeff_cv_v,
                 indeces        = indeces,
                 residuals      = as.vector(y - X %*% est$estimate[1:coeff_f_n]),
                 parameters     = est$estimate,
                 N              = ccn,
                 N_total        = length(y),
                 ineff          = ineff,
                 ineff_name     = if (ineff == -1) "production" else "cost",
                 data           = list(y = y,
                                       X = X,
                                       CM = CM,
                                       CV_u = CV_u,
                                       CV_v = CV_v),
                 call           = list(ineff = ineff,
                                       intercept = intercept,
                                       dist = dist,
                                       spec = spec,
                                       model_spec = model_spec,
                                       sv = sv,
                                       structure = structure),
                 loglik         = est$maximum,
                 hessian        = est$hessian,
                 lmfit          = lmfit,
                 opt            = est)

  attr(result$data, "missingness") <- missingness
  if (missingness) attr(result$data, "cc") <- cc

  class(result) <- c("SFA")

  return(result)
}


# FORMULA METHOD --------------------------------------------------------

#' @rdname SFA
#' @export
SFA.formula <- function(formula,
                        data = NULL,
                        cm = ~1,
                        cv_u = ~1,
                        cv_v = ~1,
                        form = c("production", "cost"),
                        deb = FALSE,
                        ...){

  form <- match.arg(form)
  ineff <- if (form == "production") -1L else 1L

  # extract the response variable
  y <- data[, all.vars(formula)[1]]

  # extract the explanatory variables
  temp <- options("na.action")
  options(na.action = "na.pass")
  fMat <- model.matrix(formula, data = data)
  if (!is.null(cm)   & class(cm)   == "formula") {cmMat  <- model.matrix(cm,   data = data)} else cmMat  <- NULL
  if (!is.null(cv_u) & class(cv_u) == "formula") {cvuMat <- model.matrix(cv_u, data = data)} else cvuMat <- NULL
  if (!is.null(cv_v) & class(cv_v) == "formula") {cvvMat <- model.matrix(cv_v, data = data)} else cvvMat <- NULL
  # options(na.action = temp)

  if (deb) {
    print(head(fMat)); print(head(cmMat)); print(head(cvuMat)); print(head(cvvMat))
  }

  sfa.fit(y = y,
          X = fMat,
          ineff = ineff,
          CM = cmMat,
          CV_u = cvuMat,
          CV_v = cvvMat,
          intercept = list(f = F, cm = F, cv_u = F, cv_v = F),
          deb = deb,
          # intercept = list(f = attr(terms(formula), "intercept"),
          #                  cm = attr(terms(cm), "intercept"),
          #                  cv_u = attr(terms(cv_u), "intercept"),
          #                  cv_v = attr(terms(cv_v), "intercept")),
          ...)
}

#' @rdname SFA
#' @export
SFA.list <- function(formulas,
                     data = NULL,
                     ...) {
  if (length(names(formulas)) != length(formulas))
    names(formulas) <- c("f", "cm", "cv_u", "cv_v")[1:length(formulas)]

  do.call(SFA.formula, args = list(formulas, data = data, ...))
}


#' @rdname SFA
#' @export
SFA <- function(object, ...) {
  UseMethod("SFA")
}

# SUMMARY METHOD --------------------------------------------------------
#' test statistics for SFA model
#' @param object object of class SFA
#' @export
summary.SFA <- function(object) {

  coef_sd <- sqrt(-diag(solve(object$hessian, tol = 100*.Machine$double.xmin)))
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

  indeces <- object$indeces
  coeff_f_table    <- coef_table[              1  : indeces[1], , drop = F]
  coeff_cv_u_table <- coef_table[(indeces[1] + 1) : indeces[2], , drop = F]
  coeff_cv_v_table <- coef_table[(indeces[2] + 1) : indeces[3], , drop = F]
  coeff_cm_table   <- if (!is.null(object$coeff_cm))
    coef_table[(indeces[3] + 1) : indeces[4], , drop = F] else NULL

  ans <- list(call = object$call,
              N = object$N,
              loglik = object$loglik,
              parameters = object$parameters,
              coefficients = coef_table,
              coefficients_frontier = coeff_f_table,
              coefficients_cm = coeff_cm_table,
              coefficients_cv_u = coeff_cv_u_table,
              coefficients_cv_v = coeff_cv_v_table)

  class(ans) <- "summary.SFA"

  return(ans)
}

#' Print summary for SFA objects.
#' @param object object of class SFA
#' @details Prints table of summary statistics for SFA model fitted by \code{sfa.fit} funtion.
#' @export
print.summary.SFA <- function(object) {

  head_cv_u <- matrix(NA, 1, ncol(object$coefficients_frontier))
  head_cv_v <- head_cv_u
  if (!is.null(object$coefficients_cm)) head_mu <- head_cv_u

  row.names(head_cv_u) <- "log(sigma2_u):"
  row.names(head_cv_v) <- "log(sigma2_v):"
  if (!is.null(object$coefficients_cm)) row.names(head_mu) <- "mu:"

  row.names(object$coefficients_cv_u) <- paste("  ", row.names(object$coefficients_cv_u))
  row.names(object$coefficients_cv_v) <- paste("  ", row.names(object$coefficients_cv_v))
  if (!is.null(object$coefficients_cm)) row.names(object$coefficients_cm)   <- paste("  ", row.names(object$coefficients_cm))

  outtable <- rbind(object$coefficients_frontier,
                    if (!is.null(object$coefficients_cm)) head_mu else NULL,
                    if (!is.null(object$coefficients_cm)) round(object$coefficients_cm, 3) else NULL,
                    head_cv_u,
                    object$coefficients_cv_u,
                    head_cv_v,
                    object$coefficients_cv_v)

  cat("Stochastic frontier model",
      "=========================", sep = "\n")

  cat(if (object$call$structure == "panel") "Panel" else "Cross-section", "data.\n")
  cat("Total observations:", object$N, if (object$call$structure == "panel") "in cross-section." else NULL, "\n")
  cat("Log-likelihood:", object$loglik, "\n")
  cat("\n", sep = "")
  print(round(outtable, digits = 4), digits = 4, nsmall = 3, na.print = "", quote = F)
  cat("\n", sep = "")
}

#' @export
print.SFA <- function(object) {
  print(summary.SFA(object))
}
