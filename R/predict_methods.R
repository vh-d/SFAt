# PREDICT FUNCTIONS --------------------------------------------------------
#' predict function for SFA models
#' @param object    SFA object
#' @param estimator estimator of inefficiency term
#' @export
predict.SFA <- function(object,
                        newdata = NULL,
                        type = c("frontier", "efficiency", "inefficiency"),
                        estimator = c("JLMS", "ME", "BC")) {

  type <- match.arg(type)
  estimator <- match.arg(estimator)

  ret <- switch(type,
                frontier     = predictFrontier(object, newdata),
                efficiency   = efficiency.SFA(object,       estimator = estimator),
                inefficiency = inefficiencyTerm.SFA(object, estimator = estimator)
  )

  return(ret)
}

#' compute efficiency for SFA models
#' @rdname predict.SFA
#' @details
#' Currently, pedict methods assume that frontier model on \code{object} is specified in (trans)log form.
#' @export
efficiency.SFA <- function(object, estimator = "JLMS") {
  ineff_vec <- inefficiencyTerm.SFA(object, estimator)
  return(exp(object$ineff * ineff_vec))
}

#' predict function for SFA models
#' @rdname predict.SFA
#' @export
inefficiencyTerm.SFA <- function(object, estimator) {
  spec <- object$call$spec
  dist <- object$call$dist
  structure <- object$call$structure

  u_fn_call <- paste0("u", "_",
                      structure, "_",
                      dist,
                      if (is.null(spec)) NULL else paste0("_", spec))

  out <- do.call(u_fn_call, list(object = object, estimator = estimator))

  return(out)
}


predictFrontier <- function(object, newdata) {

  if (is.null(newdata)) {
    X <- object$data$X
  } else {
    X <- cbind(if (object$intercept) 1 else NULL,
               newdata)
  }

  return(as.vector(X %*% object$coeff_frontier))
}


# S3 methods -----------

#' @export
efficiency <- function(object, ...){
  UseMethod("efficiency")
}

#' @export
inefficiencyTerm <- function(object, ...){
  UseMethod("inefficiencyTerm")
}