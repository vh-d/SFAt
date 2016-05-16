# PREDICT FUNCTIONS --------------------------------------------------------
#' predict function for SFA models
#' @export
predict.SFA <- function(object,
                        newdata = NULL,
                        type = c("frontier", "efficiency", "inefficiency"),
                        estimator = c("JLMS", "ME", "BC")) {

  type <- match.arg(type)
  estimator <- match.arg(estimator)

  ret <- switch(type,
                frontier = predictFrontier(object, newdata),
                efficiency = efficiency.SFA(object, estimator = estimator),
                inefficiency = inefficiencyTerm.SFA(object,  estimator = estimator)
  )

  return(ret)
}

#' compute efficiency for SFA models
#' @export
efficiency.SFA <- function(object, estimator = "JLMS") {

  ineff_vec <- inefficiencyTerm.SFA(object, estimator)

  if (estimator == "BC") {
    return(exp(-ineff_vec))
  } else {
    pred_vec <- predict.SFA(object)
    return((pred_vec + object$ineff * ineff_vec) / pred_vec)
  }
}

#' predict function for SFA models
#' @export
inefficiencyTerm.SFA <- function(object, estimator) {
  spec <- object$call$spec
  dist <- object$call$dist
  structure <- object$call$structure

  u_fn_call <- paste0("u", "_",
                      structure, "_",
                      dist,
                      if (is.null(spec)) NULL else paste0("_", spec))

  # print(u_fn_call)
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

  return(as.vector(X %*% object$coeff))
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