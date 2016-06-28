gradient <- function(fn, par) {
  n <- length(par)
  hh <- matrix(0, n, n)
  diag(hh) <- .Machine$double.eps^(1/3)

  res <- function(params,
                  indeces,
                  y, X,
                  CV_u,
                  CV_v,
                  CM,
                  ineff,
                  minmax,
                  deb){
    sapply(
      1:n,
      function(i) {
        ( do.call(what = fn, args = list(params + hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb)) -
          do.call(what = fn, args = list(params - hh[i, ], indeces, y, X, CV_u, CV_v, CM, ineff, minmax, deb))) / (2 * .Machine$double.eps^(1/3))
      }
    )
  }

  if (require(compiler)) res <- cmpfun(res)

  return(res)
}