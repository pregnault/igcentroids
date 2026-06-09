# Header ----

# Desc.: Miscellanous functions used in package igcentroids.


## Binary search ----


#' Binary search algorithm
#'
#' An R implementation of the binary search of the zero of a
#' continous strictly monotone real function.
#'
#' @param lb A single numeric value, the lower bound of the interval to seek in the zero point.
#' @param ub A single numeric value, the upper bound of the interval to seek in the zero point.
#' @param fun A continuous monotone real function defined on the interval
#' `[lb, up]` such that `fun(lb)*fun(up) < 0`.
#' @param tol The precision of the approximation of the zero point to be attained.
#'
#' @returns A list with three components :
#'  - `sol`, the approximation of the zero point for function `fun` ;
#'  - `prec`, the precision of the approximation, of the form `(ub - lb)/2^(iter+1)` ;
#'  - `iter`, the number of iterations of the binary search while processing the approximation.
#' @export
#'
#' @examples
#' tmpfun <- function(x) return(sqrt(x)-2)
#' tmpub <- 9
#' tmplb <- 1
#' binarysearch(lb = tmplb, ub = tmpub, fun = tmpfun)
#'
binarysearch <- function(lb, ub, fun, tol = 1e-8) {
  if (fun(lb) * fun(ub) > 0 | fun(lb) > fun(ub) ) {
    stop("Function fun must be strictly monotone,
          with 0 belonging between the images by fun of
          lower and upper bound lb and up.")
  }
  if (isTRUE(all.equal(fun(lb), 0, tolerance = tol))) {
    return(list(sol = lb, prec = 0, iter = 0))
  }
  if (isTRUE(all.equal(fun(ub), 0, tolerance = tol))) {
    return(list(sol = ub, prec = 0, iter = 0))
  }
  if (fun(lb) < 0 && fun(ub) > 0) {
    lba <- lb
    uba <- ub
  }
  else {
    lba <- ub
    uba <- lb
  }
  # From there, fun is assumed to be a continusous strictly
  # increasing function with fun(lba) <0 and fun(uba) >0
  cpt <- 0
  while (uba - lba > 2 * tol) {
    cpt <- cpt + 1
    x <- (uba + lba) / 2
    if (fun(x) > 0) {uba <- x} else {lba <- x}
  }
  res <- (uba + lba) / 2
  err <- (uba - lba) / 2
  return(list(sol = res, prec = err, iter = cpt))
}
