# Header ----

# Desc.: Miscellanous functions used in package igcentroids.


# Binary search ----


#' Binary search algorithm
#'
#' An **R** implementation of the binary search of the zero of a
#' continous strictly monotone real function.
#'
#' @param lb A single numeric value, the lower bound of the interval to seek in the zero point.
#' @param ub A single numeric value, the upper bound of the interval to seek in the zero point.
#' @param fun A continuous monotone real function defined on the interval
#' `[lb, ub]` such that `fun(lb)*fun(ub) < 0`.
#' @param tol The precision of the approximation of the zero point to be attained.
#' Default value is `1e-8`.
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

# Support functions for computing jensen-shannon centroids of categorical distributions ----

## Gradient of the neg-entropy ----
grad <- function(theta) {
  return(log(theta / (1-sum(theta))))
}

## Inverse of the gradient of the neg-entropy ----
invgrad <- function(eta) {
  return(exp(eta)/(1+sum(exp(eta))))
}

# Integration into flexclust workflow ----

#' @name fcwintegration
#' @rdname fcwintegration
#'
#' @title Wrapped-up functions for seemless-integration in `{flexclust}` workflow.
#' @description
#' Functions `preproc_igc()`, `jdivflex()`, `jcentlex()`, `jsdivflex()`and `jscentlex()`
#' are intended to ease `igcentroids` integration into clustering workflow with `{flexclust}`.
#' @param p An object of class `igmatcatdistn`.
#' @param x A positive matrix.
#' @param weight A positive vector of length `n`.
#' The arithmetic mean is balanced by weighting each row of `p` with
#' the corresponding element of `weight`.
#' @param base A positive value, the base in which logarithms are computed.
#' @param method A character string specifying
#' the method to be used for computing the approximation of Jeffreys centroid. See [jcent.igcmatcatdistn()] for details.
#' @param centers A positive matrix whose rows are the (temporary) centers of clusters.
#' @seealso [flexclust::kccaFamily()]
#' @returns Objects suitable for defining a clustering family with `flexclust::kccaFamily()`.
NULL

#' @rdname fcwintegration
#' @export
#' @examples
#' # Create an object of class igcmatcatdistn
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
#' matp <- rbind(p, q, r)
#' matptest <- igcmatcatdistn(p = matp)
#' # Avoid class so as to return an object of type matrix
#' undo_igc(matptest)
undo_igc <- function(p) {
  if (! any(grepl(pattern = "igc", class(p)))) {
    stop("p must be an object of class igccatdist, igccatdistn, igcmatcatdist or igcmatcatdistn.")
  }
  res <- p
  attr(res, "class") <- NULL
  attr(res, "mass") <- NULL
  return(res)
}

#' @rdname fcwintegration
#' @export
#' @examples
#' # Compute Jeffreys divergences from a matrix as returned by undo_igc()
#' center <- matrix(rep(1/5, times = 5), nrow = 1)
#' jdivflex(x = matp, centers = center)
jdivflex <- function(x, centers, base = exp(1)) {
  jdiv(p = igcmatcatdistn(x), q = igcmatcatdistn(centers), base = base)
}

#' @rdname fcwintegration
#' @export
#' @examples
#' # Compute Jeffreys centroid from a matrix as returned by undo_igc()
#' jcentflex(x = matp)
jcentflex <- function(x, weight = NULL, method = "gbic", base = exp(1)) {
  jcent(p = igcmatcatdistn(x), weight = weight, method = method, base = base)
}

#' @rdname fcwintegration
#' @export
#' @examples
#' # Compute Jensen-Shannon divergences from a matrix as returned by undo_igc()
#' jsdivflex(x = matp, centers = center)
jsdivflex <- function(x, centers, base = exp(1)) {
  jsdiv(p = igcmatcatdistn(x), q = igcmatcatdistn(centers), base = base)
}

#' @rdname fcwintegration
#' @export
#' @examples
#' # Compute Jensen-Shannon centroid from a matrix as returned by undo_igc()
#' jscentflex(x = matp)
jscentflex <- function(x, weight = NULL, method = "gbic", base = exp(1)) {
  jscent(p = igcmatcatdistn(x), weight = weight, method = method, base = base)
}
