# Header ----

## Desc.: Implementation of S3 class igccatdist, i.e.,
## finite support, positive measures.

# Constructors ----

## Classes igccatdist and igccatdistn ----

new_igccatdist <- function(p) {
  stopifnot(is.vector(p), is.numeric(p), length(p) > 0, all(p > 0))
  structure(p, class = "igccatdist", mass = sum(p))
}

new_igccatdistn <- function(p) {
  stopifnot(is.vector(p), is.numeric(p), length(p) > 0, all(p > 0))
  structure(p / sum(p), class = c("igccatdistn", "igccatdist"), mass = sum(p))
}

## Class igcmatcatdist ----

new_igcmatcatdist <- function(p) {
  stopifnot(is.matrix(p), is.numeric(p), all(dim(p)) > 0, all(p > 0))
  structure(p, class = "igcmatcatdist", mass = apply(p, 1, sum))
}

new_igcmatcatdistn <- function(p) {
  stopifnot(is.matrix(p), is.numeric(p), all(dim(p)) > 0, all(p > 0))
  structure(t(apply(p, 1, function(row) row / sum(row))),
            class = c("igcmatcatdistn", "igcmatcatdist"),
            mass = rep(1, times = nrow(p)))
}

# Validators ----

## Classes igccatdist and igccatdistn ----

validate_igccatdist <- function(p) {
  mass_p <- attr(p, "mass")
  if (!(is.numeric(p) && length(p) > 0 && all(p > 0))) {
    stop("Objet of class igccatdist must be a non-zero length
         numeric vector with positive values.",
         call. = FALSE)
  }
  if (!(is.numeric(mass_p) && length(mass_p) == 1 && mass_p > 0)) {
    stop("Attribute mass of  class igccatdist must be a positive value.",
         call. = FALSE)
  }
  return(p)
}

validate_igccatdistn <- function(p) {
  mass_p <- attr(p, "mass")
  if (!(is.numeric(p) && length(p) > 0 && all(p > 0))) {
    stop("Objet of class igccatdist must be a non-zero length
         numeric vector with positive values.",
         call. = FALSE)
  }
  if (!(is.numeric(mass_p) && length(mass_p) == 1 && isTRUE(all.equal(mass_p, 1)))) {
    stop("Attribute mass of  class igccatdistn must equal 1.", call. = FALSE)
  }
  return(p)
}

## Class igcmatcatdist ----

validate_igcmatcatdist <- function(p) {
  mass_p <- attr(p, "mass")
  if (!(is.matrix(p) && is.numeric(p) && all(dim(p) > 0) && all(p > 0))) {
    stop("Objet of class igcmatcatdist must be a positive matrix.", call. = FALSE)
  }
  if (!(is.numeric(mass_p) && length(mass_p) == dim(p)[1])) {
    stop("Attribute mass of class igcmatcatdist must be a vector of positive values,
         whose length equals the number of rows of the object.",
         call. = FALSE)
  }
  return(p)
}

validate_igcmatcatdistn <- function(p) {
  mass_p <- attr(p, "mass")
  if (!(is.matrix(p) && is.numeric(p) && all(dim(p) > 0) && all(p > 0))) {
    stop("Objet of class igcmatcatdist must be a positive matrix.", call. = FALSE)
  }
  if (!(is.numeric(mass_p) && length(mass_p) == dim(p)[1] &&
          all(all(isTRUE(all.equal(mass_p, rep(1L, times = nrow(p)))))))) {
    stop("Attribute mass of class igcmatcatdistn must be a
         constant vector of ones, whose length equals
         the number of rows of the object.",
         call. = FALSE)
  }
  return(p)
}
# Helpers ----

## Classes igccatdist and igccatdistn ----

#' @name catdist
#' @rdname catdist
#'
#' @title Classes igccatdist, igccatdistn, igcmatcatdist, igcmatcatdistn
#' for positive (normalized) vectors and matrices.
#' @param p A numeric vector (for `igcatdist`, `igcatdistn`) or
#' matrix (`igmatcatdist`, `igmatcatdistn`) of positive values.
#'
#' @returns An object of the eponyme class.
NULL

#' @rdname catdist
#' @export
#' @examples
#' igccatdist(c(0.3, 0.1, 0.4, 0.8))
igccatdist <- function(p = numeric()) {
  p <- as.numeric(p)
  validate_igccatdist(new_igccatdist(p))
}

#' @rdname catdist
#' @export
#' @examples
#' p <- igccatdistn(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdistn(rep(1/5, times = 5))
#' p+q

igccatdistn <- function(p = numeric()) {
  p <- as.numeric(p)
  validate_igccatdistn(new_igccatdistn(p))
}

## Class igcmatcatdist ----


#' @rdname catdist
#' @export
#' @examples
#' igcmatcatdist(matrix(c(dbinom(0:4, size = 4, prob = 0.5), 1:5),
#'                      nrow = 2, byrow = TRUE))

igcmatcatdist <- function(p = numeric()) {
  validate_igcmatcatdist(new_igcmatcatdist(p))
}

#' @rdname catdist
#' @export
#' @examples
#' matp <- matrix(c(dbinom(0:4, size = 4, prob = 0.5),
#'                  rep(1/5, times = 5)),
#'                nrow = 2, byrow = TRUE)
#' matptest <- igcmatcatdistn(p = matp)

igcmatcatdistn <- function(p = numeric()) {
  res <- validate_igcmatcatdistn(new_igcmatcatdistn(p))
  attributes(res)$mass <- rep(1, time = nrow(p))
  return(res)
}


# Methods for classes igccatdist and igccatdistn ----

## Methods for computing divergences ----

#' @export
kldiv.igccatdist <- function(p, q, base = exp(1)) {
  if (! ("igccatdist" %in% class(p) && "igccatdist" %in% class(q) && length(p) == length(q))) {
    stop("Arguments p and q must be igccatdist objects of same size")
  }
  res <- sum(p * log(p / q, base = base) - p + q)
  return(res)
}

#' @export
jdiv.igccatdist <- function(p, q, base = exp(1)) {
  if (! ("igccatdist" %in% class(p) && "igccatdist" %in% class(q) && length(p) == length(q))) {
    stop("Arguments p and q must be igccatdist objects of same size")
  }
  res <- sum((p - q) * log(p / q, base = base))
  return(res)
}

## Methods for computing centroids ----

#' @export
arithcent.igcmatcatdist <- function(p, weight = NULL) {
  np <- dim(p)[1] # Number of points
  if (is.null(weight)) {
    weight <- rep(1, times = np)
    }
  if (! (is.numeric(weight) && length(weight) == np && all(weight > 0))) {
    stop("Argument weight must be a positive vector of length,
         the number of rows of the igcmatcatdist object p.")
  }
  res <- apply(p, 2, function(col) sum(col * weight / sum(weight)))
  return(igccatdist(res))
}

#' @export
geomcent.igcmatcatdist <- function(p, weight = NULL, base = exp(1)) {
  np <- dim(p)[1] # Number of points
  if (is.null(weight)) {
    weight <- rep(1, times = np)
  }
  if (! (is.numeric(weight) && length(weight) == np && all(weight > 0))) {
    stop("Argument weight must be a positive vector of length,
         the number of rows of the igcmatcatdist object p.")
  }
  res <- apply(p, 2, function(col) exp(sum(log(col, base = base) * weight / sum(weight))))
  return(igccatdist(res))
}

#' @export
#' @importFrom lamW lambertW0
jcent.igcmatcatdist <- function(p, weight = NULL, base = exp(1), lambda = 0, ...) {
  a <- arithcent.igcmatcatdist(p, weight = weight)
  g <- geomcent.igcmatcatdist(p, weight = weight, base = base)
  res <- a * log(base) / lambertW0(a / g * exp(log(base) * (lambda + 1)))
  return(igccatdist(res))
}

#' @export
#' @param method For `igcmatcatdistn` objects only. A character string specifying
#' the method to be used for computing the approximation of Jeffreys centroid.
#' Possible values are:
#'  - "mag" for the mean of arithmetic and normalized geometric centroids;
#'  - "njc" for the normalized Jeffreys centroid of positive vectors;
#'  - "jfr" for the Jeffreys-Fisher-Rao centroid, i.e., the Fisher-Rao mean of arithmetic
#' and normalized geometric centroids;
#'  - "bsa" for a single-value binary search alogrithm based on Lagrangian optimization;
#'  - "gbic" for the Gauss-Bregman inductive centroid.
#'  Default is "bsa". See Details section.
#' @param tol For `igcmatcatdistn` objects only. A positive value,
#' the tolerance in the approximation of Jeffreys centroid. Default is 1e-8.
jcent.igcmatcatdistn <- function(p, weight = NULL, base = exp(1),
                                 method = "bsa", tol = 1e-8, ...) {
  # agm method: arithmetic mean of arithmetic and normalized geometric centroids
  if (method == "mag") {
    a <- arithcent.igcmatcatdist(p, weight = weight)
    g <- geomcent.igcmatcatdist(p, weight = weight, base = base)
    res <- (a + g / sum(g)) / 2
  }
  # njc method: normalized jeffreys centroid of positive measures catdist
  if (method == "njc") {
    j <- jcent.igcmatcatdist(p, weight = weight, base = base)
    res <- j / sum(j)
    attributes(res)$mass <- sum(res)
  }
  # bsa method: binary search approximation
  if (method == "bsa") {
    a <- arithcent.igcmatcatdist(p, weight = weight)
    g <- geomcent.igcmatcatdist(p, weight = weight, base = base)
    # Lower bound for lagrangian multiplier:
    lb <- max(log(exp(a) * (g / sum(g)))) - 1
    # Upper bound:
    ub <- 0
    objfun <- function(lambda) {
      j <- jcent.igcmatcatdist(p, weight, base, lambda)
      return(-sum(j) + 1)
    }
    bsres <- binarysearch(lb, ub, objfun, tol = tol)
    lsol <- bsres$sol
    res <- jcent.igcmatcatdist(p, weight, base, lsol)
  }
  if (method == "gbic") {
    a <- arithcent.igcmatcatdist(p, weight = weight)
    g <- geomcent.igcmatcatdist(p, weight = weight, base = base)
    g <- g / sum(g)
    iter <- 0
    while (sum(abs(a - g)) > 2 * tol && iter < 1000) {
      a <- (a + g) / 2
      g <- sqrt(a * g)
      g <- g / sum(g)
      iter <- iter + 1
    }
    res <- a
    attributes(res)$mass <- sum(res)
  }
  if (method == "jfr") {
    a <- arithcent.igcmatcatdist(p, weight = weight)
    g <- geomcent.igcmatcatdist(p, weight = weight, base = base)
    g <- g / sum(g)
    res <- (sqrt(a) + sqrt(g))^2
    res <- res / sum(res)
    attributes(res)$mass <- sum(res)
  }
  return(igccatdistn(res))
}
