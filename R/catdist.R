# Header ----

## Desc.: Implementation of S3 class igccatdist, i.e.,
## finite support, positive measures.

# Constructors ----

## Classes igccatdist0, igccatdist and igccatdistn ----


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
  if(ncol(p) == 1) {
    res <- structure(p/p,
                     class = c("igcmatcatdistn", "igcmatcatdist"),
                     mass = rep(1, times = nrow(p)))
  }
  else {
    res <- structure(t(apply(p, 1, function(row) row / sum(row))),
                     class = c("igcmatcatdistn", "igcmatcatdist"),
                     mass = rep(1, times = nrow(p)))
  }
  return(res)
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

#' @export
jsdiv.igccatdist <- function(p, q, base = exp(1)) {
  if (! ("igccatdist" %in% class(p) && "igccatdist" %in% class(q) && length(p) == length(q))) {
    stop("Arguments p and q must be igccatdist objects of same size")
  }
  return((kldiv(p, (p+q)/2, base = base) + kldiv(q, (p+q)/2, base = base))*.5)
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


#' Compute Jeffreys centroids of a set of normalized categorical distributions.
#'
#' Method for the generic `jcent()` and the class `igcmatcatdistn`.
#'
#' @export
#' @param p An object of classes `igcmatcatdist` or `igcmatcatdistn`, i.e.,
#' a positive or stochastic matrix (for the former, rows sum to one).
#' @param weight A vector of weights for the rows of `p`.
#' @param base A positive value, the base of the logarithm to be used
#' to compute Jeffreys divergence.
#' @param method A character string specifying
#' the method to be used for computing the approximation of Jeffreys centroid.
#' Possible values are:
#'  - "mag" for the mean of arithmetic and normalized geometric centroids;
#'  - "njc" for the normalized Jeffreys centroid of positive vectors;
#'  - "jfr" for the Jeffreys-Fisher-Rao centroid, i.e., the Fisher-Rao mean of arithmetic
#' and normalized geometric centroids;
#'  - "bsa" for a single-value binary search alogrithm based on Lagrangian optimization;
#'  - "gbic" for the Gauss-Bregman inductive centroid.
#' Default is "bsa". See Details section.
#' @param tol For `igcmatcatdistn` objects only. A positive value,
#' the tolerance in the approximation of Jeffreys centroid. Default is 1e-8.
#' @param ... Additionnal arguments for other methods.
#' @seealso Generic function [jcent()]
#' @references
#'
#' Nielsen, F. (2013). Jeffreys Centroids: A Closed-From Expresion for
#' Positive Histograms and a Guaranteed Tight Approximation for
#' Frequency Histograms.
#' *IEEE Signal Processing Letters*, **20**(7), 657--660.
#' [doi:10.1109/LSP.2013.2260538](https://ieeexplore.ieee.org/document/6509930)
#'
#' Nielsen, F. (2024). Fast Proxy Centers for the Jeffreys Centroid:
#' The Jeffreys-Fisher-Rado Center and the Gauss-Bregman Inductive Center.
#' *Entropy*, **26**, 1008-1037.
#' [doi:10.3390/e26121008](https://www.mdpi.com/1099-4300/26/12/1008)
#'
#' @examples
#' u <- rep(1/3, 3)
#' prob <- c(0.9, 0.05, 0.05)
#' p <- matrix(c(u, prob),
#'              nrow = 2, ncol =3, byrow = TRUE)
#' p <- igcmatcatdistn(p)
#' meths <- c("mag", "njc", "jfr", "bsa", "gbic")
#' t(sapply(meths, function(m) jcent(p, method = m)))

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
    while (sum(abs(a - g)) > 2 * tol && iter < 10000) {
      na <- (a + g) / 2
      ng <- sqrt(a * g)
      ng <- ng / sum(ng)
      a <- na
      g <- ng
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

#' @export
jscent.igcmatcatdistn <- function(p, weight = NULL, base = exp(1), tol = 1e-8, ...) {
  np <- dim(p)[1] # Number of points
  D <- dim(p)[2] # Number of states of the distributions
  if (is.null(weight)) {
    weight <- rep(1, times = np)
  }
  if (! (is.numeric(weight) && length(weight) == np && all(weight > 0))) {
    stop("Argument weight must be a positive vector of length
         the number of rows of the igcmatcatdist object p.")
  }
  if (! all.equal(weight, rep(1, times = np))) {
    stop("Only uniform weights are supported for now.")
  }
  # Initialize CCCP algorithm with the arithmetic centroid
  ptemp <- arithcent(p, weight=weight)
  # Compute the loss function for this initial value
  obj_actu <- sum(apply(p, 1, function(cd) jsdiv(igccatdist(cd), ptemp)))
  # CCCP algorithm : iterate procedure until convergence
  repeat {
    theta_actu <- ptemp[1:(D-1)] # Mixture coordinates
    # First step : gradient of the loss function at the current value
    grads <- apply(p, 1, function(cd) {
      theta <- cd[1:(D-1)]
      return(grad((theta + theta_actu)/2))
    })
    # Second step : Inverse of the gradient at the arithmetic mean of the gradients
    theta_new <- invgrad(apply(grads, 1, function(col) mean(col)))
    # Convert back to catdist.
    ptemp <- igccatdistn(c(theta_new, 1- sum(theta_new)))
    if (sum(abs(theta_new - theta_actu)) < tol) {break}
  }
  return(ptemp)
}



# Plot methods for categorical distributions ----

#' Scatter plot of normalized categorical  distributions
#'
#' Methods `plot()` and `points()` for objects of classes `igccatdistn` and `igcmatcatdistn`.
#' The type of the chart produced depends on the size of distributions to be drawn.
#' More precisely, two-state distributions are drawn as points on a unit segment,
#' three-state distributions are drawn as points in the equilateral triangle.
#' Distributions with more than three states are drawn as polylines.
#'
#' @param x An object of class igcmatcatdistn, i.e., a matrice with positive entries and
#' rows sum to 1.
#' @param y NULL. Argument present only for compatibility with generic function plot.
#' @param add Logical. Indicates whether the normalized categorical distributions should
#' be drawn in a new chart (FALSE) or should be added to the active chart (TRUE).
#' @param pch An integer value specifying the shape point
#' @param ... Further arguments to be passed to default methods for `plot()` and `points()`;
#' see `?par` for details.
#' @importFrom graphics axis lines plot.new plot.window points segments
#'
#' @returns NULL
#' @export
#'
#' @examples
#' # Two-state  categorical distributions plot
#' u <- c(0.5, 0.5)
#' p1 <- c(0.3, 0.7)
#' p2 <- c(0.4, 0.6)
#' x <- igcmatcatdistn(matrix(c(u, p1, p2),
#'                     nrow = 3, ncol = 2, byrow = TRUE))
#' plot(x)
#' # Three-state categorical distributions plot
#' u <- rep(1/3, 3)
#' prob <- c(0.99, 0.005, 0.005)
#' x <- matrix(c(u, prob),
#'             nrow = 2, ncol =3, byrow = TRUE)
#' x <- igcmatcatdistn(x)
#' plot(x)
#' # Many-state categorical distributions plot
#' set.seed(745)
#' var1 <- factor(sample(x = LETTERS[1:4], size = 50L,
#'                    replace = TRUE, prob = c(0.2, 0.3, 0.4, 0.1)),
#'             levels = LETTERS[1:4], ordered = TRUE)
#' p1 <- table(var1)
#' var2 <- factor(sample(x = LETTERS[1:4], size = 50L,
#'                    replace = TRUE, prob = c(0.4, 0.2, 0.1, 0.3)),
#'             levels = LETTERS[1:4], ordered = TRUE)
#' p2 <- table(var2)
#' p <- igcmatcatdistn(matrix(c(p1/sum(p1), p2/sum(p2)),
#'                            nrow = 2, ncol = 4, byrow = TRUE))
#' plot(p)
#' axis(side = 1, at = 1:4, labels = LETTERS[1:4])
#'
plot.igcmatcatdistn <- function(x, y = NULL, add = FALSE, pch = 20, ...) {
  if (ncol(x) < 2) {
    stop("Plot method for igcmatcatdistn available only for probability distributions with two states or more. Input has ", ncol(x), " states.")
  }
  if (ncol(x) == 2) {
    # New plot with a segment representing the set of two-state p.d.
    if (!add) {
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(-0.1, 0.1), asp = 1)
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0,
               lty = 1, lwd = 2, col = "gray")
      segments(x0=c(0, 1), y0 = c(-0.05, -0.05), x1 = c(0, 1), y1 = c(0.05, 0.05),
               lty = 1, lwd = 2, col = "gray")
    }
    # remove structure from igcmatcatdistn object so as to use default points() method
    coords <- x
    attr(coords, "class")  <- NULL
    points(x = coords[,1], y = rep(0, times = nrow(x)), pch = pch, ...)
    typeplot <- "1D segment plot"
    coordsplot <- coords[,1]
  }
  if (ncol(x) == 3) {
    if (!add) {
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(0, sqrt(3)/2), asp = 1)
      segments(x0 = c(0, 1, 1/2), y0 = c(0, 0, sqrt(3)/2), x1 = c(1, 1/2, 0), y1 = c(0, sqrt(3)/2, 0),
               lty = 1, lwd = 2)
    }
    plan2triangle <- matrix(c(1, 0, 1/2, sqrt(3)/2), nrow = 2,ncol = 2, byrow = FALSE)
    coords_simplex <- plan2triangle %*% t(x[, 1:2, drop = FALSE])
    points(t(coords_simplex),
           pch = pch, ...)
    typeplot <- "Equilateral triangle plot"
    coordsplot <- t(coords_simplex)
  }
  if (ncol(x) > 3) {
    if (!add) {
      plot.new()
      plot.window(xlim = c(1,ncol(x)), ylim = c(0, max(x)))
      axis(side = 2)
      apply(x, 1, lines, ...)
    }
    typeplot <- "Polyline plot"
    coordsplot <- NULL
  }
  invisible(list(type = typeplot, coords = coordsplot))
}

#' @rdname plot.igcmatcatdistn
#' @export
#' @examples
#' u <- rep(1/3, 3)
#' prob <- c(0.99, 0.005, 0.005)
#' x <- matrix(c(u, prob),
#'             nrow = 2, ncol =3, byrow = TRUE)
#' x <- igcmatcatdistn(x)
#' plot(x)
#' y <- matrix(c(0.7, 0.15, 0.15,
#'               0.5, 0.25, 0.25),
#'               nrow = 2, ncol = 3, byrow = TRUE)
#' y <- igcmatcatdistn(y)
#' points(y, pch = 3)

points.igcmatcatdistn <- function(x, y = NULL, ...) {
  plot.igcmatcatdistn(x, add = TRUE, ...)
}

#' @rdname plot.igcmatcatdistn
#' @export
#' @examples
#' u <- rep(1/3, 3)
#' prob <- c(0.99, 0.005, 0.005)
#' x <- matrix(c(u, prob),
#'             nrow = 2, ncol =3, byrow = TRUE)
#' x <- igcmatcatdistn(x)
#' plot(x)
#' acent <- arithcent(p = x)
#' acent <- igccatdistn(acent)
#' gcent <- geomcent(p = x)
#' gcent <- igccatdistn(gcent/sum(gcent))
#' points(x = acent, pch = 3, col = "red")
#' points(x = gcent, pch = 3, col = "blue")
points.igccatdistn <- function(x, y = NULL, ...) {
  dim(x) <- c(1,length(x))
  x <- igcmatcatdistn(x)
  points.igcmatcatdistn(x, ...)
}


