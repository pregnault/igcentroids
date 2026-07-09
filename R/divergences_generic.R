# Generics for divergences, including Kullback-Leibler, Jeffreys.

#' Kullback-Leibler Divergence
#'
#' Compute the Kullback-Leibler divergence of a positive measure with respect to another dominating the first one.
#'
#' @param p,q Two **R** objects of same class. Currently,
#' there is a method for `igccatdist` objects.
#' @param base A positive value, the base in which logarithms are computed.
#'
#' @returns The value of the Kullback-Leibler divergence from `p` to `q`.
#' @export
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' kldiv(p, q)
kldiv <- function(p, q , base = exp(1)) {
  UseMethod("kldiv")
}

#' Jeffreys Divergence
#'
#' Compute the Jeffreys divergence between two positive mesasures with same support i.e.,
#' the symmetrization of Kullback-Leibler divergence.
#'
#' @param p,q Two **R** objects of same class. Currently, there is a method for
#' `igccatdist` objects.
#' @param base A positive value specifying the base of logarithm used for
#' computing Jeffreys divergence. Default is exp(1) (natural logarithm).
#'
#' @returns The value of the Jeffreys divergence (or symmetrized Kullback-Leibler)
#' from `p` to `q`.
#' @seealso [kldiv()], [jsdiv()]
#' @export
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' jdiv(p, q)
jdiv <- function(p, q, base = exp(1)) {
  UseMethod("jdiv")
}

#' Jensen-Shannon divergence
#'
#' Compute the Jensen-Shannon divergence between two probability with same support,
#' i.e., the average of Kullback-Leibler divergences of probability distributions
#' with respect to the mixture of them.
#'
#' @param p,q Two **R** objects of same class. Currently, there is a method for
#' `igccatdist` objects.
#' @param base A positive value specifying the base of logarithm used for
#' computing Jeffreys divergence. Default is exp(1) (natural logarithm).
#'
#' @returns The value of (symmetric) the Jensen-Shannon divergence (or capacity discrimination)
#' between `p` and `q`.
#' @seealso [kldiv()], [jdiv()]
#' @export
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' jsdiv(p, q)

jsdiv <- function(p, q, base =exp(1)) {
  UseMethod("jsdiv")
}


#' @title Arithmetic centroid of n objects
#' @description
#' `arithcent()` is an S3 generic for computing the arithmetic mean of R objects,
#' with methods for the following classes:
#'
#' `r doclisting::methods_list('arithcent')`
#'
#' @param p An **R** object gathering a number `n` of objects of same nature.
#' Currently, there is a method for `igcmatcatdist` matrix, i.e.,
#' a matrix of categorical distribution; see.
#' @param weight A positive vector of length `n`.
#'
#' @returns The arithmetic centroid of the input objects `p` weighted by `weight`.
#' @seealso [geomcent()], [jcent()], [jscent()]
#' @export
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' matp <- matrix(c(p, q), nrow = 2)
#' matptest <- igcmatcatdist(p = matp)
#' arithcent(p = matptest)
arithcent <- function(p, weight = NULL) {
  UseMethod("arithcent")
}

#' @title Geometrix centroid of n objects
#' @description
#' `geomcent()` is an S3 generic for computing the geometric mean of R objects,
#' with methods for the following classes:
#'
#' `r doclisting::methods_list('geomcent')`
#'
#' @param p An **R** object gathering a number `n` of objects of same nature.
#' Currently, there is a method for `igcmatcatdist` matrix, i.e.,
#' a matrix of categorical distribution.
#' @param weight A positive vector of length `n`.
#' @param base A positive value specifying the base of logarithm used
#' for computing Jeffreys divergence. Default is exp(1) (natural logarithm).
#'
#' @returns The geometric centroid of the input objects `p` weighted by `weight`.
#' @seealso [arithcent()], [jcent()], [jscent()]
#' @export
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' matp <- matrix(c(p, q), nrow = 2)
#' matptest <- igcmatcatdist(p = matp)
#' geomcent(p = matptest)
geomcent <- function(p, weight = NULL, base = exp(1)) {
  UseMethod("geomcent")
}

#' @title Jeffreys centroid of positive vectors and
#' categorical probability distributions.
#' @description
#' `jcent()` is an S3 generic for computing the Jeffreys centroid of R objects,
#' with methods for the following classes:
#'
#' `r doclisting::methods_list('jcent')`
#'
#' @param p An object of classes `igcmatcatdist` or `igcmatcatdistn`, i.e.,
#' a positive or stochastic matrix (for the former, rows sum to one).
#' @param weight A vector of weights for the rows of `p`.
#' @param base A positive value, the base of the logarithm to be used
#' to compute Jeffreys divergence.
#' @param ... Additionnal arguments for other methods.
#'
#' @returns A positive vector (object of class `igccatdist`) if `p`
#' is of class `igcmatcatdist` or a categorical probability distribution
#'  (object of class `igccatdistn`), approximation of the Jeffreys centroid
#'  of the rows of `p`.
#'
#' @details
#' Additional details...
#'
#' @seealso [arithcent()], [geomcent()], [jscent()], [lamW::lambertW0()]
#' @export
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

jcent <- function(p, weight = NULL, base = exp(1), ...) {
  UseMethod("jcent")
}

#' @title Jensen-Shannon centroid of categorical probability distributions.
#' @description
#' `jscent()` is an S3 generic for computing the Jensen-Shannon centroid of R objects,
#' with methods for the following classes:
#'
#' `r doclisting::methods_list('jscent')`
#'
#' @param p An object of class `igcmatcatdistn`, i.e.,
#' a stochastic matrix (positive matrix with row sums equal to one).
#' @param weight A vector of weights for the rows of `p`.
#' Default is NULL, yielding uniform weights on the rows of `p`.
#' @param base A positive value, the base of the logarithm to be used
#' to compute Jensen-Shannon divergence. Default is `exp(1)`.
#' @param ... Additionnal arguments for other methods.
#'
#' @returns A positive numeric vector summing to one  (object of class `igccatdistn`),
#'  i.e., a categorical probability distribution,
#'  approximation of the Jensen-Shannon centroid of the rows of `p`.
#'
#' @details
#' Additional details...
#'
#' @seealso [arithcent()], [geomcent()], [jcent()]
#' @export
#' @references
#'
#' Nielsen, F. (2020). On a Generalization of the Jensen-Shannon Divergence
#' and the Jensen-Shannon Centroid.
#' *Entropy*, **22**, 221--244.
#' [doi:10.3390/e22020221](https://www.mdpi.com/1099-4300/22/2/221)
#'
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
#' matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
#' matptest <- igcmatcatdistn(p = matp)
#' jscent(matptest)
jscent <- function(p, weight = NULL, base = exp(1), ...) {
  UseMethod("jscent")
}
