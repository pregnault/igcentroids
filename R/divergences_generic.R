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
#' @details
#' For any finite-mass positive measures \eqn{p} and \eqn{q} on a measurable space,
#' such that \eqn{p} is dominated by \eqn{q} and both are dominated by some measure \eqn{\mu},
#' the Kullback-Leibler divergence \eqn{\mathbb{K}(p:q)} from \eqn{p} to \eqn{q} is given by
#' \deqn{\mathbb{K}(p:q) = \int p \log(p/q) \textrm{d} \mu - \int p \textrm{d} \mu +  \int q \textrm{d} \mu.}
#'
#' If \eqn{p} and \eqn{q} are probability distributions on a finite set \eqn{\mathcal{X}},
#' the Kullback-Leibler divergence from \eqn{p} to \eqn{q} reduces to
#' \deqn{\mathbb{K}(p:q) = \sum_{x \in \mathcal{X}} p(x) \log \frac{p(x)}{q(x)},}
#' with the conventions \eqn{0 \log \frac{0}{0} = 0}, \eqn{0 \log \frac{0}{q} = 0} and
#' \eqn{p \log \frac{p}{0} = \infty}.
#' @examples
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' kldiv(p, q)
kldiv <- function(p, q , base = exp(1)) {
  UseMethod("kldiv")
}

#' @title Jeffreys Divergence
#'
#' @description
#' `jdiv()` is an S3 generic for computing the Jeffreys divergence --
#' or symmetrized Kullback-Leibler divergence -- between
#' R objects, with methods for the following classes:
#' `r doclisting::methods_list('jdiv')`
#'
#' @param p,q Two **R** objects of same class. Currently, there are methods for
#' `igccatdist` and `igcmatcatdist` objects.
#' @param base A positive value specifying the base of logarithm used for
#' computing Jeffreys divergence. Default is exp(1) (natural logarithm).
#'
#' @returns
#'  - If `p` and `q` are `igccatdist` objects, returns the value of the Jeffreys divergence
#'  (or symmetrized Kullback-Leibler) between `p` and `q`;
#'  - if `p` and `q` are `igcmatcatdist` objects, returns a matrix of `nrow(p)` rows and
#'  `nrow(q)` columns whose coefficients are the Jeffreys divergence
#'  from rows of `p` to rows of `q`.
#' @details
#' The Jeffreys divergence or symmetric Kullback-Leibler divergence
#' is a symmetric and bounded discrimination measure between positive measures defined on
#' a measurable space.
#' It is formally defined, for \eqn{p} and \eqn{q} two positive measures with same support, by
#' \deqn{\mathbb{J}(p,q) = \mathbb{K}\left(p: q \right) +  \mathbb{K}\left(q:p \right),}
#' where \eqn{\mathbb{K}} denotes the Kullback-Leibler divergence; see [kldiv()].

#'
#' @seealso [kldiv()], [jsdiv()], [catdist]
#' @export
#'
#' @examples
#' # Compute the Jeffreys divergence
#' # between two probability distributions with common support.
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' jdiv(p, q)
#' # Compute the (squared) matrix of Jeffreys divergences between
#' # the elements of set of probabilily distributions
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
#' matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
#' matptest <- igcmatcatdistn(p = matp)
#' jdiv(p = matptest)
#' # Compute the matrix of Jeffreys divergences
#' # between two sets of probability distributions
#' matp2 <- matrix(c(dbinom(0:4, size = 5, prob = 0.3),
#'                       c(0.1, 0.1, 0.2, 0.2, 0.4)),
#'                     nrow = 2, byrow = TRUE)
#' matptest2 <- igcmatcatdist(matp2)
#' jdiv(p = matptest, q = matptest2)
jdiv <- function(p, q, base = exp(1)) {
  UseMethod("jdiv")
}

#' @title Jensen-Shannon divergence
#'
#' @description
#' `jsdiv()` is an S3 generic for computing the Jensen-Shannon divergence between
#' R objects, with methods for the following classes:
#' `r doclisting::methods_list('jdiv')`
#'
#' @param p,q Two **R** objects of same class. Currently, there are methods for
#' `igccatdist` and `igcmatcatdist` objects.
#' @param base A positive value specifying the base of logarithm used for
#' computing Jensen-Shannon divergence. Default is exp(1) (natural logarithm).
#'
#' @returns
#'  - If `p` and `q` are `igccatdist` objects, returns the value of the
#'  Jensen-Shannon divergence between `p` and `q`;
#'  - if `p` and `q` are `igcmatcatdist` objects, returns a matrix of `nrow(p)` rows and
#'  `nrow(q)` columns whose coefficients are the Jensen-Shannon divergences
#'  between rows of `p` and rows of `q`.
#' @details
#' The Jensen-Shannon divergence, also known as capacitory discrimination,
#' is a symmetric and bounded discrimination measure between finite-mass positive measures defined on a
#' measurable space.
#' It is formally defined, for \eqn{p} and \eqn{q} two finite-mass positive measures, by
#' \deqn{\mathbb{JS}(p,q) = \frac{1}{2} \left( \mathbb{K}\left(p: \frac{p+q}{2}  \right) +  \mathbb{K}\left(q: \frac{p+q}{2}  \right)\right),}
#' where \eqn{\mathbb{K}} denotes the Kullback-Leibler divergence; see [kldiv()].
#'
#' For probability distributions (pds), it is linked to Shannon entropy \eqn{\mathbb{S}} through the relation
#' \deqn{\mathbb{JS}(p,q) = \mathbb{S} \left( \frac{p+q}{2} \right) - \frac{\mathbb{S}(p) + \mathbb{S}(q)}{2}.}
#' Properties of Jensen-Shannon divergence include:
#'  - definition for pds with distinct supports:
#'  \eqn{\mathbb{JS}(p,q)} is well-defined wether or not \eqn{p} and \eqn{q} have same support,
#'  - symmetry: \eqn{\mathbb{JS}(p,q) = \mathbb{JS}(q,p)},
#'  - boundedness: \eqn{\mathbb{JS}(p,q) \leq \log(2)}, with equality whenever
#'  supports of \eqn{p} and \eqn{q} are disjoint.
#'
#' @seealso [kldiv()], [jdiv()], [catdist].
#' @export
#'
#' @examples
#' # Compute the Jensen-Shannon (JS) divergence
#' # between two probability distributions with common support.
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' jsdiv(p, q)
#' # Compute the (squared) matrix of JS divergences between
#' # the elements of set of probabilily distributions
#' p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
#' q <- igccatdist(rep(1/5, times = 5))
#' r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
#' matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
#' matptest <- igcmatcatdistn(p = matp)
#' jsdiv(p = matptest)
#' # Compute the matrix of JS divergences between
#' # two sets of probability distributions
#' matp2 <- matrix(c(dbinom(0:4, size = 5, prob = 0.3),
#'                       c(0.1, 0.1, 0.2, 0.2, 0.4)),
#'                     nrow = 2, byrow = TRUE)
#' matptest2 <- igcmatcatdist(matp2)
#' jsdiv(p = matptest, q = matptest2)
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
#' a matrix of categorical distribution; see [catdist].
#' @param weight A positive vector of length `n`.
#' The arithmetic mean is balanced by weighting each row of `p` with
#' the corresponding element of `weight`.
#'
#' @returns The arithmetic centroid of the input objects `p` weighted by `weight`.
#' @seealso [geomcent()], [jcent()], [jscent()], [catdist].
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

#' @title Geometric centroid of n objects
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
#' The geometric mean is balanced by weighting each row of `p` with
#' the corresponding element of `weight`.
#' @param base A positive value specifying the base of logarithm used
#' for computing Jeffreys divergence. Default is exp(1) (natural logarithm).
#'
#' @returns The geometric centroid of the input objects `p` weighted by `weight`.
#' @seealso [arithcent()], [jcent()], [jscent()], [catdist].
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
#' If input `p` is a positive matrix (object of class `igcmatcatdist`),
#' jcent()` returns the positive vector minimizing the weighted mean of
#' Jeffreys divergences with the rows of `p`.
#' More precisely, by denoting \eqn{p_1, \cdots, p_n} the rows of `p` and
#' \eqn{w_1, \cdots, w_n} the components of the vector `weight`, the output of jcent()` is
#' the positive vector \eqn{p^*} satisfying
#' \deqn{p^* \in \textrm{Argmin} \left\{ \sum_{i=1}^n w_i \mathbb{J}(p_i:q) \ | \  q \textrm{ pos. vector} \right\},}
#' where \eqn{\mathbb{J}((p_i:q)} denotes the Jeffreys divergence between \eqn{p_i} and \deqn{q};
#' see [jdiv()] for details.
#' The minimizer \eqn{p^*} has a closed-from expression given by
#' \deqn{p^*(x) = \frac{a(x)}{W( \exp( \log(b) a(x)/g(x)))},}
#' where \eqn{W} denotes the Lambert function,
#' \eqn{a = \sum w_i p_i} is the arithmetic centroid of the vectors \eqn{p_i},
#' \eqn{g = \prod p_i^{w_i}} is their geometric centroid and
#' \eqn{b >0} is the logarithm basis used to compute the Jeffreys divergence.
#'
#' If input `p` is a stochastic matrix (object of class `igcmatcatdistn`),
#' the output is the probability distribution minimizing the weighted mean of
#' Jeffreys divergences with the rows of `p`; formally,
#' \deqn{p^* \in \textrm{Argmin} \left\{ \sum_{i=1}^n w_i \mathbb{J}(p_i:q) |  q >= 0 \sum_{x} q(x) =1 \right\}.}
#' No closed-form expression exists for \eqn{p^*} is this case.
#' Several methods are implemented so as to compute numerical approximations;
#' see Details section in [jcent.igcmatcatdistn()].
#'
#' @seealso [arithcent()], [geomcent()], [jscent()], [lamW::lambertW0()], [jcent.igcmatcatdistn()]
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
#' @param tol For `igcmatcatdistn` objects only. A positive value,
#' the tolerance in the optimization procedure; see Details Section. Default is 1e-8.
#' @param ... Additionnal arguments for other methods.
#'
#' @returns A positive numeric vector summing to one  (object of class `igccatdistn`),
#'  i.e., a categorical probability distribution,
#'  approximation of the Jensen-Shannon centroid of the rows of `p`.
#'
#' @details
#' If input `p` is a stochastic matrix (object of class `igcmatcatdistn`),
#' the output is the probability distribution minimizing the weighted mean of
#' Jensen-Shannon divergences with the rows of `p`. More precisely,
#' let \eqn{p_1, \cdots, p_n} denote the rows of `p` and
#' \eqn{w_1, \cdots, w_n} the components of the vector `weight`.
#' The output of `jscent()` is a numerical approximation of the
#' probability distribution \eqn{p^*} (object of class `igccatdistn`) satisfying
#' \deqn{p^* \in \textrm{Argmin} \left\{ \sum_{i=1}^n w_i \ \mathbb{JS}(p_i:q) \ | \  q >= 0,  \sum_{x} q(x) =1 \right\},}
#' where \eqn{\mathbb{JS}(p_i:q)} denotes the Jensen-Shannon divergence between \eqn{p_i} and \eqn{q};
#' see [jsdiv()] for details.
#'
#' This numerical approximation is obtained by the ConCave-Convexe Procedure (CCCP) introduced by
#' Yuille and Rangarajan (2002), as stated in Nielsen (2020).
#' It is an iterative optimization procedure for objective functions that can be written as
#' the sum of two functions, one concave \eqn{B}, the other convex \eqn{A}.
#' Starting from an initial pd \eqn{p^{(0)}} (the arithmetic centroid of rows of `p`),
#' a sequence of pds \eqn{p^{(1)}, \cdots, p^{(r)}} is computed recursively
#' through the formula
#' \deqn{p^{(t+1)} = (\nabla B)^{-1} \left( \nabla A (p^{(t)})\right), \quad t =1, \cdots, r,}
#' where \eqn{\nabla A} and \eqn{\nabla B} denote gradients of functions \eqn{A} and  \eqn{B}.
#' The algorithm stops as soon as
#' the norm-1 between \eqn{p^{(r)}} and \eqn{p^{(r-1)}} is less
#' than the tolerance argument `tol` and outputs \eqn{p^{(r)}}.
#'
#' @seealso [arithcent()], [geomcent()], [jcent()], [catdist].
#' @export
#' @references
#'
#' Yuille, A.L., Rangarajan, A. (2001). The concave-convex procedure (CCCP).
#' In Proceedings of the Neural Information Processing Systems 2001,
#' Vancouver, BC, Canada, 9--14, December 2001, 1033-1040.
#' [doi:10.5555/2980539.2980672](https://dl.acm.org/doi/10.5555/2980539.2980672)
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
jscent <- function(p, weight = NULL, base = exp(1), tol = 1e-8, ...) {
  UseMethod("jscent")
}
