# Jensen-Shannon centroid of categorical probability distributions.

`jscent()` is an S3 generic for computing the Jensen-Shannon centroid of
R objects, with methods for the following classes:

- `igcmatcatdistn`

## Usage

``` r
jscent(p, weight = NULL, base = exp(1), tol = 1e-08, ...)
```

## Arguments

- p:

  An object of class `igcmatcatdistn`, i.e., a stochastic matrix
  (positive matrix with row sums equal to one).

- weight:

  A vector of weights for the rows of `p`. Default is NULL, yielding
  uniform weights on the rows of `p`.

- base:

  A positive value, the base of the logarithm to be used to compute
  Jensen-Shannon divergence. Default is `exp(1)`.

- tol:

  For `igcmatcatdistn` objects only. A positive value, the tolerance in
  the optimization procedure; see Details Section. Default is 1e-8.

- ...:

  Additionnal arguments for other methods.

## Value

A positive numeric vector summing to one (object of class
`igccatdistn`), i.e., a categorical probability distribution,
approximation of the Jensen-Shannon centroid of the rows of `p`.

## Details

If input `p` is a stochastic matrix (object of class `igcmatcatdistn`),
the output is the probability distribution minimizing the weighted mean
of Jensen-Shannon divergences with the rows of `p`. More precisely, let
\\p_1, \cdots, p_n\\ denote the rows of `p` and \\w_1, \cdots, w_n\\ the
components of the vector `weight`. The output of `jscent()` is a
numerical approximation of the probability distribution \\p^\*\\ (object
of class `igccatdistn`) satisfying \$\$p^\* \in \textrm{Argmin} \left\\
\sum\_{i=1}^n w_i \\ \mathbb{JS}(p_i:q) \\ \| \\ q \>= 0, \sum\_{x} q(x)
=1 \right\\,\$\$ where \\\mathbb{JS}(p_i:q)\\ denotes the Jensen-Shannon
divergence between \\p_i\\ and \\q\\; see
[`jsdiv()`](https://pregnault.github.io/igcentroids/reference/jsdiv.md)
for details.

This numerical approximation is obtained by the ConCave-Convexe
Procedure (CCCP) introduced by Yuille and Rangarajan (2002), as stated
in Nielsen (2020). It is an iterative optimization procedure for
objective functions that can be written as the sum of two functions, one
concave \\B\\, the other convex \\A\\. Starting from an initial pd
\\p^{(0)}\\ (the arithmetic centroid of rows of `p`), a sequence of pds
\\p^{(1)}, \cdots, p^{(r)}\\ is computed recursively through the formula
\$\$p^{(t+1)} = (\nabla B)^{-1} \left( \nabla A (p^{(t)})\right), \quad
t =1, \cdots, r,\$\$ where \\\nabla A\\ and \\\nabla B\\ denote
gradients of functions \\A\\ and \\B\\. The algorithm stops as soon as
the norm-1 between \\p^{(r)}\\ and \\p^{(r-1)}\\ is less than the
tolerance argument `tol` and outputs \\p^{(r)}\\.

## References

Yuille, A.L., Rangarajan, A. (2001). The concave-convex procedure
(CCCP). In Proceedings of the Neural Information Processing Systems
2001, Vancouver, BC, Canada, 9–14, December 2001, 1033-1040.
[doi:10.5555/2980539.2980672](https://dl.acm.org/doi/10.5555/2980539.2980672)

Nielsen, F. (2020). On a Generalization of the Jensen-Shannon Divergence
and the Jensen-Shannon Centroid. *Entropy*, **22**, 221–244.
[doi:10.3390/e22020221](https://www.mdpi.com/1099-4300/22/2/221)

## See also

[`arithcent()`](https://pregnault.github.io/igcentroids/reference/arithcent.md),
[`geomcent()`](https://pregnault.github.io/igcentroids/reference/geomcent.md),
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md),
[catdist](https://pregnault.github.io/igcentroids/reference/catdist.md).

## Examples

``` r
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
matptest <- igcmatcatdistn(p = matp)
jscent(matptest)
#> [1] 0.1821431 0.2907138 0.2210574 0.1863583 0.1197274
#> attr(,"class")
#> [1] "igccatdistn" "igccatdist" 
#> attr(,"mass")
#> [1] 1
```
