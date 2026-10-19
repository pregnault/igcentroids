# Compute Jeffreys centroids of a set of normalized categorical distributions.

Method for the generic
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md)
and the class `igcmatcatdistn`.

## Usage

``` r
# S3 method for class 'igcmatcatdistn'
jcent(p, weight = NULL, base = exp(1), method = "gbic", tol = 1e-08, ...)
```

## Arguments

- p:

  An object of classes `igcmatcatdist` or `igcmatcatdistn`, i.e., a
  positive or stochastic matrix (for the former, rows sum to one).

- weight:

  A vector of weights for the rows of `p`.

- base:

  A positive value, the base of the logarithm to be used to compute
  Jeffreys divergence.

- method:

  A character string specifying the method to be used for computing the
  approximation of Jeffreys centroid. Possible values are:

  - "mag" for the mean of arithmetic and normalized geometric centroids;

  - "njc" for the normalized Jeffreys centroid of positive vectors;

  - "jfr" for the Jeffreys-Fisher-Rao centroid, i.e., the Fisher-Rao
    mean of arithmetic and normalized geometric centroids;

  - "bsa" for a single-value binary search alogrithm based on Lagrangian
    optimization;

  - "gbic" for the Gauss-Bregman inductive centroid. Default is "gbic".
    See Details section.

- tol:

  For `igcmatcatdistn` objects only. A positive value, the tolerance in
  the approximation of Jeffreys centroid. Default is 1e-8.

- ...:

  Additionnal arguments for other methods.

## Details

Let \\p_1, \cdots, p_n\\ denote the rows of `p` and \\w_1, \cdots, w_n\\
the components of the vector `weight`. The output of
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md)
is a numerical approximation of the probability distribution \\p^\*\\
(object of class `igccatdistn`) satisfying \$\$p^\* \in \textrm{Argmin}
\left\\ \sum\_{i=1}^n w_i \mathbb{J}(p_i:q) \\ \| \\ q \>= 0 \sum\_{x}
q(x) =1 \right\\,\$\$ where \\\mathbb{J}((p_i:q)\\ denotes the Jeffreys
divergence between \\p_i\\ and \\q\\; see
[`jdiv()`](https://pregnault.github.io/igcentroids/reference/jdiv.md)
for details. Nielsen (2013) establishes using Lagrangian multipliers
method that \\p^\*\\ is of the form \$\$p^\*(x) = \frac{a(x)}{W( \exp(
\log(b) (\lambda + 1) a(x)/g(x)))},\$\$ where \\W\\ denotes the Lambert
function, \\a = \sum w_i p_i\\ is the arithmetic centroid of the vectors
\\p_i\\, \\g = \prod p_i^{w_i}\\ is their geometric centroid, \\b \>0\\
is the logarithm basis used to compute the Jeffreys divergence and
\\\lambda\\ is the unique real number such that \$\$\mathbb{K}(p^\* \|
g) = - \lambda.\$\$ There is no known closed-form expression for
\\\lambda\\.

The argument `method` identifies the method to be used to compute a
numerical approximation of \\\lambda\\ or directly of \\p^\*\\. More
precisely, if:

- `method = "bsa"`, a numerical approximation of \\\lambda\\ satisfying
  the last equation is computed thanks to a binary search. The argument
  `tol` refers to the numerical tolerance for this approximation. Then,
  the corresponding probability distribution \\p^\*\\ is computed. Due
  to the computational cost of evaluating the Lambert function, this
  method can be timec-consuming, particularly when the state space
  (i.e., the number of columns of `p`) is large;

- `method = "gbic"`, an iterative procedure called Gauss-Bregman
  inductive center is performed; see Nielsen (2024). It consists in
  iteratively computing arithmetic and geometric means of centroids
  until convergence. It is by far computationally less expansive than
  "bsa", with a similar quality of approximation. Its excellent
  performances make it the default method used by
  [`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md);

- `method = "mag"`, the mean \\(a+g)/2\\ of the arithmetic and geometric
  centroids is returned; it corresponds to the first iteration of
  "gbic";

- `method = "jfr"`, the Fisher-Rao midpoint of the arithmetic and
  geometric centroids is returned, i.e., the normalized version of
  \\(\sqrt{a} +\sqrt{g})^2\\;

- `method = "njc"`, the normalized version of the positive vector
  Jeffreys centroid is returned.

## References

Nielsen, F. (2013). Jeffreys Centroids: A Closed-From Expresion for
Positive Histograms and a Guaranteed Tight Approximation for Frequency
Histograms. *IEEE Signal Processing Letters*, **20**(7), 657–660.
[doi:10.1109/LSP.2013.2260538](https://ieeexplore.ieee.org/document/6509930)

Nielsen, F. (2024). Fast Proxy Centers for the Jeffreys Centroid: The
Jeffreys-Fisher-Rado Center and the Gauss-Bregman Inductive Center.
*Entropy*, **26**, 1008-1037.
[doi:10.3390/e26121008](https://www.mdpi.com/1099-4300/26/12/1008)

## See also

Generic function
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md),
[`jscent()`](https://pregnault.github.io/igcentroids/reference/jscent.md)
for computing Jensen-Shannon centroids

## Examples

``` r
u <- rep(1/3, 3)
prob <- c(0.9, 0.05, 0.05)
p <- matrix(c(u, prob),
             nrow = 2, ncol =3, byrow = TRUE)
p <- igcmatcatdistn(p)
meths <- c("mag", "njc", "jfr", "bsa", "gbic")
t(sapply(meths, function(m) jcent(p, method = m)))
#>           [,1]      [,2]      [,3]
#> mag  0.6481447 0.1759276 0.1759276
#> njc  0.6468378 0.1765811 0.1765811
#> jfr  0.6484677 0.1757661 0.1757661
#> bsa  0.6484927 0.1757537 0.1757537
#> gbic 0.6484681 0.1757659 0.1757659
```
