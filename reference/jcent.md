# Jeffreys centroid of positive vectors and categorical probability distributions.

`jcent()` is an S3 generic for computing the Jeffreys centroid of R
objects, with methods for the following classes:

- `igcmatcatdist`

- [`igcmatcatdistn`](https://pregnault.github.io/igcentroids/reference/jcent.igcmatcatdistn.md)

## Usage

``` r
jcent(p, weight = NULL, base = exp(1), ...)
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

- ...:

  Additionnal arguments for other methods.

## Value

A positive vector (object of class `igccatdist`) if `p` is of class
`igcmatcatdist` or a categorical probability distribution (object of
class `igccatdistn`), approximation of the Jeffreys centroid of the rows
of `p`.

## Details

If input `p` is a positive matrix (object of class `igcmatcatdist`),
jcent()`returns the positive vector minimizing the weighted mean of Jeffreys divergences with the rows of`p`. More precisely, by denoting \eqn{p_1, \cdots, p_n} the rows of `p`and \eqn{w_1, \cdots, w_n} the components of the vector`weight`, the output of jcent()`
is the positive vector \\p^\*\\ satisfying \$\$p^\* \in \textrm{Argmin}
\left\\ \sum\_{i=1}^n w_i \mathbb{J}(p_i:q) \\ \| \\ q \textrm{ pos.
vector} \right\\,\$\$ where \\\mathbb{J}((p_i:q)\\ denotes the Jeffreys
divergence between \\p_i\\ and \$\$q\$\$; see
[`jdiv()`](https://pregnault.github.io/igcentroids/reference/jdiv.md)
for details. The minimizer \\p^\*\\ has a closed-from expression given
by \$\$p^\*(x) = \frac{a(x)}{W( \exp( \log(b) a(x)/g(x)))},\$\$ where
\\W\\ denotes the Lambert function, \\a = \sum w_i p_i\\ is the
arithmetic centroid of the vectors \\p_i\\, \\g = \prod p_i^{w_i}\\ is
their geometric centroid and \\b \>0\\ is the logarithm basis used to
compute the Jeffreys divergence.

If input `p` is a stochastic matrix (object of class `igcmatcatdistn`),
the output is the probability distribution minimizing the weighted mean
of Jeffreys divergences with the rows of `p`; formally, \$\$p^\* \in
\textrm{Argmin} \left\\ \sum\_{i=1}^n w_i \mathbb{J}(p_i:q) \| q \>= 0
\sum\_{x} q(x) =1 \right\\.\$\$ No closed-form expression exists for
\\p^\*\\ is this case. Several methods are implemented so as to compute
numerical approximations; see Details section in
[`jcent.igcmatcatdistn()`](https://pregnault.github.io/igcentroids/reference/jcent.igcmatcatdistn.md).

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

[`arithcent()`](https://pregnault.github.io/igcentroids/reference/arithcent.md),
[`geomcent()`](https://pregnault.github.io/igcentroids/reference/geomcent.md),
[`jscent()`](https://pregnault.github.io/igcentroids/reference/jscent.md),
[`lamW::lambertW0()`](https://rdrr.io/pkg/lamW/man/lamW.html),
[`jcent.igcmatcatdistn()`](https://pregnault.github.io/igcentroids/reference/jcent.igcmatcatdistn.md)

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
