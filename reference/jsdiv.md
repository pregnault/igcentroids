# Jensen-Shannon divergence

`jsdiv()` is an S3 generic for computing the Jensen-Shannon divergence
between R objects, with methods for the following classes:

- `igccatdist`

- `igcmatcatdist`

## Usage

``` r
jsdiv(p, q, base = exp(1))
```

## Arguments

- p, q:

  Two **R** objects of same class. Currently, there are methods for
  `igccatdist` and `igcmatcatdist` objects.

- base:

  A positive value specifying the base of logarithm used for computing
  Jensen-Shannon divergence. Default is exp(1) (natural logarithm).

## Value

- If `p` and `q` are `igccatdist` objects, returns the value of the
  Jensen-Shannon divergence between `p` and `q`;

- if `p` and `q` are `igcmatcatdist` objects, returns a matrix of
  `nrow(p)` rows and `nrow(q)` columns whose coefficients are the
  Jensen-Shannon divergences between rows of `p` and rows of `q`.

## Details

The Jensen-Shannon divergence, also known as capacitory discrimination,
is a symmetric and bounded discrimination measure between finite-mass
positive measures defined on a measurable space. It is formally defined,
for \\p\\ and \\q\\ two finite-mass positive measures, by
\$\$\mathbb{JS}(p,q) = \frac{1}{2} \left( \mathbb{K}\left(p:
\frac{p+q}{2} \right) + \mathbb{K}\left(q: \frac{p+q}{2}
\right)\right),\$\$ where \\\mathbb{K}\\ denotes the Kullback-Leibler
divergence; see
[`kldiv()`](https://pregnault.github.io/igcentroids/reference/kldiv.md).

For probability distributions (pds), it is linked to Shannon entropy
\\\mathbb{S}\\ through the relation \$\$\mathbb{JS}(p,q) = \mathbb{S}
\left( \frac{p+q}{2} \right) - \frac{\mathbb{S}(p) +
\mathbb{S}(q)}{2}.\$\$ Properties of Jensen-Shannon divergence include:

- definition for pds with distinct supports: \\\mathbb{JS}(p,q)\\ is
  well-defined wether or not \\p\\ and \\q\\ have same support,

- symmetry: \\\mathbb{JS}(p,q) = \mathbb{JS}(q,p)\\,

- boundedness: \\\mathbb{JS}(p,q) \leq \log(2)\\, with equality whenever
  supports of \\p\\ and \\q\\ are disjoint.

## See also

[`kldiv()`](https://pregnault.github.io/igcentroids/reference/kldiv.md),
[`jdiv()`](https://pregnault.github.io/igcentroids/reference/jdiv.md),
[catdist](https://pregnault.github.io/igcentroids/reference/catdist.md).

## Examples

``` r
# Compute the Jensen-Shannon (JS) divergence
# between two probability distributions with common support.
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
jsdiv(p, q)
#> [1] 0.05418392
# Compute the (squared) matrix of JS divergences between
# the elements of set of probabilily distributions
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
matptest <- igcmatcatdistn(p = matp)
jsdiv(p = matptest)
#>            [,1]       [,2]       [,3]
#> [1,] 0.00000000 0.05418392 0.11222490
#> [2,] 0.05418392 0.00000000 0.04750864
#> [3,] 0.11222490 0.04750864 0.00000000
# Compute the matrix of JS divergences between
# two sets of probability distributions
matp2 <- matrix(c(dbinom(0:4, size = 5, prob = 0.3),
                      c(0.1, 0.1, 0.2, 0.2, 0.4)),
                    nrow = 2, byrow = TRUE)
matptest2 <- igcmatcatdist(matp2)
jsdiv(p = matptest, q = matptest2)
#>            [,1]       [,2]
#> [1,] 0.03165236 0.10241485
#> [2,] 0.05793705 0.03397981
#> [3,] 0.04961128 0.13952469
```
