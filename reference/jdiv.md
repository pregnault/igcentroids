# Jeffreys Divergence

`jdiv()` is an S3 generic for computing the Jeffreys divergence – or
symmetrized Kullback-Leibler divergence – between R objects, with
methods for the following classes:

- `igccatdist`

- `igcmatcatdist`

## Usage

``` r
jdiv(p, q, base = exp(1))
```

## Arguments

- p, q:

  Two **R** objects of same class. Currently, there are methods for
  `igccatdist` and `igcmatcatdist` objects.

- base:

  A positive value specifying the base of logarithm used for computing
  Jeffreys divergence. Default is exp(1) (natural logarithm).

## Value

- If `p` and `q` are `igccatdist` objects, returns the value of the
  Jeffreys divergence (or symmetrized Kullback-Leibler) between `p` and
  `q`;

- if `p` and `q` are `igcmatcatdist` objects, returns a matrix of
  `nrow(p)` rows and `nrow(q)` columns whose coefficients are the
  Jeffreys divergence from rows of `p` to rows of `q`.

## Details

The Jeffreys divergence or symmetric Kullback-Leibler divergence is a
symmetric and bounded discrimination measure between positive measures
defined on a measurable space. It is formally defined, for \\p\\ and
\\q\\ two positive measures with same support, by \$\$\mathbb{J}(p,q) =
\mathbb{K}\left(p: q \right) + \mathbb{K}\left(q:p \right),\$\$ where
\\\mathbb{K}\\ denotes the Kullback-Leibler divergence; see
[`kldiv()`](https://pregnault.github.io/igcentroids/reference/kldiv.md).

## See also

[`kldiv()`](https://pregnault.github.io/igcentroids/reference/kldiv.md),
[`jsdiv()`](https://pregnault.github.io/igcentroids/reference/jsdiv.md),
[catdist](https://pregnault.github.io/igcentroids/reference/catdist.md)

## Examples

``` r
# Compute the Jeffreys divergence
# between two probability distributions with common support.
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
jdiv(p, q)
#> [1] 0.4521873
# Compute the (squared) matrix of Jeffreys divergences between
# the elements of set of probabilily distributions
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
matp <- matrix(c(p, q, r), nrow = 3, byrow = TRUE)
matptest <- igcmatcatdistn(p = matp)
jdiv(p = matptest)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.4521873 0.9615984
#> [2,] 0.4521873 0.0000000 0.3871201
#> [3,] 0.9615984 0.3871201 0.0000000
# Compute the matrix of Jeffreys divergences
# between two sets of probability distributions
matp2 <- matrix(c(dbinom(0:4, size = 5, prob = 0.3),
                      c(0.1, 0.1, 0.2, 0.2, 0.4)),
                    nrow = 2, byrow = TRUE)
matptest2 <- igcmatcatdist(matp2)
jdiv(p = matptest, q = matptest2)
#>           [,1]      [,2]
#> [1,] 0.2594416 0.9027330
#> [2,] 0.5102637 0.2772589
#> [3,] 0.4152278 1.1901285
```
