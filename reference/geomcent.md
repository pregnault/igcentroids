# Geometric centroid of n objects

`geomcent()` is an S3 generic for computing the geometric mean of R
objects, with methods for the following classes:

- `igcmatcatdist`

## Usage

``` r
geomcent(p, weight = NULL, base = exp(1))
```

## Arguments

- p:

  An **R** object gathering a number `n` of objects of same nature.
  Currently, there is a method for `igcmatcatdist` matrix, i.e., a
  matrix of categorical distribution.

- weight:

  A positive vector of length `n`. The geometric mean is balanced by
  weighting each row of `p` with the corresponding element of `weight`.

- base:

  A positive value specifying the base of logarithm used for computing
  Jeffreys divergence. Default is exp(1) (natural logarithm).

## Value

The geometric centroid of the input objects `p` weighted by `weight`.

## See also

[`arithcent()`](https://pregnault.github.io/igcentroids/reference/arithcent.md),
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md),
[`jscent()`](https://pregnault.github.io/igcentroids/reference/jscent.md),
[catdist](https://pregnault.github.io/igcentroids/reference/catdist.md).

## Examples

``` r
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
matp <- matrix(c(p, q), nrow = 2)
matptest <- igcmatcatdist(p = matp)
geomcent(p = matptest)
#> [1] 0.1250000 0.3061862 0.1118034 0.2000000 0.2000000
#> attr(,"class")
#> [1] "igccatdist"
#> attr(,"mass")
#> [1] 0.9429896
```
