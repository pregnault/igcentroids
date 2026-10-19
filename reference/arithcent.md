# Arithmetic centroid of n objects

`arithcent()` is an S3 generic for computing the arithmetic mean of R
objects, with methods for the following classes:

- `igcmatcatdist`

## Usage

``` r
arithcent(p, weight = NULL)
```

## Arguments

- p:

  An **R** object gathering a number `n` of objects of same nature.
  Currently, there is a method for `igcmatcatdist` matrix, i.e., a
  matrix of categorical distribution; see
  [catdist](https://pregnault.github.io/igcentroids/reference/catdist.md).

- weight:

  A positive vector of length `n`. The arithmetic mean is balanced by
  weighting each row of `p` with the corresponding element of `weight`.

## Value

The arithmetic centroid of the input objects `p` weighted by `weight`.

## See also

[`geomcent()`](https://pregnault.github.io/igcentroids/reference/geomcent.md),
[`jcent()`](https://pregnault.github.io/igcentroids/reference/jcent.md),
[`jscent()`](https://pregnault.github.io/igcentroids/reference/jscent.md),
[catdist](https://pregnault.github.io/igcentroids/reference/catdist.md).

## Examples

``` r
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
matp <- matrix(c(p, q), nrow = 2)
matptest <- igcmatcatdist(p = matp)
arithcent(p = matptest)
#> [1] 0.15625 0.31250 0.13125 0.20000 0.20000
#> attr(,"class")
#> [1] "igccatdist"
#> attr(,"mass")
#> [1] 1
```
