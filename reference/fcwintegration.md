# Wrapped-up functions for seemless-integration in `{flexclust}` workflow.

Functions `preproc_igc()`, `jdivflex()`, `jcentlex()`, `jsdivflex()`and
`jscentlex()` are intended to ease `igcentroids` integration into
clustering workflow with `{flexclust}`.

## Usage

``` r
undo_igc(p)

jdivflex(x, centers, base = exp(1))

jcentflex(x, weight = NULL, method = "gbic", base = exp(1))

jsdivflex(x, centers, base = exp(1))

jscentflex(x, weight = NULL, method = "gbic", base = exp(1))
```

## Arguments

- p:

  An object of class `igmatcatdistn`.

- x:

  A positive matrix.

- centers:

  A positive matrix whose rows are the (temporary) centers of clusters.

- base:

  A positive value, the base in which logarithms are computed.

- weight:

  A positive vector of length `n`. The arithmetic mean is balanced by
  weighting each row of `p` with the corresponding element of `weight`.

- method:

  A character string specifying the method to be used for computing the
  approximation of Jeffreys centroid. See
  [`jcent.igcmatcatdistn()`](https://pregnault.github.io/igcentroids/reference/jcent.igcmatcatdistn.md)
  for details.

## Value

Objects suitable for defining a clustering family with
[`flexclust::kccaFamily()`](https://rdrr.io/pkg/flexclust/man/kcca.html).

## See also

[`flexclust::kccaFamily()`](https://rdrr.io/pkg/flexclust/man/kcca.html)

## Examples

``` r
# Create an object of class igcmatcatdistn
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
r <- igccatdist(c(0.3, 0.4, 0.1, 0.1, 0.1))
matp <- rbind(p, q, r)
matptest <- igcmatcatdistn(p = matp)
# Avoid class so as to return an object of type matrix
undo_igc(matptest)
#>     [,1] [,2]  [,3] [,4]   [,5]
#> p 0.0625 0.25 0.375 0.25 0.0625
#> q 0.2000 0.20 0.200 0.20 0.2000
#> r 0.3000 0.40 0.100 0.10 0.1000
# Compute Jeffreys divergences from a matrix as returned by undo_igc()
center <- matrix(rep(1/5, times = 5), nrow = 1)
jdivflex(x = matp, centers = center)
#>           [,1] [,2]      [,3]
#> [1,] 0.4521873    0 0.3871201
# Compute Jeffreys centroid from a matrix as returned by undo_igc()
jcentflex(x = matp)
#> [1] 0.1798936 0.2922401 0.2211160 0.1865516 0.1201987
#> attr(,"class")
#> [1] "igccatdistn" "igccatdist" 
#> attr(,"mass")
#> [1] 1
# Compute Jensen-Shannon divergences from a matrix as returned by undo_igc()
jsdivflex(x = matp, centers = center)
#>            [,1] [,2]       [,3]
#> [1,] 0.05418392    0 0.04750864
# Compute Jensen-Shannon centroid from a matrix as returned by undo_igc()
jscentflex(x = matp)
#> [1] 0.1821431 0.2907138 0.2210574 0.1863583 0.1197274
#> attr(,"class")
#> [1] "igccatdistn" "igccatdist" 
#> attr(,"mass")
#> [1] 1
```
