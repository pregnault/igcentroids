# Classes igccatdist, igccatdistn, igcmatcatdist, igcmatcatdistn for positive (normalized) vectors and matrices.

`igccatdist`, `igccatdistn`, `igcmatcatdist`, `igcmatcatdistn` are S3
classes for, respectively, positive vectors, probability distributions,
positive matrices and stochastic matrices (positive matrices whose rows
sum to one).

## Usage

``` r
igccatdist(p = numeric())

igccatdistn(p = numeric())

igcmatcatdist(p = numeric())

igcmatcatdistn(p = numeric())
```

## Arguments

- p:

  A numeric vector (for `igcatdist`, `igcatdistn`) or matrix
  (`igmatcatdist`, `igmatcatdistn`) of positive values.

## Value

An object of the eponyme class.

## Examples

``` r
igccatdist(c(0.3, 0.1, 0.4, 0.8))
#> [1] 0.3 0.1 0.4 0.8
#> attr(,"class")
#> [1] "igccatdist"
#> attr(,"mass")
#> [1] 1.6
p <- igccatdistn(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdistn(rep(1/5, times = 5))
p+q
#> [1] 0.2625 0.4500 0.5750 0.4500 0.2625
#> attr(,"class")
#> [1] "igccatdistn" "igccatdist" 
#> attr(,"mass")
#> [1] 1
igcmatcatdist(matrix(c(dbinom(0:4, size = 4, prob = 0.5), 1:5),
                     nrow = 2, byrow = TRUE))
#>        [,1] [,2]  [,3] [,4]   [,5]
#> [1,] 0.0625 0.25 0.375 0.25 0.0625
#> [2,] 1.0000 2.00 3.000 4.00 5.0000
#> attr(,"class")
#> [1] "igcmatcatdist"
#> attr(,"mass")
#> [1]  1 15
matp <- matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                 rep(1/5, times = 5)),
               nrow = 2, byrow = TRUE)
matptest <- igcmatcatdistn(p = matp)
```
