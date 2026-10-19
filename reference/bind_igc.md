# Bind positive matrices and vectors

`bind_igc()` behaves for positive matrices and vectors as
[`base::rbind()`](https://rdrr.io/r/base/cbind.html) does for matrices:
it binds the matrices if the numbers of columns are the same for all
arguments. If all arugments are stochastic matrices or categorical
distributions, the result is a stochastic matrix. Otherwise, it returns
a positive matrix.

## Usage

``` r
bind_igc(...)
```

## Arguments

- ...:

  The positive matrices and vectors to be binded.

## Examples

``` r
# Bind positive matrices
A <- igcmatcatdist(matrix(1:6, nrow = 2))
B <- igcmatcatdistn(matrix(c(dbinom(0:2, size = 2, prob = 0.2),
dbinom(0:2, size = 2, prob = 0.5)),
nrow = 2, byrow = TRUE))
bind_igc(A,B)
#>      [,1] [,2] [,3]
#> [1,] 1.00 3.00 5.00
#> [2,] 2.00 4.00 6.00
#> [3,] 0.64 0.32 0.04
#> [4,] 0.25 0.50 0.25
#> attr(,"class")
#> [1] "igcmatcatdist"
#> attr(,"mass")
#> [1]  9 12  1  1
# Bind positive matrices and vectors
C <- igccatdist(7:9)
D <- igccatdistn(dbinom(0:2, size = 2, prob = 0.7))
bind_igc(A, B, C, D)
#>      [,1] [,2] [,3]
#> [1,] 1.00 3.00 5.00
#> [2,] 2.00 4.00 6.00
#> [3,] 0.64 0.32 0.04
#> [4,] 0.25 0.50 0.25
#> [5,] 7.00 8.00 9.00
#> [6,] 0.09 0.42 0.49
#> attr(,"class")
#> [1] "igcmatcatdist"
#> attr(,"mass")
#> [1]  9 12  1  1 24  1
# Returns a stochastic matrix when arguments are stochastic matrices or categorical distributions
bind_igc(B, D)
#>      [,1] [,2] [,3]
#> [1,] 0.64 0.32 0.04
#> [2,] 0.25 0.50 0.25
#> [3,] 0.09 0.42 0.49
#> attr(,"class")
#> [1] "igcmatcatdistn" "igcmatcatdist" 
#> attr(,"mass")
#> [1] 1 1 1
```
