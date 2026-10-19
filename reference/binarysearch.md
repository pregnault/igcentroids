# Binary search algorithm

An **R** implementation of the binary search of the zero of a continous
strictly monotone real function.

## Usage

``` r
binarysearch(lb, ub, fun, tol = 1e-08)
```

## Arguments

- lb:

  A single numeric value, the lower bound of the interval to seek in the
  zero point.

- ub:

  A single numeric value, the upper bound of the interval to seek in the
  zero point.

- fun:

  A continuous monotone real function defined on the interval `[lb, ub]`
  such that `fun(lb)*fun(ub) < 0`.

- tol:

  The precision of the approximation of the zero point to be attained.
  Default value is `1e-8`.

## Value

A list with three components :

- `sol`, the approximation of the zero point for function `fun` ;

- `prec`, the precision of the approximation, of the form
  `(ub - lb)/2^(iter+1)` ;

- `iter`, the number of iterations of the binary search while processing
  the approximation.

## Examples

``` r
tmpfun <- function(x) return(sqrt(x)-2)
tmpub <- 9
tmplb <- 1
binarysearch(lb = tmplb, ub = tmpub, fun = tmpfun)
#> $sol
#> [1] 4
#> 
#> $prec
#> [1] 7.450581e-09
#> 
#> $iter
#> [1] 29
#> 
```
