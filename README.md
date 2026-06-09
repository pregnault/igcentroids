
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `igcentroids`: An R package for computing Informational Geometrical centroids of postive vectors and categorical probability distributions

<!-- badges: start -->

![](https://img.shields.io/badge/lifecycle-maturing-blue.svg) [![Project
Status: Active – The project has reached a stable, usable state and is
being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/pregnault/igcentroids/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pregnault/igcentroids/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`{igcentroids}` provides an implementation for R statistical software of
several methods for computing exact or approximated centroids of
positive vectors and probability distributions with finite support.

## Installation

You can install the development version of `{igcentroids}` like so:

``` r
# Install package remotes if not already installed
if (! "remotes" %in% installed.packages()[,"Package"]) {install.packages("remotes")}
# Install package samplelim from its GitHub repo
remotes::install_github("https://github.com/pregnault/igcentroids")
```

## Examples

Here is a first example, illustrating the five methods implemented for
computing approximated centroids of probability distributions.

``` r
library(igcentroids)
## basic example code
u <- rep(1/3, 3)
prob <- c(0.9, 0.05, 0.05)
p <- matrix(c(u, prob),
           nrow = 2, ncol =3, byrow = TRUE)
p <- igcmatcatdistn(p)
meths <- c("mag", "njc", "jfr", "bsa", "gbic")
t(sapply(meths, function(m) jcent(p, method = m)))
```

## Credits

The algorithms have been designed by Frank NIELSEN; see Nielsen (2013)
and Nielsen (2024).

The core implementation of those algorithms and R packaging have been
performed by Philippe REGNAULT.

## Licensing

You may redistribute or modify the software under the GNU General Public
Licence as published by Free Software Foundation, either version 3 of
the License, or any later version. It is distributed in the hope that it
may be helpful to the community, but WITHOUT ANY WARRANTY.

## References

Nielsen, F. (2013). Jeffreys Centroids: A Closed-From Expresion for
Positive Histograms and a Guaranteed Tight Approximation for Frequency
Histograms. *IEEE Signal Processing Letters*, **20**(7), 657–660. DOI:
<https://doi.org/10.1109/LSP.2013.2260538>

Nielsen, F. (2024). Fast Proxy Centers for the Jeffreys Centroid: The
Jeffreys-Fisher-Rado Center and the Gauss-Bregman Inductive Center.
*Entropy*, **26**, 1008-1037. DOI: <https://doi.org/10.3390/e26121008>
