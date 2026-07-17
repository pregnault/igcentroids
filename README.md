
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `igcentroids`: An R package for computing information-geometric centroids of probability distributions

<!-- badges: start -->

![](https://img.shields.io/badge/lifecycle-maturing-blue.svg) [![Project
Status: Active – The project has reached a stable, usable state and is
being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/pregnault/igcentroids/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pregnault/igcentroids/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`{igcentroids}` provides an implementation for R statistical software of
several methods for computing exact or approximated centroids of
probability distributions based on Jeffreys and Jensen-Shannon
divergences.

## Installation

You can install the development version of `{igcentroids}` like so:

``` r
# Install package remotes if not already installed
if (! "remotes" %in% installed.packages()[,"Package"]) {
  install.packages("remotes")
}
# Install package samplelim from its GitHub repo
remotes::install_github("https://github.com/pregnault/igcentroids")
```

## Purpose

`{igcentroids}` provides information-geometric tools for clustering
probability distributions. More precisely, it introduces S3 classes and
methods for:

- computing Kullback-Leibler, Jeffreys and Jensen-Shannon divergences of
  :

  - positive vectors;
  - probability distributions (pds) with finite support;

- computing centroids of a set of such objects, i.e., an object of the
  same mathematical nature minimizing the (weighted) mean of the
  divergences ;

- integrating these computations in a clustering workflow performed with
  the package `{flexclust}`.

## How to use

Here is a toy example, illustrating a clustering workflow of a random
set of three-state pds using Jeffreys divergence and related centroids.

We first generate a random set of 32 two-state pds.

``` r
library(igcentroids) # Load package
# Function for randomly generating points in the 2D-simplex
runifsimplex <- function(n = 1L) {
  x <- runif(n)
  y <- runif(n, min = 0, max = 1-x)
  res <- igcmatcatdistn(matrix(c(x,y, 1-x-y), ncol = 3, byrow = FALSE))
  return(res)
}
set.seed(527) # Set seed of PNRG
samp <- runifsimplex(n = 32L) # Generate 32 two-state prob. dist.
sampsimp <- igcmatcatdistn(samp) # Convert into class igcmatcatdist
```

We can visualize these pds as points of the 2D-simplex, identified to
the equilateral triangle.

``` r
plot(sampsimp)
```

![](README_files/figure-gfm/vizsimplex-1.png)<!-- -->

The matrix of the Jeffreys divergences between these pds is obtained as
follows.

``` r
matjdiv <- jdiv(sampsimp)
# Display only the five first rows and columns
matjdiv[1:5, 1:5]
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]
    ## [1,] 0.0000000 0.4298333 0.6262441 0.4615052 1.6558647
    ## [2,] 0.4298333 0.0000000 0.3150661 1.2050240 0.9587228
    ## [3,] 0.6262441 0.3150661 0.0000000 2.1362794 0.1565056
    ## [4,] 0.4615052 1.2050240 2.1362794 0.0000000 4.2208317
    ## [5,] 1.6558647 0.9587228 0.1565056 4.2208317 0.0000000

Then, the Jeffreys centroid of these pds.

``` r
jcent(sampsimp)
```

    ## [1] 0.5068580 0.2079931 0.2851488
    ## attr(,"class")
    ## [1] "igccatdistn" "igccatdist" 
    ## attr(,"mass")
    ## [1] 1

`{igcentroids}` provides functions for an easy integration of these
methods for clustering tasks using the package `{flexclust}`, as
illustrated below.

``` r
library("flexclust")
# Set jclust family: clustering based on Jeffreys divergence and centroids
jfam <- kccaFamily(name = "jclust",
                   dist = jdivflex,
                   cent = jcentflex,
                   preproc = undo_igc)
# Operate clustering
res_clust <- kcca(x = sampsimp, k= 4, family = jfam)
# Visualize clusters
plot(sampsimp, col = res_clust@cluster+1)
```

![](README_files/figure-gfm/clusteringintegration-1.png)<!-- -->

We can compute the Jeffreys divergences between the pds and the centroid
of their clusters.

``` r
cents <- igcmatcatdist(res_clust@centers) # Centroids
clusters <- split(1:nrow(sampsimp), res_clust@cluster)
jdiv2cents <- list(4)
for (i in 1:4) {
  jdiv2cents[[i]] <- jdiv(p = igcmatcatdistn(sampsimp[clusters[[i]],]), 
                          q = cents)
}
jdiv2cents
```

    ## [[1]]
    ##             [,1]      [,2]      [,3]      [,4]
    ##  [1,] 0.04179681 1.4627622 1.1562806 1.1215317
    ##  [2,] 0.13513820 1.7277522 2.4780415 1.3170095
    ##  [3,] 0.53523277 1.4404718 4.7757247 2.6307397
    ##  [4,] 0.18807317 1.5849733 2.8979431 1.6023330
    ##  [5,] 0.05086308 0.9447104 1.6021936 1.6647139
    ##  [6,] 0.15196904 0.5591045 2.2893058 2.2933588
    ##  [7,] 0.10399919 1.3179992 1.0700267 1.3256607
    ##  [8,] 0.08402536 0.8653833 1.5972523 1.7847049
    ##  [9,] 0.25892925 1.9605293 0.6047796 0.9368975
    ## [10,] 0.33862482 1.2621630 0.9428206 1.7266241
    ## [11,] 0.05337381 1.2726543 2.2856323 1.5286908
    ## [12,] 0.17427924 2.0451277 0.7215748 0.7551340
    ## [13,] 0.06298114 1.8954377 1.7132045 0.9090367
    ## [14,] 0.34775212 2.9823803 1.6943987 0.4620123
    ## 
    ## [[2]]
    ##           [,1]        [,2]     [,3]     [,4]
    ## [1,] 1.1762538 0.007547095 4.521785 4.870659
    ## [2,] 1.2544066 0.033535451 5.298650 5.113569
    ## [3,] 1.8335925 0.042377553 4.412256 5.916831
    ## [4,] 1.7485681 0.019194765 5.303618 5.951605
    ## [5,] 0.4653694 0.233081325 3.387479 3.293830
    ## [6,] 1.3210977 0.031782600 3.763081 4.955200
    ## [7,] 1.7726012 0.020084774 5.160121 5.966074
    ## [8,] 2.1064831 0.050871191 5.169548 6.498406
    ## 
    ## [[3]]
    ##          [,1]     [,2]         [,3]      [,4]
    ## [1,] 1.576205 5.116070 0.1514593751 0.3379225
    ## [2,] 1.635318 4.452160 0.0009273484 0.9430653
    ## [3,] 4.246695 7.499680 0.6236321861 2.7844979
    ## [4,] 1.253974 3.698391 0.0372771886 0.9635677
    ## [5,] 0.797882 3.081353 0.1464568790 0.7678379
    ## [6,] 1.196562 4.011474 0.0398383290 0.6123780
    ## 
    ## [[4]]
    ##           [,1]     [,2]      [,3]        [,4]
    ## [1,] 0.5119198 3.457708 0.8365206 0.151698243
    ## [2,] 1.2399195 5.115600 1.0657271 0.009254796
    ## [3,] 3.4750676 9.046286 0.8596527 0.344943580
    ## [4,] 0.9565685 4.545494 1.3616577 0.085255848

Learn more about the possibilities offered by igcentroids by reading the
vignettes of the package.

## Credits

The algorithms have been designed by Frank NIELSEN; see Nielsen (2013),
Nieslen (2020) and Nielsen (2024).

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

Nielsen, F. (2020). On a Generalization of the Jensen-Shannon Divergence
and the Jensen-Shannon Centroid. *Entropy*, **22**, 221–244. DOI:
<https://www.mdpi.com/1099-4300/22/2/221>

Nielsen, F. (2024). Fast Proxy Centers for the Jeffreys Centroid: The
Jeffreys-Fisher-Rao Center and the Gauss-Bregman Inductive Center.
*Entropy*, **26**, 1008-1037. DOI: <https://doi.org/10.3390/e26121008>
