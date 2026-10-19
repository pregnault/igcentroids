# Kullback-Leibler Divergence

Compute the Kullback-Leibler divergence of a positive measure with
respect to another dominating the first one.

## Usage

``` r
kldiv(p, q, base = exp(1))
```

## Arguments

- p, q:

  Two **R** objects of same class. Currently, there is a method for
  `igccatdist` objects.

- base:

  A positive value, the base in which logarithms are computed.

## Value

The value of the Kullback-Leibler divergence from `p` to `q`.

## Details

For any finite-mass positive measures \\p\\ and \\q\\ on a measurable
space, such that \\p\\ is dominated by \\q\\ and both are dominated by
some measure \\\mu\\, the Kullback-Leibler divergence
\\\mathbb{K}(p:q)\\ from \\p\\ to \\q\\ is given by \$\$\mathbb{K}(p:q)
= \int p \log(p/q) \textrm{d} \mu - \int p \textrm{d} \mu + \int q
\textrm{d} \mu.\$\$

If \\p\\ and \\q\\ are probability distributions on a finite set
\\\mathcal{X}\\, the Kullback-Leibler divergence from \\p\\ to \\q\\
reduces to \$\$\mathbb{K}(p:q) = \sum\_{x \in \mathcal{X}} p(x) \log
\frac{p(x)}{q(x)},\$\$ with the conventions \\0 \log \frac{0}{0} = 0\\,
\\0 \log \frac{0}{q} = 0\\ and \\p \log \frac{p}{0} = \infty\\.

## Examples

``` r
p <- igccatdist(dbinom(0:4, size = 4, prob = 0.5))
q <- igccatdist(rep(1/5, times = 5))
kldiv(p, q)
#> [1] 0.2019062
```
