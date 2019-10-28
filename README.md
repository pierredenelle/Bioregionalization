README
================

[![Travis build
status](https://travis-ci.org/pierredenelle/Bioregionalization.svg?branch=master)](https://travis-ci.org/pierredenelle/Bioregionalization)

`Bioregionalization` is a package designed to compare several methods of
clustering based on multivariate or network
approaches.

## 1 Installation

``` r
devtools::install_github("pierredenelle/Bioregionalization", build_vignettes = TRUE)
```

## 2 Tutorial

A tutorial vignette showing the main steps of the processing can be
visualised with the following command line:

``` r
vignette("tutorial", package = "Bioregionalization")
```

    ## Warning: vignette 'tutorial' not found

## 3 Dependencies

`Bioregionalization` depends on `dplyr`, `ecodist`, `reshape2`,
`DescTools`, `ade4`, `cluster`, `sf`, `ggplot2`, `Rcpp`, `SMUT` and
`igraph`.
