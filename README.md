
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JM.landslide: Joint Modeling of Landslide Counts and Sizes

<!-- badges: start -->
<!-- badges: end -->

The goal of `JM.landslide` package is to demonstrate how to jointly
model landslide counts and sizes using a sub-asymptotic distribution for
the landslide sizes, following the joint modeling approach of [Yadav et
al. (2024)](https://arxiv.org/abs/2404.09156) and [Yadav et
al. (2023)](https://doi.org/10.1093/jrsssc/qlad077). The package
performs simulation-based MCMC inference to simultaneously carry out
inference and prediction of landslide counts and sizes. Additionally,
the package provides estimates of landslide hazards and susceptibility
maps over the study regions, key components often sought in landslide
literature. Moreover, the package also generates outputs that provide
estimates and summaries of fixed and random effects, along with the
fitted counts and sizes.

## Installation

You may install the development version of package using:

``` r
devtools::install_github("yadavrishikesh/JM.landslide")
```

## Vignette

See the
[Vignette](https://yadavrishikesh.github.io/JM.landslide/JM.landslide.html)
for detailed description of the package.

## References

1.  Yadav, R., Lombardo, L., & Huser, R. (2024). *Statistics of extremes
    for natural hazards: landslides and earthquakes*. arXiv preprint
    [arXiv:2404.09156](https://arxiv.org/abs/2404.09156).
2.  Yadav, R., Huser, R., Opitz, T., & Lombardo, L. (2024) Rishikesh
    Yadav, Raphael Huser, Thomas Opitz, and Luigi Lombardo (2023).
    *Joint modeling of landslide counts and sizes using spatial marked
    point processes with sub-asymptotic mark distributions.*
    `Journal of the Royal Statistical Society Series C (JRSSC): Applied Statistics, qlad077`.
    <https://doi.org/10.1093/jrsssc/qlad077>
