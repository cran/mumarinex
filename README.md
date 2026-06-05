
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mumarinex

<!-- badges: start -->

<img src="man/figures/logo.png" align="right" height="120" alt="" />
<!-- badges: end -->

The ‘mumarinex’ package aims to facilitate the computation of the
multivariate marine recovery index, providing functions for data
visualization and ecological diagnostics of marine ecosystems. The
computational details are described in the original publication (please,
see below for reference).

## Installation

You can install the development version of mumarinex from CRAN or from
[GitHub](https://github.com/), respectively with:

``` r
install.packages("mumarinex") # CRAN
devtools::install_github("Nathan-Chauvel/mumarinex") # GitHub
```

## Functions of the mumarinex package

| Name | Short description |
|----|----|
| `decomplot()` | Generates a graphical representation of MUMARINEX sub-indices. |
| `diagnostic_tool()` | Identifies, for each sub-index, the species/taxa that most contribute to its variation.. |
| `mumarinex()` | Computes the MUMARINEX index and its sub-indices (SCSR, CBCS, and SPI). |
| `mum_GUI()` | Launch the online Graphical User Interface to compute the MUMARINEX index and its sub-indices. |

## Reference

Chauvel, N., Grall, J., Thiébaut, E., Houbin, C., Pezy, J.-P., 2026. A
general-purpose multivariate marine recovery index (MUMARINEX) for
quantifying the influence of human activities on benthic habitat
ecological status. Ecological Indicators 188, 115002.
<https://doi.org/10.1016/j.ecolind.2026.115002>
