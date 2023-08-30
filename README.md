
# gg4way

<!-- badges: start -->

<!-- badges: end -->

4way plots enable a comparison of two contrasts of differential gene expression. The gg4way package creates 4way plots using the ggplot2 framework and supports popular Bioconductor objects. The package also provides information about the correlation and genes of interest.

## Installation

gg4way can be installed from [Bioconductor](https://bioconductor.org/packages/gg4way):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
  }

BiocManager::install("gg4way")
```

To install the development version directly from GitHub:

```r
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
  }

remotes::install_github("ben-laufer/gg4way")
```

## Examples

See the package vignette.

