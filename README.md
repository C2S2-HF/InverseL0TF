
# A Splicing Approach to $\ell_0$ Trend Filtering Using Inverse Transformation

<!-- badges: start -->
<!-- badges: end -->

## Introduction

This package provides an efficient solution for $\ell_0$ Trend Filtering, avoiding the traditional methods of using Lagrange duality or ADMM algorithms. It employ a splicing approach that minimizes L0-regularized sparse approximation by transforming the $\ell_0$ Trend Filtering problem. 

## R Package Installation

`L0TFinv` can be installed from Github as follows:

``` r
library(devtools)
install_github("C2S2-HF/InverseL0TF", repos = NULL, type = "source")
```

Alternatively, you can run the following code in R to install `L0TFinv` after downloading [L0TFinv_1.0.1.tar.gz](./L0TFinv_1.0.1.tar.gz).

```r
install.packages("Your_download_path/L0TFinv_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Usage

For a tutorial, please refer to [L0TFinv's Vignette](./vignettes). 



