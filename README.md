
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPCWJK

<!-- badges: start -->

<!-- badges: end -->

The goal of IPCWJK is to …

## Installation

You can install the development version of IPCWJK like so:

``` r
# With remotes
remotes::install_github("IDEN-Project-UAS-Darmstadt/IPCWJK")
# With pak (recommended for speed)
pak::pkg_install("IDEN-Project-UAS-Darmstadt/IPCWJK")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(IPCWJK)
## basic example code
```

# Development

Restore the development environment with:

``` r
renv::restore()
```

then use everything available in the `devtools` package to develop the
package.

``` r
library(devtools)
document() # to update documentation and roxygen functionality
load_all() # to load the package functions for development
build_readme() # to update the README
test() # to run tests
check() # to check the package
covr::package_coverage() # to check code coverage
styler::style_pkg() # to style the code
lint() # to check the code for linting issues
```
