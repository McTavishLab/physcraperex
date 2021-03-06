
# Physcraper examples

<!-- badges: start -->
<!-- badges: end -->

The goal of `physcraperex` R package is to present tools for visualizing and analysing results from Physcraper, a python toolkit for updating existing phylogenetic knowledge with new molecular data.

Physcraper's documentation is available at [https://physcraper.readthedocs.io](https://physcraper.readthedocs.io/en/latest/index.html).

## Installing physcraperex

For now, you can only install the development version of the `physcraperex` package with 

``` r
devtools::install_github("McTavishLab/physcraperex")
```

For instructions on how to install Physcraper go to the [documentation website](https://physcraper.readthedocs.io/en/latest/install.html).

<!--You can install the released version of physcraperex from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("physcraperex")
```
-->

## Examples


Browse the Articles section of the [`physcraperex` website](https://mctavishlab.github.io/physcraperex/) to find various examples showing usage of available visualization functions.

<!--
This is a basic example which shows you how to solve a common problem:

``` r
library(physcraperex)
## basic example code
```
-->


## Building this website

To build this website, we used the R package `pkgdown`. 
We ran once:

```
usethis::use_pkgdown()
```

And then the following every time we wanted to update the package:

```
pkgdown::build_site()
```
