
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DISCLAIMER

This package/application is provided without warranty of any kind.

While every effort has been made to ensure that that routines provide
correct and accurate results, all responsibility for the use of the
results belongs to the user.

# Public Domain

This repository constitutes a work of the United States Government and
is not subject to domestic copyright protection under 17 USC § 105. This
repository is in the public domain within the United States, and
copyright and related rights in the work worldwide are waived through
the [CC0 1.0 Universal public domain
dedication](https://creativecommons.org/publicdomain/zero/1.0/). All
contributions to this repository will be released under the CC0
dedication. By submitting a pull request you are agreeing to comply with
this waiver of copyright interest.

# PooledInfRate

The PooledInfRate package provides functions to estimate proportions
from pooled or group testing data.

## Installation

You can install PooledInfRate from GitLab:

``` r
# To install from GitLab CDC
# library(devtools) # need to install this if you do not have it
# This will build on the fly, so you need
# RTools (available from CRAN)
# and need to set
# Sys.setenv(R_QPDF="true") # for the build
# then 
devtools::install_git("https://git.cdc.gov/bkb5/PooledInfRate")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PooledInfRate)
## basic example code
x <- c(1,0,0,0)
m <- c(50,25,10,5)
pooledBin(x,m)
#>            P        Lower        Upper 
#> 0.0125791125 0.0007261606 0.0782849595
```
