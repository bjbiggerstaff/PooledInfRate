
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PooledInfRate

<!-- badges: start -->
<!-- badges: end -->

The PooledInfRate package provides functions to estimate proportions
from pooled or group testing data.

## Installation

You can install PooledInfRate from GitLab:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
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
#>      P  Lower  Upper 
#> 0.0126 0.0007 0.0783
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
