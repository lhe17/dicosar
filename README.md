# DICOSAR v1.0.8

## Overview

*dicosar* is an R package for testing the equality of Pearson
correlation coefficients (PCCs) between two groups.

## Installation

### Most recent version

To install the latest version from github:

``` r
install.packages("devtools")
library(devtools)
install_github("lhe17/dicosar")
```

The installation requires Rcpp-1.0.7 and has been tested on R-4.1.0.

Please contact <liang.he@duke.edu> for more information.

## Basic usage

For the illustration, we first generate an example data set containing
100 samples per group. In this example, we assume that the two variables
follow a standard normal distribution and they are independent.

``` r
library(dicosar)
#> 
#> Attaching package: 'dicosar'
#> The following object is masked from 'package:base':
#> 
#>     seq
n = 100
x1 = matrix(rnorm(2*n),nrow=n)
x2 = matrix(rnorm(2*n),nrow=n)
```

Using the two matrices `x1` and `x2` as the input data, the following
script tests whether the PCCs are the same between the two groups.

``` r
p = dicosar(x1,x2)
p
#> $r1
#> [1] -0.01703903
#> 
#> $r2
#> [1] -0.1509276
#> 
#> $p
#> [1] 0.3879636
#> 
#> $p_delta
#> [1] NA
#> 
#> $code
#> [1] 1
```

The output shows the PCCs in the two groups, and the p-value for testing
the equality of the PCCs. The function will return NA if the sample size
is too small, or there are no enough different levels in the data, or
the algorithm does not converge. The detailed information can be checked
through `code`. Please see the manual for the definition of `code`.
DICOSAR also provides a p-value using the Delta method. This p-value can
be obtained using `delta=TRUE` as follows.

``` r
p = dicosar(x1,x2,delta=TRUE)
p
#> $r1
#> [1] -0.01703903
#> 
#> $r2
#> [1] -0.1509276
#> 
#> $p
#> [1] 0.3879636
#> 
#> $p_delta
#> [1] 0.3767165
#> 
#> $code
#> [1] 1
```

Here, `p_delta` is the p-value based on the Delta method. As shown in
our manuscript, the Delta method works quite well under a large sample
size, and is very fast. If the input matrices have more than two
variables (columns), DICOSAR will compute the p-value for each pair of
the variables and return three matrices for the pairwise PCCs and
p-values, respectively. Here, we show an example using four variables.

``` r
x1 = matrix(rnorm(4*n),nrow=n)
x2 = matrix(rnorm(4*n),nrow=n)
p = dicosar(x1,x2)
p
#> $r1
#>      [,1]       [,2]        [,3]        [,4]
#> [1,]   NA 0.05903522 -0.02908886 -0.07837627
#> [2,]   NA         NA -0.17183379 -0.10759844
#> [3,]   NA         NA          NA  0.01196397
#> [4,]   NA         NA          NA          NA
#> 
#> $r2
#>      [,1]        [,2]        [,3]         [,4]
#> [1,]   NA -0.08758217 -0.31576718 -0.258570221
#> [2,]   NA          NA -0.08022411 -0.009324102
#> [3,]   NA          NA          NA  0.098498614
#> [4,]   NA          NA          NA           NA
#> 
#> $p_single
#>      [,1]      [,2]       [,3]      [,4]
#> [1,]   NA 0.2881095 0.03831703 0.2246889
#> [2,]   NA        NA 0.48915950 0.4752875
#> [3,]   NA        NA         NA 0.6097066
#> [4,]   NA        NA         NA        NA
#> 
#> $p_global
#> [1] 0.1723248
```

The function also returns a p-value `p_global` of a global test for the
equality of these two correlation matrices using a Cauchy combination
test (CCT). The CCT is accurate when the number of aggregated p-values
is not large. If this number is very large, further assumptions about
the structure of the correlations are required by this test for
obtaining an accurate p-value.
