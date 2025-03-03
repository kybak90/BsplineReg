
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BsplineReg

<!-- badges: start -->
<!-- badges: end -->

The goal of BsplineReg is to fit a regression spline estimator based on
a set of data $\{(x_i, y_i)\}_{i=1}^n$.

## Installation

You can install the development version of BsplineReg from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("kybak90/BsplineReg")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BsplineReg)
set.seed(923)  
n = 30
x_values = sort(runif(n, 0, 1))
y_values = sin(2 * pi * x_values) + cos(4 * pi * x_values) + rnorm(n, sd = 0.2)

# spline degree specification
degree = 3

# knot generation
num_interior_knots = 5
interior_knots = knots_quantile(x_values, num_interior_knots)

# model fitting
model = fit_spline(x_values, y_values, interior_knots, degree)

grid_x = seq(min(x_values), max(x_values), length.out = 100)
plot_spline(x_values, y_values, model, grid_x)
```

<img src="man/figures/README-example-1.png" width="100%" />
