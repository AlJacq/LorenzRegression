
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LorenzRegression

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlJacq/LorenzRegression/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AlJacq/LorenzRegression/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The LorenzRegression package proposes a toolbox to estimate, produce
inference on and interpret Lorenz regressions. These regressions are
used to determine the explanatory power of a set of covariates on the
inequality of a response variable.

## Installation

You can install the released version of LorenzRegression from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LorenzRegression")
```

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AlJacq/LorenzRegression")
```

## What’s new

### Version 2.0.0 (2024-08-01)

- **Breaking Changes**:
  - The package structure has been standardized and builds upon existing
    libraries. Please review the updated function documentation for
    current usage.
  - `Lorenz.Reg`: The function structure has changed. It now acts as a
    wrapper for the fitting functions `Lorenz.GA`, `Lorenz.FABS`and
    `Lorenz.SCADFABS`. It returns an object of class `''LR''`
    (non-penalized regression) or `''PLR''` (penalized regression). Each
    class has a set of designated methods.
  - `Lorenz.boot`: This function performs bootstrap calculations for
    objects of class `''LR''` or `''PLR''`. It uses the `boot` function
    from the `boot` package. For penalized regression, it also computes
    an out-of-bag score for tuning parameter selection. The function
    returns the updated object with bootstrap results and adds the class
    `''LR_boot''` (non-penalized regression) or `''PLR_boot''`
    (penalized regression).
  - `PLR.CV`: This function performs cross-validation for objects of
    class `''PLR''`. It computes a cross-validation score for tuning
    parameter selection. The folds are constructed using `vfold_cv` from
    the `rsample` package. The function returns the updated object with
    cross-validation results and adds the class `''PLR_cv''`.
  - Method availability for classes `''LR''` and `''PLR''` is now
    documented in the `Lorenz.Reg` help page. Each method also has its
    own help page.
- **New Features**:
  - The `grid.arg` and `grid.value` arguments in `Lorenz.Reg` allow
    users to specify one tuning parameter for penalized regression and
    construct a grid for it. Fitting is repeated for each grid value,
    and optimal values are determined using available methods (among
    BIC, bootstrap and cross-validation).
  - `diagnostic.PLR` provides diagnostic information for penalized
    Lorenz regression.
