# simpredsel: Simple predictor selection

<!-- badges: start -->
<!-- badges: end -->

`simpredsel` is a small R package intended to help in the development of scales consisting of simple sum scores of items, which are potentially predictive for a given binary criterion variable. Such scales are common in clinical research. The central functions are

- `sim_pred_sel()` for finding predictors in a set that contribute to the predictive validity of the sum score and

- `mc_crossvalidation_sps()` for cross-validating predictor selections with `sim_pred_sel()`. 

For purposes of comparison, there is also a cross-validation function that uses stepwise ordinary or logistic regression for finding predictors, `mc_crossvalidation_regression()`.

## Installation

You can install `simpredsel` from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("mihock/simpredsel")
```

Under Windows may need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for this.
