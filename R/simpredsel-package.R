#' @details The functions are intended to help in the development of scales consisting of simple sum scores of items, which are potentially predictive for a given binary criterion variable. Such scales are common in clinical research. The central functions are [sim_pred_sel] for finding predictors in a set that contribute to the predictive validity of the sum score and [mc_crossvalidation_sps] for cross-validating predictor selections with `sim_pred_sel`. For comparison, there is also a cross-validation function that uses stepwise ordinary or logistic regression for finding predictors, see [mc_crossvalidation_regression].

"_PACKAGE"

## usethis namespace: start
#' @importFrom graphics hist lines
#' @importFrom stats cor var as.formula coef glm predict reformulate update lm step
## usethis namespace: end
NULL
