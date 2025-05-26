#' @details The functions are intended to help in the development of scales consisting of simple sum scores of items, which are potentially predictive for a given binary criterion variable. Such scales are common in clinical research. The central functions are [sim_pred_sel] for finding predictors in a set that contribute to the predictive validity of the sum score and [mc_crossvalidation_sps] for cross-validating predictor selections with `sim_pred_sel`.
#'
#' For purposes of comparison, there is also a cross-validation function that uses stepwise ordinary or logistic regression for finding predictors, see [mc_crossvalidation_regression]. The two cross-validation functions use Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections.
#'
#' @importFrom MASS mvrnorm
#' @examples
#' # Generate a toy data set from a correlation matrix with 10 predictors
#' # (V1 to V10) and a criterion variable (crit).
#'
#' # Construct correlation matrix ----------------------------------------
#' n <- 200
#' k <- 10 # number of items
#' S <- diag(k)
#' S[upper.tri(S)][1:10] <- 0.3
#' S[lower.tri(S)] <- t(S)[lower.tri(S)]
#' S
#'
#' # Generate data ----------------------------------------
#' set.seed(123)
#' dat <- as.data.frame(MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = S))
#' # Binary criterion
#' dat$crit <- ifelse(rowMeans(dat[, 1:3],) + rnorm(n, 0, 1) > 0, 1, 0)
#' round(cor(dat), 2)
#'
#' # Analysis ----------------------------------------
#'
#' # The following would do one run:
#' # sim_pred_sel(dat, "crit")
#'
#' # Cross-validation runs multiple (by default 100 ) times and yields
#' # summary statistics and plots of the distribution of the interesting
#' # statistics.
#'
#' # Simple predictor selection:
#' res_mc <- mc_crossvalidation_sps(dat, criterion = "crit")
#' res_mc
#' plot(res_mc)
#'
#' # Selection via 'stats::step' for logistic regression:
#' res_mcr <- mc_crossvalidation_regression(dat, criterion = "crit")
#' res_mcr
#' plot(res_mcr)

"_PACKAGE"

## usethis namespace: start
#' @importFrom graphics hist lines
#' @importFrom stats cor var as.formula coef glm predict reformulate update lm step
## usethis namespace: end
NULL
