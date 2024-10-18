#' Simple Predictor Selection
#'
#' The package provides functions for simple (naive) predictor selection and cross-validation.
#'
#' @name simpredsel
#'
#' @importFrom stats cor var
### #' @importFrom splitTools partition
### #' @importFrom pROC auc
NULL

#' Compute AUCs for one or more predictors
#'
#' Compute AUCs with the first component of `x` as response and each of the remainings components
#' as predictors. The first component must be binary (coded as 0, 1).
#'
#' @param x a data frame
#' @returns a vector of AUC values
#'
#' @export
compute_aucs <- function(x) {
    stopifnot(is.data.frame(x))
    myauc <- numeric(ncol(x)-1)
    for (i in 2:ncol(x)) {
        myauc[i - 1] <- pROC::auc(x[[1]], x[[i]], direction = "<", quiet = TRUE)
    }
    names(myauc) <- names(x)[2:ncol(x)]
    myauc
}

#' Select predictors that contribute to the predictive validity of a sum score of the predictors
#'
#' Forward selection of predictors that maximize the association of a sum score of the selected predictors with a criterion. All predictors should be coded in the same direction as the criterion, so that the (true) correlation is positive.
#'
#' @param predictors data frame of predictor variables
#' @param criterion a numeric vector representing the criterion (dependent variable)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion should be binary (coded 0 for absence, 1 for presence of a feature).
#' @param only_positive consider only predictors with nonnegative correlations with the criterion?
#' @param show_progress show progress?
#'
#' @returns A list containing the components
#'
#' - `k` (the number of predictors used for the final sum score),
#' - `assoc` (the association of the final sum score with the criterion),
#' - `sel_pred_names` (the names of the selected predictors),
#' - `final_sum_score` (the final sum score of the predictors)
#'
#' @export
forward_selection_preds <- function(predictors, criterion, assoc_measure = c("auc", "cor"), only_positive = TRUE, show_progress = TRUE) {
    # requires function `compute_aucs` (which requires pROC)
    assoc_measure <- match.arg(assoc_measure)
    stopifnot(is.data.frame(predictors), is.atomic(criterion))
    # Remove predictors that have no variance.
    # Should that be done *before* the CV iterations?
    n_pred <- ncol(predictors)
    item_vars <- diag(var(predictors))
    predictors <- predictors[which(item_vars > 0)]
    k0 <- n_pred - ncol(predictors)
    if (k0 > 0) {
        message(k0, " predictors removed because of zero variance")
        n_pred <- ncol(predictors)
    }
    # Remove predictors with nonpositive correlation
    item_corrs <- cor(predictors, criterion)[,1]
    if (only_positive) {
        predictors <- predictors[which(item_corrs > 0)]
        k0 <- n_pred - ncol(predictors)
        if (k0 > 0) {
            message(k0, " predictors removed because of nonpositive correlation with the criterion")
        }
    }
    # Main loop
    best_pred <- 0
    old_best_assoc <- 0
    sel_pred_names <- character() # names of selected predictors
    updated_preds <- predictors
    remaining_preds <- predictors
    for (i in 1:ncol(predictors)) {
        # Add the values of the best predictor from the last step (i-1)
        # to all predictors. Thus, 'updated_preds' are the new
        # candidates for the best sum score.
        updated_preds <- updated_preds + best_pred
        # Compute associations of all remaining predictors with the criterion.
        if (assoc_measure == "auc") {
            assocs <- compute_aucs(cbind(criterion, updated_preds))
        } else {
            assocs <- cor(updated_preds, criterion)[,1]
        }
        # Select the best predictor from the remaining predictors.
        best_pred_idx <- which.max(assocs)
        best_pred <- remaining_preds[, best_pred_idx]
        best_assoc <- assocs[best_pred_idx]
        best_pred_name <- names(assocs)[best_pred_idx]

        if (best_assoc < old_best_assoc) {
            output <- list(k = i - 1,
                assoc = old_best_assoc,
                sel_pred_names = sel_pred_names,
                final_sum_score = final_sum_score)
            return(output)
        } else {
            if (i == 1) {
                sel_pred_names <- best_pred_name
            } else {
                sel_pred_names <- c(sel_pred_names, best_pred_name)
            }
        }
        final_sum_score <- updated_preds[, best_pred_idx]
        remaining_preds[best_pred_idx] <- NULL
        updated_preds[best_pred_idx] <- NULL
        old_best_assoc <- best_assoc
        if (show_progress) {
            if (i == 1) {
                cat("No | Variable | Assoc/sum score | Correlation/single predictor | \n")
            }
            cat(i, best_pred_name, best_assoc, item_corrs[best_pred_name], "\n", sep = " | ")
        }
    }
}

#' Monte Carlo Cross-Validation
#'
#' Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections with `forward_selection_preds`.
#'
#' @param x data frame with predictors
#' @param criterion criterion (to be predicted variable)
#' @param n number of Monte Carlo runs (i.e., training/validation samples drawn)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion should be binary (coded for 0, 1 for presence of a feature).
#'
#' @returns A list containing the components
#'
#' - `assoc` (associations in the training set),
#' - `assoc_valid` (associations in the validation set), and
#' - `k` (number of predictors identified in each run).
#'
#' The components are vectors with each value representing the result of one Monte Carlo run.
#'
#' @export
mc_crossvalidation <- function(x, criterion, n = 100L, assoc_measure = c("auc", "cor")) {
    # Needs
    #  - library splitTools
    #  - function forward_selection_preds
    assoc_measure <- match.arg(assoc_measure)
    stopifnot(is.data.frame(x), is.atomic(criterion), is.numeric(n), n > 0)

    k <- numeric(n)
    assoc <- numeric(n)
    assoc_valid <- numeric(n)
    for (i in 1:n) {
        mysplit <- splitTools::partition(criterion,
            p = c(train = 0.5, valid = 0.5))
        train_x <- x[mysplit$train,]
        criterion_train <- criterion[mysplit$train]
        valid_x <- x[mysplit$valid,]
        criterion_valid <- criterion[mysplit$valid]
        res <- forward_selection_preds(train_x, criterion_train,
            assoc_measure = assoc_measure)
        k[i] <- res$k
        assoc[i] <- res$assoc
        y <- rowSums(valid_x[res$sel_pred_names])
        if (assoc_measure == "auc") {
            assoc_valid[i] <- pROC::auc(criterion_valid, y, direction = "<",
                quiet = TRUE)
        } else {
            assoc_valid[i] <- cor(criterion_valid, y)
        }
    }
    return(list(assoc = assoc, assoc_valid = assoc_valid, k = k))
}
