#' Compute AUCs for one or more predictors
#'
#' Compute AUCs (areas under the receiver operating characteristic curve) for the first component of `x` as response and each of its remainings components as predictors. The first component must be binary (coded as 0 for absence, and 1 for presence of a feature). The computation is done via [pROC::auc()].
#'
#' @param x a data frame
#' @returns a vector of AUC values
#'
#' @export
compute_aucs <- function(x) {
    stopifnot(is.data.frame(x), all(x[[1]] %in% c(0, 1)))
    myauc <- numeric(ncol(x)-1)
    for (i in 2:ncol(x)) {
        myauc[i - 1] <- pROC::auc(x[[1]], x[[i]], direction = "<", quiet = TRUE)
    }
    names(myauc) <- names(x)[2:ncol(x)]
    myauc
}

#' Simple Predictor Selection
#'
#' Find predictors in a set that contribute to the predictive validity of the sum score of the selected predictors.
#'
#' @details The search for predictors is done by forward selection. The procedure starts with selecting the predictor with the highest validity (i.e., highest association with the criterion) and then successively adds predictors to the selected set that maximize the *incremental* validity of the sum score.
#'
#' Association may be measured by the area under the ROC curve (AUC) or by the correlation coefficient. All predictors should be coded in the same direction as the criterion, so that the correlation between the predictor and the criterion is positive. When the association is measured by the AUC, the criterion should be binary, with its values coded 0 or 1 (indicator coding).
#'
#' By default, predictors with nonpositive correlations with the criterion are removed from the predictor set. Predictors with no variance are also removed.
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion must be binary (coded 0 for absence, 1 for presence of a feature).
#' @param only_positive consider only predictors with positive correlations with the criterion?
#' @param show_progress show progress?
#'
#' @returns A list containing the components
#'
#' - `k` (the number of predictors used for the final sum score),
#' - `assoc` (the association of the final sum score with the criterion),
#' - `sel_pred_names` (the names of the selected predictors),
#' - `final_sum_score` (the final sum score of the selected predictors)
#'
#' If no valid predictors are found NULL is returned.
#'
#' @export
sim_pred_sel <- function(x, criterion, assoc_measure = c("auc", "cor"), only_positive = TRUE, show_progress = TRUE) {
    assoc_measure <- match.arg(assoc_measure)
    stopifnot(is.data.frame(x),
        is.character(criterion), criterion %in% names(x))
    crit_var <- x[[criterion]]
    if (assoc_measure == "auc") {
        stopifnot(all(crit_var %in% c(0, 1)))
    }
    x[[criterion]] <- NULL
    # Remove predictors that have no variance.
    # Should that be done *before* the CV iterations?
    n_pred <- ncol(x)
    item_vars <- diag(var(x))
    x <- x[which(item_vars > 0)]
    k0 <- n_pred - ncol(x)
    if (k0 > 0) {
        message(k0, " predictor(s) removed because of zero variance")
        n_pred <- ncol(x)
    }
    # Remove predictors with nonpositive correlation
    item_corrs <- cor(x, crit_var)[, 1]
    if (only_positive) {
        x <- x[which(item_corrs > 0)]
        k0 <- n_pred - ncol(x)
        if (k0 > 0) {
            message(k0, " predictor(s) removed because of nonpositive correlation with the criterion")
        }
    }
    # Do we still have predictors?
    if (ncol(x) == 0) {
        warning("no predictors with positive variance/correlation remaining")
        return(NULL)
    }

    # Main loop
    best_pred <- 0
    old_best_assoc <- 0
    sel_pred_names <- character() # names of selected predictors
    updated_preds <- remaining_preds <- x
    for (i in 1:ncol(x)) {
        # Add the values of the best predictor from the last step (i-1)
        # to all predictors. Thus, 'updated_preds' are the new
        # candidates for the best sum score.
        updated_preds <- updated_preds + best_pred
        # Compute associations of all remaining predictors with the criterion.
        if (assoc_measure == "auc") {
            assocs <- compute_aucs(cbind(crit_var, updated_preds))
        } else {
            assocs <- cor(updated_preds, crit_var)[,1]
        }
        # Select the best predictor from the remaining predictors.
        best_pred_idx <- which.max(assocs)
        best_pred <- remaining_preds[, best_pred_idx]
        best_assoc <- assocs[best_pred_idx]
        best_pred_name <- names(assocs)[best_pred_idx]
        # Return if there is no validity increment
        if (best_assoc <= old_best_assoc) {
            names(old_best_assoc) <- "assoc"
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
        remaining_preds[best_pred_idx] <- updated_preds[best_pred_idx] <- NULL
        old_best_assoc <- best_assoc
        if (show_progress) {
            mydf <- data.frame(i = i, variable = best_pred_name, assoc = best_assoc,
                cor_single = item_corrs[best_pred_name])
            print(mydf)
        }
    }
    ##### ??? old_best_assoc?
    # Output / all predictors processed
    names(best_assoc) <- "assoc"
    output <- list(k = length(sel_pred_names),
        assoc = best_assoc,
        sel_pred_names = sel_pred_names,
        final_sum_score = final_sum_score)
    return(output)
}

#' Monte Carlo Cross-Validation for Simple Predictor Selection
#'
#' Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections with [sim_pred_sel].
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param n number of Monte Carlo runs (i.e., training/validation samples drawn)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion should be binary (coded for 0, 1 for presence of a feature).
#' @param only_positive consider only predictors with positive correlations with the criterion?
#' @param show_progress show progress?
#'
#' @returns A list containing the components
#'
#' - `assoc` (associations in the training set),
#' - `assoc_valid` (associations in the validation set), and
#' - `k` (number of predictors identified in each run).
#'
#' The components are vectors with each value representing the result of one Monte Carlo run.
#'
#' @references Xu, Q.-S., & Liang, Y.-Z. (2001). Monte Carlo cross validation. *Chemometrics and Intelligent Laboratory Systems*, *56*(1), 1–11. https://doi.org/10.1016/S0169-7439(00)00122-2
#'
#' @export
mc_crossvalidation_sps <- function(x, criterion, n = 100L, assoc_measure = c("auc", "cor"), only_positive = TRUE, show_progress = FALSE) {
    assoc_measure <- match.arg(assoc_measure)
    stopifnot(is.data.frame(x),
        is.character(criterion), criterion %in% names(x),
        is.numeric(n), n > 0)
    k <- numeric(n)
    assoc_train <- numeric(n)
    assoc_valid <- numeric(n)
    for (i in 1:n) {
        mysplit <- splitTools::partition(x[[criterion]],
            p = c(train = 0.5, valid = 0.5))
        train_x <- x[mysplit$train,]
        valid_x <- x[mysplit$valid,]
        res <- sim_pred_sel(train_x, criterion,
            assoc_measure = assoc_measure,
            only_positive = only_positive,
            show_progress = show_progress)
        if (is.null(res$k)) {
            warning("No results for predictor selection in MC iteration step ", i)
            k[i] <- 0
            assoc_train[i] <- NA
            y <- NA
            assoc_valid[i] <- NA
        } else {
            k[i] <- res$k
            assoc_train[i] <- res$assoc
            y <- rowSums(valid_x[res$sel_pred_names])
            if (assoc_measure == "auc") {
                stopifnot(all(valid_x[[criterion]] %in% c(0, 1)))
                assoc_valid[i] <- pROC::auc(valid_x[[criterion]], y, direction = "<",
                    quiet = TRUE)
            } else {
                assoc_valid[i] <- cor(valid_x[[criterion]], y)
            }
        }
    }
    return(list(assoc_measure = assoc_measure, assoc_train = assoc_train, assoc_valid = assoc_valid, k = k))
}

#' Monte Carlo Cross-Validation for Predictor Selection with Linear or Logistic Regression
#'
#' Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections with [MASS::stepAIC]. The use of ordinary or logistic regression depends on the criterion values. Logistic regression is used when the criterion values are binary and coded with 0 and 1. Otherwise, ordinary regression is used.
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param n number of Monte Carlo runs (i.e., training/validation samples drawn)
#' @param only_positive keep only predictors with positive regression coefficients during the initial training run in the model?
#' @param plot plot a histogram of the distribution of the associations for the validation sample?
#' @param show_progress show progress?
#'
#' @returns A list containing the components
#'
#' - `assoc_train` (associations in the training set),
#' - `assoc_valid` (associations in the validation set), and
#' - `k` (number of predictors identified in each run).
#'
#' The components are vectors with each value representing the result of one Monte Carlo run.
#'
#' @references Xu, Q.-S., & Liang, Y.-Z. (2001). Monte Carlo cross validation. *Chemometrics and Intelligent Laboratory Systems*, *56*(1), 1–11. https://doi.org/10.1016/S0169-7439(00)00122-2
#'
#' @export
mc_crossvalidation_regression <- function(x, criterion, n = 100L, only_positive = TRUE, plot = TRUE, show_progress = TRUE) {
    stopifnot(is.data.frame(x),
        is.character(criterion), criterion %in% names(x),
        is.numeric(n), n > 0)
    logistic <- FALSE
    if (all(x[[criterion]] %in% c(0, 1))) {
        logistic <- TRUE
    }
    if (logistic) {
        message("Using logistic regression.")
    } else {
        message("Using ordinary regression.")
    }

    k <- numeric(n)
    assoc_train <- numeric(n)
    assoc_valid <- numeric(n)
    for (i in 1:n) {
        if (show_progress) cat(i, ".", sep = "")
        # Create training and validation data frames
        mysplit <- splitTools::partition(x[[criterion]],
            p = c(train = 0.5, valid = 0.5))
        train_x <- x[mysplit$train,]
        valid_x <- x[mysplit$valid,]

        # Training
        myformula <- reformulate(names(train_x)[names(train_x) != criterion],
            response = criterion)
        if (logistic) {
            res <- glm(formula = myformula, family = "binomial", data = train_x)
        } else {
            res <- lm(formula = myformula, data = train_x)
        }
        final <- MASS::stepAIC(res, trace = 0)
        if (only_positive) {
            # Identify variables with negative coefficients
            negative_vars <- names(coef(final))[coef(final) < 0]
            # Exclude the intercept
            negative_vars <- negative_vars[negative_vars != "(Intercept)"]
            # Update the model by removing variables with negative coefficients
            if (length(negative_vars) > 0) {
                # Notice that the model may now contain *new* negative coefficients
                final <- update(final, as.formula(paste(". ~ .",
                    paste(negative_vars, collapse = " - "), sep = " - ")))
            }
        }

        # Compute assoc for training data
        fitted <- predict(final)
        if (logistic) {
            assoc_train[i] <- pROC::auc(train_x[[criterion]], fitted, direction = "<",
                quiet = TRUE)
        } else {
            assoc_train[i] <- cor(train_x[[criterion]], fitted)
            #####
            ####### include R2 from model
            # su <- summary(final)
            # R2 <- su$r.squared
            # adjR2 <- su$adj.r.squared
        }
        k[i] <- length(names(coef(final))) - 1

        # Validation
        fitted <- predict(final, newdata = valid_x)
        if (logistic) {
            assoc_valid[i] <- pROC::auc(valid_x[[criterion]], fitted, direction = "<",
                quiet = TRUE)
        } else {
            assoc_valid[i] <- cor(valid_x[[criterion]], fitted)
            #####
            ####### include R2 from model
            # su <- summary(final)
            # R2 <- su$r.squared
            # adjR2 <- su$adj.r.squared
        }
        if (is.null(k[i]) || k[i] == 0) {
            warning("k (number of predictors) is 0 in MC iteration step ", i)
            k[i] <- 0
            assoc_train[i] <- NA
            assoc_valid[i] <- NA
        }
    }
    if (show_progress) cat("\n")
    if (plot) {
        if (logistic) {
            cnt <- max( hist(assoc_valid, main = "AUC for validation data", xlab = "AUC")$counts )
            m <- mean(assoc_valid)
            lines(c(m, m), c(-0.5,cnt), col = "red")
        } else {
            cnt <- max( hist(assoc_valid, main = "Correlation for validation data", xlab = "Correlation")$counts )
            m <- mean(assoc_valid)
            lines(c(m, m), c(-0.5,cnt), col = "red")
        }
    }
    return(list(method = ifelse(logistic, "logistic regression", "ordinary regression"),
        assoc_train = assoc_train, assoc_valid = assoc_valid, k = k))
}

