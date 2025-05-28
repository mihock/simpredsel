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
#' Find predictors in a set that contribute to the criterion validity of the sum score of the selected predictors.
#'
#' @details The search for predictors is done by forward selection. The procedure starts with selecting the predictor with the highest validity (i.e., highest association with the criterion) and then successively adds predictors to the selected set that maximize the *incremental* validity of the sum score.
#'
#' Association may be measured by the area under the ROC curve (AUC) or by the correlation coefficient. All predictors should be coded in the same direction as the criterion, so that the correlations between the predictors and the criterion are positive. When the association is measured by the AUC, the criterion should be binary, with its values coded 0 or 1 (indicator coding).
#'
#' Predictors with zero variance are removed from the predictor set before testing incremental validity. By default, predictors with nonpositive correlation with the criterion are also removed in advance. The latter may be included, however, this does normally not matter because predictors with negative correlations with the criterion rarely contribute to the validity of a sum score. (This would be a kind of suppressor effect.)
#'
#' For an example, see [simpredsel].
#'
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion must be binary (coded 0 for absence, 1 for presence of a feature).
#' @param only_positive consider only predictors with positive correlations with the criterion?
#' @param delta_val minimum validity increment that a chosen predictor must reach
#' @param show_progress show progress?
#' @param show_messages show messages?
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
sim_pred_sel <- function(x, criterion, assoc_measure = c("auc", "cor"), only_positive = TRUE, delta_val = 0, show_progress = TRUE, show_messages = TRUE) {
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
        if (show_messages) {
            message(k0, " predictor(s) removed because of zero variance")
        }
        n_pred <- ncol(x)
    }
    # Remove predictors with nonpositive correlation
    item_corrs <- cor(x, crit_var)[, 1]
    if (only_positive) {
        x <- x[which(item_corrs > 0)]
        if (show_messages) {
            k0 <- n_pred - ncol(x)
            if ( k0 > 0 ) {
                message(k0, " predictor(s) removed because of nonpositive correlation with the criterion")
            }
        }
    }
    # Do we still have predictors?
    if (ncol(x) == 0) {
        warning("no predictors with positive variance/correlation remaining")
        return(NULL)
    }

    # Main loop
    final_sum_score <- -1
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
        # Return if there is no validity increment of at least delta_val
        if (best_assoc <= (old_best_assoc + delta_val)) {
            names(old_best_assoc) <- "assoc"
            output <- list(
                k = i - 1,
                assoc = old_best_assoc,
                sel_pred_names = sel_pred_names,
                final_sum_score = final_sum_score
            )
            class(output) <- c("sim_pred_sel", class(output))
            invisible(output)
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
            myDf <- data.frame(i = i,
                variable = best_pred_name,
                cor_single = item_corrs[best_pred_name],
                assoc = best_assoc)
            if (i == 1) {
                progressDf <- myDf
            } else {
                progressDf <- rbind(progressDf, myDf)
            }
        }
    }
    if (show_progress) {
        print(knitr::kable(progressDf, format = "simple", digits = 3,
          row.names = FALSE,
          col.names = c("i", "Item", "Correlation/Item", "Association/Sum score")))
        cat("\n")
    }
    # Output / all predictors processed
    names(best_assoc) <- "assoc"
    output <- list(
        k = length(sel_pred_names),
        assoc = best_assoc,
        sel_pred_names = sel_pred_names,
        final_sum_score = final_sum_score
    )
    class(output) <- c("sim_pred_sel", class(output))
    return(output)
    # invisible(output)
}

### Not working anymore?
### @method print sim_pred_sel
### @noRd
#' @export
print.sim_pred_sel <- function(x, ...) {
    cat("k = ", x$k, ", association = ", round(x$assoc, 3), "\n", sep = "")
}

#' Monte Carlo Cross-Validation for Simple Predictor Selection
#'
#' Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections with [sim_pred_sel].
#'
#' For an example, see [simpredsel].
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param n number of Monte Carlo runs (i.e., training/validation samples drawn)
#' @param assoc_measure type of association measure. May be `auc` (area under the ROC) or `cor` (correlation). For the former, the criterion should be binary (coded for 0 absence, 1 for presence of a feature).
#' @param only_positive consider only predictors with positive correlations with the criterion?
#' @param delta_val minimum validity increment that a chosen predictor must reach
#' @param show_progress show progress from [sim_pred_sel]?
#' @param show_messages show messages from [sim_pred_sel]?
#'
#' @returns A list containing the components
#'
#' - `k` (number of predictors identified in each run),
#' - `assoc_measure` (type of association measure),
#' - `assoc` (associations in the training set),
#' - `assoc_valid` (associations in the validation set).
#'
#' The components are vectors with each value representing the result of one Monte Carlo run.
#'
#' @references Xu, Q.-S., & Liang, Y.-Z. (2001). Monte Carlo cross validation. *Chemometrics and Intelligent Laboratory Systems*, *56*(1), 1–11. https://doi.org/10.1016/S0169-7439(00)00122-2
#'
#' @export
mc_crossvalidation_sps <- function(x, criterion, n = 100L, assoc_measure = c("auc", "cor"), only_positive = TRUE, delta_val = 0, show_progress = FALSE, show_messages = FALSE) {
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
            delta_val = delta_val,
            show_progress = show_progress,
            show_messages = show_messages)
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
    output <- list(k = k,
        assoc_measure = assoc_measure,
        assoc_train = assoc_train,
        assoc_valid = assoc_valid)
    class(output) <- c("mc_crossvalidation", class(output))
    return(output)
}

#' Monte Carlo Cross-Validation for Predictor Selection with Linear or Logistic Regression
#'
#' Stratified Monte Carlo (repeated random sub-sampling) cross-validation for predictor selections with [stats::step]. The use of ordinary or logistic regression depends on the criterion values. Logistic regression is used when the criterion values are binary and coded with 0 and 1. Otherwise, ordinary regression is used.
#'
#' For an example, see [simpredsel].
#'
#' @param x data frame containing predictors and criterion
#' @param criterion character string specifying the criterion (must be in `x`)
#' @param n number of Monte Carlo runs (i.e., training/validation samples drawn)
#' @param only_positive keep only predictors with positive regression coefficients in the training run in the model?
#' @param show_progress show progress?
#'
#' @returns A list containing the components
#'
#' - `k` (number of predictors identified in each run),
#' - `method` (logistic or ordinary regression),
#' - `assoc_train` (associations in the training set),
#' - `assoc_valid` (associations in the validation set),
#' - `nonpos_vars_excluded` (number of variables excluded due to nonpositive regression coefficients).
#'
#' The components are vectors with each value representing the result of one Monte Carlo run. In the case of logistic regression, associations are measured by the area under the ROC curve (AUC), otherwise associations are measured by the correlation coefficient.
#'
#' @references Xu, Q.-S., & Liang, Y.-Z. (2001). Monte Carlo cross validation. *Chemometrics and Intelligent Laboratory Systems*, *56*(1), 1–11. https://doi.org/10.1016/S0169-7439(00)00122-2
#'
#' @export
mc_crossvalidation_regression <- function(x, criterion, n = 100L, only_positive = TRUE,  show_progress = TRUE) {
    stopifnot(is.data.frame(x),
        is.character(criterion), criterion %in% names(x),
        is.numeric(n), n > 0)
    logistic <- FALSE
    if (all(x[[criterion]] %in% c(0, 1))) {
        logistic <- TRUE
    }

    k <- numeric(n)
    assoc_train <- numeric(n)
    assoc_valid <- numeric(n)
    nonpos_vars_excluded <- numeric(n) ######
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
        final <- stats::step(res, trace = 0)
        if (only_positive) {
            # Identify variables with nonpositive coefficients
            nonpositive_vars <- names(coef(final))[coef(final) <= 0]
            # Exclude the intercept
            nonpositive_vars <- nonpositive_vars[nonpositive_vars != "(Intercept)"]
            nonpos_vars_excluded[i] <- length(nonpositive_vars)

            # Update the model by removing variables with nonpositive coefficients
            if (nonpos_vars_excluded[i] > 0) {
                # Notice that the model may now contain *new* nonpositive coefficients
                final <- update(final, as.formula(paste(". ~ .",
                    paste(nonpositive_vars, collapse = " - "), sep = " - ")))
            }
        }

        # Compute assoc for training data
        fitted <- predict(final)
        if (logistic) {
            assoc_train[i] <- pROC::auc(train_x[[criterion]], fitted, direction = "<",
                quiet = TRUE)
        } else {
            assoc_train[i] <- cor(train_x[[criterion]], fitted)
        }
        k[i] <- length(names(coef(final))) - 1

        # Validation
        fitted <- predict(final, newdata = valid_x)
        if (logistic) {
            assoc_valid[i] <- pROC::auc(valid_x[[criterion]], fitted, direction = "<",
                quiet = TRUE)
        } else {
            assoc_valid[i] <- cor(valid_x[[criterion]], fitted)
        }
        if (is.null(k[i]) || k[i] == 0) {
            warning("k (number of predictors) is 0 in MC iteration step ", i)
            k[i] <- 0
            assoc_train[i] <- NA ######
            assoc_valid[i] <- NA ######
        }
    }
    if (show_progress) cat("\n")
    output <- list(
        method = ifelse(logistic,
            "logistic regression (association measure: AUC)",
            "ordinary regression (association measure: correlation)"),
        k = k,
        assoc_train = assoc_train,
        assoc_valid = assoc_valid,
        nonpos_vars_excluded = nonpos_vars_excluded)
    class(output) <- c("mc_crossvalidation", class(output))
    return(output)
}

#' @export
print.mc_crossvalidation <- function(x, ...) {
  if ("assoc_measure" %in% names(x)) {
    # x stems from mc_crossvalidation_sps
    cat("Measure: ", x$assoc_measure, "\n")
  } else {
    # x stems from mc_crossvalidation_regression
    cat("Method: ", x$method, "\n")
  }
  cat(
    "Mean k                = ", round(mean(x$k, na.rm = TRUE), 3), "\n",
    "SD k                  = ", round(stats::sd(x$k, na.rm = TRUE), 3), "\n",
    "Mean Assoc/Training   = ", round(mean(x$assoc_train, na.rm = TRUE), 3), "\n",
    "SD Assoc/Training     = ", round(stats::sd(x$assoc_train, na.rm = TRUE), 3), "\n",
    "Mean Assoc/Validation = ", round(mean(x$assoc_valid, na.rm = TRUE), 3), "\n",
    "SD Assoc/Validation   = ", round(stats::sd(x$assoc_valid, na.rm = TRUE), 3), "\n",
    sep = "")
  if ("method" %in% names(x)) {
    cat("Non-positive predictors excluded: ", min(x$nonpos_vars_excluded), " - ",
      max(x$nonpos_vars_excluded), ", M = ", mean(x$nonpos_vars_excluded), "\n",
      sep = "")
  }
}

#' Plot Results From a Cross-Validation
#'
#' Plot results from a cross-validation with [mc_crossvalidation_sps] or [mc_crossvalidation_regression]: a barplot for the distribution of the number of items in the training sets and two histograms for the distributions of the association measures in the training and the validation set.
#'
#' @param x results from one of the cross-validation functions
#' @param ... may include `binwidth` (binwidth for plot of histograms of the associations, default: 0.01), `fill` (fill color for bars, default: "steelblue"), `col_mean` (color of lines indicating means, default: "black")
#'
#' @export
#' @importFrom rlang .data
#' @importFrom ggplot2 aes facet_wrap geom_histogram geom_vline ggplot labs
#' @importFrom dplyr summarise group_by
#'
plot.mc_crossvalidation <- function(x, ...) {
  arglist <- list(...)
  if (is.null(arglist$binwidth)) arglist$binwidth <- 0.01
  if (is.null(arglist$fill)) arglist$fill <- "steelblue"
  if (is.null(arglist$col_mean)) arglist$col_mean <- "black"

  if ("assoc_measure" %in% names(x)) {
    # x stems from mc_crossvalidation_sps
    xlab <- ifelse(x$assoc_measure == "auc", "AUC", "Correlation")
  } else {
    # x stems from mc_crossvalidation_regression
    xlab <- ifelse(substr(x$method, 1, 4) == "logi", "AUC", "Correlation")
  }

  # Barplot for k
  graphics::barplot(table(x$k), xlab  = "Number of Selected Predictors",
    ylab = "Frequency")

  # Histograms using ggplot2
  n <- length(x$assoc_train)
  mydat <- data.frame(source = c(rep("Training", n), rep("Validation", n)),
    assoc = c(x$assoc_train, x$assoc_valid))
  p <- ggplot(mydat, aes(x = .data$assoc )) +
    geom_histogram(binwidth = arglist$binwidth, col = "black", fill = arglist$fill) +
    labs(x = xlab, y = "Frequency") +
    facet_wrap(~ source, nrow = 1) ###### 'source' okay?
  means <- dplyr::group_by(mydat, source) |> dplyr::summarise(m = mean(.data$assoc)) ###### 'source' okay?
  p <- p + geom_vline(aes(xintercept = .data$m), means, col = arglist$col_mean, linetype = "dotted")
  print(p)

  message("crossvalidation: added two summary plots (see plot history)")
}

# #' @export
# plot.mc_crossvalidation_regression <- function(x, ...) {
#     ### The next line is the only difference to plot.mc_crossvalidation_sps.
#     xlab <- ifelse(substr(x$method, 1, 4) == "logi", "AUC", "Correlation")
#     graphics::barplot(table(x$k), xlab  = "Number of Selected Predictors", ylab = "requency")
#     h <- graphics::hist(x$assoc_train, main = "Training", xlab = xlab)
#     m <- mean(x$assoc_train, na.rm = TRUE)
#     lines(x = c(m, m), y = c(0, max(h$counts)), col = "red")
#     h <- graphics::hist(x$assoc_valid, main = "Validation", xlab = xlab)
#     m <- mean(x$assoc_valid, na.rm = TRUE)
#     lines(x = c(m, m), y = c(0, max(h$counts)), col = "red")
#     message("crossvalidation: added three summary plots (see plot history)")
# }

