#' IPCW Model Class
#'
#' Constructs an object of class \code{ipcwmodel}.
#'
#' @description
#' \loadmathjax
#' The resulting object contains the fitted model, jackknife resamples,
#' prediction function, training Brier score, and metadata about the model.
#'
#' @details
#' Models can be used by calling the \code{predict(model, newdata)} method.
#' It returns a data frame with columns for the prediction, lower and upper
#' confidence intervals, and standard error.
#'
#' The confidence intervals and standard errors can be computed using either
#' the naive approach (using the number of jackknife models minus one as the
#' denominator) or the robust approach (using the IPCW weights). This
#' is controlled by the \code{naive} argument in the \code{predict} method.
#'
#' The confidence intervals then are computed using the Wald
#' approach (using a z-score, default is 1.96 for 95% CI).
#'
#' @param model_name Character. The name of the model.
#' @param full_model The fitted model object for the full dataset.
#' @param jackknife_models List of fitted models, each omitting one observation
#'   (the jackknife resamples).
#' @param tau Numeric. The time horizon at which survival is estimated.
#' @param predict Function. The prediction function for the fitted model.
#' @param train_brier Numeric. The Brier score for the fitted model on the
#'   training data.
#' @param time_var Character. The name of the time variable in the data.
#' @param status_var Character. The name of the status variable in the data.
#' @param training_vars Character vector. The names of the covariates used for
#'   model fitting.
#' @param w Numeric vector. The IPCW weights used for model fitting.
#' @param additional_information List. Additional information to be stored in
#'   the model object (default is an empty list).
#' @return An object of class \code{ipcwmodel} containing the model details.
#' @importFrom Rdpack reprompt
#' @import mathjaxr
#' @references
#' \insertAllCited{}
#' @family ipcwbaseclass
#' @export
ipcwmodel <- function(model_name, full_model, jackknife_models, tau,
                      predict, train_brier,
                      time_var, status_var,
                      training_vars,
                      w,
                      additional_information = list()) {
  stopifnot(is.character(model_name) && length(model_name) == 1)
  stopifnot(is.numeric(train_brier) && length(train_brier) == 1)
  stopifnot(is.numeric(tau) && length(tau) == 1)
  stopifnot(is.character(time_var) && length(time_var) == 1)
  stopifnot(is.character(status_var) && length(status_var) == 1)
  stopifnot(is.character(training_vars) && length(training_vars) >= 0)
  stopifnot(is.numeric(w) && length(w) > 0)
  stopifnot(length(w) == length(jackknife_models))
  stopifnot(is.list(additional_information))
  structure(
    list(
      model_name = model_name,
      full_model = full_model,
      jackknife_models = jackknife_models,
      predict = predict,
      train_brier = train_brier,
      tau = tau,
      time_var = time_var,
      status_var = status_var,
      training_vars = training_vars,
      w = w,
      additional_information = additional_information
    ),
    class = "ipcwmodel"
  )
}

#' @importFrom utils capture.output
#' @export
print.ipcwmodel <- function(x, ...) {
  cat("IPCW Model\n")
  cat("-----------\n")
  cat("Model name: ", x$model_name, "\n")
  cat("Tau (time horizon): ", x$tau, "\n")
  cat("Time variable: ", x$time_var, "\n")
  cat("Status variable: ", x$status_var, "\n")
  cat("Training variables: ", paste(x$training_vars, collapse = ", "), "\n")
  cat("Number of training samples: ", length(x$w), "\n")
  cat("Number of unusable training samples: ", sum(x$w == 0), "\n")
  w_notnorm <- x$w * length(x$w)
  cat("Number of (effective) training samples: ", sum(w_notnorm), "\n")
  cat("Train Brier score: ", round(x$train_brier, 3), "\n")
  if (length(x$additional_information) > 0) {
    cat("Additional information:\n")
    for (name in names(x$additional_information)) {
      cat(
        "  -", name, ":",
        paste(capture.output(print(x$additional_information[[name]])),
          collapse = "\n    "
        ), "\n"
      )
    }
  }
  invisible(x)
}

#' @export
#' @param object An object of class \code{ipcwmodel}.
#' @param newdata A data frame containing the covariates for which predictions
#'   are to be made.
#' @param naive Logical. If \code{TRUE}, use the naive jackknife variance
#'   estimator. If \code{FALSE}, use the IPCW-weighted estimator.
#' @param z Numeric. The z-score to use for the confidence interval.
#'   Default is 1.96, corresponding to a 95% confidence interval.
#' @param ... Additional arguments (currently ignored).
#' @importFrom stats pnorm
#' @describeIn ipcwmodel Predict method for \code{ipcwmodel} objects.
predict.ipcwmodel <- function(object, newdata, naive = FALSE,
                              z = 1.96, ...) {
  if (missing(newdata)) {
    stop("Argument 'newdata' is required for prediction.")
  }
  if (!is.data.frame(newdata)) {
    stop("Argument 'newdata' must be a data frame.")
  }
  # Check 'naive' argument
  if (!is.logical(naive) || length(naive) != 1 || is.na(naive)) {
    stop("Argument 'naive' must be a single logical value (TRUE or FALSE).")
  }
  # Check 'z' argument
  if (!is.numeric(z) || length(z) != 1 || is.na(z) || z <= 0) {
    stop("Argument 'z' must be a single positive numeric value.")
  }

  if (!all(object$training_vars %in% names(newdata))) {
    missing_vars <- setdiff(object$training_vars, names(newdata))
    stop(
      paste(
        "The following variables are missing in 'newdata':",
        paste(missing_vars, collapse = ", ")
      )
    )
  }
  n_preds <- nrow(newdata)
  n_models <- length(object$jackknife_models)

  pred_fun <- object$predict
  pred <- pred_fun(object$full_model, newdata)

  stopifnot(is.numeric(pred), length(pred) == n_preds)
  preds <- matrix(0, nrow = n_preds, ncol = n_models)
  for (i in seq_along(object$jackknife_models)) {
    preds[, i] <- pred_fun(object$jackknife_models[[i]], newdata)
  }

  jk_dev <- preds - matrix(pred, nrow = n_preds, ncol = n_models, byrow = FALSE)
  jk_dev <- jk_dev^2

  if (naive) {
    weights <- (n_models - 1) / (n_models)
    weights <- rep(weights, each = n_models)
  } else {
    weights <- ifelse(object$w == 0, 0, 1 - object$w)
  }
  se <- sqrt(rowSums(weights * jk_dev))
  wald_logit(pred, z, se)
}
