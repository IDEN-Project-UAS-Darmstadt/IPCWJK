#' IPCW XGBoost Binary Classifier
#'
#' @description
#' \loadmathjax
#' Fits a binary classification model using XGBoost with IPCW for
#' right-censored survival data. Hyperparameter tuning and Jackknife
#' model training are performed.
#'
#' @details
#' Training is performed using the `xgboost` package
#' \insertCite{xgboost}{IPCWJK} based on the `"binary:logistic"` objective.
#'
#' Hyperparameter tuning is done using three (`nfold`) fold cross-validation
#' with a grid of parameters. The best parameters are selected based on the
#' minimum test log loss over 100 (`nrounds`) rounds with early stopping
#' (10 rounds, `early_stopping_rounds`).
#' Note that the tested hyperparameters are
#' based on our simulation and will not be useful for all datasets.
#'
#' The tested hyperparameters include:
#' * `booster`: `"gbtree"` or `"gblinear"`.
#' * `eta`: Learning rate, tested as `1 / 10^(0:5)`.
#' * For`booster="gblinear"`:
#'   * `max_depth`: Maximum depth of the tree, tested as `c(12, 6, 3, 1)`.
#'
#' With the best parameters, the model is trained on the full dataset.
#'
#' @inheritParams ipcw_weights
#' @param grid Data frame. Grid of hyperparameters to test in cross-validation.
#'   The default is the return of `ipcw_xgboost_default_grid()``.
#' @param nrounds Integer. Maximum number of boosting rounds for XGBoost
#'   training and cross-validation (default is 100).
#' @param early_stopping_rounds Integer. Number of rounds with no improvement
#'   to trigger early stopping during cross-validation (default is 10).
#' @param nfold Integer. Number of folds for cross-validation (default is 3).
#' @param verbose Integer. Verbosity level for XGBoost training and
#'   cross-validation (default is 0).
#' @param nthread Integer. Number of threads to use for XGBoost training
#'   (default is 1).
#' @return An object of class [ipcwmodel].
#' @seealso [ipcw_weights()] for the underlying implementation of the weights
#' and [IPCWJK] for more information.
#' @import xgboost
#' @importFrom Rdpack reprompt
#' @import mathjaxr
#' @references
#' \insertAllCited{}
#' @family IPCW models
#' @examples
#' library(survival)
#' tau <- 100
#' df <- veteran[, c("time", "status", "trt")]
#' newdata <- data.frame(trt = c(1, 2))
#'
#' fit <- ipcw_xgboost(df,
#'   tau = tau, time_var = "time",
#'   status_var = "status"
#' )
#' predict(fit, newdata)
#' @export
ipcw_xgboost <- function(
    data, tau, time_var = "t", status_var = "delta",
    verbose = 0, grid = ipcw_xgboost_default_grid(),
    nrounds = 100, early_stopping_rounds = 10,
    nfold = 3,
    nthread = 1) {
  n <- nrow(data)
  w <- ipcw_weights(data, tau, time_var = time_var, status_var = status_var)
  # Normalize weights
  w <- w / sum(w)
  # This invalid for the censored observations before tau,
  # but their weights are 0 anyway
  y <- 1 * (data[[time_var]] > tau)
  training_vars <- setdiff(names(data), c(time_var, status_var))
  dtrain <- as.matrix(data[, training_vars])
  dtrain <- xgb.DMatrix(
    data = dtrain, weight = w, nthread = nthread,
    label = y
  )

  if (!is.data.frame(grid)) {
    stop("The 'grid' argument must be a data.frame.")
  }
  params_totest <- grid
  cvfun <- function(...) {
    withCallingHandlers(
      xgb.cv(
        data = dtrain, early_stopping_rounds = early_stopping_rounds,
        nrounds = nrounds, objective = "binary:logistic",
        nfold = nfold, verbose = verbose, nthread = nthread, ...
      ),
      warning = function(w) {
        if (grepl("NaNs produced", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  cvres <- lapply(seq.int(nrow(params_totest)), function(i) {
    current_params <- as.list(params_totest[i, ])
    current_params <- current_params[!is.na(unlist(current_params))]
    do.call(
      cvfun, current_params
    )
  })
  besttestscores <- sapply(
    cvres,
    function(res) min(res$evaluation_log$test_logloss_mean)
  )
  touse <- cvres[[which.min(besttestscores)]]
  if (verbose > 0) {
    cat("Best parameters:\n")
    print(touse)
  }

  params_ignore <- c("weight", "nthread", "silent")
  params_passed <- touse$params[!names(touse$params) %in% params_ignore]
  full_model <- xgb.train(
    data = dtrain, params = params_passed,
    nrounds = touse$best_iteration, nthread = nthread
  )

  train_brier <- sum(w * (predict(full_model, dtrain) - y)**2)
  jk_models <- list()

  for (i in c(1:n)) {
    wadj <- (n - 1) * w[-i] / (sum(w[-i]))
    dtraini <- xgb.DMatrix(
      data = as.matrix(data[-i, training_vars]),
      label = y[-i], weight = wadj, nthread = nthread
    )
    model_wjk <- xgb.train(
      data = dtraini, params = params_passed,
      nrounds = touse$best_iteration, nthread = nthread
    )
    jk_models[[i]] <- model_wjk
  }
  predfun <- function(model, x) {
    mat <- as.matrix(x[, training_vars])
    predict(model, mat)
  }
  ipcwmodel(
    model_name = "IPCW XGBoost",
    full_model = full_model,
    jackknife_models = jk_models,
    predict = predfun,
    train_brier = train_brier,
    tau = tau,
    time_var = time_var,
    status_var = status_var,
    training_vars = training_vars,
    w = w,
    additional_information = list(params_passed = params_passed)
  )
}
#' @describeIn ipcw_xgboost Returns a default grid of hyperparameters.
#' @export
ipcw_xgboost_default_grid <- function() {
  etas <- 1 / 10^(c(0:5))
  params_totest <- expand.grid(
    booster = "gbtree", eta = etas,
    max_depth = c(12, 6, 3, 1), stringsAsFactors = FALSE
  )
  params_totest <- rbind(params_totest, expand.grid(
    booster = "gblinear", eta = etas, max_depth = NA, stringsAsFactors = FALSE
  ))
  params_totest
}
