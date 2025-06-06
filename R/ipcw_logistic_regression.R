#' IPCW Logistic Regression
#'
#' @description
#' \loadmathjax
#' Fits a logistic regression model with IPCW for
#' right-censored survival data. Jackknife model training is performed.
#'
#' @details
#' Training is performed using the `glm`
#' using a quasibinomial family to account for the weights.
#'
#' @inheritParams ipcw_xgboost
#' @inherit ipcw_xgboost return
#' @seealso [ipcw_weights()] for the underlying implementation of the weights
#' and [IPCWJK] for more information.
#' @importFrom stats glm as.formula predict
#' @importFrom Rdpack reprompt
#' @import mathjaxr
#' @family IPCW models
#' @examples
#' library(survival)
#' tau <- 100
#' df <- veteran[, c("time", "status", "trt")]
#' newdata <- data.frame(trt = c(1, 2))
#'
#' fit <- ipcw_logistic_regression(df,
#'   tau = tau, time_var = "time",
#'   status_var = "status"
#' )
#' predict(fit, newdata)
#' @export
ipcw_logistic_regression <- function(
    data, tau, time_var = "t", status_var = "delta") {
  n <- nrow(data)
  w <- ipcw_weights(data, tau, time_var = time_var, status_var = status_var)
  w <- w / sum(w)
  y <- 1 * (data[[time_var]] > tau)

  training_vars <- setdiff(names(data), c(time_var, status_var))
  dtrain <- as.data.frame(cbind(data[, training_vars], y = y))
  colnames(dtrain) <- c(training_vars, "y")
  formula <- as.formula(paste("y ~", paste(training_vars, collapse = "+")))
  full_model <- glm(formula,
    data = dtrain, weights = w,
    family = "quasibinomial"
  )

  train_pred <- predict(full_model, dtrain, type = "response")
  train_brier <- sum(w * (train_pred - y)**2)
  jk_models <- list()

  for (i in c(1:n)) {
    wadj <- (n - 1) * w[-i] / (sum(w[-i]))
    dtraini <- as.data.frame(cbind(data[-i, training_vars], y = y[-i]))
    colnames(dtraini) <- c(training_vars, "y")
    model_wjk <- glm(formula,
      data = dtraini, weights = wadj,
      family = "quasibinomial"
    )
    jk_models[[i]] <- model_wjk
  }
  predfun <- function(model, x) {
    mat <- as.data.frame(x[, training_vars])
    colnames(mat) <- training_vars
    predict(model, mat, type = "response")
  }
  ipcwmodel(
    model_name = "IPCW Logistic Regression",
    full_model = full_model,
    jackknife_models = jk_models,
    predict = predfun,
    train_brier = train_brier,
    tau = tau,
    time_var = time_var,
    status_var = status_var,
    training_vars = training_vars,
    w = w,
    additional_information = list()
  )
}
