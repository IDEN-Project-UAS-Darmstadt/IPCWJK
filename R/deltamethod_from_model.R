#' Predict the standard error from a fitted model using the delta method
#'
#' @description
#' \loadmathjax
#' Constructs a prediction function that estimates the survival probability
#' at a specified time horizon \mjeqn{|tau}{tau} from a fitted model,
#' and computes its standard error and Wald confidence interval
#' using the delta method.
#'
#' @details
#' The function supports models of from [survival::survreg()]
# \insertCite{survival-package}{IPCWJK} with
#' log-logistic distribution and models of class `binreg`
#' (such as those fitted by [mets::logitIPCW()])
#' \insertCite{Blanche2023,mets1,mets2}{IPCWJK}. For `binreg` models,
#' the function can use either the naive or the robust variance estimator,
#' depending on the value of the `naive` argument.
#'
#' @param model A fitted model object. Supported types are
#'   [survival::survreg()] with log-logistic distribution and
#'  [mets::logitIPCW()].
#' @inheritParams ipcw_weights
#' @param naive Logical. If `TRUE`, use the naive variance
#'   estimator for `binreg` models. If `FALSE` (default),
#'   use the robust variance estimator.
#' @inherit deltamethod_pred_function return
#' @seealso [deltamethod_pred_function()] for the underlying implementation and
#' [IPCWJK] for confidence interval calculation
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#' @importFrom stats vcov coef deriv
#' @references
#' \insertAllCited{}
#' @examples
#' library(survival)
#' tau <- 100
#' df <- veteran[, c("time", "status", "trt")]
#' newdata <- data.frame(trt = c(1, 2))
#'
#' # Fit a log-logistic survival model
#' survreg_fit <- survreg(Surv(time, status) ~ trt,
#'   data = df,
#'   dist = "loglogistic"
#' )
#' pred_fun <- deltamethod_from_model(survreg_fit, tau = tau)
#' pred_fun(newdata)
#'
#' # Fit a logitIPCW model
#' library(mets)
#' logipcw_fit <- logitIPCW(Event(time, status) ~ trt, time = tau, data = df)
#' predfun_logit <- deltamethod_from_model(logipcw_fit, tau = tau)
#' pred_fun(newdata)
#'
#' @export
deltamethod_from_model <- function(model, tau, naive = FALSE) {
  if (inherits(model, "survreg") && model$dist == "loglogistic") {
    # Covariance matrix uses the log scale
    lscale <- log(model$scale)
    coefs <- c(coef(model), lscale = lscale)
    vcoef_cov <- vcov(model)
    colnames(vcoef_cov)[length(coefs)] <- "lscale"
    rownames(vcoef_cov)[length(coefs)] <- "lscale"
    f <- "1 / (1 + (tau / exp(LP))^(1 / exp(lscale)))"
    predfun <- deltamethod_pred_function(f,
      coefs,
      additional_coefs = c("lscale"),
      coef_cov = vcoef_cov,
      fixed_vars = c(tau = tau)
    )
  } else if (inherits(model, "binreg") && ("naive.var" %in% names(model))) {
    if (model$time != tau) {
      stop("The 'tau' argument must match the time used in the model.")
    }
    # Only this model has both the corrected and naive variance covariance
    if (naive) {
      vcoef_cov <- model$naive.var
    } else {
      vcoef_cov <- vcov(model)
    }
    rownames(vcoef_cov) <- colnames(vcoef_cov) <- names(coef(model))

    predfun <- deltamethod_pred_function("1-1 / (1 + exp(-(LP)))",
      coef(model),
      additional_coefs = character(),
      coef_cov = vcoef_cov,
      fixed_vars = c(tau = tau)
    )
  } else {
    stop("Unsupported model type. Only 'survreg' with 'loglogistic'
          distribution and 'logitIPCW' models are supported.")
  }
  return(predfun)
}
