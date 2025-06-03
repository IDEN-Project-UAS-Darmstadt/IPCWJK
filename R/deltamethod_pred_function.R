#' Apply the delta method based on a prediction formula
#'
#' @description
#' \loadmathjax
#' Constructs a prediction function that computes model predictions and
#' their standard errors and Wald confidence interval using the delta method
#' for a user-specified prediction formula.
#'
#' @details
#' This function is intended for use with fitted
#' regression models, including those for survival analysis, and allows
#' flexible specification of the linear predictor, additional coefficients,
#' and fixed covariate values. The resulting function can be applied to
#' new data to obtain both predictions and delta method-based standard
#' errors.
#'
#' This function is used internally in the package to provide
#' delta method-based standard errors for as a comparison.
#'
#' For a prediction function
#' \mjeqn{g(\theta)}{g(theta)}, the standard error is approximated as
#' \mjdeqn{\sqrt{\nabla g(\theta)^\top \Sigma \nabla g(\theta)}}{
#' sqrt(grad g(theta)^T Sigma grad g(theta))}
#' where \mjeqn{\Sigma}{Sigma} is the covariance matrix of the estimated
#' coefficients.
#'
#' The input \code{prediction_str} must contain one occurrence of
#' \code{"LP"}, which will be replaced by the linear predictor constructed
#' from the provided coefficients and covariate values.
#'
#' @param prediction_str Character. A string specifying the prediction formula,
#'   which must contain one occurrence of \code{"LP"} to be replaced by
#'   the linear predictor constructed from the coefficients and covariate
#'   values.
#' @param coefs Named numeric vector. The estimated coefficients from the fitted
#'   model. The names must correspond to the covariate names used in the model.
#' @param coef_cov Square numeric matrix. The covariance matrix of the estimated
#'   coefficients. The row and column names must match the names of
#'   \code{coefs}.
#' @param additional_coefs Character vector. Names of coefficients in
#'   \code{coefs} that are not part of the linear predictor but are required in
#'   the prediction formula (default is an empty character vector).
#' @param fixed_vars List. Named list of fixed values for variables used in the
#'   prediction formula but not present in the new data (default is an empty
#'   list).
#' @return A function that takes a data frame of covariates and returns a
#'   data frame with columns for the prediction, lower and upper Wald
#'   confidence intervals, and standard error. The function can also take
#'   an optional argument \code{z} for the z-score used in confidence
#'   interval calculation (default is 1.96 for 95% CI).
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
#' @family helpers
#' @examples
#' coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "sex" = -0.2)
#' coef_cov <- diag(c(0.01, 0.0025, 0.0025))
#' rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
#' pred_fun <- deltamethod_pred_function(
#'   prediction_str = "1 / (1 + exp(-(LP)))",
#'   coefs = coefs,
#'   coef_cov = coef_cov
#' )
#' newdata <- data.frame(age = c(50, 60), sex = c(1, 0))
#' pred_fun(newdata)
#' @export
deltamethod_pred_function <- function(
    prediction_str,
    coefs,
    coef_cov,
    additional_coefs = character(),
    fixed_vars = list()) {
  coef_names <- names(coefs)
  if (!all(additional_coefs %in% coef_names)) {
    stop("additional_coefs must be a subset of the names of coefs")
  }
  if (!all(coef_names %in% rownames(coef_cov))) {
    stop("coef_cov must have the same rownames as the names of coefs")
  }
  if (!all(coef_names %in% colnames(coef_cov))) {
    stop("coef_cov must have the same colnames as the names of coefs")
  }
  if (!all(rownames(coef_cov) == colnames(coef_cov))) {
    stop("coef_cov must be a square matrix with the same row
    and column names as coefs")
  }
  if (!all(rownames(coef_cov) %in% coef_names)) {
    stop("coef_cov must have the same names as the names of coefs")
  }
  if (length(grep("LP", prediction_str)) != 1) {
    stop("prediction_str must contain LP")
  }

  lpcoefs_sel <- !coef_names %in% additional_coefs
  lpcoefs <- coef_names[lpcoefs_sel]
  lpcoefs_form <- paste0("coef__", seq(0, length(lpcoefs) - 1))
  # skip the intercept
  lpcovar <- paste0("covariate__", seq(1, length(lpcoefs) - 1))
  # Add intercept multiplication
  lpcovar <- c("1", lpcovar)
  lpstring <- paste0(lpcoefs_form, " * ", lpcovar, collapse = " + ")
  prediction_f <- as.formula(paste0("~", gsub("LP", lpstring, prediction_str)))

  toanalyze <- c(lpcoefs_form, additional_coefs)
  toanalyze_old <- c(lpcoefs, additional_coefs)
  prediction_f_deriv <- deriv(prediction_f, toanalyze)


  # resort the cov matrix
  coef_cov_new <- coef_cov[toanalyze_old, toanalyze_old]
  colnames(coef_cov_new) <- toanalyze
  rownames(coef_cov_new) <- toanalyze

  prediction_fun <- function(xrow, z) {
    # skip the intercept, use the names in the formula
    newvals <- as.list(xrow[lpcoefs[2:length(lpcoefs)]])
    names(newvals) <- lpcovar[2:length(lpcovar)]
    trained_params <- coefs[toanalyze_old]
    names(trained_params) <- toanalyze
    envir <- as.list(c(newvals, fixed_vars, trained_params))
    envir <- list2env(envir)
    prediction <- eval(prediction_f_deriv, envir = envir)
    prediction_grad <- attr(prediction, "gradient")
    attr(prediction, "gradient") <- NULL
    ses <- sqrt(diag(prediction_grad %*% coef_cov_new %*% t(prediction_grad)))
    c("prediction" = prediction, "se" = ses)
  }

  prediction_fun_df <- function(df, z = 1.96) {
    preds <- apply(df, 1, prediction_fun, z = z)
    preds <- as.data.frame(t(preds))
    wald_logit(preds$prediction, z, preds$se)
  }
  return(prediction_fun_df)
}
