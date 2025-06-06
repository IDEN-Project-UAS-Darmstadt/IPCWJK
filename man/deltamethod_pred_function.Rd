% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltamethod_pred_function.R
\name{deltamethod_pred_function}
\alias{deltamethod_pred_function}
\title{Predict the standard error from a prediction formula using the delta method}
\usage{
deltamethod_pred_function(
  prediction_str,
  coefs,
  coef_cov,
  additional_coefs = character(),
  fixed_vars = numeric(),
  logit = FALSE
)
}
\arguments{
\item{prediction_str}{Character. A string specifying the prediction formula,
which must contain one occurrence of \code{"LP"} to be replaced by
the linear predictor constructed from the coefficients and covariate
values.}

\item{coefs}{Named numeric vector. The estimated coefficients from the fitted
model. The names must correspond to the covariate names used in the model.}

\item{coef_cov}{Square numeric matrix. The covariance matrix of the estimated
coefficients in \code{coefs}. The row and column names
must match the names of \code{coefs}.}

\item{additional_coefs}{Character vector. Names of coefficients in
\code{coefs} that are not part of the linear predictor but are required in
the prediction formula (default is an empty character vector). An
example is the estimated scale from a parametric survival model.}

\item{fixed_vars}{Named numeric vector. Fixed values for variables used in
the prediction formula but not present in the new data (default is an empty
numeric vector).}

\item{logit}{Logical. \code{TRUE}, the delta method is used on the
logit scale. This ensures CIs are between 0 and 1 (default is
\code{TRUE}).}
}
\value{
A \verb{function(df, z = 1.96)} that takes a data frame of covariates
\code{df} and returns a data frame with columns for the
prediction \code{prediction}, lower \code{lower} and upper \code{upper} Wald
confidence intervals, and standard error \code{se}. The function can also take
an optional argument \code{z} for the z-score used in confidence
interval calculation (default is 1.96 for 95\% confidence intervals).
}
\description{
\loadmathjax
Constructs a prediction function that computes model predictions,
their standard errors and Wald confidence interval using the delta method
for a user-specified prediction formula.
}
\details{
This function is intended for use with fitted
regression models, including those for survival analysis, and allows
flexible specification of the linear predictor, additional coefficients,
and fixed covariate values. The resulting function can be applied to
new data to obtain predictions, delta method-based standard
errors and confidence intervals.

This function is used internally in the package by
\code{\link[=deltamethod_from_model]{deltamethod_from_model()}}.

For a prediction function
\mjeqn{g(\theta)}{g(theta)}, the standard error is approximated as
\mjdeqn{\sqrt{\nabla g(\theta)^\top \Sigma \nabla g(\theta)}}{
sqrt(grad g(theta)^T Sigma grad g(theta))}
where \mjeqn{\Sigma}{Sigma} is the covariance matrix of the estimated
coefficients.

The input \code{prediction_str} must contain one occurrence of
\code{"LP"}, which will be replaced by the linear predictor constructed
from the provided coefficients and covariate values.

With \code{logit=TRUE} the confidence intervals can be calculated on a logit
scale, see \link{IPCWJK} for more information.
}
\examples{
coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "sex" = -0.2)
coef_cov <- diag(c(0.01, 0.0025, 0.0025))
rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
pred_fun <- deltamethod_pred_function(
  prediction_str = "1 / (1 + exp(-(LP)))",
  coefs = coefs,
  coef_cov = coef_cov
)
newdata <- data.frame(age = c(50, 60), sex = c(1, 0))
pred_fun(newdata)
}
