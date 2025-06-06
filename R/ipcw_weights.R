#' Use a Kaplan-Meier estimator to estimate IPCW weights
#'
#' @description
#' \loadmathjax
#' Computes inverse probability of censoring weights (IPCW) for right-censored
#' survival data using the [pec::ipcw()] function based on a Kaplan-Meier
#' estimate of the censoring distribution \mjeqn{C}{C}.
#'
#' @details
#' Weighting observations by their Inverse Probability of Censoring has been
#' proposed to prevent bias introduced by removing censored individuals
#' \insertCite{Vock2016,Ginestet2021,Blanche2023}{IPCWJK} and improve model
#' performance.
#'
#' \mjdeqn{
#' \tilde{w}_i =
#' \left\lbrace
#' \begin{array}{ll}
#'  0, & c_i < \tau \wedge t^*_i \cr
#'  \frac{1}{P(C > \tau \mid X = x_i)}, & \tau < c_i \wedge t^*_i \cr
#'  \frac{1}{P(C > t_i \mid X = x_i)}, & t^*_i < c_i \wedge \tau
#' \end{array}
#' \right.
#' }{w_i_tilde = 0, if c_i < tau and t*_i; 1 / P(C > tau | X = x_i),
#' if tau < c_i and t*_i; 1 / P(C > t_i | X = x_i), if t*_i < c_i and tau }
#'
#' The function uses the [pec::ipcw()] \insertCite{pec}{IPCWJK} function to
#' compute the IPCW weights. They are not normalized (do **not** sum
#' to one).
#'
#' @param data A data frame containing the survival data. Must include columns
#'   for the observed time and event indicator.
#' @param tau Numeric scalar. The time horizon at which the survival
#'   probability is to be estimated.
#' @param time_var Character. The name of the variable in `data`
#'   representing the observed time to event or censoring.
#'   Default is `"t"`.
#' @param status_var Character. The name of the variable in `data`
#'   representing the event indicator (1 if event occurred, 0 if censored).
#'   Default is `"delta"`.
#' @return A numeric vector of IPCW weights, ordered as in the original
#'  `data`. The weights are not normalized (do not sum to one).
#' @seealso [IPCWJK] for more information and [ipcwmodel] for implementations
#'  of ready to use models.
#' @importFrom pec ipcw
#' @importFrom prodlim Hist
#' @importFrom stats as.formula
#' @importFrom Rdpack reprompt
#' @import mathjaxr
#' @references
#' \insertAllCited{}
#' @examples
#' data <- data.frame(
#'   t = c(5, 8, 12, 15, 20),
#'   delta = c(1, 0, 1, 0, 1)
#' )
#' w <- ipcw_weights(data, tau = 10)
#' @export
ipcw_weights <- function(data, tau, time_var = "t", status_var = "delta") {
  if (!is.numeric(tau) || length(tau) != 1) {
    stop("'tau' must be a single numeric value.")
  }
  # All other values of tau are valid, but we require it to be positive
  if (tau <= 0) {
    stop("'tau' must be a positive numeric value.")
  }
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!is.character(time_var) || !is.character(status_var)) {
    stop("'time_var' and 'status_var' must be character strings.")
  }
  if (length(time_var) != 1 || length(status_var) != 1) {
    stop("'time_var' and 'status_var' must be single character strings.")
  }
  if (!all(c(time_var, status_var) %in% names(data))) {
    stop(paste(
      "Data frame must contain columns:", time_var, "and",
      status_var
    ))
  }
  if (any(is.na(data[[time_var]])) || any(is.na(data[[status_var]]))) {
    stop("Input data contains NA values in time or status variables.")
  }
  if (!is.numeric(data[[time_var]]) || !is.numeric(data[[status_var]])) {
    stop("Time and status variables must be numeric.")
  }
  if (!all(data[[status_var]] %in% c(0, 1))) {
    stop("Status variable must be binary (0 or 1).")
  }
  if (!all(data[[time_var]] >= 0)) {
    stop("Time variable must be non-negative.")
  }
  if (nrow(data) == 0) {
    stop("Input data frame is empty.")
  }
  if (sum(data[[status_var]]) == nrow(data)) {
    # If all subjects are events, we cannot compute IPCW weights
    warning("All subjects have events. IPCW weights cannot be computed.")
    return(rep(1, nrow(data)))
  }

  # Use time, status and then row number to order the data to address ties
  ordered <- order(data[[time_var]], -data[[status_var]], seq_len(nrow(data)))
  data <- data[ordered, ]

  # If the time and status variables cause problems for pec, we will
  # let pec handle the error messages
  w <- pec::ipcw(
    as.formula(paste0("Hist(", time_var, ", ", status_var, ") ~ 1")),
    data = data,
    times = unique(data[[time_var]]),
    subjectTimes = data[[time_var]],
    method = "marginal"
  )
  wtau <- pec::ipcw(
    as.formula(paste0("Hist(", time_var, ", ", status_var, ") ~ 1")),
    data = data,
    times = tau,
    what = "IPCW.times",
    method = "marginal"
  )$IPCW.times

  w <- w$IPCW.subjectTimes
  w[data[[time_var]] > tau] <- rep(wtau, sum(data[[time_var]] > tau))
  w <- 1 / w
  w[!data[[status_var]] & data[[time_var]] < tau] <- 0

  if (all(w == 0)) {
    warning("All weights are zero.
    This may indicate that all subjects were censored before tau.")
  }

  w <- w[order(ordered)]
  return(w)
}
