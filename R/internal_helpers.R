wald <- function(pred, z, se, logit) {
  if (logit) {
    pred_logit <- log(pred / (1 - pred))
    z_logit <- z * se / (pred * (1 - pred))
    lower <- exp(pred_logit - z_logit)
    upper <- exp(pred_logit + z_logit)
    lower <- lower / (1 + lower)
    upper <- upper / (1 + upper)
  } else {
    upper <- pred + z * se
    lower <- pred - z * se
  }
  data.frame(
    prediction = pred,
    lower = lower,
    upper = upper,
    se = se
  )
}
