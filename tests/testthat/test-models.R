library(flexsurv)

lambda_c <- 0.10 # ~ 25% censoring rate
b0 <- 1
b1 <- -0.5
b2 <- 0.5
shape <- 1
tau <- 5
n <- 100

x1 <- rbinom(n, 1, 0.5)
x1 <- x1 - 1 * (x1 == 0)
x2 <- rbinom(n, 1, 0.5)
scalevec <- exp(b0 + x1 * b1 + x2 * b2)
t.star <- rllogis(n, shape = shape, scale = scalevec)
cens <- rexp(n, lambda_c)
delta <- 1 * (t.star <= cens)
t <- pmin(t.star, cens)
df <- data.frame(t, delta, x1, x2)
test_data <- data.frame(x1 = c(1, 1), x2 = c(1, 0))

models <- list(
  ipcw_xgboost = ipcw_xgboost,
  ipcw_logistic_regression = ipcw_logistic_regression
)


for (model_name in names(models)) {
  model <- models[[model_name]]
  test_that(paste("Model", model_name, "works"), {
    fit <- model(df, tau = tau, time_var = "t", status_var = "delta")
    expect_s3_class(fit, "ipcwmodel")
    preds <- predict(fit, test_data)
    expect_equal(nrow(preds), nrow(test_data))
  })
  fit <- model(df, tau = tau, time_var = "t", status_var = "delta")
  for (wald in c(TRUE, FALSE)) {
    for (naive in c(TRUE, FALSE)) {
      test_that(paste(
        "Model", model_name, "predicts with wald =",
        wald, "and naive =", naive
      ), {
        preds <- predict(fit, test_data, wald = wald, naive = naive)
        expect_s3_class(preds, "data.frame")
        expect_true(all(predcols %in% colnames(preds)))
        expect_equal(nrow(preds), nrow(test_data))
        expect_true(all(sapply(preds, is.numeric)))
        expect_true(all(apply(preds, 1, function(row) {
          row["lower"] <= row["prediction"] &&
            row["prediction"] <= row["upper"]
        })))
        expect_true(all(preds$se >= 0))
      })
    }
  }
}
