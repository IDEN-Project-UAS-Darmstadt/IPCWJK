library(survival)
library(mets)

df <- data.frame(
  time = c(5, 10, 15, 20, 25, 30),
  status = c(1, 1, 0, 1, 0, 1),
  trt = c(1, 1, 2, 2, 1, 2)
)
test_data <- data.frame(trt = c(1, 2))
tau <- 12

test_that("deltamethod_from_model works for survreg loglogistic", {
  survreg_fit <- survreg(Surv(time, status) ~ trt,
    data = df,
    dist = "loglogistic"
  )
  pred_fun <- deltamethod_from_model(survreg_fit, tau = tau)
  expect_type(pred_fun, "closure")

  res <- pred_fun(test_data)
  expect_true(is.data.frame(res) || is.matrix(res))
  expect_true(all(predcols %in% colnames(res)))
})

test_that("deltamethod_from_model works", {
  logipcw_fit <- logitIPCW(Event(time, status) ~ trt, time = tau, data = df)
  # Robust
  pred_fun_robust <- deltamethod_from_model(logipcw_fit,
    tau = tau,
    naive = FALSE
  )
  expect_type(pred_fun_robust, "closure")
  res_robust <- pred_fun_robust(test_data)
  expect_true(is.data.frame(res_robust) || is.matrix(res_robust))
  expect_true(all(predcols %in% colnames(res_robust)))
  # Naive
  pred_fun_naive <- deltamethod_from_model(logipcw_fit,
    tau = tau,
    naive = TRUE
  )
  expect_type(pred_fun_naive, "closure")
  res_naive <- pred_fun_naive(test_data)
  expect_true(is.data.frame(res_naive) || is.matrix(res_naive))
  expect_true(all(predcols %in% colnames(res_naive)))
})

test_that("deltamethod_from_model errors for unsupported model types", {
  lm_fit <- lm(mpg ~ cyl, data = mtcars)
  expect_error(
    deltamethod_from_model(lm_fit, tau = 5),
    regexp = "Unsupported|supported"
  )
})

test_that("deltamethod_from_model errors if tau does not match model$time", {
  logipcw_fit <- logitIPCW(Event(time, status) ~ trt, time = tau, data = df)
  expect_error(
    deltamethod_from_model(logipcw_fit, tau = 10),
    regexp = "tau.*match|match.*tau"
  )
})

test_that("deltamethod_from_model errors for wrong argument types", {
  expect_error(
    deltamethod_from_model("not_a_model", tau = 5),
    regexp = "Unsupported|supported"
  )
})
