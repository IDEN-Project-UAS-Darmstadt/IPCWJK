test_that("deltamethod_pred_function computes correct predictions and SE", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "sex" = -0.2)
  coef_cov <- diag(c(0.01, 0.0025, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP)))",
    coefs = coefs,
    coef_cov = coef_cov
  )
  newdata <- data.frame(age = c(50, 60), sex = c(1, 0))
  preds <- pred_fun(newdata)
  expect_s3_class(preds, "data.frame")
  expect_true(all(c("prediction", "lower", "upper", "se") %in% colnames(preds)))
  expect_equal(nrow(preds), 2)
  expect_true(all(preds$se > 0))
})

test_that("deltamethod_pred_function uses z argument for CI", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1)
  coef_cov <- diag(c(0.01, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP)))",
    coefs = coefs,
    coef_cov = coef_cov
  )
  newdata <- data.frame(age = 40)
  preds1 <- pred_fun(newdata, z = 1.96)
  preds2 <- pred_fun(newdata, z = 2.58)
  preds2 <- abs(preds2$upper - preds2$prediction)
  preds1 <- abs(preds1$upper - preds1$prediction)
  expect_true(preds1 < preds2)
})

test_that("deltamethod_pred_function errors for invalid additional_coefs", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1)
  coef_cov <- diag(c(0.01, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  expect_error(
    deltamethod_pred_function(
      prediction_str = "1 / (1 + exp(-(LP)))",
      coefs = coefs,
      coef_cov = coef_cov,
      additional_coefs = "not_a_coef"
    ),
    regexp = "additional_coefs"
  )
})

test_that("deltamethod_pred_function errors for coef_cov with wrong names", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1)
  coef_cov <- diag(c(0.01, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- c("wrong", "names")
  expect_error(
    deltamethod_pred_function(
      prediction_str = "1 / (1 + exp(-(LP)))",
      coefs = coefs,
      coef_cov = coef_cov
    ),
    regexp = "coef_cov"
  )
})

test_that("deltamethod_pred_function errors for prediction_str with no LP", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1)
  coef_cov <- diag(c(0.01, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  expect_error(
    deltamethod_pred_function(
      prediction_str = "1 / (1 + exp(-(NO)))",
      coefs = coefs,
      coef_cov = coef_cov
    ),
    regexp = "LP"
  )
})

test_that("deltamethod_pred_function uses fixed_vars correctly", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "sex" = -0.2)
  coef_cov <- diag(c(0.01, 0.0025, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP + grad)))",
    coefs = coefs,
    coef_cov = coef_cov,
    fixed_vars = c(grad = 0.3)
  )
  newdata <- data.frame(age = 40, sex = 1)
  preds <- pred_fun(newdata)
  expect_true(is.numeric(preds$prediction))
  lp <- 0.5 + 0.1 * 40 - 0.2 * 1
  expect_equal(preds$prediction, 1 / (1 + exp(-(lp + 0.3))))
})

test_that("deltamethod_pred_function uses additional_coefs correctly", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "sex" = -0.2, grad = 0.05)
  coef_cov <- diag(c(0.01, 0.0025, 0.0025, 1))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP)))+ grad",
    coefs = coefs,
    coef_cov = coef_cov,
    additional_coefs = "grad"
  )
  newdata <- data.frame(age = 40, sex = 1)
  preds <- pred_fun(newdata)
  expect_true(is.numeric(preds$prediction))
  lp <- (0.5 + 0.1 * 40 - 0.2 * 1)
  expect_equal(preds$prediction, 1 / (1 + exp(-lp)) + 0.05)
})

test_that("deltamethod_pred_function output is a data.frame", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1)
  coef_cov <- diag(c(0.01, 0.0025))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP)))",
    coefs = coefs,
    coef_cov = coef_cov
  )
  newdata <- data.frame(age = c(30, 40, 50))
  preds <- pred_fun(newdata)
  expect_s3_class(preds, "data.frame")
  expect_equal(nrow(preds), 3)
  expect_true(all(c("prediction", "lower", "upper", "se") %in% colnames(preds)))
})

test_that("deltamethod_pred_function handles additional_coefs", {
  coefs <- c("(Intercept)" = 0.5, "age" = 0.1, "extra" = 0.2)
  coef_cov <- diag(c(0.01, 0.0025, 0.001))
  rownames(coef_cov) <- colnames(coef_cov) <- names(coefs)
  pred_fun <- deltamethod_pred_function(
    prediction_str = "1 / (1 + exp(-(LP + extra)))",
    coefs = coefs,
    coef_cov = coef_cov,
    additional_coefs = "extra"
  )
  newdata <- data.frame(age = 40)
  preds <- pred_fun(newdata)
  expect_true(is.numeric(preds$prediction))
})
