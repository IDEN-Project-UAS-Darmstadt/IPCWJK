# R

test_that("ipcw_weights computes correct weights and handles censoring", {
  data <- data.frame(
    t = c(5, 8, 12, 15, 20),
    delta = c(1, 0, 1, 0, 1)
  )
  tau <- 10
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, nrow(data))
  expect_true(all(w >= 0))
  # Censored before tau (row 2) should be zero
  expect_equal(w[2], 0)
  # Events/censored after tau (rows 3,4,5) should have same weight
  expect_equal(w[3], w[4])
  expect_equal(w[4], w[5])
})

test_that("ipcw_weights errors for invalid tau", {
  data <- data.frame(t = 1:3, delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = "a"), "numeric")
  expect_error(ipcw_weights(data, tau = c(1, 2)), "single")
  expect_error(ipcw_weights(data, tau = 0), "positive")
  expect_error(ipcw_weights(data, tau = -1), "positive")
})

test_that("ipcw_weights errors for invalid data", {
  expect_error(ipcw_weights(list(a = 1), tau = 1), "data frame")
})

test_that("ipcw_weights errors for invalid time_var/status_var", {
  data <- data.frame(t = 1:3, delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = 1, time_var = 1), "character")
  expect_error(ipcw_weights(data, tau = 1, status_var = 1), "character")
  expect_error(ipcw_weights(data, tau = 1, time_var = "foo"), "columns")
  expect_error(ipcw_weights(data, tau = 1, status_var = "bar"), "columns")
})

test_that("ipcw_weights errors for NA in time/status", {
  data <- data.frame(t = c(1, NA, 3), delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = 1), "NA")
  data2 <- data.frame(t = c(1, 2, 3), delta = c(1, NA, 1))
  expect_error(ipcw_weights(data2, tau = 1), "NA")
})

test_that("ipcw_weights errors for non-numeric time/status", {
  data <- data.frame(t = c("a", "b", "c"), delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = 1), "numeric")
  dat2 <- data.frame(t = c(1, 2, 3), delta = c("x", "y", "z"))
  expect_error(ipcw_weights(dat2, tau = 1), "numeric")
})

test_that("ipcw_weights errors for non-binary status", {
  data <- data.frame(t = 1:3, delta = c(1, 2, 1))
  expect_error(ipcw_weights(data, tau = 1), "binary")
})

test_that("ipcw_weights errors for negative time", {
  data <- data.frame(t = c(-1, 2, 3), delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = 1), "non-negative")
})

test_that("ipcw_weights warns if all weights are zero", {
  data <- data.frame(t = c(1, 2, 3), delta = c(0, 0, 0))
  tau <- 5
  expect_warning(ipcw_weights(data, tau), "zero")
})

test_that("ipcw_weights handles tau equal to observed time", {
  data <- data.frame(
    t = c(5, 10, 15),
    delta = c(1, 0, 1)
  )
  tau <- 10
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights handles all events (no censoring)", {
  data <- data.frame(
    t = c(2, 4, 6),
    delta = c(1, 1, 1)
  )
  tau <- 5
  expect_warning(
    w <- ipcw_weights(data, tau),
    "events"
  )
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w == 1)) # All weights should be 1 since no censoring
})

test_that("ipcw_weights handles all censored after tau", {
  data <- data.frame(
    t = c(12, 15, 20),
    delta = c(0, 0, 0)
  )
  tau <- 10
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights handles duplicate times with different status", {
  data <- data.frame(
    t = c(5, 5, 10, 10),
    delta = c(1, 0, 1, 0)
  )
  tau <- 7
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 4)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights errors for empty data frame", {
  data <- data.frame(t = numeric(0), delta = numeric(0))
  expect_error(ipcw_weights(data, tau = 1), "empty", ignore.case = TRUE)
})

test_that("ipcw_weights handles tau greater than all observed times", {
  data <- data.frame(t = c(2, 4, 6), delta = c(1, 0, 1))
  tau <- 10
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights handles tau less than all observed times", {
  data <- data.frame(t = c(12, 14, 16), delta = c(1, 0, 1))
  tau <- 5
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights errors for factor time_var/status_var", {
  data <- data.frame(
    t = factor(c(1, 2, 3)),
    delta = c(1, 0, 1)
  )
  expect_error(ipcw_weights(data, tau = 1), "numeric")
  data2 <- data.frame(
    t = c(1, 2, 3),
    delta = factor(c(1, 0, 1))
  )
  expect_error(ipcw_weights(data2, tau = 1), "numeric")
})

test_that("ipcw_weights works with non-default column names", {
  data <- data.frame(time = c(1, 2, 3), status = c(1, 0, 1))
  tau <- 2
  w <- ipcw_weights(data, tau, time_var = "time", status_var = "status")
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights works with extra columns in data", {
  data <- data.frame(
    t = c(1, 2, 3),
    delta = c(1, 0, 1),
    x = c(10, 20, 30)
  )
  tau <- 2
  w <- ipcw_weights(data, tau)
  expect_type(w, "double")
  expect_length(w, 3)
  expect_true(all(w >= 0))
})

test_that("ipcw_weights errors for NA tau", {
  data <- data.frame(t = c(1, 2, 3), delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = NA), "numeric")
})

test_that("ipcw_weights errors for Inf or NaN in time/status", {
  data <- data.frame(t = c(1, Inf, 3), delta = c(1, 0, 1))
  expect_error(ipcw_weights(data, tau = 1), "NA|NaN|Inf")
  data2 <- data.frame(t = c(1, 2, 3), delta = c(1, NaN, 1))
  expect_error(ipcw_weights(data2, tau = 1), "NA|NaN|Inf")
})
