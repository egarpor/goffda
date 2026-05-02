## Constructor + S3 methods for the "flm_est" class
## Covers the four input branches: X (functional/scalar) x Y (functional/scalar)

make_funct_funct <- function(seed = 12345, n = 40, m = 31) {

  set.seed(seed)
  t <- seq(0, 1, l = m)
  X <- r_ou(n = n, t = t, sigma = 1)
  E <- r_ou(n = n, t = t, sigma = 0.3)
  Y <- 2 * X + E
  list(X = X, Y = Y, t = t)

}


make_funct_scalar <- function(seed = 12345, n = 40, m = 31) {

  set.seed(seed)
  t <- seq(0, 1, l = m)
  X <- fdata_cen(r_ou(n = n, t = t, sigma = 1))
  beta <- r_ou(n = 1, t = t, sigma = 0.3, x0 = 1)
  Y <- drop(inprod_fdata(X_fdata1 = X, X_fdata2 = beta)) + rnorm(n, sd = 0.1)
  list(X = X, Y = Y, t = t)

}


make_scalar_funct <- function(seed = 12345, n = 40, m = 31) {

  set.seed(seed)
  t <- seq(0, 1, l = m)
  X <- rnorm(n)
  beta <- r_ou(n = 1, t = t, sigma = 0.3, x0 = 1)
  beta$data <- matrix(beta$data, nrow = n, ncol = ncol(beta$data),
                      byrow = TRUE)
  E <- r_ou(n = n, t = t, sigma = 0.2)
  Y <- X * beta + E
  list(X = X, Y = Y, t = t)

}


make_scalar_scalar <- function(seed = 12345, n = 40) {

  set.seed(seed)
  X <- rnorm(n)
  Y <- 2 * X + rnorm(n, sd = 0.1)
  list(X = X, Y = Y)

}


test_that("flm_est (functional X, functional Y) returns an flm_est object", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr_l1s")

  expect_s3_class(fit, "flm_est")
  expect_true(!isTRUE(fit$scalar_X))
  expect_true(!isTRUE(fit$scalar_Y))
  expect_true(fda.usc::is.fdata(fit$Y_hat))
  expect_true(fda.usc::is.fdata(fit$residuals))

})


test_that("flm_est handles the four predictor/response branches", {

  expect_s3_class(flm_est(X = make_funct_funct()$X,
                          Y = make_funct_funct()$Y,
                          est_method = "fpcr"),
                  "flm_est")
  expect_s3_class(flm_est(X = make_funct_scalar()$X,
                          Y = make_funct_scalar()$Y,
                          est_method = "fpcr"),
                  "flm_est")
  expect_s3_class(flm_est(X = make_scalar_funct()$X,
                          Y = make_scalar_funct()$Y,
                          est_method = "fpcr"),
                  "flm_est")
  expect_s3_class(flm_est(X = make_scalar_scalar()$X,
                          Y = make_scalar_scalar()$Y,
                          est_method = "fpcr"),
                  "flm_est")

})


test_that("predict.flm_est matches fitted() on training data (4 branches)", {

  ## Functional X, functional Y
  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")
  pred <- predict(fit, newdata = d$X)
  expect_true(fda.usc::is.fdata(pred))
  expect_lt(max(abs(pred$data - fit$Y_hat$data)), 1e-6)

  ## Functional X, scalar Y
  d <- make_funct_scalar()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")
  pred <- predict(fit, newdata = d$X)
  expect_type(pred, "double")
  expect_length(pred, length(d$Y))
  expect_lt(max(abs(pred - fit$Y_hat)), 1e-6)

  ## Scalar X, functional Y
  d <- make_scalar_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")
  pred <- predict(fit, newdata = d$X)
  expect_true(fda.usc::is.fdata(pred))
  expect_lt(max(abs(pred$data - fit$Y_hat$data)), 1e-6)

  ## Scalar X, scalar Y
  d <- make_scalar_scalar()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")
  pred <- predict(fit, newdata = d$X)
  expect_type(pred, "double")
  expect_length(pred, length(d$Y))
  expect_lt(max(abs(pred - fit$Y_hat)), 1e-6)

})


test_that("predict.flm_est without newdata equals fitted()", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")

  expect_identical(predict(fit), fitted(fit))

})


test_that("predict.flm_est validates newdata", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")

  # Wrong type for functional predictor
  expect_error(predict(fit, newdata = c(1, 2, 3)), regexp = "fdata")

  # Scalar fit but functional newdata
  d2 <- make_scalar_scalar()
  fit2 <- flm_est(X = d2$X, Y = d2$Y, est_method = "fpcr")
  expect_error(predict(fit2, newdata = d$X), regexp = "numeric")

})


test_that("fitted/residuals error if compute_residuals = FALSE", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr",
                 compute_residuals = FALSE)

  expect_error(fitted(fit), regexp = "compute_residuals")
  expect_error(residuals(fit), regexp = "compute_residuals")

})


test_that("coef returns the requested representation", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")

  expect_true(is.matrix(coef(fit, type = "matrix")))
  expect_true(is.matrix(coef(fit, type = "scores")))

})


test_that("summary.flm_est returns a sensible R-squared", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")
  s <- summary(fit)

  expect_s3_class(s, "summary.flm_est")
  expect_true(s$r_squared >= 0 && s$r_squared <= 1 + 1e-8)

})


test_that("print and plot smoke tests on flm_est objects", {

  d <- make_funct_funct()
  fit <- flm_est(X = d$X, Y = d$Y, est_method = "fpcr")

  expect_output(print(fit), regexp = "flm_est")
  expect_output(print(summary(fit)), regexp = "Summary")

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(fit))

})
