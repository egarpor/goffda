## Constructor + S3 methods for the "cv_glmnet" class

make_xy <- function(seed = 123, n = 50, p = 4, q = 2) {

  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n)
  beta <- matrix(rnorm(p * q), nrow = p)
  y <- x %*% beta + matrix(rnorm(n * q, sd = 0.1), nrow = n)
  list(x = x, y = y)

}


test_that("cv_glmnet returns an object with the expected class and fields", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge")

  expect_s3_class(cv, "cv_glmnet")
  expect_named(cv, c("beta_hat", "intercept", "lambda", "cv", "fit"),
               ignore.order = TRUE)
  expect_equal(dim(cv$beta_hat), c(4, 2))
  expect_equal(dim(cv$intercept), c(1, 2))
  expect_type(cv$lambda, "double")
  expect_length(cv$lambda, 1)

})


test_that("cv_glmnet skips CV when lambda is supplied", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge", lambda = 0.1)

  expect_null(cv$cv)
  expect_equal(cv$lambda, 0.1)

})


test_that("predict.cv_glmnet matches manual prediction", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge")

  newx <- matrix(rnorm(20), nrow = 5, ncol = 4)
  yhat <- predict(cv, newx = newx)

  manual <- newx %*% as.matrix(cv$beta_hat)
  manual <- sweep(manual, 2, as.numeric(cv$intercept), "+")
  expect_equal(yhat, manual, tolerance = 1e-12)

})


test_that("predict.cv_glmnet validates dimensions", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge")
  expect_error(predict(cv, newx = matrix(0, nrow = 3, ncol = 99)),
               regexp = "columns")

})


test_that("coef.cv_glmnet shapes match intercept toggle", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge")

  with_int <- coef(cv, intercept = TRUE)
  no_int <- coef(cv, intercept = FALSE)

  expect_equal(nrow(with_int), nrow(no_int) + 1)
  expect_equal(ncol(with_int), ncol(no_int))

})


test_that("print and plot smoke tests on cv_glmnet objects", {

  d <- make_xy()
  cv <- cv_glmnet(x = d$x, y = d$y, alpha = "ridge")

  expect_output(print(cv), regexp = "cv_glmnet")

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(cv))

})
