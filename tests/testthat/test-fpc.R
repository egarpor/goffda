## Constructor + S3 methods for the "fpc" class

make_fdata <- function(seed = 1, n = 30, m = 21) {

  set.seed(seed)
  r_ou(n = n, t = seq(0, 1, l = m), sigma = 1)

}


test_that("fpc() returns an object with all expected fields and class", {

  X <- make_fdata()
  fp <- fpc(X_fdata = X, n_fpc = 5)

  expect_s3_class(fp, "fpc")
  expect_named(fp, c("d", "rotation", "scores", "l", "equispaced"),
               ignore.order = TRUE)
  expect_length(fp$d, 5)
  expect_equal(dim(fp$scores), c(30, 5))
  expect_true(fda.usc::is.fdata(fp$rotation))
  expect_equal(length(fp$l), 5)
  expect_type(fp$equispaced, "logical")

})


test_that("predict.fpc round-trips on the training data", {

  X <- make_fdata()
  fp <- fpc(X_fdata = X, n_fpc = 5)

  scores_new <- predict(fp, newdata = X)
  expect_equal(dim(scores_new), c(30, 5))
  # Round-trip: scores recomputed from centered training data must match
  expect_lt(max(abs(scores_new - fp$scores)), 1e-8)

})


test_that("predict.fpc validates inputs", {

  X <- make_fdata()
  fp <- fpc(X_fdata = X, n_fpc = 3)

  # Wrong newdata type
  expect_error(predict(fp, newdata = matrix(0, 5, 21)),
               regexp = "fdata")
  # argvals mismatch
  X2 <- r_ou(n = 5, t = seq(0, 1, l = 11))
  expect_error(predict(fp, newdata = X2), regexp = "argvals")
  # object of wrong class
  fp_unclassed <- unclass(fp)
  expect_error(predict.fpc(fp_unclassed, newdata = X), regexp = "fpc")

})


test_that("summary.fpc produces consistent variance proportions", {

  X <- make_fdata()
  fp <- fpc(X_fdata = X, n_fpc = 5)
  s <- summary(fp)

  expect_s3_class(s, "summary.fpc")
  expect_equal(s$n_fpc, 5)
  expect_equal(rownames(s$importance),
               c("Standard deviation", "Proportion of Variance",
                 "Cumulative Proportion"))
  cum_prop <- s$importance["Cumulative Proportion", ]
  expect_true(all(diff(cum_prop) >= -1e-12))
  expect_lte(max(cum_prop), 1 + 1e-12)

})


test_that("print and plot smoke tests on fpc objects", {

  X <- make_fdata()
  fp <- fpc(X_fdata = X, n_fpc = 5)

  expect_output(print(fp), regexp = "Functional principal components")
  expect_output(print(summary(fp)), regexp = "Importance")

  # plot should run without error (output suppressed)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(fp))
  expect_silent(plot(fp, ind = 1:2))
  expect_error(plot(fp, ind = 100))

})
