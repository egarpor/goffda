## Constructor + S3 methods for the "flm_test" class

make_data <- function(seed = 1, n = 40, m = 21) {

  set.seed(seed)
  t <- seq(0, 1, l = m)
  X <- r_ou(n = n, t = t, sigma = 1)
  Y <- r_ou(n = n, t = t, sigma = 0.5)
  list(X = X, Y = Y, t = t)

}


test_that("flm_test inherits from htest and exposes the new class", {

  d <- make_data()
  out <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE)

  expect_s3_class(out, "flm_test")
  expect_s3_class(out, "htest")
  expect_named(out,
               c("statistic", "p.value", "boot_statistics", "method",
                 "parameter", "p", "q", "fit_flm", "boot_lambda",
                 "boot_p", "boot_E_scores", "data.name"),
               ignore.order = TRUE)
  expect_true(is.numeric(out$statistic))
  expect_true(out$p.value >= 0 && out$p.value <= 1)

})


test_that("flm_test deprecates plot_dens / plot_proc with a warning", {

  d <- make_data()

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_warning(
    flm_test(X = d$X, Y = d$Y, B = 20, verbose = FALSE, plot_dens = TRUE),
    regexp = "deprecated"
  )

})


test_that("print on flm_test uses htest formatting", {

  d <- make_data()
  out <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE)

  expect_output(print(out), regexp = "p-value")

})


test_that("summary.flm_test reports bootstrap moments", {

  d <- make_data()
  out <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE)
  s <- summary(out)

  expect_s3_class(s, "summary.flm_test")
  expect_true(!is.null(s$boot_summary))
  expect_named(s$boot_summary, c("mean", "sd", "median"))
  expect_true(is.finite(s$boot_summary[["sd"]]))

})


test_that("plot.flm_test runs for type 'dens' and 'proc'", {

  d <- make_data()
  out <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE)

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(out, type = "dens"))
  expect_silent(plot(out, type = "proc", plot_max_p = 1, plot_max_q = 1))

})


test_that("plot.flm_test errors when the requested data was pruned", {

  d <- make_data()
  out <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE,
                  save_boot_stats = FALSE)

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(out, type = "dens"), regexp = "save_boot_stats")

  out2 <- flm_test(X = d$X, Y = d$Y, B = 30, verbose = FALSE,
                   save_fit_flm = FALSE)
  expect_error(plot(out2, type = "proc"), regexp = "save_fit_flm")

})
