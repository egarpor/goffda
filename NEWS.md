# goffda 0.2.0

* New S3 classes and methods to expose an idiomatic R API:
  * `"fpc"` (already assigned by `fpc()`): adds `print`, `summary`,
    `plot`, and `predict` methods. `predict.fpc()` projects new
    `fdata` curves onto the stored eigenfunctions.
  * `"cv_glmnet"`: `cv_glmnet()` now returns an object of this class with
    `print`, `coef`, `predict`, and `plot` methods.
  * `"flm_est"`: `flm_est()` now returns an object of this class with
    `print`, `summary`, `coef`, `fitted`, `residuals`, `predict`, and
    `plot` methods. `predict.flm_est()` covers the four
    functional/scalar combinations of `X` and `Y`.
  * `c("flm_test", "htest")`: `flm_test()` is now a subclass of
    `htest`, so the standard `print` keeps working. New
    `summary.flm_test` and `plot.flm_test` are provided; `plot()` accepts
    `type = "dens"`, `"proc"`, or `"both"`.
* `flm_test()` now stores `boot_E_scores` (a 3-D array of bootstrap
  residual FPC scores) in the returned object, used by
  `plot.flm_test()` to draw the residual marked empirical processes.
  This is `NA` when `save_boot_stats = FALSE`.
* The arguments `plot_dens` and `plot_proc` of `flm_test()` are
  deprecated and now default to `FALSE`. Setting them to `TRUE`
  triggers a deprecation warning. Plotting should be performed via
  `plot()` on the returned object. The arguments will be removed in a
  future release.
* `flm_est()` now records `scalar_X` and `scalar_Y` flags in the
  returned object to support dispatch in `predict.flm_est()`.
* New unit tests under `tests/testthat/`.

# goffda 0.0.5

* Initial version

# goffda 0.0.6

* Comply with _R_CLASS_MATRIX_ARRAY_=true.
* Roxygen 7.0.1.

# goffda 0.0.7

* Added option refit_lambda in `flm_test`.
* Added Lee et al. (2020) test.
* Extended applications.
* Fixes in documentation.

# goffda 0.1.0

* Support forthcoming `Rcpp`'s STRICT_R_HEADERS.
* Update bibliographical details for main reference.
* Rename `r_gof2019_flmfr` to `r_gof2021_flmfr`.

# goffda 0.1.1

* Fix in HTML version of manual.

# goffda 0.1.2

* Drop C++11 requirement to adhere to new CRAN policies.
* Drop `personList()` and `citEntry()`.
