# CLAUDE.md

Guidance for Claude Code sessions working on this repository.

## Package overview

**goffda** is a CRAN package implementing goodness-of-fit tests for
functional linear models with functional/scalar response and
functional/scalar predictor. It is the software companion of
García-Portugués, Álvarez-Liébana, Álvarez-Pérez and González-Manteiga
(2021), *Scandinavian Journal of Statistics* 48(2):502–528,
[doi:10.1111/sjos.12486](https://doi.org/10.1111/sjos.12486).

Current development cycle (branch
`claude/prepare-r-journal-submission-gpLmM`) is preparing the package
for an **R Journal submission** with version `0.2.0`.

## Layout

| Path | Purpose |
|---|---|
| `R/fpc.R` | FPC constructor + S3 methods (`print`, `summary`, `plot`, `predict`) at the bottom |
| `R/cv_glmnet.R` | Regularised linear fit + S3 methods (`print`, `coef`, `predict`, `plot`) |
| `R/flm_est.R` | FLM estimator + S3 methods (`print`, `summary`, `coef`, `fitted`, `residuals`, `predict`, `plot`); `predict` covers all four functional/scalar branches |
| `R/flm_test.R` | Goodness-of-fit test + helper `plot_flm_test_internal` + S3 methods (`summary`, `plot`); inherits `print` from `htest` |
| `R/auxiliary.R`, `R/quadrature.R`, `R/scenarios-fr.R`, `R/data.R` | Supporting code, simulators, datasets |
| `src/flm_stat.cpp` | Test statistic in C++ via Rcpp/RcppArmadillo |
| `man/*.Rd` | Roxygen-generated docs (do not edit by hand; run `devtools::document()`) |
| `tests/testthat/` | Unit tests, one file per source file |
| `data/`, `data-raw/`, `application/` | Datasets and replication scripts (excluded from build via `.Rbuildignore`) |
| `.claude/plans/` | Design/plan documents (see `s3-classes.md` for the v0.2.0 work) |

## Conventions

### S3 organisation
- Methods live **at the end of their constructor's file**, not in
  separate `R/s3-*.R` files (tidyverse / r-lib convention).
- All methods are documented under a shared roxygen page named
  `<constructor>-S3` via `@rdname`.
- `print.summary.X` shares the same page as the rest of the methods
  for class `X` so `R CMD check` does not flag missing docs.

### Style
- Two-space indent, snake_case for functions, blank lines between
  logical blocks (matches the existing code base).
- Roxygen with `RoxygenNote: 7.2.3` and `Roxygen: list(old_usage =
  TRUE)`. Do not change those without intent — the diff is noisy.
- Default to no comments; only add them when the *why* is non-obvious.
- Keep S3 method behaviour additive: returning lists with new class
  attributes must not break existing `$` / `[[` access patterns
  (there are 27+ internal accesses in `flm_est` and many more in
  `flm_test` that rely on this).

### Deprecation policy
- `plot_dens` and `plot_proc` of `flm_test()` are deprecated since
  0.2.0. They default to `FALSE` and emit `warning()` if `TRUE`.
  Removal scheduled for 0.3.0. Do **not** introduce new uses.

## Workflow

After modifying `R/*.R`, regenerate documentation and NAMESPACE so
roxygen and the manual files stay in sync:

```r
devtools::document()         # regenerate NAMESPACE + man/*.Rd
devtools::test()             # run testthat suite
devtools::check(args = "--as-cran")
```

CI runs the same `R CMD check` matrix
(`macos-latest`, `windows-latest`,
`ubuntu-latest` with `release`, `devel`, `oldrel-1`) via
`.github/workflows/R-CMD-check.yaml`.

## Roadmap (towards R Journal submission)

Done in 0.2.0 (commit `bd1fc8e` on the current branch):
- S3 classes for `fpc`, `cv_glmnet`, `flm_est`, `flm_test`.
- `predict.flm_est` for all four predictor/response combinations.
- Refactor of inline plotting in `flm_test` into `plot.flm_test`.
- testthat suite with one file per source file.
- README showcase of the new API; `BugReports` now points to
  `/issues`.

Open work (in roughly this order):
1. Run `devtools::document()` locally and commit any roxygen-generated
   diff in `man/*.Rd` and `NAMESPACE` (the manual edits in this
   commit are roxygen-equivalent but a clean regeneration is the
   reference).
2. Run `devtools::check(args = "--as-cran")` and resolve any
   NOTE/WARNING.
3. Test coverage with `covr` and a `test-coverage.yaml` workflow plus
   a Codecov badge in the README.
4. `pkgdown` site published via GitHub Pages.
5. Restore the `lint.yaml` workflow that was historically present in
   the commit history but is missing from the current
   `.github/workflows/`.
6. Manuscript using `rticles::rjournal_article`, placed under
   `paper/` with `^paper$` added to `.Rbuildignore`.
7. Vignette covering the FLM workflow (estimate → test → predict),
   sized so it can be reused as the JSS-style example in the paper.

## Branch policy

Develop on `claude/prepare-r-journal-submission-gpLmM`. Push to that
same branch (`git push -u origin
claude/prepare-r-journal-submission-gpLmM`). Do not push to `master`
without explicit user approval.

## Resources

- `.claude/plans/s3-classes.md` — full design doc for the v0.2.0 S3
  refactor (constructors, methods, file-by-file changes, test plan,
  end-to-end verification checklist).
- `NEWS.md` — user-facing change log.
- `cran-comments.md` — to be filled in for the next CRAN submission.
