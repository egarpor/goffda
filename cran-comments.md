
# Update 0.0.6

## Test environments

* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## Test results

No errors. 

One PREPERROR on `Debian Linux, R-devel, GCC ASAN/UBSAN`:

> PREPERROR ERROR: dependency ‘glmnet’ is not available for package ‘goffda’

but `glmnet` is a massively used package likely to be available soon on such platform.

One note on `Fedora Linux, R-devel, clang, gfortran`:

> Possibly mis-spelled words in DESCRIPTION:
> García (20:6)
> González (20:59)
> Liébana (20:32)
> Manteiga (20:68)
> Portugués (20:13)
> Pérez (20:49)
> Álvarez (20:24, 20:41)

but words are properly spelled.

