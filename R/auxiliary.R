

#' @title Auxiliary functions for the \pkg{goffda} package
#'
#' @description Auxiliary functions required for the methods
#' implemented in the \pkg{goffda} package, as enhancements of the auxiliary
#' functions \code{\link[fda.usc]{fdata.cen}} and
#' \code{\link[fda.usc]{func.mean}} from the
#' \code{\link[fda.usc]{fda.usc-package}}.
#'
#' @param X_fdata sample of functional data as an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @param mean_X functional mean of \code{X_fdata}.
#' @param X_fdata1,X_fdata2 samples of functional
#' data as \code{\link[fda.usc]{fdata}} objects of lengths \eqn{n_1}
#' and \eqn{n_2}, respectively. Sample sizes can be different.
#' @inheritParams quadrature
#' @param as_matrix flag to indicate if \code{inprod_fdata} returns a matrix
#' or the vector of its lower triangular part in column-major order.
#' Defaults to \code{TRUE}.
#' @param verbose whether to show or not information about the
#' \code{inprod_fdata} procedure.
#' @details
#' \itemize{
#'   \item{\code{func_mean}: computes the functional mean of
#'   \code{X_fdata}.}
#'   \item{\code{fdata_cen}: centers the
#'   functional data \code{X_fdata}.}
#'   \item{\code{inprod_fdata(X_fdata1)}:  computes as a row vector the
#'   elements of the lower triangular part of the inner products matrix
#'   (\code{X_fdata} vs \code{X_fdata}). If \code{as_matrix = TRUE}, the
#'   matrix of inner products is given.}
#'   \item{\code{inprod_fdata(X_fdata1, X_fdata2)}: computes the matrix of
#'   inner products (\code{as_matrix = TRUE} is forced) between \code{X_fdata1}
#'   and \code{X_fdata2}.}
#' }
#' @examples
#' ## fdata_cen() vs fda.usc::fdata_cen()
#'
#' data(phoneme, package = "fda.usc")
#' mlearn <- phoneme$learn[1:10, ]
#' plot(fda.usc::fdata.cen(mlearn)$Xcen)
#' plot(fdata_cen(mlearn))
#'
#' ## inprod_fdata() vs fda.usc::inprod.fdata()
#'
#' # inprod_fdata between mlearn and mlearn: as a row vector
#'
#' A <- fda.usc::inprod.fdata(fdata1 = mlearn)
#' A[upper.tri(A, diag = TRUE)]
#' inprod_fdata(X_fdata1 = mlearn, int_rule = "trapezoid", as_matrix = FALSE)
#'
#' # inprod_fdata between mlearn and mlearn: as a matrix
#'
#' A <- fda.usc::inprod.fdata(fdata1 = mlearn)
#' A
#' inprod_fdata(X_fdata1 = mlearn, int_rule = "trapezoid", as_matrix = TRUE)
#'
#' # inprod_fdata between mlearn and mlearn2: as a matrix
#'
#' mlearn2 <- phoneme$learn[11:30, ]
#' A <- fda.usc::inprod.fdata(fdata1 = mlearn, fdata2 = mlearn2)
#' A
#' B <- inprod_fdata(X_fdata1 = mlearn, X_fdata2 = mlearn2,
#' int_rule = "trapezoid", as_matrix = TRUE)
#' B
#' \donttest{
#' ## Efficiency comparisons
#'
#' microbenchmark::microbenchmark(fda.usc::fdata.cen(mlearn), fdata_cen(mlearn),
#'                                times = 1e3, control = list(warmup = 20))
#'
#' microbenchmark::microbenchmark(fda.usc::inprod.fdata(fdata1 = mlearn),
#'                                inprod_fdata(X_fdata1 = mlearn,
#'                                as_matrix = FALSE), times = 1e3,
#'                                control = list(warmup = 20))
#' }
#' @author Code iterated by Eduardo García-Portugués, Gonzalo Álvarez-Pérez,
#' and Javier Álvarez-Liébana from the \code{\link[fda.usc]{fda.usc-package}}
#' originals.
#' @export
#' @keywords internal
#' @name fda.usc_efic


#' @rdname fda.usc_efic
#' @export
fdata_cen <- function(X_fdata, mean_X = func_mean(X_fdata)) {

  # Check fdata
  stopifnot(fda.usc::is.fdata(X_fdata))

  # Data centering employing an implicit column recycling
  X_fdata[["data"]] <- t(t(X_fdata[["data"]]) - drop(mean_X[["data"]]))
  return(X_fdata)

}


#' @rdname fda.usc_efic
#' @export
func_mean <- function(X_fdata) {

  # Check fdata
  stopifnot(fda.usc::is.fdata(X_fdata))

  # Functional mean
  X_fdata[["data"]] <- rbind(colMeans(X_fdata[["data"]], na.rm = TRUE))
  return(X_fdata)

}


#' @rdname fda.usc_efic
#' @export
inprod_fdata <- function(X_fdata1, X_fdata2 = NULL, int_rule = "trapezoid",
                         as_matrix = TRUE, verbose = FALSE) {

  # Check fdata
  stopifnot(fda.usc::is.fdata(X_fdata1))

  # Basic info of X_fdata1
  t1 <- X_fdata1[["argvals"]]
  n1 <- dim(X_fdata1[["data"]])[1]

  # Check if X_fdata1 is equispaced
  eps <- sqrt(.Machine[["double.eps"]])
  equispaced <- all(abs(diff(t1, differences = 2)) < eps)

  # Check if X_fdata2 is NULL: X_fdata1 vs X_fdata1, we compute the matrix of
  # inner products between X_fdata1 and itself
  if (is.null(X_fdata2)) {

    inprods <- numeric(n1 * (n1 + 1) / 2)
    inprods[1] <- integral1D(X_fdata1[["data"]][1, ]^2, t = t1,
                             int_rule = int_rule, equispaced = equispaced)
    k <- 2
    for (i in 2:n1) {

      # Inner product in L2 spaces by using integral1D
      inprods[k:(k + i - 1)] <- apply(t(t(X_fdata1[["data"]][1:i, ]) *
                                          X_fdata1[["data"]][i, ]), 1,
                                      FUN = "integral1D", t1, int_rule,
                                      equispaced)
      k <- k + i

    }

    # If as_matrix = TRUE, we construct the symmetric matrix
    if (!as_matrix) {

      return(inprods)

    } else {

      mat_inprods <- matrix(0, n1, n1)
      mat_inprods[1, 1] <- integral1D(X_fdata1[["data"]][1, ]^2, t = t1,
                                      int_rule = int_rule,
                                      equispaced = equispaced)
      for (i in 2:n1) {

        # Inner product in L2 spaces by using integral1D
        mat_inprods[i, 1:i] <- apply(t(t(X_fdata1[["data"]][1:i, ]) *
                                         X_fdata1[["data"]][i, ]), 1,
                                     FUN = "integral1D", t1, int_rule,
                                     equispaced)

      }

      # We complete the upper part of the matrix
      mat_inprods[upper.tri(mat_inprods, diag = FALSE)] <-
        mat_inprods[lower.tri(mat_inprods, diag = FALSE)]
      return(mat_inprods)

    }

  } else {

    # Check fdata
    stopifnot(fda.usc::is.fdata(X_fdata2))

    # as_matrix = TRUE is forced
    as_matrix <- TRUE
    if (verbose) {

      message("Two \"fdata\" objects have been introduced: ",
              "as_matrix is forced to be TRUE")

    }

    # Basic info of X_fdata2
    t2 <- X_fdata2[["argvals"]]
    n2 <- dim(X_fdata2[["data"]])[1]

    # Check if X_fdata1 and X_fdata2 are valued in the same grid
    if (length(t1) != length(t2)) {

      stop("Functional data 1 and 2 must be valued in the same grid")

    }

    # If X_fdata2 is available the whole (non symmetric) matrix must be computed
    mat_inprods <- matrix(0, n1, n2)
    for (i in 1:n1) {

      mat_inprods[i, ] <- apply(t(t(X_fdata2[["data"]]) *
                                     X_fdata1[["data"]][i, ]), 1,
                                 FUN = "integral1D", t1, int_rule, equispaced)

    }
    return(mat_inprods)

  }

}
