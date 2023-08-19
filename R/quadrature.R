

#' @title Quadrature rules for the \pkg{goffda} package
#'
#' @description Quadrature rules for unidimensional and bidimensional
#' functions, as enhancements of \code{\link[fda.usc]{int.simpson}} from the
#' \code{\link[fda.usc]{fda.usc-package}}.
#'
#' @param fx a vector of length \code{length(t)} with the evaluation of
#' a univariate function at \code{t}.
#' @param fxy a matrix of size \code{c(length(s), length(t))} with the
#' evaluation of a bivariate function at the bivariate grid formed by
#' \code{s} and \code{t}.
#' @param s,t vectors with grid points where functions are valued.
#' @param int_rule quadrature rule for approximating the definite
#' unidimensional integral: trapezoidal rule (\code{int_rule = "trapezoid"})
#' and extended Simpson rule (\code{int_rule = "Simpson"}) are available.
#' Defaults to \code{"trapezoid"}.
#' @param equispaced flag to indicate if \code{X_fdata$data} is valued in
#' an equispaced grid or not. Defaults to \code{FALSE}.
#' @param equispaced_x,equispaced_y flags to indicate if \code{fxy} is valued
#' in equispaced grids or not, at each one of dimensions. Both default to
#' \code{FALSE}.
#' @param verbose flag to show information about the procedures. Defaults
#' to \code{FALSE}.
#' @return
#' \itemize{
#'   \item{\code{w_integral1D}: a vector of length \code{t} with the weights
#'   required for the quadrature rule \code{int_rule}.}
#'   \item{\code{integral1D}: a scalar with the approximation of the univariate
#'   integral.}
#'   \item{\code{integral2D}: a scalar with the approximation of the bivariate
#'   integral.}
#' }
#' @examples
#' ## Numerical integral of 1-D functions
#'
#' # Equispaced grid points
#' t1 <- seq(0, 1, l = 201)
#' t2 <- seq(pi / 4, 3 * pi / 2, l = 201)
#' fx1 <- 2 * (t1^3) - t1^2 + 5 * t1 - 2 # True integral is equal to 2/3
#' fx2 <- sin(sqrt(t2)) # True integral is equal to 3.673555
#' int_fx1_trap <- integral1D(fx = fx1, t = t1, int_rule = "trapezoid",
#'                            equispaced = TRUE)
#' int_fx1_Simp <- integral1D(fx = fx1, t = t1, int_rule = "Simpson",
#'                            equispaced = TRUE)
#'
#' int_fx2_trap <- integral1D(fx = fx2, t = t2, int_rule = "trapezoid",
#'                            equispaced = TRUE)
#' int_fx2_Simp <- integral1D(fx = fx2, t = t2, int_rule = "Simpson",
#'                            equispaced = TRUE)
#'
#' # Check if the true integrals is approximated properly
#' abs(int_fx1_trap - 2/3) / (2/3)
#' abs(int_fx1_Simp - 2/3) / (2/3)
#' abs(int_fx2_trap - 3.673555) / 3.673555
#' abs(int_fx2_Simp - 3.673555) / 3.673555
#'
#' # Non equispaced grid points
#' t <- c(seq(0, 0.3, l = 50), seq(0.31, 0.6, l = 150),
#'        seq(0.61, 1, l = 100))
#' fx <- 2 * (t^3) - t^2 + 5 * t - 2
#' int_fx_trap <- integral1D(fx = fx, t = t, int_rule = "trapezoid",
#'                           equispaced = FALSE)
#' int_fx_Simp <- integral1D(fx = fx, t = t, int_rule = "Simpson",
#'                           equispaced = FALSE)
#'
#' # Check if the true integral is approximated properly
#' abs(int_fx_trap - 2/3) / (2/3)
#' abs(int_fx_Simp - 2/3) / (2/3)
#'
#' ## Numerical integral of 2-dimensional functions
#'
#' # Equispaced grid points
#' s <- seq(0, 2, l = 101)
#' t <- seq(1, 5, l = 151)
#' fxy <- outer(s^2, t^3) # True integral is equal to 416
#' int_fxy_trap <- integral2D(fxy = fxy, s = s, t = t, int_rule = "trapezoid",
#'                            equispaced_x = TRUE, equispaced_y = TRUE)
#' int_fxy_Simp <- integral2D(fxy = fxy, s = s, t = t, int_rule = "Simpson",
#'                            equispaced_x = TRUE, equispaced_y = TRUE)
#'
#' # Check if the true integral is approximated properly
#' abs(int_fxy_trap - 416) / 416
#' abs(int_fxy_Simp - 416) / 416
#'
#' # Non equispaced grid points
#' s <- c(seq(0, 0.3, l = 150), seq(0.31, 1.6, l = 100),
#'        seq(1.61, 2, l = 250))
#' t <- c(seq(1, 2.6, l = 170), seq(2.61, 4, l = 100),
#'        seq(4.01, 5, l = 140))
#' fxy <- outer(s^2, t^3)
#'
#' int_fxy_trap <- integral2D(fxy = fxy, s = s, t = t, int_rule = "trapezoid",
#'                            equispaced_x = FALSE, equispaced_y = FALSE)
#'
#' # Check if the true integral is approximated properly
#' abs(int_fxy_trap - 416) / 416
#' @author Code iterated by Javier Álvarez-Liébana and Eduardo García-Portugués
#' from the \code{\link[fda.usc]{fda.usc-package}} originals.
#' @references
#' Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery B. P. (1997).
#' Numerical Recipes in Fortran 77: the art of scientific computing (2nd ed).
#' @keywords internal
#' @name quadrature


#' @rdname quadrature
#' @export
integral1D <- function(fx, t, int_rule = "trapezoid", equispaced = FALSE,
                       verbose = FALSE) {

  # Check if grid points has a properly length
  if (length(fx) != length(t)) {

    stop(paste("The number of evaluations of the function must be the same",
               "as the number of grid points"))

  }

  # Check if NA
  if (any(is.na(fx)) && verbose) {

    warning("NA values have been found")

  }

  # Vector of weights
  weights <- w_integral1D(t = t, int_rule = int_rule,
                          equispaced = equispaced, verbose = verbose)

  # An interpolation (in an equispaced and finer grid) is required when
  # Simpson's rule is used and equispaced = FALSE (either when number of grid
  # points is less than 7)
  if (length(weights) != length(t)) {

    fx <- spline(x = t, y = fx, n = length(weights))[["y"]]

  }

  # Numerical integral as the weighted sum of the evaluated function in
  # the grid points
  return(sum(weights * fx))

}


#' @rdname quadrature
#' @export
integral2D <- function(fxy, s, t, int_rule = "trapezoid", equispaced_x = FALSE,
                       equispaced_y = FALSE, verbose = FALSE) {

  # Check if fxy represents a surface given by a 2D matrix
  if (!is.matrix(fxy)) {

    stop("Surface must be provided by a matrix")

  }

  # Only trapezoidal rule has been implemented for bidimensional functions with
  # non equispaced grids
  if (int_rule == "Simpson" && (equispaced_x + equispaced_y) != 2) {

    stop(paste("Only trapezoidal rule (int_rule = \"trapezoid\") is allowed",
                "for 2D functions and non equispaced grids"))

  }

  # Check if grid points has a properly length
  lx <- length(s)
  ly <- length(t)
  dim_fxy <- dim(fxy)
  if (dim_fxy[1] != lx || dim_fxy[2] != ly) {

    stop(paste("The number of evaluations of the function must",
               "be the same as the number of grid points"))

  }

  # Allow for integral1D as particular case
  if (lx == 1) {

    return(integral1D(fx = fxy[1, ], t = t, int_rule = int_rule,
                      equispaced = equispaced_y, verbose = verbose))

  } else if (ly == 1) {

    return(integral1D(fx = fxy[, 1], t = s, int_rule = int_rule,
                      equispaced = equispaced_x, verbose = verbose))

  } else {

    # Vector of weights for each of the grid intervals
    weights_x <- w_integral1D(t = s, int_rule = int_rule,
                              equispaced = equispaced_x, verbose = verbose)
    weights_y <- w_integral1D(t = t, int_rule = int_rule,
                              equispaced = equispaced_y, verbose = verbose)

    # Numerical integral as a double weighted sum of function valued in the grid
    # points for each one of the dimensions
    matrix_weights <- matrix(rep(weights_x, length(weights_y)),
                             nrow = length(weights_x)) *
      t(matrix(rep(weights_y, length(weights_x)), nrow = length(weights_y)))

    # Numerical integral as the weighted sum of the evaluated function in
    # the grid points
    return(sum(matrix_weights * fxy))

  }

}


#' @rdname quadrature
#' @export
w_integral1D <- function(t, int_rule = "trapezoid", equispaced = FALSE,
                         verbose = FALSE) {

  # Check if int_rule has been implemented
  if (!(int_rule %in% c("trapezoid", "Simpson"))) {

    stop("Only trapezoidal and extended Simpson's rules are allowed")

  }

  # Check if a minimum number of grid points have been introduced, but
  # allow length one to account for "scalar integration" = multiplication
  lx <- length(t)
  if (lx == 1) {

    return(1)

  }
  if (int_rule == "trapezoid" && lx < 2) {

    stop(paste("Trapezoidal rule requires at least 2 grid points,",
               "but the length of t is", lx))

  } else if (int_rule == "Simpson" && lx < 7) {

    stop(paste("The extended Simpson's rule requires at least 7 grid points,",
               "but the length of t is", lx))

  }

  # Simpson's rule requires equispaced grid points
  if (!equispaced && int_rule == "Simpson") {

    t <- seq(t[1], t[lx], l = lx) # We consider an equispaced grid
    lx <- length(t)
    if (verbose) {

      message(paste("Non equispaced grid points. An interpolation is",
                    "performed for applying the extended Simpsons's rule"))

    }

  }

  # Discretization step
  h <- (t[lx] - t[1]) / (lx - 1)

  # Weights for Simpson's rule (equispace grid has been forced)
  if (int_rule == "Simpson") {

    weights <- h * c(3 / 8, 7 / 6, 23 / 24, rep(1, lx - 6),
                     23 / 24, 7 / 6, 3 / 8)

  } else {

    # Weights for trapezoidal rule depending on value in equispaced
    weights <- switch(equispaced + 1,
                      c(0.5 * (t[2] - t[1]), t[3:lx] - t[3:lx - 1],
                        0.5 * (t[lx] - t[lx - 1])),
                      h * c(1 / 2, rep(1, lx - 2), 1 / 2))

  }

  # Output
  return(weights)

}
