

#' @title Computation of functional principal components
#'
#' @description Computation of Functional Principal Components (FPC) for
#' equispaced and non equispaced functional data.
#'
#' @inheritParams fda.usc_efic
#' @inheritParams quadrature
#' @param centered flag to indicate if \code{X_fdata} is centered or not.
#' Defaults to \code{FALSE}.
#' @param n_fpc number of FPC to be computed. If \code{n_fpc > n}, \code{n_fpc}
#' is set to \code{n}. Defaults to \code{3}.
#' @param verbose whether to show or not information about the \code{fpc}
#' procedure. Defaults to \code{FALSE}.
#' @return An \code{"fpc"} object containing the following elements:
#' \item{d}{standard deviations of the FPC (i.e., square roots of
#' eigenvalues of the empirical autocovariance estimator).}
#' \item{rotation}{orthonormal eigenfunctions (loadings or functional
#' principal components), as an \code{\link[fda.usc]{fdata}} class object.}
#' \item{scores}{rotated samples: inner products.
#' between \code{X_fdata} and eigenfunctions in\cr \code{rotation}.}
#' \item{l}{vector of index of FPC.}
#' \item{equispaced}{\code{equispaced} flag.}
#' @details The FPC are obtained by performing the single value decomposition
#' \deqn{\mathbf{X} \mathbf{W}^{1/2} =
#' \mathbf{U} \mathbf{D} (\mathbf{V}' \mathbf{W}^{1/2})}{
#' X W^{1/2} = U D (V' W^{1/2})}
#' where \eqn{\mathbf{X}}{X} is the matrix of discretized functional data,
#' \eqn{\mathbf{W}}{W} is a diagonal matrix of weights (computed by
#' \code{w_integral1D} according to \code{int_rule}), \eqn{\mathbf{D}}{D}
#' is the diagonal matrix with singular values (standard deviations of FPC),
#' \eqn{\mathbf{U}}{U} is a matrix whose columns contain the left singular
#' vectors, and \eqn{\mathbf{V}}{V} is a matrix whose columns contain the
#' right singular vectors (FPC).
#' @examples
#' ## Computing FPC for equispaced data
#'
#' # Sample data
#' X_fdata1 <- r_ou(n = 200, t = seq(2, 4, l = 201))
#'
#' # FPC with trapezoid rule
#' X_fpc1 <- fpc(X_fdata = X_fdata1, n_fpc = 50, equispaced = TRUE,
#'               int_rule = "trapezoid")
#'
#' # FPC with Simpsons's rule
#' X_fpc2 <- fpc(X_fdata = X_fdata1, n_fpc = 50, equispaced = TRUE,
#'                int_rule = "Simpson")
#'
#' # Check if FPC are orthonormal
#' norms1 <- rep(0, length(X_fpc1$l))
#' for (i in X_fpc1$l) {
#'
#'   norms1[i] <- integral1D(fx = X_fpc1$rotation$data[i, ]^2,
#'                           t = X_fdata1$argvals)
#'
#' }
#' norms2 <- rep(0, length(X_fpc2$l))
#' for (i in X_fpc2$l) {
#'
#'   norms2[i] <- integral1D(fx = X_fpc2$rotation$data[i, ]^2,
#'                           t = X_fdata1$argvals)
#'
#' }
#'
#' ## Computing FPC for non equispaced data
#'
#' # Sample data
#' X_fdata2 <- r_ou(n = 200, t = c(seq(0, 0.5, l = 201), seq(0.51, 1, l = 301)))
#'
#' # FPC with trapezoid rule
#' X_fpc3 <- fpc(X_fdata = X_fdata2, n_fpc = 5, int_rule = "trapezoid",
#'               equispaced = FALSE)
#'
#' # Check if FPC are orthonormal
#' norms3 <- rep(0, length(X_fpc3$l))
#' for (i in X_fpc3$l) {
#'
#'   norms3[i] <- integral1D(fx = X_fpc3$rotation$data[i, ]^2,
#'                           t = X_fdata2$argvals)
#'
#' }
#' \donttest{
#' ## Efficiency comparisons
#'
#' # fpc() vs. fda.usc::fdata2pc()
#' data(phoneme, package = "fda.usc")
#' mlearn <- phoneme$learn[1:10, ]
#' res1 <- fda.usc::fdata2pc(mlearn, ncomp = 3)
#' res2 <- fpc(X_fdata = mlearn, n_fpc = 3)
#' plot(res1$x[, 1:3], col = 1)
#' points(res2$scores, col = 2)
#'
#' microbenchmark::microbenchmark(fda.usc::fdata2pc(mlearn, ncomp = 3),
#'                                fpc(X_fdata = mlearn, n_fpc = 3), times = 1e3,
#'                                control = list(warmup = 20))
#' }
#' @author Javier Álvarez-Liébana and Gonzalo Álvarez-Pérez.
#' @references
#' Jolliffe, I. T. (2002). Principal Component Analysis. Springer-Verlag,
#' New York.
#' @export
fpc <- function(X_fdata, n_fpc = 3, centered = FALSE, int_rule = "trapezoid",
                equispaced = FALSE, verbose = FALSE) {

  # Check fdata
  stopifnot(fda.usc::is.fdata(X_fdata))

  # Basic info of X_fdata
  t <- X_fdata[["argvals"]]
  n <- dim(X_fdata[["data"]])[1]
  lx <- ncol(X_fdata[["data"]])

  # Check n_fpc
  if (n_fpc < 1) {

    stop("Number of components (n_fpc) must be greater or equal than one")

  }

  # Centering X if not centered
  if (!centered) {

    X_fdata <- fdata_cen(X_fdata = X_fdata)

  }

  # Check if a vector of grid points is provided
  if (lx == 1) {

    if (verbose) {

      message("Number of grid points is equal to one, X_fdata is scalar:",
              "data is just centered and n_fpc is forced to be one")

    }

    # Output as "fpc" object
    out <- list("d" = 1,
                "rotation" = fda.usc::fdata(mdata = cbind(1), argvals = t,
                                            rangeval = X_fdata[["rangeval"]]),
                "scores" = X_fdata[["data"]], "l" = 1, "equispaced" = TRUE)
    class(out) <- "fpc"
    return(out)

  } else {

    # Check if n_fpc is less than maximum
    if (n_fpc > min(n, lx)) {

      n_fpc <- min(n, lx)
      if (verbose) {

        message("Number of components n_fpc must be smaller than sample size",
                " and smaller than the number of grid points since an SVD",
                " is implemented: n_fpc = ", min(n, lx), " is enforced")

      }

    }

    # Components indexes
    l <- 1:n_fpc

    # weights is a diagonal matrix storing, as diagonal entries, the
    # square roots of weights provided by w_integral1D
    weights <- sqrt(w_integral1D(t = t, int_rule = int_rule,
                                 equispaced = equispaced))
    mat_weights <- diag(weights)

    # FPC are computed by performing a Singular Value Decomposition
    # X * weights^{1/2} = U * D * (V' * weights^{1/2}), more efficient than
    # solving the eigenequation as usual
    eigenres <- svd(X_fdata[["data"]] %*% mat_weights)

    # Functional Principal Components: eigenfunctions must be multiplied by the
    # inverse of the square root of weights to obtain orthonormal eigenfunctions
    v <- (1 / weights) * eigenres$v
    d <- eigenres$d # Square root of eigenvalues (standard deviation)

    # Data projected into FPC, recycling the weights already computed: weights^2
    # appears since eigenfunctions are divided by weights
    scores <- t((weights^2) * t(X_fdata[["data"]])) %*% v[, l]

    # Rotation as fdata
    rotation <- fda.usc::fdata(mdata = t(v[, l, drop = FALSE]), argvals = t,
                               rangeval = X_fdata[["rangeval"]])

    # Output as "fpc" object
    out <- list("d" = d[l], "rotation" = rotation, "scores" = scores, "l" = l,
                "equispaced" = equispaced)
    class(out) <- "fpc"
    return(out)

  }

}


#' @title Utilities for functional principal components
#'
#' @description Computation of coefficients and reconstructions based on
#' Functional Principal Components (FPC). The function \code{fpc_coefs} allows
#' to project a functional data sample into a basis of FPC; the reconstruction
#' of the sample from its projections and the FPC is done with
#' \code{fpc_to_fdata}. The functions \code{beta_fpc_coefs} and
#' \code{fpc_to_beta} do analogous operations but for the
#' \link[=flm_est]{bivariate kernel} \eqn{\beta} and the tensor product
#' of two FPC bases.
#'
#' @inheritParams fpc
#' @inheritParams quadrature
#' @param beta a matrix containing the bivariate kernel \eqn{\beta} evaluated
#' on a grid. Must be of size \code{c(length(X_fpc$rotation$argvals),
#' length(Y_fpc$rotation$argvals))}.
#' @param X_fpc,Y_fpc \code{"fpc"} objects as resulted from calling
#' \code{\link{fpc}}.
#' @param ind_X_fpc,ind_Y_fpc vectors giving the FPC indexes for whom the
#' coefficients are computed. Their lengths must be smaller than the number of
#' FPC in \code{X_fpc} and \code{Y_fpc}, respectively. Default to \code{1:3}.
#' @param coefs a vector of coefficients to combine linearly the FPC. Its
#' length must be smaller than the number of FPC in \code{X_fpc}.
#' @param beta_coefs a matrix of coefficients to combine linearly the tensor
#' products of FPC. Its size must be smaller than the number of FPC in
#' \code{X_fpc} and \code{Y_fpc}.
#' @param ind_coefs,ind_coefs_X,ind_coefs_Y indexes of FPC to associate to the
#' provided coefficients. By default, from the first FPC to the sizes of
#' \code{coefs} or \code{beta_coefs}.
#' @return
#' \item{fpc_coefs}{a vector of the same length as \code{coefs} containing
#' the coefficients of \code{X_fdata} in the FPC of \code{X_fpc}.}
#' \item{beta_fpc_coefs}{a matrix of the same size as \code{beta_coefs}
#' containing the coefficients of \eqn{\beta} in the tensor product of
#' the FPC in \code{X_fpc} and \code{Y_fpc}.}
#' \item{fpc_to_fdata}{an \code{\link[fda.usc]{fdata}} object of the same
#' type as \code{X_fpc$rotation}.}
#' \item{fpc_to_beta}{a matrix with the reconstructed kernel and size\cr
#' \code{c(length(X_fpc$rotation$argvals), length(Y_fpc$rotation$argvals))}.}
#' @examples
#' ## Compute FPC coefficients and reconstruct data
#'
#' # Sample data
#' X_fdata <- r_ou(n = 200, t = seq(2, 4, l = 201))
#'
#' # Compute FPC
#' X_fpc <- fpc(X_fdata = X_fdata, n_fpc = 50)
#'
#' # FPC coefficients are the same if the data is centered
#' fpc_coefs(X_fdata = fdata_cen(X_fdata), X_fpc = X_fpc)[1:4, ]
#' X_fpc$scores[1:4, 1:3]
#'
#' # Reconstruct the first two curves for an increasing number of FPC
#' plot(X_fdata[1:2, ], col = 1)
#' n_fpc <- c(2, 5, 10, 25, 50)
#' for (j in 1:5) {
#'   lines(fpc_to_fdata(X_fpc = X_fpc,
#'                      coefs = X_fpc$scores[, 1:n_fpc[j]])[1:2, ], col = j + 1)
#' }
#'
#' ## Project and reconstruct beta
#'
#' # Surface
#' beta_fun <- function(s, t) sin(6 * pi * s) + cos(6 * pi * t)
#' s <- seq(0, 1, l = 101)
#' t <- seq(0, 1, l = 201)
#' beta_surf <- outer(s, t, FUN = beta_fun)
#'
#' # Functional data as zero-mean Gaussian process with exponential variogram
#' X_fdata <- fda.usc::rproc2fdata(n = 100, t = s, sigma = "vexponential",
#'                                 list = list(scale = 2.5))
#' Y_fdata <- flm_term(X_fdata = X_fdata, beta = beta_surf, t = t) +
#'   r_ou(n = 100, t = t, sigma = sqrt(0.075) * 2)
#'
#' # FPC
#' X_fpc <- fpc(X_fdata = X_fdata, n_fpc = 50)
#' Y_fpc <- fpc(X_fdata = Y_fdata, n_fpc = 50)
#'
#' # Coefficients
#' beta_coefs <- beta_fpc_coefs(beta = beta_surf, X_fpc = X_fpc, Y_fpc = Y_fpc,
#'                              ind_X_fpc = 1:50, ind_Y_fpc = 1:50)
#'
#' # Reconstruction
#' beta_surf1 <- fpc_to_beta(beta_coefs = beta_coefs[1:2, 1:5],
#'                           X_fpc = X_fpc, Y_fpc = Y_fpc)
#' beta_surf2 <- fpc_to_beta(beta_coefs = beta_coefs[1:15, 1:10],
#'                           X_fpc = X_fpc, Y_fpc = Y_fpc)
#' beta_surf3 <- fpc_to_beta(beta_coefs = beta_coefs[1:50, 1:50],
#'                           X_fpc = X_fpc, Y_fpc = Y_fpc)
#'
#' # Show reconstructions
#' old_par <- par(mfrow = c(2, 2))
#' col <- viridisLite::viridis(20)
#' image(s, t, beta_surf, col = col, zlim = c(-2.5, 2.5), main = "Original")
#' image(s, t, beta_surf1, col = col, zlim = c(-2.5, 2.5), main = "2 x 5")
#' image(s, t, beta_surf2, col = col, zlim = c(-2.5, 2.5), main = "15 x 10")
#' image(s, t, beta_surf3, col = col, zlim = c(-2.5, 2.5), main = "50 x 50")
#' par(old_par)
#' @author Eduardo García-Portugués.
#' @references
#' Jolliffe, I. T. (2002). Principal Component Analysis. Springer-Verlag,
#' New York.
#' @name fpc_utils


#' @rdname fpc_utils
#' @export
fpc_coefs <- function(X_fdata, X_fpc, ind_X_fpc = 1:3, int_rule = "trapezoid") {

  # Sample size
  n <- length(X_fdata)

  # Number of FPC for the coefs
  n_fpc <- length(X_fpc[["rotation"]])
  if (max(ind_X_fpc) > n_fpc) {

    stop("ind_X_fpc must contain up to", n_fpc, "components")

  }

  # Check X_fpc
  if (!(class(X_fpc) %in% "fpc")) {

    stop(paste("X_fpc must be an \"fpc\" class object",
               "as returned by the fpc function"))

  }

  # Check if the lengths of argvals coincide
  l_argvals <- length(X_fdata[["argvals"]])
  if (l_argvals != length(X_fpc[["rotation"]][["argvals"]])) {

    stop("The argvals of X_fdata and X_fpc must coincide")

  }

  # Coefficients matrix
  coefs <-
    X_fdata[["data"]][rep(1:n, times = length(ind_X_fpc)), , drop = FALSE] *
    X_fpc[["rotation"]][["data"]][rep(ind_X_fpc, each = n), , drop = FALSE]
  if (l_argvals > 1) {

    coefs <- matrix(apply(coefs, 1, FUN = integral1D,
                          t = X_fdata[["argvals"]], int_rule = int_rule,
                          equispaced = X_fpc[["equispaced"]]), nrow = n)

  }
  return(coefs)

}


#' @rdname fpc_utils
#' @export
beta_fpc_coefs <- function(beta, X_fpc, Y_fpc, ind_X_fpc = 1:3, ind_Y_fpc = 1:3,
                           int_rule = "trapezoid") {

  # Check X_fpc and Y_fpc
  if (!(class(X_fpc) %in% "fpc")) {

    stop(paste("X_fpc must be an \"fpc\" class object",
               "as returned by the fpc function"))

  } else if (!(class(Y_fpc) %in% "fpc")) {

    stop(paste("Y_fpc must be an \"fpc\" class object",
               "as returned by the fpc function"))

  }

  # Check if beta is a matrix
  if (!is.matrix(beta)) {

    stop("beta must be a matrix")

  }

  # Check the dimension of beta with the argvals of X_fdata and Y_fdata
  s <- X_fpc[["rotation"]][["argvals"]]
  t <- Y_fpc[["rotation"]][["argvals"]]
  if (nrow(beta) != length(s) | ncol(beta) != length(t)) {

    stop(paste("The dimensions of beta do not coincide with the argvals",
               "of X_fpc and Y_fpc"))

  }

  # Compute coefficients
  equispaced_x <- X_fpc[["equispaced"]]
  equispaced_y <- Y_fpc[["equispaced"]]
  beta_coefs <- matrix(nrow = length(ind_X_fpc), ncol = length(ind_Y_fpc))
  for (i in seq_along(ind_X_fpc)) {
    for (j in seq_along(ind_Y_fpc)) {

      # Compute ij-th coefficient
      psi_ij <- tcrossprod(x = X_fpc[["rotation"]][["data"]][ind_X_fpc[i], ],
                           y = Y_fpc[["rotation"]][["data"]][ind_Y_fpc[j], ])
      beta_coefs[i, j] <- integral2D(fxy = beta * psi_ij, s = s, t = t,
                                     int_rule = int_rule,
                                     equispaced_x = equispaced_x,
                                     equispaced_y = equispaced_y)

    }
  }

  # Coefficients matrix
  return(beta_coefs)

}


#' @rdname fpc_utils
#' @export
fpc_to_fdata <- function(coefs, X_fpc, ind_coefs = 1:ncol(coefs)) {

  # coefs as matrix
  if (!is.matrix(coefs)) {

    coefs <- cbind(coefs)

  }
  n_coefs <- ncol(coefs)

  # Check size of coefs and ind_coefs
  n_fpc <- nrow(X_fpc[["rotation"]][["data"]])
  if (n_coefs > n_fpc) {

    stop(paste0("The number of coefs (", n_coefs, ") must be smaller or equal",
               " than the number of FPC (", n_fpc, ")"))

  } else if (length(ind_coefs) != n_coefs) {

    stop("The lengths of coefs and ind_coefs do not match")

  }

  # Construct fdata
  X_fpc[["rotation"]][["data"]] <-
    coefs %*% X_fpc[["rotation"]][["data"]][ind_coefs, , drop = FALSE]
  return(X_fpc[["rotation"]])

}


#' @rdname fpc_utils
#' @export
fpc_to_beta <- function(beta_coefs, X_fpc, Y_fpc,
                        ind_coefs_X = 1:nrow(beta_coefs),
                        ind_coefs_Y = 1:ncol(beta_coefs)) {

  # Check beta_coefs
  n_fpc_X <- nrow(X_fpc[["rotation"]][["data"]])
  n_fpc_Y <- nrow(Y_fpc[["rotation"]][["data"]])
  max_ind_coefs_X <- max(ind_coefs_X)
  max_ind_coefs_Y <- min(ind_coefs_Y)
  if (!is.matrix(beta_coefs)) {

    stop("beta_coefs must be a matrix")

  } else {

    if (nrow(beta_coefs) != length(ind_coefs_X) |
        ncol(beta_coefs) != length(ind_coefs_Y)) {

      stop(paste("The dimensions of beta_coefs do not match with",
                 "ind_coefs_X or ind_coefs_Y"))

    } else if (n_fpc_X < max_ind_coefs_X) {

      stop(paste0("The number of FPC in ind_coefs_X (", max_ind_coefs_X,
                  ") must be smaller or equal than the number of FPC in X (",
                  n_fpc_X, ")"))

    } else if (n_fpc_Y < max_ind_coefs_Y) {

      stop(paste0("The number of FPC in ind_coefs_Y (", max_ind_coefs_Y,
                  ") must be smaller or equal than the number of FPC in Y (",
                  n_fpc_Y, ")"))

    }

  }

  # Construct beta
  beta <- t(X_fpc[["rotation"]][["data"]][ind_coefs_X, , drop = FALSE]) %*%
    beta_coefs %*% Y_fpc[["rotation"]][["data"]][ind_coefs_Y, , drop = FALSE]
  return(beta)

}
