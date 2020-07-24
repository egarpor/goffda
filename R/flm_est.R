

#' @title Estimation of functional linear models
#'
#' @description Estimation of the linear operator relating a
#' functional predictor \eqn{X} with a functional response \eqn{Y} in the
#' linear model
#' \deqn{Y(t) = \int_a^b \beta(s, t) X(s)\,\mathrm{d}s + \varepsilon(t),}{
#' Y(t) = \int_a^b \beta(s, t) X(s) ds + \epsilon(t),}
#' where \eqn{X} is a random variable in the Hilbert space of
#' square-integrable functions in \eqn{[a, b]}, \eqn{L^2([a, b])},
#' \eqn{Y} and \eqn{\varepsilon}{\epsilon} are random variables
#' in \eqn{L^2([c, d])}, and \eqn{s \in [a, b]} and \eqn{t \in [c, d]}.
#'
#' The linear, Hilbert--Schmidt, integral operator is parametrized by
#' the bivariate kernel \eqn{\beta \in L^2([a, b]) \otimes
#' L^2([c, d])}. Its estimation is done through the truncated expansion
#' of \eqn{\beta} in the tensor product of the data-driven
#' bases of the Functional Principal Components (FPC) of
#' \eqn{X} and \eqn{Y}, and through the fitting of the resulting multivariate
#' linear model. The FPC basis for \eqn{X} is truncated in \eqn{p}
#' components, while the FPC basis for \eqn{Y} is truncated in \eqn{q}
#' components. Automatic selection of \eqn{p} and \eqn{q} is detailed below.
#'
#' The particular cases in which either \eqn{X} or \eqn{Y} are
#' \emph{constant} functions give either a scalar predictor or response.
#' The simple linear model arises if both \eqn{X} and \eqn{Y} are scalar,
#' for which \eqn{\beta} is a constant.
#'
#' @param X,Y samples of functional/scalar predictors and functional/scalar
#' response. Either \code{\link[fda.usc]{fdata}} objects (for functional
#' variables) or vectors of length \code{n} (for scalar variables).
#' @param est_method either \code{"fpcr"} (Functional Principal Components
#' Regression; FPCR), \code{"fpcr_l2"} (FPCR with ridge penalty),
#' \code{"fpcr_l1"} (FPCR with lasso penalty) or \code{"fpcr_l1s"}
#' (FPCR with lasso-selected FPC). If \code{X} is scalar, \code{flm_est}
#' only considers \code{"fpcr"} as estimation method. See details below.
#' Defaults to \code{"fpcr_l1s"}.
#' @param p,q index vectors indicating the specific FPC to be
#' considered for the truncated bases expansions of \code{X} and \code{Y},
#' respectively. If a single number for \code{p} is provided, then
#' \code{p <- 1:max(p)} internally (analogously for \code{q}) and the first
#' \code{max(p)} FPC are considered. If \code{NULL} (default), then a
#' data-driven selection of \code{p} and \code{q} is done. See details below.
#' @param thre_p,thre_q thresholds for the \emph{proportion} of variance
#' that is explained, \emph{at least}, by the first \eqn{p} and \eqn{q} FPC
#' of \code{X} and \code{Y}, respectively. These thresholds are employed
#' for an (initial) automatic selection of \eqn{p} and \eqn{q}.
#' Default to \code{0.99}. \code{thre_p} (\code{thre_q}) is ignored if
#' \code{p} (\code{q}) is provided.
#' @param lambda regularization parameter \eqn{\lambda} for the estimation
#' methods \code{"fpcr_l2"}, \code{"fpcr_l1"}, and \code{"fpcr_l1s"}. If
#' \code{NULL} (default), it is chosen with \code{\link[goffda]{cv_glmnet}}.
#' @param X_fpc,Y_fpc FPC decompositions of \code{X} and \code{Y}, as
#' returned by \code{\link[goffda]{fpc}}. Computed if not provided.
#' @param compute_residuals whether to compute the fitted values \code{Y_hat}
#' and its \code{Y_hat_scores}, and the \code{residuals} of the fit
#' and its \code{residuals_scores}. Defaults to \code{TRUE}.
#' @param centered flag to indicate if \code{X} and \code{Y} have been
#' centered or not, in order to avoid their recentering. Defaults to
#' \code{FALSE}.
#' @inheritParams fpc
#' @param cv_verbose flag to display information about the estimation procedure
#' (passed to \code{\link{cv_glmnet}}). Defaults to \code{FALSE}.
#' @param ... further parameters to be passed to \code{\link{cv_glmnet}}
#' (and then to \code{\link[glmnet]{cv.glmnet}}) such as \code{cv_1se},
#' \code{cv_nlambda} or \code{cv_parallel}, among others.
#' @return A list with the following entries:
#' \item{Beta_hat}{estimated \eqn{\beta}, a matrix with values
#' \eqn{\hat\beta (s, t)} evaluated at the grid points for \code{X}
#' and \code{Y}. Its size is \code{c(length(X$argvals), length(Y$argvals))}.}
#' \item{Beta_hat_scores}{the matrix of coefficients of \code{Beta_hat}
#' (resulting from projecting it into the tensor basis of \code{X_fpc} and
#' \code{Y_fpc}), with dimension \code{c(p_hat, q_thre)}.}
#' \item{H_hat}{hat matrix of the associated fitted multivariate
#' linear model, a matrix of size \code{c(n, n)}. \code{NULL} if
#' \code{est_method = "fpcr_l1"}, since lasso estimation does not provide it
#' explicitly.}
#' \item{p_thre}{index vector indicating the FPC of \code{X}
#' considered for estimating the model. Chosen by \code{thre_p} or equal
#' to \code{p} if given.}
#' \item{p_hat}{index vector of the FPC considered by the methods
#' \code{"fpcr_l1"} and \code{"fpcr_l1s"} methods after further selection
#' of the FPC considered in \code{p_thre}. For methods \code{"fpcr"} and
#' \code{"fpcr_l2"}, \code{p_hat} equals \code{p_thre}.}
#' \item{q_thre}{index vector indicating the FPC of \code{Y}
#' considered for estimating the model. Chosen by \code{thre_q} or equal
#' to \code{q} if given. Note that zeroing by lasso procedure only affects
#' in the rows.}
#' \item{est_method}{the estimation method employed.}
#' \item{Y_hat}{fitted values, either an \code{\link[fda.usc]{fdata}}
#' object or a vector, depending on \code{Y}.}
#' \item{Y_hat_scores}{the matrix of coefficients of \code{Y_hat}, with
#' dimension \code{c(n, q_thre)}.}
#' \item{residuals}{residuals of the fitted model, either an
#' \code{\link[fda.usc]{fdata}} object or a vector, depending on \code{Y}.}
#' \item{residuals_scores}{the matrix of coefficients of
#' \code{residuals}, with dimension \code{c(n, q_thre)}.}
#' \item{X_fpc, Y_fpc}{FPC of \code{X} and \code{Y}, as
#' returned by \code{\link{fpc}} with \code{n_fpc = n}.}
#' \item{lambda}{regularization parameter \eqn{\lambda} used for the
#' estimation methods \code{"fpcr_l2"}, \code{"fpcr_l1"}, and
#' \code{"fpcr_l1s"}.}
#' \item{cv}{cross-validation object returned by
#' \code{\link{cv_glmnet}}.}
#' @details
#' \code{flm_est} deals seamlessly with either functional or scalar inputs
#' for the predictor and response. In the case of scalar inputs, the
#' corresponding dimension-related arguments (\code{p}, \code{q},
#' \code{thre_p} or \code{thre_q}) will be ignored as in these cases either
#' \eqn{p = 1} or \eqn{q = 1}.
#'
#' The function translates the functional linear model into a multivariate
#' model with multivariate response and then estimates the
#' \eqn{p \times q}{p x q} matrix of coefficients of \eqn{\beta} in the
#' tensor basis of the FPC of \code{X} and \code{Y}. The following estimation
#' methods are implemented:
#' \itemize{
#'   \item \code{"fpcr"}: Functional Principal Components Regression (FPCR);
#'   see details in Ramsay and Silverman (2005).
#'   \item \code{"fpcr_l2"}: FPCR, with ridge penalty on the associated
#'   multivariate linear model.
#'   \item \code{"fpcr_l1"}: FPCR, with lasso penalty on the associated
#'   multivariate linear model.
#'   \item \code{"fpcr_l1s"}: FPCR, with FPC selected by lasso regression
#'   on the associated multivariate linear model.
#' }
#' The last three methods are explained in García-Portugués et al. (2019).
#'
#' The \eqn{p} FPC of \code{X} and \eqn{q} FPC of \code{Y} are determined
#' as follows:
#' \itemize{
#'   \item If \code{p = NULL}, then \code{p} is set as
#'   \code{p_thre <- 1:j_thre}, where \code{j_thre} is the \eqn{j}-th FPC of
#'   \code{X} for which the cumulated proportion of explained variance is
#'   greater than \code{thre_p}. If \code{p != NULL}, then \code{p_thre <- p}.
#'   \item If \code{q = NULL}, then the same procedure is followed with
#'   \code{thre_q}, resulting \code{q_thre}.
#' }
#' Once \code{p_thre} and \code{q_thre} have been obtained, the methods
#' \code{"fpcr_l1"} and \code{"fpcr_l1s"} perform a second selection
#' of the FPC that are effectively considered in the estimation of \eqn{\beta}.
#' This subset of FPC (of \code{p_thre}) is encoded in \code{p_hat}. No further
#' selection of FPC is done for the methods \code{"fpcr"} and \code{"fpcr_l2"}.
#'
#' The flag \code{compute_residuals} controls if \code{Y_hat},
#' \code{Y_hat_scores}, \code{residuals}, and \code{residuals_scores} are
#' computed. If \code{FALSE}, they are set to \code{NULL}. \code{Y_hat} equals
#' \eqn{\hat Y_i(t) = \int_a^b \hat \beta(s, t) X_i(s) \,\mathrm{d}s}{
#' \hat Y_i(t) = \int_a^b \hat \beta(s, t) X_i(s) ds} and \code{residuals}
#' stands for \eqn{\hat \varepsilon_i(t) = Y_i(t) - \hat Y_i(t)}{
#' \hat \epsilon_i(t) = Y_i(t) - \hat Y_i(t)}, both for
#' \eqn{i = 1, \ldots, n}. \code{Y_hat_scores} and\cr \code{residuals_scores}
#' are the \eqn{n\times q}{n x q} matrices of coefficients (or scores) of these
#' functions in the FPC of \code{Y}.
#'
#' Missing values on \code{X} and \code{Y} are automatically removed.
#' @examples
#' ## Quick example of functional response and functional predictor
#'
#' # Generate data
#' set.seed(12345)
#' n <- 50
#' X_fdata <- r_ou(n = n, t = seq(0, 1, l = 201), sigma = 2)
#' epsilon <- r_ou(n = n, t = seq(0, 1, l = 201), sigma = 0.5)
#' Y_fdata <- 2 * X_fdata + epsilon
#'
#' # Lasso-selection FPCR (p and q are estimated)
#' flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s")
#' \donttest{
#' ## Functional response and functional predictor
#'
#' # Generate data
#' set.seed(12345)
#' n <- 50
#' X_fdata <- r_ou(n = n, t = seq(0, 1, l = 201), sigma = 2)
#' epsilon <- r_ou(n = n, t = seq(0, 1, l = 201), sigma = 0.5)
#' Y_fdata <- 2 * X_fdata + epsilon
#'
#' # FPCR (p and q are estimated)
#' fpcr_1 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr")
#' fpcr_1$Beta_hat_scores
#' fpcr_1$p_thre
#' fpcr_1$q_thre
#'
#' # FPCR (p and q are provided)
#' fpcr_2 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr",
#'                   p = c(1, 5, 2, 7), q = 2:1)
#' fpcr_2$Beta_hat_scores
#' fpcr_2$p_thre
#' fpcr_2$q_thre
#'
#' # Ridge FPCR (p and q are estimated)
#' l2_1 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2")
#' l2_1$Beta_hat_scores
#' l2_1$p_hat
#'
#' # Ridge FPCR (p and q are provided)
#' l2_2 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2",
#'                 p = c(1, 5, 2, 7), q = 2:1)
#' l2_2$Beta_hat_scores
#' l2_2$p_hat
#'
#' # Lasso FPCR (p and q are estimated)
#' l1_1 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1")
#' l1_1$Beta_hat_scores
#' l1_1$p_thre
#' l1_1$p_hat
#'
#' # Lasso estimator (p and q are provided)
#' l1_2 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1",
#'                 p = c(1, 5, 2, 7), q = 2:1)
#' l1_2$Beta_hat_scores
#' l1_2$p_thre
#' l1_2$p_hat
#'
#' # Lasso-selection FPCR (p and q are estimated)
#' l1s_1 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s")
#' l1s_1$Beta_hat_scores
#' l1s_1$p_thre
#' l1s_1$p_hat
#'
#' # Lasso-selection FPCR (p and q are provided)
#' l1s_2 <- flm_est(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s",
#'                  p = c(1, 5, 2, 7), q = 1:4)
#' l1s_2$Beta_hat_scores
#' l1s_2$p_thre
#' l1s_2$p_hat
#'
#' ## Scalar response
#'
#' # Generate data
#' set.seed(12345)
#' n <- 50
#' beta <- r_ou(n = 1, t = seq(0, 1, l = 201), sigma = 0.5, x0 = 3)
#' X_fdata <- fdata_cen(r_ou(n = n, t = seq(0, 1, l = 201), sigma = 2))
#' epsilon <- rnorm(n, sd = 0.25)
#' Y <- drop(inprod_fdata(X_fdata1 = X_fdata, X_fdata2 = beta)) + epsilon
#'
#' # FPCR
#' fpcr_4 <- flm_est(X = X_fdata, Y = Y, est_method = "fpcr")
#' fpcr_4$p_hat
#'
#' # Ridge FPCR
#' l2_4 <- flm_est(X = X_fdata, Y = Y, est_method = "fpcr_l2")
#' l2_4$p_hat
#'
#' # Lasso FPCR
#' l1_4 <- flm_est(X = X_fdata, Y = Y, est_method = "fpcr_l1")
#' l1_4$p_hat
#'
#' # Lasso-selection FPCR
#' l1s_4 <- flm_est(X = X_fdata, Y = Y, est_method = "fpcr_l1s")
#' l1s_4$p_hat
#'
#' ## Scalar predictor
#'
#' # Generate data
#' set.seed(12345)
#' n <- 50
#' X <- rnorm(n)
#' epsilon <- r_ou(n = n, t = seq(0, 1, l = 201), sigma = 0.5)
#' beta <- r_ou(n = 1, t = seq(0, 1, l = 201), sigma = 0.5, x0 = 3)
#' beta$data <- matrix(beta$data, nrow = n, ncol = ncol(beta$data),
#'                     byrow = TRUE)
#' Y_fdata <- beta * X + epsilon
#'
#' # FPCR
#' fpcr_4 <- flm_est(X = X, Y = Y_fdata, est_method = "fpcr")
#' plot(beta, col = 2)
#' lines(beta$argvals, drop(fpcr_4$Beta_hat))
#' }
#' @author Eduardo García-Portugués and Javier Álvarez-Liébana.
#' @references
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' González-Manteiga, W. (2019). A goodness-of-fit test for the functional
#' linear model with functional response. \emph{arXiv:1909.07686}.
#' \url{https://arxiv.org/abs/1909.07686}
#'
#' Ramsay, J. and Silverman, B. W. (2005). \emph{Functional Data Analysis}.
#' Springer-Verlag, New York.
#' @export
flm_est <- function(X, Y, est_method = "fpcr_l1s",
                    p = NULL, q = NULL, thre_p = 0.99, thre_q = 0.99,
                    lambda = NULL, X_fpc = NULL, Y_fpc = NULL,
                    compute_residuals = TRUE, centered = FALSE,
                    int_rule = "trapezoid", cv_verbose = FALSE, ...) {

  # Check if the estimation method is implemented
  if (!(est_method %in% c("fpcr", "fpcr_l1", "fpcr_l2", "fpcr_l1s"))) {

    stop(paste("Estimation method", est_method, "not implemented"))

  }

  # Check if X_fpc and Y_fpc make sense
  if (!(class(X_fpc) %in% c("fpc", "NULL"))) {

    stop(paste("If X_fpc is non-NULL, then it must be an \"fpc\" class object",
               "as returned by the fpc function"))

  }
  if (!(class(Y_fpc) %in% c("fpc", "NULL"))) {

    stop(paste("If Y_fpc is non-NULL, then it must be an \"fpc\" class object",
               "as returned by the fpc function"))

  }

  # Check that p and q are positive (real p or q are truncated to integers)
  if (any(p <= 0)) {

    stop("p must be a vector of positive integers greater than zero")

  }
  if (any(q <= 0)) {

    stop("q must be a vector of positive integers greater than zero")

  }

  # Check thresholds (threshold = 1 is OK)
  if (!is.numeric(thre_p) | (thre_p <= 0 | thre_p > 1)) {

    stop(paste("Parameter thre_p must be a real numeric parameter greater",
               "than 0 and less than 1. Value of 1 is allowed."))

  }
  if (!is.numeric(thre_q) | (thre_q <= 0 | thre_q > 1)) {

    stop(paste("Parameter thre_q must be a real numeric parameter greater",
               "than 0 and less than 1. Value of 1 is allowed."))

  }

  # Scalar predictor or response? Transform into an fdata object
  # and keep track of the different nature of the input
  scalar_X <- is.numeric(X)
  scalar_Y <- is.numeric(Y)
  if (scalar_X) {

    # p is known in this case
    p <- 1

    # As fdata
    X <- fda.usc::fdata(mdata = cbind(X), argvals = 0)

    # The only sensible method is ordinary least squares
    if (est_method != "fpcr") {

      warning(paste("Setting est_method to \"fpcr\", as is the only",
                    "method available for scalar predictor"))
      est_method <- "fpcr"

    }

  }
  if (scalar_Y) {

    # q is known in this case
    q <- 1

    # As fdata
    Y <- fda.usc::fdata(mdata = cbind(Y), argvals = 0)

  }

  # Check sample sizes
  if (nrow(X[["data"]]) != nrow(Y[["data"]])) {

    stop("The sample sizes of X and Y do not match")

  }

  # Missing data on X? Drop the pairs (Xi, Yi) with missing values on Xi
  ind_X_NA <- !stats::complete.cases(X[["data"]])
  if (any(ind_X_NA)) {

    n_NA <- sum(ind_X_NA)
    if (cv_verbose) {

      message("Missing values in X, dropping ", n_NA,
              " curves with NAs (new sample size n = ",
              nrow(X[["data"]]) - n_NA, ")")

    }
    X <- X[!ind_X_NA]
    Y <- Y[!ind_X_NA]

  }

  # Missing at random situation on Y?
  ind_Y_NA <- !stats::complete.cases(Y[["data"]])
  if (any(ind_Y_NA)) {

    n_NA <- sum(ind_Y_NA)
    if (cv_verbose) {

      message("\nMissing values in Y, dropping ", n_NA,
              " responses with NAs (new sample size n = ",
              length(Y) - n_NA, ")")

    }
    X <- X[!ind_Y_NA]
    Y <- Y[!ind_Y_NA]

  }

  # Sample size after removing possible NA's
  n <- nrow(X[["data"]])

  # Center the data if it was not done before
  if (!centered) {

    X <- fdata_cen(X_fdata = X)
    Y <- fdata_cen(X_fdata = Y)

  }

  # Check if data is equispaced for FPC
  eps <- sqrt(.Machine[["double.eps"]])
  equispaced_X <- all(abs(diff(X[["argvals"]], differences = 2)) < eps)
  equispaced_Y <- all(abs(diff(Y[["argvals"]], differences = 2)) < eps)

  # FPC of X and Y (skip them if any is scalar: the only eigenfunction is just
  # the identity function valued in a single grid point)
  if (scalar_X) {

    X_fpc <- list("scores" = cbind(X[["data"]]),
                  "rotation" = fda.usc::fdata(mdata = cbind(1),
                                              argvals = X[["argvals"]]))

  } else {

    if (is.null(X_fpc)) {

      X_fpc <- fpc(X_fdata = X, n_fpc = n, centered = TRUE,
                   int_rule = int_rule, equispaced = equispaced_X)

    }

  }
  if (scalar_Y) {

    Y_fpc <- list("scores" = cbind(Y[["data"]]),
                  "rotation" = fda.usc::fdata(mdata = cbind(1),
                                              argvals = Y[["argvals"]]))

  } else {

    if (is.null(Y_fpc)) {

      Y_fpc <- fpc(X_fdata = Y, n_fpc = n, centered = TRUE,
                   int_rule = int_rule, equispaced = equispaced_Y)

    }

  }

  # Data-driven truncations? These will be the estimated p's and q's for
  # "fpcr" and "fpcr_l2", but not for "fpcr_l1" and "fpcr_l1s"
  if (is.null(p)) {

    s <- cumsum(X_fpc[["d"]]^2)
    p_thre <- 1:which(s / s[length(s)] > thre_p)[1]

  } else {

    # Set p as vector
    if (length(p) == 1) {

      p <- 1:p

    }
    p_thre <- p

  }
  if (is.null(q)) {

    s <- cumsum(Y_fpc[["d"]]^2)
    q_thre <- 1:which(s / s[length(s)] > thre_q)[1]

  } else {

    # Set q as vector
    if (length(q) == 1) {

      q <- 1:q

    }
    q_thre <- q

  }

  # Carry out the estimation of the multivariate linear model

  # FPCR
  if (est_method == "fpcr") {

    # Compute estimator and hat matrix
    H_hat <- solve(a = crossprod(X_fpc[["scores"]][, p_thre, drop = FALSE]),
                   b = t(X_fpc[["scores"]][, p_thre]))
    Beta_hat_scores <- H_hat %*% Y_fpc[["scores"]][, q_thre, drop = FALSE]
    H_hat <- X_fpc[["scores"]][, p_thre, drop = FALSE] %*% H_hat

    # Does not change the estimated dimension p
    p_hat <- p_thre

    # Cross-validation object
    cv <- NULL

  # FPC + ridge regression
  } else if (est_method == "fpcr_l2") {

    # Find the optimal lambda for ridge regression (if lambda = NULL)
    if (is.null(lambda)) {

      cv <- cv_glmnet(x = X_fpc[["scores"]][, p_thre, drop = FALSE],
                      y = Y_fpc[["scores"]][, q_thre, drop = FALSE],
                      alpha = "ridge", cv_verbose = cv_verbose,
                      intercept = FALSE, ...)
      lambda <- cv[["lambda"]]

    } else{

      cv <- NULL

    }

    # Compute estimator and hat matrix. For that, we have to change the scale
    # of lambda from the standardized scale to the unstandardized scale (see
    # examples in ?cv_glmnet)
    var_x <- apply(X_fpc[["scores"]][, p_thre, drop = FALSE], 2, var) *
      (n - 1) / n
    sd_y <- ifelse(length(q_thre) == 1,
                   sd(Y_fpc[["scores"]][, q_thre, drop = FALSE]) *
                     sqrt((n - 1) / n),
                   1)
    H_hat <- solve(a = crossprod(X_fpc[["scores"]][, p_thre, drop = FALSE]) +
                     diag(x = n * var_x / sd_y * lambda, nrow = length(p_thre)),
                   b = t(X_fpc[["scores"]][, p_thre, drop = FALSE]))
    Beta_hat_scores <- H_hat %*% Y_fpc[["scores"]][, q_thre, drop = FALSE]
    H_hat <- X_fpc[["scores"]][, p_thre, drop = FALSE] %*% H_hat

    # Does not change the estimated dimension p
    p_hat <- p_thre

  # FPC + lasso regression
  } else if (est_method == "fpcr_l1") {

    # Fit the lasso and find the optimal lambda (if lambda = NULL)
    cv <- cv_glmnet(x = X_fpc[["scores"]][, p_thre, drop = FALSE],
                    y = Y_fpc[["scores"]][, q_thre, drop = FALSE],
                    alpha = "lasso", lambda = lambda,
                    intercept = FALSE, cv_verbose = cv_verbose, ...)
    lambda <- cv[["lambda"]]

    # What are the non-null coefficients? Updates p_thre with the feature
    # selection performed by lasso regression (zeroing by rows)
    p_hat <- p_thre[which(apply(abs(cv[["beta_hat"]]) > 0, 1, any))]

    # p_hat can have length zero, causing further problems when computing
    # Beta_hat_scores and Beta_hat. To avoid this, set p_hat = 1 and
    # Beta_hat_scores = 0 when computing the estimator
    if (length(p_hat) == 0) {

      p_hat <- 1
      Beta_hat_scores <- rbind(rep(0, length(q_thre)))

    } else {

      Beta_hat_scores <- as.matrix(cv[["beta_hat"]][p_hat, q_thre,
                                                    drop = FALSE])

    }

    # Explicit hat matrix cannot be computed for the lasso regression
    H_hat <- NULL

  # FPCR + lasso regression + reduced linear model
  } else if (est_method == "fpcr_l1s") {

    # Fit the lasso and find the optimal lambda (if lambda = NULL)
    cv <- cv_glmnet(x = X_fpc[["scores"]][, p_thre, drop = FALSE],
                    y = Y_fpc[["scores"]][, q_thre, drop = FALSE],
                    alpha = "lasso", lambda = lambda,
                    intercept = FALSE, cv_verbose = cv_verbose, ...)
    lambda <- cv[["lambda"]]

    # What are the non-null coefficients? Updates p_thre with the feature
    # selection performed by lasso regression (zeroing by rows)
    p_hat <- p_thre[which(apply(abs(cv[["beta_hat"]]) > 0, 1, any))]

    # p_hat can have length zero, causing further problems when computing
    # Beta_hat_scores and Beta_hat. To avoid this, set p_hat = 1 and
    # Beta_hat_scores = 0 when computing the estimator
    if (length(p_hat) == 0) {

      p_hat <- 1
      Beta_hat_scores <- rbind(rep(0, length(q_thre)))
      H_hat <- matrix(0, nrow = n, ncol = n)

    } else {

      H_hat <- solve(a = crossprod(X_fpc[["scores"]][, p_hat, drop = FALSE]),
                     b = t(X_fpc[["scores"]][, p_hat, drop = FALSE]))
      Beta_hat_scores <- H_hat %*% Y_fpc[["scores"]][, q_thre, drop = FALSE]
      H_hat <- X_fpc[["scores"]][, p_hat, drop = FALSE] %*% H_hat

    }

  }

  # Computing estimated surface
  Beta_hat <- fpc_to_beta(beta_coefs = Beta_hat_scores, X_fpc = X_fpc,
                          Y_fpc = Y_fpc, ind_coefs_X = p_hat,
                          ind_coefs_Y = q_thre)

  # Compute fitted values and residuals
  if (compute_residuals) {

    # Fitted values of the FPC coefficients of the response (used for
    # defining the functional fitted values)
    Y_hat_scores <- X_fpc[["scores"]][, p_hat, drop = FALSE] %*%
      Beta_hat_scores

    # Functional fitted values (same object as Y)
    Y_hat <- fpc_to_fdata(coefs = Y_hat_scores, X_fpc = Y_fpc,
                          ind_coefs = q_thre)
    if (scalar_Y) {

      Y_hat <- drop(Y_hat[["data"]])

    }

    # Residuals FPC coefficients
    residuals_scores <- Y_fpc[["scores"]][, q_thre, drop = FALSE] - Y_hat_scores

    # Functional residuals (the difference works for fdata and vector objects)
    residuals <- Y - Y_hat
    if (scalar_Y) {

      residuals <- drop(residuals)

    }

  } else {

    # Empty outputs
    Y_hat <- NULL
    Y_hat_scores <- NULL
    residuals <- NULL
    residuals_scores <- NULL

  }

  # Final object
  return(list("Beta_hat" = Beta_hat, "Beta_hat_scores" = Beta_hat_scores,
              "H_hat" = H_hat, "p_thre" = p_thre, "q_thre" = q_thre,
              "p_hat" = p_hat, "est_method" = est_method,
              "Y_hat" = Y_hat, "Y_hat_scores" = Y_hat_scores,
              "residuals" = residuals, "residuals_scores" = residuals_scores,
              "X_fpc" = X_fpc, "Y_fpc" = Y_fpc, "lambda" = lambda, "cv" = cv))

}
