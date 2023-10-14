

#' @title Goodness-of-fit test for functional linear models
#'
#' @description Goodness-of-fit test of a functional linear model with
#' functional response \eqn{Y \in L^2([c, d])} and functional predictor
#' \eqn{X \in L^2([a, b])}, where \eqn{L^2([a, b])} is the Hilbert space of
#' square-integrable functions in \eqn{[a, b]}.
#'
#' The goodness-of-fit test checks the \emph{linearity} of the regression model
#' \eqn{m:L^2([a, b])\rightarrow L^2([c, d])} that relates \eqn{Y} and \eqn{X}
#' by
#' \deqn{Y(t) = m(X) + \varepsilon(t),}{Y(t) = m(X)(t) + \epsilon(t),}
#' where \eqn{\varepsilon}{\epsilon} is a random variable in
#' \eqn{L^2([c, d])} and \eqn{t \in [c, d]}. The check is formalized as the
#' test of the composite hypothesis
#' \deqn{H_0: m \in \{m_\beta : \beta \in L^2([a, b]) \otimes L^2([c, d])\},}{
#' H_0: m \in {m_\beta : \beta \in L^2([a, b]) \otimes L^2([c, d])},}
#' where
#' \deqn{m_\beta(X(s))(t) = \int_a^b \beta(s, t) X(s)\,\mathrm{d}s}{
#' m_\beta(X(s))(t) = \int_a^b \beta(s, t) X(s) ds}
#' is the linear, Hilbert--Schmidt, integral operator parametrized by
#' the bivariate kernel \eqn{\beta}. Its estimation is done by the
#' truncated expansion of \eqn{\beta} in the tensor product of the
#' data-driven bases of \emph{Functional Principal Components} (FPC) of
#' \eqn{X} and \eqn{Y}. The FPC basis for \eqn{X} is truncated in \eqn{p}
#' components, while the FPC basis for \eqn{Y} is truncated in \eqn{q}
#' components.
#'
#' The particular cases in which either \eqn{X} or \eqn{Y} are
#' \emph{constant} functions give either a scalar predictor or response.
#' The simple linear model arises if both \eqn{X} and \eqn{Y} are scalar,
#' for which \eqn{\beta} is a constant.
#'
#' @param X,Y samples of functional/scalar predictors and functional/scalar
#' response. Either \code{\link[fda.usc]{fdata}} objects (for functional
#' variables) or vectors of length \code{n} (for scalar variables).
#' @param beta0 if provided (defaults to \code{NULL}), the \emph{simple} null
#' hypothesis \eqn{H_0: m = m_{\beta_0}} is tested. \code{beta0} must be a
#' matrix of size\cr \code{c(length(X$argvals), length(Y$argvals))}. If \code{X}
#' or \code{Y} are scalar, \code{beta0} can be also an
#' \code{\link[fda.usc]{fdata}} object, with the same \code{argvals} as
#' \code{X} or \code{Y}. Can also be a constant (understood as a shorthand for
#' a matrix with \emph{all} its entries equal to the constant).
#' @param B number of bootstrap replicates. Defaults to \code{500}.
#' @inheritParams flm_est
#' @param p,q either index vectors indicating the specific FPC to be
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
#' @param boot_scores flag to indicate if the bootstrap shall be applied to the
#' scores of the residuals, rather than to the functional residuals. This
#' improves the computational expediency notably. Defaults to \code{TRUE}.
#' @param verbose flag to show information about the testing progress. Defaults
#' to \code{TRUE}.
#' @param plot_dens flag to indicate if a kernel density estimation of the
#' bootstrap statistics shall be plotted. Defaults to \code{TRUE}.
#' @param plot_proc whether to display a graphical tool to identify the
#' degree of departure from the null hypothesis. If \code{TRUE} (default),
#' the residual marked empirical process, projected in several FPC directions
#' of \code{X} and \code{Y}, is shown, together with bootstrap analogues.
#' The FPC directions are ones selected at the estimation stage.
#' @param plot_max_procs maximum number of bootstrapped processes to plot in
#' the graphical tool. Set as the minimum of \code{plot_max_procs} and \code{B}.
#' Defaults to \code{100}.
#' @param plot_max_p,plot_max_q maximum number of FPC directions to be
#' considered in the graphical tool. They limit the resulting plot to be at
#' most of size \code{c(plot_max_p, plot_max_q)}. Default to \code{2}.
#' @param save_fit_flm,save_boot_stats flag to return \code{fit_flm} and
#' \code{boot_*}. If \code{FALSE}, these memory-expensive objects
#' are set to \code{NA}. Default to \code{TRUE}.
#' @param refit_lambda flag to reselect \eqn{lambda} in each bootstrap
#' replicate, incorporating its variability in the bootstrap calibration.
#' Much more time consumig. Defaults to \code{FALSE}.
#' @return An object of the \code{htest} class with the following elements:
#' \item{statistic}{test statistic.}
#' \item{p.value}{\eqn{p}-value of the test.}
#' \item{boot_statistics}{the bootstrapped test statistics, a vector
#' of length \code{B}.}
#' \item{method}{information on the type of test performed.}
#' \item{parameter}{a vector with the dimensions \eqn{p} and \eqn{q}
#' considered in the test statistic. These are the lengths of the outputs
#' \code{p} and \code{q}.}
#' \item{p}{the index of the FPC considered for \code{X}.}
#' \item{q}{the index of the FPC considered for \code{Y}.}
#' \item{fit_flm}{the output resulted from calling \code{\link{flm_est}}.}
#' \item{boot_lambda}{bootstrapped \eqn{lambda}.}
#' \item{boot_p}{a list with the bootstrapped indexes of the FPC considered
#' for \code{X}.}
#' \item{data.name}{name of the value of \code{data}.}
#' @details
#' The function implements the bootstrap-based goodness-of-fit test for
#' the functional linear model with functional/scalar response and
#' functional/scalar predictor, as described in Algorithm 1 in
#' García-Portugués et al. (2021). The specifics are detailed there.
#'
#' By default \code{cv_1se = TRUE} for \code{\link{cv_glmnet}} is
#' considered, unless it is changed via \code{...}. This is the recommended
#' choice for conducting the goodness-of-fit test based on regularized
#' estimators, as the oversmoothed estimate of the regression model under the
#' null hypothesis notably facilitates the calibration of the test (see
#' García-Portugués et al., 2021).
#'
#' The graphical tool obtained with \code{plot_proc = TRUE} is based on
#' an extension of the tool described in García-Portugués et al. (2014).
#'
#' Repeated observations on \code{X} are internally removed, as otherwise they
#' would cause \code{NaN}s in \code{Adot}. Missing values on \code{X} and
#' \code{Y} are also automatically removed.
#' @examples
#' ## Quick example for functional response and predictor
#'
#' # Generate data under H0
#' n <- 100
#' set.seed(987654321)
#' X_fdata <- r_ou(n = n, t = seq(0, 1, l = 101), sigma = 2)
#' epsilon <- r_ou(n = n, t = seq(0, 1, l = 101), sigma = 0.5)
#' Y_fdata <- epsilon
#'
#' # Test the FLMFR
#' flm_test(X = X_fdata, Y = Y_fdata)
#'
#' # Simple hypothesis
#' flm_test(X = X_fdata, Y = Y_fdata, beta0 = 0)
#'
#' # Generate data under H1
#' n <- 100
#' set.seed(987654321)
#' sample_frm_fr <- r_frm_fr(n = n, scenario = 3, s = seq(0, 1, l = 101),
#'                           t = seq(0, 1, l = 101), nonlinear = "quadratic")
#' X_fdata <- sample_frm_fr[["X_fdata"]]
#' Y_fdata <- sample_frm_fr[["Y_fdata"]]
#'
#' # Test the FLMFR
#' flm_test(X = X_fdata, Y = Y_fdata)
#' \donttest{
#' ## Functional response and predictor
#'
#' # Generate data under H0
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' t <- seq(0, 1, l = 201)
#' X_fdata <- r_ou(n = n, t = t, sigma = 2)
#' epsilon <- r_ou(n = n, t = t, sigma = 0.5)
#' Y_fdata <- epsilon
#'
#' # With boot_scores = TRUE
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s",
#'          boot_scores = FALSE, B = B)
#'
#' # Simple hypothesis
#' flm_test(X = X_fdata, Y = Y_fdata, beta0 = 2, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, beta0 = 0, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, beta0 = 0, est_method = "fpcr_l1s", B = B)
#'
#' # Generate data under H1
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' sample_frm_fr <- r_frm_fr(n = n, scenario = 3, s = t, t = t,
#'                           nonlinear = "quadratic")
#' X_fdata <- sample_frm_fr$X_fdata
#' Y_fdata <- sample_frm_fr$Y_fdata
#'
#' # With boot_scores = TRUE
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2", B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l2",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y_fdata, est_method = "fpcr_l1s",
#'          boot_scores = FALSE, B = B)
#'
#' ## Scalar response and functional predictor
#'
#' # Generate data under H0
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' t <- seq(0, 1, l = 201)
#' X_fdata <- r_ou(n = n, t = t, sigma = 2)
#' beta <- r_ou(n = 1, t = t, sigma = 0.5, x0 = 2)
#' epsilon <- rnorm(n = n)
#' Y <- drop(inprod_fdata(X_fdata1 = X_fdata, X_fdata2 = beta) + epsilon)
#'
#' # With boot_scores = TRUE
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l2", B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1s", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l2",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1s",
#'          boot_scores = FALSE, B = B)
#'
#' # Simple hypothesis
#' flm_test(X = X_fdata, Y = Y, beta0 = beta, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y, beta0 = 0, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y, beta0 = 0, est_method = "fpcr_l1s", B = B)
#'
#' # Generate data under H1
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' X_fdata <- r_ou(n = n, t = t, sigma = 2)
#' beta <- r_ou(n = 1, t = t, sigma = 0.5)
#' epsilon <- rnorm(n = n)
#' Y <- drop(exp(inprod_fdata(X_fdata1 = X_fdata^2, X_fdata2 = beta)) + epsilon)
#'
#' # With boot_scores = TRUE
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr", B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l2", B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1s", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l2",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1",
#'          boot_scores = FALSE, B = B)
#' flm_test(X = X_fdata, Y = Y, est_method = "fpcr_l1s",
#'          boot_scores = FALSE, B = B)
#'
#' ## Functional response and scalar predictor
#'
#' # Generate data under H0
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' X <- rnorm(n)
#' t <- seq(0, 1, l = 201)
#' beta <- r_ou(n = 1, t = t, sigma = 0.5, x0 = 3)
#' beta$data <- matrix(beta$data, nrow = n, ncol = ncol(beta$data),
#'                     byrow = TRUE)
#' epsilon <- r_ou(n = n, t = t, sigma = 0.5)
#' Y_fdata <- X * beta + epsilon
#'
#' # With boot_scores = TRUE
#' flm_test(X = X, Y = Y_fdata, est_method = "fpcr", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X, Y = Y_fdata, est_method = "fpcr", boot_scores = FALSE, B = B)
#'
#' # Simple hypothesis
#' flm_test(X = X, Y = Y_fdata, beta0 = beta[1], est_method = "fpcr", B = B)
#' flm_test(X = X, Y = Y_fdata, beta0 = 0, est_method = "fpcr", B = B)
#'
#' # Generate data under H1
#' n <- 50
#' B <- 100
#' set.seed(987654321)
#' X <- rexp(n)
#' beta <- r_ou(n = 1, t = t, sigma = 0.5, x0 = 3)
#' beta$data <- matrix(beta$data, nrow = n, ncol = ncol(beta$data),
#'                     byrow = TRUE)
#' epsilon <- r_ou(n = n, t = t, sigma = 0.5)
#' Y_fdata <- log(X * beta) + epsilon
#'
#' # With boot_scores = TRUE
#' flm_test(X = X, Y = Y_fdata, est_method = "fpcr", B = B)
#'
#' # With boot_scores = FALSE
#' flm_test(X = X, Y = Y_fdata, est_method = "fpcr", boot_scores = FALSE, B = B)
#' }
#' @author Eduardo García-Portugués.
#' @references
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' Gonzalez-Manteiga, W. (2021). A goodness-of-fit test for the functional
#' linear model with functional response. \emph{Scandinavian Journal of
#' Statistics}, 48(2):502--528. \doi{10.1111/sjos.12486}
#'
#' García-Portugués, E., González-Manteiga, W. and Febrero-Bande, M. (2014). A
#' goodness-of-fit test for the functional linear model with scalar response.
#' \emph{Journal of Computational and Graphical Statistics}, 23(3):761--778.
#' \doi{10.1080/10618600.2013.812519}
#' @export
flm_test <- function(X, Y, beta0 = NULL, B = 500, est_method = "fpcr",
                     p = NULL, q = NULL, thre_p = 0.99, thre_q = 0.99,
                     lambda = NULL, boot_scores = TRUE, verbose = TRUE,
                     plot_dens = TRUE, plot_proc = TRUE, plot_max_procs = 100,
                     plot_max_p = 2, plot_max_q = 2, save_fit_flm = TRUE,
                     save_boot_stats = TRUE, int_rule = "trapezoid",
                     refit_lambda = FALSE, ...) {

  ## Preprocessing

  # Detect non-implemented bootstrap resampling
  if (boot_scores && (est_method == "fpcr_l1")) {

    stop(paste("Bootstrap resampling of the residuals FPC coefficients is",
               "not implemented for est_method = \"fpcr_l1\""))

  }

  # Data name
  data_name <- paste(deparse(substitute(Y)), deparse(substitute(X)),
                     sep = " ~ ")

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

  # Repeated observations on X? Drop the pairs (Xi, Yi) with duplicated
  # values on Xi or thei will cause NaNs in Adot()
  ind_X_dup <- duplicated(X[["data"]])
  if (any(ind_X_dup)) {

    n_dup <- sum(ind_X_dup)
    if (verbose) {

      message("Repeated observations in X, dropping ", n_dup,
              " repeated curves (new sample size n = ",
              nrow(X[["data"]]) - n_dup, ")")

    }
    X <- X[!ind_X_dup]
    Y <- Y[!ind_X_dup]

  }

  # Missing data on X? Drop the pairs (Xi, Yi) with missing values on Xi
  ind_X_NA <- !stats::complete.cases(X[["data"]])
  if (any(ind_X_NA)) {

    n_NA <- sum(ind_X_NA)
    if (verbose) {

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
    if (verbose) {

      message("\nMissing values in Y, dropping ", n_NA,
              " responses with NAs (new sample size n = ",
              length(Y) - n_NA, ")")

    }
    X <- X[!ind_Y_NA]
    Y <- Y[!ind_Y_NA]

  }

  # Sample size after removing possible NA's
  n <- nrow(X[["data"]])

  # Center the data
  X <- fdata_cen(X_fdata = X)
  Y <- fdata_cen(X_fdata = Y)

  ## Estimation of beta

  # Display progress
  if (verbose) {

    message("Estimating beta... ", appendLF = FALSE)

  }

  # Fit the model
  fit_flm <- flm_est(X = X, Y = Y, est_method = est_method, p = p, q = q,
                     thre_p = thre_p, thre_q = thre_q, lambda = lambda,
                     X_fpc = NULL, Y_fpc = NULL, centered = TRUE,
                     int_rule = int_rule, cv_verbose = FALSE, ...)

  # Estimated FPC indexes
  p_thre <- fit_flm[["p_thre"]]
  q_thre <- fit_flm[["q_thre"]]
  p_hat <- fit_flm[["p_hat"]]

  # Dimensions of the predictor and response
  p <- length(p_hat)
  q <- length(q_thre)

  # Composite hypothesis
  if (is.null(beta0)) {

    # Method
    meth <- paste0("Goodness-of-fit test for the functional linear model ",
                  "(", ifelse(scalar_X, "scalar", "functional"),
                  " predictor and ", ifelse(scalar_Y, "scalar", "functional"),
                  " response)")

    # Residuals FPC coefficients
    E_hat_scores <- fit_flm[["residuals_scores"]]

    # Specific objects required for each type of bootstrap
    if (boot_scores) {

      # Projection matrix
      P_hat <- diag(1, nrow = n, ncol = n) - fit_flm[["H_hat"]]

      # Prediction FPC coefficients
      Y_hat_scores <- fit_flm[["Y_hat_scores"]]

    } else {

      # Functional predictions
      Y_hat <- fit_flm[["Y_hat"]]

      # Functional residuals
      E_hat <- fit_flm[["residuals"]]

    }

  # Simple hypothesis
  } else {

    # Method
    meth <- paste0("Goodness-of-fit test for the functional linear model ",
                   "(", ifelse(scalar_X, "scalar", "functional"),
                   " predictor and ", ifelse(scalar_Y, "scalar", "functional"),
                   " response; simple hypothesis)")

    # Check if matrix or fdata
    if (!is.numeric(beta0)) {

      if (fda.usc::is.fdata(beta0)) {

        if (scalar_X || scalar_Y) {

          l_beta0 <- length(beta0)
          if (l_beta0 != 1) {

            stop("beta0 must be a single fdata, but has length ", l_beta0)

          } else {

            beta0 <- beta0[["data"]]
            if (scalar_Y) {

              beta0 <- t(beta0)

            }

          }

        } else {

          stop("beta0 must be a matrix")

        }

      } else {

        stop("beta0 must be a matrix or an fdata object")

      }

    }

    # If a single number is passed (beta0 = 0, e.g.), then expand it to a vector
    if (length(beta0) == 1) {

      beta0 <- rep(beta0, ifelse(scalar_X, 1, length(X[["argvals"]])) *
                     ifelse(scalar_Y, 1, length(Y[["argvals"]])))

    }

    # Set dimensions if necessary
    if (is.null(dim(beta0))) {

      if (scalar_Y) {

        beta0 <- cbind(beta0)

      } else if (scalar_X) {

        beta0 <- rbind(beta0)

      } else if (length(beta0) == length(X$argvals) * length(Y$argvals)) {

        beta0 <- matrix(beta0, nrow = length(X$argvals),
                        ncol = length(Y$argvals))

      }

    }

    # Check adequate dimensions
    if (length(X$argvals) != nrow(beta0) ||
        length(Y$argvals) != ncol(beta0)) {

      stop("beta0 must be a matrix of size ", length(X$argvals), " x ",
           length(Y$argvals))

    }

    # Express beta0 in the FPC bases of X and Y
    beta0_coefs <- beta_fpc_coefs(beta = beta0, X_fpc = fit_flm[["X_fpc"]],
                                  Y_fpc = fit_flm[["Y_fpc"]],
                                  ind_X_fpc = p_hat, ind_Y_fpc = q_thre,
                                  int_rule = int_rule)

    # Prediction FPC coefficients
    Y_hat_scores <- fit_flm[["X_fpc"]][["scores"]][, p_hat, drop = FALSE] %*%
      beta0_coefs

    # Residuals FPC coefficients
    E_hat_scores <- fit_flm[["Y_fpc"]][["scores"]][, q_thre, drop = FALSE] -
      Y_hat_scores

    # If bootstrap is done on the functional errors
    if (!boot_scores) {

      # Functional fitted values (same object as Y)
      Y_hat <- fpc_to_fdata(coefs = Y_hat_scores, X_fpc = fit_flm[["Y_fpc"]],
                            ind_coefs = q_thre)
      if (scalar_Y) {

        Y_hat <- drop(Y_hat)

      }

      # Functional residuals
      E_hat <- Y - Y_hat

    }

  }

  ## Statistic

  # Display progress
  if (verbose) {

    message("Done.\nComputing statistic... ", appendLF = FALSE)

  }

  # Adot matrix for using in the original and bootstraped statistics. Observe
  # that the inner products returned by inprod_fdata(X) form the vector of the
  # upper triangular of the inner product matrix that stacks the columns
  Adot_vec <- Adot(X = fit_flm[["X_fpc"]][["scores"]])

  # Original statistic. Constant computed only if refit_lambda = TRUE as in
  # that case there are different normalizing constants depending on the
  # selected p
  orig_stat <- flm_stat(E = E_hat_scores, p = p, Adot_vec = Adot_vec,
                        constant = refit_lambda)

  ## Bootstrap

  if (verbose) {

    message("Done.\nBootstrap calibration...\n", appendLF = FALSE)
    pb <- txtProgressBar(style = 3)

  }

  # Method employed
  est_method_star <- switch(est_method, "fpcr_l1s" = "fpcr", est_method)
  p_star <- switch(est_method, "fpcr_l1s" = p_hat, p_thre)

  # Objects employed in the bootstrap resampling
  phi <- (1 + sqrt(5)) / 2
  prob <- (phi + 2) / 5
  boot_stats <- rep(NA, B)
  boot_lambda <- rep(ifelse(is.null(fit_flm[["lambda"]]),
                            NA, fit_flm[["lambda"]]), B)
  boot_p_hat <- rep(list(p_hat), B)

  # Array with the estimated bootstrap errors, saved for the processes plot
  # done plot_max_procs bootstrap processes
  plot_max_procs <- min(plot_max_procs, B)
  E_star_hat_scores <-
    array(dim = c(n, q, ifelse(plot_proc, plot_max_procs, 1)))

  # Composite hypothesis
  if (is.null(beta0)) {

    # Bootstrap resampling
    for (i in 1:B) {

      # Save the bootstrap scores? Only plot_max_procs of them in a cyclic way
      j <- ifelse(plot_proc, (i - 1) %% plot_max_procs + 1, 1)

      # Golden section binary variable
      V <- sample(x = c(1 - phi, phi), prob = c(prob, 1 - prob),
                  size = n, replace = TRUE)

      # Bootstrap on the FPC coefficients of the functional residuals
      if (boot_scores) {

        # Perturb FPC coefficients of the residuals. Each row of E_hat_scores
        # is multiplied by the *same* Vi (there is an implicit recycling)
        E_star_scores <- E_hat_scores * V

        # Obtain new bootstrap observations of the FPC coefficients
        Y_star_scores <- Y_hat_scores + E_star_scores

        # Recenter the bootstrap observations before refitting the model,
        # imitating what it is done in the original fitting by flm_est
        Y_star_scores <- t(t(Y_star_scores) - colMeans(Y_star_scores))

        # Keep lambda fixed in the bootstrap (much faster)
        if (!refit_lambda) {

          # Avoid refitting the model by computing the residuals of refitted
          # model as projections of Y_star_scores
          E_star_hat_scores[, , j] <- P_hat %*% Y_star_scores

        # Refit lambda to capture the variability on its selection
        } else {

          # Reconstruct the functional response from bootstrapped scores,
          # useful for recalling flm_est using an fdata response
          Y_star <- fpc_to_fdata(coefs = Y_star_scores,
                                 X_fpc = fit_flm[["Y_fpc"]])

          # Refit model searching for the optimal lambda and using p_thre
          # (instead of p_hat)
          fit_flm_star <- flm_est(X = X, Y = Y_star, est_method = est_method,
                                  p = p_thre, q = q_thre, lambda = lambda,
                                  X_fpc = fit_flm[["X_fpc"]], Y_fpc = NULL,
                                  centered = TRUE, int_rule = int_rule,
                                  cv_verbose = FALSE, ...)

          # Residuals FPC coefficients
          E_star_hat_scores[, , j] <- fit_flm_star[["residuals_scores"]]

          # Bootstrapped lambda and p_hat
          boot_lambda[i] <- ifelse(is.null(fit_flm_star[["lambda"]]),
                                   NA, fit_flm_star[["lambda"]])
          boot_p_hat[[i]] <- fit_flm_star[["p_hat"]]

        }

      # Bootstrap on the functional residuals
      } else {

        # Perturb functional residuals: each functional observation is
        # multiplied by one Vi
        E_star <- E_hat * V

        # Obtain new bootstrap functional observations
        Y_star <- Y_hat + E_star

        # Center Y_star to avoid the internal recentering in flm_est (because
        # it will recenter also X, which is already centered)
        Y_star <- fdata_cen(X_fdata = Y_star)

        # Keep lambda fixed in the bootstrap (much faster)
        if (!refit_lambda) {

          # Refit the model using the same tuning parameters (lambda, p_thre and
          # q_thre) and reusing the FPC of X. Observe that p_hat is *not* used
          # because passing p = p_hat if the estimation method is "fpcr_l1" or
          # "fpcr_l1s" will result in a different fit
          E_star_hat_scores[, , j] <- flm_est(X = X, Y = Y_star,
                                              est_method = est_method_star,
                                              p = p_star, q = q_thre,
                                              lambda = fit_flm[["lambda"]],
                                              X_fpc = fit_flm[["X_fpc"]],
                                              Y_fpc = NULL, centered = TRUE,
                                              int_rule = int_rule,
                                              cv_verbose = FALSE,
                                              ...)[["residuals_scores"]]

        # Refit lambda to capture the variability on its selection
        } else {

          # Refit model searching for the optimal lambda and using p_thre
          # (instead of p_hat)
          fit_flm_star <- flm_est(X = X, Y = Y_star, est_method = est_method,
                                  p = p_thre, q = q_thre, lambda = lambda,
                                  X_fpc = fit_flm[["X_fpc"]], Y_fpc = NULL,
                                  centered = TRUE, int_rule = int_rule,
                                  cv_verbose = FALSE, ...)

          # Residuals FPC coefficients
          E_star_hat_scores[, , j] <- fit_flm_star[["residuals_scores"]]

          # Bootstrapped lambda and p_hat
          boot_lambda[i] <- ifelse(is.null(fit_flm_star[["lambda"]]),
                                   NA, fit_flm_star[["lambda"]])
          boot_p_hat[[i]] <- fit_flm_star[["p_hat"]]

        }

      }

      # Bootstrap statistic
      boot_stats[i] <- flm_stat(E = as.matrix(E_star_hat_scores[, , j]),
                                p = p, Adot_vec = Adot_vec,
                                constant = refit_lambda)

      # Display progress
      if (verbose) {

        setTxtProgressBar(pb, i / B)

      }

    }

  # Simple hypothesis
  } else {

    # Bootstrap resampling
    for (i in 1:B) {

      # Save the bootstrap scores? Only plot_max_procs of them in a cyclic way
      j <- ifelse(plot_proc, (i - 1) %% plot_max_procs + 1, 1)

      # Golden section binary variable
      V <- sample(x = c(1 - phi, phi), prob = c(prob, 1 - prob),
                  size = n, replace = TRUE)

      # Bootstrap on the FPC coefficients of the functional residuals
      if (boot_scores) {

        # Perturb FPC coefficients of the residuals. Each row of E_hat_scores
        # is multiplied by the *same* Vi (there is an implicit recycling)
        E_star_hat_scores[, , j] <- E_hat_scores * V

        # Recenter the bootstrap scores
        E_star_hat_scores[, , j] <-
          t(t(E_star_hat_scores[, , j]) -
              colMeans(as.matrix(E_star_hat_scores[, , j])))

      # Bootstrap on the functional residuals
      } else {

        # Perturb functional residuals: each functional observation is
        # multiplied by one Vi
        E_star <- E_hat * V

        # Obtain new bootstrap functional observations
        Y_star <- Y_hat + E_star

        # Center Y_star
        Y_star <- fdata_cen(X_fdata = Y_star)

        # Responses FPC coefficients
        Y_star_scores <- fpc_coefs(X_fdata = Y_star, X_fpc = fit_flm[["Y_fpc"]],
                                   ind_X_fpc = q_thre)

        # Residuals FPC coefficients
        E_star_hat_scores[, , j] <- Y_star_scores - Y_hat_scores

      }

      # Bootstrap statistic
      boot_stats[i] <- flm_stat(E = as.matrix(E_star_hat_scores[, , j]),
                                p = p, Adot_vec = Adot_vec, constant = FALSE)

      # Display progress
      if (verbose) {

        setTxtProgressBar(pb, i / B)

      }

    }

  }

  ## p-value and final result

  # Add constants to statistics
  const <- ifelse(refit_lambda, 1, 2 * pi^(0.5 * (p + q) - 1) /
                    (q * gamma(0.5 * p) * gamma(0.5 * q) * n^2))
  orig_stat <- const * orig_stat
  boot_stats <- const * boot_stats

  # Approximation of the p-value by MC
  p_value <- mean(orig_stat < boot_stats)

  # Plot the density of the bootstrapped statistics and the value of
  # the statistic
  if (plot_dens) {

    plot(ks::kde(x = boot_stats, positive = TRUE),
         xlim = range(c(boot_stats, orig_stat)),
         main = paste("p-value:", p_value))
    rug(boot_stats)
    abline(v = orig_stat, col = 2)

  }

  # Plot the original and bootstrapped processes projected on FPC directions
  if (plot_proc) {

    # Function that computes the residual marked projected empirical process
    # based on two-side projections Rn(u, gamma_X, gamma_Y). It does so for
    # gamma_X and gamma_Y being a FPC of X or Y, respectively, such that the
    # projections of X become the scores of X and the projected residuals are
    # the scores of the residuals on the FPC of Y
    Rn <- function(X_scores, E_hat_scores, ind_X_fpc, ind_Y_fpc) {

      # Projections on gamma_X and gamma_Y as scores
      proj_X_gamma <- X_scores[, ind_X_fpc]
      proj_E_gamma <- E_hat_scores[, ind_Y_fpc]

      # Rn(u, gamma_X, gamma_Y)
      ord <- order(proj_X_gamma)
      y <- cumsum(proj_E_gamma[ord]) / sqrt(nrow(E_hat_scores))
      stepfun(x = proj_X_gamma[ord], y = c(0, y))

    }

    # Plotting function with original and bootstrapped processes
    plot_Rn <- function(X_scores, E_hat_scores, E_star_hat_scores,
                        ind_X_fpc, ind_Y_fpc, main = "") {

      # Rn processes, both for original and bootstrapped data
      Rn_processes <- c(
        Rn(X_scores = X_scores, E_hat_scores = E_hat_scores,
           ind_X_fpc = ind_X_fpc, ind_Y_fpc = ind_Y_fpc),
        lapply(seq_len(dim(E_star_hat_scores)[3]), function(i) {
          Rn(X_scores = X_scores,
             E_hat_scores = as.matrix(E_star_hat_scores[, , i]),
             ind_X_fpc = ind_X_fpc, ind_Y_fpc = ind_Y_fpc)})
        )

      # Create and decorate plot for the original process
      ylim <- range(lapply(Rn_processes, function(x) x(knots(x))))
      xlab <- substitute(paste(symbol("\xe1"), list(Chi, hat(Psi)[i]),
                               symbol("\xf1")), list(i = ind_X_fpc))
      ylab <- substitute(R[n](u, hat(Psi)[i], hat(Phi)[j]),
                         list(i = ind_X_fpc, j = ind_Y_fpc))
      plot(Rn_processes[[1]], lwd = 2, pch = NA, ylim = ylim, main = main,
           xlab = "", ylab = "")
      mtext(text = xlab, side = 1, line = 3, cex = 0.85)
      mtext(text = ylab, side = 2, line = 2.5, cex = 0.85)
      rug(knots(Rn_processes[[1]]))

      # Add bootstrap processes
      sapply(seq_len(dim(E_star_hat_scores)[3]), function(i) {
        plot(Rn_processes[[i + 1]], add = TRUE,
             col = gray(0.75, alpha = 0.75), pch = NA)
      })

      # Replot original process
      plot(Rn_processes[[1]], add = TRUE, lwd = 2, pch = NA)

    }

    # Reset par() for the user
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # Produce the pairs plot
    p_max <- min(length(p_hat), plot_max_p)
    q_max <- min(length(q_thre), plot_max_q)
    par(mfrow = c(p_max, q_max),
        mar = c(4, 4, 3, 2) + 0.1)
    for (i in p_hat[1:p_max]) {
      for (j in q_thre[1:q_max]) {

        plot_Rn(X_scores = fit_flm[["X_fpc"]][["scores"]],
                E_hat_scores = E_hat_scores,
                E_star_hat_scores = E_star_hat_scores,
                ind_X_fpc = i, ind_Y_fpc = j,
                main = substitute(Chi * " FPC " * i * ", " * Y * " FPC " * j,
                                  list(i = i, j = j)))

      }
    }
    par(mfrow = c(1, 1))

  }

  # Memory-expensive parts
  if (!save_fit_flm) {

    fit_flm <- NA

  }
  if (!save_boot_stats) {

    boot_stats <- NA
    boot_lambda <- NA
    boot_p_hat <- NA

  }

  # Return htest object
  result <- structure(list(statistic = c("statistic" = orig_stat),
                           p.value = p_value, boot_statistics = boot_stats,
                           method = meth, parameter = c("p" = p, "q" = q),
                           p = p_hat, q = q_thre, fit_flm = fit_flm,
                           boot_lambda = boot_lambda, boot_p = boot_p_hat,
                           data.name = data_name))
  class(result) <- "htest"
  return(result)

}
