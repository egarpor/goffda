

#' @title Sampling functional regression models with functional responses
#'
#' @description Simulation of a Functional Regression Model with Functional
#' Response (FRMFR) comprised of an additive mix of a linear and nonlinear
#' terms:
#' \deqn{Y(t) = \int_a^b X(s) \beta(s,t) ds + \Delta(X)(t) + \varepsilon(t),}{
#' Y(t) = \int_a^b X(s) \beta(s,t) ds + \Delta(X)(t) + \epsilon(t),}
#' where \eqn{X} is a random variable in the Hilbert space of
#' square-integrable functions in \eqn{[a, b]}, \eqn{L^2([a, b])},
#' \eqn{\beta} is the bivariate kernel of the FRMFR,
#' \eqn{\varepsilon}{\epsilon} is a random variable in \eqn{L^2([c, d])},
#' and \eqn{\Delta(X)}{\Delta(X)} is a nonlinear term.
#'
#' In particular, the scenarios considered in García-Portugués et al. (2019)
#' can be easily simulated.
#'
#' @param n sample size, only required when \code{scenario} is given.
#' @param scenario an index from \code{1} to \code{3} (default) denoting
#' one of the scenarios (S1, S2 or S3) simulated in
#' García-Portugués et al. (2019) (see details below). If
#' \code{scenario = NULL}, \code{X_fdata}, \code{error_fdata}, and \code{beta}
#' have to be provided. Otherwise, \code{X_fdata}, \code{error_fdata}, and
#' \code{beta} will be ignored.
#' @param X_fdata sample of functional covariates \eqn{X(s)} as
#' \code{\link[fda.usc]{fdata}} objects of length \code{n}, with \eqn{s} in
#' \eqn{[a, b]}. Defaults to \code{NULL}.
#' @param error_fdata sample of functional errors \eqn{\varepsilon(t)} as
#' \code{\link[fda.usc]{fdata}} objects of length \code{n}, with \eqn{t} in
#' \eqn{[c, d]}. If \code{concurrent = TRUE}, \code{X_fdata} and
#' \code{error_fdata} must be valued in the same grid. Defaults to \code{NULL}.
#' @param beta matrix containing the values \eqn{\beta(s, t)}, for each grid
#' point \eqn{s} in \eqn{[a, b]} and \eqn{t} in \eqn{[c, d]}. If
#' \code{concurrent = TRUE} (see details below), a row/column vector
#' must be introduced, valued in the same grid as \code{error_fdata}.
#' If \code{beta = NULL} (default), \code{scenario != NULL} is required.
#' @param std_error standard deviation of the random variables
#' involved in the generation of the functional error \code{error_fdata}.
#' Defaults to \code{0.15}.
#' @param s,t grid points. If \code{X_fdata}, \code{error_fdata} and
#' \code{beta} are provided, \code{s} and \code{t} are ignored. Default to
#' \code{s = seq(0, 1, l = 101)} and \code{t = seq(0, 1, l = 101)},
#' respectively.
#' @param n_fpc number of components to be considered for the generation of
#' functional variables. Defaults to \code{50}.
#' @param nonlinear nonlinear term. Either a character string (\code{"exp"},
#' \code{"quadratic"} or \code{"sin"}) or an \code{\link[fda.usc]{fdata}}
#' object of length \code{n}, valued in the same grid as \code{error_fdata}.
#' If \code{nonlinear = NULL} (default), the nonlinear part is set to zero.
#' @param concurrent flag to consider a concurrent FLRFR (degenerate case).
#' Defaults to \code{FALSE}.
#' @inheritParams quadrature
#' @inheritParams flm_term
#' @param verbose flag to display information about the sampling procedure.
#' Defaults to \code{FALSE}.
#' @param ... further parameters passed to
#' \code{\link[goffda]{r_cm2013_flmfr}}, \code{\link[goffda]{r_gof2019_flmfr}}
#' and\cr \code{\link[goffda]{r_ik2018_flmfr}}, depending on the
#' chosen \code{scenario}.
#' @return A list with the following elements:
#' \item{\code{X_fdata}}{functional covariates, an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#' \item{\code{Y_fdata}}{functional responses, an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#' \item{\code{error_fdata}}{functional errors, an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#' \item{\code{beta}}{either the matrix with \eqn{\beta(s, t)} evaluated at
#' the \code{argvals} of \code{X_fdata} and \code{Y_fdata} (if
#' \code{concurrent = FALSE}) or a vector with \eqn{\beta(t)}
#' evaluated at the \code{argvals} of \code{X_fdata} (if
#' \code{concurrent = TRUE}).}
#' \item{\code{nl_dev}}{nonlinear term, an \code{\link[fda.usc]{fdata}}
#' object of length \code{n}.}
#' @details
#' \itemize{
#'   \item{\code{r_frm_fr} samples the above regression model,
#'   where the nonlinear term \eqn{\Delta(X)} is computed by \code{nl_dev}.
#'   Functional covariates, errors, and \eqn{\beta} are generated
#'   automatically from the scenarios in García-Portugués et al. (2019) when
#'   \code{scenario != NULL} (see the documentation of
#'   \code{\link{r_gof2019_flmfr}}). If \code{scenario = NULL},
#'   covariates, errors and \eqn{\beta} must be provided.
#'
#'   When \code{concurrent = TRUE}, the concurrent FRMFR
#'   \deqn{Y(t) = X(t) \beta(t) +
#'   \Delta(X)(t) + \varepsilon(t)}{Y(t) = X(t) * \beta(t)
#'   + \Delta(X)(t) + \epsilon(t)}
#'   is considered.}
#'   \item{\code{nl_dev} computes a nonlinear deviation
#'   \eqn{\Delta(X)}{\Delta(X)}:
#'   \eqn{\exp(\sqrt{X(a + (t - c) ((b - a) / (d - c)))})}{
#'   exp(\sqrt(X(a + (t - c) * ((b - a) / (d - c)))))}
#'   (for \code{"exp"}),
#'   \eqn{(X^2 (a + (t - c) ((b - a) / (d - c))) - 1)}{
#'   (X^2 * (a + (t - c) * ((b - a) / (d - c))) - 1)}
#'   (\code{"quadratic"}) or
#'   \eqn{(\sin(2\pi t) - \cos(2 \pi t)) \| X \|^2}{
#'   (sin(2 \pi t) - cos(2 \pi t)) * ||X||^2}
#'   (\code{"sin"}). Also, \eqn{\Delta(X)} can be manually set as an
#'   \code{\link[fda.usc]{fdata}} object of length \code{n} and valued in
#'   the same grid as \code{error_fdata}.}
#' }
#' @examples
#' ## Generate samples for the three scenarios
#'
#' # Equispaced grids and Simpson's rule
#'
#' s <- seq(0, 1, l = 101)
#' samp <- list()
#' old_par <- par(mfrow = c(3, 5))
#' for (i in 1:3) {
#'   samp[[i]] <- r_frm_fr(n = 100, scenario = i, s = s, t = s,
#'                         int_rule = "Simpson")
#'   plot(samp[[i]]$X_fdata)
#'   plot(samp[[i]]$error_fdata)
#'   plot(samp[[i]]$Y_fdata)
#'   plot(samp[[i]]$nl_dev)
#'   image(x = s, y = s, z = samp[[i]]$beta, col = viridisLite::viridis(20))
#' }
#' par(old_par)
#'
#' ## Linear term as a concurrent model
#'
#' # The grids must be have the same number of grid points for a given
#' # nonlinear term and a given beta function
#'
#' s <- seq(1, 2, l = 101)
#' t <- seq(0, 1, l = 101)
#' samp_c_1 <- r_frm_fr(n = 100, scenario = 3, beta = sin(t) - exp(t),
#'                      s = s, t = t, nonlinear = fda.usc::fdata(mdata =
#'                        t(matrix(rep(sin(t), 100), nrow = length(t))),
#'                        argvals = t),
#'                      concurrent = TRUE)
#' old_par <- par(mfrow = c(3, 2))
#' plot(samp_c_1$X_fdata)
#' plot(samp_c_1$error_fdata)
#' plot(samp_c_1$Y_fdata)
#' plot(samp_c_1$nl_dev)
#' plot(samp_c_1$beta)
#' par(old_par)
#'
#' ## Sample for given X_fdata, error_fdata, and beta
#'
#' # Non equispaced grids with sinusoidal nonlinear term and intensity 0.5
#' s <- c(seq(0, 0.5, l = 50), seq(0.51, 1, l = 101))
#' t <- seq(2, 4, len = 151)
#' X_fdata <- r_ou(n = 100, t = s, alpha = 2, sigma = 4, x0 = 1:100)
#' error_fdata <- r_ou(n = 100, t = t, alpha = 1, sigma = 1, x0 = 1:100)
#' beta <- r_gof2019_flmfr(n = 100, s = s, t = t)$beta
#' samp_Xeps <- r_frm_fr(scenario = NULL, X_fdata = X_fdata,
#'                       error_fdata = error_fdata, beta = beta,
#'                       nonlinear = "exp", int_rule = "trapezoid")
#' old_par <- par(mfrow = c(3, 2))
#' plot(samp_Xeps$X_fdata)
#' plot(samp_Xeps$error_fdata)
#' plot(samp_Xeps$Y_fdata)
#' plot(samp_Xeps$nl_dev)
#' image(x = s, y = t, z = beta, col = viridisLite::viridis(20))
#' par(old_par)
#' @author Javier Álvarez-Liébana.
#' @references
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' Gonzalez-Manteiga, W. (2019). A goodness-of-fit test
#' linear model with functional response. \emph{arXiv:1909.07686}.
#' \url{https://arxiv.org/abs/1909.07686}
#' @name sim-frmfr


#' @rdname sim-frmfr
#' @export
r_frm_fr <- function(n, scenario = 3, X_fdata = NULL, error_fdata = NULL,
                     beta = NULL, s = seq(0, 1, l = 101),
                     t = seq(0, 1, l = 101), std_error = 0.15, nonlinear = NULL,
                     concurrent = FALSE, int_rule = "trapezoid",
                     n_fpc = 50, verbose = FALSE, ...) {

  # If scenario != NULL, some of the scenarios in
  # García-Portugués et al. (2019) will be generated
  if (!is.null(scenario)) {

    # If scenario is not NULL, X_fdata = NULL, error_fdata = NULL and
    # beta  = NULL are forced
    if (verbose & (!is.null(X_fdata) | !is.null(error_fdata) |
                   !is.null(beta))) {

      message("Scenario encoded as scenario = ", scenario,
              ". X_fdata, beta and error_fdata will be ignored")

    }
    X_fdata <- error_fdata <- beta <- NULL

    # Check if the scenario is implemented
    if (!(scenario %in% 1:3)) {

      stop(paste("Scenario encoded by", scenario, "is not implemented. Only",
                 "scenarios 1 (Crambes and Mas, 2013),",
                 "2 (Garcia-Portugues et al., 2019) and",
                 "3 (Imaizumi and Kato, 2018) are directly implemented"))

    }

    # Check if s and t are intervals
    if (!is.vector(s) | length(s) < 1) {

      stop("s must be a vector with length greater than zero")

    }
    if (!is.vector(t) | length(t) < 1) {

      stop("t must be a vector with length greater than zero")

    }

  } else { # If scenario = NULL, X_fdata, error_fdata and beta are required

    # Check if functional variables are properly provided
    if (!fda.usc::is.fdata(X_fdata) | !fda.usc::is.fdata(error_fdata)) {

      stop(paste("scenario = NULL: \"fdata\" objects must be provided",
                 "as X_fdata and error_fdata, and beta must be a matrix"))

    }

    # Basic info of data
    n <- nrow(X_fdata[["data"]])
    s <- X_fdata[["argvals"]]
    t <- error_fdata[["argvals"]]

    # Check sample sizes
    if (n != nrow(error_fdata[["data"]])) {

      stop("The sample sizes of X_fdata and error_fdata do not match")

    }

    # Check concurrent model
    if (concurrent) {

      # Check if grids are the same
      if (length(s) != length(t)) {

        stop(paste("If concurrent model is generated, grid intervals of",
                   "X_fdata and error_fdata must be the same number of points"))

      }

      # Check if beta is a vector with properly length
      if (!is.vector(beta) |
          length(beta) != length(t)) {

        stop(paste("If concurrent model is generated, beta must be a vector",
                   "with the same length as error_fdata"))

      }

    } else {

      # Check if beta has properly dimensions
      if (!(nrow(beta) == length(s) &
            ncol(beta) == length(t) & is.matrix(beta))) {

        stop(paste("scenario = NULL: beta must be a matrix whose",
                   "dimensions are c(length(X_fdata$argvals),",
                   "length(error_fdata$argvals))"))


      }

    }

  }

  # Sampling adopting one of the scenarios considered in
  # García-Portugués et al. (2019)
  if (!is.null(scenario)) {

    # S1: based on the data simluated by Crambes and Mas 2013 (Section 3)
    if (scenario == 1) {

      # Sampling functional covariates and functional errors
      cm2013 <- r_cm2013_flmfr(n = n, s = s, t = t, std_error = std_error,
                               n_fpc = n_fpc, concurrent = concurrent)
      X_fdata <- cm2013[["X_fdata"]]
      error_fdata <- cm2013[["error_fdata"]]

      # beta(s,t) (or beta(t) for the concurrent model) valued in grid points
      beta <- cm2013[["beta"]]

    }

    # S2: example of García-Portugués et al. 2019 (Sections 2 and 4)
    if (scenario == 2) {

      # Sampling functional covariates and functional errors
      gof2019 <- r_gof2019_flmfr(n = n, s = s, t = t, std_error = std_error,
                                 concurrent = concurrent)
      X_fdata <- gof2019[["X_fdata"]]
      error_fdata <- gof2019[["error_fdata"]]

      # beta(s,t) (or beta(t) for the concurrent model) valued in grid points
      beta <- gof2019[["beta"]]

    }

    # S3: based on the data simulated by Imaizumi and Kato 2018 (Section 4)
    if (scenario == 3) {

      # Sampling functional covariates and functional errors
      ik2018 <- r_ik2018_flmfr(n = n, s = s, t = t, std_error = std_error,
                               n_fpc = n_fpc, concurrent = concurrent, ...)
      X_fdata <- ik2018[["X_fdata"]]
      error_fdata <- ik2018[["error_fdata"]]

      # beta(s,t) (or beta(t) for the concurrent model) valued in grid points
      beta <- ik2018[["beta"]]

    }

  }

  # Check if X_fdata is equispaced
  eps <- sqrt(.Machine[["double.eps"]])
  equispaced_x <- all(abs(diff(s, differences = 2)) < eps)

  # Sampling the linear term generated by X_fdata, error_fdata and beta
  linear_fdata <- flm_term(X_fdata = X_fdata, beta = beta, t = t,
                           int_rule = int_rule, equispaced = equispaced_x,
                           concurrent = concurrent)

  # Sampling deviations from the linearity
  # Check for nonlinear term
  if (!is.null(nonlinear) & !is.character(nonlinear)) {

    if (!fda.usc::is.fdata(nonlinear)) {

      stop(paste("If nonlinear is not NULL neither \"exp\", \"quadratic\"",
                 "or \"sin\", it must be an \"fdata\" class object"))

    } else {

      # Check if dimensions match
      if (nrow(nonlinear[["data"]]) != n |
               ncol(nonlinear[["data"]]) != length(t)) {

        stop(paste("If nonlinear is not NULL neither \"exp\", \"quadratic\"",
                   "or \"sin\", it must be an \"fdata\" class object with",
                   "the same dimensions as error_fdata"))

      } else {

        # Nonlinear term is given as input
        d_nl_dev <- nonlinear

      }

    }

  } else {

    # Nonlinear term providing by nl_dev function. If nonlinear = NULL,
    # d_nl_dev just stores a null fdata
    d_nl_dev <- switch(is.null(nonlinear) + 1, nl_dev(X_fdata = X_fdata, t = t,
                                                      nonlinear = nonlinear,
                                                      int_rule = int_rule,
                                                      equispaced = equispaced_x,
                                                      verbose = verbose),
                       fda.usc::fdata(mdata = matrix(0, nrow = n,
                                                     ncol = length(t)),
                                      argvals = t))

  }

  # Generative model
  Y_fdata <- linear_fdata + d_nl_dev + error_fdata

  # Output
  return(list("X_fdata" = X_fdata, "Y_fdata" = Y_fdata,
              "error_fdata" = error_fdata, "beta" = beta,
              "linear_fdata" = linear_fdata, "nl_dev" = d_nl_dev))

}


#' @rdname sim-frmfr
#' @export
nl_dev <- function(X_fdata, t = seq(0, 1, l = 101), nonlinear = NULL,
                   int_rule = "trapezoid", equispaced = equispaced,
                   verbose = FALSE) {

  # Check X_fdata
  if (!fda.usc::is.fdata(X_fdata)) {

    stop("X_fdata must be an \"fdata\" class object")

  }

  # Check grid points
  if (!is.vector(t) | length(t) < 1) {

    stop("Grid points in t must be a vector of length greater than zero")

  }

  # Check nonlinear term
  if (is.null(nonlinear)) {

    if (verbose) {

      message("Nonlinear = NULL: nonlinear term will be zero")

    }

    # If nonlinear = NULL, we just return a null fdata
    return(fda.usc::fdata(mdata = matrix(0, nrow(X_fdata[["data"]]), length(t)),
                          argvals = t))

  } else if (!(nonlinear %in% c("exp", "quadratic", "sin"))) {

    stop(paste("Only exponential, quadratic or sinusoidal nonlinear terms",
               "have been implemented"))

  } else {

    # If lengths of intervals is not the same, just rescaling is not possible
    # An interpolation is required
    if (nonlinear != "sin" & (length(X_fdata[["argvals"]]) != length(t))) {

      X_fdata_int <- fda.usc::fdata(mdata =
                                    matrix(0, nrow = nrow(X_fdata[["data"]]),
                                           ncol = length(t)), argvals = t)
      for (i in 1:nrow(X_fdata[["data"]])) {

        X_fdata_int[["data"]][i, ] <- spline(x = X_fdata[["argvals"]],
                                             y = X_fdata[["data"]][i, ],
                                             n = length(t))[["y"]]

      }

      if (verbose) {

        message("When exponential or quadratic nonlinear terms are considered,",
                " X_fdata must be valued in a grid with the same number of",
                " points as t: an interpolation will be performed")

      }
    } else {

      X_fdata_int <- X_fdata

    }

    # Output nonlinear term
    nonlinear_term <- switch((nonlinear == "quadratic") +
                               2 * (nonlinear == "sin") +  1,
                      fda.usc::fdata(mdata = exp(sqrt(X_fdata_int[["data"]])),
                                     argvals = t),
                      fda.usc::fdata(mdata = X_fdata_int[["data"]]^2 - 1,
                                     argvals = t),
                      fda.usc::fdata(mdata = outer(apply(X_fdata[["data"]]^2, 1,
                                                         FUN = "integral1D",
                                                         X_fdata[["argvals"]],
                                                         int_rule, equispaced),
                                                   sin(2 * pi * t) -
                                                     cos(2 * pi * t), "*"),
                                     argvals = t))

  }

  # Output
  return(nonlinear_term)

}


#' @title Covariate, error, and kernel of a functional linear model
#' with functional response
#'
#' @description Simulation of \eqn{X}, a random variable in the Hilbert space
#' of square-integrable functions in \eqn{[a, b]}, \eqn{L^2([a, b])}, and
#' \eqn{\varepsilon}{\epsilon}, a random variable in \eqn{L^2([c, d])}.
#' Together with the bivariate kernel \eqn{\beta}, they are the necessary
#' elements for sampling a Functional Linear Model with  Functional Response
#' (FLMFR):
#' \deqn{Y(t) = \int_a^b X(s) \beta(s,t) ds + \varepsilon(t).}{
#' Y(t) = \int_a^b X(s) \beta(s,t) ds + \epsilon(t).}
#'
#' The next functions sample \eqn{X} and \eqn{\varepsilon}{\epsilon}, and
#' construct \eqn{\beta}, using different proposals in the literature:
#' \itemize{
#'   \item{\code{r_cm2013_flmfr} is based on the numerical example given in
#'   Section 3 of Crambes and Mas (2013). Termed as S1 in Section 2 of
#'   García-Portugués et al. (2019).}
#'   \item{\code{r_ik2018_flmfr} is based on the numerical example given in
#'   Section 4 of Imaizumi and Kato (2018), but zeroing the first Functional
#'   Principal Components (FPC) coefficients of \eqn{\beta} (so the first FPC
#'   are not adequate for estimation). S3 in Section 2 of
#'   García-Portugués et al. (2019).}
#'   \item{\code{r_gof2019_flmfr} gives a numerical example in Section 2
#'   of García-Portugués et al. (2019), denoted therein as S2.}
#' }
#'
#' @param s,t grid points where functional covariates and responses are valued,
#' respectively.
#' @inheritParams r_ou
#' @inheritParams sim-frmfr
#' @param parameters vector of parameters, only required for
#' \code{r_ik2018_flmfr}. Defaults to\cr \code{c(1.75, 0.8, 2.4, 0.25)}.
#' @param n_fpc number of FPC to be taken into account for the data generation.
#' Must be greater than \code{4} when \code{r_ik2018_flmfr} is applied, since
#' the first \eqn{4} FPC are null. Defaults to \code{50}.
#' @param concurrent flag to consider a concurrent FLMFR (degenerate case).
#' Defaults to \code{FALSE}.
#' @return A list with the following elements:
#' \item{\code{X_fdata}}{functional covariates, an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#' \item{\code{error_fdata}}{functional errors, an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#' \item{\code{beta}}{either the matrix with \eqn{\beta(s, t)} evaluated at
#' the \code{argvals} of \code{X_fdata} and \code{Y_fdata} (if
#' \code{concurrent = FALSE}) or a vector with \eqn{\beta(t)}
#' evaluated at the \code{argvals} of \code{X_fdata} (if
#' \code{concurrent = TRUE}).}
#' @details
#' Descriptions of the processes \eqn{X} and \eqn{\varepsilon}{\epsilon},
#' and of \eqn{\beta} can be seen in the references.
#' @examples
#' # FLMFR based on Imaizumi and Kato (2018) adopting different Hilbert spaces
#' s <- seq(0, 1, l = 201)
#' t <- seq(2, 4, l = 301)
#' r_ik2018 <- r_ik2018_flmfr(n = 50, s = s, t = t, std_error = 1.5,
#'                            parameters = c(1.75, 0.8, 2.4, 0.25), n_fpc = 50)
#' plot(r_ik2018$X_fdata)
#' plot(r_ik2018$error_fdata)
#' image(x = s, y = t, z = r_ik2018$beta, col = viridisLite::viridis(20))
#'
#' # FLMFR based on Cardot and Mas (2013) adopting different Hilbert spaces
#' r_cm2013 <- r_cm2013_flmfr(n = 50, s = s, t = t, std_error = 0.15,
#'                            n_fpc = 50)
#' plot(r_cm2013$X_fdata)
#' plot(r_cm2013$error_fdata)
#' image(x = s, y = t, z = r_cm2013$beta, col = viridisLite::viridis(20))
#'
#' # FLMFR in García-Portugués et al. (2019) adopting different Hilbert spaces
#' r_gof2019 <- r_gof2019_flmfr(n = 50, s = s, t = t, std_error = 0.35,
#'                              concurrent = FALSE)
#' plot(r_gof2019$X_fdata)
#' plot(r_gof2019$error_fdata)
#' image(x = s, y = t, z = r_gof2019$beta, col = viridisLite::viridis(20))
#'
#' # Concurrent model in García-Portugués et al. (2019)
#' r_gof2019 <- r_gof2019_flmfr(n = 50, s = s, t = s, std_error = 0.35,
#'                              concurrent = TRUE)
#' plot(r_gof2019$X_fdata)
#' plot(r_gof2019$error_fdata)
#' plot(r_gof2019$beta)
#' @author Javier Álvarez-Liébana.
#' @references
#' Cardot, H. and Mas, A. (2013). Asymptotics of prediction in functional linear
#' regression with functional outputs. \emph{Bernoulli}, 19(5B):2627--2651.
#' \url{https://doi.org/10.3150/12-BEJ469}
#'
#' Imaizumi, M. and Kato, K. (2018). PCA-based estimation for functional linear
#' regression with functional responses. \emph{Journal of Multivariate
#' Analysis}, 163:15--36. \url{https://doi.org/10.1016/j.jmva.2017.10.001}
#'
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' Gonzalez-Manteiga, W. (2019). A goodness-of-fit test
#' linear model with functional response. \emph{arXiv:1909.07686}.
#' \url{https://arxiv.org/abs/1909.07686}
#' @name elem-flmfr


#' @rdname elem-flmfr
#' @export
r_cm2013_flmfr <- function(n, s = seq(0, 1, len = 101),
                           t = seq(0, 1, len = 101), std_error = 0.15,
                           n_fpc = 50, concurrent = FALSE) {

  # Check for standard deviation of errors
  if (std_error < 0) {

    stop("Standard deviation of error must be positive")

  }

  # Check for sample size
  if (n < 1) {

    stop("Sample size n must be greater than zero")

  }

  # Check n_fpc
  if (n_fpc < 1) {

    stop("Number of components (n_fpc) must be greater or equal than one")

  }

  # Bivariate kernel function beta(s,t) as an uniform surface for non
  # concurrent model and sqrt(|sin(pi * t) - cos(pi * t)|) for concurrent models
  beta <- switch(concurrent + 1,  outer((s - s[1])^2, (t - t[1])^2, "+"),
                 sqrt(abs(sin(pi * t) - cos(pi * t))))

  # Functional covariate based on Crambes and Mas (2013), in terms of the
  # Karhunen Loeve expansion
  X_fdata <- fda.usc::fdata((t(replicate(n, 1 / ((pi^2) * (1:n_fpc - 0.5)^2))) *
                               matrix(rnorm(n * n_fpc, mean = 0, sd = 2),
                                      nrow = n)) %*%
                              (sqrt(2) * sin((1:n_fpc - 0.5) * pi *
                                             t(matrix(rep(s, n_fpc),
                                                      nrow = length(s))))),
                              argvals = s)

  # Functional error as Brownian motion with sigma equal to std_error
  error_fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian",
                                      par.list = list("scale" = std_error^2))

  # Output
  return(list("X_fdata" = X_fdata, "error_fdata" = error_fdata, "beta" = beta))

}


#' @rdname elem-flmfr
#' @export
r_ik2018_flmfr <- function(n, s = seq(0, 1, l = 101), t = seq(0, 1, l = 101),
                           std_error = 1.5,
                           parameters = c(1.75, 0.8, 2.4, 0.25),
                           n_fpc = 50, concurrent = FALSE) {

  # Check for standard deviation of errors
  if (std_error < 0) {

    stop("Standard deviation of error must be positive")

  }

  # Check for sample size
  if (n < 1) {

    stop("Sample size n must be greater than zero")

  }

  # Check n_fpc
  if (n_fpc < 4) {

    stop(paste("Number of functional components must be greater than 4 when",
               "the example based on Imaizumi and Kato (2018) is implemented"))

  }

  # X_fdata generated as a basis expansion in terms of uniform
  # distributions Uk ~ U(-sqrt(3),sqrt(3))
  Psi <- rbind(rep(1, length(s)), sqrt(2) *
                 cos((1:(n_fpc - 1)) * pi * t(matrix(rep(s, n_fpc - 1),
                                                     nrow = length(s)))))
  X_fdata <- fda.usc::fdata((t(replicate(n, (1:n_fpc)^(-parameters[1]))) *
                               matrix(runif(n * n_fpc, min = -sqrt(5),
                                            max = sqrt(5)), nrow = n)) %*% Psi,
                            argvals = s)

  # Functional errors as a basis expansin in terms of gaussian distributions
  Phi <- rbind(rep(1, length(t)),
               sqrt(2) * cos((1:(n_fpc - 1)) * pi *
                               t(matrix(rep(t, n_fpc - 1), nrow = length(t)))))

  error_fdata <- fda.usc::fdata((t(replicate(n, (1:n_fpc)^(-parameters[2]))) *
                                   matrix(rnorm(n * n_fpc, mean = 0,
                                                sd = std_error), nrow = n)) %*%
                                  Phi, argvals = t)

  # Bivariate kernel function beta(s,t)

  # Coefficients of beta based on the beta proposed by Imaizumi and Kato (2018),
  # but zeroing the first 4x4 coefficients
  coef_beta <- matrix(0, nrow = n_fpc, ncol = n_fpc)
  null_rows <- 4
  coef_beta[(null_rows + 1):n_fpc, (null_rows + 1):n_fpc] <-
    (6 * (-1)^(outer((null_rows + 1):n_fpc, (null_rows + 1):n_fpc, "+"))) *
    outer((1:(n_fpc - null_rows))^(-parameters[3]),
          (1:(n_fpc - null_rows))^(-parameters[4]), "*")

  # Bivariate kernel function beta(s,t) as surface in terms of the generated
  # coefficients for non concurrent model and (t - 0.5)^3 for concurrent models
  beta <- switch(concurrent + 1,  t(Psi) %*% coef_beta %*% Phi,
                 (t - 0.5)^3)

  # Output
  return(list("X_fdata" = X_fdata, "error_fdata" = error_fdata, "beta" = beta))

}


#' @rdname elem-flmfr
#' @export
r_gof2019_flmfr <- function(n, s = seq(0, 1, len = 101),
                            t = seq(0, 1, len = 101), std_error = 0.35,
                            concurrent = FALSE) {

  # Check for standard deviation of errors
  if (std_error < 0) {

    stop("Standard deviation of error must be positive")

  }

  # Check for sample size
  if (n < 1) {

    stop("Sample size n must be greater than zero")

  }

  # Bivariate kernel function beta(s,t) as an egg carton shape for non
  # concurrent model and log(t - t[1] - 0.5) for concurrent models
  beta <- switch(concurrent + 1, 2 * outer(sin(6 * pi * (s - s[1])),
                                           cos(6 * pi * (t - t[1])), "+"),
                 log(t - t[1] + 0.5))

  # Functional covariate as zero-mean Gaussian process with exponential
  # variogram (scale = 6)
  X_fdata <- fda.usc::rproc2fdata(n = n, t = s, sigma = "vexponential",
                                  par.list = list("scale" = 6))

  # Functional error as OU process with unitary drift and stationary std
  # equal to std_error
  error_fdata <- r_ou(n = n, t = t, sigma = sqrt(2) * std_error)

  # Output
  return(list("X_fdata" = X_fdata, "error_fdata" = error_fdata, "beta" = beta))

}


#' @title Functional linear model term with bivariate kernel
#'
#' @description Computation of the functional linear term
#' \deqn{\int_a^b \beta(s, t) X(s)\,\mathrm{d}s,}{\int_a^b \beta(s, t) X(s) ds,}
#' of a Functional Linear Model with Functional Response (FLMFR), where
#' \eqn{X} is a random variable in the Hilbert space of
#' square-integrable functions in \eqn{[a, b]}, \eqn{L^2([a, b])},
#' \eqn{\beta} is the bivariate kernel of the FLMFR, and
#' \eqn{\varepsilon}{\epsilon} is a random variable in \eqn{L^2([c, d])}.
#'
#' @param X_fdata sample of functional data as an
#' \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @param beta matrix containing the values  \eqn{\beta(s, t)},
#' for each grid point \eqn{s} in  \eqn{[a, b]} and \eqn{t} in \eqn{[c, d]}. If
#' \code{concurrent = TRUE}, a row/column vector must be introduced, valued in
#' the same grid as \code{error_fdata}, with the same length as
#' \code{length(X_fdata$argvals)}.
#' @inheritParams quadrature
#' @inheritParams elem-flmfr
#' @return Functional linear model term as the integral (in \code{s}) between
#' \code{X_fdata} and \code{beta}, as an \code{\link[fda.usc]{fdata}} object of
#' length \code{n}.
#' @examples
#' ## Generate a sample of functional responses via FLMFR
#'
#' # Bivariate kernel beta(s,t) as an egg carton shape
#' s <- seq(0, 1, l = 101)
#' t <- seq(0, 1, l = 201)
#' beta <- outer(s, t, FUN = function(s, t) sin(6 * pi * s) + cos(6 * pi * t))
#'
#' # Functional data as zero-mean Gaussian process with exponential variogram
#' X_fdata <- fda.usc::rproc2fdata(n = 50, t = s, sigma = "vexponential",
#'                                 list = list(scale = 2.5))
#'
#' # Functional error as an OU process with variance 0.075
#' sigma <- sqrt(0.075) * 2
#' error_fdata <- r_ou(n = 50, t = t, sigma = sigma)
#' Y_fdata <- flm_term(X_fdata = X_fdata, beta = beta, t = t) + error_fdata
#' plot(Y_fdata)
#'
#' ## Generate a sample of functional responses via concurrent model
#'
#' # Function beta(t)
#' s <- seq(1, 3, l = 201)
#' t <- seq(2, 5, l = 201)
#' beta <- sin(pi * t) + cos(pi * t)
#'
#' # Functional data as zero-mean Gaussian process with exponential variogram
#' X_fdata <- fda.usc::rproc2fdata(n = 50, t = s, sigma = "vexponential",
#'                                 list = list(scale = 2.5))
#'
#' # Functional error as an OU process with variance 0.075
#' sigma <- sqrt(0.075) * 2
#' error_fdata <- r_ou(n = 50, t = t, sigma = sigma)
#' Y_fdata <- flm_term(X_fdata = X_fdata, beta = beta, t = t,
#'                     concurrent = TRUE) + error_fdata
#' plot(Y_fdata)
#' @author Javier Álvarez-Liébana.
#' @references
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' González-Manteiga, W. (2019). A goodness-of-fit test for the functional
#' linear model with functional response. \emph{arXiv:1909.07686}.
#' \url{https://arxiv.org/abs/1909.07686}
#' @export
flm_term <- function(X_fdata, beta, t, int_rule = "trapezoid",
                     equispaced = NULL, concurrent = FALSE) {

  # Basic info of X_fdata and beta
  n <- dim(X_fdata[["data"]])[1]
  s <- X_fdata[["argvals"]]

  # Check if X_fdata, error_fdata and beta have proper dimensions
  if (!concurrent) {

    # Check dimensions
    dim_beta <- dim(beta)
    if (length(s) != dim_beta[1]) {

      stop(paste("Matrix beta must have as number of rows the number",
                 "of grid points where functional regressors are valued"))

    }
    if (length(t) != dim_beta[2]) {

      stop(paste("Matrix beta must have as number of columns the number",
                 "of grid points stored in t"))

    }

  } else {

    # As vector
    beta <- as.vector(beta)

    # Check that dimensions of beta match lengths of s and t
    if ((length(s) != length(t)) | (length(t) != length(beta))) {

      stop(paste("When concurrent model is considered, X_fdata and beta",
                 "must be valued in the same grid interval t"))

    }

  }

  # Concurrent model is just generated: different intervals with the same
  # number of grid points are allowed
  if (concurrent) {

    return(fda.usc::fdata(mdata = t(t(X_fdata[["data"]]) * beta),
                          argvals = t))
  } else {

    # Check if X_fdata is equispaced
    if (is.null(equispaced)) {

      # Check if X_fdata is equispaced
      eps <- sqrt(.Machine[["double.eps"]])
      equispaced <- all(abs(diff(s, differences = 2)) < eps)

    }

    # Y_fdata is computed, for each sample element i = 1:n, as the linear
    # operator given by the integral of X(s) * beta(s,t)ds, such that this
    # operator is valued in t. Loop is avoided by using "apply" function

    # Output
    return(matrix(apply(t(beta)[rep(1:dim_beta[2], each = n), ] *
                          X_fdata[["data"]][rep(1:n, times = dim_beta[2]), ],
                        1, FUN = integral1D, t = s,
                        int_rule = int_rule,
                        equispaced = equispaced), nrow = n))

  }

}


#' @title Simulation of an Ornstein--Uhlenbeck process
#'
#' @description Simulation of trajectories of the Ornstein--Uhlenbeck process
#' \eqn{\{X_t\}}{{X_t}}. The process is the solution to the stochastic
#' differential equation
#' \deqn{\mathrm{d}X_t = \alpha (X_t - \mu)\mathrm{d}t + \sigma \mathrm{d}W_t,
#' }{dX_t = \alpha (X_t - \mu) dt + \sigma dW_t,}
#' whose stationary distribution is \eqn{N(\mu, \sigma^2 / (2 \alpha))}, for
#' \eqn{\alpha, \sigma > 0} and \eqn{\mu \in R}.
#'
#' Given an initial point \eqn{x_0} and the evaluation times
#' \eqn{t_1, \ldots, t_m}, a sample trajectory \eqn{X_{t_1}, \ldots, X_{t_m}}
#' can be obtained by sampling the joint Gaussian distribution of
#' \eqn{(X_{t_1}, \ldots, X_{t_m})}.
#'
#' @param n number of trajectories to sample.
#' @param t evaluation times for the trajectories, a vector.
#' @param alpha strength of the drift, a positive scalar.
#' @param mu mean of the process, a scalar.
#' @param sigma diffusion coefficient, a positive scalar.
#' @param x0 a vector of length \code{n} giving the initial
#' values of the Ornstein--Uhlenbeck trajectories. By default, \code{n}
#' points are sampled from the stationary distribution. If a single scalar
#' is passed, then the same \code{x0} is employed for all the trajectories.
#' @return Random trajectories, an \code{\link[fda.usc]{fdata}} object of
#' length \code{n} and \code{t} as \code{argvals}.
#' @examples
#' # Same initial point
#' plot(r_ou(n = 20, x0 = 5), col = viridisLite::viridis(20))
#'
#' # Different initial points
#' plot(r_ou(n = 100, alpha = 2, sigma = 4, x0 = 1:100),
#'      col = viridisLite::viridis(100))
#' @author Eduardo García-Portugués.
#' @export
r_ou <- function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1,
                 x0 = rnorm(n, mean = mu, sd = sigma / sqrt(2 * alpha))) {

  # Check for drift
  if (alpha < 0) {

    stop("Strength of the drift must be positive")

  }

  # Check for diffusion coefficient
  if (sigma < 0) {

    stop("Diffusion coefficient must be positive")

  }

  # Check for the initial conditions
  lx0 <- length(x0)
  if (lx0 == 1) {

    x0 <- rep(x0, n)

  } else {

    if (lx0 != floor(n)) {

      stop(paste0("Length of vector of initial conditions x0 is ", lx0,
                  ". It must be 1 or the same as the sample size n = ", n))

    }

  }

  # Time-varying covariance matrix
  St <- (sigma^2) / (2 * alpha) * outer(t, t, function(s, t) {
    exp(alpha * (2 * pmin(s, t) - (s + t))) - exp(-alpha * (s + t))
  })

  # Cholesky decomposition cannot be computed in t = 0: we remove the
  # associated rows and columns, compute the Cholesky decomposition, and then
  # include zero vectors
  lt <- length(t)
  t_zero <- which(t == 0)
  if (length(t_zero) > 0) {

    R <- matrix(0, nrow = lt, ncol = lt)
    R[-t_zero, -t_zero] <- chol(St[-t_zero, -t_zero])

  } else {

    R <- chol(St)

  }

  # Sample paths of an OU process from a multivariate N(0, St), using the
  # Cholesky decomposition of St
  OU_paths <- matrix(rnorm(n = n * lt), nrow = n, ncol = lt,
                     byrow = TRUE) %*% R +
    outer(x0, t, function(x0, t) mu + (x0 - mu) * exp(-alpha * t))

  # As fdata object
  return(fda.usc::fdata(mdata = OU_paths, argvals = t))

}
