

#' @title Fitting of regularized linear models
#'
#' @description Convenience function for fitting multivariate linear models
#' with multivariate response by relying on \code{\link[glmnet]{cv.glmnet}}
#' from the \code{\link[glmnet]{glmnet-package}}. The function fits the
#' multivariate linear model
#' \deqn{\mathbf{Y} = \mathbf{X}\mathbf{B} + \mathbf{E},}{
#' Y = XB + E,}
#' where \eqn{\mathbf{X}}{X} is a \eqn{p}-dimensional vector,
#' \eqn{\mathbf{Y}}{Y} and \eqn{\mathbf{E}}{E} are two
#' \eqn{q}-dimensional vectors, and \eqn{\mathbf{B}}{B} is a
#' \eqn{p\times q}{p x q} matrix.
#'
#' If \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} are \emph{centered}
#' (i.e., have zero-mean columns), the function estimates \eqn{\mathbf{B}}{B}
#' by solving, for the sample
#' \eqn{(\mathbf{X}_1, \mathbf{Y}_1), \ldots, (\mathbf{X}_n, \mathbf{Y}_n)}{
#' (X_1, Y_1), \ldots, (X_n, Y_n)}, the elastic-net optimization problem
#' \deqn{
#' \min_{\mathbf{B}\in R^{q \times p}}
#' \frac{1}{2n}\sum_{i=1}^n
#' \|\mathbf{Y}_i-\mathbf{X}_i\mathbf{B}\|^2 +
#' \lambda\left[(1-\alpha)\|\mathbf{B}\|_\mathrm{F}^2 / 2 +
#' \alpha \sum_{j=1}^p \|\mathbf{B}_j\|_2\right],
#' }{
#' \min_{B\in R^{q \times p}}
#' \frac{1}{2n} \sum_{i=1}^n
#' ||Y_i - X_i B||_F^2 +
#' \lambda[(1 - \alpha) ||B||_F^2 / 2 +
#' \alpha \sum_{j=1}^p ||B_j||_2],
#' }
#' where \eqn{\|\mathbf{B}\|_\mathrm{F}}{||B||_F} stands for
#' the Frobenious norm of the matrix \eqn{\mathbf{B}}{B} and
#' \eqn{\|\mathbf{B}_j\|_2}{||B_j||_2} for the Euclidean norm
#' of the \eqn{j}-th \emph{row} of \eqn{\mathbf{B}}{B}. The choice
#' \eqn{\alpha = 0} in the elastic-net penalization corresponds to ridge
#' regression, whereas \eqn{\alpha = 1} yields a lasso-type estimator.
#' The unpenalized least-squares estimator is obtained with \eqn{\lambda = 0}.
#'
#' @param x input matrix of size \code{c(n, p)}, or a vector of length
#' \code{n}.
#' @param y response matrix of size \code{c(n, q)}, or a vector of
#' length \code{n}.
#' @param alpha elastic net mixing argument in \code{\link[glmnet]{glmnet}},
#' with \eqn{0 \le \alpha \le 1}. Alternatively, a character string indicating
#' whether the \code{"ridge"} (\eqn{\alpha = 0}) or \code{"lasso"}
#' (\eqn{\alpha = 1}) fit is to be performed.
#' @param lambda scalar giving the regularization parameter \eqn{\lambda}. If
#' \code{NULL} (default), the optimal \eqn{\lambda} is searched by
#' cross-validation. If \code{lambda} is provided, then cross-validation is
#' skipped and the fit is performed for the given \code{lambda}.
#' @param intercept flag passed to the \code{intercept} argument in
#' \code{\link[glmnet]{glmnet}} to indicate if the intercept should be fitted
#' (default; does not assume that the data is centered) or set to zero
#' (the optimization problem above is solved as-is). Defaults to \code{TRUE}.
#' @param thresh convergence threshold of the coordinate descending algorithm,
#' passed to the \code{thresh} argument in \code{\link[glmnet]{glmnet}}.
#' Defaults to \code{1e-10}.
#' @param cv_1se shall the \emph{optimal} lambda be the \code{lambda.1se}, as
#' returned by \code{\link[glmnet]{cv.glmnet}}? This favors sparser fits. If
#' \code{FALSE}, then the optimal lambda is \code{lambda.min}, the minimizer
#' of the cross-validation loss. Defaults to \code{TRUE}.
#' @param cv_nlambda the length of the sequence of \eqn{\lambda} values,
#' passed to the \code{nlambda} argument in \code{\link[glmnet]{cv.glmnet}}
#' for the cross-validation search. Defaults to \code{50}.
#' @param cv_folds number of folds to perform cross-validation. If
#' \code{NULL} (the default), then\cr \code{cv_folds <- n} internally,
#' that is, leave-one-out cross-validation is performed. Passed to the
#' \code{nfolds} argument in \code{\link[glmnet]{cv.glmnet}}.
#' @param cv_grouped passed to the \code{grouped} argument in
#' \code{\link[glmnet]{cv.glmnet}}. Defaults to \code{FALSE}.
#' @param cv_lambda passed to the \code{lambda} argument in
#' \code{\link[glmnet]{cv.glmnet}}. Defaults to\cr
#' \code{10^seq(2, -3, length.out = cv_nlambda)}.
#' @param cv_second flag to perform a second cross-validation search if the
#' optimal \eqn{\lambda} was found at the extremes of the first \eqn{\lambda}
#' sequence (indicating that the minimum may not be reliable). Defaults to
#' \code{TRUE}.
#' @param cv_tol_second tolerance for performing a second search if
#' \code{second = TRUE}. If the minimum is found at the
#' \code{100 * cv_tol_second} lower/upper percentile of search interval, then
#' the search interval is expanded for a second search. Defaults to
#' \code{0.025}.
#' @param cv_log10_exp expansion of the \eqn{\lambda} sequence if the minimum
#' is found close to its \emph{upper} extreme. If that is the case, the
#' sequence for the is set as \code{10^(log10(lambda_min) + cv_log10_exp)},
#' where \code{lambda_min} is the minimum obtained in the first sequence. If
#' the minimum is found close to the lower extreme of the sequence, then
#' \code{-rev(cv_log10_exp)} is considered. Defaults to \code{c(-0.5, 5)}.
#' @param cv_thresh convergence threshold used during cross-validation in
#' \code{\link[glmnet]{cv.glmnet}}. Defaults to \code{1e-5}.
#' @param cv_parallel passed to the \code{parallel} argument in
#' \code{\link[glmnet]{cv.glmnet}}. Defaults to \code{FALSE}.
#' @param cv_verbose flag to display information about the cross-validation
#' search with plots and messages. More useful if \code{second = TRUE}.
#' Defaults to \code{FALSE}.
#' @param ... further parameters to be passed to \code{\link[glmnet]{glmnet}}
#' to perform the final model fit.
#' @return
#' A list with the following entries:
#' \item{beta_hat}{the estimated \eqn{\mathbf{B}}{B},
#' a matrix of size \code{c(p, q)}.}
#' \item{lambda}{the optimal \eqn{\lambda} obtained by cross-validation and
#' according to \code{cv_1se}.}
#' \item{cv}{if \code{lambda = NULL}, the result of the cross-validation
#' search for the optimal \eqn{\lambda}. Otherwise, \code{NULL}.}
#' \item{fit}{the \code{glmnet} fit, computed with \code{thresh} threshold
#' and with an automatically chosen \eqn{\lambda} sequence.}
#' @details
#' If \eqn{\alpha = 1}, then the lasso-type fit shrinks to zero,
#' \emph{simultaneously}, all the elements of certain rows of
#' \eqn{\mathbf{B}}{B}, thus
#' helping the selection of the \eqn{p} most influential variables in
#' \eqn{\mathbf{X}}{X} for explaining/predicting \eqn{\mathbf{Y}}{Y}.
#'
#' The function first performs a cross-validation search for the optimal
#' \eqn{\lambda} if \code{lambda = NULL} (using \code{cv_thresh} to control
#' the convergence threshold). After the optimal penalization parameter is
#' determined, a second fit (now with convergence threshold \code{thresh})
#' using the default \eqn{\lambda} sequence in \code{\link[glmnet]{glmnet}}
#' is performed. The final estimate is obtained via
#' \code{\link[glmnet]{predict.glmnet}} from the optimal \eqn{\lambda}
#' determined in the first step.
#'
#' Due to its cross-validatory nature, \code{cv_glmnet} can be
#' computationally demanding. Approaches for reducing the computation time
#' include: considering a smaller number of folds than \code{n}, such as
#' \code{cv_folds = 10} (but will lead to random partitions of the
#' data); decrease the tolerance of the coordinate descending
#' algorithm by increasing \code{cv_thresh}; reducing the number of
#' candidate \eqn{\lambda} values with \code{nlambda}; setting
#' \code{second = FALSE} to avoid a second cross-validation; or considering
#' \code{cv_parallel = TRUE} to use a parallel backend (must be registered
#' before hand; see examples).
#'
#' By default, the \eqn{\lambda} sequence is used with \emph{standardized}
#' \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} (both divided by their
#' columnwise variances); see \code{\link[glmnet]{glmnet}} and the
#' \code{standardized} argument. Therefore, the optimal selected \eqn{\lambda}
#' value assumes standardization and must be used with care if the variables
#' are not standardized. For example, when using the ridge analytical
#' solution, a prior change of scale that depends on \eqn{q} needs to be done.
#' See the examples for the details.
#' @examples
#' ## Quick example for multivariate linear model with multivariate response
#'
#' # Simulate data
#' n <- 100
#' p <- 10; q <- 5
#' set.seed(123456)
#' x <- matrix(rnorm(n * p, sd = rep(1:p, each = n)), nrow = n, ncol = p)
#' e <- matrix(rnorm(n * q, sd = rep(q:1, each = n)), nrow = n, ncol = q)
#' beta <- matrix(((1:p - 1) / p)^2, nrow = p, ncol = q)
#' y <- x %*% beta + e
#'
#' # Fit lasso (model with intercept, the default)
#' cv_glmnet(x = x, y = y, alpha = "lasso", cv_verbose = TRUE)$beta_hat
#' \donttest{
#' ## Multivariate linear model with multivariate response
#'
#' # Simulate data
#' n <- 100
#' p <- 10; q <- 5
#' set.seed(123456)
#' x <- matrix(rnorm(n * p, sd = rep(1:p, each = n)), nrow = n, ncol = p)
#' e <- matrix(rnorm(n * q, sd = rep(q:1, each = n)), nrow = n, ncol = q)
#' beta <- matrix(((1:p - 1) / p)^2, nrow = p, ncol = q)
#' y <- x %*% beta + e
#'
#' # Fit ridge
#' cv0 <- cv_glmnet(x = x, y = y, alpha = "ridge", intercept = FALSE,
#'                  cv_verbose = TRUE)
#' cv0$beta_hat
#'
#' # Same fit for the chosen lambda
#' cv_glmnet(x = x, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE)$beta_hat
#'
#' # Fit lasso (model with intercept, the default)
#' cv1 <- cv_glmnet(x = x, y = y, alpha = "lasso", cv_verbose = TRUE)
#' cv1$beta_hat
#'
#' # Use cv_1se = FALSE
#' cv1_min <- cv_glmnet(x = x, y = y, alpha = "lasso", cv_verbose = TRUE,
#'                      cv_1se = FALSE)
#'
#' # Compare it with ridge analytical solution. Observe the factor in lambda,
#' # necessary since lambda is searched for the standardized data. Note also
#' # that, differently to the case q = 1, no standardarization with respect to
#' # y happens
#' sd_x <- apply(x, 2, function(x) sd(x)) * sqrt((n - 1) / n)
#' cv_glmnet(x = x, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           thresh = 1e-20, intercept = FALSE)$beta_hat
#' solve(crossprod(x) + diag(cv0$lambda * sd_x^2 * n)) %*% t(x) %*% y
#'
#' # If x is standardized, the match between glmnet and usual ridge
#' # analytical expression does not require scaling of lambda
#' x_std <- scale(x, scale = sd_x, center = TRUE)
#' cv_glmnet(x = x_std, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE, thresh = 1e-20)$beta_hat
#' solve(crossprod(x_std) + diag(rep(cv0$lambda * n, p))) %*% t(x_std) %*% y
#'
#' ## Simple linear model
#'
#' # Simulate data
#' n <- 100
#' p <- 1; q <- 1
#' set.seed(123456)
#' x <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' e <- matrix(rnorm(n * q), nrow = n, ncol = q)
#' beta <- 2
#' y <- x * beta + e
#'
#' # Fit by ridge (model with intercept, the default)
#' cv0 <- cv_glmnet(x = x, y = y, alpha = "ridge", cv_verbose = TRUE)
#' cv0$beta_hat
#' cv0$intercept
#'
#' # Comparison with linear model with intercept
#' lm(y ~ 1 + x)$coefficients
#'
#' # Fit by ridge (model without intercept)
#' cv0 <- cv_glmnet(x = x, y = y, alpha = "ridge", cv_verbose = TRUE,
#'                  intercept = FALSE)
#' cv0$beta_hat
#'
#' # Comparison with linear model without intercept
#' lm(y ~ 0 + x)$coefficients
#'
#' # Same fit for the chosen lambda (and without intercept)
#' cv_glmnet(x = x, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE)$beta_hat
#'
#' # Same for lasso (model with intercept, the default)
#' cv1 <- cv_glmnet(x = x, y = y, alpha = "lasso")
#' cv1$beta_hat
#'
#' ## Multivariate linear model (p = 3, q = 1)
#'
#' # Simulate data
#' n <- 50
#' p <- 10; q <- 1
#' set.seed(123456)
#' x <- matrix(rnorm(n * p, mean = 1, sd = rep(1:p, each = n)),
#'             nrow = n, ncol = p)
#' e <- matrix(rnorm(n * q), nrow = n, ncol = q)
#' beta <- ((1:p - 1) / p)^2
#' y <- x %*% beta + e
#'
#' # Fit ridge (model without intercept)
#' cv0 <- cv_glmnet(x = x, y = y, alpha = "ridge", intercept = FALSE,
#'                  cv_verbose = TRUE)
#' cv0$beta_hat
#'
#' # Same fit for the chosen lambda
#' cv_glmnet(x = x, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE)$beta_hat
#'
#' # Compare it with ridge analytical solution. Observe the factor in lambda,
#' # necessary since lambda is searched for the standardized data
#' sd_x <- apply(x, 2, function(x) sd(x)) * sqrt((n - 1) / n)
#' sd_y <- sd(y) * sqrt((n - 1) / n)
#' cv_glmnet(x = x, y = y, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE, thresh = 1e-20)$beta_hat
#' solve(crossprod(x) + diag(cv0$lambda * sd_x^2 / sd_y * n)) %*% t(x) %*% y
#'
#' # If x and y are standardized, the match between glmnet and usual ridge
#' # analytical expression does not require scaling of lambda
#' x_std <- scale(x, scale = sd_x, center = TRUE)
#' y_std <- scale(y, scale = sd_y, center = TRUE)
#' cv_glmnet(x = x_std, y = y_std, alpha = "ridge", lambda = cv0$lambda,
#'           intercept = FALSE, thresh = 1e-20)$beta_hat
#' solve(crossprod(x_std) + diag(rep(cv0$lambda * n, p))) %*% t(x_std) %*% y_std
#'
#' # Fit lasso (model with intercept, the default)
#' cv1 <- cv_glmnet(x = x, y = y, alpha = "lasso", cv_verbose = TRUE)
#' cv1$beta_hat
#'
#' # ## Parallelization
#' #
#' # # Parallel
#' # doMC::registerDoMC(cores = 2)
#' # microbenchmark::microbenchmark(
#' # cv_glmnet(x = x, y = y, nlambda = 100, cv_parallel = TRUE),
#' # cv_glmnet(x = x, y = y, nlambda = 100, cv_parallel = FALSE),
#' # times = 10)
#' }
#' @author Eduardo García-Portugués. Initial contributions by Gonzalo
#' Álvarez-Pérez.
#' @references
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. \emph{Journal of
#' Statistical Software}, 33(1):1--22. \doi{10.18637/jss.v033.i01}
#' @export
cv_glmnet <- function(x, y, alpha = c("lasso", "ridge")[1], lambda = NULL,
                      intercept = TRUE, thresh = 1e-10, cv_1se = TRUE,
                      cv_nlambda = 50, cv_folds = NULL, cv_grouped = FALSE,
                      cv_lambda = 10^seq(2, -3, length.out = cv_nlambda),
                      cv_second = TRUE, cv_tol_second = 0.025,
                      cv_log10_exp = c(-0.5, 3), cv_thresh = 1e-5,
                      cv_parallel = FALSE, cv_verbose = FALSE, ...) {

  # Check for x and y as matrices
  if (is.null(dim(x))) {

    x <- cbind(x)

  }
  if (is.null(dim(y))) {

    y <- cbind(y)

  }

  # Sample size
  n <- nrow(x)

  # Default cv_folds adopted in the cross-validation procedure: n for
  # leave-one-out cross-validation
  cv_folds <- ifelse(is.null(cv_folds), n, cv_folds)

  # Ridge or lasso regularization?
  alpha <- ifelse(is.character(alpha),
                  switch(alpha, "ridge" = 0, "lasso" = 1), alpha)

  # Distinguish the type of family in glmnet::cv.glmnet to avoid errors
  # in edge cases
  family <- ifelse(ncol(y) == 1, "gaussian", "mgaussian")

  # What if ncol(x) == 1?
  if (ncol(x) == 1) {

    x <- cbind(x, 0)
    p <- 1

  } else {

    p <- ncol(x)

  }

  # Search for the optimal lambda
  if (is.null(lambda)) {

    # Caution! glmnet::cv.glmnet affects the seed, as a call to sample()
    # is performed for determining the folds of the data if foldid is not
    # provided (even if the deterministic leave-one-out partition is used).
    # To avoid undesirable replication side effects, the seed is reset when
    # exiting the function
    if (exists(".Random.seed")) {

      old <- .Random.seed
      on.exit({.Random.seed <<- old})

    }

    # First fit with a given lambda sequence (cv_lambda) to avoid letting
    # glmnet::cv.glmnet choosing the lambda sequence automatically (as of
    # version 2.0-16, this can generate errors). The ... are related only to
    # arguments for glmnet::glmnet
    weights <- rep(1, n)
    cv <- glmnet::cv.glmnet(x = x, y = y, weights = weights,
                            offset = NULL, family = family, alpha = alpha,
                            lambda = cv_lambda, nfolds = cv_folds,
                            type.measure = "mse", grouped = cv_grouped,
                            keep = FALSE, parallel = cv_parallel,
                            intercept = intercept, thresh = cv_thresh, ...)

    # Compute where the minimum is reached, as the *upper* percentile of the
    # sequence. Beware: cv$cvm is sorted in a decreasing way!
    ind_lambda_min <- which(cv$lambda ==
                              ifelse(cv_1se, cv$lambda.1se, cv$lambda.min))
    prc_lambda_min <- ind_lambda_min / (cv_nlambda + 1)

    # Save for message()
    cv_min <- cv$cvm[ind_lambda_min]
    lambda_min <- cv$lambda[ind_lambda_min]

    # Avoid problems with glmnet::plot.cv.glmnet (cv$nzero stores the number
    # of non-zero coefficients at each lambda)
    if (length(cv$nzero) == 1) {

      cv$nzero <- rep(NA, cv_nlambda)

    }

    # Reset par() for the user
    if (cv_verbose) {

      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))

    }

    # Check if it is required a second call because the minimum is found at
    # one extreme of the sequence (quantified as either the 100 * cv_tol_second
    # largest or smallest values)
    if (cv_second & (prc_lambda_min <= cv_tol_second |
                     prc_lambda_min >= (1 - cv_tol_second))) {

      # Save for message()
      cv_min_old <- cv_min
      lambda_min_old <- lambda_min
      log_lambda_range_old <- log10(range(cv$lambda))

      # Expand sequence on the lower part ($lambda is sorted decreasingly!),
      # otherwise the sequence is expanded on the upper part
      if (prc_lambda_min >= (1 - cv_tol_second)) {

        cv_log10_exp <- -rev(cv_log10_exp)

      }
      log_lambda_range <- log10(lambda_min_old) + cv_log10_exp

      # Decreasing sequence of lambdas
      cv_lambda <- 10^seq(log_lambda_range[2], log_lambda_range[1],
                          length.out = cv_nlambda / 2)

      # Second fit with expanded lambda sequence and same options
      cv_old <- cv
      cv <- glmnet::cv.glmnet(x = x, y = y, weights = weights,
                              offset = NULL, family = family, alpha = alpha,
                              lambda = cv_lambda, nfolds = cv_folds,
                              type.measure = "mse", grouped = cv_grouped,
                              keep = FALSE, parallel = cv_parallel,
                              intercept = intercept, thresh = cv_thresh, ...)

      # Save for message()
      ind_lambda_min <- which(cv$lambda ==
                                ifelse(cv_1se, cv$lambda.1se, cv$lambda.min))
      cv_min <- cv$cvm[ind_lambda_min]
      lambda_min <- cv$lambda[ind_lambda_min]

      # Avoid problems with glmnet::plot.cv.glmnet
      if (length(cv$nzero) == 1) {

        cv$nzero <- rep(NA, cv_nlambda / 2)

      }

      # Display information
      if (cv_verbose) {

        message(paste0(
          "Optimum (", sprintf("%.3f", cv_min_old), ") at ",
          sprintf("%.3f", log10(lambda_min_old)),
          " is close to the boundary of (",
          sprintf("%.3f", log_lambda_range_old[1]), ", ",
          sprintf("%.3f", log_lambda_range_old[2]),
          "), expanded the search to (",
          sprintf("%.3f", log_lambda_range[1]),
          ", ", sprintf("%.3f", log_lambda_range[2]),
          ") with new optimum (",
          sprintf("%.3f", cv_min), ") at ",
          sprintf("%.3f", log10(lambda_min))
        ))

        # Plot the first and second cross-validation search, in log10 scale
        par(mfrow = c(1, 2))
        plot(cv_old, axes = FALSE, xlab = "log10(lambda)")
        labs <- pretty(log10(cv_old$lambda))
        axis(1, at = log(10^labs), labels = labs); axis(2); box()
        plot(cv, axes = FALSE, xlab = "log10(lambda)")
        labs <- pretty(log10(cv$lambda))
        axis(1, at = log(10^labs), labels = labs); axis(2); box()

      }

    } else {

      if (cv_verbose) {

        message(paste0(
          "Optimum (", sprintf("%.3f", cv_min), ") at ",
          sprintf("%.3f", log10(lambda_min)), " found in (",
          sprintf("%.3f", log10(min(cv$lambda))), ", ",
          sprintf("%.3f", log10(max(cv$lambda))), ")"
        ))

        par(mfrow = c(1, 1))
        plot(cv, axes = FALSE, xlab = "log10(lambda)")
        labs <- pretty(log10(cv$lambda))
        axis(1, at = log(10^labs), labels = labs); axis(2); box()

      }

    }

    # Optimal regularization penalty parameter
    lambda <- lambda_min

    # Known lambda: optimal choice of lambda is not performed
  } else {

    # Check lambda
    if (length(lambda) > 1) {

      stop("lambda (if provided) must be a single number")

    }

    # cv object for coherency
    cv <- NULL

  }

  # Fit glmnet::glmnet with thresh accuracy and with default lambda sequence
  fit <- glmnet::glmnet(x = x, y = y, family = family, alpha = alpha,
                        intercept = intercept, thresh = thresh, ...)

  # Prediction for the supplied lambda (exact = TRUE is required since
  # the desired lambda may not be in the lambda sequence). This implies a
  # refitting, hence why all the required arguments are passed
  beta_hat <- predict(fit, type = "coefficients", s = lambda, exact = TRUE,
                      x = x, y = y, alpha = alpha, family = family,
                      thresh = thresh, intercept = intercept, ...)
  beta_hat <- switch(family, "gaussian" = beta_hat,
                     "mgaussian" = do.call(cbind, args = as.list(beta_hat)))

  # The intercept is returned in the first row even if intercept = FALSE
  # (in this case it is set to zero)
  a0 <- beta_hat[1, , drop = FALSE]

  # We exclude the first row of beta_hat
  beta_hat <- beta_hat[2:(p + 1), , drop = FALSE]

  # Result
  return(list("beta_hat" = beta_hat, "intercept" = a0,
              "lambda" = lambda, "cv" = cv, "fit" = fit))

}
