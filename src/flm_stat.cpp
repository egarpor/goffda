#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

//' @title Projected Cramér--von Mises test statistic for the
//' goodness-of-fit test of functional linear models
//'
//' @description Computation of the Projected Cramér--von Mises (PCvM) test
//' statistic and its associated \eqn{\mathbf{A}_\bullet}{A_\bullet} matrix.
//' For a sample of functional covariates \eqn{X_1, \ldots, X_n}, the test
//' statistic is computed from
//' \eqn{\mathbf{x}_{1,p}, \ldots, \mathbf{x}_{n,p}}{
//' x_{1, p}, \ldots, x_{n, p}},
//' the coefficients (scores) of the sample in a \eqn{p}-truncated basis
//' expansion, such as Functional Principal Components (FPC).
//'
//' The PCvM statistic is defined as
//' \deqn{\mathrm{PCvM}_{n,p,q} =
//' c \cdot \mathrm{tr}(\hat\mathbf{E}_q' \mathbf{A}_\bullet \hat\mathbf{E}_q)}{
//' PCvM_{n, p, q} = c * tr(\hat E_q' A_\bullet \hat E_q)
//' }
//' where
//' \deqn{c = 2 \pi^{(p + q) / 2 - 1} / (q \Gamma(p / 2) \Gamma(q / 2) n^2),}
//' \eqn{\hat\mathbf{E}_q}{\hat E_q} is the \eqn{n \times q}{n x q}
//' matrix of multivariate residuals, and
//' \eqn{\mathbf{A}_\bullet}{A_\bullet} is a \eqn{n \times n}{n x n}
//' matrix whose \eqn{ij}-th element is
//' \eqn{\sum_{r = 1}^n A_{ijr}}, for \eqn{A_{ijr}} depending on
//' \eqn{(\mathbf{x}_{i,p}, \mathbf{x}_{j,p}, \mathbf{x}_{r,p})}{
//' (x_{i, p}, x_{j, p}, x_{r, p})}. Its exact expression can be seen in
//' Escanciano (2006) and García-Portugués et al. (2021).
//'
//' @param E the matrix of multivariate residuals, with dimension
//' \code{c(n, q)}. A vector if \eqn{q = 1}.
//' @param p dimension of the covariates space. Must be a positive integer.
//' @param Adot_vec output from \code{\link{Adot}}. A vector of length
//' \code{n * (n - 1) / 2 + 1}. This corresponds to the most expensive
//' computation in the test statistic.
//' @param constant whether to include the constant of the PCvM test
//' statistic, \eqn{c}, in its computation. Defaults to \code{TRUE}.
//' @param X a matrix of size \code{c(n, p)} containing the coefficients
//' (scores) of the functional data in a \code{p}-truncated \emph{orthonormal}
//' basis expansion, such as FPC. Must not contain repeated rows.
//' @return
//' \itemize{
//'   \item{\code{flm_stat}: the value of the test statistic, a scalar.}
//'   \item{\code{A_dot}: a vector of length \code{n * (n - 1) / 2 + 1}.
//'   The first entry contains the common diagonal element of
//'   \eqn{\mathbf{A}_\bullet}{A_\bullet}. The remaining entries are the upper
//'   triangular matrix (excluding the diagonal) of
//'   \eqn{\mathbf{A}_\bullet}{A_\bullet}, stacked by columns.}
//' }
//' @author Eduardo García-Portugués.
//' @details
//' \code{Adot} assumes that \code{X} does not have repeated rows or otherwise
//' \code{NaN}s will be present in the result. If \code{X} has repeated rows,
//' \code{Adot} will throw a warning.
//'
//' The implementation of the PCvM test statistic for scalar response is
//' addressed in García-Portugués et al. (2014), whereas García-Portugués et al.
//' (2021) presents its multivariate extension and shows that
//' \eqn{\mathbf{A}_\bullet}{A_\bullet} induces a weighted quadratic norm (if
//' there are no repetitions in the sample). The PCvM statistic is rooted in
//' the proposal by Escanciano (2006).
//'
//' Both \code{flm_stat} and \code{A_dot} are coded in C++.
//' @references
//' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
//' Gonzalez-Manteiga, W. (2021). A goodness-of-fit test for the functional
//' linear model with functional response. \emph{Scandinavian Journal of
//' Statistics}, 48(2):502--528. \doi{10.1111/sjos.12486}
//'
//' Escanciano, J. C. (2006) A consistent diagnostic test for regression
//' models using projections. \emph{Econometric Theory}, 22(6):1030–-1051.
//' \doi{10.1017/S0266466606060506}
//'
//' García-Portugués, E., González-Manteiga, W. and Febrero-Bande, M. (2014). A
//' goodness-of-fit test for the functional linear model with scalar response.
//' \emph{Journal of Computational and Graphical Statistics}, 23(3):761--778.
//' \doi{10.1080/10618600.2013.812519}
//' @examples
//' ## flm_stat
//'
//' # Generate data
//' n <- 200
//' q <- 2
//' p <- 3
//' E <- matrix(rnorm(n * q), nrow = n, ncol = q)
//' X_fdata <- r_ou(n = n, t = seq(0, 1, l = 101))
//'
//' # Compute FPC
//' X_fpc <- fpc(X_fdata)
//'
//' # Adot
//' Adot_vec <- Adot(X = X_fpc[["scores"]])
//'
//' # Check equality
//' constant <- n^(-2) * 2 * pi^((p / 2) - 1) / gamma(p / 2)
//' constant * .Fortran("pcvm_statistic", n = as.integer(n),
//'                     Adot_vec = Adot_vec, residuals = E[, 2],
//'                     statistic = 0)$statistic
//' flm_stat(E = E[, 2, drop = FALSE], p = p, Adot_vec = Adot_vec,
//'          constant = FALSE)
//'
//' ## Adot
//'
//' # Generate data
//' n <- 200
//' X_fdata <- r_ou(n = n, t = seq(0, 1, l = 101))
//'
//' # Compute FPC
//' X_fpc <- fpc(X_fdata)
//'
//' # Using inprod_fdata and Adot
//' Adot_vec <- Adot(X = X_fpc[["scores"]])
//'
//' # Check with fda.usc::Adot with adequate inprod
//' head(drop(Adot_vec))
//' head(fda.usc::Adot(X_fdata))
//'
//' # Obtention of the entire Adot matrix
//' Ad <- diag(rep(Adot_vec[1], n))
//' Ad[upper.tri(Ad, diag = FALSE)] <- Adot_vec[-1]
//' head(Ad <- t(Ad) + Ad - diag(diag(Ad)))
//'
//' # Positive definite
//' eigen(Ad)$values
//'
//' # # Warning if X contains repeated observations
//' # Adot(X = rbind(1:3, 1:3, 3:5))
//' \donttest{
//' # Comparison with .Fortran("adot", PACKAGE = "fda.usc")
//' n <- as.integer(n)
//' a <- as.double(rep(0, (n * (n - 1) / 2 + 1)))
//' inprod <- X_fpc[["scores"]] %*% t(X_fpc[["scores"]])
//' inprod <- inprod[upper.tri(inprod, diag = TRUE)]
//' X <- X_fpc[["scores"]]
//' microbenchmark::microbenchmark(
//'   .Fortran("adot", n = n, inprod = inprod, Adot_vec = a,
//'            PACKAGE = "fda.usc"),
//'   Adot(X = X),
//'   times = 50, control = list(warmup = 10))
//' }
//' @export
// [[Rcpp::export]]
double flm_stat(arma::mat E, int p, arma::vec Adot_vec, bool constant = true) {

  // n and q
  arma::uword n = E.n_rows;
  arma::uword q = E.n_cols;

  // Check n from the length of Adot_vec ((n^2 - n) / 2 + 1)
  if (n != 0.5 * (1 + std::sqrt(1 + 8.0 * (Adot_vec.n_elem - 1)))) {

    stop("Dimensions of Adot and E do not coincide");

  }

  // Constant
  double c = 1;
  if (constant) {

    c = (2 * std::pow(M_PI, 0.5 * (p + q) - 1)) /
      (q * R::gammafn(0.5 * p) * R::gammafn(0.5 * q) * n * n);

  }

  // Transpose E to access columns instead of rows
  E = E.t();

  // Use the double sum expression, as in this expression the symmetry of
  // Adot can be exploited
  double sym_sum = 0;
  for (arma::uword i = 1; i < n; i++) {
    for (arma::uword j = 1; j <= i; j++) {

      // Sum for the symmetric part
      arma::uword ind_ij = i * (i - 1) / 2 + j;
      sym_sum += arma::dot(E.unsafe_col(i), E.unsafe_col(j - 1)) *
        Adot_vec[ind_ij];

    }

    // Check for user interruptions 1 out of 10 times
    double test = n * 0.1;
    test -= floor(test);
    if (test == 0) {

      checkUserInterrupt();

    }

  }

  // Statistic computed as the sum of the diagonal and the symmetric part
  double stat = Adot_vec[0] * arma::trace(E * E.t()) + 2 * sym_sum;

  // Result
  return c * stat;

}


//' @rdname flm_stat
//' @export
// [[Rcpp::export]]
arma::vec Adot(arma::mat X) {

  // Get n
  arma::uword n = X.n_rows;

  // Compute the n x n matrix of inner products of X
  arma::mat inprod_full = arma::trimatu(X * X.t(), 0);

  // Index of the upper tridiagonal matrix, by columns, and including
  // the diagonal
  arma::uvec ind_tri = arma::find(arma::trimatu(arma::ones(n, n), 0));

  // Vector half of inprod_full (containing the diagonal as well)
  arma::vec inprod = inprod_full.elem(ind_tri);

  // Adot
  arma::vec Adot_vec = arma::zeros(n * (n - 1) / 2 + 1);

  // The first element of Adot_vec is the common diagonal element
  Adot_vec[0] = M_PI * (n + 1);

  // Loop on the rest of the elements are the lower triangle matrix of Adot
  for (arma::uword i = 2; i <= n; i++) {
    for (arma::uword j = 1; j < i; j++) {

      // Sum on the r index
      double sum_r = 0;
      for (arma::uword r = 1; r <= n; r++) {

        // From the definition of Aijr0
        if ((i == r) | (j == r)) {

          // Sum variable
          sum_r += M_PI;

        } else {

          // Auxiliar variables for the indexes
          arma::uword aux_i = i * (i - 1) / 2;
          arma::uword aux_j = j * (j - 1) / 2;
          arma::uword aux_r = r * (r - 1) / 2;

          // Indexes
          arma::uword ij = aux_i + j;
          arma::uword ii = aux_i + i;
          arma::uword jj = aux_j + j;
          arma::uword rr = aux_r + r;

          // Lower triangular part of the matrix
          arma::uword ir = 0;
          if (i > r) {

            ir = aux_i + r;

          // Upper triangular part
          } else {

            ir = aux_r + i;

          }

          // Lower triangular part of the matrix
          arma::uword rj = 0;
          if (r > j) {

            rj = aux_r + j;

          // Upper triangular part
          } else {

            rj = aux_j + r;

          }

          // Symmetry
          arma::uword jr = rj;

          // Computation of the quotient
          double quo = inprod[ij - 1] - inprod[ir - 1] - inprod[rj - 1] +
            inprod[rr - 1];
          quo /=
            std::sqrt((inprod[ii - 1] - 2 * inprod[ir - 1] + inprod[rr - 1]) *
              (inprod[jj - 1] - 2 * inprod[jr - 1] + inprod[rr - 1]));

          // Avoid numerical problems on acos
          quo = std::min(std::max(quo, -1.0), 1.0);

          // Sum
          sum_r += std::acos(-quo);

        }

      }

      // Enter the ij-th element of A
      arma::uword ind_ij = (i - 1) * (i - 2) / 2 + j;
      Adot_vec[ind_ij] = sum_r;

    }

    // Check for user interruptions 1 out 10 times
    double test = n * 0.1;
    test -= std::floor(test);
    if (test == 0) {

      checkUserInterrupt();

    }

  }

  // Check if there are NaNs
  if (!Adot_vec.is_finite()) {

    warning("Non-finite values in Adot: check X for repeated rows");

  }

  // Result
  return Adot_vec;

}
