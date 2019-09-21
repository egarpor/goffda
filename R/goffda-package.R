

#' @title \pkg{goffda} -- Goodness-of-Fit Tests for Functional Data
#'
#' @description Implementation of several goodness-of-fit tests for functional
#' data. Currently, mostly related with the functional linear model with
#' functional/scalar response and functional/scalar predictor. The package
#' allows for the replication of the data applications considered in
#' García-Portugués, Álvarez-Liébana, Álvarez-Pérez and González-Manteiga
#' (2019) <arXiv:1909.07686>.
#'
#' @author Eduardo García-Portugués and Javier Álvarez-Liébana.
#' @references
#' García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
#' González-Manteiga, W. (2019). A goodness-of-fit test for the functional
#' linear model with functional response. \emph{arXiv:1909.07686}.
#' \url{https://arxiv.org/abs/1909.07686}
#'
#' García-Portugués, E., González-Manteiga, W. and Febrero-Bande, M. (2014). A
#' goodness-of-fit test for the functional linear model with scalar response.
#' \emph{Journal of Computational and Graphical Statistics}, 23(3):761--778.
#' \url{http://doi.org/10.1080/10618600.2013.812519}
#' @docType package
#' @name goffda-package
#' @import graphics stats Rcpp
#' @importFrom fda.usc plot.fdata
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom grDevices gray
#' @useDynLib goffda
#' @aliases goffda goffda-package
NULL
