

#' @title Ontario temperature and electricity consumption during 2010--2014
#'
#' @description Real dataset employed Benatia et al. (2017). Contains the
#' hourly electricity consumption and air temperature curves in the province
#' of Ontario (Canada). It features a set of daily curves during the summer
#' months of 2010--2014.
#'
#' @docType data
#' @format A list with the following entries:
#' \describe{
#'   \item{temp}{an \code{\link[fda.usc]{fdata}} with 368 smoothed
#'   daily temperature (in Celsius degrees) curves of the Ontario province,
#'   discretized on 73 equispaced grid points on \eqn{[-24, 48]}
#'   (see examples).}
#'   \item{elec}{an \code{\link[fda.usc]{fdata}} with the daily
#'   electricity consumption (in gigawatts) curves of the Ontario province.
#'   Discretized on 25 equispaced grid points on \eqn{[0, 24]}.}
#'   \item{df}{a dataframe with time metadata for each curve:
#'   \itemize{
#'     \item{\code{date}: the date of the observation, a \code{\link{POSIXct}}
#'     object.}
#'     \item{\code{weekday}: the weekday of the observation.}
#'   }}
#' }
#' @details
#' The summer months correspond to June 1st to September 15th. Weekend days and
#' holidays are disregarded.
#'
#' The smoothed temperature curves are constructed by a weighted average of
#' the temperatures of 41 Ontarian cities that is afterwards smoothed with
#' a local polynomial regression. The curves correspond to a 3-days window
#' of the temperature (see examples). The temperature is standardized such
#' that its original minimum, 6 ºC, is subtracted.
#'
#' The electricity consumption curves are discretized on the interval
#' \eqn{[0, 24]}. That means that the last observation of the
#' \eqn{i}-th curve is the same as the first observation of the
#' \eqn{(i + 1)}-th curve \emph{if} the curves correspond to consecutive days.
#'
#' See more details about the construction of the dataset in Benatia et al.
#' (2017).
#' @source
#' The dataset comes from the companion data to Benatia et al. (2017), which
#' was retrieved from the \href{https://www.davidbenatia.com/publication/}{
#' first author's website}. The source of the electricity consumption data is
#' the \href{https://www.ieso.ca/}{System operator's website}. The source
#' of the preprocessed temperature values is the
#' \href{https://climat.meteo.gc.ca/}{Environment Canada's website}.
#' @author Data gathered and processed by David Benatia, Marine Carrasco, and
#' Jean-Pierre Florens. Javier Álvarez-Liébana and Eduardo García-Portugués
#' imported the dataset and added temporal metadata.
#' @references
#' Benatia, D., Carrasco, M. and Florens, J. P. (2017) Functional linear
#' regression with functional response. \emph{Journal of Econometrics},
#' 201(2):269--291. \doi{10.1016/j.jeconom.2017.08.008}
#' @examples
#' ## Show data
#'
#' # Load data
#' data("ontario")
#'
#' # Plot
#' old_par <- par(mfrow = c(1, 2))
#' plot(ontario$temp)
#' plot(ontario$elec)
#' par(old_par)
#'
#' # Observe the 3-day windows for each observation
#' plot(ontario$temp$argvals, ontario$temp$data[2, ], type = "o",
#'      xlim = c(-48, 72), ylim = c(7, 13), xlab = "Hours",
#'      ylab = "Electricity consumption", pch = 16)
#' points(ontario$temp$argvals - 24, ontario$temp$data[1, ], col = 3, pch = 2)
#' points(ontario$temp$argvals + 24, ontario$temp$data[3, ], col = 2, cex = 1.5)
#' abline(v = 24 * -2:3, lty = 2)
#' legend("top", legend = c("Curve 1", "Curve 2", "Curve 3"), col = c(3, 1, 2),
#'        pt.cex = c(1, 1, 1.5), pch = c(2, 16, 1))
#'
#' # If the days are not consecutive, then the electricity consumptions at the
#' # end of one day and the beginning of the next do not match
#' head(abs(ontario$elec$data[-368, 25] - ontario$elec$data[-1, 1]))
#' head(diff(ontario$df$date))
#' \donttest{
#' ## Test the linear model with functional response and predictor
#'
#' (comp_flmfr <- flm_test(X = ontario$temp, Y = ontario$elec,
#'                         est_method = "fpcr_l1s"))
#' (simp_flmfr <- flm_test(X = ontario$temp, Y = ontario$elec,
#'                         beta0 = 0, est_method = "fpcr_l1s"))
#'
#' # Visualize estimation
#' filled.contour(x = ontario$temp$argvals, y = ontario$elec$argvals,
#'                z = comp_flmfr$fit_flm$Beta_hat,
#'                color.palette = viridisLite::viridis, nlevels = 20)
#' }
"ontario"


#' @title AEMET daily temperatures during 1974--2013
#'
#' @description Series of daily temperatures of 73 Spanish weather stations
#' during the 40-year period 1974--2013.
#'
#' @docType data
#' @format A list with the following entries:
#' \describe{
#'   \item{temp}{an \code{\link[fda.usc]{fdata}} with 2892 temperature
#'   (in Celsius degrees) curves, discretized on 365 equispaced
#'   grid points (days) on \eqn{[0.5, 364.5]}. Each curve corresponds to the
#'   yearly records of a weather station.}
#'   \item{df}{a dataframe with metadata for each curve:
#'   \itemize{
#'     \item{\code{ind}: identifier of the weather station.}
#'     \item{\code{name}: name of the weather station.}
#'     \item{\code{year}: year of the observation.}
#'   }}
#' }
#' @details
#' For consistency with the \code{\link[fda.usc]{fda.usc-package}}'s
#' \code{\link[fda.usc]{aemet}} dataset, the names and identifiers of the 73
#' weather stations are the same as in that dataset. Only a minor fix has been
#' applied to the "A CORUÑA/ALVEDRO" station, whose identifier was the same
#' as the "A CORUÑA" station, \code{"1387"}. The former was set to
#' \code{"1387A"}.
#'
#' Information about the province, altitude, longitude, and latitude of
#' each weather station can be retrieved in \code{df} from the
#' \code{\link[fda.usc]{fda.usc-package}}'s \code{\link[fda.usc]{aemet}}
#' dataset.
#'
#' The dataset is a curated version of a larger database of 115 stations. It
#' excludes stations with inconsistent records or that were relocated, closed,
#' or opened during the 40-year period. There are 9 stations with missing
#' years. The total of missing years is 28.
#'
#' In leap years, the daily-average temperature is computed as the average of
#' February 28th and 29th.
#' @source The data was retrieved from the FTP of the
#' \href{http://www.aemet.es/es/portada/}{Meteorological State Agency of Spain
#' (AEMET)} in 2014 using a processing script by the authors of the
#' \code{\link[fda.usc]{fda.usc-package}}.
#' @author
#' Original data processing scripts by Manuel Febrero-Bande and Manuel Oviedo
#' de la Fuente. Adaptations by Eduardo García-Portugués.
#' @references
#' Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical Computing
#' in Functional Data Analysis: The R Package fda.usc. \emph{Journal of
#' Statistical Software}, 51(4):1--28. \doi{10.18637/jss.v051.i04}
#' @examples
#' ## Data splitting
#'
#' # Load raw data
#' data("aemet_temp")
#'
#' # Partition the dataset in the first and last 20 years
#' with(aemet_temp, {
#'   ind_pred <- which((1974 <= df$year) & (df$year <= 1993))
#'   ind_resp <- which((1994 <= df$year) & (df$year <= 2013))
#'   aemet_temp_pred <<- list("df" = df[ind_pred, ], "temp" = temp[ind_pred])
#'   aemet_temp_resp <<- list("df" = df[ind_resp, ], "temp" = temp[ind_resp])
#' })
#'
#' # Average the temperature on each period
#' mean_aemet <- function(x) {
#'   m <- tapply(X = 1:nrow(x$temp$data), INDEX = x$df$ind,
#'               FUN = function(i) colMeans(x$temp$data[i, , drop = FALSE],
#'                                          na.rm = TRUE))
#'  x$temp$data <- do.call(rbind, m)
#'  return(x$temp)
#' }
#'
#' # Build predictor and response fdatas
#' aemet_temp_pred <- mean_aemet(aemet_temp_pred)
#' aemet_temp_resp <- mean_aemet(aemet_temp_resp)
#'
#' # Plot
#' old_par <- par(mfrow = c(1, 2))
#' plot(aemet_temp_pred)
#' plot(aemet_temp_resp)
#' par(old_par)
#'
#' # Average daily temperatures
#' day_avg_pred <- func_mean(aemet_temp_pred)
#' day_avg_resp <- func_mean(aemet_temp_resp)
#'
#' # Average yearly temperatures
#' avg_year_pred <- rowMeans(aemet_temp_pred$data)
#' avg_year_resp <- rowMeans(aemet_temp_resp$data)
#' \donttest{
#' ## Test the linear model with functional response and predictor
#'
#' (comp_flmfr <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
#'                         est_method = "fpcr_l1s"))
#' beta0 <- diag(rep(1, length(aemet_temp_pred$argvals)))
#' (simp_flmfr <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
#'                         beta0 = beta0, est_method = "fpcr_l1s"))
#'
#' # Visualize estimation
#' filled.contour(x = aemet_temp_pred$argvals, y = aemet_temp_resp$argvals,
#'                z = comp_flmfr$fit_flm$Beta_hat,
#'                color.palette = viridisLite::viridis, nlevels = 20)
#'
#' ## Test the linear model with scalar response and functional predictor
#'
#' (comp_flmsr <- flm_test(X = aemet_temp_pred, Y = avg_year_resp,
#'                         est_method = "fpcr_l1s"))
#' (simp_flmsr <- flm_test(X = aemet_temp_pred, Y = avg_year_resp,
#'                         beta0 = 1 / 365, est_method = "fpcr_l1s"))
#'
#' # Visualize estimation
#' plot(aemet_temp_pred$argvals, comp_flmsr$fit_flm$Beta_hat, type = "l",
#'      ylim = c(0, 30 / 365))
#' abline(h = 1 / 365, col = 2)
#'
#' ## Test the linear model with functional response and scalar predictor
#'
#' (comp_frsp <- flm_test(X = avg_year_pred, Y = aemet_temp_resp))
#' (simp_frsp <- flm_test(X = avg_year_pred, Y = aemet_temp_resp, beta0 = 1))
#'
#' ## Test the linear model with scalar response and predictor
#'
#' (comp_srsp <- flm_test(X = avg_year_pred, Y = avg_year_resp))
#' (simp_srsp <- flm_test(X = avg_year_pred, Y = avg_year_resp, beta0 = 1))
#' }
"aemet_temp"
