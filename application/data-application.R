
# Clean workspace
rm(list = ls())

# Load packages
library(fda.usc)
library(goffda)
library(lubridate)
library(viridisLite)
(version_goffda <- packageVersion("goffda")) # Originally run with 0.0.5

### Required functions and constants

# Bootstrap replicates
B <- 1e4

# Save figures?
save_fig <- TRUE

# Improve visualizaation of saved figures
if (save_fig) {

  pdf <- function(..., cex.spe = 1.85) {

    grDevices::pdf(...)
    par(cex.axis = cex.spe, cex.main = cex.spe, cex.lab = cex.spe,
        mar = c(5.5, 4.5, 3.5, 1.5) + 0.1)

  }

}

# Test for no functional effects by Kokoszka et al. (2008)
kokoszka_test <- function(X_fpc, Y_fpc, thre_p = 0.99, thre_q = 0.99) {

  # Sample size
  n <- nrow(X_fpc[["scores"]])
  stopifnot(n == nrow(Y_fpc[["scores"]]))

  # Cuttofs for p and q
  p <- min(which((cumsum(X_fpc$d^2) / sum(X_fpc$d^2)) > thre_p))
  q <- min(which((cumsum(Y_fpc$d^2) / sum(Y_fpc$d^2)) > thre_q))

  # Statistic
  stat <- 0
  for (j in 1:p) {
    for (k in 1:q) {

      stat <- stat + 1 / (X_fpc[["d"]][j] * Y_fpc[["d"]][k])^2 *
        mean(X_fpc[["scores"]][, j] * Y_fpc[["scores"]][, k])^2

    }
  }
  stat <- n^3 * stat

  # p-value
  p_value <- pchisq(stat, p * q, lower.tail = FALSE)

  # Return htest object
  meth <- "Kokoszka et al. (2008) test for significance"
  result <- structure(list(statistic = c("statistic" = stat),
                           p.value = p_value, method = meth,
                           parameter = c("p" = p , "q" = q),
                           p = p, q = q, data.name = "Y ~ X"))
  class(result) <- "htest"
  return(result)

}

# Add custom month axis
months_axis <- function(side = 1, cex.spe = 1.85, text = "Day", line = 4, ...) {

  axis(side = side, at = pmin(0.5 + c(0, cumsum(days_in_month(1:12))), 364.5),
       labels = c(month.abb, "Jan"), las = 2, ...)
  mtext(text = text, side = side, line = line, cex = cex.spe, las = 0)

}

# Convenience function to get the "nonoverlapping curves" (understood as those
# curves measured on disjoint intervals) on the ontario dataset. day_start
# measures when the first nonoverlapping curve is considered (first day of
# available data, second, or third), the choice affecting the number of
# nonoverlapping observations
ontario_nonoverlap_curves <- function(day_start = c(1, 2, 3)[1]) {

  # Avoid function misusage
  stopifnot(day_start %in% 1:3 & length(day_start) == 1)

  # Split the five years measured in ontario
  split_years <- lapply(2010:2014, function(year) year(ontario$df$date) == year)

  # Ensure that the difference between observations is larger than 2 days
  # (since the observation for each day also includes the observation of
  # the next and previous day)
  index_years <- lapply(split_years, function(ind_year) {

    # Use the year days
    ydays <- yday(ontario$df$date[ind_year])
    ydays <- ydays - min(ydays) + 1

    # There are a lot of irregularities on the sequence od days due to weekends
    # and holidays, better use a sequential search
    ind <- day_start # The first day considered, on each year
    k <- 1
    repeat {
      j <- which(ydays > (ydays[ind[k]] + 2))[1] # First nonoverlaping curve
      if (is.na(j)) break # Stop once no more nonoverlapping curves are found
      ind <- c(ind, j)
      k <- k + 1 
    }

    # Return the positions of nonoverlapping curves
    return(which(ind_year)[ind])

  })

  # Concatenate the indexes for all years
  return(unlist(index_years))

}

### Ontario dataset

## Data

# Load dataset
data("ontario")

# Plot
if (save_fig) pdf(file = "ontario_temp.pdf", width = 14, height = 7)
plot(ontario$temp, main = "Ontario temperature (2010-2014)", axes = FALSE)
axis(1, at = seq(-24, 48, by = 6)); axis(2); box()
if (save_fig) dev.off()
if (save_fig) pdf(file = "ontario_elec.pdf", width = 7, height = 7)
plot(ontario$elec, main = "Ontario electricity consumption (2010-2014)",
     axes = FALSE)
axis(1, at = seq(0, 24, by = 6)); axis(2, at = seq(12, 26, by = 2)); box()
if (save_fig) dev.off()

## Test for the FLMFR

# With the data used in Benatia et al. (2017)
set.seed(123456789)
(flmfr_ontario_1 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr_l1s", save_fit_flm = TRUE,
                             save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_ontario_2 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr", save_fit_flm = FALSE,
                             save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_ontario_3 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr_l2", save_fit_flm = FALSE,
                             save_boot_stats = FALSE))
# Emphatic rejections

# With iid-ish data (curves without overlapping recordings)
iid <- ontario_nonoverlap_curves(day_start = 3)
set.seed(123456789)
(flmfr_ontario_iid_1 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr_l1s",
                                 save_fit_flm = TRUE, save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_ontario_iid_2 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr",
                                 save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_ontario_iid_3 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr_l2",
                                 save_fit_flm = FALSE, save_boot_stats = FALSE))
# Emphatic rejections (same for day_start = 2 or day_start = 3)

# Visualize estimate
if (save_fig) pdf(file = "ontario_beta.pdf", width = 15.1, height = 6)
lev <- seq(-0.03, 0.07, by = 0.01)
filled.contour(x = ontario$temp$argvals, y = ontario$elec$argvals,
               z = flmfr_ontario_1$fit_flm$Beta_hat, color.palette = viridis,
               levels = lev, asp = 1,
               plot.axes = {
                 axis(1, at = seq(-24, 48, by = 6))
                 axis(2, at = seq(0, 24, by = 6))
                 box()
                 abline(v = seq(-48, 48, by = 24))
                 abline(a = 0, b = 1)
               }, key.axes = {
                 axis(4, at = lev, cex.axis = 1.5)
               }, xlab = "Temperature", ylab = "Electricity consumption")
if (save_fig) dev.off()

## Test for significance

# With the data used in Benatia et al. (2017)
set.seed(123456789)
(noeff_ontario_1 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr_l1s", beta0 = 0,
                             save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_ontario_2 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr", beta0 = 0,
                             save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_ontario_3 <- flm_test(X = ontario$temp, Y = ontario$elec, B = B,
                             est_method = "fpcr_l2", beta0 = 0,
                             save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(kok_ontario <- kokoszka_test(X_fpc = flmfr_ontario_1$fit_flm$X_fpc,
                              Y_fpc = flmfr_ontario_1$fit_flm$Y_fpc))
# Emphatic rejections

# With iid-ish data (curves without overlapping recordings)
set.seed(123456789)
(noeff_ontario_iid_1 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr_l1s", beta0 = 0,
                                 save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_ontario_iid_2 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr", beta0 = 0,
                                 save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_ontario_iid_3 <- flm_test(X = ontario$temp[iid], Y = ontario$elec[iid],
                                 B = B, est_method = "fpcr_l2", beta0 = 0,
                                 save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(kok_ontario_iid <- kokoszka_test(X_fpc = flmfr_ontario_iid_1$fit_flm$X_fpc,
                                  Y_fpc = flmfr_ontario_iid_1$fit_flm$Y_fpc))
# Emphatic rejections (same for day_start = 2 or day_start = 3)

### AEMET temperatures

## Data

# Load dataset
data("aemet_temp")

# Partition the dataset in the first and last 20 years
with(aemet_temp, {
  ind_pred <- which((1974 <= df$year) & (df$year <= 1993))
  ind_resp <- which((1994 <= df$year) & (df$year <= 2013))
  aemet_temp_pred <<- list("df" = df[ind_pred, ], "temp" = temp[ind_pred])
  aemet_temp_resp <<- list("df" = df[ind_resp, ], "temp" = temp[ind_resp])
})

# Average the temperature on each period
mean_aemet <- function(x) {
  m <- tapply(X = 1:nrow(x$temp$data), INDEX = x$df$ind,
              FUN = function(i) colMeans(x$temp$data[i, , drop = FALSE],
                                         na.rm = TRUE))
  x$temp$data <- do.call(rbind, m)
  return(x$temp)
}

# Build predictor and response fdatas
aemet_temp_pred <- mean_aemet(aemet_temp_pred)
aemet_temp_resp <- mean_aemet(aemet_temp_resp)

# Smooth the data
aemet_temp_pred_smooth <-
  fda.usc::optim.np(fdataobj = aemet_temp_pred)$fdata.est
aemet_temp_resp_smooth <-
  fda.usc::optim.np(fdataobj = aemet_temp_resp)$fdata.est

# Plot raw data
if (save_fig) pdf(file = "aemet_pred.pdf", width = 7, height = 7)
plot(aemet_temp_pred, main = "AEMET temperature (1974-1993)",
     axes = FALSE, ylim = c(-2.5, 30), xlab = "")
months_axis(); axis(2); box()
if (save_fig) dev.off()
if (save_fig) pdf(file = "aemet_resp.pdf", width = 7, height = 7)
plot(aemet_temp_resp, main = "AEMET temperature (1994-2013)",
     axes = FALSE, ylim = c(-2.5, 30), xlab = "")
months_axis(); axis(2); box()
if (save_fig) dev.off()

# Plot smoothed data
plot(aemet_temp_pred_smooth, main = "AEMET temperature (1974-1993)",
     axes = FALSE, ylim = c(-2.5, 30), xlab = "")
months_axis(); axis(2); box()
plot(aemet_temp_resp_smooth, main = "AEMET temperature (1994-2013)",
     axes = FALSE, ylim = c(-2.5, 30), xlab = "")
months_axis(); axis(2); box()

# Average temperatures on both periods
if (save_fig) pdf("aemet_means.pdf", width = 7, height = 7)
plot(func_mean(aemet_temp_pred), main = "AEMET average temperature",
     axes = FALSE, ylim = c(-2.5, 30), xlab = "")
plot(func_mean(aemet_temp_resp), col = 2, add = TRUE)
months_axis(); axis(2); box()
legend("topleft", legend = c("1974-1993", "1994-2013"), col = 1:2, lwd = 2,
       cex = 1.5)
if (save_fig) dev.off()

# Smoothed average temperatures on both periods
plot(func_mean(aemet_temp_pred_smooth),
     main = "AEMET smoothed average temperature", axes = FALSE,
     ylim = c(-2.5, 30), xlab = "")
plot(func_mean(aemet_temp_resp_smooth), col = 2, add = TRUE)
months_axis(); axis(2); box()
legend("topleft", legend = c("1974-1993", "1994-2013"), col = 1:2, lwd = 2,
       cex = 1.5)

## Test for FLMFR

# Raw data
set.seed(123456789)
(flmfr_aemet_1 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr_l1s",
                           save_fit_flm = TRUE, save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_aemet_2 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr", save_fit_flm = FALSE,
                           save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_aemet_3 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr_l2", save_fit_flm = FALSE,
                           save_boot_stats = FALSE))

# Smooth data
set.seed(123456789)
(flmfr_aemet_smooth_1 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr_l1s",
                                  save_fit_flm = TRUE,
                                  save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_aemet_smooth_2 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr",
                                  save_fit_flm = FALSE,
                                  save_boot_stats = FALSE))
set.seed(123456789)
(flmfr_aemet_smooth_3 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr_l2",
                                  save_fit_flm = FALSE,
                                  save_boot_stats = FALSE))
# Clear no rejections in all of them

# Visualize raw estimate
if (save_fig) pdf(file = "aemet_beta.pdf", width = 7.8, height = 7,
                  cex.spe = 1.4)
cols <- c("blue", "white", "red")
lev <- seq(-0.022, 0.022, by = 0.004)
out <- 1:2
n_lev <- length(lev)
lev <- lev[-out]
filled.contour(x = aemet_temp_pred$argvals, y = aemet_temp_resp$argvals,
               z = flmfr_aemet_1$fit_flm$Beta_hat,
               col = colorRampPalette(colors = cols)(n_lev - 1)[-out],
               levels = lev, asp = 1,
               plot.axes = {
                 months_axis(side = 1, text = "Temperature in 1974-1993",
                             cex.spe = 1.5, line = 3.75)
                 months_axis(side = 2, text = "Temperature in 1994-2013",
                             cex.spe = 1.5, line = 3.5)
                 box()
                 for (k in seq(-365, 365, by = 365)) {
                   abline(a = k, b = 1)
                   abline(a = k - 90, b = 1, lty = 2)
                   abline(a = k + 90, b = 1, lty = 2)
                 }
               }, key.axes = {
                 axis(4, at = lev, cex.axis = 1.25)
               }, xlab = "", ylab = "")
if (save_fig) dev.off()

# Visualize smooth estimate
cols <- c("blue", "white", "red")
lev <- seq(-0.011, 0.011, by = 0.002)
out <- 1:3
n_lev <- length(lev)
lev <- lev[-out]
filled.contour(x = aemet_temp_pred$argvals, y = aemet_temp_resp$argvals,
               z = flmfr_aemet_smooth_1$fit_flm$Beta_hat,
               col = colorRampPalette(colors = cols)(n_lev - 1)[-out],
               levels = lev, asp = 1,
               plot.axes = {
                 months_axis(side = 1)
                 months_axis(side = 2)
                 box()
                 for (k in seq(-365, 365, by = 365)) {
                   abline(a = k, b = 1)
                   abline(a = k - 90, b = 1, lty = 2)
                   abline(a = k + 90, b = 1, lty = 2)
                 }
               }, key.axes = {
                 axis(4, at = lev, cex = 1.25)
               }, xlab = "Temperature in 1974-1993",
               ylab = "Temperature in 1994-2013")

## Test for significance

# Raw data
set.seed(123456789)
(noeff_aemet_1 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr_l1s", beta0 = 0,
                           save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_aemet_2 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr", beta0 = 0,
                           save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(noeff_aemet_3 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                           B = B, est_method = "fpcr_l2", beta0 = 0,
                           save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(kok_aemet <- kokoszka_test(X_fpc = flmfr_aemet_1$fit_flm$X_fpc,
                            Y_fpc = flmfr_aemet_1$fit_flm$Y_fpc))

# Smooth data
set.seed(123456789)
(noeff_aemet_smooth_1 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr_l1s",
                                  beta0 = 0, save_fit_flm = FALSE,
                                  save_boot_stats = FALSE))
set.seed(123456789)
(noeff_aemet_smooth_2 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr",
                                  beta0 = 0, save_fit_flm = FALSE,
                                  save_boot_stats = FALSE))
set.seed(123456789)
(noeff_aemet_smooth_3 <- flm_test(X = aemet_temp_pred_smooth,
                                  Y = aemet_temp_resp_smooth,
                                  B = B, est_method = "fpcr_l2",
                                  beta0 = 0, save_fit_flm = FALSE,
                                  save_boot_stats = FALSE))
set.seed(123456789)
(kok_aemet_smooth <- kokoszka_test(X_fpc = flmfr_aemet_smooth_1$fit_flm$X_fpc,
                                   Y_fpc = flmfr_aemet_smooth_1$fit_flm$Y_fpc))
# Emphatic rejections in all of them

## Stationarity: beta0 = diag(1)

# Raw data
beta0 <- matrix(0, nrow = length(aemet_temp_pred$argvals),
                ncol = length(aemet_temp_pred$argvals))
mgcv::sdiag(beta0, k = 0) <- 1
set.seed(123456789)
(stat_0_aemet_1 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l1s", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_0_aemet_2 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_0_aemet_3 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l2", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))

# Smooth data
set.seed(123456789)
(stat_0_aemet_smooth_1 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l1s",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
set.seed(123456789)
(stat_0_aemet_smooth_2 <- flm_test(X = aemet_temp_pred_smooth,
                                 Y = aemet_temp_resp_smooth,
                                 B = B, est_method = "fpcr",
                                 beta0 = beta0, save_fit_flm = FALSE,
                                 save_boot_stats = FALSE))
set.seed(123456789)
(stat_0_aemet_smooth_3 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l2",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
# Rejections

## Stationarity: beta0 = c = mean(beta_hat)

# Raw data
beta0 <- mean(flmfr_aemet_1$fit_flm$Beta_hat)
set.seed(123456789)
(stat_c_aemet_1 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l1s", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_c_aemet_2 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_c_aemet_3 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l2", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))

# Smooth data
set.seed(123456789)
beta0 <- mean(flmfr_aemet_smooth_1$fit_flm$Beta_hat)
(stat_c_aemet_smooth_1 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l1s",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
set.seed(123456789)
(stat_c_aemet_smooth_2 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
set.seed(123456789)
(stat_c_aemet_smooth_3 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l2",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
# Rejections

## Stationarity: beta0 = periodic diagonal averaging

# Raw data
beta0 <- flmfr_aemet_1$fit_flm$Beta_hat
n <- nrow(beta0)
for (k in 0:(n - 1)) {

  m <- mean(c(mgcv::sdiag(beta0, k = k), mgcv::sdiag(beta0, k = -(n - k))),
            na.rm = TRUE)
  if (k < n - 1) mgcv::sdiag(beta0, k = k) <- m
  if (k > 0) mgcv::sdiag(beta0, k = -(n - k)) <- m

}
set.seed(123456789)
(stat_d_aemet_1 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l1s", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_d_aemet_2 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))
set.seed(123456789)
(stat_d_aemet_3 <- flm_test(X = aemet_temp_pred, Y = aemet_temp_resp,
                            B = B, est_method = "fpcr_l2", beta0 = beta0,
                            save_fit_flm = FALSE, save_boot_stats = FALSE))

# Smooth data
beta0 <- flmfr_aemet_smooth_1$fit_flm$Beta_hat
n <- nrow(beta0)
for (k in 0:(n - 1)) {

  m <- mean(c(mgcv::sdiag(beta0, k = k), mgcv::sdiag(beta0, k = -(n - k))),
            na.rm = TRUE)
  if (k < n - 1) mgcv::sdiag(beta0, k = k) <- m
  if (k > 0) mgcv::sdiag(beta0, k = -(n - k)) <- m

}
set.seed(123456789)
beta0 <- mean(flmfr_aemet_smooth_1$fit_flm$Beta_hat)
(stat_d_aemet_smooth_1 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l1s",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
set.seed(123456789)
(stat_d_aemet_smooth_2 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
set.seed(123456789)
(stat_d_aemet_smooth_3 <- flm_test(X = aemet_temp_pred_smooth,
                                   Y = aemet_temp_resp_smooth,
                                   B = B, est_method = "fpcr_l2",
                                   beta0 = beta0, save_fit_flm = FALSE,
                                   save_boot_stats = FALSE))
# Rejections

# Load RData for the following chunk, which is time consuming and not so
# easy to replicate due to the specificities of the installation of
# "fdapss_1.0.zip"
load(file = "anova-pss-application.RData")

# ## ANOVAs for stationarity
#
# # Raw data
# set.seed(123456789)
# group <- as.factor(rep(1:2, each = length(aemet_temp_pred)))
# (anv_aemet <- fda.usc::anova.RPm(object = c(aemet_temp_pred,
#                                            aemet_temp_resp),
#                                  formula = ~group,
#                                  data.fac = data.frame(group), nboot = B))
#
# # Smooth data
# set.seed(123456789)
# group <- as.factor(rep(1:2, each = length(aemet_temp_pred_smooth)))
# (anv_aemet_smooth <- fda.usc::anova.RPm(object = c(aemet_temp_pred_smooth,
#                                                    aemet_temp_resp_smooth),
#                                         formula = ~group,
#                                         data.fac = data.frame(group),
#                                         nboot = B))
#
# ### Test for significance with Patilea et al. (2016) test
#
# # The "fdapss_1.0.zip" was kindly provided by César Sánchez Sellero. It
# # contains Windows executables. Installation works for R version < 3.5 and
# # a Windows machine
# install.packages(pkg = "fdapss_1.0.zip", repos = NULL)
# library(fdapss)
#
# ## Ontario
#
# set.seed(123456789)
# (pat_ontario <- pss.test(X = ontario$temp$data, Y = ontario$elec$data,
#                          nb = B, pr_pca = 0.99, nq = 50, vec_alpha = 1,
#                          ch = 1))
#
# ## AEMET
#
# set.seed(123456789)
# (pat_aemet <- pss.test(X = aemet_temp_pred$data, Y = aemet_temp_resp$data,
#                        nb = B, pr_pca = 0.99, nq = 50, vec_alpha = 1,
#                        ch = 1))
#
# set.seed(123456789)
# (pat_aemet_smooth <- pss.test(X = aemet_temp_pred_smooth$data,
#                               Y = aemet_temp_resp_smooth$data,
#                               nb = B, pr_pca = 0.99, nq = 50, vec_alpha = 1,
#                               ch = 1))
