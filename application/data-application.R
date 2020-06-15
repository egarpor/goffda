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
  
  pdf <- function(..., cex.spe = 1.5) {
    
    grDevices::pdf(...)
    par(cex.axis = cex.spe, cex.main = cex.spe, cex.lab = cex.spe)
    
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

# Test for no functional effects by Lee et al. (2020)
FMDD_test <- function(X_fdata, Y_fdata, int_rule = "trapezoid", B = 1e3) {
  
  # Check if X_fdata and Y_fdata are functional data as fdata
  stopifnot(fda.usc::is.fdata(X_fdata))
  stopifnot(fda.usc::is.fdata(Y_fdata))
  
  # Sample size
  n <- nrow(X_fdata[["data"]])
  
  # Check if sample sizes are matched
  if (n != nrow(Y_fdata[["data"]])) {
    stop("X_fdata and Y_fdata samples must have the same number of elements")
  }
  
  # Grid points in which our functional samples are valued
  s <- X_fdata[["argvals"]]
  t <- Y_fdata[["argvals"]]
  
  # Check if grids are equispaced
  eps <- sqrt(.Machine[["double.eps"]])
  equispaced_x <- all(abs(diff(s, differences = 2)) < eps)
  equispaced_y <- all(abs(diff(t, differences = 2)) < eps)
  
  # Building matrices A and B
  a_ij <- b_ij <- A_mat <- B_mat <- matrix(0, n, n)
  for (i in 1:n) {

    a_ij[i, ] <- sqrt(apply(t(t(X_fdata[["data"]]) - X_fdata[["data"]][i, ])^2,
                            1, FUN = "integral1D", s, int_rule, equispaced_x))
    b_ij[i, ] <- apply(t(t(Y_fdata[["data"]]) - Y_fdata[["data"]][i, ])^2,
                       1, FUN = "integral1D", t, int_rule, equispaced_y)/2
  }
  
  a_i_dot <- (1 / (n - 2)) * rowSums(a_ij)
  a_dot_j <- (1 / (n - 2)) * colSums(a_ij)
  a_dot_dot <- (1 / ((n - 2) * (n - 1))) * sum(sum(a_ij))
  
  b_i_dot <- (1 / (n - 2)) * rowSums(b_ij)
  b_dot_j <- (1 / (n - 2)) * colSums(b_ij)
  b_dot_dot <- (1 / ((n - 2) * (n - 1))) * sum(sum(b_ij))
  
  # A and B matrices: null elements in the diagonal
  A_mat <- a_ij - matrix(rep(a_i_dot, n), nrow = n) -
    t(matrix(rep(a_dot_j, n), nrow = n)) + a_dot_dot
  B_mat <- b_ij - matrix(rep(b_i_dot, n), nrow = n) -
    t(matrix(rep(b_dot_j, n), nrow = n)) + b_dot_dot
  diag(B_mat) <- diag(A_mat) <- 0
  
  # FMDD_n statistic (sample version)
  FMDD_n <- (1 / (n * (n - 3))) * sum(sum(A_mat * B_mat))
  FMDD_orig_stat <- n * FMDD_n
  
  # Bootstrapped statistics
  phi <- (1 + sqrt(5))/2
  prob <- (phi + 2)/5
  FMDD_star <- FMDD_boot_stats <- numeric(B)
  for (b in 1:B) {
    
    V <- sample(x = c(1 - phi, phi), prob = c(prob, 1 - prob), size = n,
                replace = TRUE)
    
    FMDD_star[b] <- (1 / (n * (n - 3))) * sum(sum(outer(V, V, "*") *
                                                    A_mat * B_mat))
    FMDD_boot_stats[b] <- n * FMDD_star[b]
    
  }
  
  # p-value
  p_value <- mean(FMDD_orig_stat <= FMDD_boot_stats)
  
  # Output
  result <-
    structure(list(statistic = c(statistic = FMDD_orig_stat), p.value = p_value,
                   boot_statistics = FMDD_boot_stats,
                   method =
                     paste0("FMDDD (Functional martingale difference ",
                            "divergence)-based FLMFR significance test"),
                   A_mat = A_mat, a_ij = a_ij, a_i_dot = a_i_dot,
                   a_dot_j = a_dot_j, a_dot_dot = a_dot_dot,
                   B_mat = B_mat, b_ij = b_ij, b_i_dot = b_i_dot,
                   b_dot_j = b_dot_j, b_dot_dot = b_dot_dot))
  class(result) <- "htest"
  return(result)
  
}


# Add custom month axis
months_axis <- function(side = 1, ...) {
  
  axis(side = side, at = pmin(0.5 + c(0, cumsum(days_in_month(1:12))), 364.5),
       labels = c(month.abb, "Jan"), las = 2, ...)
  
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

# Visualize estimate
if (save_fig) pdf(file = "ontario_beta.pdf", width = 15, height = 6.03)
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
                 axis(4, at = lev)
               }, xlab = "Temperature", ylab = "Electricity consumption")
if (save_fig) dev.off()

## Tests for significance

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

# Kokoszka et al.
set.seed(123456789)
(kok_ontario <- kokoszka_test(X_fpc = flmfr_ontario_1$fit_flm$X_fpc,
                              Y_fpc = flmfr_ontario_1$fit_flm$Y_fpc))

# Lee et al.: statistic = 230564, p-value < 2.2e-16
set.seed(123456789)
(lee_ontario <-
    FMDD_test(X_fdata = ontario$temp, Y_fdata = ontario$elec, B = B))



# Emphatic rejections

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
aemet_temp_pred_smooth <- fda.usc::min.np(fdataobj = aemet_temp_pred)$fdata.est
aemet_temp_resp_smooth <- fda.usc::min.np(fdataobj = aemet_temp_resp)$fdata.est

# Plot raw data
if (save_fig) pdf(file = "aemet_pred.pdf", width = 7, height = 7)
plot(aemet_temp_pred, main = "AEMET temperature (1974-1993)",
     axes = FALSE, ylim = c(-2.5, 30))
months_axis(); axis(2); box()
if (save_fig) dev.off()
if (save_fig) pdf(file = "aemet_resp.pdf", width = 7.1, height = 7.1)
plot(aemet_temp_resp, main = "AEMET temperature (1994-2013)",
     axes = FALSE, ylim = c(-2.5, 30))
months_axis(); axis(2); box()
if (save_fig) dev.off()

# Plot smoothed data
plot(aemet_temp_pred_smooth, main = "AEMET temperature (1974-1993)",
     axes = FALSE, ylim = c(-2.5, 30))
months_axis(); axis(2); box()
plot(aemet_temp_resp_smooth, main = "AEMET temperature (1994-2013)",
     axes = FALSE, ylim = c(-2.5, 30))
months_axis(); axis(2); box()

# Average temperatures on both periods
if (save_fig) pdf("aemet_means.pdf", width = 7.1, height = 7.1)
plot(func_mean(aemet_temp_pred), main = "AEMET average temperature",
     axes = FALSE, ylim = c(-2.5, 30))
plot(func_mean(aemet_temp_resp), col = 2, add = TRUE)
months_axis(); axis(2); box()
legend("topleft", legend = c("1974-1993", "1994-2013"), col = 1:2, lwd = 2)
if (save_fig) dev.off()

# Smoothed average temperatures on both periods
plot(func_mean(aemet_temp_pred_smooth), main = "AEMET smoothed average temperature",
     axes = FALSE, ylim = c(-2.5, 30))
plot(func_mean(aemet_temp_resp_smooth), col = 2, add = TRUE)
months_axis(); axis(2); box()
legend("topleft", legend = c("1974-1993", "1994-2013"), col = 1:2, lwd = 2)

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
if (save_fig) pdf(file = "aemet_beta.pdf", width = 7.6, height = 7,
                  cex.spe = 1.2)
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

# Kokoszka et al.
set.seed(123456789)
(kok_aemet <- kokoszka_test(X_fpc = flmfr_aemet_1$fit_flm$X_fpc,
                            Y_fpc = flmfr_aemet_1$fit_flm$Y_fpc))

# Lee et al.: statistic = 9112344, p-value < 2.2e-16
set.seed(123456789)
(lee_aemet <-
    FMDD_test(X_fdata = aemet_temp_pred, Y_fdata = aemet_temp_resp, B = B))



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

# Kokozska et al.
set.seed(123456789)
(kok_aemet_smooth <- kokoszka_test(X_fpc = flmfr_aemet_smooth_1$fit_flm$X_fpc,
                                   Y_fpc = flmfr_aemet_smooth_1$fit_flm$Y_fpc))

# Lee et al.
set.seed(123456789)
(lee_aemet_smooth <- FMDD_test(X_fdata = aemet_temp_pred_smooth,
                               Y_fdata = aemet_temp_resp_smooth, B = B))



# Emphatic rejections in all of them

## Stationarity: beta0 = diag(1)

# Raw data
beta0 <- matrix(0, nrow = length(aemet_temp_pred$argvals),
                ncol = length(aemet_temp_pred$argvals))
sdiag(beta0, k = 0) <- 1
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
  
  m <- mean(c(sdiag(beta0, k = k), sdiag(beta0, k = -(n - k))), na.rm = TRUE)
  if (k < n - 1) sdiag(beta0, k = k) <- m
  if (k > 0) sdiag(beta0, k = -(n - k)) <- m
  
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
  
  m <- mean(c(sdiag(beta0, k = k), sdiag(beta0, k = -(n - k))), na.rm = TRUE)
  if (k < n - 1) sdiag(beta0, k = k) <- m
  if (k > 0) sdiag(beta0, k = -(n - k)) <- m
  
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