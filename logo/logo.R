# Generates the goffda hex sticker (logo/logo.png + man/figures/logo.png).
#
# Motif: a sheaf of AEMET seasonal temperature curves (the functional data)
# over the bilinear functional-linear-model coefficient surface beta(s, t)
# (the linear model), in the package's blue-white-red diverging palette.
#
# Everything is drawn directly in hexSticker's hexagon coordinate system
# (pointy-top hexagon centred at (1, 1)) and clipped to the hexagon, so the
# content fills the whole hex and the hex outline itself does the cropping.
#
# Run from the package root:  Rscript logo/logo.R
# Output: logo/logo.png (master) and man/figures/logo.png (shipped mirror).

library(hexSticker)
library(ggplot2)
library(showtext)
showtext_auto()

## ---- Shared logo standard (identical across egarpor packages) -------------
# goffda draws in the hexagon coordinate system (geom_pkgname/geom_url), so it
# keeps showtext for the Aller font; the shared constants below match the other
# egarpor package logos.
FONT    <- "Aller_Rg"   # typeface for the package name and the GitHub URL
P_SIZE  <- 31.2         # package-name size (shared across packages)
U_SIZE  <- 9.0          # GitHub URL size (large enough to read)
U_X     <- 1.00         # GitHub URL position: along the lower-right hex edge
U_Y     <- 0.08
U_ANGLE <- 30
H_SIZE  <- 1.5          # hexagon border thickness
DPI     <- 600

## Hexagon geometry (matches hexSticker::geom_hexagon) ----------------------
HX <- sqrt(3) / 2                                # half-width; vertices y in [0, 2]
in_hex <- function(x, y) {                       # point inside the hexagon?
  u <- x - 1; v <- y - 1
  abs(u) <= HX + 1e-9 & abs(v) + abs(u) / sqrt(3) <= 1 + 1e-9
}
map_x <- function(f) (1 - HX) + f * (2 * HX)     # unit fraction -> hex x span

## 1. Curve sheaf: AEMET seasonal arches ------------------------------------

data("aemet_temp", package = "goffda")
Y   <- aemet_temp$temp$data                      # 2892 x 365 (has some NA)
grd <- aemet_temp$temp$argvals                   # day-of-year grid

complete <- which(rowSums(is.na(Y)) == 0)
set.seed(20260721)
idx <- sort(sample(complete, 26))

# Smooth each daily record into a clean seasonal curve (functional data).
Ys <- t(apply(Y[idx, ], 1, function(y) smooth.spline(grd, y, df = 8)$y))

warmth <- rowMeans(Ys)                           # blue (cool) -> red (warm)
band   <- c(0.42, 1.22)                          # vertical span of the sheaf
                                                 # (kept low so "goffda" breathes)
fx <- (grd - min(grd)) / diff(range(grd))
fy <- (Ys - min(Ys)) / diff(range(Ys))
curves <- do.call(rbind, lapply(seq_along(idx), function(k) data.frame(
  id = k,
  x  = map_x(fx),
  y  = band[1] + fy[k, ] * diff(band),
  w  = warmth[k]
)))
curves <- curves[in_hex(curves$x, curves$y), ]   # clip to the hexagon

## 2. beta(s, t) field: diagonal red/blue kernel that fills the hexagon -----

gx <- seq(1 - HX, 1 + HX, length.out = 240)
gy <- seq(0, 2, length.out = 280)
bs <- expand.grid(x = gx, y = gy)
bump <- function(x, y, x0, y0, r) exp(-((x - x0)^2 + (y - y0)^2) / (2 * r^2))
bs$z <- bump(bs$x, bs$y, 0.62, 0.70, 0.62) -     # red positive lobe (lower-left)
        bump(bs$x, bs$y, 1.38, 1.30, 0.62)       # blue negative lobe (upper-right)
bs <- bs[in_hex(bs$x, bs$y), ]                   # clip to the hexagon

## 3. Assemble the sticker in hex coordinates -------------------------------

pal_curve <- c("#2166AC", "#4393C3", "grey80", "#D6604D", "#B2182B")

sticker <- ggplot() +
  geom_hexagon(size = H_SIZE, fill = "#FBFBFB", color = NA) +       # base fill
  geom_tile(data = bs, aes(x, y, fill = z), alpha = 0.85) +         # beta field
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, guide = "none") +
  geom_line(data = curves, aes(x, y, group = id, colour = w),       # curves
            linewidth = 0.5, alpha = 0.85, lineend = "round") +
  scale_colour_gradientn(colours = pal_curve, guide = "none") +
  geom_hexagon(size = H_SIZE, fill = NA, color = "#2166AC") +       # border on top
  geom_pkgname("goffda", x = 1, y = 1.52, color = "#B2182B",
               family = FONT, size = P_SIZE) +
  geom_url("github.com/egarpor/goffda", x = U_X, y = U_Y,
           color = "#2166AC", family = FONT, size = U_SIZE, angle = U_ANGLE) +
  theme_sticker(size = H_SIZE)

dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
save_sticker("logo/logo.png", sticker, dpi = DPI)
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)

message("Wrote logo/logo.png and man/figures/logo.png")
