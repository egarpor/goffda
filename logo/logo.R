
# Hex sticker for the 'goffda' package: a sheaf of AEMET seasonal temperature
# curves (the functional data) over the bilinear functional-linear-model
# coefficient surface beta(s, t), in a blue-white-red diverging palette. Drawn
# directly in hexSticker's hexagon coordinates and clipped to the hexagon.

library(hexSticker)
library(ggplot2)
library(showtext) # geom_pkgname/geom_url need showtext for the Aller font
showtext_auto()

# Shared logo standards
font <- "Aller_Rg"
name_size <- 31.2
url_size <- 9.0
url_x <- 1.00
url_y <- 0.08
url_angle <- 30
hex_border <- 1.5
dpi <- 600

# Hexagon geometry (matches hexSticker::geom_hexagon)
hex_half <- sqrt(3) / 2 # half-width; vertices y in [0, 2]
in_hex <- function(x, y) {
  u <- x - 1; v <- y - 1
  abs(u) <= hex_half + 1e-9 & abs(v) + abs(u) / sqrt(3) <= 1 + 1e-9
}
map_x <- function(f) (1 - hex_half) + f * (2 * hex_half) # fraction -> hex x

# Curve sheaf: AEMET seasonal arches
data("aemet_temp", package = "goffda")
temp_raw <- aemet_temp$temp$data # 2892 x 365 (has some NA)
grd <- aemet_temp$temp$argvals # day-of-year grid

# Pick 26 complete records at random
complete <- which(rowSums(is.na(temp_raw)) == 0)
set.seed(20260721)
idx <- sort(sample(complete, 26))

# Smooth each daily record into a clean seasonal curve.
temp_smooth <- t(apply(temp_raw[idx, ], 1,
                       function(y) smooth.spline(grd, y, df = 8)$y))

# Map the smoothed curves into the hexagon as a coloured sheaf
warmth <- rowMeans(temp_smooth) # blue (cool) -> red (warm)
band <- c(0.42, 1.22) # vertical span of the sheaf (kept low)
fx <- (grd - min(grd)) / diff(range(grd))
fy <- (temp_smooth - min(temp_smooth)) / diff(range(temp_smooth))
curves <- do.call(rbind, lapply(seq_along(idx), function(k) data.frame(
  id = k,
  x = map_x(fx),
  y = band[1] + fy[k, ] * diff(band),
  w = warmth[k]
)))
curves <- curves[in_hex(curves$x, curves$y), ] # clip to the hexagon

# beta(s, t) field: diagonal red/blue kernel that fills the hexagon
gx <- seq(1 - hex_half, 1 + hex_half, length.out = 240)
gy <- seq(0, 2, length.out = 280)
bs <- expand.grid(x = gx, y = gy)
bump <- function(x, y, x0, y0, r) exp(-((x - x0)^2 + (y - y0)^2) / (2 * r^2))
bs$z <- bump(bs$x, bs$y, 0.62, 0.70, 0.62) - # red lobe (lower-left)
        bump(bs$x, bs$y, 1.38, 1.30, 0.62) # blue lobe (upper-right)
bs <- bs[in_hex(bs$x, bs$y), ] # clip to the hexagon

# Curve colour ramp (cool -> warm)
pal_curve <- c("#2166AC", "#4393C3", "grey80", "#D6604D", "#B2182B")

# Assemble the sticker in hex coordinates
hex_plot <- ggplot() +
  geom_hexagon(size = hex_border, fill = "#FBFBFB", color = NA) +
  geom_tile(data = bs, aes(x, y, fill = z), alpha = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, guide = "none") +
  geom_line(data = curves, aes(x, y, group = id, colour = w),
            linewidth = 0.5, alpha = 0.85, lineend = "round") +
  scale_colour_gradientn(colours = pal_curve, guide = "none") +
  geom_hexagon(size = hex_border, fill = NA, color = "#2166AC") +
  geom_pkgname("goffda", x = 1, y = 1.52, color = "#B2182B",
               family = font, size = name_size) +
  geom_url("github.com/egarpor/goffda", x = url_x, y = url_y,
           color = "#2166AC", family = font, size = url_size,
           angle = url_angle) +
  theme_sticker(size = hex_border)

# Render and write the sticker to logo/ and man/figures/
dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
save_sticker("logo/logo.png", hex_plot, dpi = dpi)
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)
