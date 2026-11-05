# make_icon.R: PaleoMapR icon, a globe with an ancient continent, a shallow sea and a fossil-locality pin ####
#
# Run in tools/: writes launcher/paleomapr.ico, the favicons in pkgdown/favicon and man/figures/logo.png.

blob <- function(cx, cy, r, amp, phase, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  rr <- r * (1 + amp[1] * cos(2 * t + phase[1]) + amp[2] * cos(3 * t + phase[2]) + amp[3] * cos(5 * t + phase[3]))
  list(x = cx + rr * cos(t), y = cy + rr * sin(t))
}

draw_icon <- function(file, size) {
  png(file, width = size, height = size, bg = "transparent")
  op <- par(mar = c(0, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i", asp = 1)
  t <- seq(0, 2 * pi, length.out = 300)
  # ocean: concentric circles from deep to lighter blue
  pal <- grDevices::colorRampPalette(c("#0f3d66", "#2f7fc0"))(30)
  for (k in 1:30) {
    r <- 0.47 * (1 - (k - 1) / 34)
    polygon(0.5 + r * cos(t) - 0.004 * k, 0.5 + r * sin(t) + 0.004 * k, col = pal[k], border = NA)
  }
  # shelf (light blue) and continent (green to brown)
  shelf <- blob(0.47, 0.5, 0.27, c(0.12, 0.08, 0.05), c(0.4, 1.1, 2.0))
  polygon(shelf$x, shelf$y, col = "#7cc3e8", border = NA)
  land <- blob(0.45, 0.52, 0.21, c(0.15, 0.1, 0.06), c(0.4, 1.1, 2.0))
  polygon(land$x, land$y, col = "#7fae5a", border = NA)
  hill <- blob(0.40, 0.58, 0.09, c(0.2, 0.1, 0.05), c(1, 2, 3))
  polygon(hill$x, hill$y, col = "#b08d57", border = NA)
  # inland sea joined to the ocean by a seaway (the editor's correction)
  sea <- blob(0.52, 0.44, 0.06, c(0.15, 0.05, 0.05), c(0, 1, 2))
  polygon(sea$x, sea$y, col = "#7cc3e8", border = NA)
  lw <- max(1, size / 40)
  lines(c(0.55, 0.62, 0.72), c(0.42, 0.37, 0.34), col = "#7cc3e8", lwd = lw * 2.2, lend = 1)
  # globe rim
  lines(0.5 + 0.47 * cos(t), 0.5 + 0.47 * sin(t), col = "#0b2a47", lwd = max(1, size / 64))
  # fossil-locality pin
  px <- 0.68
  py <- 0.62
  polygon(c(px - 0.07, px + 0.07, px), c(py, py, py - 0.17), col = "#d9442b", border = NA)
  polygon(px + 0.085 * cos(t), py + 0.085 * sin(t), col = "#d9442b", border = NA)
  polygon(px + 0.035 * cos(t), py + 0.035 * sin(t), col = "white", border = NA)
}

dir.create("icon_tmp", showWarnings = FALSE)
sizes <- c(256, 128, 64, 48, 32, 16)
for (s in sizes) {
  draw_icon(file.path("icon_tmp", sprintf("icon_%d.png", s)), s)
}
library(magick)
im <- image_read(file.path("icon_tmp", sprintf("icon_%d.png", sizes)))
image_write(im, "../launcher/paleomapr.ico", format = "ico")
dir.create("../pkgdown/favicon", recursive = TRUE, showWarnings = FALSE)
dir.create("../man/figures", recursive = TRUE, showWarnings = FALSE)
image_write(image_read(file.path("icon_tmp", "icon_16.png")), "../pkgdown/favicon/favicon.ico", format = "ico")
draw_icon("../pkgdown/favicon/favicon-96x96.png", 96)
draw_icon("../pkgdown/favicon/apple-touch-icon.png", 180)
draw_icon("../pkgdown/favicon/web-app-manifest-192x192.png", 192)
draw_icon("../pkgdown/favicon/web-app-manifest-512x512.png", 512)
draw_icon("../man/figures/logo.png", 240)
unlink("icon_tmp", recursive = TRUE)
