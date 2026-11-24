# Four seasons on one page ####
#
# Figures with the four seasonal means (DJF, MAM, JJA, SON) of a climate layer side by side and one colour scale
# (as in De Vleeschouwer et al. 2014, Global and Planetary Change 120, 65-80). Used by make_maps(seasons4 = TRUE)
# and the Figure tab's "Four seasons" option.

#' one page with the four seasons (DJF, MAM, JJA, SON) of a layer, one colour scale
draw_season_page <- function(age, model, g, cfg_a, cs, ocean_pct, corrected, borders, modern_coast, indicators,
                             recs, user_points, sites, center, equator, cfg, climates, detail = "auto",
                             ice_front = TRUE, bar_pos = "bottom") {
  par(fig = .pf(c(0, 1, 0, 1)), mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  c1 <- climates[[1]]
  text(0.02, 0.985, sprintf("%s Ma — %s", format(age), stage_label(age, cfg)), adj = c(0, 1), font = 2, cex = 1.15)
  region <- if (!is.null(cfg_a$view)) {
    paste0(" · region: ", cfg_a$view$label)
  } else {
    paste0(" · ", PROJ_LABEL[[.proj_name(cfg_a)]], " projection")
  }
  text(0.02, 0.945, sprintf("%s (Scotese & Wright 2018), PALEOMAP frame · ocean %.1f%% of surface%s",
      if (corrected) { "Fossil-corrected PaleoDEM" } else { "PaleoDEM" }, ocean_pct, region),
      adj = c(0, 1), cex = 0.68, col = "#52514e")
  if (!is.null(c1$f)) {
    .caption_text(sub(", (December-February|annual mean)", "", climate_caption(c1)), 0.92, 0.6)
  }
  boxes <- list(c(0.01, 0.495, 0.49, 0.87), c(0.505, 0.99, 0.49, 0.87), c(0.01, 0.495, 0.11, 0.48), c(0.505, 0.99, 0.11, 0.48))
  labs <- c(djf = "DJF", mam = "MAM", jja = "JJA", son = "SON")
  for (k in seq_along(climates)) {
    par(fig = .pf(boxes[[k]]), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    draw_map_panel(age, model, g$z, g$lat, g$lon, cfg_a, corrected, borders, modern_coast, indicators, cs, recs,
        user_points = user_points, sites = sites, legend = k == 1, width_px = round(0.48 * 11 * cfg$map_dpi),
        equator = equator, detail = detail, ice_front = ice_front, climate = climates[[k]])
    u <- par("usr")
    text(u[1] + 0.02 * diff(u[1:2]), u[4] - 0.03 * diff(u[3:4]), labs[[climates[[k]]$season]], adj = c(0, 1), font = 2, cex = 0.9)
  }
  # the colour bar names the four seasons; CLIMATE_SEASONS lives in the editor's own environment
  env <- environment(draw_season_page)
  old <- get("CLIMATE_SEASONS", envir = env)
  on.exit(assign("CLIMATE_SEASONS", old, envir = env), add = TRUE)
  cs4 <- old
  cs4[["seasons4"]] <- "seasonal means (DJF, MAM, JJA, SON)"
  assign("CLIMATE_SEASONS", cs4, envir = env)
  cb <- c1
  cb$season <- "seasons4"
  par(fig = .pf(c(0.25, 0.75, 0.065, 0.083)), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  climate_colour_bar(1, cb)
}

#' the four seasonal climates of a figure with one common colour scale
.season_climates <- function(make_one) {
  seasons <- c("djf", "mam", "jja", "son")
  cl <- vector("list", 4)
  for (k in 1:4) {
    cl[[k]] <- make_one(seasons[k])
  }
  lay <- cl[[1]]$layer
  if (isTRUE(lay %in% c("temperature_change", "precipitation_change"))) {
    lims <- numeric(4)
    for (k in 1:4) {
      lims[k] <- cl[[k]]$lim %||% CHANGE_LIMITS[[lay]]
    }
    for (k in 1:4) {
      cl[[k]]$lim <- max(lims)
    }
  }
  if (!is.null(cl[[1]]$scale) && identical(cl[[1]]$scale$type, "div")) {
    lims <- numeric(4)
    for (k in 1:4) {
      lims[k] <- cl[[k]]$scale$lim %||% 0
    }
    for (k in 1:4) {
      cl[[k]]$scale$lim <- max(lims)
    }
  }
  cl
}
