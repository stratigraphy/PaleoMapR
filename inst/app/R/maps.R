`%||%` <- function(a, b) if (is.null(a)) b else a
# =============================================================================
#  Topographic palaeogeographic maps (Mollweide), in R
#   - hypsometric tints + bathymetry, hillshading
#   - corrected 0 m coastline, PaleoCoastlines shoreline (dashed)
#   - present-day countries (borders black, coastlines white; dotted = crust the
#     plate model adds later or approximate position), rotated with pyGPlates
#   - climate indicators (reefs, evaporites, glacial deposits) and ice front,
#     incl. your own points and edits
#   - localities (points/localities.csv: add / remove / show freely)
#  Output: one PNG per age and a PDF built from the PNGs (thin lines stay visible
#  in every PDF viewer).
# =============================================================================
library(sf); library(png)

# ---------------------------------------------------------------- colours ----
.sea_stops  <- list(z = c(-4000, -3000, -1400, -600, -200, 0),
                    col = c("#0d2b4e", "#163f6c", "#2c6aa0", "#4f93c4", "#86bfe0", "#c4e3f1"))
.land_stops <- list(z = c(0, 150, 400, 800, 1300, 1900, 2500, 3000),
                    col = c("#9cc489", "#c4d696", "#e3d9a0", "#d6b277", "#b98a56", "#956a46", "#b9a89a", "#f4f1ec"))
.ramp <- function(st, z) {
  f <- (pmin(pmax(z, min(st$z)), max(st$z)) - min(st$z)) / (max(st$z) - min(st$z))
  pos <- (st$z - min(st$z)) / (max(st$z) - min(st$z))
  rgbm <- col2rgb(st$col) / 255
  sapply(1:3, function(k) approx(pos, rgbm[k, ], xout = f)$y)
}
#' z (matrix) -> RGB array [n, p, 3]
colourise <- function(z) {
  v <- as.vector(z); m <- .ramp(.sea_stops, v); land <- v >= 0
  m[land, ] <- .ramp(.land_stops, v[land])
  array(m, c(nrow(z), ncol(z), 3))
}
hillshade <- function(z, lat, az = 315, alt = 40, exag_land = 60, exag_sea = 8) {
  d <- abs(lat[2] - lat[1]) * pi / 180 * 6371000
  dx <- d * pmin(pmax(cos(lat * pi / 180), 0.05), 1)
  zz <- z * ifelse(z >= 0, exag_land, exag_sea)
  n <- nrow(z)
  gy <- rbind(zz[2, ] - zz[1, ], (zz[3:n, ] - zz[1:(n - 2), ]) / 2, zz[n, ] - zz[n - 1, ]) / d
  gx <- (shift(zz, 0, 1) - shift(zz, 0, -1)) / 2 / dx
  slope <- atan(sqrt(gx^2 + gy^2)); aspect <- atan2(-gx, gy)
  a <- az * pi / 180; b <- alt * pi / 180
  hs <- sin(b) * cos(slope) + cos(b) * sin(slope) * cos(a - aspect)
  pmin(pmax(hs, 0), 1)
}
render_rgb <- function(z, lat) {
  rgb <- colourise(z); hs <- hillshade(z, lat)
  k <- ifelse(z >= 0, 0.75, 0.45)
  for (c in 1:3) rgb[, , c] <- pmin(pmax(rgb[, , c] * (1 - k + k * (0.35 + 0.9 * hs)), 0), 1)
  rgb
}

#' the Map-editor grid of an age (plate model, colours with relief shading, 0 m line), kept in memory for the last
#' `keep` maps, so switching back to a map is instant; the key includes the file's time and size, so a rebuilt grid
#' is read again. The colours are kept as bytes (5 MB per map instead of 39).
.EDITOR_GRIDS <- new.env()
editor_grid <- function(file, age, model, label = "", cfg = SILUCARB, keep = 8) {
  inf <- file.info(file)
  key <- paste(normalizePath(file, winslash = "/", mustWork = FALSE), model, format(as.numeric(inf$mtime), digits = 15), inf$size)
  items <- .EDITOR_GRIDS$items %||% list()
  ord <- .EDITOR_GRIDS$order %||% character(0)
  g <- items[[key]]
  if (is.null(g)) {
    dem <- read_dem_nc(file)
    m <- cached_grid_to_model(age, dem, model, cfg)
    img <- render_rgb(m$z, m$lat)
    zz <- t(cbind(m$z, m$z[, 1]))
    coast <- grDevices::contourLines(c(m$lon, m$lon[1] + 360), m$lat, zz, levels = 0)
    g <- list(img = array(as.raw(round(255 * pmin(pmax(img, 0), 1))), dim(img)), z = m$z, lat = m$lat, lon = m$lon, coast = coast)
    items[[key]] <- g
  }
  ord <- c(setdiff(ord, key), key)
  if (length(ord) > keep) {
    drop <- ord[seq_len(length(ord) - keep)]
    items[drop] <- NULL
    ord <- setdiff(ord, drop)
  }
  .EDITOR_GRIDS$items <- items
  .EDITOR_GRIDS$order <- ord
  g$img <- array(as.integer(g$img) / 255, dim(g$img))
  g$label <- label
  g
}

#' the 0 m line of an editor grid (from editor_grid) as one x/y pair with NA breaks
editor_coast_xy <- function(coast) {
  if (!length(coast)) {
    return(list(x = numeric(0), y = numeric(0)))
  }
  xs <- vector("list", length(coast))
  ys <- vector("list", length(coast))
  for (i in seq_along(coast)) {
    xs[[i]] <- c(coast[[i]]$x, NA)
    ys[[i]] <- c(coast[[i]]$y, NA)
  }
  list(x = unlist(xs), y = unlist(ys))
}

# ------------------------------------------------------------- projection ----
#' map projections (make_maps(projection = ...), map_projection in R/config.R, or the app)
#'   whole-globe maps: name -> PROJ projection
PROJ_GLOBAL <- c(mollweide = "moll", robinson = "robin", equal_earth = "eqearth", winkel_tripel = "wintri",
                 hammer = "hammer", eckert4 = "eck4", plate_carree = "eqc")
#'   globe / polar views (a round map)
PROJ_DISC <- c("orthographic", "south_polar", "north_polar")
PROJ_LABEL <- c(mollweide = "Mollweide", robinson = "Robinson", equal_earth = "Equal Earth",
                winkel_tripel = "Winkel tripel", hammer = "Hammer", eckert4 = "Eckert IV",
                plate_carree = "Plate Carr\u00e9e (equirectangular)", orthographic = "orthographic (globe)",
                south_polar = "south polar (Lambert equal-area)", north_polar = "north polar (Lambert equal-area)")
.proj_name <- function(cfg) { p <- cfg$map_projection %||% "mollweide"
  if (!p %in% c(names(PROJ_GLOBAL), PROJ_DISC))
    stop("projection must be one of: ", paste(c(names(PROJ_GLOBAL), PROJ_DISC), collapse = ", "), call. = FALSE)
  p }
map_crs <- function(cfg = SILUCARB) {
  if (!is.null(cfg$view)) return(cfg$view$crs)
  sprintf("+proj=%s +lon_0=%g +datum=WGS84 +units=m +no_defs", PROJ_GLOBAL[[.proj_name(cfg)]], cfg$central_lon)
}
#' outline of a whole-globe map: the two edge meridians (x, y in map units, right edge first)
.global_outline <- function(cfg) {
  la <- seq(-90, 90, length.out = 361)
  sf_project("EPSG:4326", map_crs(cfg), rbind(cbind(cfg$central_lon + 179.9999, la),
                                             cbind(cfg$central_lon - 179.9999, rev(la))), keep = TRUE, warn = FALSE)
}
.EXT <- new.env()
.global_ext <- function(cfg) {
  k <- .proj_name(cfg)
  if (is.null(.EXT[[k]])) { e <- .global_outline(list(map_projection = k, central_lon = 0))
    .EXT[[k]] <- list(x = c(-1, 1) * max(abs(e[, 1]), na.rm = TRUE), y = c(-1, 1) * max(abs(e[, 2]), na.rm = TRUE),
                      ex = e[1:361, 1], ey = e[1:361, 2]) }
  .EXT[[k]]
}
view_xlim <- function(cfg) if (is.null(cfg$view)) .global_ext(cfg)$x else cfg$view$xlim
view_ylim <- function(cfg) if (is.null(cfg$view)) .global_ext(cfg)$y else cfg$view$ylim
.is_disc <- function(cfg) !is.null(cfg$view) && !is.null(cfg$view$disc)

#' palaeo lat/lon of a map centre: NULL -> (0, central_lon); one number -> (0, that longitude);
#' a locality name or c(lat, lon) present-day -> its palaeo-position at the age
center_latlon <- function(center, age, model, cfg = SILUCARB) {
  if (is.null(center)) return(c(0, cfg$central_lon))
  if (is.numeric(center) && length(center) == 1) return(c(0, center))
  if (is.character(center)) {
    loc <- read_localities(cfg); k <- match(center, loc$name)
    if (is.na(k)) stop("unknown locality '", center, "'")
    center <- c(loc$lat[k], loc$lon[k])
  }
  p <- gp_reconstruct(center[1], center[2], age, model)
  if (!is.finite(p$plat)) { warning("centre point not on a plate at ", age, " Ma; using default"); return(c(0, cfg$central_lon)) }
  c(p$plat, p$plon)
}

#' round views: an orthographic globe centred on `center`, or a polar map (Lambert
#' azimuthal equal-area) reaching to latitude cfg$polar_edge_lat (0 = the equator)
projection_view <- function(center, age, model, cfg = SILUCARB) {
  p <- .proj_name(cfg); if (!p %in% PROJ_DISC) return(NULL)
  R <- 6371000
  if (p == "orthographic") {
    c0 <- center_latlon(center, age, model, cfg)
    return(list(crs = sprintf("+proj=ortho +lat_0=%g +lon_0=%g +R=%d +units=m +no_defs", c0[1], c0[2], R),
                xlim = c(-R, R), ylim = c(-R, R), lat0 = c0[1], lon0 = c0[2], disc = R, far_deg = 89.5,
                label = sprintf("globe centred on %.0f\u00b0%s, %.0f\u00b0%s", abs(c0[1]), if (c0[1] < 0) "S" else "N",
                                abs(c0[2]), if (c0[2] < 0) "W" else "E")))
  }
  south <- p == "south_polar"; e <- cfg$polar_edge_lat %||% 0
  d <- if (south) 90 + e else 90 - e                     # angular radius of the map from the pole
  d <- min(max(d, 10), 175); r <- 2 * R * sin(d / 2 * pi / 180)
  lon0 <- map_center(center, age, model, cfg)
  list(crs = sprintf("+proj=laea +lat_0=%d +lon_0=%g +R=%d +units=m +no_defs", if (south) -90 else 90, lon0, R),
       xlim = c(-r, r), ylim = c(-r, r), lat0 = if (south) -90 else 90, lon0 = lon0, disc = r, far_deg = d,
       label = { el <- if (south) -90 + d else 90 - d
                 if (abs(el) < 0.5) "reaching to the equator" else sprintf("reaching to %.0f\u00b0%s", abs(el), if (el < 0) "S" else "N") })
}

#' a zoomed-in regional view (Lambert azimuthal equal-area centred on the region)
#' region = c(lat_min, lat_max, lon_min, lon_max)  palaeo-coordinates of the map age, or
#' region = list(center = "Huy" | c(lat, lon) present-day, width_km = 4000, height_km = 3000)
#'          (follows the locality through time)
#' region = list(center_palaeo = c(lat, lon), width_km, height_km)  a fixed palaeo-position (the editor)
.deg_ns <- function(v) sprintf("%g\u00b0%s", abs(v), if (v < 0) "S" else if (v > 0) "N" else "")
.deg_ew <- function(v) { v <- ((v + 180) %% 360) - 180; sprintf("%g\u00b0%s", abs(v), if (v < 0) "W" else if (v > 0) "E" else "") }
make_view <- function(region, age, model, cfg = SILUCARB) {
  if (is.null(region)) return(NULL)
  if (is.numeric(region) && length(region) == 4) {
    la <- sort(region[1:2]); lo1 <- region[3]; lo2 <- region[4]; if (lo2 < lo1) lo2 <- lo2 + 360
    c_la <- mean(la); c_lo <- ((mean(c(lo1, lo2)) + 180) %% 360) - 180
    crs <- sprintf("+proj=laea +lat_0=%g +lon_0=%g +datum=WGS84 +units=m +no_defs", c_la, c_lo)
    b <- rbind(cbind(seq(lo1, lo2, length.out = 60), la[1]), cbind(lo2, seq(la[1], la[2], length.out = 60)),
               cbind(seq(lo2, lo1, length.out = 60), la[2]), cbind(lo1, seq(la[2], la[1], length.out = 60)))
    b[, 1] <- ((b[, 1] + 180) %% 360) - 180
    xy <- sf_project("EPSG:4326", crs, b, warn = FALSE)
    return(list(crs = crs, xlim = range(xy[, 1], na.rm = TRUE), ylim = range(xy[, 2], na.rm = TRUE),
                lat0 = c_la, lon0 = c_lo, label = sprintf("%s\u2013%s, %s\u2013%s", .deg_ns(la[1]), .deg_ns(la[2]), .deg_ew(region[3]), .deg_ew(region[4]))))
  }
  ctr <- region$center
  if (!is.null(region$center_palaeo)) {             # centre given in palaeo-coordinates of this map
    p <- list(plat = region$center_palaeo[1], plon = region$center_palaeo[2])
  } else {
    if (is.character(ctr)) {
      loc <- read_localities(cfg); k <- match(ctr, loc$name)
      if (is.na(k)) stop("unknown locality '", ctr, "'")
      ctr <- c(loc$lat[k], loc$lon[k])
    }
    p <- gp_reconstruct(ctr[1], ctr[2], age, model)
    if (!is.finite(p$plat)) stop("region centre is not on a plate at ", age, " Ma")
  }
  w <- if (is.null(region$width_km)) 4000 else region$width_km
  h <- if (is.null(region$height_km)) 0.6 * w else region$height_km
  list(crs = sprintf("+proj=laea +lat_0=%g +lon_0=%g +datum=WGS84 +units=m +no_defs", p$plat, p$plon),
       xlim = c(-w, w) * 500, ylim = c(-h, h) * 500, lat0 = p$plat, lon0 = p$plon,
       label = sprintf("%s \u00d7 %s km around %s", format(round(w), big.mark = ","), format(round(h), big.mark = ","),
                       if (is.character(region$center)) region$center else if (!is.null(region$center_palaeo))
                         sprintf("%.1f\u00b0%s, %.1f\u00b0%s", abs(p$plat), if (p$plat < 0) "S" else "N", abs(p$plon), if (p$plon < 0) "W" else "E")
                       else "the chosen point"))
}
#' regional views: hide what lies more than ~100 deg from the centre (far side of the globe)
.far <- function(lon, lat, cfg) {
  if (is.null(cfg$view)) return(rep(FALSE, length(lon)))
  v <- cfg$view
  cosd <- sin(lat * pi / 180) * sin(v$lat0 * pi / 180) + cos(lat * pi / 180) * cos(v$lat0 * pi / 180) * cos((lon - v$lon0) * pi / 180)
  !is.na(cosd) & cosd < cos((v$far_deg %||% 100) * pi / 180)
}

#' RGB lat-lon image -> raster array in Mollweide (width px)
#' fast path of project_rgb() for the Figure-tab preview (equirectangular view around lon0 with R = 6371 km, as set by
#' figure_preview): each output pixel takes the bilinear mean of the four nearest grid cells (longitudes wrap around),
#' instead of a terra re-projection; ~10x faster. Cells beyond the poles are transparent.
.project_rgb_eqc <- function(rgb, lat, lon, width, cfg) {
  v <- cfg$view
  k <- 6371000 * pi / 180
  xl <- v$xlim
  yl <- v$ylim
  h <- max(10, round(width * diff(yl) / diff(xl)))
  xc <- seq(xl[1], xl[2], length.out = width + 1)
  xc <- (xc[-1] + xc[-length(xc)]) / 2
  yc <- seq(yl[2], yl[1], length.out = h + 1)
  yc <- (yc[-1] + yc[-length(yc)]) / 2
  dl <- abs(lat[2] - lat[1])
  nr <- length(lat)
  nc <- length(lon)
  fi <- (yc / k - lat[1]) / dl + 1                    # fractional row (lat ascending) of each output row
  fj <- (((v$lon0 + xc / k - lon[1]) / dl) %% nc) + 1  # fractional column, wrapped
  i0 <- floor(fi)
  wi <- fi - i0
  inside <- fi >= 0.5 & fi <= nr + 0.5
  i0 <- pmin(pmax(i0, 1), nr)
  i1 <- pmin(i0 + 1, nr)
  j0 <- floor(fj)
  wj <- fj - j0
  j0 <- ((j0 - 1) %% nc) + 1
  j1 <- (j0 %% nc) + 1
  WJ <- matrix(wj, h, width, byrow = TRUE)
  a <- array(1, c(h, width, 4))
  for (ch in 1:3) {
    m <- rgb[, , ch]
    top <- m[i0, j0, drop = FALSE] * (1 - WJ) + m[i0, j1, drop = FALSE] * WJ
    bot <- m[i1, j0, drop = FALSE] * (1 - WJ) + m[i1, j1, drop = FALSE] * WJ
    a[, , ch] <- top * (1 - wi) + bot * wi           # wi recycles down the rows
  }
  a[, , 4] <- matrix(ifelse(inside, 1, 0), h, width)
  a
}

project_rgb <- function(rgb, lat, lon, width, cfg = SILUCARB, wrap = TRUE) {
  if (isTRUE(cfg$view$fast_eqc)) {
    return(.project_rgb_eqc(rgb, lat, lon, width, cfg))
  }
  dl <- abs(lat[2] - lat[1]); e <- if (wrap) 1 else 0
  p <- dim(rgb)[2]; if (wrap) rgb <- rgb[, c(p, seq_len(p), 1), , drop = FALSE]   # wrap one column: no gap at 180 deg
  src <- rast(nrows = length(lat), ncols = length(lon) + 2 * e, nlyrs = 3,
              xmin = lon[1] - (1 + 2 * e) * dl / 2, xmax = lon[length(lon)] + (1 + 2 * e) * dl / 2,
              ymin = lat[1] - dl / 2, ymax = lat[length(lat)] + dl / 2, crs = "EPSG:4326")
  for (k in 1:3) values(src[[k]]) <- as.vector(t(rgb[nrow(rgb):1, , k]))
  xl <- view_xlim(cfg); yl <- view_ylim(cfg)
  tmpl <- rast(nrows = max(10, round(width * diff(yl) / diff(xl))), ncols = width, xmin = xl[1], xmax = xl[2],
               ymin = yl[1], ymax = yl[2], crs = map_crs(cfg), nlyrs = 3)
  pr <- suppressWarnings(project(src, tmpl, method = "bilinear"))   # cells beyond the poles: harmless
  a <- as.array(pr)
  # keep only what lies inside the map outline (ellipse, Robinson outline, globe disc ...)
  xc <- seq(xl[1], xl[2], length.out = dim(a)[2] + 1); xc <- (xc[-1] + xc[-length(xc)]) / 2
  yc <- seq(yl[2], yl[1], length.out = dim(a)[1] + 1); yc <- (yc[-1] + yc[-length(yc)]) / 2
  out <- if (is.null(cfg$view)) {
    E <- .global_ext(cfg); o <- order(E$ey); xe <- approx(E$ey[o], abs(E$ex[o]), yc, ties = max, rule = 2)$y
    outer(xe, abs(xc), "<")
  } else if (.is_disc(cfg)) outer(yc^2, xc^2, "+") > cfg$view$disc^2 else FALSE
  for (k in 1:3) { ak <- a[, , k]; ak[out] <- NA; a[, , k] <- ak }
  alpha <- ifelse(is.na(a[, , 1]), 0, 1)
  a[is.na(a)] <- 1
  array(c(a, alpha), c(dim(a)[1:2], 4))
}

#' longitude relative to the central meridian, with points lying exactly on the map edge
#' (e.g. a grid column at central_lon + 180) moved to the side of their neighbours, so a
#' line along the edge is never drawn across the whole map
.rel_lon <- function(lon, cfg) {
  rel <- ((lon - cfg$central_lon + 180) %% 360) - 180
  e <- !is.na(rel) & abs(rel) > 179.999
  if (any(e)) {
    s <- ifelse(e | is.na(rel), NA, sign(rel))
    fwd <- s; bwd <- s
    for (i in seq_along(s)[-1]) if (is.na(fwd[i]) && !is.na(rel[i]) && !is.na(rel[i - 1])) fwd[i] <- fwd[i - 1]
    for (i in rev(seq_along(s))[-1]) if (is.na(bwd[i]) && !is.na(rel[i]) && !is.na(rel[i + 1])) bwd[i] <- bwd[i + 1]
    sg <- ifelse(is.na(fwd), bwd, fwd); sg[is.na(sg)] <- -1
    rel[e] <- 179.999 * sg[e]
  }
  rel
}

#' lon/lat -> map x/y, with NA breaks where a line crosses the map edge
proj_lines <- function(lon, lat, cfg = SILUCARB) {
  rel <- .rel_lon(lon, cfg); lon <- cfg$central_lon + rel
  br <- c(FALSE, abs(diff(rel)) > 180) & !.is_disc(cfg)
  xy <- sf_project("EPSG:4326", map_crs(cfg), cbind(lon, lat), keep = TRUE, warn = FALSE)
  far <- .far(lon, lat, cfg); xy[far, ] <- NA
  if (any(br)) {
    idx <- cumsum(br); out <- NULL
    for (s in unique(idx)) out <- rbind(out, xy[idx == s, , drop = FALSE], c(NA, NA))
    xy <- out
  }
  xy
}
draw_ll <- function(lon, lat, ..., cfg = SILUCARB) { xy <- proj_lines(lon, lat, cfg); lines(xy[, 1], xy[, 2], ...) }
#' many lat/lon lines (list of data.frames with lat, lon) in one projected call
draw_ll_list <- function(L, cfg = SILUCARB, ...) {
  if (!length(L)) return(invisible())
  lon <- unlist(lapply(L, function(p) c(p$lon, NA))); lat <- unlist(lapply(L, function(p) c(p$lat, NA)))
  rel <- .rel_lon(lon, cfg); lon <- cfg$central_lon + rel
  jump <- c(FALSE, abs(diff(rel)) > 180); jump[is.na(jump)] <- FALSE
  if (.is_disc(cfg)) jump[] <- FALSE                # round maps have no edge meridian
  # insert a break (NA) before every point that jumps across the map edge
  k <- which(jump); if (length(k)) {
    o <- order(c(seq_along(lon), k - 0.5)); lon <- c(lon, rep(NA, length(k)))[o]; lat <- c(lat, rep(NA, length(k)))[o] }
  ok <- !is.na(lon); xy <- matrix(NA_real_, length(lon), 2)
  xy[ok, ] <- sf_project("EPSG:4326", map_crs(cfg), cbind(lon[ok], lat[ok]), keep = TRUE, warn = FALSE)
  xy[ok & .far(lon, lat, cfg), ] <- NA
  lines(xy[, 1], xy[, 2], ...)
}

map_frame <- function(cfg = SILUCARB, ...) {
  if (.is_disc(cfg)) {
    a <- list(...); t <- seq(0, 2 * pi, length.out = 721)
    polygon(cfg$view$disc * cos(t), cfg$view$disc * sin(t), border = a$border, lwd = a$lwd); return(invisible())
  }
  if (!is.null(cfg$view)) {
    a <- list(...); rect(cfg$view$xlim[1], cfg$view$ylim[1], cfg$view$xlim[2], cfg$view$ylim[2],
                         border = a$border, lwd = a$lwd); return(invisible())
  }
  e <- .global_outline(cfg)
  polygon(e[, 1], e[, 2], ...)
}
draw_graticule <- function(cfg = SILUCARB) {
  if (!is.null(cfg$view) && !.is_disc(cfg)) return(draw_graticule_regional(cfg))
  la <- seq(-90, 90, 1); lo <- seq(-180, 180, 1)
  if (.is_disc(cfg)) {                             # globe / polar: 30-degree graticule
    for (m in seq(-180, 150, 30)) draw_ll(rep(m, length(la)), la, col = adjustcolor("white", 0.55), lwd = 0.5, cfg = cfg)
    for (p in seq(-60, 60, 30)) draw_ll(lo, rep(p, length(lo)), col = adjustcolor("white", 0.55), lwd = 0.5, cfg = cfg)
    return(invisible())
  }
  for (m in seq(-180, 180, 60)) draw_ll(rep(m, length(la)), la, col = adjustcolor("white", 0.55), lwd = 0.5, cfg = cfg)
  for (p in seq(-60, 60, 30)) draw_ll(lo, rep(p, length(lo)), col = adjustcolor("white", 0.55), lwd = 0.5, cfg = cfg)
}
#' regional graticule (5/10/20 deg depending on the extent) with labels on the frame
draw_graticule_regional <- function(cfg) {
  xl <- cfg$view$xlim; yl <- cfg$view$ylim
  step <- { w <- diff(xl) / 1000; if (w < 2500) 5 else if (w < 7000) 10 else 20 }
  col <- adjustcolor("white", 0.55); la <- seq(-90, 90, 0.5); lo <- seq(-180, 180, 0.5)
  inside <- function(xy) !is.na(xy[, 1]) & xy[, 1] >= xl[1] & xy[, 1] <= xl[2] & xy[, 2] >= yl[1] & xy[, 2] <= yl[2]
  lab <- list()
  for (m in seq(-180, 180 - step, step)) {
    xy <- proj_lines(rep(m, length(la)), la, cfg); lines(xy, col = col, lwd = 0.5)
    k <- which(inside(xy)); if (length(k)) { b <- k[which.min(xy[k, 2])]
      if (xy[b, 2] < yl[1] + 0.05 * diff(yl))
        lab[[length(lab) + 1]] <- list(side = 1, at = xy[b, 1], txt = sprintf("%d\u00b0%s", abs(m), if (m < 0) "W" else if (m > 0) "E" else "")) }
  }
  for (q in seq(-90 + step, 90 - step, step)) {
    xy <- proj_lines(lo, rep(q, length(lo)), cfg); lines(xy, col = col, lwd = 0.5)
    k <- which(inside(xy)); if (length(k)) { b <- k[which.min(xy[k, 1])]
      if (xy[b, 1] < xl[1] + 0.05 * diff(xl))
        lab[[length(lab) + 1]] <- list(side = 2, at = xy[b, 2], txt = sprintf("%d\u00b0%s", abs(q), if (q < 0) "S" else if (q > 0) "N" else "")) }
  }
  # labels in the margin last (drawing in the margin resets the clipping region), then clip again
  for (l in lab) mtext(l$txt, side = l$side, at = l$at, line = 0.1, las = 1, cex = 0.55, col = "#52514e")
  clip(xl[1], xl[2], yl[1], yl[2])
}
#' a label on a soft background box: ONE text object and ONE box per label (easy to edit
#' in the vector PDF; no stacked copies of the same text)
halo_text <- function(x, y, labels, col = "white", halo = "#222222", cex = 0.75, font = 2, adj = 0, r = 0.06) {
  ax <- adj[1]; ay <- if (length(adj) > 1) adj[2] else 0.5
  w <- strwidth(labels, cex = cex, font = font); h <- strheight("Mg", cex = cex, font = font)
  px <- 0.25 * strwidth("M", cex = cex, font = font); py <- 0.22 * h
  rect(x - ax * w - px, y - ay * h - py, x + (1 - ax) * w + px, y + (1 - ay) * h + py,
       col = adjustcolor(halo, if (halo %in% c("white", "#ffffff")) 0.75 else 0.5), border = NA)
  text(x, y, labels, col = col, cex = cex, font = font, adj = c(ax, ay))
}

# ----------------------------------------------------------- data layers ----
#' Natural Earth shapefile (downloaded once)
ne_file <- function(kind = c("borders", "coast"), cfg = SILUCARB) {
  kind <- match.arg(kind); res <- cfg$ne_resolution
  name <- if (kind == "borders") sprintf("ne_%s_admin_0_boundary_lines_land", res) else sprintf("ne_%s_coastline", res)
  f <- file.path(cfg$dir_ne, paste0(name, ".shp"))
  if (!file.exists(f)) {
    dir.create(cfg$dir_ne, showWarnings = FALSE, recursive = TRUE)
    url <- sprintf("https://naciscdn.org/naturalearth/%s/%s/%s.zip", res,
                   if (kind == "borders") "cultural" else "physical", name)
    zf <- file.path(cfg$dir_ne, paste0(name, ".zip"))
    download.file(url, zf, mode = "wb", quiet = TRUE); unzip(zf, exdir = cfg$dir_ne)
  }
  f
}

#' PaleoCoastlines CS outlines (lat/lon line list), cut where they follow the antimeridian
cs_lines <- function(age, cfg = SILUCARB) {
  f <- file.path(cfg$dir_coast, "CS", paste0(age, "Ma_CS_v7.shp"))
  if (!file.exists(f)) return(NULL)
  g <- st_read(f, quiet = TRUE)
  co <- st_coordinates(st_cast(st_geometry(g), "MULTILINESTRING"))
  grp <- interaction(co[, "L1"], co[, "L2"], drop = TRUE)
  out <- list()
  for (gi in levels(grp)) {
    p <- co[grp == gi, c("X", "Y"), drop = FALSE]
    keep <- abs(p[, 1]) < 179.5
    br <- c(FALSE, !keep[-1] | !keep[-nrow(p)] | abs(diff(p[, 1])) > 180)
    seg <- cumsum(br)
    for (s in unique(seg)) { q <- p[seg == s & keep, , drop = FALSE]; if (nrow(q) > 1) out[[length(out) + 1]] <- q }
  }
  out
}

#' PALEOMAP topography moved into another plate model (display product)
grid_to_model <- function(age, z, lat, lon, target, res = 0.25) {
  if (target == "scotese_and_wright2018") return(list(z = z, lat = lat, lon = lon))
  LON <- matrix(lon, nrow(z), ncol(z), byrow = TRUE); LAT <- matrix(lat, nrow(z), ncol(z))
  m <- z >= -3000
  tr <- gp_transfer(LAT[m], LON[m], age, "scotese_and_wright2018", target)
  ok <- is.finite(tr$lat)
  glat <- seq(-90 + res / 2, 90, by = res); glon <- seq(-180 + res / 2, 180, by = res)
  i <- pmin(pmax(floor((tr$lat[ok] + 90) / res) + 1, 1), length(glat))
  j <- pmin(pmax(floor((tr$lon[ok] + 180) / res) + 1, 1), length(glon))
  key <- (j - 1) * length(glat) + i
  acc <- tapply(z[m][ok], key, mean)
  g <- matrix(NA_real_, length(glat), length(glon)); g[as.integer(names(acc))] <- acc
  for (it in 1:5) {                           # fill small gaps with neighbour means
    nas <- is.na(g); if (!any(nas)) break
    s <- matrix(0, nrow(g), ncol(g)); c_ <- s
    for (di in -1:1) for (dj in -1:1) { v <- shift(g, di, dj, fill = NA); ok2 <- !is.na(v); s[ok2] <- s[ok2] + v[ok2]; c_ <- c_ + ok2 }
    g[nas & c_ > 0] <- (s / pmax(c_, 1))[nas & c_ > 0]
  }
  g[is.na(g)] <- -3480
  g <- seafloor_profile(g, glat)
  list(z = smooth_transferred(g, glat), lat = glat, lon = glon)
}
#' Moving plates apart in another plate model leaves land and shelf directly against
#' deep-ocean fill (gaps between plates, crust that does not exist in that model):
#' seafloor_profile() (grid_tools.R) raises the sea floor there to a shelf + slope profile.
#' grid_to_model with a disk cache (re-used until the source grid changes)
cached_grid_to_model <- function(age, dem, model, cfg = SILUCARB) {
  if (model == "scotese_and_wright2018") return(list(z = dem$z, lat = dem$lat, lon = dem$lon))
  dir <- file.path(cfg$dir_results, "cache"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  inf <- file.info(dem$file)
  f <- file.path(dir, sprintf("%s_%s_%s_v2.rds", model, tools::file_path_sans_ext(basename(dem$file)),
                              format(as.numeric(inf$mtime) * 1000 + inf$size, scientific = FALSE)))
  if (file.exists(f)) return(readRDS(f))
  old <- list.files(dir, paste0("^", model, "_", tools::file_path_sans_ext(basename(dem$file)), "_"), full.names = TRUE)
  unlink(old)
  g <- grid_to_model(age, dem$z, dem$lat, dem$lon, model); saveRDS(g, f); g
}
#' smooth the blocky coasts left by moving a regular grid between plate models
smooth_transferred <- function(g, la, sigma_km = 40) {
  s <- smooth_km((g < 0) * 1, la, sigma_km, nbands = 40) - 0.5
  gz <- smooth_km(g, la, sigma_km, nbands = 40)
  ramp <- -80 * s
  out <- ifelse(s > 0, pmin(gz, -1), pmax(gz, 1))
  band <- abs(s) < 0.45
  out <- ifelse(band, ifelse(s > 0, pmin(out, ramp), pmin(pmax(out, 0.5), pmax(ramp, 0.5))), out)
  matrix(out, nrow(g))
}

IND_STYLE <- list(reef = list(pch = 21, bg = "#e87ba4", col = "#3b0d22", cex = 0.85, label = "Reef"),
                  evaporite = list(pch = 22, bg = "#eda100", col = "#3d2800", cex = 0.85, label = "Evaporite"),
                  glacial = list(pch = 23, bg = "white", col = "#0b2a4a", cex = 1.0, label = "Glacial deposit"))
USER_STYLE <- list(pch = 24, bg = "#7fffd4", col = "#003b2e", cex = 0.8)

# ------------------------------------------------------- detail / furniture ----
#' cut the grid to the area of the view and refine it by `detail` (bilinear); "auto"
#' refines regional views so a grid cell is ~1.5 output pixels, up to cfg$map_refine_max
#' (set by map_quality), and global maps by cfg$map_refine_global. A number forces that factor.
refine_for_view <- function(z, lat, lon, cfg, width_px, detail = "auto") {
  keep <- list(z = z, lat = lat, lon = lon, wrap = TRUE)
  dl <- abs(lat[2] - lat[1])
  if (is.null(cfg$view)) {
    k <- if (is.numeric(detail)) round(detail) else (cfg$map_refine_global %||% 1)
    if (k <= 1) return(keep)
    r <- rast(z[nrow(z):1, ], extent = ext(lon[1] - dl / 2, lon[length(lon)] + dl / 2, lat[1] - dl / 2, lat[length(lat)] + dl / 2))
    r <- disagg(r, k, method = "bilinear"); m <- as.matrix(r, wide = TRUE)[nrow(r):1, ]
    return(list(z = m, lat = seq(lat[1] - dl / 2 + dl / (2 * k), by = dl / k, length.out = nrow(m)),
                lon = seq(lon[1] - dl / 2 + dl / (2 * k), by = dl / k, length.out = ncol(m)), wrap = TRUE))
  }
  v <- cfg$view
  gx <- seq(v$xlim[1], v$xlim[2], length.out = 40); gy <- seq(v$ylim[1], v$ylim[2], length.out = 40)
  ll <- sf_project(v$crs, "EPSG:4326", as.matrix(expand.grid(gx, gy)), keep = TRUE, warn = FALSE)
  ll <- ll[is.finite(ll[, 1]) & is.finite(ll[, 2]), , drop = FALSE]
  if (!nrow(ll)) return(keep)
  la_r <- range(ll[, 2]) + c(-2, 2); rel <- ((ll[, 1] - v$lon0 + 180) %% 360) - 180
  lo_r <- v$lon0 + range(rel) + c(-3, 3)
  if (diff(lo_r) > 300 || la_r[1] < -88 || la_r[2] > 88 || lo_r[1] < -180 || lo_r[2] > 180) return(keep)
  ii <- which(lat >= la_r[1] & lat <= la_r[2]); jj <- which(lon >= lo_r[1] & lon <= lo_r[2])
  if (length(ii) < 4 || length(jj) < 4) return(keep)
  km_px <- diff(v$xlim) / 1000 / width_px
  k <- if (is.numeric(detail)) round(detail) else min(cfg$map_refine_max %||% 2, max(1, ceiling(dl * KM_PER_DEG / (1.5 * km_px))))
  zc <- z[ii, jj, drop = FALSE]
  if (k > 1) {
    r <- rast(zc[nrow(zc):1, ], extent = ext(lon[jj[1]] - dl / 2, lon[max(jj)] + dl / 2, lat[ii[1]] - dl / 2, lat[max(ii)] + dl / 2))
    r <- disagg(r, k, method = "bilinear"); zc <- as.matrix(r, wide = TRUE)[nrow(r):1, ]
  }
  list(z = zc, lat = seq(lat[ii[1]] - dl / 2 + dl / (2 * k), by = dl / k, length.out = nrow(zc)),
       lon = seq(lon[jj[1]] - dl / 2 + dl / (2 * k), by = dl / k, length.out = ncol(zc)), wrap = FALSE)
}

.inv <- function(cfg, x, y) sf_project(map_crs(cfg), "EPSG:4326", cbind(x, y), keep = TRUE, warn = FALSE)
.fwd <- function(cfg, lon, lat) sf_project("EPSG:4326", map_crs(cfg), cbind(lon, lat), keep = TRUE, warn = FALSE)
#' point at distance d_km and bearing (deg) from (lat, lon) on the sphere
.dest <- function(lat, lon, d_km, bearing) {
  p1 <- lat * pi / 180; l1 <- lon * pi / 180; b <- bearing * pi / 180; dr <- d_km / 6371
  p2 <- asin(sin(p1) * cos(dr) + cos(p1) * sin(dr) * cos(b))
  l2 <- l1 + atan2(sin(b) * sin(dr) * cos(p1), cos(dr) - sin(p1) * sin(p2))
  c(p2 * 180 / pi, ((l2 * 180 / pi + 180) %% 360) - 180)
}

#' north arrow. Regional maps: placed in the top-right corner, pointing to the (palaeo)
#' north pole at that spot. Global maps: placed outside the ellipse, top right, pointing up
#' (north is straight up along the central meridian of the Mollweide map)
draw_north_arrow <- function(cfg) {
  if (.is_disc(cfg) && abs(cfg$view$lat0) > 89) return(invisible())   # polar maps: north is along every meridian
  xl <- view_xlim(cfg); yl <- view_ylim(cfg); W <- diff(xl); H <- diff(yl)
  if (is.null(cfg$view)) {
    P <- c(xl[1] + 0.965 * W, yl[1] + 0.86 * H); d <- c(0, 1)
  } else {
    P <- c(xl[2] - 0.06 * W, yl[2] - 0.16 * H)
    Q <- if (.is_disc(cfg)) 0.9 * cfg$view$disc * P / sqrt(sum(P^2)) else P   # round maps: north at the rim
    ll <- .inv(cfg, Q[1], Q[2]); if (!is.finite(ll[1, 2]) || abs(ll[1, 2]) > 88) return(invisible())
    q <- .fwd(cfg, ll[1, 1], ll[1, 2] + 0.5); d <- q[1, ] - Q; d <- d / sqrt(sum(d^2))
  }
  L <- 0.08 * H; tip <- P + d * L; nrm <- c(-d[2], d[1]) * 0.28 * L
  polygon(rbind(tip, P + nrm, P - d * 0.15 * L), col = "#1b1b1b", border = "#1b1b1b", lwd = 0.6)
  polygon(rbind(tip, P - nrm, P - d * 0.15 * L), col = "white", border = "#1b1b1b", lwd = 0.6)
  halo_text(tip[1] + d[1] * 0.2 * L, tip[2] + d[2] * 0.2 * L, "N", adj = c(0.5, 0), cex = 0.8,
            col = "#1b1b1b", halo = "white")
}

#' scale bar. Regional maps: bottom right, exact at its position. Global maps: outside the
#' ellipse, bottom right; Mollweide changes scale across the map, so the bar gives the
#' east-west scale at the map centre (equator, central meridian) and says so
draw_scale_bar <- function(cfg) {
  xl <- view_xlim(cfg); yl <- view_ylim(cfg); W <- diff(xl); H <- diff(yl)
  if (is.null(cfg$view)) {
    P <- c(xl[1] + 0.80 * W, yl[1] + 0.035 * H); ll <- cbind(cfg$central_lon, 0); P0 <- .fwd(cfg, ll[1, 1], 0)[1, ]
    frac <- 0.15
  } else if (.is_disc(cfg)) {                     # round maps: scale at the centre, bar in the corner
    P <- c(xl[2] - 0.27 * W, yl[1] + 0.03 * H); P0 <- c(0, 0)
    ll <- cbind(cfg$view$lon0, if (abs(cfg$view$lat0) > 89) sign(cfg$view$lat0) * 89.9 else cfg$view$lat0); frac <- 0.15
  } else {
    P <- c(xl[2] - 0.30 * W, yl[1] + 0.06 * H); ll <- .inv(cfg, P[1], P[2]); P0 <- P
    if (!is.finite(ll[1, 2])) return(invisible()); frac <- 0.2
  }
  e <- .dest(ll[1, 2], ll[1, 1], 100, 90); q <- .fwd(cfg, e[2], e[1])
  per_km <- sqrt(sum((q[1, ] - P0)^2)) / 100                 # map units per km there
  raw <- frac * W / per_km; p10 <- 10^floor(log10(raw)); cand <- c(1, 2, 5, 10) * p10
  nice <- cand[which.min(abs(log(cand / raw)))]
  len <- nice * per_km; h <- 0.012 * H
  xs <- P[1] + len * (0:4) / 4
  for (k in 1:4) rect(xs[k], P[2], xs[k + 1], P[2] + h, col = if (k %% 2) "#1b1b1b" else "white", border = "#1b1b1b", lwd = 0.6)
  lab <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  halo_text(xs[c(1, 3, 5)], P[2] + 2.6 * h, c("0", lab(nice / 2), paste(lab(nice), "km")), adj = 0.5, cex = 0.6,
            col = "#1b1b1b", halo = "white", font = 1)
  if (is.null(cfg$view) || .is_disc(cfg))
    halo_text(P[1] + len / 2, P[2] - 1.8 * h, if (is.null(cfg$view)) "scale at the equator, map centre" else "scale at the map centre", adj = 0.5, cex = 0.5,
              col = "#1b1b1b", halo = "white", font = 1)
}

# ------------------------------------------------------------- one panel ----
#' draw one map panel into the current plot region (already sized)
#' the colours of a map layer on a grid (zs: heights for shading, z: heights; lat, lon: grid),
#' with sea ice and ice sheets; also returns the sea-ice fields and the ice for the outlines
.layer_rgb <- function(climate, zs, z, lat, lon) {
  hs <- list(z = zs)
  rgb_img <- if (identical(climate$layer, "precipitation"))
    render_precip_rgb(hs$z, lat, precip_on_grid(climate$f, lat, lon)) else
    if (identical(climate$layer, "temperature")) render_temp_rgb(hs$z, lat, .temp_on_map(climate, z, lat, lon)) else
    if (isTRUE(climate$layer %in% names(SCALAR_LAYERS)) && !is.null(climate$scale)) {
      L <- climate$scale; vv <- .scalar_on_grid(climate$f[[L$field]], lat, lon) * .scalar_factor(L, climate$season)
      render_scalar_rgb(render_rgb(hs$z, lat), hs$z, lat, vv, L) } else
    if (isTRUE(climate$layer %in% c("precipitation_change", "temperature_change")))
      render_change_rgb(hs$z, lat, change_on_grid(climate$f$change[[if (climate$layer == "temperature_change") "t0" else "precip"]], lat, lon),
                        climate$layer, climate$lim %||% CHANGE_LIMITS[[climate$layer]]) else render_rgb(hs$z, lat)
  if (!is.null(climate$sea_ice)) {                   # permanent / seasonal sea ice
    si_min <- field_on_grid(climate$sea_ice$min, lat, lon); si_max <- field_on_grid(climate$sea_ice$max, lat, lon)
    rgb_img <- overlay_sea_ice_rgb(rgb_img, hs$z, si_min, si_max)
  }
  ic <- NULL
  if (!is.null(climate$ice)) {                       # ice sheets of an ice run
    ic <- ice_on_latlon(climate$ice$state, lat, lon)
    # ice painted white on the topography / thickness map; on climate layers only its outline is drawn
    if (!is.null(ic) && (is.null(climate$layer) || climate$layer %in% c("topography", "ice_thickness", "ice_smb")))
      rgb_img <- overlay_ice_rgb(rgb_img, hs$z, lat, ic, climate$layer)
  }
  list(rgb = rgb_img, si_min = if (exists("si_min", inherits = FALSE)) si_min, si_max = if (exists("si_max", inherits = FALSE)) si_max, ic = ic)
}

draw_map_panel <- function(age, model, z, lat, lon, cfg = SILUCARB, corrected = TRUE, borders = TRUE,
                           modern_coast = FALSE, indicators = TRUE, cs = NULL, recs = NULL,
                           user_points = TRUE, sites = NULL, legend = TRUE, width_px = 2300, equator = TRUE,
                           detail = "auto", north_arrow = FALSE, scale_bar = FALSE, ice_front = TRUE,
                           climate = NULL, rgb_cache = NULL) {
  ind_kinds <- if (isTRUE(indicators)) names(IND_STYLE) else if (isFALSE(indicators)) character(0) else
    intersect(indicators, names(IND_STYLE))
  xl <- view_xlim(cfg); yl <- view_ylim(cfg)
  plot.window(xlim = xl, ylim = yl, asp = 1, xaxs = "i", yaxs = "i")
  if (!is.null(cfg$view)) clip(xl[1], xl[2], yl[1], yl[2])      # nothing outside the regional frame
  # higher detail: the grid is cut to the view and refined (bilinear) before shading, so
  # zoomed-in maps get smooth relief and coastlines instead of enlarged 0.2-degree cells
  zs <- z                                         # shading grid (land relief optionally smoothed)
  if ((cfg$map_relief_smooth_km %||% 0) > 0 && nrow(z) == length(lat) && abs(diff(range(lon))) > 300) {
    sm <- smooth_km(pmax(z, 0), lat, cfg$map_relief_smooth_km)
    zs <- ifelse(z > 0, pmax(sm, pmin(z, 1)), z)  # land stays land, sea and coastline untouched
  }
  z0 <- z; lat0 <- lat; lon0 <- lon                  # the global grid (for the climate arrows)
  if (isTRUE(cfg$view$fast_eqc)) {                 # Figure-tab preview (plain lon/lat view): shade the whole grid once,
    hr <- list(z = z, lat = lat, lon = lon, wrap = TRUE)   # .project_rgb_eqc() then cuts out the window
    hs <- hr
  } else {
    hr <- refine_for_view(z, lat, lon, cfg, width_px, detail)
    hs <- if (identical(zs, z)) hr else refine_for_view(zs, lat, lon, cfg, width_px, detail)
  }
  z <- hr$z; lat <- hr$lat; lon <- hr$lon
  ckey <- sprintf("rgb_%d_%d_%.5f_%.5f", nrow(z), ncol(z), lat[1], lon[1])
  if (!is.null(rgb_cache) && exists(ckey, envir = rgb_cache, inherits = FALSE)) {
    lr <- get(ckey, envir = rgb_cache)              # same map and layer as the last draw: re-use the colours and shading
  } else {
    lr <- .layer_rgb(climate, hs$z, z, lat, lon)
    if (!is.null(rgb_cache)) {
      assign(ckey, lr, envir = rgb_cache)
    }
  }
  rgb_img <- lr$rgb; si_min <- lr$si_min; si_max <- lr$si_max; ic <- lr$ic
  img <- project_rgb(rgb_img, lat, lon, width_px, cfg, wrap = hr$wrap)
  rasterImage(img, xl[1], yl[1], xl[2], yl[2], interpolate = TRUE)
  if (!is.null(cs)) draw_ll_list(lapply(cs, function(p) data.frame(lon = p[, 1], lat = p[, 2])), cfg,
                                  col = adjustcolor("#5a1e10", 0.85), lwd = 0.6, lty = "22")
  if (corrected) {                                  # 0 m line of the displayed grid
    if (hr$wrap) { zz <- cbind(z, z[, 1]); ll <- c(lon, lon[1] + 360) } else { zz <- z; ll <- lon }
    cl <- contourLines(ll, lat, t(zz), levels = 0)
    draw_ll_list(lapply(cl, function(q) data.frame(lon = ((q$x + 180) %% 360) - 180, lat = q$y)), cfg,
                 col = "#2b2b2b", lwd = 0.6)
  }
  draw_set <- function(L, col, a_old, a_young) {
    y <- vapply(L, function(p) isTRUE(attr(p, "young")), TRUE)
    draw_ll_list(L[!y], cfg, col = adjustcolor(col, a_old), lwd = 0.65)
    draw_ll_list(L[y], cfg, col = adjustcolor(col, a_young), lwd = 0.65, lty = 3)
  }
  if (modern_coast || borders) draw_set(gp_modern_lines(ne_file("coast", cfg), model, age), "white", 0.95, 0.75)
  if (borders) draw_set(gp_modern_lines(ne_file("borders", cfg), model, age), "#1f1f1f", 0.85, 0.65)
  if (!is.null(climate$f) && isTRUE(climate$flag)) draw_climate_flag(climate$f$flag, cfg)
  if (!is.null(climate$sea_ice)) draw_sea_ice_lines(si_min, si_max, z, lat, lon, cfg, wrap = hr$wrap)
  if (!is.null(climate$rivers)) draw_rivers(climate$rivers, cfg)
  if (!is.null(ic) && !identical(climate$layer, "ice_smb")) draw_ice_lines(ic, lat, lon, cfg, wrap = hr$wrap)
  clim_sp <- if (!is.null(climate$f)) {
    fa <- climate$f
    if (isTRUE(climate$layer %in% c("precipitation_change", "temperature_change"))) {   # arrows: wind change
      fa$wind_u <- fa$change$wind_u; fa$wind_v <- fa$change$wind_v; fa$wind_ref <- 2 }
    draw_climate_arrows(fa, cfg, z0, lat0, lon0, isTRUE(climate$wind), isTRUE(climate$currents))
  }
  draw_graticule(cfg)
  if (isTRUE(equator)) {                         # palaeo-equator of this map
    lo <- seq(-180, 180, 0.5)
    draw_ll(lo, rep(0, length(lo)), col = "#ffd54a", lwd = 1.6, lty = "44", lend = "butt", cfg = cfg)
  }
  map_frame(cfg, border = "#3b3a37", lwd = 1.1)
  # localities (points/localities.csv); sites = FALSE hides them, or give names
  loc <- read_localities(cfg)
  st_ <- if (is.null(sites) || isTRUE(sites)) loc[loc$show, ] else if (isFALSE(sites)) loc[0, ] else loc[loc$name %in% sites, ]
  if (nrow(st_)) {
  s <- gp_reconstruct(st_$lat, st_$lon, rep(age, nrow(st_)), model)
  okp <- is.finite(s$plat); s <- s[okp, , drop = FALSE]; st_ <- st_[okp, , drop = FALSE]
  }
  if (nrow(st_)) {
  xy <- sf_project("EPSG:4326", map_crs(cfg), cbind(s$plon, s$plat), keep = TRUE, warn = FALSE)
  vis <- is.finite(xy[, 1]) & is.finite(xy[, 2]) & !.far(s$plon, s$plat, cfg)   # e.g. far side of a globe
  xy <- xy[vis, , drop = FALSE]; s <- s[vis, , drop = FALSE]; st_ <- st_[vis, , drop = FALSE]
  }
  if (nrow(st_)) {
  points(xy, pch = 21, bg = "#e34948", col = "white", cex = 1.2, lwd = 1.5)
  rel <- ((s$plon - cfg$central_lon + 180) %% 360) - 180
  right <- xy[, 1] > xl[1] + 0.78 * diff(xl)       # labels near the right edge go to the left of the point
  off <- 0.012 * diff(xl)
  # labels of nearby localities are moved apart vertically so they do not overlap
  ly <- xy[, 2]; dy <- 0.028 * diff(yl); o <- order(-ly)
  for (a in seq_along(o)[-1]) for (b in seq_len(a - 1)) {
    i <- o[a]; j <- o[b]
    if (right[i] == right[j] && abs(xy[i, 1] - xy[j, 1]) < 0.12 * diff(xl) && abs(ly[i] - ly[j]) < dy) ly[i] <- ly[j] - dy
  }
  moved <- abs(ly - xy[, 2]) > 0
  if (any(moved)) segments(xy[moved, 1], xy[moved, 2], xy[moved, 1] + ifelse(right[moved], -off, off) * 0.8, ly[moved],
                           col = "#222222", lwd = 0.6)
  labs <- sprintf("%s %.0f\u00b0%s", st_$name, abs(s$plat), ifelse(s$plat < 0, "S", "N"))
  for (i in seq_along(labs))                         # one call per label: text() uses only one adj
    halo_text(xy[i, 1] + if (right[i]) -off else off, ly[i], labs[i], adj = if (right[i]) 1 else 0)
  }
  handles <- NULL
  if (isTRUE(equator)) handles <- data.frame(label = "Equator", pch = NA, bg = NA, col = "#e0b100", lty = 3, lwd = 1.6)
  if (length(ind_kinds) || isTRUE(ice_front)) {
    d <- records_for_age(age, model, purpose = "indicator", recs = recs, cfg = cfg)
    d <- d[!duplicated(paste(d$indicator, round(d$lat, 4), round(d$lon, 4))), ]   # one symbol per site
    for (k in ind_kinds) {
      e <- d[d$indicator == k, ]; st <- IND_STYLE[[k]]
      if (nrow(e)) points(sf_project("EPSG:4326", map_crs(cfg), cbind(e$plon, e$plat), keep = TRUE, warn = FALSE),
                          pch = st$pch, bg = st$bg, col = st$col, cex = st$cex, lwd = 0.6)
      if (nrow(e)) handles <- rbind(handles, data.frame(label = sprintf("%s (%d)", st$label, nrow(e)), pch = st$pch,
                                           bg = st$bg, col = st$col, lty = NA, lwd = 0.6))
    }
    g <- d[d$indicator == "glacial", ]
    if (nrow(g) && isTRUE(ice_front)) {                 # ice front: own toggle (also without glacial symbols)
      front <- g$plat[which.min(abs(g$plat))]; lo <- seq(-180, 180, 1)
      draw_ll(lo, rep(front, length(lo)), col = "#6fd0ff", lwd = 2, lty = "62", cfg = cfg)   # one line object
      handles <- rbind(handles, data.frame(label = sprintf("Ice front from glacial deposits (%.0f\u00b0%s)",
                       abs(front), ifelse(front < 0, "S", "N")), pch = NA, bg = NA, col = "#6fd0ff", lty = 2, lwd = 2))
    }
  }
  if (isTRUE(user_points)) {                           # your own points (any setting)
    u <- recs[recs$source == "user", ]
    if (!is.null(u) && nrow(u)) {
      uu <- records_for_age(age, model, purpose = "dem", recs = u, cfg = cfg)
      if (nrow(uu)) points(sf_project("EPSG:4326", map_crs(cfg), cbind(uu$plon, uu$plat), keep = TRUE, warn = FALSE),
                           pch = USER_STYLE$pch, bg = USER_STYLE$bg, col = USER_STYLE$col, cex = USER_STYLE$cex, lwd = 0.7)
      if (nrow(uu)) handles <- rbind(handles, data.frame(label = sprintf("Your points (%d)", nrow(uu)), pch = USER_STYLE$pch,
                                           bg = USER_STYLE$bg, col = USER_STYLE$col, lty = NA, lwd = 0.7))
    }
  }
  if (!is.null(clim_sp)) draw_climate_key(cfg, clim_sp, isTRUE(climate$wind), isTRUE(climate$currents))
  if (isTRUE(north_arrow)) draw_north_arrow(cfg)
  if (isTRUE(scale_bar)) draw_scale_bar(cfg)
  if (legend && !is.null(handles))
    legend(if (is.null(cfg$view)) "bottomleft" else xl[1], if (is.null(cfg$view)) NULL else yl[1],
           xjust = 0, yjust = 0, legend = handles$label, pch = handles$pch, pt.bg = handles$bg, col = ifelse(handles$col == "white", "#0b2a4a", handles$col),
           lty = handles$lty, lwd = handles$lwd, bg = adjustcolor("white", 0.9), box.col = "#c9c8c2", cex = 0.62,
           pt.cex = 0.9, inset = c(0, 0.0))
  invisible(handles)
}

#' elevation / depth colour bar; side: 1 = horizontal, labels below; 3 = horizontal, labels above;
#' 2 = vertical, labels on the left; 4 = vertical, labels on the right
colour_bar <- function(side = 1) {
  zb <- seq(-4000, 3000, by = 25); cols <- colourise(matrix(zb[-1] - 12.5, 1))   # vector colour boxes
  cc <- rgb(cols[1, , 1], cols[1, , 2], cols[1, , 3]); n <- length(zb)
  ticks <- c(-4000, -3000, -1400, -200, 0, 400, 1300, 2500, 3000)
  lab <- "Elevation / depth (m, relative to the map's sea level)"
  if (side %in% c(1, 3)) {
    plot.window(xlim = c(-4000, 3000), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(zb[-n], 0, zb[-1], 1, col = cc, border = NA)
    axis(side, at = ticks, labels = format(ticks, big.mark = ","), cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    abline(v = 0, col = "#222222", lwd = 1)
    mtext(lab, side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = c(-4000, 3000), xaxs = "i", yaxs = "i")
    rect(0, zb[-n], 1, zb[-1], col = cc, border = NA)
    axis(side, at = ticks, labels = format(ticks, big.mark = ",", trim = TRUE), las = 1, cex.axis = 0.6, lwd = 0,
         lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    abline(h = 0, col = "#222222", lwd = 1)
    mtext("Elevation / depth (m)", side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}

#' page layout for the colour bar position ("bottom", "top", "left", "right"): M is the map
#' area the page uses with the bar at the bottom, B that bottom bar. Returns the new map area,
#' a function that moves any map figure region from M into it, the bar's figure region and side.
.bar_layout <- function(pos, M, B) {
  pos <- match.arg(pos %||% "bottom", c("bottom", "top", "left", "right"))
  M2 <- switch(pos, bottom = M, top = c(M[1], M[2], 0.04, M[4] - 0.09),
               right = c(M[1], M[2] - 0.1, 0.04, M[4]), left = c(M[1] + 0.1, M[2], 0.04, M[4]))
  h <- M2[4] - M2[3]
  bar <- switch(pos, bottom = B, top = c(0.25, 0.75, M[4] - 0.08, M[4] - 0.055),
                right = c(M2[2] + 0.02, M2[2] + 0.038, M2[3] + 0.2 * h, M2[4] - 0.2 * h),
                left = c(M2[1] - 0.038, M2[1] - 0.02, M2[3] + 0.2 * h, M2[4] - 0.2 * h))
  tf <- function(f) c(M2[1] + (f[1:2] - M[1]) / (M[2] - M[1]) * (M2[2] - M2[1]),
                      M2[3] + (f[3:4] - M[3]) / (M[4] - M[3]) * (M2[4] - M2[3]))
  list(map = M2, tf = tf, bar = bar, side = c(bottom = 1, top = 3, left = 2, right = 4)[[pos]])
}

# -------------------------------------------------------------- map pages ----
#' map centre for one age and plate model.
#' center: NULL (cfg$central_lon), a palaeo-longitude (one number), a site name from
#' points/localities.csv (e.g. "Huy": the map follows that locality), or c(lat, lon) present-day
#' coordinates of any locality (followed through time)
map_center <- function(center, age, model, cfg = SILUCARB) {
  if (is.null(center)) return(cfg$central_lon)
  if (is.numeric(center) && length(center) == 1) return(center)
  if (is.character(center)) {
    loc <- read_localities(cfg); k <- match(center, loc$name)
    if (is.na(k)) stop("unknown locality '", center, "'; localities: ", paste(loc$name, collapse = ", "))
    center <- c(loc$lat[k], loc$lon[k])
  }
  p <- gp_reconstruct(center[1], center[2], age, model)
  if (!is.finite(p$plon)) { warning("centre point not on a plate at ", age, " Ma; using default"); return(cfg$central_lon) }
  round(p$plon, 1)
}
.region_tag <- function(region) {
  if (is.null(region)) return("")
  if (is.numeric(region)) return(sprintf("_region%+.0f%+.0f%+.0f%+.0f", region[1], region[2], region[3], region[4]))
  cc <- if (!is.null(region$center_palaeo)) sprintf("p%.1fN%.1fE", region$center_palaeo[1], region$center_palaeo[2]) else
        if (is.character(region$center)) region$center else sprintf("%gN%gE", region$center[1], region$center[2])
  sprintf("_region-%s-%gx%gkm", cc, round(region$width_km %||% 4000),
          round(region$height_km %||% (0.6 * (region$width_km %||% 4000))))
}
.center_tag <- function(center) {
  if (is.null(center)) return("")
  if (is.character(center)) return(paste0("_centre-", center))
  if (length(center) == 1) return(sprintf("_centre%g", center))
  sprintf("_centre%gN%gE", center[1], center[2])
}

#' render the per-age maps (PNG + one PDF)
#' indicators:  TRUE / FALSE, or a selection, e.g. c("reef", "glacial")
#' user_points: show your own points (TRUE / FALSE)
#' sites:       NULL = localities marked show in points/localities.csv; FALSE = none;
#'              or locality names, e.g. c("Huy", "Prague")
#' shoreline:   draw the original PaleoCoastlines shoreline (dashed; Kocsis & Scotese 2021)
#' equator:     draw the palaeo-equator (yellow dashed line)
#' ice_front:   draw the ice front (latitude of the glacial deposit nearest the equator)
#' center:      NULL, a longitude, a site name ("Huy") or c(lat, lon) - see map_center()
#' ages:        any PaleoDEM ages, 0-540 Ma (available_ages() lists them all)
#' north_arrow, scale_bar: TRUE / FALSE
#' projection:  "mollweide" (default, map_projection in R/config.R), "robinson", "equal_earth",
#'              "winkel_tripel", "hammer", "eckert4", "plate_carree", or a round map:
#'              "orthographic" (globe centred on `center`), "south_polar", "north_polar"
#'              (polar_edge_lat: latitude of the polar map's edge, 0 = equator, e.g. 30 for a
#'              south polar map reaching to 30N). A `region` map always uses its own projection.
#' quality:     "draft", "standard" (default, from map_quality in R/config.R) or "high"
#' dpi:         overrides the quality's resolution; detail: a number forces the refinement factor
#' overview:    with a region: also draw the whole-earth map (in `projection`) with a box around
#'              the region, joined to the zoomed-in map by two lines
#' overview_corner: where the zoomed-in map sits, partly over the whole-earth map: "auto" (the
#'              corner away from the box), "bottomright", "bottomleft", "topright", "topleft"
#' overview_ratio: size of the zoom-in relative to the whole-earth map (1 = same size,
#'              2 = twice as large, 0.7 = smaller); overview_overlap: how much the two overlap
#'              (fraction of the smaller one; 0 = corners touch, negative = a gap)
#' colour_bar_pos: where the elevation / depth colour bar goes: "bottom" (default, colour_bar_pos
#'              in R/config.R), "top", "left" or "right"
#' layer:       "topography" (default), "precipitation", "temperature" (near-surface air, deg C),
#'              or for an experiment "precipitation_change" / "temperature_change"
#'              (experiment - reference, or experiment - `compare`)
#' climate_method: "local" (default; HadCM3L corrected locally for your geography), "model"
#'              (HadCM3L as run) or "rerun" (+ SPEEDY response; needs run_speedy() first)
#' climate_flag: outline areas where the corrected geography differs much from the model's
#' experiment:  climate_experiment(...) (orbit, CO2 in PAL, Sun age); needs run_speedy() first.
#'              orbit, co2 are a short form: an orbit preset and a CO2 factor relative to HadCM3L
#' compare:     another climate_experiment(): change layers then show experiment - compare
#' map:         a map loaded with import_map() (or climate_map()) instead of the editor's map of
#'              the age; ages is then the map's age
#' wind, currents: TRUE draws 10 m wind / surface-current arrows
#' season:      "ann" (annual mean), "djf" (Dec-Feb) or "jja" (Jun-Aug)
#'              climate: HadCM3L runs of Valdes et al. (2021), see R/climate.R; PALEOMAP frame only
#' pdf:         "vector" (default): lines, points, text, legend and colour bar are vector
#'              objects you can edit (e.g. in Illustrator / Inkscape); only the shaded
#'              topography is an embedded image. "raster": pages are images (as before).
#' short file-name tag of an experiment
.exp_file_tag <- function(x) {
  if (.is_reference(x)) return("reference")
  sprintf("e%.3f_o%.1f_p%.0f_%s%s%s_%s", x$ecc, x$obl, x$peri,
          if (!is.null(x$co2_pal)) sprintf("co2-%gPAL", x$co2_pal) else sprintf("co2x%g", x$co2 %||% 1),
          if (!is.null(x$sun_age)) sprintf("_sun%gMa", x$sun_age) else "",
          if (abs((x$gradient %||% 1) - 1) > 1e-3) sprintf("_grad%g", x$gradient) else "", x$length %||% "standard")
}

#' the climate of a figure (layer, fields, sea ice, ice sheets) for one age: used by make_maps()
#' and by the editor's interactive Figure map
.figure_climate <- function(age, model, layer, season, climate_method, experiment, compare, map, corrected, cfg,
                            climate_flag, wind, currents, rivers, sea_ice, ice, sn, sl_ice, o2_pal, lapse_rate, use_clim,
                            diff_with = NULL) {
  clim <- NULL
  if (use_clim) {
    fcl <- if (corrected) climate_for_map(map %||% age, season, climate_method, cfg, climate_flag,
                                          experiment, compare)
           else { f0 <- climate_fields(age, season, cfg); f0$temp_href <- orig_href(age, cfg); f0 }
    lay <- layer
    if (lay %in% c("precipitation_change", "temperature_change") && is.null(fcl$change)) {
      message("No finished orbit / CO2 / Sun experiment for ", age, " Ma: drawing ",
              if (lay == "temperature_change") "temperature" else "precipitation", " instead of the change.")
      lay <- if (lay == "temperature_change") "temperature" else "precipitation"
    }
    clim <- list(layer = lay, wind = isTRUE(wind), currents = isTRUE(currents), season = season, f = fcl,
                 flag = isTRUE(climate_flag) && corrected)
    if (lay == "o2_surface") {
      fcl$o2_surface <- .fld(o2_saturation(.to_1deg(fcl$sst_ocean), .to_1deg(fcl$salinity), o2_pal)); clim$f <- fcl
      clim$extra <- sprintf("O2 at saturation (Garcia & Gordon 1992) for the sea-surface temperature and salinity, atmospheric O2 %g PAL; deep-water O2 is not modelled", o2_pal)
    }
    if (lay == "surface_density") {
      fcl$surface_density <- .fld(rho_sw(pmax(.to_1deg(fcl$sst_ocean), -1.8), .to_1deg(fcl$salinity)) - 1000); clim$f <- fcl
      clim$extra <- "surface density (simplified TEOS-10) of the sea-surface temperature (incl. SPEEDY change) and the HadCM3L salinity; the densest winter water can sink and ventilate the deep sea (DJF for the north, JJA for the south)"
    }
    if (lay == "w_ocean") clim$extra <- paste0("vertical velocity: HadCM3L", if (isTRUE(fcl$ekman)) " + change of Ekman pumping by the SPEEDY wind stress" else "")
    if (lay %in% c("salinity", "mld")) clim$extra <- "ocean of the HadCM3L run (not changed by the map or the experiment)"
    if (lay %in% c("pme", "runoff")) clim$extra <- "runoff on land where the model had none, and its change, from P - E (long-term water balance)"
    if (lay %in% c("freshwater", "river_input", "salinity_source") || isTRUE(rivers)) {
      fcl <- tryCatch(freshwater_fields(fcl, map$file %||% find_dem(if (corrected) cfg$dir_corrected else cfg$dir_dem, age)),
                      error = function(e) { message("rivers: ", conditionMessage(e)); fcl })
      clim$f <- fcl; clim$rivers <- if (isTRUE(rivers)) fcl$rivers
      if (lay %in% c("freshwater", "river_input")) clim$extra <- "runoff routed downhill over the map to the coast and spread ~300 km from the mouths; P - E over the sea"
      if (lay == "salinity_source") clim$extra <- sprintf(paste(
        "HadCM3L surface salinity compared with its ocean mean (\u00b1%g psu) and the local freshwater flux (P - E + rivers, \u00b1%g mm/yr):",
        "where the local flux cannot explain it, the water was brought by currents, mixing, upwelling or sea ice"), SALINITY_DS, SALINITY_FLUX)
    }
    if (lay == "temperature") {
      clim$lapse <- lapse_rate
      clim$extra <- if (isTRUE(lapse_rate > 0)) sprintf("temperature brought to the 0.2\u00b0 map heights%s with a lapse rate of %g K/km",
                                                        if (!is.null(ice)) " and the ice surface" else "", lapse_rate)
                    else "temperature at the climate model's smooth terrain (no lapse-rate correction)"
    }
    if (lay %in% names(SCALAR_LAYERS)) clim$scale <- scalar_scale(lay, fcl[[SCALAR_LAYERS[[lay]]$field]], season)
    if (lay %in% c("precipitation_change", "temperature_change"))
      clim$lim <- change_limit(fcl$change[[if (lay == "temperature_change") "t0" else "precip"]], lay)
    if (!is.null(diff_with)) clim <- .difference_climate(clim, diff_with, age, season, climate_method, corrected, cfg, o2_pal, map)
  }
  if (isTRUE(sea_ice) && model == "scotese_and_wright2018") {   # permanent / seasonal sea ice
    ist <- if (!is.null(ice)) ice_snapshot(ice$id, ice$t_kyr, cfg)
    si <- tryCatch(sea_ice_for_map(map %||% age, climate_method, experiment, cfg, ist),
                   error = function(e) { message("sea ice: ", conditionMessage(e)); NULL })
    if (!is.null(si)) { if (is.null(clim)) clim <- list(layer = layer, wind = FALSE, currents = FALSE, season = season); clim$sea_ice <- si }
  }
  if (!is.null(ice)) {                               # ice sheets of an ice run (a snapshot)
    if (!is.null(sn)) {
      if (is.null(clim)) clim <- list(layer = layer, wind = FALSE, currents = FALSE, season = season)
      else if (layer %in% c("ice_thickness", "ice_smb")) clim$layer <- layer
      clim$ice <- list(state = sn, label = ice$label, coupled = isTRUE(ice$coupled), lower_sea = sl_ice < 0)
    }
  }
  clim
}

make_maps <- function(ages = SILUCARB$ages, model = "scotese_and_wright2018", corrected = TRUE,
                      borders = TRUE, modern_coast = FALSE, indicators = TRUE, user_points = TRUE,
                      sites = NULL, center = NULL, shoreline = FALSE, region = NULL, overview = FALSE,
                      overview_corner = "auto", overview_ratio = 1.25, overview_overlap = 0.15,
                      overview_show = character(0),
                      equator = TRUE, ice_front = TRUE, colour_bar_pos = cfg$colour_bar_pos,
                      layer = c("topography", "precipitation", "temperature", "precipitation_change", "temperature_change",
                                "ice_thickness", "ice_smb", "sea_ice", "pme", "runoff", "freshwater", "river_input", "salinity_source",
                                "salinity", "w_ocean", "mld", "o2_surface", "surface_density"),
                      wind = FALSE, currents = FALSE, season = "ann",
                      climate_method = c("local", "model", "rerun"), climate_flag = TRUE, orbit = "present", co2 = 1,
                      sensitivity = CLIMATE_SENSITIVITY, experiment = NULL, compare = NULL, map = NULL, ice = NULL, sea_ice = FALSE, o2_pal = 1, lapse_rate = 6.5,
                      rivers = FALSE, diff_with = NULL, seasons4 = FALSE,
                      north_arrow = FALSE, scale_bar = FALSE, projection = cfg$map_projection, polar_edge_lat = cfg$polar_edge_lat,
                      quality = cfg$map_quality, dpi = NULL, detail = "auto",
                      comparison = is.null(region), comparison_age = if (425 %in% ages) 425 else ages[1],
                      pdf = c("vector", "raster"), cfg = SILUCARB) {
  pdf <- match.arg(pdf); layer <- match.arg(layer); season <- match.arg(season, names(CLIMATE_SEASONS))
  climate_method <- match.arg(climate_method)
  experiment <- experiment %||% climate_experiment(orbit, co2, sensitivity = sensitivity)
  if (!is.null(map)) { map <- climate_map(map, cfg = cfg); ages <- map$age; corrected <- TRUE }
  if (layer %in% c("ice_thickness", "ice_smb") && is.null(ice)) {
    message("Ice layers need an ice run (ice = list(id = ..., t_kyr = ...)): drawing the topography.")
    layer <- "topography"
  }
  if (layer == "sea_ice") sea_ice <- TRUE
  use_clim <- !layer %in% c("topography", "ice_thickness", "ice_smb", "sea_ice") || isTRUE(wind) || isTRUE(currents) || isTRUE(rivers)
  if (use_clim && model != "scotese_and_wright2018") {
    message("Precipitation, wind and currents come from climate-model runs on the PALEOMAP maps; ",
            "they are only drawn for model = \"scotese_and_wright2018\".")
    use_clim <- FALSE; layer <- "topography"; wind <- FALSE; currents <- FALSE
  }
  cfg$map_projection <- projection %||% "mollweide"; cfg$polar_edge_lat <- polar_edge_lat %||% 0
  invisible(.proj_name(cfg))                             # checks the name
  q <- cfg$MAP_QUALITY[[quality]]
  if (is.null(q)) stop("quality must be one of: ", paste(names(cfg$MAP_QUALITY), collapse = ", "))
  cfg$map_dpi <- if (is.null(dpi)) q$dpi else dpi
  cfg$map_refine_max <- q$refine_max; cfg$map_refine_global <- q$refine_global
  cfg$map_relief_smooth_km <- q$relief_smooth_km %||% 0
  dir.create(cfg$dir_maps, showWarnings = FALSE, recursive = TRUE)
  recs <- all_records(cfg)
  has_ind <- !isFALSE(indicators) && length(indicators)
  suffix <- paste0(if (corrected) "_corrected" else "", if (model != "scotese_and_wright2018") paste0("_", model) else "",
                   if (has_ind) "_indicators" else "", if (isTRUE(user_points)) "_userpoints" else "",
                   if (borders) "_borders" else "", if (modern_coast && !borders) "_moderncoast" else "",
                   if (shoreline) "_shoreline" else "", .region_tag(region),
                   .center_tag(center), if ((is.null(region) || overview) && cfg$map_projection != "mollweide") paste0("_", cfg$map_projection) else "",
                   if (!is.null(region) && overview) "_overview" else "",
                   if (use_clim || isTRUE(sea_ice)) paste0(if (layer == "precipitation") "_precip" else "", if (isTRUE(wind)) "_wind" else "",
                                        if (isTRUE(currents)) "_currents" else "", if (season != "ann") paste0("_", season) else "",
                                        if (layer %in% c("temperature", names(SCALAR_LAYERS))) paste0("_", layer) else "",
                                        if (layer %in% c("precipitation_change", "temperature_change")) paste0("_", layer) else "",
                                        if (climate_method != "local") paste0("_", climate_method) else "",
                                        if (!.is_reference(experiment)) paste0("_", .exp_file_tag(experiment)) else "",
                                        if (!is.null(compare)) paste0("_vs_", .exp_file_tag(compare)) else "") else "",
                   if (!is.null(map)) paste0("_", gsub("[^A-Za-z0-9-]+", "-", map$name)) else "",
                   if (isTRUE(sea_ice)) "_seaice" else "", if (isTRUE(rivers)) "_rivers" else "",
                   if (!is.null(diff_with)) paste0("_minus_", gsub("[^A-Za-z0-9.-]+", "-", substr(diff_with$tag %||% diff_with$label %||% "other", 1, 60))) else "",
                   if (isTRUE(seasons4)) "_4seasons" else "",
                   if (layer == "temperature" && !isTRUE(lapse_rate > 0)) "_nolapse" else
                     if (layer == "temperature" && lapse_rate != 6.5) sprintf("_lapse%g", lapse_rate) else "",
                   if (!is.null(ice)) sprintf("_%s_%s%s", if (layer == "ice_smb") "icesmb" else "ice",
                                              substr(ice$id, nchar(ice$id) - 10, nchar(ice$id)),
                                              if (isTRUE(is.finite(ice$t_kyr))) sprintf("_%gkyr", ice$t_kyr) else "") else "")
  tag <- if (model == "scotese_and_wright2018") "PALEOMAP" else model
  pages <- list(); pngs <- character(0)
  for (age in ages) {
    cfg_a <- cfg; cfg_a$central_lon <- map_center(center, age, model, cfg)
    cfg_a$view <- make_view(region, age, model, cfg) %||% projection_view(center, age, model, cfg)
    if (!is.null(cfg_a$view)) cfg_a$central_lon <- cfg_a$view$lon0
    dem <- read_dem_nc(if (!is.null(map)) map$file else find_dem(if (corrected) cfg$dir_corrected else cfg$dir_dem, age))
    sn <- if (!is.null(ice)) ice_snapshot(ice$id, ice$t_kyr, cfg) else NULL
    sl_ice <- if (!is.null(sn) && !isFALSE(ice$lower_sea) && model == "scotese_and_wright2018") min(sn$sl %||% 0, 0) else 0
    if (sl_ice < 0) dem$z <- dem$z - sl_ice              # sea lowered by the water in the ice sheets
    W <- area_matrix(dem$lat, ncol(dem$z)); ocean_pct <- 100 * sum(W[dem$z < 0]) / sum(W)
    g <- cached_grid_to_model(age, dem, model, cfg)
    cs <- if (shoreline && model == "scotese_and_wright2018") cs_lines(age, cfg) else NULL
    f <- file.path(cfg$dir_maps, sprintf("topo_%sMa_%s%s.png", format(age), tag, sub(paste0("_", model), "", suffix, fixed = TRUE)))
    ph <- 6.4 * max(1, diff(view_ylim(cfg_a)) / diff(view_xlim(cfg_a)) / 0.5)   # page height follows the map's shape
    if (!is.null(cfg_a$view)) {                     # page height so the regional map fills its area exactly
      Lm <- .bar_layout(colour_bar_pos, c(0.02, 0.98, 0.125, if (use_clim || !is.null(ice) || isTRUE(sea_ice) || layer == "sea_ice") 0.83 else 0.86), c(0.25, 0.75, 0.065, 0.09))$map
      asp <- diff(cfg_a$view$ylim) / diff(cfg_a$view$xlim)
      ph <- min(14, max(4.5, (asp * (diff(Lm[1:2]) * 11 - 0.54) + 0.2) / diff(Lm[3:4])))
    }
    if (.is_disc(cfg_a)) ph <- 10.5
    ov <- !is.null(region) && isTRUE(overview)
    clim <- .figure_climate(age, model, layer, season, climate_method, experiment, compare, map, corrected, cfg,
                            climate_flag, wind, currents, rivers, sea_ice, ice, sn, sl_ice, o2_pal, lapse_rate, use_clim, diff_with)
    if (isTRUE(seasons4) && use_clim) {                 # four seasons on one page, one colour scale
      clims <- .season_climates(function(se) .figure_climate(age, model, layer, se, climate_method, experiment, compare, map, corrected, cfg,
                                                            climate_flag, wind, currents, rivers, sea_ice, ice, sn, sl_ice, o2_pal, lapse_rate,
                                                            use_clim, diff_with))
      draw <- local({ age_ <- age; cfg_ <- cfg_a; g_ <- g; cs_ <- cs; op_ <- ocean_pct; cl_ <- clims
        function() draw_season_page(age_, model, g_, cfg_, cs_, op_, corrected, borders, modern_coast, indicators, recs,
                                    user_points, sites, center, equator, cfg, cl_, detail = detail, ice_front = ice_front) })
      ph <- if (is.null(cfg_a$view)) 7.6 else 9
      png(f, width = 11, height = ph, units = "in", res = cfg$map_dpi, type = "cairo", bg = "white")
      draw(); dev.off(); pngs <- c(pngs, f); message(f)
      pages[[length(pages) + 1]] <- list(draw = draw, h = ph)
      next
    }
    if (ov) {                                          # whole-earth map + zoom-in, joined by lines
      ph <- 8.5
      cfg_g <- cfg; cfg_g$central_lon <- map_center(center, age, model, cfg)
      cfg_g$view <- projection_view(center, age, model, cfg)
      if (!is.null(cfg_g$view)) cfg_g$central_lon <- cfg_g$view$lon0
      if (is.null(center) && is.null(cfg_g$view) && !is.null(cfg_a$view)) {   # box across the edge of the whole-earth map: centre it on the box
        ol <- .view_outline_ll(cfg_a$view)$line; ok <- is.finite(ol$lon)
        rel <- ((ol$lon[ok] - cfg_g$central_lon + 180) %% 360) - 180
        if (any(abs(diff(rel)) > 180) || any(abs(rel) > 170)) cfg_g$central_lon <- cfg_a$view$lon0
      }
    }
    draw <- local({
      age_ <- age; cfg_ <- cfg_a; g_ <- g; cs_ <- cs; op_ <- ocean_pct; cg_ <- if (ov) cfg_g else NULL; cl_ <- clim
      if (ov) function() draw_overview_page(age_, model, g_, cfg_, cg_, cs_, op_, corrected, borders, modern_coast,
                               indicators, recs, user_points, sites, center, equator, cfg, corner = overview_corner,
                               ratio = overview_ratio, overlap = overview_overlap, overview_show = overview_show,
                               detail = detail, north_arrow = north_arrow, scale_bar = scale_bar, ice_front = ice_front,
                               bar_pos = colour_bar_pos, climate = cl_)
      else function() draw_age_page(age_, model, g_, cfg_, cs_, op_, corrected, borders, modern_coast, indicators,
                               recs, user_points, sites, center, equator, cfg,
                               detail = detail, north_arrow = north_arrow, scale_bar = scale_bar, ice_front = ice_front,
                               bar_pos = colour_bar_pos, climate = cl_)
    })
    png(f, width = 11, height = ph, units = "in", res = cfg$map_dpi, type = "cairo", bg = "white")
    draw(); dev.off(); pngs <- c(pngs, f); message(f)
    pages[[length(pages) + 1]] <- list(draw = draw, h = ph)
  }
  if (comparison) {
    f <- file.path(cfg$dir_maps, sprintf("topo_%sMa_model_comparison%s.png", format(comparison_age), suffix))
    draw <- function() draw_comparison(comparison_age, corrected, borders, modern_coast, indicators, recs, cfg,
                                       user_points = user_points, sites = sites, center = center, equator = equator,
                                       ice_front = ice_front, bar_pos = colour_bar_pos)
    png(f, width = 11, height = 7.4, units = "in", res = 200, type = "cairo", bg = "white")
    draw(); dev.off(); pngs <- c(pngs, f)
    pages[[length(pages) + 1]] <- list(draw = draw, h = 7.4)
  }
  pdf_file <- file.path(cfg$dir_maps, if (all(ages >= 410 & ages <= 450)) sprintf("Silurian_topographic_maps%s.pdf", suffix) else
    sprintf("topographic_maps_%s-%sMa%s.pdf", format(min(ages)), format(max(ages)), suffix))
  if (pdf == "raster") pngs_to_pdf(pngs, pdf_file) else pages_to_vector_pdf(pages, pdf_file)
  invisible(structure(pdf_file, pngs = pngs))
}

# page layout helper: figure regions are given relative to the page; a page that is
# lower than the PDF sheet is centred vertically on it (vector PDF with mixed page sizes)
.PAGE <- new.env(); .PAGE$y <- c(0, 1)
.pf <- function(f) { y <- .PAGE$y; c(f[1], f[2], y[1] + f[3] * diff(y), y[1] + f[4] * diff(y)) }

#' one map page (title, subtitles, map, colour bar) on the current device
#' climate caption under the subtitle: wrapped to the page width (at most two lines)
.caption_text <- function(txt, y, cex) {
  wrap <- function(cx) { w <- strsplit(txt, " ", fixed = TRUE)[[1]]; lines <- character(0); cur <- ""
    for (x in w) { t2 <- if (nzchar(cur)) paste(cur, x) else x
      if (strwidth(t2, cex = cx) > 0.95 && nzchar(cur)) { lines <- c(lines, cur); cur <- x } else cur <- t2 }
    c(lines, cur) }
  lines <- wrap(cex)
  while (length(lines) > 2 && cex > 0.42) { cex <- cex - 0.04; lines <- wrap(cex) }   # smaller type for long notes
  if (length(lines) > 3) lines <- c(lines[1:2], paste(lines[-(1:2)], collapse = " "))
  dy <- if (length(lines) > 2) 0.014 else 0.019
  for (i in seq_along(lines)) text(0.02, y - (i - 1) * dy, lines[i], adj = c(0, 1), cex = cex, col = "#52514e")
}

draw_age_page <- function(age, model, g, cfg_a, cs, ocean_pct, corrected, borders, modern_coast, indicators,
                          recs, user_points, sites, center, equator, cfg,
                          detail = "auto", north_arrow = FALSE, scale_bar = FALSE, ice_front = TRUE,
                          bar_pos = "bottom", climate = NULL) {
  L <- .bar_layout(bar_pos, c(0.02, 0.98, 0.125, if (!is.null(climate$f) || !is.null(climate$ice) || !is.null(climate$sea_ice)) 0.83 else 0.86), c(0.25, 0.75, 0.065, 0.09))   # room for the climate caption
  par(mar = if (is.null(cfg_a$view)) c(0, 0, 0, 0) else c(1, 2.4, 0, 0.3), fig = .pf(L$map),
      family = "sans")
  plot.new()
  draw_map_panel(age, model, g$z, g$lat, g$lon, cfg_a, corrected, borders, modern_coast, indicators, cs, recs,
                 user_points = user_points, sites = sites, width_px = round(0.96 * 11 * cfg$map_dpi),
                 equator = equator, detail = detail, north_arrow = north_arrow, scale_bar = scale_bar,
                 ice_front = ice_front, climate = climate)
  par(fig = .pf(c(0, 1, 0, 1)), new = TRUE, mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.02, 0.985, sprintf("%s Ma \u2014 %s", format(age), stage_label(age, cfg)), adj = c(0, 1), font = 2, cex = 1.15)
  l2 <- c(if (corrected) "solid: corrected coastline incl. inferred seaways",
          if (borders) "present-day countries (black: borders, white: coastlines)" else if (modern_coast) "white: present-day coastlines")
  sub1 <- sprintf("%s (Scotese & Wright 2018), %s \u00b7 ocean %.1f%% of surface%s%s",
                  if (corrected) "Fossil-corrected PaleoDEM" else "PaleoDEM",
                  if (model == "scotese_and_wright2018") "PALEOMAP frame" else cfg$models[model], ocean_pct,
                  if (!is.null(cs)) " \u00b7 dashed: fossil-adjusted shoreline (Kocsis & Scotese 2021)" else "",
                  if (.is_disc(cfg_a)) paste0(" \u00b7 ", cfg_a$view$label) else
                  if (!is.null(cfg_a$view)) paste0(" \u00b7 region: ", cfg_a$view$label) else
                  if (!is.null(center)) sprintf(" \u00b7 centred on %s (%.0f\u00b0)",
                    if (is.character(center)) center else "chosen point", cfg_a$central_lon) else "")
  if (is.null(cfg_a$view) || .is_disc(cfg_a)) sub1 <- paste0(sub1, " \u00b7 ", PROJ_LABEL[[.proj_name(cfg_a)]], " projection")
  sub2 <- paste0(paste(l2, collapse = " \u00b7 "), if (borders || modern_coast) " (dotted = crust added later or position approximate)" else "")
  text(0.02, 0.938, sub1, adj = c(0, 1), cex = 0.68, col = "#52514e")
  if (nzchar(sub2)) text(0.02, 0.912, sub2, adj = c(0, 1), cex = 0.68, col = "#52514e")
  if (!is.null(climate$f) || !is.null(climate$ice) || !is.null(climate$sea_ice)) .caption_text(paste(c(if (!is.null(climate$f)) climate_caption(climate),
                                                                       if (!is.null(climate$ice)) ice_caption(climate$ice),
                                                                       if (!is.null(climate$sea_ice)) sea_ice_caption(climate$sea_ice)), collapse = " \u00b7 "),
                                                                  if (nzchar(sub2)) 0.889 else 0.912, 0.6)
  par(fig = .pf(L$bar), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  climate_colour_bar(L$side, climate)
}

#' outline of a regional view (its rectangle) as lon/lat, plus its four corners (bl, br, tr, tl)
.view_outline_ll <- function(v, n = 50) {
  x <- v$xlim; y <- v$ylim
  r <- rbind(cbind(seq(x[1], x[2], length.out = n), y[1]), cbind(x[2], seq(y[1], y[2], length.out = n)),
             cbind(seq(x[2], x[1], length.out = n), y[2]), cbind(x[1], seq(y[2], y[1], length.out = n)))
  ll <- sf_project(v$crs, "EPSG:4326", r, keep = TRUE, warn = FALSE)
  co <- sf_project(v$crs, "EPSG:4326", cbind(x[c(1, 2, 2, 1)], y[c(1, 1, 2, 2)]), keep = TRUE, warn = FALSE)
  list(line = data.frame(lon = ll[, 1], lat = ll[, 2]), corners = co)
}

#' one page: whole-earth map with a box around the region, and the zoomed-in map placed in a
#' corner (partly over the whole-earth map), joined by two lines
draw_overview_page <- function(age, model, g, cfg_a, cfg_g, cs, ocean_pct, corrected, borders, modern_coast,
                               indicators, recs, user_points, sites, center, equator, cfg, corner = "auto",
                               ratio = 1.25, overlap = 0.15, detail = "auto", north_arrow = FALSE, scale_bar = FALSE, ice_front = TRUE,
                               bar_pos = "bottom", climate = NULL, overview_show = character(0)) {
  # overview_show: what is also drawn on the whole-earth map: "sites" (localities + labels),
  # "user" (your points), "indicators" (reefs / evaporites / glacial deposits, as ticked for the zoom-in)
  L <- .bar_layout(bar_pos, c(0.02, 0.98, 0.115, if (!is.null(climate$f) || !is.null(climate$ice) || !is.null(climate$sea_ice)) 0.835 else 0.875), c(0.25, 0.75, 0.055, 0.08))
  A <- L$map                                          # map area on the page (fractions)
  pin <- par("din"); Wi <- diff(A[1:2]) * pin[1]; Hi <- diff(A[3:4]) * pin[2] * diff(.PAGE$y)
  # where is the region on the whole-earth map? (fractions of the map, 0-1)
  ol <- .view_outline_ll(cfg_a$view)
  cxy <- sf_project("EPSG:4326", map_crs(cfg_g), cbind(cfg_a$view$lon0, cfg_a$view$lat0), keep = TRUE, warn = FALSE)
  gx <- view_xlim(cfg_g); gy <- view_ylim(cfg_g)
  fx <- (cxy[1, 1] - gx[1]) / diff(gx); fy <- (cxy[1, 2] - gy[1]) / diff(gy)
  if (!is.finite(fx)) fx <- 0.5; if (!is.finite(fy)) fy <- 0.5
  if (identical(corner, "auto") || is.null(corner))
    corner <- paste0(if (fy > 0.5) "bottom" else "top", if (fx > 0.5) "left" else "right")
  corner <- match.arg(corner, c("bottomright", "bottomleft", "topright", "topleft"))
  zr <- grepl("right", corner); zb <- grepl("bottom", corner)
  # sizes: `ratio` = linear size of the zoom-in relative to the whole-earth map (square root of
  # their area ratio); `overlap` = how far the two panels overlap, as a fraction of the smaller
  # one (0 = touching corners, negative = a gap). Both are then scaled to fill the page.
  ga <- diff(gy) / diff(gx); za <- diff(cfg_a$view$ylim) / diff(cfg_a$view$xlim)
  gw0 <- 1 / sqrt(ga); gh0 <- sqrt(ga); zw0 <- ratio / sqrt(za); zh0 <- ratio * sqrt(za)
  ovx <- overlap * min(gw0, zw0); ovy <- overlap * min(gh0, zh0)
  mz <- c(0.22, 0.5, 0.05, 0.05)                     # margin for the zoom-in's graticule labels (in)
  mzx <- mz[2] + mz[4]; mzy <- mz[1] + mz[3]
  sc <- min((Wi - mzx) / (gw0 + zw0 - ovx), (Hi - mzy) / (gh0 + zh0 - ovy))
  gw <- sc * gw0; gh <- sc * gh0; zw <- sc * zw0; zh <- sc * zh0
  totW <- max(gw, sc * (gw0 + zw0 - ovx) + mzx); totH <- max(gh, sc * (gh0 + zh0 - ovy) + mzy)
  ox <- (Wi - totW) / 2; oy <- (Hi - totH) / 2       # centre the pair on the page
  fr <- function(x0, y0, w, h) c(A[1] + x0 / Wi * diff(A[1:2]), A[1] + (x0 + w) / Wi * diff(A[1:2]),
                                  A[3] + y0 / Hi * diff(A[3:4]), A[3] + (y0 + h) / Hi * diff(A[3:4]))
  fg <- fr(if (zr) ox else ox + totW - gw, if (zb) oy + totH - gh else oy, gw, gh)
  zx0 <- if (zr) ox + totW - mz[4] - zw else ox + mz[2]; zy0 <- if (zb) oy + mz[1] else oy + totH - mz[3] - zh
  fz_out <- fr(zx0 - mz[2], zy0 - mz[1], zw + mz[2] + mz[4], zh + mz[1] + mz[3])

  # 1. whole-earth map with the region box
  par(fig = .pf(fg), mar = c(0, 0, 0, 0), family = "sans"); plot.new()
  draw_map_panel(age, model, g$z, g$lat, g$lon, cfg_g, corrected, borders, modern_coast,
                 indicators = if ("indicators" %in% overview_show) indicators else FALSE, cs = NULL,
                 recs = recs, user_points = user_points && "user" %in% overview_show,
                 sites = if ("sites" %in% overview_show) sites else FALSE, legend = FALSE,
                 width_px = round(gw * cfg$map_dpi), equator = equator, detail = 1, ice_front = FALSE,
                 climate = if (!is.null(climate)) modifyList(climate, list(wind = FALSE, currents = FALSE)))
  draw_ll_list(list(ol$line), cfg_g, col = "#d7301f", lwd = 1.6)
  bc <- sf_project("EPSG:4326", map_crs(cfg_g), ol$corners, keep = TRUE, warn = FALSE)
  bn <- cbind(grconvertX(bc[, 1], "user", "ndc"), grconvertY(bc[, 2], "user", "ndc"))
  # 2. zoom-in on a white panel (covers the whole-earth map where they overlap)
  fz <- .pf(fr(zx0, zy0, zw, zh))                    # zoom-in map area as a fraction of the device
  zn <- cbind(fz[c(1, 2, 2, 1)], fz[c(3, 3, 4, 4)])
  par(fig = pmin(pmax(fz, 0), 1), mar = c(0, 0, 0, 0), new = TRUE); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  rect(0, 0, 1, 1, col = "white", border = NA)       # only under the map; its labels may lie over the whole-earth map
  par(fig = .pf(fz_out), mai = mz, new = TRUE); plot.new()     # margins hold the graticule labels
  draw_map_panel(age, model, g$z, g$lat, g$lon, cfg_a, corrected, borders, modern_coast, indicators, cs, recs,
                 user_points = user_points, sites = sites, width_px = round(zw * cfg$map_dpi), equator = equator,
                 detail = detail, north_arrow = north_arrow, scale_bar = scale_bar, ice_front = ice_front,
                 climate = climate)
  xl <- cfg_a$view$xlim; yl <- cfg_a$view$ylim
  rect(xl[1], yl[1], xl[2], yl[2], border = "#d7301f", lwd = 1.6)
  # 3. two connecting lines between corresponding corners (the pair whose lines pass outside
  #    both boxes); a line never runs across the zoomed-in map itself
  ok <- is.finite(bn[, 1]) & is.finite(bn[, 2])
  if (sum(ok) >= 2) {
    d <- colMeans(zn) - colMeans(bn[ok, , drop = FALSE]); nrm <- c(-d[2], d[1]) / sqrt(sum(d^2))
    s_ <- (zn - matrix(colMeans(zn), 4, 2, byrow = TRUE)) %*% nrm
    s_[!ok] <- NA; k <- c(which.min(s_), which.max(s_))
    par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    e <- 1e-4; tt <- seq(0, 1, length.out = 801)
    for (i in k) {
      px <- bn[i, 1] + tt * (zn[i, 1] - bn[i, 1]); py <- bn[i, 2] + tt * (zn[i, 2] - bn[i, 2])
      inside <- px > fz[1] + e & px < fz[2] - e & py > fz[3] + e & py < fz[4] - e
      t1 <- if (any(inside)) tt[max(1, which(inside)[1] - 1)] else 1
      if (t1 > 0) segments(bn[i, 1], bn[i, 2], bn[i, 1] + t1 * (zn[i, 1] - bn[i, 1]), bn[i, 2] + t1 * (zn[i, 2] - bn[i, 2]),
                           col = "#d7301f", lwd = 1)
    }
  }
  # titles and colour bar
  par(fig = .pf(c(0, 1, 0, 1)), new = TRUE, mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.02, 0.985, sprintf("%s Ma \u2014 %s", format(age), stage_label(age, cfg)), adj = c(0, 1), font = 2, cex = 1.15)
  sub1 <- sprintf("%s (Scotese & Wright 2018), %s \u00b7 ocean %.1f%% of surface \u00b7 region: %s \u00b7 whole earth: %s projection",
                  if (corrected) "Fossil-corrected PaleoDEM" else "PaleoDEM",
                  if (model == "scotese_and_wright2018") "PALEOMAP frame" else cfg$models[model], ocean_pct,
                  cfg_a$view$label, PROJ_LABEL[[.proj_name(cfg_g)]])
  l2 <- c(if (corrected) "solid: corrected coastline incl. inferred seaways",
          if (borders) "present-day countries (black: borders, white: coastlines)" else if (modern_coast) "white: present-day coastlines")
  text(0.02, 0.945, sub1, adj = c(0, 1), cex = 0.62, col = "#52514e")
  if (length(l2)) text(0.02, 0.922, paste(l2, collapse = " \u00b7 "), adj = c(0, 1), cex = 0.62, col = "#52514e")
  if (!is.null(climate$f) || !is.null(climate$ice) || !is.null(climate$sea_ice)) .caption_text(paste(c(if (!is.null(climate$f)) climate_caption(climate),
                                                                       if (!is.null(climate$ice)) ice_caption(climate$ice),
                                                                       if (!is.null(climate$sea_ice)) sea_ice_caption(climate$sea_ice)), collapse = " \u00b7 "),
                                                                  if (length(l2)) 0.902 else 0.922, 0.56)
  par(fig = .pf(L$bar), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  climate_colour_bar(L$side, climate)
}

#' 2 x 2 panel: the same topography in four palaeomagnetic-frame plate models
draw_comparison <- function(age, corrected, borders, modern_coast, indicators, recs, cfg = SILUCARB,
                            user_points = TRUE, sites = NULL, center = NULL, equator = TRUE, ice_front = TRUE,
                            bar_pos = "bottom") {
  dem <- read_dem_nc(find_dem(if (corrected) cfg$dir_corrected else cfg$dir_dem, age))
  pos <- list(c(0.01, 0.495, 0.53, 0.93), c(0.505, 0.99, 0.53, 0.93), c(0.01, 0.495, 0.12, 0.52), c(0.505, 0.99, 0.12, 0.52))
  L <- .bar_layout(bar_pos, c(0.01, 0.99, 0.12, 0.93), c(0.25, 0.75, 0.065, 0.087)); pos <- lapply(pos, L$tf)
  for (k in seq_along(cfg$models)) {
    key <- names(cfg$models)[k]
    cfg_k <- cfg; cfg_k$central_lon <- map_center(center, age, key, cfg)
    cfg_k$view <- projection_view(center, age, key, cfg)
    if (!is.null(cfg_k$view)) cfg_k$central_lon <- cfg_k$view$lon0
    g <- cached_grid_to_model(age, dem, key, cfg)
    par(fig = .pf(pos[[k]]), new = k > 1, mar = c(0, 0, 1.2, 0)); plot.new()
    draw_map_panel(age, key, g$z, g$lat, g$lon, cfg_k, corrected, borders, modern_coast, indicators, NULL, recs,
                   user_points = user_points, sites = sites, legend = k == 1, width_px = 1100, equator = equator,
                   ice_front = ice_front)
    title(cfg$models[key], cex.main = 0.85, adj = 0, line = 0.2)
    message("  comparison: ", key)
  }
  par(fig = .pf(c(0, 1, 0, 1)), new = TRUE, mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.01, 0.985, sprintf("%s Ma (%s) \u2014 same %sPALEOMAP topography in four plate models",
                            format(age), stage_label(age, cfg), if (corrected) "fossil-corrected " else ""), adj = c(0, 1), font = 2, cex = 1.05)
  par(fig = .pf(L$bar), new = TRUE, mar = c(0, 0, 0, 0)); plot.new(); colour_bar(L$side)
}
make_comparison <- function(age, file, ...) {           # PNG of the comparison panel
  png(file, width = 11, height = 7.4, units = "in", res = 200, type = "cairo", bg = "white")
  draw_comparison(age, ...); dev.off()
}

#' editable PDF: every page drawn again on a vector device (cairo_pdf: unicode text,
#' vector lines / points / text; the shaded topography stays an embedded image)
pages_to_vector_pdf <- function(pages, pdf_file) {
  H <- max(vapply(pages, function(p) p$h, 0))
  grDevices::cairo_pdf(pdf_file, width = 11, height = H, onefile = TRUE, bg = "white")
  on.exit({ .PAGE$y <- c(0, 1) })
  for (k in seq_along(pages)) {                         # each page's first plot.new() starts a new sheet
    h <- pages[[k]]$h / H; .PAGE$y <- c((1 - h) / 2, (1 + h) / 2)
    pages[[k]]$draw()
  }
  dev.off(); message(pdf_file)
}

#' PDF made of raster pages (thin lines never disappear in PDF viewers)
pngs_to_pdf <- function(pngs, pdf_file) {
  pdf(pdf_file, width = 11, height = 7.4, onefile = TRUE)
  for (f in pngs) {
    im <- readPNG(f); asp <- dim(im)[1] / dim(im)[2]
    par(mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    h <- min(1, asp * 11 / 7.4); w <- h / (asp * 11 / 7.4)
    rasterImage(im, (1 - w) / 2, (1 - h) / 2, (1 + w) / 2, (1 + h) / 2, interpolate = TRUE)
  }
  dev.off(); message(pdf_file)
}

#' near-surface temperature on the map grid, with the lapse-rate correction to the map heights
#' (and the ice surface of an ice run, the sea lowered by the ice)
.temp_on_map <- function(climate, z, lat, lon) {
  T <- field_on_grid(climate$f$temp, lat, lon)
  lapse <- climate$lapse %||% 6.5
  if (!isTRUE(lapse > 0) || is.null(climate$f$temp_href)) return(T)
  surf <- pmax(z, 0)
  if (!is.null(climate$ice)) {
    ic <- ice_on_latlon(climate$ice$state, lat, lon)
    surf <- ifelse(ic$H > 5, pmax(z - ic$w + ic$H, 0), surf)
  }
  temp_downscale(T, field_on_grid(climate$f$temp_href, lat, lon), surf, lapse)
}

# ---- interactive Figure map (editor) ---------------------------------------------------
#' everything a figure needs except the view: map, grid, climate layer, sea ice, ice sheets,
#' records (the same as make_maps() uses, for one age); for the editor's interactive Figure map
figure_preview_setup <- function(age, model = "scotese_and_wright2018", corrected = TRUE, layer = "topography",
                                 season = "ann", climate_method = "local", experiment = NULL, compare = NULL,
                                 map = NULL, ice = NULL, sea_ice = FALSE, rivers = FALSE, wind = FALSE, currents = FALSE,
                                 climate_flag = TRUE, o2_pal = 1, lapse_rate = 6.5, shoreline = FALSE, diff_with = NULL, cfg = SILUCARB) {
  if (!is.null(map)) { map <- climate_map(map, cfg = cfg); age <- map$age; corrected <- TRUE }
  if (layer %in% c("ice_thickness", "ice_smb") && is.null(ice)) layer <- "topography"
  if (layer == "sea_ice") sea_ice <- TRUE
  use_clim <- !layer %in% c("topography", "ice_thickness", "ice_smb", "sea_ice") || isTRUE(wind) || isTRUE(currents) || isTRUE(rivers)
  if (use_clim && model != "scotese_and_wright2018") { use_clim <- FALSE; layer <- "topography"; wind <- FALSE; currents <- FALSE }
  experiment <- experiment %||% climate_experiment()
  dem <- read_dem_nc(if (!is.null(map)) map$file else find_dem(if (corrected) cfg$dir_corrected else cfg$dir_dem, age))
  sn <- if (!is.null(ice)) ice_snapshot(ice$id, ice$t_kyr, cfg) else NULL
  sl_ice <- if (!is.null(sn) && !isFALSE(ice$lower_sea) && model == "scotese_and_wright2018") min(sn$sl %||% 0, 0) else 0
  if (sl_ice < 0) dem$z <- dem$z - sl_ice
  g <- cached_grid_to_model(age, dem, model, cfg)
  clim <- .figure_climate(age, model, layer, season, match.arg(climate_method, c("local", "model", "rerun")), experiment, compare,
                          map, corrected, cfg, climate_flag, wind, currents, rivers, sea_ice, ice, sn, sl_ice, o2_pal, lapse_rate, use_clim, diff_with)
  list(age = age, model = model, corrected = corrected, g = g, clim = clim, cache = new.env(),
       cs = if (isTRUE(shoreline) && model == "scotese_and_wright2018") cs_lines(age, cfg) else NULL,
       recs = all_records(cfg))
}

#' draw the Figure map for a window (cx, cy: centre lon / lat; w, h: width / height in degrees)
#' in an equirectangular view, with all layers and overlays as in the rendered figure.
#' Returns the plot coordinates to convert mouse positions (usr, plt, k = metres per degree, cx).
figure_preview <- function(S, window = NULL, borders = TRUE, modern_coast = FALSE, indicators = TRUE,
                           user_points = TRUE, sites = NULL, equator = TRUE, ice_front = TRUE,
                           north_arrow = FALSE, scale_bar = FALSE, width_px = 1100, cfg = SILUCARB) {
  w <- window %||% list(cx = cfg$central_lon, cy = 0, w = 360, h = 180)
  k <- 6371000 * pi / 180
  cfg$view <- list(crs = sprintf("+proj=eqc +lon_0=%g +R=6371000 +units=m +no_defs", w$cx),
                   xlim = c(-w$w, w$w) / 2 * k, ylim = (w$cy + c(-w$h, w$h) / 2) * k,
                   lat0 = w$cy, lon0 = w$cx, far_deg = 180.5,
                   label = sprintf("%.1f\u00b0 \u00d7 %.1f\u00b0 around %.1f\u00b0, %.1f\u00b0", w$w, w$h, w$cy, w$cx))
  cfg$central_lon <- w$cx; cfg$map_refine_max <- 1; cfg$map_relief_smooth_km <- 0
  cfg$view$fast_eqc <- TRUE                          # plain lon/lat view: fast raster path, shading cached in S$cache
  g <- S$g; dl <- abs(g$lat[2] - g$lat[1])
  st <- max(1, floor(1.5 * w$w / width_px / dl))      # coarser grid for wide windows (fast)
  if (st > 1) { ii <- seq(1, nrow(g$z), st); jj <- seq(1, ncol(g$z), st); g <- list(z = g$z[ii, jj], lat = g$lat[ii], lon = g$lon[jj]) }
  fm <- c(0, 1, 0.1, 1)                              # map above, colour bar below
  par(fig = fm, mar = c(0, 0, 0, 0)); plot.new()
  draw_map_panel(S$age, S$model, g$z, g$lat, g$lon, cfg, S$corrected, borders, modern_coast, indicators, S$cs, S$recs,
                 user_points = user_points, sites = sites, legend = TRUE, width_px = width_px, equator = equator,
                 detail = 1, north_arrow = north_arrow, scale_bar = scale_bar, ice_front = ice_front, climate = S$clim,
                 rgb_cache = S$cache)
  usr <- par("usr"); pl <- par("plt")
  plt <- c(fm[1] + pl[1:2] * (fm[2] - fm[1]), fm[3] + pl[3:4] * (fm[4] - fm[3]))   # plot region as fractions of the image
  par(fig = c(0.25, 0.75, 0.045, 0.075), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  try(climate_colour_bar(1, S$clim %||% list(layer = "topography")), silent = TRUE)
  invisible(list(usr = usr, plt = plt, k = k, cx = w$cx))
}
