# =============================================================================
#  Climate layers: precipitation, 10 m wind and surface ocean currents.
#
#  Source: the HadCM3L climate-model simulations of Valdes et al. (2021, Climate of the
#  Past 17, 1483-1506, doi:10.5194/cp-17-1483-2021; model: Valdes et al. 2017, GMD 10, 3715), "scotese_02"
#  series (Foster et al. 2017 CO2; runs texqe = 541 Ma ... texpa1 = 0 Ma; verified 2026-10-20 against the
#  BRIDGE run list .../html_bridge/scotese_02.html and Zenodo 10.5281/zenodo.8237751), run on the same
#  PALEOMAP palaeogeographies as the PaleoDEMs. Output from the BRIDGE group, University of
#  Bristol (https://www.paleo.bristol.ac.uk/ummodel/scripts/papers/Valdes_et_al_2021.html).
#  Per age only ~3-6 MB is downloaded once (data/climate/); nothing is simulated here, so
#  this runs on any desktop PC.
#
#  Fitting to the corrected map (hybrid):
#   - precipitation and wind are interpolated smoothly (they exist over land and sea);
#   - currents are shown only where the corrected map is sea. Where the corrected map has sea
#     but the model had land (e.g. inferred seaways), the current is ESTIMATED from the model's
#     wind: ~3 % of the wind speed, turned 45 deg (Ekman; right in the north, left in the
#     south) and kept parallel to nearby coasts. These arrows are drawn lighter.
# =============================================================================

# experiment id -> model age (Ma) and CO2 (ppm), scotese_02 series (CO2 after Foster et al. 2017;
# from ExperimentInfo.xlsx of the PhanDA HadCM3L priors, Zenodo 10.5281/zenodo.8237751)
PAL_PPM <- 280                                   # 1 PAL = pre-industrial CO2
CLIMATE_EXPS <- data.frame(
  id = c("texqe","texqd","texqc","texqb","texqa","teXPz","teXPy","teXPx","teXPw","teXPv","teXPu","teXPt","teXPs",
         "teXPr","teXPq","teXPp","teXPo","teXPn","teXPm","teXPl","teXPk","teXPj","teXPi","teXPh","teXPg","teXPf",
         "teXPe","teXPd","teXPc","teXPb","teXPa","teXpz","teXpy","teXpx","teXpw","teXpv","teXpu","teXpt","teXps",
         "teXpr","teXpq","teXpp","teXpo","teXpn","teXpm","teXpl","teXpk","teXpj","teXpi","teXph","teXpg","teXpf",
         "teXpe","teXpd","teXpc","teXpb","teXpa","texPz","texPy","texPx","texPw2","texPv1","texPu1","texPt1",
         "texPs1","texPr1","texPq1","texPp1","texPo1","texPn1","texPm1","texPl1","texPk1","texPj1","texPi1",
         "texPh2","texPg1","texPf1","texPe2","texPd2","texPc2","texPb2","texPa2","texpz1","texpy1","texpx2",
         "texpw2","texpv1","texpu1","texpt2","texps2","texpr2","texpq","texpp","texpo1","texpn1","texpm1",
         "texpl1","texpk2","texpj2","texpi1","texph1","texpg1","texpf","texpe","texpd","texpc","texpb","texpa1"),
  age = c(541,535,530,525,520,515,510,505,498.8,495.5,491.8,485.4,481.6,475,470,465,460,455.7,449.1,444.5,
          441.2,436,430.4,425.2,421.1,415,410,405,400,395,390.5,385.2,380,375,370,365.6,358.9,354,349,344,
          338.8,333,330.9,327,319.2,314.6,311.1,305.4,301.3,297,292.6,286.8,280,275,268.7,265.1,262.5,256,252,
          244.6,239.5,233.6,232,227,222.4,217.8,213.2,204.9,201.3,196,190.8,186.8,178.4,172.2,168.2,164.8,
          160.4,154.7,148.6,145,142.4,136.4,131.2,127.2,121.8,115.8,111,107,102.6,97.2,91.9,86.7,80.8,75,69,
          66,60.6,55.8,52.2,44.5,39.5,35.9,31,25.6,19.5,14.9,10.5,4,0),
  co2 = c(3374,3290,3220,3150,3080,3010,2940,2870,2776.39,2711.48,2638.69,2512.79,2438.03,2308.2,2209.84,
          2111.48,2013.11,1928.52,1798.69,1708.2,1643.28,1575.99,1531.16,1489.54,1456.72,1579.21,2102.11,
          1689.11,1731,1296.84,1093.36,1377.28,1269.38,1052.02,810.78,805.89,587.2,588.87,644.64,564.74,232.91,
          262.84,296.31,316.52,328.44,327.16,338.13,358.05,392.99,444.7,495.35,625.67,690.31,555.97,310.63,
          320.82,352.02,810.9,878.74,418.58,1033.59,1491.97,1613.52,2059.04,1810.4,1481.01,1223.32,1503.18,
          1728.9,1783.74,1180.68,948.81,986.22,1046.43,1018.89,867.55,785.38,802.03,721.3,713.01,831.71,
          1020.29,896.15,897.93,986.19,1102.96,1148.07,1028.16,1007.8,707.39,465.65,590.28,667.16,559.27,
          261.65,228.84,334.83,570.28,735.76,751.41,795.57,901.42,763.85,502.09,338.28,310.19,299.4,298.33,
          276.01),
  stringsAsFactors = FALSE)
CLIMATE_SEASONS <- c(ann = "annual mean", djf = "December-February", mam = "March-May", jja = "June-August",
                     son = "September-November", jan = "January", feb = "February", mar = "March", apr = "April",
                     may = "May", jun = "June", jul = "July", aug = "August", sep = "September", oct = "October",
                     nov = "November", dec = "December")
#' months of a season name
season_months <- function(season) switch(season, ann = 1:12, djf = c(12, 1, 2), mam = 3:5, jja = 6:8, son = 9:11,
                                         match(season, tolower(month.abb)))
CLIMATE_URL <- "https://www.paleo.bristol.ac.uk/ummodel/data/%s/climate/%s"

#' download with up to 3 tries (the server sometimes drops a connection)
.fetch_retry <- function(url, dest, tries = 3) {
  for (k in seq_len(tries)) { if (file.exists(dest) || isTRUE(fetch(url, dest))) return(invisible(TRUE)); Sys.sleep(2) }
  invisible(file.exists(dest))
}

#' the model run closest in age
climate_exp <- function(age) CLIMATE_EXPS[which.min(abs(CLIMATE_EXPS$age - age)), ]

#' download the files for one age and season (once; ~1.3 MB atmosphere + 0.5 MB ocean)
climate_files <- function(age, season = "ann", cfg = SILUCARB) {
  e <- climate_exp(age); d <- file.path(cfg$dir_data, "climate", e$id)
  f <- c(atm = file.path(d, sprintf("%sa.pdcl%s.nc", e$id, season)),
         ocn = file.path(d, sprintf("%so.pfcl%s.nc", e$id, season)))
  for (k in names(f)) .fetch_retry(sprintf(CLIMATE_URL, e$id, basename(f[[k]])), f[[k]])
  if (!all(file.exists(f))) stop("climate data for ", age, " Ma could not be downloaded (see the message above)", call. = FALSE)
  f
}

#' one variable as a [lat x lon] matrix, latitude increasing, longitude -180..180
.read_clim_var <- function(file, var) {
  n <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(n))
  v <- n$var[[var]]; z <- ncdf4::ncvar_get(n, var)            # [lon x lat]
  lon <- as.vector(v$dim[[1]]$vals); lat <- as.vector(v$dim[[2]]$vals)
  z[abs(z) > 1e19] <- NA
  m <- t(z)                                                   # [lat x lon]
  o <- order(lat); m <- m[o, , drop = FALSE]; lat <- lat[o]
  lon <- ((lon + 180) %% 360) - 180; p <- order(lon)
  list(z = m[, p, drop = FALSE], lat = lat, lon = lon[p])
}

#' the climate fields for a map age: precipitation (mm/day), 10 m wind (m/s) and surface
#' current (cm/s, NA on model land); cached in memory
.CLIM <- new.env()
climate_fields <- function(age, season = "ann", cfg = SILUCARB) {
  season <- match.arg(season, names(CLIMATE_SEASONS))
  key <- paste(climate_exp(age)$id, season)
  if (!is.null(.CLIM[[key]])) return(.CLIM[[key]])
  f <- climate_files(age, season, cfg); e <- climate_exp(age)
  pr <- .read_clim_var(f[["atm"]], "precip_mm_srf"); pr$z <- pr$z * 86400   # kg m-2 s-1 -> mm/day
  out <- list(exp = e$id, model_age = e$age, season = season, hadcm3l_co2 = e$co2,
              precip = pr, temp = local({ t <- .read_clim_var(f[["atm"]], "temp_mm_1_5m"); t$z <- t$z - 273.15; t }),
              wind_u = .read_clim_var(f[["atm"]], "u_mm_10m"), wind_v = .read_clim_var(f[["atm"]], "v_mm_10m"),
              curr_u = .read_clim_var(f[["ocn"]], "ucurrTot_mm_dpth"), curr_v = .read_clim_var(f[["ocn"]], "vcurrTot_mm_dpth"))
  # water balance: precipitation - evaporation (evaporation from the latent heat flux) and runoff
  ev <- .read_clim_var(f[["atm"]], "lh_mm_srf"); ev$z <- ev$z / 2.501e6 * 86400
  out$evap <- ev; out$pme <- pr; out$pme$z <- pr$z - ev$z
  ro <- .read_clim_var(f[["atm"]], "srfRunoff_mm_srf"); rs <- .read_clim_var(f[["atm"]], "subsrfRunoff_mm_srf")
  out$runoff <- ro; out$runoff$z <- (ifelse(is.finite(ro$z), ro$z, 0) + ifelse(is.finite(rs$z), rs$z, 0)) * 86400
  out$runoff$z[!is.finite(ro$z) & !is.finite(rs$z)] <- NA
  # ocean (HadCM3L as run): salinity at 5 m (psu), vertical velocity at 10 m (m/day, + = up), mixed layer (m)
  fillsea <- function(F) { if (exists(".fill_from")) F$z <- .fill_from(F$z, is.finite(F$z), iters = 60); F }  # into new seaways
  sal <- .read_clim_var(f[["ocn"]], "salinity_mm_dpth"); sal$z <- sal$z * 1000 + 35; out$salinity <- fillsea(sal)
  w <- .read_clim_var(f[["ocn"]], "W_mm_dpth"); w$z <- w$z * 864; out$w_ocean <- fillsea(w)
  out$mld <- fillsea(.read_clim_var(f[["ocn"]], "mixLyrDpth_mm_uo"))
  sst <- .read_clim_var(f[["ocn"]], "temp_mm_dpth"); out$sst_ocean <- fillsea(sst)          # C, 5 m
  .CLIM[[key]] <- out
  out
}

#' bilinear interpolation of a global field at points (wraps in longitude, NA stays NA)
clim_interp <- function(F, qlat, qlon) {
  lat <- F$lat; lon <- F$lon; z <- F$z; nl <- length(lon)
  lonx <- c(lon, lon[1] + 360); zx <- cbind(z, z[, 1])
  ql <- ((qlon - lon[1]) %% 360) + lon[1]
  j <- pmin(findInterval(ql, lonx, all.inside = TRUE), nl)
  i <- findInterval(pmin(pmax(qlat, lat[1]), lat[length(lat)]), lat, all.inside = TRUE)
  tx <- (ql - lonx[j]) / (lonx[j + 1] - lonx[j])
  ty <- pmin(pmax((qlat - lat[i]) / (lat[i + 1] - lat[i]), 0), 1)
  z00 <- zx[cbind(i, j)]; z01 <- zx[cbind(i, j + 1)]; z10 <- zx[cbind(i + 1, j)]; z11 <- zx[cbind(i + 1, j + 1)]
  (1 - ty) * ((1 - tx) * z00 + tx * z01) + ty * ((1 - tx) * z10 + tx * z11)
}

#' precipitation on a lat-lon grid (e.g. the map grid), smoothed so no model cells show
precip_on_grid <- function(cl, lat, lon) {
  g <- expand.grid(lon = lon, lat = lat)
  m <- matrix(clim_interp(cl$precip, g$lat, g$lon), length(lat), length(lon), byrow = TRUE)
  m[is.na(m)] <- 0
  pmax(m, 0)
}

# ---- colours ------------------------------------------------------------------------
.PRECIP_STOPS <- list(p = c(0, 0.5, 1, 2, 3, 5, 8, 12),
                      col = c("#8c510a", "#c7883a", "#e9cf8a", "#e8f0c0", "#9fd39a", "#3fae8e", "#1f78b4", "#283a8f"))
precip_colour <- function(p) {
  x <- pmin(pmax(p, 0), 12); s <- .PRECIP_STOPS; pos <- seq(0, 1, length.out = length(s$p))
  f <- approx(s$p, pos, xout = x)$y; rgbm <- col2rgb(s$col) / 255
  sapply(1:3, function(k) approx(pos, rgbm[k, ], xout = f)$y)
}
#' RGB image of precipitation with relief shading on land; the sea is drawn paler
render_precip_rgb <- function(z, lat, pr) {
  m <- precip_colour(as.vector(pr)); rgb <- array(m, c(nrow(z), ncol(z), 3))
  hs <- hillshade(z, lat); land <- z >= 0
  for (c in 1:3) {
    ch <- rgb[, , c]
    ch[land] <- pmin(pmax(ch[land] * (0.55 + 0.6 * hs[land]), 0), 1)
    ch[!land] <- 0.55 * ch[!land] + 0.45 * 0.93
    rgb[, , c] <- ch
  }
  rgb
}
colour_bar_precip <- function(side = 1, season = "ann") {
  s <- .PRECIP_STOPS; k <- length(s$p) - 1
  pb <- seq(0, k, length.out = 281)                         # bar position: one unit per colour stop
  val <- approx(seq(0, k), s$p, xout = (pb[-1] + pb[-length(pb)]) / 2)$y
  cols <- precip_colour(val); cc <- rgb(cols[, 1], cols[, 2], cols[, 3]); n <- length(pb)
  labs <- c(as.character(s$p[-length(s$p)]), "12+")
  if (side %in% c(1, 3)) {
    plot.window(xlim = c(0, k), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(pb[-n], 0, pb[-1], 1, col = cc, border = NA)
    axis(side, at = 0:k, labels = labs, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    mtext(sprintf("Precipitation, %s (mm/day)", CLIMATE_SEASONS[[season]]), side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = c(0, k), xaxs = "i", yaxs = "i")
    rect(0, pb[-n], 1, pb[-1], col = cc, border = NA)
    axis(side, at = 0:k, labels = labs, las = 1, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    mtext("Precipitation (mm/day)", side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}

# ---- arrows -----------------------------------------------------------------------------
#' arrow positions: a regular grid in map units, inside the map outline
.arrow_points <- function(cfg, n_across = 34) {
  xl <- view_xlim(cfg); yl <- view_ylim(cfg); s <- diff(xl) / n_across
  g <- expand.grid(x = seq(xl[1] + s / 2, xl[2], by = s), y = seq(yl[1] + s / 2, yl[2], by = s))
  if (is.null(cfg$view)) {
    E <- .global_ext(cfg); o <- order(E$ey)
    xe <- approx(E$ey[o], abs(E$ex[o]), g$y, ties = max, rule = 2)$y
    g <- g[abs(g$x) < xe * 0.985, ]
  } else if (.is_disc(cfg)) g <- g[g$x^2 + g$y^2 < (0.985 * cfg$view$disc)^2, ]
  ll <- .inv(cfg, g$x, g$y); ok <- is.finite(ll[, 1]) & is.finite(ll[, 2])
  if (!is.null(cfg$view)) ok <- ok & !.far(ll[, 1], ll[, 2], cfg)
  list(x = g$x[ok], y = g$y[ok], lon = ll[ok, 1], lat = ll[ok, 2], spacing = s)
}

#' draw vectors (u east, v north) at points; the arrow length is proportional to speed
#' (ref = speed drawn one grid spacing long)
draw_vectors <- function(P, u, v, cfg, ref, col, lwd = 0.8) {
  sp <- sqrt(u^2 + v^2); ok <- is.finite(sp) & sp > 0.02 * ref
  if (!any(ok)) return(invisible())
  lon <- P$lon[ok]; lat <- P$lat[ok]; ue <- u[ok] / sp[ok]; vn <- v[ok] / sp[ok]
  d <- 20 / 111.2                                         # a 20 km step along the vector
  lat2 <- pmin(pmax(lat + vn * d, -89.9), 89.9); lon2 <- lon + ue * d / pmax(cos(lat * pi / 180), 0.05)
  q <- .fwd(cfg, lon2, lat2); dx <- q[, 1] - P$x[ok]; dy <- q[, 2] - P$y[ok]
  nn <- sqrt(dx^2 + dy^2); good <- is.finite(nn) & nn > 0
  L <- P$spacing * 0.95 * pmin(sqrt(sp[ok] / ref), 1.4)   # length ~ square root of speed (slow flow stays visible)
  x0 <- P$x[ok] - 0.5 * L * dx / nn; y0 <- P$y[ok] - 0.5 * L * dy / nn
  x1 <- P$x[ok] + 0.5 * L * dx / nn; y1 <- P$y[ok] + 0.5 * L * dy / nn
  g <- good & L > 0.12 * P$spacing
  arrows(x0[g], y0[g], x1[g], y1[g], length = 0.035, angle = 25, col = col, lwd = lwd)
}

#' surface currents fitted to the corrected land/sea: model currents where the model has
#' sea, a wind-driven estimate in sea the model did not have; nothing on corrected land
currents_at <- function(cl, P, z, lat, lon) {
  # corrected land/sea at the points and a land fraction smoothed over ~150 km (coast direction)
  dl <- abs(lat[2] - lat[1])
  ii <- pmin(pmax(round((P$lat - lat[1]) / dl) + 1, 1), length(lat))
  jj <- pmin(pmax(round((((P$lon - lon[1]) %% 360)) / dl) + 1, 1), length(lon))
  sea <- z[cbind(ii, jj)] < 0
  u <- clim_interp(cl$curr_u, P$lat, P$lon); v <- clim_interp(cl$curr_v, P$lat, P$lon)
  est <- sea & !(is.finite(u) & is.finite(v))
  if (any(est)) {
    wu <- clim_interp(cl$wind_u, P$lat[est], P$lon[est]); wv <- clim_interp(cl$wind_v, P$lat[est], P$lon[est])
    a <- -pi / 4 * pmax(-1, pmin(1, P$lat[est] / 10))      # 45 deg, right in the N, left in the S, 0 at the equator
    eu <- 3 * (wu * cos(a) - wv * sin(a)); ev <- 3 * (wu * sin(a) + wv * cos(a))   # 3 % of the wind, in cm/s
    # keep them parallel to the coast: remove the part flowing onto land
    lf <- .land_frac_cache(z, lat)
    gy <- (lf[pmin(ii[est] + 2, nrow(lf)) + nrow(lf) * (jj[est] - 1)] - lf[pmax(ii[est] - 2, 1) + nrow(lf) * (jj[est] - 1)])
    jr <- ((jj[est] + 1) %% ncol(lf)) + 1; jl <- ((jj[est] - 3) %% ncol(lf)) + 1
    gx <- (lf[ii[est] + nrow(lf) * (jr - 1)] - lf[ii[est] + nrow(lf) * (jl - 1)]) / pmax(cos(P$lat[est] * pi / 180), 0.1)
    gn <- sqrt(gx^2 + gy^2); nx <- ifelse(gn > 0, gx / gn, 0); ny <- ifelse(gn > 0, gy / gn, 0)
    on <- eu * nx + ev * ny; on <- pmax(on, 0)
    u[est] <- eu - on * nx; v[est] <- ev - on * ny
  }
  u[!sea] <- NA; v[!sea] <- NA
  list(u = u, v = v, est = est & sea)
}
.LF <- new.env()
.land_frac_cache <- function(z, lat) {
  key <- paste(dim(z), sum(z >= 0), collapse = "_")
  if (is.null(.LF[[key]])) .LF[[key]] <- smooth_km((z >= 0) * 1, lat, 120)
  .LF[[key]]
}

#' climate overlays on a map panel (after the topography / precipitation image)
#' returns legend rows and the reference-arrow key text
draw_climate_arrows <- function(cl, cfg, z, lat, lon, wind = FALSE, currents = FALSE) {
  if (!wind && !currents) return(invisible(NULL))
  P <- .arrow_points(cfg, if (is.null(cfg$view)) 40 else 26)
  if (!length(P$x)) return(invisible(NULL))
  if (currents) {
    cu <- currents_at(cl, P, z, lat, lon)
    draw_vectors(P, ifelse(cu$est, NA, cu$u), ifelse(cu$est, NA, cu$v), cfg, ref = 20, col = "#c2185b", lwd = 0.9)
    draw_vectors(P, ifelse(cu$est, cu$u, NA), ifelse(cu$est, cu$v, NA), cfg, ref = 20, col = adjustcolor("#c2185b", 0.45), lwd = 0.9)
  }
  if (wind) {
    if (currents) { P2 <- P; P2$x <- P$x + 0.5 * P$spacing; P2$y <- P$y + 0.5 * P$spacing   # between the current arrows
      ll <- .inv(cfg, P2$x, P2$y); P2$lon <- ll[, 1]; P2$lat <- ll[, 2] } else P2 <- P
    wu <- clim_interp(cl$wind_u, P2$lat, P2$lon); wv <- clim_interp(cl$wind_v, P2$lat, P2$lon)
    draw_vectors(P2, wu, wv, cfg, ref = cl$wind_ref %||% 8, col = "#1b1b1b", lwd = 0.7)
  }
  invisible(structure(P$spacing, wind_ref = cl$wind_ref %||% 8))
}

#' key with reference arrows (top-left of the map)
draw_climate_key <- function(cfg, spacing, wind, currents, any_est = TRUE) {
  wref <- attr(spacing, "wind_ref") %||% 8
  if (!wind && !currents) return(invisible())
  xl <- view_xlim(cfg); yl <- view_ylim(cfg); W <- diff(xl); H <- diff(yl)
  x0 <- xl[1] + 0.012 * W; y <- yl[2] - (if (is.null(cfg$view)) 0.07 else 0.035) * H; L <- spacing * 0.95; dy <- 0.035 * H
  rows <- c(if (wind) "wind", if (currents) "cur", if (currents && any_est) "est")
  rect(x0 - 0.006 * W, y - (length(rows) - 0.4) * dy, x0 + L + 0.2 * W, y + 0.6 * dy,
       col = adjustcolor("white", 0.85), border = "#c9c8c2", lwd = 0.6)
  for (r in rows) {
    col <- switch(r, wind = "#1b1b1b", cur = "#c2185b", est = adjustcolor("#c2185b", 0.45))
    arrows(x0, y, x0 + L, y, length = 0.035, angle = 25, col = col, lwd = 0.9)
    text(x0 + L + 0.008 * W, y, switch(r, wind = if (wref == 8) "wind (10 m), 8 m/s (length ~ \u221aspeed)" else sprintf("wind change, %g m/s (length ~ \u221aspeed)", wref), cur = "surface current, 20 cm/s",
                                        est = "current estimated from wind (sea not in model)"),
         adj = c(0, 0.5), cex = 0.55, col = "#222222")
    y <- y - dy
  }
}

#' one-line source note for the map subtitle
climate_caption <- function(climate) {
  f <- climate$f
  how <- switch(f$method %||% "model",
    model = "as run",
    local = "with a local correction for the corrected geography",
    rerun = sprintf("+ SPEEDY response to the corrected geography%s (%s yr)",
                    if (!is.null(f$experiment)) paste0(" and ", f$experiment$label) else "", format(f$speedy_years)))
  chg <- if (climate$layer %in% c("precipitation_change", "temperature_change") && !is.null(f$change))
    paste0(" \u00b7 colours: ", if (!is.null(f$compare)) paste0((f$experiment$label %||% "reference"), " minus ", f$compare$label)
           else "experiment - reference (HadCM3L orbit, CO2 and Sun)", ", SPEEDY") else ""
  sprintf("climate: HadCM3L, BRIDGE series scotese_02 (Valdes et al. 2021, Clim. Past 17, 1483), run %s, %s Ma (CO2 %.0f ppm), %s, %s%s%s%s%s", f$exp, format(f$model_age),
          f$hadcm3l_co2 %||% NA, CLIMATE_SEASONS[[f$season]], how, chg,
          if (!is.null(f$flag)) " \u00b7 dashed: geography differs much from the model's (least reliable)" else "",
          if (!is.null(climate$extra)) paste0(" \u00b7 ", climate$extra) else "",
          if (!is.null(climate$rivers)) paste0(" \u00b7 ", river_caption(climate$rivers)) else "")
}

# ---- near-surface air temperature ------------------------------------------------------
# white at exactly 0 degC: blues below freezing, yellows - oranges - reds above
.TEMP_STOPS <- list(t = c(-30, -20, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40),
                    col = c("#2c1e6e", "#3b4cc0", "#6a9bd8", "#bcdcf0", "#ffffff", "#fff4c0", "#fde68a",
                            "#fbc561", "#f6a253", "#e8693a", "#c93127", "#981c2a", "#5c0f24"))
temp_colour <- function(t) {
  s <- .TEMP_STOPS; x <- pmin(pmax(t, min(s$t)), max(s$t)); pos <- seq(0, 1, length.out = length(s$t))
  f <- approx(s$t, pos, xout = x)$y; rgbm <- col2rgb(s$col) / 255
  sapply(1:3, function(k) approx(pos, rgbm[k, ], xout = f)$y)
}
render_temp_rgb <- function(z, lat, tt) {
  rgb <- array(temp_colour(as.vector(tt)), c(nrow(z), ncol(z), 3)); hs <- hillshade(z, lat); land <- z >= 0
  for (c in 1:3) { ch <- rgb[, , c]; ch[land] <- pmin(pmax(ch[land] * (0.6 + 0.5 * hs[land]), 0), 1)
    ch[!land] <- 0.8 * ch[!land] + 0.2; rgb[, , c] <- ch }        # sea a little paler (0 degC stays white)
  rgb
}
field_on_grid <- function(F, lat, lon) {
  g <- expand.grid(lon = lon, lat = lat)
  m <- matrix(clim_interp(F, g$lat, g$lon), length(lat), length(lon), byrow = TRUE); m[is.na(m)] <- mean(m, na.rm = TRUE); m
}
colour_bar_temp <- function(side = 1, season = "ann") {
  s <- .TEMP_STOPS; v <- seq(min(s$t), max(s$t), length.out = 281); m <- temp_colour(v)
  cc <- rgb(m[, 1], m[, 2], m[, 3]); n <- length(v); at <- s$t
  lab <- sprintf("Near-surface air temperature, %s (\u00b0C)", CLIMATE_SEASONS[[season]])
  if (side %in% c(1, 3)) {
    plot.window(xlim = range(v), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(v[-n], 0, v[-1], 1, col = cc[-1], border = NA)
    axis(side, at = at, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    mtext(lab, side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = range(v), xaxs = "i", yaxs = "i")
    rect(0, v[-n], 1, v[-1], col = cc[-1], border = NA)
    axis(side, at = at, las = 1, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    mtext("Temperature (\u00b0C)", side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}

#' the colour bar of a climate layer (or the elevation bar)
climate_colour_bar <- function(side, climate) {
  if (isTRUE(climate$layer %in% names(SCALAR_LAYERS)) && !is.null(climate$scale)) colour_bar_scalar(side, climate$layer, climate$season, climate$scale) else
  if (isTRUE(climate$layer %in% c("ice_thickness", "ice_smb")) && !is.null(climate$ice)) colour_bar_ice(side, climate$layer) else
  if (identical(climate$layer, "precipitation")) colour_bar_precip(side, climate$season) else
  if (identical(climate$layer, "temperature")) colour_bar_temp(side, climate$season) else
  if (isTRUE(climate$layer %in% c("precipitation_change", "temperature_change")))
    colour_bar_change(side, climate$layer, climate$season, climate$lim %||% CHANGE_LIMITS[[climate$layer]],
                      what = if (!is.null(climate$f$compare)) "difference" else "change") else
  colour_bar(side)
}

# ---- change maps (experiment - reference) ---------------------------------------------
.diverging <- function(x, lim, neg = c("#8c510a", "#d8b365", "#f6e8c3"), pos = c("#c7eae5", "#5ab4ac", "#01665e")) {
  cols <- col2rgb(c(neg, "#f7f7f7", pos)) / 255; pos_ <- seq(-1, 1, length.out = ncol(cols))
  f <- pmax(pmin(x / lim, 1), -1)
  sapply(1:3, function(k) approx(pos_, cols[k, ], xout = f)$y)
}
CHANGE_LIMITS <- c(precipitation_change = 3, temperature_change = 5)
#' a round colour-scale limit for a change field (98 % of the values inside)
change_limit <- function(F, layer) {
  if (is.null(F)) return(CHANGE_LIMITS[[layer]])
  q <- stats::quantile(abs(F$z), 0.98, na.rm = TRUE)
  if (!is.finite(q) || q <= 0) return(CHANGE_LIMITS[[layer]])
  nice <- c(0.2, 0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30); nice[which(nice >= q)[1]] %||% 30
}
render_change_rgb <- function(z, lat, v, layer, lim = CHANGE_LIMITS[[layer]]) {
  m <- if (layer == "temperature_change") .diverging(as.vector(v), lim,
         neg = c("#2166ac", "#67a9cf", "#d1e5f0"), pos = c("#fddbc7", "#ef8a62", "#b2182b"))
       else .diverging(as.vector(v), lim)
  rgb <- array(m, c(nrow(z), ncol(z), 3)); hs <- hillshade(z, lat); land <- z >= 0
  for (c in 1:3) { ch <- rgb[, , c]; ch[land] <- pmin(pmax(ch[land] * (0.6 + 0.5 * hs[land]), 0), 1); rgb[, , c] <- ch }
  rgb
}
change_on_grid <- function(F, lat, lon) {
  g <- expand.grid(lon = lon, lat = lat)
  m <- matrix(clim_interp(F, g$lat, g$lon), length(lat), length(lon), byrow = TRUE); m[is.na(m)] <- 0; m
}
colour_bar_change <- function(side = 1, layer, season = "ann", lim = CHANGE_LIMITS[[layer]], what = "change") {
  v <- seq(-lim, lim, length.out = 201)
  m <- if (layer == "temperature_change") .diverging(v, lim, neg = c("#2166ac", "#67a9cf", "#d1e5f0"), pos = c("#fddbc7", "#ef8a62", "#b2182b"))
       else .diverging(v, lim)
  cc <- rgb(m[, 1], m[, 2], m[, 3]); n <- length(v)
  lab <- sprintf("%s %s, %s (%s)", if (layer == "temperature_change") "Near-surface temperature" else "Precipitation", what,
                 CLIMATE_SEASONS[[season]], if (layer == "temperature_change") "K" else "mm/day")
  at <- pretty(c(-lim, lim)); at <- at[abs(at) <= lim]
  if (side %in% c(1, 3)) {
    plot.window(xlim = c(-lim, lim), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(v[-n], 0, v[-1], 1, col = cc[-1], border = NA)
    axis(side, at = at, labels = ifelse(at > 0, paste0("+", at), at), cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    mtext(lab, side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = c(-lim, lim), xaxs = "i", yaxs = "i")
    rect(0, v[-n], 1, v[-1], col = cc[-1], border = NA)
    axis(side, at = at, labels = ifelse(at > 0, paste0("+", at), at), las = 1, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    mtext(if (layer == "temperature_change") "Temperature change (K)" else "Precipitation change (mm/day)", side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}

#' outline of the areas where the geography differs much from the climate model's
draw_climate_flag <- function(F, cfg) {
  if (is.null(F)) return(invisible())
  zz <- cbind(F$z, F$z[, 1]); ll <- c(F$lon, F$lon[1] + 360)
  cl <- contourLines(ll, F$lat, t(zz), levels = 0.5)
  if (length(cl)) draw_ll_list(lapply(cl, function(q) data.frame(lon = ((q$x + 180) %% 360) - 180, lat = q$y)), cfg,
                               col = "#3a3a3a", lwd = 0.9, lty = "33")
}

# ---- water balance and ocean layers ------------------------------------------------------
#' map layers of single fields: where they are drawn (land / sea), colours, units
SCALAR_LAYERS <- list(
  pme = list(field = "pme", title = "Precipitation \u2212 evaporation", type = "div", mask = "all",
             neg = c("#8c510a", "#d8b365", "#f6e8c3"), pos = c("#c7eae5", "#5ab4ac", "#01665e")),
  runoff = list(field = "runoff", title = "Total runoff (surface + subsurface)", type = "seq", mask = "land",
                stops = c(0, 25, 50, 100, 200, 400, 800, 1600, 3200), stops_day = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8),
                col = c("#fbf8ef", "#eaf2d7", "#c9e6c4", "#94cfbf", "#5cb2c6", "#3a8cc0", "#2b5fa8", "#24357f", "#1b1a4f")),
  salinity = list(field = "salinity", title = "Sea-surface salinity (5 m, HadCM3L)", type = "seq", mask = "sea", units = "psu",
                  stops = c(30, 32, 33, 34, 34.5, 35, 35.5, 36, 37, 38),
                  col = c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081", "#041c3f")),
  w_ocean = list(field = "w_ocean", title = "Ocean vertical velocity at 10 m (+ = upwelling)", type = "div", mask = "sea", units = "m/day",
                 neg = c("#2166ac", "#67a9cf", "#d1e5f0"), pos = c("#fddbc7", "#ef8a62", "#b2182b")),
  surface_density = list(field = "surface_density", title = "Sea-surface density \u03c3 (season; where deep water can form)", type = "seq", mask = "sea",
                         units = "kg/m\u00b3 - 1000", stops = c(20, 22, 23, 24, 25, 25.5, 26, 26.5, 27, 27.5, 28),
                         col = c("#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000", "#4a0000", "#200000")),
  o2_surface = list(field = "o2_surface", title = "Surface O2 at saturation", type = "seq", mask = "sea", units = "\u00b5mol/kg",
                    stops = c(100, 150, 175, 200, 225, 250, 275, 300, 350),
                    col = c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac")),
  freshwater = list(field = "freshwater", title = "Freshwater into the sea (P \u2212 E + rivers; + = freshening)", type = "div", mask = "sea",
                    neg = c("#8c510a", "#d8b365", "#f6e8c3"), pos = c("#d1e5f0", "#67a9cf", "#2166ac")),
  river_input = list(field = "river_input", title = "River water into the sea (spread ~300 km from the mouths)", type = "seq", mask = "sea",
                     stops = c(0, 10, 30, 100, 300, 1000, 3000, 10000, 30000),
                     col = c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")),
  salinity_source = list(field = "salinity_source", title = "Why the surface water is fresher / saltier", type = "cat", mask = "sea",
                         units = "", codes = 1:6,
                         labels = c("fresher: river water", "fresher: rain (P > E)", "fresher: currents / mixing",
                                    "near the ocean mean", "saltier: currents / mixing", "saltier: evaporation (E > P)"),
                         col = c("#08519c", "#6baed6", "#9e9ac8", "#eeeeee", "#f1b6da", "#b35806")),
  mld = list(field = "mld", title = "Ocean mixed-layer depth (HadCM3L)", type = "seq", mask = "sea", units = "m",
             stops = c(0, 20, 40, 60, 100, 150, 250, 400, 700),
             col = c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")))
.scalar_units <- function(L, season) L$units %||% (if (season == "ann") "mm/yr" else "mm/day")
.cat_colour <- function(v, L) {
  m <- t(col2rgb(L$col) / 255); i <- match(round(v), L$codes)
  out <- matrix(NA_real_, length(v), 3); ok <- !is.na(i); out[ok, ] <- m[i[ok], , drop = FALSE]; out
}
.scalar_factor <- function(L, season) if (is.null(L$units) && season == "ann") 365 else 1
.seq_colour <- function(v, L) {
  pos <- seq(0, 1, length.out = length(L$stops)); x <- pmin(pmax(v, min(L$stops)), max(L$stops))
  f <- approx(L$stops, pos, xout = x)$y; m <- col2rgb(L$col) / 255
  sapply(1:3, function(k) approx(pos, m[k, ], xout = f)$y)
}
#' the limit / stops of a layer for a field (auto for diverging layers)
scalar_scale <- function(layer, F, season, z = NULL) {
  L <- SCALAR_LAYERS[[layer]]; fac <- .scalar_factor(L, season)
  if (L$type == "cat") return(L)
  if (L$type == "seq") { if (is.null(L$units) && season != "ann") L$stops <- L$stops_day %||% (L$stops / 365); return(L) }
  q <- stats::quantile(abs(F$z * fac), 0.98, na.rm = TRUE)
  nice <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000)
  L$lim <- nice[which(nice >= q)[1]] %||% max(nice); L
}
render_scalar_rgb <- function(base, z, lat, v, L) {
  m <- if (L$type == "div") .diverging(as.vector(v), L$lim, neg = L$neg, pos = L$pos) else
       if (L$type == "cat") .cat_colour(as.vector(v), L) else .seq_colour(as.vector(v), L)
  sel <- switch(L$mask, land = z >= 0, sea = z < 0, all = rep(TRUE, length(z)))
  sel <- sel & is.finite(as.vector(v)); hs <- hillshade(z, lat); land <- z >= 0
  for (c in 1:3) { ch <- base[, , c]; mc <- m[, c]
    mc[land] <- pmin(pmax(mc[land] * (0.6 + 0.5 * hs[land]), 0), 1); ch[sel] <- mc[sel]; base[, , c] <- ch }
  base
}
colour_bar_scalar <- function(side = 1, layer, season = "ann", L) {
  if (L$type == "cat") {                              # classes: boxes with labels
    k <- length(L$codes); lab <- sprintf("%s, %s", L$title, CLIMATE_SEASONS[[season]])
    if (side %in% c(1, 3)) {
      plot.window(xlim = c(0, k), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      rect(0:(k - 1) + 0.04, 0, 1:k - 0.04, 1, col = L$col, border = "#9a9a9a", lwd = 0.4)
      mtext(sub(": ", ":\n", L$labels), side = side, at = 1:k - 0.5, line = if (side == 1) 0.9 else 0.3, cex = 0.45, padj = if (side == 1) 0 else 1)
      mtext(lab, side = side, line = 2.1, cex = 0.62, col = "#52514e")
    } else {
      plot.window(xlim = c(0, 1), ylim = c(0, k), xaxs = "i", yaxs = "i")
      rect(0, 0:(k - 1) + 0.06, 1, 1:k - 0.06, col = L$col, border = "#9a9a9a", lwd = 0.4)
      axis(side, at = 1:k - 0.5, labels = sub(": ", ":\n", L$labels), las = 1, cex.axis = 0.5, lwd = 0, lwd.ticks = 0, mgp = c(1, 0.3, 0))
    }
    return(invisible())
  }
  u <- .scalar_units(L, season)
  lab <- sprintf("%s, %s (%s)", L$title, CLIMATE_SEASONS[[season]], u)
  if (L$type == "div") { v <- seq(-L$lim, L$lim, length.out = 201); m <- .diverging(v, L$lim, neg = L$neg, pos = L$pos)
    at <- pretty(c(-L$lim, L$lim)); at <- at[abs(at) <= L$lim]; xl <- c(-L$lim, L$lim); labs <- ifelse(at > 0, paste0("+", at), at)
  } else { k <- length(L$stops) - 1; pb <- seq(0, k, length.out = 281); v <- pb
    vv <- approx(seq(0, k), L$stops, xout = pb)$y; m <- .seq_colour(vv, L); at <- 0:k; xl <- c(0, k)
    labs <- format(signif(L$stops, 3), trim = TRUE, drop0trailing = TRUE) }
  cc <- rgb(m[, 1], m[, 2], m[, 3]); n <- length(v)
  if (side %in% c(1, 3)) {
    plot.window(xlim = xl, ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(v[-n], 0, v[-1], 1, col = cc[-1], border = NA)
    axis(side, at = at, labels = labs, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    mtext(lab, side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = xl, xaxs = "i", yaxs = "i")
    rect(0, v[-n], 1, v[-1], col = cc[-1], border = NA)
    axis(side, at = at, labels = labs, las = 1, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    mtext(sprintf("%s (%s)", sub(" \\(.*", "", L$title), u), side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}
#' a field on a map grid, keeping gaps (NA) where the field has none (e.g. the ocean on land)
.scalar_on_grid <- function(F, lat, lon) {
  g <- expand.grid(lon = lon, lat = lat)
  matrix(clim_interp(F, g$lat, g$lon), length(lat), length(lon), byrow = TRUE)
}
