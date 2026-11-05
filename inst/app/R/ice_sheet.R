# =============================================================================
#  Ice sheets: build-up and melting under a climate run (Climate model tab, option 4).
#
#  Surface mass balance (positive-degree-day model):
#    - monthly near-surface temperature T and precipitation P of the climate run (HadCM3L +
#      SPEEDY response, as on the maps), brought to the ice-sheet grid and to the actual surface
#      height with a lapse rate (T) and less precipitation in colder air (P, ~7 % per K);
#    - snow = the part of the precipitation that falls below +1 C (normal distribution of the
#      daily temperature around the monthly mean, standard deviation sigma);
#    - positive degree days from the same distribution (Calov & Greve 2005);
#    - melt: 3 mm w.e. per degree day for snow, 8 mm for ice; up to 60 % of the snowfall
#      refreezes (Reeh 1991).
#  Ice flow: shallow-ice approximation (Glen's law n = 3) plus basal sliding (Weertman type,
#  u_b = C tau_b^3) where the bed is slippery: on soft marine sediments (bed below `soft_below`
#  m) under ice thicker than `thaw_thick` m (a thawed base), little sliding on hard rock; on an
#  equal-area grid
#  around each pole (Lambert azimuthal, default 40 km), explicit time steps; ice that would
#  float is removed (no ice shelves); the bed sinks under the ice (local isostasy, 3 kyr).
#  Coupling: every `couple_kyr` the climate model SPEEDY is run again with the ice sheet (its
#  height and albedo) and the change of the climate by the ice is added (asynchronous coupling).
#  Sea level: the ice volume lowers the sea (sea-level equivalent over the ocean area), so ice
#  can ground on the exposed shelves; maps can be drawn with the lowered sea level.
#  Equilibrium: ablation at the margins, calving where ice would float, less snowfall in the
#  colder and higher interior and the sinking bed limit the growth; the run stops early when
#  the volume changes less than eq_change (0.5 %) over eq_window_kyr (5 kyr).
#  Limitations: no ice shelves, no ice temperature (one flow factor A, the thawed base is set by
#  thickness), the ocean temperature does not respond to the ice in the coupling runs.
# =============================================================================

ICE_ALBEDO <- 0.65
RHO_I <- 910; RHO_W <- 1028; RHO_M <- 3300; G_ACC <- 9.81; R_EARTH <- 6371000
ICE_DEFAULTS <- list(duration_kyr = 50, couple_kyr = 10, res_km = 40, edge_lat = 30,
                     pdd_snow = 3, pdd_ice = 8, refreeze = 0.6, sigma = 5, lapse = 6.5, precip_k = 0.07,
                     snow_threshold = 1, A = 1e-16, tau_iso = 3000,
                     init = "none", init_lat = 60, init_thick = 2000, init_run = NA, save_kyr = 1,
                     bed = "mixed", slide_soft = 1e-13, slide_hard = 1e-15, thaw_thick = 500, soft_below = 200,
                     sea_level = TRUE, stop_equilibrium = TRUE, eq_window_kyr = 5, eq_change = 0.005,
                     max_sle = NA, g_min = 0.3, g_step = 0.02)

#' settings of an ice-sheet run (see ICE_DEFAULTS); init: "none" (ice-free start), "cap" (ice of
#' init_thick m on all land poleward of init_lat, in isostatic balance: for melting experiments)
#' or "run" (the end state of another ice run, init_run = its id)
ice_settings <- function(...) {
  s <- modifyList(ICE_DEFAULTS, list(...))
  s$duration_kyr <- max(1, s$duration_kyr); s$couple_kyr <- max(0, s$couple_kyr)
  s
}
.ice_tag <- function(st) {
  k <- st[c("duration_kyr", "couple_kyr", "res_km", "edge_lat", "pdd_snow", "pdd_ice", "refreeze", "sigma",
            "lapse", "precip_k", "snow_threshold", "A", "tau_iso", "init", "init_lat", "init_thick", "init_run",
            "bed", "slide_soft", "slide_hard", "thaw_thick", "soft_below", "sea_level", "stop_equilibrium")]
  if (is.finite(st$max_sle %||% NA)) k$max_sle <- c(st$max_sle, st$g_min %||% 0.3, st$g_step %||% 0.02)   # (old tags unchanged)
  tf <- tempfile(); on.exit(unlink(tf)); writeLines(deparse(k), tf)
  paste0("ice", substr(unname(tools::md5sum(tf)), 1, 8))
}
ice_dir <- function(ice_id, cfg = SILUCARB) file.path(cfg$dir_results, "speedy", "ice", ice_id)

# ---- grids -----------------------------------------------------------------------------
#' equal-area grid around a pole (Lambert azimuthal), reaching to edge_lat
ice_grid <- function(pole = c("south", "north"), res_km = 40, edge_lat = 35) {
  pole <- match.arg(pole); res <- res_km * 1000
  rmax <- 2 * R_EARTH * sin((90 - edge_lat) / 2 * pi / 180)
  n <- ceiling(rmax / res); x <- seq(-n, n) * res; y <- x
  X <- matrix(x, length(y), length(x), byrow = TRUE); Y <- matrix(y, length(y), length(x))
  r <- sqrt(X^2 + Y^2); cc <- 2 * asin(pmin(r / (2 * R_EARTH), 1)) * 180 / pi
  lat <- if (pole == "south") -(90 - cc) else 90 - cc
  lon <- (if (pole == "south") atan2(X, Y) else atan2(X, -Y)) * 180 / pi
  list(pole = pole, res = res, x = x, y = y, lat = lat, lon = lon, inside = r <= rmax, rmax = rmax,
       edge_lat = edge_lat)
}
.laea_xy <- function(lat, lon, pole) {
  cc <- if (pole == "south") 90 + lat else 90 - lat
  r <- 2 * R_EARTH * sin(cc / 2 * pi / 180); l <- lon * pi / 180
  list(x = r * sin(l), y = if (pole == "south") r * cos(l) else -r * cos(l))
}
#' bilinear value of a grid matrix (rows = y, cols = x) at points; NA outside
.grid_bilin <- function(M, gx, gy, qx, qy) {
  fx <- (qx - gx[1]) / (gx[2] - gx[1]) + 1; fy <- (qy - gy[1]) / (gy[2] - gy[1]) + 1
  ok <- fx >= 1 & fx <= length(gx) & fy >= 1 & fy <= length(gy)
  j <- pmin(floor(fx), length(gx) - 1); i <- pmin(floor(fy), length(gy) - 1); tx <- fx - j; ty <- fy - i
  j[!ok] <- 1; i[!ok] <- 1
  v <- (1 - ty) * ((1 - tx) * M[cbind(i, j)] + tx * M[cbind(i, j + 1)]) + ty * ((1 - tx) * M[cbind(i + 1, j)] + tx * M[cbind(i + 1, j + 1)])
  v[!ok] <- NA; v
}
#' a field list(z, lat, lon) sampled on an ice grid
.on_grid <- function(F, g) matrix(clim_interp(F, as.vector(g$lat), as.vector(g$lon)), nrow(g$lat), ncol(g$lat))

#' SPEEDY input-oriented matrix (96 lon x 48 lat north->south) as a field list(z, lat, lon)
.speedy_in_field <- function(M) {
  lon <- ((SPEEDY_LON + 180) %% 360) - 180; o <- order(lon)
  list(z = t(M[, 48:1])[, o], lat = rev(speedy_lats), lon = lon[o])
}
#' SPEEDY output-oriented matrix (96 x 48 south->north) as a field
.speedy_out_field <- function(M) {
  lon <- ((SPEEDY_LON + 180) %% 360) - 180; o <- order(lon)
  list(z = t(M)[, o], lat = rev(speedy_lats), lon = lon[o])
}

# ---- climate for the ice sheet ------------------------------------------------------------
#' monthly near-surface temperature (C) and precipitation (mm/day) of a climate run on the
#' HadCM3L grid (as on the maps): HadCM3L + SPEEDY response
climate_run_monthly <- function(F, cfg = SILUCARB) {
  age <- F$setup$map$age; mf <- climate_month_files(age, cfg)
  resp <- function(k, m) .speedy_out_field(F$map_effect[, , m, k] + F$change[, , m, k])
  lapply(1:12, function(m) {
    pr <- .read_clim_var(mf[m], "precip_mm_srf"); pr$z <- pr$z * 86400
    tt <- .read_clim_var(mf[m], "temp_mm_1_5m"); tt$z <- tt$z - 273.15
    dp <- resp("precip", m); dt <- resp("t0", m)
    g <- expand.grid(lon = tt$lon, lat = tt$lat)
    add <- function(Fh, D) { Fh$z <- Fh$z + matrix(clim_interp(D, g$lat, g$lon), length(tt$lat), length(tt$lon), byrow = TRUE); Fh }
    pr2 <- add(pr, dp); pr2$z <- pmax(pr2$z, 0)
    list(t = add(tt, dt), p = pr2)
  })
}

#' surface mass balance (m ice per year) on a grid: T, P lists of 12 matrices (C, mm/day)
ice_smb <- function(Tm, Pm, st) {
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31); s <- st$sigma
  acc <- 0; pdd <- 0
  for (m in 1:12) {
    T <- Tm[[m]]
    acc <- acc + Pm[[m]] * days[m] * pnorm((st$snow_threshold - T) / s)
    pdd <- pdd + days[m] * (s / sqrt(2 * pi) * exp(-T^2 / (2 * s^2)) + T * pnorm(T / s))
  }
  ms <- pmin(st$pdd_snow * pdd, acc)                     # snow melt (mm w.e.)
  mi <- st$pdd_ice * pmax(pdd - ms / st$pdd_snow, 0)     # then ice melt
  refr <- pmin(st$refreeze * acc, ms)
  smb <- acc - (ms + mi - refr)                          # mm w.e. / yr
  list(smb = smb / 1000 * 1000 / RHO_I, acc = acc / 1000, melt = (ms + mi - refr) / 1000, pdd = pdd)
}

#' monthly T and P on an ice grid at the surface height s (m), from the climate (cl), the
#' reference height of the climate (href, m, on the grid) and an extra change by the ice (fb)
.ice_forcing <- function(clg, s, href, st) {
  dz <- pmax(s, 0) - href
  Tm <- lapply(1:12, function(m) clg$t[[m]] - st$lapse / 1000 * dz)
  Pm <- lapply(1:12, function(m) clg$p[[m]] * pmin(pmax(exp(-st$precip_k * st$lapse / 1000 * dz), 0.1), 2))
  list(t = Tm, p = Pm)
}

# ---- ice flow ------------------------------------------------------------------------------
#' shallow-ice flow, explicit, `years` long; returns the new state of one domain
ice_flow <- function(D, years, st, smb_fun, smb_every = 100, sl = 0) {
  H <- D$H; b0 <- D$b0; w <- D$w; dx <- D$g$res; inside <- D$g$inside
  n <- 3; Gam <- 2 * st$A * (RHO_I * G_ACC)^n / (n + 2)
  ny <- nrow(H); nx <- ncol(H); t <- 0; t_smb <- -Inf; smb <- D$smb
  cb <- D$cb %||% (0 * H)                                # sliding coefficient (m / yr / Pa^3)
  thaw_h <- st$thaw_thick %||% 500
  sh <- function(M, di, dj) { i <- pmin(pmax(seq_len(ny) + di, 1), ny); j <- pmin(pmax(seq_len(nx) + dj, 1), nx); M[i, j] }
  G3 <- (RHO_I * G_ACC)^3
  while (t < years - 1e-9) {
    b <- b0 - w; s <- b + H
    if (t - t_smb >= smb_every) { smb <- smb_fun(s - sl); t_smb <- t }   # height above the (lowered) sea
    cbx <- cb * (H > thaw_h)                            # sliding only on a thawed base
    # staggered diffusivities (x faces between j and j+1, y faces between i and i+1)
    sE <- sh(s, 0, 1); HE <- sh(H, 0, 1)
    dsdx <- (sE - s) / dx
    dsdy_x <- ((sh(s, 1, 0) - sh(s, -1, 0)) + (sh(sE, 1, 0) - sh(sE, -1, 0))) / (4 * dx)
    Hx <- (H + HE) / 2; g2 <- dsdx^2 + dsdy_x^2
    Dx <- Gam * Hx^(n + 2) * g2^((n - 1) / 2) + (cbx + sh(cbx, 0, 1)) / 2 * G3 * Hx^4 * g2
    Dx[(dsdx > 0 & HE <= 0) | (dsdx < 0 & H <= 0)] <- 0   # no flow out of an ice-free cell
    sN <- sh(s, 1, 0); HN <- sh(H, 1, 0)
    dsdy <- (sN - s) / dx
    dsdx_y <- ((sh(s, 0, 1) - sh(s, 0, -1)) + (sh(sN, 0, 1) - sh(sN, 0, -1))) / (4 * dx)
    Hy <- (H + HN) / 2; g2 <- dsdy^2 + dsdx_y^2
    Dy <- Gam * Hy^(n + 2) * g2^((n - 1) / 2) + (cbx + sh(cbx, 1, 0)) / 2 * G3 * Hy^4 * g2
    Dy[(dsdy > 0 & HN <= 0) | (dsdy < 0 & H <= 0)] <- 0
    Dmax <- max(Dx, Dy, 1)
    dt <- min(10, 0.2 * dx^2 / Dmax, years - t)
    qx <- -Dx * dsdx; qy <- -Dy * dsdy                   # flux through the east / north face
    qx[, nx] <- 0; qy[ny, ] <- 0
    div <- (qx - sh(qx, 0, -1) * (col(H) > 1)) / dx + (qy - sh(qy, -1, 0) * (row(H) > 1)) / dx
    H <- pmax(H + dt * (smb - div), 0)
    float <- b < sl & H * RHO_I < (sl - b) * RHO_W       # would float at the (lowered) sea level: calves
    H[float | !inside] <- 0
    w <- w + dt * (RHO_I / RHO_M * H - w) / st$tau_iso   # bed sinks under the ice
    t <- t + dt
  }
  D$H <- H; D$w <- w; D$smb <- smb; D
}
#' sliding coefficient of the bed: "rough" (hard bed, little sliding), "slippery" (soft, wet bed
#' everywhere) or "mixed" (slippery on marine sediments, i.e. where the bed lies below
#' `soft_below` m, rough on higher bedrock)
.bed_slide <- function(b0, st) {
  switch(st$bed %||% "mixed",
         rough = 0 * b0 + st$slide_hard, slippery = 0 * b0 + st$slide_soft,
         mixed = ifelse(b0 < st$soft_below, st$slide_soft, st$slide_hard))
}

# ---- run state ----------------------------------------------------------------------------
#' the domains (south / north) that have land poleward of the edge latitude
.ice_domains <- function(dem, st) {
  out <- list()
  for (pole in c("south", "north")) {
    sel <- if (pole == "south") dem$lat <= -st$edge_lat else dem$lat >= st$edge_lat
    if (sum(dem$z[sel, ] > 0) > 50) out[[pole]] <- ice_grid(pole, st$res_km, st$edge_lat)
  }
  out
}
.ocean_area <- function(dem, sl = 0) {
  w <- cos(dem$lat * pi / 180); sum(w * rowSums(dem$z < sl)) / (sum(w) * ncol(dem$z)) * 4 * pi * R_EARTH^2
}
#' sea level (m, negative = lowered) for the ice volume of a state: the water locked in the ice
#' sheets, spread over the ocean (whose area shrinks as the sea falls); no water-load isostasy
.sea_level_of <- function(state, dem) {
  V <- sum(vapply(state$D, function(D) sum(D$H) * D$g$res^2, 0)) * RHO_I / RHO_W
  if (V <= 0) return(0)
  sl <- -V / .ocean_area(dem)
  for (k in 1:3) sl <- -V / .ocean_area(dem, sl)       # smaller ocean when lower
  sl
}
#' initial state of an ice run
ice_init <- function(map, st, cfg = SILUCARB) {
  dem <- read_dem_nc(map$file); F0 <- list(z = dem$z, lat = dem$lat, lon = dem$lon)
  doms <- .ice_domains(dem, st)
  if (!length(doms)) stop("no land poleward of ", st$edge_lat, " degrees: nothing for ice sheets", call. = FALSE)
  D <- lapply(doms, function(g) {
    b0 <- .on_grid(F0, g); b0[!g$inside] <- -4000
    # work only on the part of the grid with land (and shallow sea around it), plus a margin
    keep <- which(b0 > -1000 & g$inside, arr.ind = TRUE)
    ri <- max(1, min(keep[, 1]) - 3):min(nrow(b0), max(keep[, 1]) + 3)
    ci <- max(1, min(keep[, 2]) - 3):min(ncol(b0), max(keep[, 2]) + 3)
    g$x <- g$x[ci]; g$y <- g$y[ri]; g$lat <- g$lat[ri, ci, drop = FALSE]; g$lon <- g$lon[ri, ci, drop = FALSE]
    g$inside <- g$inside[ri, ci, drop = FALSE]; g$inside[c(1, nrow(g$inside)), ] <- FALSE; g$inside[, c(1, ncol(g$inside))] <- FALSE
    b0 <- b0[ri, ci, drop = FALSE]
    H <- matrix(0, nrow(b0), ncol(b0))
    if (st$init == "cap") H[b0 > 0 & abs(g$lat) >= st$init_lat & g$inside] <- st$init_thick
    list(g = g, b0 = b0, H = H, w = RHO_I / RHO_M * H, smb = 0 * H, cb = .bed_slide(b0, st))
  })
  if (st$init == "run" && !is.na(st$init_run)) {
    prev <- tryCatch(readRDS(file.path(ice_dir(st$init_run, cfg), "state.rds")), error = function(e) NULL)
    if (is.null(prev)) stop("the ice run to start from has no saved state", call. = FALSE)
    for (p in intersect(names(D), names(prev$D))) {
      P <- prev$D[[p]]; g <- D[[p]]$g
      xy <- list(x = rep(g$x, each = length(g$y)), y = rep(g$y, length(g$x)))
      H <- .grid_bilin(P$H, P$g$x, P$g$y, xy$x, xy$y); H[!is.finite(H)] <- 0
      Wp <- .grid_bilin(P$w, P$g$x, P$g$y, xy$x, xy$y); Wp[!is.finite(Wp)] <- 0
      D[[p]]$H <- matrix(H, length(g$y)); D[[p]]$w <- matrix(Wp, length(g$y))
    }
  }
  st0 <- list(t = 0, k = 0, phase = "ice", D = D, ocean_area = .ocean_area(dem), fb = NULL, href = NULL, job = NULL, sl = 0,
              dem_z = NULL)
  if (isTRUE(st$sea_level)) st0$sl <- .sea_level_of(st0, dem)
  st0
}

#' the ice thickness and bed sinking of the state on the map grid (for SPEEDY and the maps)
ice_on_latlon <- function(state, lat, lon) {
  LAT <- matrix(lat, length(lat), length(lon)); LON <- matrix(lon, length(lat), length(lon), byrow = TRUE)
  H <- matrix(0, length(lat), length(lon)); W <- H; SMB <- matrix(NA_real_, length(lat), length(lon))
  for (D in state$D) {
    g <- D$g; sel <- if (g$pole == "south") LAT <= -g$edge_lat + 0.5 else LAT >= g$edge_lat - 0.5
    if (!any(sel)) next
    xy <- .laea_xy(LAT[sel], LON[sel], g$pole)
    h <- .grid_bilin(D$H, g$x, g$y, xy$x, xy$y); ww <- .grid_bilin(D$w, g$x, g$y, xy$x, xy$y)
    sm <- .grid_bilin(D$smb, g$x, g$y, xy$x, xy$y)
    H[sel] <- ifelse(is.finite(h), h, 0); W[sel] <- ifelse(is.finite(ww), ww, 0); SMB[sel] <- sm
  }
  list(H = H, w = W, smb = SMB)
}

#' totals of a state: volume (km3), area (km2), sea-level equivalent (m)
ice_totals <- function(state) {
  out <- do.call(rbind, lapply(names(state$D), function(p) {
    D <- state$D[[p]]; a <- D$g$res^2; ice <- D$H > 1
    V <- sum(D$H[ice]) * a
    data.frame(t_kyr = state$t / 1000, domain = p, volume_km3 = V / 1e9, mass_Gt = V * RHO_I / 1e12, area_km2 = sum(ice) * a / 1e6,
               sle_m = V * RHO_I / RHO_W / state$ocean_area, max_thick_m = max(D$H), sea_level_m = state$sl %||% 0)
  }))
  if (!is.null(state$g_ice)) out$gradient <- state$g_ice      # runs with an ice limit
  out
}
#' the climate of the run on the ice grids, with the SPEEDY change by the ice (state$fb)
.ice_climate_on_grids <- function(state, clm, href0) {
  lapply(state$D, function(D) {
    tt <- lapply(1:12, function(m) { v <- .on_grid(clm[[m]]$t, D$g)
      if (!is.null(state$fb)) v <- v + .on_grid(state$fb$t[[m]], D$g); v })
    if (!is.null(state$g_ice) && abs(state$g_ice - (state$g_run %||% 1)) > 1e-6)   # ice limit: extra gradient change
      tt <- lapply(1:12, function(m) {
        a <- gradient_anomaly(state$age, state$g_ice, GRID1$lat, m) - gradient_anomaly(state$age, state$g_run %||% 1, GRID1$lat, m)
        tt[[m]] + .on_grid(list(z = matrix(a, 180, 360), lat = GRID1$lat, lon = GRID1$lon), D$g) })
    pp <- lapply(1:12, function(m) { v <- .on_grid(clm[[m]]$p, D$g)
      if (!is.null(state$fb)) v <- v + .on_grid(state$fb$p[[m]], D$g); pmax(v, 0) })
    list(t = tt, p = pp, href = .on_grid(state$href %||% href0, D$g))
  })
}

#' run the ice sheets for `years` under the current climate (state$fb: change by the ice)
ice_segment <- function(state, clm, href0, st, years, ts_file = NULL, snap_dir = NULL, dem = NULL) {
  clg <- .ice_climate_on_grids(state, clm, href0)
  step <- min(years, st$save_kyr * 1000); done <- 0
  if (!is.null(ts_file) && !file.exists(ts_file))      # the start
    utils::write.table(ice_totals(state), ts_file, sep = ",", row.names = FALSE)
  while (done < years - 1e-6) {
    dt <- min(step, years - done)
    .heartbeat()                                        # runner lock (and this run's busy lock) stay fresh
    for (p in names(state$D)) {
      cg <- clg[[p]]
      smb_fun <- function(s) { f <- .ice_forcing(cg, s, cg$href, st); ice_smb(f$t, f$p, st)$smb }
      state$D[[p]] <- ice_flow(state$D[[p]], dt, st, smb_fun, sl = if (isTRUE(st$sea_level)) state$sl %||% 0 else 0)
    }
    done <- done + dt; state$t <- state$t + dt
    if (isTRUE(st$sea_level) && !is.null(dem)) state$sl <- .sea_level_of(state, dem)
    tot <- ice_totals(state)
    if (!is.null(ts_file)) utils::write.table(tot, ts_file, sep = ",", row.names = FALSE,
                                              append = file.exists(ts_file), col.names = !file.exists(ts_file))
    if (!is.null(snap_dir)) {
      gf <- file.path(dirname(snap_dir), "grids.rds")
      if (!file.exists(gf)) saveRDS(list(grids = lapply(state$D, `[[`, "g"), ocean_area = state$ocean_area), gf)
      snap <- lapply(state$D, function(D) list(H = round(D$H), w = round(D$w), smb = round(D$smb, 3)))
      attr(snap, "sea_level") <- state$sl %||% 0
      saveRDS(snap, file.path(snap_dir, sprintf("ice_%07.1fkyr.rds", state$t / 1000)))
    }
    message(sprintf("  ice %.1f kyr: %s; sea level %.1f m", state$t / 1000,
                    paste(sprintf("%s %.3g km3 (%.2f m SLE)", tot$domain, tot$volume_km3, tot$sle_m), collapse = ", "), state$sl %||% 0))
    state$vol_hist <- rbind(state$vol_hist, data.frame(t = state$t, V = sum(tot$volume_km3)))
    # ice limit: above max_sle the equator-to-pole gradient of the ice model's climate is lowered step by step
    # ice limit: the factor follows the ice volume AND its trend, so it does not keep lowering the factor
    # while the ice is already shrinking (that over-correction melted whole ice sheets):
    #   above the limit and still growing      -> lower the factor by a step (more steps the further above)
    #   above the limit but shrinking          -> keep it
    #   below the limit and shrinking          -> raise it by a step (stop the melt near the limit)
    #   below half the limit and not growing   -> raise it by a step (back towards the run's own factor)
    if (is.finite(st$max_sle %||% NA) && !is.null(state$g_ice)) {
      sle <- sum(tot$sle_m); g0 <- state$g_ice; ds <- sle - (state$sle_prev %||% sle); state$sle_prev <- sle
      stp <- st$g_step %||% 0.02; lim <- st$max_sle
      if (sle > lim && ds > 0) state$g_ice <- g0 - stp * min(3, 1 + (sle - lim) / max(lim, 1))
      else if (lim > 0 && sle <= lim && ds < -0.01) state$g_ice <- g0 + stp
      else if (lim > 0 && sle < 0.5 * lim && ds <= 0.01) state$g_ice <- g0 + stp
      state$g_ice <- min(state$g_run %||% 1, max(st$g_min %||% 0.3, state$g_ice))
      if (abs(state$g_ice - g0) > 1e-9) {
        clg <- .ice_climate_on_grids(state, clm, href0)
        state$vol_hist <- tail(state$vol_hist, 1)        # equilibrium is judged again from here
        message(sprintf("  ice %.2f m SLE (limit %.1f, %+.2f m/kyr): gradient factor of the ice climate now %.2f", sle, lim, ds, state$g_ice))
      }
    }
    # equilibrium (with climate updates: judged over at least one update interval)
    ste <- st; ste$eq_window_kyr <- max(st$eq_window_kyr, st$couple_kyr)
    if (isTRUE(st$stop_equilibrium) && .ice_equilibrium(state$vol_hist, ste)) { state$equilibrium <- state$t; break }
  }
  state
}
#' equilibrium: the ice volume changed less than eq_change (relative) over the last eq_window_kyr
.ice_equilibrium <- function(h, st) {
  if (is.null(h) || nrow(h) < 2) return(FALSE)
  t1 <- max(h$t); old <- h[h$t <= t1 - st$eq_window_kyr * 1000 + 1e-6, ]
  if (!nrow(old)) return(FALSE)
  V0 <- old$V[nrow(old)]; V1 <- h$V[nrow(h)]
  if (V1 < 1 && V0 < 1) return(TRUE)                        # stays ice-free
  abs(V1 - V0) / max(V1, V0) < st$eq_change
}

# ---- the coupled run (one step at a time, for the background runner) -----------------------
#' the SPEEDY runs of the ice coupling: "iceref" (no ice) and "ice<k>" (with the ice sheet of
#' coupling step k), with the experiment's orbit / CO2 / Sun and prescribed sea surface
.ice_job <- function(map, x, ice_id, k, geo, cfg) {
  age <- map$age; xr <- .resolve(x, age); L <- SPEEDY_LENGTHS[[x$length %||% "standard"]]$reference
  tag <- sprintf("%s_%s", ice_id, if (k == 0) "iceref" else paste0("ice", k))     # ice_id starts with the map tag
  list(tag = tag, dem = map$file, sea_mode = 1, x = xr, spinup = L[["spinup"]], years = L[["years"]], shift = xr$shift,
       geo = geo, age = age, dir = speedy_dir(age, tag, cfg))
}

#' one step of an ice run (map, experiment, ice settings): runs the ice model up to the next
#' coupling time, or starts / checks the SPEEDY runs of the coupling. Returns the state.
ice_continue_one <- function(row, cfg = SILUCARB) {
  s <- run_setup(row, cfg); ice_id <- row$ice_id; dirI <- ice_dir(ice_id, cfg)
  # only one process works on an ice run at a time (busy.lock, refreshed every model kyr)
  bl <- file.path(dirI, "busy.lock")
  if (file.exists(bl) && difftime(Sys.time(), file.info(bl)$mtime, units = "mins") < RUNNER_STALE_MIN) {
    who <- suppressWarnings(as.integer(tryCatch(readLines(bl, warn = FALSE)[2], error = function(e) NA)))
    if (!identical(who, Sys.getpid())) return(invisible(NULL))
  }
  writeLines(c(format(Sys.time()), Sys.getpid()), bl)
  on.exit(unlink(bl), add = TRUE)
  hb0 <- getOption("silucarb.heartbeat")
  options(silucarb.heartbeat = function() { if (is.function(hb0)) hb0(); writeLines(c(format(Sys.time()), Sys.getpid()), bl) })
  on.exit(options(silucarb.heartbeat = hb0), add = TRUE)
  st <- modifyList(ICE_DEFAULTS, readRDS(file.path(dirI, "settings.rds")))
  jobs <- speedy_jobs(s$map, s$x, cfg)
  if (any(vapply(jobs, .job_status, "") != "done")) return(invisible(NULL))   # climate run first
  fs <- file.path(dirI, "state.rds")
  state <- if (file.exists(fs)) readRDS(fs) else ice_init(s$map, st, cfg)
  if (state$phase %in% c("done", "failed")) return(invisible(state))
  if (is.finite(st$max_sle %||% NA) && is.null(state$g_ice)) {    # ice limit: starts at the run's gradient factor
    state$g_run <- s$x$gradient %||% 1; state$g_ice <- state$g_run; state$age <- s$map$age }
  dir.create(file.path(dirI, "snapshots"), showWarnings = FALSE, recursive = TRUE)
  geo0 <- speedy_geography(s$map$file); href0 <- .speedy_in_field(geo0$orog)
  if (state$phase == "waiting") {
    jr <- .ice_job(s$map, s$x, ice_id, 0, NULL, cfg); jk <- .ice_job(s$map, s$x, ice_id, state$k, NULL, cfg)
    sr <- .job_status(jr); sk <- .job_status(jk)
    if ("failed" %in% c(sr, sk)) { state$phase <- "failed"; saveRDS(state, fs); return(invisible(state)) }
    if (sr == "stopped" || sk == "stopped") { state$phase <- "couple"; saveRDS(state, fs); return(invisible(state)) }
    if (sr != "done" || sk != "done") return(invisible(state))
    a <- .read_means(jr$dir); b <- .read_means(jk$dir)
    state$fb <- list(t = lapply(1:12, function(m) .speedy_out_field(smooth_sp(b[, , m, "t0"] - a[, , m, "t0"]))),
                     p = lapply(1:12, function(m) .speedy_out_field(smooth_sp(b[, , m, "precip"] - a[, , m, "precip"]))))
    state$href <- .speedy_in_field(matrix(readBin(file.path(jk$dir, "orog.bin"), "numeric", size = 4, n = 96 * 48, endian = "little"), 96, 48))
    state$phase <- "ice"
  }
  if (state$phase == "couple") {                       # start the SPEEDY runs with the ice sheet
    if (.n_running(cfg) + 2 > cfg$workers && .n_running(cfg) > 0) return(invisible(state))
    dem <- read_dem_nc(s$map$file); ic <- ice_on_latlon(state, dem$lat, dem$lon)
    geo <- speedy_geography(list(z = dem$z - ic$w - (state$sl %||% 0), lat = dem$lat, lon = dem$lon), ice = ic$H)   # heights above the lowered sea
    forcing <- speedy_forcing(s$map$age, cfg)
    jr <- .ice_job(s$map, s$x, ice_id, 0, geo0, cfg)
    li <- land_ice_flags_forcing(s$map$age, jr$shift, s$map$file, cfg)$flags        # sea ice <-> land ice
    if (!.job_status(jr) %in% c("done", "running")) {
      fr <- shift_forcing(forcing, jr$shift); fr$icec <- cap_perennial(fr$icec, speedy_lats, li)
      speedy_launch(jr, fr, NULL, cfg)
    }
    jk <- .ice_job(s$map, s$x, ice_id, state$k, geo, cfg)
    ia <- vapply(c("south", "north"), function(p) { D <- state$D[[p]]; if (is.null(D)) 0 else sum(D$H > 10) * D$g$res^2 / 1e6 }, 0)
    lk <- li | land_ice_flags(ia)                      # an ice sheet allows permanent sea ice
    fk <- shift_forcing(forcing, jk$shift); fk$icec <- cap_perennial(fk$icec, speedy_lats, lk)
    speedy_launch(jk, fk, NULL, cfg)
    state$phase <- "waiting"; saveRDS(state, fs)
    return(invisible(state))
  }
  # phase "ice": run the ice sheets to the next coupling time (or the end)
  clm <- climate_run_monthly(climate_run_fields(row, cfg), cfg)
  dur <- st$duration_kyr * 1000
  nxt <- if (st$couple_kyr > 0) min(dur, (floor(state$t / (st$couple_kyr * 1000) + 1e-9) + 1) * st$couple_kyr * 1000) else dur
  if (st$couple_kyr == 0) nxt <- min(dur, state$t + 5000)   # no coupling: save every 5 kyr of model time
  dem <- read_dem_nc(s$map$file)
  state <- ice_segment(state, clm, href0, st, nxt - state$t, file.path(dirI, "timeseries.csv"), file.path(dirI, "snapshots"), dem)
  if (state$t >= dur - 1e-6 || !is.null(state$equilibrium)) state$phase <- "done"
  else if (st$couple_kyr > 0) { state$k <- state$k + 1; state$phase <- "couple" }
  saveRDS(state, fs)
  invisible(state)
}
#' light smoothing of a SPEEDY output field (noise of short runs)
smooth_sp <- function(M) { f <- .speedy_out_field(M); z <- smooth_km(f$z, f$lat, 250, nbands = 24); lon <- ((SPEEDY_LON + 180) %% 360) - 180
  o <- order(lon); back <- matrix(0, 96, 48); back[o, ] <- t(z); back }

#' state of an ice run: list(status, progress, t_kyr)
ice_status <- function(ice_id, cfg = SILUCARB) {
  dirI <- ice_dir(ice_id, cfg); st <- tryCatch(modifyList(ICE_DEFAULTS, readRDS(file.path(dirI, "settings.rds"))), error = function(e) NULL)
  fs <- file.path(dirI, "state.rds")
  if (is.null(st)) return(list(status = "missing", progress = 0, t_kyr = 0))
  if (!file.exists(fs)) return(list(status = "waiting", progress = 0, t_kyr = 0))
  s <- readRDS(fs)
  list(status = switch(s$phase, done = if (!is.null(s$equilibrium)) "done (equilibrium)" else "done", failed = "failed", "running"),
       progress = if (s$phase == "done") 100 else round(100 * s$t / (st$duration_kyr * 1000)), t_kyr = s$t / 1000,
       gradient = s$g_ice)
}
#' time series of an ice run (volume, area, sea-level equivalent per domain)
ice_timeseries <- function(ice_id, cfg = SILUCARB) {
  f <- file.path(ice_dir(ice_id, cfg), "timeseries.csv")
  if (!file.exists(f)) return(NULL)
  ts <- utils::read.csv(f)
  ts <- ts[!duplicated(ts[c("t_kyr", "domain")], fromLast = TRUE), , drop = FALSE]   # a step written twice: keep the last
  ts[order(ts$t_kyr, ts$domain), , drop = FALSE]
}
#' the snapshot times (kyr) of an ice run
ice_snapshots <- function(ice_id, cfg = SILUCARB) {
  f <- list.files(file.path(ice_dir(ice_id, cfg), "snapshots"), "^ice_.*kyr\\.rds$")
  sort(as.numeric(sub("^ice_([0-9.]+)kyr\\.rds$", "\\1", f)))
}
#' an ice snapshot (kyr; NULL = the last) as a state-like list for ice_on_latlon()
ice_snapshot <- function(ice_id, t_kyr = NULL, cfg = SILUCARB) {
  tt <- ice_snapshots(ice_id, cfg); if (!length(tt)) return(NULL)
  t_kyr <- if (is.null(t_kyr) || !isTRUE(is.finite(t_kyr))) max(tt) else tt[which.min(abs(tt - t_kyr))]
  snap <- readRDS(file.path(ice_dir(ice_id, cfg), "snapshots", sprintf("ice_%07.1fkyr.rds", t_kyr)))
  st <- readRDS(file.path(ice_dir(ice_id, cfg), "settings.rds"))
  gr <- readRDS(file.path(ice_dir(ice_id, cfg), "grids.rds"))
  D <- lapply(names(snap), function(p) c(list(g = gr$grids[[p]]), snap[[p]])); names(D) <- names(snap)
  list(D = D, t = t_kyr * 1000, ocean_area = gr$ocean_area, sl = attr(snap, "sea_level") %||% 0)
}

#' plot of the ice volume (sea-level equivalent) and area through time
plot_ice_timeseries <- function(ice_id, cfg = SILUCARB, main = NULL) {
  ts <- ice_timeseries(ice_id, cfg); if (is.null(ts)) { plot.new(); text(0.5, 0.5, "no ice results yet"); return(invisible()) }
  op <- par(mfrow = c(1, 2), mar = c(3.4, 4.2, 2.2, 3.2), mgp = c(2.1, 0.6, 0), cex = 0.85); on.exit(par(op))
  cols <- c(south = "#2166ac", north = "#b2182b")
  tot <- aggregate(cbind(sle_m, area_km2) ~ t_kyr, ts, sum)
  eq <- tryCatch(readRDS(file.path(ice_dir(ice_id, cfg), "state.rds"))$equilibrium, error = function(e) NULL)
  plot(tot$t_kyr, tot$sle_m, type = "l", lwd = 2, xlab = "Time (kyr)", ylab = "Ice volume (m sea-level equiv.)",
       main = main %||% "Ice volume", font.main = 1, las = 1, ylim = range(0, tot$sle_m, na.rm = TRUE))
  for (p in unique(ts$domain)) lines(ts$t_kyr[ts$domain == p], ts$sle_m[ts$domain == p], col = cols[p], lty = 2)
  legend("topleft", bty = "n", cex = 0.8, legend = c("total", paste(unique(ts$domain), "pole")), lty = c(1, rep(2, length(unique(ts$domain)))),
         col = c("black", cols[unique(ts$domain)]), lwd = c(2, 1, 1))
  if (!is.null(eq)) { abline(v = eq / 1000, col = "#2b8a3e", lty = 3); text(eq / 1000, par("usr")[3] + 0.05 * diff(par("usr")[3:4]), sprintf("equilibrium at %.0f kyr", eq / 1000), adj = c(1.05, 0), cex = 0.7, col = "#2b8a3e") }
  st <- tryCatch(modifyList(ICE_DEFAULTS, readRDS(file.path(ice_dir(ice_id, cfg), "settings.rds"))), error = function(e) ICE_DEFAULTS)
  if (is.finite(st$max_sle %||% NA)) {
    abline(h = st$max_sle, col = "#d7301f", lty = 2, lwd = 1.5)
    if (!is.null(ts$gradient)) { g <- aggregate(gradient ~ t_kyr, ts, function(v) v[1])
      mtext(sprintf("limit %.1f m SLE (red dashed); gradient factor of the ice climate (orange, right axis): start %.2f, lowest %.2f, end %.2f",
                    st$max_sle, g$gradient[1], min(g$gradient), tail(g$gradient, 1)), side = 3, line = 0.1, cex = 0.6, col = "#d7301f")
      u <- par("usr"); gl <- c(min(0.3, g$gradient), 1)
      lines(g$t_kyr, u[3] + (g$gradient - gl[1]) / diff(gl) * (u[4] - u[3]), col = "#e6550d", lwd = 1.3)
      at <- pretty(gl); at <- at[at >= gl[1] & at <= gl[2]]
      axis(4, at = u[3] + (at - gl[1]) / diff(gl) * (u[4] - u[3]), labels = at, col.axis = "#e6550d", col = "#e6550d", las = 1, cex.axis = 0.75) }
  }
  if (!is.null(ts$sea_level_m)) { sl <- aggregate(sea_level_m ~ t_kyr, ts, function(v) v[1])
    mtext(sprintf("sea level %.1f m at the end", tail(sl$sea_level_m, 1)), side = 1, line = 2.2, adj = 1, cex = 0.7, col = "#555555") }
  plot(tot$t_kyr, tot$area_km2 / 1e6, type = "l", lwd = 2, xlab = "Time (kyr)", ylab = "Ice area (million km2)",
       main = "Ice area", font.main = 1, las = 1, ylim = range(0, tot$area_km2 / 1e6, na.rm = TRUE))
  for (p in unique(ts$domain)) lines(ts$t_kyr[ts$domain == p], ts$area_km2[ts$domain == p] / 1e6, col = cols[p], lty = 2)
}

#' export an ice run as NetCDF (each snapshot: thickness, surface mass balance, bed sinking;
#' equal-area polar grids with their latitude / longitude)
export_ice_run <- function(ice_id, file, cfg = SILUCARB) {
  tt <- ice_snapshots(ice_id, cfg); if (!length(tt)) stop("no ice results yet", call. = FALSE)
  st <- readRDS(file.path(ice_dir(ice_id, cfg), "settings.rds"))
  first <- ice_snapshot(ice_id, tt[1], cfg)
  defs <- list(); dims <- list()
  for (p in names(first$D)) {
    g <- first$D[[p]]$g
    dx <- ncdf4::ncdim_def(paste0("x_", p), "m", g$x); dy <- ncdf4::ncdim_def(paste0("y_", p), "m", g$y)
    dt <- ncdf4::ncdim_def("time", "kyr since the start", tt)
    dims[[p]] <- list(dx, dy, dt)
    defs[[paste0("lat_", p)]] <- ncdf4::ncvar_def(paste0("lat_", p), "degrees_north", list(dx, dy), 1e20, "palaeo-latitude", prec = "float")
    defs[[paste0("lon_", p)]] <- ncdf4::ncvar_def(paste0("lon_", p), "degrees_east", list(dx, dy), 1e20, "palaeo-longitude", prec = "float")
    defs[[paste0("thk_", p)]] <- ncdf4::ncvar_def(paste0("thk_", p), "m", dims[[p]], 1e20, "ice thickness", prec = "float")
    defs[[paste0("smb_", p)]] <- ncdf4::ncvar_def(paste0("smb_", p), "m ice / yr", dims[[p]], 1e20, "surface mass balance", prec = "float")
    defs[[paste0("bedsink_", p)]] <- ncdf4::ncvar_def(paste0("bedsink_", p), "m", dims[[p]], 1e20, "isostatic sinking of the bed", prec = "float")
  }
  nc <- ncdf4::nc_create(file, defs, force_v4 = TRUE); on.exit(ncdf4::nc_close(nc))
  for (p in names(first$D)) {
    g <- first$D[[p]]$g
    ncdf4::ncvar_put(nc, defs[[paste0("lat_", p)]], t(g$lat)); ncdf4::ncvar_put(nc, defs[[paste0("lon_", p)]], t(g$lon))
  }
  for (k in seq_along(tt)) {
    sn <- ice_snapshot(ice_id, tt[k], cfg)
    for (p in names(sn$D)) {
      for (v in c("H", "smb", "w")) {
        nm <- c(H = "thk_", smb = "smb_", w = "bedsink_")[[v]]
        ncdf4::ncvar_put(nc, defs[[paste0(nm, p)]], t(sn$D[[p]][[v]]), start = c(1, 1, k), count = c(-1, -1, 1))
      }
    }
  }
  for (k in names(st)) if (!is.null(st[[k]]) && length(st[[k]]) == 1 && !is.na(st[[k]])) ncdf4::ncatt_put(nc, 0, paste0("ice_", k), st[[k]])
  ncdf4::ncatt_put(nc, 0, "method", "PDD surface mass balance (Reeh 1991; Calov & Greve 2005) + shallow-ice flow, Lambert azimuthal equal-area polar grids; asynchronous coupling with SPEEDY")
  invisible(file)
}

# ---- drawing ----------------------------------------------------------------------------------
ICE_THICK_MAX <- 4000; ICE_SMB_LIM <- 2
.ice_col <- function(H) {                              # white -> pale blue with thickness
  f <- pmin(pmax(H / ICE_THICK_MAX, 0), 1)
  cbind(0.97 - 0.30 * f, 0.985 - 0.17 * f, 1.0 - 0.05 * f)
}
#' put the ice sheet (layer "ice_thickness") or its surface mass balance ("ice_smb") on a map image
overlay_ice_rgb <- function(rgb, z, lat, ic, layer) {
  if (identical(layer, "ice_smb")) {
    m <- .diverging(as.vector(ic$smb), ICE_SMB_LIM, neg = c("#b2182b", "#ef8a62", "#fddbc7"), pos = c("#d1e5f0", "#67a9cf", "#2166ac"))
    sel <- (z >= 0 | ic$H > 5) & is.finite(as.vector(ic$smb))
    for (c in 1:3) { ch <- rgb[, , c]; ch[sel] <- m[sel, c]; rgb[, , c] <- ch }
    return(rgb)
  }
  ice <- ic$H > 5; if (!any(ice)) return(rgb)
  hs <- hillshade(z - ic$w + ic$H, lat, exag_land = 25)
  m <- .ice_col(as.vector(ic$H))
  for (c in 1:3) { ch <- rgb[, , c]; ch[ice] <- pmin(pmax(m[ice, c] * (0.75 + 0.35 * hs[ice]), 0), 1); rgb[, , c] <- ch }
  rgb
}
#' ice margin and thickness contours (every 1000 m)
draw_ice_lines <- function(ic, lat, lon, cfg, wrap = TRUE) {
  H <- ic$H; if (wrap) { H <- cbind(H, H[, 1]); lon <- c(lon, lon[1] + 360) }
  for (lev in c(10, 1000, 2000, 3000, 4000)) {
    cl <- contourLines(lon, lat, t(H), levels = lev); if (!length(cl)) next
    draw_ll_list(lapply(cl, function(q) data.frame(lon = ((q$x + 180) %% 360) - 180, lat = q$y)), cfg,
                 col = if (lev == 10) "#08306b" else adjustcolor("#2c5f8a", 0.6), lwd = if (lev == 10) 0.9 else 0.4)
  }
}
colour_bar_ice <- function(side = 1, layer) {
  if (identical(layer, "ice_smb")) {
    lim <- ICE_SMB_LIM; v <- seq(-lim, lim, length.out = 201)
    m <- .diverging(v, lim, neg = c("#b2182b", "#ef8a62", "#fddbc7"), pos = c("#d1e5f0", "#67a9cf", "#2166ac"))
    lab <- "Surface mass balance (m ice per year; blue = build-up, red = melt)"; at <- pretty(c(-lim, lim))
  } else {
    v <- seq(0, ICE_THICK_MAX, length.out = 201); m <- .ice_col(v); lab <- "Ice thickness (m)"; at <- seq(0, ICE_THICK_MAX, 1000)
  }
  cc <- rgb(m[, 1], m[, 2], m[, 3]); n <- length(v)
  if (side %in% c(1, 3)) {
    plot.window(xlim = range(v), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(v[-n], 0, v[-1], 1, col = cc[-1], border = NA); box(lwd = 0.5, col = "#999999")
    axis(side, at = at, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.2, 0))
    mtext(lab, side = side, line = 1.2, cex = 0.62, col = "#52514e")
  } else {
    plot.window(xlim = c(0, 1), ylim = range(v), xaxs = "i", yaxs = "i")
    rect(0, v[-n], 1, v[-1], col = cc[-1], border = NA); box(lwd = 0.5, col = "#999999")
    axis(side, at = at, las = 1, cex.axis = 0.6, lwd = 0, lwd.ticks = 0.6, tcl = -0.3, mgp = c(1, 0.3, 0))
    mtext(if (identical(layer, "ice_smb")) "Mass balance (m/yr)" else "Ice thickness (m)", side = side, line = 2.3, cex = 0.62, col = "#52514e")
  }
}
ice_caption <- function(ice) {
  tot <- ice_totals(ice$state)
  sprintf("ice sheets: %s, %.1f kyr: %.3g km3 = %.3g Gt of ice = %.1f m of sea level%s (PDD mass balance + shallow-ice flow with basal sliding%s)",
          ice$label %||% "", ice$state$t / 1000, sum(tot$volume_km3), sum(tot$mass_Gt), sum(tot$sle_m),
          if (isTRUE(ice$lower_sea)) sprintf("; map drawn with the sea %.0f m lower", -(ice$state$sl %||% 0)) else "",
          if (isTRUE(ice$coupled)) ", climate updated with the ice by SPEEDY" else "")
}

# ---- sea ice and land ice ------------------------------------------------------------------
#' Sea ice forms where the sea surface reaches the freezing point of sea water (-1.8 C, not 0 C;
#' SPEEDY: sea ice where the prescribed or slab sea-surface temperature is below -1.8 C).
#' Link between sea ice and land ice, for a hemisphere with land around the pole (>= 15 % land
#' poleward of 60 degrees): permanent (perennial) sea ice is only kept if that hemisphere has
#' land ice - an ice sheet of an ice run, or land where snow survives the summer
#' (positive mass balance by the degree-day model) over at least `min_area_km2`. Otherwise the
#' sea ice of each cell is reduced to at most `cap` in its month with the least ice, so it
#' becomes seasonal. Applied to the experiment runs of SPEEDY (prescribed sea ice and, with the
#' slab ocean, no new sea ice), the ice-sheet coupling runs, and the sea-ice maps / export.
SEA_ICE_LAND_LINK <- list(on = TRUE, min_area_km2 = 5e4, cap = 0.1, polar_land = 0.15)

#' area (km2) of land with a positive surface mass balance, per hemisphere; Tm, Pm: lists of
#' 12 monthly fields list(z, lat, lon) (C, mm/day) that refer to the height field href (m)
land_ice_area <- function(map_file, Tm, Pm, href, st = ice_settings()) {
  g <- geo_1deg_cached(map_file); lat <- GRID1$lat
  G <- expand.grid(lon = GRID1$lon, lat = lat)
  on1 <- function(F) matrix(clim_interp(F, G$lat, G$lon), 180, 360, byrow = TRUE)
  hl <- ifelse(g$lf > 0, g$h / pmax(g$lf, 1e-6), 0)
  dz <- hl - on1(href)
  T <- lapply(1:12, function(m) on1(Tm[[m]]) - st$lapse / 1000 * dz)
  P <- lapply(1:12, function(m) pmax(on1(Pm[[m]]), 0) * pmin(pmax(exp(-st$precip_k * st$lapse / 1000 * dz), 0.1), 2))
  smb <- ice_smb(T, P, st)$smb
  area <- 111.2^2 * outer(cos(lat * pi / 180), rep(1, 360)) * g$lf
  pos <- is.finite(smb) & smb > 0 & g$lf > 0.3
  c(south = sum(area[pos & lat < 0]), north = sum(area[pos & lat > 0]))
}
land_ice_flags <- function(areas, map_file = NULL) {
  if (!isTRUE(SEA_ICE_LAND_LINK$on)) return(c(south = TRUE, north = TRUE))
  f <- areas >= SEA_ICE_LAND_LINK$min_area_km2
  # a pole without land (open polar ocean): sea ice follows the temperature alone
  if (!is.null(map_file)) f <- f | .polar_land(map_file) < SEA_ICE_LAND_LINK$polar_land
  f
}
#' fraction of land poleward of 60 degrees, per hemisphere
.polar_land <- function(map_file) {
  g <- geo_1deg_cached(map_file); w <- cos(GRID1$lat * pi / 180)
  frac <- function(sel) sum(w[sel] * rowSums(g$lf[sel, , drop = FALSE])) / (sum(w[sel]) * 360)
  c(south = frac(GRID1$lat <= -60), north = frac(GRID1$lat >= 60))
}
#' monthly sea ice (array lon x lat x 12 or list of 12 lat x lon matrices) with permanent ice
#' removed in hemispheres without land ice; lat: latitudes along the lat dimension
cap_perennial <- function(Z, lat, flags, lat_dim = 2) {
  if (all(flags)) return(Z)
  is_list <- is.list(Z)
  A <- if (is_list) simplify2array(lapply(Z, function(F) F$z)) else Z      # list: lat x lon x 12
  ld <- if (is_list) 1 else lat_dim
  dims <- dim(A); cells <- matrix(A, ncol = 12)
  latc <- if (ld == 1) rep(lat, dims[2]) else rep(lat, each = dims[1])
  hem <- (latc < 0 & !flags[["south"]]) | (latc > 0 & !flags[["north"]])
  if (any(hem)) {
    mmin <- max.col(-cells[hem, , drop = FALSE], ties.method = "first")
    idx <- cbind(which(hem), mmin)
    cells[idx] <- pmin(cells[idx], SEA_ICE_LAND_LINK$cap)
  }
  A <- array(cells, dims)
  if (is_list) { for (m in 1:12) Z[[m]]$z <- A[, , m]; Z } else A
}
#' the land-ice flags for an experiment's SPEEDY forcing: HadCM3L monthly climate of the age,
#' with the sea-surface shift dT (land 1.3 x), on the map; cached
.LIF <- new.env()
land_ice_flags_forcing <- function(age, shift, map_file, cfg = SILUCARB) {
  dT <- .shift_dT(shift); g <- .shift_grad(shift)
  key <- paste(age, round(dT, 3), round(g, 3), .map_tag(map_file))
  if (!is.null(.LIF[[key]])) return(.LIF[[key]])
  mf <- climate_month_files(age, cfg)
  pat <- function(lat) 0.75 + 0.75 * sin(lat * pi / 180)^2
  Tm <- lapply(seq_along(mf), function(m) { t <- .read_clim_var(mf[m], "temp_mm_1_5m")
    t$z <- t$z - 273.15 + 1.3 * dT * pat(t$lat) + gradient_anomaly(age, g, t$lat, m, cfg); t })
  Pm <- lapply(mf, function(f) { p <- .read_clim_var(f, "precip_mm_srf"); p$z <- p$z * 86400; p })
  href <- .speedy_in_field(speedy_geography(find_dem(cfg$dir_dem, age))$orog)
  a <- land_ice_area(map_file, Tm, Pm, href)
  .LIF[[key]] <- list(flags = land_ice_flags(a, map_file), area = a)
  .LIF[[key]]
}


#' the largest equator-to-pole gradient factor (closest to HadCM3L) at which the land where snow
#' survives the summer stays at or below max_area_km2 (0 = ice-free start), for a map and an
#' experiment (its CO2 / Sun / sensitivity). Uses the HadCM3L climate + the CO2 / Sun shift +
#' the gradient change and the degree-day model, without SPEEDY runs (a quick first estimate).
#' Returns the factor, the area at it (with the margin) and a note.
#' margin_K: the search is made this much colder (K) as a safety margin: the SPEEDY response to the
#' map, the ice-elevation feedback and the coarse 1-degree heights can make it colder than estimated.
gradient_for_ice <- function(map, x = NULL, max_area_km2 = 0, lower = 0.3, margin_K = 2, cfg = SILUCARB, tol = 0.01) {
  m <- climate_map(map, cfg = cfg); xr <- .resolve(x %||% climate_experiment(), m$age)
  area <- function(g) sum(land_ice_flags_forcing(m$age, c(dT = xr$dT - margin_K / 1.3, grad = g), m$file, cfg)$area)
  a1 <- area(1)
  if (a1 <= max_area_km2) return(list(gradient = 1, area_km2 = a1, note = "no change needed"))
  if (area(lower) > max_area_km2) return(list(gradient = lower, area_km2 = area(lower), note = sprintf("even %.2f is not enough", lower)))
  lo <- lower; hi <- 1                                   # area(lo) <= target < area(hi)
  while (hi - lo > tol) { mid <- (lo + hi) / 2; if (area(mid) <= max_area_km2) lo <- mid else hi <- mid }
  g <- floor(lo * 100) / 100
  list(gradient = g, area_km2 = area(g), note = "")
}
