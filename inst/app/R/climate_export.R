# =============================================================================
#  Export of a climate run (Climate model tab) as NetCDF: every field SPEEDY writes, for all
#  12 months, plus the HadCM3L-based fields that the maps show.
# =============================================================================

#' monthly fields of a finished climate run (a row of climate_runs() or its id):
#' list(run, reference, change, map_effect) of arrays (96 lon, 48 lat S->N, 12, variables)
climate_run_fields <- function(run, cfg = SILUCARB) {
  if (is.character(run)) { r <- climate_runs(cfg, status = FALSE); run <- r[r$id == run, , drop = FALSE] }
  if (!nrow(run)) stop("unknown climate run", call. = FALSE)
  s <- run_setup(run, cfg); jobs <- speedy_jobs(s$map, s$x, cfg)
  if (any(vapply(jobs, .job_status, "") != "done")) stop("this climate run is not finished yet", call. = FALSE)
  rd <- lapply(jobs, function(j) .read_means(j$dir))
  vars <- Reduce(intersect, lapply(rd, function(a) dimnames(a)[[4]]))
  rd <- lapply(rd, function(a) a[, , , vars, drop = FALSE])
  change <- 0 * rd$corr
  if (!is.null(rd$forcing)) change <- change + rd$forcing - rd$corr
  if (!is.null(rd$exp)) change <- change + rd$exp - rd$slab
  list(setup = s, jobs = jobs, vars = vars, reference = rd$corr, change = change, run = rd$corr + change,
       map_effect = rd$corr - rd$orig, experiment = length(jobs) >= 3)
}

#' write a climate run to a NetCDF file (one file, all fields, 12 months)
#'   <var>             SPEEDY climate of the run (reference + change), 96 x 48 Gaussian grid
#'   <var>_change      change by the orbit / CO2 / Sun experiment (run - reference)  [experiments]
#'   <var>_map_effect  your map - the original PALEOMAP map (SPEEDY)
#'   precip_map, tas_map, u10_map, v10_map: HadCM3L + SPEEDY response as drawn on the maps
#'                     (HadCM3L grid, 3.75 x 2.5 deg)
export_climate_run <- function(run, file, cfg = SILUCARB, o2_pal = 1, lapse_rate = 6.5) {
  F <- climate_run_fields(run, cfg); s <- F$setup; age <- s$map$age; xr <- .resolve(s$x, age)
  lat <- rev(speedy_lats); lon <- SPEEDY_LON
  d_lon <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  d_lat <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  d_mon <- ncdf4::ncdim_def("month", "month of the year", 1:12)
  meta <- SPEEDY_VARS[match(F$vars, SPEEDY_VARS$name), ]
  vdefs <- list()
  add <- function(name, units, long, dims = list(d_lon, d_lat, d_mon))
    vdefs[[name]] <<- ncdf4::ncvar_def(name, units, dims, missval = 1e20, longname = long, prec = "float")
  for (i in seq_len(nrow(meta))) {
    add(meta$name[i], meta$units[i], paste0(meta$long_name[i], " (SPEEDY, this run)"))
    if (F$experiment) add(paste0(meta$name[i], "_change"), meta$units[i], paste0(meta$long_name[i], ": change by the experiment (run - reference)"))
    add(paste0(meta$name[i], "_map_effect"), meta$units[i], paste0(meta$long_name[i], ": your map - original PALEOMAP map"))
  }
  # HadCM3L + SPEEDY response, as on the maps
  mf <- climate_month_files(age, cfg)
  g0 <- .read_clim_var(mf[1], "temp_mm_1_5m")
  dh_lon <- ncdf4::ncdim_def("lon_map", "degrees_east", g0$lon); dh_lat <- ncdf4::ncdim_def("lat_map", "degrees_north", g0$lat)
  hd <- list(dh_lon, dh_lat, d_mon)
  add("precip_map", "mm/day", "precipitation: HadCM3L + SPEEDY response (as on the maps)", hd)
  add("tas_map", "degC", "near-surface (1.5 m) air temperature: HadCM3L + SPEEDY response", hd)
  add("u10_map", "m/s", "eastward 10 m wind: HadCM3L + SPEEDY response", hd)
  add("v10_map", "m/s", "northward 10 m wind: HadCM3L + SPEEDY response", hd)
  add("pme_map", "mm/day", "precipitation - evaporation: HadCM3L + SPEEDY response", hd)
  add("runoff_map", "mm/day", "total runoff (surface + subsurface; land): HadCM3L + SPEEDY change of P - E", hd)
  hd1 <- list(dh_lon, dh_lat)
  add("pme_map_annual", "mm/yr", "precipitation - evaporation, annual total: HadCM3L + SPEEDY response", hd1)
  add("runoff_map_annual", "mm/yr", "total runoff, annual total: HadCM3L + SPEEDY change of P - E", hd1)
  # ocean of the HadCM3L run (not changed by SPEEDY; currents on the velocity grid)
  mo <- climate_month_files_ocean(age, cfg)
  u0 <- .read_clim_var(mo[1], "ucurrTot_mm_dpth")
  du_lon <- ncdf4::ncdim_def("lon_uv", "degrees_east", u0$lon); du_lat <- ncdf4::ncdim_def("lat_uv", "degrees_north", u0$lat)
  d_seas <- ncdf4::ncdim_def("season", "1 annual, 2 DJF, 3 MAM, 4 JJA, 5 SON", 1:5)
  add("salinity_5m", "psu", "sea-surface salinity at 5 m (HadCM3L run)", hd)
  add("salinity_5m_seasonal", "psu", "sea-surface salinity at 5 m: annual, DJF, MAM, JJA, SON means (HadCM3L)", list(dh_lon, dh_lat, d_seas))
  add("w_10m", "m/day", "ocean vertical velocity at 10 m, + = upward (HadCM3L)", hd)
  add("mixed_layer_depth", "m", "ocean mixed-layer depth (HadCM3L)", hd)
  add("o2_sat_surface", "umol/kg", sprintf("surface O2 at saturation (Garcia & Gordon 1992), SST incl. SPEEDY change, atmospheric O2 %g PAL", o2_pal), hd)
  add("w_ekman_change", "m/day", "change of Ekman pumping by the SPEEDY wind-stress change (+ = more upwelling)", list(d_lon, d_lat, d_mon))
  add("sea_ice_map", "fraction", "sea-ice fraction: HadCM3L + SPEEDY change, permanent ice only with land ice (land-ice link)", hd)
  # near-surface temperature brought to the map heights (0.5 deg) with the lapse rate
  dm <- read_dem_nc(s$map$file); blk <- .dem_block(dm, 0.5)
  d5_lon <- ncdf4::ncdim_def("lon_05", "degrees_east", blk$lon); d5_lat <- ncdf4::ncdim_def("lat_05", "degrees_north", blk$lat)
  add("surface_height_05", "m", "mean land surface height of the map in 0.5 deg cells (0 at sea)", list(d5_lon, d5_lat))
  add("tas_downscaled", "degC", sprintf("near-surface air temperature at the map heights (0.5 deg): tas_map - %g K/km x (height - SPEEDY T30 terrain of the map)", lapse_rate),
      list(d5_lon, d5_lat, d_mon))
  add("u_ocean_5m", "cm/s", "eastward ocean current at 5 m (HadCM3L)", list(du_lon, du_lat, d_mon))
  add("v_ocean_5m", "cm/s", "northward ocean current at 5 m (HadCM3L)", list(du_lon, du_lat, d_mon))
  # rivers and freshwater into the sea (1 deg, annual)
  d1_lon <- ncdf4::ncdim_def("lon_1deg", "degrees_east", GRID1$lon); d1_lat <- ncdf4::ncdim_def("lat_1deg", "degrees_north", GRID1$lat)
  h1 <- list(d1_lon, d1_lat)
  add("river_flow", "km3/yr", "river flow: annual runoff routed downhill over the map to the coast (land, 1 deg)", h1)
  add("river_discharge_mouth", "km3/yr", "river water reaching the sea, per coastal sea cell (river mouths)", h1)
  add("river_input_annual", "mm/yr", "river water into the sea, spread ~300 km (Gaussian sigma 150 km) around the mouths", h1)
  add("freshwater_into_sea_annual", "mm/yr", "freshwater into the sea: P - E over the sea + river water (+ = freshening)", h1)
  add("salinity_source", "class", paste("why the annual surface water is fresher / saltier than the ocean mean: 1 river water, 2 rain (P > E),",
                                        "3 fresher by currents/mixing, 4 near the mean (+-0.3 psu), 5 saltier by currents/mixing, 6 evaporation (E > P)"), h1)
  nc <- ncdf4::nc_create(file, vdefs, force_v4 = TRUE); on.exit(ncdf4::nc_close(nc))
  put <- function(name, a) ncdf4::ncvar_put(nc, vdefs[[name]], ifelse(is.finite(a), a, 1e20))
  for (v in F$vars) {
    put(v, F$run[, , , v])
    if (F$experiment) put(paste0(v, "_change"), F$change[, , , v])
    put(paste0(v, "_map_effect"), F$map_effect[, , , v])
  }
  # SPEEDY response on the HadCM3L grid (wraps in longitude)
  sp_lon <- ((lon + 180) %% 360) - 180; o <- order(sp_lon)
  resp <- function(k, m) list(z = t(F$map_effect[, , m, k] + F$change[, , m, k])[, o], lat = lat, lon = sp_lon[o])
  g <- expand.grid(lon = g0$lon, lat = g0$lat)
  on_h <- function(Fl) matrix(clim_interp(Fl, g$lat, g$lon), length(g0$lon), length(g0$lat))
  out <- list(precip_map = array(NA, c(length(g0$lon), length(g0$lat), 12)))
  out$tas_map <- out$u10_map <- out$v10_map <- out$pme_map <- out$runoff_map <- out$precip_map
  out$salinity_5m <- out$w_10m <- out$mixed_layer_depth <- out$o2_sat_surface <- out$precip_map
  out$u_ocean_5m <- out$v_ocean_5m <- array(NA, c(length(u0$lon), length(u0$lat), 12))
  g_uv <- NULL
  for (m in 1:12) {
    pr <- .read_clim_var(mf[m], "precip_mm_srf"); pr$z <- pr$z * 86400
    tt <- .read_clim_var(mf[m], "temp_mm_1_5m"); tt$z <- tt$z - 273.15
    uu <- .read_clim_var(mf[m], "u_mm_10m"); vv <- .read_clim_var(mf[m], "v_mm_10m")
    out$precip_map[, , m] <- pmax(on_h(pr) + on_h(resp("precip", m)), 0)
    out$tas_map[, , m] <- on_h(tt) + on_h(resp("t0", m))
    out$u10_map[, , m] <- on_h(uu) + on_h(resp("u0", m))
    out$v10_map[, , m] <- on_h(vv) + on_h(resp("v0", m))
    ev <- .read_clim_var(mf[m], "lh_mm_srf"); ev$z <- ev$z / 2.501e6 * 86400
    dpe <- on_h(resp("precip", m)) - (if ("evap" %in% F$vars) on_h(resp("evap", m)) else 0)
    out$pme_map[, , m] <- on_h(pr) - on_h(ev) + dpe
    ro <- .read_clim_var(mf[m], "srfRunoff_mm_srf"); rs <- .read_clim_var(mf[m], "subsrfRunoff_mm_srf")
    ro$z <- (ifelse(is.finite(ro$z), ro$z, 0) + ifelse(is.finite(rs$z), rs$z, 0)) * 86400; ro$z[!is.finite(rs$z) & !is.finite(.read_clim_var(mf[m], "srfRunoff_mm_srf")$z)] <- NA
    out$runoff_map[, , m] <- ifelse(is.finite(t(ro$z)), pmax(t(ro$z) + dpe, 0), NA)
    sal <- .read_clim_var(mo[m], "salinity_mm_dpth"); out$salinity_5m[, , m] <- t(sal$z * 1000 + 35)
    w <- .read_clim_var(mo[m], "W_mm_dpth"); out$w_10m[, , m] <- t(w$z * 864)
    out$mixed_layer_depth[, , m] <- t(.read_clim_var(mo[m], "mixLyrDpth_mm_uo")$z)
    to <- .read_clim_var(mo[m], "temp_mm_dpth")
    dsst <- if ("sst" %in% F$vars) on_h(resp("sst", m)) else on_h(resp("t0", m))
    out$o2_sat_surface[, , m] <- o2_saturation(t(to$z) + dsst, t(sal$z * 1000 + 35), o2_pal)
    out$u_ocean_5m[, , m] <- t(.read_clim_var(mo[m], "ucurrTot_mm_dpth")$z); out$v_ocean_5m[, , m] <- t(.read_clim_var(mo[m], "vcurrTot_mm_dpth")$z)
  }
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  out$pme_map_annual <- apply(out$pme_map, 1:2, function(v) sum(v * days))
  out$runoff_map_annual <- apply(out$runoff_map, 1:2, function(v) sum(v * days))
  ss <- list(1:12, c(12, 1, 2), 3:5, 6:8, 9:11)
  out$salinity_5m_seasonal <- array(NA, c(length(g0$lon), length(g0$lat), 5))
  for (k in 1:5) out$salinity_5m_seasonal[, , k] <- apply(out$salinity_5m[, , ss[[k]], drop = FALSE], 1:2, mean)
  if (all(c("taux", "tauy") %in% F$vars)) {
    dtx <- F$map_effect[, , , "taux"] + F$change[, , , "taux"]; dty <- F$map_effect[, , , "tauy"] + F$change[, , , "tauy"]
    out$w_ekman_change <- array(NA, c(96, 48, 12)); for (m in 1:12) out$w_ekman_change[, , m] <- ekman_w(dtx[, , m], dty[, , m])
  }
  si <- tryCatch(sea_ice_for_map(s$map, "rerun", s$x, cfg, monthly = TRUE), error = function(e) NULL)
  if (!is.null(si$monthly)) { out$sea_ice_map <- array(NA, c(length(g0$lon), length(g0$lat), 12)); for (m in 1:12) out$sea_ice_map[, , m] <- t(si$monthly[[m]]$z) }
  href <- .speedy_in_field(speedy_geography(s$map$file)$orog)
  G5 <- expand.grid(lon = blk$lon, lat = blk$lat)
  h5 <- matrix(clim_interp(href, G5$lat, G5$lon), length(blk$lon), length(blk$lat))
  out$surface_height_05 <- t(blk$h)
  out$tas_downscaled <- array(NA, c(length(blk$lon), length(blk$lat), 12))
  for (m in 1:12) {
    Fm <- list(z = t(out$tas_map[, , m]), lat = g0$lat, lon = g0$lon)
    t5 <- matrix(clim_interp(Fm, G5$lat, G5$lon), length(blk$lon), length(blk$lat))
    out$tas_downscaled[, , m] <- temp_downscale(t5, h5, t(blk$h), lapse_rate)
  }
  fw <- tryCatch({
    hl <- function(a) list(z = t(a), lat = g0$lat, lon = g0$lon)
    freshwater_fields(list(runoff = hl(out$runoff_map_annual / 365.25), pme = hl(out$pme_map_annual / 365.25),
                           salinity = hl(out$salinity_5m_seasonal[, , 1])), s$map$file)
  }, error = function(e) { message("rivers: ", conditionMessage(e)); NULL })
  if (!is.null(fw)) {
    out$river_flow <- t(ifelse(fw$rivers$route$land, fw$rivers$Q$acc, NA)); out$river_discharge_mouth <- t(fw$rivers$Q$mouth)
    out$river_input_annual <- t(fw$river_input$z * 365.25); out$freshwater_into_sea_annual <- t(fw$freshwater$z * 365.25)
    out$salinity_source <- t(fw$salinity_source$z)
  }
  for (k in names(out)) put(k, out[[k]])
  e <- climate_exp(age)
  att <- list(title = paste("Climate run:", s$map$name, "-", xr$label),
              map = s$map$name, map_file = s$map$file, map_age_Ma = age,
              experiment = if (.is_reference(s$x, age)) "reference (HadCM3L orbit, CO2 and Sun)" else xr$label,
              eccentricity = xr$ecc, obliquity_deg = xr$obl, longitude_of_perihelion_deg = xr$peri,
              co2_ppm = xr$co2_ppm, co2_PAL = xr$co2_pal_eff, sun_age_Ma = xr$sun, solar_factor = solar_factor_at(xr$sun),
              climate_sensitivity_K = xr$sensitivity, run_length = s$x$length,
              sst_shift_K = xr$dT, equator_to_pole_gradient_factor = xr$gradient %||% 1, hadcm3l_run = e$id, hadcm3l_age_Ma = e$age, hadcm3l_co2_ppm = e$co2,
              spinup_and_averaged_years = paste(vapply(F$jobs, function(j) sprintf("%s %d+%d", basename(j$dir), j$spinup, j$years), ""), collapse = "; "),
              method = paste("SPEEDY (Molteni 2003; Kucharski et al. 2013) runs on the map; run = reference + change;",
                             "change = (CO2/Sun run - reference) + (orbit slab run - slab run with today's orbit);",
                             "*_map fields = HadCM3L (Valdes et al. 2021) + SPEEDY response (map effect + change)."),
              grid = "SPEEDY T30 Gaussian grid, latitudes south to north; *_map fields on the HadCM3L grid",
              created = format(Sys.time(), "%Y-%m-%d %H:%M"), source = "paleogeographic map editor")
  for (k in names(att)) ncdf4::ncatt_put(nc, 0, k, att[[k]])
  invisible(file)
}

# ---- sea ice ----------------------------------------------------------------------------
#' monthly SPEEDY response of one variable for a map and experiment: (your map - original
#' map) + (experiment - reference); array 96 x 48 (south->north) x 12, or NULL
.speedy_monthly_response <- function(map, x, var, cfg = SILUCARB) {
  jobs <- tryCatch(speedy_jobs(map, x, cfg), error = function(e) NULL)
  if (is.null(jobs) || any(vapply(jobs, .job_status, "") != "done")) return(NULL)
  rd <- lapply(jobs, function(j) { a <- .read_means(j$dir); if (var %in% dimnames(a)[[4]]) a[, , , var] else NULL })
  if (any(vapply(rd, is.null, TRUE))) return(NULL)
  d <- rd$corr - rd$orig
  if (!is.null(rd$forcing)) d <- d + rd$forcing - rd$corr
  if (!is.null(rd$exp)) d <- d + rd$exp - rd$slab
  d
}
#' sea ice of a map: annual minimum (permanent ice) and maximum of the monthly sea-ice
#' fraction - HadCM3L, plus the SPEEDY change for a climate run (method "rerun")
sea_ice_for_map <- function(map, method = "local", experiment = NULL, cfg = SILUCARB, ice_state = NULL, monthly = FALSE) {
  age <- if (is.list(map)) map$age else map
  mf <- climate_month_files(age, cfg)
  m <- tryCatch(climate_map(map, cfg = cfg), error = function(e) NULL)
  resp <- if (method == "rerun" && !is.null(m)) .speedy_monthly_response(m, experiment, "sea_ice", cfg)
  mon <- lapply(1:12, function(k) {
    x <- .read_clim_var(mf[k], "iceconc_mm_srf"); x$z[!is.finite(x$z)] <- 0
    if (!is.null(resp)) {
      g <- expand.grid(lon = x$lon, lat = x$lat)
      x$z <- x$z + matrix(clim_interp(.speedy_out_field(resp[, , k]), g$lat, g$lon), length(x$lat), length(x$lon), byrow = TRUE)
    }
    x$z <- pmin(pmax(x$z, 0), 1); x
  })
  # link to land ice: permanent sea ice only in a hemisphere with land ice (see SEA_ICE_LAND_LINK)
  link <- NULL
  if (isTRUE(SEA_ICE_LAND_LINK$on) && !is.null(m)) {
    rt <- if (method == "rerun") .speedy_monthly_response(m, experiment, "t0", cfg)
    rp <- if (method == "rerun") .speedy_monthly_response(m, experiment, "precip", cfg)
    if (!is.null(rt) && !is.null(rp)) {
      add <- function(F, R) { g <- expand.grid(lon = F$lon, lat = F$lat)
        F$z <- F$z + matrix(clim_interp(.speedy_out_field(R), g$lat, g$lon), length(F$lat), length(F$lon), byrow = TRUE); F }
      Tm <- lapply(1:12, function(k) { t <- .read_clim_var(mf[k], "temp_mm_1_5m"); t$z <- t$z - 273.15; add(t, rt[, , k]) })
      Pm <- lapply(1:12, function(k) { p <- .read_clim_var(mf[k], "precip_mm_srf"); p$z <- p$z * 86400; add(p, rp[, , k]) })
      area <- land_ice_area(m$file, Tm, Pm, .speedy_in_field(speedy_geography(m$file)$orog))
    } else area <- land_ice_flags_forcing(age, 0, m$file, cfg)$area
    if (!is.null(ice_state)) {                         # ice sheets of an ice run count as land ice
      ia <- vapply(c("south", "north"), function(p) { D <- ice_state$D[[p]]; if (is.null(D)) 0 else sum(D$H > 10) * D$g$res^2 / 1e6 }, 0)
      area <- pmax(area, ia)
    }
    link <- list(area = area, flags = land_ice_flags(area, m$file))
    mon <- cap_perennial(mon, mon[[1]]$lat, link$flags)
  }
  mn <- mon[[1]]; mx <- mon[[1]]
  mn$z <- Reduce(pmin, lapply(mon, `[[`, "z")); mx$z <- Reduce(pmax, lapply(mon, `[[`, "z"))
  # smooth the coarse model cells (~300 km) so the ice edges are not blocky
  mn$z <- pmin(pmax(smooth_km(mn$z, mn$lat, 300, nbands = 24), 0), 1); mx$z <- pmin(pmax(smooth_km(mx$z, mx$lat, 300, nbands = 24), 0), 1)
  list(min = mn, max = mx, changed = !is.null(resp), link = link, monthly = if (monthly) mon)
}
SEA_ICE_EDGE <- 0.15
#' sea ice on a map image: permanent (all months >= 15 %) white, seasonal paler
overlay_sea_ice_rgb <- function(rgb, z, si_min, si_max) {
  sea <- z < 0; perm <- sea & si_min >= SEA_ICE_EDGE; seas <- sea & si_max >= SEA_ICE_EDGE & !perm
  white <- c(0.94, 0.96, 0.99); pale <- c(0.80, 0.88, 0.95)
  for (c in 1:3) { ch <- rgb[, , c]
    ch[perm] <- white[c]; ch[seas] <- 0.2 * ch[seas] + 0.8 * pale[c]; rgb[, , c] <- ch }
  rgb
}
draw_sea_ice_lines <- function(si_min, si_max, z, lat, lon, cfg, wrap = TRUE) {
  for (k in 1:2) {
    M <- if (k == 1) si_min else si_max; M[z >= 0] <- NA      # edges at sea only
    if (wrap) { M <- cbind(M, M[, 1]); lo <- c(lon, lon[1] + 360) } else lo <- lon
    cl <- contourLines(lo, lat, t(M), levels = SEA_ICE_EDGE); if (!length(cl)) next
    draw_ll_list(lapply(cl, function(q) data.frame(lon = ((q$x + 180) %% 360) - 180, lat = q$y)), cfg,
                 col = "#2f5577", lwd = if (k == 1) 1.1 else 0.8, lty = if (k == 1) 1 else 2)
  }
}
sea_ice_caption <- function(si) {
  lk <- si$link
  note <- if (!is.null(lk) && !all(lk$flags)) paste0("; no land ice in the ", paste(c(south = "southern", north = "northern")[names(lk$flags)[!lk$flags]], collapse = " and "),
                                                    " hemisphere, so no permanent sea ice there (land-ice link)") else ""
  paste0("sea ice (HadCM3L", if (isTRUE(si$changed)) " + SPEEDY change" else "",
         "): white / solid line = permanent (all months \u2265 15 %), pale / dashed = seasonal", note)
}

#' the 12 monthly HadCM3L ocean files of the run for an age (downloaded once, ~6 MB)
climate_month_files_ocean <- function(age, cfg = SILUCARB) {
  e <- climate_exp(age); d <- file.path(cfg$dir_data, "climate", e$id)
  f <- file.path(d, sprintf("%so.pfcl%s.nc", e$id, .MONTHS))
  for (k in seq_along(f)) .fetch_retry(sprintf(CLIMATE_URL, e$id, basename(f[k])), f[k])
  if (!all(file.exists(f))) stop("monthly ocean files for ", age, " Ma could not be downloaded", call. = FALSE)
  f
}

# ---- Ekman pumping and surface oxygen --------------------------------------------------------
#' Ekman pumping (m/day, + = upwelling) from the surface stress on the ATMOSPHERE as SPEEDY writes
#' it (taux, tauy: 96 x 48, south->north); the ocean feels the opposite stress. Not computed within
#' 5 degrees of the equator (f -> 0).
ekman_w <- function(taux, tauy) {
  lat <- rev(speedy_lats); phi <- lat * pi / 180; f <- 2 * 7.292e-5 * sin(phi)
  R <- 6.371e6; dlon <- 3.75 * pi / 180
  tx <- -taux; ty <- -tauy                                   # stress on the ocean
  F <- matrix(f, 96, 48, byrow = TRUE); F[abs(F) < 2 * 7.292e-5 * sin(5 * pi / 180)] <- NA
  ax <- ty / F; ay <- tx / F
  dxdl <- (ax[c(2:96, 1), ] - ax[c(96, 1:95), ]) / (2 * dlon)                       # d/dlambda
  cphi <- matrix(cos(phi), 96, 48, byrow = TRUE)
  dyd <- matrix(NA, 96, 48)
  for (j in 2:47) dyd[, j] <- (ay[, j + 1] - ay[, j - 1]) / (phi[j + 1] - phi[j - 1])
  w <- (dxdl / (R * cphi) - dyd / R) / 1025                  # m/s
  w * 86400
}
#' change of the Ekman pumping for a climate run (map + experiment), months `mo`, on SPEEDY's grid
ekman_change <- function(map, x, mo, cfg = SILUCARB) {
  tx <- .speedy_monthly_response(map, x, "taux", cfg); ty <- .speedy_monthly_response(map, x, "tauy", cfg)
  if (is.null(tx) || is.null(ty)) return(NULL)
  w <- ekman_w(apply(tx[, , mo, drop = FALSE], 1:2, mean), apply(ty[, , mo, drop = FALSE], 1:2, mean))
  w[!is.finite(w)] <- 0
  .speedy_out_field(w)
}
#' dissolved O2 at saturation with the atmosphere (umol/kg; Garcia & Gordon 1992, Benson & Krause
#' data), for temperature t (C), salinity S (psu) and atmospheric O2 in PAL (1 = 20.95 %)
o2_saturation <- function(t, S, o2_pal = 1) {
  Ts <- log((298.15 - t) / (273.15 + pmax(t, -2)))
  lnC <- 5.80871 + 3.20291 * Ts + 4.17887 * Ts^2 + 5.10006 * Ts^3 - 9.86643e-2 * Ts^4 + 3.80369 * Ts^5 +
    S * (-7.01577e-3 - 7.70028e-3 * Ts - 1.13864e-2 * Ts^2 - 9.51519e-3 * Ts^3) - 2.75915e-7 * S^2
  exp(lnC) * o2_pal
}

#' mean land surface height (m, sea = 0) of a map in cells of `res` degrees
.dem_block <- function(dem, res = 0.5) {
  nlat <- round(180 / res); nlon <- round(360 / res)
  ii <- pmin(pmax(floor((dem$lat + 90) / res) + 1, 1), nlat)
  jj <- pmin(pmax(floor(((dem$lon + 180) %% 360) / res) + 1, 1), nlon)
  agg <- function(M) { a <- rowsum(M, ii, reorder = TRUE); t(rowsum(t(a), jj, reorder = TRUE)) }
  h <- agg(pmax(dem$z, 0)) / agg(matrix(1, nrow(dem$z), ncol(dem$z)))
  list(h = h, lat = -90 + (seq_len(nlat) - 0.5) * res, lon = -180 + (seq_len(nlon) - 0.5) * res)
}
