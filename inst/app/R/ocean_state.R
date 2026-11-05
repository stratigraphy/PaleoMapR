# =============================================================================
#  Ocean states of the Silurian: Jeppsson's P and S episodes
#
#  Jeppsson (1990, J. Geol. Soc. 147; 1998, N.Y. State Mus. Bull. 491): the Silurian ocean
#  switched between
#    P (primo) episodes: deep water formed at high latitudes by cooling (thermal mode); an
#       oxygenated deep sea, upwelling of nutrients, high productivity, marls and shales;
#    S (secundo) episodes: deep water formed in warm, evaporative low-latitude (shelf) seas
#       (haline mode); poorly ventilated deep sea, nutrient-poor surface, carbonates and reefs.
#  Bickert et al. (1997) linked humid / arid periods to estuarine / anti-estuarine circulation
#  of the epicontinental seas and gave the d13C the opposite sign; Cramer & Saltzman (2005,
#  2007), Loydell (2007), Trotter et al. (2016) and Fryda et al. (2021) tie the positive d13C
#  excursions to cooling and glaciation. Pohl et al. (2021, Nat. Geosci., cGENIE) show how
#  cooling (lower CO2) can REDUCE high-latitude deep-water formation: seasonal sea ice and its
#  melt water stratify the polar sea, and the deep sea becomes poorly ventilated while the
#  shelves are better oxygenated. So the same two physical modes are linked to high CO2 by
#  some authors and to low CO2 by others.
#
#  Approach here: the state is not set by a CO2 label but diagnosed from the climate:
#   - which surface water is densest in winter (polar water: thermal; warm saline shelf
#     water: haline), from the sea-surface temperature, salinity and freshwater of the climate
#     runs (HadCM3L + SPEEDY), and
#   - a small ocean box model (Stommel 1961 type with several possible sources: south polar,
#     north polar and low-latitude shelf seas, an open low-latitude surface box and a deep
#     box; temperature restored, salinity driven by the freshwater flux, phosphate and O2) that
#     finds which sources ventilate the deep sea, how strongly, and the deep O2 - from two
#     starts (polar and shelf), so two possible states (bistability, hysteresis) show.
#   P-like = deep water mainly from the poles; S-like = mainly from the shelf seas; stagnant =
#   weak overturning. Sea ice: `seaice_effect` -1 (melt-water lid, Pohl et al. 2021) ... +1
#   (brine rejection makes polar water denser, as around Antarctica today).
#  The CO2 x gradient-factor plane: CO2 and the equator-to-pole gradient factor both change the
#  polar winter temperature, so they trade off. The regime map scans both (sea surface shifted as
#  in the SPEEDY runs; freshwater scaled ~7 %/K, Held & Soden 2006) and draws the constraints
#  that separate them: tropical SST from conodont d18O (target band) and the land-ice line.
# =============================================================================

#' melt-water lid calibrated so that the global P/S switch lies at 6-8 PAL (435 Ma, SPEEDY baltic_switch_1
#' climates, gradient 1: switch at ~6.2 PAL with the high-e orbit, ~8.0 PAL with the low-e orbit)
LID_CALIBRATED <- list(seaice_effect = -1, seaice_psu = 0.95)
#' polar fresh-water corrections (Sv) so that, as in HadCM3L at 435 Ma (teXPj) and cGENIE at 440 Ma, the
#' south polar sea forms the deep water and the north polar sea is fresh and does not sink (box polar
#' salinity minus open sea: S -1.8, N -4.0 psu; HadCM3L -2.5 / -4.4). Without them the box model sinks
#' in the north and its south polar sea never sinks.
SALT_CALIBRATED <- list(fwcorr_polarS_Sv = 0.25, fwcorr_polarN_Sv = 3)
#' ice-sheet melt water: the box model leaves the P state at ~0.16 (4 PAL) / 0.18 Sv (8 PAL), cGENIE at
#' 0.18 / 0.22 Sv; melt water is given in cGENIE-equivalent Sv and scaled by this factor
MELT_EFF <- 0.85
OCEAN_ICESHEET <- c(SALT_CALIBRATED, list(melt_eff = MELT_EFF))
#' melt water (cGENIE-equivalent Sv) at a CO2, linear in log2(CO2) between the values at 8 and 4 PAL
melt_at_co2 <- function(co2, at8 = 0.15, at4 = 0.40) pmax(0, at8 + (at4 - at8) * log2(8 / co2))
OCEAN_PARAMS <- list(polar_lat = 55, shelf_depth = 200, h_polar = 250, h_open = 150, h_shelf_max = 150,
                     k_sink = 1,          # Sv per (kg/m3) per million km2 of source area
                     mix = 0.25,          # Sv per million km2: exchange of each source box with the open box
                     mix_v = 10,          # Sv: vertical mixing open box <-> deep box
                     mix_pd = 1,          # Sv per million km2: convective mixing polar box <-> deep box
                     mix_shelf = 0.1,     # Sv per million km2: exchange of the (restricted) evaporative shelf seas with the open sea
                     seaice_effect = 0, seaice_psu = 0.5, seaice_o2 = 0.3,
                     seaice_conv = 1,
                     seaice_auto = FALSE, h_ice = 1.5, h_ml = 50, ice_salinity = 5,   # TRUE: the lid's freshening from the seasonal sea ice (h_ice m of ice into an h_ml m mixed layer) instead of seaice_psu
                     ice_edge = FALSE,    # TRUE: polar seas sink only where open in winter (ice edge); the ice-covered part neither sinks nor convects
                     tau_open = 1, tau_polar = 20, tau_shelf = 0.5,     # yr: export of the surface phosphate
                     f_deep = 0.4,        # share of the export remineralised in the deep box (the rest in the surface layer)
                     p_mean = 2.1, o2_p = 170, o2_pal = 1, cc = 0.07,
                     years = 15000, q_min = 1, S0 = 35,
                     deep_layers = 1,     # 2: intermediate (300-1500 m) + abyssal layer: a lighter source can ventilate the intermediate layer while a denser one fills the abyss
                     mix_ia = 10, f_int = 0.67,
                     melt_polarS_Sv = 0, melt_eff = 1,   # ice-sheet melt water into the south polar sea (Sv) and its efficiency (calibration vs cGENIE)
                     fwcorr_polarS_Sv = 0, fwcorr_polarN_Sv = 0, mix_polarN = 1,   # fresh-water corrections of the polar seas (Sv), see SALT_CALIBRATED   # ice-sheet melt water into the south polar sea (Sv) and its efficiency (calibration vs cGENIE)
                     fw_polarS = 1, fw_polarN = 1,  # scaling of the fresh water into the south / north polar sea (sensitivity)
                     south_basins = 1, basin_min_km2 = 1e6,   # >1: the south polar sea as that many basins
                     river_polarS = 1,    # scaling of the river water into the south polar sea (basins)
                     shelf_regions = 4,   # low-latitude (shelf) sources as that many regions (connected pieces); 0 = one shelf box
                     named_regions = c("Baltic", "Laurentia"), region_km = 700, region_depth = 500,   # basins of their own (LOWLAT_REGIONS): all sea within region_km of their sites, shallower than region_depth
                     river_regions = c(Baltic = 1, Laurentia = 1))   # scaling of the river water into the named basins (0 = no runoff: is the basin stratified by runoff?)

#' sea-water density (kg/m3) at the surface, simplified TEOS-10 (Vallis 2017; Roquet et al. 2015)
rho_sw <- function(T, S) 1027 * (1 - 1.67e-4 * (T - 10) - 0.5 * 1.0e-5 * (T - 10)^2 + 7.8e-4 * (S - 35))

# ---- inputs from a climate run ---------------------------------------------------------------
#' the ocean-surface fields of a climate run on the 1-degree grid: monthly SST (C) and sea ice
#' (HadCM3L + SPEEDY response, as on the maps), annual freshwater into the sea (P - E + rivers,
#' mm/day), HadCM3L surface salinity, sea depth (m) and land fraction of the map
ocean_inputs <- function(run, cfg = SILUCARB) {
  if (is.character(run)) { r <- climate_runs(cfg, status = FALSE); run <- r[r$id == run, , drop = FALSE] }
  if (!nrow(run)) stop("unknown climate run", call. = FALSE)
  s <- run_setup(run, cfg); m <- s$map; x <- s$x; age <- m$age
  cl <- suppressMessages(climate_for_map(m, "ann", "rerun", cfg, flag = FALSE, experiment = if (.is_reference(x, age)) NULL else x))
  cl <- freshwater_fields(cl, m$file)
  mf <- climate_month_files(age, cfg)
  r_sst <- .speedy_monthly_response(m, if (.is_reference(x, age)) NULL else x, "sst", cfg)
  r_ice <- .speedy_monthly_response(m, if (.is_reference(x, age)) NULL else x, "sea_ice", cfg)
  sst <- ice <- array(NA_real_, c(180, 360, 12))
  for (k in 1:12) {
    t <- .read_clim_var(mf[k], "temp_mm_srf"); t$z <- t$z - 273.15
    ic <- .read_clim_var(mf[k], "iceconc_mm_srf"); ic$z[!is.finite(ic$z)] <- 0
    a <- .to_1deg(t); b <- .to_1deg(ic)
    if (!is.null(r_sst)) a <- a + .to_1deg(.speedy_out_field(r_sst[, , k]))
    if (!is.null(r_ice)) b <- b + .to_1deg(.speedy_out_field(r_ice[, , k]))
    sst[, , k] <- a; ice[, , k] <- pmin(pmax(b, 0), 1)
  }
  dem <- read_dem_nc(m$file)
  ii <- pmin(pmax(floor(dem$lat + 90) + 1, 1), 180); jj <- pmin(pmax(floor(((dem$lon + 180) %% 360)) + 1, 1), 360)
  agg <- function(M) { a <- rowsum(M, ii, reorder = TRUE); t(rowsum(t(a), jj, reorder = TRUE)) }
  nsea <- agg((dem$z < 0) * 1); depth <- agg(pmin(dem$z, 0)) / pmax(nsea, 1); n <- agg(matrix(1, nrow(dem$z), ncol(dem$z)))
  xr <- .resolve(x, age)
  list(age = age, run = run$id, label = xr$label, map_file = m$file, dT = xr$dT, gradient = xr$gradient %||% 1,
       sensitivity = x$sensitivity, co2_pal = xr$co2_pal_eff, sun_age = x$sun_age,
       sst = sst, ice = ice, fw = .to_1deg(cl$freshwater), sal = .to_1deg(cl$salinity), w = .to_1deg(cl$w_ocean), riv = .to_1deg(cl$river_input), depth = -depth, lf = 1 - nsea / n,
       river_km3 = cl$rivers$total_km3) |> ocean_inputs_regions()
}

# ---- named low-latitude basins ---------------------------------------------------------------
#' present-day positions of sites in basins that get a box of their own; reconstructed to the map age
#' with PALEOMAP (Scotese & Wright 2018), as the PaleoDEM maps
LOWLAT_REGIONS <- list(
  Baltic = data.frame(name = c("Gotland", "Saaremaa", "Tallinn", "Riga", "Lithuania", "Kaliningrad", "Scania"),
                      lat = c(57.5, 58.4, 59.4, 56.95, 55.3, 54.7, 55.6), lon = c(18.5, 22.5, 24.7, 24.1, 23.9, 20.5, 13.4)),
  Laurentia = data.frame(name = c("Anticosti", "Niagara", "Iowa", "Copenhagen Canyon"),
                         lat = c(49.5, 43.1, 42.0, 39.6), lon = c(-63.0, -79.0, -93.5, -116.1)))
#' add the palaeo-positions of the LOWLAT_REGIONS sites to ocean inputs (needs pyGPlates)
ocean_inputs_regions <- function(inp) {
  inp$regions <- tryCatch(lapply(LOWLAT_REGIONS, function(p) { r <- gp_reconstruct(p$lat, p$lon, rep(inp$age, nrow(p)), "scotese_and_wright2018")
    cbind(p, plat = r$plat, plon = r$plon) }), error = function(e) { message("reconstructing the regional sites failed: ", conditionMessage(e)); NULL })
  inp
}

# ---- boxes -------------------------------------------------------------------------------
#' the cells of each box: south / north polar sea (poleward of polar_lat), shelf seas (shallower
#' than shelf_depth, equatorward), open sea (the rest, equatorward)
.ocean_cells <- function(inp, par) {
  lat <- matrix(GRID1$lat, 180, 360); A <- area_matrix(GRID1$lat, 360); A <- A / sum(A) * 5.1007e14   # m2
  sea <- inp$lf < 0.5 & is.finite(inp$sst[, , 1]) & is.finite(inp$fw)
  dep <- inp$depth; dep[!is.finite(dep)] <- 0
  box <- matrix(NA_character_, 180, 360)
  box[sea & lat <= -par$polar_lat] <- "polarS"; box[sea & lat >= par$polar_lat] <- "polarN"
  eq <- sea & abs(lat) < par$polar_lat
  # shelf seas that lose water (E > P + rivers: anti-estuarine, Bickert's 'arid') can make dense brine;
  # humid shelf seas (estuarine) export fresher water and belong to the open surface box
  fw <- inp$fw; fw[!is.finite(fw)] <- 0
  box[eq & dep < par$shelf_depth & fw < 0] <- "shelf"; box[eq & (dep >= par$shelf_depth | fw >= 0)] <- "open"
  hum <- which(eq & dep < par$shelf_depth & fw >= 0)
  # low-latitude sources as regions: named basins (all their sea, humid or arid) and connected pieces of the
  # remaining evaporative shelf seas
  rid <- matrix(NA_character_, 180, 360); linfo <- NULL
  if (par$shelf_depth > 0 && isTRUE((par$shelf_regions %||% 0) > 0)) {
    lonm <- matrix(GRID1$lon, 180, 360, byrow = TRUE)
    for (nm in intersect(par$named_regions %||% character(0), names(inp$regions %||% list()))) {
      P <- inp$regions[[nm]]; P <- P[is.finite(P$plat), ]; if (!nrow(P)) next
      d <- Reduce(pmin, lapply(seq_len(nrow(P)), function(i) .gc_km(lat, lonm, P$plat[i], P$plon[i])))
      sel <- eq & d <= par$region_km & dep < par$region_depth & is.na(rid)
      if (!any(sel)) next
      rid[sel] <- nm; box[sel] <- "shelf"
      linfo <- rbind(linfo, data.frame(source = paste0("shelf_", nm), label = sprintf("%s basin", nm), sites = paste(P$name, collapse = ", "), stringsAsFactors = FALSE))
    }
    rest <- box == "shelf" & is.na(rid); rest[is.na(rest)] <- FALSE
    if (any(rest)) {
      cid <- .connected_ids(rest)
      area <- tapply(A[rest] / 1e12, cid[rest], sum); keep <- as.integer(names(sort(area, decreasing = TRUE)))[seq_len(min(par$shelf_regions, length(area)))]
      keep <- keep[area[as.character(keep)] >= 0.2 | keep == keep[1]]
      kc <- which(rest & cid %in% keep); kr <- (kc - 1) %% 180 + 1; kcc <- (kc - 1) %/% 180 + 1
      for (c1 in which(rest & !cid %in% keep)) { r <- (c1 - 1) %% 180 + 1; cc <- (c1 - 1) %/% 180 + 1
        dd <- (kr - r)^2 + (((kcc - cc + 180) %% 360) - 180)^2 * cos(GRID1$lat[r] * pi / 180)^2; cid[c1] <- cid[kc[which.min(dd)]] }
      j <- 0
      for (k0 in keep) { j <- j + 1; sel <- rest & cid == k0
        cy <- mean(lat[sel]); cx <- atan2(mean(sin(lonm[sel] * pi / 180)), mean(cos(lonm[sel] * pi / 180))) * 180 / pi
        nm <- sprintf("R%d", j); rid[sel] <- nm
        linfo <- rbind(linfo, data.frame(source = paste0("shelf_", nm), label = sprintf("shelf seas %s (%.0f\u00b0%s, %s)", nm, abs(cy), if (cy < 0) "S" else "N", .lonlab(cx)), sites = "", stringsAsFactors = FALSE)) }
    }
  }
  out <- lapply(c(polarS = "polarS", polarN = "polarN", shelf = "shelf", open = "open"), function(b) {
    k <- which(box == b)
    list(k = k, row = (k - 1) %% 180 + 1, A = A[k], lat = lat[k], depth = dep[k],
         sst = matrix(inp$sst, ncol = 12)[k, , drop = FALSE], ice = matrix(inp$ice, ncol = 12)[k, , drop = FALSE],
         fw = inp$fw[k], sal = inp$sal[k], riv = if (!is.null(inp$riv)) inp$riv[k] else NA * k) })
  hum <- setdiff(hum, which(!is.na(rid)))
  attr(out, "humid_shelf") <- c(A = sum(A[hum]), F = sum(A[hum] * fw[hum]) / 1000 / 86400)
  if (!is.null(linfo)) {
    C <- out$shelf; id <- rid[C$k]
    for (i in seq_len(nrow(linfo))) { nm <- sub("^shelf_", "", linfo$source[i]); sel <- which(id == nm)
      if (length(sel)) out[[linfo$source[i]]] <- lapply(C, function(v) if (is.matrix(v)) v[sel, , drop = FALSE] else v[sel]) }
    linfo$area_Mkm2 <- vapply(linfo$source, function(b) if (is.null(out[[b]])) 0 else sum(out[[b]]$A) / 1e12, 0)
    linfo$lat_centre <- vapply(linfo$source, function(b) if (is.null(out[[b]])) NA else mean(out[[b]]$lat), 0)
    attr(out, "lowlat") <- linfo[linfo$area_Mkm2 > 0, ]
  }
  # the south polar sea split into basins (sources of their own)
  if (isTRUE((par$south_basins %||% 1) > 1) && length(out$polarS$k)) {
    bid <- .south_basin_ids(box == "polarS", inp$lf, A, par$south_basins, par$polar_lat, par$basin_min_km2 %||% 1e6)
    C <- out$polarS; id <- bid[C$k]; info <- attr(bid, "info")
    for (i in sort(unique(id))) { sel <- which(id == i)
      out[[sprintf("polarS_%d", i)]] <- lapply(C, function(v) if (is.matrix(v)) v[sel, , drop = FALSE] else v[sel]) }
    attr(out, "basins") <- info
  }
  out
}
#' basins of the south polar sea: the sea poleward of polar_lat is split at the n meridians where land
#' reaches furthest into it (peninsulas, island chains: natural barriers) and by connectivity; pieces
#' smaller than min_km2 join their nearest basin. Returns a 180 x 360 matrix of basin numbers (west to
#' east) with attr "info" (basin, lon range, centre, area)
.south_basin_ids <- function(M, lf, A, n, polar_lat, min_km2 = 1e6) {
  lat <- GRID1$lat; lon <- GRID1$lon; M[is.na(M)] <- FALSE
  band <- lat <= -polar_lat + 0 & lat >= -polar_lat - 25
  landlon <- colMeans(lf[band, , drop = FALSE] >= 0.5)
  sm <- stats::filter(c(tail(landlon, 7), landlon, head(landlon, 7)), rep(1 / 15, 15))[8:367]   # smoothed, wrapping
  cut <- integer(0); w <- as.numeric(sm); sep <- max(25, 180 / n)
  for (j in seq_len(n)) { m <- which.max(w); if (!is.finite(w[m]) || w[m] < 0) break; cut <- c(cut, m)
    d <- abs(((seq_len(360) - m + 180) %% 360) - 180); w[d < sep] <- -Inf }
  cut <- sort(cut)
  sector <- function(j) { if (length(cut) < 2) return(1L); k <- findInterval(j, cut); if (k == 0) length(cut) else k }
  # connected pieces within each sector (8-neighbours, longitude wraps)
  id <- matrix(0L, 180, 360); sec <- matrix(vapply(rep(1:360, each = 180), sector, 1L), 180, 360)
  cells <- which(M); nid <- 0L
  for (c0 in cells) { if (id[c0] > 0) next
    nid <- nid + 1L; queue <- c0; id[c0] <- nid
    while (length(queue)) { c1 <- queue[1]; queue <- queue[-1]
      r <- (c1 - 1) %% 180 + 1; cc <- (c1 - 1) %/% 180 + 1
      for (dr in -1:1) for (dc in -1:1) { r2 <- r + dr; if (r2 < 1 || r2 > 180) next; c2 <- ((cc + dc - 1) %% 360) + 1; k2 <- r2 + (c2 - 1) * 180
        if (M[k2] && id[k2] == 0L && sec[k2] == sec[c1]) { id[k2] <- nid; queue <- c(queue, k2) } } } }
  # small pieces join the nearest large one
  area <- tapply(A[cells] / 1e6, id[cells], sum); big <- as.integer(names(area)[area >= min_km2])
  if (!length(big)) big <- as.integer(names(which.max(area)))
  bigc <- cells[id[cells] %in% big]; br <- (bigc - 1) %% 180 + 1; bc <- (bigc - 1) %/% 180 + 1
  for (c1 in cells[!id[cells] %in% big]) { r <- (c1 - 1) %% 180 + 1; cc <- (c1 - 1) %/% 180 + 1
    d <- (br - r)^2 + (((bc - cc + 180) %% 360) - 180)^2 * cos(lat[r] * pi / 180)^2; id[c1] <- id[bigc[which.min(d)]] }
  # number west to east by the circular mean longitude
  cl <- tapply(lon[(cells - 1) %/% 180 + 1], id[cells], function(x) atan2(mean(sin(x * pi / 180)), mean(cos(x * pi / 180))) * 180 / pi)
  o <- order(cl); newid <- match(as.integer(names(cl)), as.integer(names(cl))[o]); map <- setNames(newid, names(cl))
  id[cells] <- map[as.character(id[cells])]
  info <- do.call(rbind, lapply(sort(unique(id[cells])), function(i) { k <- cells[id[cells] == i]; x <- lon[(k - 1) %/% 180 + 1]
    xs <- sort(unique(x)); gaps <- diff(c(xs, xs[1] + 360)); g <- which.max(gaps)
    from <- xs[(g %% length(xs)) + 1]; to <- xs[g]
    data.frame(basin = sprintf("polarS_%d", i), lon_from = from, lon_to = to, lon_centre = cl[[names(map)[map == i][1]]],
               lat_centre = mean(lat[(k - 1) %% 180 + 1]), area_Mkm2 = sum(A[k]) / 1e12,
               label = sprintf("S%d (%s to %s)", i, .lonlab(from - 0.5), .lonlab(to + 0.5))) }))
  attr(id, "info") <- info; id
}
.lonlab <- function(x) { x <- ((x + 180) %% 360) - 180; sprintf("%.0f\u00b0%s", abs(x), if (x < 0) "W" else "E") }
#' box values (vectors over configurations) for SST shifts dmon (list over configs of 12 x 180
#' matrices: K by month and latitude row) and freshwater scaling
.ocean_box_values <- function(cells, par, dmon = list(NULL), sea_level = 0) {
  tf <- -1.8
  out <- list()
  for (b in names(cells)) {
    C <- cells[[b]]; nb <- length(C$k)
    A <- C$A

    res <- lapply(seq_along(dmon), function(j) {
      D <- if (is.null(dmon[[j]])) matrix(0, nb, 12) else t(dmon[[j]][, C$row, drop = FALSE])
      raw <- C$sst + D
      Tc <- pmax(raw, tf)
      ic <- pmax(C$ice * pmax(0, 1 - 0.25 * pmax(D, 0)), pmin(1, pmax(0, (tf - (pmax(C$sst, tf) + D)) / 2)))
      Aj <- C$A
      sl <- if (length(sea_level) >= j) sea_level[j] else sea_level[1]
      if (startsWith(b, "shelf") && sl < 0) Aj <- C$A * (C$depth > -sl)  # shelf falls dry when the sea is lower
      polar <- startsWith(b, "polar"); mo <- if (startsWith(b, "polarS")) 6:8 else if (b == "polarN") c(12, 1, 2) else 1:12
      shf <- startsWith(b, "shelf"); wm <- if (shf) ifelse(C$lat < 0, 1, 0) else NULL
      Tw <- if (shf) rowMeans(Tc[, 6:8, drop = FALSE]) * wm + rowMeans(Tc[, c(12, 1, 2), drop = FALSE]) * (1 - wm) else rowMeans(Tc[, mo, drop = FALSE])
      fi <- if (shf) rowMeans(ic[, 6:8, drop = FALSE]) * wm + rowMeans(ic[, c(12, 1, 2), drop = FALSE]) * (1 - wm) else rowMeans(ic[, mo, drop = FALSE])
      Aall <- Aj; fice <- sum(Aall * fi) / max(sum(Aall), 1)
      ms <- if (startsWith(b, "polarS")) c(12, 1, 2) else 6:8                          # summer months (polar seas)
      fsum <- if (polar) sum(Aall * rowMeans(ic[, ms, drop = FALSE])) / max(sum(Aall), 1) else 0   # summer sea ice (perennial part)
      # ice-edge mode: a polar sea sinks only where it is open in winter (sea ice < 50 %), i.e. at the
      # ice edge; the ice-covered part does not sink or convect, its fresh water goes to the open part
      if (isTRUE(par$ice_edge) && polar) Aj <- Aj * (fi < 0.5)
      W <- Aj / max(sum(Aj), 1); if (sum(Aj) <= 0) W <- Aall / max(sum(Aall), 1)
      dTa <- rowMeans(D)
      Rv <- if (all(is.na(C$riv))) NA else sum(Aall * pmax(C$riv, 0), na.rm = TRUE) / 1000 / 86400 * (1 + par$cc * mean(dTa))
      rsc <- if (startsWith(b, "polarS")) par$river_polarS %||% 1 else if (startsWith(b, "shelf_") && sub("shelf_", "", b) %in% names(par$river_regions)) par$river_regions[[sub("shelf_", "", b)]] else 1
      rs <- if (is.finite(Rv)) (1 - rsc) * Rv else 0                                       # less (more) river water
      c(A = sum(Aj), T = sum(W * Tw), Tann = sum(W * rowMeans(Tc)), ice = fice, ice_summer = fsum, dTann = sum(W * dTa),
        F = (sum(Aall * C$fw * (1 + par$cc * dTa)) / 1000 / 86400 - rs) * # m3/s, + = freshening
            (if (startsWith(b, "polarS")) par$fw_polarS %||% 1 else if (b == "polarN") par$fw_polarN %||% 1 else 1),
        R = Rv * rsc,                                                     # of which river water (m3/s)
        Sobs = sum(W * C$sal, na.rm = TRUE) / max(sum(W[is.finite(C$sal)]), 1e-12),
        depth = if (sum(Aj) > 0) sum(W * pmax(C$depth, 10)) else 50)
    })
    out[[b]] <- do.call(rbind, res)
  }
  out
}

# ---- box model ---------------------------------------------------------------------------
#  Boxes: south polar, north polar, evaporative shelf seas (the three possible deep-water
#  sources), open surface sea, deep sea. Surface temperatures follow the climate (winter for the
#  sources); salinity follows the freshwater flux and the exchanges. A source sinks when its water
#  is denser than the deep water: q = k_sink * area * (rho_source - rho_deep) (Stommel 1961 type);
#  the deep water upwells into the open box, which feeds the sources. Mixing: sources <-> open box,
#  polar boxes <-> deep box (convection), open <-> deep (vertical). Phosphate is exported from the
#  surface boxes and partly remineralised in the deep box, which uses O2 (170 mol per mol P).
SRC <- c("polarS", "polarN", "shelf")
#' the deep-water sources of a set of box values: the south polar sea (or its basins polarS_1...),
#' the north polar sea and the evaporative shelf seas
.src_names <- function(bv) { sb <- sort(grep("^polarS_[0-9]+$", names(bv), value = TRUE)); sh <- grep("^shelf_", names(bv), value = TRUE)
  sh <- c(sh[!grepl("^shelf_R[0-9]+$", sh)], sort(grep("^shelf_R[0-9]+$", sh, value = TRUE)))
  c(if (length(sb)) sb else "polarS", "polarN", if (length(sh)) sh else "shelf") }
.box_setup <- function(bv, par) {
  Sv <- 1e6 * 3.15576e7
  src <- .src_names(bv); ns <- length(src); io <- ns + 1; pol <- startsWith(src, "polar"); ish <- which(src == "shelf")
  N <- nrow(bv$open)
  col <- function(v) matrix(vapply(c(src, "open"), function(b) bv[[b]][, v], numeric(N)), N, io)
  A <- col("A"); colnames(A) <- c(src, "open")
  ish <- which(!pol)
  h <- matrix(par$h_polar, N, io); h[, ish] <- pmin(col("depth")[, ish, drop = FALSE], par$h_shelf_max); h[, io] <- par$h_open
  V <- pmax(A * h, 1e12); VD <- pmax(rowSums(A) * 3700 - rowSums(V), 1e16)
  Tsrc <- col("T")[, seq_len(ns), drop = FALSE]; To <- bv$open[, "Tann"]
  fice <- col("ice")[, seq_len(ns), drop = FALSE]; fice[, !pol] <- 0
  Fw <- col("F") * 3.15576e7
  # melt water of an ice sheet on the south polar continent (Sv, e.g. Gondwana), shared over the south
  # polar box(es) by area; times melt_eff (box model vs cGENIE, see MELT_CALIBRATED)
  # plus a constant fresh-water correction of each polar sea (Sv; flux correction, calibrated so that the
  # polar surface salinities match HadCM3L at 435 Ma, see SALT_CALIBRATED)
  addfw <- function(Fw, Sv, boxes) { if (!any(Sv != 0) || !length(boxes)) return(Fw); Ab <- A[, boxes, drop = FALSE]
    Fw[, boxes] <- Fw[, boxes] + rep_len(Sv, N) * 1e6 * 3.15576e7 * Ab / pmax(rowSums(Ab), 1); Fw }
  ms <- (par$melt_polarS_Sv %||% 0) * (par$melt_eff %||% 1)
  Fw <- addfw(Fw, ms + (par$fwcorr_polarS_Sv %||% 0), which(startsWith(src, "polarS")))
  Fw <- addfw(Fw, par$fwcorr_polarN_Sv %||% 0, which(src == "polarN"))
  Fw[, io] <- Fw[, io] - rowSums(Fw)                                                     # P - E + R of the whole sea = 0
  As <- A[, seq_len(ns), drop = FALSE] / 1e12
  P1 <- matrix(pol, N, ns, byrow = TRUE)
  kq <- par$k_sink * As * Sv
  mx <- As * ifelse(P1, par$mix, par$mix_shelf) * Sv
  mx[, src == "polarN"] <- mx[, src == "polarN"] * (par$mix_polarN %||% 1)    # exchange of the north polar sea with the open sea (x)
  # polar convection; a melt-water lid (seaice_effect < 0) shuts it in proportion to the winter sea ice
  # (Pohl et al. 2021: seasonal sea ice and its melt water stratify the polar sea and cut its ventilation)
  lid <- pmin(1, (par$seaice_conv %||% 1) * max(-par$seaice_effect, 0) * fice)
  if (isTRUE(par$ice_edge)) lid <- 0 * lid                                          # the ice-covered part is already left out
  mpd <- par$mix_pd * As * (1 - lid) * P1 * Sv
  sie <- par$seaice_effect * par$seaice_psu * fice * P1
  if (isTRUE(par$seaice_auto)) {
    # freshening computed from the seasonal sea ice: the ice that forms each winter and melts each summer
    # (winter minus summer cover, thickness h_ice) mixed into the polar mixed layer (h_ml); ice keeps
    # ice_salinity of the sea water's salt. Full seasonal cover, 1.5 m into 50 m: about 0.9 psu.
    fs <- tryCatch(col("ice_summer")[, seq_len(ns), drop = FALSE], error = function(e) 0 * fice); fs[, !pol] <- 0
    psu <- par$S0 * (1 - (par$ice_salinity %||% 5) / par$S0) * (par$h_ice %||% 1.5) / (par$h_ml %||% 50)
    sie <- par$seaice_effect * psu * pmax(fice - fs, 0) * P1
  }
  ofac <- 1 - par$seaice_o2 * fice * P1                                         # less gas exchange under sea ice
  osat <- o2_saturation(Tsrc, par$S0, par$o2_pal) * ofac                         # (updated with the box salinity in the runs)
  osat_o <- o2_saturation(To, par$S0, par$o2_pal)
  tau <- matrix(c(ifelse(pol, par$tau_polar, par$tau_shelf), par$tau_open), N, io, byrow = TRUE)
  rate <- max((3 * kq + mx + mpd) / V[, seq_len(ns), drop = FALSE], (par$mix_v * Sv + 3 * rowSums(kq) + rowSums(mx)) / V[, io], (par$mix_v * Sv + 3 * rowSums(kq)) / VD * 10)
  list(N = N, src = src, ns = ns, io = io, pol = pol, ish = ish, south = startsWith(src, "polarS"),
       A = A, V = V, VD = VD, Tsrc = Tsrc, To = To, fice = fice, Fw = Fw, kq = kq, mx = mx, mpd = mpd,
       mv = par$mix_v * Sv, sie = sie, osat = osat, osat_o = osat_o, ofac = ofac, o2_pal = par$o2_pal, tau = tau, dt = min(1, 0.4 / rate), Sv = Sv)
}
.box_init <- function(bv, par, start = c("thermal", "haline"), k = .box_setup(bv, par)) {
  start <- match.arg(start); N <- k$N
  S <- matrix(par$S0, N, k$io); P <- matrix(par$p_mean, N, k$io)
  if (start == "thermal") { TD <- apply(k$Tsrc[, k$pol, drop = FALSE], 1, min); SD <- rep(par$S0, N) }
  else { S[, k$ish] <- par$S0 + 2; TD <- rowSums(k$Tsrc[, k$ish, drop = FALSE] * k$A[, k$ish, drop = FALSE]) / pmax(rowSums(k$A[, k$ish, drop = FALSE]), 1); SD <- rep(par$S0 + 1, N) }
  list(S = S, P = P, TD = TD, SD = SD, PD = rep_len(par$p_mean, N), OD = rep(150, N), q = matrix(0, N, k$ns), EPt = rep(0, N), deficit = rep(FALSE, N))
}
.box_run <- function(st, k, par, years, keep = FALSE) {
  nstep <- max(1, ceiling(years / k$dt)); dt <- years / nstep
  S <- st$S; P <- st$P; TD <- st$TD; SD <- st$SD; PD <- st$PD; OD <- st$OD
  is <- seq_len(k$ns); io <- k$io
  hist <- if (keep) list() else NULL
  for (it in seq_len(nstep)) {
    rD <- rho_sw(TD, SD)
    q <- k$kq * pmin(pmax(rho_sw(k$Tsrc, S[, is, drop = FALSE] + k$sie) - rD, 0), 3); qt <- rowSums(q)
    ex <- q + k$mx; So <- S[, io]; Po <- P[, io]
    if (it %% 20 == 1) { k$osat <- o2_saturation(k$Tsrc, S[, is, drop = FALSE], k$o2_pal) * k$ofac; k$osat_o <- o2_saturation(k$To, So, k$o2_pal) }   # warm, saline water holds less O2
    dS <- cbind((ex * (So - S[, is, drop = FALSE]) + k$mpd * (SD - S[, is, drop = FALSE]) - par$S0 * k$Fw[, is, drop = FALSE]) / k$V[, is, drop = FALSE],
                (qt * (SD - So) + rowSums(k$mx * (S[, is, drop = FALSE] - So)) + k$mv * (SD - So) - par$S0 * k$Fw[, io]) / k$V[, io])
    dP <- cbind((ex * (Po - P[, is, drop = FALSE]) + k$mpd * (PD - P[, is, drop = FALSE])) / k$V[, is, drop = FALSE],
                (qt * (PD - Po) + rowSums(k$mx * (P[, is, drop = FALSE] - Po)) + k$mv * (PD - Po)) / k$V[, io])
    Pn <- pmax(P + dt * dP, 0); Pe <- Pn / (1 + dt * par$f_deep / k$tau)               # export (implicit)
    EPt <- rowSums(k$V * (Pn - Pe)) / dt
    dSD <- (rowSums(q * S[, is, drop = FALSE]) - qt * SD + k$mv * (So - SD) + rowSums(k$mpd * (S[, is, drop = FALSE] - SD))) / k$VD
    dTD <- (rowSums(q * k$Tsrc) - qt * TD + k$mv * (k$To - TD) + rowSums(k$mpd * (k$Tsrc - TD))) / k$VD
    dPD <- (rowSums(q * P[, is, drop = FALSE]) - qt * PD + k$mv * (Po - PD) + rowSums(k$mpd * (P[, is, drop = FALSE] - PD)) + EPt) / k$VD
    dOD <- (rowSums(q * k$osat) - qt * OD + k$mv * (k$osat_o - OD) + rowSums(k$mpd * (k$osat - OD)) - par$o2_p * EPt) / k$VD
    S <- S + dt * dS; P <- Pe; SD <- SD + dt * dSD; TD <- TD + dt * dTD; PD <- PD + dt * dPD; OD <- OD + dt * dOD
    deficit <- OD < 0; OD <- pmax(OD, 0)
    if (keep && it %% max(1, round(50 / dt)) == 0)
      hist[[length(hist) + 1]] <- data.frame(year = it * dt, q_polarS = sum(q[1, k$south]) / k$Sv, q_polarN = q[1, k$src == "polarN"] / k$Sv,
                                             q_shelf = sum(q[1, k$ish]) / k$Sv, O2_deep = OD[1], T_deep = TD[1], S_deep = SD[1])
  }
  out <- list(S = S, P = P, TD = TD, SD = SD, PD = PD, OD = OD, q = q / k$Sv, EPt = EPt, deficit = deficit)
  if (keep) attr(out, "history") <- do.call(rbind, hist)
  out
}
#' per-source columns (q_, rho_, S_, layer_, ...) and the south polar sea as a whole (sum of its basins)
.src_columns <- function(k, fo, rs, S, layer = NULL, sink = NULL, conv = NULL) {
  out <- list()
  for (j in seq_len(k$ns)) { b <- k$src[j]
    out[[paste0("q_", b)]] <- fo[, j]; out[[paste0("rho_", b)]] <- rs[, j]; out[[paste0("S_", b)]] <- S[, j]
    if (!is.null(layer)) out[[paste0("layer_", b)]] <- layer[, j]
    if (!is.null(sink)) out[[paste0("sink_", b)]] <- sink[, j]
    if (!is.null(conv)) out[[paste0("conv_", b)]] <- conv[, j] }
  agg <- function(sel, nm) { w <- k$A[, seq_len(k$ns), drop = FALSE][, sel, drop = FALSE]; ws <- pmax(rowSums(w), 1)
    out[[paste0("q_", nm)]] <<- rowSums(fo[, sel, drop = FALSE]); out[[paste0("rho_", nm)]] <<- rowSums(rs[, sel, drop = FALSE] * w) / ws
    out[[paste0("S_", nm)]] <<- rowSums(S[, sel, drop = FALSE] * w) / ws
    if (!is.null(layer)) { L <- layer[, sel, drop = FALSE]; F <- fo[, sel, drop = FALSE]
      out[[paste0("layer_", nm)]] <<- vapply(seq_len(k$N), function(i) if (max(F[i, ]) > 0.05) L[i, which.max(F[i, ])] else "none", "") }
    if (!is.null(sink)) out[[paste0("sink_", nm)]] <<- rowSums(sink[, sel, drop = FALSE])
    if (!is.null(conv)) out[[paste0("conv_", nm)]] <<- rowSums(conv[, sel, drop = FALSE]) }
  if (length(k$ish) > 1 || k$src[k$ish[1]] != "shelf") agg(k$ish, "shelf")
  if (sum(k$south) > 1) {                                                      # the south polar sea as a whole
    w <- k$A[, seq_len(k$ns), drop = FALSE][, k$south, drop = FALSE]; ws <- pmax(rowSums(w), 1)
    out$q_polarS <- rowSums(fo[, k$south, drop = FALSE]); out$rho_polarS <- rowSums(rs[, k$south, drop = FALSE] * w) / ws
    out$S_polarS <- rowSums(S[, k$south, drop = FALSE] * w) / ws
    if (!is.null(layer)) { L <- layer[, k$south, drop = FALSE]; F <- fo[, k$south, drop = FALSE]
      out$layer_polarS <- vapply(seq_len(k$N), function(i) if (max(F[i, ]) > 0.05) L[i, which.max(F[i, ])] else "none", "") }
    if (!is.null(sink)) out$sink_polarS <- rowSums(sink[, k$south, drop = FALSE])
    if (!is.null(conv)) out$conv_polarS <- rowSums(conv[, k$south, drop = FALSE])
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
.box_diag <- function(st, k, par) {
  q <- st$q; qt <- rowSums(q); rD <- rho_sw(st$TD, st$SD)
  rs <- rho_sw(k$Tsrc, st$S[, seq_len(k$ns), drop = FALSE] + k$sie) - 1000
  pw <- k$A[, seq_len(k$ns), drop = FALSE] * matrix(k$pol, k$N, k$ns, byrow = TRUE)
  r <- data.frame(q_total = qt, share_polar = ifelse(qt > 0, rowSums(q[, k$pol, drop = FALSE]) / qt, NA), share_shelf = ifelse(qt > 0, rowSums(q[, k$ish, drop = FALSE]) / qt, NA),
                  T_deep = st$TD, S_deep = st$SD, O2_deep = st$OD, P_deep = st$PD, P_surface = st$P[, k$io],
                  export_Tmol = st$EPt * 1025e-6 / 1e12, anoxic = st$deficit, rho_deep = rD - 1000, S_open = st$S[, k$io],
                  T_polar_winter = rowSums(k$Tsrc * pw) / pmax(rowSums(pw), 1), shelf_area = rowSums(k$A[, k$ish, drop = FALSE]) / 1e12,
                  vent_age_yr = k$VD / (qt * k$Sv + rowSums(k$mpd) + k$mv))                     # deep box volume / renewal
  r <- cbind(.src_columns(k, q, rs, st$S[, seq_len(k$ns), drop = FALSE]), r)
  r$state <- .ocean_class(r, par); r
}
#' steady states of the box model for box values bv (from .ocean_box_values; one row per
#' configuration), from a polar ("thermal") or a shelf ("haline") start
ocean_boxes <- function(bv, par = OCEAN_PARAMS, years = par$years, start = c("thermal", "haline"), keep = FALSE) {
  if (isTRUE((par$deep_layers %||% 1) >= 2)) return(.ocean_boxes2(bv, par, years, match.arg(start)))
  k <- .box_setup(bv, par); st <- .box_init(bv, par, start, k); st <- .box_run(st, k, par, years, keep)
  out <- .box_diag(st, k, par); if (keep) attr(out, "history") <- attr(st, "history"); out
}
#' two deep layers: intermediate (300-1500 m) and abyss. A source sinks into the abyss when it is
#' denser than the abyssal water, else into the intermediate layer when denser than that; the abyss
#' upwells into the intermediate layer, which upwells into the open surface box. Polar convection
#' mixes a polar box with the layer it reaches. So two sources can sink at once (e.g. a southern
#' intermediate and a northern abyssal cell, as in Pohl et al. 2021 Fig. 4i).
.ocean_boxes2 <- function(bv, par, years, start = "thermal", st = NULL, return_state = FALSE) {
  k <- .box_setup(bv, par); N <- k$N; Sv <- k$Sv; is <- seq_len(k$ns); io <- k$io
  VI <- pmin(rowSums(k$A) * 1200, 0.4 * k$VD); VA <- k$VD - VI
  mia <- par$mix_ia * Sv; fI <- par$f_int
  if (is.null(st)) {
    S <- matrix(par$S0, N, io); P <- matrix(par$p_mean, N, io)
    if (start == "thermal") { TA <- apply(k$Tsrc[, k$pol, drop = FALSE], 1, min); SA <- rep(par$S0, N) } else {
      S[, k$ish] <- par$S0 + 2; TA <- rowSums(k$Tsrc[, k$ish, drop = FALSE] * k$A[, k$ish, drop = FALSE]) / pmax(rowSums(k$A[, k$ish, drop = FALSE]), 1); SA <- rep(par$S0 + 1, N) }
    st <- list(S = S, P = P, TI = (TA + k$To) / 2, SI = rep(par$S0, N), PI = rep_len(par$p_mean, N), OI = rep(150, N),
               TA = TA, SA = SA, PA = rep_len(par$p_mean, N), OA = rep(150, N))
  }
  S <- st$S; P <- st$P; TI <- st$TI; SI <- st$SI; PI <- st$PI; OI <- st$OI; TA <- st$TA; SA <- st$SA; PA <- st$PA; OA <- st$OA
  rate <- max(k$dt^-1 * 0.4, (mia + 3 * rowSums(k$kq)) / VI * 10)
  dt <- min(1, 0.4 / rate); nstep <- max(1, ceiling(years / dt)); dt <- years / nstep
  acc_n <- 0; acc_qA <- acc_qI <- acc_mA <- acc_mI <- 0 * k$kq
  for (it in seq_len(nstep)) {
    rI <- rho_sw(TI, SI); rA <- rho_sw(TA, SA); rs <- rho_sw(k$Tsrc, S[, is, drop = FALSE] + k$sie)
    toA <- pmin(pmax((rs - rA) / 0.04 + 0.5, 0), 1)                         # share going to the abyss (smooth)
    qA <- k$kq * pmin(pmax(rs - rA, 0), 3); qI <- k$kq * pmin(pmax(rs - rI, 0), 3) * (1 - toA)
    qAt <- rowSums(qA); qIt <- rowSums(qI); qt <- qAt + qIt
    # polar convection only where the surface water is (nearly) as dense as the water below it
    cf <- pmin(pmax((rs - (toA * rA + (1 - toA) * rI) + 0.05) / 0.1, 0), 1)
    mA <- k$mpd * toA * cf; mI <- k$mpd * (1 - toA) * cf
    if (it > nstep * 0.9) { acc_n <- acc_n + 1; acc_qA <- acc_qA + qA; acc_qI <- acc_qI + qI; acc_mA <- acc_mA + mA; acc_mI <- acc_mI + mI }
    ex <- qA + qI + k$mx; So <- S[, io]; Po <- P[, io]; Ss <- S[, is, drop = FALSE]; Ps <- P[, is, drop = FALSE]
    if (it %% 20 == 1) { k$osat <- o2_saturation(k$Tsrc, Ss, k$o2_pal) * k$ofac; k$osat_o <- o2_saturation(k$To, So, k$o2_pal) }
    dS <- cbind((ex * (So - Ss) + mI * (SI - Ss) + mA * (SA - Ss) - par$S0 * k$Fw[, is, drop = FALSE]) / k$V[, is, drop = FALSE],
                (qt * (SI - So) + rowSums(k$mx * (Ss - So)) + k$mv * (SI - So) - par$S0 * k$Fw[, io]) / k$V[, io])
    dP <- cbind((ex * (Po - Ps) + mI * (PI - Ps) + mA * (PA - Ps)) / k$V[, is, drop = FALSE],
                (qt * (PI - Po) + rowSums(k$mx * (Ps - Po)) + k$mv * (PI - Po)) / k$V[, io])
    Pn <- pmax(P + dt * dP, 0); Pe <- Pn / (1 + dt * par$f_deep / k$tau); EPt <- rowSums(k$V * (Pn - Pe)) / dt
    lay <- function(XS, XI, XA, Xo) list(
      I = (rowSums(qI * XS) + qAt * XA - qt * XI + k$mv * (Xo - XI) + mia * (XA - XI) + rowSums(mI * (XS - XI))) / VI,
      A = (rowSums(qA * XS) - qAt * XA + mia * (XI - XA) + rowSums(mA * (XS - XA))) / VA)
    dSL <- lay(Ss, SI, SA, So); dTL <- lay(k$Tsrc, TI, TA, k$To)
    dPL <- lay(Ps, PI, PA, Po); dOL <- lay(k$osat, OI, OA, k$osat_o)
    S <- S + dt * dS; P <- Pe
    SI <- SI + dt * dSL$I; SA <- SA + dt * dSL$A; TI <- TI + dt * dTL$I; TA <- TA + dt * dTL$A
    PI <- PI + dt * (dPL$I + fI * EPt / VI); PA <- PA + dt * (dPL$A + (1 - fI) * EPt / VA)
    OI <- pmax(OI + dt * (dOL$I - par$o2_p * fI * EPt / VI), 0); OA <- pmax(OA + dt * (dOL$A - par$o2_p * (1 - fI) * EPt / VA), 0)
  }
  st <- list(S = S, P = P, TI = TI, SI = SI, PI = PI, OI = OI, TA = TA, SA = SA, PA = PA, OA = OA)
  if (return_state) return(st)
  if (acc_n > 0) { qA <- acc_qA / acc_n; qI <- acc_qI / acc_n; mA <- acc_mA / acc_n; mI <- acc_mI / acc_n; qAt <- rowSums(qA); qIt <- rowSums(qI) }   # mean of the last 10 %
  q <- (qA + qI) / Sv; rs <- rho_sw(k$Tsrc, S[, is, drop = FALSE] + k$sie) - 1000
  conv <- (mA + mI) / Sv                                                    # convective renewal (Sv) per source
  layer <- matrix(ifelse(qA + mA + qI + mI < 0.05 * Sv, "none", ifelse(qA + mA >= qI + mI, "abyss", "intermediate")), N, k$ns)
  fo <- q + conv; ft <- rowSums(fo)                                        # deep-water formation = sinking + convective renewal
  pw <- k$A[, is, drop = FALSE] * matrix(k$pol, N, k$ns, byrow = TRUE)
  r <- data.frame(q_total = ft, share_polar = ifelse(ft > 0, rowSums(fo[, k$pol, drop = FALSE]) / ft, NA), share_shelf = ifelse(ft > 0, rowSums(fo[, k$ish, drop = FALSE]) / ft, NA),
                  T_deep = TA, S_deep = SA, O2_deep = OA, P_deep = PA, P_surface = P[, io], export_Tmol = EPt * 1025e-6 / 1e12, anoxic = OA < 1,
                  rho_deep = rho_sw(TA, SA) - 1000, S_open = S[, io],
                  T_polar_winter = rowSums(k$Tsrc * pw) / pmax(rowSums(pw), 1), shelf_area = rowSums(k$A[, k$ish, drop = FALSE]) / 1e12,
                  vent_age_yr = VA / (qAt + mia + rowSums(mA)),
                  T_int = TI, S_int = SI, O2_int = OI, rho_int = rho_sw(TI, SI) - 1000, q_abyss = qAt / Sv, q_int = qIt / Sv, stringsAsFactors = FALSE)
  r <- cbind(.src_columns(k, fo, rs, S[, is, drop = FALSE], layer, q, conv), r)
  r$state <- .ocean_class(r, par)
  r
}
.ocean_class <- function(r, par) {
  ifelse(r$q_total < par$q_min, "stagnant",
         ifelse(r$share_polar >= 2 / 3, "P-like (thermal)", ifelse(r$share_shelf >= 2 / 3, "S-like (haline)", "mixed")))
}
OCEAN_STATE_COLS <- c("P-like (thermal)" = "#2c7fb8", "mixed" = "#9e9ac8", "S-like (haline)" = "#e6550d", "stagnant" = "#969696")

#' the ocean state of one climate run: box values and the steady states from both starts
ocean_state <- function(run, par = list(), cfg = SILUCARB, inputs = NULL) {
  par <- modifyList(OCEAN_PARAMS, par)
  inp <- inputs %||% ocean_inputs(run, cfg)
  cells <- .ocean_cells(inp, par); bv <- .ocean_box_values(cells, par)
  a <- ocean_boxes(bv, par, start = "thermal", keep = TRUE); b <- ocean_boxes(bv, par, start = "haline", keep = TRUE)
  list(inputs = inp, par = par, boxes = bv, thermal = a, haline = b, bistable = a$state != b$state,
       label = inp$label, humid_shelf = attr(cells, "humid_shelf"), summary = ocean_state_text(bv, a, b, attr(cells, "humid_shelf")))
}
ocean_state_text <- function(bv, a, b, hum = NULL) {
  f <- function(r) sprintf("%s: deep water %.1f Sv (south pole %.1f, north pole %.1f, shelf seas %.1f), deep sea %.1f \u00b0C, %.2f psu, O2 %.0f \u00b5mol/kg%s",
                           r$state, r$q_total, r$q_polarS, r$q_polarN, r$q_shelf, r$T_deep, r$S_deep, r$O2_deep, if (isTRUE(r$anoxic)) " (anoxic: O2 used up)" else "")
  bx <- sprintf("Winter surface water: south polar %.1f \u00b0C (%.0f %% sea ice, %.1f Mkm\u00b2), north polar %.1f \u00b0C (%.0f %% ice, %.1f Mkm\u00b2), evaporative (anti-estuarine) shelf seas %.1f \u00b0C (%.1f Mkm\u00b2), open sea %.1f \u00b0C. Freshwater into the boxes (mSv): S %.0f, N %.0f, shelf %.0f, open %.0f.",
                bv$polarS[1, "T"], 100 * bv$polarS[1, "ice"], bv$polarS[1, "A"] / 1e12, bv$polarN[1, "T"], 100 * bv$polarN[1, "ice"], bv$polarN[1, "A"] / 1e12,
                bv$shelf[1, "T"], bv$shelf[1, "A"] / 1e12, bv$open[1, "Tann"], bv$polarS[1, "F"] / 1e3, bv$polarN[1, "F"] / 1e3, bv$shelf[1, "F"] / 1e3, bv$open[1, "F"] / 1e3)
  if (!is.null(hum)) bx <- paste0(bx, sprintf(" Humid (estuarine) shelf seas: %.1f Mkm\u00b2, %.0f mSv fresh water (part of the open box).", hum[["A"]] / 1e12, hum[["F"]] / 1e3))
  paste0(bx, "\nFrom a polar start \u2192 ", f(a), ".\nFrom a shelf start \u2192 ", f(b), ".",
         if (a$state != b$state) "\nTwo states are possible here (bistable): which one the ocean is in depends on its history." else "")
}

# ---- the CO2 x gradient plane --------------------------------------------------------------------
#' SST shift (K, 12 months x 180 latitude rows) from the base run to CO2 co2_pal and gradient g:
#' as shift_forcing() does for SPEEDY (global shift with polar amplification + gradient change)
.ocean_dmon <- function(inp, co2_pal, g) {
  x <- .resolve(climate_experiment(co2_pal = co2_pal, sun_age = inp$sun_age, sensitivity = inp$sensitivity, gradient = g), inp$age)
  pat <- 0.75 + 0.75 * sin(GRID1$lat * pi / 180)^2
  t(vapply(1:12, function(m) (x$dT - inp$dT) * pat + gradient_anomaly(inp$age, g, GRID1$lat, m) - gradient_anomaly(inp$age, inp$gradient, GRID1$lat, m),
           numeric(180)))
}
#' regime map: ocean state over CO2 (PAL) x gradient factor, from both starts; with the mean
#' annual tropical SST (0-20 deg, for the conodont d18O constraint) and, optionally, the area of
#' land where snow survives the summer (land-ice line)
ocean_regime_map <- function(inp, co2 = 2^seq(log2(2), log2(16), length.out = 15), gradient = seq(0.5, 1.1, 0.05),
                             par = list(), ice_line = FALSE, progress = NULL, cfg = SILUCARB) {
  par <- modifyList(OCEAN_PARAMS, par)
  G <- expand.grid(co2 = co2, gradient = gradient)
  cells <- .ocean_cells(inp, par)
  dm <- lapply(seq_len(nrow(G)), function(i) .ocean_dmon(inp, G$co2[i], G$gradient[i]))
  if (is.function(progress)) progress(0.2, "sea surface of each setting")
  bv <- .ocean_box_values(cells, par, dm)
  # tropical SST (annual, sea 20S-20N)
  lat <- GRID1$lat; tr <- abs(lat) <= 20; A <- area_matrix(lat, 360); sea <- inp$lf < 0.5 & is.finite(inp$sst[, , 1])
  ann <- apply(pmax(inp$sst, -1.8), 1:2, mean)
  G$sst_tropics <- vapply(dm, function(D) { Z <- ann + matrix(colMeans(D), 180, 360); sel <- sea & matrix(tr, 180, 360)
    sum((A * Z)[sel]) / sum(A[sel]) }, 0)
  G$sst_polar_winter <- (bv$polarS[, "T"] * bv$polarS[, "A"] + bv$polarN[, "T"] * bv$polarN[, "A"]) / pmax(bv$polarS[, "A"] + bv$polarN[, "A"], 1)
  if (is.function(progress)) progress(0.4, "box model, polar start")
  a <- ocean_boxes(bv, par, start = "thermal")
  if (is.function(progress)) progress(0.7, "box model, shelf start")
  b <- ocean_boxes(bv, par, start = "haline")
  G$state_thermal <- a$state; G$state_haline <- b$state; G$bistable <- a$state != b$state
  G$O2_thermal <- a$O2_deep; G$O2_haline <- b$O2_deep; G$q_thermal <- a$q_total; G$q_haline <- b$q_total
  G$share_polar_thermal <- a$share_polar; G$share_polar_haline <- b$share_polar
  if (ice_line) {
    if (is.function(progress)) progress(0.85, "land-ice line")
    G$ice_area <- vapply(seq_len(nrow(G)), function(i) {
      x <- .resolve(climate_experiment(co2_pal = G$co2[i], sun_age = inp$sun_age, sensitivity = inp$sensitivity, gradient = G$gradient[i]), inp$age)
      sum(land_ice_flags_forcing(inp$age, c(dT = x$dT, grad = G$gradient[i]), inp$map_file, cfg)$area) }, 0)
  }
  structure(list(grid = G, co2 = co2, gradient = gradient, inputs = inp[c("age", "label", "co2_pal", "gradient", "run")], par = par),
            class = "ocean_regime")
}
#' figure of the regime map: (a) the state from a polar start (dots: a second state is possible
#' from a shelf start) with the constraints; (b) deep O2
plot_ocean_regime <- function(rm, sst_band = NULL, which = c("thermal", "haline")) {
  which <- match.arg(which); G <- rm$grid
  st <- G[[paste0("state_", which)]]; lv <- names(OCEAN_STATE_COLS)
  Z <- matrix(match(st, lv), length(rm$co2)); X <- log2(rm$co2); Y <- rm$gradient
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1), mgp = c(2.4, 0.6, 0), cex = 0.82, las = 1); on.exit(par(op))
  dx <- diff(X)[1] / 2; dy <- diff(Y)[1] / 2
  image(X, Y, Z, breaks = seq(0.5, 4.5, 1), col = OCEAN_STATE_COLS, xaxt = "n", xlab = "CO2 (PAL)", ylab = "Equator-to-pole gradient factor",
        main = sprintf("Ocean state (%s start)", if (which == "thermal") "polar" else "shelf"), font.main = 1, cex.main = 0.95)
  axis(1, at = X, labels = signif(rm$co2, 2), cex.axis = 0.8)
  bi <- matrix(G$bistable, length(X))
  if (any(bi)) points(rep(X, length(Y))[bi], rep(Y, each = length(X))[bi], pch = 20, cex = 0.5, col = "white")
  Ts <- matrix(G$sst_tropics, length(X))
  contour(X, Y, Ts, add = TRUE, col = "#7a0177", labcex = 0.6, lwd = 0.7, nlevels = 8)
  if (length(sst_band) == 2 && all(is.finite(sst_band))) {
    contour(X, Y, Ts, levels = sst_band, add = TRUE, col = "#7a0177", lwd = 2.2, drawlabels = FALSE)
    .shade_between(X, Y, Ts, sst_band)
  }
  if (!is.null(G$ice_area)) contour(X, Y, matrix(G$ice_area, length(X)), levels = 5e4, add = TRUE, col = "#08306b", lwd = 2, lty = 2, drawlabels = FALSE)
  b <- rm$inputs; points(log2(b$co2_pal), b$gradient, pch = 4, cex = 1.3, lwd = 2)
  legend("bottomleft", bty = "o", bg = adjustcolor("white", 0.85), cex = 0.66,
         fill = c(OCEAN_STATE_COLS, NA, NA, NA, NA), border = c(rep("black", 4), NA, NA, NA, NA),
         legend = c(lv, "\u2022 second state from a shelf start", "purple: tropical SST (\u00b0C); thick = proxy band",
                    if (!is.null(G$ice_area)) "blue dashed: land ice starts (colder side)" else "", "\u00d7 the climate run"))
  O <- matrix(G[[paste0("O2_", which)]], length(X))
  image(X, Y, O, col = grDevices::hcl.colors(30, "RdYlBu"), xaxt = "n", xlab = "CO2 (PAL)", ylab = "", main = "Deep-sea O2 (\u00b5mol/kg)",
        font.main = 1, cex.main = 0.95, zlim = c(0, max(O, 10, na.rm = TRUE)))
  axis(1, at = X, labels = signif(rm$co2, 2), cex.axis = 0.8)
  contour(X, Y, O, add = TRUE, labcex = 0.6, col = adjustcolor("black", 0.5), levels = c(0, 25, 50, 100, 150, 200, 250, 300))
  if (length(sst_band) == 2 && all(is.finite(sst_band))) contour(X, Y, Ts, levels = sst_band, add = TRUE, col = "#7a0177", lwd = 2, drawlabels = FALSE)
  points(log2(b$co2_pal), b$gradient, pch = 4, cex = 1.3, lwd = 2)
  mtext(sprintf("%s Ma; %s", format(b$age), b$label), 1, 3.3, cex = 0.55, adj = 1, col = "grey30")
}
.shade_between <- function(X, Y, Z, band) {
  ok <- Z >= band[1] & Z <= band[2]
  if (!any(ok)) return(invisible())
  xx <- rep(X, length(Y))[ok]; yy <- rep(Y, each = length(X))[ok]
  rect(xx - diff(X)[1] / 2, yy - diff(Y)[1] / 2, xx + diff(X)[1] / 2, yy + diff(Y)[1] / 2, col = adjustcolor("#7a0177", 0.12), border = NA)
}
#' the gradient factor that gives a tropical SST (annual mean, 0-20 deg) at each CO2: the line
#' through the CO2 x gradient plane that a conodont d18O temperature allows
gradient_for_tropics <- function(inp, sst_target, co2 = 2^seq(1, 4, 0.25), lower = 0.3, upper = 1.3) {
  lat <- GRID1$lat; tr <- abs(lat) <= 20; A <- area_matrix(lat, 360); sea <- inp$lf < 0.5 & is.finite(inp$sst[, , 1])
  ann <- apply(pmax(inp$sst, -1.8), 1:2, mean); sel <- sea & matrix(tr, 180, 360)
  f <- function(c2, g) { D <- .ocean_dmon(inp, c2, g); Z <- ann + matrix(colMeans(D), 180, 360); sum((A * Z)[sel]) / sum(A[sel]) - sst_target }
  data.frame(co2 = co2, gradient = vapply(co2, function(c2) {
    lo <- f(c2, lower); hi <- f(c2, upper); if (sign(lo) == sign(hi)) return(NA_real_)
    stats::uniroot(function(g) f(c2, g), c(lower, upper), tol = 1e-3)$root }, 0))
}
#' hysteresis: CO2 lowered step by step and raised again (gradient factor fixed); the box model
#' continues from the state of the step before, so the two branches show where the ocean flips
ocean_hysteresis <- function(inp, co2 = 2^seq(log2(16), log2(2), length.out = 13), gradient = inp$gradient, par = list(), years_step = 4000) {
  par <- modifyList(OCEAN_PARAMS, par); cells <- .ocean_cells(inp, par)
  path <- c(co2, rev(co2)[-1]); dirn <- c(rep("CO2 falling", length(co2)), rep("CO2 rising", length(co2) - 1))
  bv <- .ocean_box_values(cells, par, lapply(path, function(c2) .ocean_dmon(inp, c2, gradient)))
  b1 <- lapply(bv, function(m) m[1, , drop = FALSE])
  k1 <- .box_setup(b1, par); st <- .box_run(.box_init(b1, par, "thermal", k1), k1, par, par$years)
  out <- vector("list", length(path))
  for (i in seq_along(path)) {
    bi <- lapply(bv, function(m) m[i, , drop = FALSE])
    ki <- .box_setup(bi, par); st <- .box_run(st, ki, par, years_step)
    out[[i]] <- cbind(co2 = path[i], direction = dirn[i], .box_diag(st, ki, par))
  }
  do.call(rbind, out)
}
plot_ocean_hysteresis <- function(h) {
  op <- par(mar = c(4.2, 4.4, 2.4, 1), mgp = c(2.4, 0.6, 0), cex = 0.85, las = 1); on.exit(par(op))
  dn <- h$direction == "CO2 falling"
  plot(h$co2, h$share_polar, log = "x", type = "n", ylim = c(0, 1), xlab = "CO2 (PAL)", ylab = "Share of deep water from the poles",
       main = "Hysteresis: CO2 down (blue) and up (red)", font.main = 1, cex.main = 0.9)
  lines(h$co2[dn], h$share_polar[dn], type = "b", col = "#2c7fb8", pch = 19); lines(h$co2[!dn], h$share_polar[!dn], type = "b", col = "#d7301f", pch = 17)
  abline(h = c(1 / 3, 2 / 3), lty = 3, col = "grey50"); text(max(h$co2), c(0.85, 0.15), c("P-like", "S-like"), adj = 1, cex = 0.75, col = "grey30")
}

#' schematic of the boxes of one ocean state: boxes with winter T, S and density; arrows ~ flow
plot_ocean_state <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); r <- os[[which]]; bv <- os$boxes
  op <- par(mar = c(1, 1, 2.4, 1), cex = 0.8); on.exit(par(op))
  plot.new(); plot.window(c(0, 10), c(0, 7))
  title(sprintf("%s \u2013 %s start", r$state, if (which == "thermal") "polar" else "shelf"), font.main = 1, cex.main = 1)
  bx <- list(polarS = c(0.2, 4.6, 2.4, 6.6), shelf = c(2.8, 4.6, 4.9, 6.6), open = c(5.3, 4.6, 7.4, 6.6), polarN = c(7.8, 4.6, 9.8, 6.6))
  lab <- c(polarS = "South polar sea", shelf = "Evaporative shelf seas", open = "Open sea", polarN = "North polar sea")
  for (b in names(bx)) {
    v <- bx[[b]]; S <- r[[paste0("S_", b)]] %||% NA; Tb <- if (b == "open") bv$open[1, "Tann"] else bv[[b]][1, "T"]
    rho <- if (b == "open") rho_sw(Tb, S) - 1000 else r[[paste0("rho_", b)]]
    rect(v[1], v[2], v[3], v[4], col = if (b %in% c("polarS", "polarN")) "#deebf7" else if (b == "shelf") "#fee6ce" else "#f0f0f0", border = "grey40")
    text(mean(v[c(1, 3)]), v[4] - 0.3, lab[[b]], font = 2, cex = 0.8)
    text(mean(v[c(1, 3)]), v[2] + 0.9, sprintf("%.1f \u00b0C%s\n%.2f psu\n\u03c3 %.2f", Tb, if (b == "open") " (ann.)" else " (winter)", S, rho), cex = 0.72)
    if (b %in% c("polarS", "polarN")) text(mean(v[c(1, 3)]), v[2] + 0.2, sprintf("sea ice %.0f %%", 100 * bv[[b]][1, "ice"]), cex = 0.65, col = "#08519c")
  }
  rect(0.2, 0.4, 9.8, 3.6, col = "#e7e1ef", border = "grey40")
  text(5, 2.6, sprintf("Deep sea: %.1f \u00b0C, %.2f psu, \u03c3 %.2f", r$T_deep, r$S_deep, r$rho_deep), font = 2, cex = 0.8)
  text(5, 1.9, sprintf("O2 %.0f \u00b5mol/kg%s; phosphate %.2f \u00b5mol/kg", r$O2_deep, if (isTRUE(r$anoxic)) " (anoxic)" else "", r$P_deep), cex = 0.75)
  for (b in c("polarS", "shelf", "polarN")) { q <- r[[paste0("q_", b)]]; v <- bx[[b]]
    if (q > 0.05) arrows(mean(v[c(1, 3)]), v[2], mean(v[c(1, 3)]), 3.6, lwd = 1 + 2 * sqrt(q), col = "#3f007d", length = 0.1)
    text(mean(v[c(1, 3)]) + 0.15, 4.1, sprintf("%.1f Sv", q), adj = 0, cex = 0.7, col = "#3f007d") }
  arrows(6.35, 3.6, 6.35, 4.6, lwd = 1 + 2 * sqrt(r$q_total), col = "#636363", length = 0.1)
  text(6.5, 4.1, "upwelling", adj = 0, cex = 0.65, col = "#636363")
  mtext(os$label, 1, -0.2, cex = 0.6, col = "grey30")
}

# ---- through time: the ocean state along an orbit series -----------------------------------------
#' deep-water source and deep O2 through time: the box model run transiently along an orbit series,
#' its surface boxes following the emulated climate of an orbital ensemble (SST, sea ice, P - E),
#' and optionally the sea level of an orbital ice run (shelves fall dry)
ocean_series <- function(ens, astro, par = list(), orbital_ice = NULL, step_kyr = 1, cfg = SILUCARB, progress = NULL) {
  par <- modifyList(OCEAN_PARAMS, par)
  if (is.character(ens)) ens <- ensemble_load(ens, cfg)
  emu <- ensemble_emulator(ens, cfg = cfg)
  st <- ensemble_status(ens, cfg); done <- which(st$status == "done")
  d <- ens$design[done, ]; ctr <- done[which.min(abs(d$obl - mean(range(d$obl))) + d$ecc * 10)]
  inp <- ocean_inputs(ens$ids[ctr], cfg)
  R0 <- emulate_fields(emu, ens$design$obl[ctr], ens$design$ecc[ctr] * sin(ens$design$peri[ctr] * pi / 180), ens$design$ecc[ctr] * cos(ens$design$peri[ctr] * pi / 180))
  cells <- .ocean_cells(inp, par)
  sl <- if (!is.null(orbital_ice)) orbital_ice_timeseries(orbital_ice, cfg) else NULL
  tt <- seq(min(astro$t_kyr), max(astro$t_kyr), by = step_kyr)
  ap <- function(v) stats::approx(astro$t_kyr, v, tt, rule = 2)$y
  orb <- data.frame(t_kyr = tt, ecc = ap(astro$ecc), obl = ap(astro$obl), esin = ap(astro$esin), ecos = ap(astro$ecos))
  if (!is.null(astro$age_ma)) orb$age_ma <- ap(astro$age_ma)
  orb$sea_level <- if (!is.null(sl)) stats::approx(sl$t_kyr, sl$sea_level_m, tt, rule = 2)$y else 0
  # box values per step: emulated change of SST / sea ice / P - E relative to the central run
  rows <- vector("list", length(tt))
  for (i in seq_along(tt)) {
    R <- emulate_fields(emu, orb$obl[i], orb$esin[i], orb$ecos[i])
    inp2 <- inp
    for (k in 1:12) {
      inp2$sst[, , k] <- inp$sst[, , k] + .to_1deg(.speedy_out_field(R[, , k, "sst"] - R0[, , k, "sst"]))
      inp2$ice[, , k] <- pmin(pmax(inp$ice[, , k] + .to_1deg(.speedy_out_field(R[, , k, "sea_ice"] - R0[, , k, "sea_ice"])), 0), 1)
    }
    dpe <- apply((R[, , , "precip"] - R[, , , "evap"]) - (R0[, , , "precip"] - R0[, , , "evap"]), 1:2, mean)
    inp2$fw <- inp$fw + .to_1deg(.speedy_out_field(dpe))
    c2 <- .ocean_cells(inp2, par)
    bv <- .ocean_box_values(c2, par, list(NULL), sea_level = orb$sea_level[i])
    rows[[i]] <- bv
    if (is.function(progress) && i %% 20 == 0) progress(0.6 * i / length(tt), sprintf("climate of %d of %d kyr", i, length(tt)))
  }
  bvall <- lapply(names(rows[[1]]), function(b) do.call(rbind, lapply(rows, `[[`, b))); names(bvall) <- names(rows[[1]])
  # transient box model: one configuration, forcing changes every step
  r <- .ocean_transient(bvall, par, step_kyr * 1000, progress)
  cbind(orb, r)
}
.ocean_transient <- function(bv, par, step_yr, progress = NULL) {
  n <- nrow(bv$open); out <- vector("list", n)
  # spin-up on the first step from the polar start, then each step continues from the last state
  b1 <- lapply(bv, function(m) m[1, , drop = FALSE])
  k1 <- .box_setup(b1, par); st <- .box_run(.box_init(b1, par, "thermal", k1), k1, par, par$years)
  for (i in seq_len(n)) {
    bi <- lapply(bv, function(m) m[i, , drop = FALSE])
    ki <- .box_setup(bi, par); st <- .box_run(st, ki, par, step_yr)
    out[[i]] <- .box_diag(st, ki, par)
    if (is.function(progress) && i %% 50 == 0) progress(0.6 + 0.4 * i / n, sprintf("box model %d of %d kyr", i, n))
  }
  do.call(rbind, out)
}
plot_ocean_series <- function(os, main = NULL) {
  deep <- !is.null(os$age_ma); x <- if (deep) os$age_ma else os$t_kyr / 1000; xl <- if (deep) rev(range(x)) else range(x)
  op <- par(mfrow = c(4, 1), mar = c(0.4, 4.6, 0.4, 4.2), oma = c(3.4, 0, 2, 0), mgp = c(2.6, 0.6, 0), cex = 0.8, las = 1); on.exit(par(op))
  plot(x, os$ecc, type = "l", xlim = xl, xaxt = "n", ylab = "Eccentricity", col = "#555555")
  par(new = TRUE); plot(x, os$obl, type = "l", xlim = xl, axes = FALSE, xlab = "", ylab = "", col = "#31a354", lwd = 0.6); axis(4, col.axis = "#31a354")
  plot(x, os$q_polarS + os$q_polarN, type = "l", xlim = xl, xaxt = "n", ylab = "Deep water (Sv)", col = "#2c7fb8", ylim = range(0, os$q_total))
  lines(x, os$q_shelf, col = "#e6550d"); lines(x, os$q_total, col = "black", lwd = 0.6, lty = 3)
  legend("topleft", bty = "n", horiz = TRUE, cex = 0.7, lty = c(1, 1, 3), col = c("#2c7fb8", "#e6550d", "black"), legend = c("from the poles", "from the shelf seas", "total"))
  plot(x, os$O2_deep, type = "l", xlim = xl, xaxt = "n", ylab = "Deep O2 (\u00b5mol/kg)", col = "#7a0177", ylim = range(0, os$O2_deep))
  cols <- OCEAN_STATE_COLS[os$state]
  plot(x, os$share_polar, type = "n", xlim = xl, ylim = c(0, 1), ylab = "Share from the poles")
  segments(x[-length(x)], 0, x[-1], 0.08, col = cols[-length(cols)], lwd = 6, lend = 1)
  lines(x, os$share_polar, col = "#2c7fb8"); abline(h = c(1 / 3, 2 / 3), lty = 3, col = "grey60")
  if (any(os$sea_level != 0)) { par(new = TRUE); plot(x, os$sea_level, type = "l", xlim = xl, axes = FALSE, xlab = "", ylab = "", col = "#08306b", lty = 2); axis(4, col.axis = "#08306b"); mtext("sea level (m)", 4, 2.6, cex = 0.7, col = "#08306b", las = 0) }
  mtext(if (deep) "Age (Ma)" else "Time (Myr, relative)", 1, 2.2, cex = 0.8)
  mtext(main %||% "Ocean state through time: P-like (blue bar) / S-like (orange) / mixed / stagnant", 3, 0.5, outer = TRUE, cex = 0.8)
}


# ---- where the water sinks and rises: map and pole-to-pole section (from the model) ---------------
#' winter surface density of every sea cell of a box-model state and its share of the sinking:
#' winter SST of the run (each hemisphere's winter), salinity of the run shifted by the box model's
#' salinity change of the cell's box, + the sea-ice effect in the polar boxes. A source box sinks
#' q (Sv) in the model; that flow is shared over its cells by area x density excess over the deep
#' water (only cells denser than the deep water sink). Returns lat, lon, row, box, area_km2,
#' sigma (kg/m3 - 1000), excess, sv, ice (winter sea ice), w (wind-driven vertical velocity, m/day).
ocean_cell_state <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); r <- os[[which]]; par <- os$par; inp <- os$inputs
  cells <- .ocean_cells(inp, par); tf <- -1.8; out <- list()
  A1 <- area_matrix(GRID1$lat, 360) / sum(area_matrix(GRID1$lat, 360)) * 5.1007e8
  bnames <- names(cells); if (any(grepl("^polarS_", bnames))) bnames <- setdiff(bnames, "polarS")
  if (any(grepl("^shelf_", bnames))) bnames <- setdiff(bnames, "shelf")
  for (b in bnames) {
    C <- cells[[b]]; if (!length(C$k)) next
    Tc <- pmax(C$sst, tf); s <- C$lat < 0
    Tw <- ifelse(s, rowMeans(Tc[, 6:8, drop = FALSE]), rowMeans(Tc[, c(12, 1, 2), drop = FALSE]))
    fi <- ifelse(s, rowMeans(C$ice[, 6:8, drop = FALSE]), rowMeans(C$ice[, c(12, 1, 2), drop = FALSE]))
    sal <- ifelse(is.finite(C$sal), C$sal, par$S0); W <- C$A / sum(C$A)
    Sbox <- r[[paste0("S_", b)]]; dS <- if (is.null(Sbox)) 0 else Sbox - sum(W * sal)
    Seff <- sal + dS + if (startsWith(b, "polar")) par$seaice_effect * par$seaice_psu * fi else 0
    tgt <- if (!is.null(r$rho_int) && identical(r[[paste0("layer_", b)]], "intermediate")) r$rho_int else r$rho_deep
    sig <- rho_sw(Tw, Seff) - 1000; ex <- sig - tgt
    q <- if (b == "open") 0 else r[[paste0("q_", b)]]
    sv <- 0 * ex
    if (isTRUE(q > 0)) { w <- C$A * pmax(ex, 0); if (sum(w) <= 0) w <- C$A * exp(5 * (ex - max(ex))); sv <- q * w / sum(w) }
    out[[b]] <- data.frame(lat = C$lat, lon = GRID1$lon[(C$k - 1) %/% 180 + 1], row = C$row, box = b, area_km2 = A1[C$k],
                           sigma = sig, excess = ex, sv = sv, ice = fi, w = if (!is.null(inp$w)) inp$w[C$k] else NA)
  }
  do.call(rbind, out)
}

#' map of the circulation of an ocean state (all from the model): background = vertical velocity at
#' 10 m of the run (HadCM3L + SPEEDY Ekman change; red up, blue down); shading = where the box model
#' sinks its deep water (Sv per million km2; purple polar seas, orange evaporative shelf seas),
#' triangles = total per source at its centre; winter sea ice; evaporative shelf seas; palaeo-coordinates
plot_ocean_circulation_map <- function(os, which = c("thermal", "haline"), wlim = NULL, notes = TRUE, main = NULL) {
  which <- match.arg(which); inp <- os$inputs; r <- os[[which]]; par <- os$par
  lon <- GRID1$lon; lat <- GRID1$lat; land <- inp$lf >= 0.5
  op <- par(mar = if (notes) c(5.4, 3.6, 2.6, 1) else c(2, 3, 1.6, 0.5), mgp = c(2, 0.5, 0), cex = if (notes) 0.85 else 0.7, las = 1); on.exit(par(op))
  plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i", xlab = "", ylab = "Palaeo-latitude", asp = 1, axes = FALSE)
  axis(1, at = seq(-180, 180, 60), labels = paste0(seq(-180, 180, 60), "\u00b0")); axis(2, at = seq(-90, 90, 30), labels = paste0(seq(-90, 90, 30), "\u00b0"))
  lim <- NA
  if (!is.null(inp$w)) {
    w <- inp$w; w[land] <- NA
    lim <- wlim %||% max(0.1, signif(stats::quantile(abs(w), 0.97, na.rm = TRUE), 1))
    cols <- grDevices::colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7", "#ef8a62", "#b2182b"))(21)
    image(lon, lat, t(pmin(pmax(w, -lim), lim)), breaks = seq(-lim, lim, length.out = 22), col = cols, add = TRUE, useRaster = FALSE)
  } else rect(-180, -90, 180, 90, col = "#f0f4f8", border = NA)
  cs <- ocean_cell_state(os, which)
  M <- matrix(0, 180, 360); M[cbind(cs$row, match(cs$lon, lon))[startsWith(cs$box, "shelf"), , drop = FALSE]] <- 1
  contour(lon, lat, t(M), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#e6550d", lwd = 0.8)
  I <- matrix(0, 180, 360); I[cbind(cs$row, match(cs$lon, lon))] <- cs$ice
  .hatch_cells(lon, lat, I >= 0.5, col = "white")
  contour(lon, lat, t(I), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#4a6fa5", lwd = 0.9, lty = 2)
  L <- matrix(NA, 180, 360); L[land] <- 1
  image(lon, lat, t(L), col = "#cfc6b4", add = TRUE, useRaster = FALSE)
  contour(lon, lat, t(inp$lf), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#5b5246", lwd = 0.7)
  abline(h = c(-1, 1) * par$polar_lat, lty = 3, col = "grey40", lwd = 0.6); abline(h = 0, col = "grey60", lwd = 0.4)
  d <- cs[cs$sv > 0, ]
  for (bx in unique(d$box)) {
    e <- d[d$box == bx, ]; F <- matrix(NA, 180, 360); F[cbind(e$row, match(e$lon, lon))] <- e$sv / e$area_km2 * 1e6
    col0 <- .src_col(bx); mx <- max(F, na.rm = TRUE)
    image(lon, lat, t(F), breaks = seq(0, mx * 1.0001, length.out = 7), col = adjustcolor(grDevices::colorRampPalette(c("white", col0))(8)[3:8], 0.85), add = TRUE, useRaster = FALSE)
    q <- sum(e$sv); cy <- stats::weighted.mean(e$lat, e$sv)
    tl <- tapply(e$sv, round(e$lon / 10) * 10, sum); cx <- as.numeric(names(tl)[which.max(tl)])
    points(cx, cy, pch = 25, cex = 1.2 + 0.45 * sqrt(q), bg = col0, col = "white", lwd = 1)
    text(cx, cy + if (cy > 0) -8 else 8, sprintf("%.1f Sv sinks", q), cex = 0.75, font = 2, col = col0)
  }
  for (nb in grep("^shelf_[A-Za-z]+$", unique(cs$box), value = TRUE)) if (!grepl("^shelf_R[0-9]", nb)) {   # named basins
    Mb <- matrix(0, 180, 360); e <- cs[cs$box == nb, ]; Mb[cbind(e$row, match(e$lon, lon))] <- 1
    contour(lon, lat, t(Mb), levels = 0.5, add = TRUE, drawlabels = FALSE, col = .src_col(nb), lwd = 1.6)
    xl <- min(e$lon) - 2; yl <- mean(range(e$lat)); w <- strwidth(sprintf("%s basin", sub("shelf_", "", nb)), cex = 0.75)
    rect(xl - w - 1, yl - 3, xl + 1, yl + 3, col = adjustcolor("white", 0.8), border = NA)
    text(xl, yl, sprintf("%s basin", sub("shelf_", "", nb)), col = .src_col(nb), font = 2, cex = 0.75, adj = c(1, 0.5)) }
  box()
  title(main %||% sprintf("Where water sinks and rises \u2013 %s (%s start)", r$state, if (which == "thermal") "polar" else "shelf"), font.main = 1, cex.main = 0.95)
  if (!notes) return(invisible())
  mtext(sprintf("Shading: where the box model forms its deep water (cells denser than the deep sea, \u03c3 > %.2f; share ~ area \u00d7 density excess): purple = north polar sea, green = south polar sea, orange = evaporative shelf seas; triangles = total per source.", r$rho_deep), 1, 2.0, cex = 0.6, adj = 0)
  mtext(sprintf("Background: vertical velocity at 10 m of the run%s (HadCM3L + SPEEDY Ekman change): red = upwelling, blue = downwelling. Return of the %.1f Sv deep water: upwelling over the open sea.",
                if (is.finite(lim)) sprintf(" (\u00b1%g m/day)", lim) else "", r$q_total), 1, 2.8, cex = 0.6, adj = 0)
  mtext(sprintf("Orange line: evaporative (anti-estuarine) shelf seas; white hatching, dashed: winter sea ice \u2265 50 %%; dotted: polar boxes (%g\u00b0). %s Ma, %s.", par$polar_lat, format(inp$age), inp$label), 1, 3.6, cex = 0.6, adj = 0)
}
.hatch_cells <- function(lon, lat, M, col = "white", step = 2) {
  k <- which(M, arr.ind = TRUE); if (!nrow(k)) return(invisible())
  k <- k[(k[, 1] + k[, 2]) %% step == 0, , drop = FALSE]
  segments(lon[k[, 2]] - 0.5, lat[k[, 1]] - 0.5, lon[k[, 2]] + 0.5, lat[k[, 1]] + 0.5, col = col, lwd = 0.6)
}

#' pole-to-pole section of an ocean state, built from the model:
#'  top: winter surface density by latitude (area mean and densest 5 % of the sea cells) against the
#'       deep-water density - where the curve is above the line the water can sink (the model's rule);
#'  bottom: zonal-mean sea depth and land of the map, zonal-mean winter sea ice; sinking at the
#'       sinking-weighted latitude of each source (width ~ Sv); upwelling of the same total shared
#'       over the latitudes where the run's open sea upwells (zonal sum of w > 0); deep sea coloured
#'       by the model's deep O2. Sizes are the model's; only the depth axis is schematic (not to scale).
plot_ocean_section <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); r <- os[[which]]; par <- os$par; inp <- os$inputs
  cs <- ocean_cell_state(os, which); lat <- GRID1$lat; pl <- par$polar_lat
  # zonal profiles (1-degree rows)
  zrow <- function(v, f = function(x) mean(x)) vapply(1:180, function(i) { x <- v[cs$row == i]; if (length(x)) f(x) else NA }, 0)
  sig_m <- vapply(1:180, function(i) { k <- cs$row == i; if (any(k)) sum(cs$sigma[k] * cs$area_km2[k]) / sum(cs$area_km2[k]) else NA }, 0)
  sig_x <- zrow(cs$sigma, function(x) stats::quantile(x, 0.95, names = FALSE))
  ice <- zrow(cs$ice)
  landf <- rowMeans(inp$lf >= 0.5)
  dep <- vapply(1:180, function(i) { z <- inp$depth[i, ]; z <- z[is.finite(z) & inp$lf[i, ] < 0.5]; if (length(z)) mean(z) else NA }, 0)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2, 2), heights = c(1, 2.4))
  # --- top: density criterion
  par(mar = c(0.5, 4.2, 3, 4.2), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
  yl <- range(c(sig_m, sig_x, r$rho_deep), na.rm = TRUE); yl <- yl + c(-0.2, 0.45 * diff(yl))
  plot(lat, sig_m, type = "n", xlim = c(-90, 90), ylim = yl, xaxs = "i", xaxt = "n", ylab = "Winter \u03c3 (kg/m\u00b3)")
  # what the box model uses: one winter density per box (thick); the run's spatial detail (thin)
  shl <- cs[startsWith(cs$box, "shelf"), ]; shc <- if (nrow(shl)) stats::weighted.mean(shl$lat, pmax(shl$sv, 1e-9) * shl$area_km2) else 20
  bxs <- list(polarS = c(-90, -pl), polarN = c(pl, 90), shelf = shc + c(-8, 8), open = c(-pl, pl))
  bsig <- unlist(lapply(list(polarS = r$rho_polarS, polarN = r$rho_polarN, shelf = r$rho_shelf,
                             open = rho_sw(os$boxes$open[1, "Tann"], r$S_open) - 1000), function(v) unname(as.numeric(v)[1])))
  bq <- unlist(lapply(list(polarS = r$q_polarS, polarN = r$q_polarN, shelf = r$q_shelf, open = 0), function(v) unname(as.numeric(v)[1])))
  bcol <- c(polarS = "#3f007d", polarN = "#3f007d", shelf = "#e6550d", open = "#636363")
  for (b in names(bxs)) if (bq[[b]] > 0.05) rect(bxs[[b]][1], yl[1], bxs[[b]][2], yl[2], col = adjustcolor(bcol[[b]], 0.12), border = NA)
  lines(lat, sig_m, col = "#6baed6", lwd = 0.9); lines(lat, sig_x, col = "#fdae6b", lwd = 0.8, lty = 2)
  abline(h = r$rho_deep, col = "black", lwd = 1.4, lty = 5); abline(v = c(-pl, pl), lty = 3, col = "grey50")
  segments(bxs$open[1], bsig[["open"]], bxs$open[2], bsig[["open"]], col = bcol[["open"]], lwd = 3)
  for (b in c("polarS", "polarN", "shelf")) if (is.finite(bsig[[b]]) && any(startsWith(cs$box, b))) {
    segments(bxs[[b]][1], bsig[[b]], bxs[[b]][2], bsig[[b]], col = bcol[[b]], lwd = 4)
    text(mean(bxs[[b]]), bsig[[b]], sprintf("%.2f", bsig[[b]]), pos = 3, cex = 0.65, col = bcol[[b]], offset = 0.3) }
  legend("topleft", bty = "n", ncol = 3, cex = 0.62, lty = c(1, 1, 1, 5, 1, 2), lwd = c(4, 4, 3, 1.4, 0.9, 0.8),
         col = c("#3f007d", "#e6550d", "#636363", "black", "#6baed6", "#fdae6b"),
         legend = c("polar box", "shelf box", "open box", sprintf("deep water %.2f", r$rho_deep), "run: zonal mean", "run: densest 5 %"))
  title(sprintf("%s \u2013 pole-to-pole section (%s start)", r$state, if (which == "thermal") "polar" else "shelf"), font.main = 2, cex.main = 1.05)
  # --- bottom: section
  par(mar = c(4.6, 4.2, 0.5, 4.2))
  .ocean_section_panel(os, which, cs)
  mtext("Top: thick = winter surface density of each box as used by the box model (shaded = denser than the deep water, so it sinks); thin = the run's zonal pattern inside the boxes.", 1, 2.3, cex = 0.58, col = "grey30")
  mtext(sprintf("Bottom: sea floor and land = zonal means of the map; white = winter sea ice (zonal mean); arrows = box-model flows at their model latitudes (width ~ Sv); deep colour = model O2. %s Ma; %s.",
                format(inp$age), inp$label), 1, 3.0, cex = 0.58, col = "grey30")
}
#' the section panel (sea floor, land, sea ice, sinking, upwelling, deep O2) of plot_ocean_section()
.ocean_section_panel <- function(os, which, cs = ocean_cell_state(os, which), key = TRUE, notes = TRUE, cex_lab = 1) {
  r <- os[[which]]; par <- os$par; inp <- os$inputs; lat <- GRID1$lat; pl <- par$polar_lat
  ice <- vapply(1:180, function(i) { x <- cs$ice[cs$row == i]; if (length(x)) mean(x) else NA }, 0)
  landf <- rowMeans(inp$lf >= 0.5)
  dep <- vapply(1:180, function(i) { z <- inp$depth[i, ]; z <- z[is.finite(z) & inp$lf[i, ] < 0.5]; if (length(z)) mean(z) else NA }, 0)
  zmax <- 4500; ys <- function(z) -pmin(z, zmax)                                     # depth (m, positive) -> y
  plot(NA, xlim = c(-90, 90), ylim = c(-zmax - 400, 1300), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(1, at = seq(-90, 90, 30), labels = c("South Pole", "60\u00b0S", "30\u00b0S", "Equator", "30\u00b0N", "60\u00b0N", "North Pole"), cex.axis = 0.8)
  axis(2, at = c(0, -1000, -2000, -3000, -4000), labels = c(0, 1, 2, 3, 4), cex.axis = 0.8); mtext("Depth (km, zonal mean)", 2, 2.4, cex = 0.8, las = 0)
  oc <- grDevices::colorRampPalette(c("#7f0000", "#d7301f", "#fdae61", "#ffffbf", "#abd9e9", "#4575b4", "#313695"))(100)
  o2c <- oc[max(1, min(100, round(r$O2_deep / 3.2) + 1))]
  ok <- is.finite(dep)
  # deep sea (O2 colour) and surface layer
  dd <- pmin(dep, zmax); dd[!ok] <- 0
  polygon(c(lat, rev(lat)), c(pmax(ys(dd), -300), rev(ys(dd))), col = adjustcolor(o2c, 0.6), border = NA)
  two <- !is.null(r$O2_int)
  if (two) {                                                                            # intermediate layer (300-1500 m)
    o2i <- oc[max(1, min(100, round(r$O2_int / 3.2) + 1))]
    polygon(c(lat, rev(lat)), c(pmax(ys(dd), -300), rev(pmax(ys(dd), -1500))), col = adjustcolor(o2i, 0.6), border = NA)
    lines(lat[dd > 1500], rep(-1500, sum(dd > 1500)), lty = 3, col = "grey40")
    text(-88, -1450, sprintf("intermediate: O2 %.0f", r$O2_int), adj = c(0, 0), cex = 0.6 * cex_lab, col = "grey20")
    text(-88, -1600, sprintf("abyss: O2 %.0f", r$O2_deep), adj = c(0, 1), cex = 0.6 * cex_lab, col = "grey20")
  }
  polygon(c(lat, rev(lat)), c(rep(0, 180), rev(pmax(ys(dd), -300))), col = "#eaf2fa", border = NA)
  lines(lat, ys(dd), col = "#5b5246", lwd = 1.2)                                     # sea floor (zonal mean)
  polygon(c(-90, lat, 90), c(-zmax - 400, ys(dd), -zmax - 400), col = "#d9d0bf", border = NA)
  # land of the map: height ~ land fraction of the latitude
  polygon(c(lat, rev(lat)), c(rep(0, 180), rev(landf * 900)), col = "#cfc6b4", border = "#5b5246", lwd = 0.6)
  # winter sea ice
  si <- ifelse(is.finite(ice), ice, 0)
  if (any(si > 0.05)) polygon(c(lat, rev(lat)), c(rep(0, 180), rev(si * 160)), col = "white", border = "#4a6fa5", lwd = 0.6)
  if (par$seaice_effect < 0 && any(si > 0.05)) {                                   # melt-water lid: fresh layer under the ice
    polygon(c(lat, rev(lat)), c(rep(0, 180), rev(-si * 220)), col = adjustcolor("#6baed6", 0.7), border = NA)
    k <- which.max(si); text(lat[k], -260, "melt-water lid", cex = 0.62, col = "#08519c", font = 3, adj = c(0.5, 1)) }
  abline(v = c(-pl, pl), lty = 3, col = "grey50")
  # sinking
  lw <- function(q) 1 + 7 * sqrt(q / max(r$q_total, 1))
  sh <- unique(grep("^shelf", cs$box, value = TRUE))
  src <- list(); for (b in c("polarS", "polarN", sh)) { e <- cs[(if (b == "polarS") startsWith(cs$box, "polarS") else cs$box == b) & cs$sv > 0, ]
    if (nrow(e) && sum(e$sv) > 0.05) src[[b]] <- c(lat = stats::weighted.mean(e$lat, e$sv), q = sum(e$sv)) }
  for (b in names(src)) { x <- src[[b]][["lat"]]; q <- src[[b]][["q"]]; col0 <- .src_col(b)
    zb <- min(zmax - 500, max(dd[max(1, round(x + 90.5) + c(-3:3))], na.rm = TRUE) - 300)
    lay <- if (two) r[[paste0("layer_", b)]] else "abyss"
    yend <- if (identical(lay, "intermediate")) -1400 else -max(zb, 1500)
    arrows(x, -40, x, yend, col = col0, lwd = lw(q), length = 0.12)
    text(x + 1.5, -900, sprintf("%.1f Sv%s", q, if (two) paste0("\n", if (identical(lay, "intermediate")) "intermediate" else "abyssal") else ""),
         adj = 0, col = col0, font = 2, cex = 0.85 * cex_lab) }
  for (b in c("polarS", "polarN")) if (is.null(src[[b]])) {
    x <- if (b == "polarS") -(pl + 90) / 2 else (pl + 90) / 2
    text(x, -1300, if (notes) "no deep-water\nformation\n(lighter than\nthe deep water)" else "no sinking", cex = 0.68 * cex_lab, col = "#525252", font = 3) }
  # upwelling: shared by the run's open-sea upwelling, 10-degree bands
  if (r$q_total > 0.05 && any(is.finite(cs$w))) {
    o <- cs[cs$box == "open" & is.finite(cs$w) & cs$w > 0, ]
    up <- tapply(o$w * o$area_km2, floor(o$lat / 10) * 10 + 5, sum); up <- up / sum(up)
    up <- sort(up, decreasing = TRUE); up <- up[up >= 0.08]; up <- up / sum(up)
    for (k in seq_along(up)) { x <- as.numeric(names(up)[k]); q <- r$q_total * up[[k]]
      arrows(x, -zmax + 700, x, -320, col = "#636363", lwd = lw(q), length = 0.1)
      text(x + 1.5, -zmax + 1000, sprintf("%.1f", q), adj = 0, cex = 0.7, col = "#525252") }
    xs_src <- vapply(src, function(v) v[["lat"]], 0)
    for (x in xs_src) { u <- as.numeric(names(up))[which.min(abs(as.numeric(names(up)) - x))]
      if (abs(u - x) > 4) arrows(x + sign(u - x) * 3, -zmax + 500, u - sign(u - x) * 3, -zmax + 500, col = "#636363", lwd = 1.6, length = 0.08)
      if (abs(u - x) > 4) arrows(u + sign(x - u) * 3, -150, x - sign(x - u) * 3, -150, col = "#969696", lwd = 1.2, length = 0.08) }
    if (notes) text(-88, -zmax + 250, sprintf("upwelling of the %.1f Sv, shared by the run's open-sea upwelling (Sv per 10\u00b0 band)", r$q_total), adj = 0, cex = 0.62, col = "#525252")
  } else text(0, -2200, "almost no overturning: stagnant deep sea", cex = 0.9, col = "#525252")
  text(0, -zmax - 200, sprintf("Deep sea: %.1f \u00b0C, %.2f psu, O2 %.0f \u00b5mol/kg%s, phosphate %.2f \u00b5mol/kg, renewal ~%.0f yr", r$T_deep, r$S_deep,
                               r$O2_deep, if (r$O2_deep < 5) " (anoxic)" else "", r$P_deep, r$vent_age_yr %||% NA), cex = 0.78 * cex_lab, font = 2)
  if (!key) return(invisible())
  # O2 colour key
  par(xpd = NA); kx <- seq(91.5, 95, length.out = 2); ky <- seq(-4000, -1000, length.out = 101)
  for (k in 1:100) rect(kx[1], ky[k], kx[2], ky[k + 1], col = oc[k], border = NA)
  rect(kx[1], -4000, kx[2], -1000, border = "grey40"); text(kx[2] + 0.6, c(-4000, -2500, -1000), c(0, 160, 320), adj = 0, cex = 0.6)
  text(kx[1], -800, "deep O2", adj = 0, cex = 0.6); points(mean(kx), -4000 + 3000 * min(r$O2_deep, 320) / 320, pch = 17, cex = 0.9)
  par(xpd = FALSE)
}

# ---- Pohl et al. (2021) scenarios, simplified ------------------------------------------------------
#  Pohl, A., Lu, Z., ... Ridgwell, A. (2021) Vertical decoupling in Late Ordovician anoxia due to
#  reorganization of ocean circulation. Nature Geoscience 14, 868-873. cGENIE (3D ocean, EMBM
#  atmosphere, sea ice), 445 Ma (Scotese & Wright 2018), solar constant 3.65 % lower, circular
#  orbit, atmospheric O2 x0.4 (8.38 %), 10 kyr to equilibrium; ensemble CO2 x5-x24 PI x PO4 x0.2-x1
#  (x1 = 2.159 umol/kg). Fig. 4: x5, x8.5, x24 CO2 at x0.4 PO4 - overturning, ventilation age,
#  deep O2, sea ice: as CO2 falls, seasonal sea ice over the North Pole and its melt water shut the
#  northern deep-water formation. Scenario 1: x24 -> x8.5 CO2 at x0.4 PO4 (cooling alone expands
#  seafloor anoxia); scenario 2: x7 -> x5 CO2 with PO4 x0.2 -> x0.6.
#  Simplified here: the climate of one run (HadCM3L + SPEEDY) is shifted to each CO2 as in the
#  SPEEDY forcing (sea surface, sea ice where it freezes, freshwater 7 %/K), and the box model gives
#  the deep-water sources, overturning, deep O2 and renewal time. Not reproduced: the 3D circulation
#  and the percentage of anoxic seafloor (one deep box), the wind (kept from the run).
POHL2021 <- list(co2 = c(5, 8.5, 10, 24), po4_rel = 0.4, po4_modern = 2.159, o2_pal = 0.4,
                 grid_co2 = c(5, 6, 7, 8.5, 10, 24), grid_po4 = c(0.2, 0.4, 0.6, 0.8, 1.0),
                 # set-up closer to their cGENIE: two deep layers (both poles can form deep water at once),
                 # flat ocean floor (no shelf seas), melt-water lid; the south polar sea as basins (split at land
                 # barriers), each with its own runoff, temperature and sea ice
                 par = list(deep_layers = 2, shelf_depth = 0, seaice_effect = -1, seaice_psu = 2, seaice_conv = 1, fw_polarS = 1,
                            south_basins = 4, years = 10000))

#' the inputs of a run shifted to another CO2 (PAL) and gradient factor, as in the SPEEDY forcing
ocean_inputs_shift <- function(inp, co2_pal, gradient = inp$gradient, cc = OCEAN_PARAMS$cc) {
  D <- .ocean_dmon(inp, co2_pal, gradient); tf <- -1.8; out <- inp
  for (k in 1:12) {
    Dm <- matrix(D[k, ], 180, 360); raw <- inp$sst[, , k]
    out$sst[, , k] <- raw + Dm
    out$ice[, , k] <- pmax(inp$ice[, , k] * pmax(0, 1 - 0.25 * pmax(Dm, 0)), pmin(1, pmax(0, (tf - (pmax(raw, tf) + Dm)) / 2)))
  }
  out$fw <- inp$fw * (1 + cc * matrix(colMeans(D), 180, 360))
  x <- .resolve(climate_experiment(co2_pal = co2_pal, sun_age = inp$sun_age, sensitivity = inp$sensitivity, gradient = gradient), inp$age)
  out$dT <- x$dT; out$gradient <- gradient; out$co2_pal <- x$co2_pal_eff
  out$label <- sprintf("%s shifted to CO2 %g PAL%s (sea surface, sea ice, freshwater; wind of the run)", inp$label, co2_pal,
                       if (abs(gradient - inp$gradient) > 1e-6) sprintf(", gradient factor %.2f", gradient) else "")
  out
}
#' the three climate states of Pohl et al. (2021, Fig. 4): ocean states at CO2 x5, x8.5, x24 PI with
#' atmospheric O2 x0.4 and phosphate x0.4 (melt-water lid on: the mechanism they found)
pohl2021_scenarios <- function(inp, co2 = POHL2021$co2, gradient = inp$gradient, par = list()) {
  par <- modifyList(modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), list(o2_pal = POHL2021$o2_pal, p_mean = POHL2021$po4_rel * POHL2021$po4_modern)), par)
  co2 <- sort(co2, decreasing = TRUE)
  sts <- lapply(co2, function(c2) { i2 <- ocean_inputs_shift(inp, c2, gradient); os <- ocean_state(NULL, par, inputs = i2); os$co2 <- c2; os })
  tab <- do.call(rbind, lapply(sts, function(os) { r <- os$thermal
    data.frame(co2_PI = os$co2, state = r$state, south_polar_Sv = r$q_polarS, south_layer = r$layer_polarS %||% NA,
               north_polar_Sv = r$q_polarN, north_layer = r$layer_polarN %||% NA, shelf_Sv = r$q_shelf,
               total_Sv = r$q_total, deep_T = r$T_deep, deep_O2 = r$O2_deep, intermediate_O2 = r$O2_int %||% NA, renewal_yr = r$vent_age_yr,
               sea_ice_N = os$boxes$polarN[1, "ice"], sea_ice_S = os$boxes$polarS[1, "ice"], bistable = os$bistable) }))
  list(states = sts, table = tab, par = par)
}
#' figure as Pohl et al. (2021) Fig. 4: one row per CO2; pole-to-pole section (overturning), map of
#' deep-water formation and up- / downwelling, annual mean sea ice
plot_pohl2021 <- function(sc) {
  n <- length(sc$states)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(seq_len(3 * n), n, 3, byrow = TRUE), widths = c(1.25, 1.25, 0.9))
  for (os in sc$states) {
    r <- os$thermal; inp <- os$inputs
    par(mar = c(2.2, 3.2, 2, 2.2), mgp = c(1.8, 0.5, 0), cex = 0.62, las = 1)
    .ocean_section_panel(os, "thermal", key = FALSE, notes = FALSE, cex_lab = 0.9)
    title(sprintf("\u00d7%g CO2: south %.1f Sv, north %.1f Sv%s; abyssal O2 %.0f \u00b5mol/kg, renewal ~%.0f yr", os$co2, r$q_polarS, r$q_polarN,
                  if (r$q_shelf > 0.05) sprintf(", shelf %.1f Sv", r$q_shelf) else "", r$O2_deep, r$vent_age_yr), font.main = 2, cex.main = 1)
    plot_ocean_circulation_map(os, notes = FALSE, main = sprintf("\u00d7%g CO2: deep-water formation (shaded), up- / downwelling", os$co2))
    # annual mean sea ice
    par(mar = c(2, 2.6, 1.6, 0.5), cex = 0.62)
    I <- apply(inp$ice, 1:2, mean); I[inp$lf >= 0.5] <- NA
    image(GRID1$lon, GRID1$lat, t(I), zlim = c(0, 1), col = grDevices::hcl.colors(20, "Blues 3", rev = TRUE), xlab = "", ylab = "", asp = 1, axes = FALSE)
    L <- matrix(NA, 180, 360); L[inp$lf >= 0.5] <- 1; image(GRID1$lon, GRID1$lat, t(L), col = "#cfc6b4", add = TRUE)
    contour(GRID1$lon, GRID1$lat, t(inp$lf), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#5b5246", lwd = 0.5)
    axis(1, at = seq(-180, 180, 90)); axis(2, at = seq(-90, 90, 45)); box()
    title(sprintf("Sea ice (annual mean): N %.0f %%, S %.0f %% of the polar seas in winter", 100 * os$boxes$polarN[1, "ice"], 100 * os$boxes$polarS[1, "ice"]), font.main = 1, cex.main = 0.95)
  }
}
#' Pohl et al. (2021) Fig. 2 analogue: deep O2 of the box model over CO2 (x PI) x phosphate (x modern),
#' with their scenario 1 (x24 -> x8.5 CO2 at x0.4 PO4) and 2 (x7 -> x5 CO2 with x0.2 -> x0.6 PO4)
pohl2021_grid <- function(inp, co2 = POHL2021$grid_co2, po4 = POHL2021$grid_po4, gradient = inp$gradient, par = list()) {
  par <- modifyList(modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), list(o2_pal = POHL2021$o2_pal)), par)
  G <- expand.grid(co2 = co2, po4 = po4); cells <- .ocean_cells(inp, par)
  bv <- .ocean_box_values(cells, par, lapply(G$co2, function(c2) .ocean_dmon(inp, c2, gradient)))
  p2 <- par; p2$p_mean <- G$po4 * POHL2021$po4_modern
  r <- ocean_boxes(bv, p2)
  cbind(G, O2_deep = r$O2_deep, state = r$state, q_total = r$q_total, anoxic = r$O2_deep < 5)
}
plot_pohl2021_grid <- function(g) {
  co2 <- sort(unique(g$co2)); po4 <- sort(unique(g$po4))
  Z <- matrix(g$O2_deep[order(g$po4, g$co2)], length(co2))
  op <- par(mar = c(4.2, 4.6, 2.6, 5), mgp = c(2.4, 0.6, 0), cex = 0.85, las = 1); on.exit(par(op))
  cols <- grDevices::colorRampPalette(c("#7f0000", "#d7301f", "#fdae61", "#ffffbf", "#abd9e9", "#4575b4"))(40)
  zl <- c(0, max(Z, 50))
  image(seq_along(co2), seq_along(po4), Z, col = cols, zlim = zl, axes = FALSE, xlab = "CO2 (\u00d7 pre-industrial = PAL)", ylab = "Phosphate (\u00d7 modern 2.159 \u00b5mol/kg)",
        main = "Deep-sea O2 (\u00b5mol/kg) over CO2 and phosphate (as Pohl et al. 2021, Fig. 2)", font.main = 1, cex.main = 0.9)
  axis(1, at = seq_along(co2), labels = paste0("\u00d7", co2)); axis(2, at = seq_along(po4), labels = po4); box()
  for (i in seq_along(co2)) for (j in seq_along(po4)) { k <- g$co2 == co2[i] & g$po4 == po4[j]
    text(i, j, sprintf("%.0f%s", g$O2_deep[k], if (g$anoxic[k]) "\nanoxic" else ""), cex = 0.75, font = if (g$anoxic[k]) 2 else 1)
    text(i, j - 0.32, substr(g$state[k], 1, 1), cex = 0.6, col = "grey20") }
  ix <- function(v) match(v, co2); jy <- function(v) match(v, po4)
  if (all(c(24, 8.5) %in% co2) && 0.4 %in% po4) { arrows(ix(24) - 0.15, jy(0.4) + 0.12, ix(8.5) + 0.15, jy(0.4) + 0.12, lwd = 2, length = 0.1); text(mean(c(ix(24), ix(8.5))), jy(0.4) + 0.22, "1", font = 2) }
  if (all(c(7, 5) %in% co2) && all(c(0.2, 0.6) %in% po4)) { arrows(ix(7), jy(0.2) + 0.15, ix(5) + 0.1, jy(0.6) - 0.15, lwd = 2, length = 0.1); text(mean(c(ix(7), ix(5))) + 0.15, jy(0.4), "2", font = 2) }
  mtext("Letters: P-like / S-like / mixed / stagnant. Pohl et al. show the % of anoxic seafloor; the box model has one deep box, so its deep O2 is shown.", 1, 3.1, cex = 0.62, col = "grey30")
}

#' deep-water formation at each pole over a fine CO2 scan (the transitions of Pohl et al. 2021):
#' both poles -> south off -> north off and south on again, as CO2 falls
pohl2021_scan <- function(inp, co2 = 2^seq(log2(4), log2(28), length.out = 25), gradient = inp$gradient, par = list()) {
  par <- modifyList(modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), list(o2_pal = POHL2021$o2_pal, p_mean = POHL2021$po4_rel * POHL2021$po4_modern)), par)
  cells <- .ocean_cells(inp, par)
  bv <- .ocean_box_values(cells, par, lapply(co2, function(c2) .ocean_dmon(inp, c2, gradient)))
  r <- ocean_boxes(bv, par)
  out <- data.frame(co2 = co2, south = r$q_polarS, north = r$q_polarN, shelf = r$q_shelf, south_layer = r$layer_polarS %||% NA, north_layer = r$layer_polarN %||% NA,
                    O2_abyss = r$O2_deep, O2_int = r$O2_int %||% NA, ice_S = bv$polarS[, "ice"], ice_N = bv$polarN[, "ice"], renewal_yr = r$vent_age_yr)
  for (b in grep("^polarS_[0-9]+$", names(bv), value = TRUE)) { out[[b]] <- r[[paste0("q_", b)]]; out[[paste0("layer_", b)]] <- r[[paste0("layer_", b)]] }
  attr(out, "basins") <- attr(cells, "basins")
  out
}
plot_pohl2021_scan <- function(sc) {
  op <- par(mfrow = c(2, 1), mar = c(0.6, 4.4, 2.4, 4.4), mgp = c(2.4, 0.6, 0), cex = 0.85, las = 1); on.exit(par(op))
  x <- sc$co2; ym <- max(c(sc$south, sc$north, 1), na.rm = TRUE)
  plot(x, sc$north, log = "x", type = "n", ylim = c(0, ym * 1.12), xaxt = "n", xlab = "", ylab = "Deep-water formation (Sv)",
       main = "Deep-water formation at each pole as CO2 falls (right to left), as in Pohl et al. (2021)", font.main = 1, cex.main = 0.9)
  pch <- function(l) ifelse(l == "intermediate", 1, ifelse(l == "abyss", 19, NA))
  lines(x, sc$north, col = "#3f007d", lwd = 2); points(x, sc$north, pch = pch(sc$north_layer), col = "#3f007d")
  lines(x, sc$south, col = "#1b9e77", lwd = 2); points(x, sc$south, pch = pch(sc$south_layer), col = "#1b9e77")
  bn <- grep("^polarS_[0-9]+$", names(sc), value = TRUE); bi <- attr(sc, "basins")
  for (b in bn) { lines(x, sc[[b]], col = .src_col(b), lwd = 1, lty = 2); points(x, sc[[b]], pch = pch(sc[[paste0("layer_", b)]]), col = .src_col(b), cex = 0.6) }
  if (any(sc$shelf > 0.05)) lines(x, sc$shelf, col = "#e6550d", lwd = 1.2, lty = 2)
  abline(v = c(5, 8.5, 24), lty = 3, col = "grey60"); text(c(5, 8.5, 24), ym * 1.1, c("\u00d75", "\u00d78.5", "\u00d724"), cex = 0.7, col = "grey40")
  lb <- if (length(bn) && !is.null(bi)) bi$label[match(bn, bi$basin)] else character(0)
  legend("topleft", bty = "n", cex = 0.66, lwd = c(2, 2, rep(1, length(bn)), 1.2, NA, NA), lty = c(1, 1, rep(2, length(bn)), 2, NA, NA),
         pch = c(NA, NA, rep(NA, length(bn)), NA, 19, 1), col = c("#3f007d", "#1b9e77", vapply(bn, .src_col, ""), "#e6550d", "grey30", "grey30"),
         legend = c("north polar sea", "south polar sea (all)", lb, "shelf seas", "fills the abyss", "fills the intermediate layer"))
  par(mar = c(4, 4.4, 0.6, 4.4))
  plot(x, sc$O2_abyss, log = "x", type = "l", lwd = 2, col = "#7a0177", ylim = c(0, max(c(sc$O2_abyss, sc$O2_int), 10, na.rm = TRUE)), xlab = "CO2 (\u00d7 pre-industrial)", ylab = "O2 (\u00b5mol/kg)", xaxt = "n")
  if (any(is.finite(sc$O2_int))) lines(x, sc$O2_int, col = "#c51b8a", lty = 2)
  axis(1, at = c(4, 5, 6, 8.5, 10, 14, 20, 24, 28), labels = c(4, 5, 6, 8.5, 10, 14, 20, 24, 28))
  par(new = TRUE); plot(x, sc$ice_N, log = "x", type = "l", col = "#4a6fa5", lty = 3, axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
  lines(x, sc$ice_S, col = "#66c2a5", lty = 3); axis(4, col.axis = "#4a6fa5"); mtext("winter sea ice (fraction)", 4, 2.6, cex = 0.8, las = 0, col = "#4a6fa5")
  legend("topleft", bty = "n", cex = 0.72, lty = c(1, 2, 3, 3), col = c("#7a0177", "#c51b8a", "#4a6fa5", "#66c2a5"),
         legend = c("O2 abyss", "O2 intermediate layer", "sea ice north", "sea ice south"))
}

#' colour of a deep-water source (south basins: shades of green)
.src_col <- function(b) {
  if (startsWith(b, "polarS_")) { i <- as.integer(sub("polarS_", "", b)); return(c("#1b9e77", "#66a61e", "#00441b", "#41ae76", "#7fbc41", "#006d2c")[(i - 1) %% 6 + 1]) }
  if (b == "shelf_Baltic") return("#b2182b"); if (b == "shelf_Laurentia") return("#d01c8b")
  if (startsWith(b, "shelf_R")) { i <- as.integer(sub("shelf_R", "", b)); return(c("#e6550d", "#fd8d3c", "#a63603", "#fdae6b", "#d94801", "#8c2d04")[(i - 1) %% 6 + 1]) }
  c(shelf = "#e6550d", polarS = "#1b9e77", polarN = "#3f007d", open = "#636363")[[b]]
}

# ---- the south polar sea as basins ---------------------------------------------------------------
#' table of the south polar basins at one CO2: area, fresh water (P - E + rivers; mm/yr over the basin
#' and mSv), of which river water, winter temperature and sea ice
south_basin_table <- function(inp, co2 = inp$co2_pal, par = list(), gradient = inp$gradient) {
  par <- modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), par)
  cells <- .ocean_cells(inp, par); bi <- attr(cells, "basins"); if (is.null(bi)) stop("set south_basins > 1", call. = FALSE)
  bv <- .ocean_box_values(cells, par, list(.ocean_dmon(inp, co2, gradient)))
  d <- do.call(rbind, lapply(bi$basin, function(b) { v <- bv[[b]][1, ]
    data.frame(basin = b, label = bi$label[bi$basin == b], area_Mkm2 = v[["A"]] / 1e12, freshwater_mSv = v[["F"]] / 1e3,
               freshwater_mm_yr = v[["F"]] / max(v[["A"]], 1) * 3.15576e7 * 1000, rivers_mSv = v[["R"]] / 1e3,
               winter_T = v[["T"]], winter_ice = v[["ice"]]) }))
  d$co2 <- co2; d
}
#' formation of each south basin (and the north) against the river water into the south polar sea
#' (x0 ... x2 of the run) at a few CO2 levels: which basins switch on or off with more or less runoff
south_basin_runoff_scan <- function(inp, co2 = c(5, 8.5, 24), river = seq(0, 2, 0.25), par = list(), gradient = inp$gradient) {
  par <- modifyList(modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), list(o2_pal = POHL2021$o2_pal, p_mean = POHL2021$po4_rel * POHL2021$po4_modern)), par)
  out <- list()
  for (rv in river) { p2 <- par; p2$river_polarS <- rv
    cells <- .ocean_cells(inp, p2); bv <- .ocean_box_values(cells, p2, lapply(co2, function(c2) .ocean_dmon(inp, c2, gradient)))
    r <- ocean_boxes(bv, p2)
    for (b in c(grep("^polarS_[0-9]+$", names(bv), value = TRUE), "polarN"))
      out[[length(out) + 1]] <- data.frame(co2 = co2, river = rv, source = b, formation = r[[paste0("q_", b)]], layer = r[[paste0("layer_", b)]] %||% NA,
                                           O2_abyss = r$O2_deep, stringsAsFactors = FALSE)
    bi <- attr(cells, "basins") }
  d <- do.call(rbind, out); attr(d, "basins") <- bi; d
}
plot_south_basin_runoff <- function(d) {
  co2 <- sort(unique(d$co2), decreasing = TRUE); bi <- attr(d, "basins")
  op <- par(mfrow = c(1, length(co2)), mar = c(4, 4, 2.6, 0.8), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1); on.exit(par(op))
  ym <- max(d$formation, 1, na.rm = TRUE)
  for (c2 in co2) { e <- d[d$co2 == c2, ]
    plot(NA, xlim = range(e$river), ylim = c(0, ym * 1.05), xlab = "River water into the south polar sea (\u00d7 run)", ylab = "Deep-water formation (Sv)",
         main = sprintf("\u00d7%g CO2", c2), font.main = 2, cex.main = 1)
    abline(v = 1, lty = 3, col = "grey60")
    for (b in unique(e$source)) { f <- e[e$source == b, ]
      lines(f$river, f$formation, col = .src_col(b), lwd = if (b == "polarN") 2 else 1.6, lty = if (b == "polarN") 1 else 2)
      points(f$river, f$formation, pch = ifelse(f$layer == "abyss", 19, ifelse(f$layer == "intermediate", 1, NA)), col = .src_col(b), cex = 0.7) }
    if (c2 == co2[1]) legend("topright", bty = "n", cex = 0.7, lwd = 1.6, lty = c(1, rep(2, nrow(bi))), col = c("#3f007d", vapply(bi$basin, .src_col, "")),
                             legend = c("north polar sea", bi$label))
  }
}
#' map of the south polar basins: their area, fresh water (mm/yr; mSv of river water), winter temperature
#' and sea ice, and - for an ocean state - their deep-water formation
plot_south_basins <- function(inp, par = list(), co2 = inp$co2_pal, os = NULL, which = "thermal") {
  par <- modifyList(modifyList(OCEAN_PARAMS, POHL2021$par), par)
  cells <- .ocean_cells(inp, par); bi <- attr(cells, "basins"); if (is.null(bi)) stop("set south_basins > 1", call. = FALSE)
  tb <- south_basin_table(inp, co2, par)
  lon <- GRID1$lon; lat <- GRID1$lat; ylim <- c(-90, -par$polar_lat + 22)
  op <- par(mar = c(4.6, 3.6, 2.4, 1), mgp = c(2, 0.5, 0), cex = 0.85, las = 1); on.exit(par(op))
  plot(NA, xlim = c(-180, 180), ylim = ylim, xaxs = "i", yaxs = "i", asp = 1, axes = FALSE, xlab = "", ylab = "Palaeo-latitude")
  axis(1, at = seq(-180, 180, 60), labels = paste0(seq(-180, 180, 60), "\u00b0")); axis(2, at = seq(-90, ylim[2], 15))
  rect(-180, -90, 180, ylim[2], col = "#eef4fa", border = NA)
  for (b in bi$basin) { C <- cells[[b]]; Z <- matrix(NA, 180, 360); Z[C$k] <- 1
    image(lon, lat, t(Z), col = adjustcolor(.src_col(b), 0.55), add = TRUE, useRaster = FALSE) }
  if (!is.null(inp$riv)) { R <- inp$riv; R[inp$lf >= 0.5 | R <= 0.3] <- NA
    contour(lon, lat, t(ifelse(is.na(R), 0, R)), levels = c(1, 3), add = TRUE, drawlabels = FALSE, col = "#08519c", lwd = 0.7, lty = 1) }
  L <- matrix(NA, 180, 360); L[inp$lf >= 0.5] <- 1; image(lon, lat, t(L), col = "#cfc6b4", add = TRUE, useRaster = FALSE)
  contour(lon, lat, t(inp$lf), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#5b5246", lwd = 0.7)
  abline(h = -par$polar_lat, lty = 3, col = "grey40")
  for (i in seq_len(nrow(tb))) { b <- tb$basin[i]; x <- bi$lon_centre[bi$basin == b]; y <- bi$lat_centre[bi$basin == b]
    f <- if (!is.null(os)) os[[which]][[paste0("q_", b)]] else NA
    txt <- sprintf("%s\n%.1f Mkm\u00b2\nfresh water %.0f mm/yr\n(rivers %.0f mSv)\nwinter %.1f \u00b0C, ice %.0f %%%s", sub(" \\(.*", "", tb$label[i]), tb$area_Mkm2[i],
                   tb$freshwater_mm_yr[i], tb$rivers_mSv[i], tb$winter_T[i], 100 * tb$winter_ice[i], if (is.finite(f)) sprintf("\nsinks %.1f Sv", f) else "")
    text(x, y + 8, txt, cex = 0.62, font = 2, col = "grey10") }
  box()
  title(sprintf("South polar basins at \u00d7%g CO2 (split at the land barriers reaching furthest into the polar sea)", co2), font.main = 1, cex.main = 0.9)
  mtext("Blue lines: river plumes (\u2265 1 and 3 mm/day of river water). Fresh water = P - E + rivers over each basin (more fresh water per area = lighter surface water).", 1, 2.6, cex = 0.62, adj = 0)
}

#' great-circle distance (km) from grid cells to a point
.gc_km <- function(lat, lon, plat, plon) {
  p1 <- lat * pi / 180; p2 <- plat * pi / 180; dl <- (lon - plon) * pi / 180
  6371 * acos(pmin(1, pmax(-1, sin(p1) * sin(p2) + cos(p1) * cos(p2) * cos(dl))))
}
#' connected pieces (8 neighbours, longitude wraps) of a logical 180 x 360 mask
.connected_ids <- function(M) {
  id <- matrix(0L, 180, 360); nid <- 0L
  for (c0 in which(M)) { if (id[c0] > 0) next
    nid <- nid + 1L; queue <- c0; id[c0] <- nid
    while (length(queue)) { c1 <- queue[1]; queue <- queue[-1]
      r <- (c1 - 1) %% 180 + 1; cc <- (c1 - 1) %/% 180 + 1
      for (dr in -1:1) for (dc in -1:1) { r2 <- r + dr; if (r2 < 1 || r2 > 180) next; c2 <- ((cc + dc - 1) %% 360) + 1; k2 <- r2 + (c2 - 1) * 180
        if (M[k2] && id[k2] == 0L) { id[k2] <- nid; queue <- c(queue, k2) } } } }
  id
}

# ---- contributions of the deep-water sources ---------------------------------------------------
#' table of the deep-water sources of an ocean state: area, winter temperature and (model) salinity, density,
#' formation (Sv) and share, the layer it fills, O2 at saturation of its water (warm, saline water holds less),
#' and the O2 it carries into the deep sea (Tmol/yr) and its share
ocean_source_table <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); r <- os[[which]]; par <- os$par
  cells <- .ocean_cells(os$inputs, par); bv <- os$boxes; src <- .src_names(bv)
  lab <- c(polarS = "south polar sea", polarN = "north polar sea", shelf = "evaporative shelf seas")
  bi <- attr(cells, "basins"); li <- attr(cells, "lowlat")
  if (!is.null(bi)) lab[bi$basin] <- paste("south polar", bi$label)
  if (!is.null(li)) lab[li$source] <- li$label
  Sv <- 1e6 * 3.15576e7
  d <- do.call(rbind, lapply(src, function(b) { v <- bv[[b]][1, ]; q <- r[[paste0("q_", b)]] %||% 0; S <- r[[paste0("S_", b)]] %||% NA
    fi <- if (startsWith(b, "polar")) v[["ice"]] else 0
    o2 <- o2_saturation(v[["T"]], S, par$o2_pal) * (1 - if (startsWith(b, "polar")) par$seaice_o2 * fi else 0)
    data.frame(source = b, label = lab[[b]] %||% b, area_Mkm2 = v[["A"]] / 1e12, winter_T = v[["T"]], salinity = S,
               sigma = r[[paste0("rho_", b)]] %||% NA, freshwater_mm_yr = v[["F"]] / max(v[["A"]], 1) * 3.15576e7 * 1000,
               formation_Sv = q, layer = r[[paste0("layer_", b)]] %||% (if (q > 0.05) "deep" else "none"), O2_sat = o2,
               O2_flux_Tmol_yr = q * Sv * 1025 * o2 * 1e-6 / 1e12,
               salinity_to_sink = 35 + ((min(r$rho_deep, r$rho_int %||% Inf) + 1000) / 1027 - 1 + 1.67e-4 * (v[["T"]] - 10) + 0.5e-5 * (v[["T"]] - 10)^2) / 7.8e-4,
               stringsAsFactors = FALSE) }))
  d$share_formation <- d$formation_Sv / max(sum(d$formation_Sv), 1e-9)
  d$share_O2 <- d$O2_flux_Tmol_yr / max(sum(d$O2_flux_Tmol_yr), 1e-9)
  d$O2_per_Sv_vs_polar <- d$O2_sat / max(d$O2_sat[startsWith(d$source, "polar") & d$formation_Sv > 0.05], d$O2_sat[startsWith(d$source, "polar")], na.rm = TRUE)
  d
}
#' figure of the contributions: share of the deep-water formation and of the O2 carried down, per source,
#' with each source's O2 at saturation
plot_ocean_sources <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); d <- ocean_source_table(os, which); d <- d[order(-d$formation_Sv, d$source), ]
  op <- par(mar = c(4.2, 13, 2.6, 6), mgp = c(2.3, 0.6, 0), cex = 0.82, las = 1); on.exit(par(op))
  n <- nrow(d); y <- rev(seq_len(n)); cols <- vapply(d$source, .src_col, "")
  plot(NA, xlim = c(0, max(100 * c(d$share_formation, d$share_O2), 10) * 1.05), ylim = c(0.4, n + 0.6), yaxt = "n", xlab = "Share of the deep-sea renewal (%)", ylab = "",
       main = sprintf("Deep-water sources \u2013 %s (%s start): %.1f Sv in all", os[[which]]$state, if (which == "thermal") "polar" else "shelf", sum(d$formation_Sv)), font.main = 1, cex.main = 0.95)
  rect(0, y - 0.38, 100 * d$share_formation, y + 0.02, col = cols, border = NA)
  rect(0, y - 0.02, 100 * d$share_O2, y + 0.38, col = adjustcolor(cols, 0.4), border = cols)
  axis(2, at = y, labels = sprintf("%s\n%.1f Mkm\u00b2, %.1f \u00b0C, %.1f psu", d$label, d$area_Mkm2, d$winter_T, d$salinity), cex.axis = 0.7, tick = FALSE)
  par(xpd = TRUE)
  text(par("usr")[2], y, sprintf("%.1f Sv\nO2 sat. %.0f", d$formation_Sv, d$O2_sat), adj = -0.1, cex = 0.68)
  par(xpd = FALSE)
  legend("bottomright", bty = "n", cex = 0.7, fill = c("grey40", adjustcolor("grey40", 0.4)), border = c(NA, "grey40"),
         legend = c("share of the formation (Sv)", "share of the O2 carried down"))
  mtext("O2 sat.: O2 (\u00b5mol/kg) of the sinking water at saturation (Garcia & Gordon 1992) for its winter temperature and box salinity - warm, saline water carries less.", 1, 3, cex = 0.6, adj = 1)
}

# ---- behaviour of a basin: stratified by runoff, or a downwelling centre? ---------------------------
#' the behaviour of a low-latitude basin (named basin, e.g. "shelf_Baltic", or a shelf region) in an
#' ocean state: its fresh-water budget (P - E and river water, mm/yr over the basin; share of the area
#' that gains fresh water), its salinity and winter density in the box model against the open sea and the
#' deep layers, the wind-driven vertical velocity, and the same box model without its river water.
#' Classes: deep-water source (downwelling centre) / anti-estuarine (dense outflow, lighter than the deep
#' water) / stratified by river runoff (estuarine only because of its rivers) / stratified by rain and runoff
basin_behaviour <- function(os, basin = "shelf_Baltic", which = c("thermal", "haline")) {
  which <- match.arg(which); r <- os[[which]]; par <- os$par; inp <- os$inputs
  cells <- .ocean_cells(inp, par); C <- cells[[basin]]
  if (is.null(C)) stop("no box '", basin, "' (named basins: ", paste(grep("^shelf_", names(cells), value = TRUE), collapse = ", "), ")", call. = FALSE)
  A <- C$A; mm <- function(x) sum(A * x, na.rm = TRUE) / sum(A) * 365.25
  riv <- ifelse(is.finite(C$riv), pmax(C$riv, 0), 0); pme <- C$fw - riv
  bv <- os$boxes[[basin]][1, ]
  q <- r[[paste0("q_", basin)]] %||% 0; Sb <- r[[paste0("S_", basin)]]; sb <- r[[paste0("rho_", basin)]]
  s_open <- rho_sw(os$boxes$open[1, "Tann"], r$S_open) - 1000
  s_int <- r$rho_int %||% NA; s_deep <- r$rho_deep; s_target <- min(s_deep, s_int, na.rm = TRUE)
  # the same without the basin's river water
  nm <- sub("shelf_", "", basin); p0 <- par
  if (startsWith(basin, "shelf_") && nm %in% names(par$river_regions)) p0$river_regions[[nm]] <- 0 else if (startsWith(basin, "shelf_")) p0$river_regions <- c(par$river_regions, setNames(0, nm))
  bv0 <- .ocean_box_values(.ocean_cells(inp, p0), p0); r0 <- ocean_boxes(bv0, p0, start = which)
  w <- if (!is.null(inp$w)) inp$w[C$k] else NA
  out <- list(basin = basin, area_Mkm2 = sum(A) / 1e12, lat = range(C$lat), lon = range(GRID1$lon[(C$k - 1) %/% 180 + 1]),
              P_minus_E_mm_yr = mm(pme), rivers_mm_yr = mm(riv), net_mm_yr = mm(C$fw), rivers_mSv = sum(A * riv) / 1000 / 86400 / 1e3,
              net_mSv = bv[["F"]] / 1e3, gaining_area = sum(A[C$fw > 0]) / sum(A), runoff_gaining_area = sum(A[C$fw > 0 & pme <= 0]) / sum(A),
              winter_T = bv[["T"]], S_box = Sb, sigma_box = sb, sigma_open = s_open, sigma_int = s_int, sigma_deep = s_deep,
              S_to_sink = 35 + ((s_target + 1000) / 1027 - 1 + 1.67e-4 * (bv[["T"]] - 10) + 0.5e-5 * (bv[["T"]] - 10)^2) / 7.8e-4,
              formation_Sv = q, layer = r[[paste0("layer_", basin)]] %||% NA,
              S_box_norivers = r0[[paste0("S_", basin)]], sigma_box_norivers = r0[[paste0("rho_", basin)]], formation_norivers_Sv = r0[[paste0("q_", basin)]] %||% 0,
              w_mean_m_day = mean(w, na.rm = TRUE), downwelling_area = sum(A[is.finite(w) & w < 0]) / sum(A))
  F <- out$net_mm_yr
  out$river_offset <- if (out$P_minus_E_mm_yr < 0) out$rivers_mm_yr / -out$P_minus_E_mm_yr else NA
  out$runoff_note <- if (is.finite(out$river_offset)) sprintf("its rivers offset %.0f %% of its net evaporation (P - E %.0f, rivers +%.0f mm/yr)", 100 * out$river_offset, out$P_minus_E_mm_yr, out$rivers_mm_yr) else
    sprintf("P - E %+.0f, rivers +%.0f mm/yr", out$P_minus_E_mm_yr, out$rivers_mm_yr)
  out$class <- if (q > 0.05) sprintf("deep-water source (downwelling centre): %.1f Sv into the %s", q, if (is.na(out$layer)) "deep sea" else out$layer) else
    if (F > 0 && out$P_minus_E_mm_yr <= 0) "stratified by river runoff: estuarine only because of its rivers (evaporative without them)" else
    if (F > 0) "stratified by rain and runoff: estuarine even without its rivers" else
    if (sb > s_open) sprintf("anti-estuarine (evaporative): its saltier water sinks below the open-sea surface water and flows out at depth, but is lighter than the deep water (would need %.1f psu)", out$S_to_sink) else
    "anti-estuarine, but too warm to be denser than the open-sea surface water"
  out$class_norivers <- if (out$formation_norivers_Sv > 0.05) sprintf("without its rivers: deep-water source (%.1f Sv)", out$formation_norivers_Sv) else
    sprintf("without its rivers: no deep-water formation (%.2f psu, \u03c3 %.2f)", out$S_box_norivers, out$sigma_box_norivers)
  out$ekman <- if (is.finite(out$w_mean_m_day)) sprintf("wind-driven %s (%.2f m/day on average; %.0f %% of the basin downwelling)", if (out$w_mean_m_day < 0) "downwelling" else "upwelling",
                                                        out$w_mean_m_day, 100 * out$downwelling_area) else ""
  out
}
#' behaviour of all low-latitude boxes (named basins and shelf regions) as a table
basin_behaviour_table <- function(os, which = c("thermal", "haline")) {
  which <- match.arg(which); cells <- .ocean_cells(os$inputs, os$par); li <- attr(cells, "lowlat")
  bs <- intersect(li$source %||% character(0), names(cells))
  do.call(rbind, lapply(bs, function(b) { x <- basin_behaviour(os, b, which)
    data.frame(basin = li$label[li$source == b], area_Mkm2 = x$area_Mkm2, P_minus_E = x$P_minus_E_mm_yr, rivers = x$rivers_mm_yr, net = x$net_mm_yr,
               gaining_area = x$gaining_area, winter_T = x$winter_T, S = x$S_box, sigma = x$sigma_box, S_to_sink = x$S_to_sink,
               formation_Sv = x$formation_Sv, S_norivers = x$S_box_norivers, formation_norivers_Sv = x$formation_norivers_Sv,
               w_m_day = x$w_mean_m_day, river_offset = x$river_offset, class = x$class, stringsAsFactors = FALSE) }))
}
#' figure of a basin's behaviour: map of its fresh-water flux (P - E + rivers) with river plumes and the
#' box outline; its budget; and its winter density against the open sea and the deep layers, with and
#' without its river water
plot_basin_behaviour <- function(os, basin = "shelf_Baltic", which = c("thermal", "haline")) {
  which <- match.arg(which); x <- basin_behaviour(os, basin, which); inp <- os$inputs
  lon <- GRID1$lon; lat <- GRID1$lat
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(c(1, 1, 2, 3), 2, 2), widths = c(1.35, 1)); par(oma = c(4.4, 0, 0, 0))
  cl <- mean(x$lat); co <- mean(x$lon); xl <- co + c(-40, 40); yl <- pmax(pmin(cl + c(-25, 25), 90), -90)
  par(mar = c(3, 3.4, 2.6, 1), mgp = c(1.9, 0.5, 0), cex = 0.8, las = 1)
  plot(NA, xlim = xl, ylim = yl, asp = 1, xaxs = "i", yaxs = "i", xlab = "", ylab = "Palaeo-latitude")
  F <- inp$fw * 365.25; F[inp$lf >= 0.5] <- NA; lim <- 1500
  cols <- grDevices::colorRampPalette(c("#8c510a", "#d8b365", "#f6e8c3", "#f5f5f5", "#c7eae5", "#5ab4ac", "#01665e"))(21)
  image(lon, lat, t(pmin(pmax(F, -lim), lim)), breaks = seq(-lim, lim, length.out = 22), col = cols, add = TRUE, useRaster = FALSE)
  if (!is.null(inp$riv)) contour(lon, lat, t(ifelse(is.finite(inp$riv), inp$riv, 0)), levels = c(1, 3, 10), add = TRUE, drawlabels = FALSE, col = "#08519c", lwd = c(0.6, 1, 1.4))
  L <- matrix(NA, 180, 360); L[inp$lf >= 0.5] <- 1; image(lon, lat, t(L), col = "#cfc6b4", add = TRUE, useRaster = FALSE)
  contour(lon, lat, t(inp$lf), levels = 0.5, add = TRUE, drawlabels = FALSE, col = "#5b5246", lwd = 0.7)
  C <- .ocean_cells(inp, os$par)[[basin]]; Mb <- matrix(0, 180, 360); Mb[C$k] <- 1
  contour(lon, lat, t(Mb), levels = 0.5, add = TRUE, drawlabels = FALSE, col = .src_col(basin), lwd = 2)
  P <- inp$regions[[sub("shelf_", "", basin)]]; if (!is.null(P)) { points(P$plon, P$plat, pch = 21, bg = "white", cex = 0.8); text(P$plon, P$plat, P$name, pos = 4, cex = 0.55) }
  box(); title(sprintf("%s: fresh water into the sea (P - E + rivers, mm/yr)", sub("shelf_", "", basin)), font.main = 1, cex.main = 0.9)
  # budget
  par(mar = c(3, 4.2, 2.6, 1))
  v <- c(`P - E` = x$P_minus_E_mm_yr, rivers = x$rivers_mm_yr, net = x$net_mm_yr)
  bp <- barplot(v, col = c("#8c510a", "#2171b5", if (x$net_mm_yr > 0) "#01665e" else "#8c510a"), ylab = "mm/yr over the basin", ylim = range(0, v) * 1.25,
                main = sprintf("Budget: %.0f %% of the basin gains fresh water", 100 * x$gaining_area), font.main = 1, cex.main = 0.9)
  abline(h = 0); text(bp, v, sprintf("%.0f", v), pos = ifelse(v >= 0, 3, 1), cex = 0.75)
  # density ladder
  par(mar = c(5.6, 4.2, 2.6, 1))
  d <- unlist(lapply(list(`basin` = x$sigma_box, `basin, no rivers` = x$sigma_box_norivers, `open sea surface` = x$sigma_open, intermediate = x$sigma_int, `deep water` = x$sigma_deep),
                     function(v) unname(as.numeric(v)[1])))
  d <- d[is.finite(d)]
  plot(seq_along(d), d, pch = 19, cex = 1.4, col = c(.src_col(basin), adjustcolor(.src_col(basin), 0.5), "grey50", "#9e9ac8", "#3f007d")[seq_along(d)],
       xaxt = "n", xlim = c(0.5, length(d) + 0.5), ylim = range(d) + c(-0.5, 0.5), xlab = "", ylab = "Winter \u03c3 (kg/m\u00b3)", main = "Is it dense enough to sink?", font.main = 1, cex.main = 0.9)
  axis(1, at = seq_along(d), labels = names(d), cex.axis = 0.7, las = 2); abline(h = x$sigma_deep, lty = 2, col = "#3f007d")
  text(seq_along(d), d, sprintf("%.2f", d), pos = 4, cex = 0.7)
  tl <- strwrap(paste0(sub("shelf_", "", basin), " basin: ", x$class, "; ", x$runoff_note, ". ", x$class_norivers, ". ", x$ekman, "."), 170)
  for (i in seq_along(tl)) mtext(tl[i], 1, 0.2 + 0.85 * (i - 1), outer = TRUE, cex = 0.62, adj = 0.01, font = if (i == 1) 2 else 1)
  mtext("Map: brown = net evaporation (anti-estuarine), green = net fresh water (estuarine); blue lines: river water 1, 3, 10 mm/day; outline = the basin box; dots = its sites.",
        1, 0.2 + 0.85 * length(tl), outer = TRUE, cex = 0.6, adj = 0.01, col = "grey30")
}

#' the behaviour of a basin over CO2 (sea surface shifted to each CO2, as in the regime map): budget,
#' salinity, formation, with and without its rivers (box model from a polar start)
basin_behaviour_scan <- function(inp, basin = "shelf_Baltic", co2 = c(2, 4, 6, 8.5, 12, 24), par = list()) {
  par <- modifyList(OCEAN_PARAMS, par); nm <- sub("shelf_", "", basin)
  p0 <- par; p0$river_regions <- c(par$river_regions[setdiff(names(par$river_regions), nm)], setNames(0, nm))
  dm <- lapply(co2, function(c2) .ocean_dmon(inp, c2, inp$gradient))
  bv <- .ocean_box_values(.ocean_cells(inp, par), par, dm); r <- ocean_boxes(bv, par)
  bv0 <- .ocean_box_values(.ocean_cells(inp, p0), p0, dm); r0 <- ocean_boxes(bv0, p0)
  B <- bv[[basin]]; if (is.null(B)) stop("no box '", basin, "'", call. = FALSE)
  mmyr <- function(v) v / B[, "A"] * 3.15576e7 * 1000
  tgt <- pmin(r$rho_deep, r$rho_int %||% Inf)
  data.frame(co2 = co2, state = r$state, winter_T = B[, "T"], P_minus_E = mmyr(B[, "F"] - B[, "R"]), rivers = mmyr(B[, "R"]), net = mmyr(B[, "F"]),
             S = r[[paste0("S_", basin)]], S_to_sink = 35 + ((tgt + 1000) / 1027 - 1 + 1.67e-4 * (B[, "T"] - 10) + 0.5e-5 * (B[, "T"] - 10)^2) / 7.8e-4,
             formation_Sv = r[[paste0("q_", basin)]], S_norivers = r0[[paste0("S_", basin)]], formation_norivers_Sv = r0[[paste0("q_", basin)]])
}

# ---- scenarios for one basin: when is it a downwelling centre, when (partly) stratified? -------------------
BASIN_SCENARIOS <- list(
  base = list(label = "as the run", par = list()),
  lid = list(label = "melt-water lid (P/S switch at ~9 PAL)", par = list(seaice_effect = -1, seaice_psu = 2)),
  lid_cal = list(label = "melt-water lid 0.95 psu (calibrated: P/S switch at 6-8 PAL in the 435 Ma SPEEDY runs)", par = list(seaice_effect = -1, seaice_psu = 0.95)),
  warm_poles = list(label = "warmer poles (gradient factor 0.7)", par = list(), gradient = 0.7),
  restricted = list(label = "restricted exchange with the open sea (\u00d70.5)", par = list(mix_shelf = 0.05)),
  open = list(label = "open exchange with the open sea (\u00d72)", par = list(mix_shelf = 0.2)))
#' a basin over CO2 x river water for a set of scenarios and box sizes. Classes per setting:
#'   downwelling centre (sinks >= 0.02 Sv), stratified (net fresh-water gain: estuarine), partly stratified
#'   (>= 25 % of the basin under a fresh-water lid, i.e. gains fresh water there), anti-estuarine (not sinking)
basin_scenarios <- function(inp, basin = "shelf_Baltic", co2 = c(2, 4, 6, 8.5, 12, 24), river = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3),
                            scenarios = BASIN_SCENARIOS, radius = c(300, 700), par = list(), progress = NULL) {
  nm <- sub("shelf_", "", basin); out <- list(); n <- 0; ntot <- length(scenarios) * length(radius) * length(river)
  for (sn in names(scenarios)) for (rk in radius) {
    sc <- scenarios[[sn]]; g <- sc$gradient %||% inp$gradient
    p <- modifyList(modifyList(OCEAN_PARAMS, par), c(sc$par, list(region_km = rk)))
    dm <- lapply(co2, function(c2) .ocean_dmon(inp, c2, g))
    cells <- .ocean_cells(inp, p); C <- cells[[basin]]; if (is.null(C)) next
    riv <- ifelse(is.finite(C$riv), pmax(C$riv, 0), 0); pme <- C$fw - riv
    dTa <- vapply(dm, function(D) mean(colMeans(D)[C$row]), 0)                # basin-mean warming for the 7 %/K scaling
    for (rv in river) {
      p$river_regions <- c(p$river_regions[setdiff(names(p$river_regions), nm)], setNames(rv, nm))
      bv <- .ocean_box_values(cells, p, dm); r <- ocean_boxes(bv, p); B <- bv[[basin]]
      lid_area <- vapply(seq_along(co2), function(j) { f <- (pme + rv * riv) * (1 + p$cc * dTa[j]); sum(C$A[f > 0]) / sum(C$A) }, 0)
      net <- B[, "F"] / B[, "A"] * 3.15576e7 * 1000; q <- r[[paste0("q_", basin)]]
      tgt <- pmin(r$rho_deep, r$rho_int %||% Inf)
      cls <- ifelse(q >= 0.02, ifelse(lid_area >= 0.25, "downwelling centre, partly stratified", "downwelling centre"),
                    ifelse(net > 0, "stratified (estuarine)", ifelse(lid_area >= 0.25, "partly stratified", "anti-estuarine, not sinking")))
      out[[length(out) + 1]] <- data.frame(scenario = sn, scenario_label = sc$label, radius_km = rk, co2 = co2, river = rv, area_Mkm2 = sum(C$A) / 1e12,
                                           net_mm_yr = net, freshlid_area = lid_area, winter_T = B[, "T"], S = r[[paste0("S_", basin)]],
                                           S_to_sink = 35 + ((tgt + 1000) / 1027 - 1 + 1.67e-4 * (B[, "T"] - 10) + 0.5e-5 * (B[, "T"] - 10)^2) / 7.8e-4,
                                           formation_Sv = q, ocean_state = r$state, class = cls, stringsAsFactors = FALSE)
      n <- n + 1; if (is.function(progress)) progress(n / ntot, sprintf("%s, %g km, rivers x%g", sn, rk, rv))
    }
  }
  d <- do.call(rbind, out); attr(d, "basin") <- basin; d
}
BASIN_CLASS_COLS <- c("downwelling centre" = "#3f007d", "downwelling centre, partly stratified" = "#807dba",
                      "anti-estuarine, not sinking" = "#fdae6b", "partly stratified" = "#6baed6", "stratified (estuarine)" = "#08519c")
#' figure: one panel per scenario and box size, class over CO2 x river water; numbers: sinking (Sv) or the salinity
#' still missing to sink (psu)
plot_basin_scenarios <- function(d) {
  sn <- unique(d$scenario); rk <- sort(unique(d$radius_km)); co2 <- sort(unique(d$co2)); rv <- sort(unique(d$river))
  op <- par(mfrow = c(length(rk), length(sn)), mar = c(3.4, 3.4, 2.6, 0.6), oma = c(4.2, 0, 2, 0), mgp = c(2, 0.5, 0), cex = 0.62, las = 1); on.exit(par(op))
  lv <- names(BASIN_CLASS_COLS)
  for (k in rk) for (s in sn) { e <- d[d$scenario == s & d$radius_km == k, ]
    Z <- matrix(match(e$class[order(e$river, e$co2)], lv), length(co2))
    image(seq_along(co2), seq_along(rv), Z, breaks = seq(0.5, length(lv) + 0.5, 1), col = BASIN_CLASS_COLS, axes = FALSE,
          xlab = "CO2 (PAL)", ylab = "River water (\u00d7 run)", main = sprintf("%s\n%g km (%.1f Mkm\u00b2)", e$scenario_label[1], k, e$area_Mkm2[1]), font.main = 1, cex.main = 0.95)
    axis(1, at = seq_along(co2), labels = co2); axis(2, at = seq_along(rv), labels = rv); box()
    abline(h = match(1, rv), lty = 3, col = "white")
    for (i in seq_len(nrow(e))) { x <- match(e$co2[i], co2); y <- match(e$river[i], rv)
      text(x, y, if (e$formation_Sv[i] >= 0.02) sprintf("%.2f", e$formation_Sv[i]) else { ds <- e$S_to_sink[i] - e$S[i]; if (ds > 0.05) sprintf("-%.1f", ds) else if (e$formation_Sv[i] > 0.001) "<0.02" else "0" }, cex = 0.75,
           col = if (grepl("downwelling|estuarine", e$class[i]) && !grepl("anti", e$class[i])) "white" else "black") }
  }
  mtext(sprintf("%s: when is it a downwelling centre, when (partly) stratified?", sub("shelf_", "", attr(d, "basin") %||% "")), 3, 0.4, outer = TRUE, font = 2, cex = 0.8)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", inset = c(0, 0.035), horiz = TRUE, bty = "n", fill = BASIN_CLASS_COLS, border = "grey40", legend = lv, cex = 0.95, xpd = NA)
  par(op); par(oma = c(4.2, 0, 2, 0))
  mtext("Numbers: sinking (Sv), or the salinity still missing to sink (psu, negative). Dotted line: the run's river water. Partly stratified: \u2265 25 % of the basin gains fresh water (fresh-water lid).",
        1, 2.2, outer = TRUE, cex = 0.62)
}
