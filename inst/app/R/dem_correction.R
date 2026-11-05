# =============================================================================
#  STEP 6 (R) - Fossil-corrected PaleoDEMs
#
#  The PaleoDEMs (Scotese & Wright 2018) put 9-27 % of Silurian marine records on
#  land. correct_dem() builds a corrected grid for one map age:
#   1. land/sea from the PBDB-adjusted PaleoCoastlines (Kocsis & Scotese 2021)
#   2. local flooding around every marine / evaporite record still on land
#   3. seaways: inland marine areas joined to the ocean along the lowest route,
#      gulf/funnel shaped (widen, never narrow, towards the ocean), meandering
#   4. clean-up: no cliffs (slope limits), islets cut off by seaways drowned,
#      stray pools removed, no spits/barriers, smooth 0 m coastline
#  All records come from records_for_age() (database records + your points and
#  edits), so adding / removing points changes the corrected maps.
#
#  correct_dem() is pure R (no Python) so ages can run in parallel.
# =============================================================================

# ------------------------------------------------------------------ helpers ----
.cos_clip <- function(lat, lo = 0.05) pmin(pmax(cos(pi / 180 * lat), lo), 1)

#' grid (row, col) positions -> local km coordinates (east-west scaled by cos lat)
.to_km <- function(i, j, lat) {
  dl <- abs(lat[2] - lat[1]); n <- length(lat)
  cbind(i * dl * KM_PER_DEG, j * dl * KM_PER_DEG * .cos_clip(lat[pmin(pmax(round(i), 1), n)]))
}
.unwrap <- function(x, p) { if (length(x) < 2) return(x); d <- diff(x); d <- d - p * round(d / p); c(x[1], x[1] + cumsum(d)) }
.wrapj <- function(j, p) ((j - 1) %% p) + 1

#' smooth random field (sd 1) that is isotropic on the sphere (sum of random
#' plane waves in 3-D): no north-south streaks at high latitude. 300-1500 km.
fractal_noise <- function(lat, lon, seed) {
  set.seed(seed)
  LA <- matrix(pi / 180 * lat, length(lat), length(lon)); LO <- matrix(pi / 180 * lon, length(lat), length(lon), byrow = TRUE)
  X <- cos(LA) * cos(LO); Y <- cos(LA) * sin(LO); Z <- sin(LA); rm(LA, LO)
  f <- matrix(0, length(lat), length(lon))
  for (s in list(c(300, 0.5), c(700, 1), c(1500, 1))) {
    k <- 2 * pi * 6371 / s[1]
    for (q in 1:40) {
      v <- rnorm(3); v <- v / sqrt(sum(v^2))
      f <- f + (s[2] / sqrt(40)) * cos(k * (v[1] * X + v[2] * Y + v[3] * Z) + runif(1, 0, 2 * pi))
    }
  }
  (f - mean(f)) / sd(as.vector(f))
}

#' land/sea mask of the PaleoCoastlines CS polygons (TRUE = land) on the DEM grid
cs_land_mask <- function(age, lat, lon, cfg = SILUCARB) {
  dl <- abs(lat[2] - lat[1])
  f <- file.path(cfg$dir_coast, "CS", paste0(age, "Ma_CS_v7.shp"))
  if (!file.exists(f)) return(NULL)                 # no PaleoCoastlines for this age
  v <- vect(f)
  r <- rast(nrows = length(lat), ncols = length(lon), xmin = lon[1] - dl / 2, xmax = lon[length(lon)] + dl / 2,
            ymin = lat[1] - dl / 2, ymax = lat[length(lat)] + dl / 2, crs = "EPSG:4326")
  m <- as.matrix(rasterize(v, r, field = 1, background = 0), wide = TRUE)
  m[nrow(m):1, ] == 1
}

# --------------------------------------------------------------- clean-up ----
#' remove cliffs around everything the correction changed: land <= 1 m + slope x
#' distance to new sea; sea floor >= shelf (-30 m - 15 m/km x d) then continental
#' slope (25 m/km) from the changed coasts. Land only lowered, sea only raised.
relax_banks <- function(zc, z0, lat, cfg) {
  W <- area_matrix(lat, ncol(zc))
  newsea <- zc < 0 & z0 >= 0; newland <- zc >= 0 & z0 < 0
  d_land <- dist_km(newsea, lat)
  land <- zc >= 0; cap <- 1 + cfg$MAX_BANK_SLOPE * d_land
  low <- land & zc > cap
  out <- zc; out[low] <- pmax(cap[low], 1)
  d_sea <- dist_km(newsea | newland, lat)
  d1 <- (1400 - 30) / cfg$MAX_SEA_SLOPE
  floor_ <- ifelse(d_sea <= d1, -30 - cfg$MAX_SEA_SLOPE * d_sea, -1400 - cfg$MAX_DEEP_SLOPE * (d_sea - d1))
  floor_ <- smooth_km(pmax(floor_, -6000), lat, 25, nbands = 45)   # smooth: no cell-by-cell steps
  up <- out < 0 & out < floor_
  out[up] <- pmin(floor_[up], -1)
  list(z = out, stats = list(
    land_lowered_Mkm2 = round(sum(W[low]), 2),
    land_lowered_mean_m = round(sum(((zc - out) * W)[low]) / max(sum(W[low]), 1e-9)),
    seafloor_raised_Mkm2 = round(sum(W[up]), 2),
    seafloor_raised_mean_m = round(sum(((out - zc) * W)[up]) / max(sum(W[up]), 1e-9))))
}

#' cone slope limiter on land: z(x) <= z(y) + smax * dist(x, y) for neighbours,
#' only lowering land (never below +1 m) in and around the changed zone; the zone
#' grows with the lowering so a cascade never stops at a fixed edge (no new cliff).
#' Exact Jacobi iteration, evaluated only where a neighbour changed.
limit_land_slopes <- function(zc, z0, lat, smax, iters = 1000) {
  n <- nrow(zc); p <- ncol(zc); dy <- abs(lat[2] - lat[1]) * KM_PER_DEG
  dx <- dy * .cos_clip(lat, 0.02)
  z <- zc
  zone <- edt(abs(zc - z0) < 1) <= 6
  zone <- zone | dilate(abs(z - z0) > 1e-6, 2)
  nb <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  neigh <- function(idx, di, dj) {
    i <- (idx - 1) %% n + 1 + di; j <- ((idx - 1) %/% n + dj) %% p
    ifelse(i >= 1 & i <= n, i + j * n, NA_integer_)
  }
  ring <- function(idx, r) {
    out <- idx
    for (di in -r:r) for (dj in -r:r) if (di || dj) out <- c(out, neigh(idx, di, dj))
    unique(out[!is.na(out)])
  }
  active <- which(zone & z >= 0)
  for (it in seq_len(iters)) {
    if (!length(active)) break
    i_act <- (active - 1) %% n + 1
    cap <- rep(Inf, length(active))
    for (k in 1:8) {
      q <- neigh(active, nb[k, 1], nb[k, 2])
      d <- sqrt((dy * abs(nb[k, 1]))^2 + (dx[i_act] * abs(nb[k, 2]))^2)
      v <- ifelse(is.na(q), Inf, z[ifelse(is.na(q), 1L, q)] + smax * d)
      cap <- pmin(cap, v)
    }
    newv <- pmax(cap, 1)
    low <- z[active] >= 0 & z[active] > newv + 1e-6
    if (!any(low)) break
    lowered <- active[low]
    z[lowered] <- newv[low]
    grow <- ring(lowered, 2); zone[grow] <- TRUE
    act <- ring(lowered, 1); act <- unique(c(act, grow))
    active <- act[zone[act] & z[act] >= 0]
  }
  z
}

#' small land patches cut off by new seas become shallow sea (islands that existed
#' before the seaways and islets with an unambiguous land record are kept)
drown_islets <- function(zc, z0, lat, lon, cfg, protect = NULL, zref = NULL) {
  W <- area_matrix(lat, ncol(zc))
  lab <- label_wrap(zc >= 0, 8)
  if (!any(lab > 0)) return(list(z = zc, n = 0L, area = 0))
  nl <- max(lab)
  size <- numeric(nl); s <- tapply(as.vector(W), as.vector(lab), sum)
  ids <- as.integer(names(s)); size[ids[ids > 0]] <- s[ids > 0] * 1e6
  newsea <- zc < 0 & z0 >= 0
  nbm <- dilate(newsea, 1)
  touch <- logical(nl); tl <- unique(lab[nbm & lab > 0]); touch[tl] <- TRUE
  small <- size < cfg$ISLET_MAX_KM2 & touch
  if (!is.null(zref) && any(small)) {
    lr <- label_wrap(zref >= 0, 8)
    size_ref <- numeric(max(lr, 1)); s2 <- tapply(as.vector(W), as.vector(lr), sum)
    id2 <- as.integer(names(s2)); size_ref[id2[id2 > 0]] <- s2[id2 > 0] * 1e6
    sel <- small[pmax(lab, 1)] & lab > 0 & lr > 0
    pr <- unique(data.frame(k = lab[sel], par = lr[sel]))
    for (k in which(small)) {
      par <- pr$par[pr$k == k]
      if (!length(par) || max(size_ref[par]) < 2 * size[k]) small[k] <- FALSE
    }
  }
  if (!is.null(protect) && nrow(protect)) {
    l <- lab[grid_index(lat, lon, protect$plat, protect$plon)]; small[l[l > 0]] <- FALSE
  }
  drown <- lab > 0 & small[pmax(lab, 1)]
  out <- zc; out[drown] <- cfg$SEAWAY_DEPTH
  list(z = out, n = sum(small), area = sum(W[drown]))
}

#' new sea cells in small pools not connected to the ocean and holding no marine
#' record go back to low land (+1 m)
drop_stray_pools <- function(zc, z0, lat, lon, recs) {
  W <- area_matrix(lat, ncol(zc))
  lab <- label_wrap(zc < 0, 4)
  s <- tapply(as.vector(W), as.vector(lab), sum); ids <- as.integer(names(s))
  ocean <- ids[ids > 0][which.max(s[ids > 0])]
  keep <- c(ocean, lab[grid_index(lat, lon, recs$plat, recs$plon)])
  stray <- lab > 0 & !(lab %in% keep) & z0 >= 0
  stray <- matrix(stray, nrow(zc))
  out <- zc; out[stray] <- 1
  list(z = out, area = sum(W[stray]))
}

#' coastline rounding in the changed zone (binary smoothing of land/sea)
smooth_coastline <- function(zc, z0, lat, keep_sea, keep_land, sigma_km) {
  zone <- edt(abs(zc - z0) < 1) <= 6
  sm <- smooth_km((zc < 0) * 1, lat, sigma_km, nbands = 30)
  new_sea <- ifelse(zone, sm > 0.5, zc < 0)
  new_sea <- (new_sea | keep_sea) & !keep_land
  out <- zc
  to_sea <- new_sea & zc >= 0; to_land <- !new_sea & zc < 0
  out[to_sea] <- -10; out[to_land] <- 1
  list(z = out, to_sea = sum(to_sea), to_land = sum(to_land))
}

#' spits / barrier ridges: land that does not survive a morphological opening
#' (radius SPIT_HALF_KM) near new sea becomes shallow sea
remove_spits <- function(zc, z0, lat, keep_land, cfg) {
  land <- zc >= 0
  core <- land & dist_km(!land, lat) > cfg$SPIT_HALF_KM
  opened <- land & dist_km(core, lat) <= cfg$SPIT_HALF_KM
  near_new <- dist_km(zc < 0 & z0 >= 0, lat) <= cfg$SPIT_ZONE_KM
  thin <- land & !opened & near_new & !keep_land
  out <- zc; out[thin] <- -10
  list(z = out, n = sum(thin))
}

#' elevation ramp across the corrected coasts (no cell staircase in the 0 m line)
smooth_coast_elevation <- function(zc, z0, lat, keep_sea, keep_land, sigma_km) {
  sm <- smooth_km((zc < 0) * 1, lat, sigma_km, nbands = 30)
  sfield <- sm - 0.5
  coast <- edt(zc < 0) <= 2 & edt(zc >= 0) <= 2
  zone <- edt(abs(zc - z0) < 1) <= 4 & coast
  ramp <- -80 * sfield
  ramp <- ifelse(ramp >= 0, pmax(ramp, 0.5), pmin(ramp, -0.5))
  flip <- sign(ramp) != sign(zc)
  out <- ifelse(zone, ifelse(flip, ramp, ifelse(zc < 0, zc, pmin(zc, ramp))), zc)
  out <- ifelse(keep_sea, pmin(out, -5), out)
  out <- ifelse(keep_land, pmax(out, 5), out)
  out <- matrix(out, nrow(zc))
  list(z = out, flips = sum(sign(out) != sign(zc)))
}

#' cells within frac x designed half-width of any seaway route
seaway_core_mask <- function(dims, lat, paths, frac = 1) {
  n <- dims[1]; p <- dims[2]; dl <- abs(lat[2] - lat[1]); mask <- matrix(FALSE, n, p)
  for (P in paths) {
    jc <- .unwrap(P[, 2], p); ii <- round(P[, 1])
    r <- ceiling(frac * max(P[, 3]) / (dl * KM_PER_DEG * 0.3)) + 2
    I <- max(min(ii) - r, 1):min(max(ii) + r, n); J <- (floor(min(jc)) - r):(ceiling(max(jc)) + r)
    GI <- rep(I, times = length(J)); GJ <- rep(J, each = length(I))
    nn <- nn2(.to_km(P[, 1], jc, lat), .to_km(GI, GJ, lat), k = 1)
    inside <- nn$nn.dists[, 1] <= frac * P[nn$nn.idx[, 1], 3]
    mask[cbind(GI[inside], .wrapj(GJ[inside], p))] <- TRUE
  }
  mask
}

#' final step: smooth 0 m line everywhere (true-km gaussian of land/sea; stronger
#' near the corrections), seaway cores protected with a soft score, elevation ramp
final_coast_smooth <- function(zc, lat, keep_sea, keep_land, cfg, zone_dist = NULL, soft_sea = NULL, soft_land = NULL) {
  wz <- if (is.null(zone_dist)) 0 else pmin(pmax(1 - zone_dist / cfg$FINAL_COAST_ZONE_KM, 0), 1)
  blend <- function(a) (1 - wz) * smooth_km(a, lat, cfg$FINAL_COAST_SIGMA_KM) +
    wz * smooth_km(a, lat, cfg$FINAL_COAST_SIGMA_CORR_KM)
  # land made from land points is left out here; its own smooth outline is added below
  sea_ind <- if (!is.null(soft_land) && any(soft_land)) (zc < 0 | soft_land) * 1 else (zc < 0) * 1
  sfield <- blend(sea_ind) - 0.5
  if (!is.null(soft_sea)) {
    c_ <- blend(soft_sea * 1)
    score <- ifelse(c_ < 0.3, (c_ - 0.3) * (0.5 / 0.3), (c_ - 0.3) * (0.5 / 0.7))
    sfield <- pmax(sfield, score)
  }
  if (!is.null(soft_land) && any(soft_land)) {
    # land from land points: its own (already smooth) outline, lightly smoothed, wins over the sea
    # several grid cells wide, so the 0 m line is interpolated between cells (no staircase)
    sl <- smooth_km(soft_land * 1, lat, max(cfg$FINAL_COAST_SIGMA_KM, 0.4 * cfg$LAND_POINT_RADIUS_KM),
                    nbands = 45) - 0.5
    sfield <- pmin(sfield, -sl)
  }
  want_sea <- (sfield > 0 | keep_sea) & !keep_land
  ramp <- -80 * sfield
  ramp <- ifelse(want_sea, pmin(ramp, -0.5), pmax(ramp, 0.5))
  flip <- want_sea != (zc < 0)
  band <- abs(sfield) < 0.45 | flip
  out <- zc
  out[flip] <- ramp[flip]
  keep <- band & !flip
  out[keep] <- pmin(zc[keep], ramp[keep])
  if (!is.null(soft_land) && any(soft_land)) {
    # around land made from land points (often in deep water) the sea side also follows the
    # ramp: a shallow apron, so the 0 m line is interpolated smoothly instead of cell by cell
    # continuous sea-floor profile around that land: the gentle ramp 80 m x sl across the
    # smoothed outline (|sl| < 0.45), then a shelf at MAX_SEA_SLOPE m/km down to -1400 m and
    # MAX_DEEP_SLOPE below - the sea floor is only raised, never lowered
    band_l <- sl > -0.45
    d_b <- dist_km(band_l, lat)
    SEA_RAMP <- 400                                  # m per unit sl on the sea side (~ -180 m at the band edge)
    d1 <- (1400 - 0.45 * SEA_RAMP) / cfg$MAX_SEA_SLOPE
    shelf <- ifelse(band_l, SEA_RAMP * sl,
             ifelse(d_b <= d1, -0.45 * SEA_RAMP - cfg$MAX_SEA_SLOPE * d_b, -1400 - cfg$MAX_DEEP_SLOPE * (d_b - d1)))
    shelf <- smooth_km(pmax(shelf, -4500), lat, 30, nbands = 45)   # no cell-by-cell steps offshore
    shelf[band_l] <- pmax(SEA_RAMP * sl[band_l], shelf[band_l])   # exact ramp at the coast itself
    apron <- want_sea & d_b < d1 + 3000 / cfg$MAX_DEEP_SLOPE
    out[apron] <- pmax(out[apron], pmin(shelf[apron], -0.5))
    sea_b <- want_sea & band_l                                 # at the coast: exactly the ramp, both sides
    out[sea_b] <- pmin(pmax(SEA_RAMP * sl[sea_b], shelf[sea_b]), -0.5)
    land_b <- !want_sea & band_l & sl < 0.45
    out[land_b] <- pmax(pmin(SEA_RAMP * sl[land_b], 60), 0.5)   # same slope both sides: 0 m line exactly at sl = 0
  }
  out[keep_sea] <- pmin(out[keep_sea], -5)
  out[keep_land] <- pmax(out[keep_land], 5)
  list(z = out, flips = sum(flip))
}

# --------------------------------------------------------- deeper marine ----
#' weight field (0-1) around 'marine_deep' records: 1 within MARINE_DEEP_RADIUS_KM, then a
#' smooth S-shaped (smoothstep) taper to 0 over MARINE_DEEP_TAPER_KM
deep_weight <- function(deep, lat, lon, cfg) {
  w <- matrix(0, length(lat), length(lon))
  if (!nrow(deep)) return(w)
  tg <- matrix(FALSE, length(lat), length(lon)); tg[grid_index(lat, lon, deep$plat, deep$plon)] <- TRUE
  R <- cfg$MARINE_DEEP_RADIUS_KM; Tk <- cfg$MARINE_DEEP_TAPER_KM
  d <- dist_km(tg, lat, max_km = R + Tk + 50)
  u <- pmin(pmax((d - R) / Tk, 0), 1)
  1 - u * u * (3 - 2 * u)
}

#' land that is land again in the result (flooded on the way by the PaleoCoastlines,
#' then made land again by land points, their bridges or the coast smoothing) gets back
#' the original PaleoDEM relief where the PaleoDEM had land: min(z0, bank ramp from the
#' final coast), only raising. Without this such land stays a +1...70 m lowland strip, e.g.
#' a dip between land points and a mountain range. Only around the land points
#' (full within LAND_LINK_KM of their land, fading out by 1.5 x that distance).
restore_relief <- function(zc, z0, lat, cfg, near_of) {
  land <- zc >= 0
  L <- max(cfg$LAND_LINK_KM, 2 * cfg$LAND_POINT_RADIUS_KM)
  w <- pmin(pmax(1 - (dist_km(near_of, lat, max_km = 1.6 * L) - L) / (0.5 * L), 0), 1)
  sel <- land & z0 > zc + 1 & w > 0
  if (!any(sel)) return(list(z = zc, area_Mkm2 = 0))
  ramp <- 1 + cfg$MAX_BANK_SLOPE * dist_km(!land, lat, max_km = 1500)   # no cliffs at the final coast
  tgt <- pmin(z0, ramp)
  sel <- sel & tgt > zc
  zc[sel] <- zc[sel] + w[sel] * (tgt[sel] - zc[sel])
  W <- area_matrix(lat, ncol(zc))
  list(z = zc, area_Mkm2 = round(sum(W[sel]), 3))
}

# ------------------------------------------------------------- land points ----
.hav_km <- function(la1, lo1, la2, lo2) {
  a <- sin((la2 - la1) * pi / 360)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin((lo2 - lo1) * pi / 360)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
#' points along the great circle between two points (every ~step_km)
.gc_line <- function(la1, lo1, la2, lo2, step_km = 10) {
  d <- .hav_km(la1, lo1, la2, lo2); k <- max(2, ceiling(d / step_km) + 1)
  v1 <- c(cos(la1 * pi / 180) * cos(lo1 * pi / 180), cos(la1 * pi / 180) * sin(lo1 * pi / 180), sin(la1 * pi / 180))
  v2 <- c(cos(la2 * pi / 180) * cos(lo2 * pi / 180), cos(la2 * pi / 180) * sin(lo2 * pi / 180), sin(la2 * pi / 180))
  f <- seq(0, 1, length.out = k); v <- outer(1 - f, v1) + outer(f, v2); v <- v / sqrt(rowSums(v^2))
  cbind(lat = asin(v[, 3]) * 180 / pi, lon = atan2(v[, 2], v[, 1]) * 180 / pi)
}

#' Land kept around land records, as connected, rounded land instead of separate
#' pixel islets:
#'  - land points closer together than LAND_LINK_KM are joined by land bridges
#'    (minimum spanning tree), and each group is joined to the nearest existing land
#'    within LAND_LINK_KM;
#'  - a bridge is not drawn if it would pass within LAND_MARINE_CLEAR_KM of a marine
#'    or evaporite record (the sea there is documented);
#'  - the land is everything within LAND_POINT_RADIUS_KM of points and bridges
#'    (rounded outline), minus cells closer to a marine record than to the land skeleton;
#'  - elevation: a low dome, LAND_POINT_ELEV at the coast rising inland.
land_point_mask <- function(hard, me, zref, lat, lon, cfg) {
  n <- length(lat); p <- length(lon)
  pts <- hard[, c("plat", "plon")]
  clear_ok <- function(L) {                        # bridge far enough from marine records?
    if (!nrow(me)) return(TRUE)
    for (k in seq_len(nrow(L))) {
      near <- abs(me$plat - L[k, 1]) < cfg$LAND_MARINE_CLEAR_KM / 100
      if (any(near) && min(.hav_km(L[k, 1], L[k, 2], me$plat[near], me$plon[near])) < cfg$LAND_MARINE_CLEAR_KM)
        return(FALSE)
    }
    TRUE
  }
  lines_ <- list()
  # minimum spanning tree (Prim) over the land points
  m <- nrow(pts)
  if (m > 1) {
    D <- outer(seq_len(m), seq_len(m), function(i, j) .hav_km(pts$plat[i], pts$plon[i], pts$plat[j], pts$plon[j]))
    intree <- c(TRUE, rep(FALSE, m - 1)); best <- D[1, ]; from <- rep(1, m)
    for (step in seq_len(m - 1)) {
      cand <- which(!intree); k <- cand[which.min(best[cand])]
      if (best[k] <= cfg$LAND_LINK_KM && best[k] > 0) {
        L <- .gc_line(pts$plat[from[k]], pts$plon[from[k]], pts$plat[k], pts$plon[k])
        if (clear_ok(L)) lines_[[length(lines_) + 1]] <- L
      }
      intree[k] <- TRUE
      upd <- !intree & D[k, ] < best; best[upd] <- D[k, upd]; from[upd] <- k
    }
  }
  # each point: bridge to the nearest existing land within LAND_LINK_KM
  land <- zref >= 0
  for (k in seq_len(m)) {
    dl <- abs(lat[2] - lat[1]); r_deg <- cfg$LAND_LINK_KM / 111.2
    ii <- which(abs(lat - pts$plat[k]) <= r_deg)
    jj <- which(abs(((lon - pts$plon[k] + 180) %% 360) - 180) <= r_deg / max(cos(pts$plat[k] * pi / 180), 0.05))
    sub <- land[ii, jj, drop = FALSE]
    if (!any(sub)) next
    w <- which(sub, arr.ind = TRUE)
    d <- .hav_km(pts$plat[k], pts$plon[k], lat[ii[w[, 1]]], lon[jj[w[, 2]]])
    if (min(d) < 1e-6 || min(d) > cfg$LAND_LINK_KM) next
    b <- which.min(d)
    L <- .gc_line(pts$plat[k], pts$plon[k], lat[ii[w[b, 1]]], lon[jj[w[b, 2]]])
    if (clear_ok(L)) lines_[[length(lines_) + 1]] <- L
  }
  skel <- matrix(FALSE, n, p)
  skel[grid_index(lat, lon, pts$plat, pts$plon)] <- TRUE
  for (L in lines_) skel[grid_index(lat, lon, L[, 1], L[, 2])] <- TRUE
  d_s <- dist_km(skel, lat)
  free <- matrix(TRUE, n, p)                        # cells nearer a marine record stay sea
  if (nrow(me)) {
    tm <- matrix(FALSE, n, p); tm[grid_index(lat, lon, me$plat, me$plon)] <- TRUE
    d_m <- dist_km(tm, lat)
    free <- d_m > pmin(d_s, cfg$LAND_MARINE_CLEAR_KM)     # round sea notch (<= 40 km) around marine records
  }
  Rk <- cfg$LAND_POINT_RADIUS_KM
  mask <- d_s <= Rk & free
  # morphological closing (radius Rk): fills the notches between neighbouring disks and
  # bridges, so a group of points becomes one compact landmass rather than a chain of lobes
  dil <- dist_km(mask, lat) <= Rk
  mask <- !(dist_km(!dil, lat) <= Rk) & free
  # smooth outline: gaussian of the land/sea indicator re-thresholded (true km, per latitude band)
  sig <- max(Rk / 2, 20)
  mask <- ((smooth_km(mask * 1, lat, sig, nbands = 45) > 0.5) | (skel & free)) & free
  d_in <- dist_km(!mask, lat)                      # km from the new coast, inland
  core <- (mask & d_in >= 0.5 * Rk) | (skel & mask)   # kept land even under the final coast smoothing
  # smooth dome: LAND_POINT_ELEV at the coast rising ~60 m over one radius inland
  elev <- cfg$LAND_POINT_ELEV + 60 * pmin(1, d_in / Rk)
  list(mask = mask, core = core, elev = elev, n_bridges = length(lines_))
}

# ---------------------------------------------------------------- seaways ----
#' sea areas that hold a marine / evaporite record but are not connected to the
#' ocean get a seaway: least-cost route over the land surface (cost 1 + z/H_COST,
#' with meander noise) to the ocean, carved as a funnel that never narrows
#' towards the mouth. Appends the routes to st$paths.
connect_to_ocean <- function(zc, z0, lat, lon, recs, cfg, st) {
  n <- nrow(zc); p <- ncol(zc); dl <- abs(lat[2] - lat[1])
  lab <- label_wrap(zc < 0, 4)
  W <- area_matrix(lat, p)
  s <- tapply(as.vector(W), as.vector(lab), sum); ids <- as.integer(names(s))
  sizes <- numeric(max(lab)); sizes[ids[ids > 0]] <- s[ids > 0]
  ocean <- which.max(sizes)
  rl <- unique(lab[grid_index(lat, lon, recs$plat, recs$plon)])
  targets <- sort(rl[rl != 0 & rl != ocean])
  log <- list()
  if (!length(targets)) return(list(z = zc, log = log))
  pad <- p %/% 4
  wrapc <- function(a) cbind(a[, (p - pad + 1):p, drop = FALSE], a, a[, 1:pad, drop = FALSE])
  zp <- wrapc(zc); labp <- wrapc(lab)
  coslat <- .cos_clip(lat)
  if (is.null(st$noise)) {
    lonp <- c(lon[(p - pad + 1):p] - 360, lon, lon[1:pad] + 360)
    st$noise <- fractal_noise(lat, lonp, cfg$NOISE_SEED)
    st$noise2 <- fractal_noise(lat, lonp, cfg$NOISE_SEED + 1)
  }
  cost <- (1 + pmax(zp, 0) / cfg$H_COST) * (0.5 + 0.5 * coslat) * exp(cfg$MEANDER * st$noise)
  if (!is.null(st$hard) && any(st$hard)) cost[wrapc(st$hard)] <- cost[wrapc(st$hard)] * 1000   # land records: avoid
  cost[labp == ocean] <- 0
  acc <- cost_distance(cost)
  P2 <- ncol(zp)
  depress <- matrix(0, n, P2); core <- matrix(FALSE, n, P2); coredepth <- matrix(Inf, n, P2)
  for (t in targets) {
    cells <- which(labp == t)
    best <- cells[which.min(acc[cells])]
    bi <- (best - 1) %% n + 1; bj <- (best - 1) %/% n + 1
    path <- trace_path(acc, cost, bi, bj)
    path <- path[nrow(path):1, , drop = FALSE] * 1.0          # ocean end first
    if (nrow(path) > 7) {
      sm <- cbind(gauss1d(path[, 1], 3), gauss1d(path[, 2], 3))
      sm[1, ] <- path[1, ]; sm[nrow(sm), ] <- path[nrow(path), ]
      path <- sm
    }
    if (nrow(path) < 2) path <- rbind(path, path)
    seg <- sqrt(rowSums(diff(path)^2)); nseg <- pmax(ceiling(seg * 3), 1)
    dense <- do.call(rbind, lapply(seq_len(nrow(path) - 1), function(k) {
      f <- (seq_len(nseg[k]) - 1) / nseg[k]
      cbind(path[k, 1] + f * (path[k + 1, 1] - path[k, 1]), path[k, 2] + f * (path[k + 1, 2] - path[k, 2]))
    }))
    dense <- rbind(dense, path[nrow(path), ])
    ii <- pmin(pmax(round(dense[, 1]), 1), n); jj <- pmin(pmax(round(dense[, 2]), 1), P2)
    zline <- zp[cbind(ii, jj)]
    zpath <- gauss1d(pmax(zline, 0), 15)
    onland <- zline >= 0
    stepkm <- sqrt((diff(dense[, 1]) * dl)^2 + (diff(dense[, 2]) * dl * coslat[ii[-1]])^2) * KM_PER_DEG
    len_km <- sum(stepkm[onland[-1]])
    log[[length(log) + 1]] <- data.frame(water_body = t, area_Mkm2 = round(sizes[t], 3),
      seaway_length_km = round(len_km), max_elev_crossed_m = round(if (any(onland)) max(zline[onland]) else 0, 1),
      lat = round(lat[bi], 2))
    # funnel geometry: dense[1, ] at the ocean, last row inland
    Pk <- .to_km(dense[, 1], dense[, 2], lat)
    s_km <- c(0, cumsum(sqrt(rowSums(diff(Pk)^2))))
    L <- max(s_km[length(s_km)], 1); frac <- 1 - s_km / L
    w_raw <- cfg$SEAWAY_HALF_KM_INLAND * pmin(pmax(1 + 0.25 * st$noise2[cbind(ii, jj)], 0.7), 1.4)
    w_mono <- rev(cummax(rev(w_raw)))
    w_mouth <- min(max(L / (2 * cfg$SEAWAY_ASPECT), cfg$SEAWAY_MOUTH_MIN * cfg$SEAWAY_HALF_KM_INLAND), cfg$SEAWAY_HALF_KM_MAX)
    w <- w_mono * (1 + (w_mouth / cfg$SEAWAY_HALF_KM_INLAND - 1) * frac^cfg$SEAWAY_SHAPE)
    depth <- cfg$SEAWAY_DEPTH - 30 * frac
    reach <- max(w) + 3 * cfg$SEAWAY_BANK_KM
    rc <- ceiling(reach / (dl * KM_PER_DEG * 0.3))
    I <- max(min(ii) - rc, 1):min(max(ii) + rc, n); J <- max(min(jj) - rc, 1):min(max(jj) + rc, P2)
    GI <- rep(I, times = length(J)); GJ <- rep(J, each = length(I))
    nn <- nn2(Pk, .to_km(GI, GJ, lat), k = 1, searchtype = "radius", radius = reach)
    idx <- nn$nn.idx[, 1]; dkm <- nn$nn.dists[, 1]; ok <- idx > 0
    st$paths[[length(st$paths) + 1]] <- cbind(dense[, 1], .wrapj(dense[, 2] - pad, p), w)
    incore <- ok; incore[ok] <- dkm[ok] <= w[idx[ok]]
    od <- ok & !incore
    D <- numeric(length(GI))
    D[od] <- (zpath[idx[od]] + 40) * exp(-((dkm[od] - w[idx[od]]) / cfg$SEAWAY_BANK_KM)^2)
    cellsW <- cbind(GI, GJ)
    depress[cellsW] <- pmax(depress[cellsW], D)
    cd <- rep(Inf, length(GI)); cd[incore] <- depth[idx[incore]]
    coredepth[cellsW] <- pmin(coredepth[cellsW], cd)
    core[cellsW] <- core[cellsW] | incore
  }
  fold <- function(a, f) {
    out <- a[, (pad + 1):(pad + p)]
    out[, (p - pad + 1):p] <- f(out[, (p - pad + 1):p], a[, 1:pad])
    out[, 1:pad] <- f(out[, 1:pad], a[, (pad + p + 1):(2 * pad + p)])
    out
  }
  dep <- fold(depress, pmax); cor <- fold(core, `|`); cdep <- fold(coredepth, pmin)
  land <- zc >= 0
  new <- zc; new[land] <- pmax(zc[land] - dep[land], 1)
  new[cor] <- pmin(new[cor], cdep[cor])
  list(z = new, log = log)
}

#' make every newly cut channel at least as wide downstream as anywhere upstream
enforce_widening <- function(zc, lat, paths, zref, cfg, tol_km = cfg$WIDEN_TOL_KM, rounds = 3) {
  n <- nrow(zc); p <- ncol(zc); z <- zc; widened <- 0L
  for (round_ in seq_len(rounds)) {
    changed <- FALSE
    for (P in paths) {
      ii <- pmin(pmax(round(P[, 1]), 1), n); jj <- .wrapj(round(P[, 2]), p)
      new_part <- zref[cbind(ii, jj)] >= 0
      if (sum(new_part) < 2) next
      jc <- .unwrap(P[, 2], p); r <- 60
      I <- max(min(ii) - r, 1):min(max(ii) + r, n); J <- (floor(min(jc)) - r):(ceiling(max(jc)) + r)
      GI <- rep(I, times = length(J)); GJ <- rep(J, each = length(I)); GJw <- .wrapj(GJ, p)
      G <- .to_km(GI, GJ, lat)
      landw <- z[cbind(GI, GJw)] >= 0
      if (!any(landw)) next
      Q <- .to_km(P[, 1], jc, lat)
      k_new <- rev(which(new_part))                       # inland -> mouth
      hw <- nn2(G[landw, , drop = FALSE], Q[k_new, , drop = FALSE], k = 1)$nn.dists[, 1]
      hw <- pmin(hw, max(P[, 3]))
      run <- cummax(hw)
      for (q in which(hw < run - tol_km)) {
        k <- k_new[q]; rr <- run[q]
        near <- (G[, 1] - Q[k, 1])^2 + (G[, 2] - Q[k, 2])^2 <= rr^2
        cells <- cbind(GI[near], GJw[near])
        sel <- z[cells] >= 0
        if (any(sel)) { z[cells[sel, , drop = FALSE]] <- cfg$SEAWAY_DEPTH; changed <- TRUE; widened <- widened + sum(sel) }
      }
    }
    if (!changed) break
  }
  list(z = z, widened = widened)
}

#' measured half-width along every seaway and the largest narrowing towards the mouth
width_check <- function(zc, lat, paths, zref) {
  n <- nrow(zc); p <- ncol(zc); sea <- zc < 0
  do.call(rbind, lapply(paths, function(P) {
    ii <- pmin(pmax(round(P[, 1]), 1), n); jj <- .wrapj(round(P[, 2]), p)
    new_part <- zref[cbind(ii, jj)] >= 0
    allm <- all(sea[cbind(ii, jj)])
    if (sum(new_part) < 2) return(data.frame(max_narrowing_km = 0, route_all_marine = allm,
                                             halfwidth_inland_km = NA, halfwidth_mouth_km = NA))
    jc <- .unwrap(P[, 2], p); r <- 40
    I <- max(min(ii) - r, 1):min(max(ii) + r, n); J <- (floor(min(jc)) - r):(ceiling(max(jc)) + r)
    GI <- rep(I, times = length(J)); GJ <- rep(J, each = length(I))
    landw <- !sea[cbind(GI, .wrapj(GJ, p))]
    hw <- if (!any(landw)) rep(Inf, sum(new_part)) else          # open water all around
      nn2(.to_km(GI[landw], GJ[landw], lat), .to_km(P[new_part, 1], jc[new_part], lat), k = 1)$nn.dists[, 1]
    wloc <- rev(P[new_part, 3]); itm <- rev(hw)
    run <- cummax(pmin(itm, wloc))
    narrowing <- pmax(pmin(run, wloc) - itm, 0)
    data.frame(max_narrowing_km = round(max(narrowing), 1), route_all_marine = allm,
               halfwidth_inland_km = round(mean(head(itm, 5)), 1), halfwidth_mouth_km = round(mean(tail(itm, 5)), 1))
  }))
}

# ------------------------------------------------------------------- main ----
#' corrected PaleoDEM for one map age.
#' recs: records_for_age(age, purpose = "dem") (columns plat, plon, setting)
correct_dem <- function(age, recs, cfg = SILUCARB, write = TRUE, verbose = TRUE) {
  t_start <- Sys.time()
  say <- function(...) if (verbose) message(sprintf("[%g Ma %5.0fs] ", age, as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ...)
  dem <- read_dem_nc(find_dem(cfg$dir_dem, age))
  z <- dem$z; lat <- dem$lat; lon <- dem$lon; z0 <- z; n <- nrow(z); p <- ncol(z); dl <- abs(lat[2] - lat[1])
  W <- area_matrix(lat, p)
  st <- new.env(); st$paths <- list()
  mar <- recs[recs$setting %in% c("marine", "marine_deep"), ]; eva <- recs[recs$setting == "evaporite", ]
  deep <- recs[recs$setting == "marine_deep", ]
  ter <- recs[recs$setting == "terrestrial", ]; me <- rbind(mar, eva)
  misfit <- function(zz) c(
    if (nrow(mar)) 100 * mean(sample_grid(zz, lat, lon, mar$plat, mar$plon) >= 0) else NA,
    if (nrow(eva)) 100 * mean(sample_grid(zz, lat, lon, eva$plat, eva$plon) > 0) else NA,
    if (nrow(ter)) 100 * mean(sample_grid(zz, lat, lon, ter$plat, ter$plon) < 0) else NA)
  before <- misfit(z)

  # ---- tier 1: PaleoCoastlines land/sea
  land_cs <- cs_land_mask(age, lat, lon, cfg)
  if (is.null(land_cs)) { land_cs <- z >= 0; say("no PaleoCoastlines for this age: land/sea from the PaleoDEM") }
  flood <- z >= 0 & !land_cs; raise_ <- z < 0 & land_cs & z > -200
  target <- z
  target[flood] <- pmin(cfg$SEA_DEPTH, -0.15 * z[flood] - 5)
  target[raise_] <- cfg$LAND_NEW

  # ---- tier 2: local flooding at remaining marine / evaporite records
  carve <- matrix(0, n, p); cdepth <- matrix(0, n, p)
  nr <- ceiling(2 * cfg$R_POINT / dl) + 1
  for (set in list(list(d = mar, depth = cfg$SEA_DEPTH), list(d = eva, depth = cfg$EVAP_DEPTH))) {
    if (!nrow(set$d)) next
    zz <- sample_grid(target, lat, lon, set$d$plat, set$d$plon)
    for (k in seq_len(nrow(set$d))) {
      if ((set$depth == cfg$SEA_DEPTH && zz[k] < 0) || (set$depth == cfg$EVAP_DEPTH && zz[k] <= 0)) next
      la <- set$d$plat[k]; lo <- set$d$plon[k]
      i0 <- round((la - lat[1]) / dl) + 1; j0 <- round((lo - lon[1]) / dl)
      ii <- (i0 - nr):(i0 + nr); ii <- ii[ii >= 1 & ii <= n]; jj <- ((j0 - nr):(j0 + nr)) %% p + 1
      sla <- matrix(lat[ii], length(ii), length(jj)); slo <- matrix(lon[jj], length(ii), length(jj), byrow = TRUE)
      dlon <- (slo - lo + 180) %% 360 - 180
      dist <- sqrt((sla - la)^2 + (dlon * cos(pi / 180 * la))^2)
      wgt <- ifelse(dist <= cfg$R_POINT, 1, ifelse(dist <= 2 * cfg$R_POINT,
                    0.5 * (1 + cos(pi * (dist - cfg$R_POINT) / cfg$R_POINT)), 0))
      blk <- carve[ii, jj, drop = FALSE]; dep <- cdepth[ii, jj, drop = FALSE]
      upd <- wgt > blk; dep[upd] <- set$depth; blk[upd] <- wgt[upd]
      carve[ii, jj] <- blk; cdepth[ii, jj] <- dep
    }
  }
  corem <- carve >= 1
  target <- ifelse(carve > 0, (1 - carve) * target + carve * pmin(target, cdepth), target)
  sel <- corem & target >= 0; target[sel] <- cdepth[sel]
  if (nrow(ter)) {                         # protect terrestrial records outside marine cores
    gi <- grid_index(lat, lon, ter$plat, ter$plon)
    sel <- target[gi] < 0 & !corem[gi]; target[gi[sel, , drop = FALSE]] <- 5
  }
  # ---- land records: a small island / land patch that every later step keeps
  #      (your terrestrial points always; database terrestrial records when no marine
  #      record lies within 1 deg). Cells closer to a marine record stay free for the sea.
  ter_clean <- ter
  if (nrow(ter) && nrow(mar)) {
    dmin <- vapply(seq_len(nrow(ter)), function(k) min(sqrt((mar$plat - ter$plat[k])^2 +
      (((mar$plon - ter$plon[k] + 180) %% 360 - 180) * cos(pi / 180 * ter$plat[k]))^2)), 0)
    ter_clean <- ter[dmin > 1, ]
  }
  hard <- rbind(ter[ter$source == "user", ], if (isTRUE(cfg$HARD_LAND_DATABASE)) ter_clean[ter_clean$source != "user", ])
  hard_mask <- matrix(FALSE, n, p); hard_elev <- matrix(-Inf, n, p); hard_core <- hard_mask
  if (nrow(hard)) {
    lm <- land_point_mask(hard, me, target, lat, lon, cfg)
    hard_mask <- lm$mask; hard_elev <- lm$elev; hard_core <- lm$core
  }
  st$hard <- hard_mask
  # ---- deeper-shelf marine records: lower the sea floor to MARINE_DEEP_DEPTH around them
  #      (cosine taper to the surroundings; only lowers sea floor, never turns land into deep sea)
  deep_w <- deep_weight(deep, lat, lon, cfg)
  D <- function(zz) {
    if (!any(deep_w > 0)) return(zz)
    sea <- zz < 0
    # how much to lower each sea cell (never raise, never touch land) ...
    delta <- ifelse(sea & deep_w > 0, pmin(0, deep_w * (cfg$MARINE_DEEP_DEPTH - zz)), 0)
    # ... but not steeper than a normal shelf next to a nearby coast, and leave the
    # coastal ramp itself alone (keeps the 0 m line smooth)
    d_land <- dist_km(!sea, lat)
    shelf_floor <- -1 - cfg$MAX_SEA_SLOPE * d_land
    delta <- pmax(delta, pmin(0, shelf_floor - zz)) * pmin(pmax((d_land - 20) / 60, 0), 1)
    # smooth the lowering itself, so the deeper area blends into the untouched shelf
    delta <- pmin(smooth_km(delta, lat, cfg$MARINE_DEEP_SMOOTH_KM, nbands = 45), 0)
    zz[sea] <- pmin(zz[sea] + delta[sea], -0.5)
    zz
  }
  target <- D(target)
  H <- function(zz) { if (any(hard_mask)) zz[hard_mask] <- pmax(zz[hard_mask], hard_elev[hard_mask]); zz }
  # after the final coast smoothing only the inner part is forced (the coast itself follows the smooth ramp)
  Hc <- function(zz) { if (any(hard_core)) zz[hard_core] <- pmax(zz[hard_core], hard_elev[hard_core]); zz }
  target <- H(target)
  changed <- abs(target - z0) > 1e-6
  sm <- smooth_cells(target, 1)
  near <- edt(!changed) <= 2
  zc <- ifelse(near & !changed, sm, target)
  zc <- ifelse(!changed & sign(zc) != sign(z0), z0, zc)
  zc <- matrix(zc, n)
  say("tiers 1-2 done")

  # ---- tier 3: seaways
  z_before_seaways <- zc
  r <- connect_to_ocean(zc, z0, lat, lon, me, cfg, st); zc <- H(r$z); seaways <- r$log
  chg <- edt(abs(zc - z0) < 1) <= 3
  sl_before <- max_slopes(zc, lat, chg)
  protect_pts <- rbind(ter_clean, hard[!hard$record_id %in% ter_clean$record_id, ])
  n_isl <- 0L; a_isl <- 0; stray_area <- 0
  isl <- function(zz) { r <- drown_islets(zz, z0, lat, lon, cfg, protect_pts, z_before_seaways)
                        n_isl <<- n_isl + r$n; a_isl <<- a_isl + r$area; H(r$z) }
  zc <- isl(zc)
  r <- relax_banks(zc, z0, lat, cfg); zc <- H(r$z); relax <- r$stats
  zc <- isl(zc)
  r <- drop_stray_pools(zc, z0, lat, lon, me); zc <- H(r$z); stray_area <- stray_area + r$area
  for (round_ in 1:3) {                    # repair lost connections
    r <- connect_to_ocean(zc, z0, lat, lon, me, cfg, st); zc <- H(r$z)
    if (!length(r$log)) break
    seaways <- c(seaways, r$log)
    zc <- isl(zc)
    zc <- H(relax_banks(zc, z0, lat, cfg)$z)
    r <- drop_stray_pools(zc, z0, lat, lon, me); zc <- H(r$z); stray_area <- stray_area + r$area
  }
  sl_orig <- max_slopes(z0, lat, matrix(TRUE, n, p))
  paths <- st$paths
  say(length(paths), " seaways")

  # ---- geometry rules: widening, coast smoothing, spits, slopes
  r <- enforce_widening(zc, lat, paths, z_before_seaways, cfg); zc <- H(r$z); widened <- r$widened
  zc <- H(relax_banks(zc, z0, lat, cfg)$z)
  zc <- H(limit_land_slopes(zc, z0, lat, cfg$MAX_LAND_SLOPE_FINAL))
  zc <- isl(zc)
  keep_sea <- matrix(FALSE, n, p); keep_land <- matrix(FALSE, n, p)
  keep_sea[grid_index(lat, lon, me$plat, me$plon)] <- TRUE
  for (P in paths) {
    pi_ <- pmin(pmax(round(P[, 1]), 1), n); pj_ <- .wrapj(round(P[, 2]), p)
    keep_sea[cbind(pi_, pj_)] <- TRUE
    keep_sea[cbind(pmin(pi_ + 1, n), pj_)] <- TRUE; keep_sea[cbind(pi_, .wrapj(pj_ + 1, p))] <- TRUE
  }
  if (nrow(ter_clean)) keep_land[grid_index(lat, lon, ter_clean$plat, ter_clean$plon)] <- TRUE
  keep_land_db <- keep_land
  keep_land <- keep_land | hard_mask
  r <- smooth_coastline(zc, z0, lat, keep_sea, keep_land, cfg$COAST_SMOOTH_KM); zc <- H(r$z)
  sm_to_sea <- r$to_sea; sm_to_land <- r$to_land
  r <- enforce_widening(zc, lat, paths, z_before_seaways, cfg); zc <- H(r$z); widened <- widened + r$widened
  r <- remove_spits(zc, z0, lat, keep_land, cfg); zc <- H(r$z); spit_cells <- r$n
  r <- enforce_widening(zc, lat, paths, z_before_seaways, cfg); zc <- H(r$z); widened <- widened + r$widened
  zc <- H(relax_banks(zc, z0, lat, cfg)$z)
  zc <- H(limit_land_slopes(zc, z0, lat, cfg$MAX_LAND_SLOPE_FINAL))
  zc <- isl(zc)
  r <- drop_stray_pools(zc, z0, lat, lon, me); zc <- H(r$z); stray_area <- stray_area + r$area
  r <- smooth_coast_elevation(zc, z0, lat, keep_sea, keep_land, cfg$COAST_ELEV_SMOOTH_KM)
  zc <- H(r$z); coast_flips <- r$flips
  r <- enforce_widening(zc, lat, paths, z_before_seaways, cfg); zc <- H(r$z); widened <- widened + r$widened
  zc <- H(relax_banks(zc, z0, lat, cfg)$z)
  zc <- H(limit_land_slopes(zc, z0, lat, cfg$MAX_LAND_SLOPE_FINAL))
  r <- connect_to_ocean(zc, z0, lat, lon, me, cfg, st); zc <- H(r$z); seaways <- c(seaways, r$log)
  paths <- st$paths
  r <- drop_stray_pools(zc, z0, lat, lon, me); zc <- H(r$z); stray_area <- stray_area + r$area
  zc <- isl(zc)
  say("geometry rules done")

  # ---- final: smooth 0 m line everywhere, then sign-preserving slope limits
  core_mask <- seaway_core_mask(dim(zc), lat, paths)
  zone_dist <- dist_km(abs(zc - z0) >= 1, lat)
  final_flips <- 0L
  for (round_ in 1:3) {
    r <- final_coast_smooth(zc, lat, keep_sea, (keep_land_db & !core_mask) | hard_core, cfg, zone_dist, core_mask,
                            soft_land = hard_mask)
    zc <- Hc(r$z); final_flips <- final_flips + r$flips
    st$paths <- list()
    r <- connect_to_ocean(zc, z0, lat, lon, me, cfg, st); zc <- Hc(r$z)
    if (!length(r$log)) break
    seaways <- c(seaways, r$log); paths <- c(paths, st$paths)
    core_mask <- seaway_core_mask(dim(zc), lat, paths)
  }
  zc <- Hc(relax_banks(zc, z0, lat, cfg)$z)
  zc <- Hc(limit_land_slopes(zc, z0, lat, cfg$MAX_LAND_SLOPE_FINAL))
  restored <- 0
  if (any(hard_mask)) { r <- restore_relief(zc, z0, lat, cfg, hard_mask); zc <- r$z; restored <- r$area_Mkm2 }
  # shelf edges of the corrected areas (flooded shelves next to deep ocean): smooth shelf +
  # continental-slope profile instead of a step, within 400 km of anything that was changed
  near_chg <- dist_km(abs(zc - z0) >= 1, lat)
  zc <- seafloor_profile(zc, lat, cfg, where = near_chg <= 400)
  # gentle blur of the sea floor in and around the corrected areas: flat flooded shelves,
  # carved patches and seaways blend into the original sea floor (no blocky steps);
  # sign-preserving (sea stays sea, land untouched)
  sea <- zc < 0; blur <- smooth_km(ifelse(sea, zc, 0), lat, cfg$SEAFLOOR_BLUR_KM, nbands = 45)
  wb <- pmin(pmax(1 - (near_chg - 100) / 200, 0), 1)          # full within 100 km, fading out by 300 km
  wb <- wb * pmin(pmax((dist_km(!sea, lat) - 30) / 50, 0), 1)  # keep the coastal ramp (smooth 0 m line) untouched
  sel <- sea & wb > 0
  zc[sel] <- pmin((1 - wb[sel]) * zc[sel] + wb[sel] * blur[sel], -0.5)
  zc <- D(zc)                                       # deeper-shelf records keep their depth
  st$paths <- list()
  still <- connect_to_ocean(zc, z0, lat, lon, me, cfg, st)$log      # dry run: should be empty
  chg <- edt(abs(zc - z0) < 1) <= 3
  sl_after <- max_slopes(zc, lat, chg)
  wc <- width_check(zc, lat, paths, z_before_seaways)
  after <- misfit(zc)
  sw <- if (length(seaways)) do.call(rbind, seaways) else NULL
  say("done: ", length(seaways), " seaways, ", length(still), " still isolated")

  stats <- data.frame(map_age = age, n_marine = nrow(mar), n_evaporite = nrow(eva), n_terrestrial = nrow(ter),
    n_user_points = sum(recs$source == "user"), n_land_points_kept = nrow(hard), n_marine_deep = nrow(deep),
    land_points_in_sea_after = if (nrow(hard)) sum(sample_grid(zc, lat, lon, hard$plat, hard$plon) < 0) else 0L,
    marine_on_land_before_pct = round(before[1], 1), marine_on_land_after_pct = round(after[1], 1),
    evaporite_above_sl_before_pct = round(before[2], 1), evaporite_above_sl_after_pct = round(after[2], 1),
    terrestrial_in_sea_before_pct = round(before[3], 1), terrestrial_in_sea_after_pct = round(after[3], 1),
    flooded_by_coastlines_Mkm2 = round(sum(W[flood]), 2), raised_by_coastlines_Mkm2 = round(sum(W[raise_]), 2),
    flooded_by_points_Mkm2 = round(sum(W[zc < 0 & z0 >= 0 & !flood]), 2),
    relax, relief_restored_Mkm2 = restored, stray_pools_removed_Mkm2 = round(stray_area, 3), islets_drowned = n_isl,
    islets_drowned_Mkm2 = round(a_isl, 3),
    slope_p99_original_all = round(sl_orig[1], 1), slope_max_changed_before_relax = round(sl_before[2], 1),
    slope_p99_changed_after = round(sl_after[1], 1), slope_max_changed_after = round(sl_after[2], 1),
    seaways_added = length(seaways), marine_water_bodies_still_isolated = length(still),
    cells_widened = widened, coast_smoothing_cells_to_sea = sm_to_sea, coast_smoothing_cells_to_land = sm_to_land,
    spit_cells_removed = spit_cells, coast_smoothing_sign_flips = coast_flips, final_smoothing_flips = final_flips,
    seaway_routes_fully_marine = if (is.null(wc)) 0 else sum(wc$route_all_marine),
    seaway_max_narrowing_km = if (is.null(wc)) 0 else max(wc$max_narrowing_km),
    seaway_total_km = if (is.null(sw)) 0 else sum(sw$seaway_length_km),
    ocean_DEM_Mkm2 = round(sum(W[z0 < 0]), 2), ocean_corrected_Mkm2 = round(sum(W[zc < 0]), 2),
    shallow_corrected_Mkm2 = round(sum(W[zc < 0 & zc >= -1400]), 2),
    runtime_s = round(as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  if (!is.null(sw) && !is.null(wc) && nrow(sw) == nrow(wc)) sw <- cbind(map_age = age, sw, wc)
  else if (!is.null(sw)) sw <- cbind(map_age = age, sw)
  if (write) {
    dir.create(cfg$dir_corrected, showWarnings = FALSE, recursive = TRUE)
    write_dem_nc(dem, zc, file.path(cfg$dir_corrected, basename(dem$file)),
      history = paste0("Fossil-corrected (paleogeographic map editor, R): land/sea from PaleoCoastlines v7.2 CS; ",
                       "flooding at marine / evaporite records (PBDB, Macrostrat, user points, with user edits); ",
                       "least-cost funnel seaways to the ocean; slope limits; islets/spits removed; 0 m line smoothed. ",
                       format(Sys.time())))
  }
  list(z = zc, z0 = z0, lat = lat, lon = lon, stats = stats, seaways = sw, paths = paths,
       records = recs, z_before_seaways = z_before_seaways)
}
