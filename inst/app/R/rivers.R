# =============================================================================
#  Rivers and freshwater into the sea
#
#  Runoff (surface + subsurface, mm/day) is routed downhill over the map on a 1-degree grid to
#  the coast: depressions are filled ("priority flood", Barnes et al. 2014) so every land cell
#  drains to the sea, then each cell passes its water to the neighbour it was reached from.
#  The river water that reaches the sea at a mouth is spread over the sea within ~300 km
#  (Gaussian, sigma 150 km) - a simple stand-in for a river plume.
#
#  freshwater into the sea = P - E over the sea + river water        (mm/day; > 0 = freshening)
#  salinity source: where the HadCM3L surface salinity is lower / higher than the ocean mean,
#  whether the local freshwater flux explains it (rivers, rain, evaporation) or the water must
#  have been brought by currents / mixing / upwelling / sea ice (local flux small or opposite).
# =============================================================================

RIVER_SPREAD_KM <- 150          # sigma of the plume (~300 km wide)
RIVER_MIN_KM3 <- 50             # rivers drawn from this discharge (km3/yr)
RIVER_MOUTH_KM3 <- 100          # river mouths marked from this discharge (km3/yr)

# ---- routing ---------------------------------------------------------------------------
.ROUTE <- new.env()
#' drainage on the 1-degree grid of a map file: for every land cell the cell it drains to
#' (`to`, cell index; an ocean cell = the mouth) and the order to accumulate in (`order`)
river_routing <- function(file) {
  inf <- file.info(file); key <- paste(file, inf$mtime, inf$size)
  if (!is.null(.ROUTE[[key]])) return(.ROUTE[[key]])
  g <- geo_1deg_cached(file); nr <- 180; nc <- 360; N <- nr * nc
  land <- g$lf >= 0.5
  h <- ifelse(land, g$h / pmax(g$lf, 1e-6), 0)                     # mean height of the land part
  # neighbours (8), longitude wraps
  rr <- rep(1:nr, nc); cc <- rep(1:nc, each = nr)
  nb <- matrix(NA_integer_, N, 8); k <- 0
  for (dr in -1:1) for (dc in -1:1) if (dr != 0 || dc != 0) {
    k <- k + 1; r2 <- rr + dr; c2 <- ((cc + dc - 1) %% nc) + 1
    nb[, k] <- ifelse(r2 >= 1 & r2 <= nr, r2 + (c2 - 1) * nr, NA_integer_)
  }
  to <- integer(N); done <- !as.vector(land); hv <- as.vector(h)
  # binary heap of (priority, cell)
  hk <- numeric(N); hc <- integer(N); n <- 0L
  push <- function(p, c) {
    n <<- n + 1L; i <- n
    while (i > 1L) { j <- i %/% 2L; if (hk[j] <= p) break; hk[i] <<- hk[j]; hc[i] <<- hc[j]; i <- j }
    hk[i] <<- p; hc[i] <<- c
  }
  pop <- function() {
    top <- c(hk[1], hc[1]); p <- hk[n]; c <- hc[n]; n <<- n - 1L; i <- 1L
    while (TRUE) { j <- 2L * i; if (j > n) break
      if (j < n && hk[j + 1L] < hk[j]) j <- j + 1L
      if (hk[j] >= p) break; hk[i] <<- hk[j]; hc[i] <<- hc[j]; i <- j }
    if (n >= 1L) { hk[i] <<- p; hc[i] <<- c }
    top
  }
  # seeds: coastal land cells drain to their lowest... nearest ocean neighbour
  coast <- which(!done & apply(nb, 1, function(v) any(!is.na(v) & !as.vector(land)[v])))
  for (c in coast) {
    v <- nb[c, ]; v <- v[!is.na(v) & !as.vector(land)[v]]
    to[c] <- v[1]; done[c] <- TRUE; push(hv[c], c)
  }
  ord <- integer(sum(land)); k <- 0L
  while (n > 0L) {
    t <- pop(); p <- t[1]; c <- as.integer(t[2]); k <- k + 1L; ord[k] <- c
    for (v in nb[c, ]) if (!is.na(v) && !done[v]) {
      done[v] <- TRUE; to[v] <- c; push(max(hv[v], p), v)
    }
  }
  ord <- ord[seq_len(k)]
  r <- list(to = to, order = ord, land = land, lat = GRID1$lat, lon = GRID1$lon)
  .ROUTE[[key]] <- r
  r
}

#' river discharge from a runoff field (mm/day on the 1-degree grid, land): accumulated flow in
#' every land cell and the water reaching the sea per ocean cell (km3/yr)
river_discharge <- function(route, runoff) {
  ro <- runoff; ro[!is.finite(ro) | !route$land] <- 0; ro <- pmax(ro, 0)
  A <- area_matrix(GRID1$lat, 360) / sum(area_matrix(GRID1$lat, 360)) * 5.1007e8   # km2 per cell
  acc <- as.vector(ro * 365.25 * 1e-6 * A)                                    # km3/yr made in each cell
  mouth <- numeric(length(acc)); land <- as.vector(route$land); to <- route$to
  for (c in rev(route$order)) {
    d <- to[c]
    if (land[d]) acc[d] <- acc[d] + acc[c] else mouth[d] <- mouth[d] + acc[c]
  }
  acc[!land] <- 0
  list(acc = matrix(acc, 180, 360), mouth = matrix(mouth, 180, 360), area = A)
}

#' river water spread over the sea around the mouths (mm/day of fresh water, sea cells)
river_plume <- function(Q, land, spread_km = RIVER_SPREAD_KM) {
  sea <- !land; A <- Q$area
  fd <- Q$mouth / A * 1e6 / 365.25                                             # mm/day at the mouth cell
  S <- smooth_km(fd, GRID1$lat, spread_km); M <- smooth_km(sea * 1, GRID1$lat, spread_km)
  out <- ifelse(sea, S / pmax(M, 1e-3), 0)
  tot_in <- sum(Q$mouth); tot_out <- sum(out * A) * 365.25 / 1e6
  if (tot_out > 0) out <- out * tot_in / tot_out                               # conserve the water
  out[!sea] <- NA
  out
}

SALINITY_CLASSES <- data.frame(
  code = 1:6,
  label = c("fresher: river water", "fresher: rain (P > E)", "fresher: brought by currents / mixing",
            "near the ocean mean", "saltier: brought by currents / mixing", "saltier: evaporation (E > P)"),
  col = c("#08519c", "#6baed6", "#9e9ac8", "#eeeeee", "#f1b6da", "#b35806"), stringsAsFactors = FALSE)
SALINITY_DS <- 0.3              # psu: departure from the ocean mean taken as fresher / saltier
SALINITY_FLUX <- 100            # mm/yr: local freshwater flux taken as able to explain it

#' freshwater fields for a climate (climate_for_map() result) on a map file: adds
#' river_input, freshwater (mm/day), salinity_source (class 1-6) and rivers (for drawing)
freshwater_fields <- function(cl, file) {
  rt <- river_routing(file); land <- rt$land
  ro <- .to_1deg(cl$runoff); pme <- .to_1deg(cl$pme)
  miss <- land & !is.finite(ro); ro[miss] <- pmax(pme[miss], 0)              # land the model had as sea
  Q <- river_discharge(rt, ro)
  riv <- river_plume(Q, land)
  pme_sea <- pme; pme_sea[land] <- NA; pme_sea[!land & !is.finite(pme_sea)] <- 0
  fw <- pme_sea + riv
  cl$river_input <- .fld(riv); cl$freshwater <- .fld(fw)
  if (!is.null(cl$salinity)) {
    S <- .to_1deg(cl$salinity); S[land] <- NA
    W <- area_matrix(GRID1$lat, 360); ok <- is.finite(S)
    dS <- S - sum(W[ok] * S[ok]) / sum(W[ok])
    F <- fw * 365.25; share <- riv / pmax(riv + pmax(pme_sea, 0), 1e-9)
    cls <- matrix(4L, 180, 360)
    cls[dS < -SALINITY_DS & F > SALINITY_FLUX & share >= 0.5] <- 1L
    cls[dS < -SALINITY_DS & F > SALINITY_FLUX & share < 0.5] <- 2L
    cls[dS < -SALINITY_DS & F <= SALINITY_FLUX] <- 3L
    cls[dS > SALINITY_DS & F >= -SALINITY_FLUX] <- 5L
    cls[dS > SALINITY_DS & F < -SALINITY_FLUX] <- 6L
    cls[!ok] <- NA
    cl$salinity_source <- .fld(cls)
  }
  cl$rivers <- list(route = rt, Q = Q, total_km3 = sum(Q$mouth))
  cl
}

#' rivers (flow lines) and river mouths (circles, area ~ discharge) on a map
draw_rivers <- function(rv, cfg, qmin = RIVER_MIN_KM3, col = "#1f5fa8", qmouth = RIVER_MOUTH_KM3) {
  if (is.null(rv)) return(invisible())
  acc <- as.vector(rv$Q$acc); to <- rv$route$to; land <- as.vector(rv$route$land)
  sel <- which(land & acc >= qmin)
  if (length(sel)) {
    r <- (sel - 1) %% 180 + 1; c <- (sel - 1) %/% 180 + 1
    d <- to[sel]; r2 <- (d - 1) %% 180 + 1; c2 <- (d - 1) %/% 180 + 1
    lon1 <- GRID1$lon[c]; lon2 <- GRID1$lon[c2]; lon2 <- lon1 + ((lon2 - lon1 + 180) %% 360) - 180
    cls <- cut(acc[sel], c(qmin, 100, 500, 2000, Inf), labels = FALSE, include.lowest = TRUE)
    for (k in unique(cls)) { i <- which(cls == k)
      draw_ll_list(lapply(i, function(j) data.frame(lon = c(lon1[j], lon2[j]), lat = c(GRID1$lat[r[j]], GRID1$lat[r2[j]]))),
                   cfg, col = adjustcolor(col, 0.85), lwd = c(0.5, 0.9, 1.5, 2.3)[k], lend = "round") }
  }
  m <- which(as.vector(rv$Q$mouth) >= qmouth)
  if (length(m)) {
    r <- (m - 1) %% 180 + 1; c <- (m - 1) %/% 180 + 1; q <- as.vector(rv$Q$mouth)[m]
    xy <- sf_project("EPSG:4326", map_crs(cfg), cbind(cfg$central_lon + ((GRID1$lon[c] - cfg$central_lon + 180) %% 360) - 180, GRID1$lat[r]),
                     keep = TRUE, warn = FALSE)
    vis <- is.finite(xy[, 1]) & !.far(GRID1$lon[c], GRID1$lat[r], cfg)
    o <- order(-q)
    o <- o[vis[o]]
    symbols(xy[o, 1], xy[o, 2], circles = 0.004 * diff(view_xlim(cfg)) * sqrt(q[o] / 1000), inches = FALSE, add = TRUE,
            bg = adjustcolor("#9ecae1", 0.55), fg = adjustcolor("#08306b", 0.9), lwd = 0.6)
  }
  invisible()
}
river_caption <- function(rv, qmin = RIVER_MIN_KM3)
  sprintf("rivers: runoff routed downhill on the map (1°), lines ≥ %g km³/yr, circles at mouths ≥ %g km³/yr (area ~ discharge), %.0f km³/yr into the sea in total",
          qmin, RIVER_MOUTH_KM3, rv$total_km3)
