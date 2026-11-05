# =============================================================================
#  Grid tools for global lat-lon grids stored as matrices z[lat, lon]
#  (row 1 = southernmost latitude, columns = longitudes, periodic in longitude;
#  the duplicated 180 deg column of the PaleoDEMs is removed on reading).
#  Heavy lifting is done by compiled code: terra (distance, patches, costDist),
#  base R FFT and stats::filter, RANN (nearest neighbours).
# =============================================================================
library(terra); library(ncdf4); library(RANN)
terraOptions(progress = 0)

EARTH_MKM2 <- 4 * pi * 6371.0088^2 / 1e6
KM_PER_DEG <- 111.2

#' read a PaleoDEM netCDF -> list(z, lat, lon, file, dup)  (z: lat x lon)
read_dem_nc <- function(file) {
  nc <- nc_open(file); on.exit(nc_close(nc))
  z <- ncvar_get(nc, "z"); lat <- as.numeric(ncvar_get(nc, "lat")); lon <- as.numeric(ncvar_get(nc, "lon"))
  if (nrow(z) == length(lon) && ncol(z) == length(lat)) z <- t(z)
  if (lat[1] > lat[length(lat)]) { lat <- rev(lat); z <- z[nrow(z):1, , drop = FALSE] }
  dup <- abs(lon[length(lon)] - lon[1] - 360) < 1e-6
  if (dup) { z <- z[, -ncol(z), drop = FALSE]; lon <- lon[-length(lon)] }
  storage.mode(z) <- "double"
  list(z = z, lat = lat, lon = lon, file = file, dup = dup)
}

#' write a corrected grid: copy of the source file with z replaced
write_dem_nc <- function(dem, z, out_file, history = "") {
  file.copy(dem$file, out_file, overwrite = TRUE)
  zz <- if (dem$dup) cbind(z, z[, 1]) else z
  nc <- nc_open(out_file, write = TRUE); on.exit(nc_close(nc))
  latf <- as.numeric(ncvar_get(nc, "lat"))
  if (latf[1] > latf[length(latf)]) zz <- zz[nrow(zz):1, ]
  dims <- vapply(nc$var$z$dim, function(d) d$name, "")
  if (grepl("lon|x", dims[1])) zz <- t(zz)
  ncvar_put(nc, "z", zz)
  if (nzchar(history)) ncatt_put(nc, 0, "history", history)
  invisible(out_file)
}

#' all map ages that have a PaleoDEM (0-540 Ma every ~5 Myr for PALEOMAP)
available_ages <- function(cfg = SILUCARB) {
  f <- list.files(cfg$dir_dem, "Ma\\.nc$")
  a <- as.numeric(sub(".*_([0-9.]+)Ma\\.nc$", "\\1", f))
  sort(unique(a[is.finite(a)]), decreasing = TRUE)
}
#' stage / epoch label of a map age (config table for the Silurian, else from the file name)
stage_label <- function(age, cfg = SILUCARB) {
  k <- as.character(age)
  if (!is.null(cfg$stage) && k %in% names(cfg$stage)) return(unname(cfg$stage[k]))
  f <- list.files(cfg$dir_dem, paste0("_", age, "Ma\\.nc$"))
  if (!length(f)) return("")
  lab <- sub("^Map[0-9.]+_PALEOMAP_[0-9a-z]+_", "", sub(paste0("_", age, "Ma\\.nc$"), "", f[1]))
  lab <- gsub("_", " ", lab); lab <- gsub("([a-z])([A-Z])", "\\1 \\2", lab)
  paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab)))
}
find_dem <- function(dir, age) {
  f <- list.files(dir, paste0("_", age, "Ma\\.nc$"), full.names = TRUE)
  if (!length(f) && exists("SILUCARB") && identical(normalizePath(dir, mustWork = FALSE),
                                                  normalizePath(SILUCARB$dir_dem, mustWork = FALSE)) &&
      exists("download_dems")) {
    download_dems(SILUCARB)                            # original PaleoDEMs missing: fetch them once
    f <- list.files(dir, paste0("_", age, "Ma\\.nc$"), full.names = TRUE)
  }
  if (!length(f)) stop("no PaleoDEM for ", age, " Ma in ", dir,
                       if (grepl("corrected", dir)) " - run run_dem_correction() first" else " - run download_data()",
                       call. = FALSE)
  f[1]
}

#' cell areas (Mkm2) per latitude row
area_w <- function(lat, nlon) {
  d <- abs(lat[2] - lat[1])
  (sin(pi / 180 * pmin(pmax(lat + d / 2, -90), 90)) - sin(pi / 180 * pmin(pmax(lat - d / 2, -90), 90))) / 2 / nlon * EARTH_MKM2
}
area_matrix <- function(lat, nlon) matrix(area_w(lat, nlon), length(lat), nlon)

#' nearest grid node of points (returns row, col indices)
grid_index <- function(lat, lon, plat, plon) {
  dl <- abs(lat[2] - lat[1]); nlon <- length(lon)
  i <- pmin(pmax(round((plat - lat[1]) / dl) + 1, 1), length(lat))
  j <- (round((plon - lon[1]) / dl) %% nlon) + 1
  cbind(i, j)
}
sample_grid <- function(z, lat, lon, plat, plon) z[grid_index(lat, lon, plat, plon)]

#' out[i, j] = m[i + di, j + dj]; wraps in longitude, `fill` beyond the poles
shift <- function(m, di = 0, dj = 0, fill = NA) {
  n <- nrow(m); p <- ncol(m)
  if (dj != 0) m <- m[, ((seq_len(p) - 1 + dj) %% p) + 1, drop = FALSE]
  if (di > 0) m <- rbind(m[(1 + di):n, , drop = FALSE], matrix(fill, di, p))
  if (di < 0) m <- rbind(matrix(fill, -di, p), m[1:(n + di), , drop = FALSE])
  m
}

# ---------------------------------------------------------------- smoothing ----
.gauss_kernel <- function(sig) {
  r <- max(1, ceiling(4 * sig)); k <- exp(-0.5 * ((-r:r) / sig)^2); k / sum(k)
}
#' gaussian along latitude (columns), edge values repeated ('nearest')
gauss_lat <- function(m, sig) {
  if (sig <= 0) return(m)
  k <- .gauss_kernel(sig); r <- (length(k) - 1) / 2; n <- nrow(m)
  pad <- rbind(m[rep(1, r), , drop = FALSE], m, m[rep(n, r), , drop = FALSE])
  f <- stats::filter(pad, k, sides = 2)
  matrix(f, nrow(pad))[(r + 1):(r + n), , drop = FALSE]
}
#' gaussian along longitude (rows), periodic (FFT)
gauss_lon <- function(m, sig) {
  if (sig <= 0) return(m)
  p <- ncol(m); fr <- c(0:floor(p / 2), -(ceiling(p / 2) - 1):-1) / p
  H <- exp(-2 * pi^2 * sig^2 * fr^2)
  Re(mvfft(mvfft(t(m)) * H, inverse = TRUE) / p) |> t()
}
#' gaussian with the same sigma in km in both directions on the sphere: the
#' east-west sigma in cells grows with 1/cos(lat), applied per latitude band
smooth_km <- function(a, lat, sigma_km, nbands = 45) {
  nlat <- nrow(a); nlon <- ncol(a); dl <- abs(lat[2] - lat[1])
  sy <- sigma_km / (dl * KM_PER_DEG)
  b <- gauss_lat(a, sy); out <- b
  for (rows in split(seq_len(nlat), cut(seq_len(nlat), nbands, labels = FALSE))) {
    cc <- max(cos(pi / 180 * min(abs(lat[rows]))), 0.05)
    out[rows, ] <- gauss_lon(b[rows, , drop = FALSE], min(sy / cc, nlon / 8))
  }
  out
}
smooth_cells <- function(a, sig) gauss_lon(gauss_lat(a, sig), sig)

# ---------------------------------------------------------------- distances ----
.rast <- function(m, dx = 1, dy = 1) rast(m, extent = ext(0, ncol(m) * dx, 0, nrow(m) * dy), crs = "local")
.mat <- function(r) as.matrix(r, wide = TRUE)

#' Euclidean distance in cells to the nearest TRUE cell of `target` (0 on targets);
#' longitude wrap handled by padding `pad` columns
dist_cells <- function(target, pad = 30) {
  if (!any(target)) return(matrix(Inf, nrow(target), ncol(target)))
  p <- ncol(target); pad <- min(pad, p)
  tp <- cbind(target[, (p - pad + 1):p], target, target[, 1:pad])
  r <- .rast(ifelse(tp, 1, NA))
  d <- .mat(distance(r))
  d[, (pad + 1):(pad + p), drop = FALSE]
}
#' scipy-style EDT: distance (cells) from each TRUE cell of x to the nearest FALSE cell
edt <- function(x, pad = 30) dist_cells(!x, pad)

#' distance in km to the nearest TRUE cell, with the east-west cell size of each
#' latitude band (close to the true distance at all latitudes); wraps in longitude
#' Exact separable Euclidean distance transform with different cell sizes in x and
#' y (terra::distance assumes square cells): vertical distance to the nearest
#' target in each column, then the minimum over horizontal offsets. Exact up to
#' max_km; beyond that only guaranteed to be > max_km.
dist_km <- function(target, lat, rows_per_band = 30, max_km = 700) {
  nlat <- nrow(target); p <- ncol(target); dl <- abs(lat[2] - lat[1]); dy <- dl * KM_PER_DEG
  if (!any(target)) return(matrix(Inf, nlat, p))
  g <- matrix(Inf, nlat, p)                       # rows to the nearest target in the column
  cur <- rep(Inf, p)
  for (i in seq_len(nlat)) { cur <- ifelse(target[i, ], 0, cur + 1); g[i, ] <- cur }
  cur <- rep(Inf, p)
  for (i in nlat:1) { cur <- ifelse(target[i, ], 0, cur + 1); g[i, ] <- pmin(g[i, ], cur) }
  g2 <- (g * dy)^2
  out <- matrix(Inf, nlat, p)
  for (rows in split(seq_len(nlat), ceiling(seq_len(nlat) / rows_per_band))) {
    cc <- max(cos(pi / 180 * min(abs(lat[rows]))), 0.05); dx <- dy * cc
    gb <- g2[rows, , drop = FALSE]
    if (all(!is.finite(gb))) next
    K <- min(ceiling(max_km / dx), p %/% 2)
    best <- gb
    for (k in seq_len(K)) {
      add <- (dx * k)^2
      best <- pmin(best, gb[, ((seq_len(p) - 1 + k) %% p) + 1, drop = FALSE] + add,
                   gb[, ((seq_len(p) - 1 - k) %% p) + 1, drop = FALSE] + add)
    }
    out[rows, ] <- sqrt(best)
  }
  out
}

# ---------------------------------------------------------------- labelling ----
#' connected components of TRUE cells (4 or 8 neighbours), joined across the
#' longitude seam. Returns an integer matrix (0 = background).
label_wrap <- function(mask, directions = 4) {
  if (!any(mask)) return(matrix(0L, nrow(mask), ncol(mask)))
  lab <- .mat(patches(.rast(ifelse(mask, 1, NA)), directions = directions))
  lab[is.na(lab)] <- 0; storage.mode(lab) <- "integer"
  p <- ncol(lab); n <- nrow(lab)
  pairs <- cbind(lab[, 1], lab[, p])
  if (directions == 8) pairs <- rbind(pairs, cbind(lab[-1, 1], lab[-n, p]), cbind(lab[-n, 1], lab[-1, p]))
  pairs <- unique(pairs[pairs[, 1] > 0 & pairs[, 2] > 0 & pairs[, 1] != pairs[, 2], , drop = FALSE])
  if (nrow(pairs)) {
    parent <- seq_len(max(lab))
    find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
    for (k in seq_len(nrow(pairs))) {
      a <- find(pairs[k, 1]); b <- find(pairs[k, 2])
      if (a != b) parent[max(a, b)] <- min(a, b)
    }
    root <- vapply(seq_along(parent), find, 1L)
    lab[lab > 0] <- root[lab[lab > 0]]
  }
  lab
}

#' 3x3 (or (2r+1)^2) maximum filter of a logical mask, wrapping in longitude
dilate <- function(mask, r = 1) {
  out <- mask
  for (di in -r:r) for (dj in -r:r) if (di || dj) out <- out | shift(mask, di, dj, fill = FALSE)
  out
}

# ---------------------------------------------------------- least-cost path ----
#' accumulated least cost from target cells (friction 0) over a friction matrix,
#' 8 neighbours, cost of a step = mean friction of the two cells x step length
cost_distance <- function(friction, maxiter = 100) {
  .mat(costDist(.rast(friction), target = 0, maxiter = maxiter))
}
#' trace the least-cost path from cell (i, j) down the accumulated cost surface
#' to a target cell. Returns a matrix of (row, col), from the start to the target.
trace_path <- function(acc, friction, i, j, max_steps = 20000) {
  n <- nrow(acc); p <- ncol(acc); path <- matrix(NA_integer_, max_steps, 2); k <- 1
  path[1, ] <- c(i, j)
  nb <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  len <- sqrt(rowSums(nb^2))
  while (acc[i, j] > 0 && k < max_steps) {
    ii <- i + nb[, 1]; jj <- j + nb[, 2]
    ok <- ii >= 1 & ii <= n & jj >= 1 & jj <= p
    cand <- rep(Inf, 8)
    cand[ok] <- acc[cbind(ii[ok], jj[ok])] + 0.5 * (friction[i, j] + friction[cbind(ii[ok], jj[ok])]) * len[ok]
    b <- which.min(cand)
    if (!is.finite(cand[b]) || acc[ii[b], jj[b]] >= acc[i, j]) {
      # numerical plateau: step to the lowest neighbour anyway
      v <- rep(Inf, 8); v[ok] <- acc[cbind(ii[ok], jj[ok])]; b <- which.min(v)
      if (v[b] >= acc[i, j]) break
    }
    i <- ii[b]; j <- jj[b]; k <- k + 1; path[k, ] <- c(i, j)
  }
  path[seq_len(k), , drop = FALSE]
}

#' 1-D gaussian of a vector, edge values repeated
gauss1d <- function(x, sig) as.numeric(gauss_lat(matrix(x, ncol = 1), sig))

#' 99th percentile and maximum of the cell-to-cell slope (m/km) where mask, |lat| < 80
max_slopes <- function(z, lat, mask) {
  dy <- abs(lat[2] - lat[1]) * KM_PER_DEG; n <- nrow(z); p <- ncol(z)
  dx <- dy * cos(pi / 180 * lat[-n])
  g <- pmax(abs(diff(z))[, -p] / dy, abs(t(diff(t(z))))[-n, ] / pmax(dx, 1e-6))
  m <- mask[-n, -p] & (abs(lat[-n]) < 80)
  if (!any(m)) return(c(NA, NA))
  c(unname(quantile(g[m], 0.99)), max(g[m]))
}

#' Raise (never lower) the sea floor to a shelf + continental-slope profile: at most
#' MAX_SEA_SLOPE m/km deeper per km away from shallower water down to -1400 m and
#' MAX_DEEP_SLOPE m/km below. The profile is smoothed (smooth_km) so new shelf edges have
#' no cell-by-cell steps. `where`: optional logical mask of cells that may change.
seafloor_profile <- function(g, la, cfg = SILUCARB, where = NULL, smooth_sigma_km = 25) {
  floor_ <- matrix(-Inf, nrow(g), ncol(g))
  for (c_ in c(0, -200, -600, -1000, -1400, -2000, -2600, -3200)) {
    d <- dist_km(g >= c_, la, max_km = 400)
    s_ <- if (c_ > -1400) cfg$MAX_SEA_SLOPE else cfg$MAX_DEEP_SLOPE
    floor_ <- pmax(floor_, min(c_, -1) - s_ * d)
  }
  if (smooth_sigma_km > 0) floor_ <- smooth_km(pmax(floor_, -6000), la, smooth_sigma_km, nbands = 45)
  sea <- g < 0
  if (!is.null(where)) sea <- sea & where
  g[sea] <- pmax(g[sea], pmin(floor_[sea], -1))
  g
}
