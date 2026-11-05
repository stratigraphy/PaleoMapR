# =============================================================================
#  Climate runs with SPEEDY on your own maps (orbit, CO2, Sun).
#
#  SPEEDY (speedy.f90, Kucharski, Molteni & Hatfield; Molteni 2003, Kucharski et al. 2006,
#  2013) is an atmosphere model of intermediate complexity: T30 (3.75 deg), 8 levels, about
#  2-4 minutes per model year on one desktop core. The version in tools/speedy reads and
#  writes plain binary files (no NetCDF library needed), writes monthly means and a progress
#  file, and has settings for the orbit, CO2, the Sun and a slab ocean. Windows and Linux
#  programs are included.
#
#  The climate of a map = HadCM3L run of that age (Valdes et al. 2021) + SPEEDY response:
#     orig   SPEEDY on the original PALEOMAP map  } same HadCM3L sea-surface temperatures, sea
#     corr   SPEEDY on your map                   } ice, soil moisture and snow (monthly)
#     corr - orig = response to your changes of the map (removes most of SPEEDY's own biases)
#  and for an experiment (see climate_experiment()):
#     forcing   CO2 / Sun changed: sea surface shifted by the response of a chosen climate
#               sensitivity, SPEEDY's atmosphere responds to that and to the CO2 / Sun itself
#     slab, exp orbit changed: slab-ocean runs with today's and with the chosen orbit
#  Runs are independent processes, so they go in parallel on a multicore PC; they keep running
#  in the background when the editor is closed. Results: results/speedy/.
# =============================================================================

#' 48 Gaussian latitudes of the T30 grid, north to south (as in SPEEDY's input files)
speedy_lats <- local({
  n <- 48; x <- cos(pi * (seq_len(n) - 0.25) / (n + 0.5))
  for (it in 1:100) {
    p0 <- 1; p1 <- x
    for (k in 2:n) { p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 <- p1; p1 <- p2 }
    dp <- n * (x * p1 - p0) / (x^2 - 1); dx <- p1 / dp; x <- x - dx
    if (max(abs(dx)) < 1e-14) break
  }
  asin(x) * 180 / pi
})
SPEEDY_LON <- (0:95) * 3.75

#' the SPEEDY program for this computer (Windows and Linux are included; on macOS it is
#' compiled once from tools/speedy/src if gfortran is installed)
speedy_exe <- function(cfg = SILUCARB) {
  d <- file.path(cfg$dir_tools, "speedy")
  if (.Platform$OS.type == "windows") {
    exe <- file.path(d, "windows", "speedy.exe"); nw <- file.path(d, "windows", "speedy_new.exe")
    # an update that could not be copied over the program while it ran: swap it in once no SPEEDY runs
    if (file.exists(nw)) { old <- file.path(d, "windows", "speedy_old.exe"); unlink(old)
      if (!file.exists(exe) || isTRUE(file.rename(exe, old))) { if (file.rename(nw, exe)) message("SPEEDY program updated") } }
    return(exe)
  }
  if (Sys.info()[["sysname"]] == "Linux") return(file.path(d, "linux", "speedy"))
  exe <- file.path(d, "macos", "speedy")
  if (!file.exists(exe)) {
    fc <- Sys.which("gfortran")
    if (!nzchar(fc)) stop("SPEEDY needs gfortran on macOS (e.g. 'brew install gcc'), then run again.", call. = FALSE)
    src <- file.path(d, "src"); ord <- scan(file.path(src, "order.txt"), what = "", quiet = TRUE)
    b <- file.path(tempdir(), "speedy_build"); dir.create(b, showWarnings = FALSE)
    old <- setwd(b); on.exit(setwd(old))
    system2(fc, c("-O3", "-ffast-math", "-c", shQuote(file.path(src, ord))))
    dir.create(dirname(exe), showWarnings = FALSE, recursive = TRUE)
    system2(fc, c("-O3", "-o", shQuote(exe), "*.o"))
    if (!file.exists(exe)) stop("compiling SPEEDY failed", call. = FALSE)
  }
  exe
}

# ---- boundary conditions --------------------------------------------------------------
.MONTHS <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

#' the 12 monthly HadCM3L atmosphere files of the run for this age (downloaded once, ~16 MB)
climate_month_files <- function(age, cfg = SILUCARB) {
  e <- climate_exp(age); d <- file.path(cfg$dir_data, "climate", e$id)
  f <- file.path(d, sprintf("%sa.pdcl%s.nc", e$id, .MONTHS))
  for (k in seq_along(f)) .fetch_retry(sprintf(CLIMATE_URL, e$id, basename(f[k])), f[k])
  if (!all(file.exists(f))) stop("monthly climate files for ", age, " Ma could not be downloaded", call. = FALSE)
  f
}

#' fill the cells where `keep` is FALSE with values spread from the kept cells (repeated
#' neighbour averages, wraps in longitude): e.g. sea temperature carried under former land
.fill_from <- function(z, keep, iters = 400) {
  z[!keep] <- NA; nr <- nrow(z); nc <- ncol(z)
  for (it in seq_len(iters)) {
    na <- is.na(z); if (!any(na)) break
    s <- matrix(0, nr, nc); n <- s
    for (sh in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))) {
      i <- pmin(pmax(seq_len(nr) + sh[1], 1), nr); j <- ((seq_len(nc) - 1 + sh[2]) %% nc) + 1
      v <- z[i, j]; ok <- !is.na(v); s[ok] <- s[ok] + v[ok]; n <- n + ok
    }
    upd <- na & n > 0; z[upd] <- s[upd] / n[upd]
  }
  z[is.na(z)] <- mean(z, na.rm = TRUE)
  z
}

#' HadCM3L field -> SPEEDY grid (96 x 48, longitude from 0 E, latitude north to south)
.to_speedy <- function(F) {
  g <- expand.grid(lon = ((SPEEDY_LON + 180) %% 360) - 180, lat = speedy_lats)
  matrix(clim_interp(F, g$lat, g$lon), 96, 48)
}

#' monthly forcing shared by both runs of an age (cached in the run folder)
speedy_forcing <- function(age, cfg = SILUCARB) {
  f <- climate_month_files(age, cfg)
  out <- list(sst = array(0, c(96, 48, 12)), icec = array(0, c(96, 48, 12)), stl = array(0, c(96, 48, 12)),
              swl1 = array(0, c(96, 48, 12)), swl2 = array(0, c(96, 48, 12)), snowd = array(0, c(96, 48, 12)))
  for (m in 1:12) {
    ts <- .read_clim_var(f[m], "temp_mm_srf"); ic <- .read_clim_var(f[m], "iceconc_mm_srf")
    sn <- .read_clim_var(f[m], "snowdepth_mm_srf")
    land <- is.finite(.lay_land(f[m]))
    # sea values carried under (former) land and land values carried over (new) land, so both runs share them
    sst <- ts; sst$z <- .fill_from(ts$z, !land); ice <- ic; ice$z <- pmin(pmax(.fill_from(ic$z, !land), 0), 1)
    stl <- ts; stl$z <- .fill_from(ts$z, land); sno <- sn; sno$z <- pmax(.fill_from(sn$z, land), 0)
    s1 <- .soil_layer(f[m], 1, 0.10); s2 <- .soil_layer(f[m], 2, 0.25)
    s1$z <- pmin(pmax(.fill_from(s1$z, land), 0), 0.45); s2$z <- pmin(pmax(.fill_from(s2$z, land), 0), 0.45)
    out$sst[, , m] <- .to_speedy(sst); out$icec[, , m] <- .to_speedy(ice); out$stl[, , m] <- .to_speedy(stl)
    out$swl1[, , m] <- .to_speedy(s1); out$swl2[, , m] <- .to_speedy(s2); out$snowd[, , m] <- .to_speedy(sno)
  }
  attr(out, "age") <- age
  out
}
#' HadCM3L soil layer k as volumetric water content (NA over sea), oriented like .read_clim_var
.soil_layer <- function(file, k, thick) {
  n <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(n))
  v <- n$var[["sm_mm_soil"]]; x <- ncdf4::ncvar_get(n, "sm_mm_soil")[, , k]
  lon <- as.vector(v$dim[[1]]$vals); lat <- as.vector(v$dim[[2]]$vals)
  x[abs(x) > 1e19] <- NA; m <- t(x); o <- order(lat); m <- m[o, , drop = FALSE]; lat <- lat[o]
  lon <- ((lon + 180) %% 360) - 180; p <- order(lon)
  list(z = m[, p, drop = FALSE] / 1000 / thick, lat = lat, lon = lon[p])
}
.lay_land <- function(file) .soil_layer(file, 1, 1)$z

#' land fraction and mean land height on the SPEEDY grid from a PaleoDEM file (or a list
#' z, lat, lon); ice: optional ice thickness (m) on the same grid -> ice fraction and the
#' surface height includes the ice
speedy_geography <- function(dem_file, ice = NULL) {
  d <- if (is.list(dem_file)) dem_file else read_dem_nc(dem_file)
  lb <- c(90, (speedy_lats[-1] + speedy_lats[-48]) / 2, -90)          # cell edges, north to south
  ii <- 49 - findInterval(d$lat, rev(lb), rightmost.closed = TRUE)
  ii <- pmin(pmax(ii, 1), 48)
  jj <- floor((((d$lon + 1.875) %% 360)) / 3.75) + 1
  agg <- function(m) { a <- rowsum(m, ii, reorder = TRUE); t(rowsum(t(a), jj, reorder = TRUE)) }   # 48 x 96
  n <- agg(matrix(1, nrow(d$z), ncol(d$z)))
  land <- d$z >= 0; surf <- pmax(d$z, 0)
  if (!is.null(ice)) { land <- land | ice > 10; surf <- pmax(d$z + ice, 0) }
  out <- list(lsm = t(agg(land * 1) / n), orog = t(agg(surf) / n))                                    # 96 x 48
  if (!is.null(ice)) out$ice <- t(agg((ice > 10) * 1) / n)
  out
}

.wbin <- function(x, file) { con <- file(file, "wb"); on.exit(close(con)); writeBin(as.vector(x), con, size = 4, endian = "little") }

#' sun ~3.5 % fainter in the Silurian (Gough 1981)
solar_factor_at <- function(age) 1 / (1 + 0.4 * (1 - (4.57 - age / 1000) / 4.57))

# ---- maps for climate runs -------------------------------------------------------------
#' a map for climate runs: the editor's corrected grid of an age, or a map file (NetCDF, e.g.
#' exported with export_map() on another computer and loaded with import_map())
climate_map <- function(age, file = NULL, name = NULL, cfg = SILUCARB) {
  if (is.list(age)) return(age)
  if (is.null(file)) {
    file <- tryCatch(find_dem(cfg$dir_corrected, age), error = function(e) NULL)
    if (is.null(file)) stop("no corrected grid for ", age, " Ma yet: press 'Rebuild this age' in the editor first ",
                            "(or load a map file)", call. = FALSE)
    name <- name %||% sprintf("%s Ma, editor map", format(age))
  }
  if (!file.exists(file)) stop("map file not found: ", file, call. = FALSE)
  list(age = age, file = normalizePath(file, winslash = "/"), name = name %||% basename(file), tag = .map_tag(file))
}
#' short fingerprint of a map file's content (runs are kept per map version)
.MAPTAG <- new.env()
.map_tag <- function(file) {
  i <- file.info(file); key <- paste(normalizePath(file, winslash = "/"), as.numeric(i$mtime), i$size)
  if (is.null(.MAPTAG[[key]])) {                       # from the heights only (not the file's notes)
    tf <- tempfile(); on.exit(unlink(tf)); z <- read_dem_nc(file)$z
    con <- file(tf, "wb"); writeBin(as.integer(round(z)), con); writeBin(dim(z), con); close(con)
    .MAPTAG[[key]] <- paste0("map", substr(unname(tools::md5sum(tf)), 1, 10))
  }
  .MAPTAG[[key]]
}

#' write the corrected map of an age as a self-contained NetCDF file (with its age and
#' source) that can be loaded on another computer with import_map()
export_map <- function(age, file = NULL, cfg = SILUCARB) {
  src <- find_dem(cfg$dir_corrected, age)
  dir.create(cfg$dir_exports, showWarnings = FALSE, recursive = TRUE)
  file <- file %||% file.path(cfg$dir_exports, sprintf("paleomap_%sMa_%s.nc", format(age), format(Sys.Date(), "%Y%m%d")))
  file.copy(src, file, overwrite = TRUE)
  n <- ncdf4::nc_open(file, write = TRUE); on.exit(ncdf4::nc_close(n))
  ncdf4::ncatt_put(n, 0, "paleomap_age_Ma", as.numeric(age))
  ncdf4::ncatt_put(n, 0, "paleomap_plate_model", "PALEOMAP (Scotese & Wright 2018)")
  ncdf4::ncatt_put(n, 0, "paleomap_source", "paleogeographic map editor: PaleoDEM corrected with fossil / rock records")
  ncdf4::ncatt_put(n, 0, "paleomap_created", format(Sys.time(), "%Y-%m-%d %H:%M"))
  invisible(file)
}
#' the age stored in a map file (NULL if none)
map_file_age <- function(file) {
  n <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(n))
  a <- ncdf4::ncatt_get(n, 0, "paleomap_age_Ma")
  if (isTRUE(a$hasatt)) return(as.numeric(a$value))
  m <- regmatches(basename(file), regexpr("[0-9.]+(?=Ma)", basename(file), perl = TRUE))
  if (length(m)) as.numeric(m) else NULL
}
#' copy a map file into the project (folder map_files) so climate runs can use it;
#' age: the map's age in Ma (read from the file if it was made with export_map())
import_map <- function(file, age = NULL, name = NULL, cfg = SILUCARB) {
  age <- age %||% map_file_age(file)
  if (is.null(age)) stop("give the age of this map (age = ...)", call. = FALSE)
  d <- read_dem_nc(file)                                # checks that it is a map grid
  if (is.null(d$z) || length(dim(d$z)) != 2) stop("this file is not a map grid", call. = FALSE)
  dir.create(cfg$dir_map_files, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cfg$dir_map_files, basename(name %||% file))
  if (!grepl("\\.nc$", dest)) dest <- paste0(dest, ".nc")
  if (normalizePath(file, winslash = "/", mustWork = FALSE) != normalizePath(dest, winslash = "/", mustWork = FALSE))
    file.copy(file, dest, overwrite = TRUE)
  climate_map(age, dest, name = sub("\\.nc$", "", basename(dest)), cfg = cfg)
}
#' the map files loaded so far (folder map_files), with their ages
map_files <- function(cfg = SILUCARB) {
  f <- list.files(cfg$dir_map_files, "\\.nc$", full.names = TRUE)
  if (!length(f)) return(data.frame(name = character(0), age = numeric(0), file = character(0)))
  data.frame(name = sub("\\.nc$", "", basename(f)),
             age = vapply(f, function(x) map_file_age(x) %||% NA_real_, 0), file = f, row.names = NULL)
}

# ---- experiments: orbit, CO2 and the Sun -------------------------------------------------
ORBIT_PRESENT <- c(ecc = 0.016715, obl = 23.441, peri = 102.7)
#' preset orbital configurations. peri = longitude of perihelion (precession), measured from the
#' NH spring equinox as in Berger (1978): 90 = perihelion at the December solstice -> strong
#' southern summers, 270 = at the June solstice -> strong northern summers; 102.7 today
ORBIT_PRESETS <- list(
  present   = list(label = "present-day orbit", ecc = 0.016715, obl = 23.441, peri = 102.7),
  sh_strong = list(label = "strong southern summers (e 0.05, obliquity 24.5, perihelion in Dec.)", ecc = 0.05, obl = 24.5, peri = 90),
  nh_strong = list(label = "strong northern summers (e 0.05, obliquity 24.5, perihelion in June)", ecc = 0.05, obl = 24.5, peri = 270),
  weak      = list(label = "weak seasons (circular orbit, obliquity 22.1)", ecc = 0, obl = 22.1, peri = 102.7),
  high_obl  = list(label = "high obliquity only (24.5, present eccentricity)", ecc = 0.016715, obl = 24.5, peri = 102.7),
  low_obl   = list(label = "low obliquity only (22.1, present eccentricity)", ecc = 0.016715, obl = 22.1, peri = 102.7))
CO2_EXPONENT <- 1.1      # SPEEDY CO2-band absorptivity ~ CO2^1.1 gives 3.7 W/m2 per doubling (calibrated)
#' CO2 and Sun: SPEEDY with a slab ocean is far too sensitive at the warm Silurian temperatures
#' (it does not settle), so for a CO2 or solar change the sea surface is warmed (or cooled) by
#' the equilibrium response of a chosen climate sensitivity (K per CO2 doubling, i.e. per
#' 3.7 W/m2; HadCM3 ~3.3), more at high latitudes (0.75 + 0.75 sin^2 lat, global mean 1), land
#' 1.3 x; SPEEDY's atmosphere then responds to that and to the CO2 / Sun itself. Orbital
#' experiments use the slab ocean.
CLIMATE_SENSITIVITY <- 3
#' run lengths (model years): spin-up (not averaged) and years averaged
SPEEDY_LENGTHS <- list(
  quick    = list(label = "quick (~10 min)",      reference = c(spinup = 1, years = 2), experiment = c(spinup = 2, years = 2)),
  standard = list(label = "standard (~30-40 min)", reference = c(spinup = 1, years = 3), experiment = c(spinup = 6, years = 4)),
  long     = list(label = "long (~1.5 h)",        reference = c(spinup = 2, years = 8), experiment = c(spinup = 10, years = 10)),
  decade   = list(label = "decade (10-yr means, ~1 h)", reference = c(spinup = 1, years = 10), experiment = c(spinup = 5, years = 10)))

#' HadCM3L CO2 (ppm) of the run used for a map age
hadcm3l_co2 <- function(age) climate_exp(age)$co2

#' an experiment (orbit, CO2, Sun) for climate runs
#' orbit:   a preset (names(ORBIT_PRESETS)) or list(ecc, obl, peri); ecc / obl / peri override it
#'          (ecc: eccentricity, obl: obliquity in degrees, peri: longitude of perihelion in degrees)
#' co2_pal: CO2 in PAL (1 PAL = 280 ppm); NULL = as in the HadCM3L run of the map's age
#' co2:     alternatively a factor relative to the HadCM3L run (ignored when co2_pal is given)
#' sun_age: age (Ma) for the strength of the Sun (fainter in the past, Gough 1981); NULL = the map's age
#' sensitivity: K per CO2 doubling for the CO2 / Sun response of the sea surface
#' length:  "quick", "standard" or "long"
#' gradient: equator-to-pole factor for the HadCM3L surface temperatures (1 = as HadCM3L). Each
#'          latitude's annual zonal-mean departure from the global mean is multiplied by it, so
#'          < 1 warms the poles and cools the tropics with the global mean kept (bias correction)
climate_experiment <- function(orbit = "present", co2 = 1, ecc = NULL, obl = NULL, peri = NULL,
                               sensitivity = CLIMATE_SENSITIVITY, co2_pal = NULL, sun_age = NULL,
                               length = "standard", gradient = 1) {
  o <- if (is.list(orbit)) orbit else ORBIT_PRESETS[[orbit]]
  if (is.null(o)) stop("orbit must be one of: ", paste(names(ORBIT_PRESETS), collapse = ", "), call. = FALSE)
  if (is.null(SPEEDY_LENGTHS[[length]])) stop("length must be one of: ", paste(names(SPEEDY_LENGTHS), collapse = ", "), call. = FALSE)
  e <- list(ecc = ecc %||% o$ecc, obl = obl %||% o$obl, peri = (peri %||% o$peri) %% 360, co2 = co2, co2_pal = co2_pal,
            sun_age = sun_age, sensitivity = sensitivity, length = length, gradient = gradient %||% 1)
  match_preset <- Filter(function(p) isTRUE(all.equal(c(p$ecc, p$obl, p$peri), c(e$ecc, e$obl, e$peri))), ORBIT_PRESETS)
  e$orbit_label <- if (length(match_preset)) match_preset[[1]]$label else
    sprintf("e %.3f, obliquity %.2f, perihelion %.0f", e$ecc, e$obl, e$peri)
  e$label <- e$orbit_label
  e
}
.present_orbit <- function(x) isTRUE(all.equal(c(x$ecc, x$obl, x$peri), unname(ORBIT_PRESENT), tolerance = 1e-6))
#' an experiment for a map age: CO2 factor relative to HadCM3L, Sun age, and the SST shift (K)
#' that the CO2 / Sun change brings (global mean)
.resolve <- function(x, age) {
  x <- x %||% climate_experiment()
  x$co2f <- if (!is.null(x$co2_pal)) x$co2_pal * PAL_PPM / hadcm3l_co2(age) else (x$co2 %||% 1)
  x$sun <- x$sun_age %||% age
  dF_sun <- 1361 / 4 * 0.7 * (solar_factor_at(x$sun) / solar_factor_at(age) - 1)   # W/m2, planetary albedo 0.3
  x$dT <- (x$sensitivity %||% CLIMATE_SENSITIVITY) * (log2(x$co2f) + dF_sun / 3.7)
  x$gradient <- x$gradient %||% 1
  x$forced <- abs(x$co2f - 1) > 1e-3 || abs(x$sun - age) > 0.05 || abs(x$gradient - 1) > 1e-3
  x$co2_ppm <- x$co2f * hadcm3l_co2(age); x$co2_pal_eff <- x$co2_ppm / PAL_PPM
  x$label <- sprintf("%s; CO2 %.3g PAL (%.0f ppm)%s%s", x$orbit_label, x$co2_pal_eff, x$co2_ppm,
                     if (abs(x$sun - age) > 0.05) sprintf("; Sun of %s Ma", format(x$sun)) else "",
                     if (abs(x$gradient - 1) > 1e-3) sprintf("; gradient \u00d7%.2f", x$gradient) else "")
  x$shift <- c(dT = x$dT, grad = x$gradient)
  x
}
.is_reference <- function(x, age = NULL) {
  if (is.null(x)) return(TRUE)
  if (is.null(age)) return(.present_orbit(x) && is.null(x$co2_pal) && (x$co2 %||% 1) == 1 && is.null(x$sun_age) && (x$gradient %||% 1) == 1)
  r <- .resolve(x, age); .present_orbit(r) && !r$forced
}
.forcing_tag <- function(r) sprintf("c%.3f_sun%s_s%.1f%s", r$co2f, format(r$sun), r$sensitivity %||% CLIMATE_SENSITIVITY,
                                     if (abs((r$gradient %||% 1) - 1) > 1e-3) sprintf("_g%.2f", r$gradient) else "")
.exp_tag <- function(r) {
  if (is.null(r$co2f)) stop(".exp_tag needs a resolved experiment")
  sprintf("e%.4f_o%.2f_p%.0f_%s", r$ecc, r$obl, r$peri, .forcing_tag(r))
}

#' HadCM3L zonal-mean surface temperature (degC) at a map age by month (12 x lats) and its
#' area-weighted global means by month (for the equator-to-pole gradient factor)
.ZPROF <- new.env()
hadcm3l_zonal <- function(age, lats = speedy_lats, cfg = SILUCARB) {
  key <- format(age)
  if (is.null(.ZPROF[[key]])) {
    mf <- climate_month_files(age, cfg)
    Z <- t(vapply(mf, function(f) { t <- .read_clim_var(f, "temp_mm_srf"); rowMeans(t$z, na.rm = TRUE) - 273.15 },
                  numeric(length(.read_clim_var(mf[1], "temp_mm_srf")$lat))))
    lat <- .read_clim_var(mf[1], "temp_mm_srf")$lat; w <- cos(lat * pi / 180)
    .ZPROF[[key]] <- list(lat = lat, T = unname(Z), Tg = as.vector(Z %*% w) / sum(w))
  }
  p <- .ZPROF[[key]]
  Tm <- t(apply(p$T, 1, function(v) approx(p$lat, v, lats, rule = 2)$y))
  list(lat = lats, T = Tm, Tg = p$Tg, T_ann = colMeans(Tm), Tg_ann = mean(p$Tg))
}
#' extra temperature change (K) by latitude from the gradient factor g:
#' (g - 1) (T_zonal - T_global) of month `month` (1-12), or of the annual mean (month = NULL)
gradient_anomaly <- function(age, g, lats = speedy_lats, month = NULL, cfg = SILUCARB) {
  if (is.null(g) || abs(g - 1) < 1e-6) return(rep(0, length(lats)))
  p <- hadcm3l_zonal(age, lats, cfg)
  a <- if (is.null(month)) (g - 1) * (p$T_ann - p$Tg_ann) else (g - 1) * (p$T[month, ] - p$Tg[month])
  w <- cos(lats * pi / 180); a - sum(w * a) / sum(w)              # keep the global mean exactly
}
#' text for the editor: HadCM3L's tropics / poles and what a gradient factor does to them
gradient_text <- function(age, g = 1, cfg = SILUCARB) {
  lats <- seq(-89.5, 89.5, 1); p <- hadcm3l_zonal(age, lats, cfg); a <- gradient_anomaly(age, g, lats, cfg = cfg)
  p$T <- p$T_ann; p$Tg <- p$Tg_ann
  band <- function(v, lo, hi) { k <- abs(lats) >= lo & abs(lats) <= hi; w <- cos(lats[k] * pi / 180); sum(w * v[k]) / sum(w) }
  sh <- function(v, lo, hi) { k <- lats <= -lo & lats >= -hi; mean(v[k]) }
  nh <- function(v, lo, hi) { k <- lats >= lo & lats <= hi; mean(v[k]) }
  f <- function(v) c(trop = band(v, 0, 20), s = sh(v, 70, 90), n = nh(v, 70, 90))
  t0 <- f(p$T); t1 <- f(p$T + a)
  paste0(sprintf("HadCM3L at %s Ma (annual mean surface temperature): tropics (0-20\u00b0) %.1f \u00b0C, 70-90\u00b0S %.1f \u00b0C, 70-90\u00b0N %.1f \u00b0C, global %.1f \u00b0C; equator-to-pole difference %.0f K. ",
                 format(age), t0[["trop"]], t0[["s"]], t0[["n"]], p$Tg, t0[["trop"]] - (t0[["s"]] + t0[["n"]]) / 2),
         if (abs(g - 1) > 1e-3) sprintf("With factor %.2f (before CO2 / Sun): tropics %.1f, 70-90\u00b0S %.1f, 70-90\u00b0N %.1f \u00b0C (difference %.0f K, same global mean).",
                                     g, t1[["trop"]], t1[["s"]], t1[["n"]], t1[["trop"]] - (t1[["s"]] + t1[["n"]]) / 2)
         else "Factor 1 = HadCM3L as it is.")
}
.shift_dT <- function(shift) if (is.null(shift)) 0 else unname(shift[1])
.shift_grad <- function(shift) if (length(shift) >= 2) unname(shift[2]) else 1

#' warm (cool) the prescribed sea and land temperatures by dT (global mean of the sea; see
#' CLIMATE_SENSITIVITY) and, with shift = c(dT, grad), scale the equator-to-pole contrast
shift_forcing <- function(forcing, shift, age = attr(forcing, "age")) {
  dT <- .shift_dT(shift); g <- .shift_grad(shift); grad <- abs(g - 1) > 1e-6
  if (abs(dT) < 1e-6 && !grad) return(forcing)
  pat <- matrix(0.75 + 0.75 * sin(speedy_lats * pi / 180)^2, 96, 48, byrow = TRUE)
  tf <- 273.2 - 1.8
  for (m in 1:12) {
    # gradient factor: this month's zonal departure from the global mean, scaled (less contrast
    # between latitudes and a smaller seasonal cycle towards the poles when g < 1)
    gam <- if (grad) matrix(gradient_anomaly(age, g, speedy_lats, m), 96, 48, byrow = TRUE) else 0
    d_sea <- dT * pat + gam; d_land <- 1.3 * dT * pat + gam
    forcing$sst[, , m] <- forcing$sst[, , m] + d_sea
    forcing$stl[, , m] <- forcing$stl[, , m] + d_land
    ic <- forcing$icec[, , m]; sst <- forcing$sst[, , m]
    frz <- pmin(1, pmax(0, (tf - sst) / 2))                                    # sea surface below freezing
    if (grad) {
      ic <- pmax(ic * pmax(0, 1 - 0.25 * pmax(d_sea, 0)), frz)                  # less ice when warmer, but ice where it freezes
      ic[d_sea < 0] <- pmax(ic[d_sea < 0], frz[d_sea < 0])
      ic[sst > tf + 2] <- 0
    } else {
      if (dT > 0) ic <- ic * pmax(0, 1 - 0.25 * dT * pat)                      # less sea ice when warmer
      if (dT < 0) ic <- pmax(ic, frz)                                           # more sea ice where the sea cools below freezing
    }
    forcing$icec[, , m] <- ic
  }
  forcing
}

#' write all input files of one run (x: a resolved experiment, see .resolve)
speedy_prepare <- function(dir, age, dem_file, forcing, years = 3, spinup = 1, x = NULL, sea_mode = 1, hfsea = NULL,
                           geo = NULL) {
  x <- x %||% .resolve(NULL, age)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  g <- geo %||% speedy_geography(dem_file)                             # geo: e.g. with an ice sheet
  .wbin(g$lsm, file.path(dir, "lsm.bin")); .wbin(g$orog, file.path(dir, "orog.bin"))
  bare <- age > 420; early <- age > 385                                # little or no land plants before the Devonian
  alb <- matrix(if (bare) 0.25 else if (early) 0.20 else 0.17, 96, 48)
  if (!is.null(g$ice)) alb <- alb * (1 - g$ice) + ICE_ALBEDO * g$ice  # ice-covered land
  .wbin(alb, file.path(dir, "alb.bin"))
  .wbin(matrix(if (early) 0 else 0.3, 96, 48), file.path(dir, "vegh.bin"))
  .wbin(matrix(if (bare) 0.05 else if (early) 0.2 else 0.3, 96, 48), file.path(dir, "vegl.bin"))
  for (v in names(forcing)) .wbin(forcing[[v]], file.path(dir, paste0(v, ".bin")))
  if (!is.null(hfsea)) .wbin(hfsea, file.path(dir, "hfsea.bin"))
  y0 <- 1981
  writeLines(c("&params", "nsteps_out = 99999999", "nstdia = 99999999", "/",
               "&date", sprintf("start_datetime%%year = %d", y0), "start_datetime%month = 1", "start_datetime%day = 1",
               "start_datetime%hour = 0", "start_datetime%minute = 0",
               sprintf("end_datetime%%year = %d", y0 + spinup + years), "end_datetime%month = 1", "end_datetime%day = 1",
               "end_datetime%hour = 0", "end_datetime%minute = 0", "/",
               "&paleo", sprintf("spinup_years = %d", spinup), sprintf("solar_factor = %.5f", solar_factor_at(x$sun)),
               sprintf("eccentricity = %.6f", x$ecc), sprintf("obliquity = %.4f", x$obl), sprintf("perihelion = %.3f", x$peri),
               sprintf("co2_factor = %.5f", x$co2f), sprintf("co2_exponent = %.3f", CO2_EXPONENT),
               sprintf("sea_mode = %d", sea_mode), sprintf("tdsst_days = %.1f", if (sea_mode >= 3) 1e6 else 90), sprintf("restart_days = %d", SPEEDY_RESTART_DAYS),
               sprintf("ice_growth_sh = %d", as.integer(isTRUE((x$ice_growth %||% c(south = TRUE))[["south"]]))),
               sprintf("ice_growth_nh = %d", as.integer(isTRUE((x$ice_growth %||% c(north = TRUE))[["north"]]))), "/"),
             file.path(dir, "namelist.nml"))
  writeLines(c(paste("age", age), paste("dem", normalizePath(dem_file, winslash = "/")), paste("years", years),
               paste("spinup", spinup), paste("experiment", x$label), paste("sea_mode", sea_mode),
               paste("sst_shift_K", round(x$dT %||% 0, 3)), paste("gradient_factor", x$gradient %||% 1)),
             file.path(dir, "run_info.txt"))
}

# ---- runs ---------------------------------------------------------------------------
speedy_dir <- function(age, tag, cfg = SILUCARB) file.path(cfg$dir_results, "speedy", sprintf("%sMa_%s", format(age), tag))

#' the runs needed for a map (an age, or climate_map()) and an experiment:
#'   orig     reference on the original PALEOMAP map, HadCM3L sea-surface temperatures
#'   corr     the same on your map; corr - orig = the response to your changes of the map
#'   forcing  (CO2 / Sun changed) your map with the sea surface shifted by the CO2 / Sun response
#'   slab, exp (orbit changed) slab-ocean runs with today's and with the chosen orbit, around the
#'            reference (forcing run if there is one, else corr); exp - slab = orbit effect
speedy_jobs <- function(map, x = NULL, cfg = SILUCARB) {
  m <- climate_map(map, cfg = cfg); age <- m$age
  L <- SPEEDY_LENGTHS[[(x %||% climate_experiment())$length %||% "standard"]]; r <- L$reference; e <- L$experiment
  yt <- function(v) sprintf("_y%d+%d", v[["spinup"]], v[["years"]])
  ref <- .resolve(climate_experiment(length = (x %||% climate_experiment())$length), age)
  mk <- function(tag, dem, sea_mode, xr, v, ...) list(tag = tag, dem = dem, sea_mode = sea_mode, x = xr,
                                                      spinup = v[["spinup"]], years = v[["years"]], ...)
  jobs <- list(orig = mk(paste0("orig", yt(r)), find_dem(cfg$dir_dem, age), 1, ref, r),
               corr = mk(paste0(m$tag, yt(r)), m$file, 1, ref, r))
  if (!.is_reference(x, age)) {
    xr <- .resolve(x, age); base <- jobs$corr$tag; shift <- 0
    if (xr$forced) {
      xf <- .resolve(modifyList(x, c(as.list(ORBIT_PRESENT), orbit_label = "present-day orbit")), age)          # CO2 / Sun only, today's orbit
      jobs$forcing <- mk(paste0(m$tag, "_", .forcing_tag(xr), yt(r)), m$file, 1, xf, r, shift = xr$shift)
      base <- jobs$forcing$tag; shift <- xr$shift
    }
    if (!.present_orbit(xr)) {
      xs <- .resolve(modifyList(x, c(as.list(ORBIT_PRESENT), orbit_label = "present-day orbit")), age)
      jobs$slab <- mk(paste0(base, "_slab", yt(e)), m$file, 3, xs, e, shift = shift, after = base)
      jobs$exp <- mk(paste0(m$tag, "_", .exp_tag(xr), yt(e)), m$file, 3, xr, e, shift = shift, after = base)
    }
  }
  lapply(jobs, function(j) { j$age <- age; j$dir <- speedy_dir(age, j$tag, cfg); j })
}

SPEEDY_STALL_MIN <- 20     # a run without progress for this long is taken as stopped (e.g. the PC restarted, or a virus scanner stopped it)
SPEEDY_RESTART_DAYS <- 30  # SPEEDY writes a restart file every this many model days; a stopped run continues from it
SPEEDY_MAX_RESUMES <- 25   # a run that had to be resumed this often is taken as failed (something keeps stopping it)
.has_restart <- function(dir) any(file.exists(file.path(dir, c("restart.bin", "restart.tmp"))))
.n_resumes <- function(dir) { f <- file.path(dir, "resumes.txt"); if (!file.exists(f)) return(0L)
  n <- suppressWarnings(as.integer(tryCatch(readLines(f, warn = FALSE)[1], error = function(e) NA))); if (is.na(n)) 0L else n }
.speedy_pid <- function(dir) { f <- file.path(dir, "pid.txt"); if (!file.exists(f)) return(NA_integer_)
  suppressWarnings(as.integer(tryCatch(readLines(f, warn = FALSE)[1], error = function(e) NA))) }
#' is the SPEEDY program of a run folder still alive (its process id is running speedy)?
.speedy_alive <- function(dir) {
  if (!file.exists(file.path(dir, "running.lock"))) return(FALSE)
  pid <- .speedy_pid(dir); if (is.na(pid)) return(FALSE)
  if (.Platform$OS.type == "windows") {
    out <- tryCatch(suppressWarnings(system2("tasklist", c("/FI", shQuote(sprintf("PID eq %d", pid)), "/NH"), stdout = TRUE, stderr = FALSE)), error = function(e) "")
    any(grepl("speedy", out, ignore.case = TRUE))
  } else isTRUE(tryCatch(tools::pskill(pid, 0), error = function(e) FALSE))
}
.speedy_kill <- function(dir) { pid <- .speedy_pid(dir); if (is.na(pid)) return(invisible(FALSE))
  if (.Platform$OS.type == "windows") suppressWarnings(system2("taskkill", c("/F", "/PID", pid), stdout = FALSE, stderr = FALSE)) else tools::pskill(pid)
  invisible(TRUE) }
#' how a run's SPEEDY program ended each time (written by the launcher: date, time, exit code), and
#' how often it was resumed - to find out what stops runs (virus scanner, sleep, restarts)
speedy_run_history <- function(dir) {
  ex <- file.path(dir, "exitcodes.txt"); e <- if (file.exists(ex)) tryCatch(readLines(ex, warn = FALSE), error = function(e) character(0)) else character(0)
  list(resumes = .n_resumes(dir), exits = e, alive = .speedy_alive(dir), restart_file = .has_restart(dir))
}
SPEEDY_RETRIES <- 3        # a run that stopped on a passing file problem (e.g. a locked file on Windows) is started again this often
.retry_file <- function(dir) paste0(dir, "_retries.txt")      # beside the run folder (that is emptied on a restart)
.n_retries <- function(dir) { f <- .retry_file(dir); if (!file.exists(f)) return(0L)
  n <- suppressWarnings(as.integer(tryCatch(readLines(f, warn = FALSE)[1], error = function(e) NA))); if (is.na(n)) 0L else n }
.job_status <- function(j) {
  mb <- file.path(j$dir, "means.bin")
  if (file.exists(mb) && (file.exists(file.path(j$dir, "means_vars.txt")) ||            # fully written
                          difftime(Sys.time(), file.info(mb)$mtime, units = "secs") > 60)) return("done")
  if (file.exists(file.path(j$dir, "started.txt"))) {
    if (file.exists(file.path(j$dir, "killed.txt"))) return("stopped")      # stopped with the Stop button
    lg <- file.path(j$dir, "log.txt")
    # on Windows the running program keeps log.txt locked: then it is simply still running
    txt <- if (file.exists(lg)) tryCatch(suppressWarnings(readLines(lg, warn = FALSE)), error = function(e) character(0)) else character(0)
    if (any(grepl("cannot open|too short|STOP|Error|NaN|Segmentation|error", txt, ignore.case = TRUE))) {
      # a passing file problem (file locked by another program, e.g. the editor reading progress.txt or a
      # virus scanner): start it again, up to SPEEDY_RETRIES times
      if (any(grepl("Permission denied|Cannot open file|locked|being used", txt, ignore.case = TRUE))) {
        if (.n_retries(j$dir) < SPEEDY_RETRIES) return("stopped")
      }
      return("failed")
    }
    pf <- file.path(j$dir, c("progress.txt", "started.txt")); pf <- pf[file.exists(pf)]
    last <- suppressWarnings(max(file.info(pf)$mtime, na.rm = TRUE))         # (started.txt is renewed when a run is resumed)
    if (difftime(Sys.time(), last, units = "mins") > SPEEDY_STALL_MIN) {
      if (.n_resumes(j$dir) >= SPEEDY_MAX_RESUMES) return("failed")
      return("stopped")                                                       # restarted (resumed) by speedy_start()
    }
    return("running")
  }
  "missing"
}
#' fraction done of a run (0-1)
.job_progress <- function(j) {
  st <- .job_status(j)
  if (st == "done") return(1)
  p <- file.path(j$dir, "progress.txt")
  if (!file.exists(p)) return(0)
  v <- tryCatch(suppressWarnings(scan(p, quiet = TRUE)), error = function(e) numeric(0))   # may be being written
  if (length(v) < 2) return(0)
  dm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[min(max(v[2], 1), 12)]
  dd <- if (length(v) >= 3) (v[3] - 1) / dm else 0             # day of the month (newer SPEEDY writes it daily)
  min(1, max(0, ((v[1] - 1981) * 12 + v[2] - 1 + dd) / ((j$spinup + j$years) * 12)))
}

#' state of the runs for a map and experiment: "done", "running", "missing", "failed" or "stopped"
speedy_status <- function(map, x = NULL, cfg = SILUCARB) {
  jobs <- tryCatch(speedy_jobs(map, x, cfg), error = function(e) NULL)
  if (is.null(jobs)) return(c(orig = "missing", corr = "missing"))
  sapply(jobs, .job_status)
}
#' how many SPEEDY runs are going on now (all maps and experiments)
.n_running <- function(cfg = SILUCARB) {
  d <- list.dirs(file.path(cfg$dir_results, "speedy"), recursive = FALSE)
  sum(vapply(d, function(x) .job_status(list(dir = x)) == "running", TRUE))
}

#' start the runs that are missing (in the background; returns at once). Runs that need a run
#' that is not finished yet are started on a later call (run_speedy() or the editor do this).
speedy_start <- function(maps, x = NULL, cfg = SILUCARB, max_jobs = cfg$workers) {
  if (!is.list(maps) || !is.null(maps$file)) maps <- if (is.list(maps)) list(maps) else as.list(maps)
  all_jobs <- unlist(lapply(maps, function(a) speedy_jobs(a, x, cfg)), recursive = FALSE)
  running <- .n_running(cfg)
  started <- 0; forcing <- list()
  for (j in all_jobs) {
    if (!.job_status(j) %in% c("missing", "failed", "stopped")) next
    if (running + started >= max_jobs) break
    hf <- NULL
    if (!is.null(j$after)) {                          # needs the finished reference run
      ref <- speedy_dir(j$age, j$after, cfg)
      if (!file.exists(file.path(ref, "means.bin"))) next
      hf <- .read_means(ref)[, 48:1, , 6]             # its heat flux into the sea (north to south, like the inputs)
    }
    key <- format(j$age); if (is.null(forcing[[key]])) forcing[[key]] <- speedy_forcing(j$age, cfg)
    fo <- shift_forcing(forcing[[key]], j$shift)
    if (!is.null(j$shift) || j$sea_mode >= 3) {          # experiments: sea ice linked to land ice
      li <- land_ice_flags_forcing(j$age, j$shift, j$dem, cfg)
      fo$icec <- cap_perennial(fo$icec, speedy_lats, li$flags)
      j$x$ice_growth <- li$flags
    }
    speedy_launch(j, fo, hf, cfg)
    started <- started + 1
  }
  invisible(started)
}
#' prepare and start one SPEEDY run in the background (j: a job as in speedy_jobs(); j$geo can
#' give the geography, e.g. with an ice sheet)
speedy_launch <- function(j, forcing, hf = NULL, cfg = SILUCARB) {
  exe <- normalizePath(speedy_exe(cfg))
  if (.Platform$OS.type != "windows") Sys.chmod(exe, "755")
  if (dir.exists(j$dir) && .speedy_alive(j$dir)) { .speedy_kill(j$dir); Sys.sleep(2) }       # hung: no progress but still there
  lg <- file.path(j$dir, "log.txt")
  if (file.exists(lg) && any(grepl("Permission denied|Cannot open file|locked|being used",
                                   tryCatch(readLines(lg, warn = FALSE), error = function(e) ""), ignore.case = TRUE))) {   # a restart after a file problem: count it
    writeLines(as.character(.n_retries(j$dir) + 1L), .retry_file(j$dir))
  }
  # a run that was stopped (PC asleep or restarted, virus scanner, Stop button) continues from its last
  # restart file (at most SPEEDY_RESTART_DAYS model days are lost); otherwise it starts afresh
  resume <- dir.exists(j$dir) && .has_restart(j$dir) && file.exists(file.path(j$dir, "namelist.nml")) && !file.exists(file.path(j$dir, "means.bin")) &&
    !any(grepl("Permission denied|Cannot open file|locked|being used", tryCatch(readLines(lg, warn = FALSE), error = function(e) ""), ignore.case = TRUE))
  if (resume) {
    unlink(file.path(j$dir, c("killed.txt", "running.lock")))
    n <- .n_resumes(j$dir) + 1L; writeLines(as.character(n), file.path(j$dir, "resumes.txt"))
  } else {
    unlink(j$dir, recursive = TRUE)
    speedy_prepare(j$dir, j$age, j$dem, forcing, j$years, j$spinup, j$x, j$sea_mode, hf, geo = j$geo)
  }
  writeLines(format(Sys.time()), file.path(j$dir, "started.txt"))
  old <- setwd(j$dir); on.exit(setwd(old))
  # through a small launcher that records how the program ended (exit code) in exitcodes.txt
  if (.Platform$OS.type == "windows") {
    writeLines(c("@echo off", sprintf('"%s" >> log.txt 2>&1', gsub("/", "\\\\", exe)), ">> exitcodes.txt echo %date% %time% exit %errorlevel%"), "run.cmd")
    .spawn_detached("cmd.exe /c run.cmd", j$dir)
  } else {
    system2("sh", c("-c", shQuote(sprintf('"%s" >> log.txt 2>&1; rc=$?; echo "$(date "+%%F %%T") exit $rc" >> exitcodes.txt', exe))), wait = FALSE)
  }
  message(format(Sys.time()), " SPEEDY ", if (resume) sprintf("resumed (%d)", .n_resumes(j$dir)) else "started", ": ", j$age, " Ma, ", basename(j$dir))
  invisible(TRUE)
}

# ---- list of requested climate runs ------------------------------------------------------
.reg_file <- function(cfg) file.path(cfg$dir_results, "speedy", "climate_runs.csv")
.REG_COLS <- c("id", "created", "map_name", "age", "map_file", "map_tag", "ecc", "obl", "peri", "co2_pal", "co2",
               "sun_age", "sensitivity", "length", "label", "ice_id", "gradient")
#' remember a requested run (map + experiment [+ ice sheets]), so it is continued and can be
#' plotted later; ice: ice_settings() or NULL
speedy_register <- function(map, x = NULL, cfg = SILUCARB, ice = NULL) {
  m <- climate_map(map, cfg = cfg); x <- x %||% climate_experiment(); r <- .resolve(x, m$age)
  id <- paste0(m$tag, "_", if (.is_reference(x, m$age)) "reference" else .exp_tag(r), "_", x$length)
  ice_id <- NA
  if (!is.null(ice)) {
    ice_id <- paste0(id, "_", .ice_tag(ice)); id <- ice_id
    dir.create(ice_dir(ice_id, cfg), showWarnings = FALSE, recursive = TRUE)
    saveRDS(ice, file.path(ice_dir(ice_id, cfg), "settings.rds"))
  }
  reg <- climate_runs(cfg, status = FALSE)
  if (id %in% reg$id) return(invisible(id))
  # keep a copy of the map as it was run (the editor's grid changes when an age is rebuilt)
  keep <- file.path(cfg$dir_results, "speedy", "maps", paste0(m$tag, ".nc"))
  if (!file.exists(keep)) { dir.create(dirname(keep), showWarnings = FALSE, recursive = TRUE); file.copy(m$file, keep) }
  lab <- if (.is_reference(x, m$age)) "reference (HadCM3L orbit, CO2 and Sun)" else r$label
  if (!is.null(ice)) lab <- sprintf("%s + ice sheets (%g kyr%s%s%s)", lab, ice$duration_kyr,
                                    if (ice$couple_kyr > 0) sprintf(", climate updated every %g kyr", ice$couple_kyr) else ", no climate feedback",
                                    switch(ice$init, cap = sprintf(", start: %g m ice poleward of %g\u00b0", ice$init_thick, ice$init_lat),
                                           run = ", start: end of another ice run", ""),
                                    sprintf(", %s bed%s", ice$bed %||% "mixed", if (isTRUE(ice$sea_level)) ", sea level follows the ice" else ""))
  row <- data.frame(id = id, created = format(Sys.time(), "%Y-%m-%d %H:%M"), map_name = m$name, age = m$age,
                    map_file = normalizePath(keep, winslash = "/"), map_tag = m$tag, ecc = x$ecc, obl = x$obl, peri = x$peri, co2_pal = x$co2_pal %||% NA,
                    co2 = if (is.null(x$co2_pal)) x$co2 %||% 1 else NA, sun_age = x$sun_age %||% NA,
                    sensitivity = x$sensitivity, length = x$length, label = lab, ice_id = ice_id,
                    gradient = x$gradient %||% 1)
  dir.create(dirname(.reg_file(cfg)), showWarnings = FALSE, recursive = TRUE)
  reg <- if (nrow(reg)) { if (is.null(reg$ice_id)) reg$ice_id <- NA; if (is.null(reg$gradient)) reg$gradient <- 1
                          rbind(reg[, names(row)], row) } else row
  utils::write.csv(reg, .reg_file(cfg), row.names = FALSE)
  invisible(id)
}
#' the map and experiment of a registered run (a row of climate_runs())
run_setup <- function(row, cfg = SILUCARB) {
  row <- as.list(row)
  m <- list(age = as.numeric(row$age), file = row$map_file, name = row$map_name, tag = row$map_tag)
  na <- function(v) if (is.null(v) || is.na(v)) NULL else as.numeric(v)
  x <- climate_experiment(list(ecc = as.numeric(row$ecc), obl = as.numeric(row$obl), peri = as.numeric(row$peri)),
                          co2 = na(row$co2) %||% 1, co2_pal = na(row$co2_pal), sun_age = na(row$sun_age),
                          sensitivity = as.numeric(row$sensitivity), length = row$length,
                          gradient = na(row$gradient) %||% 1)
  list(map = m, x = x)
}
#' all requested climate runs with their state and progress
climate_runs <- function(cfg = SILUCARB, status = TRUE) {
  f <- .reg_file(cfg)
  reg <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else
    as.data.frame(setNames(replicate(length(.REG_COLS), character(0), simplify = FALSE), .REG_COLS))
  if (is.null(reg$ice_id)) reg$ice_id <- rep(NA_character_, nrow(reg))
  if (is.null(reg$gradient)) reg$gradient <- rep(1, nrow(reg))
  reg$ice_id[!is.na(reg$ice_id) & !nzchar(reg$ice_id)] <- NA
  if (!status || !nrow(reg)) return(reg)
  st <- lapply(seq_len(nrow(reg)), function(i) {
    s <- run_setup(reg[i, ], cfg)
    if (!file.exists(s$map$file)) return(c(status = "map file missing", progress = 0))
    jobs <- speedy_jobs(s$map, s$x, cfg); js <- vapply(jobs, .job_status, "")
    w <- vapply(jobs, function(j) (j$spinup + j$years), 0); p <- vapply(jobs, .job_progress, 0)
    out <- c(status = if (all(js == "done")) "done" else if (any(js == "failed")) "failed" else
                 if (any(js %in% c("running"))) "running" else if (any(js == "stopped")) "stopped" else "waiting",
             progress = round(100 * sum(w * p) / sum(w)))
    if (!is.na(reg$ice_id[i])) {
      is <- ice_status(reg$ice_id[i], cfg)
      out <- if (out[["status"]] != "done") c(status = paste("climate", out[["status"]]), progress = round(as.numeric(out[["progress"]]) / 2))
             else c(status = if (is$status == "running") sprintf("ice %.1f kyr%s", is$t_kyr, if (!is.null(is$gradient)) sprintf(" (gradient %.2f)", is$gradient) else "") else if (is$status == "waiting") "ice waiting" else is$status,
                    progress = 50 + round(is$progress / 2))
    }
    out
  })
  reg$status <- vapply(st, `[[`, "", "status"); reg$progress <- as.numeric(vapply(st, `[[`, "", "progress"))
  reg
}
#' start the next steps of all requested runs that are not finished (called regularly by the
#' background runner; after a restart of the PC, stopped runs are started again)
.priority_file <- function(cfg) file.path(cfg$dir_results, "speedy", "priority.txt")
#' let these runs (ids) go first; the runner works through the registered runs in this order, then
#' through the others in the order they were registered
speedy_priority_set <- function(ids, cfg = SILUCARB) {
  f <- .priority_file(cfg); old <- if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
  writeLines(unique(c(ids, old)), f); invisible(ids)
}
#' which SPEEDY runs are going on now (folder names)
speedy_running_now <- function(cfg = SILUCARB) {
  d <- list.dirs(file.path(cfg$dir_results, "speedy"), recursive = FALSE)
  basename(d[vapply(d, function(x) .job_status(list(dir = x)) == "running", TRUE)])
}
#' restart the background runner only (e.g. to load updated code or a new queue order); SPEEDY
#' programs that are running go on undisturbed
#' end every background runner (R programs running results/speedy/runner.R), also ones that no longer
#' own the lock file - so that never two runners go on at the same time. SPEEDY programs go on.
.kill_runners <- function(cfg = SILUCARB) {
  lk <- .runner_lock(cfg)
  pid <- if (file.exists(lk)) suppressWarnings(as.integer(tryCatch(readLines(lk, warn = FALSE), error = function(e) "")[2])) else NA
  if (isTRUE(pid > 0) && pid != Sys.getpid()) try(tools::pskill(pid), silent = TRUE)
  if (.Platform$OS.type == "windows") {
    ps <- sprintf("Get-CimInstance Win32_Process | Where-Object { $_.Name -match '^(Rscript|R|Rterm)\\.exe$' -and $_.CommandLine -like '*speedy*runner.R*' -and $_.ProcessId -ne %d } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force }", Sys.getpid())
    try(suppressWarnings(system2("powershell.exe", c("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", shQuote(ps)), stdout = FALSE, stderr = FALSE)), silent = TRUE)
  } else try(suppressWarnings(system2("pkill", c("-f", shQuote(file.path(cfg$dir_results, "speedy", "runner.R"))), stdout = FALSE, stderr = FALSE)), silent = TRUE)
  unlink(lk); invisible(TRUE)
}
speedy_runner_restart <- function(cfg = SILUCARB, project_dir = getwd()) {
  .kill_runners(cfg); Sys.sleep(1)
  speedy_runner_start(cfg, project_dir, resume = TRUE)
}
speedy_continue <- function(cfg = SILUCARB) {
  reg <- climate_runs(cfg, status = FALSE); n <- 0; ice <- integer(0)
  pf <- .priority_file(cfg)
  if (file.exists(pf) && nrow(reg)) { pr <- readLines(pf, warn = FALSE); o <- match(reg$id, pr); reg <- reg[order(is.na(o), o, seq_len(nrow(reg))), , drop = FALSE] }
  # first start / resume the SPEEDY runs of every registered run, then do ice-sheet steps (which can
  # take long) - so that waiting climate runs are never held up by an ice sheet
  for (i in seq_len(nrow(reg))) {
    s <- run_setup(reg[i, ], cfg)
    if (!file.exists(s$map$file)) next
    st <- speedy_status(s$map, s$x, cfg)
    if (any(st == "failed")) next
    if (!all(st == "done")) { n <- n + speedy_start(list(s$map), s$x, cfg); next }
    if (!is.na(reg$ice_id[i])) ice <- c(ice, i)
  }
  for (i in ice) tryCatch(ice_continue_one(reg[i, ], cfg), error = function(e) message(format(Sys.time()), " ice: ", conditionMessage(e)))
  invisible(n)
}

#' run SPEEDY for these maps (ages, or climate_map() / import_map() maps) and an experiment,
#' and wait until everything is finished; several runs go in parallel (workers).
#' Example: run_speedy(425, climate_experiment(ecc = 0.05, obl = 24.5, peri = 90, co2_pal = 8, sun_age = 425))
run_speedy <- function(maps, x = NULL, cfg = SILUCARB, workers = cfg$workers, ice = NULL) {
  if (!is.list(maps) || !is.null(maps$file)) maps <- if (is.list(maps)) list(maps) else as.list(maps)
  maps <- lapply(maps, climate_map, cfg = cfg)
  ids <- vapply(maps, function(m) speedy_register(m, x, cfg, ice), "")
  repeat {
    jobs <- unlist(lapply(maps, function(a) speedy_jobs(a, x, cfg)), recursive = FALSE)
    st <- vapply(jobs, .job_status, "")
    if (all(st == "done")) break
    if (any(st == "failed")) stop("a SPEEDY run failed; see log.txt in ", file.path(cfg$dir_results, "speedy"), call. = FALSE)
    speedy_start(maps, x, cfg, workers)
    message(format(Sys.time(), "%H:%M"), "  SPEEDY: ", sum(st == "done"), " of ", length(st), " runs done")
    Sys.sleep(30)
  }
  message("SPEEDY runs finished: ", paste(vapply(maps, `[[`, "", "name"), collapse = ", "))
  if (!is.null(ice)) {                                 # ice sheets (with the coupling SPEEDY runs)
    cfg$workers <- workers
    repeat {
      r <- climate_runs(cfg); r <- r[r$id %in% ids, ]
      if (all(r$status %in% c("done", "failed"))) break
      for (i in seq_len(nrow(r))) ice_continue_one(r[i, ], cfg)
      Sys.sleep(5)
    }
    message("ice-sheet runs finished")
  }
  invisible(TRUE)
}

#' the monthly fields SPEEDY writes (means.bin; runs made before this version have only 1-7)
SPEEDY_VARS <- data.frame(
  name = c("precip", "u0", "v0", "t0", "evap", "hflux_sea", "olr", "precip_conv", "precip_ls", "snowfall",
           "ts", "sst", "sea_ice", "cloud", "cloud_strat", "ps", "q_sfc", "rh_sfc",
           "u500", "v500", "t500", "z500", "u200", "v200", "t200", "z200",
           "sw_toa_in", "sw_toa_net", "sw_sfc_down", "sw_sfc_net", "lw_sfc_down", "lw_sfc_net_up",
           "shf", "lhf", "taux", "tauy", "t0_sd", "pdd"),
  units = c("mm/day", "m/s", "m/s", "K", "mm/day", "W/m2", "W/m2", "mm/day", "mm/day", "mm/day (water)",
            "K", "K", "fraction", "fraction", "fraction", "hPa", "g/kg", "%",
            "m/s", "m/s", "K", "m", "m/s", "m/s", "K", "m",
            "W/m2", "W/m2", "W/m2", "W/m2", "W/m2", "W/m2", "W/m2", "W/m2", "N/m2", "N/m2", "K", "K"),
  long_name = c("precipitation (total)", "eastward wind, lowest level (sigma 0.95, ~950 hPa)",
                "northward wind, lowest level (sigma 0.95)", "near-surface air temperature",
                "evaporation", "net heat flux into the sea (incl. sea ice; daily)", "outgoing longwave radiation (top)",
                "convective precipitation", "large-scale precipitation", "snowfall (precipitation below -2 to +2 C)",
                "surface temperature (land skin / sea surface)", "sea-surface temperature seen by the atmosphere (incl. ice)",
                "sea-ice fraction", "total cloud cover", "stratiform cloud cover", "surface pressure",
                "specific humidity, lowest level", "relative humidity, lowest level",
                "eastward wind at sigma 0.51 (~500 hPa)", "northward wind at sigma 0.51 (~500 hPa)",
                "temperature at sigma 0.51 (~500 hPa)", "geopotential height at sigma 0.51 (~500 hPa)",
                "eastward wind at sigma 0.20 (~200 hPa)", "northward wind at sigma 0.20 (~200 hPa)",
                "temperature at sigma 0.20 (~200 hPa)", "geopotential height at sigma 0.20 (~200 hPa)",
                "incoming shortwave radiation (top)", "absorbed shortwave radiation (top, net down)",
                "downward shortwave radiation at the surface", "net shortwave radiation at the surface (down)",
                "downward longwave radiation at the surface", "net longwave radiation at the surface (up)",
                "sensible heat flux (up)", "latent heat flux (up)", "eastward surface wind stress",
                "northward surface wind stress", "standard deviation of daily-step near-surface temperature",
                "mean positive near-surface temperature (degree-days per day)"),
  stringsAsFactors = FALSE)

#' monthly means of a finished run: array (96 lon from 0 E, 48 lat south to north, 12, n variables)
#' with the variable names as the 4th dimnames (1 precipitation, 2-3 wind u, v, 4 near-surface
#' temperature, 5 evaporation, 6 heat flux into the sea, 7 OLR, ... see SPEEDY_VARS)
.read_means <- function(dir) {
  f <- file.path(dir, "means.bin"); nv <- file.size(f) / (96 * 48 * 12 * 4)
  x <- readBin(f, "numeric", size = 4, n = 96 * 48 * 12 * nv, endian = "little")
  a <- array(x, c(96, 48, 12, nv)); a[a < -9000] <- NA
  vf <- file.path(dir, "means_vars.txt")
  nm <- if (file.exists(vf)) trimws(readLines(vf, warn = FALSE)) else SPEEDY_VARS$name[seq_len(nv)]
  dimnames(a) <- list(NULL, NULL, month.abb, nm[seq_len(nv)])
  a
}

#' difference b - a of two runs as fields (precip, wind_u, wind_v, t0) for a season
.run_diff <- function(da, db, season) {
  mo <- season_months(season)
  a <- .read_means(da); b <- .read_means(db)
  lat <- rev(speedy_lats); lon <- ((SPEEDY_LON + 180) %% 360) - 180; p <- order(lon)
  fld <- function(k) { d <- apply(b[, , mo, k, drop = FALSE] - a[, , mo, k, drop = FALSE], c(1, 2), mean)
    d[!is.finite(d)] <- 0
    z <- t(d)[, p]
    # light smoothing (~ one grid cell) against the run-to-run noise of short runs
    list(z = smooth_km(z, lat, 250, nbands = 24), lat = lat, lon = lon[p]) }
  list(precip = fld(1), wind_u = fld(2), wind_v = fld(3), t0 = fld(4), evap = fld(5))
}
.DIFF_VARS <- c("precip", "wind_u", "wind_v", "t0", "evap")
.add_diff <- function(a, b, sign = 1) { if (is.null(b)) return(a); if (is.null(a)) { if (sign == 1) return(b); a <- b; for (k in .DIFF_VARS) a[[k]]$z <- 0 * a[[k]]$z }
  for (k in .DIFF_VARS) a[[k]]$z <- a[[k]]$z + sign * b[[k]]$z; a }

#' change caused by an experiment alone, on the map: (forcing run - corr) + (orbit run - slab run)
speedy_experiment_change <- function(map, season, x, cfg = SILUCARB) {
  jobs <- speedy_jobs(map, x, cfg)
  if (length(jobs) < 3 || any(vapply(jobs, .job_status, "") != "done")) return(NULL)
  d <- NULL
  if (!is.null(jobs$forcing)) d <- .run_diff(jobs$corr$dir, jobs$forcing$dir, season)
  if (!is.null(jobs$exp)) d <- .add_diff(d, .run_diff(jobs$slab$dir, jobs$exp$dir, season))
  d$label <- .resolve(x, climate_map(map, cfg = cfg)$age)$label; d
}

#' SPEEDY response for the map: (your map - original map) [+ experiment change];
#' NULL when the runs are not there (yet)
speedy_delta <- function(map, season, cfg = SILUCARB, x = NULL) {
  jobs <- tryCatch(speedy_jobs(map, x, cfg), error = function(e) NULL)
  if (is.null(jobs)) return(NULL)
  if (any(vapply(jobs, .job_status, "") != "done")) return(NULL)
  d <- .run_diff(jobs$orig$dir, jobs$corr$dir, season)
  if (length(jobs) >= 3) d <- .add_diff(d, speedy_experiment_change(map, season, x, cfg))
  d$years <- jobs[[length(jobs)]]$years
  d
}

# ---- background runner -------------------------------------------------------------------
#' a small R process that keeps the requested climate runs going (starts each next step when
#' the run before it finishes) until all are done, also when the editor and RStudio are closed.
#' Only one runner works at a time (results/speedy/runner.lock).
.runner_lock <- function(cfg) file.path(cfg$dir_results, "speedy", "runner.lock")
RUNNER_STALE_MIN <- 15    # a runner whose lock was not refreshed for this long is taken as gone
speedy_runner_active <- function(cfg = SILUCARB) {
  f <- .runner_lock(cfg)
  file.exists(f) && difftime(Sys.time(), file.info(f)$mtime, units = "mins") < RUNNER_STALE_MIN
}
#' keep the runner's lock fresh during long steps (ice sheets); called every model kyr
.heartbeat <- function() { h <- getOption("silucarb.heartbeat"); if (is.function(h)) try(h(), silent = TRUE); invisible() }
.stop_flag <- function(cfg) file.path(cfg$dir_results, "speedy", "STOPPED")
#' start the background runner (not when the runs were stopped with speedy_stop(), unless resume = TRUE)
.winpath <- function(p) gsub("/", "\\\\", normalizePath(p, winslash = "/", mustWork = FALSE))
#' start a program on Windows in its own hidden console (through Windows Script Host), so that a
#' Ctrl+C or Stop in the R console / RStudio that started it (which Windows passes on to every program
#' sharing that console - exit code -1073741510 = 0xC000013A) does not end SPEEDY or the runner
.spawn_detached <- function(cmdline, dir) {
  vbs <- file.path(dir, "launch.vbs")
  writeLines(c('Set sh = CreateObject("WScript.Shell")', sprintf('sh.CurrentDirectory = "%s"', .winpath(dir)),
               sprintf('sh.Run "%s", 0, False', gsub('"', '""', cmdline, fixed = TRUE))), vbs)
  wscript <- file.path(Sys.getenv("SystemRoot", "C:/Windows"), "System32", "wscript.exe")
  # wait = FALSE: with wait = TRUE R on Windows waits until every program started from it has ended,
  # i.e. until SPEEDY is finished - so only one run went on at a time
  ok <- file.exists(wscript) && !inherits(tryCatch(system2(wscript, c("//B", "//Nologo", shQuote(.winpath(vbs))), wait = FALSE), error = function(e) e), "error")
  if (!ok) { old <- setwd(dir); on.exit(setwd(old)); system(cmdline, wait = FALSE, invisible = TRUE) }   # fall-back: the old way
  invisible(ok)
}

speedy_runner_start <- function(cfg = SILUCARB, project_dir = getwd(), resume = FALSE) {
  if (resume) unlink(.stop_flag(cfg))
  if (file.exists(.stop_flag(cfg))) return(invisible(FALSE))
  if (speedy_runner_active(cfg)) return(invisible(FALSE))
  .kill_runners(cfg)                                   # a runner that hangs (stale lock) must not go on beside the new one
  dir.create(dirname(.runner_lock(cfg)), showWarnings = FALSE, recursive = TRUE)
  writeLines(format(Sys.time()), .runner_lock(cfg))
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  scr <- file.path(cfg$dir_results, "speedy", "runner.R")
  writeLines(c(sprintf("setwd('%s')", normalizePath(project_dir, winslash = "/")),
               "for (f in c('config', 'grid_tools', 'climate', 'climate_correct', 'speedy', 'climate_export', 'ice_sheet', 'rivers', 'ensemble', 'orbital_ice', 'ocean_state', 'baltic_sill', 'baltic_experiment', 'speedy_guard'))",
               "  source(file.path('R', paste0(f, '.R')), encoding = 'UTF-8')",
               "SILUCARB <- apply_user_settings(absolutize_config(SILUCARB, getwd()))",
               "# a Ctrl+C / Stop in the R console that started it must not end the runner",
               "repeat { done <- tryCatch({ speedy_runner_loop(SILUCARB); TRUE }, interrupt = function(i) { message(format(Sys.time()), ' interrupt ignored'); FALSE }); if (done) break }"), scr)
  lg <- file.path(cfg$dir_results, "speedy", "runner_log.txt")
  if (.Platform$OS.type == "windows") .spawn_detached(sprintf('cmd.exe /c ""%s" "%s" > "%s" 2>&1"', .winpath(rs), .winpath(scr), .winpath(lg)), dirname(scr))
  else system2(rs, shQuote(scr), wait = FALSE, stdout = lg, stderr = lg)
  invisible(TRUE)
}
speedy_runner_loop <- function(cfg = SILUCARB, every = 30) {
  if (exists("speedy_keep_awake")) try(speedy_keep_awake(cfg), silent = TRUE)        # no sleep while runs go on (Windows)
  options(silucarb.heartbeat = function() writeLines(c(format(Sys.time()), Sys.getpid()), .runner_lock(cfg)))
  repeat {
    if (file.exists(.stop_flag(cfg))) { message(format(Sys.time()), " stopped by the user"); break }
    writeLines(c(format(Sys.time()), Sys.getpid()), .runner_lock(cfg))
    cfg$workers <- tryCatch(apply_user_settings(cfg)$workers, error = function(e) cfg$workers)   # cores can be changed while it runs
    tryCatch(speedy_continue(cfg), error = function(e) message(format(Sys.time()), " ", conditionMessage(e)))
    # climate matrices and astronomically paced (orbital) ice runs
    pending <- exists("orbital_continue") && isTRUE(tryCatch(orbital_continue(cfg), error = function(e) { message(format(Sys.time()), " orbital: ", conditionMessage(e)); TRUE }))
    r <- tryCatch(climate_runs(cfg), error = function(e) NULL)
    if (is.null(r)) { Sys.sleep(every); next }
    if (!pending && (!nrow(r) || all(r$status %in% c("done", "failed", "climate failed", "map file missing")))) break
    Sys.sleep(every)
  }
  unlink(.runner_lock(cfg))
  message(format(Sys.time()), " all climate runs finished")
}

#' stop all climate runs: the background runner (and an ice-sheet step it is doing) and every
#' SPEEDY program; the runs keep their state and go on with speedy_runner_start(resume = TRUE)
speedy_stop <- function(cfg = SILUCARB) {
  dir.create(dirname(.stop_flag(cfg)), showWarnings = FALSE, recursive = TRUE)
  writeLines(format(Sys.time()), .stop_flag(cfg))
  lk <- .runner_lock(cfg)
  pid <- if (file.exists(lk)) suppressWarnings(as.integer(tryCatch(readLines(lk, warn = FALSE), error = function(e) "")[2])) else NA
  if (isTRUE(pid > 0) && pid != Sys.getpid()) tools::pskill(pid)
  .kill_runners(cfg); unlink(lk)
  unlink(list.files(file.path(cfg$dir_results, "speedy", "ice"), "^busy\\.lock$", recursive = TRUE, full.names = TRUE))
  unlink(list.files(file.path(cfg$dir_results, "speedy", "orbital_ice"), "^busy\\.lock$", recursive = TRUE, full.names = TRUE))
  d <- list.dirs(file.path(cfg$dir_results, "speedy"), recursive = FALSE)
  run <- d[vapply(d, function(x) .job_status(list(dir = x)) == "running", TRUE)]
  exe <- tryCatch(normalizePath(speedy_exe(cfg)), error = function(e) NULL)
  if (.Platform$OS.type == "windows") suppressWarnings(system2("taskkill", c("/F", "/IM", "speedy.exe"), stdout = FALSE, stderr = FALSE))
  else if (!is.null(exe)) suppressWarnings(system2("pkill", c("-f", shQuote(exe)), stdout = FALSE, stderr = FALSE))
  for (x in run) writeLines(format(Sys.time()), file.path(x, "killed.txt"))
  invisible(length(run))
}

#' what the climate runs are doing now: state "active", "halted", "stopped" (by the user) or
#' "idle" (nothing to do), with the number of SPEEDY programs running and the last progress
speedy_activity <- function(cfg = SILUCARB, reg = NULL) {
  reg <- reg %||% tryCatch(climate_runs(cfg), error = function(e) NULL)
  todo <- if (is.null(reg) || !nrow(reg)) 0 else sum(!reg$status %in% c("done", "done (equilibrium)", "failed", "climate failed", "map file missing"))
  d <- list.dirs(file.path(cfg$dir_results, "speedy"), recursive = FALSE)
  pf <- c(file.path(d, "progress.txt"), list.files(file.path(cfg$dir_results, "speedy", c("ice", "orbital_ice")), "timeseries.csv|state.rds", recursive = TRUE, full.names = TRUE))
  mt <- file.info(pf[file.exists(pf)])$mtime
  last <- if (length(mt)) max(mt) else NA
  age_min <- if (is.na(last)) Inf else as.numeric(difftime(Sys.time(), last, units = "mins"))
  n_prog <- sum(vapply(d, function(x) { p <- file.path(x, "progress.txt")
    file.exists(p) && !file.exists(file.path(x, "means.bin")) && !file.exists(file.path(x, "killed.txt")) &&
      difftime(Sys.time(), file.info(p)$mtime, units = "mins") < 5 }, TRUE))
  if (exists("orbital_ice_list")) todo <- todo + sum(vapply(orbital_ice_list(cfg), function(nm)
    !identical(tryCatch(orbital_ice_status(nm, cfg)$status, error = function(e) "done"), "done"), TRUE))
  n_ice <- sum(vapply(c(list.dirs(file.path(cfg$dir_results, "speedy", "ice"), recursive = FALSE),
                        list.dirs(file.path(cfg$dir_results, "speedy", "orbital_ice"), recursive = FALSE)), function(x) {
    f <- file.path(x, "busy.lock"); file.exists(f) && difftime(Sys.time(), file.info(f)$mtime, units = "mins") < RUNNER_STALE_MIN }, TRUE))
  runner <- speedy_runner_active(cfg); stopped <- file.exists(.stop_flag(cfg))
  state <- if (stopped && todo > 0) "stopped" else if (todo == 0 && n_prog == 0 && n_ice == 0) "idle" else
           if (n_prog > 0 || n_ice > 0 || (runner && age_min < SPEEDY_STALL_MIN)) "active" else "halted"
  list(state = state, runner = runner, n_running = n_prog, n_ice = n_ice, unfinished = todo, last = last, last_min = age_min, stopped = stopped)
}
activity_text <- function(a) {
  when <- if (is.na(a$last)) "no progress yet" else sprintf("last progress %s (%s ago)", format(a$last, "%H:%M:%S"),
            if (a$last_min < 1) "<1 min" else if (a$last_min < 120) sprintf("%.0f min", a$last_min) else sprintf("%.1f h", a$last_min / 60))
  switch(a$state,
    active = sprintf("ACTIVE \u2013 %d SPEEDY run(s) and %d ice-sheet run(s) computing, background runner %s; %d run(s) unfinished; %s.",
                     a$n_running, a$n_ice %||% 0, if (a$runner) "on" else "off", a$unfinished, when),
    halted = sprintf("HALTED \u2013 %d unfinished run(s) but nothing is computing (background runner %s); %s. Press Resume.",
                     a$unfinished, if (a$runner) "on but waiting" else "off", when),
    stopped = sprintf("STOPPED by you \u2013 %d unfinished run(s) kept; press Resume to go on; %s.", a$unfinished, when),
    idle = sprintf("IDLE \u2013 all runs finished; %s.", when))
}
