# =============================================================================
#  Climate fields fitted to YOUR corrected geography (land added / removed).
#
#  The HadCM3L runs (R/climate.R) were made on the original PALEOMAP geography. Three ways
#  to use them on the corrected map:
#   "model"  the model fields as they are (only interpolated);
#   "local"  local correction learned from the same model run: per latitude band, how
#            precipitation and near-surface wind speed depend on the local land fraction
#            (150 and 500 km), the sea fraction upwind (300-1000 km), and up- / down-slope
#            flow over the relief. The field is then shifted by
#            prediction(corrected geography) - prediction(original geography).
#            Captures continentality, moisture supply from new seaways, rain-shadows and
#            slower wind over new land; NOT changes of the large-scale circulation.
#   "rerun"  the simplified climate model SPEEDY (R/speedy.R) is run twice, on the original
#            and on the corrected geography, with the same sea-surface temperatures from the
#            HadCM3L run; their difference is added to the HadCM3L field.
#  Flag: where the corrected land fraction differs much from the model's (> 25 % of the area
#  within ~150 km), the climate fields are least reliable; these areas are outlined.
# =============================================================================

GRID1 <- list(lat = seq(-89.5, 89.5, 1), lon = seq(-179.5, 179.5, 1))

#' land fraction and mean land height (m) on the 1-degree grid
geo_1deg <- function(dem) {
  ii <- pmin(pmax(floor(dem$lat + 90) + 1, 1), 180)
  jj <- pmin(pmax(floor(((dem$lon + 180) %% 360)) + 1, 1), 360)
  agg <- function(m) { a <- rowsum(m, ii, reorder = TRUE); t(rowsum(t(a), jj, reorder = TRUE)) }
  n <- agg(matrix(1, nrow(dem$z), ncol(dem$z)))
  list(lf = agg((dem$z >= 0) * 1) / n, h = agg(pmax(dem$z, 0)) / n, lat = GRID1$lat, lon = GRID1$lon)
}
.G1 <- new.env()
geo_1deg_cached <- function(file) {
  inf <- file.info(file); key <- paste(file, inf$mtime, inf$size)
  if (is.null(.G1[[key]])) .G1[[key]] <- geo_1deg(read_dem_nc(file))
  .G1[[key]]
}
.fld <- function(z) list(z = z, lat = GRID1$lat, lon = GRID1$lon)

#' predictors on the 1-degree grid; wind (u, v in m/s on the same grid) sets "upwind" and slope flow
climate_predictors <- function(g, wu, wv) {
  lat <- GRID1$lat
  L150 <- smooth_km(g$lf, lat, 150); L500 <- smooth_km(g$lf, lat, 500)
  sp <- sqrt(wu^2 + wv^2) + 1e-6; ue <- wu / sp; vn <- wv / sp
  LAT <- matrix(lat, 180, 360); LON <- matrix(GRID1$lon, 180, 360, byrow = TRUE)
  up <- 0
  for (d in c(300, 600, 1000)) {                    # sea fraction upwind
    la <- pmin(pmax(LAT - vn * d / KM_PER_DEG, -89.5), 89.5)
    lo <- LON - ue * d / (KM_PER_DEG * pmax(cos(LAT * pi / 180), 0.1))
    up <- up + (1 - clim_interp(.fld(L150), as.vector(la), as.vector(lo))) / 3
  }
  hs <- smooth_km(g$h, lat, 100)
  dy <- KM_PER_DEG
  gy <- (rbind(hs[-1, ], hs[180, ]) - rbind(hs[1, ], hs[-180, ])) / (2 * dy)            # m per km
  gx <- (cbind(hs[, -1], hs[, 1]) - cbind(hs[, 360], hs[, -360])) / (2 * dy * pmax(cos(LAT * pi / 180), 0.1))
  flow <- (wu * gx + wv * gy) / 10                  # ~ vertical motion forced by the relief
  list(L150 = L150, L500 = L500, upwind_sea = matrix(up, 180, 360),
       upslope = pmax(flow, 0), downslope = pmax(-flow, 0), height = hs / 1000)
}

#' least-squares fits per latitude band (weights fall off over 12 deg of latitude)
.band_fit <- function(y, X, lat, centers = seq(-80, 80, 10), ridge = 0.05) {
  lat <- as.vector(lat); y <- as.vector(y)
  X <- cbind(1, X); p <- ncol(X)
  sapply(centers, function(c0) {
    w <- exp(-0.5 * ((lat - c0) / 12)^2)
    A <- crossprod(X * w, X) + diag(c(0, rep(ridge * sum(w), p - 1)))
    solve(A, crossprod(X * w, y))
  })                                                  # p x n_centers
}
.band_coef <- function(B, lat, centers = seq(-80, 80, 10))
  t(apply(B, 1, function(b) approx(centers, b, xout = pmin(pmax(lat, -80), 80))$y))   # p x length(lat)

#' fields on the 1-degree grid from a model field
.to_1deg <- function(F) {
  g <- expand.grid(lon = GRID1$lon, lat = GRID1$lat)
  matrix(clim_interp(F, g$lat, g$lon), 180, 360, byrow = TRUE)
}

#' the local statistical correction; returns precip, wind_u, wind_v on the 1-degree grid and
#' the fit quality (R2) for precipitation and wind speed
local_correction <- function(cl, g_orig, g_corr) {
  P1 <- pmax(.to_1deg(cl$precip), 0); U1 <- .to_1deg(cl$wind_u); V1 <- .to_1deg(cl$wind_v)
  U1[!is.finite(U1)] <- 0; V1[!is.finite(V1)] <- 0
  Xo <- climate_predictors(g_orig, U1, V1); Xc <- climate_predictors(g_corr, U1, V1)
  vars_p <- c("L150", "L500", "upwind_sea", "upslope", "downslope", "height")
  vars_w <- c("L150", "L500", "height")
  # training: the model's own cells (precipitation grid), predictors of the original geography
  la <- cl$precip$lat; lo <- cl$precip$lon
  tg <- expand.grid(lon = lo, lat = la); ok <- abs(tg$lat) < 88
  samp <- function(M) clim_interp(.fld(M), tg$lat[ok], tg$lon[ok])
  Xt <- sapply(vars_p, function(v) samp(Xo[[v]]))
  yp <- log(pmax(as.vector(t(cl$precip$z))[ok], 0) + 0.3)
  Bp <- .band_fit(yp, Xt, tg$lat[ok])
  wt <- sqrt(clim_interp(cl$wind_u, tg$lat[ok], tg$lon[ok])^2 + clim_interp(cl$wind_v, tg$lat[ok], tg$lon[ok])^2)
  okw <- is.finite(wt)
  Bw <- .band_fit(log(wt[okw] + 0.5), Xt[okw, vars_w, drop = FALSE], tg$lat[ok][okw])
  r2 <- function(y, X, B, lat) { cf <- .band_coef(B, lat); f <- colSums(t(cbind(1, X)) * cf); 1 - sum((y - f)^2) / sum((y - mean(y))^2) }
  # apply: change of the prediction between the two geographies
  cfp <- .band_coef(Bp, GRID1$lat); cfw <- .band_coef(Bw, GRID1$lat)
  dP <- 0; for (k in seq_along(vars_p)) dP <- dP + cfp[k + 1, ] * (Xc[[vars_p[k]]] - Xo[[vars_p[k]]])
  dW <- 0; for (k in seq_along(vars_w)) dW <- dW + cfw[k + 1, ] * (Xc[[vars_w[k]]] - Xo[[vars_w[k]]])
  dP <- pmax(pmin(dP, 2), -2); dW <- pmax(pmin(dW, 1), -1)          # at most a factor e^2 resp. e
  list(precip = .fld(pmax((P1 + 0.3) * exp(dP) - 0.3, 0)),
       wind_u = .fld(U1 * exp(dW)), wind_v = .fld(V1 * exp(dW)),
       r2 = c(precip = r2(yp, Xt, Bp, tg$lat[ok]), wind = r2(log(wt[okw] + 0.5), Xt[okw, vars_w, drop = FALSE], Bw, tg$lat[ok][okw])))
}

#' areas where the corrected land differs much from the model's: 1-degree mask (0/1)
geography_flag <- function(g_orig, g_corr, threshold = 0.25) {
  d <- abs(smooth_km(g_corr$lf, GRID1$lat, 150) - smooth_km(g_orig$lf, GRID1$lat, 150))
  (d > threshold) * 1
}

#' the climate for one map: the HadCM3L fields, optionally fitted to the corrected geography
#' map:    an age (the editor's corrected grid of that age) or climate_map() / import_map()
#' method: "model", "local" or "rerun" (falls back to "local" if no SPEEDY run is there)
#' experiment: climate_experiment() (orbit / CO2 / Sun; needs the SPEEDY runs)
#' compare:    another experiment: the change layers then show experiment - compare
climate_for_map <- function(map, season = "ann", method = c("model", "local", "rerun"), cfg = SILUCARB,
                            flag = TRUE, experiment = NULL, compare = NULL) {
  method <- match.arg(method)
  m <- if (is.list(map)) map else tryCatch(climate_map(map, cfg = cfg), error = function(e) NULL)
  age <- if (is.list(map)) map$age else map
  if (!.is_reference(experiment, age) && method != "rerun") {
    message("An orbit / CO2 / Sun experiment needs SPEEDY runs: using method = \"rerun\".")
    method <- "rerun"
  }
  cl <- climate_fields(age, season, cfg)
  f_orig <- find_dem(cfg$dir_dem, age)
  cl$method <- "model"; cl$hadcm3l_co2 <- hadcm3l_co2(age)
  cl$temp_href <- orig_href(age, cfg)                # height the temperature refers to (for the lapse rate)
  if (is.null(m)) return(cl)                         # no corrected grid: nothing to correct
  g_o <- geo_1deg_cached(f_orig); g_c <- geo_1deg_cached(m$file)
  if (flag) cl$flag <- .fld(geography_flag(g_o, g_c))
  T1 <- .to_1deg(cl$temp); T1[!is.finite(T1)] <- 0
  if (method == "rerun") {
    d <- tryCatch(speedy_delta(m, season, cfg, experiment), error = function(e) { message("  ", conditionMessage(e)); NULL })
    if (is.null(d)) {
      message("The SPEEDY runs for ", m$name, if (!.is_reference(experiment, age)) paste0(" (", .resolve(experiment, age)$label, ")"),
              " are not finished yet (run_speedy() or the editor's Climate model tab); ",
              "using the local correction instead.")
      method <- "local"; experiment <- NULL
    } else {
      P1 <- pmax(.to_1deg(cl$precip), 0); U1 <- .to_1deg(cl$wind_u); V1 <- .to_1deg(cl$wind_v)
      U1[!is.finite(U1)] <- 0; V1[!is.finite(V1)] <- 0
      cl$precip <- .fld(pmax(P1 + .to_1deg(d$precip), 0))
      cl$wind_u <- .fld(U1 + .to_1deg(d$wind_u)); cl$wind_v <- .fld(V1 + .to_1deg(d$wind_v))
      cl$temp <- .fld(T1 + .to_1deg(d$t0))
      cl <- .water_balance(cl, .to_1deg(d$precip) - if (!is.null(d$evap)) .to_1deg(d$evap) else 0)
      ek <- tryCatch(ekman_change(m, experiment, season_months(season), cfg),
                     error = function(e) NULL)
      if (!is.null(ek)) { cl$w_ocean <- .fld(.to_1deg(cl$w_ocean) + .to_1deg(ek)); cl$ekman <- TRUE }   # + change of Ekman pumping
      ds <- tryCatch(.speedy_monthly_response(m, experiment, "sst", cfg), error = function(e) NULL)
      cl$sst_ocean <- .fld(.to_1deg(cl$sst_ocean) + if (!is.null(ds)) .to_1deg(.speedy_out_field(apply(ds[, , season_months(season), drop = FALSE], 1:2, mean))) else .to_1deg(d$t0))
      cl$method <- "rerun"; cl$speedy_years <- d$years
      cl$temp_href <- .speedy_in_field(speedy_geography(m$file)$orog)   # SPEEDY saw your map at ~3.75 deg
      if (!.is_reference(experiment, age)) {
        cl$experiment <- .resolve(experiment, age)
        cl$change <- speedy_experiment_change(m, season, experiment, cfg)   # experiment - reference
      }
      if (!is.null(compare)) {                       # experiment - other experiment
        cl$compare <- if (.is_reference(compare, age)) list(label = "reference (HadCM3L orbit, CO2 and Sun)") else .resolve(compare, age)
        cc <- if (.is_reference(compare, age)) NULL else speedy_experiment_change(m, season, compare, cfg)
        if (!.is_reference(compare, age) && is.null(cc)) {
          message("The runs of the comparison experiment are not finished: showing the change against the reference.")
          cl$compare <- NULL
        } else if (!is.null(cl$change) || !is.null(cc)) cl$change <- .add_diff(cl$change, cc, -1)
      }
      return(cl)
    }
  }
  if (method == "local") {
    lc <- local_correction(cl, g_o, g_c)
    P0 <- pmax(.to_1deg(cl$precip), 0)
    cl$precip <- lc$precip; cl$wind_u <- lc$wind_u; cl$wind_v <- lc$wind_v
    cl <- .water_balance(cl, lc$precip$z - P0)
    cl$temp <- .fld(T1 - 6.5 * (smooth_km(g_c$h, GRID1$lat, 100) - smooth_km(g_o$h, GRID1$lat, 100)) / 1000)   # lapse rate for changed relief
    cl$method <- "local"; cl$r2 <- lc$r2
    cl$temp_href <- .fld(smooth_km(g_c$h, GRID1$lat, 100))
  }
  cl
}

#' P - E and runoff on the 1-degree grid, changed by dPE (mm/day, change of P - E); runoff on
#' land the model did not have is taken as P - E (long-term water balance), never below 0
.water_balance <- function(cl, dPE) {
  pme <- .to_1deg(cl$pme); pme[!is.finite(pme)] <- 0
  ro <- .to_1deg(cl$runoff)
  pme <- pme + dPE
  ro[!is.finite(ro)] <- pmax(pme[!is.finite(ro)], 0)
  cl$pme <- .fld(pme); cl$runoff <- .fld(pmax(ro + dPE, 0))
  cl
}

#' the smooth terrain the HadCM3L temperatures refer to: the original PALEOMAP map at the
#' climate-model resolution (as SPEEDY's T30 grid, ~3.75 deg); cached
.HREF <- new.env()
orig_href <- function(age, cfg = SILUCARB) {
  f <- find_dem(cfg$dir_dem, age)
  if (is.null(.HREF[[f]])) .HREF[[f]] <- .speedy_in_field(speedy_geography(f)$orog)
  .HREF[[f]]
}
#' near-surface temperature brought to the map heights with a lapse rate (K/km): T - lapse *
#' (surface - reference height); surface = land height (0 at sea) [+ ice surface]
temp_downscale <- function(T, href, surf, lapse = 6.5) T - lapse / 1000 * (surf - href)
