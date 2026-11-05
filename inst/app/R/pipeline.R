# =============================================================================
#  Pipeline: records -> corrected PaleoDEMs (parallel) -> statistics
# =============================================================================

#' records used for the DEM correction at each age (needs pyGPlates; main process)
records_by_age <- function(ages = SILUCARB$ages, cfg = SILUCARB) {
  ensure_records_for_ages(ages, cfg)          # PBDB / Macrostrat for these ages (downloaded once)
  recs <- all_records(cfg)
  setNames(lapply(ages, function(a) records_for_age(a, purpose = "dem", recs = recs, cfg = cfg)), ages)
}

#' build the fossil-corrected PaleoDEMs for `ages` (one age per worker)
run_dem_correction <- function(ages = SILUCARB$ages, cfg = SILUCARB, workers = cfg$workers) {
  dir.create(cfg$dir_results, showWarnings = FALSE, recursive = TRUE)
  message("Reconstructing records (pyGPlates) ...")
  rb <- records_by_age(ages, cfg)
  one <- function(a, rb, cfg) {
    r <- correct_dem(a, rb[[as.character(a)]], cfg)
    list(stats = r$stats, seaways = r$seaways)
  }
  workers <- max(1, min(workers, length(ages)))
  if (workers > 1) {
    message("Correcting ", length(ages), " ages on ", workers, " workers ...")
    cl <- parallel::makeCluster(workers); on.exit(parallel::stopCluster(cl))
    parallel::clusterCall(cl, function(dir_R) {
      suppressMessages({ source(file.path(dir_R, "grid_tools.R")); source(file.path(dir_R, "dem_correction.R")) }); NULL },
      cfg$dir_R)
    res <- parallel::parLapply(cl, ages, one, rb = rb, cfg = cfg)
  } else {
    res <- lapply(ages, one, rb = rb, cfg = cfg)
  }
  stats <- do.call(rbind, lapply(res, `[[`, "stats"))
  sw <- do.call(rbind, lapply(res, `[[`, "seaways"))
  # merge with earlier runs of other ages
  f1 <- file.path(cfg$dir_results, "step6_corrected_dem_stats.csv")
  f2 <- file.path(cfg$dir_results, "step6_seaways.csv")
  merge_rows <- function(old, new) {             # keep other ages from earlier runs; new columns -> NA there
    old <- old[!old$map_age %in% ages, , drop = FALSE]
    if (!nrow(old)) return(new)
    for (k in setdiff(names(new), names(old))) old[[k]] <- NA
    for (k in setdiff(names(old), names(new))) new[[k]] <- NA
    rbind(old[, names(new)], new)
  }
  if (file.exists(f1)) stats <- merge_rows(read.csv(f1), stats)
  if (file.exists(f2) && !is.null(sw)) sw <- merge_rows(read.csv(f2), sw)
  stats <- stats[order(-stats$map_age), ]
  write.csv(stats, f1, row.names = FALSE)
  if (!is.null(sw)) write.csv(sw, f2, row.names = FALSE)
  used <- do.call(rbind, lapply(names(rb), function(a) if (nrow(rb[[a]])) cbind(map_age = as.numeric(a), rb[[a]][, c(
    "record_id", "source", "name", "setting", "indicator", "lat", "lon", "plat", "plon", "edited")])))
  write.csv(used, file.path(cfg$dir_results, "step6_records_used.csv"), row.names = FALSE)
  print(stats[, c("map_age", "n_marine", "n_evaporite", "n_terrestrial", "n_user_points",
                  "marine_on_land_after_pct", "seaways_added", "marine_water_bodies_still_isolated",
                  "seaway_max_narrowing_km", "slope_max_changed_after", "ocean_DEM_Mkm2", "ocean_corrected_Mkm2")])
  invisible(stats)
}

#' check that every marine / evaporite record is in the sea of the corrected grid
verify_corrected <- function(ages = SILUCARB$ages, cfg = SILUCARB) {
  rb <- records_by_age(ages, cfg)
  do.call(rbind, lapply(ages, function(a) {
    d <- read_dem_nc(find_dem(cfg$dir_corrected, a)); r <- rb[[as.character(a)]]
    zz <- sample_grid(d$z, d$lat, d$lon, r$plat, r$plon)
    data.frame(map_age = a, marine_on_land = sum(r$setting %in% c("marine", "marine_deep") & zz >= 0),
               evaporite_above_sl = sum(r$setting == "evaporite" & zz > 0),
               terrestrial_in_sea = sum(r$setting == "terrestrial" & zz < 0), n = nrow(r))
  }))
}

#' ocean area (Mkm2) of raw and corrected grids per map age
ocean_areas <- function(ages = SILUCARB$ages, cfg = SILUCARB) {
  do.call(rbind, lapply(ages, function(a) {
    raw <- read_dem_nc(find_dem(cfg$dir_dem, a)); W <- area_matrix(raw$lat, ncol(raw$z))
    f <- list.files(cfg$dir_corrected, paste0("_", a, "Ma\\.nc$"), full.names = TRUE)
    cor <- if (length(f)) read_dem_nc(f[1])$z else NA
    data.frame(map_age = a, ocean_raw_Mkm2 = round(sum(W[raw$z < 0]), 2),
               ocean_corrected_Mkm2 = if (length(f)) round(sum(W[cor < 0]), 2) else NA,
               shelf_corrected_Mkm2 = if (length(f)) round(sum(W[cor < 0 & cor >= -1400]), 2) else NA)
  }))
}
