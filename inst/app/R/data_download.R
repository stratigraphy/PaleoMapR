# =============================================================================
#  Download the open data sets (once). Files that already exist are kept, so
#  you can also copy them in by hand.
#   PaleoDEMs 0.2 deg (Scotese & Wright 2018)       EarthByte
#   PaleoCoastlines v7.2 (Kocsis & Scotese 2021)    Zenodo 7828387
#   PBDB Silurian collections                       paleobiodb.org data1.2
#   Macrostrat units 450-410 Ma                     macrostrat.org API v2
#   Sea-level curve (van der Meer et al. 2022)      Gondwana Research suppl.
#   Natural Earth borders / coastlines              naciscdn.org (on first map)
# =============================================================================

URL_PBDB  <- "https://paleobiodb.org/data1.2/colls/list.csv?interval=Silurian&show=loc,strat,lith,env,time"
URL_MACRO <- "https://macrostrat.org/api/v2/units?age_top=410&age_bottom=450&format=json&response=long"

# A downloaded file is only kept when it is complete: it is written to "<dest>.part", checked (zip, NetCDF,
# JSON or CSV header, by extension) and then renamed. An interrupted download or an error page from a proxy is
# therefore never taken for the real file.
.download_ok <- function(path) {
  ext <- tolower(tools::file_ext(sub("\\.part$", "", path)))
  sz <- file.info(path)$size
  if (!is.finite(sz) || sz == 0) {
    return(FALSE)
  }
  head <- tryCatch(readBin(path, "raw", 16), error = function(e) { raw(0) })
  txt <- tryCatch(rawToChar(head[head != as.raw(0)]), error = function(e) { "" })
  if (grepl("^\\s*<", txt)) {
    return(FALSE)
  }
  if (ext %in% c("zip", "xlsx")) {
    return(!is.null(tryCatch(utils::unzip(path, list = TRUE), error = function(e) { NULL }, warning = function(w) { NULL })))
  }
  if (ext == "nc") {
    nc <- tryCatch(suppressWarnings(ncdf4::nc_open(path)), error = function(e) { NULL })
    if (is.null(nc)) {
      return(FALSE)
    }
    ncdf4::nc_close(nc)
    return(TRUE)
  }
  if (ext == "json") {
    return(isTRUE(tryCatch(jsonlite::validate(paste(readLines(path, warn = FALSE), collapse = "\n")),
        error = function(e) { FALSE })))
  }
  if (ext == "csv") {
    first <- tryCatch(readLines(path, n = 12, warn = FALSE), error = function(e) { character(0) })
    return(length(first) > 0 && any(grepl(",", first)))
  }
  TRUE
}

# URLs that failed in this R session, with the time: they are not tried again for 10 minutes, so an offline
# computer does not wait for every download again at each rebuild
.FETCH_FAILED <- new.env(parent = emptyenv())

#' download one file unless it exists (creates the folder; keeps only complete files)
fetch <- function(url, dest, timeout = 600, tries = 2) {
  if (file.exists(dest) && isTRUE(file.info(dest)$size > 0)) {
    return(invisible(TRUE))
  }
  t_fail <- .FETCH_FAILED[[url]]
  if (!is.null(t_fail) && difftime(Sys.time(), t_fail, units = "mins") < 10) {
    message("skipping ", basename(dest), ": its download failed a few minutes ago")
    return(invisible(FALSE))
  }
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  op <- options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(op), add = TRUE)
  part <- paste0(dest, ".part")
  message("downloading ", basename(dest), " ...")
  ok <- FALSE
  for (k in seq_len(tries)) {
    unlink(part)
    got <- tryCatch({
      utils::download.file(url, part, mode = "wb", quiet = TRUE,
          method = if (capabilities("libcurl")) { "libcurl" } else { "auto" })
      TRUE
    }, error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      FALSE
    }, warning = function(w) {
      message("  FAILED: ", conditionMessage(w))
      FALSE
    })
    if (got && .download_ok(part)) {
      ok <- file.rename(part, dest)
      if (!ok) {
        ok <- file.copy(part, dest, overwrite = TRUE)
      }
      break
    }
    if (got) {
      message("  the downloaded file is incomplete or not the expected data")
    }
  }
  unlink(part)
  if (!ok) {
    assign(url, Sys.time(), envir = .FETCH_FAILED)
    message("  download it by hand from\n  ", url, "\n  and save it as\n  ", dest)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# ---- PBDB / Macrostrat records in 50-Myr windows (any part of the Phanerozoic) ----
RECORD_WINDOW_MYR <- 50
.win_files <- function(lo, hi, cfg) list(
  pbdb = file.path(cfg$dir_data, "pbdb", sprintf("pbdb_%03d_%03d.csv", lo, hi)),
  macro = file.path(dirname(cfg$file_macro), sprintf("units_%03d_%03d.json", lo, hi)))

#' all record files present (legacy Silurian files + the window files)
record_source_files <- function(cfg = SILUCARB) {
  list(pbdb = c(if (file.exists(cfg$file_pbdb)) cfg$file_pbdb,
                list.files(file.path(cfg$dir_data, "pbdb"), "^pbdb_.*\\.csv$", full.names = TRUE)),
       macro = c(if (file.exists(cfg$file_macro)) cfg$file_macro,
                 list.files(dirname(cfg$file_macro), "^units_[0-9]{3}_[0-9]{3}\\.json$", full.names = TRUE)))
}

#' make sure PBDB and Macrostrat records cover the map ages (downloads missing 50-Myr
#' windows once, e.g. 400-450 Ma). Returns TRUE if something new was downloaded.
ensure_records_for_ages <- function(ages, cfg = SILUCARB) {
  lo <- max(0, floor((min(ages) - 15) / RECORD_WINDOW_MYR) * RECORD_WINDOW_MYR)
  hi <- min(550, ceiling((max(ages) + 15) / RECORD_WINDOW_MYR) * RECORD_WINDOW_MYR)
  new <- FALSE
  for (w in seq(lo, hi - RECORD_WINDOW_MYR, by = RECORD_WINDOW_MYR)) {
    f <- .win_files(w, w + RECORD_WINDOW_MYR, cfg)
    # the Silurian start-up files do not cover the whole 400-450 Ma window (PBDB: Silurian only), so the
    # window files are always downloaded; duplicates are removed when the record table is built
    legacy_ok <- FALSE
    if (!file.exists(f$pbdb) && !legacy_ok) {
      ok <- fetch(sprintf("https://paleobiodb.org/data1.2/colls/list.csv?max_ma=%d&min_ma=%d&timerule=overlap&show=loc,strat,lith,env,time",
                          w + RECORD_WINDOW_MYR, w), f$pbdb, timeout = 300)
      new <- new || ok
    }
    if (!file.exists(f$macro) && !legacy_ok) {
      ok <- fetch(sprintf("https://macrostrat.org/api/v2/units?age_top=%d&age_bottom=%d&format=json&response=long",
                          w, w + RECORD_WINDOW_MYR), f$macro, timeout = 300)
      new <- new || ok
    }
  }
  invisible(new)
}

#' the record files for the configured ages; called automatically when needed
download_records_data <- function(cfg = SILUCARB, ages = cfg$ages) {
  ensure_records_for_ages(ages, cfg)
  src <- record_source_files(cfg)
  if (!length(src$pbdb) || !length(src$macro)) stop("PBDB / Macrostrat data missing (see the messages above). ",
                     "Check the internet connection, or download the files by hand.", call. = FALSE)
  invisible(TRUE)
}

#' PaleoDEMs 0.2 deg (Scotese & Wright 2018): downloads the ~280 MB zip once and
#' unpacks all grids (0-540 Ma). Called automatically when grids are missing.
download_dems <- function(cfg = SILUCARB, timeout = 3600) {
  if (length(list.files(cfg$dir_dem, "Ma\\.nc$")) >= 80) return(invisible(TRUE))
  z <- file.path(cfg$dir_data, "paleoDEM_0.2Deg_grids.zip")
  message("The PaleoDEMs (0.2 deg) are not there yet: downloading them once (~280 MB, a few minutes) ...")
  if (!fetch("https://www.earthbyte.org/webdav/ftp/Data_Collections/Scotese_Wright_2018_PaleoDEM/paleoDEM_0.2Deg_grids.zip",
             z, timeout)) return(invisible(FALSE))
  lst <- tryCatch(unzip(z, list = TRUE)$Name, error = function(e) NULL)
  if (is.null(lst)) {                                   # broken / partial zip: remove, so a rerun starts fresh
    unlink(z); message("  the downloaded zip is damaged; run download_dems() again"); return(invisible(FALSE))
  }
  dir.create(cfg$dir_dem, showWarnings = FALSE, recursive = TRUE)
  unzip(z, files = lst[grepl("Ma\\.nc$", lst)], exdir = cfg$dir_dem, junkpaths = TRUE)   # all ages 0-540 Ma
  message("  PaleoDEMs unpacked to ", cfg$dir_dem)
  invisible(TRUE)
}

download_data <- function(cfg = SILUCARB, timeout = 3600) {
  # small files first, so they are there even if a large download is interrupted
  fetch(URL_PBDB, cfg$file_pbdb)                     # PBDB collections
  fetch(URL_MACRO, cfg$file_macro)                   # Macrostrat units 450-410 Ma
  fetch("https://ars.els-cdn.com/content/image/1-s2.0-S1342937X22002192-mmc1.xlsx",
        cfg$file_sealevel)                           # sea-level curve (steps 1 and 4)
  # PaleoCoastlines (shapefiles: CS = coastlines, CM = margins)
  if (!dir.exists(file.path(cfg$dir_coast, "CS"))) {
    z <- file.path(cfg$dir_coast, "shp.zip")
    if (fetch("https://zenodo.org/records/7828387/files/paleocoastlines_v7_shapefiles.zip?download=1", z, timeout)) {
      unzip(z, exdir = cfg$dir_coast)
      for (inner in list.files(cfg$dir_coast, "\\.zip$", full.names = TRUE, recursive = TRUE))
        if (normalizePath(inner) != normalizePath(z)) unzip(inner, exdir = dirname(inner))
      if (!dir.exists(file.path(cfg$dir_coast, "CS"))) {          # move CS/ and CM/ up if nested
        for (d in c("CS", "CM")) {
          src <- list.dirs(cfg$dir_coast)[basename(list.dirs(cfg$dir_coast)) == d][1]
          if (!is.na(src)) file.copy(src, cfg$dir_coast, recursive = TRUE)
        }
      }
      if (!dir.exists(file.path(cfg$dir_coast, "CS"))) {
        # a damaged zip would otherwise be reused by every later run
        unlink(z)
        message("  the PaleoCoastlines zip could not be unpacked; run the setup again")
      }
    }
  }
  download_dems(cfg, timeout)                        # PaleoDEMs 0.2 deg (~280 MB zip)
  # PaleoDEMs 1 deg (used by steps/silurian_ocean_area.R for the uncertainty range)
  d1 <- file.path(cfg$dir_data, "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2")
  if (!length(list.files(d1, "\\.nc$"))) {
    z <- file.path(cfg$dir_data, "PaleoDEMs_1deg_nc.zip")
    if (fetch("https://www.earthbyte.org/webdav/ftp/Data_Collections/Scotese_Wright_2018_PaleoDEM/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc.zip", z, timeout)) {
      dir.create(d1, showWarnings = FALSE, recursive = TRUE)
      lst <- tryCatch(unzip(z, list = TRUE)$Name, error = function(e) { NULL })
      if (is.null(lst)) {
        unlink(z)
      } else {
        unzip(z, files = lst[grepl("Ma\\.nc$", lst)], exdir = d1, junkpaths = TRUE)
      }
    }
  }
  ok <- c(PBDB = file.exists(cfg$file_pbdb), Macrostrat = file.exists(cfg$file_macro),
          PaleoCoastlines = dir.exists(file.path(cfg$dir_coast, "CS")),
          PaleoDEMs = length(list.files(cfg$dir_dem, "\\.nc$")) > 0)
  print(ok)
  if (!all(ok)) message("Missing: ", paste(names(ok)[!ok], collapse = ", "), " - see the messages above.")
  invisible(all(ok))
}
