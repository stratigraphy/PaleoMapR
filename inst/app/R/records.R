# =============================================================================
#  Records (data points) used to correct and annotate the maps
#
#  Three tables in points/:
#   records_base.csv   built from PBDB + Macrostrat by build_base_records()
#                      (regenerate, do not edit by hand)
#   user_points.csv    YOUR points: add rows (or click them in the app)
#   record_edits.csv   YOUR changes to database records: exclude / include /
#                      set_setting / set_indicator, by record_id
#
#  Columns of every record
#   record_id   PBDB:<collection_no>, MS:<unit_id> or U:<point_id>
#   setting     marine | marine_deep | evaporite | terrestrial | "" (not used for the DEM)
#               marine      -> must plot in the sea, shallow shelf (flooded / seaway if needed)
#               marine_deep -> sea, deeper part of the shelf / outer shelf (MARINE_DEEP_DEPTH m):
#                              steers the deeper-water topography; never deep ocean
#               evaporite   -> at or just below sea level (-5 m)
#               terrestrial -> kept on land where no marine record is close
#   indicator   reef | evaporite | glacial | "" (climate symbol on the maps)
#   age_young, age_old   Ma; used for map age A when
#               age_young - age_window <= A <= age_old + age_window
# =============================================================================
# jsonlite is used as jsonlite::fromJSON (not attached: it would mask shiny::validate)

USER_POINT_COLS <- c("point_id", "setting", "indicator", "lat", "lon", "coords", "model",
                     "clicked_age", "age_young", "age_old", "active", "note")
EDIT_COLS <- c("record_id", "action", "value", "note")

points_file <- function(name, cfg = SILUCARB) file.path(cfg$dir_points, name)

#' build points/records_base.csv from the PBDB and Macrostrat downloads
build_base_records <- function(cfg = SILUCARB) {
  dir.create(cfg$dir_points, showWarnings = FALSE, recursive = TRUE)
  src <- record_source_files(cfg)
  if (!length(src$pbdb) || !length(src$macro)) { download_records_data(cfg); src <- record_source_files(cfg) }
  # ---- PBDB collections (data1.2 colls/list.csv; a file may start with 'Warning:' /
  #      'Records:' lines before the header). Several age windows are merged.
  read_pbdb <- function(f) {
    first <- readLines(f, n = 10); skip <- grep('^"collection_no"', first)[1] - 1
    if (is.na(skip)) return(NULL)
    read.csv(f, skip = skip, stringsAsFactors = FALSE, colClasses = "character")
  }
  cl <- Filter(Negate(is.null), lapply(src$pbdb, read_pbdb))
  cols <- Reduce(intersect, lapply(cl, names))
  c <- do.call(rbind, lapply(cl, function(x) x[, cols]))
  c <- c[!duplicated(c$collection_no), ]
  for (k in c("lat", "lng", "max_ma", "min_ma")) c[[k]] <- as.numeric(c[[k]])
  for (k in c("environment", "lithology1", "lithology2", "lithdescript", "formation"))
    c[[k]][!is.na(c[[k]]) & c[[k]] == ""] <- NA
  c <- c[(c$max_ma - c$min_ma) <= cfg$pbdb_max_age_range, ]
  mid <- (c$max_ma + c$min_ma) / 2
  env <- tolower(ifelse(is.na(c$environment), "", c$environment))
  lith <- tolower(paste(c$lithology1, c$lithology2, c$lithdescript, sep = "|"))
  lith <- gsub("\\bna\\b", "", lith)
  setting <- ifelse(grepl(cfg$terrestrial_env, env), "terrestrial",
                    ifelse(grepl(cfg$evaporite_lith, lith), "evaporite", "marine"))
  indicator <- ifelse(grepl(cfg$reef_env, env), "reef",
                      ifelse(grepl(cfg$evaporite_lith, lith), "evaporite", ""))
  fm <- ifelse(is.na(c$formation), "", c$formation)
  glac <- grepl("Cancaniri", fm) | grepl("tillite", lith) | grepl("glacial", env)
  indicator[glac] <- "glacial"
  # glacial deposits are not marine by themselves (tillites are mostly terrestrial,
  # diamictites can be glaciomarine): an environment that is only 'glacial' gives no
  # land/sea constraint; an explicit marine/coastal environment (e.g. foreshore) still counts
  setting[glac & (env == "" | grepl("^glacial", env))] <- ""
  pb <- data.frame(record_id = paste0("PBDB:", c$collection_no), source = "PBDB",
                   name = ifelse(fm == "", c$collection_name, fm),
                   lat = c$lat, lon = c$lng, age_young = mid, age_old = mid,
                   max_ma = c$max_ma, min_ma = c$min_ma,
                   setting = setting, indicator = indicator,
                   environment = env, lithology = substr(lith, 1, 120), stringsAsFactors = FALSE)
  # ---- Macrostrat units (all projects), all downloaded age windows
  u <- do.call(base::c, lapply(src$macro, function(f) jsonlite::fromJSON(f, simplifyVector = FALSE)$success$data))
  u <- u[!duplicated(vapply(u, function(x) as.character(x$unit_id), ""))]
  ms <- do.call(rbind, lapply(u, function(x) {
    t <- x$t_age; b <- x$b_age
    if (is.null(t) || is.null(b) || b - t > cfg$macro_max_age_range) return(NULL)
    ln <- paste(vapply(x$lith, function(l) l$name, ""), collapse = " ")
    en <- paste(vapply(x$environ, function(e) e$name, ""), collapse = " ")
    cl <- vapply(x$environ, function(e) e$class, "")
    # glacial environments (tillite, diamictite, 'glacial indet.') can be terrestrial OR
    # glaciomarine: they say nothing about land or sea by themselves. Only the unit's
    # other environments (e.g. submarine fan, offshore -> marine; fluvial -> terrestrial) count.
    envn <- vapply(x$environ, function(e) e$name, "")
    cl <- cl[!grepl("glac", envn, ignore.case = TRUE)]
    evap <- grepl("evaporite|gypsum|anhydrite|halite", ln)
    setting <- if (evap) "evaporite" else if ("marine" %in% cl) "marine" else
      if ("non-marine" %in% cl) "terrestrial" else ""
    ind <- if (evap) "evaporite" else
      if (grepl("glacial", en) && grepl("diamictite|tillite", ln)) "glacial" else ""
    if (setting == "" && ind == "") return(NULL)
    nm <- if (!is.null(x$Fm) && nzchar(x$Fm)) x$Fm else x$unit_name
    data.frame(record_id = paste0("MS:", x$unit_id), source = "Macrostrat", name = nm,
               lat = x$clat, lon = x$clng, age_young = t, age_old = b, max_ma = b, min_ma = t,
               setting = setting, indicator = ind, environment = en, lithology = substr(ln, 1, 120),
               stringsAsFactors = FALSE)
  }))
  base <- rbind(pb, ms)
  base <- base[is.finite(base$lat) & is.finite(base$lon), ]
  write.csv(base, points_file("records_base.csv", cfg), row.names = FALSE)
  records_rules_stamp(cfg)
  message(nrow(pb), " PBDB collections and ", nrow(ms), " Macrostrat units -> ",
          points_file("records_base.csv", cfg))
  invisible(base)
}

#' create empty user tables (with one commented example row each) if missing
init_user_tables <- function(cfg = SILUCARB) {
  dir.create(cfg$dir_points, showWarnings = FALSE, recursive = TRUE)
  f <- points_file("user_points.csv", cfg)
  if (!file.exists(f)) {
    ex <- data.frame(point_id = "example1", setting = "marine", indicator = "",
                     lat = 50.52, lon = 5.24, coords = "present", model = "scotese_and_wright2018",
                     clicked_age = NA, age_young = 430, age_old = 440, active = FALSE,
                     note = "example row (inactive): Huy, present-day coordinates")
    write.csv(ex, f, row.names = FALSE)
  }
  f <- points_file("record_edits.csv", cfg)
  if (!file.exists(f)) {
    write.csv(data.frame(record_id = character(0), action = character(0), value = character(0),
                         note = character(0)), f, row.names = FALSE)
  }
  invisible(TRUE)
}

# Reads a points table as written by the editor, by hand or by Excel: UTF-8 with or without a byte-order mark,
# comma- or semicolon-separated (Excel in many European languages, with decimal commas). All columns are read as
# text. Stops (and so never lets the table be written back) when a required column is missing.
.read_points_csv <- function(f, required) {
  first <- readLines(f, n = 1, warn = FALSE, encoding = "UTF-8")
  if (!length(first) || !nzchar(trimws(sub("^\ufeff", "", first)))) {
    out <- as.data.frame(matrix(character(0), 0, length(required), dimnames = list(NULL, required)), stringsAsFactors = FALSE)
    return(out)
  }
  bare <- gsub('"[^"]*"', "", first)
  semi <- lengths(regmatches(bare, gregexpr(";", bare))) > lengths(regmatches(bare, gregexpr(",", bare)))
  d <- if (semi) {
    utils::read.csv2(f, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA", ""))
  } else {
    utils::read.csv(f, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA", ""))
  }
  names(d) <- trimws(sub("^\ufeff|^X\\.U\\.FEFF\\.", "", names(d)))
  miss <- setdiff(required, names(d))
  if (length(miss)) {
    stop(basename(f), ": column(s) ", paste(miss, collapse = ", "), " not found. Check the file (it must have a header ",
        "row; separators , or ;). Nothing was changed.", call. = FALSE)
  }
  attr(d, "decimal_comma") <- semi
  d
}
# text -> number, accepting a decimal comma
.as_num <- function(x) {
  suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(x)), fixed = TRUE)))
}
# text -> TRUE / FALSE in several languages (Excel writes WAAR, VRAI, WAHR ...); NA stays NA
.as_flag <- function(x) {
  v <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(v))
  out[v %in% c("TRUE", "T", "1", "YES", "Y", "WAAR", "VRAI", "WAHR", "VERO", "VERDADERO", "SANT")] <- TRUE
  out[v %in% c("FALSE", "F", "0", "NO", "N", "ONWAAR", "FAUX", "FALSCH", "FALSO", "FALSK")] <- FALSE
  out
}

read_user_points <- function(cfg = SILUCARB) {
  f <- points_file("user_points.csv", cfg)
  if (!file.exists(f)) init_user_tables(cfg)
  u <- .read_points_csv(f, c("lat", "lon"))
  for (k in setdiff(USER_POINT_COLS, names(u))) u[[k]] <- rep(NA, nrow(u))
  for (k in c("lat", "lon", "clicked_age", "age_young", "age_old")) u[[k]] <- .as_num(u[[k]])
  for (k in c("setting", "indicator", "coords", "model", "note", "point_id")) {
    u[[k]] <- trimws(ifelse(is.na(u[[k]]), "", as.character(u[[k]])))
  }
  u$setting <- tolower(u$setting)
  u$indicator <- tolower(u$indicator)
  u$coords <- tolower(u$coords)
  u$coords[u$coords %in% c("paleo", "palaeo")] <- "palaeo"
  u$coords[u$coords != "palaeo"] <- "present"
  u$model[u$model == ""] <- "scotese_and_wright2018"
  u$active <- .as_flag(u$active)
  u$active[is.na(u$active)] <- TRUE
  u$age_old[is.na(u$age_old)] <- u$age_young[is.na(u$age_old)]
  u$age_young[is.na(u$age_young)] <- u$age_old[is.na(u$age_young)]
  sw <- is.finite(u$age_young) & is.finite(u$age_old) & u$age_young > u$age_old
  if (any(sw)) {
    tmp <- u$age_young[sw]
    u$age_young[sw] <- u$age_old[sw]
    u$age_old[sw] <- tmp
  }
  bad <- !(u$setting %in% c("marine", "marine_deep", "evaporite", "terrestrial", ""))
  if (any(bad)) {
    warning("user_points.csv: setting must be marine, marine_deep, evaporite, terrestrial or empty; rows ",
        paste(which(bad) + 1, collapse = ", "), " are not used for the correction", call. = FALSE)
    u$setting[bad] <- ""
  }
  # every point needs its own ID (empty or repeated IDs would make edits hit the wrong point)
  no_id <- u$point_id == ""
  u$point_id[no_id] <- sprintf("row%d", which(no_id) + 1)
  dup <- duplicated(u$point_id)
  if (any(dup)) {
    warning("user_points.csv: repeated point_id renamed: ", paste(unique(u$point_id[dup]), collapse = ", "), call. = FALSE)
    u$point_id <- make.unique(u$point_id, sep = "_")
  }
  u[, USER_POINT_COLS]
}

read_record_edits <- function(cfg = SILUCARB) {
  f <- points_file("record_edits.csv", cfg)
  if (!file.exists(f)) init_user_tables(cfg)
  e <- .read_points_csv(f, c("record_id", "action"))
  for (k in setdiff(EDIT_COLS, names(e))) e[[k]] <- rep(NA, nrow(e))
  for (k in EDIT_COLS) {
    e[[k]] <- trimws(ifelse(is.na(e[[k]]), "", as.character(e[[k]])))
  }
  e$action <- tolower(e$action)
  bad <- e$record_id == "" | e$action == ""
  if (any(bad)) {
    warning("record_edits.csv: rows without record_id or action are ignored (rows ", paste(which(bad) + 1, collapse = ", "), ")",
        call. = FALSE)
  }
  e[!bad, EDIT_COLS]
}

# Writes a points table; a file that is open in Excel (locked on Windows) gives a clear message instead of a
# crash, and nothing is changed
.write_points_csv <- function(d, f) {
  ok <- tryCatch(suppressWarnings({
    tmp <- paste0(f, ".tmp")
    utils::write.csv(d, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (!file.rename(tmp, f)) {
      if (!file.copy(tmp, f, overwrite = TRUE)) {
        stop("locked")
      }
      unlink(tmp)
    }
    TRUE
  }), error = function(e) { FALSE })
  if (!ok) {
    unlink(paste0(f, ".tmp"))
    stop(sprintf("Could not save %s: the file is open in another program (e.g. Excel) or read-only. Close it and try again; nothing was changed.",
        basename(f)), call. = FALSE)
  }
  invisible(TRUE)
}

#' convert palaeo-coordinate user points to present day (needed to use them at
#' other ages); results are cached in the returned table (lat_present, lon_present)
user_points_present <- function(u) {
  u$lat_present <- ifelse(u$coords == "present", u$lat, NA)
  u$lon_present <- ifelse(u$coords == "present", u$lon, NA)
  pal <- which(u$coords == "palaeo" & is.finite(u$lat) & is.finite(u$clicked_age))
  lost <- u$coords == "palaeo" & !is.finite(u$clicked_age)
  if (any(lost)) {
    warning("user points in palaeo-coordinates need clicked_age (the map age they belong to); not used: ",
        paste(u$point_id[lost], collapse = ", "), call. = FALSE)
  }
  if (length(pal)) {
    grp <- split(pal, paste(u$model[pal], u$clicked_age[pal]))
    for (g in grp) {                       # one pyGPlates call per (model, age)
      p <- gp_to_present(u$lat[g], u$lon[g], u$clicked_age[g[1]], u$model[g[1]])
      u$lat_present[g] <- p$lat; u$lon_present[g] <- p$lon
    }
  }
  u
}

# The classification rules are in records.R. A checksum of that file is stored next to records_base.csv, so the
# table is rebuilt when the rules change, but not when an update of the package only gives the file a new date.
# The settings in config.R that the classification uses are part of the checksum too.
RECORD_RULE_KEYS <- c("terrestrial_env", "evaporite_lith", "reef_env", "pbdb_max_age_range", "macro_max_age_range")
.records_rules_hash <- function(src, cfg) {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(c(unname(tools::md5sum(src)), deparse(cfg[RECORD_RULE_KEYS])), tf)
  unname(tools::md5sum(tf))
}
records_rules_stamp <- function(cfg = SILUCARB) {
  src <- file.path(cfg$dir_R, "records.R")
  if (file.exists(src)) {
    writeLines(.records_rules_hash(src, cfg), paste0(points_file("records_base.csv", cfg), ".rules"))
  }
  invisible(NULL)
}
records_rules_changed <- function(f, src, cfg = SILUCARB) {
  if (!file.exists(src)) {
    return(FALSE)
  }
  stamp <- paste0(f, ".rules")
  if (!file.exists(stamp)) {
    # tables from before the checksum: compare dates as before, then keep the checksum
    changed <- file.mtime(src) > file.mtime(f)
    if (!changed) {
      writeLines(.records_rules_hash(src, cfg), stamp)
    }
    return(changed)
  }
  !identical(readLines(stamp, warn = FALSE)[1], .records_rules_hash(src, cfg))
}

#' all records after applying your edits and adding your points
#' (one row per record; coordinates are present-day)
all_records <- function(cfg = SILUCARB, include_inactive = FALSE) {
  f <- points_file("records_base.csv", cfg)
  src <- file.path(cfg$dir_R, "records.R")              # classification rules changed -> rebuild
  newest <- suppressWarnings(max(file.mtime(unlist(record_source_files(cfg)))))
  if (!file.exists(f) || records_rules_changed(f, src, cfg) ||
      (is.finite(newest) && newest > file.mtime(f))) build_base_records(cfg)
  base <- read.csv(f, stringsAsFactors = FALSE, na.strings = "NA")
  base$setting[is.na(base$setting)] <- ""; base$indicator[is.na(base$indicator)] <- ""
  base$active <- TRUE; base$edited <- ""
  base$coords <- "present"; base$model <- "scotese_and_wright2018"; base$clicked_age <- NA
  base$plat_clicked <- NA; base$plon_clicked <- NA
  # ---- edits (applied in file order: a later row wins)
  e <- read_record_edits(cfg)
  for (k in seq_len(nrow(e))) {
    i <- match(e$record_id[k], base$record_id)
    if (is.na(i)) { warning("record_edits.csv: unknown record_id ", e$record_id[k]); next }
    a <- tolower(e$action[k])
    if (a == "exclude") base$active[i] <- FALSE
    else if (a == "include") base$active[i] <- TRUE
    else if (a == "set_setting") base$setting[i] <- e$value[k]
    else if (a == "set_indicator") base$indicator[i] <- e$value[k]
    else warning("record_edits.csv: unknown action ", e$action[k])
    base$edited[i] <- trimws(paste(base$edited[i], a))
  }
  # ---- user points
  u <- read_user_points(cfg)
  if (nrow(u)) {
    u <- user_points_present(u)
    up <- data.frame(record_id = paste0("U:", u$point_id), source = "user", name = u$note,
                     lat = u$lat_present, lon = u$lon_present, age_young = u$age_young,
                     age_old = u$age_old, max_ma = u$age_old, min_ma = u$age_young,
                     setting = u$setting, indicator = u$indicator, environment = "", lithology = "",
                     active = u$active, edited = "user point", coords = u$coords, model = u$model,
                     clicked_age = u$clicked_age,
                     plat_clicked = ifelse(u$coords == "palaeo", u$lat, NA),
                     plon_clicked = ifelse(u$coords == "palaeo", u$lon, NA), stringsAsFactors = FALSE)
    base <- rbind(base, up[, names(base)])
  }
  if (!include_inactive) base <- base[base$active, ]
  base
}

.NOT_YET_SAID <- new.env(parent = emptyenv())   # said once per age, model and purpose

#' records used at map `age`, with palaeo-positions in `model`
#' purpose = "dem": records with a setting, chosen by age range (DEM correction)
#' purpose = "indicator": records with an indicator, chosen by mid-age within the window
records_for_age <- function(age, model = "scotese_and_wright2018", purpose = c("dem", "indicator"),
                            recs = NULL, cfg = SILUCARB) {
  purpose <- match.arg(purpose)
  if (is.null(recs)) recs <- all_records(cfg)
  w <- cfg$age_window
  if (purpose == "dem") {
    r <- recs[recs$setting != "" & recs$age_young - w <= age & recs$age_old + w >= age, ]
  } else {
    mid <- (recs$age_young + recs$age_old) / 2
    r <- recs[recs$indicator != "" & abs(mid - age) <= w, ]
  }
  if (!nrow(r)) return(cbind(r, plat = numeric(0), plon = numeric(0)))
  r$plat <- NA_real_; r$plon <- NA_real_
  # rows with impossible coordinates are left out (pyGPlates would stop on them)
  ok <- is.finite(r$lat) & is.finite(r$lon) & abs(r$lat) <= 90
  # (a point clicked on the palaeo-map of this age keeps its clicked position, also off any plate)
  here <- r$coords == "palaeo" & !is.na(r$clicked_age) & r$clicked_age == age & r$model == model
  bad_user <- !ok & !here & r$source == "user" & r$coords != "palaeo"
  if (any(bad_user)) {
    warning("user points with missing or impossible coordinates are not used: ",
        paste(sub("^U:", "", r$record_id[bad_user]), collapse = ", "), call. = FALSE)
  }
  r$lon <- ((r$lon + 180) %% 360) - 180
  if (any(ok)) {
    p <- gp_reconstruct(r$lat[ok], r$lon[ok], rep(age, sum(ok)), model)
    r$plat[ok] <- p$plat; r$plon[ok] <- p$plon
    if (any(p$not_yet)) {
      key <- paste(age, model, purpose, sum(p$not_yet))
      if (is.null(.NOT_YET_SAID[[key]])) {
        assign(key, TRUE, envir = .NOT_YET_SAID)
        message(sum(p$not_yet), " ", purpose, " record(s) at ", age, " Ma (", model, ") lie on crust that the plate model ",
            "adds only later; they are not used")
      }
    }
  }
  # a point clicked on a palaeo map at this age in this model keeps its exact position
  same <- r$coords == "palaeo" & !is.na(r$clicked_age) & r$clicked_age == age & r$model == model
  r$plat[same] <- r$plat_clicked[same]; r$plon[same] <- r$plon_clicked[same]
  r <- r[is.finite(r$plat), ]
  if (purpose == "dem") r <- r[!duplicated(paste(round(r$plat, 5), round(r$plon, 5), r$setting)), ]
  r
}

#' permanently delete your own point(s) from user_points.csv (point_id, or "U:<id>")
delete_user_point <- function(point_id, cfg = SILUCARB) {
  u <- read_user_points(cfg); id <- sub("^U:", "", point_id)
  .write_points_csv(u[!u$point_id %in% id, ], points_file("user_points.csv", cfg))
  invisible(sum(u$point_id %in% id))
}

#' add a point to user_points.csv (used by the app; also handy from the console)
add_user_point <- function(lat, lon, setting = "marine", indicator = "", coords = "present",
                           model = "scotese_and_wright2018", clicked_age = NA,
                           age_young, age_old = age_young, note = "", cfg = SILUCARB) {
  u <- read_user_points(cfg)
  id <- sprintf("p%04d", max(c(0, suppressWarnings(as.integer(sub("^p", "", u$point_id)))), na.rm = TRUE) + 1)
  row <- data.frame(point_id = id, setting = setting, indicator = indicator, lat = lat, lon = lon,
                    coords = coords, model = model, clicked_age = clicked_age, age_young = age_young,
                    age_old = age_old, active = TRUE, note = note, stringsAsFactors = FALSE)
  .write_points_csv(rbind(u, row), points_file("user_points.csv", cfg))
  invisible(id)
}

#' exclude (or re-include) a database record
edit_record <- function(record_id, action = c("exclude", "include", "set_setting", "set_indicator"),
                        value = "", note = "", cfg = SILUCARB) {
  action <- match.arg(action)
  if (startsWith(record_id, "U:")) {          # own point: toggle its 'active' flag instead
    u <- read_user_points(cfg); i <- match(sub("^U:", "", record_id), u$point_id)
    if (is.na(i)) {
      stop("point ", record_id, " is not in user_points.csv (was the file changed outside the editor?)", call. = FALSE)
    }
    if (action %in% c("exclude", "include")) u$active[i] <- action == "include"
    if (action == "set_setting") u$setting[i] <- value
    if (action == "set_indicator") u$indicator[i] <- value
    .write_points_csv(u, points_file("user_points.csv", cfg))
    return(invisible(TRUE))
  }
  e <- read_record_edits(cfg)
  e <- rbind(e, data.frame(record_id = record_id, action = action, value = value, note = note))
  .write_points_csv(e, points_file("record_edits.csv", cfg))
  invisible(TRUE)
}

# ----------------------------------------------------------- localities ----
#' localities drawn (labelled) on the maps and usable as map centre:
#' points/localities.csv with columns name, lat, lon (present-day), show (TRUE/FALSE)
read_localities <- function(cfg = SILUCARB) {
  f <- points_file("localities.csv", cfg)
  if (!file.exists(f)) {
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    write.csv(cbind(cfg$sites, show = TRUE), f, row.names = FALSE)
  }
  l <- .read_points_csv(f, c("name", "lat", "lon"))
  if (!"show" %in% names(l)) l$show <- NA
  l$lat <- .as_num(l$lat)
  l$lon <- .as_num(l$lon)
  l$name <- trimws(ifelse(is.na(l$name), "", l$name))
  l$show <- .as_flag(l$show)
  l$show[is.na(l$show)] <- TRUE
  l[is.finite(l$lat) & is.finite(l$lon) & abs(l$lat) <= 90 & nzchar(l$name), c("name", "lat", "lon", "show")]
}
#' write the localities; returns TRUE, or FALSE with attribute "msg" when the file
#' cannot be written (usually because localities.csv is open in Excel or read-only)
write_localities <- function(l, cfg = SILUCARB) {
  f <- points_file("localities.csv", cfg)
  ok <- tryCatch(suppressWarnings({ write.csv(l, f, row.names = FALSE); TRUE }), error = function(e) FALSE)
  if (!ok) attr(ok, "msg") <- sprintf(paste0("Could not save %s - the file is locked. Close it in Excel (or untick ",
    "'Read-only' in its Properties) and try again; nothing was changed."), f)
  ok
}
#' message for a failed locality write, or NULL when it worked
locality_write_msg <- function(ok) if (isTRUE(ok)) NULL else attr(ok, "msg") %||% "Could not save localities.csv."
#' add a locality (present-day lat/lon); replaces one with the same name. Returns TRUE/FALSE (see write_localities)
add_locality <- function(name, lat, lon, show = TRUE, cfg = SILUCARB) {
  l <- read_localities(cfg); l <- l[l$name != name, ]
  write_localities(rbind(l, data.frame(name = name, lat = lat, lon = lon, show = show)), cfg)
}
remove_locality <- function(name, cfg = SILUCARB) {
  l <- read_localities(cfg); write_localities(l[!l$name %in% name, ], cfg)
}
