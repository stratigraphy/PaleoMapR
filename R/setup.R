# One-time setup: Python, data, climate files and corrected grids ####

#' Set up the editor (once per working folder)
#'
#' Prepares a working folder in four steps. It needs internet and takes 20 to 60 minutes the first time, mostly
#' downloads (about 1 GB). Running it again is safe: finished steps and files are kept, so after an interrupted
#' download simply run it again.
#' \enumerate{
#'   \item Python with pyGPlates (plate rotations): a small Python environment is made with 'uv', which downloads
#'     Python 3.12 itself; no Python or conda installation is needed. It is kept in your R user data folder and
#'     shared by all working folders. If that fails, a conda environment "paleomap" is tried.
#'   \item The open data sets into \code{data/}: PaleoDEMs 0-540 Ma (Scotese and Wright 2018), PaleoCoastlines
#'     (Kocsis and Scotese 2021), PBDB collections, Macrostrat units and the sea-level curve of van der Meer et al.
#'     (2022); then the table of database records.
#'   \item The HadCM3L climate files (Valdes et al. 2021) of the chosen ages, about 30 MB per age.
#'   \item The fossil-corrected grids of the chosen ages (1-2 minutes per age per core).
#' }
#'
#' @param folder The working folder (see \code{\link{pm_folder}()}).
#' @param ages The map ages to prepare: \code{"config"} (the Silurian ages 445-415 Ma), \code{"all"} (all 109
#'   PaleoDEM ages, 0-540 Ma: several hours and about 3 GB of climate files), or ages in Ma such as
#'   \code{c(300, 305, 310)}, each matched to the nearest PaleoDEM age. Other ages work in the editor too: choose
#'   them and press "Rebuild this age".
#' @param workers The number of R processes for the corrected grids; by default up to 8, one less than the cores.
#' @return \code{TRUE} (invisibly) when everything is in place, \code{FALSE} when a part is missing.
#' @examples
#' \dontrun{
#' pm_setup()                                   # the Silurian ages, in pm_folder()
#' pm_setup("D:/Maps/PaleoMapR", ages = seq(250, 200, by = -5))
#' }
#' @export
pm_setup <- function(folder = pm_folder(), ages = "config", workers = NULL) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)
  message("PaleoMapR ", as.character(utils::packageVersion("PaleoMapR")), ": setting up ", folder)
  pm <- pm_load(folder)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(folder)
  cfg <- pm$SILUCARB

  # Python with pyGPlates
  message("\n== 1/4  Python with pyGPlates ==")
  py_ok <- tryCatch({
    pm$make_py_venv()
    pm_python_test(folder)
  }, error = function(e) {
    message("  uv environment FAILED: ", conditionMessage(e))
    FALSE
  })
  if (!py_ok) {
    py_ok <- pm_conda_fallback(pm, folder)
  }
  if (!py_ok) {
    stop("Python with pyGPlates could not be set up; see the messages above. ",
        "Restart R and run pm_setup() again, or see https://stratigraphy.github.io/PaleoMapR/articles/getting-started.html#problems",
        call. = FALSE)
  }
  message("pyGPlates works.")

  # Data
  message("\n== 2/4  data (PaleoDEMs, PaleoCoastlines, PBDB, Macrostrat, sea level; ~600 MB) ==")
  pm$download_data()
  pm$init_user_tables()
  invisible(pm$read_localities())
  if (!file.exists(pm$points_file("records_base.csv"))) {
    pm$build_base_records()
  }
  # the DEM zip is not needed once it is unpacked
  if (length(list.files(cfg$dir_dem, "Ma\\.nc$")) >= 80) {
    unlink(file.path(cfg$dir_data, c("paleoDEM_0.2Deg_grids.zip", "PaleoDEMs_1deg_nc.zip")))
  }
  all_ages <- pm$available_ages(cfg)
  if (identical(ages, "config")) {
    ages <- cfg$ages
  } else if (identical(ages, "all")) {
    ages <- all_ages
  } else {
    a <- as.numeric(ages)
    out <- numeric(length(a))
    for (i in seq_along(a)) {
      out[i] <- all_ages[which.min(abs(all_ages - a[i]))]
    }
    ages <- out
  }
  ages <- sort(unique(ages), decreasing = TRUE)
  message("\nAges to prepare: ", paste(ages, collapse = ", "), " Ma")

  # HadCM3L climate files
  message("\n== 3/4  HadCM3L climate (Valdes et al. 2021) for ", length(ages), " age(s) ==")
  clim_ok <- TRUE
  for (a in ages) {
    e <- pm$climate_exp(a)
    message(sprintf("  %s Ma: run %s (%s Ma)", a, e$id, e$age))
    for (s in names(pm$CLIMATE_SEASONS)) {
      ok <- tryCatch({
        pm$climate_files(a, s)
        TRUE
      }, error = function(err) {
        message("    ", conditionMessage(err))
        FALSE
      })
      clim_ok <- isTRUE(ok) && clim_ok
    }
  }

  # Corrected grids
  message("\n== 4/4  fossil-corrected grids for ", length(ages), " age(s) (~1-2 min per age per core) ==")
  todo <- ages[!pm_grid_done(cfg, ages)]
  if (length(todo)) {
    if (is.null(workers)) {
      workers <- max(1, min(length(todo), parallel::detectCores() - 1, 8))
    }
    tryCatch(pm$run_dem_correction(todo, cfg, workers = workers), error = function(e) {
      message("  ", conditionMessage(e), "\n  (you can also build them later in the editor with 'Rebuild this age')")
    })
  } else {
    message("  already there")
  }

  need <- c(PaleoDEMs = length(list.files(cfg$dir_dem, "Ma\\.nc$")) >= 80,
      PaleoCoastlines = dir.exists(file.path(cfg$dir_coast, "CS")),
      records = file.exists(pm$points_file("records_base.csv")),
      climate = clim_ok,
      corrected_grids = all(pm_grid_done(cfg, ages)))
  if (!all(need)) {
    message("\nNot complete yet: ", paste(names(need)[!need], collapse = ", "),
        ". Check the internet connection and run the setup again (finished steps and downloads are kept).")
    return(invisible(FALSE))
  }
  message("\nSetup finished. Start the editor with PaleoMapR::run_app(\"", folder, "\").")
  invisible(TRUE)
}

# Which ages have a corrected grid
pm_grid_done <- function(cfg, ages) {
  done <- logical(length(ages))
  for (i in seq_along(ages)) {
    done[i] <- length(list.files(cfg$dir_corrected, paste0("_", ages[i], "Ma\\.nc$"))) > 0
  }
  done
}

# Tests pyGPlates in a separate R process, so a failed try does not block the next one
pm_python_test <- function(folder) {
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") { "Rscript.exe" } else { "Rscript" })
  tf <- tempfile(fileext = ".R")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
      sprintf("pm <- PaleoMapR::pm_load('%s')", folder),
      "m <- pm$gp_init()",
      "p <- reticulate::import('pygplates')",
      "cat('PYGPLATES_OK', reticulate::py_config()$python, '\\n')"), tf)
  out <- suppressWarnings(system2(rs, shQuote(tf), stdout = TRUE, stderr = TRUE))
  ok <- any(grepl("PYGPLATES_OK", out))
  out <- out[!grepl("^\\s*$", out)]
  message(paste(utils::tail(out, if (ok) { 1 } else { 15 }), collapse = "\n"))
  ok
}

# Second try: a conda environment "paleomap" (miniconda is installed when needed)
pm_conda_fallback <- function(pm, folder) {
  message("trying a conda environment 'paleomap' instead ...")
  unlink(pm$py_venv_dir(), recursive = TRUE)
  Sys.setenv(CONDA_PLUGINS_AUTO_ACCEPT_TOS = "yes")
  tryCatch({
    if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) { NULL }))) {
      message("installing miniconda (once, ~5 minutes) ...")
      tryCatch(reticulate::install_miniconda(), error = function(e) { message("  ", conditionMessage(e)) })
    }
    cb <- reticulate::conda_binary()
    for (ch in paste0("https://repo.anaconda.com/pkgs/", c("main", "r", "msys2"))) {
      suppressWarnings(system2(cb, c("tos", "accept", "--override-channels", "--channel", ch), stdout = FALSE, stderr = FALSE))
    }
    if (!"paleomap" %in% reticulate::conda_list()$name) {
      reticulate::conda_create("paleomap", python_version = "3.11", forge = TRUE, channel = "conda-forge")
    }
    reticulate::py_install(pm$PY_PKGS, envname = "paleomap", pip = TRUE)
    pm_python_test(folder)
  }, error = function(e) {
    message("  conda FAILED: ", conditionMessage(e))
    FALSE
  })
}
