# Working folder and loading the editor code ####

#' The working folder of the editor
#'
#' The working folder holds the downloaded data (\code{data/}, about 1 GB after \code{\link{pm_setup}()}), the plate
#' models (\code{models/}), your points and edits (\code{points/}), the rendered maps (\code{maps/}) and the
#' statistics of the corrected grids (\code{results/}). The Windows launcher uses its own folder. From R the default
#' is the option \code{PaleoMapR.folder}, else the environment variable \code{PALEOMAPR_FOLDER}, else a folder
#' \code{PaleoMapR} in your home folder (on Windows: Documents).
#'
#' @return The path of the working folder (it is not created).
#' @examples
#' pm_folder()
#' @export
pm_folder <- function() {
  opt <- getOption("PaleoMapR.folder")
  if (!is.null(opt) && nzchar(opt)) {
    return(normalizePath(opt, winslash = "/", mustWork = FALSE))
  }
  env <- Sys.getenv("PALEOMAPR_FOLDER")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(path.expand("~"), "PaleoMapR"), winslash = "/", mustWork = FALSE)
}

# The installed editor code: inst/app with app.R, R/ (the scripts) and python/
pm_code_dir <- function() {
  d <- system.file("app", package = "PaleoMapR")
  if (!nzchar(d)) {
    stop("The editor code was not found; reinstall PaleoMapR.", call. = FALSE)
  }
  normalizePath(d, winslash = "/")
}

#' Load the functions behind the editor
#'
#' Sources the editor scripts into an environment and points their settings (\code{SILUCARB}) to a working
#' folder, so that the steps of the editor can be run from an R script. The working directory is restored
#' afterwards.
#'
#' @param folder The working folder (see \code{\link{pm_folder}()}). It is created if needed.
#' @param envir The environment to load into; a new environment by default.
#' @return The environment (invisibly) with all functions and the settings list \code{SILUCARB}.
#' @examples
#' \dontrun{
#' pm <- pm_load("D:/Maps/PaleoMapR")
#' pm$run_dem_correction(425)                                # corrected grid for 425 Ma
#' pm$make_maps(425, layer = "precipitation", wind = TRUE)   # PNG and PDF in maps/
#' }
#' @export
pm_load <- function(folder = pm_folder(), envir = new.env(parent = globalenv())) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)
  code <- pm_code_dir()
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(folder)
  assign("PROJECT_DIR", folder, envir = envir)
  assign("CODE_DIR", code, envir = envir)
  files <- sort(list.files(file.path(code, "R"), "\\.R$", full.names = TRUE))
  for (f in files) {
    exprs <- parse(file = f, keep.source = FALSE, encoding = "UTF-8")
    for (k in seq_along(exprs)) {
      eval(exprs[[k]], envir = envir)
    }
  }
  cfg <- get("SILUCARB", envir = envir)
  cfg$dir_R <- file.path(code, "R")
  cfg$dir_python <- file.path(code, "python")
  cfg <- envir$absolutize_config(cfg, folder)
  cfg <- envir$apply_user_settings(cfg)
  assign("SILUCARB", cfg, envir = envir)
  invisible(envir)
}
