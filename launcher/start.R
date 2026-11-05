# start.R: starts PaleoMapR with this folder as working folder (started by the desktop shortcut) ####
#
# From R:  PaleoMapR::run_app("<this folder>")

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  dir <- if (length(f)) { dirname(normalizePath(sub("^--file=", "", f[1]))) } else { getwd() }
  assign("APP_DIR", normalizePath(dir, winslash = "/"), envir = globalenv())
})
lib_file <- file.path(APP_DIR, "state", "libpaths.txt")
if (file.exists(lib_file)) {
  .libPaths(c(readLines(lib_file, warn = FALSE), .libPaths()))
}
# PaleoMapR.launcher: lets "Update and restart" hand over to update.R in this folder
options(PaleoMapR.launcher = APP_DIR)
if (!file.exists(file.path(APP_DIR, "state", "setup_ok.txt"))) {
  cat("The data setup is not finished: the editor downloads the basic data first (about 400 MB).\n",
      "Running SETUP.bat again prepares everything, including the climate files and corrected grids.\n\n", sep = "")
}
PaleoMapR::run_app(folder = APP_DIR, port = 7793L, launch.browser = TRUE, stop_on_close = TRUE)
