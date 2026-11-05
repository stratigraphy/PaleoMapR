# update.R: started by "Update and restart" in PaleoMapR ####
#
# Waits until the editor has stopped (a package must not be replaced while R uses it), installs the newest version
# from GitHub and starts the editor again. Arguments: GitHub repositories, e.g. stratigraphy/PaleoMapR. The log is
# written to state/update_log.txt.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  dir <- if (length(f)) { dirname(normalizePath(sub("^--file=", "", f[1]))) } else { getwd() }
  assign("APP_DIR", dir, envir = globalenv())
})
PORT <- 7793L
repos <- commandArgs(trailingOnly = TRUE)
repos <- repos[grepl("^stratigraphy/[A-Za-z0-9.]+$", repos)]
if (!length(repos)) {
  repos <- "stratigraphy/PaleoMapR"
}
lib_file <- file.path(APP_DIR, "state", "libpaths.txt")
if (file.exists(lib_file)) {
  .libPaths(c(readLines(lib_file, warn = FALSE), .libPaths()))
}
log_file <- file.path(APP_DIR, "state", "update_log.txt")
cat(sprintf("Update, %s\n", format(Sys.time())), file = log_file)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_file, append = TRUE)
}

# wait (at most a minute) until the editor no longer answers on its port
port_open <- function() {
  con <- suppressWarnings(try(socketConnection("127.0.0.1", PORT, open = "r", blocking = TRUE, timeout = 1), silent = TRUE))
  if (inherits(con, "try-error")) {
    return(FALSE)
  }
  close(con)
  TRUE
}
say("Waiting for PaleoMapR to close ...")
k <- 0
while (port_open() && k < 60) {
  Sys.sleep(1)
  k <- k + 1
}
Sys.sleep(2)

ok <- TRUE
options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- .libPaths()[1]
for (r in repos) {
  say("Installing the newest ", basename(r), " from GitHub ...")
  res <- try(remotes::install_github(r, lib = lib, upgrade = "never", dependencies = NA), silent = TRUE)
  if (inherits(res, "try-error")) {
    ok <- FALSE
    say("FAILED: ", conditionMessage(attr(res, "condition")))
  }
}
# a fresh check is needed after an update
unlink(file.path(tools::R_user_dir("PaleoMapR", "cache"), "update_check.rds"))

vbs <- shQuote(normalizePath(file.path(APP_DIR, "start_app.vbs"), winslash = "\\"), type = "cmd")
if (ok) {
  say("Done. Starting PaleoMapR ...")
  system2("wscript", vbs, wait = FALSE)
  Sys.sleep(3)
} else {
  say("The update did not complete; PaleoMapR still has the previous version. Log: ", log_file)
  cat("Press Enter to close this window and start PaleoMapR ...\n")
  invisible(readLines(file("stdin"), n = 1L))
  system2("wscript", vbs, wait = FALSE)
}
