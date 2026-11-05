# Update check ####
#
# PaleoMapR is installed from GitHub. Once a day (in a background process) the installed commit is compared with
# the newest commit on GitHub; when that is newer, the editor offers to update. With the desktop launcher the
# update is installed after the editor has closed (a package must not be replaced while R uses it) and the editor
# is started again; without the launcher the R code for the update is shown.

pm_github_repo <- function() {
  "stratigraphy/PaleoMapR"
}

# Where the installed package came from: "github" (with its commit), "cran", "local" (a development install from a
# source folder) or "none"
pm_pkg_origin <- function(pkg = "PaleoMapR") {
  d <- suppressWarnings(utils::packageDescription(pkg))
  if (!is.list(d)) {
    return(list(origin = "none", sha = NA_character_, version = NA_character_))
  }
  if (!is.null(d$RemoteSha) && identical(tolower(as.character(d$RemoteType)), "github")) {
    return(list(origin = "github", sha = d$RemoteSha, version = d$Version))
  }
  if (identical(d$Repository, "CRAN")) {
    return(list(origin = "cran", sha = NA_character_, version = d$Version))
  }
  list(origin = "local", sha = NA_character_, version = d$Version)
}

# The newest commit of the default branch of a GitHub repository (NA when GitHub cannot be reached)
pm_github_head <- function(repo, timeout = 8) {
  old <- options(timeout = timeout)
  on.exit(options(old), add = TRUE)
  txt <- tryCatch(suppressWarnings(readLines(sprintf("https://api.github.com/repos/%s/commits/HEAD", repo), warn = FALSE)),
      error = function(e) { NULL })
  if (is.null(txt)) {
    return(NA_character_)
  }
  txt <- paste(txt, collapse = "")
  m <- regmatches(txt, regexpr("\"sha\" *: *\"[0-9a-f]{40}\"", txt))
  if (!length(m)) {
    return(NA_character_)
  }
  sub("^.*\"([0-9a-f]{40})\"$", "\\1", m)
}

#' Check GitHub for a newer version of PaleoMapR
#'
#' Compares the installed commit with the newest commit on GitHub. The editor does this once a day by itself and
#' offers "Update and restart".
#'
#' @param force Check now instead of using the result of the last check (kept for a day).
#' @return A one-row data frame: package, repository, origin of the installation, version, installed and newest
#'   commit, and whether an update is available.
#' @examples
#' \dontrun{
#' pm_update_status(force = TRUE)
#' }
#' @export
pm_update_status <- function(force = FALSE) {
  cache <- file.path(tools::R_user_dir("PaleoMapR", "cache"), "update_check.rds")
  if (!force && file.exists(cache)) {
    old <- try(readRDS(cache), silent = TRUE)
    if (!inherits(old, "try-error") && difftime(Sys.time(), old$time, units = "hours") < 24) {
      return(old$status)
    }
  }
  repo <- pm_github_repo()
  o <- pm_pkg_origin("PaleoMapR")
  head <- if (identical(o$origin, "local")) { NA_character_ } else { pm_github_head(repo) }
  upd <- !is.na(head) && (o$origin %in% c("cran", "none") || (identical(o$origin, "github") && !identical(o$sha, head)))
  status <- data.frame(package = "PaleoMapR", repo = repo, origin = o$origin, version = o$version,
      installed = o$sha, latest = head, update = upd, stringsAsFactors = FALSE)
  if (!is.na(head)) {
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    try(saveRDS(list(time = Sys.time(), status = status), cache), silent = TRUE)
  }
  status
}

pm_update_code <- function(status) {
  c("# close the editor first, then in a fresh R session:",
      sprintf("remotes::install_github(\"%s\", upgrade = \"never\")", status$repo[1]))
}

# Starts the launcher's updater (it waits until the editor has stopped, installs, and starts it again)
pm_start_updater <- function(status) {
  dir <- getOption("PaleoMapR.launcher")
  upd <- if (is.null(dir)) { "" } else { file.path(dir, "update.R") }
  if (.Platform$OS.type != "windows" || !nzchar(upd) || !file.exists(upd)) {
    return(FALSE)
  }
  rscript <- normalizePath(file.path(R.home("bin"), "Rscript.exe"), winslash = "\\", mustWork = FALSE)
  args <- c("/c", "start", "\"PaleoMapR update\"", shQuote(rscript, type = "cmd"),
      shQuote(normalizePath(upd, winslash = "\\", mustWork = FALSE), type = "cmd"), status$repo[status$update])
  res <- try(system2("cmd", args, wait = FALSE), silent = TRUE)
  !inherits(res, "try-error")
}
