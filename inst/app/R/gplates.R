# =============================================================================
#  pyGPlates through reticulate: plate rotations only.
#  Functions return data.frames / lists of plain R vectors.
# =============================================================================
library(reticulate)

.gp <- new.env()
PY_PKGS <- c("pygplates>=1.0", "plate-model-manager", "numpy")

#' the Python environment that PaleoMapR::pm_setup() makes (a plain virtual environment, created with
#' 'uv'; no conda). It lives in your R user data folder, so all copies of the editor share it.
py_venv_dir <- function() file.path(tools::R_user_dir("paleomapeditor", "data"), "venv")
py_venv_python <- function(v = py_venv_dir())
  if (.Platform$OS.type == "windows") file.path(v, "Scripts", "python.exe") else file.path(v, "bin", "python")

#' the 'uv' program (downloads Python and packages): the copy reticulate keeps, or a
#' standalone copy downloaded once from GitHub
uv_path <- function() {
  b <- tryCatch(reticulate:::uv_binary(), error = function(e) NULL)
  if (length(b) && nzchar(b) && file.exists(b)) return(b)
  win <- .Platform$OS.type == "windows"
  d <- file.path(tools::R_user_dir("paleomapeditor", "data"), "uv")
  f <- list.files(d, if (win) "^uv\\.exe$" else "^uv$", recursive = TRUE, full.names = TRUE)
  if (length(f)) return(f[1])
  arch <- if (grepl("arm|aarch64", R.version$arch)) "aarch64" else "x86_64"
  asset <- if (win) "uv-x86_64-pc-windows-msvc.zip" else if (Sys.info()[["sysname"]] == "Darwin")
    sprintf("uv-%s-apple-darwin.tar.gz", arch) else sprintf("uv-%s-unknown-linux-gnu.tar.gz", arch)
  dir.create(d, recursive = TRUE, showWarnings = FALSE); z <- file.path(d, asset)
  download.file(paste0("https://github.com/astral-sh/uv/releases/latest/download/", asset), z, mode = "wb", quiet = TRUE)
  if (win) unzip(z, exdir = d) else untar(z, exdir = d)
  f <- list.files(d, if (win) "^uv\\.exe$" else "^uv$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("could not get the 'uv' program", call. = FALSE)
  if (!win) Sys.chmod(f[1], "755")
  f[1]
}

#' make (or update) the virtual environment with pyGPlates: uv downloads Python 3.12 itself
make_py_venv <- function() {
  uv <- uv_path(); v <- py_venv_dir(); py <- py_venv_python(v)
  if (!file.exists(py)) {
    dir.create(dirname(v), recursive = TRUE, showWarnings = FALSE)
    out <- system2(uv, c("venv", shQuote(v), "--python", "3.12"), stdout = TRUE, stderr = TRUE)
    if (!file.exists(py)) stop("uv venv failed:\n", paste(tail(out, 10), collapse = "\n"), call. = FALSE)
  }
  out <- system2(uv, c("pip", "install", "--python", shQuote(py), shQuote(PY_PKGS)), stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("installing the Python packages failed:\n", paste(tail(out, 10), collapse = "\n"), call. = FALSE)
  invisible(py)
}

#' which Python to use: the environment made by pm_setup(); otherwise a conda environment
#' "paleomap" (or "silucarb", from older versions); otherwise reticulate's own managed Python
gp_choose_python <- function(cfg = SILUCARB) {
  py <- py_venv_python()
  if (file.exists(py)) return(invisible(use_python(py, required = TRUE)))
  envs <- tryCatch(conda_list()$name, error = function(e) character(0))
  for (e in unique(c(cfg$conda_env, "silucarb"))) if (e %in% envs) return(invisible(use_condaenv(e, required = TRUE)))
  if (tryCatch(virtualenv_exists(cfg$conda_env), error = function(e) FALSE))
    return(invisible(use_virtualenv(cfg$conda_env, required = TRUE)))
  if (utils::packageVersion("reticulate") >= "1.41")
    return(invisible(py_require(PY_PKGS, python_version = ">=3.10,<3.13")))
  message("No Python environment for pyGPlates found; using the default Python. Run PaleoMapR::pm_setup() if pyGPlates is missing.")
}

gp_init <- function(cfg = SILUCARB) {
  if (!is.null(.gp$mod)) return(invisible(.gp$mod))
  if (!py_available(initialize = FALSE)) gp_choose_python(cfg)
  .gp$mod <- import_from_path("paleomap_gplates", path = cfg$dir_python, convert = TRUE)
  .gp$mod$init(normalizePath(cfg$dir_models, mustWork = FALSE))
  invisible(.gp$mod)
}

#' present-day points -> palaeo-positions (each point at its own age)
gp_reconstruct <- function(lat, lon, ages, model = "scotese_and_wright2018") {
  if (!length(lat)) return(data.frame(plat = numeric(0), plon = numeric(0), plate_id = integer(0)))
  r <- gp_init()$reconstruct(as.numeric(lat), as.numeric(lon), as.numeric(ages), model)
  data.frame(plat = as.numeric(r$lat), plon = as.numeric(r$lon), plate_id = as.integer(r$plate_id))
}

#' palaeo-positions at `age` -> present day (for points you click on a palaeo map)
gp_to_present <- function(plat, plon, age, model = "scotese_and_wright2018") {
  if (!length(plat)) return(data.frame(lat = numeric(0), lon = numeric(0), plate_id = integer(0)))
  r <- gp_init()$to_present(as.numeric(plat), as.numeric(plon), as.numeric(age), model)
  data.frame(lat = as.numeric(r$lat), lon = as.numeric(r$lon), plate_id = as.integer(r$plate_id))
}

#' palaeo-positions in one plate model -> the same crust in another model, same age
gp_transfer <- function(plat, plon, age, source = "scotese_and_wright2018", target, valid_only = TRUE) {
  r <- gp_init()$transfer(as.numeric(plat), as.numeric(plon), as.numeric(age), source, target, valid_only)
  data.frame(lat = as.numeric(r$lat), lon = as.numeric(r$lon))
}

#' present-day lines (shapefile) rotated to `age`: list of data.frames (lat, lon) with attr young
gp_modern_lines <- function(shapefile, model, age) {
  key <- paste(shapefile, model, age)
  if (!is.null(.gp$lines[[key]])) return(.gp$lines[[key]])
  r <- gp_init()$modern_lines(normalizePath(shapefile), model, as.numeric(age))
  out <- lapply(r, function(p) structure(data.frame(lat = as.numeric(p$lat), lon = as.numeric(p$lon)),
                                         young = isTRUE(p$young)))
  if (is.null(.gp$lines)) .gp$lines <- list()
  .gp$lines[[key]] <- out
  out
}

gp_plate_info <- function(lat, lon, model = "scotese_and_wright2018") {
  r <- gp_init()$plate_info(as.numeric(lat), as.numeric(lon), model)
  data.frame(plate_id = as.integer(r$plate_id), name = unlist(r$name), valid_from = as.numeric(r$valid_from))
}
