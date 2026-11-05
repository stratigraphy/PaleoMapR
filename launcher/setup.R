# setup.R: installs PaleoMapR, makes the desktop shortcut and prepares the data (run it with SETUP.bat) ####
#
# 1. Installs the newest PaleoMapR from GitHub with the CRAN packages it needs.
# 2. Writes the hidden launcher start_app.vbs and the desktop shortcut "PaleoMapR".
# 3. Runs PaleoMapR::pm_setup() with this folder as working folder: Python with pyGPlates, the open data sets
#    (about 1 GB), the HadCM3L climate files and the corrected grids of the Silurian ages (20-60 minutes).
# Arguments of SETUP.bat: "all" (all 109 map ages, several hours), map ages such as 300 305 310, or "nodata"
# (install and shortcut only). Everything is also written to state/setup_log.txt.
# Developers: a file local_source.txt in this folder with the source folder of PaleoMapR installs that instead of
# the GitHub version.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  dir <- if (length(f)) { dirname(normalizePath(sub("^--file=", "", f[1]))) } else { getwd() }
  assign("APP_DIR", dir, envir = globalenv())
})

APP_TITLE <- "PaleoMapR"
PORT <- 7793L
REPO <- "stratigraphy/PaleoMapR"
CRAN <- "https://cloud.r-project.org"

state_dir <- file.path(APP_DIR, "state")
dir.create(state_dir, showWarnings = FALSE)
log_file <- file.path(state_dir, "setup_log.txt")
cat(sprintf("Setup of %s, %s\n", APP_TITLE, format(Sys.time())), file = log_file)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_file, append = TRUE)
}
win <- function(p) {
  normalizePath(p, winslash = "\\", mustWork = FALSE)
}
run_ps <- function(lines) {
  f <- tempfile(fileext = ".ps1")
  on.exit(unlink(f), add = TRUE)
  writeLines(lines, f)
  suppressWarnings(system2("powershell", c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", shQuote(f)),
      stdout = TRUE, stderr = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
ages <- "config"
do_data <- TRUE
if (length(args)) {
  if (any(tolower(args) == "nodata")) {
    do_data <- FALSE
  } else if (any(tolower(args) == "all")) {
    ages <- "all"
  } else {
    num <- suppressWarnings(as.numeric(args))
    if (any(is.finite(num))) {
      ages <- num[is.finite(num)]
    }
  }
}

# A user library when the R library is not writable (R installed for all users without admin rights) ####
lib <- .libPaths()[1]
if (file.access(lib, 2) != 0) {
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  say("Packages go to the user library ", lib)
}

# Packages ####
options(repos = c(CRAN = CRAN))
if (!requireNamespace("remotes", quietly = TRUE)) {
  say("Installing remotes from CRAN ...")
  utils::install.packages("remotes", lib = lib)
}
if (isNamespaceLoaded("PaleoMapR")) {
  unloadNamespace("PaleoMapR")
}
local_src <- file.path(APP_DIR, "local_source.txt")
if (file.exists(local_src)) {
  src <- trimws(readLines(local_src, warn = FALSE))
  src <- src[nzchar(src)]
  for (s in src) {
    say("Installing the packages PaleoMapR needs, then PaleoMapR from the source folder ", s, " ...")
    remotes::install_deps(s, lib = lib, upgrade = "never")
    utils::install.packages(s, lib = lib, repos = NULL, type = "source")
  }
} else {
  say("Installing the newest PaleoMapR from GitHub (", REPO, ") and the packages it needs (a few minutes the first time) ...")
  res <- try(remotes::install_github(REPO, lib = lib, upgrade = "never", dependencies = NA), silent = TRUE)
  if (inherits(res, "try-error")) {
    say("FAILED: ", conditionMessage(attr(res, "condition")))
  }
}
if (!requireNamespace("PaleoMapR", quietly = TRUE)) {
  say("PaleoMapR could not be installed; see the messages above and ", win(log_file), ".")
  quit(save = "no", status = 1)
}
say("PaleoMapR ", as.character(utils::packageVersion("PaleoMapR")), ", shiny ", as.character(utils::packageVersion("shiny")),
    ", terra ", as.character(utils::packageVersion("terra")), ", reticulate ", as.character(utils::packageVersion("reticulate")))

writeLines(.libPaths(), file.path(state_dir, "libpaths.txt"))
rscript <- win(file.path(R.home("bin"), "Rscript.exe"))
say("R: ", R.version.string, " (", rscript, ")")

# Launcher: opens the page if the editor already runs; else starts R without a window, or with a window while the
# data setup is not finished (the first start then downloads the basic data and shows its progress) ####
vbs <- file.path(APP_DIR, "start_app.vbs")
start_r <- win(file.path(APP_DIR, "start.R"))
ok_file <- win(file.path(state_dir, "setup_ok.txt"))
con <- file(vbs, "wb")
writeLines(c(
    'Set sh = CreateObject("WScript.Shell")',
    'Set fs = CreateObject("Scripting.FileSystemObject")',
    sprintf('rc = sh.Run("cmd /c netstat -ano | findstr /R /C:"":%d .*LISTENING"" >nul", 0, True)', PORT),
    "If rc = 0 Then",
    sprintf('  sh.Run "http://127.0.0.1:%d", 1, False', PORT),
    sprintf('ElseIf fs.FileExists("%s") Then', ok_file),
    sprintf('  sh.Run """%s"" ""%s""", 0, False', rscript, start_r),
    "Else",
    sprintf('  sh.Run """%s"" ""%s""", 1, False', rscript, start_r),
    "End If"), con, sep = "\r\n")
close(con)

# Desktop shortcut ####
say(paste(run_ps(c(
    "$desktop = [Environment]::GetFolderPath('Desktop')",
    "$w = New-Object -ComObject WScript.Shell",
    sprintf("$s = $w.CreateShortcut((Join-Path $desktop '%s.lnk'))", APP_TITLE),
    "$s.TargetPath = 'wscript.exe'",
    sprintf("$s.Arguments = '\"%s\"'", win(vbs)),
    sprintf("$s.WorkingDirectory = '%s'", win(APP_DIR)),
    sprintf("$s.IconLocation = '%s,0'", win(file.path(APP_DIR, "paleomapr.ico"))),
    "$s.Description = 'PaleoMapR: fossil-corrected paleogeographic maps'",
    "$s.Save()",
    sprintf("Write-Output 'Desktop shortcut written: %s.lnk'", APP_TITLE))), collapse = "\n"))

# Data, climate files and corrected grids (this folder is the working folder) ####
if (do_data) {
  say("\nPreparing the data in ", win(APP_DIR), " (20-60 minutes the first time; running SETUP.bat again continues where it stopped) ...")
  # the messages of pm_setup() also go to the log
  to_log <- function(m) {
    cat(conditionMessage(m), file = log_file, append = TRUE)
  }
  ok <- tryCatch(withCallingHandlers(PaleoMapR::pm_setup(APP_DIR, ages = ages), message = to_log), error = function(e) {
    say("FAILED: ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) {
    writeLines(format(Sys.time()), file.path(state_dir, "setup_ok.txt"))
    say("Done. Start PaleoMapR with the desktop icon. Log: ", win(log_file))
  } else {
    say("The data setup is not complete. Run SETUP.bat again (finished parts are kept). Log: ", win(log_file))
  }
} else {
  say("Done (without the data). The first start downloads the basic data. Log: ", win(log_file))
}
