test_that("pm_folder follows the option and the environment variable", {
  oldop <- options(PaleoMapR.folder = file.path(tempdir(), "pm_a"))
  on.exit(options(oldop), add = TRUE)
  expect_equal(basename(pm_folder()), "pm_a")
  options(PaleoMapR.folder = NULL)
  withr_env <- Sys.getenv("PALEOMAPR_FOLDER", unset = NA)
  Sys.setenv(PALEOMAPR_FOLDER = file.path(tempdir(), "pm_b"))
  on.exit(if (is.na(withr_env)) { Sys.unsetenv("PALEOMAPR_FOLDER") } else { Sys.setenv(PALEOMAPR_FOLDER = withr_env) },
      add = TRUE)
  expect_equal(basename(pm_folder()), "pm_b")
})

test_that("the editor scripts and the app parse", {
  code <- system.file("app", package = "PaleoMapR")
  files <- c(list.files(file.path(code, "R"), "\\.R$", full.names = TRUE), file.path(code, "app.R"))
  expect_true(length(files) >= 16)
  for (f in files) {
    expect_no_error(parse(file = f, keep.source = FALSE, encoding = "UTF-8"))
  }
  expect_true(file.exists(file.path(code, "python", "paleomap_gplates.py")))
})

test_that("pm_load points the settings to the working folder and the code to the package", {
  skip_on_cran()
  folder <- file.path(tempdir(), "pm_load_test")
  wd <- getwd()
  pm <- pm_load(folder)
  expect_equal(getwd(), wd)
  cfg <- pm$SILUCARB
  root <- normalizePath(folder, winslash = "/")
  expect_equal(cfg$dir_points, file.path(root, "points"))
  expect_equal(cfg$dir_dem, file.path(root, "data/paleoDEM_0.2Deg_grids"))
  expect_equal(cfg$dir_R, file.path(system.file("app", package = "PaleoMapR"), "R"))
  expect_true(file.exists(file.path(cfg$dir_python, "paleomap_gplates.py")))
  for (fn in c("run_dem_correction", "make_maps", "download_data", "gp_init", "all_records", "climate_files")) {
    expect_true(is.function(pm[[fn]]), info = fn)
  }
  # the user tables are written into the working folder
  pm$init_user_tables()
  invisible(pm$read_localities())
  expect_true(file.exists(file.path(root, "points", "user_points.csv")))
  expect_true(file.exists(file.path(root, "points", "localities.csv")))
})

test_that("the update helpers report the origin and the update code", {
  o <- PaleoMapR:::pm_pkg_origin("PaleoMapR")
  expect_true(o$origin %in% c("github", "cran", "local", "none"))
  st <- data.frame(package = "PaleoMapR", repo = "stratigraphy/PaleoMapR", update = TRUE, stringsAsFactors = FALSE)
  code <- PaleoMapR:::pm_update_code(st)
  expect_match(code[2], "remotes::install_github\\(\"stratigraphy/PaleoMapR\"")
})

test_that("the About tab names the version and the data sources", {
  html <- as.character(PaleoMapR:::pm_about_tab())
  expect_match(html, as.character(utils::packageVersion("PaleoMapR")), fixed = TRUE)
  expect_match(html, "10.5194/cp-17-1483-2021", fixed = TRUE)
})
