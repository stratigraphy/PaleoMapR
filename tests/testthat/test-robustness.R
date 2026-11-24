pm_test_env <- function() {
  skip_on_cran()
  folder <- file.path(tempdir(), paste0("pm_rob_", sample.int(1e6, 1)))
  pm <- pm_load(folder)
  pm$init_user_tables()
  pm
}

test_that("points tables saved by Excel (BOM, semicolons, decimal commas, local TRUE) are read", {
  pm <- pm_test_env()
  f <- pm$points_file("user_points.csv")
  con <- file(f, "wb")
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeLines(c("point_id;setting;indicator;lat;lon;coords;model;clicked_age;age_young;age_old;active;note",
      "a1;Marine;;50,52;5,24;present;scotese_and_wright2018;;440;430;VRAI;Huy",
      "a1;terrestrial;;51;6;present;scotese_and_wright2018;;425;425;WAAR;same id"), con)
  close(con)
  u <- suppressWarnings(pm$read_user_points())
  expect_equal(u$point_id, c("a1", "a1_1"))
  expect_equal(u$setting, c("marine", "terrestrial"))
  expect_equal(u$lat, c(50.52, 51))
  expect_equal(u$age_young[1], 430)
  expect_equal(u$age_old[1], 440)
  expect_true(all(u$active))
})

test_that("a table without the needed columns is never written back", {
  pm <- pm_test_env()
  f <- pm$points_file("user_points.csv")
  writeLines(c("a;b", "1;2"), f)
  expect_error(pm$add_user_point(1, 1, age_young = 425), "not found")
  expect_equal(readLines(f), c("a;b", "1;2"))
})

test_that("empty and incomplete edit tables are read", {
  pm <- pm_test_env()
  f <- pm$points_file("record_edits.csv")
  writeLines(character(0), f)
  expect_equal(nrow(pm$read_record_edits()), 0)
  writeLines(c("record_id,action,value,note", "PBDB:1,,,", "PBDB:2,Exclude,,"), f)
  e <- suppressWarnings(pm$read_record_edits())
  expect_equal(e$record_id, "PBDB:2")
  expect_equal(e$action, "exclude")
})

test_that("grid distances are great-circle distances at every longitude and latitude", {
  pm <- pm_test_env()
  lat <- seq(-90, 90, 0.2)
  km_per_deg <- pm$KM_PER_DEG
  for (la in c(0, 30, 60)) {
    i <- (la + 90) / 0.2 + 1
    for (j in c(10, 1000, 2700)) {
      p <- pm$.to_km(c(i, i + 1), c(j, j), lat)
      expect_equal(sqrt(sum(diff(p)^2)), 0.2 * km_per_deg, tolerance = 1e-3)
      q <- pm$.to_km(c(i, i), c(j, j + 1), lat)
      expect_equal(sqrt(sum(diff(q)^2)), 0.2 * km_per_deg * cos(la * pi / 180), tolerance = 1e-3)
    }
  }
})

test_that("colour ramps work for no values, file names are safe, coarse fields extend over gaps", {
  pm <- pm_test_env()
  m <- pm$colourise(matrix(-100, 3, 4))
  expect_true(is.numeric(m))
  expect_equal(dim(m), c(3, 4, 3))
  expect_equal(pm$.file_safe("Gotland/Visby: \"N\""), "Gotland-Visby-N")
  z <- matrix(c(NA, 1, 2, NA, NA, NA), 2, 3)
  ze <- pm$.extend_field(z, 1)
  expect_true(all(is.finite(ze[1, ])))
})

test_that("downloads are kept only when complete", {
  pm <- pm_test_env()
  d <- tempfile()
  dir.create(d)
  html <- file.path(d, "a.zip.part")
  writeLines("<html>proxy error</html>", html)
  expect_false(pm$.download_ok(html))
  csv <- file.path(d, "b.csv")
  writeLines(c("\"collection_no\",\"lat\"", "1,2"), csv)
  expect_true(pm$.download_ok(csv))
})
