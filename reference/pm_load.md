# Load the functions behind the editor

Sources the editor scripts into an environment and points their settings
(`SILUCARB`) to a working folder, so that the steps of the editor can be
run from an R script. The working directory is restored afterwards.

## Usage

``` r
pm_load(folder = pm_folder(), envir = new.env(parent = globalenv()))
```

## Arguments

- folder:

  The working folder (see
  [`pm_folder()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_folder.md)).
  It is created if needed.

- envir:

  The environment to load into; a new environment by default.

## Value

The environment (invisibly) with all functions and the settings list
`SILUCARB`.

## Examples

``` r
if (FALSE) { # \dontrun{
pm <- pm_load("D:/Maps/PaleoMapR")
pm$run_dem_correction(425)                                # corrected grid for 425 Ma
pm$make_maps(425, layer = "precipitation", wind = TRUE)   # PNG and PDF in maps/
} # }
```
