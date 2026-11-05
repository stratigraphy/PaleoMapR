# The working folder of the editor

The working folder holds the downloaded data (`data/`, about 1 GB after
[`pm_setup()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_setup.md)),
the plate models (`models/`), your points and edits (`points/`), the
rendered maps (`maps/`) and the statistics of the corrected grids
(`results/`). The Windows launcher uses its own folder. From R the
default is the option `PaleoMapR.folder`, else the environment variable
`PALEOMAPR_FOLDER`, else a folder `PaleoMapR` in your home folder (on
Windows: Documents).

## Usage

``` r
pm_folder()
```

## Value

The path of the working folder (it is not created).

## Examples

``` r
pm_folder()
#> [1] "/home/runner/PaleoMapR"
```
