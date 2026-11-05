# Set up the editor (once per working folder)

Prepares a working folder in four steps. It needs internet and takes 20
to 60 minutes the first time, mostly downloads (about 1 GB). Running it
again is safe: finished steps and files are kept, so after an
interrupted download simply run it again.

1.  Python with pyGPlates (plate rotations): a small Python environment
    is made with 'uv', which downloads Python 3.12 itself; no Python or
    conda installation is needed. It is kept in your R user data folder
    and shared by all working folders. If that fails, a conda
    environment "paleomap" is tried.

2.  The open data sets into `data/`: PaleoDEMs 0-540 Ma (Scotese and
    Wright 2018), PaleoCoastlines (Kocsis and Scotese 2021), PBDB
    collections, Macrostrat units and the sea-level curve of van der
    Meer et al. (2022); then the table of database records.

3.  The HadCM3L climate files (Valdes et al. 2021) of the chosen ages,
    about 30 MB per age.

4.  The fossil-corrected grids of the chosen ages (1-2 minutes per age
    per core).

## Usage

``` r
pm_setup(folder = pm_folder(), ages = "config", workers = NULL)
```

## Arguments

- folder:

  The working folder (see
  [`pm_folder()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_folder.md)).

- ages:

  The map ages to prepare: `"config"` (the Silurian ages 445-415 Ma),
  `"all"` (all 109 PaleoDEM ages, 0-540 Ma: several hours and about 3 GB
  of climate files), or ages in Ma such as `c(300, 305, 310)`, each
  matched to the nearest PaleoDEM age. Other ages work in the editor
  too: choose them and press "Rebuild this age".

- workers:

  The number of R processes for the corrected grids; by default up to 8,
  one less than the cores.

## Value

`TRUE` (invisibly) when everything is in place, `FALSE` when a part is
missing.

## Examples

``` r
if (FALSE) { # \dontrun{
pm_setup()                                   # the Silurian ages, in pm_folder()
pm_setup("D:/Maps/PaleoMapR", ages = seq(250, 200, by = -5))
} # }
```
