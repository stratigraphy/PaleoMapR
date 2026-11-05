# Check GitHub for a newer version of PaleoMapR

Compares the installed commit with the newest commit on GitHub. The
editor does this once a day by itself and offers "Update and restart".

## Usage

``` r
pm_update_status(force = FALSE)
```

## Arguments

- force:

  Check now instead of using the result of the last check (kept for a
  day).

## Value

A one-row data frame: package, repository, origin of the installation,
version, installed and newest commit, and whether an update is
available.

## Examples

``` r
if (FALSE) { # \dontrun{
pm_update_status(force = TRUE)
} # }
```
