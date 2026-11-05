# Start the paleogeographic map editor

Opens the editor in your web browser. It runs on your own computer; your
points, edits and maps are plain files in the working folder. Run
[`pm_setup()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_setup.md)
once first: without it, the first start downloads the basic data (about
400 MB) before the window opens, and climate layers and corrected grids
are made when you first need them.

## Usage

``` r
run_app(
  folder = pm_folder(),
  port = NULL,
  launch.browser = TRUE,
  stop_on_close = FALSE
)
```

## Arguments

- folder:

  The working folder (see
  [`pm_folder()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_folder.md)).

- port:

  The port of the local web server; `NULL` picks a free one.

- launch.browser:

  Open the editor in the web browser.

- stop_on_close:

  Stop R when the browser tab is closed (used by the desktop launcher).

## Value

Nothing; returns when the editor is stopped.

## Examples

``` r
if (FALSE) { # \dontrun{
run_app()                        # working folder pm_folder()
run_app("D:/Maps/PaleoMapR")
} # }
```
