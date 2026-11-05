# Contributing to PaleoMapR

PaleoMapR is in alpha, so the most useful contribution right now is
feedback from making real maps with it.

## Reporting a problem

Open a [bug
report](https://github.com/stratigraphy/PaleoMapR/issues/new?template=bug_report.yml)
and include:

- what you did (tab, map age, plate model, button or settings) and what
  you expected;
- the message in the status line or the R console; for a failed figure,
  `maps/last_render_error.txt`; for the Windows setup,
  `state/setup_log.txt`;
- if the problem depends on your points, `points/user_points.csv` and
  `points/record_edits.csv`.

Running [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) in R
and pasting the output helps for installation problems.

## Suggesting a change

Open a [feature
request](https://github.com/stratigraphy/PaleoMapR/issues/new?template=feature_request.yml),
ideally with a reference to the data set or method, or an example map.

## Code

Pull requests are welcome; please open an issue first for anything
larger than a small fix. The editor code is in `inst/app/`: `app.R` (the
Shiny app) and the scripts in `inst/app/R/` (settings in `config.R`, the
correction in `dem_correction.R`, maps in `maps.R`, climate in
`climate*.R`); `R/` holds the package functions that start and set up
the editor. Code style for new code: base R, braces on every `if`/`for`
body, no apply family or tidyverse,
[`par()`](https://rdrr.io/r/graphics/par.html) and
[`options()`](https://rdrr.io/r/base/options.html) restored with
[`on.exit()`](https://rdrr.io/r/base/on.exit.html).
