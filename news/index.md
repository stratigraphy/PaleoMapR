# Changelog

## PaleoMapR 0.1.0 (alpha)

First public release, for testing and feedback. Buttons, defaults and
file names may still change.

- The paleogeographic map editor (standalone version of 2026-10-26) as
  an R package:
  [`run_app()`](https://stratigraphy.github.io/PaleoMapR/reference/run_app.md)
  starts it,
  [`pm_setup()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_setup.md)
  prepares a working folder (Python with pyGPlates, the open data sets,
  HadCM3L climate files and the corrected grids of the chosen ages),
  [`pm_load()`](https://stratigraphy.github.io/PaleoMapR/reference/pm_load.md)
  gives its functions for scripts.
- Windows launcher (`PaleoMapR-launcher.zip` on the release page):
  SETUP.bat installs everything, prepares the data in the launcher
  folder and puts a desktop icon.
- Map editor: every PaleoDEM age 0–540 Ma, four plate models, PBDB and
  Macrostrat records, your own points (marine, deeper shelf, evaporite,
  terrestrial, reef, glacial), excluded records, undo, rebuilding of the
  fossil-corrected grid with statistics.
- Figure tab: many projections, whole map / zoomed-in / composite,
  localities, symbols, scale bar and north arrow, HadCM3L climate layers
  (as run, or with a local correction to your map), wind and current
  arrows, rivers, sea ice.
- Faster map switching and Figure redraws (cached editor maps and relief
  shading).
- About tab with the version, the works to cite, a daily update check
  with “Update and restart”, and a button that closes the editor.
- The table of database records is rebuilt when its classification rules
  change (checksum), not after every update of the package.
