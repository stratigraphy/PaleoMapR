# PaleoMapR 0.1.1

Bug fixes and speed-ups after a full review of the editor.

## Editor (app)

* An error in a button or click (for example a CSV file in `points/` that is open in Excel) is now shown in the
  status line; before, it ended the session and, with the launcher, stopped R without a message.
* Undo restores `user_points.csv` and `record_edits.csv` exactly as they were before the last change made in the
  editor (up to 50 steps). Before, it removed the last row of a table, which could delete the wrong point.
* "Delete selected points" and "Remove selected edits" remove the rows you selected by their content, also when the
  files changed meanwhile.
* After a drag or zoom, clicks and wheel steps are ignored until the new map is shown; a quick click could land in
  the wrong place. A sideways trackpad swipe no longer zooms in.
* The Figure map is updated after adding points, edits or rebuilding a grid.
* Ticking symbols, north arrow or scale bar no longer recomputes the whole Figure map.
* A locality with a latitude outside -90 to 90 is refused (it used to break the editor map).
* A point whose age range does not include the map age gets a warning when it is added.
* The download buttons appear once a figure has been rendered.
* A zoomed-in or composite figure of an area wider than 180 degrees of longitude is refused with a message.
* Rebuilding reports ages that failed instead of stopping; the other ages are kept.

## Maps

* Winkel tripel maps work (they always failed) and use the proper standard parallel (50.467 degrees).
* Relief shading is lit from the north-west, as intended (it was lit from the south-west).
* "Four seasons on one page" works (it stopped with "could not find function .season_climates").
* Globe and polar overviews of a composite turn towards the zoomed-in area, so its box is visible.
* Symbols and points on the far side of a globe or outside a polar map are no longer drawn.
* Zoomed-in maps of open sea (no land) no longer fail; small windows keep their graticule labels.
* Graphics devices are always closed, also after an error; a PDF or PNG that is open in a viewer gives a clear
  message. Without cairo, the standard PNG and PDF devices are used.
* Locality names with / : * ? and similar characters work as map centre (file names are made safe).
* Single-age renders get their own PDF name (rendering 430 Ma no longer overwrites the 425 Ma PDF); polar maps
  with a different edge latitude and zoom windows less than 1 degree apart get different file names.
* Rendering is about twice as fast: each page's image is computed once for the PNG and re-used for the PDF.
* The Figure preview's scale bar refers to the map centre.

## Climate layers

* Climate files are stored per run and model age (`data/climate/<run>_<age>Ma/`). Run IDs that differ only in
  upper / lower case (e.g. teXPh 425 Ma and teXph 297 Ma) shared one folder on Windows and macOS, so another age's
  climate could be shown. Existing files in the old `data/climate/<run>/` folders are not used; the new folders are
  filled when needed (about 30 MB per age).
* Downloaded files are checked before they are kept; an interrupted download or a proxy error page is no longer used
  as data.
* Local correction: runoff on new land no longer counts the precipitation change twice.
* Temperature: both climate options use the model's own terrain as reference and your lapse rate; lapse rate 0 now
  switches the correction off completely, and the two options agree where the map is unchanged.
* "Why the sea is fresher / saltier" uses the class of the nearest model cell (no blended, invented classes).
* No blank strips along the coasts in the river, freshwater, salinity-source and runoff ("as run") layers.
* Mixed-layer depth and sea ice extend into new seaways from the nearby sea (they were 0 there).
* With "Original PaleoDEM", sea ice is linked to the land ice of that map, not of the corrected one.

## Corrected grids

* Seaway widths and banks are measured with true distances on the sphere. The old east-west distance depended on
  the grid column, so north-south seaways far from 180 degrees came out too wide. Rebuild your ages to apply this.
* The flooded patch around a marine record is round also at high latitude (it was cut off east-west above ~45 deg).
* Records on crust that the plate model adds only after the map age (e.g. later accreted terranes) are left out;
  their number is reported.
* PBDB and Macrostrat records of the whole 400-450 Ma window are downloaded (the start-up Silurian file missed the
  latest Ordovician and Early Devonian, e.g. 415 Ma had no PBDB records).
* Points tables written by Excel are read correctly: byte-order mark, semicolons and decimal commas, TRUE/FALSE in
  other languages (WAAR, VRAI, ...), swapped ages, repeated IDs. A table that lacks its columns is never written back.
* User points with impossible coordinates are skipped with a warning instead of stopping every rebuild.
* `records_base.csv` is rebuilt when the classification settings in `config.R` change.
* Rebuilding no longer resets R's random numbers, and a failing age no longer discards the others' statistics;
  `step6_records_used.csv` and `step6_seaways.csv` keep the other ages.
* Downloads that failed are not retried for 10 minutes, so an offline computer does not wait at every rebuild.

## Launcher and setup

* The launcher works in folders with accents or apostrophes in their path (UTF-16 VBS, UTF-8 PowerShell script with
  byte-order mark, quoted paths), and the "already running" check works in every Windows language.
* After an R upgrade, the desktop icon says that SETUP.bat has to be run again instead of failing silently.
* The editor started by the icon writes its messages to `state/app_log.txt`.
* SETUP.bat reports arguments it does not understand.
* The update check reads today's result directly instead of starting an R process at every page load.
* The editor scripts are no longer loaded twice at start.

# PaleoMapR 0.1.0 (alpha)

First public release, for testing and feedback. Buttons, defaults and file names may still change.

* The paleogeographic map editor (standalone version of 2026-10-26) as an R package: `run_app()` starts it,
  `pm_setup()` prepares a working folder (Python with pyGPlates, the open data sets, HadCM3L climate files and the
  corrected grids of the chosen ages), `pm_load()` gives its functions for scripts.
* Windows launcher (`PaleoMapR-launcher.zip` on the release page): SETUP.bat installs everything, prepares the
  data in the launcher folder and puts a desktop icon.
* Map editor: every PaleoDEM age 0–540 Ma, four plate models, PBDB and Macrostrat records, your own points
  (marine, deeper shelf, evaporite, terrestrial, reef, glacial), excluded records, undo, rebuilding of the
  fossil-corrected grid with statistics.
* Figure tab: many projections, whole map / zoomed-in / composite, localities, symbols, scale bar and north arrow,
  HadCM3L climate layers (as run, or with a local correction to your map), wind and current arrows, rivers, sea ice.
* Faster map switching and Figure redraws (cached editor maps and relief shading).
* About tab with the version, the works to cite, a daily update check with "Update and restart", and a button
  that closes the editor.
* The table of database records is rebuilt when its classification rules change (checksum), not after every
  update of the package.
