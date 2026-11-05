PaleoMapR - fossil-corrected paleogeographic maps (Windows launcher)
=====================================================================

Install
  1. Install R from https://cloud.r-project.org/bin/windows/base/ (default options are fine).
  2. Unzip this somewhere permanent with a short, simple path and about 2 GB free, e.g. into Documents\Apps
     (not in Downloads, and preferably not in a OneDrive or Dropbox folder); this gives the folder PaleoMapR.
     This folder is also your working folder: the data, your points and your maps are kept in it.
  3. Double-click SETUP.bat. It installs PaleoMapR and the packages it needs, puts the shortcut "PaleoMapR" on
     the desktop, and then prepares the data: Python with pyGPlates, the open data sets (about 1 GB), the
     HadCM3L climate files and the corrected grids of the Silurian ages (445-415 Ma). This takes 20-60 minutes
     the first time. If a download fails, run SETUP.bat again: finished parts are kept.
     Other ages: run SETUP.bat from a command window with the ages, e.g.  SETUP.bat 300 305 310  or
     SETUP.bat all  (all 109 ages, 0-540 Ma: several hours). Any age also works in the editor itself:
     choose it and press "Rebuild this age".

Start
  Double-click "PaleoMapR" on the desktop. The editor opens in your browser at http://127.0.0.1:7793; it runs
  on your own computer. Closing the browser tab (or About > Close the editor) stops R.

What is in this folder
  data\      downloaded data sets, HadCM3L climate files and the corrected grids
  models\    plate models (downloaded when first used)
  points\    your work: user_points.csv, record_edits.csv, localities.csv, settings.csv (plain CSV files)
  maps\      rendered maps (PNG and PDF)
  results\   statistics of the corrected grids
  state\     logs of the setup and of updates

Updates
  Once a day PaleoMapR checks GitHub for a newer version and offers "Update and restart". About > Check for
  updates checks at once. Running SETUP.bat again also updates. Your data and points are kept.

Problems, questions, ideas
  https://github.com/stratigraphy/PaleoMapR/issues   (state\setup_log.txt and state\update_log.txt help)
  Manual and tutorial: https://stratigraphy.github.io/PaleoMapR/

Remove
  Delete the desktop shortcut and this folder; in R: remove.packages("PaleoMapR"). The Python environment is in
  %APPDATA%\R\data\R\paleomapeditor\venv.
