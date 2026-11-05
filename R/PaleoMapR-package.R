#' PaleoMapR: fossil-corrected paleogeographic maps in a point-and-click editor
#'
#' Start the editor with \code{\link{run_app}()}. The first time, \code{\link{pm_setup}()} downloads the open data
#' sets (PaleoDEMs, PaleoCoastlines, PBDB, Macrostrat, sea level, HadCM3L climate) into a working folder and builds
#' the corrected grids of the Silurian ages. \code{\link{pm_load}()} gives the functions behind the editor for use in
#' an R script.
#'
#' @keywords internal
#' @importFrom DT datatable
#' @importFrom RANN nn2
#' @importFrom jsonlite fromJSON
#' @importFrom ncdf4 nc_open
#' @importFrom png readPNG
#' @importFrom readxl read_excel
#' @importFrom reticulate py_available
#' @importFrom sf st_read
#' @importFrom terra rast
#' @importFrom callr r_bg
"_PACKAGE"
