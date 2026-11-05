# Start the editor ####

#' Start the paleogeographic map editor
#'
#' Opens the editor in your web browser. It runs on your own computer; your points, edits and maps are plain files
#' in the working folder. Run \code{\link{pm_setup}()} once first: without it, the first start downloads the basic
#' data (about 400 MB) before the window opens, and climate layers and corrected grids are made when you first need
#' them.
#'
#' @param folder The working folder (see \code{\link{pm_folder}()}).
#' @param port The port of the local web server; \code{NULL} picks a free one.
#' @param launch.browser Open the editor in the web browser.
#' @param stop_on_close Stop R when the browser tab is closed (used by the desktop launcher).
#' @return Nothing; returns when the editor is stopped.
#' @examples
#' \dontrun{
#' run_app()                        # working folder pm_folder()
#' run_app("D:/Maps/PaleoMapR")
#' }
#' @export
run_app <- function(folder = pm_folder(), port = NULL, launch.browser = TRUE, stop_on_close = FALSE) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)
  oldop <- options(PaleoMapR.folder = folder, PaleoMapR.stop_on_close = isTRUE(stop_on_close))
  on.exit(options(oldop), add = TRUE)
  message("PaleoMapR ", as.character(utils::packageVersion("PaleoMapR")), ", working folder ", folder)
  args <- list(appDir = pm_code_dir(), launch.browser = launch.browser)
  if (!is.null(port)) {
    args$port <- as.integer(port)
  }
  do.call(shiny::runApp, args)
  invisible(NULL)
}
