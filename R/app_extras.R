# About tab, update check and stopping, used by inst/app/app.R ####

.pm_state <- new.env(parent = emptyenv())
.pm_state$sessions <- 0L

pm_links <- function() {
  list(site = "https://stratigraphy.github.io/PaleoMapR/",
      repo = "https://github.com/stratigraphy/PaleoMapR",
      issues = "https://github.com/stratigraphy/PaleoMapR/issues/new/choose",
      discussions = "https://github.com/stratigraphy/PaleoMapR/discussions")
}

# The works to cite: label, reference, link
pm_data_citations <- function() {
  data.frame(
      what = c("Climate (HadCM3L)", "Climate model", "CO2 of the climate runs", "PaleoDEMs and PALEOMAP plate model",
          "PaleoCoastlines", "Fossil records", "Rock units", "Sea level", "Plate rotations", "Other plate models (if used)",
          "Present-day borders and coastlines"),
      ref = c("Valdes, P. J., Scotese, C. R. and Lunt, D. J. (2021). Deep ocean temperatures through time. Climate of the Past 17, 1483-1506.",
          "Valdes, P. J. et al. (2017). The BRIDGE HadCM3 family of climate models: HadCM3@Bristol v1.0. Geoscientific Model Development 10, 3715-3743.",
          "Foster, G. L., Royer, D. L. and Lunt, D. J. (2017). Future climate forcing potentially without precedent in the last 420 million years. Nature Communications 8, 14845.",
          "Scotese, C. R. and Wright, N. (2018). PALEOMAP Paleodigital Elevation Models (PaleoDEMS) for the Phanerozoic. PALEOMAP Project.",
          "Kocsis, \u00c1. T. and Scotese, C. R. (2021). Mapping paleocoastlines and continental flooding during the Phanerozoic. Earth-Science Reviews 213, 103463.",
          "The Paleobiology Database (data service API 1.2).",
          "Peters, S. E., Husson, J. M. and Czaplewski, J. (2018). Macrostrat: a platform for geological data integration and deep-time Earth crust research. Geochemistry, Geophysics, Geosystems 19, 1393-1409.",
          "van der Meer, D. G. et al. (2022). Long-term Phanerozoic global mean sea level: insights from strontium isotope variations and estimates of continental glaciation. Gondwana Research 111, 103-121.",
          "M\u00fcller, R. D. et al. (2018). GPlates: building a virtual Earth through deep time. Geochemistry, Geophysics, Geosystems 19, 2243-2261.",
          "Merdith, A. S. et al. (2021). Extending full-plate tectonic models into deep time. Earth-Science Reviews 214, 103477. Torsvik, T. H. and Cocks, L. R. M. (2017). Earth History and Palaeogeography. Cambridge University Press.",
          "Natural Earth (public domain)."),
      link = c("https://doi.org/10.5194/cp-17-1483-2021", "https://doi.org/10.5194/gmd-10-3715-2017",
          "https://doi.org/10.1038/ncomms14845", "https://www.earthbyte.org/paleodem-resource-scotese-and-wright-2018/",
          "https://doi.org/10.1016/j.earscirev.2020.103463", "https://paleobiodb.org",
          "https://doi.org/10.1029/2018GC007467", "https://doi.org/10.1016/j.gr.2022.07.014",
          "https://doi.org/10.1029/2018GC007584", "https://doi.org/10.1016/j.earscirev.2020.103477",
          "https://www.naturalearthdata.com"),
      stringsAsFactors = FALSE)
}

pm_version <- function() {
  as.character(utils::packageVersion("PaleoMapR"))
}

# The About tab
pm_about_tab <- function() {
  a <- function(href, txt) {
    shiny::tags$a(href = href, target = "_blank", txt, .noWS = "after")
  }
  L <- pm_links()
  cit <- pm_data_citations()
  items <- vector("list", nrow(cit))
  for (i in seq_len(nrow(cit))) {
    items[[i]] <- shiny::tags$li(shiny::tags$b(cit$what[i], .noWS = "after"), ": ", cit$ref[i], " ", a(cit$link[i], cit$link[i]))
  }
  shiny::tabPanel("About", icon = shiny::icon("circle-info"), value = "about",
      shiny::tags$head(shiny::tags$script(shiny::HTML(paste0(
          "Shiny.addCustomMessageHandler('pm_closed', function(m) {",
          "  document.body.innerHTML = '<div style=\"font-family:sans-serif;max-width:560px;margin:80px auto;font-size:15px\">' +",
          "    '<h3>' + m.title + '</h3><p>' + m.text + '</p></div>';",
          "});")))),
      shiny::fluidRow(shiny::column(8, offset = 1,
          shiny::h3(sprintf("PaleoMapR %s", pm_version())),
          shiny::p("Paleogeographic map editor by Michiel Arts: checks PaleoDEMs against the fossil and rock record,",
              "corrects them, and draws publication maps with climate layers. Everything you do is stored as plain",
              "files in the working folder:"),
          shiny::tags$pre(getOption("PaleoMapR.folder", getwd())),
          shiny::p("Manual and tutorial: ", a(L$site, L$site), shiny::br(),
              "Source code: ", a(L$repo, L$repo), shiny::br(),
              "Something wrong or missing? ", a(L$issues, "Report it"), " or ", a(L$discussions, "ask in Discussions"), "."),
          shiny::div(style = "margin:12px 0",
              shiny::actionButton("pm_check_update", "Check for updates", icon = shiny::icon("rotate")),
              shiny::actionButton("pm_quit", "Close the editor (stops R)", icon = shiny::icon("power-off"))),
          shiny::h4("How to cite"),
          shiny::p("Please cite PaleoMapR (see ", shiny::tags$code("citation(\"PaleoMapR\")"),
              ") and the data behind your maps:"),
          shiny::tags$ul(items),
          shiny::p("The 'local correction' of the climate to your map is this editor's own statistical adjustment,",
              "not HadCM3L output; say which climate option you used in your figure captions."))))
}

# Update check, "Update and restart", "Close" and stopping R when the last browser tab is closed
pm_server_extras <- function(input, session) {
  if (isTRUE(getOption("PaleoMapR.stop_on_close"))) {
    .pm_state$sessions <- .pm_state$sessions + 1L
    session$onSessionEnded(function() {
      .pm_state$sessions <- .pm_state$sessions - 1L
      # wait a few seconds: a page reload ends the session and starts a new one
      later::later(function() {
        if (.pm_state$sessions <= 0L) {
          shiny::stopApp()
        }
      }, 6)
    })
  }

  upd <- shiny::reactiveValues(job = NULL, status = NULL, force = FALSE)
  start_check <- function(force = FALSE) {
    if (!is.null(shiny::isolate(upd$job))) {
      return(invisible(NULL))
    }
    upd$force <- force
    job <- try(callr::r_bg(function(force) { PaleoMapR::pm_update_status(force) }, args = list(force = force)),
        silent = TRUE)
    upd$job <- if (inherits(job, "try-error")) { NULL } else { job }
  }
  if (isTRUE(getOption("PaleoMapR.check_updates", TRUE))) {
    start_check()
  }
  shiny::observe({
    j <- upd$job
    if (is.null(j)) {
      return()
    }
    if (j$is_alive()) {
      shiny::invalidateLater(1000)
      return()
    }
    st <- try(j$get_result(), silent = TRUE)
    upd$job <- NULL
    force <- shiny::isolate(upd$force)
    if (inherits(st, "try-error")) {
      if (force) {
        shiny::showNotification("The update check failed.", id = "pm_upd", type = "error")
      }
      return()
    }
    upd$status <- st
    pm_update_note(st, force)
  })
  shiny::observeEvent(input$pm_check_update, {
    shiny::showNotification("Checking GitHub for updates ...", id = "pm_upd", duration = 4)
    start_check(force = TRUE)
  })
  shiny::observeEvent(input$pm_upd_go, {
    st <- upd$status
    if (is.null(st) || !any(st$update)) {
      return()
    }
    if (pm_start_updater(st)) {
      session$sendCustomMessage("pm_closed", list(title = "PaleoMapR is being updated",
          text = "A window shows the progress; the editor opens again in a new tab when the update is installed."))
      shiny::stopApp()
    } else {
      shiny::showModal(shiny::modalDialog(title = "Update PaleoMapR", easyClose = TRUE,
          shiny::tags$p("Close the editor, then run this in a fresh R session and start the editor again:"),
          shiny::tags$pre(paste(pm_update_code(st), collapse = "\n"))))
    }
  })
  shiny::observeEvent(input$pm_quit, {
    session$sendCustomMessage("pm_closed", list(title = "The editor is closed",
        text = "R has stopped. Your points, edits and maps are saved in the working folder. Start the editor again with its desktop icon or run_app()."))
    shiny::stopApp()
  })
  invisible(NULL)
}

pm_update_note <- function(st, force = FALSE) {
  if (any(st$update)) {
    shiny::showNotification(id = "pm_upd", duration = NULL, type = "warning",
        ui = shiny::tags$span("A newer version of PaleoMapR is available on GitHub."),
        action = shiny::actionLink("pm_upd_go", "Update and restart"))
    return(invisible(TRUE))
  }
  if (force) {
    o <- st$origin[1]
    txt <- if (identical(o, "local")) {
      "development install (not checked)"
    } else if (is.na(st$latest[1])) {
      "GitHub could not be reached"
    } else {
      "up to date"
    }
    shiny::showNotification(sprintf("PaleoMapR %s: %s", st$version[1], txt), id = "pm_upd", type = "message", duration = 8)
  }
  invisible(FALSE)
}
