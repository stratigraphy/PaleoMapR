# Paleogeographic map editor (Shiny): Map editor, Figure and About tabs ####
#
# Started by PaleoMapR::run_app(folder), which sets the option PaleoMapR.folder (the working folder with data/,
# points/, maps/ ...). The editor scripts in R/ next to this file are loaded into this app's environment.
#
# - Click on the palaeo-map to add a point (marine, evaporite, terrestrial, reef, glacial) at the selected map
#   age, or to exclude / re-include the nearest database record or your own point.
# - Hold and drag to move the map; mouse wheel or the +/- buttons to zoom.
# - Everything is written to points/user_points.csv and points/record_edits.csv, so you can also edit those
#   files by hand.
# - 'Rebuild this age' re-runs the DEM correction (about 1 minute per age); 'Render map' draws the map (PNG +
#   PDF) for this age.
library(shiny)
library(DT)
PROJECT_DIR <- normalizePath(getOption("PaleoMapR.folder", getwd()), winslash = "/", mustWork = FALSE)
PaleoMapR::pm_load(PROJECT_DIR, envir = environment())
# shiny restores the working directory when the app stops
setwd(PROJECT_DIR)
init_user_tables()
invisible(read_localities())
# first start without pm_setup(): fetch the data here (finished downloads are kept)
if (!dir.exists(file.path(SILUCARB$dir_coast, "CS"))) download_data()
download_dems()                                   # ~280 MB zip, once
if (!file.exists(points_file("records_base.csv"))) build_base_records()

# every PaleoDEM age (Phanerozoic, 0-540 Ma) with its stage name, in chronological order (oldest first)
AGE_CHOICES <- local({
  a <- sort(unique(c(available_ages(SILUCARB), SILUCARB$ages)), decreasing = TRUE)
  setNames(as.character(a), sprintf("%s Ma \u2013 %s", vapply(a, function(v) format(v), ""),
                                    vapply(a, stage_label, "", cfg = SILUCARB)))
})

SETTING_STYLE <- data.frame(setting = c("marine", "marine_deep", "evaporite", "terrestrial", ""),
                            pch = c(21, 21, 22, 24, 23), bg = c("#2f7fd8", "#0b2a6e", "#eda100", "#a0522d", "#ffffff"),
                            stringsAsFactors = FALSE)
MODES <- c("Add marine point (shallow shelf)" = "add_marine",
           "Add marine point (deeper shelf)" = "add_marine_deep", "Add evaporite point" = "add_evaporite",
           "Add terrestrial point" = "add_terrestrial", "Add reef (marine + reef symbol)" = "add_reef",
           "Add glacial deposit" = "add_glacial",
           "Exclude / re-include nearest record" = "toggle", "Inspect nearest record" = "inspect",
           "Add locality here (name = note)" = "add_locality",
           "Delete nearest of YOUR points (permanent)" = "delete_point")


MAP_JS <- "
// drag to pan, mouse wheel to zoom, click: for the editor map (inputs map_*) and the Figure map (fig_map_*)
function bindNav(id) {
  var el = document.getElementById(id); if (!el) return;
  var down = null, moved = false;
  el.addEventListener('mousedown', function(e) {
    if (e.button !== 0) return;
    var img = el.querySelector('img'); if (!img) return;
    down = {x: e.clientX, y: e.clientY, img: img, r: img.getBoundingClientRect()}; moved = false;
    e.preventDefault();
  });
  window.addEventListener('mousemove', function(e) {
    if (!down) return;
    var dx = e.clientX - down.x, dy = e.clientY - down.y;
    if (Math.abs(dx) + Math.abs(dy) > 5) moved = true;
    if (moved) { down.img.style.transform = 'translate(' + dx + 'px,' + dy + 'px)'; el.style.cursor = 'grabbing'; }
  });
  window.addEventListener('mouseup', function(e) {
    if (!down) return;
    var r = down.r, dx = e.clientX - down.x, dy = e.clientY - down.y;
    if (moved) Shiny.setInputValue(id + '_pan', {dx: dx / r.width, dy: dy / r.height, n: Math.random()});
    else Shiny.setInputValue(id + '_click2', {fx: (e.clientX - r.left) / r.width, fy: (e.clientY - r.top) / r.height, n: Math.random()});
    down = null; el.style.cursor = '';
  });
  var acc = 0, last = null, timer = null;
  el.addEventListener('wheel', function(e) {
    var img = el.querySelector('img'); if (!img) return;
    e.preventDefault();
    var r = img.getBoundingClientRect();
    acc += e.deltaY; last = {fx: (e.clientX - r.left) / r.width, fy: (e.clientY - r.top) / r.height};
    clearTimeout(timer);
    timer = setTimeout(function() {
      Shiny.setInputValue(id + '_wheel', {steps: acc > 0 ? -1 : 1, fx: last.fx, fy: last.fy, n: Math.random()});
      acc = 0;
    }, 180);
  }, {passive: false});
}
$(document).on('shiny:connected', function() { bindNav('map'); bindNav('fig_map'); });
"
ui <- navbarPage(paste("PaleoMapR", APP_VERSION), id = "main", collapsible = TRUE,
  windowTitle = "PaleoMapR",
  header = tagList(tags$script(HTML(MAP_JS)),
    tags$head(tags$style(HTML("body{font-size:13px} .well{padding:10px} #map{cursor:crosshair; overflow:hidden} #fig_map{cursor:grab; overflow:hidden; border:1px solid #ddd}
                               #fig_preview img{max-width:100%; height:auto; border:1px solid #ddd} pre{white-space:pre-wrap}
                               summary{display:list-item; cursor:pointer; color:#2a6496; font-weight:600; margin:6px 0} summary:hover{text-decoration:underline}
                               details{border-left:3px solid #ccc; padding-left:8px; margin:6px 0}")))),
  # ---------------------------------------------------------------- map editor --
  tabPanel("Map editor", icon = icon("pen"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("age", "Map age (Ma)", choices = AGE_CHOICES, selected = 425),
      selectInput("model", "Plate model of the map", choices = setNames(names(SILUCARB$models), SILUCARB$models),
                  selected = "scotese_and_wright2018"),
      radioButtons("grid", "Background grid", c("Corrected" = "corrected", "Original PaleoDEM" = "raw"), inline = TRUE),
      radioButtons("mode", "Click action", MODES, selected = "inspect"),
      fluidRow(column(6, numericInput("age_young", "Point age: young (Ma)", 425, step = 1)),
               column(6, numericInput("age_old", "old (Ma)", 425, step = 1))),
      helpText("A point is used at every map age within its age range \u00b1",
               SILUCARB$age_window, "Myr."),
      textInput("note", "Note for new points / edits", ""),
      checkboxGroupInput("show", "Show", c("marine", "deeper marine" = "marine_deep", "evaporite", "terrestrial", "indicators", "excluded",
                                           "country borders" = "borders", "coastlines"),
                         selected = c("marine", "marine_deep", "evaporite", "terrestrial", "indicators", "excluded", "borders", "coastlines"),
                         inline = TRUE),
      helpText("Present-day borders (black) and coastlines (white) are rotated to the map age with the selected",
               "plate model; dotted = crust added later or approximate position."),
      hr(),
      actionButton("undo", "Undo last change", icon = icon("rotate-left")),
      br(), br(),
      actionButton("rebuild", "Rebuild this age", class = "btn-primary", icon = icon("hammer")),
      actionButton("render", "Render map", icon = icon("map")),
      helpText("'Render map' uses the settings of the Figure tab."),
      br(), br(),
      actionButton("rebuild_all", "Rebuild the Silurian ages (445\u2013415 Ma)", icon = icon("layer-group")),
      helpText("Rebuilding writes the corrected grid to", SILUCARB$dir_corrected,
               "and the statistics to", SILUCARB$dir_results),
      hr(),
      h5("Land points (used when rebuilding)"),
      fluidRow(column(6, numericInput("land_radius", "Land radius (km)", SILUCARB$LAND_POINT_RADIUS_KM, 5, 1000, 5)),
               column(6, numericInput("land_link", "Join within (km)", SILUCARB$LAND_LINK_KM, 0, 3000, 25))),
      helpText("Land kept around each land point, and the distance within which land points are joined",
               "into one landmass. Saved in points/settings.csv.")
    ),
    mainPanel(width = 9,
      div(style = "margin-bottom:4px",
          actionButton("zoom_in", NULL, icon = icon("magnifying-glass-plus"), title = "Zoom in"),
          actionButton("zoom_out", NULL, icon = icon("magnifying-glass-minus"), title = "Zoom out"),
          actionButton("zoom_reset", NULL, icon = icon("earth-europe"), title = "Whole map"),
          actionButton("pan_w", NULL, icon = icon("arrow-left")), actionButton("pan_e", NULL, icon = icon("arrow-right")),
          actionButton("pan_n", NULL, icon = icon("arrow-up")), actionButton("pan_s", NULL, icon = icon("arrow-down")),
          span(style = "color:#666; margin-left:8px",
               "Hold and drag to move the map \u00b7 mouse wheel or +/\u2212 to zoom \u00b7 click for the chosen action")),
      plotOutput("map", height = "600px"),
      verbatimTextOutput("status"),
      tabsetPanel(
        tabPanel("Nearest records", DTOutput("near")),
        tabPanel("Your points", DTOutput("upoints"), actionButton("del_points", "Delete selected points")),
        tabPanel("Your edits", DTOutput("uedits"), actionButton("del_edits", "Remove selected edits")),
        tabPanel("Localities", DTOutput("loctab"),
                 fluidRow(column(4, textInput("loc_name", "Name")),
                          column(3, numericInput("loc_lat", "Present-day lat", NA, -90, 90)),
                          column(3, numericInput("loc_lon", "Present-day lon", NA, -180, 180)),
                          column(2, br(), actionButton("loc_add", "Add"))),
                 actionButton("loc_del", "Remove selected localities")),
        tabPanel("Statistics", tableOutput("stats"))
      )
    )
  )),
  tabPanel("Figure", icon = icon("image"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Map and climate"),
      selectInput("fig_layer", "Map colours", c("topography" = "topography", "precipitation" = "precipitation",
                                                 "temperature (near-surface air)" = "temperature",
                                                 "precipitation \u2212 evaporation" = "pme",
                                                 "runoff (total, land)" = "runoff",
                                                 "freshwater into the sea (P \u2212 E + rivers)" = "freshwater",
                                                 "river water into the sea" = "river_input",
                                                 "why the sea is fresher / saltier (rivers, rain, evaporation, currents)" = "salinity_source",
                                                 "sea-surface salinity (5 m)" = "salinity",
                                                 "ocean up- / downwelling (10 m)" = "w_ocean",
                                                 "ocean mixed-layer depth" = "mld",
                                                 "surface O2 (saturation)" = "o2_surface",
                                                 "sea-surface density (where deep water can form)" = "surface_density",
                                                 "sea ice: permanent and seasonal" = "sea_ice")),
      radioButtons("fig_kind", "Figure layout",
                   c("Whole map" = "whole", "Zoomed-in area only" = "zoom",
                     "Composite: whole map + zoomed-in area (joined by lines)" = "composite"), selected = "whole"),
      helpText("For a zoomed-in or composite figure, move and zoom the Figure map (right) to the area you want, as in the",
               "Map editor (or choose another way under 'Area to render' below). The Figure map shows the chosen layer and",
               "everything ticked under 'What is drawn', as in the rendered figure."),
      selectInput("season", "Season / month (climate)", setNames(names(CLIMATE_SEASONS), CLIMATE_SEASONS)),
      checkboxInput("fig_seasons4", "Four seasons (DJF, MAM, JJA, SON) on one page, one colour scale", FALSE),
      conditionalPanel("input.fig_layer == 'temperature'",
        numericInput("lapse_rate", "Lapse rate for the map heights (K per km; 0 = off)", 6.5, 0, 10, 0.5),
        helpText("The model temperatures belong to its smooth ~3.75\u00b0 terrain; they are brought to the 0.2\u00b0 map heights",
                 "(and the ice surface of an ice run). 6.5 K/km: free atmosphere; ~5 K/km near the surface over land; up to 8 K/km over ice.")),
      conditionalPanel("input.fig_layer == 'o2_surface'",
        numericInput("o2_pal", "Atmospheric O2 (PAL; 1 = 20.95 %)", 1, 0.1, 2, 0.05),
        helpText("Silurian atmospheric O2 is uncertain: estimates range from well below to about present levels (e.g. Lenton et al. 2018; Krause et al. 2018).")),
      conditionalPanel("input.fig_layer == 'freshwater' || input.fig_layer == 'river_input' || input.fig_layer == 'salinity_source'",
        helpText("Runoff is routed downhill over the map (1\u00b0 grid, depressions filled) to the coast; the river water is spread",
                 "~300 km around each mouth. Freshwater into the sea = P \u2212 E + river water (blue = freshening, brown = net evaporation).",
                 "'Why fresher / saltier' compares the HadCM3L surface salinity with its ocean mean (\u00b10.3 psu) and the local",
                 "freshwater flux (\u00b1100 mm/yr): where the local flux does not explain it, the water was brought by currents,",
                 "mixing, upwelling or sea ice. Tick 'rivers and river mouths' to draw the rivers on any layer (e.g. salinity).")),
      helpText("P \u2212 E and runoff are shown in mm per year for the annual mean, mm per day otherwise. Salinity, mixed layer and",
               "currents are from the HadCM3L ocean; up- / downwelling is the HadCM3L ocean's."),
      radioButtons("climate_method", "Climate",
                   c("HadCM3L + local correction for your map (fast)" = "local", "HadCM3L as run" = "model"), selected = "local"),
      helpText("The map of the age selected in the Map editor tab is drawn."),
      tags$details(tags$summary("Where the climate comes from"),
        helpText(HTML(CLIMATE_SOURCE_HTML))),
      hr(),
      h5("What is drawn"),
      checkboxGroupInput("map_show", NULL,
                         c("reefs" = "reef", "evaporites" = "evaporite", "glacial deposits" = "glacial",
                           "your points" = "user", "country borders + coastlines" = "borders",
                           "equator" = "equator", "ice front (from glacial deposits)" = "ice_front",
                           "original PaleoCoastlines shoreline (dashed)" = "shoreline",
                           "north arrow" = "north_arrow", "scale bar" = "scale_bar",
                           "wind arrows (10 m)" = "wind", "ocean current arrows" = "currents",
                           "sea ice (permanent / seasonal) on every layer" = "sea_ice",
                           "rivers and river mouths (runoff routed over the map)" = "rivers",
                           "outline where your geography differs from the climate model's" = "climate_flag"),
                         selected = c("reef", "evaporite", "glacial", "user", "borders", "equator", "ice_front")),
      selectInput("projection", "Map projection",
                  c("Mollweide (equal-area)" = "mollweide", "Robinson" = "robinson", "Equal Earth (equal-area)" = "equal_earth",
                    "Winkel tripel" = "winkel_tripel", "Hammer (equal-area)" = "hammer", "Eckert IV (equal-area)" = "eckert4",
                    "Plate Carr\u00e9e (lat-lon grid)" = "plate_carree", "Orthographic globe (centred on the map centre)" = "orthographic",
                    "South polar (equal-area)" = "south_polar", "North polar (equal-area)" = "north_polar"),
                  selected = SILUCARB$map_projection),
      conditionalPanel("input.projection == 'south_polar' || input.projection == 'north_polar'",
        numericInput("polar_edge", "Polar map reaches to latitude (0 = equator; + = north)", SILUCARB$polar_edge_lat, -60, 60, 5)),
      helpText("The projection is used for whole-map renders; a zoomed-in render uses its own equal-area view."),
      radioButtons("colour_bar_pos", "Elevation / depth colour bar",
                   c("bottom" = "bottom", "top" = "top", "left" = "left", "right" = "right"),
                   selected = SILUCARB$colour_bar_pos, inline = TRUE),
      radioButtons("quality", "Rendered map quality",
                   c("Draft (fast)" = "draft", "Standard" = "standard", "High resolution (slow)" = "high"),
                   selected = SILUCARB$map_quality, inline = TRUE),
      helpText("Only affects the rendered PNG / PDF, not the editor view. High: 300 dpi and zoomed-in",
               "maps refined up to 8x - takes longer and uses more memory. Saved in points/settings.csv."),
      selectizeInput("localities", "Localities on the map", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "none")),
      tags$details(tags$summary(tags$b("Add a locality")),
        textInput("fig_loc_name", "Name", placeholder = "e.g. Gotland"),
        fluidRow(column(6, numericInput("fig_loc_lat", "Latitude (\u00b0N)", NA, -90, 90, 0.01)),
                 column(6, numericInput("fig_loc_lon", "Longitude (\u00b0E)", NA, -180, 360, 0.01))),
        radioButtons("fig_loc_kind", NULL, c("present-day coordinates" = "present", "palaeo-coordinates at the map age" = "paleo"), "present"),
        actionButton("fig_loc_add", "Add locality", icon = icon("location-dot")),
        helpText("South and west are negative. Palaeo-coordinates are taken at the map age and plate model of the Map editor and stored as present-day coordinates, so the locality moves with its plate at other ages. A locality with the same name is replaced. Remove localities in the Localities tab (Map editor).")),
      helpText("Add or remove localities here, in the Localities tab, or click 'Add locality here' on the editor map."),
      selectInput("center", "Centre the map on", c("default (150\u00b0W)" = "", "longitude below" = "lon")),
      checkboxInput("render_zoom", "Render only a zoomed-in area", FALSE),
      conditionalPanel("input.render_zoom",
        radioButtons("zoom_size", "Area to render",
                     c("Area shown in the Figure map (drag / zoom it)" = "fig", "Same as the Map editor window" = "window",
                       "Set width and height (km)" = "custom"), selected = "fig"),
        textOutput("zoom_window_txt"),
        conditionalPanel("input.zoom_size == 'custom'",
          fluidRow(column(6, numericInput("zoom_w_km", "Width (km)", 3000, 100, 20000, 100)),
                   column(6, numericInput("zoom_h_km", "Height (km)", 2000, 100, 20000, 100))),
          helpText("Centred on the 'Centre the map on' locality when one is chosen, otherwise on the middle",
                   "of the editor window; drawn in an equal-area projection.")),
        checkboxInput("overview", "Also show the whole earth, with a box around this area joined to the zoom-in", FALSE),
        conditionalPanel("input.overview",
          selectInput("overview_corner", "Corner of the zoom-in (partly over the whole-earth map)",
                      c("automatic (away from the box)" = "auto", "bottom right" = "bottomright", "bottom left" = "bottomleft",
                        "top right" = "topright", "top left" = "topleft")),
          fluidRow(column(6, numericInput("overview_ratio", "Zoom-in size \u00d7 whole earth", 1.25, 0.1, 4, 0.25)),
                   column(6, numericInput("overview_overlap", "Overlap (0 = corners touch)", 0.15, -0.5, 0.9, 0.05))),
          checkboxGroupInput("overview_show", "Also draw on the whole-earth map",
                             c("localities (with names)" = "sites", "your points" = "user",
                               "reefs / evaporites / glacial deposits" = "indicators"), inline = TRUE),
          helpText("Only what is ticked under 'What is drawn' / 'Localities' is drawn; the legend stays with the zoom-in.",
                   "Size: 1 = both maps the same size, 2 = zoom-in twice as large, 0.7 = smaller.",
                   "The whole-earth map uses the projection chosen above."))),
      conditionalPanel("input.center == 'lon'", numericInput("center_lon", "Central palaeo-longitude", 0, -180, 180)),
      hr(),
      actionButton("fig_render", "Render figure", class = "btn-primary", icon = icon("image"))
    ),
    mainPanel(width = 9,
      div(style = "margin-bottom:4px",
          actionButton("fz_in", NULL, icon = icon("magnifying-glass-plus"), title = "Zoom in"),
          actionButton("fz_out", NULL, icon = icon("magnifying-glass-minus"), title = "Zoom out"),
          actionButton("fz_reset", NULL, icon = icon("earth-europe"), title = "Whole map"),
          actionButton("fz_w", NULL, icon = icon("arrow-left")), actionButton("fz_e", NULL, icon = icon("arrow-right")),
          actionButton("fz_n", NULL, icon = icon("arrow-up")), actionButton("fz_s", NULL, icon = icon("arrow-down")),
          span(style = "color:#666; margin-left:8px",
               "Hold and drag to move the map \u00b7 mouse wheel or +/\u2212 to zoom \u00b7 the area shown is what a zoomed-in or composite figure renders")),
      plotOutput("fig_map", height = "560px"),
      div(style = "color:#555; margin:3px 0 8px 0", textOutput("fig_win_txt", inline = TRUE)),
      verbatimTextOutput("fig_msg"),
      div(downloadButton("fig_pdf", "Download PDF"), downloadButton("fig_png", "Download PNG")),
      br(),
      div(id = "fig_preview", imageOutput("fig_img", height = "auto"))
    )
  )),
  PaleoMapR:::pm_about_tab()
)

server <- function(input, output, session) {
  PaleoMapR:::pm_server_extras(input, session)
  rv <- reactiveValues(version = 0, zoom = NULL, last = NULL, msg = "Ready.", gridv = 0, locv = 0, panv = 0)

  # ---- localities -------------------------------------------------------------
  locs <- reactive({ rv$locv; read_localities() })
  observe({
    l <- locs()
    updateSelectizeInput(session, "localities", choices = l$name, selected = l$name[l$show])
    updateSelectInput(session, "center", choices = c("default (150\u00b0W)" = "", setNames(l$name, l$name),
                                                     "longitude below" = "lon"),
                      selected = isolate(if (is.null(input$center)) "" else input$center))
  })
  observeEvent(input$localities, {                 # remember which localities are shown
    l <- read_localities(); new <- l$name %in% input$localities
    if (!identical(new, l$show)) { l$show <- new; e <- locality_write_msg(write_localities(l)); if (!is.null(e)) rv$msg <- e }
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  output$loctab <- renderDT(datatable(locs(), rownames = FALSE, options = list(pageLength = 8, dom = "tip")))
  observeEvent(input$loc_add, {
    req(nzchar(input$loc_name), is.finite(input$loc_lat), is.finite(input$loc_lon))
    e <- locality_write_msg(add_locality(input$loc_name, input$loc_lat, input$loc_lon))
    if (!is.null(e)) { rv$msg <- e; return() }
    rv$locv <- rv$locv + 1
    rv$msg <- paste("Locality added:", input$loc_name)
  })
  observeEvent(input$fig_loc_add, {                # add a locality from the Figure tab
    say <- function(m) { rv$msg <- m; fig$msg <- m }
    nm <- trimws(input$fig_loc_name %||% ""); la <- input$fig_loc_lat; lo <- input$fig_loc_lon
    if (!nzchar(nm)) return(say("Give the locality a name."))
    if (!is.numeric(la) || !is.numeric(lo) || !is.finite(la) || !is.finite(lo) || abs(la) > 90) return(say("Give a latitude (-90 to 90) and a longitude."))
    lo <- ((lo + 180) %% 360) - 180
    note <- ""
    if (identical(input$fig_loc_kind, "paleo")) {
      age <- as.numeric(input$age)
      pp <- tryCatch(gp_to_present(la, lo, age, input$model), error = function(e) NULL)
      if (is.null(pp) || !is.finite(pp$lat[1])) return(say(sprintf("%.2f\u00b0N, %.2f\u00b0E is not on a plate at %s Ma in this plate model; give present-day coordinates instead.", la, lo, age)))
      note <- sprintf(" (palaeo %.2f, %.2f at %s Ma)", la, lo, age); la <- pp$lat[1]; lo <- pp$lon[1]
    }
    had <- nm %in% locs()$name
    e <- locality_write_msg(add_locality(nm, round(la, 4), round(lo, 4)))
    if (!is.null(e)) return(say(e))
    rv$locv <- rv$locv + 1
    updateTextInput(session, "fig_loc_name", value = "")
    updateNumericInput(session, "fig_loc_lat", value = NA); updateNumericInput(session, "fig_loc_lon", value = NA)
    say(sprintf("Locality '%s' %s: present-day %.2f\u00b0N, %.2f\u00b0E%s. It is shown on the maps (render the figure again).", nm, if (had) "replaced" else "added", la, lo, note))
  })
  observeEvent(input$loc_del, {
    s <- input$loctab_rows_selected; req(length(s) > 0)
    e <- locality_write_msg(remove_locality(locs()$name[s]))
    if (!is.null(e)) { rv$msg <- e; return() }
    rv$locv <- rv$locv + 1
    rv$msg <- sprintf("Removed %d locality(ies).", length(s))
  })

  # ---- land-point settings (saved to points/settings.csv) ----------------------
  observeEvent(input$colour_bar_pos, {
    if (!identical(input$colour_bar_pos, apply_user_settings(SILUCARB)$colour_bar_pos))
      save_user_setting("colour_bar_pos", input$colour_bar_pos, SILUCARB)
  }, ignoreInit = TRUE)
  observeEvent(input$projection, {
    if (!identical(input$projection, apply_user_settings(SILUCARB)$map_projection))
      save_user_setting("map_projection", input$projection, SILUCARB)
  }, ignoreInit = TRUE)
  observeEvent(input$polar_edge, {
    req(is.finite(input$polar_edge))
    if (input$polar_edge != apply_user_settings(SILUCARB)$polar_edge_lat) save_user_setting("polar_edge_lat", input$polar_edge, SILUCARB)
  }, ignoreInit = TRUE)
  observeEvent(input$quality, {
    if (!identical(input$quality, apply_user_settings(SILUCARB)$map_quality))
      save_user_setting("map_quality", input$quality, SILUCARB)
  }, ignoreInit = TRUE)
  observeEvent(input$land_radius, {
    req(is.finite(input$land_radius), input$land_radius > 0)
    if (input$land_radius != apply_user_settings(SILUCARB)$LAND_POINT_RADIUS_KM)
      save_user_setting("LAND_POINT_RADIUS_KM", input$land_radius, SILUCARB)
  }, ignoreInit = TRUE)
  observeEvent(input$land_link, {
    req(is.finite(input$land_link), input$land_link >= 0)
    if (input$land_link != apply_user_settings(SILUCARB)$LAND_LINK_KM)
      save_user_setting("LAND_LINK_KM", input$land_link, SILUCARB)
  }, ignoreInit = TRUE)
  bump <- function(msg) { rv$version <- rv$version + 1; rv$msg <- msg }
  observeEvent(input$age, { updateNumericInput(session, "age_young", value = as.numeric(input$age))
                            updateNumericInput(session, "age_old", value = as.numeric(input$age)) })

  # ---- data -----------------------------------------------------------------
  recs_all <- reactive({ rv$version; all_records(include_inactive = TRUE) })
  shown <- reactive({
    age <- as.numeric(input$age); r <- recs_all()
    a <- records_for_age(age, input$model, "dem", recs = r)
    b <- records_for_age(age, input$model, "indicator", recs = r)
    d <- rbind(a, b[!b$record_id %in% a$record_id, names(a)])
    d
  })
  grid <- reactive({
    rv$gridv; age <- as.numeric(input$age)
    dir <- if (input$grid == "corrected" && length(list.files(SILUCARB$dir_corrected, paste0("_", age, "Ma\\.nc$"))))
      SILUCARB$dir_corrected else SILUCARB$dir_dem
    shiny::validate(shiny::need(length(list.files(dir, paste0("_", age, "Ma\\.nc$"))) > 0,
                  paste0("No PaleoDEM for ", age, " Ma in ", dir, ". Run download_data() in the console ",
                         "(or download paleoDEM_0.2Deg_grids.zip from EarthByte by hand and unzip it there).")))
    editor_grid(find_dem(dir, age), age, input$model, label = if (identical(dir, SILUCARB$dir_corrected)) { "corrected" } else { "original" },
                cfg = SILUCARB)                     # kept in memory: switching back to a map is instant
  })

  # ---- present-day borders / coastlines rotated to the map age ---------------
  modern <- function(kind) {
    L <- tryCatch(gp_modern_lines(ne_file(kind), input$model, as.numeric(input$age)),
                  error = function(e) { rv$msg <- paste("Could not load present-day lines:", conditionMessage(e)); list() })
    L
  }
  draw_modern <- function(L, col, lwd) {
    if (!length(L)) return(invisible())
    young <- vapply(L, function(p) isTRUE(attr(p, "young")), TRUE)
    for (grp in list(list(sel = !young, lty = 1, a = 0.95), list(sel = young, lty = 3, a = 0.75))) {
      if (!any(grp$sel)) next
      xy <- do.call(rbind, lapply(L[grp$sel], function(p) {
        seg <- cumsum(c(FALSE, abs(diff(p$lon)) > 180))       # break where a line crosses 180 deg
        do.call(rbind, lapply(split(seq_along(p$lon), seg), function(k) rbind(cbind(p$lon[k], p$lat[k]), c(NA, NA))))
      }))
      lines(xy[, 1], xy[, 2], col = adjustcolor(col, grp$a), lwd = lwd, lty = grp$lty)
    }
  }

  # ---- map ------------------------------------------------------------------
  loc_cache <- new.env()
  output$map <- renderPlot({
    g <- grid(); d <- shown()
    xl <- if (is.null(rv$zoom)) c(-180, 180) else rv$zoom$x
    yl <- if (is.null(rv$zoom)) c(-90, 90) else rv$zoom$y
    par(mar = c(2.5, 2.5, 1.5, 0.5), mgp = c(1.4, 0.4, 0), cex.axis = 0.8)
    rv$panv
    plot(NA, xlim = xl, ylim = yl, xaxs = "i", yaxs = "i", asp = 1, xlab = "palaeo-longitude", ylab = "palaeo-latitude",
         main = sprintf("%s Ma \u2013 %s grid \u2013 %s", input$age, g$label, SILUCARB$models[input$model]), cex.main = 0.9)
    geo$usr <- par("usr"); geo$plt <- par("plt")      # for converting mouse positions to map coordinates
    if (!identical(isolate(rv$usr), geo$usr)) rv$usr <- geo$usr   # visible area (for 'render the editor area')
    dl <- abs(g$lat[2] - g$lat[1])
    rasterImage(g$img[nrow(g$img):1, , ], min(g$lon) - dl / 2, min(g$lat) - dl / 2, max(g$lon) + dl / 2, max(g$lat) + dl / 2,
                interpolate = FALSE)
    cxy <- editor_coast_xy(g$coast)                 # 0 m line, computed once per map
    lines(cxy$x, cxy$y, lwd = 0.7, col = "#222222")
    abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60), col = adjustcolor("white", 0.4), lwd = 0.5)
    if ("equator" %in% input$map_show) abline(h = 0, col = "#ffd54a", lwd = 1.3, lty = 2)
    if ("coastlines" %in% input$show) draw_modern(modern("coast"), "#ffffff", 0.9)
    if ("borders" %in% input$show) draw_modern(modern("borders"), "#000000", 0.8)
    if (!nrow(d)) return()
    act <- d[d$active, ]; exc <- d[!d$active, ]
    for (k in seq_len(nrow(SETTING_STYLE))) {
      s <- SETTING_STYLE[k, ]
      if (s$setting != "" && !s$setting %in% input$show) next
      e <- act[act$setting == s$setting & act$source != "user", ]
      if (s$setting == "" && !"indicators" %in% input$show) next
      points(e$plon, e$plat, pch = s$pch, bg = s$bg, col = "#111111", cex = 0.8, lwd = 0.4)
    }
    if ("indicators" %in% input$show) {
      e <- act[act$indicator != "", ]
      points(e$plon, e$plat, pch = ifelse(e$indicator == "reef", 1, ifelse(e$indicator == "evaporite", 0, 5)),
             col = ifelse(e$indicator == "reef", "#d6337a", ifelse(e$indicator == "evaporite", "#b36b00", "#ffffff")),
             cex = 1.4, lwd = 1.2)
    }
    if ("excluded" %in% input$show && nrow(exc)) points(exc$plon, exc$plat, pch = 4, col = "#ff2020", cex = 1, lwd = 1.5)
    u <- act[act$source == "user", ]
    if (nrow(u)) points(u$plon, u$plat, pch = 24, bg = "#7fffd4", col = "#003b2e", cex = 1.3, lwd = 1)
    if (!is.null(rv$last)) points(rv$last[1], rv$last[2], pch = 3, col = "magenta", cex = 2, lwd = 2)
    l <- locs(); l <- l[l$name %in% input$localities, ]
    if (nrow(l)) {
      lkey <- paste(input$age, input$model, paste(l$name, l$lat, l$lon, collapse = ";"))
      pl <- loc_cache[[lkey]]
      if (is.null(pl)) {                            # rotate the localities once per map (pyGPlates), not at every redraw
        pl <- gp_reconstruct(l$lat, l$lon, rep(as.numeric(input$age), nrow(l)), input$model)
        assign(lkey, pl, envir = loc_cache)
      }
      points(pl$plon, pl$plat, pch = 21, bg = "#e34948", col = "white", cex = 1.4, lwd = 1.5)
      text(pl$plon, pl$plat, l$name, pos = 4, cex = 0.8, font = 2, col = "white")
    }
    legend("bottomleft", c("marine (shallow)", "marine (deeper shelf)", "evaporite", "terrestrial",
                           "reef / evaporite / glacial symbol", "excluded", "your point", "present-day border",
                           "present-day coastline (white)"),
           pch = c(21, 21, 22, 24, 1, 4, 24, NA, NA), pt.bg = c("#2f7fd8", "#0b2a6e", "#eda100", "#a0522d", NA, NA, "#7fffd4", NA, NA),
           lty = c(NA, NA, NA, NA, NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, NA, NA, NA, NA, 1, 1.2),
           col = c(1, 1, 1, 1, "#d6337a", "#ff2020", "#003b2e", "#000000", "#9a9a9a"),
           bg = adjustcolor("white", 0.85), cex = 0.75)
  })
  # ---- navigation: drag to pan, wheel / buttons to zoom ------------------------
  geo <- new.env()
  to_map <- function(fx, fy) {                     # image fraction -> palaeo lon / lat
    u <- geo$usr; p <- geo$plt
    list(x = u[1] + (fx - p[1]) / (p[2] - p[1]) * (u[2] - u[1]),
         y = u[3] + ((1 - fy) - p[3]) / (p[4] - p[3]) * (u[4] - u[3]))
  }
  # the area visible in the editor (palaeo-coordinates): centre and true size in km
  window_box <- function(any_zoom = FALSE) {
    u <- rv$usr
    if ((is.null(rv$zoom) && !any_zoom) || is.null(u)) return(NULL)
    x <- pmin(pmax(u[1:2], -180), 180); y <- pmin(pmax(u[3:4], -90), 90)
    la <- mean(y); kmd <- 6371 * pi / 180
    list(center_palaeo = c(la, mean(x)), width_km = round(diff(x) * kmd * cos(la * pi / 180)),
         height_km = round(diff(y) * kmd))
  }
  output$zoom_window_txt <- renderText({
    input$render_zoom; fdrawn()
    if (input$zoom_size %in% c("fig", "box")) return(if (is.null(fig_region())) "Zoom in on the Figure map (mouse wheel or +) to choose the area." else
                                                      "The area shown in the Figure map is rendered (see below the map).")
    wb <- window_box()
    if (identical(input$zoom_size, "custom") && !is.null(input$center) && !input$center %in% c("", "lon"))
      return(sprintf("The area is centred on %s ('Centre the map on').", input$center))
    if (is.null(wb)) { if (identical(input$zoom_size, "custom")) "The area is centred on the middle of the editor window."
                       else "Zoom in on the editor map first (mouse wheel or +)." } else
      sprintf("Editor window: %s \u00d7 %s km around %.1f\u00b0%s, %.1f\u00b0%s", format(wb$width_km, big.mark = ","),
              format(wb$height_km, big.mark = ","), abs(wb$center_palaeo[1]), if (wb$center_palaeo[1] < 0) "S" else "N",
              abs(wb$center_palaeo[2]), if (wb$center_palaeo[2] < 0) "W" else "E")
  })
  set_view <- function(x, y) {                     # keep the window inside the globe
    w <- min(diff(x), 360); h <- min(diff(y), 180)
    x <- mean(x) + c(-w, w) / 2; y <- mean(y) + c(-h, h) / 2
    if (x[1] < -180) x <- x - (x[1] + 180); if (x[2] > 180) x <- x - (x[2] - 180)
    if (y[1] < -90) y <- y - (y[1] + 90); if (y[2] > 90) y <- y - (y[2] - 90)
    rv$zoom <- if (w >= 359.9 && h >= 179.9) NULL else list(x = x, y = y)
    rv$panv <- rv$panv + 1
  }
  zoom_by <- function(f, cx = NULL, cy = NULL) {
    req(geo$usr); u <- geo$usr
    if (is.null(cx)) { cx <- mean(u[1:2]); cy <- mean(u[3:4]) }
    set_view(cx + (u[1:2] - cx) / f, cy + (u[3:4] - cy) / f)
  }
  pan_by <- function(fx, fy) {                     # fractions of the plot width / height
    req(geo$usr); u <- geo$usr; p <- geo$plt
    sx <- fx / (p[2] - p[1]) * (u[2] - u[1]); sy <- fy / (p[4] - p[3]) * (u[4] - u[3])
    set_view(u[1:2] - sx, u[3:4] + sy)
  }
  observeEvent(input$map_pan, pan_by(input$map_pan$dx, input$map_pan$dy))
  observeEvent(input$map_wheel, {
    w <- input$map_wheel; c_ <- to_map(w$fx, w$fy)
    zoom_by(if (w$steps > 0) 1.5 else 1 / 1.5, c_$x, c_$y)
  })
  observeEvent(input$zoom_in, zoom_by(1.5))
  observeEvent(input$zoom_out, zoom_by(1 / 1.5))
  observeEvent(input$zoom_reset, { rv$zoom <- NULL; rv$panv <- rv$panv + 1 })
  observeEvent(input$pan_w, pan_by(-0.25, 0)); observeEvent(input$pan_e, pan_by(0.25, 0))
  observeEvent(input$pan_n, pan_by(0, -0.25)); observeEvent(input$pan_s, pan_by(0, 0.25))

  nearest <- function(x, y, d) {
    if (!nrow(d)) return(NULL)
    dd <- sqrt((d$plat - y)^2 + ((((d$plon - x) + 180) %% 360 - 180) * cos(y * pi / 180))^2)
    span <- if (is.null(rv$zoom)) 360 else diff(rv$zoom$x)
    o <- order(dd); o <- o[dd[o] <= max(0.6, span / 60)]
    if (length(o)) d[o, ] else NULL
  }

  # ---- clicks ---------------------------------------------------------------
  observeEvent(input$map_click2, {
    req(geo$usr); c_ <- to_map(input$map_click2$fx, input$map_click2$fy)
    x <- c_$x; y <- c_$y; age <- as.numeric(input$age)
    if (x < -180 || x > 180 || y < -90 || y > 90) return()
    rv$last <- c(x, y)
    m <- input$mode
    if (m == "delete_point") {
      u <- shown(); u <- u[u$source == "user", ]
      near <- nearest(x, y, u)
      if (is.null(near)) { rv$msg <- "None of your points near the click."; return() }
      r <- near[1, ]
      showModal(modalDialog(title = "Delete point permanently?",
        sprintf("%s (%s%s, %s\u2013%s Ma) will be removed from user_points.csv.", r$record_id, r$setting,
                if (nzchar(r$indicator)) paste0(", ", r$indicator) else "", r$age_young, r$age_old),
        footer = tagList(modalButton("Cancel"), actionButton("confirm_delete", "Delete", class = "btn-danger"))))
      rv$to_delete <- r$record_id
      return()
    }
    if (m == "add_locality") {
      nm <- trimws(input$note)
      if (!nzchar(nm)) { rv$msg <- "Type the locality name in the 'Note' box first, then click on the map."; return() }
      pp <- gp_to_present(y, x, age, input$model)
      if (!is.finite(pp$lat)) { rv$msg <- "That spot is not on a plate at this age; add the locality with present-day coordinates in the Localities tab."; return() }
      e <- locality_write_msg(add_locality(nm, round(pp$lat, 4), round(pp$lon, 4)))
      if (!is.null(e)) { rv$msg <- e; return() }
      rv$locv <- rv$locv + 1
      rv$msg <- sprintf("Locality '%s' added (present-day %.2f, %.2f).", nm, pp$lat, pp$lon)
      return()
    }
    if (startsWith(m, "add_")) {
      kind <- sub("add_", "", m)
      setting <- switch(kind, marine = "marine", marine_deep = "marine_deep", evaporite = "evaporite",
                        terrestrial = "terrestrial", reef = "marine", glacial = "")
      indicator <- switch(kind, reef = "reef", glacial = "glacial", evaporite = "evaporite", "")
      ay <- min(input$age_young, input$age_old); ao <- max(input$age_young, input$age_old)
      id <- add_user_point(lat = y, lon = x, setting = setting, indicator = indicator, coords = "palaeo",
                           model = input$model, clicked_age = age, age_young = ay, age_old = ao, note = input$note)
      bump(sprintf("Added %s point %s at %.2f, %.2f (palaeo, %s Ma). Rebuild this age to update the grid.",
                   kind, id, y, x, age))
    } else {
      near <- nearest(x, y, shown())
      if (is.null(near)) { rv$msg <- "No record near the click (zoom in for precise picking)."; return() }
      r <- near[1, ]
      if (m == "toggle") {
        act <- if (r$active) "exclude" else "include"
        edit_record(r$record_id, act, note = input$note)
        bump(sprintf("%sd %s (%s, %s). Rebuild this age to update the grid.", tools::toTitleCase(act),
                     r$record_id, r$setting, r$name))
      } else rv$msg <- sprintf("%s: %s | setting %s | indicator %s | %s\u2013%s Ma | %s",
                               r$record_id, r$name, r$setting, r$indicator, r$age_young, r$age_old,
                               if (r$active) "active" else "EXCLUDED")
    }
  })

  output$near <- renderDT({
    req(rv$last); n <- nearest(rv$last[1], rv$last[2], shown())
    if (is.null(n)) return(datatable(data.frame(info = "no records near the last click")))
    datatable(head(n[, c("record_id", "source", "name", "setting", "indicator", "age_young", "age_old",
                         "active", "edited", "lat", "lon", "plat", "plon")], 25), rownames = FALSE,
              options = list(pageLength = 5, dom = "tip", scrollX = TRUE)) |>
      formatRound(c("lat", "lon", "plat", "plon"), 2)
  })

  # ---- tables of your changes -----------------------------------------------
  output$upoints <- renderDT({ rv$version; datatable(read_user_points(), rownames = FALSE,
                                                     options = list(pageLength = 8, scrollX = TRUE)) })
  output$uedits <- renderDT({ rv$version; datatable(read_record_edits(), rownames = FALSE,
                                                    options = list(pageLength = 8, scrollX = TRUE)) })
  observeEvent(input$confirm_delete, {
    removeModal(); req(rv$to_delete)
    delete_user_point(rv$to_delete); bump(paste("Deleted", rv$to_delete, "permanently. Rebuild this age to update the grid."))
    rv$to_delete <- NULL
  })
  observeEvent(input$del_points, {
    s <- input$upoints_rows_selected; req(length(s) > 0)
    u <- read_user_points(); write.csv(u[-s, ], points_file("user_points.csv"), row.names = FALSE)
    bump(sprintf("Deleted %d point(s).", length(s)))
  })
  observeEvent(input$del_edits, {
    s <- input$uedits_rows_selected; req(length(s) > 0)
    e <- read_record_edits(); write.csv(e[-s, ], points_file("record_edits.csv"), row.names = FALSE)
    bump(sprintf("Removed %d edit(s).", length(s)))
  })
  observeEvent(input$undo, {
    # undo = remove the most recent row of whichever table changed last
    fu <- points_file("user_points.csv"); fe <- points_file("record_edits.csv")
    tu <- file.mtime(fu); te <- file.mtime(fe)
    if (!is.na(te) && (is.na(tu) || te >= tu) && nrow(read_record_edits())) {
      e <- read_record_edits(); write.csv(e[-nrow(e), ], fe, row.names = FALSE); bump("Undid the last edit.")
    } else if (nrow(read_user_points())) {
      u <- read_user_points(); write.csv(u[-nrow(u), ], fu, row.names = FALSE); bump("Removed the last point.")
    }
  })

  # ---- rebuild / render -----------------------------------------------------
  observeEvent(input$rebuild, {
    age <- as.numeric(input$age)
    withProgress(message = sprintf("Rebuilding %d Ma (about 1 minute)...", age), value = 0.3, {
      st <- run_dem_correction(age, cfg = apply_user_settings(SILUCARB), workers = 1)
    })
    rv$gridv <- rv$gridv + 1
    s <- st[st$map_age == age, ]
    rv$msg <- sprintf("%d Ma rebuilt: %d marine / %d evaporite / %d terrestrial records, %d seaways, %d isolated, ocean %.1f Mkm2 (raw %.1f).",
                      age, s$n_marine, s$n_evaporite, s$n_terrestrial, s$seaways_added,
                      s$marine_water_bodies_still_isolated, s$ocean_corrected_Mkm2, s$ocean_DEM_Mkm2)
  })
  observeEvent(input$rebuild_all, {
    withProgress(message = "Rebuilding all ages (parallel)...", value = 0.3, { run_dem_correction(cfg = apply_user_settings(SILUCARB)) })
    rv$gridv <- rv$gridv + 1; rv$msg <- "All ages rebuilt."
  })
  observeEvent(input$fig_kind, {                    # layout -> the zoom-in and overview settings
    updateCheckboxInput(session, "render_zoom", value = input$fig_kind != "whole")
    updateCheckboxInput(session, "overview", value = input$fig_kind == "composite")
  }, ignoreInit = TRUE)
  observeEvent(list(input$render_zoom, input$overview), {
    k <- if (!isTRUE(input$render_zoom)) "whole" else if (isTRUE(input$overview)) "composite" else "zoom"
    if (!identical(k, input$fig_kind)) updateRadioButtons(session, "fig_kind", selected = k)
  }, ignoreInit = TRUE)
  # ---- Figure tab: interactive map (same layers and overlays as the rendered figure) ----
  figv <- reactiveValues(cx = NULL, cy = 0, w = 360, h = 180, v = 0)
  fgeo <- new.env(); fdrawn <- reactiveVal(0)
  fig_ctr <- reactive(if (input$center == "") NULL else if (input$center == "lon") input$center_lon else input$center)
  fig_setup <- reactive({
    src <- fig_source(function(t) fig$msg <- t); req(src)
    ms <- input$map_show; run <- src$run
    withProgress(message = "Preparing the Figure map ...", value = 0.4,
      tryCatch(figure_preview_setup(src$age, model = if (is.null(run)) input$model else "scotese_and_wright2018",
                                    corrected = src$corr, layer = input$fig_layer %||% "topography", season = input$season %||% "ann",
                                    climate_method = input$climate_method %||% "local",
                                    experiment = if (!is.null(run)) run$x else NULL, diff_with = src$cmp,
                                    map = if (!is.null(run)) run$map else NULL, ice = src$ice_arg,
                                    sea_ice = "sea_ice" %in% ms, rivers = "rivers" %in% ms, wind = "wind" %in% ms,
                                    currents = "currents" %in% ms, climate_flag = "climate_flag" %in% ms,
                                    o2_pal = if (isTRUE(is.finite(input$o2_pal))) input$o2_pal else 1,
                                    lapse_rate = if (isTRUE(is.finite(input$lapse_rate))) input$lapse_rate else 6.5,
                                    shoreline = "shoreline" %in% ms, cfg = SILUCARB),
               error = function(e) { fig$msg <- paste("Figure map:", conditionMessage(e)); NULL }))
  })
  fig_window <- function() {
    cx <- figv$cx %||% tryCatch(map_center(fig_ctr(), as.numeric(input$age), input$model, SILUCARB), error = function(e) SILUCARB$central_lon)
    list(cx = cx, cy = figv$cy, w = figv$w, h = figv$h)
  }
  output$fig_map <- renderPlot({
    S <- fig_setup(); req(S); figv$v
    ms <- input$map_show
    ind <- intersect(ms, c("reef", "evaporite", "glacial")); if (!length(ind)) ind <- FALSE
    locsel <- input$localities; if (!length(locsel)) locsel <- FALSE
    r <- tryCatch(figure_preview(S, fig_window(), borders = "borders" %in% ms, indicators = ind, user_points = "user" %in% ms,
                                 sites = locsel, equator = "equator" %in% ms, ice_front = "ice_front" %in% ms,
                                 north_arrow = "north_arrow" %in% ms, scale_bar = "scale_bar" %in% ms,
                                 width_px = session$clientData$output_fig_map_width %||% 1000, cfg = SILUCARB),
                  error = function(e) { plot.new(); text(0.5, 0.5, paste("Figure map:", conditionMessage(e))); NULL })
    fgeo$r <- r
    fdrawn(isolate(fdrawn()) + 1)                    # the shown area is known now (for the text below the map)
  })
  fig_to_ll <- function(fx, fy) {                  # image fraction -> palaeo lon / lat
    r <- fgeo$r; req(r); u <- r$usr; p <- r$plt
    x <- u[1] + (fx - p[1]) / (p[2] - p[1]) * (u[2] - u[1]); y <- u[3] + ((1 - fy) - p[3]) / (p[4] - p[3]) * (u[4] - u[3])
    c(lon = r$cx + x / r$k, lat = y / r$k)
  }
  fig_visible <- function() {                      # the area shown (lon / lat ranges)
    r <- fgeo$r; if (is.null(r)) return(NULL)
    list(lon = r$cx + r$usr[1:2] / r$k, lat = pmin(pmax(r$usr[3:4] / r$k, -90), 90))
  }
  fig_set <- function(cx, cy, w, h) {
    w <- min(max(w, 1), 360); h <- min(max(h, 0.5), 180)
    if (w >= 359.9 && h >= 179.9) { cy <- 0 }
    cy <- min(max(cy, -90 + h / 2), 90 - h / 2)
    figv$cx <- ((cx + 180) %% 360) - 180; figv$cy <- cy; figv$w <- w; figv$h <- h; figv$v <- figv$v + 1
  }
  fig_zoom <- function(f, p = NULL) {
    w <- fig_window(); vis <- fig_visible(); if (!is.null(vis)) { w$w <- min(360, diff(vis$lon)); w$h <- diff(vis$lat) }
    if (is.null(p)) p <- c(lon = w$cx, lat = w$cy)
    dlon <- ((w$cx - p[["lon"]] + 180) %% 360) - 180
    fig_set(p[["lon"]] + dlon / f, p[["lat"]] + (w$cy - p[["lat"]]) / f, w$w / f, w$h / f)
  }
  fig_pan <- function(fx, fy) {                    # fractions of the image width / height
    w <- fig_window(); vis <- fig_visible(); if (is.null(vis)) return()
    fig_set(w$cx - fx * diff(vis$lon), w$cy + fy * diff(vis$lat), w$w, w$h)
  }
  observeEvent(input$fig_map_pan, fig_pan(input$fig_map_pan$dx, input$fig_map_pan$dy))
  observeEvent(input$fig_map_wheel, { k <- input$fig_map_wheel; fig_zoom(if (k$steps > 0) 1.5 else 1 / 1.5, fig_to_ll(k$fx, k$fy)) })
  observeEvent(input$fz_in, fig_zoom(1.5)); observeEvent(input$fz_out, fig_zoom(1 / 1.5))
  observeEvent(input$fz_reset, { figv$cx <- NULL; figv$cy <- 0; figv$w <- 360; figv$h <- 180; figv$v <- figv$v + 1 })
  observeEvent(input$fz_w, fig_pan(0.25, 0)); observeEvent(input$fz_e, fig_pan(-0.25, 0))
  observeEvent(input$fz_n, fig_pan(0, 0.25)); observeEvent(input$fz_s, fig_pan(0, -0.25))
  observeEvent(input$center, { figv$cx <- NULL; figv$v <- figv$v + 1 }, ignoreInit = TRUE)
  #' the area shown as a region for make_maps (NULL when the whole map is shown)
  fig_region <- function() {
    vis <- fig_visible(); if (is.null(vis) || figv$w >= 359.9) return(NULL)
    lo <- ((vis$lon + 180) %% 360) - 180
    c(round(vis$lat[1], 2), round(vis$lat[2], 2), round(lo[1], 2), round(lo[2], 2))
  }
  output$fig_win_txt <- renderText({
    fdrawn(); vis <- fig_visible()
    if (is.null(vis) || figv$w >= 359.9) return("Whole map. Zoom in to choose the area of a zoomed-in or composite figure.")
    la <- mean(vis$lat); kmd <- 6371 * pi / 180
    sprintf("Area shown: %s\u2013%s, %s\u2013%s (about %s \u00d7 %s km). A zoomed-in or composite figure renders this area (equal-area view centred on it).",
            .deg_ns(round(vis$lat[1], 1)), .deg_ns(round(vis$lat[2], 1)), .deg_ew(round(vis$lon[1], 1)), .deg_ew(round(vis$lon[2], 1)),
            format(round(diff(vis$lon) * kmd * cos(la * pi / 180)), big.mark = ","), format(round(diff(vis$lat) * kmd), big.mark = ","))
  })
  # the map of the figure (NULL + a message when not possible)
  fig_source <- function(say) {
    age <- as.numeric(input$age)
    run <- NULL; cmp <- NULL; ice_arg <- NULL           # (climate runs and ice sheets: full Silucarb version only)
    corr <- input$grid == "corrected"
    if (is.null(run) && corr && !length(list.files(SILUCARB$dir_corrected, paste0("_", age, "Ma\\.nc$")))) {
      say(paste0("No corrected grid for ", age, " Ma yet: press 'Rebuild this age' first, ",
                 "or choose 'Original PaleoDEM' as background grid.")); return(NULL)
    }
    list(age = age, run = run, cmp = cmp, ice_arg = ice_arg, corr = corr)
  }
  # ---- figure: one function for 'Render map' (editor) and 'Render figure' (Figure tab) ----
  fig <- reactiveValues(msg = "Choose the settings and press 'Render figure'.", pdf = NULL, png = NULL)
  render_figure <- function() {
    say <- function(t) { rv$msg <- t; fig$msg <- t }
    age <- as.numeric(input$age)
    ms <- input$map_show
    ind <- intersect(ms, c("reef", "evaporite", "glacial")); if (!length(ind)) ind <- FALSE
    ctr <- if (input$center == "") NULL else if (input$center == "lon") input$center_lon else input$center
    locsel <- input$localities; if (!length(locsel)) locsel <- FALSE
    src <- fig_source(say); if (is.null(src)) return()
    age <- src$age; run <- src$run; cmp <- src$cmp; ice_arg <- src$ice_arg; corr <- src$corr
    reg <- NULL
    if (isTRUE(input$render_zoom)) {
      if (input$zoom_size %in% c("fig", "box")) {
        reg <- fig_region()
        if (is.null(reg)) { say("Zoom in on the Figure map first (mouse wheel or +) to choose the area, or choose another 'Area to render'."); return() }
      } else if (identical(input$zoom_size, "custom")) {
        w <- input$zoom_w_km; h <- input$zoom_h_km
        if (!isTRUE(is.finite(w) && is.finite(h) && w > 0 && h > 0)) { say("Give a width and height in km."); return() }
        loc_ctr <- !is.null(input$center) && !input$center %in% c("", "lon")
        if (loc_ctr) {
          reg <- list(center = input$center, width_km = w, height_km = h)        # follows the locality
        } else {
          wb <- window_box(any_zoom = TRUE)
          if (is.null(wb)) { say("The editor map is not drawn yet; wait a moment and try again."); return() }
          reg <- list(center_palaeo = wb$center_palaeo, width_km = w, height_km = h)
        }
      } else {
        reg <- window_box()
        if (is.null(reg)) {
          say(paste("Zoom in on the editor map first (mouse wheel or +), or choose 'Set width and height (km)'",
                    "to render an area around the 'Centre the map on' locality.")); return()
        }
      }
    }
    num <- function(v, d) if (isTRUE(is.finite(v))) v else d
    f <- NULL
    trace_ <- NULL                                  # where a rendering error happened (for the message + log)
    withProgress(message = "Rendering map...", value = 0.3, {
      f <- tryCatch(withCallingHandlers(make_maps(age, model = if (is.null(run)) input$model else "scotese_and_wright2018",
                     corrected = corr, comparison = FALSE, indicators = ind,
                     region = reg, user_points = "user" %in% ms, sites = locsel, borders = "borders" %in% ms,
                     shoreline = "shoreline" %in% ms, equator = "equator" %in% ms, center = ctr,
                     north_arrow = "north_arrow" %in% ms, scale_bar = "scale_bar" %in% ms,
                     ice_front = "ice_front" %in% ms, projection = input$projection,
                     overview = !is.null(reg) && isTRUE(input$overview),
                     overview_corner = input$overview_corner %||% "auto",
                     overview_ratio = max(0.1, num(input$overview_ratio, 1.25)),
                     overview_overlap = num(input$overview_overlap, 0.15),
                     overview_show = input$overview_show %||% character(0),
                     polar_edge_lat = num(input$polar_edge, 0),
                     colour_bar_pos = input$colour_bar_pos %||% "bottom",
                     layer = input$fig_layer %||% "topography", wind = "wind" %in% ms, currents = "currents" %in% ms,
                     season = input$season %||% "ann", climate_method = input$climate_method %||% "local",
                     climate_flag = "climate_flag" %in% ms,
                     map = if (!is.null(run)) run$map else NULL,
                     experiment = if (!is.null(run)) run$x else NULL, diff_with = cmp, ice = ice_arg,
                     seasons4 = isTRUE(input$fig_seasons4),
                     sea_ice = "sea_ice" %in% ms, rivers = "rivers" %in% ms, o2_pal = if (isTRUE(is.finite(input$o2_pal))) input$o2_pal else 1,
                     lapse_rate = if (isTRUE(is.finite(input$lapse_rate))) input$lapse_rate else 6.5,
                     quality = input$quality, cfg = SILUCARB),
                    error = function(e) trace_ <<- vapply(sys.calls(), function(cl) paste(deparse(cl, nlines = 2), collapse = " "), "")),
                    error = function(e) {
                      tr <- trace_ %||% character(0)
                      i0 <- tail(grep("^make_maps\\(", tr), 1); i1 <- grep("^\\.handleSimpleError\\(", tr)[1]
                      if (length(i0)) tr <- tr[i0:(if (is.na(i1)) length(tr) else i1 - 1)]
                      tr <- tr[!grepl("^(tryCatch|withCallingHandlers|doTryCatch|tryCatchList|tryCatchOne|withVisible|eval|stop|h)\\(", tr)]
                      where <- if (length(tr)) sub("\\(.*", "", tail(tr, 1)) else "?"
                      log <- file.path(SILUCARB$dir_maps, "last_render_error.txt")
                      try(writeLines(c(format(Sys.time()), paste("Error:", conditionMessage(e)), "", "Calls (outer to inner):",
                                       substr(tr, 1, 300), "", "Inputs:",
                                       capture.output(str(lapply(reactiveValuesToList(input)[intersect(sort(c("age", "model", "grid", "center", "center_lon",
                                         "localities", "map_show", "projection", "fig_layer", "fig_kind", "render_zoom", "zoom_size", "season",
                                         "fig_seasons4", "climate_method", "fig_run", "fig_compare", "quality", "overview", "colour_bar_pos",
                                         "polar_edge", "overview_ratio", "overview_overlap", "overview_corner", "overview_show")), names(input))], function(v) v)))), log), silent = TRUE)
                      say(sprintf("Rendering failed in %s(): %s  (details: maps/last_render_error.txt)", where, conditionMessage(e))); NULL })
    })
    if (!is.null(f)) { say(paste("Map written:", f)); fig$pdf <- as.character(f); fig$png <- tail(attr(f, "pngs"), 1); fig$v <- (fig$v %||% 0) + 1 }   # v: redraw even when the file name is unchanged
  }
  observeEvent(input$render, render_figure())
  observeEvent(input$fig_render, render_figure())
  output$fig_msg <- renderText(fig$msg)
  output$fig_img <- renderImage({
    fig$v; req(fig$png, file.exists(fig$png))
    list(src = fig$png, contentType = "image/png", width = "100%", alt = "rendered map")
  }, deleteFile = FALSE)
  output$fig_pdf <- downloadHandler(filename = function() basename(fig$pdf %||% "map.pdf"),
                                    content = function(file) file.copy(fig$pdf, file))
  output$fig_png <- downloadHandler(filename = function() basename(fig$png %||% "map.png"),
                                    content = function(file) file.copy(fig$png, file))
  output$stats <- renderTable({
    rv$gridv; f <- file.path(SILUCARB$dir_results, "step6_corrected_dem_stats.csv")
    if (!file.exists(f)) return(NULL)
    read.csv(f)[, c("map_age", "n_marine", "n_evaporite", "n_terrestrial", "n_user_points", "seaways_added",
                    "marine_water_bodies_still_isolated", "seaway_max_narrowing_km", "slope_max_changed_after",
                    "ocean_DEM_Mkm2", "ocean_corrected_Mkm2")]
  })
  output$status <- renderText(rv$msg)
}

shinyApp(ui, server)
