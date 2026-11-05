if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
# =============================================================================
#  Paleogeographic map editor - settings. Edit here; every other script reads these values.
# =============================================================================

SILUCARB <- list(
  # ---- folders (relative to the project folder, i.e. where run_workflow.R is) ----
  dir_dem       = "data/paleoDEM_0.2Deg_grids",     # PaleoDEMs (Scotese & Wright 2018), 0.2 deg
  dir_coast     = "data/paleocoastlines",           # PaleoCoastlines v7.2 (holds CS/ and CM/)
  dir_corrected = "data/paleoDEM_0.2Deg_corrected", # output: fossil-corrected PaleoDEMs
  dir_ne        = "data/naturalearth",              # present-day borders / coastlines
  dir_models    = "models",                         # plate models (plate-model-manager)
  dir_results   = "results",
  dir_maps      = "maps",
  dir_points    = "points",                         # editable point tables
  dir_R         = "R",
  dir_python    = "python",
  dir_tools     = "tools",                          # SPEEDY climate model (option "rerun")
  dir_map_files = "map_files",                      # map files loaded for climate runs (import_map)
  dir_exports   = "exports",                        # maps exported for use elsewhere (export_map)
  dir_reference = "reference",                      # reference data (CO2 curve of Foster et al. 2017)
  dir_data      = "data",
  file_pbdb     = "data/pbdb_silurian_colls.csv",
  file_sealevel = "data/vanderMeer2022_SuppTable1.xlsx",
  file_macro    = "data/macrostrat/units_global_410_450.json",

  # ---- python environment for pyGPlates (created by PaleoMapR::pm_setup()) ----
  conda_env = "paleomap",   # (an older environment called "silucarb" is used too)

  # ---- map ages ----
  ages = c(445, 440, 435, 430, 425, 420, 415),
  age_window = 2.5,   # Myr: a record is used for map age A if A lies within its
                      # age range widened by this much (PBDB: mid-age +- 2.5 Myr)

  # ---- record selection from the databases ----
  pbdb_max_age_range = 12,   # Myr: PBDB collections with a wider age range are not used
  macro_max_age_range = Inf, # Myr: Macrostrat units spanning more than this are not used (e.g. 30)
  terrestrial_env = "fluvial|lacustrine|terrestrial|delta plain|floodplain|eolian|channel|alluvial|swamp|cave|karst",
  evaporite_lith  = "evapor|gypsum|anhydrite|halite",
  reef_env        = "reef|buildup|bioherm",

  # ---- DEM correction (step 6) ----
  R_POINT               = 0.55,   # deg, fully flooded core around a marine point (~60 km)
  SEA_DEPTH             = -30,    # m, new shallow-marine depth ('marine' records)
  MARINE_DEEP_DEPTH     = -400,   # m, 'marine_deep' records: deeper shelf (not deep ocean)
  MARINE_DEEP_RADIUS_KM = 60,     # km, full depth within this distance of the point ...
  MARINE_DEEP_TAPER_KM  = 200,    # km, ... then a smooth S-shaped rise to the untouched shelf over this distance
  MARINE_DEEP_SMOOTH_KM = 40,     # km, extra gaussian smoothing of the lowering (no edges)
  EVAP_DEPTH            = -5,     # m, depth at evaporite points
  LAND_NEW              = 5,      # m, CS land where the DEM had sea
  SEAWAY_DEPTH          = -20,    # m, seaway depth at the inland end (deepens 30 m to the mouth)
  SEAWAY_HALF_KM_INLAND = 60,     # km, seaway half-width at its inland end
  SEAWAY_ASPECT         = 3,      # length : mouth width of a seaway (gulf-like 2-5)
  SEAWAY_MOUTH_MIN      = 3,      # mouth half-width >= this x inland half-width
  SEAWAY_HALF_KM_MAX    = 400,    # km, largest mouth half-width
  SEAWAY_SHAPE          = 0.6,    # widening profile exponent (<1 opens quickly)
  SEAWAY_BANK_KM        = 60,     # km, e-folding width of the lowered banks
  MEANDER               = 0.35,   # strength of the meander noise in the route cost
  NOISE_SEED            = 7,
  H_COST                = 25,     # m, elevation scale of the route cost
  ISLET_MAX_KM2         = 50000,  # km2, islets cut off by seaways that are drowned
  MAX_BANK_SLOPE        = 6,      # m/km, land next to new sea
  MAX_SEA_SLOPE         = 15,     # m/km, shelf next to changed coasts (to -1400 m)
  MAX_DEEP_SLOPE        = 25,     # m/km, continental slope below -1400 m
  MAX_LAND_SLOPE_FINAL  = 9,      # m/km, cone slope limiter on land in the changed zone
  COAST_SMOOTH_KM       = 35,     # km, coastline smoothing in the changed zone
  SPIT_HALF_KM          = 50,     # km, land strips thinner than ~2x this near new sea are removed
  SPIT_ZONE_KM          = 150,    # km, ... within this distance of new sea
  COAST_ELEV_SMOOTH_KM  = 25,     # km, elevation ramp across corrected coasts
  FINAL_COAST_SIGMA_KM      = 40,  # km, final 0 m line smoothing everywhere
  FINAL_COAST_SIGMA_CORR_KM = 100, # km, ... in and near the corrected areas
  FINAL_COAST_ZONE_KM       = 300, # km, extent of 'near the corrected areas'
  WIDEN_TOL_KM          = 8,      # km, tolerated narrowing before widening kicks in
  SEAFLOOR_BLUR_KM      = 30,     # km, final blur of the sea floor in and near corrected areas (no steps)
  # land records (terrestrial setting): kept as land through every step
  LAND_POINT_RADIUS_KM  = 50,     # km, land kept around each land record (and half-width of land bridges)
  LAND_POINT_ELEV       = 10,     # m, elevation at the edge of that land (low dome, +60 m in the centre)
  LAND_LINK_KM          = 300,    # km, land points closer than this are joined into one landmass, and
                                  #     joined to existing land closer than this
  LAND_MARINE_CLEAR_KM  = 40,     # km, no land bridge closer than this to a marine / evaporite record
  HARD_LAND_DATABASE    = TRUE,   # also database terrestrial records with no marine record within 1 deg
                                  # (your own terrestrial points are always kept)

  # ---- maps ----
  map_projection = "mollweide",   # mollweide, robinson, equal_earth, winkel_tripel, hammer, eckert4,
                                  # plate_carree, orthographic (globe), south_polar, north_polar
  colour_bar_pos = "bottom",      # elevation / depth colour bar: bottom, top, left or right
  polar_edge_lat = 0,             # polar maps reach to this latitude (0 = equator; 30 = to 30N for south_polar)
  central_lon   = -150,           # Mollweide central meridian (cut at ~30E keeps Gondwana whole)
  # map quality used by make_maps() and the editor's "Render map" button. The editor's own
  # preview always stays light; this only affects the rendered PNG / PDF files.
  #   "draft"    150 dpi, grid not refined                      (fastest)
  #   "standard" 220 dpi, zoomed-in maps refined up to 2x        (default)
  #   "high"     300 dpi, zoomed-in maps refined up to 8x, whole-globe maps 2x, and the land relief
  #              shading smoothed over relief_smooth_km (removes the 1-degree steps of the source
  #              PaleoDEM; display only - coastlines and the grids are not changed)  (slow)
  map_quality   = "standard",
  MAP_QUALITY   = list(draft    = list(dpi = 150, refine_max = 1, refine_global = 1, relief_smooth_km = 0),
                       standard = list(dpi = 220, refine_max = 2, refine_global = 1, relief_smooth_km = 0),
                       high     = list(dpi = 300, refine_max = 8, refine_global = 2, relief_smooth_km = 35)),
  map_dpi       = 220,            # set from map_quality by make_maps(); dpi = ... overrides it
  ne_resolution = "50m",          # Natural Earth detail: 10m, 50m, 110m
  # localities shown on the maps: edit points/localities.csv (or the app's Localities tab);
  # this is only the starting list written the first time
  sites = data.frame(name = c("Huy", "Kosov"), lat = c(50.52, 49.94), lon = c(5.24, 14.06)),
  stage = c("445" = "Late Ordovician (Katian–Hirnantian)", "440" = "Llandovery (Rhuddanian–Aeronian)",
            "435" = "Llandovery (Telychian)", "430" = "Wenlock", "425" = "Ludlow (Gorstian)",
            "420" = "Ludlow (Ludfordian) – Přídolí", "415" = "Early Devonian (Lochkovian)"),
  models = c(scotese_and_wright2018 = "PALEOMAP (Scotese & Wright 2018)",
             merdith2021 = "Merdith et al. 2021",
             torsvikcocks2017 = "Torsvik & Cocks 2017 (palaeomag. frame)",
             golonka = "Golonka (EarthByte)"),

  # ---- parallel processing (DEM correction, one age per worker) ----
  workers = max(1, min(7, parallel::detectCores() - 1)),
  keep_awake = TRUE                # keep Windows awake while climate runs go on in the background
)

#' make every dir_* / file_* setting an absolute path (so the scripts also work
#' when the working directory changes, e.g. inside the Shiny app)
absolutize_config <- function(cfg, root = getwd()) {
  for (k in grep("^(dir|file)_", names(cfg), value = TRUE))
    if (!grepl("^(/|[A-Za-z]:)", cfg[[k]])) cfg[[k]] <- normalizePath(file.path(root, cfg[[k]]), winslash = "/", mustWork = FALSE)
  cfg
}
SILUCARB <- absolutize_config(SILUCARB)

#' settings you change in the app are stored in points/settings.csv and override the
#' values above (delete a row, or the file, to go back to the defaults here)
USER_SETTING_KEYS <- c("workers", "keep_awake", "LAND_POINT_RADIUS_KM", "LAND_POINT_ELEV", "LAND_LINK_KM", "LAND_MARINE_CLEAR_KM",
                       "MARINE_DEEP_DEPTH", "MARINE_DEEP_RADIUS_KM", "MARINE_DEEP_TAPER_KM",
                       "HARD_LAND_DATABASE", "R_POINT", "FINAL_COAST_SIGMA_KM", "FINAL_COAST_SIGMA_CORR_KM",
                       "map_quality", "map_projection", "polar_edge_lat",
                       "colour_bar_pos")
apply_user_settings <- function(cfg = SILUCARB) {
  f <- file.path(cfg$dir_points, "settings.csv")
  if (!file.exists(f)) return(cfg)
  s <- read.csv(f, stringsAsFactors = FALSE)
  for (k in seq_len(nrow(s))) if (s$key[k] %in% USER_SETTING_KEYS) {
    v <- s$value[k]; old <- cfg[[s$key[k]]]
    cfg[[s$key[k]]] <- if (is.logical(old)) as.logical(v) else if (is.character(old)) as.character(v) else as.numeric(v)
  }
  cfg
}
save_user_setting <- function(key, value, cfg = SILUCARB) {
  stopifnot(key %in% USER_SETTING_KEYS)
  f <- file.path(cfg$dir_points, "settings.csv")
  s <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else data.frame(key = character(0), value = character(0))
  s <- s[s$key != key, , drop = FALSE]
  s <- rbind(s, data.frame(key = key, value = as.character(value)))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  write.csv(s, f, row.names = FALSE)
  invisible(apply_user_settings(cfg))
}
SILUCARB <- apply_user_settings(SILUCARB)
