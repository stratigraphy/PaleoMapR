# =============================================================================
#  Version and data sources of the paleogeographic map editor (shown in the app,
#  README.md and CITATION.md).
# =============================================================================
APP_VERSION <- tryCatch(as.character(utils::packageVersion("PaleoMapR")), error = function(e) { "dev" })

#' where the HadCM3L climate fields come from (verified 20 Oct 2026 against the BRIDGE run
#' list of the series and the PhanDA model-prior record on Zenodo)
CLIMATE_SOURCE_HTML <- paste0(
  "The climate layers come from the <b>HadCM3L</b> coupled atmosphere&ndash;ocean simulations ",
  "(HadCM3BL-M2.1aD) of <b>Valdes, Scotese &amp; Lunt (2021)</b>, series <b>scotese_02</b>: 109 runs, ",
  "541&ndash;0 Ma, on the PALEOMAP palaeogeographies of Scotese &amp; Wright (2018) (the same maps as the ",
  "PaleoDEMs) with CO<sub>2</sub> after Foster, Royer &amp; Lunt (2017). Run IDs <i>texqe</i> (541 Ma) to ",
  "<i>texpa1</i> (0 Ma); e.g. <i>teXPj</i> 436 Ma, <i>teXPh</i> 425 Ma, <i>teXPg</i> 421 Ma. For each map age ",
  "the run closest in age is used; its climatological means (atmosphere <i>.pdcl</i>, ocean <i>.pfcl</i>) are ",
  "downloaded from the BRIDGE server of the University of Bristol into <code>data/climate/</code>.<br>",
  "Please cite: Valdes, P. J., Scotese, C. R. and Lunt, D. J. (2021) Deep ocean temperatures through time. ",
  "<i>Climate of the Past</i> 17, 1483&ndash;1506, ",
  "<a href='https://doi.org/10.5194/cp-17-1483-2021' target='_blank'>doi:10.5194/cp-17-1483-2021</a>; ",
  "model: Valdes et al. (2017) <i>Geosci. Model Dev.</i> 10, 3715&ndash;3743, ",
  "<a href='https://doi.org/10.5194/gmd-10-3715-2017' target='_blank'>doi:10.5194/gmd-10-3715-2017</a>. ",
  "Run list: <a href='https://www.paleo.bristol.ac.uk/ummodel/scripts/html_bridge/scotese_02.html' target='_blank'>",
  "BRIDGE scotese_02</a>. 'Local correction' is this editor's own statistical adjustment to your map, not part ",
  "of the HadCM3L output.")
