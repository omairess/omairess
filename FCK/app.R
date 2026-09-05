# ==============================================================================
# F*CK — Functional data analysis, Circadian regression, K-means clustering
# ==============================================================================
# One integrated Shiny app merging two previously separate apps:
#
#   WaPaa1_3.R  "Functional Data Analysis Suite"
#                 -> fPCA / time-warped PCA, functional ANOVA (between and
#                    repeated measures), post-hoc tests, functional clustering
#   CIRCAREG.R  "Functional Regression Suite"
#                 -> function-on-scalar regression, scalar-on-function
#                    regression, harmonic (cosinor) regression, pairwise
#                    comparisons of circadian parameters
#
# The two apps had four steps in common — file import, variable selection,
# smoothing, and smoothing diagnostics — implemented twice, with the same
# input ids and (for smoothing) very nearly the same code.  In this app those
# four steps exist ONCE.  Every analysis downstream of them reads the same
# `values$data`, `values$smooth_data`, `values$fd_obj`, `values$covariates`
# and `values$time_numeric`, so a curve is imported and smoothed one way and
# then analysed by any of the seven analysis tabs.
#
# Everything downstream of smoothing is the original code, carried across by
# line range (see tools/port_fck.py and PORTING_NOTES.md), so each analysis
# still computes and prints exactly what it did in its own app.
#
# Layout
#   app.R            this file: packages, UI assembly, server assembly
#   ui/              one file per tab; each defines a `ui_tab_*` object
#   server/          server body, split by section, sourced with local = TRUE
#                    into ONE server environment — so every section sees the
#                    same `values`, the same helpers, and the same `input`,
#                    exactly as in the two single-file originals.
# ==============================================================================

# --- Packages ---------------------------------------------------------------
# Required: the app cannot start without these.
required_packages <- c("shiny", "shinydashboard", "shinyWidgets", "fda", "mgcv",
                       "plotly", "DT", "dplyr", "tidyr", "ggplot2", "cluster",
                       "readxl")

# Optional: each backs ONE feature.  In the two source apps a failed install of
# any of these stopped the whole app from starting.  Merging makes that blast
# radius much larger (a missing reticulate would have taken the cosinor tabs
# down with it), so these are loaded if present and reported if not; the
# feature that needs one reports the missing package itself when used.
optional_packages <- c(
  gridExtra    = "multi-panel plot export",
  viridis      = "continuous colour scales",
  # rmfanova was listed as optional for repeated-measures fANOVA, but the app
  # never called its API correctly and now does not call it at all (P1.b). Left
  # out rather than installed for a code path that does not exist.
  fda.usc      = "functional k-means clustering",
  reticulate   = "DCF (density-core-finding) clustering via Python",
  refund       = "scalar-on-function regression (refund::pfr)",
  minpack.lm   = "robust exponential-saturation cosinor fits"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

missing_optional <- character(0)
for (pkg in names(optional_packages)) {
  ok <- suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))
  if (!ok) {
    ok <- tryCatch({
      install.packages(pkg, dependencies = TRUE)
      suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))
    }, error = function(e) FALSE)
  }
  if (!isTRUE(ok)) missing_optional <- c(missing_optional, pkg)
}
if (length(missing_optional)) {
  message("\n--- F*CK: optional packages not available ---")
  for (pkg in missing_optional)
    message(sprintf("  %-12s -> disables: %s", pkg, optional_packages[[pkg]]))
  message("Everything else works. Install with install.packages(\"<pkg>\").\n")
}

# --- Source the UI and remember the server files ----------------------------
# The UI files run now (they build tab objects).  The server files must run
# INSIDE the server function, once per session, so they are only listed here.
# source() re-encodes to the native locale, which mangles (or refuses) the UTF-8
# in these files when R is started in a non-UTF-8 locale. parse(encoding=) does
# not, so every file in this app is loaded through here.
fck_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir)
  invisible(NULL)
}

fck_dir <- function(sub) {
  d <- sub
  if (!dir.exists(d)) d <- file.path("FCK", sub)          # run from repo root
  if (!dir.exists(d)) stop("F*CK: cannot find the '", sub, "' directory. ",
                           "Start the app with shiny::runApp(\"FCK\").")
  d
}

for (f in sort(list.files(fck_dir("ui"), full.names = TRUE, pattern = "[.]R$")))
  fck_source(f, envir = globalenv())

FCK_SERVER_FILES <- sort(list.files(fck_dir("server"), full.names = TRUE,
                                    pattern = "[.]R$"))

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  dashboardHeader(title = "F*CK", titleWidth = 300),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "sidebar",

      # ---- shared pipeline: done once, feeds every analysis below ----------
      menuItem("Data Import", tabName = "import", icon = icon("upload")),
      menuItem("Data Preprocessing/Smoothing", tabName = "preprocess",
               icon = icon("cogs")),
      menuItem("Smoothing Diagnostics", tabName = "smooth_diag",
               icon = icon("chart-area")),

      # ---- F: functional data analysis (from WaPaa) -----------------------
      menuItem("fPCA/time-warped PCA Settings", tabName = "settings",
               icon = icon("sliders-h")),
      menuItem("Functional PCA Results", tabName = "results",
               icon = icon("chart-line")),
      menuItem("Functional ANOVA", tabName = "fanova", icon = icon("chart-bar")),
      menuItem("fANOVA: post-hoc tests", tabName = "pairwise",
               icon = icon("exchange-alt")),

      # ---- C: circadian / functional regression (from CIRCAREG) -----------
      menuItem("Function-on-Scalar (FoSR)", tabName = "fosr",
               icon = icon("chart-line")),
      menuItem("Scalar-on-Function (SoFR)", tabName = "sofr",
               icon = icon("chart-bar")),
      menuItem("Harmonic Regression", tabName = "harmonic", icon = icon("sync")),
      menuItem("Cosinor: pairwise tests", tabName = "harm_pairwise",
               icon = icon("not-equal")),

      # ---- K: clustering (from WaPaa) -------------------------------------
      menuItem("Functional Clustering", tabName = "kmeans",
               icon = icon("project-diagram")),

      menuItem("Data Export", tabName = "export", icon = icon("download"))
    )
  ),

  dashboardBody(
    ui_theme_css,
    tabItems(
      ui_tab_import,
      ui_tab_preprocess,
      ui_tab_smooth_diag,
      ui_tab_settings,
      ui_tab_results,
      ui_tab_fanova,
      ui_tab_posthoc,
      ui_tab_fosr,
      ui_tab_sofr,
      ui_tab_harmonic,
      ui_tab_cosinor_pairwise,
      ui_tab_clustering,
      ui_tab_export
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
# Every server/*.R file is sourced with local = TRUE into THIS function's
# environment.  That is what makes the merge work without rewriting 14 000
# lines of analysis code: helpers defined in one file (the clock-time helpers,
# the cosinor fitters, the warping functions) are visible to every other file,
# exactly as when each app was a single script.
server <- function(input, output, session) {
  cat("===== F*CK SERVER STARTED =====\n")
  for (f in FCK_SERVER_FILES) fck_source(f, envir = environment())
  cat("===== F*CK SERVER SETUP COMPLETE =====\n")
}

shinyApp(ui = ui, server = server)
