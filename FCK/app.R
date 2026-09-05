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
  # rmfanova is optional again, and for the first time the app actually calls
  # its documented API (P2.6): the GLOBAL repeated-measures test, which the
  # app's own pointwise procedure does not provide. Note that rmfanova declares
  # refund in Imports but never calls it, so an install can fail for a
  # dependency nothing needs.
  rmfanova     = "global repeated-measures functional ANOVA test",
  fda.usc      = "functional k-means clustering",
  reticulate   = "DCF (density-core-finding) clustering via Python",
  minpack.lm   = "robust exponential-saturation cosinor fits"
)

# AUDIT (P2.2): this block used to call install.packages() at startup -- for the
# required packages unconditionally, and for each missing optional one inside a
# tryCatch. A statistical application must not rewrite its own library while
# launching. Three reasons, in order of how much they matter here:
#
#   1. It silently changes the analysis. Whatever CRAN holds on the day you
#      press run becomes the estimator. An analysis re-run next year against
#      newer fda or mgcv can produce different numbers with nothing in the
#      project recording why.
#   2. It cannot work where it is most needed. A deployed app, a locked-down
#      institutional machine or an offline analysis box has no writable library
#      and no network, so the "fallback" fails anyway -- after a long pause.
#   3. It hides the real problem. A missing package is a five-second fix once
#      you are told which one; discovering it through a half-installed
#      dependency tree is not.
#
# The app now checks and stops with the exact command to run. Pin the
# environment with renv (see renv.lock in the project root) so the library the
# analysis ran against is recorded with the analysis.
missing_required <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required)) {
  stop("F*CK cannot start: required packages are missing.\n",
       "  ", paste(missing_required, collapse = ", "), "\n\n",
       "Install them with:\n",
       "  install.packages(c(",
       paste(sprintf('"%s"', missing_required), collapse = ", "), "))\n\n",
       "Or restore the recorded environment with:  renv::restore()",
       call. = FALSE)
}
invisible(lapply(required_packages, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

missing_optional <- names(optional_packages)[
  !vapply(names(optional_packages), requireNamespace, logical(1), quietly = TRUE)]
invisible(lapply(setdiff(names(optional_packages), missing_optional), function(p)
  suppressPackageStartupMessages(
    suppressWarnings(require(p, character.only = TRUE, quietly = TRUE)))))

if (length(missing_optional)) {
  message("\n--- F*CK: optional packages not available ---")
  for (pkg in missing_optional)
    message(sprintf("  %-12s -> disables: %s", pkg, optional_packages[[pkg]]))
  message("Everything else works. Install with:")
  message(sprintf("  install.packages(c(%s))\n",
                  paste(sprintf('"%s"', missing_optional), collapse = ", ")))
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
