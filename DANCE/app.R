# ==============================================================================
# DANCE — Direct Analysis of multivariate Nonlinear and Circadian Effects in R
# ==============================================================================
# One integrated Shiny app merging two previously separate apps:
#
#   WaPaa1_3.R  "Functional Data Analysis Suite"
#                 -> fPCA / time-warped PCA, functional ANOVA (between and
#                    repeated measures), post-hoc tests, functional clustering
#   CIRCAREG.R  "Functional Regression Suite"
#                 -> function-on-scalar regression, harmonic (cosinor)
#                    regression, pairwise comparisons of circadian parameters
#
# The two apps had four steps in common — file import, variable selection,
# smoothing, and smoothing diagnostics — implemented twice, with the same
# input ids and (for smoothing) very nearly the same code.  In this app those
# four steps exist ONCE.  Every analysis downstream of them reads the same
# `values$data`, `values$smooth_data`, `values$fd_obj`, `values$covariates`
# and `values$time_numeric`, so a curve is imported and smoothed one way and
# then analysed by any of the analysis tabs.
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
  minpack.lm   = "robust exponential-saturation cosinor fits",
  # The mixed-design tab. mgcv is already required (the FoSR GAM branch needs
  # it), so only lme4 is added here -- and it is optional for the same reason
  # as the rest: a missing lme4 must cost you the mixed cosinor, not the app.
  lme4         = "mixed (between x within) cosinor",
  # P21 phase 4. The trajectory framework in tab 6 rests on these, and none of
  # them was declared when it was built -- so the app started cleanly, the
  # fitted curves drew, and the first sign of a missing emmeans was a terse
  # refusal inside a panel the user had already waited for. Each is optional for
  # the same reason as the rest (a missing one must cost you a feature, not the
  # app), but it has to be REPORTED at startup, which is what declaring it here
  # does. What each one actually costs when absent:
  emmeans      = "per-cell amplitude and acrophase, and every contrast built on them",
  pbkrtest     = "Kenward-Roger degrees of freedom (inference falls back to Satterthwaite)",
  lmerTest     = "Satterthwaite degrees of freedom (inference falls back to a likelihood-ratio test)",
  glmmTMB      = "residual correlation structures (AR(1) and continuous-time OU)"
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
# The app now checks and OFFERS TO INSTALL, with consent, and stops with the
# exact command when it cannot.
#
# AMENDED (P21 phase 4). "Stop with the command" was right about the three
# hazards and wrong about the remedy: being told to install four packages, one
# refusal at a time, as you discover each feature that needs one, is not a
# five-second fix -- it is four restarts. What the reasons above actually argue
# against is installing SILENTLY and AUTOMATICALLY, not installing at all. So
# the app now asks:
#
#   interactive session   names the packages, says what each one costs you if
#                         absent, asks once, and installs only on a yes. The
#                         installed versions are printed, so what changed is on
#                         the record rather than inferred later.
#   anything else         deployed app, Rscript, CI, or DANCE_NO_INSTALL=1:
#                         behaves exactly as before -- no install, stop with
#                         the command. Reason 2 above is untouched: a machine
#                         that cannot install still fails immediately and says
#                         why, instead of pausing to try.
#
# That keeps all three properties -- nothing silent, nothing attempted where it
# cannot work, the missing package named before anything happens -- while
# costing one keystroke instead of a hunt.
#
# AUDIT (P4.9): this note used to point the reader at a lockfile in the project
# root. THERE IS NO SUCH FILE, and a reviewer was right to say so. (The old
# sentence is deliberately not quoted here. Three times in this audit a fix has
# reproduced the exact string it was removing, inside the comment explaining the
# removal, which makes the grep guard that proves the removal pass for the wrong
# reason. Describe the removed text; do not repeat it.) Pointing at a
# lockfile that does not exist is the same category of defect as the rest of
# this audit -- a claim of a guarantee that is not there. A lockfile records a
# library that EXISTS, and it can only be written on a machine where the app's
# packages are actually installed; the audit container is not that machine (see
# tools/renv_bootstrap.R, which refuses to snapshot an incomplete library
# rather than writing a fiction). So: pin the environment yourself, once, on
# the machine you analyse on --
#
#   Rscript tools/renv_bootstrap.R      # writes renv.lock; commit it
#   renv::restore()                     # later, or elsewhere
#
# -- and until you have done that, this project is NOT environment-pinned.
# Ask, install on a yes, and report what landed. Returns the packages still
# missing afterwards. Never installs without an answer, and never asks where an
# answer cannot be given.
dance_offer_install <- function(pkgs, why = NULL) {
  if (!length(pkgs)) return(character(0))
  can_ask <- interactive() && !identical(Sys.getenv("DANCE_NO_INSTALL"), "1")
  if (!can_ask) return(pkgs)
  message("\nDANCE needs ", length(pkgs), " package(s) that are not installed:")
  for (p in pkgs)
    message("  ", p, if (!is.null(why[[p]])) paste0("  -- ", why[[p]]) else "")
  message("\nThey will be installed into: ", .libPaths()[1])
  ans <- tolower(trimws(readline("Install them now? [y/N] ")))
  if (!ans %in% c("y", "yes")) {
    message("Not installing. Re-run when you are ready, or:  install.packages(c(",
            paste(sprintf('"%s"', pkgs), collapse = ", "), "))")
    return(pkgs)
  }
  utils::install.packages(pkgs)
  still <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  got <- setdiff(pkgs, still)
  if (length(got))
    message("Installed: ", paste(sprintf("%s %s", got,
            vapply(got, function(p) as.character(utils::packageVersion(p)), character(1))),
            collapse = ", "))
  if (length(still)) message("Still missing: ", paste(still, collapse = ", "))
  still
}

missing_required <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required))
  missing_required <- dance_offer_install(missing_required)
if (length(missing_required)) {
  stop("DANCE cannot start: required packages are missing.\n",
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
# THE SAME OFFER, ALL AT ONCE. Discovering four optional packages one refusal at
# a time -- each inside a panel you had already waited for -- is the failure
# that prompted this. They are offered together, here, each named with the
# feature it costs, before any analysis starts.
if (length(missing_optional))
  missing_optional <- dance_offer_install(missing_optional, optional_packages)
invisible(lapply(setdiff(names(optional_packages), missing_optional), function(p)
  suppressPackageStartupMessages(
    suppressWarnings(require(p, character.only = TRUE, quietly = TRUE)))))

if (length(missing_optional)) {
  message("\n--- DANCE: optional packages not available ---")
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
dance_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir)
  invisible(NULL)
}

dance_dir <- function(sub) {
  d <- sub
  if (!dir.exists(d)) d <- file.path("DANCE", sub)          # run from repo root
  if (!dir.exists(d)) stop("DANCE: cannot find the '", sub, "' directory. ",
                           "Start the app with shiny::runApp(\"DANCE\").")
  d
}

for (f in sort(list.files(dance_dir("ui"), full.names = TRUE, pattern = "[.]R$")))
  dance_source(f, envir = globalenv())

DANCE_SERVER_FILES <- sort(list.files(dance_dir("server"), full.names = TRUE,
                                    pattern = "[.]R$"))

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  dashboardHeader(title = "DANCE", titleWidth = 300),

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
  cat("===== DANCE SERVER STARTED =====\n")
  for (f in DANCE_SERVER_FILES) dance_source(f, envir = environment())
  cat("===== DANCE SERVER SETUP COMPLETE =====\n")
}

shinyApp(ui = ui, server = server)
