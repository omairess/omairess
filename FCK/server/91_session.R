# ==============================================================================
# server/91_session.R — save and restore a whole session   (hand-written)
#
# Neither source app could be reopened. Close the browser and the import, the
# variable selection, the smoothing and every fitted model were gone; the CSV
# exports gave you the numbers but nothing that could be loaded back.
#
# WHAT IS SAVED
#   * the entire values bus: raw data, analysis matrix, covariates and group
#     variables, the smoothed curves and the fd object, and every fitted model
#     (fPCA, warping, fANOVA, post-hoc, clustering, FoSR, cosinor)
#   * the analysis-defining settings, so a restored session can be re-run and
#     not merely re-read
#   * the exact package versions, because an fda or mgcv update can move the
#     numbers a saved model is sitting next to
#
# WHAT IS NOT RESTORED
#   Only the settings listed in RESTORE_INPUTS are pushed back into the
#   widgets. Cosmetic controls (tick spacing, which subject a plot shows) come
#   back at their defaults. Restoring every input generically means guessing
#   each widget's message format, which fails silently on the ones it gets
#   wrong — worse than a documented, checked subset.
# ==============================================================================

# Settings that DEFINE an analysis, restored with the proper update*Input call.
# Anything not here is cosmetic; the results themselves come back via `values`.
RESTORE_INPUTS <- list(
  numeric  = c("n_basis", "n_basis_manual", "smooth_factor", "min_bound",
               "max_bound", "n_components", "n_permutations", "alpha_level",
               "harmonic_period", "n_harmonics", "n_clusters", "n_boot"),
  checkbox = c("is_cyclic", "use_real_time", "constrain_bounds",
               "use_bootstrap", "cv_stratified"),
  select   = c("smooth_method", "pca_type", "fanova_design", "fanova_test_type",
               "fanova_data_source", "fanova_group_var", "harmonic_time_var",
               "harmonic_trend_type", "harmonic_group_var", "reg_method",
               "cluster_method", "pairwise_correction", "hp_param",
               "hp_correction", "data_format"),
  text     = c("harmonic_manual_times")
)

FCK_SESSION_FORMAT <- 1L

# shiny does not export %||%, and base R only gained it in 4.4.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

fck_package_versions <- function() {
  pkgs <- c("shiny", "shinydashboard", "shinyWidgets", "fda", "mgcv", "plotly",
            "DT", "dplyr", "tidyr", "ggplot2", "cluster", "readxl",
            "rmfanova", "fda.usc", "reticulate", "minpack.lm")
  vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                    error = function(e) "not installed"),
         character(1))
}

session_note <- reactiveVal(
  "No session loaded. Saving writes the data, the smoothing and every result you have run.")

# --- Save --------------------------------------------------------------------
output$save_session <- downloadHandler(
  filename = function() sprintf("fck_session_%s.rds", format(Sys.time(), "%Y%m%d_%H%M%S")),
  content = function(file) {
    # reactiveValuesToList gives the whole bus in one go, so a field added to
    # 00_state.R later is saved without anyone remembering to list it here.
    state <- shiny::reactiveValuesToList(values)

    settings <- list()
    for (nm in unlist(RESTORE_INPUTS, use.names = FALSE)) {
      v <- tryCatch(input[[nm]], error = function(e) NULL)
      if (!is.null(v)) settings[[nm]] <- v
    }

    # Every input, for the record. Not restored — read it to see what was set.
    all_inputs <- tryCatch(
      shiny::reactiveValuesToList(input),
      error = function(e) list())
    all_inputs <- all_inputs[vapply(all_inputs, function(v)
      is.atomic(v) && length(v) <= 50, logical(1))]

    saveRDS(list(
      format       = FCK_SESSION_FORMAT,
      app          = "F*CK",
      saved_at     = Sys.time(),
      r_version    = R.version.string,
      packages     = fck_package_versions(),
      values       = state,
      settings     = settings,
      all_inputs   = all_inputs
    ), file)
  }
)

# --- Restore -----------------------------------------------------------------
observeEvent(input$load_session, {
  req(input$load_session)

  saved <- tryCatch(readRDS(input$load_session$datapath),
                    error = function(e) e)
  if (inherits(saved, "error")) {
    showNotification(paste("Could not read that file:", conditionMessage(saved)),
                     type = "error", duration = 10)
    return()
  }
  if (!is.list(saved) || is.null(saved$format) || is.null(saved$values)) {
    showNotification("That .rds is not an F*CK session file.", type = "error", duration = 10)
    return()
  }
  if (saved$format > FCK_SESSION_FORMAT) {
    showNotification(
      sprintf("This session was written by a newer version of the app (format %d > %d). Update before opening it.",
              saved$format, FCK_SESSION_FORMAT),
      type = "error", duration = 12)
    return()
  }

  for (nm in names(saved$values)) values[[nm]] <- saved$values[[nm]]

  restored <- character(0)
  push <- function(names_vec, fn) {
    for (nm in intersect(names_vec, names(saved$settings))) {
      tryCatch({
        fn(session, nm, value = saved$settings[[nm]])
        restored <<- c(restored, nm)
      }, error = function(e) NULL)
    }
  }
  push(RESTORE_INPUTS$numeric,  updateNumericInput)
  push(RESTORE_INPUTS$checkbox, updateCheckboxInput)
  push(RESTORE_INPUTS$select,   updateSelectInput)
  push(RESTORE_INPUTS$text,     updateTextInput)

  # Versions the results were produced under vs. the ones running now. A model
  # object restored beside a different fda is exactly the case worth flagging.
  now  <- fck_package_versions()
  then <- saved$packages
  drift <- character(0)
  if (!is.null(then)) {
    shared <- intersect(names(then), names(now))
    drift <- shared[then[shared] != now[shared] &
                    then[shared] != "not installed" &
                    now[shared]  != "not installed"]
  }

  ran <- c(
    if (!is.null(values$smooth_data))       "smoothing",
    if (!is.null(values$pca_results))       "functional PCA",
    if (!is.null(values$warping_results))   "time warping",
    if (!is.null(values$fanova_results))    "functional ANOVA",
    if (!is.null(values$pairwise_results))  "post-hoc tests",
    if (!is.null(values$clustering_results))"clustering",
    if (!is.null(values$reg_model))         "FoSR",
    if (!is.null(values$harmonic_model))    "cosinor")

  session_note(paste(c(
    sprintf("Loaded session saved %s", format(saved$saved_at, "%Y-%m-%d %H:%M:%S")),
    sprintf("  saved under: %s", saved$r_version %||% "unknown R version"),
    sprintf("  data:        %s",
            if (is.null(values$data)) "none"
            else sprintf("%d subjects x %d time points", nrow(values$data), ncol(values$data))),
    sprintf("  results:     %s", if (length(ran)) paste(ran, collapse = ", ") else "none"),
    sprintf("  settings restored: %d of %d",
            length(restored), length(unlist(RESTORE_INPUTS, use.names = FALSE))),
    "",
    if (length(drift)) c(
      "PACKAGE VERSIONS HAVE CHANGED since these results were computed:",
      sprintf("  %-12s saved %-10s now %s", drift, then[drift], now[drift]),
      "Re-run the analyses if you intend to report these numbers.")
    else "Package versions match the ones these results were computed with.",
    "",
    "Cosmetic display controls are at their defaults; every result above is the loaded one."
  ), collapse = "\n"))

  showNotification(
    sprintf("Session restored (%s).",
            if (length(ran)) paste(ran, collapse = ", ") else "data only"),
    type = "message", duration = 8)
})

output$session_status <- renderText(session_note())
