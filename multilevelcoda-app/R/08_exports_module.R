# ============================================================================
# 08_exports_module.R — pipeline log, reproducible script, data + results
# (house rules 6, 7, 8)
# ============================================================================

exportsModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(12,
        shiny::downloadButton(ns("script"), "Reproducible R script (.R)"),
        shiny::downloadButton(ns("log"), "Pipeline log (.txt)"),
        shiny::downloadButton(ns("ana_xlsx"), "Analysis data (.xlsx)"),
        shiny::downloadButton(ns("ana_csv"), "Analysis data (.csv)"),
        shiny::downloadButton(ns("diary_csv"), "All nights + sleep parameters (.csv)"),
        shiny::downloadButton(ns("bundle"), "Everything as R objects (.rds)"))),
    shiny::h4("Pipeline (in order)"),
    shiny::verbatimTextOutput(ns("log_view")),
    shiny::h4("Reproducible R script"),
    shiny::helpText("Stand-alone: it contains the import and diary functions, so",
                    "it runs in a fresh R session next to your data files."),
    shiny::verbatimTextOutput(ns("script_view"))
  )
}

exportsModuleServer <- function(id, rec, data_bus, comp, model, sub, helpers_path) {
  shiny::moduleServer(id, function(input, output, session) {
    script <- shiny::reactive(build_script(rec$steps(), helpers_path))
    logtxt <- shiny::reactive(build_pipeline_log(rec$steps()))
    output$script_view <- shiny::renderText(script())
    output$log_view <- shiny::renderText(logtxt())

    stamp <- function() format(Sys.time(), "%Y%m%d_%H%M")
    output$script <- shiny::downloadHandler(
      function() paste0("multilevelcoda_analysis_", stamp(), ".R"),
      function(f) writeLines(script(), f, useBytes = TRUE))
    output$log <- shiny::downloadHandler(
      function() paste0("pipeline_", stamp(), ".txt"),
      function(f) writeLines(logtxt(), f, useBytes = TRUE))
    output$ana_csv <- shiny::downloadHandler(
      function() paste0("analysis_data_", stamp(), ".csv"),
      function(f) utils::write.csv(comp$ana(), f, row.names = FALSE))
    output$ana_xlsx <- shiny::downloadHandler(
      function() paste0("analysis_data_", stamp(), ".xlsx"),
      function(f) writexl::write_xlsx(list(analysis = comp$ana(),
                                           nights = as.data.frame(comp$diary())), f))
    output$diary_csv <- shiny::downloadHandler(
      function() paste0("nights_", stamp(), ".csv"),
      function(f) utils::write.csv(comp$diary(), f, row.names = FALSE))
    output$bundle <- shiny::downloadHandler(
      function() paste0("multilevelcoda_results_", stamp(), ".rds"),
      function(f) {
        fit <- model$fit()
        s <- sub()
        saveRDS(list(
          raw = data_bus$data(), nights = comp$diary(), analysis_data = comp$ana(),
          parts = comp$parts(), total = comp$total(),
          complr = if (!is.null(fit)) fit$complr,
          brmcoda = fit,
          substitution = s,
          substitution_summary = if (!is.null(s)) as.data.frame(summary(s, digits = "asis")),
          pipeline = logtxt(), script = script(),
          versions = pkg_versions()), f)
      })
  })
}
