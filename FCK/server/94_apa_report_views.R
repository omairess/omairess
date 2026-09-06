# ==============================================================================
# server/94_apa_report_views.R — the Shiny wiring for the publication report
#
# The report's LOGIC is in server/93_apa_report.R and is pure: it takes `values`
# and `input` and returns character vectors. It lives in its own file so it can
# be sourced and unit-tested outside a Shiny session (tests/testthat/
# test-apa-report.R does exactly that); an `output$x <- render...` at the top
# level of that file would make it unsourceable.
#
# Nothing here is cached. The report is regenerated on every preview and every
# download, from `values` as they stand, so there is no stored copy that could
# go stale against the screen.
# ==============================================================================
  # The report is regenerated on every download and on every preview, from
  # `values` as they stand. There is no cached copy that could go stale against
  # the screen.
  fck_report_title <- function(input) {
    ttl <- input$apa_report_title
    if (is.null(ttl) || !nzchar(trimws(ttl))) NULL else trimws(ttl)
  }

  output$apa_report_preview <- renderText({
    md <- tryCatch(fck_apa_report(values, input, fck_report_title(input)),
                   error = function(e)
                     paste("The report could not be generated:", conditionMessage(e)))
    paste(md, collapse = "\n")
  })

  output$download_apa_md <- downloadHandler(
    filename = function() paste0("analysis_report_", Sys.Date(), ".md"),
    content  = function(file)
      writeLines(fck_apa_report(values, input, fck_report_title(input)), file,
                 useBytes = TRUE)
  )

  output$download_apa_html <- downloadHandler(
    filename = function() paste0("analysis_report_", Sys.Date(), ".html"),
    content  = function(file) {
      md <- fck_apa_report(values, input, fck_report_title(input))
      writeLines(fck_apa_html(md, fck_report_title(input) %||%
                                 "Functional data analysis: methods and results"),
                 file, useBytes = TRUE)
    }
  )
