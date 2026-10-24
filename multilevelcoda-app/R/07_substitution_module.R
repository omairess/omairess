# ============================================================================
# 07_substitution_module.R — compositional isotemporal substitution
#
# multilevelcoda::substitution(): the expected difference in the outcome when
# `delta` minutes move from one part to another (one-to-one) or from one part
# to all others proportionally (one-to-all), around the grand-mean (or
# cluster-mean) composition, separately for between- and within-person
# reallocations (Le et al., 2025; Dumuid et al., 2019).
# ============================================================================

substitutionModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 4,
      shiny::textInput(ns("delta"), "Minutes to reallocate (comma-separated)", "15, 30, 60"),
      shiny::uiOutput(ns("max_delta")),
      shiny::uiOutput(ns("level_ui")),
      shiny::radioButtons(ns("type"), "Reallocation",
        c("One part to one other part" = "one-to-one",
          "One part to all others (proportionally)" = "one-to-all")),
      shiny::radioButtons(ns("ref"), "Reference composition",
        c("Grand mean (sample compositional mean)" = "grandmean",
          "Each person's own mean, averaged (slower)" = "clustermean")),
      shiny::actionButton(ns("run"), "Run substitution", class = "btn-primary"),
      shiny::hr(),
      partColoursUI(ns),
      shiny::sliderInput(ns("base_size"), "Text size", 8, 24, 12)
    ),
    shiny::mainPanel(
      width = 8,
      shiny::uiOutput(ns("status")),
      shiny::plotOutput(ns("plot"), height = "560px"),
      plotExportUI(ns, "sub", w = 3000, h = 2000),
      shiny::h4("Estimates (posterior mean, SD, 95% credible interval)"),
      DT::DTOutput(ns("table"))
    )
  )
}

substitutionModuleServer <- function(id, model, rec) {
  shiny::moduleServer(id, function(input, output, session) {
    output$level_ui <- shiny::renderUI({
      if (!isTRUE(model$multilevel())) {
        return(shiny::helpText("Single-level model: substitution at the aggregate level."))
      }
      shiny::checkboxGroupInput(session$ns("level"), "Level",
        c("Between-person" = "between", "Within-person" = "within"),
        selected = c("between", "within"))
    })

    deltas <- shiny::reactive({
      d <- suppressWarnings(as.numeric(strsplit(input$delta, "[,; ]+")[[1]]))
      d <- unique(abs(round(d[!is.na(d) & d != 0])))
      shiny::validate(shiny::need(length(d) > 0, "Give at least one number of minutes."))
      as.integer(d)
    })

    level_arg <- shiny::reactive({
      f <- model$fit()
      shiny::req(f)
      v <- multilevelcoda::get_variables(f)
      if (!isTRUE(model$multilevel()) || identical(v$fixef_type, "aggregate")) "aggregate"
      else input$level %||% c("between", "within")
    })

    sub <- shiny::reactiveVal(NULL)
    sub_error <- shiny::reactiveVal(NULL)
    shiny::observeEvent(model$fit(), {
      sub(NULL); sub_error(NULL); rec$clear_from("substitution")
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # Minutes cannot be taken from a part that holds fewer minutes than that
    # at the reference composition: show the compositional mean of each part.
    output$max_delta <- shiny::renderUI({
      f <- model$fit()
      shiny::req(f)
      out <- f$complr$output[[1]]
      centre <- comp_center(as.data.frame(f$complr$datain)[, out$parts, drop = FALSE],
                            out$total)
      shiny::helpText(sprintf(
        "Compositional mean: %s. Moves much larger than the smallest part (about %.1f min, %s) are impossible; the exact limit is reported if exceeded.",
        paste(sprintf("%s %.0f", names(centre), centre), collapse = ", "),
        min(centre), names(centre)[which.min(centre)]))
    })

    output$status <- shiny::renderUI({
      if (is.null(model$fit()))
        return(shiny::div(class = "alert alert-info", "Fit a model on the Model tab first."))
      if (!identical(model$role(), "predictor"))
        return(shiny::div(class = "alert alert-info",
          "Substitution needs the composition as a predictor (Model tab)."))
      if (!is.null(sub_error()))
        return(shiny::div(class = "alert alert-danger",
                          paste("Substitution failed:", sub_error())))
      if (is.null(sub())) return(shiny::helpText("Choose the settings and press Run."))
      NULL
    })

    shiny::observeEvent(input$run, {
      f <- model$fit()
      shiny::req(f, identical(model$role(), "predictor"))
      lv <- level_arg()
      res <- tryCatch(
        shiny::withProgress(message = "Computing substitution effects...", value = 0.3, {
          multilevelcoda::substitution(object = f, delta = deltas(), level = lv,
                                       ref = input$ref, type = input$type)
        }), error = function(e) e)
      if (inherits(res, "error")) {
        sub(NULL)
        sub_error(conditionMessage(res))
        rec$clear_from("substitution")
        return()
      }
      sub_error(NULL)
      sub(res)
      rec$record("substitution",
        sprintf("Substitution analysis: %s reallocations of %s min at the %s level(s), reference = %s",
                input$type, paste(deltas(), collapse = "/"), paste(lv, collapse = " and "),
                input$ref),
        c(sprintf("sub <- substitution(object = fit, delta = %s, level = %s,\n  ref = %s, type = %s)",
                  r_lit(deltas()), r_lit(lv), r_lit(input$ref), r_lit(input$type)),
          "print(summary(sub))"))
    })

    parts <- shiny::reactive({ s <- sub(); shiny::req(s); s$parts })
    render_part_colours(output, session, parts)
    cols_r <- shiny::reactive(read_part_colours(input, parts()))

    sub_plot <- shiny::reactive({
      s <- sub()
      shiny::req(s)
      plot_substitution(s, cols_r(), input$base_size)
    })
    output$plot <- shiny::renderPlot(sub_plot())
    plotExportServer(input, output, "sub", sub_plot, "substitution")

    shiny::observe({
      shiny::req(sub())
      rec$record("sub_plot", "Plotted substitution estimates (pastel colours chosen in the app)",
        c(sprintf("sub_colours <- %s", colours_literal(cols_r())),
          sprintf("print(plot_substitution(sub, sub_colours, base_size = %s))", input$base_size)))
    }, priority = 10)

    output$table <- DT::renderDT({
      s <- sub()
      shiny::req(s)
      tab <- as.data.frame(summary(s, digits = 2))
      DT::datatable(tab, rownames = FALSE, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE))
    })

    sub
  })
}
