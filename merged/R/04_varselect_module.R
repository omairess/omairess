# ============================================================================
# 04_varselect_module.R — shared variable selection (house rule 3)
#
# DECISION: standardized on PsychoNetrix's pickerInput pattern (the only one
#   of the three source apps compliant with rule 3); BootSON's dual-listbox
#   Add/Remove UI and Dagger's sortable drag-and-drop are dropped.
# DECISION: each analysis tab gets its OWN instance of this module (its own
#   node set), but all instances read the same data_bus$wide() — tabs may
#   analyse different variable subsets of the same canonical frame.
#   The selection is recorded per-tab (varselect phase) AND baked as an
#   explicit c("V1","V2",...) vector into that tab's analysis fragment, so a
#   later selection change cannot silently redefine an already-run analysis.
# ============================================================================

varselectUI <- function(id, label = "Variables to analyse") {
  ns <- shiny::NS(id)
  # TODO(stage3-picker): replace with
  #   shinyWidgets::pickerInput(ns("vars"), label, choices = NULL,
  #     multiple = TRUE,
  #     options = shinyWidgets::pickerOptions(actionsBox = TRUE,
  #                                           liveSearch = TRUE))
  # plus an optional group-variable selectizeInput(ns("group_var")) with an
  # empty "(none)" default, used by tabs that support split/multi-group runs.
  shiny::tagList(
    shiny::selectInput(shiny::NS(id)("vars"), label,
                       choices = NULL, multiple = TRUE),
    shiny::selectInput(ns("group_var"), "Grouping variable (optional)",
                       choices = c("(none)" = ""))
  )
}

varselectServer <- function(id, data_bus, rec, rec_prefix) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(data_bus$wide(), {
      vars <- names(data_bus$wide())
      shiny::updateSelectInput(session, "vars", choices = vars)
      shiny::updateSelectInput(session, "group_var",
                               choices = c("(none)" = "", vars))
    })

    selected <- shiny::reactive({
      shiny::req(input$vars)
      shiny::validate(shiny::need(length(input$vars) >= 2,
                                  "Select at least 2 variables."))
      rec_upsert(
        rec, paste0(rec_prefix, "_varselect"), "varselect",
        description = sprintf("[%s] Selected %d variables: %s.",
                              rec_prefix, length(input$vars),
                              paste(input$vars, collapse = ", ")),
        code = sprintf("%s_vars <- c(%s)", rec_prefix,
                       paste0('"', input$vars, '"', collapse = ", "))
      )
      input$vars
    })

    group_var <- shiny::reactive({
      if (is.null(input$group_var) || input$group_var == "") NULL
      else input$group_var
    })

    list(vars = selected, group_var = group_var)
  })
}

# Helper every tab uses to bake its subset explicitly into a code fragment.
vars_literal <- function(vars) {
  paste0("c(", paste0('"', vars, '"', collapse = ", "), ")")
}
