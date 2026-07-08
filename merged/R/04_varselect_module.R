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

# show_group = FALSE for tabs that provide their OWN grouping control (the
# NCT tab), so the module doesn't render a second, redundant one.
varselectUI <- function(id, label = "Variables to analyse", show_group = TRUE) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shinyWidgets::pickerInput(
      ns("vars"), label, choices = NULL, multiple = TRUE,
      options = shinyWidgets::pickerOptions(
        actionsBox = TRUE, liveSearch = TRUE,
        title = "Select variables...", selectedTextFormat = "count > 4")
    ),
    # Data transform sits BETWEEN variable selection and the grouping
    # variable in every module (user request): selected -> transformed ->
    # (optionally grouped) -> analysed.
    transformUI(ns("transform")),
    if (show_group) shinyWidgets::pickerInput(
      ns("group_var"), "Grouping variable (optional)",
      choices = c("(none)" = ""), multiple = FALSE,
      options = shinyWidgets::pickerOptions(liveSearch = TRUE)
    )
  )
}

varselectServer <- function(id, data_bus, rec, rec_prefix) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(data_bus$wide(), {
      vars <- names(data_bus$wide())
      shinyWidgets::updatePickerInput(session, "vars", choices = vars)
      shinyWidgets::updatePickerInput(session, "group_var",
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

    transform <- shiny::reactive(input$transform %||% "none")

    list(vars = selected, group_var = group_var, transform = transform)
  })
}

# Helper every tab uses to bake its subset explicitly into a code fragment.
vars_literal <- function(vars) {
  paste0("c(", paste0('"', vars, '"', collapse = ", "), ")")
}
