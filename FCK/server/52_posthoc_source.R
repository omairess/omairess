# ==============================================================================
# server/52_posthoc_source.R — the post-hoc "what am I comparing?" control
#
# The resolution logic is in server/06_helpers_posthoc.R; this is the UI that
# drives it and the note that reports what it resolved to.
# ==============================================================================

# --- the control on the post-hoc tab -----------------------------------------
output$pairwise_source_ui <- renderUI({
  req(values$fanova_results)
  res <- values$fanova_results
  vars <- values$selected_group_vars
  omnibus_var <- res$group_var %||% "the fANOVA variable"
  omnibus_design <- if (identical(res$design, "within")) "within-subjects" else "between-subjects"

  tagList(
    h4("What is being compared"),
    radioButtons("posthoc_source", NULL,
      choices = c(
        setNames("fanova", sprintf("Follow the functional ANOVA — '%s', %s",
                                   omnibus_var, omnibus_design)),
        "Choose a variable and design here" = "custom"),
      selected = "fanova"),
    conditionalPanel(
      condition = "input.posthoc_source == 'custom'",
      div(style = "border-left: 3px solid #fab219; padding-left: 12px; margin-bottom: 10px;",
        helpText(HTML(
          "Comparisons on a variable the omnibus test was not run on are a
           <b>new family of tests</b>, not post-hoc ones. The correction below
           covers this family only.")),
        radioButtons("posthoc_design", "Design:",
          choices = c("Between-subjects (independent)" = "between",
                      "Within-subjects (paired)"       = "within"),
          selected = if (identical(res$design, "within")) "within" else "between",
          inline = TRUE),
        conditionalPanel(
          condition = "input.posthoc_design == 'between'",
          selectInput("posthoc_group_var", "Grouping variable:",
                      choices = vars, selected = vars[1])
        ),
        conditionalPanel(
          condition = "input.posthoc_design == 'within'",
          selectInput("posthoc_subject_var", "Subject ID variable:",
                      choices = vars, selected = vars[1]),
          selectInput("posthoc_rm_var", "Repeated-measures factor:",
                      choices = vars,
                      selected = if (length(vars) > 1) vars[2] else vars[1]),
          helpText("Each subject must appear at more than one level of the factor",
                   "for a paired comparison to mean anything.")
        )
      )
    ),
    verbatimTextOutput("posthoc_source_note"),
    hr()
  )
})

output$posthoc_source_note <- renderText({
  req(values$fanova_results)
  spec <- fck_posthoc_spec(input, values)
  if (!isTRUE(spec$ok)) return(paste("Cannot run yet:", spec$message))
  paste0(spec$description,
         if (!isTRUE(spec$matches_omnibus))
           "\n\nNOTE: this does not match the omnibus test. Report it as a separate family of comparisons."
         else "")
})
