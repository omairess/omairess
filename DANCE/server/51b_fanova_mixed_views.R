# ==============================================================================
# server/51b_fanova_mixed_views.R — the mixed pickers and readout for the
# Functional ANOVA tab
# ==============================================================================
# Kept out of server/50_fanova.R so that file stays what it is (the one-way
# kernels and their wiring), and out of the kernel files so those stay pure.
# ==============================================================================

output$fanova_mixed_subject_ui <- renderUI({
  if (is.null(values$subject_ids))
    return(helpText(HTML(
      "<b>No participant identifier.</b> A mixed design needs to know which rows
       belong to the same person; pick an ID column on the Data Import tab.")))
  helpText(HTML(sprintf(
    "Participant identifier: <b>%d</b> participants across <b>%d</b> rows.",
    length(unique(values$subject_ids)), length(values$subject_ids))))
})

output$fanova_mixed_between_ui <- renderUI({
  ch <- dance_mixed_factor_choices(values)
  if (!length(ch)) return(helpText("No categorical covariate with 2-12 levels was found."))
  selectInput("fanova_mixed_between", "Between-subject factor:", choices = ch)
})

output$fanova_mixed_within_ui <- renderUI({
  ch <- dance_mixed_factor_choices(values)
  if (!length(ch)) return(NULL)
  selectInput("fanova_mixed_within", "Within-subject factor (repeated):", choices = ch,
              selected = if (length(ch) > 1) ch[2] else ch[1])
})

# The long frame, built from the fANOVA tab's own pickers. Same shape as the
# Mixed tab's, and deliberately the same helper, so the two cannot drift.
dance_fanova_mixed_frame <- reactive({
  if (is.null(values$data) || is.null(values$subject_ids)) return(NULL)
  bn <- input$fanova_mixed_between; wn <- input$fanova_mixed_within
  if (is.null(bn) || is.null(wn) || !nzchar(bn) || !nzchar(wn)) return(NULL)
  if (is.null(values$covariates) || !all(c(bn, wn) %in% names(values$covariates))) return(NULL)
  if (length(values$subject_ids) != nrow(values$data)) return(NULL)
  axis <- dance_smoothing_axis(
    list(use_real_time = isTRUE(input$fanova_mixed_real_time), is_cyclic = FALSE), values)
  out <- tryCatch(
    dance_mixed_long(values$data, axis$t_full, subject = values$subject_ids,
                     between = values$covariates[[bn]], within = values$covariates[[wn]]),
    error = function(e) NULL)
  if (!is.null(out)) attr(out, "time_values") <- axis$t_full
  out
})
