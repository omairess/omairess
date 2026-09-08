# ==============================================================================
# server/55_mixed.R — wiring for the mixed-design tab
# ==============================================================================
# The statistics live in server/06_helpers_mixed.R, which is pure. This file
# reads the controls, builds the long form, calls a kernel and renders it.
# ==============================================================================

# --- which columns can play which role ---------------------------------------
# The between and within factors come from the scalar covariates; the subject
# identifier is whatever the import step captured (it is what makes the design
# mixed, so it is not optional here).
output$mixed_subject_ui <- renderUI({
  if (is.null(values$subject_ids))
    return(helpText(HTML(
      "<b>No participant identifier.</b> A mixed design needs to know which rows
       belong to the same person. Pick an ID column on the Data Import tab.")))
  helpText(HTML(sprintf(
    "Participant identifier: <b>%d</b> distinct participants across <b>%d</b> rows.",
    length(unique(values$subject_ids)), length(values$subject_ids))))
})

dance_mixed_factor_choices <- function(values) {
  out <- character(0)
  if (!is.null(values$covariates) && ncol(values$covariates)) {
    for (nm in names(values$covariates)) {
      v <- values$covariates[[nm]]
      k <- length(unique(v[!is.na(v)]))
      # a usable factor: at least 2 levels, and few enough that a cell can be
      # occupied. A 27-level "ID" column is not a design factor.
      if (k >= 2 && k <= 12) out <- c(out, nm)
    }
  }
  out
}

output$mixed_between_ui <- renderUI({
  ch <- dance_mixed_factor_choices(values)
  if (!length(ch)) return(helpText("No categorical covariate with 2-12 levels was found."))
  selectInput("mixed_between", "Between-subject factor:", choices = ch)
})

output$mixed_within_ui <- renderUI({
  ch <- dance_mixed_factor_choices(values)
  if (!length(ch)) return(NULL)
  selectInput("mixed_within", "Within-subject factor (repeated):", choices = ch,
              selected = if (length(ch) > 1) ch[2] else ch[1])
})

# --- a live description of the design, before anything is fitted --------------
output$mixed_design_summary <- renderUI({
  d <- dance_mixed_frame()
  if (is.null(d)) return(NULL)
  bad <- dance_mixed_check(d)
  b   <- dance_mixed_balance(d)
  cells <- b$cells
  txt <- sprintf("%d observations, %d participants (%d with every level of the within factor).",
                 b$n_obs, b$n_subjects, b$n_complete)
  if (length(bad))
    return(div(style = "color:#b71c1c",
               HTML(paste0("<b>Not a mixed design as selected.</b><br>",
                           paste(bad, collapse = "<br>")))))
  div(HTML(paste0("<b>Design looks valid.</b><br>", txt, "<br>",
                  paste(sprintf("%s x %s: %d", cells$within, cells$between, cells$Freq),
                        collapse = "<br>"))))
})

# --- the long form, shared by the summary and the fit ------------------------
dance_mixed_frame <- reactive({
  if (is.null(values$data) || is.null(values$subject_ids)) return(NULL)
  bn <- input$mixed_between; wn <- input$mixed_within
  if (is.null(bn) || is.null(wn) || !nzchar(bn) || !nzchar(wn)) return(NULL)
  if (is.null(values$covariates) ||
      !all(c(bn, wn) %in% names(values$covariates))) return(NULL)
  if (length(values$subject_ids) != nrow(values$data)) return(NULL)

  # The time axis. Real elapsed hours when asked for and available, else the
  # column index -- the same rule the smoothing module uses, through the same
  # helper, so the two cannot drift (the lesson of P9.3/P10.2/P11.3).
  axis <- dance_smoothing_axis(
    list(use_real_time = isTRUE(input$mixed_real_time),
         is_cyclic = FALSE), values)
  out <- tryCatch(
    dance_mixed_long(values$data, axis$t_full, subject = values$subject_ids,
                     between = values$covariates[[bn]],
                     within  = values$covariates[[wn]]),
    error = function(e) NULL)
  # P18.2: the exported script has to rebuild this exact axis, so the values
  # travel with the frame rather than being recomputed from a description of
  # them -- the same reason the registration kernels are emitted rather than
  # re-implemented.
  if (!is.null(out)) attr(out, "time_values") <- axis$t_full
  out
})

# --- run ---------------------------------------------------------------------
observeEvent(input$run_mixed, {
  d <- dance_mixed_frame()
  if (is.null(d)) {
    showNotification("Load data, choose a participant identifier on the Import tab, and pick both factors first.",
                     type = "error", duration = 8)
    return()
  }
  bad <- dance_mixed_check(d)
  if (length(bad)) {
    showNotification(paste(bad, collapse = " "), type = "error", duration = 15)
    return()
  }

  withProgress(message = "Fitting the mixed model...", value = 0.3, {
    if (identical(input$mixed_analysis, "cosinor")) {
      res <- dance_mixed_cosinor(d,
                                 period = suppressWarnings(as.numeric(input$mixed_period %||% 24)),
                                 n_harmonics = input$mixed_harmonics %||% 1)
      res$kind <- "cosinor"
    } else {
      res <- dance_mixed_fanova(d,
                               k_time = input$mixed_k_time %||% 12,
                               k_subject = input$mixed_k_subject %||% 6)
      res$kind <- "fanova"
      if (isTRUE(res$ok)) res$curves <- dance_mixed_fanova_curves(res)
    }
  })

  if (!isTRUE(res$ok)) {
    showNotification(res$message %||% "The mixed model could not be fitted.",
                     type = "error", duration = 15)
    values$mixed_results <- NULL
    return()
  }
  res$between_name <- input$mixed_between
  res$within_name  <- input$mixed_within
  res$time_values  <- attr(d, "time_values")
  res$time_axis    <- if (isTRUE(input$mixed_real_time)) "real elapsed time" else "column index"
  values$mixed_results <- res
  showNotification("Mixed model fitted.", type = "message", duration = 4)
})

# --- readout -----------------------------------------------------------------
output$mixed_results <- renderPrint({
  res <- values$mixed_results
  if (is.null(res)) {
    cat("Run a mixed analysis to see results.\n\n")
    cat("This tab is for a design with one BETWEEN-subject factor and one\n")
    cat("WITHIN-subject factor. The Functional ANOVA tab handles one or the\n")
    cat("other, not both, and carries no interaction term.\n")
    return(invisible(NULL))
  }
  f2 <- function(x) if (is.finite(x)) sprintf("%.2f", x) else "--"
  f3 <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "--"
  pf <- function(p) if (!is.finite(p)) "--" else if (p < .001) "< .001" else sub("^0", "", sprintf("%.3f", p))

  b <- res$balance
  cat("Design\n======\n")
  cat(sprintf("  between: %s   within: %s   subjects: %d   observations: %d\n",
              res$between_name, res$within_name, b$n_subjects, res$n_obs))
  cat(sprintf("  time axis: %s\n", res$time_axis))
  if (b$n_partial > 0)
    cat(sprintf("  %d participant(s) do not have every level of the within factor.\n",
                b$n_partial))
  cat("\n")

  if (identical(res$kind, "fanova")) {
    cat("Mixed functional model\n======================\n")
    cat("  ", res$formula_full, "\n\n", sep = "")
    if (!is.null(res$s_table)) {
      cat("Smooth terms (approximate)\n")
      st <- res$s_table
      for (i in seq_len(nrow(st)))
        cat(sprintf("  %-28s edf %7s   F %8s   p %s\n", rownames(st)[i],
                    f2(st[i, "edf"]), f2(st[i, "F"]), pf(st[i, "p-value"])))
      cat("\n")
    }
    if (!is.null(res$p_table)) {
      cat("Parametric terms (mean level)\n")
      pt <- res$p_table
      for (i in seq_len(nrow(pt)))
        cat(sprintf("  %-28s est %8s   SE %7s   t %7s   p %s\n", rownames(pt)[i],
                    f2(pt[i, "Estimate"]), f2(pt[i, "Std. Error"]),
                    f2(pt[i, "t value"]), pf(pt[i, "Pr(>|t|)"])))
      cat("\n")
    }
    cat(sprintf("Deviance explained: %s%%\n", f2(100 * res$dev_expl)))
    if (is.finite(res$aic_delta)) {
      cat(sprintf("Interaction in the SHAPE (%s-based model comparison)\n",
                  res$aic_basis %||% "ML"))
      cat(sprintf("  AIC %s with the interaction, %s without; delta %+.1f\n",
                  f2(res$aic_full), f2(res$aic_additive), res$aic_delta))
      # P18.5: reported as model-comparison evidence, not as a hypothesis test.
      # A delta is a weight of evidence; calling it a decision at 2 units dresses
      # a continuous quantity as a verdict.
      cat(sprintf("  %s\n", if (res$aic_delta > 2)
        "the comparison favours letting each cell have its own temporal shape"
        else if (res$aic_delta < -2)
        "the comparison favours the additive model: no support for cell-specific shapes"
        else "the two models are within 2 AIC: this comparison does not separate them"))
      cat("  This is evidence for one model over another, not a test of a null\n")
      cat("  hypothesis, and no p-value should be quoted from it.\n")
    }
    cat("\nWhat this does not establish\n----------------------------\n")
    cat("  The p-values are APPROXIMATE. The smoothing parameters were estimated\n")
    cat("  from these data and the tests condition on those estimates; this is not\n")
    cat("  the exact permutation guarantee the one-way fANOVA gives. The subject\n")
    cat("  term is a random functional effect, so the cell smooths are population\n")
    cat("  curves and no single participant's curve is claimed. A smooth term's\n")
    cat("  p-value tests whether that cell's curve is flat, NOT whether two cells\n")
    cat("  differ -- the interaction comparison above is what addresses that.\n")

  } else {
    cat("Mixed cosinor\n=============\n")
    cat("  ", res$formula, "\n", sep = "")
    cat(sprintf("  period %s, %d harmonic(s); random rhythm per subject: %s\n\n",
                format(res$period), res$n_harmonics,
                if (isTRUE(res$random_rhythm)) "yes" else
                  "NO -- the full model did not converge, so only a random MESOR was fitted"))
    cl <- res$cells
    cat(sprintf("  %-14s %-10s %8s %10s %10s\n", res$within_name, res$between_name,
                "MESOR", "amplitude", "acrophase"))
    for (i in seq_len(nrow(cl)))
      cat(sprintf("  %-14s %-10s %8s %10s %10s\n", cl$within[i], cl$between[i],
                  f2(cl$mesor[i]), f2(cl$amplitude_1[i]), f2(cl$acrophase_1[i])))
    cat("\n  Acrophase is in time units from the start of the observation window.\n")
    if (!is.null(res$tests)) {
      cat("\nLikelihood-ratio tests on the (cosine, sine) pair\n")
      for (i in seq_len(nrow(res$tests)))
        cat(sprintf("  %-40s chi2(%d) = %s, p %s\n", res$tests$term[i],
                    res$tests$df[i], f2(res$tests$chisq[i]), pf(res$tests$p[i])))
      cat("\n  Each test drops a cosine/sine pair, so it asks whether the rhythm\n")
      cat("  differs in amplitude OR phase -- a 2-df question per harmonic. A test\n")
      cat("  of one factor also drops its higher-order terms, so it is the effect\n")
      cat("  of that factor overall, not conditional on the interaction.\n")
    }
    if (isTRUE(res$singular))
      cat("\n  WARNING: the random-effects fit is SINGULAR. A variance component is\n",
          "  estimated at zero; the rhythm parameters are still usable but the\n",
          "  random structure is not supported by this sample.\n", sep = "")
    cat("\nWhat this does not establish\n----------------------------\n")
    cat("  No standard errors are quoted for amplitude or acrophase. Amplitude is\n")
    cat("  a norm and acrophase an angle, both nonlinear in the coefficients, and\n")
    cat("  a delta-method interval on a phase near the period boundary misleads.\n")
    cat("  The tests above are on the coefficient pair, which is what the model\n")
    cat("  can test exactly. The period was FIXED, not estimated, so everything\n")
    cat("  here is conditional on that choice.\n")
  }
})

# --- plot --------------------------------------------------------------------
output$mixed_plot <- renderPlotly({
  res <- values$mixed_results
  if (is.null(res)) return(plot_ly(type = "scatter", mode = "lines") %>%
                             layout(title = "Run a mixed analysis first"))
  if (identical(res$kind, "fanova") && !is.null(res$curves)) {
    cv <- res$curves
    p <- plot_ly()
    # dance_group_colors() is keyed on the LEVEL NAME and returns a named vector;
    # its own documentation says to look up by name and never by position, so
    # that filtering one cell out cannot repaint the others.
    cols <- dance_group_colors(levels(cv$cell))
    for (cl in levels(cv$cell)) {
      s <- cv[cv$cell == cl, ]
      col <- cols[[cl]]
      p <- p %>%
        add_ribbons(data = s, x = ~t, ymin = ~fit - 1.96 * se, ymax = ~fit + 1.96 * se,
                    name = paste(cl, "95% CI"), line = list(width = 0),
                    fillcolor = paste0(col, "33"), showlegend = FALSE) %>%
        add_lines(data = s, x = ~t, y = ~fit, name = cl,
                  line = list(color = col, width = 3))
    }
    return(p %>% layout(
      title = "Fitted cell curves (population, subject effect excluded)",
      xaxis = list(title = res$time_axis), yaxis = list(title = "value")))
  }
  if (identical(res$kind, "cosinor")) {
    cl <- res$cells
    lab <- paste(cl$within, cl$between, sep = " x ")
    return(plot_ly(x = lab, y = cl$amplitude_1, type = "bar",
                   name = "amplitude",
                   marker = list(color = unname(dance_group_colors(lab)[lab]))) %>%
             layout(title = "First-harmonic amplitude per cell",
                    xaxis = list(title = ""), yaxis = list(title = "amplitude")))
  }
  plot_ly(type = "scatter", mode = "lines") %>% layout(title = "No plot for this result")
})
