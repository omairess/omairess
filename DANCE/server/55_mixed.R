# ==============================================================================
# server/55_mixed.R — the shared views for mixed-design results
# ==============================================================================
# AUDIT (P19). This began as a standalone tab. The user's point was right: a
# separate tab meant a bare-bones one, because every option the Functional ANOVA
# and Cosinor tabs already had would have to be built again. The analysis now
# lives inside those two modules -- a third "Design type" in fANOVA and a third
# "Approach" in the Cosinor tab -- and inherits their permutation count, alpha,
# correction, period, harmonics and group selection.
#
# What is left here is the pair of outputs both hosts render, plus the factor
# picker they share. The statistics are in server/06_helpers_mixed.R (model) and
# server/07_helpers_mixed_perm.R (exact permutation).
# ==============================================================================

dance_mixed_factor_choices <- function(values) {
  out <- character(0)
  if (!is.null(values$covariates) && ncol(values$covariates)) {
    for (nm in names(values$covariates)) {
      v <- values$covariates[[nm]]
      k <- length(unique(v[!is.na(v)]))
      # a usable design factor: at least 2 levels, few enough that a cell can be
      # occupied. A 27-level participant code is not a design factor.
      if (k >= 2 && k <= 12) out <- c(out, nm)
    }
  }
  out
}

output$mixed_results <- renderPrint({
  dance_mixed_readout(values$mixed_results)
})

# --- plot --------------------------------------------------------------------
output$mixed_plot <- renderPlotly({
  res <- values$mixed_results
  if (is.null(res)) return(plot_ly(type = "scatter", mode = "lines") %>%
                             layout(title = "Run a mixed analysis first"))
  if (identical(res$kind, "permutation")) {
    lab <- c(within = "within-subject", between = "between-subject", interaction = "interaction")
    p <- plot_ly(); cols <- dance_group_colors(names(lab))
    for (nm in names(lab)) {
      r <- res[[nm]]; if (is.null(r)) next
      p <- p %>% add_lines(x = res$time, y = r$statistic, name = lab[[nm]],
                           line = list(color = cols[[nm]], width = 2.5))
    }
    return(p %>% layout(
      title = "Pointwise permutation statistic by effect",
      xaxis = list(title = res$time_axis %||% "time"),
      yaxis = list(title = "between-group / between-condition sum of squares")))
  }
  if ((identical(res$kind, "fanova") || identical(res$kind, "model")) && !is.null(res$curves)) {
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
