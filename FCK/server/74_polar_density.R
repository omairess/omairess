# ==============================================================================
# server/74_polar_density.R — circular density of acrophases, day at the top
#
# The existing polar plot shows one marker per subject at (amplitude, acrophase).
# That answers "where is each subject?" but not "where does the sample sit?" —
# with 85 rows the markers overplot and the shape of the distribution is lost.
# This is the same circle with a von Mises kernel density on it, oriented as a
# clock face: noon at the top, midnight at the bottom, so the upper half is
# daytime and the lower half is night.
#
# Maths in server/07_helpers_circular.R; this file is the plot and its controls.
# ==============================================================================

# Categorical slots 1-3 of the reference palette, in fixed order, never cycled.
# Overlapping density rings are an all-pairs comparison, and only the first
# three slots clear the all-pairs gates (worst CVD dE 9.2, normal-vision 24.0).
# Past three, line dash carries identity as well as hue — see the note the plot
# shows for itself.
FCK_DENSITY_COLORS <- c("#2a78d6", "#eb6834", "#1baf7a",
                        "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
FCK_DENSITY_DASH   <- c("solid", "dash", "dot", "dashdot",
                        "longdash", "longdashdot", "solid", "dash")

# The acrophases to plot, plus the weights and grouping that go with them.
fck_density_data <- function(input, values) {
  mod <- values$harmonic_model
  if (is.null(mod) || is.null(mod$individual_params)) return(NULL)
  params <- mod$individual_params
  period <- if (!is.null(mod$period)) mod$period else 24

  h <- if (!is.null(input$density_harmonic)) as.integer(input$density_harmonic) else 1L
  h <- max(1L, min(h, mod$n_harmonics %||% 1L))

  acro_col <- paste0("acrophase_time_", h)
  amp_col  <- paste0("amplitude_", h)
  if (!(acro_col %in% names(params))) return(NULL)

  hours <- as.numeric(params[[acro_col]])
  amps  <- if (amp_col %in% names(params)) as.numeric(params[[amp_col]]) else rep(1, length(hours))

  grp <- NULL
  gv <- input$harmonic_group_var
  if (!is.null(gv) && nzchar(gv) && !identical(gv, "_none_") &&
      !is.null(values$covariates) && gv %in% names(values$covariates) &&
      nrow(values$covariates) == length(hours)) {
    grp <- droplevels(as.factor(values$covariates[[gv]]))
  }

  keep <- is.finite(hours)
  list(hours = hours[keep], amps = amps[keep],
       group = if (is.null(grp)) NULL else grp[keep],
       period = period, harmonic = h, group_var = gv)
}

output$density_controls_ui <- renderUI({
  req(values$harmonic_model)
  mod <- values$harmonic_model
  nh <- mod$n_harmonics %||% 1L
  tagList(
    if (nh > 1) selectInput("density_harmonic", "Harmonic:",
                            choices = stats::setNames(seq_len(nh),
                                                      paste0("H", seq_len(nh))),
                            selected = 1),
    radioButtons("density_bw_mode", "Kernel bandwidth:",
                 choices = c("Automatic (Taylor plug-in rule)" = "auto",
                             "Set it myself" = "manual"),
                 selected = "auto"),
    conditionalPanel(
      condition = "input.density_bw_mode == 'manual'",
      sliderInput("density_bw", "Bandwidth (hours):",
                  min = 0.25, max = 6, value = 1.5, step = 0.25)
    ),
    checkboxInput("density_weight_amp", "Weight subjects by amplitude", FALSE),
    helpText("Unweighted, a barely-detectable rhythm counts as much as a strong",
             "one, even though its acrophase is mostly noise."),
    hr(),
    checkboxInput("density_show_points", "Show individual acrophases", TRUE),
    checkboxInput("density_show_mean", "Show mean direction (length = R-bar)", TRUE),
    checkboxInput("density_show_night", "Shade the night sector", TRUE),
    conditionalPanel(
      condition = "input.density_show_night == true",
      fluidRow(
        column(6, numericInput("density_dusk", "Night from (h):", 18, min = 0, max = 24, step = 0.5)),
        column(6, numericInput("density_dawn", "Night to (h):",    6, min = 0, max = 24, step = 0.5))
      )
    )
  )
})

output$harmonic_density_plot <- renderPlotly({
  req(values$harmonic_model)
  dd <- fck_density_data(input, values)
  validate(need(!is.null(dd) && length(dd$hours) >= 2,
                "Run the harmonic regression first: at least two fitted acrophases are needed."))

  period <- dd$period
  weights <- if (isTRUE(input$density_weight_amp)) dd$amps else NULL

  groups <- if (is.null(dd$group)) list(all = seq_along(dd$hours)) else split(seq_along(dd$hours), dd$group)
  groups <- groups[vapply(groups, length, 1L) >= 2]
  validate(need(length(groups) > 0, "No group has two or more fitted acrophases."))

  # One bandwidth for every group, so the curves are comparable. Chosen from the
  # pooled sample when automatic — a per-group bandwidth would make a diffuse
  # group look as peaked as a tight one.
  bw <- if (identical(input$density_bw_mode %||% "auto", "manual") && !is.null(input$density_bw)) {
    input$density_bw
  } else {
    fck_default_bandwidth(dd$hours, period, weights)
  }

  dens <- lapply(groups, function(idx)
    fck_circular_density(dd$hours[idx], bw, period, weights = weights[idx]))
  max_d <- max(vapply(dens, function(d) if (is.null(d)) 0 else max(d$density), 0), na.rm = TRUE)
  validate(need(is.finite(max_d) && max_d > 0, "Density could not be computed."))

  p <- plot_ly()

  # --- night sector, behind everything ---------------------------------------
  if (isTRUE(input$density_show_night)) {
    for (arc in fck_night_arcs(input$density_dusk %||% 18,
                               input$density_dawn %||% 6, period)) {
      hrs <- seq(arc[1], arc[2], length.out = 60)
      p <- add_trace(p, type = "scatterpolar", mode = "lines",
                     r = rep(1.12, length(hrs)),
                     theta = fck_hour_to_theta(hrs, period),
                     fill = "toself", fillcolor = "rgba(11,11,11,0.07)",
                     line = list(color = "rgba(0,0,0,0)"),
                     hoverinfo = "skip", showlegend = FALSE)
    }
  }

  # --- one density ring per group -------------------------------------------
  nm <- names(groups)
  for (i in seq_along(groups)) {
    d <- dens[[i]]; if (is.null(d)) next
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    p <- add_trace(p, type = "scatterpolar", mode = "lines",
                   r = d$density / max_d,
                   theta = fck_hour_to_theta(d$hours, period),
                   name = if (is.null(dd$group)) "Density" else nm[i],
                   fill = if (length(groups) == 1) "toself" else "none",
                   fillcolor = if (length(groups) == 1) "rgba(42,120,214,0.18)" else NULL,
                   line = list(color = col, width = 2,
                               dash = FCK_DENSITY_DASH[((i - 1) %% length(FCK_DENSITY_DASH)) + 1]),
                   hoverinfo = "text",
                   text = sprintf("%s%02d:%02.0f<br>density %.3f",
                                  if (is.null(dd$group)) "" else paste0(nm[i], "<br>"),
                                  floor(d$hours %% period), (d$hours %% 1) * 60,
                                  d$density))
  }

  # --- individual acrophases, on the rim ------------------------------------
  if (isTRUE(input$density_show_points)) {
    for (i in seq_along(groups)) {
      idx <- groups[[i]]
      col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
      p <- add_trace(p, type = "scatterpolar", mode = "markers",
                     r = rep(1.06, length(idx)),
                     theta = fck_hour_to_theta(dd$hours[idx], period),
                     name = paste(nm[i], "(subjects)"), showlegend = FALSE,
                     marker = list(color = col, size = 8, opacity = 0.75,
                                   symbol = "line-ns-open",
                                   line = list(color = col, width = 2)),
                     hoverinfo = "text",
                     text = sprintf("%02d:%02.0f%s",
                                    floor(dd$hours[idx] %% period),
                                    (dd$hours[idx] %% 1) * 60,
                                    if (isTRUE(input$density_weight_amp))
                                      sprintf("<br>amplitude %.2f", dd$amps[idx]) else ""))
    }
  }

  # --- mean direction, length = R-bar ---------------------------------------
  if (isTRUE(input$density_show_mean)) {
    for (i in seq_along(groups)) {
      idx <- groups[[i]]
      s <- fck_circular_summary(dd$hours[idx], period, weights[idx])
      if (is.null(s)) next
      col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
      p <- add_trace(p, type = "scatterpolar", mode = "lines+markers",
                     r = c(0, s$r_bar), theta = rep(fck_hour_to_theta(s$mean_hour, period), 2),
                     name = paste(nm[i], "mean"), showlegend = FALSE,
                     line = list(color = col, width = 3),
                     marker = list(color = col, size = c(1, 9)),
                     hoverinfo = "text",
                     text = sprintf("%smean %02d:%02.0f<br>R-bar %.2f (n = %d)",
                                    if (is.null(dd$group)) "" else paste0(nm[i], " "),
                                    floor(s$mean_hour), (s$mean_hour %% 1) * 60,
                                    s$r_bar, s$n))
    }
  }

  tick_h <- seq(0, period - 1, by = if (period == 24) 2 else max(1, round(period / 12)))
  plotly::layout(p,
    polar = list(
      bgcolor = "#fcfcfb",
      radialaxis = list(range = c(0, 1.16), showticklabels = FALSE,
                        ticks = "", gridcolor = "rgba(11,11,11,0.10)"),
      angularaxis = list(
        tickmode = "array",
        tickvals = fck_hour_to_theta(tick_h, period),
        ticktext = if (period == 24) sprintf("%02d:00", tick_h) else sprintf("%.1f", tick_h),
        direction = "counterclockwise", rotation = 0,
        gridcolor = "rgba(11,11,11,0.10)")),
    showlegend = length(groups) > 1,
    legend = list(orientation = "h", y = -0.08),
    margin = list(t = 40, b = 60))
})

output$harmonic_density_note <- renderText({
  req(values$harmonic_model)
  dd <- fck_density_data(input, values)
  if (is.null(dd) || length(dd$hours) < 2) return("")
  weights <- if (isTRUE(input$density_weight_amp)) dd$amps else NULL
  bw <- if (identical(input$density_bw_mode %||% "auto", "manual") && !is.null(input$density_bw))
    input$density_bw else fck_default_bandwidth(dd$hours, dd$period, weights)
  s <- fck_circular_summary(dd$hours, dd$period, weights)
  n_groups <- if (is.null(dd$group)) 1L else nlevels(dd$group)

  paste0(
    sprintf("Von Mises kernel, bandwidth %.2f h (%s); n = %d.\n",
            bw, if (identical(input$density_bw_mode %||% "auto", "manual"))
              "set by hand" else "Taylor plug-in rule", s$n),
    sprintf("Circular mean %02d:%02.0f, R-bar %.2f — %s.\n",
            floor(s$mean_hour), (s$mean_hour %% 1) * 60, s$r_bar,
            if (s$r_bar < 0.1) "essentially uniform, the mean direction means nothing"
            else if (s$r_bar < 0.45) "dispersed; treat the mean direction cautiously"
            else if (s$r_bar < 0.7) "moderately concentrated" else "tightly concentrated"),
    if (isTRUE(input$density_weight_amp))
      "Weighted by amplitude: stronger rhythms pull the curve more.\n"
    else "Unweighted: every subject counts the same regardless of rhythm strength.\n",
    if (!isTRUE(all.equal(dd$period, 24)))
      sprintf("NOTE: the period is %g h, not 24 — the angle is phase within the period, and the day/night halves do not correspond to clock time.\n", dd$period)
    else "",
    if (n_groups > 3)
      sprintf("NOTE: %d groups. Only the first three colours are separable for every colour-vision type when curves overlap, so line style also distinguishes them; consider comparing fewer groups at a time.\n", n_groups)
    else "")
})
