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
  mod <- values$harmonic_model
  nh <- mod$n_harmonics %||% 1L
  has_clock <- !is.null(values$time_clock) && !is.null(values$smooth_data)
  tagList(
    radioButtons("density_what", "Radius shows:",
      choices = c("Harmonic regression fit" = "fit",
                  "Density of fitted acrophases" = "acrophase",
                  "Signal averaged over the clock" = "profile"),
      selected = if (is.null(mod) && has_clock) "profile" else "fit"),
    conditionalPanel(
      condition = "input.density_what == 'fit'",
      helpText(HTML("The fitted curves from <b>1. Fitted Curves</b>, wrapped onto
                     the clock: same coefficients, same band. Show/hide the band
                     with the CI checkbox on that tab."))
    ),
    conditionalPanel(
      condition = "input.density_what == 'profile'",
      helpText("The smoothed curves themselves, averaged across subjects at each",
               "clock time. Needs clock times parsed at import."),
      checkboxInput("density_avg_days",
                    "Average repeated days together (otherwise one ring per day)", TRUE),
      checkboxInput("density_band", "Shade +/- 1 SD across subjects", TRUE)
    ),
    conditionalPanel(
      condition = "input.density_what == 'acrophase'",
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
             "one, even though its acrophase is mostly noise.")
    ),
    hr(),
    checkboxInput("density_fill", "Fill the shape", TRUE),
    radioButtons("density_radial_from", "Radius measured from:",
                 choices = c("The smallest value plotted" = "range",
                             "Zero" = "zero"),
                 selected = "range"),
    helpText(HTML("From the smallest value, the shape is stretched to fill the
                   circle and small differences are easy to see. <b>From zero the
                   radius is proportional to the value</b>, which is the honest
                   scale but flattens a rhythm whose MESOR is far from zero.")),
    sliderInput("density_inner", "Inner radius (0 = the shape may reach the centre):",
                min = 0, max = 0.6, value = 0.15, step = 0.05),
    helpText(HTML("An inner radius keeps a low-density stretch of the clock from
                   collapsing to a point, which is what turns a ring into petals.
                   <b>It also stops the radius being proportional to the value</b>,
                   so read the shape, not the area. Set it to 0 for a faithful
                   radial scale.")),
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

# The fitted cosinor curves from the "Fitted Curves" tab, on a clock face.
#
# Same coefficients, same CI construction, so this IS that plot in polar
# coordinates rather than a second opinion about the same data. The band is
# the app's own approximation -- the curve rescaled about the MESOR by the
# relative standard error of the first harmonic's amplitude, which is what the
# fit plot draws and calls a 95% CI. It is not a pointwise interval, and the
# note says so.
fck_fit_rings <- function(input, values) {
  mod <- values$harmonic_model
  if (is.null(mod)) return(NULL)
  period <- mod$period %||% 24
  nh     <- mod$n_harmonics %||% 1L
  trend  <- mod$trend_type %||% "none"

  # One full turn of the clock. A cosinor is defined everywhere, so this is
  # exact for the rhythm even when the recording covered less than a period.
  t <- seq(0, period, length.out = 361)

  band <- function(pred, mesor, amp_mean, amp_sd, n) {
    if (!isTRUE(input$harmonic_show_ci)) return(NULL)
    if (!is.finite(amp_mean) || amp_mean == 0 || !is.finite(amp_sd) || !is.finite(n) || n < 2)
      return(NULL)
    amp_se <- amp_sd / sqrt(n)
    up <- 1 + 1.96 * amp_se / amp_mean
    lo <- max(0, 1 - 1.96 * amp_se / amp_mean)
    list(lower = mesor + (pred - mesor) * lo,
         upper = mesor + (pred - mesor) * up)
  }

  out <- list()
  if (!is.null(mod$group_fits) && length(mod$group_fits) >= 1) {
    for (g in names(mod$group_fits)) {
      gf <- mod$group_fits[[g]]
      if (is.null(gf$mean_coefs)) next
      pred <- fck_rhythm_from_coefs(gf$mean_coefs, t, period, nh, trend)
      if (!any(is.finite(pred))) next
      b <- band(pred, gf$mean_mesor, gf$mean_amplitudes[1], gf$sd_amplitudes[1], gf$n)
      out[[g]] <- list(hours = t, pred = pred, lower = b$lower, upper = b$upper,
                       n = gf$n)
    }
  }
  if (!length(out)) {
    coefs <- mod$pop_mean_fit$mean_coefs
    params <- mod$individual_params
    if (is.null(coefs) && !is.null(params)) {
      # Same population cosinor the app forms elsewhere: mean of the
      # coefficients, never a mean of acrophases (those are circular).
      coefs <- c(mean(params$mesor, na.rm = TRUE),
                 rep(NA_real_, switch(as.character(trend), "none" = 0, "linear" = 1,
                                      "log" = 1, "exp_sat" = 2, 0)),
                 unlist(lapply(seq_len(nh), function(h) c(
                   mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                   mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE)))))
    }
    if (is.null(coefs)) return(NULL)
    pred <- fck_rhythm_from_coefs(coefs, t, period, nh, trend)
    if (!any(is.finite(pred))) return(NULL)
    b <- if (!is.null(params)) band(pred, coefs[1],
                                    mean(params$amplitude_1, na.rm = TRUE),
                                    stats::sd(params$amplitude_1, na.rm = TRUE),
                                    sum(!is.na(params$amplitude_1))) else NULL
    out[["Population mean"]] <- list(hours = t, pred = pred,
                                     lower = b$lower, upper = b$upper,
                                     n = if (!is.null(params)) nrow(params) else NA_integer_)
  }
  if (!length(out)) NULL else list(rings = out, period = period, trend = trend, nh = nh)
}

# The signal itself over the clock: mean of the smoothed curves at each clock
# time, per day and per group. Returns NULL when clock times were not parsed.
fck_profile_rings <- function(input, values) {
  if (is.null(values$smooth_data) || is.null(values$time_clock)) return(NULL)
  Y <- values$smooth_data
  cum <- fck_cumulative_hours(values$time_labels)
  if (is.null(cum) || length(cum) != ncol(Y)) return(NULL)
  w <- fck_wrap_to_clock(cum, 24)

  gv <- input$harmonic_group_var
  grp <- if (!is.null(gv) && nzchar(gv) && !identical(gv, "_none_") &&
             !is.null(values$covariates) && gv %in% names(values$covariates) &&
             nrow(values$covariates) == nrow(Y)) {
    droplevels(as.factor(values$covariates[[gv]]))
  } else NULL

  day_sets <- if (isTRUE(input$density_avg_days %||% TRUE)) {
    list(`all days` = seq_len(ncol(Y)))
  } else split(seq_len(ncol(Y)), paste("day", w$day))

  rows <- if (is.null(grp)) list(all = seq_len(nrow(Y))) else split(seq_len(nrow(Y)), grp)

  out <- list()
  for (rn in names(rows)) for (dn in names(day_sets)) {
    cols <- day_sets[[dn]]
    if (length(cols) < 3) next
    sub <- Y[rows[[rn]], cols, drop = FALSE]
    hrs <- w$clock[cols]
    # several columns can share a clock time once days are averaged
    m  <- tapply(seq_along(hrs), hrs, function(k) mean(sub[, k, drop = FALSE], na.rm = TRUE))
    sd_ <- tapply(seq_along(hrs), hrs, function(k) stats::sd(as.vector(sub[, k, drop = FALSE]), na.rm = TRUE))
    h  <- as.numeric(names(m))
    label <- paste(c(if (!is.null(grp)) rn, if (length(day_sets) > 1) dn), collapse = " - ")
    if (!nzchar(label)) label <- "Mean"
    out[[label]] <- list(hours = h, mean = as.numeric(m), sd = as.numeric(sd_))
  }
  if (!length(out)) NULL else out
}

output$harmonic_density_plot <- renderPlotly({
  mode <- input$density_what %||% "fit"
  inner <- input$density_inner %||% 0.15
  do_fill <- isTRUE(input$density_fill %||% TRUE)
  period <- 24

  # radius scaled into [inner, 1] so a low stretch does not collapse to a point
  from_zero <- identical(input$density_radial_from %||% "range", "zero")
  scale_r <- function(v, lo, hi) {
    if (from_zero) lo <- min(0, lo, na.rm = TRUE)
    if (!is.finite(hi) || hi <= lo) return(rep((inner + 1) / 2, length(v)))
    inner + (1 - inner) * pmax(0, pmin(1, (v - lo) / (hi - lo)))
  }

  rings <- list(); pts <- NULL; means <- list(); bands <- list()

  if (identical(mode, "fit")) {
    fr <- fck_fit_rings(input, values)
    validate(need(!is.null(fr),
      "Run the harmonic regression first — this is its fitted curves on a clock face."))
    period <- fr$period
    vals <- unlist(lapply(fr$rings, function(r) c(r$pred, r$lower, r$upper)))
    lo <- min(vals, na.rm = TRUE); hi <- max(vals, na.rm = TRUE)
    for (nm2 in names(fr$rings)) {
      rr <- fr$rings[[nm2]]
      cr <- fck_close_ring(rr$hours, rr$pred)
      if (is.null(cr)) next
      rings[[nm2]] <- list(hours = cr$hours, r = scale_r(cr$values, lo, hi),
                           raw = cr$values, unit = "fitted")
      if (!is.null(rr$lower) && !is.null(rr$upper)) {
        b <- fck_band_ring(rr$hours, rr$lower, rr$upper)
        if (!is.null(b)) bands[[nm2]] <- list(hours = b$hours, r = scale_r(b$values, lo, hi))
      }
    }
    validate(need(length(rings) > 0, "The fitted curves could not be drawn."))

  } else if (identical(mode, "profile")) {
    pr <- fck_profile_rings(input, values)
    validate(need(!is.null(pr),
      "The signal profile needs smoothed curves and clock times parsed from the column names. Apply smoothing, or use the acrophase density."))
    lo <- min(unlist(lapply(pr, function(p) c(p$mean, if (isTRUE(input$density_band)) p$mean - p$sd))), na.rm = TRUE)
    hi <- max(unlist(lapply(pr, function(p) c(p$mean, if (isTRUE(input$density_band)) p$mean + p$sd))), na.rm = TRUE)
    for (nm in names(pr)) {
      p <- pr[[nm]]
      r <- fck_close_ring(p$hours, p$mean)
      if (is.null(r)) next
      rings[[nm]] <- list(hours = r$hours, r = scale_r(r$values, lo, hi),
                          raw = r$values, unit = "value")
      if (isTRUE(input$density_band) && any(is.finite(p$sd))) {
        b <- fck_band_ring(p$hours, p$mean - p$sd, p$mean + p$sd)
        if (!is.null(b)) bands[[nm]] <- list(hours = b$hours, r = scale_r(b$values, lo, hi))
      }
    }
    validate(need(length(rings) > 0, "Not enough clock times to draw a ring."))

  } else {
    dd <- fck_density_data(input, values)
    validate(need(!is.null(dd) && length(dd$hours) >= 2,
      "Run the harmonic regression first: at least two fitted acrophases are needed."))
    period <- dd$period
    weights <- if (isTRUE(input$density_weight_amp)) dd$amps else NULL
    groups <- if (is.null(dd$group)) list(Density = seq_along(dd$hours)) else split(seq_along(dd$hours), dd$group)
    groups <- groups[vapply(groups, length, 1L) >= 2]
    validate(need(length(groups) > 0, "No group has two or more fitted acrophases."))

    bw <- if (identical(input$density_bw_mode %||% "auto", "manual") && !is.null(input$density_bw))
      input$density_bw else fck_default_bandwidth(dd$hours, period, weights)

    dens <- lapply(groups, function(idx)
      fck_circular_density(dd$hours[idx], bw, period, weights = weights[idx]))
    hi <- max(vapply(dens, function(d) if (is.null(d)) 0 else max(d$density), 0), na.rm = TRUE)
    validate(need(is.finite(hi) && hi > 0, "Density could not be computed."))
    for (i in seq_along(groups)) {
      d <- dens[[i]]; if (is.null(d)) next
      rings[[names(groups)[i]]] <- list(hours = d$hours, r = scale_r(d$density, 0, hi),
                                        raw = d$density, unit = "density")
      if (isTRUE(input$density_show_mean)) {
        s <- fck_circular_summary(dd$hours[groups[[i]]], period, weights[groups[[i]]])
        if (!is.null(s)) means[[names(groups)[i]]] <- s
      }
    }
    if (isTRUE(input$density_show_points)) pts <- list(groups = groups, hours = dd$hours, amps = dd$amps)
  }

  p <- plot_ly()

  # --- night sector, behind everything ---------------------------------------
  if (isTRUE(input$density_show_night)) {
    for (arc in fck_night_arcs(input$density_dusk %||% 18,
                               input$density_dawn %||% 6, period)) {
      hrs <- seq(arc[1], arc[2], length.out = 60)
      p <- add_trace(p, type = "scatterpolar", mode = "lines",
                     r = c(rep(1.10, length(hrs)), 0),
                     theta = c(fck_hour_to_theta(hrs, period), fck_hour_to_theta(hrs[1], period)),
                     fill = "toself", fillcolor = "rgba(11,11,11,0.06)",
                     line = list(color = "rgba(0,0,0,0)"),
                     hoverinfo = "skip", showlegend = FALSE)
    }
  }

  nm <- names(rings)
  # --- SD bands (profile mode), under the lines ------------------------------
  for (i in seq_along(bands)) {
    b <- bands[[nm[i]]]; if (is.null(b)) next
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    p <- add_trace(p, type = "scatterpolar", mode = "lines",
                   r = b$r, theta = fck_hour_to_theta(b$hours, period),
                   fill = "toself", fillcolor = paste0(col, "26"),
                   line = list(color = "rgba(0,0,0,0)"),
                   name = paste(nm[i], "+/- 1 SD"), showlegend = FALSE,
                   hoverinfo = "skip")
  }

  # --- the rings themselves ---------------------------------------------------
  for (i in seq_along(rings)) {
    rg <- rings[[i]]
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    p <- add_trace(p, type = "scatterpolar", mode = "lines",
                   r = rg$r, theta = fck_hour_to_theta(rg$hours, period),
                   name = nm[i],
                   fill = if (do_fill) "toself" else "none",
                   fillcolor = if (do_fill)
                     paste0(col, if (length(rings) == 1) "8C" else "40") else NULL,
                   line = list(color = col, width = 2,
                               dash = FCK_DENSITY_DASH[((i - 1) %% length(FCK_DENSITY_DASH)) + 1]),
                   hoverinfo = "text",
                   text = sprintf("%s<br>%02d:%02.0f<br>%s %.3f", nm[i],
                                  floor(rg$hours %% period), (rg$hours %% 1) * 60,
                                  rg$unit, rg$raw))
  }

  # --- individual acrophases on the rim --------------------------------------
  if (!is.null(pts)) {
    for (i in seq_along(pts$groups)) {
      idx <- pts$groups[[i]]
      col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
      p <- add_trace(p, type = "scatterpolar", mode = "markers",
                     r = rep(1.05, length(idx)),
                     theta = fck_hour_to_theta(pts$hours[idx], period),
                     showlegend = FALSE, name = "subjects",
                     marker = list(color = col, size = 9, opacity = 0.8,
                                   symbol = "line-ns-open",
                                   line = list(color = col, width = 2)),
                     hoverinfo = "text",
                     text = sprintf("%02d:%02.0f", floor(pts$hours[idx] %% period),
                                    (pts$hours[idx] %% 1) * 60))
    }
  }

  # --- mean direction, length = R-bar (acrophase mode) -----------------------
  for (i in seq_along(means)) {
    s <- means[[i]]
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    p <- add_trace(p, type = "scatterpolar", mode = "lines+markers",
                   r = c(inner, inner + (1 - inner) * s$r_bar),
                   theta = rep(fck_hour_to_theta(s$mean_hour, period), 2),
                   showlegend = FALSE, name = paste(names(means)[i], "mean"),
                   line = list(color = col, width = 3),
                   marker = list(color = col, size = c(1, 9)),
                   hoverinfo = "text",
                   text = sprintf("mean %02d:%02.0f<br>R-bar %.2f (n = %d)",
                                  floor(s$mean_hour), (s$mean_hour %% 1) * 60,
                                  s$r_bar, s$n))
  }

  tick_h <- seq(0, period - 1, by = if (period == 24) 3 else max(1, round(period / 8)))
  plotly::layout(p,
    polar = list(
      bgcolor = "#fcfcfb",
      radialaxis = list(range = c(0, 1.13), showticklabels = FALSE, ticks = "",
                        showline = FALSE, gridcolor = "rgba(11,11,11,0.08)"),
      angularaxis = list(
        tickmode = "array",
        tickvals = fck_hour_to_theta(tick_h, period),
        ticktext = if (period == 24) sprintf("%02d:00", tick_h) else sprintf("%.1f", tick_h),
        direction = "counterclockwise", rotation = 0,
        gridcolor = "rgba(11,11,11,0.12)")),
    showlegend = length(rings) > 1,
    legend = list(orientation = "h", y = -0.08),
    margin = list(t = 30, b = 60))
})

output$harmonic_density_note <- renderText({
  mode <- input$density_what %||% "fit"
  inner <- input$density_inner %||% 0.15
  from_zero <- identical(input$density_radial_from %||% "range", "zero")
  tail_note <- if (inner > 0)
    sprintf("\nInner radius %.2f%s: the radius is NOT proportional to the value, so compare shape, not area.\n",
            inner, if (from_zero) "" else ", measured from the smallest value plotted")
  else if (!from_zero)
    "\nRadius is measured from the smallest value plotted, not from zero, so it is not proportional to the value.\n"
  else "\nRadius is proportional to the value (from zero, inner radius 0).\n"

  if (identical(mode, "fit")) {
    fr <- fck_fit_rings(input, values)
    if (is.null(fr)) return("")
    peaks <- vapply(fr$rings, function(r) r$hours[which.max(r$pred)], numeric(1))
    troughs <- vapply(fr$rings, function(r) r$hours[which.min(r$pred)], numeric(1))
    return(paste0(
      sprintf("The fitted curves from tab 1, on the clock. %d harmonic%s, period %g h.\n",
              fr$nh, if (fr$nh == 1) "" else "s", fr$period),
      paste(sprintf("  %-22s peak %02d:%02.0f, trough %02d:%02.0f (n = %s)",
                    names(fr$rings), floor(peaks), (peaks %% 1) * 60,
                    floor(troughs), (troughs %% 1) * 60,
                    vapply(fr$rings, function(r) as.character(r$n), "")),
            collapse = "\n"), "\n",
      if (!identical(fr$trend, "none"))
        sprintf("The %s trend is NOT drawn: it is not periodic, so 08:00 on two different days would sit at the same angle with different values and the ring would not close. This is the rhythm; the trend is on tab 1.\n", fr$trend)
      else "",
      "A cosinor drawn in polar coordinates is a limacon: r = MESOR + amplitude*cos(angle - acrophase) traces an OFF-CENTRE RING, widest towards the acrophase. That offset is the rhythm, not a plotting artefact.\n",
      if (isTRUE(input$harmonic_show_ci))
        "The band is the approximation tab 1 draws -- the curve rescaled about the MESOR by the relative standard error of the first harmonic's amplitude. It is not a pointwise confidence interval.\n"
      else "",
      tail_note))
  }

  if (identical(mode, "profile")) {
    pr <- fck_profile_rings(input, values)
    if (is.null(pr)) return("")
    rng <- range(unlist(lapply(pr, function(p) p$mean)), na.rm = TRUE)
    peak <- lapply(pr, function(p) p$hours[which.max(p$mean)])
    return(paste0(
      sprintf("Signal averaged over the clock; %d ring%s, values %.2f to %.2f.\n",
              length(pr), if (length(pr) == 1) "" else "s", rng[1], rng[2]),
      paste(sprintf("  %s peaks at %02d:%02.0f", names(pr),
                    floor(unlist(peak)), (unlist(peak) %% 1) * 60), collapse = "\n"), "\n",
      if (isTRUE(input$density_avg_days %||% TRUE))
        "Repeated days are averaged together at each clock time.\n"
      else "One ring per day of the recording.\n",
      tail_note))
  }

  req(values$harmonic_model)
  dd <- fck_density_data(input, values)
  if (is.null(dd) || length(dd$hours) < 2) return("")
  weights <- if (isTRUE(input$density_weight_amp)) dd$amps else NULL
  bw <- if (identical(input$density_bw_mode %||% "auto", "manual") && !is.null(input$density_bw))
    input$density_bw else fck_default_bandwidth(dd$hours, dd$period, weights)
  s <- fck_circular_summary(dd$hours, dd$period, weights)
  n_groups <- if (is.null(dd$group)) 1L else nlevels(dd$group)

  # A sinusoid drawn in polar coordinates is a limacon -- r = a + b cos(theta)
  # traces an off-centre ring, not a lobed shape -- which surprises people
  # reading a cosinor on a clock for the first time.
  paste0(
    sprintf("Von Mises kernel density of the acrophases, bandwidth %.2f h (%s); n = %d.\n",
            bw, if (identical(input$density_bw_mode %||% "auto", "manual"))
              "set by hand" else "Taylor plug-in rule", s$n),
    sprintf("Circular mean %02d:%02.0f, R-bar %.2f — %s.\n",
            floor(s$mean_hour), (s$mean_hour %% 1) * 60, s$r_bar,
            if (s$r_bar < 0.1) "essentially uniform, the mean direction means nothing"
            else if (s$r_bar < 0.45) "dispersed; treat the mean direction cautiously"
            else if (s$r_bar < 0.7) "moderately concentrated" else "tightly concentrated"),
    {
      dens_all <- fck_circular_density(dd$hours, bw, dd$period, weights = weights)
      ct <- fck_density_contrast(dens_all)
      if (is.finite(ct) && ct < 1.5)
        sprintf("The ring is nearly round (density max/min %.2f): at this bandwidth the acrophases are close to uniform. That is the data, not a broken plot — narrow the bandwidth to look for finer structure.\n", ct)
      else if (is.finite(ct))
        sprintf("Density max/min %.2f.\n", ct)
      else ""
    },
    if (isTRUE(input$density_weight_amp))
      "Weighted by amplitude: stronger rhythms pull the shape more.\n"
    else "Unweighted: every subject counts the same regardless of rhythm strength.\n",
    if (!isTRUE(all.equal(dd$period, 24)))
      sprintf("NOTE: the period is %g h, not 24 — the angle is phase within the period, and the day/night halves do not correspond to clock time.\n", dd$period)
    else "",
    if (n_groups > 3)
      sprintf("NOTE: %d groups. Only the first three colours are separable for every colour-vision type when shapes overlap, so line style also distinguishes them; consider comparing fewer groups at a time.\n", n_groups)
    else "",
    tail_note)
})
