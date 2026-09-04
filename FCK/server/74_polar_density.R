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
      radioButtons("density_fit_scope", "Curve shown:",
        choices = c("The fit over the recording (trend included)" = "recording",
                    "The rhythm only (trend removed)" = "rhythm"),
        selected = "recording"),
      conditionalPanel(
        condition = "input.density_fit_scope == 'recording'",
        radioButtons("density_lap_mode", "When the recording laps the dial:",
          choices = c("One continuous spiral" = "spiral",
                      "One trace per day" = "per_lap"),
          selected = "spiral"),
        helpText(HTML("A recording longer than 24 h passes over each clock hour
                       more than once. Drawn as one path the passes overlap and
                       cannot be told apart, which is why this plot can look
                       nothing like the 2-D fit <i>even though every value is
                       identical</i> \u2014 the note below reports the largest
                       disagreement, computed live."))
      ),
      helpText(HTML("The fitted curves from <b>1. Fitted Curves</b>, wrapped onto
                     the clock: same coefficients, same band. Show/hide the band
                     with the CI checkbox on that tab.<br><br>
                     <b>The fit over the recording</b> is that plot exactly —
                     hover a clock time here and you get the value it shows at
                     the same time. Because a trend is not periodic, it is an
                     open arc covering only the hours recorded, and the gap
                     between its two ends IS the trend.<br>
                     <b>The rhythm only</b> drops the trend and keeps
                     MESOR + harmonics. That closes into a ring, but its values
                     are lower or higher than the fit plot's by exactly the
                     trend at that time. With no trend in the model the two are
                     the same curve and this choice does nothing."))
    ),
    conditionalPanel(
      condition = "input.density_what == 'profile'",
      helpText("The smoothed curves themselves, averaged across subjects at each",
               "clock time. Needs clock times parsed at import."),
      radioButtons("density_profile_mode", "A recording longer than 24 h is shown as:",
        choices = c("A spiral over every time point (nothing averaged)" = "spiral",
                    "One ring per day" = "per_day",
                    "Folded onto one dial, repeated clock times averaged" = "average"),
        selected = "spiral"),
      helpText(HTML("<b>Spiral</b> is the 2-D plot in polar coordinates: every
                     column in time order, nothing averaged, nothing dropped, and
                     the gap between successive passes over the same clock hour
                     <i>is</i> the homeostatic rise.<br>
                     <b>Averaging</b> makes 08:00 on day 1 and 08:00 on day 2 into
                     one number. Under extended wakefulness those two differ by
                     the whole trend, so the average is a level that occurred on
                     neither day. Use it only when an 'average day' is the
                     question.")),
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
    selectInput("density_ci_style", "Confidence band drawn as:",
                choices = c("Dotted lines" = "dotted",
                            "Shaded band" = "shaded",
                            "Both" = "both",
                            "Not shown" = "none"),
                selected = "dotted"),
    selectInput("density_tick_step", "Label every:",
                choices = c("hour" = 1, "2 hours" = 2, "3 hours" = 3, "6 hours" = 6),
                selected = 1),
    checkboxInput("density_radial_labels", "Show the value scale on the radius", TRUE),
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
    radioButtons("density_night_style", "Night shown as:",
                 choices = c("Gradient, deepest at solar midnight" = "gradient",
                             "Flat block" = "block",
                             "Not shown" = "none"),
                 selected = "gradient"),
    conditionalPanel(
      condition = "input.density_night_style != 'none'",
      fluidRow(
        column(6, numericInput("density_dusk", "Night from (h):", 23, min = 0, max = 24, step = 0.5)),
        column(6, numericInput("density_dawn", "Night to (h):",    7, min = 0, max = 24, step = 0.5))
      ),
      helpText("The gradient fades in at dusk, is deepest halfway through the",
               "window and fades out at dawn — night without asserting an edge.")
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

  # WHICH curve goes on the clock face. These are genuinely different curves
  # whenever the model has a trend, and the tab used to draw the second one
  # while labelling its radius "fitted response", so its numbers could not be
  # reconciled with the fit plot's.
  #
  #   "recording" -- the fit itself, trend included, evaluated on the model's
  #                  own absolute time axis (mod$time_vec: 8, 9, ... 30) and
  #                  placed on the dial at t %% 24. Agrees with the Fitted
  #                  Curves tab value-for-value. It is an OPEN arc: it covers
  #                  only the hours recorded, and where a trend makes the two
  #                  ends differ, showing that gap is the point.
  #   "rhythm"    -- MESOR + harmonics, trend removed. The periodic part alone,
  #                  which is the only thing that closes into a ring.
  #
  # With trend_type = "none" the two coincide except that one is closed.
  scope <- input$density_fit_scope %||% "recording"
  if (identical(trend, "none")) scope <- "rhythm"
  with_trend <- identical(scope, "recording")

  tv <- mod$time_vec
  if (with_trend && (is.null(tv) || !any(is.finite(tv)))) {
    with_trend <- FALSE; scope <- "rhythm"
  }

  if (with_trend) {
    t0 <- min(tv, na.rm = TRUE); t1 <- max(tv, na.rm = TRUE)
    if (!is.finite(t0) || !is.finite(t1) || t1 <= t0) {
      with_trend <- FALSE; scope <- "rhythm"
    }
  }

  if (with_trend) {
    # Absolute hours across the recording; the clock angle is t %% period, which
    # runs continuously through midnight, so no seam handling is needed.
    #
    # A recording longer than one period LAPS the dial. Drawn as one path the
    # two passes overlap and cannot be told apart -- which is why this plot can
    # look nothing like the 2-D fit even though the values are identical at
    # every point (asserted in tests/polar_agreement_test.R). lap_mode = "per_lap"
    # splits it into one trace per turn of the clock, each labelled with its day,
    # so the reader can follow the curve round instead of guessing.
    n_pts   <- if (isTRUE((max(tv, na.rm = TRUE) - min(tv, na.rm = TRUE)) > period)) 721 else 361
    t_abs   <- seq(min(tv, na.rm = TRUE), max(tv, na.rm = TRUE), length.out = n_pts)
    t_clock <- t_abs %% period
    spans   <- (max(tv, na.rm = TRUE) - min(tv, na.rm = TRUE)) >= period
  } else {
    # One full turn of the clock. A cosinor is defined everywhere, so this is
    # exact for the rhythm even when the recording covered less than a period.
    t_abs <- seq(0, period, length.out = 361)
    t_clock <- t_abs
    spans <- FALSE
  }
  t <- t_abs

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
      pred <- fck_rhythm_from_coefs(gf$mean_coefs, t, period, nh, trend,
                                    include_trend = with_trend,
                                    t_offset = mod$t_offset %||% 0)
      if (!any(is.finite(pred))) next
      b <- band(pred, gf$mean_mesor, gf$mean_amplitudes[1], gf$sd_amplitudes[1], gf$n)
      out[[g]] <- list(hours = t_clock, pred = pred, lower = b$lower, upper = b$upper,
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
    pred <- fck_rhythm_from_coefs(coefs, t, period, nh, trend,
                                  include_trend = with_trend,
                                  t_offset = mod$t_offset %||% 0)
    if (!any(is.finite(pred))) return(NULL)
    b <- if (!is.null(params)) band(pred, coefs[1],
                                    mean(params$amplitude_1, na.rm = TRUE),
                                    stats::sd(params$amplitude_1, na.rm = TRUE),
                                    sum(!is.na(params$amplitude_1))) else NULL
    out[["Population mean"]] <- list(hours = t_clock, pred = pred,
                                     lower = b$lower, upper = b$upper,
                                     n = if (!is.null(params)) nrow(params) else NA_integer_)
  }
  if (!length(out)) return(NULL)

  # ---- optionally split each arc into one trace per lap of the dial ---------
  lap_mode <- input$density_lap_mode %||% "spiral"
  if (with_trend && spans && identical(lap_mode, "per_lap")) {
    split_out <- list()
    for (nm2 in names(out)) {
      rr <- out[[nm2]]
      lap <- floor(t_abs / period)
      lap <- lap - min(lap) + 1L
      for (L in sort(unique(lap))) {
        k <- which(lap == L)
        if (length(k) < 2) next
        lbl <- if (length(unique(lap)) > 1) sprintf("%s - day %d", nm2, L) else nm2
        split_out[[lbl]] <- list(
          hours = rr$hours[k], pred = rr$pred[k],
          lower = if (!is.null(rr$lower)) rr$lower[k] else NULL,
          upper = if (!is.null(rr$upper)) rr$upper[k] else NULL,
          n = rr$n, lap = L, parent = nm2)
      }
    }
    if (length(split_out)) out <- split_out
  }

  # ---- the agreement check, computed rather than asserted -------------------
  # The tab claims to be the Fitted Curves plot in polar coordinates. That claim
  # is checkable, so it is checked here and the result is printed in the note:
  # the same coefficients evaluated the same way must give the same number.
  agreement <- NA_real_
  if (with_trend && !is.null(mod$group_fits) && length(mod$group_fits)) {
    g1 <- mod$group_fits[[1]]
    if (!is.null(g1$mean_coefs)) {
      chk <- fck_rhythm_from_coefs(g1$mean_coefs, t_abs, period, nh, trend,
                                   include_trend = TRUE, t_offset = mod$t_offset %||% 0)
      ref <- out[[1]]
      base_pred <- if (!is.null(ref$parent)) {
        # rings were split into laps; rebuild the parent's full prediction
        fck_rhythm_from_coefs(mod$group_fits[[ref$parent]]$mean_coefs, t_abs,
                              period, nh, trend, include_trend = TRUE,
                              t_offset = mod$t_offset %||% 0)
      } else ref$pred
      if (length(base_pred) == length(chk)) agreement <- max(abs(base_pred - chk))
    }
  }

  list(rings = out, period = period, trend = trend, nh = nh,
       scope = scope, closed = !with_trend, spans_period = spans,
       lap_mode = lap_mode, agreement = agreement,
       t_range = if (with_trend) range(tv, na.rm = TRUE) else NULL)
}

# The signal itself over the clock: mean of the smoothed curves at each clock
# time, per day and per group. Returns NULL when clock times were not parsed.
fck_profile_rings <- function(input, values) {
  if (is.null(values$smooth_data) || is.null(values$time_clock)) return(NULL)
  Y <- values$smooth_data
  cum <- fck_cumulative_hours(values$time_labels)
  if (is.null(cum) || length(cum) != ncol(Y)) return(NULL)
  # fck_cumulative_hours() returns ELAPSED hours from the first column (0, 1.5,
  # 3, ...), not clock times. Wrapping that onto the dial put column 1 at 00:00
  # whatever time the recording actually started, rotating the whole ring
  # backwards by the start hour -- an 08:00 start moved the 04:00 peak to 20:00.
  # The clock hour of each column is what fck_clock_hours() already parsed; cum
  # is only needed for the day index.
  clock <- fck_clock_hours(values$time_labels)
  if (is.null(clock) || length(clock) != ncol(Y)) return(NULL)
  w <- list(clock = clock, day = fck_wrap_to_clock(cum, 24)$day)

  gv <- input$harmonic_group_var
  grp <- if (!is.null(gv) && nzchar(gv) && !identical(gv, "_none_") &&
             !is.null(values$covariates) && gv %in% names(values$covariates) &&
             nrow(values$covariates) == nrow(Y)) {
    droplevels(as.factor(values$covariates[[gv]]))
  } else NULL

  # ==========================================================================
  # HOW A RECORDING LONGER THAN 24 h IS PUT ON A 24 h DIAL
  #
  # There is no single right answer, and the wrong default silently destroys
  # the thing you are looking at, so this is a choice rather than a convention.
  #
  #   "average"  fold the recording onto one dial and AVERAGE the columns that
  #              share a clock time. 08:00 on day 1 and 08:00 on day 2 become
  #              one number. That is what you want if the question is "what
  #              does an average day look like" -- and it is exactly wrong if
  #              the two days differ, because under extended wakefulness they
  #              differ by the whole homeostatic rise. Averaging then reports a
  #              level that occurs on neither day.
  #   "spiral"   plot EVERY column in time order as one continuous path that
  #              laps the dial. Nothing is averaged and nothing is dropped:
  #              this is the 2-D plot in polar coordinates, and the gap between
  #              successive passes over the same clock hour IS the trend.
  #   "per_day"  one ring per day, so the passes can be compared directly.
  #
  # Default is "spiral" whenever the recording actually spans more than one
  # period, and "average" when it does not (where the three coincide anyway).
  # ==========================================================================
  n_periods <- (max(cum, na.rm = TRUE) - min(cum, na.rm = TRUE)) / 24
  mode <- input$density_profile_mode %||% (if (n_periods > 1) "spiral" else "average")
  # honour the old checkbox if the new control has not been rendered yet
  if (is.null(input$density_profile_mode) && !is.null(input$density_avg_days))
    mode <- if (isTRUE(input$density_avg_days)) "average" else "per_day"

  rows <- if (is.null(grp)) list(all = seq_len(nrow(Y))) else split(seq_len(nrow(Y)), grp)

  out <- list()

  if (identical(mode, "spiral")) {
    # Column order IS time order; keep it, average across SUBJECTS only.
    ord <- order(cum)
    for (rn in names(rows)) {
      sub <- Y[rows[[rn]], ord, drop = FALSE]
      m   <- colMeans(sub, na.rm = TRUE)
      sd_ <- apply(sub, 2, stats::sd, na.rm = TRUE)
      label <- if (is.null(grp)) "Mean" else rn
      out[[label]] <- list(hours = w$clock[ord], mean = as.numeric(m),
                           sd = as.numeric(sd_), spiral = TRUE,
                           day = w$day[ord], cum = cum[ord])
    }
    return(if (!length(out)) NULL else out)
  }

  day_sets <- if (identical(mode, "average")) {
    list(`all days` = seq_len(ncol(Y)))
  } else split(seq_len(ncol(Y)), paste("day", w$day))

  for (rn in names(rows)) for (dn in names(day_sets)) {
    cols <- day_sets[[dn]]
    if (length(cols) < 3) next
    sub <- Y[rows[[rn]], cols, drop = FALSE]
    hrs <- w$clock[cols]
    # several columns can share a clock time once days are averaged
    m  <- tapply(seq_along(hrs), hrs, function(k) mean(sub[, k, drop = FALSE], na.rm = TRUE))
    sd_ <- tapply(seq_along(hrs), hrs, function(k) stats::sd(as.vector(sub[, k, drop = FALSE]), na.rm = TRUE))
    n_folded <- tapply(seq_along(hrs), hrs, length)
    h  <- as.numeric(names(m))
    label <- paste(c(if (!is.null(grp)) rn, if (length(day_sets) > 1) dn), collapse = " - ")
    if (!nzchar(label)) label <- "Mean"
    out[[label]] <- list(hours = h, mean = as.numeric(m), sd = as.numeric(sd_),
                         spiral = FALSE, n_folded = as.integer(n_folded))
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
  # The baseline is decided once per mode and kept in axis_lo/axis_hi, so the
  # radial tick labels below are guaranteed to use the same mapping as the data.
  axis_lo <- NA_real_; axis_hi <- NA_real_; axis_unit <- ""
  scale_r <- function(v, lo, hi) {
    if (!is.finite(hi) || hi <= lo) return(rep((inner + 1) / 2, length(v)))
    inner + (1 - inner) * pmax(0, pmin(1, (v - lo) / (hi - lo)))
  }
  set_axis <- function(lo, hi, unit) {
    if (from_zero) lo <- min(0, lo, na.rm = TRUE)
    axis_lo <<- lo; axis_hi <<- hi; axis_unit <<- unit
    invisible(NULL)
  }

  rings <- list(); pts <- NULL; means <- list(); bands <- list()

  if (identical(mode, "fit")) {
    fr <- fck_fit_rings(input, values)
    validate(need(!is.null(fr),
      "Run the harmonic regression first — this is its fitted curves on a clock face."))
    period <- fr$period
    vals <- unlist(lapply(fr$rings, function(r) c(r$pred, r$lower, r$upper)))
    lo <- min(vals, na.rm = TRUE); hi <- max(vals, na.rm = TRUE)
    set_axis(lo, hi, if (isTRUE(fr$closed) && !identical(fr$trend, "none"))
                       "fitted response, trend removed" else "fitted response")
    lo <- axis_lo
    fit_closed <- isTRUE(fr$closed)
    fit_scope  <- fr$scope
    fit_trend  <- fr$trend
    fit_spans  <- isTRUE(fr$spans_period)
    # A closed ring is sorted by clock hour and its first point repeated; an
    # open arc must keep the order it was generated in, because it runs through
    # midnight and sorting would draw it backwards through itself.
    shape <- if (fit_closed) fck_close_ring else fck_open_path
    ribbon <- if (fit_closed) fck_band_ring else fck_band_path
    for (nm2 in names(fr$rings)) {
      rr <- fr$rings[[nm2]]
      cr <- shape(rr$hours, rr$pred)
      if (is.null(cr)) next
      rings[[nm2]] <- list(hours = cr$hours, r = scale_r(cr$values, lo, hi),
                           raw = cr$values, unit = "fitted", closed = fit_closed)
      if (!is.null(rr$lower) && !is.null(rr$upper)) {
        b <- ribbon(rr$hours, rr$lower, rr$upper)
        eu <- shape(rr$hours, rr$upper)
        el <- shape(rr$hours, rr$lower)
        bands[[nm2]] <- list(
          hours = if (!is.null(b)) b$hours else NULL,
          r     = if (!is.null(b)) scale_r(b$values, lo, hi) else NULL,
          edges = Filter(Negate(is.null), list(
            if (!is.null(eu)) list(hours = eu$hours, r = scale_r(eu$values, lo, hi), raw = eu$values),
            if (!is.null(el)) list(hours = el$hours, r = scale_r(el$values, lo, hi), raw = el$values))))
      }
    }
    validate(need(length(rings) > 0, "The fitted curves could not be drawn."))

  } else if (identical(mode, "profile")) {
    pr <- fck_profile_rings(input, values)
    validate(need(!is.null(pr),
      "The signal profile needs smoothed curves and clock times parsed from the column names. Apply smoothing, or use the acrophase density."))
    lo <- min(unlist(lapply(pr, function(p) c(p$mean, if (isTRUE(input$density_band)) p$mean - p$sd))), na.rm = TRUE)
    hi <- max(unlist(lapply(pr, function(p) c(p$mean, if (isTRUE(input$density_band)) p$mean + p$sd))), na.rm = TRUE)
    set_axis(lo, hi, "response"); lo <- axis_lo
    for (nm in names(pr)) {
      p <- pr[[nm]]
      # A spiral keeps time order and must NOT be sorted or closed: sorting it
      # would fold day 2 back onto day 1, which is the very averaging the mode
      # exists to avoid, and closing it would join the last observation to the
      # first across a gap that was never measured.
      is_spiral <- isTRUE(p$spiral)
      shape  <- if (is_spiral) fck_open_path else fck_close_ring
      ribbon <- if (is_spiral) fck_band_path else fck_band_ring
      r <- shape(p$hours, p$mean)
      if (is.null(r)) next
      rings[[nm]] <- list(hours = r$hours, r = scale_r(r$values, lo, hi),
                          raw = r$values, unit = "value", closed = !is_spiral)
      if (isTRUE(input$density_band) && any(is.finite(p$sd))) {
        b  <- ribbon(p$hours, p$mean - p$sd, p$mean + p$sd)
        eu <- shape(p$hours, p$mean + p$sd)
        el <- shape(p$hours, p$mean - p$sd)
        bands[[nm]] <- list(
          hours = if (!is.null(b)) b$hours else NULL,
          r     = if (!is.null(b)) scale_r(b$values, lo, hi) else NULL,
          edges = Filter(Negate(is.null), list(
            if (!is.null(eu)) list(hours = eu$hours, r = scale_r(eu$values, lo, hi), raw = eu$values),
            if (!is.null(el)) list(hours = el$hours, r = scale_r(el$values, lo, hi), raw = el$values))))
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
    set_axis(0, hi, "density")
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

  # --- night, behind everything ----------------------------------------------
  night_style <- input$density_night_style %||% "gradient"
  dusk <- input$density_dusk %||% 23
  dawn <- input$density_dawn %||% 7

  if (identical(night_style, "gradient")) {
    gr <- fck_night_gradient(dusk, dawn, period, n = 180)
    if (!is.null(gr)) {
      # One barpolar trace with a colour per wedge: plotly has no angular
      # gradient, and a stack of filled polygons would seam visibly.
      # Drawn as an ANNULUS, not wedges from the centre. Wedges converge to a
      # point, so the colour is most intense at the origin -- the one place on
      # the plot that carries no time at all -- and it tints the data fills
      # worst where they overlap. Starting at the inner radius keeps the
      # intensity even along the arc and the hub clean.
      night_base <- max(inner, 0.12)
      # A background that competes with the data is a failed background. The
      # wedges are drawn at partial opacity so the curves stay the most salient
      # thing on the plot; the gradient still reads, it just stops shouting.
      p <- add_trace(p, type = "barpolar",
                     r = rep(1.10 - night_base, nrow(gr)),
                     base = night_base,
                     theta = fck_hour_to_theta(gr$hour, period),
                     width = gr$width * (360 / period),
                     marker = list(color = gr$color, line = list(width = 0),
                                   opacity = 0.55),
                     hoverinfo = "skip", showlegend = FALSE)
    }
  } else if (identical(night_style, "block")) {
    for (arc in fck_night_arcs(dusk, dawn, period)) {
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
  ci_style <- input$density_ci_style %||% "dotted"

  # --- uncertainty: shaded ribbon and/or dotted edges ------------------------
  # A ribbon is easy to read but hides the ring where two groups overlap; dotted
  # edges stay legible over a fill and over each other, which is why they are
  # the default here and the ribbon is not.
  for (i in seq_along(bands)) {
    b <- bands[[nm[i]]]; if (is.null(b)) next
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    if (ci_style %in% c("shaded", "both")) {
      p <- add_trace(p, type = "scatterpolar", mode = "lines",
                     r = b$r, theta = fck_hour_to_theta(b$hours, period),
                     fill = "toself", fillcolor = paste0(col, "26"),
                     line = list(color = "rgba(0,0,0,0)"),
                     name = paste(nm[i], "band"), showlegend = FALSE,
                     legendgroup = nm[i], hoverinfo = "skip")
    }
    if (ci_style %in% c("dotted", "both") && !is.null(b$edges)) {
      for (k in seq_along(b$edges)) {
        e <- b$edges[[k]]
        p <- add_trace(p, type = "scatterpolar", mode = "lines",
                       r = e$r, theta = fck_hour_to_theta(e$hours, period),
                       line = list(color = col, width = 1.5, dash = "dot"),
                       name = paste(nm[i], "95% CI"),
                       legendgroup = nm[i],
                       showlegend = (k == 1 && length(rings) == 1),
                       hoverinfo = "text",
                       text = sprintf("%s %s<br>%02d:%02.0f<br>%.2f", nm[i],
                                      if (k == 1) "upper" else "lower",
                                      floor(e$hours %% period), (e$hours %% 1) * 60,
                                      e$raw))
      }
    }
  }

  # --- the rings themselves ---------------------------------------------------
  # Every trace of a series carries legendgroup = that series' name. Without it
  # a legend click hides only the trace it names -- the line disappears and its
  # fill stays behind, which reads as "the curve went faint" rather than "the
  # group is off".
  # Fills first, then every line on top. Drawn ring-by-ring, the second group's
  # fill lies over the first group's line and hides exactly the comparison the
  # plot exists to make.
  if (do_fill) {
    for (i in seq_along(rings)) {
      rg <- rings[[i]]
      col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
      # A closed ring fills straight to itself. An open arc does not: 'toself'
      # would shut it with a straight chord between the two ends, inventing a
      # boundary across hours the recording never covered. Closing it back along
      # the inner radius instead makes it an annular sector, so the uncovered
      # hours read as empty rather than as a flat stretch of the curve.
      fr_hours <- rg$hours; fr_r <- rg$r
      if (!isTRUE(rg$closed %||% TRUE)) {
        fr_hours <- c(rg$hours, rev(rg$hours), rg$hours[1])
        fr_r     <- c(rg$r, rep(inner, length(rg$r)), inner)
      }
      p <- add_trace(p, type = "scatterpolar", mode = "lines",
                     r = fr_r, theta = fck_hour_to_theta(fr_hours, period),
                     fill = "toself",
                     fillcolor = paste0(col, if (length(rings) == 1) "8C" else "33"),
                     line = list(color = "rgba(0,0,0,0)"),
                     name = nm[i], showlegend = FALSE,
                     legendgroup = nm[i], hoverinfo = "skip")
    }
  }
  for (i in seq_along(rings)) {
    rg <- rings[[i]]
    col <- FCK_DENSITY_COLORS[((i - 1) %% length(FCK_DENSITY_COLORS)) + 1]
    p <- add_trace(p, type = "scatterpolar", mode = "lines",
                   r = rg$r, theta = fck_hour_to_theta(rg$hours, period),
                   name = nm[i], legendgroup = nm[i],
                   line = list(color = col, width = 2.5,
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
                     legendgroup = names(pts$groups)[i],
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
                   legendgroup = names(means)[i],
                   line = list(color = col, width = 3),
                   marker = list(color = col, size = c(1, 9)),
                   hoverinfo = "text",
                   text = sprintf("mean %02d:%02.0f<br>R-bar %.2f (n = %d)",
                                  floor(s$mean_hour), (s$mean_hour %% 1) * 60,
                                  s$r_bar, s$n))
  }

  # --- angular ticks: every hour by default ----------------------------------
  step <- suppressWarnings(as.numeric(input$density_tick_step %||% 1))
  if (!is.finite(step) || step <= 0) step <- 1
  tick_h <- seq(0, period - step, by = step)
  hour_lab <- if (period == 24) {
    # 24 labels round a circle need to be short; keep the colon only when there
    # is room for it.
    if (step >= 2) sprintf("%02d:00", tick_h) else sprintf("%02d", tick_h)
  } else sprintf("%.1f", tick_h)

  # --- radial ticks: the value scale, in the response's own units ------------
  # Placed with the SAME mapping the data went through (scale_r on axis_lo/hi),
  # so a label at 60 sits exactly where a fitted value of 60 is drawn.
  rad <- list(range = c(0, 1.13), showticklabels = FALSE, ticks = "",
              showline = FALSE, gridcolor = "rgba(11,11,11,0.08)")
  if (isTRUE(input$density_radial_labels %||% TRUE) &&
      is.finite(axis_lo) && is.finite(axis_hi) && axis_hi > axis_lo) {
    vt <- pretty(c(axis_lo, axis_hi), n = 5)
    vt <- vt[vt >= axis_lo & vt <= axis_hi]
    if (length(vt) >= 2) {
      rad <- list(
        range = c(0, 1.13),
        tickmode = "array",
        tickvals = scale_r(vt, axis_lo, axis_hi),
        ticktext = formatC(vt, format = "g", digits = 4),
        tickfont = list(size = 10, color = "#52514e"),
        angle = 112.5,          # between 12:00 and 09:00, clear of the hour ring
        tickangle = 0,
        showline = FALSE, ticks = "outside", ticklen = 3,
        gridcolor = "rgba(11,11,11,0.10)")
    }
  }

  plotly::layout(p,
    polar = list(
      bgcolor = "#fcfcfb",
      radialaxis = rad,
      angularaxis = list(
        tickmode = "array",
        tickvals = fck_hour_to_theta(tick_h, period),
        ticktext = hour_lab,
        tickfont = list(size = if (step == 1) 10 else 11),
        direction = "counterclockwise", rotation = 0,
        gridcolor = "rgba(11,11,11,0.12)")),
    annotations = if (nzchar(axis_unit)) list(list(
      text = axis_unit, showarrow = FALSE, xref = "paper", yref = "paper",
      x = 0, y = 1, xanchor = "left", yanchor = "top",
      font = list(size = 11, color = "#52514e"))) else NULL,
    showlegend = length(rings) > 1 || (ci_style %in% c("dotted", "both") && length(bands) > 0),
    legend = list(orientation = "h", y = -0.08),
    margin = list(t = 34, b = 60))
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
      if (identical(fr$trend, "none"))
        ""
      else if (isTRUE(fr$closed))
        sprintf("Showing the RHYTHM ONLY: the %s trend is not drawn. Values here differ from tab 1 by exactly the trend at that time -- they are not the fitted values. Switch to 'the fit over the recording' to read the same numbers as tab 1.\n", fr$trend)
      else
        sprintf("Showing the FIT ITSELF, %s trend included, so every value matches tab 1 at the same clock time. A trend is not periodic, so this is an open arc over the %.1f h recorded (%02d:%02.0f to %02d:%02.0f) rather than a closed ring, and the gap between its two ends is the trend across the recording.\n",
                fr$trend, diff(fr$t_range),
                floor(fr$t_range[1] %% fr$period), (fr$t_range[1] %% 1) * 60,
                floor(fr$t_range[2] %% fr$period), (fr$t_range[2] %% 1) * 60),
      if (isTRUE(fr$spans_period))
        sprintf("The recording is longer than one period (%s h), so the curve LAPS the dial: it passes over each clock hour more than once, at different values, and the gap between the passes is the trend. %s\n",
                fmt1(diff(fr$t_range)),
                if (identical(fr$lap_mode, "per_lap"))
                  "Split into one trace per day, so the passes can be told apart."
                else
                  "Drawn as one continuous spiral; switch to 'one trace per day' to tell the passes apart.")
      else "",
      # The tab claims to be tab 1 in polar coordinates. That is checkable, so
      # it is checked on every render rather than asserted in a comment.
      if (is.finite(fr$agreement))
        sprintf("Agreement with tab 1: largest difference over the whole arc = %s (identical to machine precision means the two plots ARE the same curve; if this plot still looks unlike the 2-D one, that is the lapping above, not a disagreement).\n",
                format(fr$agreement, digits = 3, scientific = TRUE))
      else "",
      if (isTRUE(fr$closed))
        "A cosinor drawn in polar coordinates is a limacon: r = MESOR + amplitude*cos(angle - acrophase) traces an OFF-CENTRE RING, widest towards the acrophase. That offset is the rhythm, not a plotting artefact.\n"
      else "",
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
      {
        pm <- input$density_profile_mode %||% "spiral"
        sp <- any(vapply(pr, function(z) isTRUE(z$spiral), logical(1)))
        if (sp)
          "Every time point, in time order, as a spiral. Nothing is averaged across days: successive passes over the same clock hour keep their own values, and the gap between them is the homeostatic rise.\n"
        else if (identical(pm, "per_day"))
          "One ring per day of the recording; the days are directly comparable.\n"
        else {
          nf <- unlist(lapply(pr, function(z) z$n_folded))
          if (!is.null(nf) && any(nf > 1))
            sprintf("Folded onto one dial: %d clock time(s) had %d observations averaged together. Under a rising homeostatic trend those observations differ by the whole rise, so the averaged value occurred on neither day. Switch to the spiral to see them separately.\n",
                    sum(nf > 1), max(nf))
          else "Folded onto one dial; no clock time was observed twice, so nothing was averaged.\n"
        }
      },
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
