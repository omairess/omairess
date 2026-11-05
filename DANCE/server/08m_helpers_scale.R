# ==============================================================================
# server/08m_helpers_scale.R -- the scale of the cosinor analysis, and the
# intervals in which nothing was sampled
#
# THE LOG SCALE is a choice the user makes on the settings panel, stated in the
# summary, the report and the exported script. When chosen, the model is fitted
# to ln(y + c): every fit, test, bound and table is on that scale. Read back on
# the original scale (exp(.) - c), the fitted curve cannot fall below -c; a
# level (MESOR, intercept, time-averaged response) is a geometric mean when
# c = 0; an amplitude A is the peak-to-trough ratio exp(2A) of its harmonic; a
# change in the trend is a ratio exp(change). Nothing here is an estimator: the
# transformation is applied to the data and read back afterwards.
#
# UNSAMPLED INTERVALS: a gap between consecutive sampling times longer than a
# threshold (the user's, or twice the median interval). The fitted curve there
# is determined by the model's form, not by data: it is drawn dashed over a
# shaded band, and a time reported inside it is flagged as extrapolated.
# CALCULATIONS_LOG.md section 24.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

dance_scale_spec <- function(type = "original", offset = 0) {
  type <- if (identical(type, "log")) "log" else "original"
  offset <- suppressWarnings(as.numeric(offset))
  list(type = type, offset = if (type == "log" && length(offset) == 1L && is.finite(offset)) offset else 0)
}
dance_scale_is_log <- function(scale) identical(scale$type, "log")

# The data on the analysis scale. Refused, with the count, when a value is not
# positive after adding c: the log of it does not exist.
dance_scale_apply <- function(Y, scale) {
  if (!dance_scale_is_log(scale)) return(list(ok = TRUE, Y = Y))
  v <- Y + scale$offset
  bad <- sum(is.finite(v) & v <= 0)
  if (bad > 0)
    return(list(ok = FALSE, n_bad = bad, message = sprintf(paste(
      "%d observed value(s) are %s 0 after adding c = %s, and have no logarithm.",
      "Choose a c that makes every value positive, or analyse on the original scale."),
      bad, "at or below", format(scale$offset))))
  list(ok = TRUE, Y = log(v))
}
dance_scale_back <- function(y, scale) if (dance_scale_is_log(scale)) exp(y) - scale$offset else y

dance_scale_label <- function(name, scale) {
  name <- if (is.null(name) || !nzchar(name)) "y" else name
  if (!dance_scale_is_log(scale)) return(name)
  if (scale$offset == 0) sprintf("ln(%s)", name) else sprintf("ln(%s + %s)", name, format(scale$offset))
}

dance_scale_statement <- function(scale, name = NULL) {
  if (!dance_scale_is_log(scale)) return("Scale: original values.")
  sprintf(paste("Scale: natural log, %s. The model, its tests, bounds and tables are on the log scale.",
                "Read back on the original scale (exp(x) - c): levels are %s, an amplitude A is the",
                "peak-to-trough ratio exp(2A) of its harmonic, a trend change is a ratio exp(change),",
                "and a difference between groups in a level is a log ratio: exp(difference) is the ratio."),
          dance_scale_label(name, scale),
          if (scale$offset == 0) "geometric means" else sprintf("exp(level) - %s", format(scale$offset)))
}

# The original-scale reading of each estimate, in its note.
dance_outcomes_backtransform <- function(o, scale, digits = 3) {
  if (!dance_scale_is_log(scale) || !isTRUE(o$ok) || is.null(o$estimates)) return(o)
  es <- o$estimates
  f <- function(x) formatC(x, format = "f", digits = digits)
  ci <- function(e, l, h, g) if (is.finite(l) && is.finite(h)) sprintf("%s [%s, %s]", f(g(e)), f(g(l)), f(g(h))) else f(g(e))
  lev <- es$outcome %in% c("intercept", "mesor", "fitted_mean")
  amp <- grepl("^amplitude[0-9]+$", es$outcome) | es$outcome %in% c("circadian")
  chg <- es$outcome %in% c("change", "sat_change")
  for (i in which(lev | amp | chg)) {
    r <- es[i, ]
    if (!is.finite(r$estimate)) next
    txt <- if (lev[i]) sprintf("original scale: %s %s", if (scale$offset == 0) "geometric mean" else sprintf("exp(x) - %s", format(scale$offset)),
                               ci(r$estimate, r$lo, r$hi, function(x) exp(x) - scale$offset))
           else if (amp[i]) sprintf("peak/trough ratio exp(2A) %s", ci(r$estimate, r$lo, r$hi, function(x) exp(2 * x)))
           else sprintf("ratio end/start exp(change) %s", ci(r$estimate, r$lo, r$hi, exp))
    es$note[i] <- if (is.na(r$note) || !nzchar(r$note)) txt else paste(r$note, txt, sep = "; ")
  }
  o$estimates <- es
  o$scale <- scale
  o
}

# ---- UNSAMPLED INTERVALS -------------------------------------------------------
# times: the sampling times on a linear axis; Y (optional): the data, so a time
# column where nobody was observed is not counted as sampled.
dance_unsampled_intervals <- function(times, Y = NULL, gap_h = NA_real_, names = NULL) {
  t <- as.numeric(times)
  nm <- if (!is.null(names) && length(names) == length(t)) as.character(names) else rep(NA_character_, length(t))
  if (!is.null(Y)) { keep <- colSums(is.finite(as.matrix(Y))) > 0; t <- t[keep]; nm <- nm[keep] }
  o <- order(t); t <- t[o]; nm <- nm[o]
  ok <- is.finite(t) & !duplicated(t); t <- t[ok]; nm <- nm[ok]
  empty <- data.frame(start = numeric(0), end = numeric(0), start_col = character(0), end_col = character(0))
  # 0 switches the marking off; blank (NA) is the default rule
  if (length(gap_h) == 1L && is.finite(gap_h) && gap_h == 0)
    return(structure(empty, threshold = Inf, rule = "marking switched off (threshold 0)"))
  if (length(t) < 3L) return(structure(empty, threshold = NA_real_, rule = "fewer than three sampling times"))
  d <- diff(t)
  user <- length(gap_h) == 1L && is.finite(gap_h) && gap_h > 0
  thr <- if (user) gap_h else 2 * stats::median(d)
  k <- which(d > thr + 1e-9)
  structure(data.frame(start = t[k], end = t[k + 1L], start_col = nm[k], end_col = nm[k + 1L],
                       stringsAsFactors = FALSE), threshold = thr,
            rule = if (user) sprintf("a gap longer than %s h (set by the user)", format(thr))
                   else sprintf("a gap longer than twice the median sampling interval (2 x %s = %s h)",
                                format(stats::median(d)), format(thr)))
}

# TRUE when time t (on the same axis) lies strictly inside a gap.
dance_in_gap <- function(t, gaps) {
  if (!nrow(gaps)) return(rep(FALSE, length(t)))
  vapply(t, function(x) is.finite(x) && any(x > gaps$start & x < gaps$end), logical(1))
}

# A time of day that recurs every `step` hours: extrapolated when no occurrence
# within the sampled window lies outside the gaps. gaps and window are on the
# clock-linear axis (hours from midnight of the first day).
dance_clock_in_gap <- function(clock_h, gaps, window, step = 24) {
  if (!is.finite(clock_h) || !nrow(gaps)) return(FALSE)
  k <- seq(floor((window[1] - clock_h) / step) - 1, ceiling((window[2] - clock_h) / step) + 1)
  occ <- clock_h + step * k
  occ <- occ[occ >= window[1] & occ <= window[2]]
  length(occ) > 0 && all(dance_in_gap(occ, gaps))
}

dance_gap_text <- function(gaps, shift = 0, columns = FALSE) {
  if (!nrow(gaps)) return("none")
  # the day is shown: 10:00 on the second day is not 10:00 on the first
  lab <- function(t) dance_clock_label(c(shift, t + shift), 24, show_day = TRUE)[-1]
  txt <- sprintf("%s-%s", lab(gaps$start), lab(gaps$end))
  if (columns && !is.null(gaps$start_col) && any(!is.na(gaps$start_col)))
    txt <- sprintf("%s (no column between '%s' and '%s')", txt, gaps$start_col, gaps$end_col)
  paste(txt, collapse = ", ")
}

# Acrophases that fall only in unsampled intervals get a note.
dance_outcomes_gap_notes <- function(o, gaps_clock, window_clock, period = 24) {
  if (!isTRUE(o$ok) || is.null(o$estimates) || !nrow(gaps_clock)) return(o)
  es <- o$estimates
  for (i in which(es$kind == "circular" & is.finite(es$estimate))) {
    h <- if (is.finite(es$harmonic[i])) es$harmonic[i] else 1L
    clk <- dance_acrophase_clock(hours = es$estimate[i], period = period, harmonic = h,
                                 clock_origin = o$clock_origin %||% 0)$hours
    if (isTRUE(dance_clock_in_gap(clk, gaps_clock, window_clock, step = period / h))) {
      txt <- "in an unsampled interval: extrapolated, not observed"
      es$note[i] <- if (is.na(es$note[i]) || !nzchar(es$note[i])) txt else paste(es$note[i], txt, sep = "; ")
    }
  }
  o$estimates <- es
  o
}

# ---- THE PLOTS -----------------------------------------------------------------
# Applied to a finished plot: every trace's y, every band, every reference line
# and annotation drawn on the y axis, and the number after "label: " in the
# hover text, all read back on the original scale together.
dance_plot_backtransform <- function(p, scale, ylab) {
  if (!dance_scale_is_log(scale)) return(p)
  b <- plotly::plotly_build(p)
  g <- function(y) { y <- suppressWarnings(as.numeric(y)); exp(y) - scale$offset }
  for (i in seq_along(b$x$data)) {
    tr <- b$x$data[[i]]
    if (is.null(tr$y) || !is.numeric(unlist(tr$y))) next
    yo <- unlist(tr$y); tr$y <- g(yo)
    tx <- tr$text
    if (is.character(tx) && length(tx) == length(yo))
      tr$text <- vapply(seq_along(tx), function(j)
        if (is.finite(tr$y[j])) sub("(<br>[^<:]+: )-?[0-9]+(\\.[0-9]+)?", sprintf("\\1%.2f", tr$y[j]), tx[j]) else tx[j], "")
    b$x$data[[i]] <- tr
  }
  sh <- b$x$layout$shapes
  for (i in seq_along(sh)) if (identical(sh[[i]]$yref %||% "y", "y")) {
    if (!is.null(sh[[i]]$y0)) sh[[i]]$y0 <- g(sh[[i]]$y0)
    if (!is.null(sh[[i]]$y1)) sh[[i]]$y1 <- g(sh[[i]]$y1)
  }
  b$x$layout$shapes <- sh
  an <- b$x$layout$annotations
  for (i in seq_along(an)) if (identical(an[[i]]$yref %||% "y", "y") && is.numeric(an[[i]]$y)) an[[i]]$y <- g(an[[i]]$y)
  b$x$layout$annotations <- an
  b$x$layout$yaxis$title <- list(text = ylab)
  b
}

# Shade each gap and draw a fitted curve dashed across it. A curve is a line
# trace whose x values are not all sampling times (raw-data lines are).
dance_plot_mark_gaps <- function(p, gaps, sample_times) {
  if (is.null(gaps) || !nrow(gaps)) return(p)
  b <- plotly::plotly_build(p)
  out <- list()
  for (tr in b$x$data) {
    out[[length(out) + 1L]] <- tr
    x <- suppressWarnings(as.numeric(unlist(tr$x)))
    mode <- tr$mode %||% ""
    fill <- tr$fill %||% "none"
    if (!grepl("lines", mode) || !identical(fill, "none") || length(x) < 3L || all(x %in% sample_times)) next
    y <- unlist(tr$y)
    ing <- dance_in_gap(x, gaps)
    if (!any(ing)) next
    # the dashed copy keeps the first point either side, so the pieces join
    edge <- ing | c(ing[-1], FALSE) | c(FALSE, ing[-length(ing)])
    solid <- tr; solid$y <- ifelse(ing, NA, y); solid$connectgaps <- FALSE
    dash <- tr; dash$y <- ifelse(edge, y, NA); dash$showlegend <- FALSE; dash$connectgaps <- FALSE
    dash$line <- utils::modifyList(tr$line %||% list(), list(dash = "dot"))
    if (is.character(tr$text) && length(tr$text) == length(y))
      dash$text <- paste0(tr$text, "<br><i>no samples here: extrapolated</i>")
    out[[length(out)]] <- solid
    out[[length(out) + 1L]] <- dash
  }
  b$x$data <- out
  sh <- b$x$layout$shapes %||% list()
  an <- b$x$layout$annotations %||% list()
  for (i in seq_len(nrow(gaps))) {
    sh[[length(sh) + 1L]] <- list(type = "rect", xref = "x", yref = "paper", x0 = gaps$start[i], x1 = gaps$end[i],
                                  y0 = 0, y1 = 1, fillcolor = "rgba(120,120,120,0.10)", line = list(width = 0),
                                  layer = "below")
    an[[length(an) + 1L]] <- list(x = (gaps$start[i] + gaps$end[i]) / 2, y = 1, xref = "x", yref = "paper",
                                  text = if (!is.null(gaps$label)) gaps$label[i] else "no samples",
                                  showarrow = FALSE, yanchor = "bottom",
                                  font = list(size = 10, color = "#777"),
                                  hovertext = attr(gaps, "hover")[i] %||% "no sampling time in this interval")
  }
  b$x$layout$shapes <- sh
  b$x$layout$annotations <- an
  b
}

# Both annotations in one step, as the app's tables and the exported script
# apply them: the original-scale reading (log runs), then the acrophases that
# fall only where nothing was sampled. gaps and time_vec are in model time;
# origin_shift is the clock time of model time 0.
dance_outcomes_annotate <- function(o, scale, gaps, origin_shift = 0, time_vec = NULL, period = 24) {
  o <- dance_outcomes_backtransform(o, scale %||% dance_scale_spec())
  if (!is.null(gaps) && nrow(gaps) && length(time_vec)) {
    sh <- as.numeric(origin_shift %||% 0)
    o <- dance_outcomes_gap_notes(o, data.frame(start = gaps$start + sh, end = gaps$end + sh),
                                  range(time_vec, na.rm = TRUE) + sh, period)
  }
  o
}
