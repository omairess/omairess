# ==============================================================================
# server/08h_helpers_twostage.R -- the TWO-STAGE comparison, as the same
# questions tab 6 asks of the mixed model
# ==============================================================================
# One cosinor per participant, then the point estimates compared between
# groups. Nothing here is new: the omnibus is a one-way MANOVA (Wilks) on the
# per-participant coefficient vectors and one-way ANOVAs on the scalars; the
# rhythm is Bingham et al.'s (1982) population-mean cosinor; the acrophase is
# Watson-Williams; the pairwise tests are Welch's t and Watson-Williams with a
# Holm family. These are the tests the removed "Cosinor: pairwise tests" module
# and the legacy tab-6 block ran, gathered into one place and given the same
# shape as the mixed-effects panels, so a reader switching approach sees the
# same layout answering the same questions from a different estimator.
#
# What two-stage CANNOT do is stated where it applies: it treats every
# per-participant estimate as exact, and it cannot respect a within-participant
# factor -- curves of one participant enter as if independent.
# ==============================================================================

# The per-curve parameter frame, with the group attached and the rows the
# summaries use: converged fits with a usable group label.
dance_ts_frame <- function(mod, group_values) {
  p <- mod$individual_params
  if (is.null(p) || !nrow(p)) return(NULL)
  g <- if (is.null(group_values)) rep("(all)", nrow(p)) else as.character(group_values[p$subject])
  p$group <- droplevels(factor(g))
  p <- p[!is.na(p$group), , drop = FALSE]
  if (!nrow(p)) return(NULL)
  # each participant's non-periodic change over the recording, from their own
  # trend parameters -- the number the "change from origin" curve ends at.
  # Stored at fit time; recomputed for a model saved before it was.
  if (!"trend_change" %in% names(p) || all(is.na(p$trend_change)))
    p$trend_change <- dance_ts_trend_change(mod, p)
  # the initial slope A/tau of a saturating rise: how fast it starts, in
  # response units per hour -- stored at fit time, derived for older models
  if (identical(mod$trend_type, "exp_sat") && all(c("A_sat", "tau") %in% names(p)) &&
      !"init_slope" %in% names(p))
    p$init_slope <- p$A_sat / p$tau
  p
}

# S(last observation) - S(first observation) for every participant, through the
# same trend function the two-stage fits use (dance_trend_value), with each
# fit's own anchoring. Zero without a trend; NA where a fit lacks a parameter.
dance_ts_trend_change <- function(mod, p) {
  tr <- as.character(mod$trend_type %||% "none")
  if (identical(tr, "none")) return(rep(0, nrow(p)))
  tv <- mod$time_vec
  if (is.null(tv) || !any(is.finite(tv))) return(rep(NA_real_, nrow(p)))
  tt <- range(tv, na.rm = TRUE)
  vapply(seq_len(nrow(p)), function(i) {
    coefs <- switch(tr,
      linear  = p$trend_linear[i],
      log     = p$trend_log[i],
      exp_sat = c(p$A_sat[i], p$tau[i]),
      NULL)
    if (is.null(coefs) || !length(coefs) || !all(is.finite(coefs))) return(NA_real_)
    fi <- tryCatch(mod$individual_fits[[p$subject[i]]], error = function(e) NULL)
    off <- fi$t_offset %||% mod$t_offset %||% 0
    s <- dance_trend_value(tr, coefs, tt, off)
    s[2] - s[1]
  }, numeric(1))
}

# The level quantities per group, as means of the per-participant estimates with
# a t interval: the intercept (each fit's constant), and with a trend the MESOR
# (each fit's midline averaged over the window, mesor_adj) and the non-periodic
# change over the recording. The pairwise differences come from
# dance_ts_pairwise(), the same Welch tests as the pairwise table.
dance_ts_level_summary <- function(df, mod, conf = 0.95) {
  has_trend <- !identical(mod$trend_type %||% "none", "none")
  qs <- c(intercept = "mesor")
  if (has_trend) qs <- c(qs, mesor = "mesor_adj", change = "trend_change")
  out <- dance_ts_group_means(df, qs, conf)
  if (isTRUE(out$ok)) out$has_trend <- has_trend
  out
}

# The saturating trend's three parameters per group, as means of the
# per-participant fits: the rise A, the time constant tau and the initial slope
# A/tau. Each participant's fit estimates A and tau jointly, so these are the
# two-stage answer to "do the groups differ in tau" -- with the caveat the
# per-participant fits carry: over a short window A and tau trade off along a
# ridge, and A/tau and the change are steadier than either alone.
dance_ts_sat_summary <- function(df, mod, conf = 0.95) {
  if (!identical(mod$trend_type, "exp_sat"))
    return(list(ok = FALSE, message = "No saturating trend in this model."))
  if (!"init_slope" %in% names(df) && all(c("A_sat", "tau") %in% names(df)))
    df$init_slope <- df$A_sat / df$tau
  qs <- c(A_sat = "A_sat", tau = "tau", init_slope = "init_slope")
  if (dance_ts_tau_held(mod)) qs <- qs[names(qs) != "tau"]
  dance_ts_group_means(df, qs, conf)
}

# The saturating trend per group from the per-participant fits: each
# participant's own A (1 - exp(-t/tau)) from the first observation, averaged
# within the group, with a pointwise band from the spread across participants
# -- the two-stage counterpart of the mixed model's per-cell trend curves.
dance_ts_sat_curves <- function(df, mod, n_time = 160, conf = 0.95) {
  tv <- mod$time_vec
  if (is.null(tv) || !all(c("A_sat", "tau") %in% names(df))) return(NULL)
  tt <- seq(0, diff(range(tv, na.rm = TRUE)), length.out = n_time)
  g <- droplevels(factor(df$group))
  q <- stats::qnorm(1 - (1 - conf) / 2)
  out <- list()
  for (lv in levels(g)) {
    i <- which(g == lv & is.finite(df$A_sat) & is.finite(df$tau) & df$tau > 0)
    if (length(i) < 2L) next
    M <- t(vapply(i, function(r) df$A_sat[r] * (1 - exp(-tt / df$tau[r])), numeric(n_time)))
    m <- colMeans(M); se <- apply(M, 2, stats::sd) / sqrt(nrow(M))
    out[[lv]] <- data.frame(cell = lv, tt = tt, fit = m, se = se, lo = m - q * se,
                            hi = m + q * se, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(table = tab, cells = names(out), conf = conf)
}

# Group means of per-participant columns with a t interval on n - 1 df.
dance_ts_group_means <- function(df, qs, conf = 0.95) {
  qs <- qs[qs %in% names(df)]
  g <- droplevels(factor(df$group))
  rows <- list()
  for (q in names(qs)) for (lv in levels(g)) {
    x <- df[[qs[[q]]]][g == lv]; x <- x[is.finite(x)]
    n <- length(x)
    m <- if (n) mean(x) else NA_real_
    se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
    tq <- if (n > 1) stats::qt(1 - (1 - conf) / 2, n - 1) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(quantity = q, cell = lv, n = n,
                                            estimate = m, se = se,
                                            lo = m - tq * se, hi = m + tq * se,
                                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(list(ok = FALSE, message = "No per-participant estimates to summarise."))
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  list(ok = TRUE, table = tab, quantities = names(qs), conf = conf)
}

# The coefficient columns a MANOVA compares: constant, trend, and every
# (cos, sin) pair -- the complete per-participant fit.
# TRUE when every participant's fit held tau at one value (Model
# Specification): tau is then a constant of the design, not an estimate, and
# has no place in a test.
dance_ts_tau_held <- function(mod)
  identical(mod$trend_type, "exp_sat") && identical(mod$tau_spec$mode, "fixed")

# How many participants' tau ran to an end of the search range -- a bound the
# optimiser was stopped at, not an estimate.
dance_ts_tau_at_bound <- function(mod, subjects) {
  hits <- vapply(subjects, function(sid) {
    f <- tryCatch(mod$individual_fits[[sid]], error = function(e) NULL)
    any(grepl("^tau \\(", f$bounds_hit %||% character(0)))
  }, logical(1))
  sum(hits)
}

dance_ts_coef_cols <- function(mod) {
  tr <- switch(mod$trend_type %||% "none",
               linear = "trend_linear", log = "trend_log",
               exp_sat = c("A_sat", if (!dance_ts_tau_held(mod)) "tau"), character(0))
  c("mesor", tr, unlist(lapply(seq_len(mod$n_harmonics %||% 1L),
                                function(h) paste0(c("beta_cos_", "beta_sin_"), h))))
}

# ---- PRIMARY: do the groups differ in the fitted curve at all? -------------
# Wilks' lambda MANOVA on the coefficient vectors. Base R, stats::manova.
dance_ts_manova <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (!length(cols)) return(list(ok = FALSE, message = "No coefficient columns."))
  Y <- as.matrix(df[, cols, drop = FALSE]); g <- droplevels(factor(df$group))
  ok <- stats::complete.cases(Y) & !is.na(g)
  Y <- Y[ok, , drop = FALSE]; g <- droplevels(g[ok])
  k <- nlevels(g)
  if (k < 2) return(list(ok = FALSE, message = "Fewer than two groups have usable fits."))
  if (any(table(g) <= ncol(Y)))
    return(list(ok = FALSE, message = sprintf(
      "Every group needs more participants than the %d coefficients compared; sizes are %s.",
      ncol(Y), paste(sprintf("%s = %d", levels(g), as.integer(table(g))), collapse = ", "))))
  if (ncol(Y) == 1L) {
    lt <- dance_group_linear_test(Y[, 1], g)
    return(list(ok = TRUE, F = lt$F, df1 = lt$df1, df2 = lt$df2, p = lt$p,
                wilks = NA_real_, n_par = 1L, n_groups = k, n = nrow(Y),
                method = "one-way ANOVA"))
  }
  fit <- tryCatch(stats::manova(Y ~ g), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, message = "The MANOVA could not be fitted."))
  st <- tryCatch(summary(fit, test = "Wilks")$stats, error = function(e) NULL)
  if (is.null(st)) return(list(ok = FALSE, message = "Wilks' lambda could not be computed."))
  list(ok = TRUE, F = unname(st[1, "approx F"]), df1 = unname(st[1, "num Df"]),
       df2 = unname(st[1, "den Df"]), p = unname(st[1, "Pr(>F)"]),
       wilks = unname(st[1, "Wilks"]), n_par = ncol(Y), n_groups = k, n = nrow(Y),
       cols = cols, method = "one-way MANOVA, Wilks' lambda")
}

# ---- SECONDARY: the components, one test each -----------------------------
dance_ts_components <- function(df, mod) {
  period <- mod$period %||% 24; K <- mod$n_harmonics %||% 1L
  rows <- list()
  push <- function(...) rows[[length(rows) + 1L]] <<- list(...)
  lt_row <- function(label, x, note = NULL) {
    lt <- dance_group_linear_test(x, df$group)
    if (is.null(lt)) push(label = label, ok = FALSE, message = "not testable (too few usable values)")
    else push(label = label, ok = TRUE, F = lt$F, df1 = lt$df1, df2 = lt$df2, p = lt$p,
              omega2 = lt$omega2, method = "one-way ANOVA", note = note)
  }
  trend <- mod$trend_type %||% "none"
  has_trend <- !identical(trend, "none")
  lab <- function(w) dance_ts_pair_label(w, period, trend, mod$dv_units)
  lt_row(lab("intercept"), df$mesor,
         if (has_trend) "the midline at the first observation, where the trend is zero" else NULL)
  if (has_trend && "mesor_adj" %in% names(df))
    lt_row(lab("mesor"), df$mesor_adj)
  if (has_trend && "trend_change" %in% names(df))
    lt_row(lab("change"), df$trend_change)
  for (tc in setdiff(dance_ts_coef_cols(mod), c("mesor", grep("^beta_", dance_ts_coef_cols(mod), value = TRUE))))
    if (tc %in% names(df)) lt_row(lab(tc), df[[tc]])
  if (identical(trend, "exp_sat") && "init_slope" %in% names(df))
    lt_row(lab("init_slope"), df$init_slope)
  for (h in seq_len(K)) {
    bc <- paste0("beta_cos_", h); bs <- paste0("beta_sin_", h)
    if (!all(c(bc, bs) %in% names(df))) next
    pop <- dance_pop_cosinor(df[[bc]], df[[bs]], df$mesor, df$group, period = period, harmonic = h)
    if (!isTRUE(pop$ok)) {
      push(label = sprintf("H%d rhythm", h), ok = FALSE, message = pop$message); next
    }
    j <- pop$joint
    if (!is.null(j) && isTRUE(j$ok))
      push(label = sprintf("H%d rhythm vector (cos, sin) jointly", h), ok = TRUE,
           F = j$F, df1 = j$df1, df2 = j$df2, p = j$p, method = "MANOVA, Wilks' lambda",
           harmonic = h, kind = "joint")
    tt <- pop$tests
    push(label = sprintf("H%d amplitude", h), ok = TRUE,
         F = tt$F[tt$parameter == "Amplitude"], df1 = pop$df1, df2 = pop$df2,
         p = tt$p[tt$parameter == "Amplitude"], method = "population-mean cosinor (Bingham 1982)",
         harmonic = h, kind = "amplitude",
         note = if (!isTRUE(pop$amplitude_interpretable))
           "not interpretable as a difference in rhythm strength: the groups' acrophases differ, or the acrophase test had no power at this separation (Bingham's caution)" else NULL)
    push(label = sprintf("H%d acrophase", h), ok = TRUE,
         F = tt$F[tt$parameter == "Acrophase"], df1 = pop$df1, df2 = pop$df2,
         p = tt$p[tt$parameter == "Acrophase"], method = "population-mean cosinor (Bingham 1982)",
         harmonic = h, kind = "acrophase",
         note = if (!isTRUE(pop$acrophase_test_supported)) pop$acrophase_test_note else NULL)
    # Watson-Williams on the per-participant phases: every participant counts
    # equally, and it has no half-cycle blind spot
    rad_col <- paste0("acrophase_rad_", h)
    rad <- if (rad_col %in% names(df)) df[[rad_col]] else
      df[[paste0("acrophase_time_", h)]] * 2 * pi / (period / h)
    sp <- split(rad, df$group); sp <- sp[vapply(sp, function(v) sum(is.finite(v)) >= 2, logical(1))]
    if (length(sp) >= 2) {
      ww <- tryCatch(dance_watson_williams_test(lapply(sp, function(v) v[is.finite(v)])),
                     error = function(e) NULL)
      if (!is.null(ww)) {
        as_ <- dance_ww_assumption(ww$r_bar, ww$kappa)
        push(label = sprintf("H%d acrophase (Watson-Williams, unweighted)", h), ok = TRUE,
             F = ww$F, df1 = ww$df1, df2 = ww$df2, p = ww$p, method = "Watson-Williams F",
             harmonic = h, kind = "ww", note = if (!isTRUE(as_$ok)) as_$msg else NULL)
      }
    }
    rows[[length(rows)]]$pop <- pop
  }
  rows
}

# The per-group rhythm, per harmonic, as the vector mean of the (cos, sin)
# pairs -- Bingham's group means, identical to the app's group_fits vector means.
dance_ts_group_rhythm <- function(df, mod, h) {
  bc <- paste0("beta_cos_", h); bs <- paste0("beta_sin_", h)
  if (!all(c(bc, bs) %in% names(df))) return(NULL)
  period <- mod$period %||% 24
  g <- droplevels(factor(df$group))
  out <- do.call(rbind, lapply(levels(g), function(lv) {
    i <- g == lv & is.finite(df[[bc]]) & is.finite(df[[bs]])
    a <- mean(df[[bc]][i]); b <- mean(df[[bs]][i])
    rad <- if (all(is.finite(c(a, b)))) df[[paste0("acrophase_rad_", h)]][i] else NA
    rr <- if (sum(i) >= 2) dance_resultants(rad[is.finite(rad)]) else NULL
    data.frame(cell = lv, n = sum(i), beta_cos = a, beta_sin = b,
               amplitude = sqrt(a^2 + b^2), acrophase_rad = atan2(b, a) %% (2 * pi),
               acrophase_time = (atan2(b, a) %% (2 * pi)) * (period / h) / (2 * pi),
               r_bar = if (!is.null(rr)) rr$r_unweighted else NA_real_,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL; out
}

# ---- PAIRWISE ---------------------------------------------------------------
# The same names as the mixed-effects vocabulary for the quantities both
# approaches have (intercept, mesor, change, amplitudeK, phaseK), so a selection
# survives switching approach; the two-stage-only ones keep their column names.
# The saturating trend has TWO per-participant parameters here, A_sat and tau,
# because each participant's fit estimates both -- unlike the mixed model, which
# holds tau fixed and estimates only the rise.
dance_ts_pair_whats <- function(mod) {
  K <- mod$n_harmonics %||% 1L
  tr <- switch(mod$trend_type %||% "none", linear = "trend_linear", log = "trend_log",
               exp_sat = c("A_sat", if (!dance_ts_tau_held(mod)) "tau", "init_slope"),
               character(0))
  out <- c("intercept", if (length(tr)) c("mesor", "change"), tr,
           sprintf("amplitude%d", seq_len(K)), sprintf("phase%d", seq_len(K)), "r_squared")
  out
}
dance_ts_pair_label <- function(what, period = 24, trend = NULL, units = NULL) {
  # only a harmonic name needs parsing; the level and trend labels stand alone
  if (grepl("^(amplitude|phase)[0-9]*$", what)) {
    pp <- dance_traj_pair_parse(what)
    return(sprintf("%s — H%d (%s h)", if (identical(pp$base, "amplitude")) "Amplitude" else "Acrophase",
                   pp$harmonic, format(round(period / pp$harmonic, 2), trim = TRUE)))
  }
  u <- if (is.null(units) || !nzchar(units)) "response units" else units
  switch(what,
    level = , intercept = if (identical(trend, "none")) "Intercept = MESOR (no trend)"
                          else "Intercept (midline at the first observation)",
    mesor = , mesor_adj = "MESOR (midline averaged over the window)",
    change = , trend_change = "Non-periodic change over the recording",
    value_at_start = "Curve value at the first observation",
    trend_linear = sprintf("Slope beta (%s per hour)", u),
    trend_log = sprintf("Log-trend beta (%s per unit of log(1 + hours))", u),
    A_sat = sprintf("Saturating rise A (%s above the intercept)", u),
    tau = "Time constant tau of the saturating rise (hours)",
    init_slope = sprintf("Initial slope A/tau (%s per hour)", u),
    r_squared = "R-squared of the participant fit", what)
}

dance_ts_pairwise <- function(df, mod, what, adjust = c("holm", "bonferroni", "none"),
                              conf = 0.95) {
  adjust <- match.arg(adjust)
  period <- mod$period %||% 24
  pp <- dance_traj_pair_parse(what)
  circular <- identical(pp$base, "phase")
  col <- if (!is.na(pp$harmonic)) {
    if (circular) paste0("acrophase_rad_", pp$harmonic) else paste0("amplitude_", pp$harmonic)
  } else switch(what, level = , intercept = "mesor", mesor = "mesor_adj",
                change = "trend_change", what)
  if (identical(col, "trend_change") && !col %in% names(df))
    df$trend_change <- dance_ts_trend_change(mod, df)
  if (!col %in% names(df)) {
    # phases may be stored as time only
    if (circular && paste0("acrophase_time_", pp$harmonic) %in% names(df)) {
      df[[col]] <- df[[paste0("acrophase_time_", pp$harmonic)]] * 2 * pi / (period / pp$harmonic)
    } else return(list(ok = FALSE, message = sprintf("'%s' is not in the per-participant results.", col)))
  }
  g <- droplevels(factor(df$group)); lv <- levels(g)
  if (length(lv) < 2) return(list(ok = FALSE, message = "Only one group."))
  eff <- if (circular) period / pp$harmonic else NA_real_
  pairs <- utils::combn(length(lv), 2)
  rows <- lapply(seq_len(ncol(pairs)), function(k) {
    a <- lv[pairs[1, k]]; b <- lv[pairs[2, k]]
    v1 <- df[[col]][g == a]; v2 <- df[[col]][g == b]
    v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
    if (length(v1) < 2 || length(v2) < 2) return(NULL)
    if (circular) {
      m1 <- dance_circular_mean(v1); m2 <- dance_circular_mean(v2)
      d <- atan2(sin(m1 - m2), cos(m1 - m2))
      ww <- tryCatch(dance_watson_williams_test(list(v1, v2)), error = function(e) NULL)
      if (is.null(ww)) return(NULL)
      as_ <- dance_ww_assumption(ww$r_bar, ww$kappa)
      data.frame(cell1 = a, cell2 = b, n1 = length(v1), n2 = length(v2),
                 estimate = d * eff / (2 * pi), se = NA_real_, lo = NA_real_, hi = NA_real_,
                 statistic = ww$F, df1 = ww$df1, df2 = ww$df2, p_raw = ww$p,
                 defined = TRUE, assumption_ok = isTRUE(as_$ok), d = NA_real_,
                 stringsAsFactors = FALSE)
    } else {
      tt <- tryCatch(stats::t.test(v1, v2, var.equal = FALSE, conf.level = conf),
                     error = function(e) NULL)
      if (is.null(tt)) return(NULL)
      sp <- sqrt(((length(v1) - 1) * stats::var(v1) + (length(v2) - 1) * stats::var(v2)) /
                   (length(v1) + length(v2) - 2))
      data.frame(cell1 = a, cell2 = b, n1 = length(v1), n2 = length(v2),
                 estimate = mean(v1) - mean(v2), se = unname(tt$stderr),
                 lo = tt$conf.int[1], hi = tt$conf.int[2],
                 statistic = unname(tt$statistic), df1 = 1, df2 = unname(tt$parameter),
                 p_raw = tt$p.value, defined = TRUE, assumption_ok = TRUE,
                 d = if (is.finite(sp) && sp > 0) (mean(v1) - mean(v2)) / sp else NA_real_,
                 stringsAsFactors = FALSE)
    }
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(list(ok = FALSE, message = "No pair has two usable values on both sides."))
  tab <- do.call(rbind, rows)
  tab$p_adj <- if (identical(adjust, "none")) tab$p_raw else stats::p.adjust(tab$p_raw, method = adjust)
  rownames(tab) <- NULL
  list(ok = TRUE, what = what, circular = circular, table = tab, adjust = adjust, conf = conf,
       joint = FALSE,
       label = dance_ts_pair_label(what, period, mod$trend_type %||% "none", mod$dv_units),
       method = if (circular) "Watson-Williams F on the per-participant phases"
                else "Welch's t on the per-participant estimates",
       unit = if (circular)
         sprintf("acrophase difference in hours on an effective period of %.4g h (shortest arc); no interval -- the test is the statement", eff)
       else paste0("difference of group means of the per-participant estimates",
                   if (!is.null(mod$dv_units) && nzchar(mod$dv_units) &&
                       (what %in% c("intercept", "level", "mesor", "change", "A_sat") ||
                        identical(pp$base, "amplitude"))) paste0(", in ", mod$dv_units)
                   else if (identical(what, "tau")) ", in hours"
                   else if (identical(what, "init_slope") && !is.null(mod$dv_units) &&
                            nzchar(mod$dv_units)) paste0(", in ", mod$dv_units, " per hour")
                   else "",
                   ", with Welch's interval and Cohen's d",
                   if (identical(what, "change"))
                     ". Each participant's change uses their own trend parameters" else ""),
       note = paste(
         sprintf("%d pairwise comparisons over %d groups, on the per-participant point estimates.",
                 nrow(tab), length(lv)),
         "Two-stage: each estimate is treated as exact, so a participant with a poorly",
         "determined rhythm counts as much as one with a precise one.",
         if (identical(adjust, "none")) "NO multiplicity adjustment."
         else sprintf("p_adj is %s-adjusted across the family of %d.", adjust, nrow(tab)),
         if (circular && any(!tab$assumption_ok))
           "At least one pair fails the Watson-Williams concentration assumption (r-bar < 0.45); that F should not be reported." else NULL))
}

# ---- THE GROUP CURVES, with a band that is a standard error ----------------
# The line is the curve of the group's mean coefficients -- the same line tab 1
# draws. The band is that line +/- t * SE(t), where SE(t) is the standard error
# across participants of their own fitted curves at time t. That replaces the
# app's earlier approximation, which rescaled the whole curve by the amplitude's
# standard error and was neither a pointwise nor a simultaneous interval.
# component -- the views the mixed approach offers, on the same scale:
#   "full"          intercept + trend + harmonics
#   "baseline_harm" intercept + harmonics, trend removed ("rhythm" is accepted)
#   "harmonics"     the harmonics alone, centred on zero
#   "trend"         the non-periodic change from the first observation:
#                   trend(t) - trend(t0), zero at the start by construction
#   "harmonicK"     harmonic K alone, centred on zero
# include_trend is the older switch (FALSE = "baseline_harm"), kept for callers.
#
# band -- "pointwise" (mean curve +/- z SE(t)) or "simultaneous": the two-stage
# counterpart of the mixed model's Scheffe band. When the curve is a fixed
# linear combination x(t)'b of the coefficients (no trend, or a linear or log
# one), each participant's b_i is one draw, and Hotelling's T^2 gives a band
# that holds for every t at once: x(t)'bbar +/- c * sqrt(x(t)' S x(t) / n),
# c^2 = d (n - 1) / (n - d) F(d, n - d), d the number of coefficients the view
# uses (Scheffe 1953 with Hotelling's T^2). sqrt(x' S x / n) is the SD of the
# participants' curves at t over sqrt(n), the pointwise SE. With a saturating
# trend every participant has their own tau, the curves are not linear in shared
# coefficients, and the band stays pointwise (simultaneous_ok = FALSE).
dance_ts_group_curves <- function(mod, group_values, times, include_trend = TRUE,
                                  conf = 0.95, component = NULL, band = "pointwise") {
  df <- dance_ts_frame(mod, group_values)
  if (is.null(df)) return(NULL)
  period <- mod$period %||% 24; K <- mod$n_harmonics %||% 1L
  trend <- mod$trend_type %||% "none"
  component <- component %||% (if (isTRUE(include_trend)) "full" else "baseline_harm")
  if (identical(component, "rhythm")) component <- "baseline_harm"
  t0 <- min(mod$time_vec, na.rm = TRUE)
  n_tr <- switch(trend, none = 0L, linear = 1L, log = 1L, exp_sat = 2L, 0L)
  kh <- dance_traj_component_harmonic(component)
  # coefficients the view uses (for the simultaneous band)
  d_view <- if (!is.na(kh)) 2L else switch(component,
    full = 1L + n_tr + 2L * K, baseline_harm = 1L + 2L * K, harmonics = 2L * K,
    trend = n_tr, 2L * K)
  linear_view <- !identical(trend, "exp_sat") || component %in% c("baseline_harm", "harmonics") ||
                 !is.na(kh)
  # one participant's (or one coefficient vector's) curve for the component
  curve_of <- function(coefs, t_offset) {
    if (!is.na(kh)) {
      off <- 1L + n_tr; w <- 2 * pi * kh / period
      return(coefs[off + 2L * kh - 1L] * cos(w * times) + coefs[off + 2L * kh] * sin(w * times))
    }
    full <- function(tt) dance_rhythm_from_coefs(coefs, tt, period, K, trend,
                                                 include_trend = TRUE, t_offset = t_offset)
    base <- function(tt) dance_rhythm_from_coefs(coefs, tt, period, K, trend,
                                                 include_trend = FALSE, t_offset = t_offset)
    switch(component,
      full          = full(times),
      baseline_harm = base(times),
      harmonics     = base(times) - coefs[1],
      trend         = (full(times) - base(times)) - (full(t0) - base(t0)),
      stop("unknown component: ", component))
  }
  fits <- mod$individual_fits
  g <- droplevels(factor(df$group))
  q <- stats::qnorm(1 - (1 - conf) / 2)
  out <- list(); curves <- list()
  simultaneous <- identical(band, "simultaneous") && linear_view
  for (lv in levels(g)) {
    rows <- which(g == lv)
    M <- do.call(rbind, lapply(rows, function(i) {
      f <- fits[[df$subject[i]]]
      if (is.null(f) || !isTRUE(f$success) || is.null(f$coefs)) return(NULL)
      pr <- tryCatch(curve_of(f$coefs, f$t_offset %||% mod$t_offset %||% 0),
                     error = function(e) NULL)
      if (is.null(pr) || !all(is.finite(pr))) NULL else pr
    }))
    if (is.null(M) || nrow(M) < 2) next
    # the group's vector-mean coefficients; with no grouping, the pooled
    # population fit the summary prints -- the same coefficients in both cases
    gf <- mod$group_fits[[lv]]
    mc <- if (!is.null(gf)) gf$mean_coefs
          else if (identical(lv, "(all)")) mod$pop_mean_fit$mean_coefs else NULL
    line <- if (!is.null(mc)) curve_of(mc, mod$t_offset %||% 0) else colMeans(M)
    se <- apply(M, 2, stats::sd) / sqrt(nrow(M))
    n <- nrow(M)
    crit <- q
    if (simultaneous) {
      # the band is centred on x(t)'bbar, the mean of the participants' curves
      line <- colMeans(M)
      crit <- if (n > d_view) sqrt(d_view * (n - 1) / (n - d_view) *
                                   stats::qf(conf, d_view, n - d_view)) else NA_real_
    }
    curves[[lv]] <- M
    out[[lv]] <- data.frame(cell = lv, t = times, fit = line, se = se,
                            lo = line - crit * se, hi = line + crit * se,
                            mean_of_curves = colMeans(M), n = n, crit = crit,
                            stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(ok = TRUE, table = tab, cells = names(out), times = times, conf = conf,
       band = if (simultaneous) "simultaneous" else "pointwise", component = component,
       curves = curves, d = d_view, simultaneous_ok = linear_view,
       note = if (simultaneous) sprintf(paste(
                "Simultaneous %.0f%% band: the mean of the participants' curves +/- c * SE(t),",
                "c from Hotelling's T^2 over the %d coefficient(s) this view uses (c = %s for",
                "n = %s) -- it covers the whole group curve at once, as the mixed model's",
                "Scheffe band does."), 100 * conf, d_view,
                paste(unique(format(signif(tab$crit, 3))), collapse = "/"),
                paste(unique(tab$n), collapse = "/"))
              else paste0(sprintf(paste("Pointwise %.0f%% band: the group's mean-coefficient curve +/- z * SE(t),",
                                  "with SE(t) the standard error across participants of their own",
                                  "fitted curves. Not simultaneous."), 100 * conf),
                          if (identical(band, "simultaneous") && !linear_view) paste(
                            " A simultaneous band is not available for this view: with a saturating",
                            "trend every participant has their own tau, so the curves are not a",
                            "linear combination of shared coefficients.") else ""))
}


# ==============================================================================
# THE TWO-STAGE FITS AS THE APP RUNS THEM
# ==============================================================================
# The per-participant fits, the row each contributes to the parameter table,
# the rule for which fits enter the summaries, the bootstrap and the nested
# model set. They used to be written out inside the module's Run handler, so
# the exported script could only refit each participant and stop there; they
# are functions now, called by the app and written verbatim into the script.

# The parameter table's columns, for an analysis with no successful fit.
# THE DIFFERENCE OF TWO GROUPS' CURVES, A - B, from dance_ts_group_curves()
# output (the same component). Pointwise: the two SEs added in quadrature, z.
# Simultaneous (curves linear in the coefficients only): the two-sample
# Hotelling T^2 band, pooled covariance, c^2 = d (N - 2) / (N - d - 1)
# F(d, N - d - 1), SE(t) = sqrt(pooled var(t) (1/nA + 1/nB)).
dance_ts_diff_curve <- function(gc, a, b, band = "pointwise", conf = 0.95) {
  if (is.null(gc) || !all(c(a, b) %in% gc$cells))
    return(list(ok = FALSE, message = "Both groups need at least two converged fits."))
  MA <- gc$curves[[a]]; MB <- gc$curves[[b]]
  nA <- nrow(MA); nB <- nrow(MB); N <- nA + nB; d <- gc$d
  diff <- colMeans(MA) - colMeans(MB)
  sim <- identical(band, "simultaneous") && isTRUE(gc$simultaneous_ok)
  if (sim) {
    vp <- ((nA - 1) * apply(MA, 2, stats::var) + (nB - 1) * apply(MB, 2, stats::var)) / (N - 2)
    se <- sqrt(vp * (1 / nA + 1 / nB))
    crit <- if (N - d - 1 > 0) sqrt(d * (N - 2) / (N - d - 1) * stats::qf(conf, d, N - d - 1)) else NA_real_
  } else {
    se <- sqrt(apply(MA, 2, stats::var) / nA + apply(MB, 2, stats::var) / nB)
    crit <- stats::qnorm(1 - (1 - conf) / 2)
  }
  lo <- diff - crit * se; hi <- diff + crit * se
  list(ok = TRUE, cell1 = a, cell2 = b, band = if (sim) "simultaneous" else "pointwise",
       table = data.frame(t = gc$times, diff = diff, se = se, lo = lo, hi = hi,
                          excludes_zero = is.finite(lo) & (lo > 0 | hi < 0)),
       note = if (sim) sprintf(paste(
                "Difference of the two groups' mean curves with a SIMULTANEOUS %.0f%% band",
                "(two-sample Hotelling T^2 over %d coefficient(s), pooled covariance; c = %s):",
                "it can be read across the whole curve as one statement."),
                100 * conf, d, format(signif(crit, 3)))
              else paste0(sprintf(paste("Difference of the two groups' mean curves with a POINTWISE",
                                 "%.0f%% band from the two standard errors added in quadrature",
                                 "(independent groups). Not simultaneous: it cannot be read across",
                                 "the whole curve as one statement."), 100 * conf),
                         if (identical(band, "simultaneous") && !isTRUE(gc$simultaneous_ok))
                           " A simultaneous band is not available with a saturating trend." else ""))
}

dance_ts_param_cols <- function(trend_type, n_harmonics) {
  cols <- c("subject", "mesor", "mesor_se", "mesor_adj", "value_at_start", "trend_change")
  cols <- c(cols, switch(as.character(trend_type),
                         linear = c("trend_linear", "trend_linear_se"),
                         log = c("trend_log", "trend_log_se"),
                         exp_sat = c("A_sat", "A_sat_se", "tau", "tau_se"),
                         character(0)))
  for (h in seq_len(n_harmonics))
    cols <- c(cols, paste0("amplitude_", h), paste0("amp_se_", h),
              paste0("acrophase_rad_", h), paste0("acrophase_time_", h),
              paste0("acro_se_time_", h), paste0("beta_cos_", h), paste0("beta_sin_", h))
  c(cols, "r_squared", "percent_rhythm", "p_value", "r_squared_S", "r_squared_C",
    "unique_S", "unique_C", "shared_SC", "percent_S", "percent_C", "percent_shared",
    "converged", "boundary_hit")
}

# One participant's row of the parameter table, from their own fit: the
# intercept, the MESOR (the midline averaged over the window), the value at
# the first time point, the non-periodic change, the trend parameters (and
# A/tau), every harmonic's amplitude and acrophase, and the fit statistics.
# `time_vec` is the model axis (t = 0 at the recording's first time point).
dance_ts_param_row <- function(i, fit_i, time_vec, period, n_harmonics, trend_type) {
  # the trend's coefficients sit after the constant (dance_trend_npar() is the
  # one count of a trend's parameters)
  coef_offset <- 1L + dance_trend_npar(trend_type)
  row_data <- list(subject = i, mesor = fit_i$mesor, mesor_se = fit_i$mesor_se)

  # the rhythm-adjusted mean and the fitted value at the first
  # observation, per subject, from this subject's own coefficients
  .tc <- switch(as.character(trend_type),
                "linear"  = c(fit_i$trend_params$trend_linear$coef %||% NA_real_),
                "log"     = c(fit_i$trend_params$trend_log$coef %||% NA_real_),
                "exp_sat" = c(fit_i$trend_params$A_sat$coef %||% NA_real_,
                              fit_i$trend_params$tau$coef %||% NA_real_),
                numeric(0))
  row_data$mesor_adj <- dance_rhythm_adjusted_mean(
    fit_i$mesor, trend_type, .tc,
    min(time_vec, na.rm = TRUE), max(time_vec, na.rm = TRUE),
    fit_i$t_offset %||% 0)
  row_data$value_at_start <- tryCatch(
    as.numeric(dance_rhythm_from_coefs(
      fit_i$coefs, min(time_vec, na.rm = TRUE), period, n_harmonics,
      trend_type, include_trend = TRUE, t_offset = fit_i$t_offset %||% 0)),
    error = function(e) NA_real_)
  # the non-periodic change over the recording: the trend at the last
  # observation minus the trend at the first, from this subject's own
  # trend parameters -- for the saturating trend, A and tau together
  row_data$trend_change <- if (identical(as.character(trend_type), "none")) 0 else
    tryCatch({
      .s <- dance_trend_value(trend_type, .tc, range(time_vec, na.rm = TRUE),
                              fit_i$t_offset %||% 0)
      if (all(is.finite(.s))) .s[2] - .s[1] else NA_real_
    }, error = function(e) NA_real_)

  # Add trend parameters based on type
  if(trend_type != "none" && !is.null(fit_i$trend_params)) {
    for(param_name in names(fit_i$trend_params)) {
      row_data[[param_name]] <- fit_i$trend_params[[param_name]]$coef
      row_data[[paste0(param_name, "_se")]] <- fit_i$trend_params[[param_name]]$se
    }
    # the initial slope of a saturating rise, A/tau: how fast it starts
    if (identical(as.character(trend_type), "exp_sat"))
      row_data$init_slope <- (fit_i$trend_params$A_sat$coef %||% NA_real_) /
                             (fit_i$trend_params$tau$coef %||% NA_real_)
  }

  for(h in 1:n_harmonics) {
    row_data[[paste0("amplitude_", h)]] <- fit_i$amplitudes[h]
    row_data[[paste0("amp_se_", h)]] <- fit_i$amp_se[h]
    row_data[[paste0("acrophase_rad_", h)]] <- fit_i$acrophases[h]
    row_data[[paste0("acrophase_time_", h)]] <- fit_i$acrophases_time[h]
    row_data[[paste0("acro_se_time_", h)]] <- fit_i$acro_se_time[h]

    # Get beta coefficients - handle both linear and nls fits
    if(trend_type == "exp_sat") {
      row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[paste0("b_cos", h)]
      row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[paste0("b_sin", h)]
    } else {
      cos_idx <- coef_offset + 2 * (h - 1) + 1
      sin_idx <- coef_offset + 2 * (h - 1) + 2
      row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[cos_idx]
      row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[sin_idx]
    }
  }
  row_data$r_squared <- fit_i$r_squared
  row_data$percent_rhythm <- fit_i$percent_rhythm
  row_data$p_value <- fit_i$p_value
  row_data$aic <- fit_i$aic
  row_data$aicc <- fit_i$aicc
  row_data$bic <- fit_i$bic
  row_data$loocv_rmse <- fit_i$loocv_rmse
  row_data$r_squared_S <- fit_i$r_squared_S
  row_data$r_squared_C <- fit_i$r_squared_C
  row_data$unique_S <- fit_i$unique_S %||% NA_real_
  row_data$unique_C <- fit_i$unique_C %||% NA_real_
  row_data$shared_SC <- fit_i$shared_SC %||% NA_real_
  row_data$percent_S <- fit_i$percent_S
  row_data$percent_C <- fit_i$percent_C
  row_data$percent_shared <- fit_i$percent_shared %||% NA_real_
  row_data$converged <- isTRUE(fit_i$converged)
  row_data$boundary_hit <- isTRUE(fit_i$boundary_hit)

  as.data.frame(row_data)
}

# Every participant's cosinor, on one time axis with one trend origin -- the
# recording's first time point -- and one tau setting. Returns the fits (one
# per row of Y; a failed fit records why), the parameter table of the fits
# that returned, and the failures.
#   progress   optional function(i, n)
dance_ts_fit_all <- function(Y, time_vec, period = 24, n_harmonics = 1, trend_type = "none",
                             use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                             amplitude_min = 0, amplitude_max = NA,
                             A_sat_min = NA, A_sat_max = NA,
                             tau_min = 0.5, tau_max = NA, tau_fixed = NA, progress = NULL) {
  n_trend_params <- dance_trend_npar(trend_type)
  fits <- list(); rows <- list(); failed <- list()
  for (i in seq_len(nrow(Y))) {
    y_i <- Y[i, ]
    n_valid <- sum(!is.na(y_i))
    fit_i <- fit_cosinor(time_vec, y_i, period = period, n_harmonics = n_harmonics,
                         trend_type = trend_type, use_bounds = use_bounds,
                         mesor_min = mesor_min, mesor_max = mesor_max,
                         amplitude_min = amplitude_min, amplitude_max = amplitude_max,
                         A_sat_min = A_sat_min, A_sat_max = A_sat_max,
                         tau_min = tau_min, tau_max = tau_max, tau_fixed = tau_fixed,
                         t_origin = min(time_vec, na.rm = TRUE))
    if (isTRUE(fit_i$success)) {
      fits[[i]] <- fit_i
      rows[[length(rows) + 1L]] <- dance_ts_param_row(i, fit_i, time_vec, period,
                                                      n_harmonics, trend_type)
    } else {
      fits[[i]] <- list(success = FALSE, message = fit_i$message, n_valid = n_valid,
                        n_required = 2 * n_harmonics + 1 + n_trend_params + 1)
      failed[[length(failed) + 1L]] <- list(subject = i, n_valid = n_valid,
                                            reason = fit_i$message)
    }
    if (is.function(progress) && i %% 10 == 0) progress(i, nrow(Y))
  }
  params <- if (length(rows)) do.call(rbind, rows) else {
    cols <- dance_ts_param_cols(trend_type, n_harmonics)
    stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
  }
  list(fits = fits, params = params, failed = failed)
}

# WHO ENTERS THE SUMMARIES. Fits that did not converge never do: there is no
# solution to average. Fits sitting on a parameter bound do when
# include_boundary is TRUE (the default, with the bounds tabled) -- the value
# is real but it is where the optimiser was stopped. If fewer than 3 fits pass,
# every returned fit is used and the audit says so (gate_relaxed).
dance_ts_gate <- function(all_params, fits, n_attempted, include_boundary = TRUE,
                          subject_ids = NULL) {
  n <- nrow(all_params)
  conv_flag  <- if ("converged" %in% names(all_params)) all_params$converged else rep(TRUE, n)
  bound_flag <- if ("boundary_hit" %in% names(all_params)) all_params$boundary_hit else rep(FALSE, n)
  conv_flag[is.na(conv_flag)] <- FALSE
  bound_flag[is.na(bound_flag)] <- FALSE
  # which bounds each returned fit sits on, per parameter
  bounds_list <- lapply(seq_len(n), function(i) {
    f <- fits[[all_params$subject[i]]]
    if (is.null(f) || is.null(f$bounds_hit)) character(0) else f$bounds_hit
  })
  ids <- function(rows) if (!is.null(subject_ids)) subject_ids[all_params$subject[rows]]
                        else all_params$subject[rows]
  fit_audit <- list(
    n_attempted = n_attempted, n_returned = n,
    n_converged = sum(conv_flag & !bound_flag), n_boundary = sum(conv_flag & bound_flag),
    n_failed = n_attempted - n, n_nonconverged = sum(!conv_flag),
    bounds = dance_bounds_summary(bounds_list, subject_ids = ids(seq_len(n))))
  keep <- conv_flag & (include_boundary | !bound_flag)
  fit_audit$include_boundary <- include_boundary
  fit_audit$n_passing <- sum(keep)
  fit_audit$gate_relaxed <- sum(keep) < 3
  if (fit_audit$gate_relaxed) keep <- rep(TRUE, n)
  fit_audit$n_summarised <- sum(keep)
  fit_audit$bounds_kept <- dance_bounds_summary(bounds_list[keep], subject_ids = ids(which(keep)))
  list(keep = keep, fit_audit = fit_audit)
}

# THE BOOTSTRAP of the pooled sample's intercept and first-harmonic amplitude
# and acrophase. The resampling unit is the PARTICIPANT, drawn with replacement
# by index (a participant drawn three times counts three times); intervals are
# percentile, the acrophase's circular. A FIXED SEED, so the same data give the
# same interval twice and the seed can be stated; the caller's random stream is
# restored afterwards.
#   progress   optional function(b, B)
dance_ts_bootstrap <- function(params, B, period = 24, seed = 20240601L, progress = NULL) {
  rows_by_subject <- split(seq_len(nrow(params)), params$subject)
  subj_ids <- names(rows_by_subject)
  n_boot_subj <- length(subj_ids)
  boot_mesor <- numeric(B); boot_amplitude <- numeric(B); boot_acrophase <- numeric(B)
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  for (b in seq_len(B)) {
    take <- dance_boot_index(n_boot_subj)
    bp <- params[unlist(rows_by_subject[subj_ids[take]], use.names = FALSE), , drop = FALSE]
    boot_mesor[b] <- mean(bp$mesor, na.rm = TRUE)
    x_b <- bp$amplitude_1 * cos(bp$acrophase_rad_1)
    y_b <- bp$amplitude_1 * sin(bp$acrophase_rad_1)
    boot_amplitude[b] <- sqrt(mean(x_b, na.rm = TRUE)^2 + mean(y_b, na.rm = TRUE)^2)
    acro_b <- atan2(mean(y_b, na.rm = TRUE), mean(x_b, na.rm = TRUE))
    if (is.finite(acro_b) && acro_b < 0) acro_b <- acro_b + 2 * pi
    boot_acrophase[b] <- acro_b
    if (is.function(progress) && b %% 50 == 0) progress(b, B)
  }
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
  acro_ci <- dance_boot_circ_ci_time(boot_acrophase, period, 1)
  list(
    mesor_ci = stats::quantile(boot_mesor, c(0.025, 0.975), na.rm = TRUE),
    amplitude_ci = stats::quantile(boot_amplitude, c(0.025, 0.975), na.rm = TRUE),
    # circular, and carrying its own width because lo > hi when it wraps
    acrophase_ci = c(`2.5%` = acro_ci$lo, `97.5%` = acro_ci$hi),
    acrophase_ci_circular = acro_ci,
    # what a linear quantile would have reported, for the readout to contrast
    # when the interval wraps
    acrophase_ci_linear = stats::quantile(phi_to_hours(boot_acrophase, period, 1),
                                          c(0.025, 0.975), na.rm = TRUE),
    resample_unit = "participant", n_resampled = n_boot_subj,
    seed = seed, ci_method = "percentile",
    boot_mesor = boot_mesor, boot_amplitude = boot_amplitude,
    boot_acrophase = boot_acrophase, B = B)
}

# THE NESTED MODEL SET: every trend crossed with the cumulative harmonic sets up
# to n_harmonics, each fitted per participant with this run's tau setting and
# origin, its AICc averaged over the participants scorable under EVERY
# candidate. Returns the Akaike table (NULL if no participant is scorable under
# all of them).
#   progress   optional function(j, n, label)
dance_ts_model_set <- function(Y, time_vec, period = 24, n_harmonics = 1, trend_type = "none",
                               tau_min = 0.5, tau_max = NA, tau_fixed = NA, progress = NULL) {
  ms_models <- dance_cosinor_model_set(n_harmonics)
  n_time <- ncol(Y)
  # participants x candidates; NA where the candidate cannot be fitted
  aicc <- matrix(NA_real_, nrow(Y), length(ms_models),
                 dimnames = list(NULL, vapply(ms_models, function(m) m$label, character(1))))
  for (j in seq_along(ms_models)) {
    m <- ms_models[[j]]
    if (is.function(progress)) progress(j, length(ms_models),
                                        sprintf("%s, %s", dance_trend_label(m$trend),
                                                dance_harmonic_label(m$k)))
    # parameters + 1: a fit needs at least one observation more than it has
    # parameters for the residual variance to be estimable
    need <- dance_model_npar(m$trend, m$k) + 1L
    if (n_time < need) next
    for (i in seq_len(nrow(Y))) {
      y_i <- as.numeric(Y[i, ])
      if (sum(!is.na(y_i)) < need) next
      # the saturating candidates use the same tau setting as the fit being
      # reported, so its row in the table is that fit
      f <- tryCatch(fit_cosinor(time_vec, y_i, period = period, n_harmonics = m$k,
                                trend_type = m$trend, tau_min = tau_min, tau_max = tau_max,
                                tau_fixed = if (identical(m$trend, "exp_sat")) tau_fixed else NA,
                                t_origin = min(time_vec, na.rm = TRUE)),
                    error = function(e) NULL)
      if (is.null(f) || !isTRUE(f$success) || !is.finite(f$aicc)) next
      aicc[i, j] <- f$aicc
    }
  }
  common <- dance_model_set_common(aicc)
  tab <- if (common$n_common > 0) dance_akaike_table(common$means) else NULL
  if (!is.null(tab)) {
    attr(tab, "per_subject_mean") <- TRUE
    attr(tab, "n_common") <- common$n_common
    attr(tab, "n_participants") <- nrow(Y)
    attr(tab, "not_fitted") <- common$not_fitted
    # the specification the run reported, so a reader need not match it by eye
    attr(tab, "selected") <- dance_model_label(trend_type, n_harmonics)
    attr(tab, "n_harmonics_max") <- n_harmonics
  }
  tab
}
