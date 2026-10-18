# ==============================================================================
# server/08j_helpers_assumptions.R -- WHAT A COSINOR RUN ASSUMED
# ==============================================================================
# Every choice a cosinor run makes that is not a setting on screen: a default,
# a convention, a rule of the estimator. They used to live in code comments and
# a few notifications, which is not where a reader of the results looks, and
# some of them decide the numbers (the time axis, the residual model, which
# participants enter a mean). This lists them for ONE run, from the fitted
# objects rather than from the inputs, so a setting changed after the run does
# not change what is reported. The summary, its download and the report all
# print this one list.
#
# Returns a data frame (area, text), one row per assumption, in reading order.
# Pure: no reactives, so it is testable without a session.
# ==============================================================================

dance_cosinor_assumptions <- function(mod, ff = NULL, df_method = "kr") {
  if (is.null(mod)) return(data.frame(area = character(0), text = character(0)))
  rows <- list()
  a <- function(area, ...) rows[[length(rows) + 1L]] <<- data.frame(
    area = area, text = paste0(...), stringsAsFactors = FALSE)
  f <- function(x, d = 3) format(signif(as.numeric(x), d))
  # a mixed run whose model did not fit has no model assumptions to list --
  # only the shared ones about the time axis and the data
  is_mixed <- identical(mod$approach %||% "two_stage", "mixed")
  mixed <- is_mixed && isTRUE(ff$ok)
  P <- mod$period %||% 24
  K <- as.integer(mod$n_harmonics %||% 1L)
  trend <- as.character(mod$trend_type %||% "none")
  tv <- mod$time_vec
  n_t <- length(tv)
  span <- if (n_t) diff(range(tv, na.rm = TRUE)) else NA_real_
  clock0 <- dance_clock_label(as.numeric(mod$origin_shift %||% 0), P, show_day = FALSE)

  # ---- the time axis -----------------------------------------------------------
  src <- mod$time_source %||% ""
  if (identical(src, "index"))
    a("Time", sprintf(paste(
      "NOT read from the data: the %d columns were taken as equally spaced over exactly",
      "one period, every %s h from 0. If the columns do not cover one period at that",
      "spacing, every time-dependent estimate is wrong; choose the clock times instead."),
      n_t, f(P / max(1L, n_t))))
  else a("Time", switch(sub(":.*$", "", src),
    shared = "the clock times parsed from the column names at import.",
    manual = "the times entered by hand.",
    covariate = sprintf("read from the column '%s'.", sub("^covariate:", "", src)),
    "as stored with the model."))
  if (!identical(src, "index"))
    a("Time", if (isTRUE(mod$wrap_applied))
      paste("The columns are taken to be in recording order: each time earlier than the one",
            "before it was read as a midnight crossing and moved 24 h later, with every time",
            "after it.")
      else "The times already ran forward, so no midnight crossing was assumed.")
  a("Time", sprintf(paste(
    "t = 0 is the first time point of the recording (%s); every trend and every harmonic",
    "is anchored there, for every participant, so the intercept is the midline at that",
    "instant -- also for a participant whose first observation is missing."), clock0))
  a("Time", sprintf(paste(
    "The period is fixed at %s h, not estimated; amplitudes and acrophases are",
    "conditional on it."), f(P)))
  a("Time", if (K > 1L) sprintf(paste(
    "The harmonics are cumulative: periods %s h (harmonics 1 to %d), none skipped."),
    paste(f(P / seq_len(K)), collapse = ", "), K)
    else sprintf("One harmonic: the fundamental, period %s h.", f(P)))
  a("Data", paste("Raw observations; a missing observation is left out, not imputed or",
                  "smoothed, and each fit uses the observations present."))

  # ---- the trend -----------------------------------------------------------------
  if (identical(trend, "log"))
    a("Trend", paste("log(1 + t), t in hours: the +1 h sets the curve's shape near t = 0 and",
                     "is a convention, not an estimate."))
  if (!identical(trend, "none"))
    a("Trend", sprintf(paste(
      "MESOR = the midline (intercept + trend) averaged over t = 0 to %s h, the first to",
      "the last time point; change = the trend at the last time point (it is 0 at t = 0)."),
      f(span)))

  if (mixed) {
    sp <- ff$spec
    # ---- the mixed-effects model -------------------------------------------------
    a("Model", paste("One linear mixed model over every observation; each design cell has its",
                     "own intercept, trend and harmonics (the basis crossed with every design",
                     "factor and their interactions)."))
    if (length(sp$covariates))
      a("Model", sprintf(paste("Covariate(s) %s enter additively, as a shift of the level",
                               "only, and are held at their mean in the cell estimates."),
                         paste(sp$covariates, collapse = ", ")))
    a("Estimation", sprintf(paste(
      "%s, lme4 through lmerTest, optimiser bobyqa (up to 200,000 evaluations)."),
      if (isTRUE(ff$REML)) "REML" else "Maximum likelihood"))
    top <- sp$re_ladder[[1]]$label %||% ff$re_label
    a("Random effects", sprintf(paste(
      "Normal, with a full covariance within each term. The structure starts at '%s' and",
      "moves down only when a fit fails or does not converge; a singular fit is kept.",
      "Used: rung %d of %d, '%s'."), top, as.integer(ff$re_rung %||% 1L),
      as.integer(ff$n_rungs %||% length(sp$re_ladder)), ff$re_label))
    ar <- ff$residual_ar
    if (isTRUE(ff$phi > 0) || !is.null(ar)) {
      a("Residuals", sprintf(paste(
        "Normal, with ONE variance for every cell and time point, and AUTOCORRELATED within",
        "each curve as a continuous-time AR(1): two residuals dt hours apart correlate",
        "phi^dt (uneven gaps and missing observations as they are), and residuals of",
        "different curves are independent. phi = %s per hour, estimated by profile",
        "likelihood (%s) on 20 values from 0 to 0.98 evenly spaced in -log(1 - phi), the",
        "maximum refined between neighbours, its 95%% interval where the log-likelihood",
        "falls 1.92; the model is fitted to the whitened data (Prais-Winsten for unequal",
        "gaps), which is the same model, and refitted by REML at phi-hat. Its tests treat",
        "phi as known. Independence is tested by the likelihood ratio, half the",
        "chi-square(1) p (phi = 0 is on the boundary)."),
        f(ar$phi %||% ff$phi), ar$method %||% "REML"))
      sf <- ff$sat
      if (!is.null(sf) && isTRUE(sf$ok))
        a("Residuals", sprintf(paste(
          "The per-cell nonlinear model has the same CAR(1) residuals (nlme corCAR1),",
          "with phi estimated within it: %s per hour."), f(sf$phi)))
    } else a("Residuals", paste(
      "Normal and INDEPENDENT, with ONE variance for every cell and time point: no",
      "autocorrelation between successive observations and no unequal variance is",
      "modelled. Serially correlated residuals make standard errors too small; the",
      "summary's diagnostics report the correlation of consecutive residuals, and",
      "Model Specification offers a continuous-time AR(1)."))
    a("Tests", sprintf(paste(
      "F tests of linear contrasts of the fixed effects with %s denominator df. A main",
      "effect averages the other factors' levels with EQUAL weights, whatever the cell",
      "sizes; pairwise tables are adjusted within each table (Holm unless changed on the",
      "comparison tab)."),
      switch(df_method, kr = "Kenward-Roger", satterthwaite = "Satterthwaite", df_method)))
    a("Intervals", paste(
      "Per-cell intervals for the intercept, MESOR and change: Wald, normal quantile.",
      "Amplitude and acrophase: simulated from the joint normal distribution of each",
      "(cos, sin) pair -- 6,000 draws per cell, 20,000 per pairwise contrast, seed",
      "20240601; an acrophase is undefined when the pair's joint 95% region contains",
      "the origin (Bingham et al., 1982). Curve bands:",
      "pointwise (normal quantile) or simultaneous (Scheffe, on the smallest",
      "Satterthwaite df), as chosen on the comparison tab. 95% throughout."))

    if (identical(trend, "exp_sat")) {
      pr <- ff$tau_profile
      if (identical(ff$tau_mode, "fixed")) {
        a("tau", sprintf(paste("Held at %s h for every cell; the rise A is estimated per cell",
                               "given that value."), f(sp$tau)))
      } else {
        rg <- ff$tau_range %||% range(pr$grid %||% NA_real_)
        a("tau", sprintf(paste(
          "Estimated ONCE for all cells by profile likelihood: maximum-likelihood fits at",
          "%d values spaced evenly in log tau from %s to %s h, all on the random structure",
          "chosen at the middle of the grid (rung %s); the maximum refined between its",
          "neighbouring grid points (to about 1%%); the 95%% interval where the",
          "log-likelihood falls 1.92 below its maximum, interpolated in log tau. An",
          "estimate at an end of the range is a bound. The model is then refitted by REML",
          "at that tau, and its tests treat it as known."),
          length(pr$grid %||% integer(0)), f(rg[1]), f(rg[2]), pr$rung %||% "?"))
        a("tau", sprintf(paste(
          "A profile that varies by less than 1.92 over the whole grid is flat: tau is then",
          "not identified and is held at the held value, %s h%s."),
          f(mod$tau_spec$value %||% sp$tau),
          if (identical(ff$tau_mode, "fixed_after_flat")) " -- which is what happened in this run" else ""))
      }
      sf <- ff$sat
      if (!is.null(sf) && isTRUE(sf$ok))
        a("Per-cell A, tau", sprintf(paste(
          "Nonlinear mixed model (nlme), maximum likelihood, every parameter per cell; tau",
          "kept inside %s to %s h as lo + (hi - lo) x logistic(theta); diagonal",
          "random-effects covariance (used: %s); covariance scaled by N/(N - p) (nlme's",
          "adjustSigma); Wald intervals and tests on %d df (participants minus",
          "between-participant cells), tau on the log scale, A/tau by the delta method;",
          "one tau against a tau per cell by a likelihood-ratio test on the same random",
          "structure. Flagged: |corr(A, log tau)| > 0.9, tau beyond the recording, tau",
          "within 1%% of an end of the range."),
          f(sf$tau_range[1]), f(sf$tau_range[2]), sf$re_label, as.integer(sf$df)))
    }
  } else if (!is_mixed) {
    # ---- two-stage -----------------------------------------------------------------
    st <- mod$dance_settings %||% list()
    held <- identical(trend, "exp_sat") && identical(mod$tau_spec$mode, "fixed")
    a("Per-participant fit", if (identical(trend, "exp_sat")) {
      if (held) sprintf(paste(
        "tau held at %s h for every participant: each fit is then linear in its other",
        "parameters, fitted by least squares with tau as a constant."), f(mod$tau_spec$value))
      else sprintf(paste(
        "Nonlinear least squares: Levenberg-Marquardt (minpack.lm, up to 300 iterations),",
        "then nls 'port', then four restarts with tau started at T x 0.1/3, 0.5/3, 2/3",
        "and 5/3 (T = the participant's last observation time, at most %s h); tau kept",
        "inside %s to %s h in every attempt. Starts: A = 1.5 x the data range and tau =",
        "T/4 for a rising series, A = the linear slope x T and tau = T/3 for a falling",
        "one, moved inside the range. A restart is accepted only if it converged with",
        "R-squared > 0."),
        f(span), f(st$tau_min), f(st$tau_max))
    } else "Ordinary least squares, closed form (the model is linear in every parameter).")
    if (isTRUE(st$use_bounds))
      a("Per-participant fit", paste(
        "Parameter bounds are on: a blank bound defaults to the participant's data range",
        "(intercept within [min - range, max + range], amplitude and A at most 2 x range);",
        "an amplitude above its maximum is refitted inside the inscribed square; the",
        "closed-form fit is kept whenever it violates no bound."))
    fa <- mod$fit_audit
    a("Inclusion", sprintf(paste(
      "Fits that did not converge are excluded; fits within 1e-6 (relative) of a",
      "parameter bound are %s.%s"),
      if (isTRUE(fa$include_boundary %||% TRUE)) "INCLUDED (Advanced)" else "excluded (Advanced)",
      if (isTRUE(fa$gate_relaxed))
        " Fewer than 3 fits passed, so ALL returned fits were used in this run." else
        " If fewer than 3 fits pass, all returned fits are used and the summary says so."))
    a("Group summaries", paste(
      "Intercept, trend, MESOR and change: arithmetic means of the participants' own",
      "values, with t intervals. Amplitude and acrophase: amplitude-weighted vector means",
      "of the (cos, sin) pairs (the population-mean cosinor). A group with fewer than 3",
      "participants is not summarised; a participant without a group label is left out",
      "of the group analyses and kept in the pooled one."))
    a("Tests", paste(
      "Every participant's estimate is treated as exact (no weighting by its precision),",
      "and the participants as independent: a within-participant factor cannot be",
      "expressed. Omnibus: one-way MANOVA, Wilks' lambda, on the coefficient vector",
      "(multivariate normal, equal covariance matrices). Components: one-way ANOVA",
      "(equal variances). Amplitude and acrophase: Bingham's population-mean cosinor;",
      "acrophase also by Watson-Williams (von Mises, common concentration; flagged when",
      "r-bar < 0.45). Pairwise: Welch's t (unequal variances) and Watson-Williams,",
      "adjusted within each table (Holm unless changed)."))
    a("Rhythm detection", paste(
      "Per participant, the F test of all harmonics given the trend; 'detected' means",
      "p < .05, with no adjustment across participants."))
    a("Cross-validation", paste(
      "Leave one observation out and refit, per participant (a nonlinear refit starts",
      "from the full-data estimates; one that fails is left out of the average)."))
    br <- mod$boot_results
    if (!is.null(br))
      a("Bootstrap", sprintf(paste(
        "Participants resampled with replacement, B = %d; percentile 95%% intervals",
        "(the acrophase's circular); first harmonic and the pooled sample only; seed %s."),
        as.integer(br$B), format(br$seed %||% "not recorded")))
    ms <- mod$model_selection
    if (!is.null(ms))
      a("Model comparison", sprintf(paste(
        "Candidates: every trend crossed with the cumulative harmonic sets up to %d, each",
        "fitted per participant with this run's tau setting; AICc averaged over the %s",
        "participant(s) scorable under every candidate; Akaike weights over this set only."),
        as.integer(attr(ms, "n_harmonics_max") %||% K),
        format(attr(ms, "n_common") %||% "same")))
    cn <- mod$conditioning
    if (!is.null(cn$tau_fixed_delta_aic))
      a("Model comparison", sprintf(paste(
        "Free tau against tau held at %s h: AIC(free) - AIC(held) per participant, on",
        "the same observations, averaged over the %d participant(s) scorable in both."),
        f(cn$tau_fixed_value), as.integer(cn$tau_fixed_n)))
  }
  do.call(rbind, rows)
}

# The list as lines of text for the summary and its download: the area once,
# then each assumption wrapped under it with a hanging indent.
dance_cosinor_assumptions_text <- function(tab, label_width = 18, width = 92) {
  if (is.null(tab) || !nrow(tab)) return(character(0))
  lab <- ifelse(c(TRUE, tab$area[-1] != tab$area[-nrow(tab)]), tab$area, "")
  pad <- strrep(" ", label_width + 1L)
  unlist(lapply(seq_len(nrow(tab)), function(i) {
    w <- strwrap(tab$text[i], width = width - label_width - 1L)
    c(sprintf(paste0("%-", label_width, "s %s"), lab[i], w[1]),
      if (length(w) > 1L) paste0(pad, w[-1]))
  }))
}
