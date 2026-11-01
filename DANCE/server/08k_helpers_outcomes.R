# ==============================================================================
# server/08k_helpers_outcomes.R -- the group comparison, organised by question
#
# Three scientific questions, and above them one omnibus question:
#
#   0. Do the groups differ AT ALL?            the complete fitted trajectory
#   1. Overall response level                  time-averaged fitted response,
#                                              mean midline (MESOR), initial midline
#   2. Circadian-like rhythmicity              per harmonic: the joint (cos, sin)
#                                              difference, amplitude, acrophase
#   3. Homeostatic-like accumulation           total trend change, average rate,
#                                              initial rate; then the
#                                              model-specific parameters (the
#                                              trend coefficient, A, tau, t50)
#
# Two tables per fit, the same shape for both estimators:
#   estimates  one row per (outcome, group): estimate and interval
#   tests      one row per test: statistic, df, p, method, and what it tests
# The comparison tab, the CSV export and the APA report all read these, so a
# number cannot be one thing on screen and another in the report.
#
# Nothing here is a new statistical procedure. The mixed model's rows are its
# existing linear functionals (Wald intervals, the trajectory model's F tests)
# and the per-cell nonlinear model's delta-method rows; the two-stage rows are
# t intervals and the existing ANOVA / MANOVA / Bingham / Watson-Williams
# tests. Bingham's population-mean confidence region for a group's amplitude
# and acrophase is computed by the same ellipse routine the per-participant fits
# use (dance_bingham_ci), with the population-mean covariance.
#
# "Primary" marks the outcome the domain is read from first; the others in a
# domain are supplementary, and the model-specific parameters are flagged as
# such -- they describe the model, not the question.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

DANCE_OUTCOME_DOMAINS <- c(omnibus = "Complete fitted trajectory (omnibus)",
                           level = "Overall response level",
                           rhythm = "Circadian-like rhythmicity",
                           accumulation = "Homeostatic-like accumulation",
                           accumulation_model = "Model-specific trend parameters")

dance_outcome_est_row <- function(domain, outcome, label, cell, estimate, lo = NA_real_,
                                  hi = NA_real_, primary = FALSE, kind = "linear",
                                  note = NA_character_, harmonic = NA_integer_) {
  data.frame(domain = domain, outcome = outcome, label = label, cell = as.character(cell),
             estimate = as.numeric(estimate), lo = as.numeric(lo), hi = as.numeric(hi),
             primary = primary, kind = kind, harmonic = as.integer(harmonic),
             note = note, stringsAsFactors = FALSE)
}

dance_outcome_test_row <- function(domain, outcome, label, statistic, df1, df2, p, method,
                                   primary = FALSE, note = NA_character_) {
  data.frame(domain = domain, outcome = outcome, label = label,
             statistic = as.numeric(statistic %||% NA_real_),
             df1 = as.numeric(df1 %||% NA_real_), df2 = as.numeric(df2 %||% NA_real_),
             p = as.numeric(p %||% NA_real_), method = as.character(method %||% ""),
             primary = primary, note = as.character(note %||% NA_character_),
             stringsAsFactors = FALSE)
}

dance_outcome_bind <- function(rows) {
  rows <- Filter(function(r) !is.null(r) && nrow(r) > 0, rows)
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# MIXED-EFFECTS
# ------------------------------------------------------------------------------
#   ff       the trajectory fit (dance_traj_fit / dance_traj_run)
#   tests    dance_traj_tests(ff, ...): the omnibus and block tests, already run
#   coefs    optional: the per-harmonic dance_traj_cell_coefs() results, named
#            H1, H2, ... (the app caches them; recomputed when absent)
dance_outcomes_mixed <- function(ff, tests = NULL, conf = 0.95, n_draw = 6000, coefs = NULL) {
  if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
  spec <- ff$spec
  trend <- spec$trend %||% "none"
  K <- length(spec$harm_terms) %/% 2L
  sat <- isTRUE(ff$sat$ok)
  u <- if (is.null(ff$dv_units) || !nzchar(ff$dv_units)) "response units" else ff$dv_units
  est <- list(); tst <- list()
  E <- function(...) est[[length(est) + 1L]] <<- dance_outcome_est_row(...)
  Tt <- function(...) tst[[length(tst) + 1L]] <<- dance_outcome_test_row(...)

  # ---- per-cell estimates of the level and accumulation quantities ----------
  ls <- dance_traj_level_summary(ff, conf)
  if (isTRUE(ls$ok)) for (i in seq_len(nrow(ls$table))) {
    r <- ls$table[i, ]; q <- r$quantity
    dom <- if (q %in% c("fitted_mean", "mesor", "intercept")) "level" else "accumulation"
    # with tau estimated per cell, the accumulation outcomes come from the
    # per-cell model below, with tau's uncertainty carried through
    if (sat && dom == "accumulation") next
    prim <- (q == "fitted_mean") || (q == "mesor" && K == 0L) || (q == "intercept" && trend == "none" && K == 0L) ||
            (q == "change")
    E(dom, q, DANCE_TRAJ_BLOCK_LABEL[[q]], r$cell, r$estimate, r$lo, r$hi, primary = prim,
      note = if (dom == "accumulation" && identical(trend, "exp_sat"))
        sprintf("given the shared tau (%s h)", format(signif(spec$tau, 4))) else NA_character_)
  }
  # the trend's own coefficient per cell (A, the slope, the log-beta)
  if (!sat && length(spec$trend_terms)) {
    fx <- dance_traj_cell_functionals(ff); bv <- dance_traj_beta(ff)
    have <- names(bv$beta)[!is.na(bv$beta)]
    z <- stats::qnorm(1 - (1 - conf) / 2)
    for (tm in spec$trend_terms) for (i in seq_len(nrow(fx$grid))) {
      rw <- fx$cells[[i]]$coef[[tm]]; nm <- intersect(names(rw), have)
      e <- sum(rw[nm] * bv$beta[nm])
      se <- sqrt(max(0, as.numeric(t(rw[nm]) %*% bv$V[nm, nm, drop = FALSE] %*% rw[nm])))
      E("accumulation_model", "trend", dance_traj_trend_coef_label(trend, spec$tau, ff$dv_units,
                                                                   estimated = isTRUE(ff$tau_estimated)),
        fx$grid$.cell[i], e, e - z * se, e + z * se)
    }
    if (identical(trend, "exp_sat")) for (cl in fx$grid$.cell) {
      E("accumulation_model", "tau", sprintf("Trend time constant tau (hours; %s, common to all cells)",
                                             if (isTRUE(ff$tau_estimated)) "profiled" else "held"),
        cl, spec$tau, kind = "fixed")
      E("accumulation_model", "t50", "Half-rise time tau ln 2 (hours)", cl, spec$tau * log(2), kind = "fixed")
    }
  }
  # the per-cell nonlinear model: tau per cell, everything through the delta method
  if (sat) {
    st <- ff$sat$table
    for (i in seq_len(nrow(st))) {
      E("accumulation", "sat_change", "Total trend change (start to end)", st$cell[i], st$change[i],
        st$change_lo[i], st$change_hi[i], primary = TRUE, kind = "delta",
        note = "tau estimated per cell; its uncertainty carried through")
      E("accumulation", "sat_avg_rate", sprintf("Average trend rate (%s per hour)", u), st$cell[i],
        st$avg_rate[i], st$avg_rate_lo[i], st$avg_rate_hi[i], kind = "delta")
      E("accumulation", "init_slope", sprintf("Initial trend rate A/tau (%s per hour)", u), st$cell[i],
        st$init_slope[i], st$init_slope_lo[i], st$init_slope_hi[i], kind = "delta")
      E("accumulation_model", "A_sat", sprintf("Asymptotic trend change A (%s)", u), st$cell[i],
        st$A[i], st$A_lo[i], st$A_hi[i])
      E("accumulation_model", "tau", "Trend time constant tau (hours)", st$cell[i], st$tau[i],
        st$tau_lo[i], st$tau_hi[i])
      E("accumulation_model", "t50", "Half-rise time tau ln 2 (hours)", st$cell[i], st$t50[i],
        st$tau_lo[i] * log(2), st$tau_hi[i] * log(2))
    }
  }
  # ---- rhythm, per harmonic: amplitude and acrophase -------------------------
  for (h in seq_len(K)) {
    co <- coefs[[sprintf("H%d", h)]] %||% tryCatch(dance_traj_cell_coefs(ff, h), error = function(e) NULL)
    if (is.null(co) || !isTRUE(co$ok)) next
    ap <- tryCatch(dance_traj_amp_phase_ci(co, conf = conf, method = "joint", n_draw = n_draw),
                   error = function(e) NULL)
    if (is.null(ap) || !isTRUE(ap$ok)) next
    for (i in seq_len(nrow(ap$table))) {
      r <- ap$table[i, ]
      E("rhythm", sprintf("amplitude%d", h), sprintf("Amplitude, H%d (%s h)", h, format(spec$period / h)),
        r$cell, r$amplitude, r$amplitude_lo, r$amplitude_hi, harmonic = h)
      E("rhythm", sprintf("phase%d", h), sprintf("Acrophase, H%d (model hours)", h), r$cell,
        r$acrophase_time, if (isTRUE(r$phase_defined)) r$acrophase_lo else NA,
        if (isTRUE(r$phase_defined)) r$acrophase_hi else NA, kind = "circular", harmonic = h,
        note = if (isTRUE(r$phase_defined)) NA_character_ else
          "not identified: the (cos, sin) region includes the origin")
    }
  }

  # ---- the tests ---------------------------------------------------------------
  if (!is.null(tests)) {
    for (nm in names(tests$effects)) {
      r <- tests$effects[[nm]]; if (is.null(r) || !isTRUE(r$ok)) next
      Tt("omnibus", "full", sprintf("Complete fitted trajectory: %s", nm), r$statistic, r$df1, r$df2, r$p,
         r$method, primary = TRUE)
    }
    blk <- function(b) { r <- tests$blocks[[b]]; if (is.null(r) || !isTRUE(r$ok)) NULL else r }
    for (b in c("fitted_mean", "mesor", "intercept")) {
      r <- blk(b); if (is.null(r)) next
      Tt("level", b, DANCE_TRAJ_BLOCK_LABEL[[b]], r$statistic, r$df1, r$df2, r$p, r$method,
         primary = (b == "fitted_mean") || (b == "mesor" && K == 0L))
    }
    r <- blk("circadian")
    if (!is.null(r))
      Tt("rhythm", "circadian", "Rhythmic component, all harmonics jointly (cos, sin)", r$statistic,
         r$df1, r$df2, r$p, r$method, primary = TRUE,
         note = "the overall harmonic-vector difference; neither an amplitude nor an acrophase test alone")
    r <- blk("trend")
    if (!is.null(r) && !sat)
      Tt("accumulation", "change", "Total trend change (start to end)", r$statistic, r$df1, r$df2, r$p,
         r$method, primary = TRUE,
         note = paste0("one trend coefficient per cell: the same test as the average",
                       if (trend %in% c("log", "exp_sat")) " and the initial" else "",
                       " trend rate and the trend coefficient",
                       if (identical(trend, "exp_sat")) " (tau shared)" else ""))
  }
  if (sat) {
    for (q in c("sat_change", "init_slope", "A_sat", "tau")) {
      r <- tryCatch(dance_traj_sat_test(ff$sat, ff, q), error = function(e) NULL)
      if (is.null(r) || !isTRUE(r$ok)) next
      Tt(if (q %in% c("A_sat", "tau")) "accumulation_model" else "accumulation", q,
         dance_traj_sat_label(q, ff$dv_units), r$statistic, r$df1, r$df2, r$p,
         "Wald F, per-cell nonlinear model (delta method)", primary = identical(q, "sat_change"),
         note = if (identical(q, "sat_change")) "the same test as the average trend rate" else
                if (identical(q, "tau")) "tested on the log scale" else NA_character_)
    }
    lr <- ff$sat$lrt
    if (isTRUE(lr$ok))
      Tt("accumulation_model", "tau_lrt", "One tau for all cells vs a tau per cell", lr$statistic,
         lr$df, NA_real_, lr$p, "likelihood-ratio chi-square (ML)")
  }
  list(ok = TRUE, approach = "mixed", estimates = dance_outcome_bind(est),
       tests = dance_outcome_bind(tst), structure = dance_model_structure(trend, K), trend = trend,
       conf = conf, units = ff$dv_units, period = spec$period,
       clock_origin = ff$clock_origin %||% 0)
}

# ------------------------------------------------------------------------------
# TWO-STAGE
# ------------------------------------------------------------------------------
# Bingham's population-mean confidence region for one group: the mean
# (cos, sin) pair, with the covariance of the participants' pairs over k, on
# F(2, k - 2) with the (k - 1)/(k - 2) factor of Bingham et al. (1982).
dance_ts_group_rhythm_ci <- function(bc, bs, period = 24, harmonic = 1, level = 0.95) {
  ok <- is.finite(bc) & is.finite(bs); bc <- bc[ok]; bs <- bs[ok]
  k <- length(bc)
  if (k < 3) return(NULL)
  S <- stats::cov(cbind(bc, bs))
  V <- S * (k - 1) / (k * (k - 2))
  dance_bingham_ci(mean(bc), mean(bs), V, n = k, n_params = 2, level = level,
                   period = period, harmonic = harmonic)
}

#   df    dance_ts_frame(mod, group_values)
dance_outcomes_ts <- function(df, mod, conf = 0.95) {
  if (is.null(df) || !nrow(df)) return(list(ok = FALSE, message = "No per-participant fits."))
  trend <- mod$trend_type %||% "none"
  K <- as.integer(mod$n_harmonics %||% 0L)
  period <- mod$period %||% 24
  u <- if (is.null(mod$dv_units) || !nzchar(mod$dv_units)) "response units" else mod$dv_units
  est <- list(); tst <- list()
  E <- function(...) est[[length(est) + 1L]] <<- dance_outcome_est_row(...)
  Tt <- function(...) tst[[length(tst) + 1L]] <<- dance_outcome_test_row(...)
  lab <- function(w) dance_ts_pair_label(w, period, trend, mod$dv_units)

  ls <- dance_ts_level_summary(df, mod, conf)
  if (isTRUE(ls$ok)) for (i in seq_len(nrow(ls$table))) {
    r <- ls$table[i, ]; q <- r$quantity
    dom <- if (q %in% c("fitted_mean", "mesor", "intercept")) "level" else "accumulation"
    prim <- (q == "fitted_mean") || (q == "mesor" && K == 0L) || (q == "change")
    E(dom, q, lab(q), r$cell, r$estimate, r$lo, r$hi, primary = prim, kind = "mean of participants")
  }
  tp <- switch(trend, linear = c(trend_linear = "trend_linear"), log = c(trend_log = "trend_log"),
               exp_sat = c(A_sat = "A_sat", if (!dance_ts_tau_held(mod)) c(tau = "tau", t50 = "t50")),
               character(0))
  if (length(tp)) {
    gm <- dance_ts_group_means(df, tp, conf)
    if (isTRUE(gm$ok)) for (i in seq_len(nrow(gm$table))) {
      r <- gm$table[i, ]
      E("accumulation_model", r$quantity, lab(r$quantity), r$cell, r$estimate, r$lo, r$hi,
        kind = "mean of participants")
    }
  }
  # tau held at one value for every participant: a constant, stated as such
  if (identical(trend, "exp_sat") && dance_ts_tau_held(mod)) {
    tau_h <- suppressWarnings(as.numeric(mod$tau_spec$value %||% NA))
    if (is.finite(tau_h)) for (cl in levels(droplevels(factor(df$group)))) {
      E("accumulation_model", "tau", "Trend time constant tau (hours; held for every participant)", cl,
        tau_h, kind = "fixed")
      E("accumulation_model", "t50", "Half-rise time tau ln 2 (hours)", cl, tau_h * log(2), kind = "fixed")
    }
  }
  for (h in seq_len(K)) {
    gr <- dance_ts_group_rhythm(df, mod, h); if (is.null(gr)) next
    for (i in seq_len(nrow(gr))) {
      sel <- as.character(df$group) == gr$cell[i]
      ci <- dance_ts_group_rhythm_ci(df[[paste0("beta_cos_", h)]][sel], df[[paste0("beta_sin_", h)]][sel],
                                     period, h, conf)
      E("rhythm", sprintf("amplitude%d", h), sprintf("Amplitude, H%d (%s h)", h, format(period / h)),
        gr$cell[i], gr$amplitude[i], ci$amplitude[1] %||% NA, ci$amplitude[2] %||% NA,
        kind = "population-mean cosinor", harmonic = h)
      E("rhythm", sprintf("phase%d", h), sprintf("Acrophase, H%d (model hours)", h), gr$cell[i],
        gr$acrophase_time[i],
        if (isTRUE(ci$identified)) ci$acrophase_hours[1] else NA,
        if (isTRUE(ci$identified)) ci$acrophase_hours[2] else NA,
        kind = "circular", harmonic = h,
        note = if (is.null(ci)) "fewer than 3 participants: no confidence region"
               else if (!isTRUE(ci$identified)) "not identified: the group's (cos, sin) region includes the origin"
               else NA_character_)
    }
  }
  # ---- tests ----
  if (nlevels(droplevels(factor(df$group))) >= 2) {
    mv <- dance_ts_manova(df, dance_ts_coef_cols(mod))
    if (isTRUE(mv$ok))
      Tt("omnibus", "full", "Complete fitted trajectory: every coefficient jointly", mv$F, mv$df1, mv$df2,
         mv$p, mv$method, primary = TRUE)
    for (r in dance_ts_components(df, mod)) {
      if (!isTRUE(r$ok) || is.na(r$domain %||% NA)) next
      Tt(r$domain, r$key %||% NA_character_, r$label, r$F, r$df1, r$df2, r$p, r$method,
         primary = isTRUE(r$primary), note = r$note %||% NA_character_)
    }
  }
  list(ok = TRUE, approach = "two-stage", estimates = dance_outcome_bind(est),
       tests = dance_outcome_bind(tst), structure = dance_model_structure(trend, K), trend = trend,
       conf = conf, units = mod$dv_units, period = period,
       clock_origin = dance_clock_origin(mod))
}

# ------------------------------------------------------------------------------
# DISPLAY: one value with its interval, in the outcome's own form -- an
# acrophase as a clock time on its harmonic's period, everything else as a
# number. Shared by the comparison tab and the APA report.
# ------------------------------------------------------------------------------
dance_outcome_value_text <- function(r, period = 24, clock_origin = 0, digits = 3) {
  if (!is.finite(r$estimate)) return("--")
  if (identical(r$kind, "circular")) {
    h <- if (is.finite(r$harmonic)) r$harmonic else 1L
    clk <- function(x) dance_clock_label(dance_acrophase_clock(hours = x, period = period, harmonic = h,
                                                               clock_origin = clock_origin)$hours,
                                         period, show_day = FALSE)
    out <- clk(r$estimate)
    if (is.finite(r$lo) && is.finite(r$hi)) out <- sprintf("%s [%s, %s]", out, clk(r$lo), clk(r$hi))
    return(out)
  }
  f <- function(x) formatC(x, format = "f", digits = digits)
  if (is.finite(r$lo) && is.finite(r$hi)) sprintf("%s [%s, %s]", f(r$estimate), f(r$lo), f(r$hi))
  else f(r$estimate)
}

# The test statistic in words: F(df1, df2) = x, or chi-square(df) = x.
dance_outcome_stat_text <- function(r) {
  if (!is.finite(r$statistic)) return("")
  if (is.finite(r$df2)) sprintf("F(%s, %s) = %.2f", format(r$df1), format(round(r$df2, 1)), r$statistic)
  else sprintf("\u03c7\u00b2(%s) = %.2f", format(r$df1), r$statistic)
}

# Both tables as one long CSV-ready frame, in domain order.
dance_outcomes_flat <- function(o) {
  if (!isTRUE(o$ok)) return(data.frame(message = o$message %||% "no outcomes"))
  ord <- names(DANCE_OUTCOME_DOMAINS)
  es <- o$estimates; ts <- o$tests
  es <- if (!is.null(es)) es[order(match(es$domain, ord)), , drop = FALSE]
  ts <- if (!is.null(ts)) ts[order(match(ts$domain, ord)), , drop = FALSE]
  a <- if (!is.null(es)) data.frame(table = "estimate", approach = o$approach,
                                    domain_label = unname(DANCE_OUTCOME_DOMAINS[es$domain]), es,
                                    statistic = NA, df1 = NA, df2 = NA, p = NA, method = NA,
                                    stringsAsFactors = FALSE)
  b <- if (!is.null(ts)) data.frame(table = "test", approach = o$approach,
                                    domain_label = unname(DANCE_OUTCOME_DOMAINS[ts$domain]),
                                    domain = ts$domain, outcome = ts$outcome, label = ts$label,
                                    cell = NA, estimate = NA, lo = NA, hi = NA, primary = ts$primary,
                                    kind = "test", harmonic = NA, note = ts$note,
                                    statistic = ts$statistic, df1 = ts$df1, df2 = ts$df2, p = ts$p,
                                    method = ts$method, stringsAsFactors = FALSE)
  out <- rbind(a, b); rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# IS TAU IDENTIFIED? One checklist per approach, and a verdict.
# ------------------------------------------------------------------------------
# Convergence is not identifiability: a fit can converge at a bound, or along a
# ridge where A and tau trade off and the data pin down only their combination
# (the early slope, the change over the window). The checks:
#   convergence    how many fits converged
#   boundary       how many tau estimates sit at an end of the search range
#   precision      the relative standard error of tau (its uncertainty)
#   dependence     the correlation of A and tau (|r| > 0.9: a ridge)
#   window         tau beyond the observation window: A is then extrapolated
#   held tau       the comparison with tau held at a value (Delta-AIC; LRT)
#   valid n        participants (or cells) contributing an interior, finite estimate
# The verdict is "identified", "weakly identified" or "not identified"; anything
# short of identified recommends a sensitivity analysis with tau held.
dance_tau_verdict <- function(rows, problems, weak) {
  status <- if (length(problems)) "not identified" else if (length(weak)) "weakly identified" else "identified"
  list(ok = TRUE, rows = rows, status = status, problems = problems, weak = weak,
       recommend_fixed = !identical(status, "identified"),
       message = switch(status,
         identified = "tau is identified by these data: converged, interior, precise, and not on an A-tau ridge.",
         paste0("tau is ", status, " by these data (", paste(c(problems, weak), collapse = "; "),
                "). Report the total change and the average rate, which the data determine, ",
                "before A or tau; and repeat the analysis with tau HELD at a plausible value ",
                "(a fixed-tau sensitivity analysis) to show the conclusions do not depend on it.")))
}

dance_tau_diagnostics_ts <- function(mod) {
  if (!identical(mod$trend_type, "exp_sat") || dance_ts_tau_held(mod))
    return(list(ok = FALSE, message = "tau is not estimated in this model."))
  p <- mod$individual_params; fits <- mod$individual_fits
  n_att <- mod$fit_audit$n_attempted %||% length(fits)
  if (is.null(p) || !nrow(p)) return(list(ok = FALSE, message = "No per-participant fits."))
  T <- diff(range(mod$time_vec, na.rm = TRUE))
  at_b <- vapply(p$subject, function(s) {
    f <- fits[[s]]; !is.null(f$bounds_hit) && any(grepl("^tau", f$bounds_hit)) }, logical(1))
  cr <- vapply(p$subject, function(s) {
    V <- fits[[s]]$vcov_full
    if (is.null(V) || !all(c("A_sat", "tau") %in% rownames(V))) return(NA_real_)
    d <- sqrt(diag(V)[c("A_sat", "tau")]); if (any(!is.finite(d) | d <= 0)) return(NA_real_)
    V["A_sat", "tau"] / prod(d) }, numeric(1))
  rse <- p$tau_se / p$tau
  conv <- if ("converged" %in% names(p)) p$converged else rep(TRUE, nrow(p))
  valid <- conv & !at_b & is.finite(p$tau) & is.finite(rse)
  cn <- mod$conditioning
  f3 <- function(x) format(signif(x, 3))
  rows <- data.frame(check = c("convergence", "boundary", "precision", "A-tau dependence", "window",
                               "tau held (Delta-AIC)", "valid n"),
    value = c(sprintf("%d of %d participants converged", sum(conv), n_att),
              sprintf("%d at an end of the tau range", sum(at_b)),
              sprintf("median relative SE of tau %s", f3(stats::median(rse[valid], na.rm = TRUE))),
              sprintf("median |corr(A, tau)| %s; %d above 0.9", f3(stats::median(abs(cr), na.rm = TRUE)),
                      sum(abs(cr) > 0.9, na.rm = TRUE)),
              sprintf("%d with tau beyond the %s h window", sum(p$tau > T, na.rm = TRUE), f3(T)),
              if (!is.null(cn$tau_fixed_delta_aic))
                sprintf("free minus held (tau = %s h): %s", f3(cn$tau_fixed_value), f3(cn$tau_fixed_delta_aic))
              else "not computed (set the held tau under Advanced)",
              sprintf("%d participants with an interior, finite tau", sum(valid))),
    stringsAsFactors = FALSE)
  problems <- c(if (sum(valid) < 3) "fewer than 3 participants have a usable tau",
                if (mean(at_b) > 0.25) "more than a quarter of the estimates sit on a bound")
  weak <- c(if (isTRUE(stats::median(abs(cr), na.rm = TRUE) > 0.9)) "A and tau lie on a ridge (median |r| > 0.9)",
            if (isTRUE(stats::median(rse[valid], na.rm = TRUE) > 0.5)) "tau is imprecise (median relative SE > 0.5)",
            if (mean(p$tau > T, na.rm = TRUE) > 0.5) "most tau estimates exceed the window, so A is extrapolated",
            if (isTRUE(cn$tau_fixed_delta_aic > 0)) "holding tau fits as well or better (Delta-AIC > 0)",
            if (any(at_b) && mean(at_b) <= 0.25) sprintf("%d estimate(s) on a bound", sum(at_b)))
  out <- dance_tau_verdict(rows, problems, weak)
  out$approach <- "two-stage"; out$corr <- cr; out$valid <- valid
  out
}

dance_tau_diagnostics_mixed <- function(ff) {
  if (!isTRUE(ff$ok) || !identical(ff$spec$trend, "exp_sat"))
    return(list(ok = FALSE, message = "No saturating trend in this model."))
  pr <- ff$tau_profile; sf <- ff$sat
  if (is.null(pr) && !isTRUE(sf$ok))
    return(list(ok = FALSE, message = "tau is held, not estimated, in this model."))
  T <- diff(range(ff$spec$data$t, na.rm = TRUE))
  f3 <- function(x) format(signif(x, 3))
  rows <- list(); problems <- character(0); weak <- character(0)
  add <- function(check, value) rows[[length(rows) + 1L]] <<- data.frame(check = check, value = value,
                                                                          stringsAsFactors = FALSE)
  if (!is.null(pr)) {
    add("profile (one tau)", if (isTRUE(pr$flat)) sprintf("flat: the log-likelihood varies by %s over the grid", f3(pr$logLik_range))
        else sprintf("tau = %s h, 95%% profile interval [%s, %s]%s", f3(pr$tau), f3(pr$ci[1]), f3(pr$ci[2]),
                     if (isTRUE(pr$ci_open)) " (open: reaches the grid's end)" else ""))
    add("boundary", if (isTRUE(pr$at_edge)) "the maximum is at an end of the tau grid" else "interior")
    add("inference", "the trajectory model's tests are conditional on the profiled tau (its uncertainty is not carried)")
    if (isTRUE(pr$flat)) problems <- c(problems, "the tau profile is flat")
    if (isTRUE(pr$at_edge)) problems <- c(problems, "the profile maximum is at the grid's edge")
    if (isTRUE(pr$ci_open)) weak <- c(weak, "the profile interval is open")
    if (is.finite(pr$tau) && pr$tau > T) weak <- c(weak, "tau exceeds the window, so A is extrapolated")
  }
  if (isTRUE(sf$ok)) {
    st <- sf$table
    add("per-cell model", sprintf("converged (%s); random structure: %s", sf$method %||% "ML", sf$re_label))
    if (isTRUE(sf$hierarchical_tau_requested))
      add("hierarchical tau", if (isTRUE(sf$hierarchical_tau))
            sprintf("tau varies between participants: SD %s h around the cell tau",
                    paste(sprintf("%s %s", st$cell, f3(st$tau_between_sd_h %||% NA)), collapse = ", "))
          else "requested, but that structure did not converge: participants share their cell's tau")
    if (isTRUE(sf$hierarchical_tau_requested) && !isTRUE(sf$hierarchical_tau))
      weak <- c(weak, "the hierarchical-tau structure did not converge")
    add("A-tau dependence", sprintf("|corr(A, log tau)| per cell: %s",
                                    paste(sprintf("%s %s", st$cell, f3(abs(st$corr_A_logtau))), collapse = ", ")))
    add("precision", sprintf("tau 95%% interval per cell: %s",
                             paste(sprintf("%s [%s, %s]", st$cell, f3(st$tau_lo), f3(st$tau_hi)), collapse = ", ")))
    if (isTRUE(sf$lrt$ok))
      add("one tau vs per cell", sprintf("LRT chi2(%d) = %s, p = %s", sf$lrt$df, f3(sf$lrt$statistic), f3(sf$lrt$p)))
    fl <- tryCatch(dance_traj_sat_flags(sf), error = function(e) character(0))
    if (any(abs(st$corr_A_logtau) > 0.9)) weak <- c(weak, "A and tau lie on a ridge in at least one cell")
    if (any(grepl("end of the search range", fl))) problems <- c(problems, "a cell's tau ran to the search range's end")
    if (any(st$tau > T)) weak <- c(weak, "a cell's tau exceeds the window, so its A is extrapolated")
  }
  out <- dance_tau_verdict(do.call(rbind, rows), unique(problems), unique(weak))
  out$approach <- "mixed"
  out
}
