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
  p
}

# The coefficient columns a MANOVA compares: constant, trend, and every
# (cos, sin) pair -- the complete per-participant fit.
dance_ts_coef_cols <- function(mod) {
  tr <- switch(mod$trend_type %||% "none",
               linear = "trend_linear", log = "trend_log",
               exp_sat = c("A_sat", "tau"), character(0))
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
  lt_row("Constant term (beta_0)", df$mesor,
         "the fitted constant; NOT a MESOR when a trend is present")
  if ("mesor_adj" %in% names(df) && !identical(mod$trend_type %||% "none", "none"))
    lt_row("MESOR (rhythm-adjusted mean over the window)", df$mesor_adj)
  for (tc in setdiff(dance_ts_coef_cols(mod), c("mesor", grep("^beta_", dance_ts_coef_cols(mod), value = TRUE))))
    if (tc %in% names(df)) lt_row(sprintf("Trend: %s", tc), df[[tc]])
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
dance_ts_pair_whats <- function(mod) {
  K <- mod$n_harmonics %||% 1L
  tr <- switch(mod$trend_type %||% "none", linear = "trend_linear", log = "trend_log",
               exp_sat = c("A_sat", "tau"), character(0))
  out <- c("level", if (length(tr)) c("mesor_adj", "value_at_start"), tr,
           sprintf("amplitude%d", seq_len(K)), sprintf("phase%d", seq_len(K)), "r_squared")
  out
}
dance_ts_pair_label <- function(what, period = 24) {
  pp <- dance_traj_pair_parse(what)
  if (!is.na(pp$harmonic))
    return(sprintf("%s — H%d (%s h)", if (identical(pp$base, "amplitude")) "Amplitude" else "Acrophase",
                   pp$harmonic, format(round(period / pp$harmonic, 2), trim = TRUE)))
  switch(what,
    level = "Constant term (beta_0)", mesor_adj = "MESOR (rhythm-adjusted mean)",
    value_at_start = "Predicted value at the first observation",
    trend_linear = "Linear trend (beta)", trend_log = "Log trend (beta)",
    A_sat = "A_sat (asymptote)", tau = "tau (time constant)",
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
  } else if (identical(what, "level")) "mesor" else what
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
       joint = FALSE, label = dance_ts_pair_label(what, period),
       method = if (circular) "Watson-Williams F on the per-participant phases"
                else "Welch's t on the per-participant estimates",
       unit = if (circular)
         sprintf("acrophase difference in hours on an effective period of %.4g h (shortest arc); no interval -- the test is the statement", eff)
       else "difference of group means of the per-participant estimates, with Welch's interval and Cohen's d",
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
dance_ts_group_curves <- function(mod, group_values, times, include_trend = TRUE,
                                  conf = 0.95) {
  df <- dance_ts_frame(mod, group_values)
  if (is.null(df)) return(NULL)
  period <- mod$period %||% 24; K <- mod$n_harmonics %||% 1L
  trend <- mod$trend_type %||% "none"
  fits <- mod$individual_fits
  g <- droplevels(factor(df$group))
  q <- stats::qnorm(1 - (1 - conf) / 2)
  out <- list()
  for (lv in levels(g)) {
    rows <- which(g == lv)
    M <- do.call(rbind, lapply(rows, function(i) {
      f <- fits[[df$subject[i]]]
      if (is.null(f) || !isTRUE(f$success) || is.null(f$coefs)) return(NULL)
      pr <- tryCatch(dance_rhythm_from_coefs(f$coefs, times, period, K, trend,
                                             include_trend = include_trend,
                                             t_offset = f$t_offset %||% mod$t_offset %||% 0),
                     error = function(e) NULL)
      if (is.null(pr) || !all(is.finite(pr))) NULL else pr
    }))
    if (is.null(M) || nrow(M) < 2) next
    # the group's vector-mean coefficients; with no grouping, the pooled
    # population fit the summary prints -- the same coefficients in both cases
    gf <- mod$group_fits[[lv]]
    mc <- if (!is.null(gf)) gf$mean_coefs
          else if (identical(lv, "(all)")) mod$pop_mean_fit$mean_coefs else NULL
    line <- if (!is.null(mc))
      dance_rhythm_from_coefs(mc, times, period, K, trend,
                              include_trend = include_trend, t_offset = mod$t_offset %||% 0)
      else colMeans(M)
    se <- apply(M, 2, stats::sd) / sqrt(nrow(M))
    out[[lv]] <- data.frame(cell = lv, t = times, fit = line, se = se,
                            lo = line - q * se, hi = line + q * se,
                            mean_of_curves = colMeans(M), n = nrow(M),
                            stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(ok = TRUE, table = tab, cells = names(out), times = times, conf = conf,
       band = "pointwise",
       note = sprintf(paste("Pointwise %.0f%% band: the group's mean-coefficient curve +/- z * SE(t),",
                            "with SE(t) the standard error across participants of their own",
                            "fitted curves. Not simultaneous."), 100 * conf))
}
