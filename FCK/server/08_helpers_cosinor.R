# ==============================================================================
# server/08_helpers_cosinor.R — shared cosinor arithmetic and formatting
#
# CHANGELOG
# ---------
# 2026-09-03  Created for the harmonic-regression audit. Everything here was
#             previously either duplicated across the report block, computed
#             wrongly, or absent:
#
#   fck_format_equation()   ONE fitted-equation renderer. The pooled and group
#                           equations were built by two separate blocks in
#                           72_harmonic.R and the pooled one silently dropped
#                           the trend for EVERY trend type. (audit 1.1)
#   fck_resultants()        Returns the unweighted resultant (what the Rayleigh
#                           test is defined on) and the amplitude-weighted one
#                           (the vector-mean estimator) as separately named
#                           fields, so they cannot be swapped again. (audit 1.2)
#   fck_rayleigh()          Rayleigh Z and p from the UNWEIGHTED resultant,
#                           with the standard small-sample correction rather
#                           than the bare exp(-Z). (audit 1.2)
#   fck_commonality()       unique_S + unique_C + shared == total R2 exactly.
#                           Replaces two overlapping marginal R2s that summed
#                           to 124.7%. (audit 1.3)
#   fck_rhythm_adjusted_mean()
#                           The genuine rhythm-adjusted mean: the time-average
#                           of the non-oscillating part over the observation
#                           window, by integration. The quantity the report used
#                           to call the MESOR is not this. (audit 1.4)
#   fck_amp_se(), fck_acro_se()
#                           Delta-method SEs that use the full covariance block,
#                           not sqrt(se_c^2 + se_s^2)/sqrt(2). (audit, extra c)
#   fck_bingham_ci()        Elliptical joint amplitude-acrophase region.
#                           (audit 2.5)
#   fck_zero_amplitude_test()
#                           Full vs trend-only F test. The old test put the
#                           WHOLE model sum of squares over the harmonics' df
#                           alone, crediting the homeostatic trend to the
#                           rhythm. (audit, extra a)
#   fmt2()/fmtn()           round-half-away-from-zero then format, so 34.975
#                           prints as 34.98 rather than 34.97. (audit 1.7)
#
# Everything is base R. No package here is allowed to be a hard dependency:
# the app must still start when fda/minpack.lm/relaimpo are missing.
# ==============================================================================


# This module uses %||%, which server/00_state.R defines for the app. Defining
# it here too (only when absent) makes the file self-contained, so it can be
# sourced on its own by a test or a script without silently failing inside a
# helper that looks unrelated.
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a


# ------------------------------------------------------------------ formatting

# R's round() uses banker's rounding on the underlying binary double, so
# round(34.975, 2) is 34.97: 34.975 is not exactly representable and the stored
# value is a hair below the tie. sprintf("%.2f", .) inherits the same problem.
# Rounding half away from zero on the scaled value fixes both.
fck_round2 <- function(x, digits = 2) {
  z <- abs(x) * 10^digits
  # nudge by one ulp-ish so a value stored just under a tie still rounds up
  z <- floor(z + 0.5 + 1e-9)
  sign(x) * z / 10^digits
}

fmtn <- function(x, digits = 2, na = "NA") {
  out <- rep(na, length(x))
  ok <- is.finite(x)
  if (any(ok)) out[ok] <- formatC(fck_round2(x[ok], digits), format = "f", digits = digits)
  out
}

fmt1 <- function(x, na = "NA") fmtn(x, 1, na)
fmt2 <- function(x, na = "NA") fmtn(x, 2, na)
fmt3 <- function(x, na = "NA") fmtn(x, 3, na)
fmt4 <- function(x, na = "NA") fmtn(x, 4, na)

# A test statistic that can be 812 or 0.004 needs a format that stays readable
# at both ends; %g with 4 significant digits does, and never silently truncates
# a large value the way "%.2f" on a rounded double can.
fmt1e <- function(x, na = "NA") {
  out <- rep(na, length(x)); ok <- is.finite(x)
  if (any(ok)) out[ok] <- formatC(x[ok], format = "g", digits = 4)
  out
}


# ------------------------------------------------------------- phase conversion

# An acrophase in radians for harmonic h occupies the EFFECTIVE period T/h, so
# the conversion to hours divides by h. Getting this wrong is invisible for the
# fundamental and doubles the error for H2.
phi_to_hours <- function(phi, period = 24, harmonic = 1) {
  phi * period / (2 * pi) / harmonic
}

hours_to_phi <- function(hours, period = 24, harmonic = 1) {
  hours * 2 * pi * harmonic / period
}

phi_to_degrees <- function(phi) {
  phi * 180 / pi
}

# Harmonic h is identified only modulo T/h. Say so next to every value.
fck_phase_convention <- function(period = 24, harmonic = 1) {
  if (harmonic <= 1) {
    sprintf("clock hours, identified modulo %g h", period)
  } else {
    sprintf("clock hours, identified modulo %g h (H%d repeats %d times per %g h; %s h is equivalent to %s h, ... )",
            period / harmonic, harmonic, harmonic, period,
            "x", paste0("x+", fmtn(period / harmonic, 0)))
  }
}


# ------------------------------------------------------------ circular measures

# TWO different resultant lengths live on the same acrophases and they are not
# interchangeable:
#
#   r_unweighted = |mean(e^{i phi})|          -- unit vectors. This is the one
#                                                the Rayleigh test is defined on
#                                                (Mardia & Jupp 2000; Berens 2009).
#   r_weighted   = |sum(A e^{i phi})| / sum(A) -- amplitude-weighted. This is the
#                                                right accompaniment to the
#                                                amplitude-weighted vector mean,
#                                                and is NOT a Rayleigh input.
#
# The audit found Z computed from the weighted one (0.824 -> Z = 886) where the
# unweighted one (0.789 -> Z = 812) was required. Returning both, named, is the
# structural fix.
fck_resultants <- function(phi, amplitude = NULL) {
  ok <- is.finite(phi)
  if (!is.null(amplitude)) ok <- ok & is.finite(amplitude)
  phi <- phi[ok]
  n <- length(phi)
  if (n < 1) return(NULL)

  cx <- mean(cos(phi)); sy <- mean(sin(phi))
  r_unweighted <- sqrt(cx^2 + sy^2)
  mean_dir_unweighted <- atan2(sy, cx) %% (2 * pi)

  if (is.null(amplitude)) {
    amplitude <- rep(1, n)
    r_weighted <- r_unweighted
    mean_dir_weighted <- mean_dir_unweighted
    mean_amplitude <- 1
  } else {
    amplitude <- amplitude[ok]
    wx <- mean(amplitude * cos(phi)); wy <- mean(amplitude * sin(phi))
    mean_amplitude <- sqrt(wx^2 + wy^2)           # the vector-mean amplitude
    mean_dir_weighted <- atan2(wy, wx) %% (2 * pi)
    denom <- mean(amplitude)
    r_weighted <- if (is.finite(denom) && denom > 0) mean_amplitude / denom else NA_real_
  }

  list(n = n,
       r_unweighted = r_unweighted,
       r_weighted = r_weighted,
       mean_dir_unweighted = mean_dir_unweighted,
       mean_dir_weighted = mean_dir_weighted,
       mean_amplitude = mean_amplitude,
       # circular SD is defined on the UNWEIGHTED resultant
       circ_sd_rad = if (r_unweighted > 0 && r_unweighted < 1)
         sqrt(-2 * log(r_unweighted)) else NA_real_)
}

# Rayleigh test of uniformity. Z = n * r^2 on UNIT vectors. The p-value uses
# the standard second-order approximation (Zar, Biostatistical Analysis, eq.
# 27.4 / Berens 2009), not the bare exp(-Z) the code used before, which is only
# the leading term.
fck_rayleigh <- function(r_unweighted, n) {
  if (!is.finite(r_unweighted) || !is.finite(n) || n < 2)
    return(list(Z = NA_real_, p = NA_real_, n = n, r = r_unweighted))
  Z <- n * r_unweighted^2
  p <- exp(sqrt(1 + 4 * n + 4 * (n^2 - (n * r_unweighted)^2)) - (1 + 2 * n))
  list(Z = Z, p = min(1, max(0, p)), n = n, r = r_unweighted)
}


# ------------------------------------------------------- variance decomposition

# Commonality analysis (Chevan & Sutherland 1991; Ray-Mukherjee et al. 2014).
#
# Marginal R2s from two overlapping predictor blocks do NOT partition variance:
# R2_S + R2_C exceeds R2_total by however much S and C share, which for a
# saturating trend against a 24 h cosine over a 22 h window is a lot (0.229 in
# the reported output, invisible because it was never computed).
#
# The partition that DOES sum to the total is
#     unique_S = R2_full - R2_C_only
#     unique_C = R2_full - R2_S_only
#     shared   = R2_full - unique_S - unique_C
# and shared may legitimately be negative (suppression); that is information,
# not an error, so it is returned rather than clamped.
fck_commonality <- function(r2_full, r2_S_only, r2_C_only) {
  if (!all(is.finite(c(r2_full, r2_S_only, r2_C_only))))
    return(list(unique_S = NA_real_, unique_C = NA_real_, shared = NA_real_,
                total = r2_full, suppression = NA))
  unique_S <- r2_full - r2_C_only
  unique_C <- r2_full - r2_S_only
  shared   <- r2_full - unique_S - unique_C
  list(unique_S = unique_S,
       unique_C = unique_C,
       shared   = shared,
       total    = r2_full,
       suppression = isTRUE(shared < 0))
}

# The three parts as percentages of the total. Kept separate from the fractions
# so a caller cannot accidentally print percentages that do not sum to 100.
fck_commonality_pct <- function(cm) {
  if (!is.finite(cm$total) || cm$total <= 0)
    return(list(unique_S = NA_real_, unique_C = NA_real_, shared = NA_real_))
  list(unique_S = 100 * cm$unique_S / cm$total,
       unique_C = 100 * cm$unique_C / cm$total,
       shared   = 100 * cm$shared   / cm$total)
}


# --------------------------------------------------------- the non-oscillating
#                                                            part and its mean

# S(t) for each trend type, evaluated on the SAME anchoring the fitters use.
# fit_cosinor_nonlinear() builds the trend on (t - t_offset) while the harmonics
# run on raw t, so t_offset must be carried explicitly and never assumed to be 0.
fck_trend_value <- function(trend_type, coefs, t, t_offset = 0) {
  trend_type <- as.character(trend_type)
  switch(trend_type,
         "none"    = rep(0, length(t)),
         "linear"  = coefs[1] * t,
         "log"     = coefs[1] * log(t - t_offset + 1),
         "exp_sat" = coefs[1] * (1 - exp(-(t - t_offset) / coefs[2])),
         rep(0, length(t)))
}

# The rhythm-adjusted mean: the time-average of the non-oscillating part
# M + S(t) across the observation window.
#
# WHY THIS EXISTS. The quantity the report called "MESOR" is the fitted constant
# M. For a model with a trend that is not the rhythm-adjusted mean in Cornelissen
# (2014)'s sense -- it is the constant of a model whose trend is anchored at the
# first observation and whose harmonics are anchored at t = 0. The harmonics
# average to (almost) zero across a whole number of periods, so the
# rhythm-adjusted mean is mean(M + S(t)) over [t_min, t_max], which is what this
# returns. For an exp_sat fit with M = 27.70, A_sat = 32.30, tau = 13.92 over
# [8, 30] anchored at 8 it is 43.8 -- not 27.70, and not the 50 you get if you
# wrongly assume the trend is anchored at midnight.
fck_rhythm_adjusted_mean <- function(mesor, trend_type, trend_coefs,
                                     t_min, t_max, t_offset = 0) {
  trend_type <- as.character(trend_type)
  if (!is.finite(mesor) || !is.finite(t_min) || !is.finite(t_max) || t_max <= t_min)
    return(NA_real_)
  if (identical(trend_type, "none") || length(trend_coefs) == 0 ||
      !all(is.finite(trend_coefs)))
    return(mesor)

  f <- function(t) fck_trend_value(trend_type, trend_coefs, t, t_offset)
  area <- tryCatch(
    stats::integrate(f, t_min, t_max, subdivisions = 500L,
                     rel.tol = 1e-8)$value,
    error = function(e) {
      # closed-form / grid fallback so a failed quadrature never returns NA
      tg <- seq(t_min, t_max, length.out = 4001)
      mean(f(tg)) * (t_max - t_min)
    })
  mesor + area / (t_max - t_min)
}

# The mean of the harmonic part over the window, which is only exactly zero for
# a whole number of periods. Reported so the reader can see how much of the
# rhythm leaks into the window mean when the recording is not a whole cycle.
fck_harmonic_window_mean <- function(beta_cos, beta_sin, period, t_min, t_max) {
  n_h <- length(beta_cos)
  if (n_h < 1 || t_max <= t_min) return(NA_real_)
  total <- 0
  for (h in seq_len(n_h)) {
    w <- 2 * pi * h / period
    total <- total +
      beta_cos[h] * (sin(w * t_max) - sin(w * t_min)) / (w * (t_max - t_min)) -
      beta_sin[h] * (cos(w * t_max) - cos(w * t_min)) / (w * (t_max - t_min))
  }
  total
}


# ---------------------------------------------------------- delta-method errors

# SE of the amplitude A = sqrt(bc^2 + bs^2) from the covariance block of
# (bc, bs). The previous code used sqrt(se_c^2 + se_s^2)/sqrt(2), which is
# neither the delta-method result nor an approximation to it: it drops the
# covariance entirely and mis-weights the two terms.
fck_amp_se <- function(beta_cos, beta_sin, V) {
  A <- sqrt(beta_cos^2 + beta_sin^2)
  if (!is.finite(A) || A <= 0 || is.null(V) || any(!is.finite(V))) return(NA_real_)
  g <- c(beta_cos, beta_sin) / A
  v <- as.numeric(t(g) %*% V %*% g)
  if (!is.finite(v) || v < 0) NA_real_ else sqrt(v)
}

# SE of the acrophase phi = atan2(bs, bc), in radians.
fck_acro_se <- function(beta_cos, beta_sin, V) {
  A2 <- beta_cos^2 + beta_sin^2
  if (!is.finite(A2) || A2 <= 0 || is.null(V) || any(!is.finite(V))) return(NA_real_)
  g <- c(-beta_sin, beta_cos) / A2
  v <- as.numeric(t(g) %*% V %*% g)
  if (!is.finite(v) || v < 0) NA_real_ else sqrt(v)
}


# ------------------------------------------- Bingham elliptical joint region

# The joint (1-alpha) confidence region for the amplitude-acrophase pair, as an
# ellipse in (beta_cos, beta_sin) space, projected onto amplitude and acrophase
# limits. Bingham, Arbogast, Cornelissen Guillaume, Lee & Halberg (1982),
# Chronobiologia 9:397-439. This is the standard cosinor reporting requirement
# and the app had nothing of the kind.
#
# Returns NULL rather than a fabricated interval when the region contains the
# origin -- when it does, the acrophase is not identified at all and quoting
# limits for it would be worse than quoting none.
fck_bingham_ci <- function(beta_cos, beta_sin, V, n, n_params, level = 0.95,
                           period = 24, harmonic = 1) {
  if (is.null(V) || any(!is.finite(V)) || !is.finite(beta_cos) || !is.finite(beta_sin))
    return(NULL)
  df2 <- n - n_params
  if (!is.finite(df2) || df2 < 1) return(NULL)

  Fc <- stats::qf(level, 2, df2)
  cval <- 2 * Fc                       # radius^2 in the metric of V

  A <- sqrt(beta_cos^2 + beta_sin^2)
  if (!is.finite(A) || A <= 0) return(NULL)

  ev <- tryCatch(eigen(V, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev) || any(ev$values <= 0)) return(NULL)

  # Sample the ellipse boundary and read the amplitude / angle extremes off it.
  th <- seq(0, 2 * pi, length.out = 721)
  axes <- sqrt(cval * ev$values)
  pts <- ev$vectors %*% rbind(axes[1] * cos(th), axes[2] * sin(th))
  bx <- beta_cos + pts[1, ]
  by <- beta_sin + pts[2, ]

  radii <- sqrt(bx^2 + by^2)
  contains_origin <- min(radii) <= 0 ||
    (sum((c(0, 0) - c(beta_cos, beta_sin)) *
           solve(V, c(0, 0) - c(beta_cos, beta_sin)))) <= cval

  amp_lo <- min(radii); amp_hi <- max(radii)

  if (contains_origin) {
    return(list(level = level, amplitude = c(amp_lo, amp_hi),
                acrophase_rad = NULL, acrophase_hours = NULL,
                identified = FALSE,
                note = "the joint region covers the origin: the acrophase is not identified at this level"))
  }

  ang <- atan2(by, bx)
  centre <- atan2(beta_sin, beta_cos)
  # unwrap around the centre so min/max are the true angular extent
  d <- ((ang - centre + pi) %% (2 * pi)) - pi
  acro <- c(centre + min(d), centre + max(d))

  list(level = level,
       amplitude = c(amp_lo, amp_hi),
       acrophase_rad = acro %% (2 * pi),
       acrophase_hours = phi_to_hours(acro, period, harmonic) %% (period / harmonic),
       identified = TRUE,
       note = NULL)
}


# ------------------------------------------------- zero-amplitude F test, fixed

# The hypothesis "no rhythm" is: all harmonic coefficients are zero, GIVEN the
# trend. So the numerator sum of squares is the difference between the
# trend-only model and the full model -- not the whole model SS.
#
# The old code computed
#     F = ((ss_total - ss_resid_full) / (2*n_harmonics)) / (ss_resid_full / df2)
# which charges the numerator with every bit of variance the homeostatic trend
# explains while paying only the harmonics' degrees of freedom for it. With a
# saturating trend accounting for ~28% of variance on its own that is a large
# upward bias, and it -- not only the FDA smoothing -- is why 95.3% of subjects
# came out "significantly rhythmic".
fck_zero_amplitude_test <- function(ss_resid_full, ss_resid_trend_only,
                                    n, n_params_full, n_harmonics) {
  df1 <- 2 * n_harmonics
  df2 <- n - n_params_full
  if (!is.finite(ss_resid_full) || !is.finite(ss_resid_trend_only) ||
      df2 < 1 || ss_resid_full <= 0)
    return(list(F = NA_real_, df1 = df1, df2 = df2, p = NA_real_))
  num <- (ss_resid_trend_only - ss_resid_full) / df1
  den <- ss_resid_full / df2
  if (!is.finite(num) || !is.finite(den) || den <= 0 || num < 0)
    return(list(F = NA_real_, df1 = df1, df2 = df2, p = NA_real_))
  Fv <- num / den
  list(F = Fv, df1 = df1, df2 = df2,
       p = stats::pf(Fv, df1, df2, lower.tail = FALSE))
}


# ------------------------------------------------------------ the ONE formatter

# Renders a fitted equation from a parameter vector. The pooled and the
# per-group equations in 72_harmonic.R were built by two separate blocks; the
# pooled one read trend parameters out of a list that never contained any, so
# it silently printed the model WITHOUT its homeostatic term while the header,
# the symbolic equation and all four group equations included it. Over the
# observed window that under-predicted by roughly 20 units everywhere.
#
# There is now exactly one of these. Both call sites use it.
fck_format_equation <- function(mesor, trend_type = "none", trend_coefs = NULL,
                                amplitudes = NULL, acrophases_rad = NULL,
                                period = 24, t_offset = 0, digits = 2) {
  trend_type <- as.character(trend_type)
  eq <- paste0("Y(t) = ", fmtn(mesor, digits))

  term <- function(v, body, dg = digits) {
    if (!is.finite(v)) return(NULL)
    paste0(if (v >= 0) " + " else " - ", fmtn(abs(v), dg), body)
  }

  toff <- if (is.finite(t_offset) && abs(t_offset) > 1e-12)
    sprintf("(t - %s)", fmtn(t_offset, 0)) else "t"

  if (identical(trend_type, "linear") && length(trend_coefs) >= 1) {
    eq <- paste0(eq, term(trend_coefs[1], "·t", 3) %||% "")
  } else if (identical(trend_type, "log") && length(trend_coefs) >= 1) {
    eq <- paste0(eq, term(trend_coefs[1], sprintf("·log(%s + 1)", toff), 3) %||% "")
  } else if (identical(trend_type, "exp_sat") && length(trend_coefs) >= 2) {
    A <- trend_coefs[1]; tau <- trend_coefs[2]
    if (is.finite(A) && is.finite(tau) && tau > 0) {
      eq <- paste0(eq, if (A >= 0) " + " else " - ",
                   fmtn(abs(A), digits),
                   sprintf("·(1 - e^(-%s/%s))", toff, fmtn(tau, 1)))
    }
  }

  n_h <- length(amplitudes)
  if (n_h > 0) {
    for (h in seq_len(n_h)) {
      A <- amplitudes[h]; phi <- acrophases_rad[h]
      if (!is.finite(A) || !is.finite(phi)) next
      eq <- paste0(eq, if (A >= 0) " + " else " - ",
                   fmtn(abs(A), digits),
                   sprintf("·cos(2π·%d·t/%s - %s)",
                           h, fmtn(period, 0), fmtn(phi, digits)))
    }
  }
  eq
}


# ------------------------------------------------------------- group bookkeeping

# The reported group sizes were 654 + 410 + 181 + 59 = 1304 against 1305 fitted.
# unique() keeps NA as a level, then which(g == NA) is empty, so the subject with
# an unmatched label fell through the "n >= 3" guard with no message at all --
# it entered every pooled statistic and no group. Groups smaller than the guard
# vanished the same way.
#
# This returns the accounting explicitly so the report can print an UNASSIGNED
# row instead of losing rows in silence.
fck_group_audit <- function(group_labels, subject_ids, min_n = 3) {
  lab <- group_labels[subject_ids]
  assigned <- !is.na(lab) & nzchar(as.character(lab))
  tab <- table(as.character(lab[assigned]))
  too_small <- names(tab)[tab < min_n]
  list(
    n_total      = length(subject_ids),
    n_unassigned = sum(!assigned),
    unassigned_ids = subject_ids[!assigned],
    levels       = names(tab),
    counts       = as.integer(tab),
    dropped_small = too_small,
    n_dropped_small = sum(tab[too_small]),
    ok = sum(!assigned) == 0 && length(too_small) == 0
  )
}


# --------------------------------------------------------- model selection table

# AIC/AICc/BIC printed as means with SDs across subjects say nothing: with no
# competing model they are constant offsets of one another (which is why the
# three SDs printed identically at 16.21). What is interpretable is Delta-AICc
# against a nested set, with Akaike weights.
fck_akaike_table <- function(aicc_by_model) {
  ok <- vapply(aicc_by_model, function(v) is.finite(v), logical(1))
  aicc_by_model <- aicc_by_model[ok]
  if (!length(aicc_by_model)) return(NULL)
  v <- unlist(aicc_by_model)
  d <- v - min(v)
  w <- exp(-d / 2); w <- w / sum(w)
  data.frame(model = names(v), AICc = as.numeric(v),
             dAICc = as.numeric(d), weight = as.numeric(w),
             stringsAsFactors = FALSE)[order(d), ]
}


# ----------------------------------------------------------------- admissibility

# The DV is never named anywhere in the report and its bounds are never checked.
# With M + A1 + A2 the trough of the reported pooled fit sits near -3 before the
# trend is added, which is structurally impossible for a non-negative scale.
fck_check_bounds <- function(fitted_values, lower = NA, upper = NA) {
  out <- list(below = 0L, above = 0L, min = NA_real_, max = NA_real_, ok = TRUE)
  fv <- fitted_values[is.finite(fitted_values)]
  if (!length(fv)) return(out)
  out$min <- min(fv); out$max <- max(fv)
  if (is.finite(lower)) out$below <- sum(fv < lower)
  if (is.finite(upper)) out$above <- sum(fv > upper)
  out$ok <- out$below == 0L && out$above == 0L
  out
}


# ------------------------------------------------- group comparison machinery

# AUDIT 2.6: the report described group differences ("H1 amplitude declines
# 25.7 -> 22.6 -> 22.5 -> 21.1", "tau declines 15.4 -> 13.2 -> 11.2 -> 10.6")
# and tested none of them, with n = 59 in the smallest band. These give the
# tests, and give them WITH effect sizes and intervals rather than p-values
# alone.

# One-way comparison of a linear parameter, with eta^2 / Cohen's d and a CI on
# the largest contrast. Returns a printable list; never stops on degenerate
# input, because a group of 1 is a thing that happens.
fck_group_linear_test <- function(x, g, conf = 0.95) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  k <- nlevels(g)
  if (k < 2 || length(x) < k + 1) return(NULL)

  means <- tapply(x, g, mean, na.rm = TRUE)
  sds   <- tapply(x, g, stats::sd, na.rm = TRUE)
  ns    <- tapply(x, g, function(v) sum(is.finite(v)))

  grand <- mean(x)
  ss_between <- sum(ns * (means - grand)^2)
  ss_total   <- sum((x - grand)^2)
  ss_within  <- ss_total - ss_between
  df1 <- k - 1; df2 <- length(x) - k
  eta2 <- if (ss_total > 0) ss_between / ss_total else NA_real_
  # omega^2 is the less optimistic of the two and is what should be quoted
  ms_within <- ss_within / df2
  omega2 <- if (ss_total + ms_within > 0)
    (ss_between - df1 * ms_within) / (ss_total + ms_within) else NA_real_
  Fv <- if (df2 > 0 && ms_within > 0) (ss_between / df1) / ms_within else NA_real_
  p  <- if (is.finite(Fv)) stats::pf(Fv, df1, df2, lower.tail = FALSE) else NA_real_

  # largest pairwise contrast, with a Welch CI and Hedges' g
  lv <- levels(g); best <- NULL
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    d <- means[lv[i]] - means[lv[j]]
    if (is.null(best) || abs(d) > abs(best$diff)) {
      s1 <- sds[lv[i]]; s2 <- sds[lv[j]]; n1 <- ns[lv[i]]; n2 <- ns[lv[j]]
      if (!is.finite(s1) || !is.finite(s2) || n1 < 2 || n2 < 2) next
      se <- sqrt(s1^2 / n1 + s2^2 / n2)
      dfw <- se^4 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
      tcrit <- stats::qt(1 - (1 - conf) / 2, dfw)
      sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
      J <- 1 - 3 / (4 * (n1 + n2) - 9)          # Hedges' small-sample correction
      best <- list(a = lv[i], b = lv[j], diff = as.numeric(d),
                   ci = as.numeric(c(d - tcrit * se, d + tcrit * se)),
                   hedges_g = as.numeric(J * d / sp))
    }
  }

  list(k = k, n = length(x), levels = lv,
       means = means, sds = sds, ns = ns,
       F = Fv, df1 = df1, df2 = df2, p = p,
       eta2 = eta2, omega2 = omega2, largest = best)
}

# The trend across ORDERED groups, which is what "declines monotonically across
# age bands" actually claims. An omnibus ANOVA does not test monotonicity; a
# linear contrast on the ordered levels does, and it is far more powerful
# against exactly this alternative.
fck_group_trend_test <- function(x, g, conf = 0.95) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  k <- nlevels(g)
  if (k < 3) return(NULL)
  ns <- tapply(x, g, length); means <- tapply(x, g, mean)
  cw <- seq_len(k) - mean(seq_len(k))            # equally spaced linear contrast
  L <- sum(cw * means)
  ss_within <- sum(tapply(x, g, function(v) sum((v - mean(v))^2)))
  df2 <- length(x) - k
  ms_within <- ss_within / df2
  se <- sqrt(ms_within * sum(cw^2 / ns))
  tv <- L / se
  tcrit <- stats::qt(1 - (1 - conf) / 2, df2)
  list(L = L, se = se, t = tv, df = df2,
       p = 2 * stats::pt(-abs(tv), df2),
       ci = c(L - tcrit * se, L + tcrit * se),
       levels = levels(g), note = "equally-spaced linear contrast on the ordered levels")
}

# Watson-Williams requires reasonably concentrated data (conventionally r-bar
# >= 0.45, and kappa >= 2 for the F approximation to hold at all). The existing
# implementation gates on this; this reports the check so it appears in the
# output rather than only steering the code.
fck_ww_assumption <- function(r_bar, kappa = NA) {
  ok_r <- is.finite(r_bar) && r_bar >= 0.45
  ok_k <- !is.finite(kappa) || kappa >= 2
  list(ok = ok_r && ok_k, r_bar = r_bar, kappa = kappa,
       msg = if (!ok_r)
         sprintf("r-bar = %s is below 0.45: the von Mises F approximation does not hold and the Watson-Williams result should not be reported.",
                 fmt3(r_bar))
       else if (!ok_k)
         sprintf("kappa = %s is below 2: the F approximation is unreliable.", fmt2(kappa))
       else
         sprintf("r-bar = %s (>= 0.45) and kappa = %s (>= 2): the concentration assumption holds.",
                 fmt3(r_bar), fmt2(kappa)))
}


# --------------------------------------------------- clock labels for a linear axis

# WHY THE MODEL RUNS ON LINEAR TIME AND THE AXIS DOES NOT
#
# The harmonic part does not need it: cos(2*pi*h*t/T) is periodic, so t = 3 and
# t = 27 give identical values and a wrapped clock axis would serve.
#
# The TREND does need it. A linear, log or saturating S(t) is not periodic, so
# 08:00 on the first day and 08:00 on the second must be different values of t
# or the two days collapse onto one point and the homeostatic rise cannot be
# estimated at all. That is why the fitter unwraps 8..23,0..6 to 8..30, and why
# it must keep doing so.
#
# So linear time is a computational requirement and clock time is a DISPLAY
# concern. These convert one to the other for axes and hover text only; nothing
# here ever feeds a fit.
fck_clock_label <- function(t, period = 24, show_day = TRUE, with_minutes = TRUE) {
  out <- rep(NA_character_, length(t))
  ok <- is.finite(t)
  if (!any(ok)) return(out)
  tt <- t[ok]
  day <- floor(tt / period)
  hr <- tt %% period
  h <- floor(hr)
  m <- round((hr - h) * 60)
  # 59.7 minutes must roll the hour, not print as ":60"
  roll <- m >= 60
  m[roll] <- 0; h[roll] <- h[roll] + 1
  wrap <- h >= period
  h[wrap] <- h[wrap] - period; day[wrap] <- day[wrap] + 1

  lab <- if (with_minutes) sprintf("%02d:%02d", h, m) else sprintf("%02d:00", h)
  if (show_day && any(day != day[1])) {
    d0 <- min(day)
    suffix <- ifelse(day == d0, "", sprintf(" (+%dd)", day - d0))
    lab <- paste0(lab, suffix)
  }
  out[ok] <- lab
  out
}

# Tick positions on the LINEAR axis, labelled in clock time. Returns NULL when
# the range is degenerate, so a caller can fall back to plotly's own ticks.
fck_clock_ticks <- function(t_range, period = 24, target_n = 12) {
  if (!all(is.finite(t_range)) || diff(t_range) <= 0) return(NULL)
  span <- diff(t_range)
  # pick a step from the usual clock-friendly set
  cand <- c(0.25, 0.5, 1, 2, 3, 4, 6, 12, 24)
  step <- cand[which.min(abs(span / cand - target_n))]
  first <- ceiling(t_range[1] / step) * step
  vals <- seq(first, t_range[2], by = step)
  if (length(vals) < 2) return(NULL)
  list(vals = vals,
       text = fck_clock_label(vals, period, show_day = FALSE,
                              with_minutes = step < 1),
       step = step)
}


# ----------------------------------------------------------- parameter bounds

# WHICH bounds a fit ended up sitting on, not merely whether it sat on one.
#
# A parameter pinned to a constraint is the optimiser being stopped by the
# constraint rather than by the data: the reported value is the edge of the
# feasible region, and its standard error is meaningless there. That is worth
# knowing per parameter -- "tau ran to its ceiling" and "the amplitude hit its
# cap" are different problems with different fixes -- and worth counting across
# subjects, because a bound that catches 70% of the sample is a badly chosen
# bound, not 70% of subjects behaving oddly.
#
# Compared on a RELATIVE scale, so a tau ceiling of 110 and a tau floor of 0.5
# are judged with the same sensitivity.
fck_bounds_hit <- function(coefs, lower = NULL, upper = NULL, tol = 1e-6) {
  nm <- names(coefs)
  if (is.null(nm) || !length(coefs)) return(character(0))
  hit <- character(0)
  near <- function(v, b) {
    if (!is.finite(v) || !is.finite(b)) return(FALSE)
    scale <- max(abs(b), abs(v), 1)
    abs(v - b) <= tol * scale
  }
  for (p in nm) {
    v <- as.numeric(coefs[[p]])
    if (!is.null(lower) && p %in% names(lower) && near(v, lower[[p]]))
      hit <- c(hit, paste0(p, " (lower)"))
    if (!is.null(upper) && p %in% names(upper) && near(v, upper[[p]]))
      hit <- c(hit, paste0(p, " (upper)"))
  }
  hit
}

# Roll the per-subject bound lists into the two tables the report prints:
# how often each individual bound was hit, and how many subjects hit 1, 2, 3+.
fck_bounds_summary <- function(bounds_list, subject_ids = NULL) {
  n <- length(bounds_list)
  if (!n) return(NULL)
  counts <- vapply(bounds_list, length, integer(1))
  all_b <- unlist(bounds_list, use.names = FALSE)

  per_bound <- if (length(all_b)) {
    tb <- sort(table(all_b), decreasing = TRUE)
    data.frame(bound = names(tb), n = as.integer(tb),
               pct = 100 * as.integer(tb) / n, stringsAsFactors = FALSE)
  } else NULL

  k <- table(factor(counts, levels = 0:max(1L, max(counts))))
  per_count <- data.frame(n_bounds = as.integer(names(k)), n_subjects = as.integer(k),
                          pct = 100 * as.integer(k) / n, stringsAsFactors = FALSE)

  multi <- which(counts >= 2)
  multi_tbl <- if (length(multi)) {
    data.frame(
      subject = if (!is.null(subject_ids)) as.character(subject_ids)[multi] else as.character(multi),
      row = multi,
      n_bounds = counts[multi],
      bounds = vapply(bounds_list[multi], paste, character(1), collapse = ", "),
      stringsAsFactors = FALSE)
  } else NULL

  list(n = n, n_any = sum(counts > 0), n_multi = length(multi),
       per_bound = per_bound, per_count = per_count, multi = multi_tbl,
       counts = counts)
}


# ------------------------------------- model time -> clock time, for DISPLAY

# THE BUG THIS EXISTS TO PREVENT
#
# time_origin = "first_observation" re-anchors the model: the fit runs on
# t = 0, 1, ... 31 while the recording actually began at 08:00. mod$time_vec
# then holds MODEL time, and anything that labels it as a clock reads t = 0 as
# midnight -- so the fit plot started at 00:00 for a recording that started at
# 08:00, and the axis appeared to invert when the toggle was flipped.
#
# The shift back is one addition and it must happen at EVERY display boundary:
# axis ticks, hover text, and the angle a polar plot puts a point at. Never in
# a fit: the model keeps its own origin, which is the whole point of the option.
fck_model_to_clock <- function(t, mod) {
  t + (mod$origin_shift %||% 0)
}

# The clock range a model axis covers, for tick generation.
fck_model_clock_range <- function(mod) {
  tv <- mod$time_vec
  if (is.null(tv) || !any(is.finite(tv))) return(NULL)
  range(fck_model_to_clock(tv, mod), na.rm = TRUE)
}


# ============================================================================
# ACROPHASE: model coordinates -> clock time
#
# THE BUG THIS EXISTS TO PREVENT
# ------------------------------
# The harmonic coefficients are fitted against the model's time axis. With
# time_origin = "first_observation" that axis is ELAPSED hours from the first
# observation, so an acrophase of 19.18 h means "19.18 h after 08:00" -- which
# is 03:11 the next morning, not 19:11 this evening. The reporting layer
# printed the elapsed number as though it were a clock time, so every acrophase
# in the app was off by the origin shift while the fitted curves, which never
# left model coordinates, were right.
#
# The conversion is one line. What makes it worth a function is that there were
# a dozen places doing phi * period / (2*pi) by hand, some with the harmonic
# divisor and some without, and none of them adding the origin. One helper, and
# no site does its own arithmetic.
#
# HARMONIC h HAS h MAXIMA PER PERIOD
# ----------------------------------
# Harmonic h repeats every period/h, so within one 24 h clock it peaks h times,
# at clock_first + k * (period/h). Reporting only the first is not wrong but it
# is half the answer: for H2 with an elapsed acrophase of 7.81 h and an origin
# of 08:00 the maxima are 15:49 AND 03:49, and which one a reader cares about
# depends entirely on their hypothesis. All of them are returned.
#
# NOT THE SAME AS THE PEAK OF THE FITTED CURVE
# --------------------------------------------
# The H1 acrophase is where the FIRST HARMONIC peaks. The fitted curve also
# contains the homeostatic trend and the higher harmonics, so its maximum sits
# somewhere else -- often hours away when the trend dominates, as it does here.
# fck_curve_peak_clock() below reports that separately; conflating the two is
# the reason a reported acrophase can look wrong against its own plot even when
# both are right.
# ============================================================================

# The clock origin a model axis is measured from. Zero unless the fit was
# re-anchored at the first observation.
fck_clock_origin <- function(mod) {
  if (is.null(mod)) return(0)
  as.numeric(mod$origin_shift %||% 0)
}

# One acrophase, from either radians or model-elapsed hours, to clock time.
#
#   phi_rad  acrophase in radians, as atan2(beta_sin, beta_cos) returns it
#   hours    OR the same thing already in model-elapsed hours on the
#            harmonic's own effective period (period / harmonic)
#
# Returns the first maximum in [0, period), every maximum in [0, period), and
# the model-elapsed value it came from, so a caller can show its working.
fck_acrophase_clock <- function(phi_rad = NULL, hours = NULL, period = 24,
                                harmonic = 1, clock_origin = 0) {
  if (is.null(hours)) {
    if (is.null(phi_rad)) return(NULL)
    # phi occupies the effective period T/h, so the divisor is h
    hours <- phi_rad * period / (2 * pi) / harmonic
  }
  ok <- is.finite(hours)
  first <- rep(NA_real_, length(hours))
  first[ok] <- (hours[ok] + clock_origin) %% period

  eff <- period / harmonic
  all_h <- if (harmonic >= 1)
    lapply(first, function(f) if (!is.finite(f)) NA_real_
           else sort(((f + eff * seq(0, harmonic - 1)) %% period))) else NULL

  list(hours = first,          # first maximum on the clock
       all_hours = all_h,      # every maximum within one period
       elapsed = hours,        # what the model reported
       effective_period = eff,
       harmonic = harmonic,
       clock_origin = clock_origin)
}

# "03:11", or "15:49 / 03:49" for a harmonic with several maxima.
fck_acrophase_label <- function(phi_rad = NULL, hours = NULL, period = 24,
                                harmonic = 1, clock_origin = 0, all = TRUE) {
  a <- fck_acrophase_clock(phi_rad, hours, period, harmonic, clock_origin)
  if (is.null(a)) return(NA_character_)
  vapply(seq_along(a$hours), function(i) {
    if (!is.finite(a$hours[i])) return(NA_character_)
    if (!all || harmonic <= 1) return(fck_clock_label(a$hours[i], period, show_day = FALSE))
    paste(fck_clock_label(a$all_hours[[i]], period, show_day = FALSE), collapse = " / ")
  }, character(1))
}

# The maximum of the WHOLE fitted curve, over the observed window, in clock
# time. This is what a reader sees peaking on the fit plot, and it is not the
# H1 acrophase whenever a trend or a second harmonic is present.
# NOTE: this calls fck_rhythm_from_coefs() from server/07_helpers_circular.R.
# The app sources 07 before 08, so it is always in scope there; a script or test
# using this function must source 07 as well.
fck_curve_peak_clock <- function(coefs, mod, n_grid = 2001) {
  tv <- mod$time_vec
  if (is.null(coefs) || is.null(tv) || !any(is.finite(tv))) return(NULL)
  period <- mod$period %||% 24
  nh <- mod$n_harmonics %||% 1L
  trend <- mod$trend_type %||% "none"
  t_grid <- seq(min(tv, na.rm = TRUE), max(tv, na.rm = TRUE), length.out = n_grid)
  y <- fck_rhythm_from_coefs(coefs, t_grid, period, nh, trend,
                             include_trend = TRUE, t_offset = mod$t_offset %||% 0)
  if (!any(is.finite(y))) return(NULL)
  o <- fck_clock_origin(mod)
  i_max <- which.max(y); i_min <- which.min(y)
  list(peak_model = t_grid[i_max],
       peak_clock = (t_grid[i_max] + o) %% period,
       peak_value = y[i_max],
       trough_model = t_grid[i_min],
       trough_clock = (t_grid[i_min] + o) %% period,
       trough_value = y[i_min],
       # a peak at the very edge of the window is the trend still climbing, not
       # a maximum the data actually contain
       peak_at_edge = i_max <= 2 || i_max >= n_grid - 1,
       trough_at_edge = i_min <= 2 || i_min >= n_grid - 1)
}


# The fitted value at a given model time, in the DV's own units.
#
# "Constant term" answers "what is beta_0", which is a coefficient. Readers
# usually want "what does the model say the response was at the start", which is
# a different number whenever the harmonics are non-zero there -- and with a
# saturating trend anchored at the first observation they differ by exactly the
# harmonic sum, since S(t_min) = 0. Reporting both stops beta_0 being read as a
# starting level it never was.
fck_value_at <- function(coefs, mod, t_model) {
  if (is.null(coefs) || !is.finite(t_model)) return(NA_real_)
  v <- fck_rhythm_from_coefs(coefs, t_model, mod$period %||% 24,
                             mod$n_harmonics %||% 1L, mod$trend_type %||% "none",
                             include_trend = TRUE, t_offset = mod$t_offset %||% 0)
  as.numeric(v[1])
}

# ==============================================================================
# THE REDUCED MODEL FOR THE ZERO-AMPLITUDE TEST (P0.3)
# ==============================================================================
# The zero-amplitude F test compares the full model against the same model with
# every harmonic coefficient set to zero. For a `linear` or `log` trend the
# reduced model is linear in its parameters, so lm() already returns the exact
# least-squares fit and nothing was wrong there.
#
# For `exp_sat` it is NOT. The reduced model is
#
#     y = mesor + A_sat * (1 - exp(-(t - t0)/tau)) + e
#
# which is nonlinear in tau. The previous code held tau at the value the FULL
# model estimated and refitted only mesor and A_sat by OLS:
#
#     b <- 1 - exp(-(time - t_offset) / trend_params$tau$coef)
#     sum(residuals(lm(y ~ b))^2)
#
# I wrote that, for numerical robustness, and documented the shortcut without
# working out what it costs. Under the null tau is a free parameter; freezing it
# leaves SSE0 un-minimised, so the numerator (SSE0 - SSE1) is inflated and F is
# biased UPWARDS. Measured on 3,000 simulated nulls on the Circaflex time grid,
# the test rejected at 14.6% against a nominal 5%, and 4.5% against a nominal
# 1%. It manufactured significance.
#
# This refits tau. If the refit fails to converge the function returns NA rather
# than silently falling back to the frozen-tau value, because a wrong F is worse
# than an absent one.
fck_reduced_exp_sat_sse <- function(y, time, t_offset, start, lower = NULL, upper = NULL) {
  ok <- is.finite(y) & is.finite(time)
  y <- y[ok]; time <- time[ok]
  if (length(y) < 4) return(NA_real_)

  st <- list(mesor = as.numeric(start[["mesor"]]),
             A_sat = as.numeric(start[["A_sat"]]),
             tau   = as.numeric(start[["tau"]]))
  if (any(!vapply(st, function(v) length(v) == 1 && is.finite(v), logical(1))))
    return(NA_real_)
  if (st$tau <= 0) st$tau <- max(1, diff(range(time)) / 4)

  pick <- function(b, nm, dflt) {
    if (is.null(b) || !nm %in% names(b)) return(dflt)
    v <- as.numeric(b[[nm]]); if (length(v) != 1 || is.na(v)) dflt else v
  }
  lo <- c(mesor = pick(lower, "mesor", -Inf),
          A_sat = pick(lower, "A_sat", -Inf),
          tau   = max(pick(lower, "tau", 1e-6), 1e-6))
  up <- c(mesor = pick(upper, "mesor", Inf),
          A_sat = pick(upper, "A_sat", Inf),
          tau   = pick(upper, "tau", Inf))
  # A start sitting exactly on a bound stalls the optimiser.
  for (nm in names(st)) {
    if (is.finite(lo[[nm]]) && st[[nm]] <= lo[[nm]])
      st[[nm]] <- lo[[nm]] + max(abs(lo[[nm]]) * 1e-3, 1e-3)
    if (is.finite(up[[nm]]) && st[[nm]] >= up[[nm]])
      st[[nm]] <- up[[nm]] - max(abs(up[[nm]]) * 1e-3, 1e-3)
  }

  resid_fn <- function(p) y - (p[1] + p[2] * (1 - exp(-(time - t_offset) / p[3])))
  fit <- tryCatch(
    minpack.lm::nls.lm(par = unlist(st), fn = resid_fn,
                       lower = lo[names(st)], upper = up[names(st)],
                       control = minpack.lm::nls.lm.control(maxiter = 200,
                                                            maxfev = 2000)),
    error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  sse <- sum(resid_fn(fit$par)^2)
  if (!is.finite(sse)) return(NA_real_) else sse
}

# ---- degenerate-fit guards (P0.6) -------------------------------------------
# A constant trajectory has SST = 0, so 1 - SSE/SST is 0/0 = NaN. A perfect fit
# has SSE = 0, so log(sigma^2) = -Inf, the Gaussian log-likelihood is +Inf and
# AIC is -Inf -- a flat subject then WINS every model comparison it enters.
# Neither case is an error in the data; both need a defined answer.
fck_r_squared <- function(ss_resid, ss_total) {
  if (!is.finite(ss_resid) || !is.finite(ss_total)) return(NA_real_)
  # No variance to explain: R^2 is undefined, not 1 and not 0.
  if (ss_total <= .Machine$double.eps) return(NA_real_)
  1 - ss_resid / ss_total
}

# Gaussian log-likelihood with a floor on sigma^2. The floor is relative to the
# data scale so it does not depend on the units, and it is only ever reached by
# a fit that reproduces the observations exactly.
fck_gaussian_loglik <- function(ss_resid, n, y_scale = NULL) {
  if (!is.finite(ss_resid) || !is.finite(n) || n < 1) return(NA_real_)
  sigma_sq <- ss_resid / n
  floor_sq <- if (!is.null(y_scale) && is.finite(y_scale) && y_scale > 0)
    (1e-8 * y_scale)^2 else 1e-16
  if (sigma_sq < floor_sq) sigma_sq <- floor_sq
  -n / 2 * (log(2 * pi) + log(sigma_sq) + 1)
}
