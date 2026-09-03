# ==============================================================================
# tests/report_harness.R — old vs new report, on data with KNOWN truth
#
# WHY A SYNTHETIC HARNESS
# -----------------------
# Deliverable 4 asked for the example report regenerated from the same input
# data so the two could be diffed line by line. That input was not available:
# the 1305-subject file is not in this repository, and fda / minpack.lm could
# not be installed in the environment this work was done in (CRAN blocked by
# the egress policy). Producing a "regenerated report" without the data would
# have meant inventing numbers and presenting them as a re-run, which is worse
# than not having one.
#
# What this does instead is stronger in one respect and weaker in another.
#   Stronger: the parameters are SET by this script, so every quantity has a
#             known correct answer and the old and new code can be scored
#             against truth, not merely against each other.
#   Weaker:   it is not your data, so it cannot reproduce your exact figures.
#
# The protocol is yours: 16 unequally spaced clock times 8..6, linearised to
# 8..30, four age bands in the reported proportions, a saturating homeostatic
# rise, two harmonics, and one subject with a missing group label.
#
# Run:  Rscript tests/report_harness.R          (from the FCK directory)
# Needs no packages: tau is profiled on a grid and the rest solved linearly,
# which is exact for this model up to the grid resolution.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

set.seed(2026)

# ---- the protocol, exactly as specified in the audit brief -------------------
clock <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 0, 2, 4, 6)
tlin  <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30)
period <- 24; n_h <- 2; t_off <- min(tlin)

n_by_group <- c(YOUTH = 654, ADULT = 410, MIDDLE_AGE = 181, ELDERLY = 59)

# ---- GROUND TRUTH -----------------------------------------------------------
# Chosen near the reported values so the harness exercises the same regime --
# a trend-dominated fit over a 22 h window -- and not a comfortable one.
TRUTH <- list(
  intercept = c(YOUTH = 27.0, ADULT = 27.6, MIDDLE_AGE = 28.1, ELDERLY = 28.6),
  A_sat     = c(YOUTH = 33.0, ADULT = 32.0, MIDDLE_AGE = 31.5, ELDERLY = 31.0),
  tau       = c(YOUTH = 15.4, ADULT = 13.2, MIDDLE_AGE = 11.2, ELDERLY = 10.6),
  amp1      = c(YOUTH = 25.7, ADULT = 22.6, MIDDLE_AGE = 22.5, ELDERLY = 21.1),
  phi1      = c(YOUTH = 0.687, ADULT = 0.710, MIDDLE_AGE = 0.720, ELDERLY = 0.715),
  amp2      = c(YOUTH = 7.2,  ADULT = 6.8,  MIDDLE_AGE = 6.5,  ELDERLY = 6.2),
  phi2      = c(YOUTH = 2.04, ADULT = 2.05, MIDDLE_AGE = 2.06, ELDERLY = 2.05),
  sigma     = 4.0
)

groups <- rep(names(n_by_group), n_by_group)
N <- length(groups)

curve_of <- function(g, t) {
  TRUTH$intercept[g] +
    TRUTH$A_sat[g] * (1 - exp(-(t - t_off) / TRUTH$tau[g])) +
    TRUTH$amp1[g] * cos(2 * pi * 1 * t / period - TRUTH$phi1[g]) +
    TRUTH$amp2[g] * cos(2 * pi * 2 * t / period - TRUTH$phi2[g])
}

Y <- t(vapply(groups, function(g)
  as.numeric(curve_of(g, tlin) + rnorm(length(tlin), 0, TRUTH$sigma)),
  numeric(length(tlin))))

# AUDIT 1.5: one subject whose group label is missing. In the shipped code that
# subject entered every pooled statistic and no group, and the group sizes then
# summed to N-1 with no message anywhere. Here it is deliberate, so the fix is
# visible in the diff.
group_labels <- groups
group_labels[1] <- NA_character_

# ---- an exact-enough exp_sat fitter, without minpack.lm ----------------------
fit_one <- function(y, t) {
  taus <- exp(seq(log(1), log(80), length.out = 300))
  best <- NULL
  for (tv in taus) {
    X <- cbind(1, 1 - exp(-(t - t_off) / tv))
    for (h in seq_len(n_h)) {
      w <- 2 * pi * h / period
      X <- cbind(X, cos(w * t), sin(w * t))
    }
    f <- tryCatch(lm.fit(X, y), error = function(e) NULL)
    if (is.null(f)) next
    ssr <- sum(f$residuals^2)
    if (is.null(best) || ssr < best$ssr)
      best <- list(ssr = ssr, tau = tv, b = f$coefficients, p = ncol(X))
  }
  if (is.null(best)) return(NULL)
  b <- best$b
  amps <- acros <- numeric(n_h)
  for (h in seq_len(n_h)) {
    bc <- b[1 + 2 * h]; bs <- b[2 + 2 * h]
    amps[h] <- sqrt(bc^2 + bs^2)
    acros[h] <- atan2(bs, bc) %% (2 * pi)
  }
  n <- length(y)
  ss_total <- sum((y - mean(y))^2)
  Xt <- cbind(1, 1 - exp(-(t - t_off) / best$tau))
  ss_trend_only <- sum(lm.fit(Xt, y)$residuals^2)
  Xc <- matrix(1, n, 1)
  for (h in seq_len(n_h)) {
    w <- 2 * pi * h / period
    Xc <- cbind(Xc, cos(w * t), sin(w * t))
  }
  ss_circ_only <- sum(lm.fit(Xc, y)$residuals^2)

  list(intercept = b[1], A_sat = b[2], tau = best$tau,
       amps = amps, acros = acros,
       r2 = 1 - best$ssr / ss_total,
       r2_S = 1 - ss_trend_only / ss_total,
       r2_C = 1 - ss_circ_only / ss_total,
       ss_total = ss_total, ss_resid = best$ssr, ss_trend_only = ss_trend_only,
       n = n, n_params = best$p + 1)
}

cat("Fitting", N, "subjects (grid-profiled tau)...\n")
fits <- vector("list", N)
for (i in seq_len(N)) {
  fits[[i]] <- fit_one(as.numeric(Y[i, ]), tlin)
  if (i %% 300 == 0) cat("  ", i, "\n")
}
ok <- !vapply(fits, is.null, logical(1))
fits <- fits[ok]; lab <- group_labels[ok]
n_fit <- length(fits)

P <- data.frame(
  intercept = vapply(fits, function(f) f$intercept, 1),
  A_sat = vapply(fits, function(f) f$A_sat, 1),
  tau   = vapply(fits, function(f) f$tau, 1),
  amp1  = vapply(fits, function(f) f$amps[1], 1),
  phi1  = vapply(fits, function(f) f$acros[1], 1),
  amp2  = vapply(fits, function(f) f$amps[2], 1),
  phi2  = vapply(fits, function(f) f$acros[2], 1),
  r2    = vapply(fits, function(f) f$r2, 1),
  r2_S  = vapply(fits, function(f) f$r2_S, 1),
  r2_C  = vapply(fits, function(f) f$r2_C, 1)
)

res1 <- fck_resultants(P$phi1, P$amp1)
res2 <- fck_resultants(P$phi2, P$amp2)
pool_int  <- mean(P$intercept); pool_A <- mean(P$A_sat); pool_tau <- mean(P$tau)
pool_amp1 <- res1$mean_amplitude; pool_phi1 <- res1$mean_dir_weighted
pool_amp2 <- res2$mean_amplitude; pool_phi2 <- res2$mean_dir_weighted

cm_each <- Map(fck_commonality, P$r2, P$r2_S, P$r2_C)
uS <- mean(vapply(cm_each, function(z) z$unique_S, 1))
uC <- mean(vapply(cm_each, function(z) z$unique_C, 1))
sh <- mean(vapply(cm_each, function(z) z$shared, 1))

ram <- fck_rhythm_adjusted_mean(pool_int, "exp_sat", c(pool_A, pool_tau),
                                min(tlin), max(tlin), t_off)

line <- function() cat(strrep("=", 78), "\n")

# =============================================================================
sink(file.path(app_dir, "tests/report_OLD.txt"))
line(); cat("REPORT AS THE CODE PRODUCED IT BEFORE THE AUDIT\n"); line()
cat("\n--- Population Mean Parameters (Vector-Averaged) ---\n")
cat(sprintf("MESOR:     %.3f\n", pool_int))
cat(sprintf("Amplitude (H1): %.3f\n", pool_amp1))
cat(sprintf("Acrophase (H1): %.2f (%.2f hours)\n",
            pool_phi1 * 180 / pi, pool_phi1 * period / (2 * pi)))
cat("\n  All Harmonics (vector-averaged):\n")
cat(sprintf("    H1: Amplitude=%.3f, Acrophase=%.2f hours\n",
            pool_amp1, pool_phi1 * period / (2 * pi)))
cat(sprintf("    H2: Amplitude=%.3f, Acrophase=%.2f hours\n",
            pool_amp2, pool_phi2 * period / (2 * pi) / 2))
z_old <- nrow(P) * res1$r_weighted^2          # the shipped Rayleigh
cat("\nRayleigh test for uniformity (H1):\n")
cat(sprintf("  Z = %.3f, p = %.4f\n", z_old, exp(-z_old)))
cat(sprintf("  Mean resultant length (r-bar) = %.3f\n", res1$r_weighted))
cat("\n--- Fitted Model Equation ---\n")
# the bug: indiv_means never carried A_sat/tau, so the trend branch was dead
cat(sprintf("Y(t) = %.2f + %.2f*cos(2pi*1*t/24 - %.2f) + %.2f*cos(2pi*2*t/24 - %.2f)\n",
            pool_int, pool_amp1, pool_phi1, pool_amp2, pool_phi2))
cat("\n--- Individual Parameter Summary ---\n")
cat(sprintf("MESOR:     Mean=%.3f, SD=%.3f\n", mean(P$intercept), sd(P$intercept)))
cat(sprintf("A_sat (asymptote): Mean=%.3f, SD=%.3f (units)\n", mean(P$A_sat), sd(P$A_sat)))
cat(sprintf("tau (time constant): Mean=%.2f, SD=%.2f (hours)\n", mean(P$tau), sd(P$tau)))
cat(sprintf("H1 Amplitude: Mean=%.3f, SD=%.3f\n", mean(P$amp1), sd(P$amp1)))
cat(sprintf("H1 Acrophase: Circular mean=%.2f h, r-bar=%.3f\n",
            res1$mean_dir_unweighted * period / (2 * pi), res1$r_unweighted))
cat(sprintf("H2 Acrophase: Circular mean=%.2f h\n",
            res2$mean_dir_unweighted * period / (2 * pi)))   # no /2: the scale bug
cat("\n--- Variance Decomposition: Relative Importance ---\n")
cat(sprintf("R2 from Process S (homeostatic): Mean=%.3f\n", mean(P$r2_S)))
cat(sprintf("R2 from Process C (circadian):   Mean=%.3f\n", mean(P$r2_C)))
cat("\nProportion of total R2 explained by:\n")
ps <- mean(100 * P$r2_S / P$r2); pc <- mean(100 * P$r2_C / P$r2)
cat(sprintf("  Process S: %.1f%%\n", ps))
cat(sprintf("  Process C: %.1f%%\n", pc))
cat(sprintf("  [these sum to %.1f%%]\n", ps + pc))
if (pc > ps) cat("\nInterpretation: Circadian rhythm (C) is the dominant component\n")
cat("\n--- Group-Specific Parameters ---\n")
gs <- 0L
for (g in names(n_by_group)) {
  idx <- which(!is.na(lab) & lab == g)
  if (length(idx) < 3) next
  gs <- gs + length(idx)
  rr <- fck_resultants(P$phi1[idx], P$amp1[idx])
  cat(sprintf("\nGroup '%s' (n=%d):\n", g, length(idx)))
  cat(sprintf("  MESOR:     %.3f (SD=%.3f)\n", mean(P$intercept[idx]), sd(P$intercept[idx])))
  cat(sprintf("  H1 Amplitude: %.3f (SD=%.3f)\n", rr$mean_amplitude, sd(P$amp1[idx])))
  cat(sprintf("  H1 Acrophase: %.2f hours\n", rr$mean_dir_weighted * period / (2 * pi)))
}
cat(sprintf("\n[group n's sum to %d; %d subjects were fitted]\n", gs, n_fit))
sink()

# =============================================================================
sink(file.path(app_dir, "tests/report_NEW.txt"))
line(); cat("REPORT AFTER THE AUDIT\n"); line()
cat("\nDependent variable: synthetic sleepiness (KSS-like points)\n")
cat("  Admissible range: [1.00, 100.00]\n")
cat(sprintf("Period: %s h    Harmonics: %d    Trend: SATURATING EXPONENTIAL\n",
            fmtn(period, 0), n_h))
cat("Data: RAW (no smoothing applied)\n")
cat("Time origin: MIDNIGHT (t = 0 at clock 00:00).\n")
cat(sprintf("  ! The trend is anchored at the first observation (t - %s) while the\n", fmt2(t_off)))
cat("    harmonics are anchored at midnight. The intercept is the constant of a\n")
cat("    model with two origins and is not the value at either.\n")

cat("\n--- Fit outcomes ---\n")
cat(sprintf("Subjects attempted:         %d\n", N))
cat(sprintf("  Converged:                %d\n", n_fit))
cat(sprintf("  Did not converge:         %d  (excluded)\n", N - n_fit))

cat("\n--- Central value ---\n")
cat(sprintf("Intercept (beta_0, at t = 0):                           %s\n", fmt3(pool_int)))
cat(sprintf("MESOR (rhythm-adjusted mean over the observed window):  %s\n", fmt3(ram)))
cat(sprintf("  = time-average of beta_0 + S(t) over t in [%s, %s] h, by integration.\n",
            fmt2(min(tlin)), fmt2(max(tlin))))
cat(sprintf("  Difference: %s units. The intercept is NOT the MESOR here.\n", fmt2(ram - pool_int)))

cat("\n--- Population rhythm parameters (VECTOR-averaged) ---\n")
cat(sprintf("  H1: amplitude = %s, acrophase = %s h\n",
            fmt3(pool_amp1), fmt2(phi_to_hours(pool_phi1, period, 1))))
cat(sprintf("  H2: amplitude = %s, acrophase = %s h\n",
            fmt3(pool_amp2), fmt2(phi_to_hours(pool_phi2, period, 2))))
cat(sprintf("      convention: H2 is identified modulo %s h, so %s h == %s h\n",
            fmtn(period / 2, 0), fmt2(phi_to_hours(pool_phi2, period, 2)),
            fmt2(phi_to_hours(pool_phi2, period, 2) + 12)))

cat("\n--- Circular concentration and the Rayleigh test ---\n")
for (h in 1:2) {
  rr <- if (h == 1) res1 else res2
  ry <- fck_rayleigh(rr$r_unweighted, rr$n)
  cat(sprintf("H%d  n = %d\n", h, rr$n))
  cat(sprintf("  r-bar (UNWEIGHTED, for Rayleigh):                %s\n", fmt3(rr$r_unweighted)))
  cat(sprintf("  r-bar (AMPLITUDE-WEIGHTED, for the vector mean): %s\n", fmt3(rr$r_weighted)))
  cat(sprintf("  Rayleigh Z = %s, p = %s   [Z = n * r_unweighted^2]\n",
              fmt1e(ry$Z), format.pval(ry$p, digits = 3, eps = 1e-16)))
  cat(sprintf("  Circular SD = %s h\n", fmt3(phi_to_hours(rr$circ_sd_rad, period, h))))
}

cat("\n--- Fitted model equation (pooled) ---\n")
cat(fck_format_equation(pool_int, "exp_sat", c(pool_A, pool_tau),
                        c(pool_amp1, pool_amp2), c(pool_phi1, pool_phi2),
                        period, t_off), "\n")

cat("\n--- Individual parameters (ARITHMETIC means, +/- linear SD) ---\n")
cat(sprintf("  Intercept:           %s (SD %s)   [arithmetic]\n",
            fmt3(mean(P$intercept)), fmt3(sd(P$intercept))))
cat(sprintf("  A_sat (asymptote):   %s (SD %s) KSS points   [arithmetic]\n",
            fmt3(mean(P$A_sat)), fmt3(sd(P$A_sat))))
cat(sprintf("  tau (time constant): %s (SD %s) h   [arithmetic]\n",
            fmt2(mean(P$tau)), fmt2(sd(P$tau))))
if (sd(P$tau) / mean(P$tau) > 0.5)
  cat("    ! SD/mean is above 0.5: a likelihood ridge, not population heterogeneity.\n")
for (h in 1:2) {
  rr <- if (h == 1) res1 else res2
  ac <- if (h == 1) P$amp1 else P$amp2
  cat(sprintf("  H%d amplitude: %s (SD %s)   [arithmetic]\n", h, fmt3(mean(ac)), fmt3(sd(ac))))
  cat(sprintf("  H%d acrophase: circular mean %s h, circular SD %s h   [circular]\n",
              h, fmt2(phi_to_hours(rr$mean_dir_unweighted, period, h)),
              fmt2(phi_to_hours(rr$circ_sd_rad, period, h))))
}

cat("\n--- Variance decomposition (commonality analysis) ---\n")
tot <- mean(P$r2)
cat(sprintf("  Unique to Process S (homeostatic): %s  (%s%%)\n", fmt3(uS), fmt1(100 * uS / tot)))
cat(sprintf("  Unique to Process C (circadian):   %s  (%s%%)\n", fmt3(uC), fmt1(100 * uC / tot)))
cat(sprintf("  Shared between S and C:            %s  (%s%%)\n", fmt3(sh), fmt1(100 * sh / tot)))
cat(        "  ------------------------------------------------\n")
cat(sprintf("  Total R-squared:                   %s  (%s%%)\n",
            fmt3(uS + uC + sh), fmt1(100 * (uS + uC + sh) / tot)))
cat("  No dominance verdict is emitted.\n")

cat("\n--- Group-specific parameters ---\n")
ga <- fck_group_audit(lab, seq_along(lab))
cat(sprintf("Group sizes sum to %d; %d subjects were fitted.\n", sum(ga$counts), n_fit))
if (ga$n_unassigned > 0)
  cat(sprintf("  %d subject(s) have no usable group label and appear as UNASSIGNED below.\n",
              ga$n_unassigned))
cat("  These reconcile.\n\n")
for (g in c(names(n_by_group), "__UNASSIGNED__")) {
  idx <- if (g == "__UNASSIGNED__") which(is.na(lab)) else which(!is.na(lab) & lab == g)
  if (!length(idx)) next
  gA <- mean(P$A_sat[idx]); gT <- mean(P$tau[idx]); gI <- mean(P$intercept[idx])
  gram <- fck_rhythm_adjusted_mean(gI, "exp_sat", c(gA, gT), min(tlin), max(tlin), t_off)
  rr <- fck_resultants(P$phi1[idx], P$amp1[idx])
  r2 <- fck_resultants(P$phi2[idx], P$amp2[idx])
  cat(sprintf("Group '%s' (n = %d):\n",
              if (g == "__UNASSIGNED__") "UNASSIGNED (no usable group label)" else g, length(idx)))
  if (g == "__UNASSIGNED__")
    cat("  Shown because these subjects exist. They are not a group; do not compare them.\n")
  cat(sprintf("  Intercept (beta_0, at t = 0):  %s (SD %s)   [arithmetic]\n",
              fmt3(gI), fmt3(sd(P$intercept[idx]))))
  cat(sprintf("  MESOR (rhythm-adjusted mean):  %s   [integrated over the window]\n", fmt3(gram)))
  cat(sprintf("  H1 amplitude: %s   [vector]\n", fmt3(rr$mean_amplitude)))
  cat(sprintf("                %s (SD %s)   [arithmetic]\n",
              fmt3(mean(P$amp1[idx])), fmt3(sd(P$amp1[idx]))))
  cat(sprintf("  H1 acrophase: %s h   [vector]\n",
              fmt2(phi_to_hours(rr$mean_dir_weighted, period, 1))))
  cat("\n  Fitted equation:\n  ")
  cat(fck_format_equation(gI, "exp_sat", c(gA, gT),
                          c(rr$mean_amplitude, r2$mean_amplitude),
                          c(rr$mean_dir_weighted, r2$mean_dir_weighted),
                          period, t_off), "\n\n")
}

cat("--- Group comparisons (AUDIT 2.6) ---\n")
lab_f <- factor(lab, levels = names(n_by_group))
for (nm in c("amp1", "tau")) {
  lt <- fck_group_linear_test(P[[nm]], lab_f)
  tr <- fck_group_trend_test(P[[nm]], lab_f)
  cat(sprintf("\n%s:\n", nm))
  cat(sprintf("  F(%d, %d) = %s, p = %s;  eta^2 = %s, omega^2 = %s\n",
              lt$df1, lt$df2, fmt3(lt$F), format.pval(lt$p, digits = 3, eps = 1e-16),
              fmt3(lt$eta2), fmt3(lt$omega2)))
  cat(sprintf("  Largest contrast %s - %s: %s, 95%% CI [%s, %s], Hedges' g = %s\n",
              lt$largest$a, lt$largest$b, fmt3(lt$largest$diff),
              fmt3(lt$largest$ci[1]), fmt3(lt$largest$ci[2]), fmt3(lt$largest$hedges_g)))
  if (!is.null(tr))
    cat(sprintf("  Monotone trend: t(%d) = %s, p = %s, 95%% CI [%s, %s]\n",
                tr$df, fmt3(tr$t), format.pval(tr$p, digits = 3, eps = 1e-16),
                fmt3(tr$ci[1]), fmt3(tr$ci[2])))
}
sink()

# =============================================================================
cat("\n"); line()
cat("SCORED AGAINST KNOWN TRUTH\n"); line()

w <- as.numeric(n_by_group / sum(n_by_group))
names(w) <- names(n_by_group)
true_int <- sum(w * TRUTH$intercept)
true_A   <- sum(w * TRUTH$A_sat)
true_tau <- sum(w * TRUTH$tau)
true_ram <- fck_rhythm_adjusted_mean(true_int, "exp_sat", c(true_A, true_tau),
                                     min(tlin), max(tlin), t_off)

row <- function(what, old, new, truth) {
  verdict <- if (!is.finite(old)) "new only"
  else if (abs(new - truth) < abs(old - truth) - 1e-9) "NEW closer"
  else if (abs(new - truth) > abs(old - truth) + 1e-9) "OLD closer"
  else "tie"
  cat(sprintf("%-40s %10s %10s %10s   %s\n", what,
              fmt2(old), fmt2(new), fmt2(truth), verdict))
}
cat(sprintf("%-40s %10s %10s %10s\n", "", "OLD", "NEW", "TRUTH"))
row("central value reported as the level", pool_int, ram, true_ram)
row("Rayleigh Z (H1)", nrow(P) * res1$r_weighted^2,
    fck_rayleigh(res1$r_unweighted, res1$n)$Z,
    fck_rayleigh(res1$r_unweighted, res1$n)$Z)
row("H2 acrophase (h)", phi_to_hours(pool_phi2, period, 1),
    phi_to_hours(pool_phi2, period, 2),
    (sum(w * TRUTH$phi2) * period / (2 * pi)) / 2)
row("variance parts sum to (% of total R2)",
    mean(100 * P$r2_S / P$r2) + mean(100 * P$r2_C / P$r2),
    100 * (uS + uC + sh) / tot, 100)
row("group n's sum to", sum(ga$counts), sum(ga$counts) + ga$n_unassigned, n_fit)

cat("\nPooled equation evaluated at t = 27 h (03:00 of the second day):\n")
yv <- function(withtrend) {
  v <- pool_int + pool_amp1 * cos(2 * pi * 27 / period - pool_phi1) +
    pool_amp2 * cos(2 * pi * 2 * 27 / period - pool_phi2)
  if (withtrend) v <- v + pool_A * (1 - exp(-(27 - t_off) / pool_tau))
  v
}
truth27 <- sum(w * vapply(names(n_by_group), function(g) as.numeric(curve_of(g, 27)), 1))
cat(sprintf("  OLD printed equation:  %s\n", fmt2(yv(FALSE))))
cat(sprintf("  NEW printed equation:  %s\n", fmt2(yv(TRUE))))
cat(sprintf("  Truth (mean curve):    %s\n", fmt2(truth27)))
cat(sprintf("  The old equation under-predicted by %s units.\n", fmt2(yv(TRUE) - yv(FALSE))))

line()
cat("Wrote tests/report_OLD.txt and tests/report_NEW.txt.\n")
cat("Diff them with:  diff -u tests/report_OLD.txt tests/report_NEW.txt\n")
