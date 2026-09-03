# ==============================================================================
# tests/polar_agreement_test.R — the polar tab must agree with the plots it
# claims to be a polar version of
#
# Two silent divergences shipped before this test existed, and neither showed
# as an error — both produced a plausible-looking circle with wrong numbers:
#
#   1. The fit ring dropped the trend term while labelling its radius "fitted
#      response". With a trend in the model, hovering 03:00 on the ring gave a
#      value the Fitted Curves tab never shows at 03:00; the two disagreed by
#      exactly the trend at that time.
#   2. The signal profile wrapped ELAPSED hours (fck_cumulative_hours(), which
#      starts at 0) onto the dial as if they were clock hours. A recording that
#      started at 08:00 was rotated 8 h backwards, so a 04:00 peak was drawn at
#      20:00 — a real result, pointing at the wrong time of day.
#
# Run with:   Rscript tests/polar_agreement_test.R      (from the FCK directory)
# Needs no packages.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
`%||%` <- function(a, b) if (is.null(a)) b else a
eval(parse(file.path(app_dir, "server/03_helpers_clock.R"), encoding = "UTF-8"))
eval(parse(file.path(app_dir, "server/07_helpers_circular.R"), encoding = "UTF-8"))

# fck_clock_hours() calls the ported WaPaa parser; the merged app supplies it.
# Here we only need labels of the form "KSS_8u30", so a small stand-in is
# enough and keeps this test free of the 18k-line server files.
extract_hour_from_colname <- function(nm) {
  m <- regmatches(nm, regexec("([0-9]+)u([0-9]+)?", nm))[[1]]
  if (length(m) >= 2 && nzchar(m[2]))
    as.numeric(m[2]) + (if (length(m) >= 3 && nzchar(m[3])) as.numeric(m[3]) / 60 else 0)
  else NA_real_
}

failures <- 0
check <- function(label, cond, detail = "") {
  if (isTRUE(cond)) cat(sprintf("ok  : %s\n", label))
  else { failures <<- failures + 1; cat(sprintf("FAIL: %s  %s\n", label, detail)) }
}
near <- function(a, b, tol = 1e-8) all(abs(a - b) < tol)

# ==============================================================================
# 1. fck_rhythm_from_coefs() vs predict_from_coefs()
# ==============================================================================
# The arithmetic of server/72_harmonic.R's predict_from_coefs(), transcribed.
# If the two ever drift apart, the fit plot and its polar twin show different
# curves for the same model, which is the bug this file exists to catch.
predict_from_coefs_ref <- function(coefs, time_vec, period, n_harmonics,
                                   trend_type = "none", t_offset = 0) {
  pred <- rep(coefs[1], length(time_vec))
  n_trend <- switch(as.character(trend_type),
                    "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
  off <- 1 + n_trend
  if (trend_type != "none") {
    pred <- pred + switch(as.character(trend_type),
      "linear"  = coefs[2] * time_vec,
      "log"     = coefs[2] * log(time_vec - t_offset + 1),
      "exp_sat" = coefs[2] * (1 - exp(-(time_vec - t_offset) / coefs[3])),
      rep(0, length(time_vec)))
  }
  for (h in seq_len(n_harmonics)) {
    omega <- 2 * pi * h / period
    pred <- pred + coefs[off + 2 * h - 1] * cos(omega * time_vec) +
                   coefs[off + 2 * h]     * sin(omega * time_vec)
  }
  pred
}

t_abs <- seq(8, 30, length.out = 45)   # an 08:00 -> 06:00 protocol

# --- no trend: the two must agree unconditionally ---------------------------
cf0 <- c(50, 12, -6)                   # mesor, beta_cos_1, beta_sin_1
check("no trend: polar equals the fit plot",
      near(fck_rhythm_from_coefs(cf0, t_abs, 24, 1, "none", include_trend = TRUE),
           predict_from_coefs_ref(cf0, t_abs, 24, 1, "none")))
check("no trend: include_trend makes no difference",
      near(fck_rhythm_from_coefs(cf0, t_abs, 24, 1, "none", include_trend = TRUE),
           fck_rhythm_from_coefs(cf0, t_abs, 24, 1, "none", include_trend = FALSE)))

# --- linear trend, the case that broke --------------------------------------
cf1 <- c(20, 1.45, 12, -6)             # mesor, slope, beta_cos_1, beta_sin_1
ref <- predict_from_coefs_ref(cf1, t_abs, 24, 1, "linear")
check("linear trend: include_trend = TRUE reproduces the fit plot",
      near(fck_rhythm_from_coefs(cf1, t_abs, 24, 1, "linear", include_trend = TRUE), ref))

no_trend <- fck_rhythm_from_coefs(cf1, t_abs, 24, 1, "linear", include_trend = FALSE)
check("linear trend: include_trend = FALSE differs by exactly the trend",
      near(ref - no_trend, cf1[2] * t_abs))
# and that difference is not cosmetic: at 03:00 of the second day it is ~39 units
i3 <- which.min(abs(t_abs - 27))
check("the gap at 03:00 next day is the whole disagreement seen on screen",
      abs((ref - no_trend)[i3] - cf1[2] * t_abs[i3]) < 1e-8 &&
        (ref - no_trend)[i3] > 30,
      sprintf("gap %.2f", (ref - no_trend)[i3]))

# --- harmonic indexing must survive the trend coefficients ------------------
# Passing trend_type = "none" to drop a trend would read the SLOPE as the first
# cosine coefficient. The periodic part must be identical either way.
check("dropping the trend does not shift the harmonic coefficients",
      near(no_trend, predict_from_coefs_ref(c(cf1[1], 0, cf1[3], cf1[4]),
                                            t_abs, 24, 1, "linear")))

# --- two harmonics and a saturating trend -----------------------------------
cf2 <- c(15, 40, 6, 10, -4, 3, 2)      # mesor, A_sat, tau, then 2 harmonics
check("exp_sat trend, 2 harmonics: matches the fit plot",
      near(fck_rhythm_from_coefs(cf2, t_abs, 24, 2, "exp_sat",
                                 include_trend = TRUE, t_offset = 8),
           predict_from_coefs_ref(cf2, t_abs, 24, 2, "exp_sat", t_offset = 8)))

# --- the periodic part really is periodic, the fit is not --------------------
check("the trend-free curve closes: value at 00:00 equals value at 24:00",
      near(fck_rhythm_from_coefs(cf1, 0, 24, 1, "linear"),
           fck_rhythm_from_coefs(cf1, 24, 24, 1, "linear")))
check("the fit does NOT close, which is why it is drawn as an open arc",
      !near(fck_rhythm_from_coefs(cf1, 8, 24, 1, "linear", include_trend = TRUE),
            fck_rhythm_from_coefs(cf1, 32, 24, 1, "linear", include_trend = TRUE)))

# ==============================================================================
# 2. open paths must not be sorted
# ==============================================================================
# An 08:00 -> 06:00 arc runs 8, 9, ... 23, 0, ... 6 in clock hours. Sorting it
# would reorder it to 0 ... 23 and draw the curve backwards through itself.
arc_h <- c(20, 22, 23, 0, 2, 4)
arc_v <- c(60, 70, 75, 80, 86, 90)
op <- fck_open_path(arc_h, arc_v)
check("fck_open_path keeps the order it was given", identical(op$hours, arc_h))
check("fck_open_path does not repeat a point to close a seam",
      length(op$hours) == length(arc_h))
cr <- fck_close_ring(arc_h, arc_v)
check("fck_close_ring, by contrast, does sort and close",
      cr$hours[1] == 0 && cr$hours[length(cr$hours)] == 0)

bp <- fck_band_path(arc_h, arc_v - 5, arc_v + 5)
check("fck_band_path runs out along the upper edge in order",
      identical(bp$hours[seq_along(arc_h)], arc_h))
check("fck_band_path returns along the lower edge reversed",
      identical(bp$values[length(arc_h) + seq_along(arc_h)], rev(arc_v - 5)))

# ==============================================================================
# 3. the signal profile must sit at real clock times
# ==============================================================================
labels <- sprintf("KSS_%du%02d", c(8:23, 0:6), 0)   # 08:00 -> 06:00, hourly
clock  <- fck_clock_hours(labels)
cum    <- fck_cumulative_hours(labels)

check("the labels parse to clock hours starting at 08:00",
      !is.null(clock) && near(clock[1], 8) && near(clock[length(clock)], 6))
check("cumulative hours are ELAPSED and start at 0, not at 08:00",
      !is.null(cum) && near(cum[1], 0) && near(cum[length(cum)], 22))
check("the two are NOT interchangeable — that was the bug",
      !near(cum %% 24, clock))
check("the offset between them is exactly the start hour",
      near((cum + clock[1]) %% 24, clock))

# A peak at the last column (04:00-ish) must not be reported in the evening.
peak_col <- which.max(seq_along(labels))          # last column = 06:00
check("the last column is 06:00 on the clock, not 22:00",
      near(clock[peak_col], 6),
      sprintf("clock %.1f, elapsed-wrapped %.1f", clock[peak_col], cum[peak_col] %% 24))

# The profile ring builder must use the clock hours, not the elapsed ones.
src <- paste(readLines(file.path(app_dir, "server/74_polar_density.R"), warn = FALSE),
             collapse = "\n")
prof <- sub(".*fck_profile_rings <- function", "", src)
prof <- sub("output\\$harmonic_density_plot.*", "", prof)
check("fck_profile_rings takes its clock from fck_clock_hours()",
      grepl("clock <- fck_clock_hours", prof, fixed = TRUE))
check("fck_profile_rings no longer wraps cum straight onto the dial",
      !grepl("w <- fck_wrap_to_clock(cum", prof, fixed = TRUE))

# ==============================================================================
# 4. the fit tab must say which of the two curves is on screen
# ==============================================================================
check("the tab offers the choice between the fit and the rhythm",
      grepl("density_fit_scope", src, fixed = TRUE))
check("the radial axis is relabelled when the trend is removed",
      grepl("fitted response, trend removed", src, fixed = TRUE))

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nPolar agreement tests passed.\n")
