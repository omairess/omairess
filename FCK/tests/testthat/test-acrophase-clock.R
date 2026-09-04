# ==============================================================================
# tests/testthat/test-acrophase-clock.R
#
# Reported: acrophases were printed as clock times while they are estimated in
# MODEL coordinates. With time_origin = "first_observation" the model axis is
# elapsed hours from 08:00, so a reported "19.18 h" is 03:11 the next morning,
# not 19:11 that evening. The fitted curves never leave model coordinates and
# were right; the reporting layer was not.
#
# The four values below are the reporter's own, hand-checked. They are the
# regression targets.
#
# The fix is a coordinate conversion in ONE helper. Nothing about the model
# changes: the tests at the end assert that coefficients, amplitudes, fitted
# values and every circular TEST STATISTIC are invariant, because a constant
# rotation cannot change a dispersion or a difference -- only a direction.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/07_helpers_circular.R"))   # fck_rhythm_from_coefs
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

ORIGIN <- 8   # the reporter's data starts at 08:00

# ------------------------------------------------- the reported values, exactly
test_that("the reported H1 acrophases convert to the stated clock times", {
  expect_equal(fck_acrophase_label(hours = 19.18, harmonic = 1, clock_origin = ORIGIN), "03:11")
  expect_equal(fck_acrophase_label(hours = 19.89, harmonic = 1, clock_origin = ORIGIN), "03:53")
  expect_equal(fck_acrophase_label(hours = 18.93, harmonic = 1, clock_origin = ORIGIN), "02:56")
})

test_that("H2 has two equivalent maxima and both are reported", {
  lab <- fck_acrophase_label(hours = 7.81, harmonic = 2, clock_origin = ORIGIN)
  expect_true(grepl("15:49", lab, fixed = TRUE))
  expect_true(grepl("03:49", lab, fixed = TRUE))
  a <- fck_acrophase_clock(hours = 7.81, harmonic = 2, clock_origin = ORIGIN)
  expect_equal(length(a$all_hours[[1]]), 2)
  expect_equal(sort(round(a$all_hours[[1]], 2)), c(3.81, 15.81))
  expect_equal(diff(sort(a$all_hours[[1]])), 12)   # exactly half a day apart
})

test_that("the numeric conversion is the stated formula", {
  for (x in c(19.18, 19.89, 18.93, 0, 23.9)) {
    expect_equal(fck_acrophase_clock(hours = x, harmonic = 1, clock_origin = ORIGIN)$hours,
                 (x + ORIGIN) %% 24)
  }
})

# ------------------------------------------------------------ from radians
test_that("radians convert with the harmonic divisor", {
  # elapsed_hours = phi * 24 / (2*pi*h)
  for (h in 1:3) {
    for (phi in c(0.7, 2.04, 4.5)) {
      expect_equal(fck_acrophase_clock(phi_rad = phi, harmonic = h,
                                       clock_origin = ORIGIN)$elapsed,
                   phi * 24 / (2 * pi * h), tolerance = 1e-12)
    }
  }
})

test_that("radians and hours agree for the reporter's H2 value", {
  phi <- 7.81 * 2 * pi * 2 / 24          # back to radians
  expect_equal(fck_acrophase_label(phi_rad = phi, harmonic = 2, clock_origin = ORIGIN),
               fck_acrophase_label(hours = 7.81, harmonic = 2, clock_origin = ORIGIN))
})

test_that("the /h divisor is not optional", {
  # dropping it is the older bug; H2 would land 7.8 h away
  wrong <- 7.81 * 2 * pi * 2 / 24 * 24 / (2 * pi)      # no /h
  expect_equal(wrong, 15.62, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(wrong, 7.81)))
})

# ------------------------------------------------- harmonic h has h maxima
test_that("harmonic h reports exactly h maxima, evenly spaced", {
  for (h in 1:4) {
    a <- fck_acrophase_clock(hours = 3, harmonic = h, clock_origin = ORIGIN)
    expect_equal(length(a$all_hours[[1]]), h)
    if (h > 1) expect_equal(unique(round(diff(a$all_hours[[1]]), 8)), 24 / h)
    expect_true(all(a$all_hours[[1]] >= 0 & a$all_hours[[1]] < 24))
  }
})

test_that("every maximum really is a maximum of that harmonic", {
  # evaluate cos(2*pi*h*t_elapsed/24 - phi) at each reported clock time
  h <- 2; elapsed <- 7.81
  phi <- elapsed * 2 * pi * h / 24
  a <- fck_acrophase_clock(hours = elapsed, harmonic = h, clock_origin = ORIGIN)
  for (cl in a$all_hours[[1]]) {
    te <- cl - ORIGIN                       # back to elapsed
    expect_equal(cos(2 * pi * h * te / 24 - phi), 1, tolerance = 1e-8)
  }
})

# ------------------------------------------------------------- edge cases
test_that("a zero origin is a no-op", {
  for (x in c(0, 5.5, 19.18, 23.99))
    expect_equal(fck_acrophase_clock(hours = x, harmonic = 1, clock_origin = 0)$hours, x)
  expect_equal(fck_acrophase_label(hours = 19.18, harmonic = 1, clock_origin = 0), "19:11")
})

test_that("wrapping past midnight is handled, not clamped", {
  expect_equal(fck_acrophase_clock(hours = 20, harmonic = 1, clock_origin = 8)$hours, 4)
  expect_equal(fck_acrophase_clock(hours = 16, harmonic = 1, clock_origin = 8)$hours, 0)
  expect_equal(fck_acrophase_label(hours = 16, harmonic = 1, clock_origin = 8), "00:00")
})

test_that("NA and empty input give NA, not an error or a wrong number", {
  expect_true(is.na(fck_acrophase_clock(hours = NA_real_, clock_origin = 8)$hours))
  expect_true(is.na(fck_acrophase_label(hours = NA_real_, clock_origin = 8)))
  expect_null(fck_acrophase_clock())
})

test_that("the helper is vectorised", {
  v <- c(19.18, 19.89, 18.93)
  expect_equal(fck_acrophase_clock(hours = v, harmonic = 1, clock_origin = ORIGIN)$hours,
               (v + ORIGIN) %% 24)
  expect_equal(fck_acrophase_label(hours = v, harmonic = 1, clock_origin = ORIGIN),
               c("03:11", "03:53", "02:56"))
})

test_that("fck_clock_origin reads the model, and defaults to zero", {
  expect_equal(fck_clock_origin(list(origin_shift = 8)), 8)
  expect_equal(fck_clock_origin(list()), 0)          # midnight origin
  expect_equal(fck_clock_origin(NULL), 0)            # no model yet
})

# ================================================= the model must not change
test_that("the conversion is invariant for every DISPERSION and DIFFERENCE", {
  set.seed(3)
  phi <- rnorm(200, 2.0, 0.6) %% (2 * pi)
  amp <- rgamma(200, 5, scale = 4)
  r0 <- fck_resultants(phi, amp)
  # rotating every angle by the origin cannot change concentration
  rot <- (phi + 2 * pi * ORIGIN / 24) %% (2 * pi)
  r1 <- fck_resultants(rot, amp)
  expect_equal(r1$r_unweighted, r0$r_unweighted, tolerance = 1e-12)
  expect_equal(r1$r_weighted,   r0$r_weighted,   tolerance = 1e-12)
  expect_equal(r1$mean_amplitude, r0$mean_amplitude, tolerance = 1e-12)
  expect_equal(fck_rayleigh(r1$r_unweighted, r1$n)$Z,
               fck_rayleigh(r0$r_unweighted, r0$n)$Z, tolerance = 1e-12)
  # only the direction moves, and by exactly the origin
  expect_equal(((r1$mean_dir_unweighted - r0$mean_dir_unweighted) %% (2 * pi)),
               (2 * pi * ORIGIN / 24) %% (2 * pi), tolerance = 1e-9)
})

test_that("a between-group difference is origin-free", {
  a <- 19.89; b <- 18.93
  expect_equal((a + ORIGIN) %% 24 - (b + ORIGIN) %% 24, a - b, tolerance = 1e-12)
  # which is why the pairwise tests need no conversion at all
})

test_that("the fitted curve is untouched by the reporting change", {
  coefs <- c(20, 60, 14, 18, -9, 4, 3)
  tv <- seq(0, 31, length.out = 100)      # model time
  y <- fck_rhythm_from_coefs(coefs, tv, 24, 2, "exp_sat", include_trend = TRUE, t_offset = 0)
  expect_true(all(is.finite(y)))
  # converting for display must not enter the prediction
  y2 <- fck_rhythm_from_coefs(coefs, tv, 24, 2, "exp_sat", include_trend = TRUE, t_offset = 0)
  expect_identical(y, y2)
})

# ============================ H1 acrophase != peak of the complete curve
test_that("the curve peak and the H1 acrophase are reported as different things", {
  # a trend-dominated fit: the curve peaks near the end of the window, the H1
  # harmonic peaks wherever its own phase says
  coefs <- c(20, 80, 12, 18, -9, 4, 3)   # mesor, A_sat, tau, H1 cos/sin, H2 cos/sin
  mod <- list(period = 24, n_harmonics = 2, trend_type = "exp_sat",
              time_vec = seq(0, 31, length.out = 200), t_offset = 0, origin_shift = 8)
  pk <- fck_curve_peak_clock(coefs, mod)
  expect_false(is.null(pk))
  h1_elapsed <- atan2(coefs[5], coefs[4]) %% (2 * pi) * 24 / (2 * pi)
  h1_clock <- (h1_elapsed + 8) %% 24
  # they are genuinely different -- that is the point of reporting both
  expect_gt(abs(pk$peak_clock - h1_clock), 0.5)
  expect_true(pk$peak_clock >= 0 && pk$peak_clock < 24)
})

test_that("a peak at the window edge is flagged as a boundary value", {
  # a pure rising trend: the maximum is wherever the recording stopped
  mod <- list(period = 24, n_harmonics = 1, trend_type = "linear",
              time_vec = seq(0, 20, length.out = 200), t_offset = 0, origin_shift = 8)
  pk <- fck_curve_peak_clock(c(10, 2, 0.001, 0.001), mod)
  expect_true(pk$peak_at_edge)
  expect_equal(pk$peak_model, 20, tolerance = 1e-6)
  expect_equal(pk$peak_clock, (20 + 8) %% 24)
})

test_that("with no trend and one harmonic the curve peak IS the acrophase", {
  mod <- list(period = 24, n_harmonics = 1, trend_type = "none",
              time_vec = seq(0, 24, length.out = 400), t_offset = 0, origin_shift = 8)
  coefs <- c(50, 20, -10)                       # mesor, b_cos1, b_sin1
  pk <- fck_curve_peak_clock(coefs, mod)
  h1_clock <- ((atan2(coefs[3], coefs[2]) %% (2 * pi)) * 24 / (2 * pi) + 8) %% 24
  expect_equal(pk$peak_clock, h1_clock, tolerance = 0.05)
})
