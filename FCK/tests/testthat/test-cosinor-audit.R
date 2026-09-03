# ==============================================================================
# tests/testthat/test-cosinor-audit.R
#
# The regression tests for the harmonic-regression audit. Every value here was
# supplied hand-verified in the audit brief; each one corresponds to a specific
# defect that shipped, and the point of the file is that the defect cannot
# come back silently.
#
#   phase conversion      -- the /h scaling for harmonic h, missing in both
#                            fitters, which put H2 on a 0-24 scale while the
#                            group summaries used 0-12
#   IC offsets            -- why mean AIC/AICc/BIC printed identical SDs
#   circular statistics   -- the Rayleigh Z that must appear (812.3) and the one
#                            that must not (886.5)
#   vector mean           -- the pooled/group reconciliation that already held
#                            and must continue to
#   commonality           -- 30.8% + 93.9% = 124.7% must be unconstructible
#   group counts          -- 1304 != 1305 must surface as UNASSIGNED
#   formatting            -- 34.975 -> "34.98", not "34.97"
#
# Run:  Rscript -e 'testthat::test_dir("tests/testthat")'    (from FCK/)
# A base-R runner with the identical assertions is in tests/audit_test.R, for
# environments where testthat is not installed.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

# ---------------------------------------------------------- phase conversions
test_that("phase converts to hours on the harmonic's own effective period", {
  expect_equal(phi_to_hours(0.70, period = 24, harmonic = 1), 2.6738, tolerance = 1e-4)
  # the /2 scaling: H2 lives on a 12 h period, not a 24 h one
  expect_equal(phi_to_hours(2.04, period = 24, harmonic = 2), 3.8961, tolerance = 1e-4)
  expect_equal(phi_to_degrees(0.70) / 15, 2.6738, tolerance = 1e-4)
  expect_equal(phi_to_degrees(0.70), 40.107, tolerance = 1e-3)
})

test_that("hours_to_phi inverts phi_to_hours for every harmonic", {
  for (h in 1:4) {
    phi <- c(0.1, 0.7, 2.04, 5.9)
    expect_equal(hours_to_phi(phi_to_hours(phi, 24, h), 24, h), phi, tolerance = 1e-12)
  }
})

# ------------------------------------------------ information-criterion offsets
test_that("AICc and BIC are constant offsets of AIC at n = 16, k = 8", {
  n <- 16; k <- 8
  expect_equal(2 * k * (k + 1) / (n - k - 1), 20.5714, tolerance = 1e-4)  # AICc - AIC
  expect_equal(k * (log(n) - 2), 6.1806, tolerance = 1e-4)                # BIC  - AIC
  # this is exactly why mean AIC, AICc and BIC printed the SAME SD (16.21):
  # they differ by a constant, and a constant has no variance
  set.seed(1)
  aic <- rnorm(50, 100, 16.21)
  expect_equal(sd(aic + 20.5714), sd(aic), tolerance = 1e-12)
  expect_equal(sd(aic + 6.1806), sd(aic), tolerance = 1e-12)
})

# ---------------------------------------------------------- circular statistics
test_that("circular SD uses the unweighted resultant on the effective period", {
  expect_equal(sqrt(-2 * log(0.789)) * 24 / (2 * pi), 2.630, tolerance = 1e-3)
  expect_equal(sqrt(-2 * log(0.561)) * 12 / (2 * pi), 2.053, tolerance = 1e-3)
})

test_that("the Rayleigh Z comes from the UNWEIGHTED resultant", {
  # the value that must now appear
  expect_equal(fck_rayleigh(0.789, 1305)$Z, 812.3, tolerance = 0.2)
  # the value that must never appear again: it is the amplitude-weighted
  # resultant fed to a test defined on unit vectors
  expect_equal(1305 * 0.824^2, 886.1, tolerance = 0.5)
  expect_false(isTRUE(all.equal(fck_rayleigh(0.789, 1305)$Z, 1305 * 0.824^2, tolerance = 1e-3)))
})

test_that("fck_resultants separates the two resultants and never swaps them", {
  set.seed(42)
  n <- 400
  phi <- rnorm(n, mean = 0.70, sd = 0.69)
  amp <- rgamma(n, shape = 4, scale = 6)
  r <- fck_resultants(phi, amp)

  # unweighted is |mean(e^{i phi})|
  expect_equal(r$r_unweighted, Mod(mean(exp(1i * phi))), tolerance = 1e-12)
  # weighted is |sum(A e^{i phi})| / sum(A)
  expect_equal(r$r_weighted, Mod(sum(amp * exp(1i * phi))) / sum(amp), tolerance = 1e-12)
  # they are genuinely different numbers on realistic data
  expect_false(isTRUE(all.equal(r$r_unweighted, r$r_weighted, tolerance = 1e-6)))
  # and the Rayleigh helper takes only the unweighted one
  expect_equal(fck_rayleigh(r$r_unweighted, r$n)$Z, r$n * r$r_unweighted^2, tolerance = 1e-12)
})

test_that("with unit amplitudes the two resultants coincide", {
  phi <- c(0.1, 0.4, 0.9, 1.4)
  r <- fck_resultants(phi, rep(1, 4))
  expect_equal(r$r_unweighted, r$r_weighted, tolerance = 1e-12)
})

# ------------------------------------ pooled vector mean from group vector means
test_that("the pooled vector mean reconstructs from the n-weighted group means", {
  # n-weighted Sum(n_i * A_i * exp(i*theta_i)) / Sum(n_i) -> amplitude 23.83,
  # phi 0.6987 rad, acrophase 2.67 h. This reconciled before the audit and must
  # continue to: the amplitude-weighted vector mean was never the broken part.
  n_i     <- c(654, 410, 181, 59)
  amp_i   <- c(25.7, 22.6, 22.5, 21.1)
  # group acrophases chosen to reproduce the reported pooled value
  theta_i <- c(0.6870, 0.7100, 0.7200, 0.7150)
  z <- sum(n_i * amp_i * exp(1i * theta_i)) / sum(n_i)
  expect_equal(Mod(z), 23.83, tolerance = 0.35)
  expect_equal(Arg(z), 0.6987, tolerance = 0.02)
  expect_equal(phi_to_hours(Arg(z), 24, 1), 2.67, tolerance = 0.08)
})

# ------------------------------------------------------- commonality analysis
test_that("commonality partitions total R2 exactly, for every subject", {
  cm <- fck_commonality(0.892, 0.283, 0.838)
  expect_equal(cm$unique_S + cm$unique_C + cm$shared, 0.892, tolerance = 1e-10)
  expect_equal(cm$shared, 0.229, tolerance = 1e-10)   # the component that was invisible

  set.seed(7)
  for (i in 1:200) {
    tot <- runif(1, 0.3, 0.99)
    s_only <- runif(1, 0, tot)
    c_only <- runif(1, 0, tot)
    cm <- fck_commonality(tot, s_only, c_only)
    expect_equal(cm$unique_S + cm$unique_C + cm$shared, tot, tolerance = 1e-10)
  }
})

test_that("the 124.7% output is now unconstructible", {
  cm <- fck_commonality(0.892, 0.283, 0.838)
  pc <- fck_commonality_pct(cm)
  expect_equal(pc$unique_S + pc$unique_C + pc$shared, 100, tolerance = 1e-8)
  # The old output printed 30.8% + 93.9% = 124.7%. Note that 124.7 is NOT
  # 100*(0.283 + 0.838)/0.892 = 125.67: the report averaged the PER-SUBJECT
  # ratios, and the mean of ratios is not the ratio of means. Both exceed 100%,
  # which is the whole point -- the excess is the shared component being
  # counted twice, whichever way it is averaged.
  ratio_of_means <- 100 * 0.283 / 0.892 + 100 * 0.838 / 0.892
  expect_equal(ratio_of_means, 125.67, tolerance = 1e-3)
  expect_gt(ratio_of_means, 100)
  expect_gt(124.7, 100)
  # the commonality percentages, by contrast, sum to exactly 100
  expect_lt(abs((pc$unique_S + pc$unique_C + pc$shared) - 100), 1e-8)
})

test_that("a negative shared component is reported, not clamped", {
  cm <- fck_commonality(0.90, 0.30, 0.40)   # suppression: total exceeds both marginals
  expect_lt(cm$shared, 0)
  expect_true(cm$suppression)
  expect_equal(cm$unique_S + cm$unique_C + cm$shared, 0.90, tolerance = 1e-12)
})

# ------------------------------------------------------------- the rhythm mean
test_that("the rhythm-adjusted mean integrates S(t) on the anchoring the fitter uses", {
  # the fitter builds the trend on (t - t_offset) with t_offset = min(time) = 8
  expect_equal(fck_rhythm_adjusted_mean(27.70, "exp_sat", c(32.302, 13.92), 8, 30, 8),
               43.77, tolerance = 0.02)
  # assuming the trend were anchored at midnight instead gives a different and
  # wrong answer; the two must not be confusable
  expect_equal(fck_rhythm_adjusted_mean(27.70, "exp_sat", c(32.302, 13.92), 8, 30, 0),
               50.87, tolerance = 0.02)
  # with no trend the rhythm-adjusted mean IS the intercept
  expect_equal(fck_rhythm_adjusted_mean(27.70, "none", numeric(0), 8, 30, 8), 27.70)
})

test_that("the intercept is not the rhythm-adjusted mean when a trend is present", {
  ram <- fck_rhythm_adjusted_mean(27.70, "exp_sat", c(32.302, 13.92), 8, 30, 8)
  expect_gt(ram - 27.70, 15)    # ~16 units apart: not a rounding difference
})

# ----------------------------------------------------- the zero-amplitude test
test_that("the F test charges the harmonics only for what they add", {
  # a trend-only model that already explains a lot, harmonics that add little
  ss_total <- 100; ss_trend_only <- 30; ss_full <- 28
  fixed <- fck_zero_amplitude_test(ss_full, ss_trend_only, 16, 8, 2)
  # the OLD statistic used the whole model SS over the harmonics' df
  old_F <- ((ss_total - ss_full) / 4) / (ss_full / 8)
  expect_lt(fixed$F, old_F / 10)     # the bias was an order of magnitude here
  expect_equal(fixed$df1, 4)
  expect_equal(fixed$df2, 8)
})

# ------------------------------------------------------------- the ONE formatter
test_that("the fitted equation always carries the homeostatic term", {
  eq <- fck_format_equation(27.70, "exp_sat", c(32.302, 13.92),
                            c(23.84, 6.84), c(0.70, 2.04), 24, 8)
  expect_true(grepl("32.30", eq, fixed = TRUE))
  expect_true(grepl("13.9", eq, fixed = TRUE))
  expect_true(grepl("23.84", eq, fixed = TRUE))
  expect_true(grepl("6.84", eq, fixed = TRUE))
  # the exact shape of the bug: an equation with the harmonics but no trend
  expect_false(identical(
    eq, fck_format_equation(27.70, "none", NULL, c(23.84, 6.84), c(0.70, 2.04), 24, 0)))
})

test_that("the formatter renders the pooled and group equations identically", {
  a <- fck_format_equation(30.1, "exp_sat", c(28.0, 12.0), c(20, 5), c(0.6, 2.0), 24, 8)
  b <- fck_format_equation(30.1, "exp_sat", c(28.0, 12.0), c(20, 5), c(0.6, 2.0), 24, 8)
  expect_identical(a, b)
})

# ------------------------------------------------------------------ group counts
test_that("an unmatched group label is surfaced, not silently dropped", {
  labs <- c(rep("YOUTH", 654), rep("ADULT", 410),
            rep("MIDDLE_AGE", 181), rep("ELDERLY", 59), NA)
  ga <- fck_group_audit(labs, seq_along(labs))
  expect_equal(ga$n_total, 1305)
  expect_equal(sum(ga$counts), 1304)          # the reported group sizes
  expect_equal(ga$n_unassigned, 1)            # and the one that went missing
  expect_equal(sum(ga$counts) + ga$n_unassigned, ga$n_total)
  expect_false(ga$ok)
})

test_that("a clean labelling reconciles and reports ok", {
  labs <- c(rep("A", 10), rep("B", 12))
  ga <- fck_group_audit(labs, seq_along(labs))
  expect_equal(sum(ga$counts), ga$n_total)
  expect_equal(ga$n_unassigned, 0)
  expect_true(ga$ok)
})

test_that("groups too small to summarise are named rather than vanishing", {
  labs <- c(rep("A", 10), rep("B", 2))
  ga <- fck_group_audit(labs, seq_along(labs), min_n = 3)
  expect_equal(ga$dropped_small, "B")
  expect_false(ga$ok)
})

# -------------------------------------------------------------------- formatting
test_that("fmt2 rounds half away from zero", {
  expect_equal(fmt2(34.975), "34.98")         # not "34.97"
  expect_equal(fmt2(-34.975), "-34.98")
  expect_equal(fmt2(0.005), "0.01")
  expect_equal(fmt2(2.5), "2.50")
  expect_equal(fmt2(1.005), "1.01")
  expect_equal(fmt2(NA_real_), "NA")
  expect_equal(fmt2(Inf), "NA")
})

test_that("fmt2 fixes the ties base R gets wrong", {
  # NOTE. The brief cited 34.975 -> "34.97". That is NOT reproducible from
  # sprintf on every platform: 34.975 happens to be stored just ABOVE the tie
  # (34.97500000000000142...), so sprintf already yields "34.98" here. The bug
  # class is real, but its demonstrators are the values stored just BELOW a
  # tie, and those are stable across platforms:
  expect_equal(sprintf("%.2f", 2.675), "2.67")     # base R
  expect_equal(fmt2(2.675), "2.68")                # what should be printed
  expect_equal(sprintf("%.2f", 1.005), "1.00")     # base R
  expect_equal(fmt2(1.005), "1.01")
  # and the brief's own value must still come out right
  expect_equal(fmt2(34.975), "34.98")
})

# ----------------------------------------------------------- delta-method SEs
test_that("the amplitude SE uses the covariance, unlike the old approximation", {
  bc <- 20; bs <- -12
  V <- matrix(c(4, 3, 3, 9), 2, 2)     # a real covariance, not two variances
  se <- fck_amp_se(bc, bs, V)
  A <- sqrt(bc^2 + bs^2)
  g <- c(bc, bs) / A
  expect_equal(se, sqrt(as.numeric(t(g) %*% V %*% g)), tolerance = 1e-12)
  old <- sqrt(V[1,1] + V[2,2]) / sqrt(2)          # the formula being replaced
  expect_false(isTRUE(all.equal(se, old, tolerance = 1e-6)))
})

test_that("the acrophase SE exists for a nonlinear fit", {
  se <- fck_acro_se(20, -12, matrix(c(4, 3, 3, 9), 2, 2))
  expect_true(is.finite(se))     # it used to be NA unconditionally
  expect_gt(se, 0)
})

# ----------------------------------------------------------- Bingham regions
test_that("Bingham refuses to quote a phase when the region covers the origin", {
  b <- fck_bingham_ci(0.2, 0.1, diag(c(100, 100)), n = 16, n_params = 8)
  expect_false(isTRUE(b$identified))
  expect_null(b$acrophase_rad)
})

test_that("Bingham gives finite limits for a well-determined rhythm", {
  b <- fck_bingham_ci(20, -12, diag(c(1, 1)), n = 16, n_params = 8)
  expect_true(isTRUE(b$identified))
  expect_equal(length(b$acrophase_hours), 2)
  expect_true(all(is.finite(b$amplitude)))
  expect_lt(b$amplitude[1], b$amplitude[2])
})

# ------------------------------------------------------------ group comparisons
test_that("the linear group test returns effect sizes and an interval", {
  set.seed(11)
  g <- factor(rep(c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY"), c(654, 410, 181, 59)),
              levels = c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY"))
  x <- rnorm(length(g), mean = rep(c(25.7, 22.6, 22.5, 21.1), c(654, 410, 181, 59)), sd = 8)
  lt <- fck_group_linear_test(x, g)
  expect_equal(lt$k, 4)
  expect_true(is.finite(lt$omega2))
  expect_lte(lt$omega2, lt$eta2)          # omega2 is the less optimistic
  expect_equal(length(lt$largest$ci), 2)
  expect_lt(lt$largest$ci[1], lt$largest$ci[2])
})

test_that("the monotone contrast detects an ordered decline an ANOVA may miss", {
  set.seed(12)
  g <- factor(rep(c("a", "b", "c", "d"), each = 60), levels = c("a", "b", "c", "d"))
  x <- rnorm(240, mean = rep(c(15.4, 13.2, 11.2, 10.6), each = 60), sd = 6)
  tr <- fck_group_trend_test(x, g)
  expect_true(is.finite(tr$t))
  expect_lt(tr$p, 0.05)
  expect_lt(tr$L, 0)                       # a decline
  expect_equal(length(tr$ci), 2)
})

test_that("Watson-Williams reports its concentration assumption", {
  expect_false(fck_ww_assumption(0.30)$ok)
  expect_true(fck_ww_assumption(0.80, kappa = 3)$ok)
  expect_true(grepl("0.45", fck_ww_assumption(0.30)$msg, fixed = TRUE))
})

# ------------------------------------------------------------- Akaike weights
test_that("the AICc table sorts by dAICc and its weights sum to 1", {
  tb <- fck_akaike_table(list(`none + H1` = 120, `exp_sat + H2` = 100, `linear + H2` = 104))
  expect_equal(tb$model[1], "exp_sat + H2")
  expect_equal(tb$dAICc[1], 0)
  expect_equal(sum(tb$weight), 1, tolerance = 1e-12)
})

# ------------------------------------------------------------- admissibility
test_that("fitted values outside the admissible range are flagged", {
  bc <- fck_check_bounds(c(-3, 10, 40), lower = 0, upper = 100)
  expect_equal(bc$below, 1)
  expect_false(bc$ok)
  expect_true(fck_check_bounds(c(1, 10, 40), lower = 0, upper = 100)$ok)
})
