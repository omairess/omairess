# ==============================================================================
# tests/testthat/test-warp-stats.R -- how much warping, and is it needed?
#   * the non-periodic linear shift moves curves the right way, by the right
#     amount, also for curves that rise across the recording;
#   * the shape-invariant shift model estimates the spread of the shifts and
#     tests it against none;
#   * the amplitude/phase decomposition is fda's, and refused where undefined;
#   * each participant's warp, and the group comparison.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
skip_if_not_installed("fda"); skip_if_not_installed("nlme")
suppressPackageStartupMessages(library(fda))
e <- new.env(parent = globalenv())
for (f in c("04_helpers_fd.R", "05_helpers_warp.R", "05b_helpers_warpstats.R", "08c_helpers_circstat.R"))
  sys.source(file.path(app_dir, "server", f), envir = e)

test_that("the non-periodic shift aligns curves (it used to move them the wrong way)", {
  tp <- seq(0, 1, length.out = 41); peak <- function(c) exp(-((tp - c) / 0.08)^2)
  fdo <- e$dance_fd_no_smoothing(rbind(peak(0.40), peak(0.50), peak(0.60)), tp)
  w <- e$linear_shift_alignment(fdo, periodic = FALSE, reference = "first", time_points = tp)
  expect_equal(w$shifts, c(0, -0.1, -0.2), tolerance = 1e-12)
  expect_equal(apply(w$registered_curves, 2, function(v) tp[which.max(v)]), rep(0.4, 3))
  # a rising curve two steps late is found two steps late (ccf() said 0)
  t18 <- 1:18; x <- seq(0, 1, length.out = 18); sg <- function(d) 10 + 50 * plogis((t18 - 9 - d) / 1.5)
  w2 <- e$linear_shift_alignment(e$dance_fd_no_smoothing(rbind(sg(0), sg(2)), x),
                                 periodic = FALSE, reference = "first", time_points = x)
  expect_equal(w2$shifts * 17, c(0, -2), tolerance = 1e-12)
  est <- e$dance_warp_estimates(w2, c("a", "b"), span = 17, unit = "h")
  expect_equal(est$estimate, c(0, 2), tolerance = 1e-12)       # positive = later
})

set.seed(2); t18 <- 1:18; grp <- rep(c("A", "B"), each = 12)
fcurve <- function(x) 10 + 50 * plogis((x - 9) / 1.5)
s_true <- rnorm(24, ifelse(grp == "B", 1.5, 0), 1)
Y <- t(sapply(1:24, function(i) rnorm(1, 0, 5) + fcurve(t18 - s_true[i]) + rnorm(18, 0, 5)))
rownames(Y) <- sprintf("S%02d", 1:24)

test_that("the shift model estimates the spread of the shifts and tests it", {
  sm <- e$dance_shift_model(Y, t18, df = 6)
  expect_true(sm$ok)
  expect_gt(sm$sd_shift, 1.0); expect_lt(sm$sd_shift, 1.8)
  expect_true(sm$sd_shift_ci[1] < sm$sd_shift && sm$sd_shift_ci[2] > sm$sd_shift)
  expect_lt(sm$lrt$p, 0.001)
  expect_equal(sm$lrt$statistic, 2 * unname(diff(rev(sm$logLik))), tolerance = 1e-8)
  expect_gt(cor(sm$shifts$shift, s_true), 0.9)
  expect_match(e$dance_shift_model_text(sm), "Is warping needed?", fixed = TRUE)
  # no timing differences: not significant (here)
  set.seed(9)
  Y0 <- t(sapply(1:24, function(i) rnorm(1, 0, 5) + fcurve(t18) + rnorm(18, 0, 8)))
  rownames(Y0) <- rownames(Y)
  sm0 <- e$dance_shift_model(Y0, t18, df = 6)
  expect_true(!isTRUE(sm0$ok) || sm0$lrt$p > 0.05)
})

test_that("the amplitude/phase decomposition is fda's, and only where it is defined", {
  x <- seq(0, 1, length.out = 18)
  fdo <- e$dance_fd_no_smoothing(Y, x)
  ws <- e$linear_shift_alignment(fdo, periodic = FALSE, time_points = x)
  expect_false(e$dance_amp_phase(fdo, ws)$ok)
  wp <- e$parametric_alignment(fdo, family = "power", time_points = x)
  ap <- e$dance_amp_phase(fdo, wp)
  expect_true(ap$ok)
  ref <- fda::AmpPhaseDecomp(fdo, wp$regfd, e$dance_refd(wp$warp_functions, x, NULL), c(0, 1))
  expect_equal(ap$RSQR, ref$RSQR); expect_equal(ap$C, ref$C)
  wc <- e$linear_shift_alignment(fdo, periodic = TRUE, time_points = x)
  apc <- e$dance_amp_phase(fdo, wc)
  expect_true(apc$ok); expect_equal(apc$C, 1)
})

test_that("each participant's warp, compared between groups", {
  x <- seq(0, 1, length.out = 18)
  w <- e$linear_shift_alignment(e$dance_fd_no_smoothing(Y, x), periodic = FALSE, time_points = x)
  est <- e$dance_warp_estimates(w, rownames(Y), span = 17, unit = "h")
  expect_gt(cor(est$estimate, s_true), 0.85)
  gt <- e$dance_warp_group_test(est$estimate, grp)
  expect_equal(gt$tests$p[1], stats::oneway.test(est$estimate ~ factor(grp))$p.value)
  expect_equal(gt$tests$p[2], stats::kruskal.test(est$estimate ~ factor(grp))$p.value)
  expect_lt(gt$tests$p[1], 0.05)
  expect_gt(diff(gt$groups$mean), 0)                             # B later
  gc <- e$dance_warp_group_test(c(0.1, 0.12, 0.09, 0.4, 0.42, 0.38), rep(c("A", "B"), each = 3),
                                circular = TRUE, cycle = 1)
  expect_identical(gc$tests$test, "Watson-Williams (circular)")
})
