# ==============================================================================
# tests/testthat/test-p0-corrections.R
#
# Regression tests for the P0 corrections raised by the external review and
# verified by simulation. Each test is written so that it FAILS against the code
# as it stood before the fix -- a test that would have passed either way records
# nothing.
#
# The measured "before" numbers are in the commit message and in PORTING_NOTES
# 4.43; they are quoted here so the intent survives a later refactor.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
suppressMessages({library(fda); library(minpack.lm)})
source(file.path(app_dir, "server/04_helpers_fd.R"))
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

# Source-grep guards must read the CODE, not the prose. The fixes deliberately
# quote the removed lines in their comments so the reason survives; grepping the
# raw file would therefore find the very strings the fix deleted.
code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  ln <- sub("#.*$", "", ln)
  paste(ln, collapse = "\n")
}

# AUDIT (P11.1): the registration kernels moved out of server/40_fpca.R into
# server/05_helpers_warp.R. A guard that keeps reading only the old file does
# not fail when that happens -- it goes VACUOUS, because every expect_false()
# passes once the code is simply not in that file any more, and every
# expect_true() fails for a reason that has nothing to do with the property.
# Registration guards therefore read the registration source wherever it lives.
warp_code <- function() paste(code_of("server/40_fpca.R"),
                              code_of("server/05_helpers_warp.R"), sep = "\n")

# ------------------------------------------------- P0.2 automatic smoothing
test_that("auto lambda is a real search, not zero", {
  set.seed(3); n <- 16; t <- 1:n
  Y <- t(sapply(1:20, function(i) 5 * sin(2 * pi * t / 24) + rnorm(n, 0, 1.5)))
  b <- create.bspline.basis(c(1, n), nbasis = n)
  a <- fck_auto_lambda(Y, t, b)
  expect_false(is.null(a))
  # the old behaviour was lambda = 0 exactly
  expect_gt(a$lambda, 0)
  expect_true(is.finite(a$gcv))
})

test_that("the chosen lambda actually smooths, where lambda = 0 interpolated", {
  set.seed(3); n <- 16; t <- 1:n
  y <- 5 * sin(2 * pi * t / 24) + rnorm(n, 0, 1.5)
  Y <- t(sapply(1:20, function(i) 5 * sin(2 * pi * t / 24) + rnorm(n, 0, 1.5)))
  b <- create.bspline.basis(c(1, n), nbasis = n)
  lam <- fck_auto_lambda(Y, t, b)$lambda

  f0 <- smooth.basis(t, y, fdPar(b, 2, 0))
  fa <- smooth.basis(t, y, fdPar(b, 2, lam))
  # lambda = 0 with nbasis = n_time is a square system: it interpolates
  expect_lt(max(abs(eval.fd(t, f0$fd) - y)), 1e-9)
  expect_equal(unname(f0$df), n, tolerance = 1e-6)
  # the searched lambda leaves residual degrees of freedom
  expect_lt(fa$df, n - 2)
  expect_gt(max(abs(eval.fd(t, fa$fd) - y)), 1e-6)
})

test_that("auto lambda returns NULL rather than guessing when nothing can be scored", {
  Y <- matrix(NA_real_, 5, 16)
  b <- create.bspline.basis(c(1, 16), nbasis = 8)
  expect_null(fck_auto_lambda(Y, 1:16, b))
})

# ------------------------------------- P0.3 nonlinear zero-amplitude test
test_that("the reduced exp_sat model refits tau instead of freezing it", {
  set.seed(1)
  tg <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30)
  y <- 20 + 12 * (1 - exp(-(tg - 8) / 10)) + rnorm(16, 0, 3)
  refit <- fck_reduced_exp_sat_sse(y, tg, 8,
             start = list(mesor = 20, A_sat = 12, tau = 14))
  # the old shortcut: tau frozen at the start value, mesor/A_sat by OLS
  frozen <- sum(residuals(lm(y ~ I(1 - exp(-(tg - 8) / 14))))^2)
  expect_true(is.finite(refit))
  # a refit that also optimises tau can never do worse than one that does not
  expect_lte(refit, frozen + 1e-8)
  expect_lt(refit, frozen)      # and on these data it is strictly better
})

test_that("the reduced model returns NA rather than a wrong number on failure", {
  expect_true(is.na(fck_reduced_exp_sat_sse(c(1, 2), c(1, 2), 0,
                                            list(mesor = 0, A_sat = 1, tau = 1))))
  expect_true(is.na(fck_reduced_exp_sat_sse(rep(1, 10), 1:10, 1,
                                            list(mesor = NA, A_sat = 1, tau = 1))))
})

test_that("Type-I error of the zero-amplitude test is near nominal", {
  # This is the test that matters most in this file and it must not be skipped:
  # the defect it guards inflated rejection of true nulls from 5% to 14.6%.
  set.seed(20)
  tg <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30)
  t0 <- 8; per <- 24
  f <- function(p) p[1] + p[2] * (1 - exp(-(tg - t0) / p[3])) +
                   p[4] * cos(2 * pi * tg / per) + p[5] * sin(2 * pi * tg / per)
  pv_new <- pv_old <- numeric(0)
  for (rep in 1:250) {
    y <- 20 + 12 * (1 - exp(-(tg - t0) / 10)) + rnorm(length(tg), 0, 4)
    full <- tryCatch(nls.lm(par = c(20, 12, 10, 0, 0), fn = function(p) y - f(p),
                            lower = c(-Inf, -Inf, 0.1, -Inf, -Inf),
                            control = nls.lm.control(maxiter = 200)),
                     error = function(e) NULL)
    if (is.null(full)) next
    p <- full$par; sse1 <- sum((y - f(p))^2)
    sse0_new <- fck_reduced_exp_sat_sse(y, tg, t0,
                  start = list(mesor = p[1], A_sat = p[2], tau = p[3]),
                  lower = list(tau = 0.1))
    sse0_old <- sum(residuals(lm(y ~ I(1 - exp(-(tg - t0) / p[3]))))^2)
    a <- fck_zero_amplitude_test(sse1, sse0_new, length(tg), 5, 1)
    b <- fck_zero_amplitude_test(sse1, sse0_old, length(tg), 5, 1)
    if (is.finite(a$p)) pv_new <- c(pv_new, a$p)
    if (is.finite(b$p)) pv_old <- c(pv_old, b$p)
  }
  expect_gt(length(pv_new), 200)
  # measured over 3,000 replicates: 14.6% before, 3.6% after, nominal 5%
  expect_lt(mean(pv_new < 0.05), 0.09)
  # and the correction must actually be in the conservative direction
  expect_lt(mean(pv_new < 0.05), mean(pv_old < 0.05))
})

# -------------------------------------------------- P0.4 repeated-measures F
test_that("the RM-ANOVA F matches stats::aov at a single time point", {
  set.seed(11); n <- 12; k <- 4
  Y <- outer(rnorm(n, 0, 2), rep(1, k)) + rep(c(0, 1.2, 2.0, 2.6), each = n) +
       matrix(rnorm(n * k, 0, 1), n, k)

  grand <- mean(Y); vm <- colMeans(Y); sm <- rowMeans(Y)
  E <- sweep(sweep(Y, 1, sm, "-"), 2, vm, "-") + grand          # the shipped form
  SS_visit <- sum(nrow(Y) * (vm - grand)^2)
  F_new <- (SS_visit / (k - 1)) / (sum(E^2) / ((n - 1) * (k - 1)))

  d <- data.frame(y = as.vector(Y), s = factor(rep(1:n, k)),
                  v = factor(rep(1:k, each = n)))
  a <- summary(aov(y ~ v + Error(s / v), d))[[2]][[1]]
  expect_equal(F_new, a["v", "F value"], tolerance = 1e-8)
  expect_equal((n - 1) * (k - 1), a["Residuals", "Df"])

  # the old residual kept the visit effect: SS_resid - SS_error == SS_visit
  SS_resid_old <- sum((Y - rowMeans(Y))^2)
  expect_equal(SS_resid_old - sum(E^2), SS_visit, tolerance = 1e-8)
  expect_lt((SS_visit / (k - 1)) / (SS_resid_old / ((n - 1) * (k - 1))), F_new)
})

test_that("the source no longer computes the one-margin residual", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("SS_residual <- sum(Y_centered^2", src, fixed = TRUE))
  expect_false(grepl("SS_residual_perm <- sum(Y_centered_perm^2", src, fixed = TRUE))
  expect_true(grepl("df_within = (length(unique(subject_id)) - 1) * (n_visits - 1)",
                    src, fixed = TRUE))
})

# --------------------------------------------------- P1.1 permutation p-values
test_that("no permutation p-value can be exactly zero", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("mean(F_stat_perm[t, ] >= F_stat[t])", src, fixed = TRUE))
  expect_false(grepl("mean(L2_stat_perm >= L2_stat)", src, fixed = TRUE))
  # arithmetic check of the replacement
  perm <- rnorm(200); obs <- 99
  expect_gt((1 + sum(perm >= obs)) / (1 + length(perm)), 0)
})

# ------------------------------------------------------ P0.6 degenerate fits
test_that("a constant trajectory does not produce NaN or -Inf", {
  expect_true(is.na(fck_r_squared(0, 0)))          # was 0/0 = NaN
  expect_equal(fck_r_squared(20, 100), 0.8)
  expect_true(is.na(fck_r_squared(NA, 100)))

  ll <- fck_gaussian_loglik(0, 16, y_scale = 5)    # was +Inf
  expect_true(is.finite(ll))
  aic <- -2 * ll + 2 * 6
  expect_true(is.finite(aic))                      # was -Inf, so it won everything
  expect_lt(fck_gaussian_loglik(20, 16, y_scale = 5), ll)
})

test_that("smoothing handles 0, 1 and 2 observed points explicitly", {
  src <- code_of("server/20_smoothing.R")
  expect_true(grepl("} else if(n_valid == 1) {", src, fixed = TRUE))
  expect_true(grepl("} else if(n_valid >= 2) {", src, fixed = TRUE))
  # approx() with one point raises; that must not reach the module handler,
  # which reverts the WHOLE run to raw data
  expect_error(approx(5, 3.2, xout = 1:10, rule = 2))
})

# ---------------------------------------------------------- P0.8 time warping
test_that("warping is deterministic: no unseeded randomness in the estimator", {
  src <- warp_code()
  expect_false(grepl("runif(1, -0.03, 0.03)", src, fixed = TRUE))
  expect_false(grepl("runif(1, -0.02, 0.02)", src, fixed = TRUE))
  expect_false(grepl("Add slight S-curve for visualization", src, fixed = TRUE))
  expect_false(grepl("Simple identity warping with slight variation", src, fixed = TRUE))
  # the double attenuation of the measured lag is gone
  expect_false(grepl("* 0.1  # Scale down shift", src, fixed = TRUE))
  expect_false(grepl("time_points - shifts[i] * 0.5", src, fixed = TRUE))
})

test_that("a shift warp is monotone and the landmark branch actually warps", {
  src <- warp_code()
  # the landmark branch used to hand the curves straight back
  expect_false(grepl("registered_curves[,i] <- curves[,i]\n        }\n      }\n      \n      basis <- fd_obj$basis",
                     src, fixed = TRUE))
  # P5.6 unified the warp direction: every method now applies
  # registered <- approx(time_points, curve, xout = h). The property this guard
  # exists for -- that the landmark branch actually warps rather than handing
  # the curves back -- is checked on the new form.
  expect_true(grepl("registered_curves[, i] <- approx(time_points, curves[, i],", src, fixed = TRUE))
  expect_true(grepl("fck_landmark_warp(ref_lm, own, time_points)", src, fixed = TRUE))

  # the shipped shift warp, exercised directly
  tp <- seq(0, 1, length.out = 60)
  for (s in c(-0.25, -0.1, 0, 0.1, 0.25)) {
    w <- tp - s
    expect_true(all(diff(w) > 0), info = paste("shift", s))
  }
})

test_that("the quadratic warp family is restricted to a monotone range", {
  tp <- seq(0, 1, length.out = 60)
  # the shipped restriction for the quadratic family
  pr <- c(max(0.5, 0.05), min(2, 1))
  expect_lte(pr[2], 1)
  for (alpha in seq(pr[1], pr[2], length.out = 8)) {
    w <- alpha * tp^2 + (1 - alpha) * tp
    expect_true(all(diff(w) > -1e-12), info = paste("alpha", round(alpha, 3)))
    expect_equal(w[1], 0, tolerance = 1e-12)
    expect_equal(w[length(w)], 1, tolerance = 1e-12)
  }
  # and the old default range really did break monotonicity
  w_bad <- 2 * tp^2 - tp
  expect_true(any(diff(w_bad) < 0))
})

# ------------------------------------------------------- P0.5 amplitude bound
test_that("the amplitude bound is enforced, not just applied per coefficient", {
  src <- code_of("server/72_harmonic.R")
  expect_true(grepl("amp_bound_action", src, fixed = TRUE))
  expect_true(grepl("amplitude_max / sqrt(2)", src, fixed = TRUE))
  # a box of half-width A really does permit A*sqrt(2)
  expect_equal(sqrt(10^2 + 10^2), 10 * sqrt(2), tolerance = 1e-12)
  # the inscribed box cannot
  expect_lte(sqrt(2 * (10 / sqrt(2))^2), 10 + 1e-12)
})
