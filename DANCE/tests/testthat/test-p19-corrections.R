# ==============================================================================
# tests/testthat/test-p19-corrections.R -- the external review of November 2026
#
#   P19.1   Curve columns read as text went into the analysis as FACTOR CODES
#           (data.matrix), with no error: "3,5", "1,0", "12,5" became 3, 1, 2.
#   P19.2   The Watson-Williams correction factor was 1 - 1/(24 kappa^2)
#           instead of 1 + 3/(8 kappa): F too small by 4% to 20%.
#   P19.3   Cyclic smoothing on the column index used a period of n - 1
#           columns (fda's default: the width of the range), so 24 hourly
#           columns were fitted with a 23-hour cycle.
#   P19.4   With real clock times the smoothing used 24 h but the fd object
#           was re-expressed with the window as its period, so the fd object
#           was not the smoothed curve.
#   P19.16  The loadings CSV used length(harmonics), which is 3 for every fd
#           object: five components exported as three, two failed.
#   P19.17  An unused smooth of the warping functions threw on fewer than 10
#           time points and the handler returned the identity silently.
#   P19.18  A periodic shift was capped at a quarter of the cycle.
#   P19.19  An amplitude below the stated minimum sent a closed-form fit to the
#           constrained optimiser, which cannot raise it, and labelled the
#           result as constrained.
#
# Each test was checked by restoring the defect and confirming it fails.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
src <- function(f) eval(parse(file.path(app_dir, "server", f), encoding = "UTF-8"),
                        envir = globalenv())
for (f in c("01_helpers_time.R", "01c_helpers_norm.R", "02e_helpers_numeric.R",
            "03_helpers_clock.R", "08_helpers_cosinor.R", "08b_helpers_popcosinor.R",
            "08c_helpers_circstat.R"))
  src(f)

# ---- P19.1 -------------------------------------------------------------------
test_that("text curve columns are read as numbers, never as factor codes", {
  d <- read.csv(text = "id;t1;t2\na;3,5;10,2\nb;1,0;7,25\nc;12,5;NA", sep = ";",
                stringsAsFactors = FALSE)
  # the defect, for the record: data.matrix() returns the factor codes
  expect_equal(unname(data.matrix(d[, 2:3])[, 1]), c(3, 1, 2))

  fx <- dance_fix_decimal_commas(d, "auto")
  expect_setequal(fx$converted, c("t1", "t2"))
  expect_equal(fx$data$t1, c(3.5, 1.0, 12.5))
  expect_equal(fx$data$t2, c(10.2, 7.25, NA))
  expect_identical(fx$data$id, d$id)                       # a label column is left alone

  # "point" switches the repair off; the matrix step then COUNTS the cells
  off <- dance_fix_decimal_commas(d, "point")
  expect_length(off$converted, 0)
  res <- dance_numeric_matrix(off$data[, c("t1", "t2")], "point")
  expect_true(all(is.na(res$matrix[, "t1"])))
  expect_equal(unname(res$n_bad), c(3L, 2L))               # the NA in t2 is a gap, not text
  expect_match(dance_numeric_matrix_message(res), "5 cells in 2 curve columns")
})

test_that("flags and gaps in a numeric column become NA and are reported", {
  d <- data.frame(t1 = c("5.1", "<LOQ", "12.3", "", "n/a"), t2 = c(1, 2, 3, 4, 5),
                  stringsAsFactors = FALSE)
  res <- dance_numeric_matrix(d)
  expect_equal(res$matrix[, "t1"], c(5.1, NA, 12.3, NA, NA))
  expect_equal(res$matrix[, "t2"], c(1, 2, 3, 4, 5))
  expect_equal(unname(res$n_bad), c(1L, 0L))               # only "<LOQ" is unreadable text
  expect_match(res$examples[["t1"]], "<LOQ", fixed = TRUE)
  expect_identical(typeof(res$matrix), "double")
  expect_null(dance_numeric_matrix_message(dance_numeric_matrix(d["t2"])))
})

test_that("a decimal-comma column that also carries a flag keeps its numbers", {
  d <- data.frame(t1 = c("3,5", "<LOQ", "12,25", ""), stringsAsFactors = FALSE)
  expect_length(dance_fix_decimal_commas(d)$converted, 0)  # the flag keeps the column text
  res <- dance_numeric_matrix(d)                           # ... and the cells are still read
  expect_equal(res$matrix[, "t1"], c(3.5, NA, 12.25, NA))
  expect_equal(unname(res$n_bad), 1L)
  expect_equal(unname(dance_numeric_matrix(d, "point")$n_bad), 3L)
})

test_that("a column that only looks comma-separated is not converted", {
  d <- data.frame(code = c("1,2", "3", "a,b"), n = c("1", "2", "3"), stringsAsFactors = FALSE)
  expect_length(dance_fix_decimal_commas(d)$converted, 0)  # "a,b" disqualifies; n has no comma
  f <- data.frame(g = factor(c("0,5", "1,5")))
  expect_equal(dance_fix_decimal_commas(f)$data$g, c(0.5, 1.5))
})

# ---- P19.2 -------------------------------------------------------------------
test_that("Watson-Williams uses the factor 1 + 3/(8 kappa)", {
  set.seed(11)
  g1 <- rnorm(20, 0.4, 0.30); g2 <- rnorm(20, 1.1, 0.30); g3 <- rnorm(20, 0.5, 0.30)
  ww <- dance_watson_williams_test(list(g1, g2, g3))
  # the statistic without its factor, rebuilt here from the definition
  R <- vapply(list(g1, g2, g3), function(x) sqrt(sum(cos(x))^2 + sum(sin(x))^2), numeric(1))
  a <- c(g1, g2, g3); Rt <- sqrt(sum(cos(a))^2 + sum(sin(a))^2); N <- 60; k <- 3
  bare <- (N - k) * (sum(R) - Rt) / ((k - 1) * (N - sum(R)))
  expect_equal(ww$F / bare, 1 + 3 / (8 * ww$kappa), tolerance = 1e-12)
  expect_gt(ww$F / bare, 1)                                # the old factor was below 1
})

test_that("Watson-Williams agrees with the circular package", {
  skip_if_not_installed("circular")
  set.seed(11)
  g1 <- rnorm(20, 0.4, 0.30); g2 <- rnorm(20, 1.1, 0.30); g3 <- rnorm(20, 0.5, 0.30)
  ref <- suppressWarnings(circular::watson.williams.test(
    circular::circular(c(g1, g2, g3)), factor(rep(1:3, each = 20))))
  ww <- dance_watson_williams_test(list(g1, g2, g3))
  expect_equal(ww$F, unname(ref$statistic), tolerance = 1e-8)
  expect_equal(ww$p, ref$p.value, tolerance = 1e-8)
})

# ---- P19.3 / P19.4 -----------------------------------------------------------
skip_if_not_installed("fda")
suppressMessages(library(fda))
src("04_helpers_fd.R")

axis_for <- function(labels, real_time, y)
  dance_smoothing_axis(list(use_real_time = real_time, is_cyclic = TRUE),
                       list(data = matrix(y, 1), time_labels = labels))

test_that("the cyclic period is one cycle, not one column short of it", {
  hourly <- sprintf("%02d:00", 0:23)
  expect_equal(dance_cyclic_period(hourly, 24, FALSE)$period, 24)
  expect_equal(dance_cyclic_period(hourly, 24, TRUE)$period, 24)
  # two-hourly columns: 12 columns to the cycle
  expect_equal(dance_cyclic_period(sprintf("%02d:00", seq(0, 22, 2)), 12, FALSE)$period, 12)
  # two days of hourly columns are two cycles of 24, not one of 48
  two_days <- sprintf("d%d_%02d:00", rep(1:2, each = 24), rep(0:23, 2))
  expect_equal(dance_cyclic_period(two_days, 48, FALSE)$period, 24)
  # no clock in the names: the columns are one cycle
  expect_equal(dance_cyclic_period(paste0("T", 101:116), 16, FALSE)$period, 16)
  # the basis carries it (it used to be 23 here)
  ax <- axis_for(hourly, FALSE, rnorm(24))
  expect_equal(dance_smoothing_basis(ax, 13, "manual")$params[1], 24)
  # an axis built by hand, without a period, gets the same rule
  ax$period <- NULL
  expect_equal(dance_smoothing_basis(ax, 13, "manual")$params[1], 24)
})

test_that("an exactly 24-hour signal is reproduced, on either axis, and by the fd object", {
  tt <- 0:23
  y <- 10 + 5 * cos(2 * pi * (tt - 3) / 24) + 2 * sin(2 * pi * tt / 24)
  for (rt in c(FALSE, TRUE)) {
    ax <- axis_for(sprintf("%02d:00", tt), rt, y)
    b <- dance_smoothing_basis(ax, 13, "manual")
    sm <- as.vector(eval.fd(ax$t_full, smooth.basis(ax$t_full, y, fdPar(b, 2, 1e-8))$fd))
    expect_lt(max(abs(sm - y)), 1e-4)                      # was 0.79 on the index axis
    # the fd object: the same cycle on 0-1, as server/20_smoothing.R builds it
    tp01 <- (ax$t_full - ax$t_rng[1]) / diff(ax$t_rng)
    b01 <- create.fourier.basis(c(0, 1), b$nbasis, period = b$params[1] / diff(ax$t_rng))
    fd01 <- smooth.basis(tp01, sm, fdPar(b01, 2, 0))$fd
    expect_lt(max(abs(as.vector(eval.fd(tp01, fd01)) - sm)), 1e-6)   # was 0.79 with real time
  }
})

test_that("the smoothing step builds its 0-1 Fourier basis with the smoothing's period", {
  # the step is dance_smooth_curves() since P20 (the app and the script call it)
  for (rt in c(FALSE, TRUE)) {
    tt <- 0:23
    y <- rbind(10 + 5 * cos(2 * pi * (tt - 3) / 24), 8 + 2 * sin(2 * pi * tt / 24))
    ax <- dance_smoothing_axis_from(sprintf("%02d:00", tt), 24, rt, cyclic = TRUE)
    b <- dance_smoothing_basis(ax, 13, "manual")
    sm <- dance_smooth_curves(y, ax, b, "manual", lambda = 1e-8)
    expect_equal(sm$fd_obj$basis$params, b$params[1] / diff(ax$t_rng))
    expect_lt(max(abs(t(eval.fd(sm$time_points_01, sm$fd_obj)) - sm$smooth_mat)), 1e-6)
  }
})

# ---- P19.16 ------------------------------------------------------------------
test_that("the loadings export writes one column per component", {
  set.seed(1); tt <- seq(0, 1, length.out = 20)
  fdo <- smooth.basis(tt, matrix(rnorm(20 * 15), 20), create.bspline.basis(c(0, 1), 8))$fd
  for (k in c(2, 5)) {
    p <- pca.fd(fdo, nharm = k)
    expect_equal(length(p$harmonics), 3)                   # why length() must not be used
    n_pc <- dance_n_harmonics(p)
    m <- matrix(as.numeric(eval.fd(seq(0, 1, length.out = 100), p$harmonics)), 100, n_pc)
    expect_equal(ncol(m), k)
    expect_equal(m[, k], as.vector(eval.fd(seq(0, 1, length.out = 100), p$harmonics[k])))
  }
  code <- paste(readLines(file.path(app_dir, "server", "90_export.R"), encoding = "UTF-8"),
                collapse = "\n")
  # the guard is on the CODE lines: comments are stripped, so describing the
  # removed expression in one cannot satisfy or defeat it
  live <- gsub("#[^\n]*", "", code)
  expect_false(grepl("length(values$pca_results$harmonics)", live, fixed = TRUE))
})

# ---- P19.17 / P19.18 ---------------------------------------------------------
src("05_helpers_warp.R")

test_that("shift registration runs on fewer than 10 time points", {
  set.seed(1)
  for (nt in c(8, 9)) {
    tt <- seq(0, 1, length.out = nt)
    sh <- seq(-0.14, 0.14, length.out = 15)
    Y <- t(sapply(sh, function(s) sin(2 * pi * (tt - s))))
    r <- linear_shift_alignment(dance_fd_no_smoothing(Y, tt), periodic = FALSE, time_points = tt)
    expect_identical(r$method, "linear_shift")             # was "identity"
    expect_null(r$failed)
    expect_gt(abs(suppressWarnings(stats::cor(sh, r$shifts))), 0.8)
  }
})

test_that("a registration that cannot run says so instead of returning a silent identity", {
  # a single evaluation point leaves the periodic estimator nothing to correlate
  tt <- seq(0, 1, length.out = 8)
  fdo <- dance_fd_no_smoothing(matrix(sin(2 * pi * tt), 3, 8, byrow = TRUE), tt)
  r <- suppressWarnings(linear_shift_alignment(fdo, periodic = TRUE, time_points = 0.5))
  expect_identical(r$method, "identity")
  expect_true(r$failed)
  expect_true(nzchar(r$error_message))
})

test_that("a periodic shift is recovered up to half a cycle", {
  tt <- seq(0, 1, length.out = 49)
  sh <- c(0, 0.1, 0.2, 0.3, 0.4, -0.3, -0.4)
  Y <- t(sapply(sh, function(s) sin(2 * pi * (tt - s)) + 0.3 * cos(4 * pi * (tt - s))))
  r <- linear_shift_alignment(dance_fd_no_smoothing(Y, tt), periodic = TRUE,
                              reference = "first", time_points = tt)
  expect_equal(abs(r$shifts), abs(sh), tolerance = 0.03)   # 0.3 and 0.4 were both 0.25
  # registered: every curve now lies on the reference
  expect_lt(max(abs(r$registered_curves - r$registered_curves[, 1])), 0.15)
  # the open interval keeps its quarter-domain limit
  r2 <- linear_shift_alignment(dance_fd_no_smoothing(Y, tt), periodic = FALSE,
                               reference = "first", time_points = tt)
  expect_lte(max(abs(r2$shifts)), 0.25 + 1e-12)
})

# ---- P19.19 ------------------------------------------------------------------
test_that("an amplitude below the stated minimum stays on the closed-form fit and is reported", {
  skip_if_not_installed("minpack.lm")
  ex <- parse(file.path(app_dir, "server", "72_harmonic.R"), encoding = "UTF-8")
  for (x in ex)
    if (is.call(x) && identical(as.character(x[[1]]), "<-") && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% c("fit_cosinor", "fit_cosinor_nonlinear"))
      eval(x, globalenv())
  set.seed(2); t <- 0:47
  y <- 10 + 3 * cos(2 * pi * (t - 3) / 24) + rnorm(48, 0, 0.3)
  free <- fit_cosinor(t, y, 24, 1, "none")
  lo <- fit_cosinor(t, y, 24, 1, "none", use_bounds = TRUE, amplitude_min = 5, amplitude_max = 10)
  expect_true(lo$success)
  expect_false(lo$bounds_active)                           # was TRUE
  expect_match(lo$fit_route, "closed-form least squares")
  expect_match(lo$fit_route, "below the stated minimum")
  expect_equal(lo$amplitude_below_min, 1L)
  expect_equal(lo$amplitudes, free$amplitudes)             # the same estimate, exactly
  expect_equal(unname(lo$coefs), unname(free$coefs))
  # nothing below the minimum: nothing reported
  ok <- fit_cosinor(t, y, 24, 1, "none", use_bounds = TRUE, amplitude_min = 1, amplitude_max = 10)
  expect_length(ok$amplitude_below_min, 0)
  expect_identical(ok$fit_route, "closed-form least squares (no bound was binding)")
  # the maximum is still a constraint, and still enforced
  hi <- fit_cosinor(t, y, 24, 1, "none", use_bounds = TRUE, amplitude_min = 0, amplitude_max = 1)
  expect_true(hi$bounds_active)
  expect_lte(hi$amplitudes[1], 1 + 1e-6)
})
