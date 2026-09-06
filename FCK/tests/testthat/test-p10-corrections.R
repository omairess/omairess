# ==============================================================================
# tests/testthat/test-p10-corrections.R
#
# A seventh review. Two findings, both of which held.
#
#   P10.1  the periodic shift estimator was off by exactly one FFT bin, and was
#          fed a duplicated endpoint. UPSTREAM of the table P9.1 had just
#          fixed -- the wrong shift is APPLIED, so the per-subject R-squared,
#          RMSE and correlation were computed on a curve that had been moved
#          when it should not have been.
#   P10.2  production and the CV still capped the basis COUNT differently, so
#          P9.3's claim that the CV fits "exactly the production basis" was not
#          universally true.
#
# The end-to-end estimator test the reviewer asked for is
# tests/periodic_shift_test.R, which drives linear_shift_alignment() on curves
# with known displacement. These are the unit checks.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"

code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  strip <- function(l) {
    ch <- strsplit(l, "", fixed = TRUE)[[1]]
    if (!length(ch)) return(l)
    inq <- ""; out <- character(0)
    for (i in seq_along(ch)) {
      c1 <- ch[i]; esc <- i > 1 && ch[i - 1] == "\\"
      if (nzchar(inq)) { if (c1 == inq && !esc) inq <- "" }
      else if (c1 %in% c("\"", "'")) inq <- c1
      else if (c1 == "#") break
      out <- c(out, c1)
    }
    paste(out, collapse = "")
  }
  paste(vapply(ln, strip, character(1), USE.NAMES = FALSE), collapse = "\n")
}

# ================================================= P10.1 the FFT bin offset ===
test_that("P10.1: which.max() is 1-based and FFT bin 1 is zero lag", {
  # the arithmetic, in isolation, so the reason is pinned and not just the fix
  set.seed(101)
  n <- 100; tp <- seq(0, 1, length.out = n)
  f <- function(t) exp(cos(2 * pi * t))
  ref <- f(tp)

  cc <- Re(fft(Conj(fft(ref - mean(ref))) * fft(ref - mean(ref)), inverse = TRUE)) / n
  # the maximum of a signal's circular autocorrelation is at ZERO lag, which R
  # stores in element 1
  expect_equal(which.max(cc), 1L)
  # so taking which.max() AS the lag is an off-by-one of exactly one grid step
  expect_equal(-which.max(cc) / n, -0.01)
  expect_equal(-(which.max(cc) - 1L) / n, 0)
})

test_that("P10.1: dropping the duplicated endpoint is required for a cycle", {
  n <- 100; tp <- seq(0, 1, length.out = n)
  # seq(0, 1, length.out = n) contains BOTH endpoints, which are the same phase
  expect_equal(tp[1], 0); expect_equal(tp[n], 1)
  # so the circular grid has n - 1 distinct phases, and the lag resolution is
  # 1/(n-1), not 1/n
  expect_equal(length(seq_len(n - 1L)), 99L)
  expect_false(isTRUE(all.equal(1 / 99, 1 / 100)))
})

test_that("P10.1: the shipped estimator drops the endpoint and un-biases the lag", {
  src <- code_of("server/40_fpca.R")
  expect_true(grepl("circ_idx <- seq_len(n_time - 1L)", src, fixed = TRUE))
  expect_true(grepl("lag_idx <- which.max(cross_corr) - 1L", src, fixed = TRUE))
  expect_true(grepl("shifts[i] <- -lag_idx / n_circ", src, fixed = TRUE))
  # the two forms it replaces
  expect_false(grepl("max_idx <- which.max(cross_corr)", src, fixed = TRUE))
  expect_false(grepl("shifts[i] <- -shift_idx / n_time", src, fixed = TRUE))
})

test_that("P10.1: the corrected estimator recovers a known periodic shift", {
  n <- 100; tp <- seq(0, 1, length.out = n)
  f <- function(t) exp(cos(2 * pi * t))
  ref <- f(tp)
  est <- function(y) {                      # as shipped
    ci <- seq_len(n - 1L); nc <- length(ci)
    yc <- y[ci]; rc <- ref[ci]
    cc <- Re(fft(Conj(fft(rc - mean(rc))) * fft(yc - mean(yc)), inverse = TRUE)) / nc
    lag <- which.max(cc) - 1L
    if (lag > nc / 2) lag <- lag - nc
    -lag / nc
  }
  expect_equal(est(f(tp)), 0)                                  # identical curves
  for (d in c(0.05, 0.10, -0.05, -0.10))
    expect_equal(-est(f(tp - d)), d, tolerance = 1.5 / (n - 1))
  # and the old form is wrong by exactly one step, at every displacement
  old <- function(y) {
    cc <- Re(fft(Conj(fft(ref - mean(ref))) * fft(y - mean(y)), inverse = TRUE)) / n
    mi <- which.max(cc)
    -(if (mi > n / 2) mi - n else mi) / n
  }
  expect_equal(old(f(tp)), -0.01)
})

# ============================================== P10.2 one basis-count rule ====
test_that("P10.2: production and the CV use one basis-count rule", {
  suppressPackageStartupMessages(library(fda))
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)

  # the rule itself, unchanged from production
  inp <- function(m, nb = 20) list(smooth_method = m, n_basis = nb, n_basis_manual = nb)
  expect_equal(fck_smoothing_nbasis(inp("auto"),   24), 20L)
  expect_equal(fck_smoothing_nbasis(inp("manual"), 24), 20L)
  expect_equal(fck_smoothing_nbasis(inp("none"),   24), 20L)
  expect_equal(fck_smoothing_nbasis(inp("auto"),   16), 16L)   # capped at n_time
  expect_equal(fck_smoothing_nbasis(inp("none"),   16), 14L)   # n_time - 2
  expect_equal(fck_smoothing_nbasis(inp("auto", 2), 24), 4L)   # floor of 4
  expect_equal(fck_smoothing_nbasis(list(), 24), 20L)          # missing input

  # the reviewer's example: 16 time points, 20 requested
  expect_equal(fck_smoothing_nbasis(inp("auto"), 16), 16L)
  expect_equal(max(4, min(20L, 16L - 2L)), 14L)                # the old CV rule
  expect_false(fck_smoothing_nbasis(inp("auto"), 16) == 14L)
})

test_that("P10.2: neither file keeps a private copy of the rule", {
  sm <- code_of("server/20_smoothing.R")
  expect_true(grepl("nb <- fck_smoothing_nbasis(input, n_time)", sm, fixed = TRUE))
  expect_false(grepl("nb <- min(20, n_time - 2)", sm, fixed = TRUE))
  expect_false(grepl("nb <- min(nb, n_time)", sm, fixed = TRUE))
  dg <- code_of("server/30_diagnostics.R")
  expect_true(grepl("nb <- fck_smoothing_nbasis(input, n_time)", dg, fixed = TRUE))
  expect_false(grepl("nb <- max(4, min(as.integer(nb_user), n_time - 2))", dg, fixed = TRUE))
  # and the basis_label no longer carries a third copy
  expect_true(grepl("cv_axis, fck_smoothing_nbasis(input, n_time),", dg, fixed = TRUE))
})
