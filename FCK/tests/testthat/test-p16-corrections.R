# ==============================================================================
# tests/testthat/test-p16-corrections.R — the performance regression
#
# Reported by the user: "the performance of the app is sensible worse.
# everything is slowed down, the analyses and the GUI." Both halves were mine,
# and both were measured rather than guessed at.
#
#   P16.1 (GUI)  The resize wrapper added at P15 did a full
#                document.querySelectorAll on EVERY shiny:value -- an event that
#                fires once per output, for all 159 of them, on every render --
#                and its ResizeObserver dispatched a GLOBAL window resize, which
#                makes Shiny recompute sizes for every bound output and
#                re-render every base-R plot. ResizeObserver also fires once
#                when observation begins, so merely creating the wrappers fired
#                one global resize per plot during page load.
#
#                Measured in Chromium on a page with 40 plots and 60 text
#                outputs, before -> after:
#                    page load           301 scans, 40 resizes  ->  2 scans, 0
#                    one full refresh    100 scans              ->  0 scans
#                    dragging one plot   100 scans, 1 resize    ->  0, 0
#                and resizing still works: plotly svg 294 -> 614 px, base-R
#                image regenerated at 840x534.
#
#   P16.2 (analyses)  P12.1 removed the 60-subject cap from fck_auto_lambda,
#                correctly -- the report claimed the whole sample and did not
#                have it. But that made every automatic smoothing, and every
#                step of the n-basis sweep, proportional to the sample. The cost
#                was never necessary: smooth.basis accepts a MATRIX of curves
#                and returns one GCV per column, and the app was calling it once
#                per subject inside a 25-point lambda grid. Subjects sharing a
#                missingness pattern share everything the fit depends on, so one
#                call per pattern gives the SAME numbers. The cap stays gone.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"

raw_of <- function(f) paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")

# The pre-P16 objective, kept as the reference implementation: the point of the
# change is that it is EXACTLY equivalent, so the old code is the oracle.
old_auto_lambda <- function(data_mat, argvals, basisobj, min_points_needed = 4,
                            log_range = c(-8, 4), n_grid = 25) {
  data_mat <- as.matrix(data_mat)
  rows <- which(rowSums(!is.na(data_mat)) >= min_points_needed)
  if (!length(rows)) return(NULL)
  mean_gcv <- function(log_lambda) {
    fdp <- fda::fdPar(basisobj, 2, 10^log_lambda)
    v <- vapply(rows, function(i) {
      ok <- !is.na(data_mat[i, ])
      f <- tryCatch(fda::smooth.basis(argvals[ok], data_mat[i, ok], fdp),
                    error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      g <- suppressWarnings(as.numeric(f$gcv))
      if (!length(g) || !is.finite(g[1])) NA_real_ else g[1]
    }, numeric(1))
    if (all(is.na(v))) return(Inf)
    mean(v, na.rm = TRUE)
  }
  grid <- seq(log_range[1], log_range[2], length.out = n_grid)
  scores <- vapply(grid, mean_gcv, numeric(1))
  if (all(!is.finite(scores))) return(NULL)
  k <- which.min(scores)
  lo <- grid[max(1, k - 1)]; hi <- grid[min(length(grid), k + 1)]
  best <- if (lo < hi)
    tryCatch(stats::optimize(mean_gcv, c(lo, hi), tol = 1e-3)$minimum,
             error = function(e) grid[k]) else grid[k]
  list(lambda = 10^best, n_used = length(rows), gcv = mean_gcv(best))
}

fd_env <- function() {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = e)
  e
}

test_that("P16.2: the vectorised GCV search returns the SAME lambda", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- fd_env()
  n_time <- 24; tp <- seq_len(n_time)
  b <- create.bspline.basis(c(1, n_time), nbasis = 12)
  mk <- function(m) t(sapply(seq_len(m), function(i)
    5 * sin(2 * pi * tp / 24 + runif(1, 0, 6)) + rnorm(n_time, 0, 1.5)))

  set.seed(11); Y <- mk(40)                       # complete: one pattern
  expect_equal(e$fck_auto_lambda(Y, tp, b)$lambda, old_auto_lambda(Y, tp, b)$lambda)

  set.seed(11); Y2 <- mk(40); Y2[, 5] <- NA       # a shared gap: still one
  expect_equal(e$fck_auto_lambda(Y2, tp, b)$lambda, old_auto_lambda(Y2, tp, b)$lambda)

  set.seed(11); Y3 <- mk(40)                      # worst case: many patterns
  for (i in 1:40) Y3[i, sample(n_time, sample(1:4, 1))] <- NA
  r_new <- e$fck_auto_lambda(Y3, tp, b); r_old <- old_auto_lambda(Y3, tp, b)
  expect_equal(r_new$lambda, r_old$lambda)
  expect_equal(r_new$n_used, r_old$n_used)
  expect_equal(r_new$gcv, r_old$gcv)
})

test_that("P16.2: subjects with too few points are still excluded", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- fd_env()
  n_time <- 24; tp <- seq_len(n_time)
  b <- create.bspline.basis(c(1, n_time), nbasis = 12)
  set.seed(5)
  Y <- t(sapply(1:30, function(i) 5 * sin(2 * pi * tp / 24) + rnorm(n_time, 0, 1)))
  Y[1, ] <- NA                      # nothing observed
  Y[2, 4:n_time] <- NA              # 3 points, below the default of 4
  r <- e$fck_auto_lambda(Y, tp, b)
  expect_equal(r$n_used, 28L)
  expect_equal(r$n_used, old_auto_lambda(Y, tp, b)$n_used)
})

test_that("P16.2: the cap removed at P12.1 has NOT come back", {
  src <- raw_of("server/04_helpers_fd.R")
  expect_false(grepl("length.out = 60", src, fixed = TRUE))
  expect_false(grepl("if (length(rows) > 60)", src, fixed = TRUE))
  # the fix is the grouping, not a shortcut on the sample
  expect_true(grepl("pat_groups <- split(rows, pat)", src, fixed = TRUE))
  expect_true(grepl("t(data_mat[idx, ok, drop = FALSE])", src, fixed = TRUE))
})

test_that("P16.2: complete data costs one smooth.basis call per lambda", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- fd_env()
  n_time <- 24; tp <- seq_len(n_time)
  b <- create.bspline.basis(c(1, n_time), nbasis = 12)
  set.seed(9)
  Y <- t(sapply(1:60, function(i) 5 * sin(2 * pi * tp / 24) + rnorm(n_time, 0, 1)))
  # 60 subjects, complete: one missingness pattern, so one group
  pat <- apply(Y, 1, function(r) paste0(as.integer(!is.na(r)), collapse = ""))
  expect_equal(length(unique(pat)), 1L)
  # and the search is fast enough to be run per basis size in the sweep
  t0 <- Sys.time(); e$fck_auto_lambda(Y, tp, b)
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 3)
})

# ====================================================== P16.1 the resize JS ===
test_that("P16.1: shiny:value enhances one element, it does not rescan", {
  js <- raw_of("ui/00_theme.R")
  expect_true(grepl("jQuery(document).on('shiny:value', function (e) {", js, fixed = TRUE))
  expect_true(grepl("if (e.target && e.target.nodeType === 1) enhance(e.target);",
                    js, fixed = TRUE))
  # the old form scanned the whole document on every output paint
  expect_false(grepl("on('shiny:value shiny:visualchange', function () {", js, fixed = TRUE))
})

test_that("P16.1: a full scan is coalesced into one animation frame", {
  js <- raw_of("ui/00_theme.R")
  expect_true(grepl("requestAnimationFrame", js, fixed = TRUE))
  expect_true(grepl("if (pending) return;", js, fixed = TRUE))
  expect_true(grepl("shiny:visualchange shown.bs.tab", js, fixed = TRUE))
})

test_that("P16.1: plotly never triggers a global resize", {
  js <- raw_of("ui/00_theme.R")
  expect_true(grepl("if (g) { try { Plotly.Plots.resize(g); } catch (e) {} return; }",
                    js, fixed = TRUE))
  expect_true(grepl("if (!isRPlot) return;", js, fixed = TRUE))
  expect_true(grepl("var isRPlot = el.classList.contains('shiny-plot-output');",
                    js, fixed = TRUE))
})

test_that("P16.1: the observer's initial callback is ignored", {
  js <- raw_of("ui/00_theme.R")
  # ResizeObserver fires once on observe(); without this, creating N wrappers
  # fired N global resizes during page load
  expect_true(grepl("if (first) { first = false; return; }", js, fixed = TRUE))
  expect_true(grepl("if (t) clearTimeout(t);", js, fixed = TRUE))
})
