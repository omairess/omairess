# ==============================================================================
# tests/testthat/test-fd-no-smoothing.R -- "no smoothing" is the data points
# joined by straight lines
#
# The fPCA family needs an fd object even when the user asks for no smoothing.
# It was a cubic B-spline with nbasis = n_time on fda's evenly spaced knots: it
# passed through every point but swung far outside the data between them, so
# the fPCA mean curve did not look like the mean of the raw data (18 hourly VAS
# sleepiness means between 10 and 65 gave a mean curve from -2 to 79). It is
# now an order-2 basis with a knot at every observation time.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"

# the raw means of the reported data set, 1h00 to 18h00
y_mean <- c(17, 13, 10, 12, 10, 12, 11, 27, 38, 50, 51, 52, 51, 59, 49, 55, 56, 65)

test_that("the old cubic representation is what produced the swings (base R B-splines)", {
  n <- length(y_mean); x <- seq(0, 1, length.out = n)
  knots <- c(0, 0, 0, seq(0, 1, length.out = n - 2), 1, 1, 1)   # fda's nbasis = n, norder = 4
  cf <- solve(splines::splineDesign(knots, x, ord = 4), y_mean)
  g <- seq(0, 1, length.out = 681)
  f <- splines::splineDesign(knots, g, ord = 4) %*% cf
  expect_lt(min(f), 0); expect_gt(max(f), 75)
  # order 2 with a knot at every observation: each function is 1 at its own
  # time point and 0 at the others, so the coefficients are the data
  expect_equal(unname(splines::splineDesign(c(0, x, 1), x, ord = 2)), diag(n))
})

test_that("no smoothing: the curves and their mean are the data joined by straight lines", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)
  set.seed(4)
  n <- length(y_mean); x <- seq(0, 1, length.out = n)
  Y <- t(vapply(1:20, function(i) pmax(0, y_mean + rnorm(n, 0, 15)), numeric(n)))
  Y[3, 5] <- NA; Y[7, c(1, 2)] <- NA; Y[9, n] <- NA
  fdo <- dance_fd_no_smoothing(Y, x)
  expect_equal(fdo$basis$nbasis, n)
  # at the observations: the data (a gap read off the line between neighbours)
  Yf <- Y
  Yf[3, 5] <- mean(Y[3, c(4, 6)]); Yf[7, 1:2] <- Y[7, 3]; Yf[9, n] <- Y[9, n - 1]
  expect_equal(unname(t(eval.fd(x, fdo))), unname(Yf))
  # between them: straight lines, so the mean curve is the line through the
  # column means and never leaves their range
  g <- seq(0, 1, length.out = 681)
  m <- as.numeric(eval.fd(g, mean.fd(fdo)))
  expect_equal(m, stats::approx(x, colMeans(Yf), g)$y, tolerance = 1e-10)
  expect_gte(min(m), min(colMeans(Yf)) - 1e-10); expect_lte(max(m), max(colMeans(Yf)) + 1e-10)
  # fPCA runs on it, and its mean is the same curve
  pc <- pca.fd(fdo, nharm = 3)
  expect_equal(as.numeric(eval.fd(g, pc$meanfd)), m, tolerance = 1e-8)
})

test_that("warped curves of an unsmoothed run go back onto straight lines", {
  # without fda: stand-ins record what the helper builds
  e <- new.env()
  e$`%||%` <- function(a, b) if (is.null(a)) b else a
  e$create.bspline.basis <- function(rangeval, breaks, norder)
    list(type = "bspline", rangeval = rangeval, nbasis = length(breaks) + norder - 2,
         params = breaks[-c(1, length(breaks))])
  e$fd <- function(coef, basisobj) list(coefs = coef, basis = basisobj)
  e$smooth.basis <- function(...) stop("a smoothing fit was attempted")
  sys.source(file.path(app_dir, "server/04_helpers_fd.R"), envir = e)
  x <- seq(0, 1, length.out = 18)
  line_basis <- e$dance_interpolating_basis(x)
  expect_true(e$dance_is_straight_line_basis(line_basis))
  expect_false(e$dance_is_straight_line_basis(list(type = "bspline", nbasis = 18, params = 1:14)))
  expect_false(e$dance_is_straight_line_basis(list(type = "fourier", nbasis = 13, params = 1)))
  R <- matrix(seq_len(18 * 3), 18, 3)            # time x curves, as the warp kernels hold them
  r <- e$dance_refd(R, x, line_basis)
  expect_equal(unname(r$coefs), unname(R))
  expect_true(e$dance_is_straight_line_basis(r$basis))
  # a smoothing basis is projected onto, as before
  expect_error(e$dance_refd(R, x, list(type = "bspline", nbasis = 12, params = 1:8)),
               "a smoothing fit was attempted")
})

test_that("time-warped curves of an unsmoothed run are straight lines through the registered values", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  e <- new.env(parent = globalenv())
  sys.source(file.path(app_dir, "server/04_helpers_fd.R"), envir = e)
  sys.source(file.path(app_dir, "server/05_helpers_warp.R"), envir = e)
  set.seed(8)
  x <- seq(0, 1, length.out = 18)
  Y <- t(vapply(1:12, function(i) 30 + 20 * sin(2 * pi * (x - runif(1, -0.05, 0.05))) + rnorm(18, 0, 3),
                numeric(18)))
  fdo <- e$dance_fd_no_smoothing(Y, x)
  w <- e$linear_shift_alignment(fdo, time_points = x)
  expect_true(e$dance_is_straight_line_basis(w$regfd$basis))
  expect_equal(unname(eval.fd(x, w$regfd)), unname(as.matrix(w$registered_curves)), tolerance = 1e-10)
})
