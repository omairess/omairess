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
  expect_equal(splines::splineDesign(c(0, x, 1), x, ord = 2), diag(n), ignore_attr = TRUE)
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
  expect_equal(t(eval.fd(x, fdo)), Yf, ignore_attr = TRUE)
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
