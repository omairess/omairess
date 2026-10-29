# Raw data drawn as lines: one path per curve through its observed values, in
# time order, broken where a design time point is missing -- never bridged.
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
e <- new.env(); sys.source(file.path(app_dir, "server/08_helpers_cosinor.R"), envir = e)

test_that("a missing value breaks the curve's line instead of bridging it", {
  t <- rep(c(0, 2, 4, 6, 8), 2); y <- c(1, 2, NA, 4, 5,  3, NA, NA, 6, 7)
  p <- e$dance_raw_paths(t, y, rep(c("a", "b"), each = 5))
  expect_equal(p$x, c(0, 2, NA, 6, 8, NA, 0, NA, 6, 8, NA))
  expect_equal(p$y, c(1, 2, NA, 4, 5, NA, 3, NA, 6, 7, NA))
})

test_that("a row missing from long data still breaks the line (the grid knows the time)", {
  p <- e$dance_raw_paths(c(0, 2, 6, 8), c(1, 2, 4, 5), rep("a", 4), grid = c(0, 2, 4, 6, 8))
  expect_equal(p$x, c(0, 2, NA, 6, 8, NA))
})

test_that("complete curves are joined point to point, in time order", {
  p <- e$dance_raw_paths(c(4, 0, 2), c(3, 1, 2), rep("a", 3))
  expect_equal(p$x, c(0, 2, 4, NA)); expect_equal(p$y, c(1, 2, 3, NA))
})
