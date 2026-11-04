# ==============================================================================
# tests/testthat/test-fanova-missing-group.R
#
# A curve whose group is missing used to enter the functional ANOVA: the
# observed F skipped it, but every permuted statistic came out undefined, so
# all pointwise p-values were missing, nothing could be significant, and the
# L2 p was (1 + 0)/(1 + 0) = 1. The curve is now dropped and counted, and the
# result is the one the complete curves give. No permutation, no p-value.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
skip_if_not_installed("fda")
suppressPackageStartupMessages(library(fda))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
e <- new.env()
for (f in c("01c_helpers_norm", "06_helpers_posthoc")) sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)
for (x in parse(file.path(app_dir, "server", "50_fanova.R"), encoding = "UTF-8"))
  if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
      as.character(x[[2]]) %in% c("perform_functional_anova", "perform_pairwise_comparisons")) eval(x, e)

set.seed(1)
tg <- seq(0, 1, length.out = 15)
g <- factor(rep(1:2, c(12, 17)))
mk <- function(i) 1 + 2.5 * exp(-(tg / 0.12)^2) + ifelse(g[i] == "1", 0.75, 0) * exp(-((tg - 0.33) / 0.06)^2) +
  rnorm(15, 0, 0.35)
Y <- sapply(seq_along(g), mk)
basis <- create.bspline.basis(c(0, 1), 9)

test_that("no valid permutation gives no p-value, not p = 1", {
  expect_true(is.na(e$dance_perm_p(0, 0)))
  expect_equal(e$dance_perm_p(c(0, 3), c(0, 9)), c(NA, 0.4))
  expect_equal(e$dance_perm_p(0, 49), 1 / 50)
})

test_that("a curve with no group is dropped and counted, and changes nothing else", {
  fd <- smooth.basis(tg, Y, basis)$fd
  fd2 <- smooth.basis(tg, cbind(Y, mk(1)), basis)$fd
  set.seed(9); r <- e$perform_functional_anova(fd, g, n_permutations = 500)
  set.seed(9); r2 <- e$perform_functional_anova(fd2, factor(c(as.character(g), NA)), n_permutations = 500)
  expect_equal(r2$n_dropped, 1L)
  expect_equal(r$n_dropped, 0L)
  expect_equal(r2$group_sizes, c(12, 17))
  expect_equal(r2$df_within, 27)
  expect_equal(r2$F_stat, r$F_stat)
  expect_equal(r2$p_value_L2, r$p_value_L2)
  expect_equal(r2$p_values_pointwise, r$p_values_pointwise)
  expect_false(anyNA(r2$p_values_pointwise))
  expect_lt(r2$p_value_L2, 0.05)
  # a blank label is missing too
  set.seed(9); r3 <- e$perform_functional_anova(fd2, factor(c(as.character(g), " ")), n_permutations = 500)
  expect_equal(r3$n_dropped, 1L)
  expect_equal(r3$p_value_L2, r$p_value_L2)
})

test_that("the shaded regions are the FDR-significant stretches", {
  sh <- e$dance_fanova_sig_shapes(c(FALSE, TRUE, TRUE, FALSE, TRUE), c(0, .25, .5, .75, 1))
  expect_length(sh, 2)
  expect_equal(c(sh[[1]]$x0, sh[[1]]$x1, sh[[2]]$x0, sh[[2]]$x1), c(.25, .5, 1, 1))
  expect_length(e$dance_fanova_sig_shapes(rep(FALSE, 3), 1:3), 0)
})
