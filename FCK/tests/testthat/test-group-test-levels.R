# ==============================================================================
# tests/testthat/test-group-test-levels.R
#
# Two user-reported defects in the harmonic tab.
#
# 1. "Error in t.test.formula: grouping factor must have exactly 2 levels"
#
#    The report chose between t-test and ANOVA on length(unique(group)) -- the
#    number of VALUES present -- while t.test.formula() branches on nlevels(),
#    the number of LEVELS the factor carries. A factor keeps its levels after
#    the rows holding them are filtered away, so a 4-level variable reduced to
#    2 surviving values took the t-test branch and t.test then saw 4 levels and
#    stopped.
#
#    It became reachable when the convergence gate started excluding
#    non-converged fits: on the Circaflex data that removes 938 of 1305
#    subjects and a sparse level can lose every member. These tests pin the
#    distinction that made it possible, so no future filter can reintroduce it.
#
# 2. The fit plot's x axis and hover reported linear model time (27.01508)
#    rather than clock time (03:01). The model must be fitted on linear time --
#    a non-periodic trend needs 08:00 on two days to be two different t -- so
#    the conversion is display-only, and these pin the formatter.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

# --------------------------------------------- the level/value distinction
test_that("a filtered factor keeps levels that no row uses any more", {
  g <- factor(c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")[c(1, 1, 2, 2, 3, 4)],
              levels = c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY"))
  keep <- g %in% c("YOUTH", "ADULT")
  g2 <- g[keep]
  # this is the trap, stated as an assertion
  expect_equal(length(unique(g2)), 2)   # what the old code counted
  expect_equal(nlevels(g2), 4)          # what t.test.formula counts
  expect_equal(nlevels(droplevels(g2)), 2)
})

test_that("an undropped factor is NOT what breaks t.test", {
  # Recorded because it is the obvious wrong diagnosis, and it is wrong.
  # t.test.formula() calls factor() on the grouping column, which drops unused
  # levels itself, so a 4-level factor holding only 2 values is fine.
  # Verified on R 4.3.3.
  set.seed(1)
  g <- factor(rep(c("A", "B", "C", "D"), each = 10),
              levels = c("A", "B", "C", "D"))
  y <- rnorm(40)
  keep <- g %in% c("A", "B")
  d <- data.frame(y = y[keep], g = g[keep])
  expect_equal(nlevels(d$g), 4)
  expect_silent(tt <- stats::t.test(y ~ g, data = d))
  expect_true(is.finite(tt$p.value))
})

test_that("a SUBSET that loses a group IS what breaks it", {
  # The real mechanism: n_groups computed once on the parent frame, then used
  # to pick the test for a subset that no longer has that many groups.
  set.seed(1)
  lv <- c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")
  g <- factor(rep(lv, c(6, 4, 2, 2)), levels = lv)
  A_sat <- rnorm(14)
  A_sat[g != "YOUTH"] <- NA          # only one group has a usable trend fit
  parent <- data.frame(A_sat = A_sat, group = g)
  expect_equal(length(unique(parent$group)), 4)
  sub <- parent[!is.na(parent$A_sat), ]
  expect_equal(length(unique(droplevels(sub$group))), 1)
  # the user's error, reproduced exactly
  expect_error(stats::t.test(A_sat ~ group, data = sub), "exactly 2 levels")
})

test_that("NA in the group column is counted as a group by unique()", {
  # A second route to the same mismatch: length(unique()) counts NA as a value,
  # so one real group plus NAs reads as two and takes the t-test branch.
  gna <- factor(c(rep("YOUTH", 6), rep(NA, 4)), levels = c("YOUTH", "ADULT"))
  expect_equal(length(unique(gna)), 2)                      # what was counted
  expect_equal(nlevels(droplevels(gna[!is.na(gna)])), 1)    # what is there
})

test_that("droplevels makes the count and the test agree for every subset", {
  set.seed(2)
  lv <- c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")
  g <- factor(sample(lv, 200, TRUE, prob = c(.5, .3, .15, .05)), levels = lv)
  y <- rnorm(200)
  # every possible surviving subset of levels
  for (k in 2:4) for (combo in combn(lv, k, simplify = FALSE)) {
    d <- data.frame(y = y[g %in% combo], g = droplevels(g[g %in% combo]))
    expect_equal(nlevels(d$g), length(unique(d$g)))
    if (nlevels(d$g) == 2) expect_silent(stats::t.test(y ~ g, data = d))
    else expect_silent(stats::aov(y ~ g, data = d))
  }
})

test_that("a subset that loses a level mid-report is the real failure mode", {
  # params_trend <- params[!is.na(params$A_sat), ] can drop a whole group even
  # when the parent frame had all four. Deciding the test from the PARENT's
  # level count is what broke; deciding it from the subset in hand does not.
  set.seed(3)
  lv <- c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")
  g <- factor(rep(lv, c(60, 40, 18, 6)), levels = lv)
  A_sat <- rnorm(124)
  A_sat[g %in% c("MIDDLE_AGE", "ELDERLY")] <- NA   # both sparse groups fail to fit
  parent <- data.frame(A_sat = A_sat, group = g)
  n_parent <- nlevels(parent$group)
  sub <- parent[!is.na(parent$A_sat), ]
  sub$group <- droplevels(sub$group)
  expect_equal(n_parent, 4)          # the parent still says 4
  expect_equal(nlevels(sub$group), 2)  # the subset is a two-group problem
  expect_silent(stats::t.test(A_sat ~ group, data = sub))
})

test_that("fewer than two surviving groups is reported, not crashed on", {
  g <- droplevels(factor(rep("YOUTH", 30), levels = c("YOUTH", "ADULT")))
  expect_equal(nlevels(g), 1)
  # fck_group_linear_test is the effect-size path and must decline politely
  expect_null(fck_group_linear_test(rnorm(30), g))
  expect_null(fck_group_trend_test(rnorm(30), g))
})

# ------------------------------------------------------- clock-time display
test_that("linear model time renders as clock time", {
  expect_equal(fck_clock_label(8, 24), "08:00")
  expect_equal(fck_clock_label(12.5, 24), "12:30")
  expect_equal(fck_clock_label(23.99, 24), "23:59")
  # The day suffix is relative to the SERIES, so it appears when the vector
  # spans more than one day -- which is how the plot calls it, passing a whole
  # trace at once. A lone value has no series to be relative to and gets none.
  expect_equal(fck_clock_label(27.01508, 24), "03:01")
  # the value from the user's hover, formatted the way the plot formats it
  tf <- seq(8, 30, length.out = 200)
  hov <- fck_clock_label(tf, 24)
  expect_equal(hov[which.min(abs(tf - 27.015))], "03:01 (+1d)")
  expect_equal(hov[1], "08:00")
  expect_equal(hov[length(hov)], "06:00 (+1d)")
})

test_that("the day marker appears only when the axis actually spans two days", {
  within_day <- fck_clock_label(c(8, 12, 18), 24)
  expect_false(any(grepl("+1d", within_day, fixed = TRUE)))
  across <- fck_clock_label(c(8, 24, 30), 24)
  expect_true(any(grepl("+1d", across, fixed = TRUE)))
  # 08:00 on day 1 and on day 2 must be distinguishable
  two <- fck_clock_label(c(8, 32), 24)
  expect_false(two[1] == two[2])
})

test_that("minutes roll into hours instead of printing :60", {
  expect_equal(fck_clock_label(8.999, 24), "09:00")     # not "08:60"
  expect_equal(fck_clock_label(23.999, 24), "00:00")   # no series, no suffix
  expect_equal(fck_clock_label(c(8, 23.999), 24)[2], "00:00 (+1d)")
  expect_false(any(grepl(":60", fck_clock_label(seq(0, 48, by = 1/601), 24), fixed = TRUE)))
})

test_that("the old label formula produced the malformed strings this replaces", {
  # paste0(t %% period, ":00") on a half-past tick
  expect_equal(paste0(12.5 %% 24, ":00"), "12.5:00")    # what was printed
  expect_equal(fck_clock_label(12.5, 24), "12:30")      # what is printed now
  # and it had no zero padding
  expect_equal(paste0(8 %% 24, ":00"), "8:00")
  expect_equal(fck_clock_label(8, 24, with_minutes = FALSE), "08:00")
})

test_that("ticks land on clock-friendly positions across the recording", {
  tk <- fck_clock_ticks(c(8, 30), 24)
  expect_equal(tk$step, 2)
  expect_equal(tk$vals[1], 8)
  expect_equal(tk$text[1], "08:00")
  expect_true("00:00" %in% tk$text)      # the midnight crossing is labelled
  expect_equal(tk$text[length(tk$text)], "06:00")
  expect_true(all(diff(tk$vals) == tk$step))
})

test_that("tick generation copes with short and degenerate ranges", {
  expect_null(fck_clock_ticks(c(5, 5), 24))
  expect_null(fck_clock_ticks(c(NA, 30), 24))
  short <- fck_clock_ticks(c(8, 11), 24)
  expect_false(is.null(short))
  expect_true(all(short$vals >= 8 & short$vals <= 11))
})

test_that("a non-24 period still labels sensibly", {
  tk <- fck_clock_ticks(c(0, 12), 12)
  expect_false(is.null(tk))
  expect_true(all(nchar(tk$text) >= 5))
})

# ------------------------------------------------ the axis is display only
test_that("clock labelling never changes the underlying linear coordinate", {
  t_lin <- c(8, 20, 24, 27.015, 30)
  lab <- fck_clock_label(t_lin, 24)
  tk <- fck_clock_ticks(range(t_lin), 24)
  # the tick POSITIONS stay on the model's own axis; only the text is clock
  expect_true(all(tk$vals >= min(t_lin) & tk$vals <= max(t_lin)))
  expect_true(is.numeric(tk$vals))
  expect_true(is.character(tk$text))
  # and the trend term, which is why linear time is required, is monotone in it
  S <- 1 - exp(-(t_lin - 8) / 18)
  expect_true(all(diff(S) > 0))
  # whereas on a wrapped axis it would not be
  S_wrapped <- 1 - exp(-((t_lin %% 24) - 8) / 18)
  expect_false(all(diff(S_wrapped) > 0))
})
