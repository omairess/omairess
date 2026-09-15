# ==============================================================================
# tests/testthat/test-emmeans-df-limit.R
#
# Reported from a real run:
#
#   Note: D.f. calculations have been disabled because the number of
#   observations exceeds 3000. To enable adjustments, add the argument
#   'pbkrtest.limit = 24776' (or larger) ...
#
# emmeans defaults pbkrtest.limit and lmerTest.limit to 3000 OBSERVATIONS and
# silently drops to asymptotic inference above that. Two things were wrong with
# leaving it alone: the decision reached the console, which nobody running a
# Shiny app reads, and the message's own advice -- raise the limit to the size
# of the data -- would have pbkrtest::vcovAdj() build matrices of order n_obs
# on a model with ~20 variance components.
#
# So the mode is now requested deliberately and reported. These tests pin that
# the switch happens at the limit, that the note travels with the numbers, and
# -- the one that matters -- that asking for it changes NOTHING numerically,
# because emmeans was already asymptotic there. A disclosure fix that moved an
# estimate would be a different and much worse change.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk_fit <- function(nsub, seed = 7) {
  set.seed(seed)
  tp <- seq(0, 23, length.out = 19)
  grp <- rep(c("g1", "g2"), each = nsub / 2)
  Y <- t(sapply(seq_len(nsub), function(i)
    10 + rnorm(1, 0, 2) + (4 + rnorm(1, 0, .8)) * cos(2 * pi * tp / 24 - 2) +
      rnorm(length(tp), 0, 1.2)))
  e$dance_traj_fit(e$dance_traj_spec(
    e$dance_traj_long(Y, tp, sprintf("S%04d", seq_len(nsub)), list(Group = grp)),
    24, 1, "none"))
}

test_that("the mode follows emmeans' own limit, not a number hard-coded here", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  small <- mk_fit(100)                       # 1900 observations
  big   <- mk_fit(200)                       # 3800 observations
  lim <- min(emmeans::get_emm_option("pbkrtest.limit"),
             emmeans::get_emm_option("lmerTest.limit"))
  ds <- e$dance_emm_df_mode(small); db <- e$dance_emm_df_mode(big)
  expect_lt(ds$n_obs, lim); expect_gt(db$n_obs, lim)
  expect_false(ds$over);    expect_true(db$over)
  expect_null(ds$mode);     expect_equal(db$mode, "asymptotic")
  expect_null(ds$note);     expect_true(nzchar(db$note))
  # the limit is READ from emmeans, so raising it moves the switch
  old <- emmeans::get_emm_option("pbkrtest.limit")
  on.exit(emmeans::emm_options(pbkrtest.limit = old, lmerTest.limit = old), add = TRUE)
  emmeans::emm_options(pbkrtest.limit = 1e6, lmerTest.limit = 1e6)
  expect_false(e$dance_emm_df_mode(big)$over)
})

test_that("asking for asymptotic changes nothing numerically above the limit", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  big <- mk_fit(200)
  specs <- ~ Group
  # what the app now does ...
  new <- e$dance_emtrends(big, specs, "c1")
  # ... against the untouched call, which emmeans was already running asymptotically
  old <- suppressMessages(emmeans::emtrends(big$model, specs = specs, var = "c1"))
  expect_equal(summary(new)$c1.trend, summary(old)$c1.trend, tolerance = 1e-12)
  expect_equal(as.matrix(stats::vcov(new)), as.matrix(stats::vcov(old)), tolerance = 1e-12)
  expect_equal(new@linfct, old@linfct, tolerance = 1e-12)
})

test_that("the console note is gone and the screen note is there instead", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  big <- mk_fit(200)
  msgs <- character(0)
  co <- withCallingHandlers(
    e$dance_traj_cell_coefs(big, 1),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(isTRUE(co$ok))
  expect_length(grep("D\\.f\\. calculations have been disabled|pbkrtest\\.limit", msgs), 0L)
  expect_equal(co$df_mode, "asymptotic")
  expect_true(grepl("asymptotic", co$df_note))
  # and it does NOT repeat the emmeans advice to raise the limit
  expect_true(grepl("not recommended", co$df_note))
})

test_that("a fit under the limit keeps the adjusted df it is entitled to", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  small <- mk_fit(100)
  co <- e$dance_traj_cell_coefs(small, 1)
  expect_equal(co$df_mode, "adjusted")
  expect_null(co$df_note)
})
