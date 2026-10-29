# ==============================================================================
# tests/testthat/test-traj-levels.R
#
# The level of a mixed-effects trajectory as three numbers per design cell:
#   intercept  the midline (intercept + trend, no harmonics) at the first
#              observation -- where every trend column is zero
#   MESOR      that midline averaged over the observed window
#   change     the trend's value at the last observation: the non-periodic
#              change over the recording, the end of the "change from origin"
#              view
# All three are exact linear functions of the fixed effects, built from the
# closed-form window summaries of the trend basis, and tested with the same
# L-matrix machinery as every other contrast. These pin the closed forms, the
# identities between the three, and the fact that the change is the trend
# coefficient's own test expressed in the response's units.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

test_that("the window summaries are the exact mean and end value of each trend column", {
  T <- 25
  num_mean <- function(f) stats::integrate(f, 0, T, rel.tol = 1e-12)$value / T
  lin <- e$dance_traj_trend_window("linear", 8, 8 + T)
  expect_equal(lin$mean[["trend_lin"]], num_mean(function(t) t), tolerance = 1e-10)
  expect_equal(lin$end[["trend_lin"]], T)
  lg <- e$dance_traj_trend_window("log", 8, 8 + T)
  expect_equal(lg$mean[["trend_log"]], num_mean(function(t) log(t + 1)), tolerance = 1e-10)
  expect_equal(lg$end[["trend_log"]], log(T + 1))
  for (tau in c(3, 18, 60)) {
    sat <- e$dance_traj_trend_window("exp_sat", 8, 8 + T, tau)
    expect_equal(sat$mean[["trend_sat"]], num_mean(function(t) 1 - exp(-t / tau)),
                 tolerance = 1e-10, info = tau)
    expect_equal(sat$end[["trend_sat"]], 1 - exp(-T / tau), info = tau)
  }
  # the saturating trend has no window without its time constant
  expect_null(e$dance_traj_trend_window("exp_sat", 8, 8 + T, NULL))
  expect_null(e$dance_traj_trend_window("linear", 8, 8))
  none <- e$dance_traj_trend_window("none", 8, 8 + T)
  expect_length(none$mean, 0); expect_length(none$end, 0)
})

test_that("the window is labelled on the clock with the day of its end", {
  w <- e$dance_window_clock(0, 25, clock_origin = 8, period = 24)
  expect_equal(w$start, "08:00")
  expect_equal(w$end, "09:00 (+1d)")
  expect_equal(w$text, "08:00 to 09:00 (+1d)")
  expect_equal(e$dance_window_clock(0, 10, 8, 24)$end, "18:00")
})

mk <- function(trend = "linear", tau = NULL, nper = 10, seed = 4) {
  set.seed(seed)
  tp <- seq(0, 25, length.out = 14)
  grp <- rep(c("A", "B", "C"), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    rise <- switch(trend, linear = (0.6 + 0.3 * (g == "B")) * tp,
                   exp_sat = (18 + 6 * (g == "C")) * (1 - exp(-tp / 15)), 0 * tp)
    30 + 3 * (g == "B") + rnorm(1, 0, 2) + rise + 6 * cos(2 * pi * (tp - 8) / 24) +
      rnorm(length(tp), 0, 1.5)
  }))
  sp <- e$dance_traj_spec(e$dance_traj_long(Y, tp, sprintf("S%03d", seq_along(grp)),
                                            list(Group = grp)), 24, 1, trend, tau = tau)
  ff <- e$dance_traj_fit(sp); ff$clock_origin <- 8; ff
}

test_that("per cell: MESOR = intercept + trend window mean, change = trend at the end", {
  skip_if_not_installed("lme4")
  for (tr in c("linear", "exp_sat")) {
    ff <- mk(tr, tau = if (tr == "exp_sat") 18 else NULL); skip_if_not(isTRUE(ff$ok))
    ls <- e$dance_traj_level_summary(ff)
    expect_identical(ls$quantities, c("intercept", "mesor", "change"))
    ce <- e$dance_traj_cell_equations(ff)
    tw <- ls$window
    for (i in seq_len(nrow(ce$table))) {
      cl <- ce$table$cell[i]; c1 <- ce$trend_coefs[[i]][[1]]
      g <- function(q) ls$table$estimate[ls$table$quantity == q & ls$table$cell == cl]
      expect_equal(g("intercept"), ce$table$intercept[i], tolerance = 1e-8, info = tr)
      expect_equal(g("mesor"), ce$table$intercept[i] + c1 * tw$mean[[1]], tolerance = 1e-8, info = tr)
      expect_equal(g("change"), c1 * tw$end[[1]], tolerance = 1e-8, info = tr)
      # and the summary's own columns agree
      expect_equal(ce$table$mesor[i], g("mesor"), tolerance = 1e-8, info = tr)
      expect_equal(ce$table$change[i], g("change"), tolerance = 1e-8, info = tr)
    }
    # the MESOR is the time average of the fitted midline (the full curve minus
    # its harmonics), by brute force on a fine grid with the trapezoid rule
    pf <- e$dance_traj_predict(ff, component = "full", n_time = 4001)
    ph <- e$dance_traj_predict(ff, component = "harmonics", n_time = 4001)
    expect_true(isTRUE(pf$ok) && isTRUE(ph$ok))
    for (cl in unique(pf$table$cell)) {
      tt <- pf$table$t[pf$table$cell == cl]
      mid <- pf$table$fit[pf$table$cell == cl] - ph$table$fit[ph$table$cell == cl]
      trap <- sum(diff(tt) * (head(mid, -1) + tail(mid, -1)) / 2) / diff(range(tt))
      expect_equal(trap, ls$table$estimate[ls$table$quantity == "mesor" & ls$table$cell == cl],
                   tolerance = 1e-6, info = paste(tr, cl))
    }
    # intervals are Wald, symmetric about the estimate
    r <- ls$table[1, ]
    expect_equal(r$hi - r$estimate, r$estimate - r$lo, tolerance = 1e-10)
    expect_equal(r$hi - r$estimate, stats::qnorm(.975) * r$se, tolerance = 1e-10)
  }
})

test_that("the change is the trend coefficient's own test, in response units", {
  skip_if_not_installed("lme4")
  ff <- mk("exp_sat", tau = 18); skip_if_not(isTRUE(ff$ok))
  chg <- e$dance_traj_contrasts(ff, "change", adjust = "none", df_method = "satterthwaite")
  trd <- e$dance_traj_contrasts(ff, "trend", adjust = "none", df_method = "satterthwaite")
  expect_true(isTRUE(chg$ok) && isTRUE(trd$ok))
  end <- 1 - exp(-25 / 18)
  expect_equal(chg$table$estimate, trd$table$estimate * end, tolerance = 1e-8)
  expect_equal(chg$table$p_raw, trd$table$p_raw, tolerance = 1e-8)
  expect_equal(chg$table$statistic, trd$table$statistic, tolerance = 1e-6)
  # the unit line says so, with the window on the clock
  expect_match(chg$unit, "08:00 to 09:00 (+1d)", fixed = TRUE)
  expect_match(chg$unit, "so its test is the test of that coefficient", fixed = TRUE)
  expect_match(chg$unit, "tau held at 18 h", fixed = TRUE)
  # A itself is offered, under its own name, and says how it relates to the change
  expect_true("trend" %in% e$dance_traj_pair_whats(ff))
  expect_identical(trd$label, "Saturating rise A (the asymptote above the intercept)")
  expect_match(trd$unit, "same test and different estimates", fixed = TRUE)
})

test_that("pairwise intercept and MESOR differences are differences of the cell values", {
  skip_if_not_installed("lme4")
  ff <- mk("linear"); skip_if_not(isTRUE(ff$ok))
  ls <- e$dance_traj_level_summary(ff)
  for (q in c("intercept", "mesor")) {
    ct <- e$dance_traj_contrasts(ff, q, adjust = "holm", df_method = "satterthwaite")
    expect_true(isTRUE(ct$ok), info = q)
    expect_false(isTRUE(ct$joint), info = q)
    for (k in seq_len(nrow(ct$table))) {
      z <- ct$table[k, ]
      v <- function(cl) ls$table$estimate[ls$table$quantity == q & ls$table$cell == cl]
      expect_equal(z$estimate, v(z$cell1) - v(z$cell2), tolerance = 1e-8, info = q)
    }
    # the omnibus over all cells is one hypothesis per quantity: G - 1 contrasts
    om <- e$dance_traj_marginal_test(ff, character(0), q, df_method = "satterthwaite")
    expect_true(isTRUE(om$ok), info = q)
    expect_equal(om$df1, 2, info = q)
  }
  # the planted intercept difference (B is 3 higher) is recovered
  ct <- e$dance_traj_contrasts(ff, "intercept", adjust = "none", df_method = "satterthwaite")
  ab <- ct$table[ct$table$cell1 == "A" & ct$table$cell2 == "B", ]
  expect_lt(ab$lo, -3 + 2.5); expect_gt(ab$hi, -3 - 2.5)
})

test_that("without a trend the intercept is the MESOR, and only it is offered", {
  skip_if_not_installed("lme4")
  ff <- mk("none"); skip_if_not(isTRUE(ff$ok))
  ls <- e$dance_traj_level_summary(ff)
  expect_identical(ls$quantities, "intercept")
  expect_identical(e$DANCE_TRAJ_BLOCK_COMPONENTS(ff$spec, "mesor"), character(0))
  expect_identical(e$DANCE_TRAJ_BLOCK_COMPONENTS(ff$spec, "change"), character(0))
  expect_identical(e$dance_traj_pair_whats(ff)[1], "intercept")
})
