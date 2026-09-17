# ==============================================================================
# tests/testthat/test-traj-cell-equations.R
#
# The summary panel and the publication report print a fitted equation per
# design cell under the mixed-effects approach. That equation must be the curve
# the fitted-curves tab and the comparison tab DRAW -- the same fixed effects,
# read once. dance_traj_cell_equations() is the one place the numbers come
# from; this checks them against dance_traj_predict() (the drawn curve) and
# dance_traj_cell_coefs() (the emmeans-extracted rhythm coefficients).
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk <- function(nper = 12, seed = 31, trend = "linear", K = 2) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  grp <- rep(c("A", "B"), each = nper)
  amp <- ifelse(grp == "A", 9, 6); ph <- ifelse(grp == "A", 15, 18)
  Y <- t(sapply(seq_along(grp), function(i)
    30 + rnorm(1, 0, 2) + 0.2 * tp +
      (amp[i] + rnorm(1, 0, 0.6)) * cos(2 * pi * (tp - ph[i]) / 24) +
      1.5 * cos(2 * pi * (tp - 4) / 12) + rnorm(length(tp), 0, 1.2)))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, K, trend)
  list(fit = e$dance_traj_fit(sp), grp = grp, tp = tp)
}

test_that("the printed cell equation IS the drawn cell curve", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk(); ff <- m$fit; skip_if_not(isTRUE(ff$ok), ff$message)
  ce <- e$dance_traj_cell_equations(ff)
  expect_true(isTRUE(ce$ok))
  expect_setequal(ce$table$cell, c("A", "B"))
  pr <- e$dance_traj_predict(ff, component = "full", n_time = 50)
  skip_if_not(isTRUE(pr$ok))
  t0 <- ff$spec$t0
  for (i in seq_len(nrow(ce$table))) {
    r <- ce$table[i, ]
    d <- pr$table[pr$table$cell == r$cell, ]
    tt <- d$t - t0
    y <- r$intercept + r$trend_lin * tt
    for (h in 1:2) {
      w <- 2 * pi * h / 24
      y <- y + r[[paste0("beta_cos_", h)]] * cos(w * tt) + r[[paste0("beta_sin_", h)]] * sin(w * tt)
    }
    expect_equal(y, d$fit, tolerance = 1e-8, info = r$cell)
    # and the amplitude/acrophase form of the same equation evaluates identically
    y2 <- r$intercept + r$trend_lin * tt
    for (h in 1:2)
      y2 <- y2 + r[[paste0("amplitude_", h)]] * cos(2 * pi * h * tt / 24 - r[[paste0("acrophase_rad_", h)]])
    expect_equal(y2, d$fit, tolerance = 1e-8, info = paste(r$cell, "polar form"))
    # the level at t0 is the curve at t0
    expect_equal(r$level_at_t0, d$fit[which.min(abs(d$t - t0))], tolerance = 1e-6)
  }
})

test_that("the cell (cos, sin) pairs agree with the emmeans extraction", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk(); ff <- m$fit; skip_if_not(isTRUE(ff$ok), ff$message)
  ce <- e$dance_traj_cell_equations(ff)
  for (h in 1:2) {
    co <- e$dance_traj_cell_coefs(ff, h); skip_if_not(isTRUE(co$ok), co$message)
    for (k in seq_along(co$cells)) {
      r <- ce$table[ce$table$cell == co$cells[k], ]
      expect_equal(r[[paste0("beta_cos_", h)]], unname(co$a[k]), tolerance = 1e-6)
      expect_equal(r[[paste0("beta_sin_", h)]], unname(co$b[k]), tolerance = 1e-6)
    }
  }
  # the planted contrast survives: A has the larger H1 amplitude, B the later phase
  A <- ce$table[ce$table$cell == "A", ]; B <- ce$table[ce$table$cell == "B", ]
  expect_gt(A$amplitude_1, B$amplitude_1)
  expect_gt(B$acrophase_time_1, A$acrophase_time_1)
  # participants n per cell is carried
  expect_equal(sort(ce$table$n_participants), c(12L, 12L))
})

test_that("with no design factor there is one 'whole sample' equation", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  set.seed(5)
  tp <- seq(0, 22, length.out = 12)
  Y <- t(sapply(1:14, function(i) 30 + rnorm(1, 0, 2) + 7 * cos(2 * pi * (tp - 5) / 24) + rnorm(12, 0, 1.2)))
  sp <- e$dance_traj_spec(e$dance_traj_long(Y, tp, sprintf("S%03d", 1:14), list()), 24, 1, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  ce <- e$dance_traj_cell_equations(ff)
  expect_equal(nrow(ce$table), 1L)
  expect_equal(ce$table$cell, "(all)")
  expect_lt(abs(ce$table$amplitude_1 - 7), 1.2)
  expect_lt(abs(ce$table$acrophase_time_1 - 5), 1.0)
  expect_null(ce$trend_coefs[[1]])
})
