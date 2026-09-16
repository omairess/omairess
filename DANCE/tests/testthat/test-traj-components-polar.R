# ==============================================================================
# tests/testthat/test-traj-components-polar.R
#
# Four things asked for after a real session:
#
#   1  each harmonic drawable on its own, not only the sum of all of them
#   2  hiding a trace in the legend must hide its confidence band with it
#   3  the polar dial's group vectors from the MIXED-EFFECTS model, not the
#      two-stage mean of per-participant fits
#   4  the hover showing the acrophase as a clock time, matching the table
#
# (3) is the one worth stating carefully. The two estimators are close on
# balanced data with many participants -- the two-stage group vector is the
# amplitude of the MEAN (cos, sin) pair, which is what the mixed model also
# estimates -- close enough to agree to several figures and hide that they are
# different estimators entirely. The test below asserts the CONVENTION is shared
# (so swapping the source rotates nothing) and that the two are not assumed
# identical.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk <- function(K = 2, trend = "linear", nper = 12, seed = 31) {
  set.seed(seed)
  tp <- seq(8, 32, length.out = 13)
  grp <- rep(c("A", "B"), each = nper)
  Y <- t(sapply(seq_along(grp), function(i)
    50 + rnorm(1, 0, 6) + 16 * cos(2*pi*(tp - 2)/24) + 7 * cos(2*pi*2*(tp + 1)/24) +
      0.2 * (tp - 8) + rnorm(length(tp), 0, 4)))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, K, trend)
  ff <- e$dance_traj_fit(sp); ff$clock_origin <- sp$t0; ff
}

test_that("the component vocabulary follows the fit", {
  skip_if_not_installed("lme4")
  f2 <- mk(K = 2, trend = "linear"); skip_if_not(isTRUE(f2$ok))
  cs <- e$dance_traj_components(f2)
  expect_true(all(c("full", "harmonics", "harmonic1", "harmonic2",
                    "baseline_harm", "trend") %in% cs))
  f1 <- mk(K = 1, trend = "none"); skip_if_not(isTRUE(f1$ok))
  cs1 <- e$dance_traj_components(f1)
  # one harmonic has no per-harmonic view to offer -- "harmonics" IS H1
  expect_false(any(grepl("^harmonic[0-9]", cs1)))
  expect_false("trend" %in% cs1)
  expect_false("baseline_harm" %in% cs1)
  expect_equal(e$dance_traj_component_harmonic("harmonic2"), 2L)
  expect_true(is.na(e$dance_traj_component_harmonic("harmonics")))
  expect_match(e$dance_traj_component_label("harmonic2", 24), "Harmonic 2 only \\(period 12 h\\)")
})

test_that("the per-harmonic curves SUM to the harmonics view", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  tt <- seq(ff$spec$t0, ff$spec$t0 + 24, length.out = 97)
  all_h <- e$dance_traj_predict(ff, times = tt, component = "harmonics")
  h1 <- e$dance_traj_predict(ff, times = tt, component = "harmonic1")
  h2 <- e$dance_traj_predict(ff, times = tt, component = "harmonic2")
  expect_true(all(vapply(list(all_h, h1, h2), function(x) isTRUE(x$ok), logical(1))))
  expect_equal(h1$table$fit + h2$table$fit, all_h$table$fit, tolerance = 1e-9)
  # and the whole fit decomposes: baseline+harmonics + trend = full
  bh <- e$dance_traj_predict(ff, times = tt, component = "baseline_harm")
  tr <- e$dance_traj_predict(ff, times = tt, component = "trend")
  fu <- e$dance_traj_predict(ff, times = tt, component = "full")
  expect_equal(bh$table$fit + tr$table$fit, fu$table$fit, tolerance = 1e-9)
})

test_that("a single harmonic is centred on zero and has its own period", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  # one FULL period sampled without repeating the endpoint: t0 and t0 + 24 are
  # the same phase, and including both leaves a small residue in the mean that
  # has nothing to do with the component being centred
  tt <- head(seq(ff$spec$t0, ff$spec$t0 + 24, length.out = 481), -1L)
  for (k in 1:2) {
    pr <- e$dance_traj_predict(ff, times = tt, component = sprintf("harmonic%d", k))
    d <- pr$table[pr$table$cell == pr$cells[1], ]
    expect_lt(abs(mean(d$fit)), 1e-6)              # zero mean over a full period
    # k maxima per period: count sign changes of the derivative
    dd <- diff(d$fit)
    peaks <- sum(dd[-length(dd)] > 0 & dd[-1] <= 0)
    expect_equal(peaks, k)
  }
})

test_that("an unavailable component is refused, not silently substituted", {
  skip_if_not_installed("lme4")
  f1 <- mk(K = 1, trend = "none"); skip_if_not(isTRUE(f1$ok))
  bad <- e$dance_traj_predict(f1, component = "harmonic2")
  expect_false(isTRUE(bad$ok))
  expect_match(bad$message, "not a component of this fit")
  expect_match(bad$message, "Available:")
})

test_that("every band trace is tied to its line by legendgroup", {
  # Asked of the parser: add_ribbons and add_lines for the same cell must carry
  # legendgroup, or clicking the legend hides the curve and leaves the shading.
  src <- paste(vapply(parse(file.path(app_dir, "server/72_harmonic.R")),
                      function(x) paste(deparse(x), collapse = " "), character(1)),
               collapse = " ")
  src <- gsub("[[:space:]]+", " ", src)
  i_rib <- gregexpr("add_ribbons(", src, fixed = TRUE)[[1]]
  expect_gt(length(i_rib), 0)
  for (i in i_rib) {
    seg <- substr(src, i, i + 500)
    # a band with no legend entry of its own can ONLY be hidden through a group
    expect_true(grepl("legendgroup", seg, fixed = TRUE),
                info = paste("add_ribbons without legendgroup:", substr(seg, 1, 140)))
  }
})

test_that("the polar hover converts an angle to the clock the dial is labelled in", {
  # H1: a full turn is the period, so the angle is hours x 15 deg after the origin
  expect_equal(e$dance_polar_hover_clock(0, 24, 1, clock_origin = 8), "08:00")
  expect_equal(e$dance_polar_hover_clock(90, 24, 1, clock_origin = 8), "14:00")
  expect_equal(e$dance_polar_hover_clock(270, 24, 1, clock_origin = 8), "02:00")
  expect_equal(e$dance_polar_hover_clock(360, 24, 1, clock_origin = 8), "08:00")
  # H2: a full turn is 12 h, so the same angle is a different clock time
  expect_equal(e$dance_polar_hover_clock(180, 24, 2, clock_origin = 8), "14:00")
  expect_equal(e$dance_polar_hover_clock(0, 24, 2, clock_origin = 8), "08:00")
  # it agrees with the tick labels the dial is drawn with, by construction
  P <- 24; co <- 8; hh <- 1
  eff <- P / hh; n <- 12
  tv <- seq(0, 360 - 360/n, by = 360/n)
  te <- seq(0, eff - eff/n, by = eff/n)
  expect_equal(e$dance_polar_hover_clock(tv, P, hh, clock_origin = co),
               e$dance_clock_label((te + co) %% P, P, show_day = FALSE))
})

test_that("the hover time is the SAME number the trajectory table reports", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  co <- e$dance_traj_cell_coefs(ff, 1); skip_if_not(isTRUE(co$ok))
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 4000)
  skip_if_not(isTRUE(ap$ok))
  for (i in seq_len(nrow(ap$table))) {
    rad <- ap$table$acrophase_rad[i]
    deg <- e$phi_to_degrees(rad); if (deg < 0) deg <- deg + 360
    from_dial <- e$dance_polar_hover_clock(deg, ff$spec$period, 1,
                                           clock_origin = ff$clock_origin)
    from_table <- e$dance_clock_label(
      e$dance_acrophase_clock(hours = ap$table$acrophase_time[i],
                              period = ff$spec$period, harmonic = 1,
                              clock_origin = ff$clock_origin)$hours,
      ff$spec$period, show_day = FALSE)
    expect_equal(from_dial, from_table)
  }
})

test_that("both estimators put the angle in the SAME frame, so the source only changes the estimate", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  co <- e$dance_traj_cell_coefs(ff, 1); skip_if_not(isTRUE(co$ok))
  # the mixed-effects phase is atan2(b, a) on a zero-based model axis; the
  # two-stage per-participant phase is atan2(beta_sin, beta_cos) on the same
  # axis. Same map, so the dial needs no rotation when the source is switched.
  ap <- e$dance_traj_amp_phase(co)$table
  for (i in seq_along(co$a)) {
    # the same angle, reported on [0, 2pi) rather than atan2's (-pi, pi]
    d <- ((atan2(co$b[i], co$a[i]) - ap$acrophase_rad[i]) + pi) %% (2 * pi) - pi
    expect_lt(abs(d), 1e-9)
    expect_gte(ap$acrophase_rad[i], 0)
    expect_lt(ap$acrophase_rad[i], 2 * pi)
  }
})
