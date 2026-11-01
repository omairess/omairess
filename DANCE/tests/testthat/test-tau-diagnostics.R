# ==============================================================================
# tests/testthat/test-tau-diagnostics.R
#
# Convergence is not identifiability. The tau checklist must call tau identified
# when the recording covers several time constants, and NOT identified when tau
# is far beyond the window (A and tau on a ridge, estimates on the bound), and
# then recommend a fixed-tau sensitivity analysis. The hierarchical-tau option
# must recover a between-participant spread of tau.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("07_helpers_circular", "08_helpers_cosinor", "08b_helpers_popcosinor", "08c_helpers_circstat",
            "08d_helpers_traj", "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband",
            "08h_helpers_twostage", "08i_helpers_trajsat", "08k_helpers_outcomes"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)
ex <- parse(file.path(app_dir, "server", "72_harmonic.R"), encoding = "UTF-8")
for (x in ex) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
                  as.character(x[[2]]) %in% c("fit_cosinor", "fit_cosinor_nonlinear")) eval(x, e)

ts_diag <- function(tp, tau, seed) {
  set.seed(seed); n <- 20
  Y <- t(sapply(1:n, function(i) 20 + rnorm(1) + 5 * (1 - exp(-tp / tau)) +
                  2 * cos(2 * pi * (tp - 8) / 24) + rnorm(length(tp), 0, 0.4)))
  fa <- e$dance_ts_fit_all(Y, tp, 24, 1, "exp_sat", tau_min = 0.5, tau_max = 200)
  mod <- list(individual_params = fa$params, individual_fits = fa$fits, period = 24, n_harmonics = 1,
              trend_type = "exp_sat", time_vec = tp, t_offset = 0, tau_spec = list(mode = "estimate"),
              fit_audit = list(n_attempted = n))
  e$dance_tau_diagnostics_ts(mod)
}

test_that("two-stage: a window of several time constants identifies tau", {
  td <- ts_diag(seq(0, 48, by = 2), 8, 1)
  expect_true(isTRUE(td$ok))
  expect_identical(td$status, "identified")
  expect_false(td$recommend_fixed)
  expect_setequal(td$rows$check, c("convergence", "boundary", "precision", "A-tau dependence", "window",
                                   "tau held (Delta-AIC)", "valid n"))
})

test_that("two-stage: tau far beyond the window is flagged, and a held tau recommended", {
  td <- ts_diag(seq(0, 12, by = 1), 60, 2)
  expect_identical(td$status, "not identified")
  expect_true(td$recommend_fixed)
  expect_match(td$message, "sensitivity analysis", fixed = TRUE)
  expect_true(any(grepl("ridge", c(td$problems, td$weak))))
})

test_that("a held tau has no identifiability question", {
  mod <- list(trend_type = "exp_sat", tau_spec = list(mode = "fixed", value = 10))
  expect_false(isTRUE(e$dance_tau_diagnostics_ts(mod)$ok))
})

test_that("mixed: hierarchical tau recovers a between-participant spread", {
  skip_if_not_installed("nlme"); skip_if_not_installed("lme4")
  set.seed(5); tp <- seq(0, 40, by = 2); n <- 40; g <- rep(c("A", "B"), each = n / 2)
  taui <- exp(log(10) + rnorm(n, 0, 0.35))
  Y <- t(sapply(1:n, function(i) 20 + rnorm(1) + 5 * (1 - exp(-tp / taui[i])) +
                  2 * cos(2 * pi * (tp - 8) / 24) + rnorm(length(tp), 0, 0.3)))
  d <- e$dance_traj_long(Y, tp, sprintf("S%02d", 1:n), list(Group = g))
  sp <- e$dance_traj_spec(d, 24, 1, "exp_sat", tau = 10)
  ff <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp)))
  ff$sat <- suppressWarnings(e$dance_traj_sat_fit(ff, tau_range = c(1, 60), hierarchical_tau = TRUE))
  skip_if_not(isTRUE(ff$sat$ok))
  expect_true(ff$sat$hierarchical_tau)
  # population tau near 10 h in each cell; a between-participant SD of a few hours
  expect_true(all(abs(ff$sat$table$tau - 10) < 3))
  expect_true(all(ff$sat$table$tau_between_sd_h > 1 & ff$sat$table$tau_between_sd_h < 8))
  td <- e$dance_tau_diagnostics_mixed(ff)
  expect_true("hierarchical tau" %in% td$rows$check)
  # the per-cell change and average rate carry tau's uncertainty: wider than zero,
  # and the average rate is the change over the window
  st <- ff$sat$table
  expect_equal(st$avg_rate, st$change / 40, tolerance = 1e-10)
  expect_true(all(st$change_hi > st$change_lo))
})
