# ==============================================================================
# tests/testthat/test-cosinor-assumptions.R
#
# Nothing a cosinor run decides may be silent. These pin the assumptions an
# audit found being made without a word -- and, where one was wrong, the fix:
#
#   * clock times parsed at import went in without unwrapping midnight, and the
#     hand-entered ones were unwrapped by one PERIOD instead of 24 h;
#   * a participant missing the first observation had their trend, and so
#     their intercept, anchored at their own first time;
#   * converged nonlinear fits were labelled "not_attempted" (an assignment
#     that went to the enclosing environment);
#   * the nested-model table averaged each candidate over whoever it fitted;
#   * and the list the summary and the report print states each of them.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
srv <- function(f) file.path(app_dir, "server", f)
for (f in c("01_helpers_time.R", "01b_kernel_contracts.R", "03_helpers_clock.R",
            "07_helpers_circular.R", "08_helpers_cosinor.R", "08c_helpers_circstat.R",
            "08j_helpers_assumptions.R"))
  if (file.exists(srv(f))) eval(parse(srv(f), encoding = "UTF-8"), envir = environment())
# the two fitters, out of the server file they live in
for (e in parse(srv("72_harmonic.R"), encoding = "UTF-8"))
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% c("fit_cosinor", "fit_cosinor_nonlinear"))
    eval(e, envir = environment())

test_that("clock times are put in recording order by 24 h, whatever the period", {
  clock <- c(8, 12, 16, 20, 0, 4, 8, 12, 15)
  u <- dance_unwrap_clock(clock)
  expect_true(u$wrapped)
  expect_equal(u$times, c(8, 12, 16, 20, 24, 28, 32, 36, 39))
  expect_true(all(diff(u$times) > 0))
  # already forward: nothing moves, and it says so
  v <- dance_unwrap_clock(c(8, 9, 10))
  expect_false(v$wrapped); expect_equal(v$times, c(8, 9, 10))
  # the Run handler unwraps every clock route through this, by 24 h -- not by
  # the period, and not only for hand-entered times
  run_src <- paste(readLines(srv("72_harmonic.R"), warn = FALSE), collapse = "\n")
  expect_equal(lengths(regmatches(run_src, gregexpr("dance_unwrap_clock\\(", run_src))), 3L)
  expect_false(grepl("time_vec[i:length(time_vec)] + period", run_src, fixed = TRUE))
})

test_that("every participant's trend starts at the recording's first time point", {
  t <- seq(0, 31, by = 1)
  # noise-free, so the intercept of a correct fit is exactly the true 10
  y_log <- 10 + 2 * log(t + 1) + 3 * cos(2 * pi * (t - 6) / 24)
  miss <- y_log; miss[1:2] <- NA
  # (an exact fit: lm warns that its summary is unreliable, which is the point)
  shared <- suppressWarnings(fit_cosinor(t, miss, 24, 1, "log", t_origin = 0))
  own <- suppressWarnings(fit_cosinor(t, miss, 24, 1, "log"))   # the old behaviour
  expect_equal(shared$t_offset, 0)
  expect_equal(unname(shared$mesor), 10, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(own$mesor), 10, tolerance = 1e-3)))

  y_sat <- 10 + 8 * (1 - exp(-t / 15)) + 3 * cos(2 * pi * (t - 6) / 24)
  ms <- y_sat; ms[1:2] <- NA
  f <- fit_cosinor(t, ms, 24, 1, "exp_sat", tau_min = 1, tau_max = 72, t_origin = 0)
  expect_true(f$success)
  expect_equal(f$t_offset, 0)
  expect_equal(unname(f$mesor), 10, tolerance = 1e-4)
  expect_equal(f$trend_params$tau$coef, 15, tolerance = 1e-3)
})

test_that("a converged nonlinear fit says it converged, and leaks nothing", {
  t <- seq(0, 31, by = 1)
  set.seed(3)
  y <- 10 + 8 * (1 - exp(-t / 15)) + 3 * cos(2 * pi * (t - 6) / 24) + rnorm(length(t), 0, 0.5)
  env <- environment()
  had <- exists("conv_status", envir = env, inherits = FALSE)
  f <- fit_cosinor(t, y, 24, 1, "exp_sat", tau_min = 1, tau_max = 72, t_origin = 0)
  expect_true(f$converged)
  expect_identical(f$convergence, "converged")
  expect_identical(exists("conv_status", envir = env, inherits = FALSE), had)
})

test_that("the nested set compares every candidate on the same participants", {
  m <- cbind(A = c(10, 12, NA, 14), B = c(11, 11, 9, 13), C = c(NA, NA, NA, NA))
  r <- dance_model_set_common(m)
  # participant 3 cannot be scored under A, so nobody's B is averaged with it
  expect_equal(r$n_common, 3L)
  expect_equal(r$means$A, mean(c(10, 12, 14)))
  expect_equal(r$means$B, mean(c(11, 11, 13)))
  # a candidate fitted to nobody is named, not allowed to empty the sample
  expect_identical(r$not_fitted, "C")
  tab <- dance_akaike_table(r$means)
  expect_setequal(tab$model, c("A", "B"))
})

# A minimal fitted mixed model: the function reads fields, it fits nothing.
mixed_mod <- function(trend = "exp_sat", tau_mode = "profile") {
  ladder <- list(list(label = "participant trend and rhythm"),
                 list(label = "participant rhythm"))
  ff <- list(ok = TRUE, REML = TRUE, re_rung = 2L, n_rungs = 2L, re_label = "participant rhythm",
             tau_mode = tau_mode, tau_range = c(2, 48),
             tau_profile = list(grid = exp(seq(log(2), log(48), length.out = 24)), rung = 2L),
             spec = list(re_ladder = ladder, covariates = character(0), tau = 18.3))
  mod <- list(approach = "mixed", period = 24, n_harmonics = 2L, trend_type = trend,
              time_vec = c(0, 2, 4, 25), origin_shift = 8, time_source = "shared",
              wrap_applied = TRUE, tau_spec = list(mode = "estimate", value = 18.3, lo = 2, hi = 48))
  list(mod = mod, ff = ff)
}

test_that("the mixed model's list states the residual model and the tau procedure", {
  x <- mixed_mod()
  at <- dance_cosinor_assumptions(x$mod, x$ff, "kr")
  txt <- paste(at$text, collapse = " ")
  expect_true(all(c("Time", "Data", "Trend", "Model", "Estimation", "Random effects",
                    "Residuals", "Tests", "Intervals", "tau") %in% at$area))
  expect_match(txt, "INDEPENDENT, with ONE variance", fixed = TRUE)
  expect_match(txt, "moved 24 h later", fixed = TRUE)
  expect_match(txt, "t = 0 is the first time point of the recording (08:00)", fixed = TRUE)
  expect_match(txt, "Kenward-Roger", fixed = TRUE)
  expect_match(txt, "EQUAL weights", fixed = TRUE)
  expect_match(txt, "seed 20240601", fixed = TRUE)
  expect_match(txt, "24 values spaced evenly in log tau from 2 to 48 h", fixed = TRUE)
  expect_match(txt, "held at the held value, 18.3 h", fixed = TRUE)
  # a flat profile in THIS run is said to have happened
  x$ff$tau_mode <- "fixed_after_flat"
  expect_match(paste(dance_cosinor_assumptions(x$mod, x$ff)$text, collapse = " "),
               "which is what happened in this run", fixed = TRUE)
  # a model that did not fit lists the shared assumptions only, not another approach's
  x$ff <- list(ok = FALSE)
  at2 <- dance_cosinor_assumptions(x$mod, x$ff)
  expect_false(any(c("Residuals", "Per-participant fit", "Tests") %in% at2$area))
  expect_true("Time" %in% at2$area)
})

test_that("the two-stage list states the fit, the inclusion rule, the tests and the seed", {
  mod <- list(approach = "two_stage", period = 24, n_harmonics = 1L, trend_type = "exp_sat",
              time_vec = c(0, 2, 4, 31), origin_shift = 8, time_source = "index",
              wrap_applied = FALSE, tau_spec = list(mode = "estimate", lo = 5, hi = 30),
              dance_settings = list(use_bounds = FALSE, tau_min = 5, tau_max = 30),
              fit_audit = list(include_boundary = TRUE, gate_relaxed = FALSE),
              boot_results = list(B = 500L, seed = 20240601L),
              model_selection = structure(data.frame(model = "x"), n_common = 40L,
                                          n_harmonics_max = 1L),
              conditioning = list(tau_fixed_delta_aic = 1.2, tau_fixed_value = 18,
                                  tau_fixed_n = 40L))
  at <- dance_cosinor_assumptions(mod)
  txt <- paste(at$text, collapse = " ")
  # the column-index time axis is called what it is
  expect_match(txt, "NOT read from the data", fixed = TRUE)
  expect_match(txt, "tau kept inside 5 to 30 h in every attempt", fixed = TRUE)
  expect_match(txt, "INCLUDED (Advanced)", fixed = TRUE)
  expect_match(txt, "Wilks' lambda", fixed = TRUE)
  expect_match(txt, "Welch's t", fixed = TRUE)
  expect_match(txt, "treated as exact", fixed = TRUE)
  expect_match(txt, "B = 500; percentile 95% intervals", fixed = TRUE)
  expect_match(txt, "seed 20240601", fixed = TRUE)
  expect_match(txt, "averaged over the 40 participant(s) scorable under every candidate", fixed = TRUE)
  expect_match(txt, "averaged over the 40 participant(s) scorable in both", fixed = TRUE)
  expect_false("Residuals" %in% at$area)
  # held tau: the fit is described as linear given tau
  mod$tau_spec <- list(mode = "fixed", value = 18)
  expect_match(paste(dance_cosinor_assumptions(mod)$text, collapse = " "),
               "tau held at 18 h for every participant", fixed = TRUE)
})

test_that("the printed list fits the summary panel's width", {
  x <- mixed_mod()
  ln <- dance_cosinor_assumptions_text(dance_cosinor_assumptions(x$mod, x$ff))
  expect_true(length(ln) > 10)
  expect_true(all(nchar(ln) <= 92))
  # the area is printed once, where its block starts
  expect_equal(sum(startsWith(ln, "Time ")), 1L)
})
