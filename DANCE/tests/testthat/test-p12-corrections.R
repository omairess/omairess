# ==============================================================================
# tests/testthat/test-p12-corrections.R — reported directly by the user
#
#   P12.1  The publication report said the smoothing parameter was "chosen by
#          generalised cross-validation over the whole sample". It was not:
#          dance_auto_lambda() capped its objective at 60 subjects, taken as a
#          systematic subsample, and returned n_used, which the report ignored.
#          The number was fine -- capped and uncapped searches agree to well
#          inside the range over which GCV is flat -- but the SENTENCE was
#          false, and that sentence is the one that goes into a paper. The cap
#          is gone and the report states the count instead of asserting a scope.
#
#   P12.2  The cosinor section was not publication-grade. Against Cornelissen
#          (2014) and Bingham et al. (1982) it was missing: the reference time
#          the acrophase is measured from, confidence regions for amplitude and
#          acrophase, and -- with a grouping variable selected and pairwise
#          tests already computed by the app -- ANY between-group result at all.
#          Three defects were found while fixing it:
#            (a) the linear (OLS) fitter never computed the Bingham joint
#                confidence region, so the default analysis had no interval on
#                either amplitude or acrophase;
#            (b) the report's effect column for the circular comparison printed
#                the difference in MEAN RESULTANT LENGTH -- a difference in
#                concentration -- next to a Watson-Williams test of MEANS;
#            (c) dance_apa_M() had no `bounded` argument, so R-squared printed
#                with a leading zero, an APA 7 (6.36) error of the same kind
#                the suite already tested for in p-values and ranges.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  ln <- sub("#.*$", "", ln)
  paste(ln, collapse = "\n")
}
rep_env <- function() {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/93_apa_report.R"), local = e)
  e
}
# The report's two-stage comparison calls the same kernels the comparison tab
# does (server/08h), so a test of that section needs them in scope.
ts_env <- function() {
  e <- new.env(parent = globalenv())
  for (f in c("server/08_helpers_cosinor.R", "server/08b_helpers_popcosinor.R",
              "server/08c_helpers_circstat.R", "server/08h_helpers_twostage.R",
              "server/93_apa_report.R"))
    source(file.path(app_dir, f), local = e)
  e
}
# A two-stage cosinor result with two groups: 30 vs 20 participants, H1 with a
# planted acrophase (hours) and amplitude per group, as the run button stores it.
mk_ts_vals <- function(acro = c(A = 16, B = 19), amp = c(A = 8, B = 9), seed = 12) {
  set.seed(seed)
  n <- c(A = 30, B = 20); g <- rep(names(n), n)
  phi <- 2 * pi * (acro[g] + rnorm(sum(n), 0, 0.4)) / 24
  A <- amp[g] + rnorm(sum(n), 0, 0.8)
  ip <- data.frame(subject = seq_len(sum(n)), mesor = 50 + rnorm(sum(n), 0, 2),
                   beta_cos_1 = A * cos(phi), beta_sin_1 = A * sin(phi),
                   amplitude_1 = A, acrophase_rad_1 = phi %% (2 * pi),
                   acrophase_time_1 = (phi %% (2 * pi)) * 24 / (2 * pi),
                   r_squared = runif(sum(n), .7, .95), p_value = 1e-4,
                   stringsAsFactors = FALSE)
  list(data = matrix(rnorm(50 * 24), 50, 24),
       covariates = data.frame(AGE = factor(g), stringsAsFactors = FALSE),
       harmonic_model = list(approach = "two_stage", period = 24, n_harmonics = 1,
                             trend_type = "none", individual_fits = list(),
                             individual_params = ip, group_var_name = "AGE",
                             time_vec = 0:23, origin_shift = 0))
}

# ============================================ P12.1 the lambda claim ==========
test_that("P12.1: the GCV search uses every scorable subject", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = e)

  set.seed(7); n_time <- 24; tp <- seq_len(n_time)
  b <- create.bspline.basis(c(1, n_time), nbasis = 10)
  m <- 80                                   # deliberately more than the old cap
  Y <- t(sapply(1:m, function(i) 5*sin(2*pi*tp/24 + runif(1,0,6)) + rnorm(n_time,0,1.5)))
  a <- e$dance_auto_lambda(Y, tp, b)
  expect_equal(a$n_used, m)                 # was 60 before P12.1

  # a subject with too few observed points still cannot be scored, so n_used and
  # the sample size legitimately differ -- which is why the report prints both
  Y2 <- Y; Y2[1, ] <- NA; Y2[2, 3:n_time] <- NA
  a2 <- e$dance_auto_lambda(Y2, tp, b)
  expect_lt(a2$n_used, m)
  expect_gt(a2$n_used, 60)
})

test_that("P12.1: the cap is gone from the source, not just from the sentence", {
  src <- code_of("server/04_helpers_fd.R")
  expect_false(grepl("rows[round(seq(1, length(rows), length.out = 60))]",
                     src, fixed = TRUE))
  expect_false(grepl("if (length(rows) > 60)", src, fixed = TRUE))
})

test_that("P12.1: the report states the count instead of claiming a scope", {
  e <- rep_env()
  src <- code_of("server/93_apa_report.R")
  expect_false(grepl("cross-validation over the whole sample", src, fixed = TRUE))
  expect_true(grepl("lambda_n_used", src, fixed = TRUE))

  vals <- list(data = matrix(rnorm(126 * 24), 126, 24),
               smooth_fit_metrics = list(
    method = "auto", n_basis = 12, lambda = 4.59e-01,
    lambda_n_used = 118L, n_subjects = 126L, time_axis = "hours"))
  md <- e$dance_apa_report(vals, list(is_cyclic = FALSE), "T")
  txt <- paste(md, collapse = "\n")
  expect_true(grepl("across all 118 of the 126 subjects", txt, fixed = TRUE))
  expect_false(grepl("over the whole sample", txt, fixed = TRUE))

  # when every subject was scored, the sentence does not invent a distinction
  vals$smooth_fit_metrics$lambda_n_used <- 126L
  txt2 <- paste(e$dance_apa_report(vals, list(is_cyclic = FALSE), "T"), collapse = "\n")
  expect_true(grepl("across all 126 subjects", txt2, fixed = TRUE))
  expect_false(grepl("of the 126", txt2, fixed = TRUE))
})

# ==================================== P12.2 publication-grade cosinor =========
test_that("P12.2a: the linear fitter computes the Bingham region too", {
  src <- code_of("server/72_harmonic.R")
  # both fitters must build it; before P12.2 only the nonlinear one did
  expect_equal(length(gregexpr("bingham[[h]] <- dance_bingham_ci", src, fixed = TRUE)[[1]]), 2L)
  expect_equal(length(gregexpr("bingham = bingham", src, fixed = TRUE)[[1]]), 2L)
})

test_that("P12.2c: dance_apa_M honours the no-leading-zero rule when bounded", {
  e <- rep_env()
  expect_equal(e$dance_apa_M(c(0.958, 0.940), 3, bounded = TRUE),
               "*M* = .949, *SD* = .013")
  # and still prints the zero for a quantity that is not bounded by 1
  expect_true(grepl("^\\*M\\* = 7\\.99", e$dance_apa_M(c(7.99, 7.99), 2)))
  # the report uses it for R-squared
  expect_true(grepl('dance_apa_M(g1("r_squared"), 3, bounded = TRUE)',
                    code_of("server/93_apa_report.R"), fixed = TRUE))
})

test_that("P12.2: the legacy pairwise module has left the report", {
  src <- code_of("server/93_apa_report.R")
  for (bad in c("hp_pairwise", "hp_acrophase", "dance_cosinor_param_label",
                "dance_cosinor_origin_label", "using_smoothed"))
    expect_false(grepl(bad, src, fixed = TRUE), info = bad)
})

test_that("P12.2: the two-stage group comparison reaches the report, as one table", {
  e <- ts_env()
  txt <- paste(e$dance_apa_report(mk_ts_vals(), list(), "T"), collapse = "\n")
  expect_true(grepl("Cosinor regression (two-stage)", txt, fixed = TRUE))
  expect_true(grepl("**Group differences (AGE).**", txt, fixed = TRUE))
  # the primary (MANOVA on the coefficient vector) and every component test
  for (row in c("Complete coefficient vector", "Constant term (beta_0)",
                "H1 rhythm vector (cos, sin) jointly", "H1 amplitude", "H1 acrophase",
                "Watson-Williams"))
    expect_true(grepl(row, txt, fixed = TRUE), info = row)
  # the two-stage caveats a chronobiology reviewer will raise
  expect_true(grepl("treated as exact", txt, fixed = TRUE))
  expect_true(grepl("von Mises", txt, fixed = TRUE))
  # and the Methods paragraph names the tests
  expect_true(grepl("one-way MANOVA (Wilks' lambda)", txt, fixed = TRUE))
})

test_that("P12.2: Bingham's amplitude/acrophase caution is CHECKED, not recited", {
  e <- ts_env()
  # acrophases three hours apart: the amplitude comparison is not interpretable
  hit  <- paste(e$dance_apa_report(mk_ts_vals(acro = c(A = 16, B = 19)), list(), "T"), collapse = "\n")
  # same acrophase: the acrophase test had power and did not reject, so it clears
  miss <- paste(e$dance_apa_report(mk_ts_vals(acro = c(A = 16, B = 16)), list(), "T"), collapse = "\n")
  expect_true(grepl("should NOT be read as a difference in rhythm strength", hit, fixed = TRUE))
  expect_false(grepl("should NOT be read as a difference in rhythm strength", miss, fixed = TRUE))
  expect_true(grepl("the caution does not apply", miss, fixed = TRUE))
})

test_that("P12.2: the methods paragraph carries what a cosinor paper needs", {
  e <- rep_env()
  vals <- list(data = matrix(rnorm(20 * 24), 20, 24),
               harmonic_model = list(
    period = 24, n_harmonics = 1, trend_type = "none", individual_fits = list(),
    dv_name = "Activity", dv_units = "counts/min",
    time_origin = "first_observation", group_var_name = "AGE"))
  txt <- paste(e$dance_apa_report(vals, list(), "T"), collapse = "\n")
  expect_true(grepl("dependent variable was Activity (counts/min)", txt, fixed = TRUE))
  expect_true(grepl("period was FIXED at 24 rather than estimated", txt, fixed = TRUE))
  expect_true(grepl("measured from the first observation in each series", txt, fixed = TRUE))
  expect_true(grepl("zero-amplitude *F* test", txt, fixed = TRUE))
  expect_true(grepl("error ellipse", txt, fixed = TRUE))
  expect_true(grepl("Bingham et al. (1982)", txt, fixed = TRUE))
})
