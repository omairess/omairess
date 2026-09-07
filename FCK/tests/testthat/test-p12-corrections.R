# ==============================================================================
# tests/testthat/test-p12-corrections.R — reported directly by the user
#
#   P12.1  The publication report said the smoothing parameter was "chosen by
#          generalised cross-validation over the whole sample". It was not:
#          fck_auto_lambda() capped its objective at 60 subjects, taken as a
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
#            (c) fck_apa_M() had no `bounded` argument, so R-squared printed
#                with a leading zero, an APA 7 (6.36) error of the same kind
#                the suite already tested for in p-values and ranges.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
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
  a <- e$fck_auto_lambda(Y, tp, b)
  expect_equal(a$n_used, m)                 # was 60 before P12.1

  # a subject with too few observed points still cannot be scored, so n_used and
  # the sample size legitimately differ -- which is why the report prints both
  Y2 <- Y; Y2[1, ] <- NA; Y2[2, 3:n_time] <- NA
  a2 <- e$fck_auto_lambda(Y2, tp, b)
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
  md <- e$fck_apa_report(vals, list(is_cyclic = FALSE), "T")
  txt <- paste(md, collapse = "\n")
  expect_true(grepl("across all 118 of the 126 subjects", txt, fixed = TRUE))
  expect_false(grepl("over the whole sample", txt, fixed = TRUE))

  # when every subject was scored, the sentence does not invent a distinction
  vals$smooth_fit_metrics$lambda_n_used <- 126L
  txt2 <- paste(e$fck_apa_report(vals, list(is_cyclic = FALSE), "T"), collapse = "\n")
  expect_true(grepl("across all 126 subjects", txt2, fixed = TRUE))
  expect_false(grepl("of the 126", txt2, fixed = TRUE))
})

# ==================================== P12.2 publication-grade cosinor =========
test_that("P12.2a: the linear fitter computes the Bingham region too", {
  src <- code_of("server/72_harmonic.R")
  # both fitters must build it; before P12.2 only the nonlinear one did
  expect_equal(length(gregexpr("bingham[[h]] <- fck_bingham_ci", src, fixed = TRUE)[[1]]), 2L)
  expect_equal(length(gregexpr("bingham = bingham", src, fixed = TRUE)[[1]]), 2L)
})

test_that("P12.2c: fck_apa_M honours the no-leading-zero rule when bounded", {
  e <- rep_env()
  expect_equal(e$fck_apa_M(c(0.958, 0.940), 3, bounded = TRUE),
               "*M* = .949, *SD* = .013")
  # and still prints the zero for a quantity that is not bounded by 1
  expect_true(grepl("^\\*M\\* = 7\\.99", e$fck_apa_M(c(7.99, 7.99), 2)))
  # the report uses it for R-squared
  expect_true(grepl('fck_apa_M(g1("r_squared"), 3, bounded = TRUE)',
                    code_of("server/93_apa_report.R"), fixed = TRUE))
})

test_that("P12.2: internal parameter ids never reach report prose", {
  e <- rep_env()
  expect_equal(e$fck_cosinor_param_label("amplitude_1"), "the amplitude of harmonic 1")
  expect_equal(e$fck_cosinor_param_label("acrophase_time_2"), "the acrophase of harmonic 2")
  expect_equal(e$fck_cosinor_param_label("mesor_adj"), "the MESOR (rhythm-adjusted mean)")
  expect_equal(e$fck_cosinor_origin_label("first_observation"),
               "the first observation in each series")
})

test_that("P12.2b: the circular effect is a difference of MEANS, not concentration", {
  e <- rep_env()
  pw <- data.frame(comparison = "A vs B", group1 = "A", group2 = "B",
                   mean1 = 16.0, sd1 = 1, n1 = 30, mean2 = 19.0, sd2 = 1, n2 = 20,
                   mean_diff = -3.0, t_stat = 750.9, df = 48,
                   p_value = 1e-9, cohens_d = -0.0001,
                   ci_lower = NA, ci_upper = NA, p_adjusted = 1e-9,
                   stringsAsFactors = FALSE)
  vals <- list(
    data = matrix(rnorm(20 * 24), 20, 24),
    harmonic_model = list(period = 24, n_harmonics = 1, trend_type = "none",
                          individual_fits = list(), group_var_name = "AGE",
                          time_origin = "first_observation"),
    hp_pairwise_all = list(list(results = pw, param = "acrophase_time_1",
                                correction = "holm")))
  txt <- paste(e$fck_apa_report(vals, list(), "T"), collapse = "\n")
  # the angular difference, in hours -- what the test is about
  expect_true(grepl("-3.00 h", txt, fixed = TRUE))
  # NOT the near-zero concentration difference that used to fill this column
  expect_false(grepl("| -0.00 |", txt, fixed = TRUE))
  expect_true(grepl("difference between the group circular means", txt, fixed = TRUE))
  # and the Watson-Williams assumptions are stated
  expect_true(grepl("von Mises", txt, fixed = TRUE))
})

test_that("P12.2: Bingham's amplitude/acrophase caution is CHECKED, not recited", {
  e <- rep_env()
  mk <- function(param, acro_verdict) {
    pw <- data.frame(comparison = "A vs B", group1 = "A", group2 = "B",
                     mean1 = 8, sd1 = 1, n1 = 30, mean2 = 9, sd2 = 1, n2 = 20,
                     mean_diff = -1, t_stat = -2.3, df = 40,
                     p_value = .02, cohens_d = -.6,
                     ci_lower = NA, ci_upper = NA, p_adjusted = .02,
                     stringsAsFactors = FALSE)
    list(data = matrix(rnorm(20 * 24), 20, 24),
         harmonic_model = list(period = 24, n_harmonics = 1, trend_type = "none",
                               individual_fits = list(), group_var_name = "AGE"),
         hp_pairwise_all = list(list(results = pw, param = param, correction = "holm")),
         hp_acrophase_differs = acro_verdict)
  }
  hit  <- paste(e$fck_apa_report(mk("amplitude_1", TRUE),  list(), "T"), collapse = "\n")
  miss <- paste(e$fck_apa_report(mk("amplitude_1", FALSE), list(), "T"), collapse = "\n")
  unk  <- paste(e$fck_apa_report(mk("amplitude_1", NULL),  list(), "T"), collapse = "\n")
  expect_true(grepl("acrophase comparison WAS significant", hit, fixed = TRUE))
  expect_true(grepl("was not significant", miss, fixed = TRUE))
  expect_true(grepl("Run the same comparison on acrophase", unk, fixed = TRUE))
  # and it is not raised when no amplitude comparison was run
  none <- paste(e$fck_apa_report(mk("mesor_adj", TRUE), list(), "T"), collapse = "\n")
  expect_false(grepl("should not be read as a difference in rhythm strength",
                     none, fixed = TRUE))
})

test_that("P12.2: every comparison the user ran reaches the report", {
  e <- rep_env()
  mk1 <- function(param) {
    pw <- data.frame(comparison = "A vs B", group1 = "A", group2 = "B",
                     mean1 = 8, sd1 = 1, n1 = 30, mean2 = 9, sd2 = 1, n2 = 20,
                     mean_diff = -1, t_stat = -2.3, df = 40, p_value = .02,
                     cohens_d = -.6, ci_lower = NA, ci_upper = NA,
                     p_adjusted = .02, stringsAsFactors = FALSE)
    list(results = pw, param = param, correction = "holm")
  }
  vals <- list(
    data = matrix(rnorm(20 * 24), 20, 24),
    harmonic_model = list(period = 24, n_harmonics = 1, trend_type = "none",
                          individual_fits = list(), group_var_name = "AGE"),
    hp_pairwise_all = list(mk1("mesor_adj"), mk1("amplitude_1"),
                           mk1("acrophase_time_1")))
  txt <- paste(e$fck_apa_report(vals, list(), "T"), collapse = "\n")
  for (lab in c("the MESOR (rhythm-adjusted mean)", "the amplitude of harmonic 1",
                "the acrophase of harmonic 1"))
    expect_true(grepl(paste0("Group differences in ", lab), txt, fixed = TRUE), info = lab)
  # and the multiplicity statement says what family was controlled
  expect_true(grepl("within each parameter, not across parameters", txt, fixed = TRUE))
  # the pairwise module keeps them rather than overwriting
  expect_true(grepl("all_pw[[param]] <- list(results = results",
                    code_of("server/73_cosinor_pairwise.R"), fixed = TRUE))
})

test_that("P12.2: the methods paragraph carries what a cosinor paper needs", {
  e <- rep_env()
  vals <- list(data = matrix(rnorm(20 * 24), 20, 24),
               harmonic_model = list(
    period = 24, n_harmonics = 1, trend_type = "none", individual_fits = list(),
    dv_name = "Activity", dv_units = "counts/min",
    time_origin = "first_observation", group_var_name = "AGE"))
  txt <- paste(e$fck_apa_report(vals, list(), "T"), collapse = "\n")
  expect_true(grepl("dependent variable was Activity (counts/min)", txt, fixed = TRUE))
  expect_true(grepl("period was FIXED at 24 rather than estimated", txt, fixed = TRUE))
  expect_true(grepl("measured from the first observation in each series", txt, fixed = TRUE))
  expect_true(grepl("zero-amplitude *F* test", txt, fixed = TRUE))
  expect_true(grepl("error ellipse", txt, fixed = TRUE))
  expect_true(grepl("Bingham et al. (1982)", txt, fixed = TRUE))
})
