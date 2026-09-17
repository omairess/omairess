# ==============================================================================
# tests/testthat/test-p13-corrections.R — the rest of the publication-grade pass
#
# The user asked for the report to be brought to publication grade and named the
# cosinor section as an EXAMPLE. Reading a full generated report end to end found
# four more things, the first of which was severe and was mine, from P12:
#
#   P13.1  `pw` held the fANOVA post-hoc object, assigned near the top of the
#          function and read much further down in Results. The cosinor Methods
#          block, added at P12, REBOUND that name. So the post-hoc Results
#          section printed
#              "0 of 0 comparisons remained significant after the -- correction
#               at alpha = --"
#          on a run whose own Methods paragraph, four lines above, correctly
#          described 6 comparisons at B = 777. Only a report containing BOTH
#          kinds of pairwise result shows it, which is why it shipped.
#
#   P13.2  the function-on-scalar caveat was printed unconditionally and
#          described the POINTWISE OLS fit -- "each time point was fitted
#          separately", "the pointwise intervals are not simultaneous bands" --
#          under a GAM fit, where the curves are penalised rather than
#          independent and the sentence above it has just said no pointwise
#          inference exists.
#
#   P13.3  missing data was reported only when a smoothing step had run, because
#          the statement was gated on values$fill_status. Cosinor needs no
#          smoothing, so import -> cosinor -> report produced a document that
#          never mentioned missingness. APA 7 asks for it either way.
#
#   P13.4  the report never said WHICH CURVES an analysis ran on. With a
#          registration in the session that is the difference between a claim
#          about amplitude alone and a claim about amplitude and phase together.
#          The stored fANOVA data_source recorded the REQUEST; the app falls back
#          to the original curves when a warped representation is unavailable, so
#          the report could have stated the opposite of what happened.
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
base_vals <- function(...) {
  c(list(data = matrix(rnorm(30 * 24), 30, 24)), list(...))
}

# fANOVA post-hoc object, shaped as perform_pairwise_comparisons returns it
mk_posthoc <- function() {
  nm <- c("A vs B", "A vs C", "B vs C")
  list(pair_names = nm, correction_method = "bonferroni", alpha = 0.05,
       n_permutations = 777,
       results = setNames(lapply(seq_along(nm), function(i)
         list(n1 = 30, n2 = 20, L2_stat = 2 + i, p_value_L2 = 0.01 * i,
              p_value_L2_adjusted = 0.03 * i, sig_global = i == 1)), nm))
}
# a two-stage cosinor result with a grouping variable, so the report's cosinor
# comparison table renders next to the fANOVA post-hoc one
ts_env <- function() {
  e <- new.env(parent = globalenv())
  for (f in c("server/08_helpers_cosinor.R", "server/08b_helpers_popcosinor.R",
              "server/08c_helpers_circstat.R", "server/08h_helpers_twostage.R",
              "server/93_apa_report.R"))
    source(file.path(app_dir, f), local = e)
  e
}
mk_cos_model <- function() {
  set.seed(13); n <- c(A = 30, B = 20); g <- rep(names(n), n)
  phi <- 2 * pi * (c(A = 16, B = 18)[g] + rnorm(sum(n), 0, 0.4)) / 24
  A <- c(A = 8, B = 9)[g] + rnorm(sum(n), 0, 0.8)
  list(covariates = data.frame(AGE = factor(g)),
       harmonic_model = list(approach = "two_stage", period = 24, n_harmonics = 1,
                             trend_type = "none", individual_fits = list(),
                             group_var_name = "AGE", time_vec = 0:23, origin_shift = 0,
                             individual_params = data.frame(
                               subject = seq_len(sum(n)), mesor = 50 + rnorm(sum(n), 0, 2),
                               beta_cos_1 = A * cos(phi), beta_sin_1 = A * sin(phi),
                               amplitude_1 = A, acrophase_rad_1 = phi %% (2 * pi),
                               acrophase_time_1 = (phi %% (2 * pi)) * 24 / (2 * pi),
                               r_squared = runif(sum(n), .7, .95), p_value = 1e-4)))
}

# ================================================ P13.1 the shadowed name =====
test_that("P13.1: a cosinor comparison result does not erase the fANOVA post-hoc", {
  e <- ts_env()
  cm <- mk_cos_model()
  vals <- base_vals(pairwise_results = mk_posthoc(),
                    covariates = cm$covariates, harmonic_model = cm$harmonic_model)
  txt <- paste(e$dance_apa_report(vals, list(), "T"), collapse = "\n")

  # the exact sentence the defect produced
  expect_false(grepl("0 of 0 comparisons remained significant", txt, fixed = TRUE))
  expect_false(grepl("after the -- correction", txt, fixed = TRUE))
  expect_false(grepl("alpha = --", txt, fixed = TRUE))
  # and what it must say instead
  expect_true(grepl("1 of 3 comparisons remained significant after the Bonferroni correction at alpha = .05",
                    txt, fixed = TRUE))
  # both post-hoc tables are present and distinct
  expect_true(grepl("Post-hoc pairwise comparisons", txt, fixed = TRUE))
  expect_true(grepl("**Group differences (AGE).**", txt, fixed = TRUE))
  expect_true(grepl("Complete coefficient vector", txt, fixed = TRUE))
})

test_that("P13.1: the fANOVA post-hoc name is not rebound anywhere", {
  src <- code_of("server/93_apa_report.R")
  # exactly one assignment to `pw`: the fANOVA object
  hits <- gregexpr("(?<![A-Za-z0-9_.])pw *<- *", src, perl = TRUE)[[1]]
  expect_equal(sum(hits > 0), 1L)
  expect_true(grepl("pw <- values$pairwise_results", src, fixed = TRUE))
})

# =========================================== P13.2 the caveat matches its fit ==
test_that("P13.2: the FoSR caveat describes the estimator that was used", {
  e <- rep_env()
  gam_vals <- base_vals(reg_model = list(
    beta.hat = matrix(0.1, 2, 20, dimnames = list(c("(Intercept)", "Age"), NULL)),
    inference = "gam-prediction-contrast", method = "GAM",
    r2_t = runif(20, .05, .3)))
  gtxt <- paste(e$dance_apa_report(gam_vals, list(), "T"), collapse = "\n")
  expect_false(grepl("Each time point was fitted separately", gtxt, fixed = TRUE))
  expect_false(grepl("pointwise intervals are not simultaneous bands", gtxt, fixed = TRUE))
  expect_true(grepl("prediction contrasts read off the fitted model", gtxt, fixed = TRUE))
  expect_true(grepl("smoothness you see is imposed by the penalty", gtxt, fixed = TRUE))

  ols_vals <- base_vals(reg_model = list(
    beta.hat = matrix(0.1, 2, 20, dimnames = list(c("(Intercept)", "Age"), NULL)),
    beta.p = matrix(0.01, 2, 20), inference = "analytic-t-fdr", method = "OLS",
    r2_t = runif(20, .05, .3)))
  otxt <- paste(e$dance_apa_report(ols_vals, list(), "T"), collapse = "\n")
  expect_true(grepl("Each time point was fitted separately", otxt, fixed = TRUE))
  expect_true(grepl("not simultaneous bands", otxt, fixed = TRUE))
})

# ================================================ P13.3 missing data ==========
test_that("P13.3: missingness is reported with or without a smoothing step", {
  e <- rep_env()
  d <- matrix(rnorm(30 * 24), 30, 24)
  clean <- paste(e$dance_apa_report(list(data = d), list(), "T"), collapse = "\n")
  expect_true(grepl("No values were missing", clean, fixed = TRUE))
  expect_true(grepl("720 participant-by-time cells", clean, fixed = TRUE))

  d2 <- d; d2[1, 1:3] <- NA; d2[5, 7] <- NA
  gappy <- paste(e$dance_apa_report(list(data = d2), list(), "T"), collapse = "\n")
  expect_true(grepl("4 (0.6%) were missing, affecting 2 of 30 curves", gappy, fixed = TRUE))

  # when the smoother HAS run, the richer split is used instead and the raw
  # count is not repeated
  fs <- rep("observed", 720); fs[1:10] <- "interpolated"; fs[11:12] <- "extrapolated"
  sm <- paste(e$dance_apa_report(list(data = d2, fill_status = fs), list(), "T"),
              collapse = "\n")
  expect_true(grepl("extrapolated beyond a curve's first or last observation", sm, fixed = TRUE))
  expect_false(grepl("No values were missing", sm, fixed = TRUE))
  expect_false(grepl("were missing, affecting", sm, fixed = TRUE))
})

# ============================== P13.4 which curves did the analysis run on ====
test_that("P13.4: the fANOVA states which curves it was computed on", {
  e <- rep_env()
  fa <- list(F_stat = rep(3, 10), p_values_adjusted = rep(.2, 10),
             eta_squared = rep(.1, 10), df_between = 3, df_within = 122,
             group_names = c("A","B"), group_sizes = c(15, 15),
             n_permutations = 200, alpha = .05, design = "between")
  wr <- list(method = "linear_shift", time_points = seq(0, 1, length.out = 50))

  on_reg <- base_vals(warping_results = wr,
                      fanova_results = c(fa, list(used_warped = TRUE,
                                                  warp_method_used = "linear_shift")))
  t1 <- paste(e$dance_apa_report(on_reg, list(), "T"), collapse = "\n")
  expect_true(grepl("computed on the REGISTERED curves", t1, fixed = TRUE))
  expect_true(grepl("difference in AMPLITUDE at aligned time", t1, fixed = TRUE))

  on_raw <- base_vals(warping_results = wr,
                      fanova_results = c(fa, list(used_warped = FALSE)))
  t2 <- paste(e$dance_apa_report(on_raw, list(), "T"), collapse = "\n")
  expect_true(grepl("computed on the UNREGISTERED curves", t2, fixed = TRUE))
  expect_true(grepl("timing, in amplitude, or in both", t2, fixed = TRUE))

  # with no registration in the session the question does not arise
  t3 <- paste(e$dance_apa_report(base_vals(fanova_results = fa), list(), "T"), collapse = "\n")
  expect_false(grepl("REGISTERED curves", t3, fixed = TRUE))
  expect_false(grepl("UNREGISTERED curves", t3, fixed = TRUE))
})

test_that("P13.4: the app records what was used, not what was asked for", {
  src <- code_of("server/50_fanova.R")
  expect_true(grepl("used_warped <- FALSE", src, fixed = TRUE))
  expect_true(grepl("values$fanova_results$used_warped <- used_warped", src, fixed = TRUE))
  # the flag is set only inside the branches that actually take the warped path
  expect_equal(length(gregexpr("used_warped <- TRUE", src, fixed = TRUE)[[1]]), 2L)
})

test_that("P13.4: landmark registration states how its landmarks were obtained", {
  e <- rep_env()
  auto <- base_vals(warping_results = list(
    method = "landmark", n_rejected = 0L,
    time_points = seq(0, 1, length.out = 50)))
  atxt <- paste(e$dance_apa_report(auto, list(), "T"), collapse = "\n")
  expect_true(grepl("detected automatically", atxt, fixed = TRUE))
  expect_true(grepl("no curve required that here", atxt, fixed = TRUE))

  man <- base_vals(warping_results = list(
    method = "landmark", n_rejected = 2L, landmarks_used = c(0.25, 0.6),
    time_points = seq(0, 1, length.out = 50)))
  mtxt <- paste(e$dance_apa_report(man, list(), "T"), collapse = "\n")
  # keeps the leading zero: a landmark is a POSITION ON THE TIME AXIS, not a
  # statistic bounded by 1, so APA 6.36 does not apply to it
  expect_true(grepl("placed by hand at 0.250, 0.600", mtxt, fixed = TRUE))
  expect_true(grepl("2 curves yielded such landmarks", mtxt, fixed = TRUE))
  expect_true(grepl("left UNREGISTERED", mtxt, fixed = TRUE))
})
