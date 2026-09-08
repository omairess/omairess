# ==============================================================================
# tests/testthat/test-p14-corrections.R — the nested model-selection set
#
# Reported by the user: the diagnostic fitted trend x harmonics as a full grid
# to three harmonics regardless of what had been selected, and higher harmonics
# should never be considered without the lower ones.
#
# Verified before changing anything: the harmonics were ALREADY cumulative.
# fit_cosinor() builds its design with `for (h in 1:n_harmonics)`, so
# n_harmonics = 2 fits cos1, sin1, cos2, sin2 and the set was properly nested.
# What was wrong was the LABEL -- "exp_sat + H2" reads as "harmonic 2 alone",
# which is a different model and not one this app can fit.
#
# Two things were genuinely wrong:
#
#   P14.1  the grid was hardcoded c("none", "linear", "exp_sat") and the UI
#          offers FOUR trends. `log` was never fitted, so a user who selected a
#          logarithmic trend was shown a Delta-AICc table that did not contain
#          their own model, with Akaike weights normalised over a set their
#          specification was excluded from. The weights were wrong, not merely
#          incomplete.
#
#   P14.2  the harmonic grid ran to 3 whatever had been selected, spending
#          compute on models the user had already ruled out and, worse,
#          spreading Akaike weight over them -- which changes the weight the
#          reported model receives.
#
# The candidate set now comes from one function, and the report finally says
# which specification was chosen and how.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  ln <- sub("#.*$", "", ln)
  paste(ln, collapse = "\n")
}
cos_env <- function() {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/08_helpers_cosinor.R"), local = e)
  e
}
labels_of <- function(ms) vapply(ms, function(m) m$label, character(1))

# ================================================= P14.2 the candidate set ====
test_that("P14.2: the set stops at the number of harmonics selected", {
  e <- cos_env()
  expect_equal(labels_of(e$dance_cosinor_model_set(1)),
               c("none + H1", "linear + H1", "log + H1", "exp_sat + H1"))
  # the user's own example: saturating exponential with 2 harmonics
  expect_equal(labels_of(e$dance_cosinor_model_set(2)),
               c("none + H1", "none + H1-H2",
                 "linear + H1", "linear + H1-H2",
                 "log + H1", "log + H1-H2",
                 "exp_sat + H1", "exp_sat + H1-H2"))
  expect_equal(length(e$dance_cosinor_model_set(3)), 12L)
  # nothing beyond the selection, ever
  for (k in 1:3) {
    ks <- vapply(e$dance_cosinor_model_set(k), function(m) m$k, integer(1))
    expect_lte(max(ks), k)
    expect_equal(sort(unique(ks)), seq_len(k))
  }
})

test_that("P14.2: no harmonic is ever offered without the ones below it", {
  e <- cos_env()
  for (k in 1:4) for (m in e$dance_cosinor_model_set(k)) {
    # the label names a cumulative set, and the design fit_cosinor() builds for
    # n_harmonics = m$k is harmonics 1..m$k, so the two agree by construction
    expect_equal(m$label, paste0(m$trend, " + ",
                                 if (m$k == 1) "H1" else paste0("H1-H", m$k)))
    expect_false(grepl("\\+ H[2-9]$", m$label))
  }
})

test_that("P14: fit_cosinor's design really is cumulative, as the labels claim", {
  src <- code_of("server/72_harmonic.R")
  # the design loop that makes "H1-H2" an honest label rather than a hope
  expect_true(grepl("for(h in 1:n_harmonics)", src, fixed = TRUE))
  expect_true(grepl('paste0("cos", h), paste0("sin", h)', src, fixed = TRUE))
})

# ================================================= P14.1 the missing trend ====
test_that("P14.1: every trend the UI offers is in the candidate set", {
  e <- cos_env()
  ui <- code_of("ui/72_harmonic.R")
  for (tt in e$DANCE_TREND_TYPES)
    expect_true(grepl(paste0('"', tt, '"'), ui, fixed = TRUE), info = tt)
  expect_true("log" %in% e$DANCE_TREND_TYPES)      # the one that was missing
  trends <- unique(vapply(e$dance_cosinor_model_set(2), function(m) m$trend, character(1)))
  expect_setequal(trends, e$DANCE_TREND_TYPES)
})

test_that("P14.1: the hardcoded grid is gone from the diagnostic", {
  src <- code_of("server/72_harmonic.R")
  expect_false(grepl('ms_trends <- c("none", "linear", "exp_sat")', src, fixed = TRUE))
  expect_false(grepl("ms_harms  <- 1:3", src, fixed = TRUE))
  expect_true(grepl("dance_cosinor_model_set(n_harmonics)", src, fixed = TRUE))
  # and the parameter count is one rule, not a switch repeated per call site
  expect_false(grepl('switch(trend_type, "none" = 0', src, fixed = TRUE))
  expect_true(grepl("dance_model_npar(", src, fixed = TRUE))
})

test_that("P14: parameter counts are right for every candidate", {
  e <- cos_env()
  # MESOR + trend + 2 per harmonic
  expect_equal(e$dance_model_npar("none", 1), 3L)
  expect_equal(e$dance_model_npar("linear", 1), 4L)
  expect_equal(e$dance_model_npar("log", 1), 4L)      # was silently counted as 0
  expect_equal(e$dance_model_npar("exp_sat", 1), 5L)
  expect_equal(e$dance_model_npar("exp_sat", 2), 7L)
  expect_equal(e$dance_trend_npar("log"), 1L)
  # a nested pair must differ by exactly the 2 parameters of one harmonic
  for (tt in e$DANCE_TREND_TYPES)
    expect_equal(e$dance_model_npar(tt, 2) - e$dance_model_npar(tt, 1), 2L)
})

# ================================================= the report says what ran ===
test_that("P14: the report carries the model set, the choice and its caveat", {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/93_apa_report.R"), local = e)
  ms <- data.frame(model = c("none + H1", "exp_sat + H1-H2"),
                   AICc = c(82.80, 93.85), dAICc = c(0, 11.05),
                   weight = c(0.505, 0.002), stringsAsFactors = FALSE)
  attr(ms, "selected") <- "exp_sat + H1-H2"
  attr(ms, "n_harmonics_max") <- 2L
  vals <- list(data = matrix(rnorm(30 * 24), 30, 24),
               harmonic_model = list(period = 24, n_harmonics = 2,
                                     trend_type = "exp_sat",
                                     individual_fits = list(),
                                     model_selection = ms))
  txt <- paste(e$dance_apa_report(vals, list(), "T"), collapse = "\n")

  expect_true(grepl("nested set of 2 candidates", txt, fixed = TRUE))
  expect_true(grepl("Harmonics are cumulative", txt, fixed = TRUE))
  expect_true(grepl("exp_sat + H1-H2 (reported)", txt, fixed = TRUE))
  # it must SAY when the reported model is not the best-supported one
  expect_true(grepl("is not the best-supported candidate", txt, fixed = TRUE))
  expect_true(grepl("11.05 AICc units behind", txt, fixed = TRUE))
  # the two caveats that make the table honest
  expect_true(grepl("conditional on the candidate set", txt, fixed = TRUE))
  expect_true(grepl("not adjusted for that selection", txt, fixed = TRUE))
  # Akaike weights are bounded by 1, so no leading zero
  expect_true(grepl("| .505 |", txt, fixed = TRUE))

  # and when the reported model IS the best, it says so instead
  attr(ms, "selected") <- "none + H1"
  vals$harmonic_model$model_selection <- ms
  t2 <- paste(e$dance_apa_report(vals, list(), "T"), collapse = "\n")
  expect_true(grepl("has the lowest AICc in this set", t2, fixed = TRUE))
  expect_false(grepl("is not the best-supported candidate", t2, fixed = TRUE))
})
