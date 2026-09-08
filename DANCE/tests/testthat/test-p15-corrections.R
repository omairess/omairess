# ==============================================================================
# tests/testthat/test-p15-corrections.R — four things the user asked for
#
#   P15.1  every plot had a fixed pixel height chosen once by whoever wrote the
#          tab; 47 of them, none resizable by the reader.
#
#   P15.2  the "Diagnostics" help text still described the pre-P14 grid --
#          "trend in {none, linear, saturating} x harmonics in {1,2,3} ... 9
#          fits per subject" -- none of which had been true since P14. Worse,
#          the nested-model checkbox and the tau field sat under one heading
#          with nothing between them, so the tau value read as if it fed the
#          Delta-AICc table. It does not: that table estimates tau FREELY in
#          every saturating-exponential cell.
#
#   P15.3  a cleared tau field yields NA, not NULL, so `%||% 18` never fired and
#          is.finite(NA) was FALSE. The free-vs-fixed check was skipped and,
#          because the readout only prints that line when a result exists, it
#          vanished with nothing said. "I left the box empty" and "this was
#          tested and found nothing" looked identical on screen.
#
#   P15.4  the fitted model's own AIC/AICc/BIC were not reported anywhere. They
#          had been removed on the argument that "with no competing model they
#          are constant offsets of one another and carry no information", which
#          was right when it was made: there is a competing set now, and the
#          free-vs-fixed tau check differences the free-tau AIC, so quoting the
#          difference while withholding both terms leaves nothing to check.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

raw_of <- function(f) paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}

# ============================================== P15.1 resizable plot boxes ====
test_that("P15.1: the theme ships the resizable wrapper, css and behaviour", {
  css <- raw_of("ui/00_theme.R")
  expect_true(grepl(".fck-resizable", css, fixed = TRUE))
  expect_true(grepl("resize: both;", css, fixed = TRUE))
  # a drag must not push a plot outside its dashboard box, nor collapse it
  expect_true(grepl("max-width: 100%;", css, fixed = TRUE))
  expect_true(grepl("min-height: 140px;", css, fixed = TRUE))
  # both output families are wrapped
  for (sel in c(".shiny-plot-output", ".html-widget-output", ".plotly.html-widget"))
    expect_true(grepl(sel, css, fixed = TRUE), info = sel)
  # plotly is resized directly; base-R plots need Shiny to recompute sizes
  expect_true(grepl("Plotly.Plots.resize", css, fixed = TRUE))
  expect_true(grepl("new Event('resize')", css, fixed = TRUE))
  # the observer fires per pixel of a drag, so the re-render must be debounced
  expect_true(grepl("setTimeout", css, fixed = TRUE))
  expect_true(grepl("clearTimeout", css, fixed = TRUE))
  # lazily-rendered tabs must be picked up when Shiny paints them
  expect_true(grepl("shiny:value", css, fixed = TRUE))
  # and wrapping must be idempotent -- scan() runs on every output paint
  expect_true(grepl("dataset.fckResizable === '1'", css, fixed = TRUE))
  expect_true(grepl("closest('.fck-resizable')", css, fixed = TRUE))
  # no new dependency: this app is not environment-pinned (see the README)
  expect_false(grepl("shinyjqui", raw_of("app.R"), fixed = TRUE))
})

# ================================================ P15.2 the diagnostics text ==
test_that("P15.2: the stale grid description is gone from the UI", {
  # Comment LINES are dropped, not "# to end of line": the UI file carries hex
  # colours and HTML, and a blanket strip would delete the very strings this
  # guard checks for, making it pass vacuously. Dropping whole comment lines
  # removes the AUDIT note that quotes the old text -- which is what this file
  # is here to catch -- while leaving every string of code intact.
  ui_lines <- readLines(file.path(app_dir, "ui/72_harmonic.R"), warn = FALSE)
  ui <- paste(ui_lines[!grepl("^\\s*#", ui_lines)], collapse = "\n")
  expect_false(grepl("harmonics \\u2208 {1,2,3}", ui, fixed = TRUE))
  expect_false(grepl("9 fits per subject", ui, fixed = TRUE))
  # it is rendered from the model set now, so it cannot go stale again
  expect_true(grepl('uiOutput("harmonic_model_selection_help")', ui, fixed = TRUE))
  expect_true(grepl("dance_cosinor_model_set(k)",
                    code_of("server/72_harmonic.R"), fixed = TRUE))
})

test_that("P15.2: the two diagnostics are presented as separate checks", {
  ui <- raw_of("ui/72_harmonic.R")
  expect_true(grepl("1. Which specification?", ui, fixed = TRUE))
  expect_true(grepl("Is \\u03c4 identified?", ui, fixed = TRUE))
  # and the tau help says explicitly that it does not feed the table
  expect_true(grepl("not part of the \\u0394AICc table", ui, fixed = TRUE))
  expect_true(grepl("estimates \\u03c4 freely", ui, fixed = TRUE))
})

test_that("P15.2: the help text is built from the same set that gets fitted", {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/08_helpers_cosinor.R"), local = e)
  # what the help text will claim, for each selection
  expect_equal(length(e$dance_cosinor_model_set(1)), 4L)
  expect_equal(length(e$dance_cosinor_model_set(2)), 8L)
  expect_equal(length(e$dance_cosinor_model_set(3)), 12L)
  labs <- vapply(e$dance_cosinor_model_set(2), function(m) m$label, character(1))
  expect_true(all(grepl("^(none|linear|log|exp_sat) \\+ H1(-H2)?$", labs)))
})

# ============================================ P15.3 an empty tau field =========
test_that("P15.3: a cleared tau field is recorded, not silently dropped", {
  src <- code_of("server/72_harmonic.R")
  # the bug: `%||%` cannot rescue NA, which is what a cleared numericInput gives
  expect_false(grepl("as.numeric(input$harmonic_tau_fixed %||% 18)", src, fixed = TRUE))
  expect_true(grepl("conditioning$tau_fixed_skipped <- TRUE", src, fixed = TRUE))
  expect_true(grepl("Free-tau vs fixed-tau: NOT RUN", src, fixed = TRUE))

  # and the arithmetic of the guard itself
  for (v in list(NULL, NA, NA_real_, "", "abc", 0, -3)) {
    tau <- suppressWarnings(as.numeric(if (is.null(v)) NA else v))
    expect_false(is.finite(tau) && tau > 0)
  }
  for (v in list(18, "18", 24.5)) {
    tau <- suppressWarnings(as.numeric(v))
    expect_true(is.finite(tau) && tau > 0)
  }
})

test_that("P15.3: the readout says what the Delta-AIC is a difference of", {
  src <- raw_of("server/72_harmonic.R")
  expect_true(grepl("Delta-AIC = AIC(free tau) - AIC(tau fixed)", src, fixed = TRUE))
  expect_true(grepl("separate check from the Delta-AICc", src, fixed = TRUE))
})

# ==================================== P15.4 criteria for the fitted model ======
test_that("P15.4: the fitted model's own AIC/AICc/BIC are reported", {
  src <- raw_of("server/72_harmonic.R")
  expect_true(grepl("Information criteria for the fitted model", src, fixed = TRUE))
  # all three, from the per-subject values the fitters already return
  expect_true(grepl('for(nmi in c("aic", "aicc", "bic"))', src, fixed = TRUE))
  # with n and k stated, so the reader can check the correction themselves
  expect_true(grepl("dance_model_npar(trend_type, nh) + 1L", src, fixed = TRUE))
  # and the part of the old argument that is still true is kept
  expect_true(grepl("they separate models, not subjects", src, fixed = TRUE))
  # tau is flagged as freely estimated, which is what the user asked about
  expect_true(grepl("tau estimated freely", src, fixed = TRUE))
})

test_that("P15.4: parameter counts used by the readout are the shared rule", {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/08_helpers_cosinor.R"), local = e)
  # exp_sat + H1 = MESOR + A_sat + tau + cos1 + sin1 = 5, plus sigma = 6
  expect_equal(e$dance_model_npar("exp_sat", 1) + 1L, 6L)
  expect_equal(e$dance_model_npar("none", 1) + 1L, 4L)
})
