# ==============================================================================
# tests/testthat/test-harmonic-ui-structure.R
#
# P21 phase 4. The course correction was that the module should stay ONE module
# with an extended design capability -- not grow a parallel workflow. These
# checks pin that shape, because it is the kind of thing that erodes one
# well-meant addition at a time.
# ==============================================================================

app_dir <- if (dir.exists("ui")) "." else if (dir.exists("../../ui")) "../.." else "DANCE"
slurp <- function(f) paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
ui  <- slurp("ui/72_harmonic.R")
srv <- slurp("server/72_harmonic.R")

# The terminology checks below are about what a USER SEES, so they run against
# the file with comment lines removed. The first draft did not, and failed on
# the comment that explains why the old wording was dropped -- a test that
# forbids naming the thing you fixed makes the fix undocumentable.
strip_comments <- function(x)
  paste(grep("^\\s*#", strsplit(x, "\n", fixed = TRUE)[[1]], value = TRUE, invert = TRUE),
        collapse = "\n")
ui_visible <- strip_comments(ui)

test_that("the study design is named explicitly, for any number of levels", {
  expect_true(grepl('h4("Study Design")', ui, fixed = TRUE))
  for (d in c('"between"', '"within"', '"mixed"'))
    expect_true(grepl(d, ui, fixed = TRUE), info = d)
  # every selector the three panels need must be rendered
  for (o in c("harmonic_between_var_ui", "harmonic_between_var_ui2",
              "harmonic_within_var_ui", "harmonic_within_var_ui2",
              "harmonic_design_readout"))
    expect_true(grepl(o, ui, fixed = TRUE) && grepl(o, srv, fixed = TRUE), info = o)
})

test_that("nothing hard-codes a two-by-two", {
  # the design terms are read from the inputs, never enumerated
  expect_true(grepl("harmonic_design_terms", srv, fixed = TRUE))
  expect_false(grepl("2 x 2 only", ui, fixed = TRUE))
  expect_false(grepl("subject:Condition", srv, fixed = TRUE))
})

test_that("the classification is shown to the user BEFORE the fit", {
  # the readout must report what the DATA say, not merely echo the selection
  expect_true(grepl("Read from your data", srv, fixed = TRUE))
  for (role in c('"between"', '"within"', '"partial"', '"constant"'))
    expect_true(grepl(role, srv, fixed = TRUE), info = role)
})

test_that("a cosinor plus a trend is not called a two-process model", {
  expect_false(grepl("Two-Process Model", ui_visible, fixed = TRUE))
  expect_false(grepl("Homeostatic Trend Model", ui_visible, fixed = TRUE))
  expect_true(grepl("Non-periodic Trend Model", ui_visible, fixed = TRUE))
  # and the explanation says what the decomposition IS
  expect_true(grepl("descriptive decomposition", ui_visible, fixed = TRUE))
  # "Process S" may still be NAMED as the shape a trend often takes, but not
  # asserted as the mechanism being fitted
  expect_false(grepl("classic Process S", ui_visible, fixed = TRUE))
})

test_that("model selection is available but does not gate an ordinary run", {
  adv <- regmatches(ui, regexpr("tags\\$details.*?Run Harmonic Regression", ui))
  expect_true(length(adv) == 1)
  # the three things that moved are behind Advanced ...
  for (ctl in c("harmonic_model_selection", "harmonic_use_bounds",
                "harmonic_include_boundary"))
    expect_true(grepl(ctl, adv, fixed = TRUE), info = ctl)
  # ... and the things a normal run needs are NOT
  head_ui <- substr(ui, 1, regexpr("tags\\$details", ui))
  for (ctl in c("harmonic_period", "n_harmonics", "harmonic_trend_type",
                "harmonic_design", "harmonic_bootstrap"))
    expect_true(grepl(ctl, head_ui, fixed = TRUE), info = ctl)
})

test_that("the tab set is unchanged apart from the one rename", {
  tabs <- regmatches(ui, gregexpr('tabPanel\\("[^"]+"', ui))[[1]]
  tabs <- sub('tabPanel\\("', "", sub('"$', "", tabs))
  # "3. Parameter Distribution" was removed in the approach overhaul: its four
  # histograms were per-participant two-stage summaries with no mixed-effects
  # analogue, and the individual table carries the same numbers.
  expect_equal(tabs, c("1. Fitted Curves", "2. Polar Plot (Acrophase)",
                       "2b. Polar Density",
                       "4. Individual Results", "5. Residual Diagnostics",
                       "6. Group / Condition Comparison"))
})

test_that("no second analysis module was added alongside this one", {
  # the trajectory backend is helper files under server/, reachable from this
  # module -- it must not acquire a sidebar entry or a tab of its own
  sb <- paste(readLines(file.path(app_dir, "ui/00_theme.R"), warn = FALSE), collapse = "\n")
  for (f in list.files(file.path(app_dir, "ui"), pattern = "[.]R$", full.names = TRUE)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl('tabName = "trajectory"', src, fixed = TRUE), info = basename(f))
    expect_false(grepl('tabName = "mixed_cosinor"', src, fixed = TRUE), info = basename(f))
  }
})


test_that("tab 6 does not offer views or rows the fitted model cannot produce", {
  # Both guards were added after running on real data, where the selector
  # offered a component that draws a flat zero line and the decomposition
  # printed one test twice under two names.
  expect_true(grepl("harmonic_traj_component <- reactive", srv, fixed = TRUE))
  expect_true(grepl("updateSelectInput(session, \"harmonic_traj_component\"", srv, fixed = TRUE))
  # the decomposition's block list is conditional on the trend
  blk <- regmatches(srv, regexpr('want <- if \\(identical\\(ff\\$spec\\$trend, "none"\\)[^\n]*\n[^\n]*\n[^\n]*', srv))
  expect_true(length(blk) == 1)
  expect_true(grepl('c\\("full", "circadian", "level"\\)', blk))
})
