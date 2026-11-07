# The residual-model control is shown greyed out under the two-stage approach:
# disabled radios that are not a Shiny input, with the reason, and the live
# control (harmonic_traj_resid) stays mixed-only, so a two-stage run cannot be
# given AR(1) residuals.
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("ui")) "." else if (dir.exists("../../ui")) "../.." else "DANCE"
test_that("two-stage shows the residual control disabled, with the reason", {
  ui <- paste(readLines(file.path(app_dir, "ui/72_harmonic.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  i <- regexpr("input.harmonic_approach == 'two_stage'\",\n              tags\\$div\\(class = \"form-group dance-resid-disabled\"", ui)
  expect_true(i > 0)
  blk <- substr(ui, i, i + 1600)
  expect_equal(lengths(regmatches(blk, gregexpr("disabled = NA", blk))), 2L)
  expect_false(grepl("radioButtons(", blk, fixed = TRUE))
  expect_match(blk, "Not available for the two-stage cosinor", fixed = TRUE)
  # the live control sits in the mixed-only panel
  j <- regexpr("condition = \"input.harmonic_approach == 'mixed'\",\n              radioButtons\\(\"harmonic_traj_resid\"", ui)
  expect_true(j > 0)
  # and the server only ever passes AR(1) residuals to a mixed run
  srv <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(srv, 'resid_model = if (identical(approach, "mixed") &&', fixed = TRUE)
})
