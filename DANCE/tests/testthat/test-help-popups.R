# Explanations sit behind an (i) (dance_help), status messages stay in plain
# view (dance_status), and every helpText() in the app is routed to the (i).
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
env <- new.env(parent = globalenv())
eval(parse(file.path(app_dir, "server/00a_help.R"), encoding = "UTF-8"), envir = env)

test_that("an explanation is hidden behind an (i); a status is shown", {
  h <- as.character(env$helpText("Some explanation of the setting."))
  expect_match(h, "dance-help-icon", fixed = TRUE)
  expect_match(h, "display:none", fixed = TRUE)
  expect_match(h, "Some explanation of the setting.", fixed = TRUE)
  s <- as.character(env$dance_status("Run the fPCA first."))
  expect_false(grepl("dance-help", s, fixed = TRUE))
  expect_match(s, "help-block", fixed = TRUE)
})

test_that("the UI loads the helpers and the pop-up script", {
  th <- paste(readLines(file.path(app_dir, "ui/00_theme.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(th, "00a_help.R", fixed = TRUE)
  expect_match(th, "window.danceHelp", fixed = TRUE)
  expect_match(th, "dance-help-modal", fixed = TRUE)
})

test_that("messages that must be seen use dance_status", {
  code <- function(f) paste(readLines(file.path(app_dir, "server", f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(code("42_fpca_anova.R"), 'dance_status("Run the fPCA (or time-warped fPCA) first.")', fixed = TRUE)
  expect_match(code("76_dlmo.R"), 'dance_status(sprintf("Hockey-stick: dlmoR', fixed = TRUE)
  expect_match(code("72_harmonic.R"), 'dance_status("No testable design terms in this model.")', fixed = TRUE)
  expect_match(code("72_harmonic.R"), 'dance_status(HTML(sprintf("<b>Warning:</b> no time is read', fixed = TRUE)
})
