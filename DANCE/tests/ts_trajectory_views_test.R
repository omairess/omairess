# ==============================================================================
# tests/ts_trajectory_views_test.R -- two-stage "Estimated trajectories"
#
# With a trend in the model, the two-stage approach offers the views the mixed
# approach does: the full trajectory, the trend alone (the non-periodic change
# from the first observation), the rhythm with the trend removed, and the
# harmonics alone. It used to offer only the first and third, so trend
# differences between groups could not be seen; and overlaying the raw data
# failed with "object 'tv' not found".
#
# Run with:   Rscript tests/ts_trajectory_views_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({ library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "shinyWidgets")) if (requireNamespace(p, quietly = TRUE)) library(p, character.only = TRUE) })
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  pickerInput <- function(inputId, label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(), ...)
    selectInput(inputId, label, choices = choices, selected = selected, multiple = multiple)
  updatePickerInput <- function(session, inputId, label = NULL, selected = NULL, choices = NULL, ...)
    updateSelectInput(session, inputId, label = label, choices = choices, selected = selected)
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }

# two groups whose saturating rise differs (A = 10 vs 25), same rhythm
set.seed(5)
clock <- c(10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 14, 16)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
grp <- factor(rep(c("0", "1"), each = 12))
raw <- t(vapply(seq_along(grp), function(i)
  15 + rnorm(1, 0, 2) + (if (grp[i] == "1") 25 else 10) * (1 - exp(-tt / 12)) +
    4 * cos(2 * pi * (lin - 18) / 24) + rnorm(length(tt), 0, 1.5), numeric(length(tt))))
rownames(raw) <- sprintf("S%02d", seq_along(grp))
server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data <- raw; values$time_labels <- sprintf("%02d:00", clock); values$time_clock <- clock
  values$covariates <- data.frame(Group = grp); values$group_labels <- grp
  values$group_variables <- list(Group = grp); values$subject_ids <- rownames(raw)
}
testServer(server, {
  session$setInputs(harmonic_approach = "two_stage", harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = "exp_sat", harmonic_traj_tau_mode = "fixed", harmonic_traj_tau = 12,
                    harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_time_var = "_shared_",
                    harmonic_design = "between", harmonic_between_var = "Group",
                    harmonic_traj_component = "full", harmonic_traj_raw = TRUE, run_harmonic = 1)
  session$flushReact()
  p <- tryCatch(output$harmonic_traj_curves, error = function(e) paste("ERROR", conditionMessage(e)))
  chk(!is.character(p) || !grepl("ERROR", p), "the full trajectory with the raw data overlaid renders",
      sprintf("overlay: %s", p))
  chk(identical(harmonic_traj_component(), "full"), "full is the default view", "the view is not full")
  session$setInputs(harmonic_traj_component = "trend"); session$flushReact()
  chk(identical(harmonic_traj_component(), "trend"), "the trend view is offered under two-stage",
      "the trend view is refused under two-stage")
  gc <- ts_group_curves("trend")
  d0 <- gc$table[gc$table$cell == "0", ]; d1 <- gc$table[gc$table$cell == "1", ]
  chk(abs(d0$fit[1]) < 1e-8 && abs(d1$fit[1]) < 1e-8,
      "the trend starts at zero at the first observation", "the trend does not start at zero")
  chk(tail(d1$fit, 1) - tail(d0$fit, 1) > 8,
      sprintf("the groups' trends differ (end of recording: %.1f vs %.1f)", tail(d0$fit, 1), tail(d1$fit, 1)),
      "the trend view does not show the groups' different rise")
  p <- tryCatch(output$harmonic_traj_curves, error = function(e) paste("ERROR", conditionMessage(e)))
  chk(!is.character(p) || !grepl("ERROR", p), "the trend view renders", sprintf("trend view: %s", p))
  for (cmp in c("baseline_harm", "harmonics")) {
    g2 <- ts_group_curves(cmp)
    chk(!is.null(g2) && all(is.finite(g2$table$fit)), sprintf("the '%s' view computes", cmp), sprintf("'%s' failed", cmp))
  }
  h <- ts_group_curves("harmonics")$table
  chk(abs(mean(h$fit[h$cell == "0"])) < 1.5, "the harmonics view is centred on zero", "the harmonics view is not centred")
  # the old name still works
  chk(identical(ts_group_curves("rhythm")$component, "baseline_harm"), "'rhythm' still means trend removed",
      "'rhythm' no longer maps")
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
