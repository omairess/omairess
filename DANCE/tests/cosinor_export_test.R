# ==============================================================================
# tests/cosinor_export_test.R -- the exported script reruns the cosinor analysis
# in R, without the app, and gets the app's numbers
#
# The cosinor section of the exported script used to refit each participant and
# stop. A mixed-effects run -- the default approach -- was not reproduced at all,
# and the section read its data as `raw_data`, an object the script never
# defined, so it could not run in the first place. This presses Run in the app,
# writes the script, runs it in a SEPARATE R process with nothing of the app
# loaded (Rscript --vanilla, the script and the data only), and compares:
#
#   mixed-effects, saturating trend, tau profiled within a range, A and tau per
#   cell: tau, the fixed effects, the per-cell intercept/MESOR/change, every
#   design test, the per-cell A, tau, A/tau;
#   two-stage, saturating trend, tau within a range, groups, bootstrap, nested
#   set: the parameter table, the inclusion audit, the group summaries, the
#   MANOVA and component tests, the bootstrap intervals, the Akaike table.
#
# What runs is the script's header, its libraries and data loading (sections 1
# and 2) and the cosinor section (10), as written. The sections between them
# rerun other analyses -- the smoothing, for one, which needs fda -- and are
# not what is checked here. The one change: `library(fda)` is skipped when fda
# is not installed (it cannot be installed in the sandbox this test was written
# in). Nothing in the cosinor section uses fda.
#
# Run with:   Rscript tests/cosinor_export_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "mgcv", "ggplot2")) if (requireNamespace(p, quietly = TRUE))
    library(p, character.only = TRUE)
})
if (requireNamespace("shinyWidgets", quietly = TRUE)) {
  suppressPackageStartupMessages(library(shinyWidgets))
} else {
  pickerInput <- function(inputId, label = NULL, choices = NULL, selected = NULL,
                          multiple = FALSE, options = list(), ...)
    selectInput(inputId, label, choices = choices, selected = selected, multiple = multiple)
  updatePickerInput <- function(session, inputId, label = NULL, selected = NULL,
                                choices = NULL, ...)
    updateSelectInput(session, inputId, label = label, choices = choices, selected = selected)
}

app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))

failures <- 0L
fail <- function(...) { cat("FAIL:", ..., "\n"); failures <<- failures + 1L }
ok   <- function(...) cat("ok  :", ..., "\n")
same <- function(label, a, b, tol = 1e-6) {
  a <- suppressWarnings(as.numeric(unlist(a))); b <- suppressWarnings(as.numeric(unlist(b)))
  if (length(a) != length(b) || !length(a)) return(fail(label, ": lengths", length(a), "vs", length(b)))
  if (!identical(is.na(a), is.na(b))) return(fail(label, ": NA pattern differs"))
  d <- if (all(is.na(a))) 0 else max(abs(a[!is.na(a)] - b[!is.na(b)]) / pmax(1, abs(a[!is.na(a)])))
  if (d <= tol) ok(sprintf("%s (%d values, max rel. diff %.2g)", label, length(a), d))
  else fail(sprintf("%s: max relative difference %.3g", label, d))
}

# ---- data: 08:00 to 15:00 the next day, two groups, a saturating rise -------
set.seed(21)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 15)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
grp <- factor(rep(c("A", "B"), each = 14))
T0 <- c(A = 10, B = 16)
# residuals autocorrelated within a curve (continuous-time AR(1), 0.5 per hour),
# so the AR(1) scenario exercises the whitened fit
Lc <- t(chol(0.5^abs(outer(tt, tt, "-"))))
raw <- t(vapply(seq_along(grp), function(i) {
  g <- as.character(grp[i])
  20 + rnorm(1, 0, 1.5) + (12 + rnorm(1, 0, 2)) * (1 - exp(-tt / T0[[g]])) +
    5 * cos(2 * pi * (lin - 16) / 24) + as.numeric(Lc %*% rnorm(length(tt), 0, 1))
}, numeric(length(tt))))
rownames(raw) <- sprintf("S%02d", seq_along(grp))
raw[3, 1:2] <- NA                   # a participant missing the first two observations
covs <- data.frame(Group = grp, stringsAsFactors = FALSE)

tmp <- tempfile("cosinor_export_"); dir.create(tmp)

server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data <- raw
  values$time_labels <- sprintf("%02d:00", clock)
  values$time_clock <- clock
  values$covariates <- covs
  values$group_labels <- grp
  values$group_variables <- list(Group = grp)
  values$subject_ids <- rownames(raw)
}

run <- function(session, approach, run_id, ...) {
  session$setInputs(harmonic_approach = approach, harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = "exp_sat", harmonic_time_var = "_shared_",
                    harmonic_dv_name = "KSS", harmonic_dv_units = "KSS points",
                    harmonic_design = "between", harmonic_between_var = "Group",
                    harmonic_df_method = "kr", harmonic_traj_adjust = "holm", ...,
                    run_harmonic = run_id)
  session$flushReact()
}

# Write the script, run it in a clean R process, return what it computed.
run_script <- function(code, tag, keep) {
  code <- strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  # sections 1-2 and 10, each from its header's rule line to the next one's
  hdr <- grep("^# [0-9]+\\. ", code)
  s3 <- hdr[grepl("^# 3\\. ", code[hdr])][1]
  s10 <- hdr[grepl("^# 10\\. ", code[hdr])][1]
  nxt <- hdr[hdr > s10][1]
  code <- c(code[seq_len(s3 - 2L)],
            code[(s10 - 1L):(if (is.na(nxt)) length(code) else nxt - 2L)])
  if (!requireNamespace("fda", quietly = TRUE))
    code <- sub("^library\\(fda\\)", "# library(fda)   [not installed in the test sandbox]", code)
  f <- file.path(tmp, paste0(tag, ".R")); writeLines(code, f)
  out <- file.path(tmp, paste0(tag, "_script.rds"))
  drv <- file.path(tmp, paste0(tag, "_driver.R"))
  writeLines(c(
    '.libPaths(c("~/Rlib", .libPaths()))',
    sprintf('dance_input_data <- readRDS("%s")', file.path(tmp, "raw.rds")),
    'pdf(NULL)',
    sprintf('source("%s", echo = FALSE)', f),
    sprintf('saveRDS(mget(c(%s), ifnotfound = list(NULL)), "%s")',
            paste0('"', keep, '"', collapse = ", "), out),
    'cat("SCRIPT-RAN-CLEAN\\n")'), drv)
  res <- system2("Rscript", c("--vanilla", shQuote(drv)), stdout = TRUE, stderr = TRUE)
  if (!any(grepl("SCRIPT-RAN-CLEAN", res, fixed = TRUE))) {
    fail(tag, ": the exported script did not run in a clean R session")
    cat(paste0("      | ", tail(res, 25), collapse = "\n"), "\n")
    return(NULL)
  }
  ok(sprintf("%s: the exported script ran in a clean R session (%d lines)", tag, length(code)))
  readRDS(out)
}
saveRDS(raw, file.path(tmp, "raw.rds"))

testServer(server, {
  # ====================================================== mixed-effects
  cat("\n-- mixed-effects: tau profiled within 2 to 48 h, A and tau per cell ----\n")
  run(session, "mixed", 1, harmonic_traj_tau_mode = "profile", harmonic_traj_tau = 18,
      harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_traj_sat_cells = TRUE)
  ff <- values$harmonic_model$traj
  if (!isTRUE(ff$ok)) fail("the app's mixed fit failed:", ff$message %||% "")
  app_tests <- dance_traj_tests_table(harmonic_traj_tests())
  sc <- run_script(generate_analysis_code(), "mixed", c("traj_fit", "traj_tests"))
  if (!is.null(sc)) {
    tf <- sc$traj_fit
    chk <- function(cond, good, bad) if (isTRUE(cond)) ok(good) else fail(bad)
    chk(identical(tf$tau_mode, ff$tau_mode),
        sprintf("tau found the same way (%s)", ff$tau_mode), "the script's tau mode differs")
    same("the trajectory model's tau", ff$spec$tau, tf$spec$tau)
    same("the fixed effects", lme4::fixef(ff$model), lme4::fixef(tf$model))
    same("the per-cell intercept, MESOR and change",
         dance_traj_level_summary(ff)$table[, c("estimate", "se")],
         dance_traj_level_summary(tf)$table[, c("estimate", "se")])
    sct <- dance_traj_tests_table(sc$traj_tests)
    chk(identical(app_tests$test, sct$test), "the same design tests, in the same order",
        "the script runs different tests")
    same("every design test's F", app_tests$F, sct$F)
    same("every design test's denominator df", app_tests$df2, sct$df2)
    same("every design test's p", app_tests$p, sct$p)
    same("A, tau and A/tau per cell",
         ff$sat$table[, c("A", "tau", "init_slope")], tf$sat$table[, c("A", "tau", "init_slope")])
    same("one tau vs per cell (likelihood ratio)", ff$sat$lrt$statistic, tf$sat$lrt$statistic)
  }

  # ====================================================== mixed, AR(1) residuals
  cat("\n-- mixed-effects with continuous-time AR(1) residuals, tau held ---------\n")
  run(session, "mixed", 3, harmonic_traj_tau_mode = "fixed", harmonic_traj_tau = 12,
      harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_traj_sat_cells = FALSE,
      harmonic_traj_resid = "car1")
  fa <- values$harmonic_model$traj
  if (!isTRUE(fa$ok)) fail("the app's AR(1) fit failed:", fa$message %||% "")
  if (isTRUE(fa$phi > 0)) ok(sprintf("the app estimated phi = %.3f per hour", fa$phi))
  else fail("phi-hat is 0: the whitened fit is not exercised")
  app_tests_a <- dance_traj_tests_table(harmonic_traj_tests())
  sc <- run_script(generate_analysis_code(), "mixed_car1", c("traj_fit", "traj_tests"))
  if (!is.null(sc)) {
    tf <- sc$traj_fit
    same("the residual autocorrelation phi", fa$residual_ar$phi, tf$residual_ar$phi)
    same("its profile interval", fa$residual_ar$ci, tf$residual_ar$ci)
    same("its test against independence", fa$residual_ar$lrt$statistic, tf$residual_ar$lrt$statistic)
    same("the fixed effects at phi-hat", lme4::fixef(fa$model), lme4::fixef(tf$model))
    sct <- dance_traj_tests_table(sc$traj_tests)
    same("every design test's F, df and p",
         app_tests_a[, c("F", "df2", "p")], sct[, c("F", "df2", "p")])
  }
  session$setInputs(harmonic_traj_resid = "independent")

  # ====================================================== two-stage
  cat("\n-- two-stage: tau within 2 to 48 h, groups, bootstrap, nested set ----\n")
  run(session, "two_stage", 2, harmonic_ts_tau_mode = "estimate", harmonic_traj_tau = 18,
      harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_bootstrap = TRUE,
      harmonic_n_boot = 150, harmonic_model_selection = TRUE)
  mod <- values$harmonic_model
  sc <- run_script(generate_analysis_code(), "two_stage",
                   c("ts_params", "ts_gate", "ts_df", "ts_model", "ts_components", "ts_boot",
                     "ts_models", "cosinor_params"))
  if (!is.null(sc)) {
    num_cols <- names(mod$individual_params)[vapply(mod$individual_params, is.numeric, logical(1))]
    same("the parameter table (every column of every summarised participant)",
         mod$individual_params[, num_cols], sc$ts_params[, num_cols])
    same("the inclusion audit",
         unlist(mod$fit_audit[c("n_attempted", "n_converged", "n_boundary", "n_nonconverged",
                                "n_failed", "n_summarised")]),
         unlist(sc$ts_gate$fit_audit[c("n_attempted", "n_converged", "n_boundary",
                                       "n_nonconverged", "n_failed", "n_summarised")]))
    ts <- harmonic_ts()
    same("the group intercept, MESOR and change",
         dance_ts_level_summary(ts$df, ts$mod)$table[, c("estimate", "lo", "hi")],
         dance_ts_level_summary(sc$ts_df, sc$ts_model)$table[, c("estimate", "lo", "hi")])
    same("A, tau and A/tau per group",
         dance_ts_sat_summary(ts$df, ts$mod)$table[, c("estimate", "lo", "hi")],
         dance_ts_sat_summary(sc$ts_df, sc$ts_model)$table[, c("estimate", "lo", "hi")])
    app_c <- dance_ts_components(ts$df, ts$mod)
    same("every component test's F",
         vapply(app_c, function(r) if (isTRUE(r$ok)) r$F else NA_real_, numeric(1)),
         vapply(sc$ts_components, function(r) if (isTRUE(r$ok)) r$F else NA_real_, numeric(1)))
    same("the bootstrap intervals (seeded)",
         mod$boot_results[c("mesor_ci", "amplitude_ci", "acrophase_ci")],
         sc$ts_boot[c("mesor_ci", "amplitude_ci", "acrophase_ci")])
    same("the nested set's AICc", mod$model_selection$AICc, sc$ts_models$AICc)
  }
})

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
