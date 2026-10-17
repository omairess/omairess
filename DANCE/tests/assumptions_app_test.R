# ==============================================================================
# tests/assumptions_app_test.R -- what a cosinor run assumes, through Run
#
# An audit of the cosinor module for choices made without a word. This presses
# Run and checks the ones that were wrong, and that every run now says what it
# assumed:
#   * clock times parsed at import are unwrapped at midnight (they went in as
#     hours of the day, so a second day came out BEFORE the first);
#   * a participant missing the first observation shares the t = 0 origin;
#   * the tau settings: an unusable range, or a held value outside the range a
#     flat profile would fall back to, is refused -- and said as it is typed;
#   * a time column of the wrong length is refused, not replaced by invented
#     equally spaced times;
#   * the bootstrap is seeded: the same Run gives the same interval;
#   * the nested-model table states the participants it averages over;
#   * the summary and the report list the run's assumptions.
#
# Run with:   Rscript tests/assumptions_app_test.R      (from the DANCE directory)
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
chk  <- function(cond, good, bad) if (isTRUE(cond)) ok(good) else fail(bad)
text_of <- function(x) paste(as.character(unlist(x)), collapse = " ")

# ---- data: 08:00 to 15:00 the next day, two groups, a linear rise ------------
set.seed(11)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 15)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
grp <- factor(rep(c("A", "B"), each = 12))
raw <- t(vapply(seq_along(grp), function(i)
  20 + rnorm(1, 0, 1.5) + 0.3 * tt + 5 * cos(2 * pi * (lin - 16) / 24) + rnorm(length(tt), 0, 1),
  numeric(length(tt))))
rownames(raw) <- sprintf("S%02d", seq_along(grp))
raw[3, 1:2] <- NA                   # one participant misses the first two observations
covs <- data.frame(Group = grp, AGE = round(runif(length(grp), 20, 60)), stringsAsFactors = FALSE)

server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data <- raw
  values$time_labels <- sprintf("%02d:00", clock)
  values$time_clock <- clock                       # hours of the day, as import parses them
  values$covariates <- covs
  values$group_labels <- grp
  values$group_variables <- list(Group = grp)
  values$subject_ids <- rownames(raw)
}

run <- function(session, approach, run_id, trend = "linear", time_var = "_shared_",
                tau_mode = "fixed", tau = 18, lo = 1, hi = 72, boot = FALSE, ms = FALSE,
                ts_mode = "estimate") {
  session$setInputs(harmonic_approach = approach, harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = trend, harmonic_traj_tau_mode = tau_mode,
                    harmonic_ts_tau_mode = ts_mode, harmonic_tau_lo = lo, harmonic_tau_hi = hi,
                    harmonic_traj_tau = tau, harmonic_traj_sat_cells = FALSE,
                    harmonic_model_selection = ms, harmonic_time_var = time_var,
                    harmonic_bootstrap = boot, harmonic_n_boot = 200,
                    harmonic_dv_name = "KSS", harmonic_dv_units = "KSS points",
                    harmonic_design = "between", harmonic_between_var = "Group",
                    harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                    run_harmonic = run_id)
  session$flushReact()
}

testServer(server, {
  # ====================================================== the clock route
  cat("\n-- shared clock times across midnight, two-stage ----------------------\n")
  run(session, "two_stage", 1)
  mod <- values$harmonic_model
  chk(!is.null(mod) && isTRUE(mod$wrap_applied) && all(diff(mod$time_vec) > 0) &&
        isTRUE(all.equal(mod$time_vec, tt)) && identical(mod$origin_shift, 8),
      "clock times parsed at import are unwrapped at midnight: t runs 0 to 31 h, origin 08:00",
      sprintf("the clock route was not unwrapped: %s", paste(mod$time_vec, collapse = ",")))
  P <- mod$individual_params
  s3 <- mod$individual_fits[[3]]
  chk(!is.null(s3) && isTRUE(s3$success) && identical(s3$t_offset %||% 0, 0),
      "the participant missing the first observations shares the t = 0 origin",
      sprintf("participant 3 anchored at %s", format(s3$t_offset)))
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("--- Assumptions and fixed settings ---", sm, fixed = TRUE) &&
        grepl("moved 24 h later", sm, fixed = TRUE) &&
        grepl("the clock times parsed from the column names at import", sm, fixed = TRUE),
      "the summary lists the assumptions, the time source and the midnight rule",
      "the summary does not list the run's assumptions")
  chk(grepl("Welch's t", sm, fixed = TRUE) && grepl("treated as exact", sm, fixed = TRUE),
      "the two-stage list names its tests and that estimates are treated as exact",
      "the two-stage tests are not stated")

  # ====================================================== the mixed model
  cat("\n-- the same, mixed-effects ---------------------------------------------\n")
  run(session, "mixed", 2)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && identical(ff$spec$t0, 0),
      "the mixed model is anchored at the recording's first time point",
      sprintf("mixed t0: %s (%s)", format(ff$spec$t0), ff$message %||% ""))
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("INDEPENDENT, with ONE variance", sm, fixed = TRUE) &&
        grepl("EQUAL weights", sm, fixed = TRUE),
      "the mixed summary states the residual model and the marginal weighting",
      "the residual model is not stated")
  rp <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
  chk(grepl("Assumptions and fixed settings of the cosinor analysis", rp, fixed = TRUE) &&
        grepl("- **Residuals.** Normal and INDEPENDENT", rp, fixed = TRUE),
      "the report's Methods list the same assumptions", "the report does not list the assumptions")
  html <- paste(dance_apa_html(dance_apa_report(values, input, "T")), collapse = "\n")
  chk(grepl("<ul>", html, fixed = TRUE) && grepl("<li><strong>Residuals.</strong>", html, fixed = TRUE),
      "and the HTML report renders them as a list", "the HTML report does not render the list")

  # ====================================================== tau settings
  cat("\n-- tau settings that cannot run are refused, and said as typed --------\n")
  before <- values$harmonic_model
  session$setInputs(harmonic_trend_type = "exp_sat", harmonic_approach = "mixed",
                    harmonic_traj_tau_mode = "profile", harmonic_tau_lo = 30,
                    harmonic_tau_hi = 60, harmonic_traj_tau = 18)
  session$flushReact()
  msg <- text_of(output$harmonic_tau_check)
  chk(grepl("outside the range", msg, fixed = TRUE),
      "a held value outside the range is flagged as it is typed", "no inline message for the held value")
  run(session, "mixed", 3, trend = "exp_sat", tau_mode = "profile", tau = 18, lo = 30, hi = 60)
  chk(identical(values$harmonic_model, before),
      "and Run refuses it: a flat profile would hold tau outside the range",
      "a held value outside the range was fitted anyway")
  session$setInputs(harmonic_tau_lo = 101, harmonic_tau_hi = 40); session$flushReact()
  msg <- text_of(output$harmonic_tau_check)
  chk(grepl("The range is not usable", msg, fixed = TRUE),
      "a range from 101 to 40 h is flagged as it is typed", "no inline message for an inverted range")
  session$setInputs(harmonic_tau_lo = 1, harmonic_tau_hi = 72); session$flushReact()
  chk(identical(text_of(output$harmonic_tau_check), ""),
      "a usable setting shows no message", "a usable setting shows a message")

  # ====================================================== a time column that is not one
  cat("\n-- a time column of the wrong length is refused -----------------------\n")
  run(session, "two_stage", 4, time_var = "AGE")
  chk(identical(values$harmonic_model, before),
      "a per-participant column is refused as the time variable, not replaced by equal spacing",
      "the run went ahead on invented times")

  # ====================================================== the bootstrap
  cat("\n-- the bootstrap is seeded ---------------------------------------------\n")
  run(session, "two_stage", 5, boot = TRUE)
  b1 <- values$harmonic_model$boot_results
  run(session, "two_stage", 6, boot = TRUE)
  b2 <- values$harmonic_model$boot_results
  chk(!is.null(b1) && identical(b1$mesor_ci, b2$mesor_ci) &&
        identical(b1$amplitude_ci, b2$amplitude_ci) && identical(b1$seed, 20240601L),
      "the same Run gives the same bootstrap interval (seed 20240601)",
      "the bootstrap interval moved between identical runs")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("percentile CIs (B = 200, participant resampling, seed 20240601)", sm, fixed = TRUE),
      "the summary states the method, B and the seed", "the bootstrap line lacks method or seed")

  # ====================================================== the nested set
  cat("\n-- the nested set says which participants it averages ----------------\n")
  run(session, "two_stage", 7, ms = TRUE)
  ms <- values$harmonic_model$model_selection
  chk(!is.null(ms) && is.finite(attr(ms, "n_common")) && attr(ms, "n_common") <= nrow(raw),
      sprintf("the table carries its common sample (%s of %d)", format(attr(ms, "n_common")), nrow(raw)),
      "the table does not record its common sample")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("scorable under every candidate", sm, fixed = TRUE),
      "the summary says the candidates are compared on the same participants",
      "the summary does not say which participants")

  # ====================================================== the exported script
  # The script's own per-participant fits must be the app's: same origin for the
  # participant missing the first observations, and a held tau held.
  cat("\n-- the exported script reproduces the app's per-participant fits ------\n")
  script_fits <- function() {
    code <- strsplit(paste(generate_analysis_code(), collapse = "\n"), "\n", fixed = TRUE)[[1]]
    a <- grep("^# 10\\. HARMONIC", code)
    b <- grep("average the \\(cos, sin\\) components", code)
    # The section's own call and fit_cosinor run from the script; helpers it
    # does not define resolve in the app (the full script emits them through
    # codetools, which tests/export_roundtrip_test.R covers end to end). What
    # is checked here is the call: the origin and a held tau.
    env <- new.env(parent = environment())
    env$raw_data <- values$data
    eval(parse(text = code[a:b]), envir = env)
    env$cosinor_params
  }
  for (cfg in list(list(id = 8, trend = "log", label = "log trend"),
                   list(id = 9, trend = "exp_sat", label = "tau held at 18 h"))) {
    run(session, "two_stage", cfg$id, trend = cfg$trend, tau = 18, ts_mode = "fixed")
    app <- values$harmonic_model$all_individual_params
    sc <- tryCatch(script_fits(), error = function(e) conditionMessage(e))
    if (is.character(sc)) { fail("the exported cosinor section did not run:", sc); next }
    m <- match(app$subject, sc$subject)
    chk(!anyNA(m) && isTRUE(all.equal(unname(app$mesor), unname(sc$mesor[m]), tolerance = 1e-6)),
        sprintf("%s: the script's intercepts equal the app's, participant 3 included", cfg$label),
        sprintf("%s: the script's intercepts differ from the app's", cfg$label))
  }
  chk(all(values$harmonic_model$all_individual_params$tau == 18),
      "and the app's run held tau at 18 h", "the held-tau run did not hold tau")
})

cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
