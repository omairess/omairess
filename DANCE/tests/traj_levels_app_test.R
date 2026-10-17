# ==============================================================================
# tests/traj_levels_app_test.R -- intercept, MESOR and the non-periodic change,
# driven through the cosinor module's own Run button
#
# Tab 6 shows the level of each trajectory as three numbers: the intercept (the
# midline at the first observation), the MESOR (that midline averaged over the
# recording) and, with a trend, the non-periodic change from the first to the
# last observation. This fits both approaches on a saturating-trend design with
# four groups and checks that
#   * the plot renders, per cell and as pairwise differences;
#   * its numbers are the pairwise table's numbers and the summary's numbers;
#   * the two-stage group MESOR is the mean of the participants' own MESORs,
#     integrated on the model axis (it was integrated on the clock axis);
#   * the mixed model holds tau at the value on screen, not a hidden default;
#   * without a trend there is one quantity, the intercept, which is the MESOR.
#
# Run with:   Rscript tests/traj_levels_app_test.R      (from the DANCE directory)
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
render <- function(label, expr) {
  v <- tryCatch(force(expr), error = function(e) structure(conditionMessage(e), class = "dance_err"))
  if (inherits(v, "dance_err")) { fail(label, "->", as.character(v)); return(NULL) }
  ok(label); v
}
text_of <- function(x) paste(as.character(unlist(x)), collapse = " ")
near <- function(a, b, tol = 1e-6) isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))

# ---- data: 08:00 to 09:00 the next day, four groups, a saturating rise ---------
set.seed(2026)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 9)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]                                  # model axis, 0 to 25 h
grp <- factor(rep(c("ADULT", "ELDERLY", "MIDDLE_AGE", "YOUTH"), each = 12))
b0  <- c(ADULT = 30, ELDERLY = 26, MIDDLE_AGE = 28, YOUTH = 34)
A   <- c(ADULT = 20, ELDERLY = 16, MIDDLE_AGE = 18, YOUTH = 22)
raw <- t(vapply(seq_along(grp), function(i) {
  g <- as.character(grp[i]); tau_i <- max(6, rnorm(1, 15, 4))
  b0[[g]] + rnorm(1, 0, 2) + A[[g]] * (1 - exp(-tt / tau_i)) +
    6 * cos(2 * pi * (lin - 16) / 24) + 2.5 * cos(4 * pi * (lin - 14) / 24) + rnorm(length(tt), 0, 1.5)
}, numeric(length(tt))))
rownames(raw) <- sprintf("S%02d", seq_along(grp))
covs <- data.frame(AGEcategory = grp, stringsAsFactors = FALSE)

server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data <- raw
  values$time_labels <- sprintf("%02d:00", clock)
  values$covariates <- covs
  values$group_labels <- grp
  values$group_variables <- list(AGEcategory = grp)
  values$subject_ids <- rownames(raw)
}

run <- function(session, approach, trend, run_id, tau = 18) {
  # the tau range's on-screen defaults: a mock session has no UI to supply them,
  # and a two-stage saturating fit estimates each participant's tau inside it
  session$setInputs(harmonic_approach = approach, harmonic_period = 24, n_harmonics = 2,
                    harmonic_trend_type = trend, harmonic_traj_tau = tau,
                    harmonic_tau_lo = 1, harmonic_tau_hi = 72,
                    harmonic_model_selection = FALSE, harmonic_time_var = "_manual_",
                    harmonic_manual_times = paste(clock, collapse = ","),
                    harmonic_dv_name = "KSS", harmonic_dv_units = "KSS points",
                    harmonic_design = "between", harmonic_between_var = "AGEcategory",
                    harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                    harmonic_traj_level_view = "cells", run_harmonic = run_id)
  session$flushReact()
}

testServer(server, {
  # ======================================================== MIXED, saturating
  cat("\n-- mixed-effects, saturating trend, tau on screen = 12 h ---------------\n")
  run(session, "mixed", "exp_sat", 1, tau = 12)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok), "the mixed model fitted", paste("the mixed model did not fit:", ff$message))
  chk(identical(ff$spec$tau, 12), "tau is the value on screen (12 h), not a hidden default",
      sprintf("tau is %s, not the 12 h on screen", format(ff$spec$tau)))
  sm <- render("the summary prints", text_of(capture.output(.print_harmonic_summary())))
  chk(grepl("tau held at 12.0 h", sm, fixed = TRUE) || grepl("held at 12.0 h", sm, fixed = TRUE),
      "the summary states the tau the fit held", "the summary does not state the held tau")

  run(session, "mixed", "exp_sat", 2, tau = 18)
  ff <- values$harmonic_model$traj
  lv <- harmonic_traj_levels()
  chk(isTRUE(lv$ok) && identical(lv$quantities, c("intercept", "mesor", "change")),
      "three level quantities with a trend: intercept, MESOR, change",
      paste("level quantities:", paste(lv$quantities, collapse = ", ")))
  # the identities, cell by cell
  tw <- lv$window
  ce <- dance_traj_cell_equations(ff)
  est <- function(q, cl) lv$table$estimate[lv$table$quantity == q & lv$table$cell == cl]
  idok <- all(vapply(seq_len(nrow(ce$table)), function(i) {
    cl <- ce$table$cell[i]; Ai <- ce$trend_coefs[[i]][[1]]
    near(est("intercept", cl), ce$table$intercept[i]) &&
      near(est("mesor", cl), ce$table$intercept[i] + Ai * tw$mean[[1]]) &&
      near(est("change", cl), Ai * tw$end[[1]])
  }, logical(1)))
  chk(idok, "MESOR = intercept + A * window mean, change = A * (1 - e^(-T/tau)), in every cell",
      "a cell's MESOR or change does not follow from its intercept and A")
  chk(near(tw$end[[1]], 1 - exp(-25 / 18)) && near(tw$mean[[1]], 1 - (18 / 25) * (1 - exp(-25 / 18))),
      "the window is 0 to 25 h on the model axis, tau 18 h",
      sprintf("the window summaries are off: end %.6f, mean %.6f", tw$end[[1]], tw$mean[[1]]))
  p1 <- render("the level plot renders per cell", output$harmonic_traj_levels_plot)
  chk(!is.null(p1) && grepl("error_y", p1, fixed = TRUE) && grepl("MESOR (08:00 to 09:00 (+1d))", p1, fixed = TRUE),
      "per-cell panels carry intervals and name the window on the clock",
      "the per-cell plot lacks intervals or the clock window")
  nt <- render("the level note renders", text_of(output$harmonic_traj_levels_note))
  chk(grepl("held at 18 h", nt, fixed = TRUE) && grepl("testing that coefficient", nt, fixed = TRUE),
      "the note says the change is the trend coefficient's test, tau held at 18 h",
      "the note does not explain the change statistic")
  session$setInputs(harmonic_traj_level_view = "pairs"); session$flushReact()
  p2 <- render("the level plot renders as pairwise differences", output$harmonic_traj_levels_plot)
  chk(!is.null(p2) && grepl("error_x", p2, fixed = TRUE), "pairwise panels carry intervals",
      "the pairwise plot lacks intervals")
  ct <- dance_traj_contrasts(ff, "mesor", adjust = "holm", df_method = "satterthwaite")
  cached <- harmonic_traj_ct("mesor", "holm", "satterthwaite")
  chk(isTRUE(ct$ok) && near(cached$table$estimate, ct$table$estimate) && near(cached$table$p_adj, ct$table$p_adj),
      sprintf("the plotted MESOR differences are the pairwise table's (%d pairs)", nrow(ct$table)),
      "the plotted MESOR differences differ from the pairwise table")
  chk(near(ct$table$estimate[1], est("mesor", ct$table$cell1[1]) - est("mesor", ct$table$cell2[1])),
      "a pairwise MESOR difference is the difference of the two cell MESORs",
      "a pairwise MESOR difference does not match the cell estimates")
  chg <- dance_traj_contrasts(ff, "change", adjust = "none", df_method = "satterthwaite")
  trd <- dance_traj_contrasts(ff, "trend", adjust = "none", df_method = "satterthwaite")
  chk(isTRUE(chg$ok) && isTRUE(trd$ok) && near(chg$table$p_raw, trd$table$p_raw, 1e-8) &&
        near(chg$table$estimate, trd$table$estimate * tw$end[[1]]),
      "the change test is the trend-coefficient test, its estimate A x (1 - e^(-T/tau))",
      "the change and the trend coefficient give different tests")
  session$setInputs(harmonic_traj_what = "change"); session$flushReact()
  pw <- render("the pairwise table renders for the change", text_of(output$harmonic_traj_pairwise))
  chk(grepl("KSS points", pw, fixed = TRUE) && grepl("so its test is the test of that coefficient", pw, fixed = TRUE),
      "the change's unit line is in KSS points and names its coefficient",
      "the change's unit line is missing its units or its coefficient")
  session$setInputs(harmonic_traj_component = "trend"); session$flushReact()
  bn <- render("the change-from-origin note renders", text_of(output$harmonic_traj_band_note))
  chk(grepl("Each curve ends at its cell's non-periodic change", bn, fixed = TRUE),
      "the change-from-origin view names the statistic it is a picture of",
      "the change-from-origin view does not point to its statistic")
  session$setInputs(harmonic_traj_component = "full"); session$flushReact()
  cm <- render("the component table renders", text_of(output$harmonic_traj_components))
  chk(grepl("Intercept (midline at the first observation)", cm, fixed = TRUE) &&
        grepl("MESOR (midline averaged over the window)", cm, fixed = TRUE) &&
        grepl("saturating rise A (KSS points above the intercept", cm, fixed = TRUE) &&
        !grepl("Baseline", cm, fixed = TRUE),
      "the component table tests the intercept and the MESOR, and names the trend with its units",
      "the component table still uses the old level rows or an unlabelled trend")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("Intercept (midline at 08:00)", sm, fixed = TRUE) &&
        grepl("MESOR (midline averaged over 08:00 to 09:00 (+1d))", sm, fixed = TRUE) &&
        grepl("Non-periodic change, 08:00 to 09:00 (+1d)", sm, fixed = TRUE) &&
        !grepl("Level at the first observation", sm, fixed = TRUE),
      "the summary prints intercept, MESOR and change per cell",
      "the summary still prints the old level line or lacks the new ones")
  r1 <- ce$table[1, ]
  chk(grepl(sprintf("MESOR (midline averaged over 08:00 to 09:00 (+1d)): %.3f", r1$mesor), sm, fixed = TRUE),
      "the summary's cell MESOR is the plotted one", "the summary's MESOR differs from the plot")

  # ======================================================== TWO-STAGE, saturating
  cat("\n-- two-stage, saturating trend -------------------------------------------\n")
  run(session, "two_stage", "exp_sat", 3)
  mod <- values$harmonic_model
  P <- mod$individual_params
  chk(!is.null(mod) && "trend_change" %in% names(P), "each participant carries their change over the recording",
      "the participant table has no trend_change")
  manual <- with(P, A_sat * ((1 - exp(-25 / tau)) - 0))
  chk(near(P$trend_change, manual, 1e-8), "trend_change = A_sat * (1 - e^(-25/tau)) per participant, from both parameters",
      "trend_change does not follow from the participant's A_sat and tau")
  g1 <- names(mod$group_fits)[1]
  gp <- P[as.character(grp[P$subject]) == g1, ]
  gf <- mod$group_fits[[g1]]
  chk(near(gf$rhythm_adjusted_mean, mean(gp$mesor_adj)),
      "the group MESOR is the mean of the members' own MESORs",
      sprintf("group MESOR %.4f is not the mean of the members' MESORs %.4f", gf$rhythm_adjusted_mean, mean(gp$mesor_adj)))
  ref <- dance_rhythm_adjusted_mean(gf$intercept, "exp_sat", gf$trend_coefs, 0, 25, 0)
  chk(near(gf$rhythm_adjusted_mean_of_mean_coefs, ref),
      "the MESOR of the mean parameters is integrated over the model window, 0 to 25 h",
      sprintf("the mean-parameter MESOR %.4f is not the 0-25 h integral %.4f",
              gf$rhythm_adjusted_mean_of_mean_coefs, ref))
  chk(near(gf$trend_change, mean(gp$trend_change)), "the group change is the mean of the members' changes",
      "the group change is not the members' mean")
  session$setInputs(harmonic_traj_level_view = "cells"); session$flushReact()
  lv <- harmonic_traj_levels()
  chk(isTRUE(lv$ok) && near(lv$table$estimate[lv$table$quantity == "mesor" & lv$table$cell == g1], gf$rhythm_adjusted_mean),
      "the plotted group MESOR is the summary's", "the plotted group MESOR is not the summary's")
  render("the two-stage level plot renders per group", output$harmonic_traj_levels_plot)
  session$setInputs(harmonic_traj_level_view = "pairs"); session$flushReact()
  p3 <- render("the two-stage level plot renders as differences", output$harmonic_traj_levels_plot)
  ts_pw <- dance_ts_pairwise(harmonic_ts()$df, mod, "change", adjust = "holm")
  d12 <- mean(P$trend_change[as.character(grp[P$subject]) == ts_pw$table$cell1[1]]) -
         mean(P$trend_change[as.character(grp[P$subject]) == ts_pw$table$cell2[1]])
  chk(isTRUE(ts_pw$ok) && near(ts_pw$table$estimate[1], d12),
      "the two-stage change difference is the difference of group means (Welch)",
      "the two-stage change difference is off")
  nt <- text_of(output$harmonic_traj_levels_note)
  chk(grepl("own rise A and time constant tau", nt, fixed = TRUE),
      "the two-stage note says the change uses both saturating parameters",
      "the two-stage note does not explain the change")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl(sprintf("MESOR (midline averaged over 08:00 to 09:00 (+1d)): %.3f", gf$rhythm_adjusted_mean), sm, fixed = TRUE),
      "the two-stage summary prints the group MESOR the tests compare",
      "the two-stage summary prints a different group MESOR")
  session$setInputs(harmonic_traj_what = "intercept"); session$flushReact()
  pw <- text_of(output$harmonic_traj_pairwise)
  chk(grepl("in KSS points", pw, fixed = TRUE), "the two-stage intercept table is in KSS points",
      "the two-stage intercept table has no units")

  # ======================================================== no trend
  cat("\n-- no trend: the intercept is the MESOR ------------------------------------\n")
  run(session, "mixed", "none", 4)
  lv <- harmonic_traj_levels()
  chk(isTRUE(lv$ok) && identical(lv$quantities, "intercept"), "mixed, no trend: one quantity",
      paste("mixed, no trend:", paste(lv$quantities, collapse = ", ")))
  session$setInputs(harmonic_traj_level_view = "cells"); session$flushReact()
  p4 <- render("mixed, no trend: the level plot renders", output$harmonic_traj_levels_plot)
  chk(!is.null(p4) && grepl("Intercept = MESOR (no trend)", p4, fixed = TRUE),
      "the single panel is named Intercept = MESOR", "the no-trend panel is misnamed")
  chk(identical(dance_traj_pair_whats(values$harmonic_model$traj)[1], "intercept") &&
        !any(c("mesor", "change", "level") %in% dance_traj_pair_whats(values$harmonic_model$traj)),
      "no MESOR, change or curve-value entry is offered without a trend",
      "the no-trend pairwise menu offers a duplicate level quantity")
  run(session, "two_stage", "none", 5)
  lv <- harmonic_traj_levels()
  chk(isTRUE(lv$ok) && identical(lv$quantities, "intercept"), "two-stage, no trend: one quantity",
      paste("two-stage, no trend:", paste(lv$quantities, collapse = ", ")))
  render("two-stage, no trend: the level plot renders", output$harmonic_traj_levels_plot)
})

cat(sprintf("\n%d failure(s)\n", failures))
if (failures > 0) quit(status = 1)
