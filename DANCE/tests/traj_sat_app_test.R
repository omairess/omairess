# ==============================================================================
# tests/traj_sat_app_test.R -- tau for the mixed cosinor, and A, tau and A/tau
# per cell, driven through the cosinor module's own Run button
#
# The saturating trend's time constant tau can be held at a value or estimated
# once for all cells by profile likelihood; on request the rise A, tau and the
# initial slope A/tau are also estimated in every cell by a nonlinear mixed
# model. This presses Run in each configuration, on four groups simulated with
# DIFFERENT time constants, and checks that
#   * a held tau is the value on screen and is said to be held;
#   * a profiled tau is estimated, carries its profile interval, and is said so;
#   * a flat profile falls back to the value on screen and says why;
#   * the per-cell model fits, its panel and pairwise tables render, tau is
#     compared as a ratio, and one-tau-vs-per-cell is a likelihood-ratio test;
#   * the change's unit line says what the factor in front of A is;
#   * two-stage carries each participant's A/tau into its tests and its panel.
#
# Run with:   Rscript tests/traj_sat_app_test.R      (from the DANCE directory)
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

# ---- data: 08:00 to 09:00 the next day, four groups with their own tau --------
set.seed(7)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 9)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
grp <- factor(rep(c("ADULT", "ELDERLY", "MIDDLE_AGE", "YOUTH"), each = 12))
A0 <- c(ADULT = 20, ELDERLY = 16, MIDDLE_AGE = 18, YOUTH = 22)
T0 <- c(ADULT = 12, ELDERLY = 20, MIDDLE_AGE = 15, YOUTH = 10)
mk <- function(rise = TRUE) t(vapply(seq_along(grp), function(i) {
  g <- as.character(grp[i])
  30 + rnorm(1, 0, 2) + (if (rise) (A0[[g]] + rnorm(1, 0, 2)) * (1 - exp(-tt / T0[[g]])) else 0) +
    6 * cos(2 * pi * (lin - 16) / 24) + 2.5 * cos(4 * pi * (lin - 14) / 24) + rnorm(length(tt), 0, 1.5)
}, numeric(length(tt))))
raw <- mk(TRUE); rownames(raw) <- sprintf("S%02d", seq_along(grp))
flat_raw <- mk(FALSE); rownames(flat_raw) <- rownames(raw)
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

run <- function(session, approach, run_id, tau_mode = "fixed", tau = 18, cells = FALSE,
                lo = 1, hi = 72, ts_mode = "estimate") {
  session$setInputs(harmonic_approach = approach, harmonic_period = 24, n_harmonics = 2,
                    harmonic_trend_type = "exp_sat", harmonic_traj_tau_mode = tau_mode,
                    harmonic_ts_tau_mode = ts_mode, harmonic_tau_lo = lo, harmonic_tau_hi = hi,
                    harmonic_traj_tau = tau, harmonic_traj_sat_cells = cells,
                    harmonic_model_selection = FALSE, harmonic_time_var = "_manual_",
                    harmonic_manual_times = paste(clock, collapse = ","),
                    harmonic_dv_name = "KSS", harmonic_dv_units = "KSS points",
                    harmonic_design = "between", harmonic_between_var = "AGEcategory",
                    harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                    harmonic_traj_level_view = "cells", harmonic_traj_sat_view = "cells",
                    run_harmonic = run_id)
  session$flushReact()
}

testServer(server, {
  # ====================================================== a held tau
  cat("\n-- mixed, tau held at 18.3 h -------------------------------------------\n")
  run(session, "mixed", 1, "fixed", 18.3)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && identical(ff$spec$tau, 18.3) && identical(ff$tau_mode, "fixed") &&
        !isTRUE(ff$tau_estimated),
      "the held tau is the value on screen, and is not marked estimated",
      sprintf("held tau: %s, mode %s", format(ff$spec$tau), ff$tau_mode %||% "?"))
  hd <- render("the tab-6 header renders", text_of(output$harmonic_traj_header))
  chk(grepl("holds τ at 18.3 h, set under Model Specification", hd, fixed = TRUE),
      "the header says tau is held, and at what", "the header does not state the held tau")
  session$setInputs(harmonic_traj_what = "change"); session$flushReact()
  pw <- render("the change pairwise table renders", text_of(output$harmonic_traj_pairwise))
  chk(grepl("1 - e^(-25 h / 18.3 h), the share of A the rise has reached by the last observation",
            pw, fixed = TRUE),
      "the change's unit line says the factor IS 1 - e^(-T/tau), with T and tau",
      "the change's unit line does not explain its factor")
  sh <- render("the saturating panel header renders", text_of(output$harmonic_traj_sat_header))
  chk(grepl("Also", sh, fixed = TRUE) && grepl("estimate A and", sh, fixed = TRUE),
      "without the per-cell model the panel says how to get it", "the panel does not say how to estimate tau per cell")
  session$setInputs(harmonic_traj_sat_view = "profile"); session$flushReact()
  pp <- render("the profile view renders a message for a held tau", output$harmonic_traj_sat_plot)
  chk(grepl("held at a value, not profiled", pp, fixed = TRUE),
      "the profile view says tau was held", "the profile view does not explain a held tau")

  # ====================================================== a profiled tau + per cell
  cat("\n-- mixed, tau profiled for all cells, and A and tau per cell ------------\n")
  t1 <- Sys.time()
  run(session, "mixed", 2, "profile", 18.3, cells = TRUE, lo = 2, hi = 48)
  cat(sprintf("   (fit, profile and per-cell model took %.0f s)\n", as.numeric(Sys.time() - t1, units = "secs")))
  ff <- values$harmonic_model$traj
  pr <- ff$tau_profile
  chk(isTRUE(ff$ok) && identical(ff$tau_mode, "profile") && isTRUE(ff$tau_estimated) &&
        !is.null(pr) && is.finite(ff$spec$tau) && ff$spec$tau >= pr$ci[1] && ff$spec$tau <= pr$ci[2],
      sprintf("tau is estimated by profile likelihood: %.2f h, inside its interval [%.2f, %.2f]",
              ff$spec$tau, pr$ci[1], pr$ci[2]),
      "the profiled tau is missing, unmarked or outside its own interval")
  chk(length(unique(pr$profile$rung)) == 1L,
      "every point of the profile was fitted on one random structure",
      "the profile mixes random structures")
  chk(all(pr$profile$tau >= 2 - 1e-9 & pr$profile$tau <= 48 + 1e-9),
      "the profile searched only the range on screen, 2 to 48 h",
      sprintf("the profile left its range: %.2f to %.2f", min(pr$profile$tau), max(pr$profile$tau)))
  hd <- text_of(output$harmonic_traj_header)
  chk(grepl("estimated for all cells by profile likelihood", hd, fixed = TRUE) &&
        grepl("conditional on that value", hd, fixed = TRUE),
      "the header says tau was profiled and the tests are conditional on it",
      "the header does not state the profiled tau")
  sm <- render("the summary prints", text_of(capture.output(.print_harmonic_summary())))
  chk(grepl("estimated for all cells by profile likelihood", sm, fixed = TRUE),
      "the summary says how tau was obtained", "the summary does not say how tau was obtained")
  session$setInputs(harmonic_traj_sat_view = "profile"); session$flushReact()
  pp <- render("the profile view renders", output$harmonic_traj_sat_plot)
  chk(!is.null(pp) && grepl("Profile log-likelihood", pp, fixed = TRUE),
      "the profile is drawn", "the profile is not drawn")

  sf <- ff$sat
  chk(isTRUE(sf$ok) && nrow(sf$table) == 4L && all(is.finite(sf$table$tau)),
      sprintf("A and tau are estimated in every cell (random structure: %s)", sf$re_label %||% "?"),
      paste("the per-cell model did not fit:", sf$message %||% ""))
  if (isTRUE(sf$ok)) {
    chk(isTRUE(sf$lrt$ok) && sf$lrt$df == 3L && sf$lrt$statistic >= 0,
        sprintf("one tau vs a tau per cell is a likelihood-ratio test: chi2(3) = %.2f, p = %.3f",
                sf$lrt$statistic, sf$lrt$p),
        "the one-tau likelihood-ratio test is missing")
    chk(all(abs(sf$table$init_slope - sf$table$A / sf$table$tau) < 1e-10),
        "A/tau is A divided by tau, cell by cell", "A/tau is not A / tau")
    chk(all(sf$table$tau_lo < sf$table$tau & sf$table$tau < sf$table$tau_hi),
        "each tau sits inside its interval", "a tau lies outside its interval")
    chk(all(sf$table$tau_lo >= 2 - 1e-9 & sf$table$tau_hi <= 48 + 1e-9),
        "every per-cell tau and its interval stay inside the range 2 to 48 h",
        "a per-cell tau or its interval left the range")
    sh <- text_of(output$harmonic_traj_sat_header)
    chk(grepl("Is one", sh, fixed = TRUE) && grepl("Likelihood-ratio test", sh, fixed = TRUE),
        "the panel header reports the one-tau test", "the panel header lacks the one-tau test")
    for (v in c("cells", "pairs", "curves")) {
      session$setInputs(harmonic_traj_sat_view = v); session$flushReact()
      pv <- render(sprintf("the saturating panel renders the '%s' view", v), output$harmonic_traj_sat_plot)
      if (identical(v, "pairs"))
        chk(!is.null(pv) && grepl('"type":"log"', pv, fixed = TRUE),
            "tau's pairwise panel is a ratio on a log axis", "tau's pairwise panel is not a ratio")
    }
    st <- render("the saturating tests render", text_of(output$harmonic_traj_sat_tests))
    chk(grepl("all cells equal", st, fixed = TRUE) && grepl("likelihood ratio", st, fixed = TRUE) &&
          grepl("delta method", st, fixed = TRUE),
        "the tests table has the Wald tests per parameter and the likelihood-ratio test",
        "the tests table is incomplete")
    ws <- harmonic_traj_whats(ff)
    chk(all(c("A_sat", "tau", "init_slope") %in% ws) &&
          match("A_sat", ws) == match("change", ws) + 1L,
        "the pairwise menu offers A, tau and A/tau right after the level quantities",
        paste("pairwise menu:", paste(ws, collapse = ", ")))
    session$setInputs(harmonic_traj_what = "tau"); session$flushReact()
    pt <- render("the tau pairwise table renders", text_of(output$harmonic_traj_pairwise))
    chk(grepl("Ratio", pt, fixed = TRUE) && grepl("tested on the log scale", pt, fixed = TRUE),
        "tau is compared as a ratio, tested on the log scale", "the tau pairwise table is not a ratio")
    session$setInputs(harmonic_traj_what = "init_slope"); session$flushReact()
    pi_ <- render("the A/tau pairwise table renders", text_of(output$harmonic_traj_pairwise))
    chk(grepl("KSS points per hour (delta method)", pi_, fixed = TRUE),
        "A/tau is compared in KSS points per hour", "the A/tau pairwise table has no units")
    sm <- text_of(capture.output(.print_harmonic_summary()))
    chk(grepl("Saturating trend per cell (nonlinear mixed model)", sm, fixed = TRUE) &&
          grepl("One tau for all cells? Likelihood-ratio", sm, fixed = TRUE),
        "the summary prints the per-cell block and the one-tau test",
        "the summary lacks the per-cell block")
    rp <- render("the APA report renders", paste(dance_apa_report(values, input, "T"), collapse = "\n"))
    chk(grepl("**Saturating trend per cell.**", rp, fixed = TRUE) &&
          grepl("estimated for all cells by profile likelihood", rp, fixed = TRUE) &&
          grepl("likelihood-ratio test", rp, fixed = TRUE),
        "the report carries the per-cell table, the profiled tau and the one-tau test",
        "the report lacks the per-cell table or the tau statement")
  }

  # ====================================================== a range that excludes the optimum
  cat("\n-- mixed, tau searched only between 30 and 60 h ------------------------\n")
  # the held value -- what a flat profile would fall back to -- has to lie in
  # the range too, or Run refuses (checked below)
  run(session, "mixed", 5, "profile", 40, lo = 30, hi = 60)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && ff$spec$tau >= 30 - 1e-9 && ff$spec$tau <= 60 + 1e-9,
      sprintf("tau stays inside a range that excludes the optimum (%.1f h)", ff$spec$tau),
      sprintf("tau left its range: %s", format(ff$spec$tau)))
  hd <- text_of(output$harmonic_traj_header)
  chk(grepl("END of the search range", hd, fixed = TRUE),
      "and an estimate at the end of the range is called a bound",
      "an estimate at the end of the range is not flagged")

  # ====================================================== settings that cannot run
  cat("\n-- a range that is not a range is refused, not fitted -----------------\n")
  before <- values$harmonic_model
  run(session, "mixed", 6, "profile", 18.3, lo = 40, hi = 10)
  chk(identical(values$harmonic_model, before),
      "a range whose lower end is above its upper end is refused",
      "an impossible range was fitted anyway")
  run(session, "mixed", 9, "profile", 18.3, lo = 30, hi = 60)
  chk(identical(values$harmonic_model, before),
      "a held value outside the range is refused: a flat profile would hold tau there",
      "a held value outside the range was accepted")

  # ====================================================== a flat profile
  cat("\n-- mixed, tau profiled on data with no saturating rise -----------------\n")
  values$data <- flat_raw
  run(session, "mixed", 3, "profile", 18.3)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && identical(ff$tau_mode, "fixed_after_flat") && identical(ff$spec$tau, 18.3),
      "a flat profile falls back to the value on screen, 18.3 h",
      sprintf("flat profile: mode %s, tau %s", ff$tau_mode %||% "?", format(ff$spec$tau)))
  hd <- text_of(output$harmonic_traj_header)
  chk(grepl("not identified by these data", hd, fixed = TRUE),
      "and the header says tau is not identified", "the flat profile is not reported")
  values$data <- raw

  # ====================================================== two-stage, tau held
  cat("\n-- two-stage, tau held at 18 h for every participant -------------------\n")
  run(session, "two_stage", 7, ts_mode = "fixed", tau = 18)
  P <- values$harmonic_model$individual_params
  chk(!is.null(P) && nrow(P) > 40 && all(P$tau == 18) && all(is.na(P$tau_se)),
      sprintf("every participant's fit holds tau at 18 h (%d fits)", nrow(P %||% data.frame())),
      "a held tau was not held in every fit")
  chk(isTRUE(all.equal(P$init_slope, P$A_sat / 18)),
      "and each A/tau is that participant's A over 18", "A/tau does not use the held tau")
  mod <- values$harmonic_model
  chk(!"tau" %in% dance_ts_coef_cols(mod) && !"tau" %in% dance_ts_pair_whats(mod),
      "a held tau is not tested: out of the MANOVA and the pairwise menu",
      "a constant tau is still offered as a test")
  cm <- text_of(output$harmonic_traj_components)
  chk(!grepl("Time constant tau", cm, fixed = TRUE),
      "the component table has no tau row when tau is held",
      "the component table tests a constant tau")
  sh <- text_of(output$harmonic_traj_sat_header)
  chk(grepl("held", sh, fixed = TRUE) && grepl("18 h", sh, fixed = TRUE),
      "the panel says tau was held, and at what", "the panel does not say tau was held")
  render("the held-tau panel renders", output$harmonic_traj_sat_plot)
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("tau held at 18.0 h for every participant", sm, fixed = TRUE) &&
        grepl("tau held at 18.0 h", sm, fixed = TRUE),
      "the summary says tau was held for every participant", "the summary does not state the held tau")

  # ====================================================== two-stage, a range
  cat("\n-- two-stage, each participant's tau within 5 to 30 h -----------------\n")
  run(session, "two_stage", 8, ts_mode = "estimate", lo = 5, hi = 30)
  P <- values$harmonic_model$individual_params
  chk(all(P$tau >= 5 - 1e-6 & P$tau <= 30 + 1e-6),
      sprintf("every participant's tau is inside 5 to 30 h (%.2f to %.2f)", min(P$tau), max(P$tau)),
      sprintf("a participant's tau left the range: %.2f to %.2f", min(P$tau), max(P$tau)))
  sh <- text_of(output$harmonic_traj_sat_header)
  chk(grepl("within 5 to 30 h", sh, fixed = TRUE),
      "the panel states the range", "the panel does not state the range")

  # ====================================================== two-stage
  cat("\n-- two-stage: each participant's A/tau ----------------------------------\n")
  run(session, "two_stage", 4)
  P <- values$harmonic_model$individual_params
  chk(!is.null(P) && "init_slope" %in% names(P) &&
        isTRUE(all.equal(P$init_slope, P$A_sat / P$tau)),
      "each participant carries their own A/tau", "the per-participant A/tau is missing")
  session$setInputs(harmonic_traj_sat_view = "cells"); session$flushReact()
  render("the two-stage saturating panel renders", output$harmonic_traj_sat_plot)
  st <- render("the two-stage saturating tests render", text_of(output$harmonic_traj_sat_tests))
  chk(grepl("Initial slope A/tau (KSS points per hour)", st, fixed = TRUE),
      "A/tau is tested between groups", "A/tau is not tested between groups")
  session$setInputs(harmonic_traj_what = "init_slope"); session$flushReact()
  pt <- render("the two-stage A/tau pairwise table renders", text_of(output$harmonic_traj_pairwise))
  chk(grepl("in KSS points per hour", pt, fixed = TRUE),
      "the two-stage A/tau pairwise is in KSS points per hour", "the two-stage A/tau pairwise has no units")
  cm <- text_of(output$harmonic_traj_components)
  chk(grepl("Initial slope A/tau", cm, fixed = TRUE),
      "the component table tests A/tau", "the component table lacks A/tau")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("Saturating trend per group (means of the per-participant fits)", sm, fixed = TRUE),
      "the two-stage summary prints A, tau and A/tau per group",
      "the two-stage summary lacks the per-group saturating block")
})

cat(sprintf("\n%d failure(s)\n", failures))
if (failures > 0) quit(status = 1)
