# ==============================================================================
# tests/traj_car1_app_test.R -- autocorrelated residuals in the mixed cosinor,
# through the module's Run button
#
# The residuals within a curve can be modelled as a continuous-time AR(1)
# (Model Specification > Residuals within a curve). This simulates curves whose
# residuals ARE autocorrelated (phi = 0.6 per hour, on an uneven 08:00 to
# 15:00 grid with missing observations) and checks that
#   * an independent-residual run says, in its diagnostics, that consecutive
#     residuals are correlated and points to the AR(1) model;
#   * an AR(1) run estimates phi with an interval that covers the truth, tests
#     it against independence, and says so in the summary, the comparison tab,
#     the assumptions list, the residual diagnostics and the report;
#   * the comparison tab's tests run on it with Kenward-Roger df;
#   * with tau profiled, tau and phi are estimated together and settle, and
#     independence is tested against the independent model at its own best tau;
#   * the per-cell nonlinear model carries the same residual model;
#   * the setting is restored with a session and written into the export.
#
# Run with:   Rscript tests/traj_car1_app_test.R      (from the DANCE directory)
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
# whitespace squashed: the summary wraps its long lines
text_of <- function(x) gsub("\\s+", " ", paste(as.character(unlist(x)), collapse = " "))

# ---- data: 08:00 to 15:00 the next day, two groups, AR(1) residuals ----------
set.seed(31)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 15)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
PHI <- 0.6
grp <- factor(rep(c("A", "B"), each = 14))
R <- PHI^abs(outer(tt, tt, "-")); Lc <- t(chol(R))
raw <- t(vapply(seq_along(grp), function(i) {
  b <- if (grp[i] == "B") 1 else 0
  20 + b + rnorm(1, 0, 1.5) + 10 * (1 - exp(-tt / 12)) +
    (4 + rnorm(1, 0, 0.8)) * cos(2 * pi * (lin - 16) / 24) + as.numeric(Lc %*% rnorm(length(tt), 0, 1.2))
}, numeric(length(tt))))
rownames(raw) <- sprintf("S%02d", seq_along(grp))
raw[4, c(2, 7)] <- NA; raw[17, 11] <- NA
covs <- data.frame(Group = grp, stringsAsFactors = FALSE)

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

run <- function(session, run_id, resid = "independent", trend = "linear", tau_mode = "fixed",
                cells = FALSE) {
  session$setInputs(harmonic_approach = "mixed", harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = trend, harmonic_traj_tau_mode = tau_mode,
                    harmonic_traj_tau = 12, harmonic_tau_lo = 2, harmonic_tau_hi = 48,
                    harmonic_traj_sat_cells = cells, harmonic_traj_resid = resid,
                    harmonic_time_var = "_shared_", harmonic_dv_name = "KSS",
                    harmonic_dv_units = "KSS points", harmonic_design = "between",
                    harmonic_between_var = "Group", harmonic_df_method = "kr",
                    harmonic_traj_adjust = "holm", harmonic_traj_what = "intercept",
                    run_harmonic = run_id)
  session$flushReact()
}

testServer(server, {
  # ====================================================== independent residuals
  cat("\n-- independent residuals: the diagnostics see the autocorrelation -----\n")
  run(session, 1)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && identical(ff$phi, 0) && is.null(ff$residual_ar),
      "the default fit has independent residuals", "the default fit is not the independent one")
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("Residuals: independent (no autocorrelation modelled)", sm, fixed = TRUE) &&
        grepl("Conditional residuals of a curve: consecutive ones correlate r =", sm, fixed = TRUE) &&
        grepl("The AR(1) residual model (Model Specification) estimates the autocorrelation", sm, fixed = TRUE),
      "the summary reports how consecutive residuals correlate and points to the AR(1) model",
      "the summary does not report the residual autocorrelation")
  lg <- dance_traj_resid_lag1(ff, "raw")
  ok(sprintf("  (consecutive conditional residuals correlate r = %.2f)", lg$r))

  # ====================================================== AR(1)
  cat("\n-- continuous-time AR(1) residuals, the generating trend (tau 12 h) ----\n")
  # with the trend the data were made with: a misfitted trend leaves smooth
  # residual structure that the AR(1) would absorb, inflating phi
  run(session, 2, resid = "car1", trend = "exp_sat")
  fa <- values$harmonic_model$traj
  ar <- fa$residual_ar
  chk(isTRUE(fa$ok) && isTRUE(fa$phi > 0) && isTRUE(ar$ok) && identical(ar$method, "REML"),
      sprintf("phi estimated by REML: %.3f per hour [%.3f, %.3f]", ar$phi, ar$ci[1], ar$ci[2]),
      sprintf("no AR(1) fit: %s", fa$message %||% ""))
  chk(ar$ci[1] <= PHI && ar$ci[2] >= PHI,
      sprintf("its interval covers the true phi of %.1f", PHI), "the interval misses the true phi")
  chk(isTRUE(ar$lrt$p < 0.001) && identical(fa$phi, ar$phi),
      sprintf("independence is rejected (chi2 = %.1f, p = %.2g) and the fit is at phi-hat",
              ar$lrt$statistic, ar$lrt$p), "the test against independence is missing or weak")
  run(session, 21, resid = "independent", trend = "exp_sat"); ff <- values$harmonic_model$traj
  lg <- dance_traj_resid_lag1(ff, "raw")
  run(session, 2, resid = "car1", trend = "exp_sat"); fa <- values$harmonic_model$traj
  ll_i <- dance_traj_loglik(ff); ll_a <- dance_traj_loglik(fa)
  chk(isTRUE(ll_a$df == ll_i$df + 1L) && ll_a$logLik > ll_i$logLik,
      sprintf("phi is counted as a parameter, and the AR(1) fits better (logLik %.1f vs %.1f)",
              ll_a$logLik, ll_i$logLik), "the AR(1) log-likelihood is not on the data's scale")
  lgn <- dance_traj_resid_lag1(fa, "normalized")
  # near zero, not exactly zero: residuals of a random-slope model carry a small
  # negative correlation of their own (fitting the random effects removes
  # curve-level structure), whitened or not
  chk(abs(lgn$r) < 0.15,
      sprintf("its normalized residuals are nearly uncorrelated (r = %.2f)", lgn$r),
      sprintf("the normalized residuals still correlate (r = %.2f)", lgn$r))
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("continuous-time AR(1) within each curve, phi =", sm, fixed = TRUE) &&
        grepl("corr(e(t), e(t')) = phi^|t - t'|", sm, fixed = TRUE) &&
        grepl("phi included", sm, fixed = TRUE) &&
        grepl("Normal, with ONE variance for every cell and time point, and AUTOCORRELATED", sm, fixed = TRUE),
      "the summary states the residual model, its equation, its parameter and its assumptions",
      "the summary does not state the AR(1) model")
  hd <- render("the comparison tab's header renders", text_of(output$harmonic_traj_header))
  chk(grepl("Residuals are autocorrelated within each curve as a continuous-time AR(1)", hd, fixed = TRUE) &&
        grepl("phi as known", hd, fixed = TRUE),
      "the header says the tests are conditional on the AR(1)", "the header does not state the residual model")
  render("the primary tests render", output$harmonic_traj_primary)
  tests <- dance_traj_tests_table(harmonic_traj_tests())
  chk(!is.null(tests) && all(grepl("Kenward-Roger", tests$method)) && all(tests$df2 > 15 & tests$df2 < 30),
      sprintf("every design test is Kenward-Roger, df %.1f to %.1f (28 participants)",
              min(tests$df2), max(tests$df2)), "the tests did not run on the AR(1) fit")
  render("the pairwise table renders", output$harmonic_traj_pairwise)
  render("the residual plot renders", output$harmonic_residual_plot)
  rs <- harmonic_residuals()
  chk(grepl("normalized residuals", rs$source, fixed = TRUE) && abs(mean(rs$fitted) - mean(raw, na.rm = TRUE)) < 2,
      "the residual plot shows normalized residuals against fitted values on the data's scale",
      "the residual plot is on the whitened scale")
  r2 <- dance_traj_r2(fa)
  chk(!is.null(r2) && r2$marginal > 0 && r2$conditional > r2$marginal && r2$conditional < 1,
      sprintf("R-squared on the data's scale (marginal %.2f, conditional %.2f)", r2$marginal, r2$conditional),
      "R-squared is off for the AR(1) fit")
  rp <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
  chk(grepl("continuous-time first-order autoregressive process", rp, fixed = TRUE) &&
        grepl("The residuals were autocorrelated within curves", rp, fixed = TRUE),
      "the report's Methods and Results describe the AR(1)", "the report omits the AR(1)")
  code <- paste(generate_analysis_code(), collapse = "\n")
  chk(grepl('residuals = "car1")', code, fixed = TRUE),
      "the exported script fits the same residual model", "the export drops the residual model")

  # ====================================================== tau profiled, jointly with phi
  cat("\n-- saturating trend, tau profiled together with phi --------------------\n")
  run(session, 3, resid = "car1", trend = "exp_sat", tau_mode = "profile")
  fj <- values$harmonic_model$traj
  aj <- fj$residual_ar
  chk(isTRUE(fj$ok) && identical(fj$tau_mode, "profile") && grepl("jointly with tau", aj$method) &&
        isTRUE(aj$settled),
      sprintf("tau %.1f h and phi %.3f estimated together by ML, settled in %d round(s)",
              fj$spec$tau, aj$phi, aj$rounds),
      sprintf("joint estimation failed: %s", fj$message %||% aj$method %||% ""))
  hd <- text_of(output$harmonic_traj_header)
  chk(grepl("estimated together, by maximum likelihood", hd, fixed = TRUE),
      "the tau statement says tau and phi were estimated together", "the tau statement is silent on phi")
  chk(grepl("The interval is computed at tau's estimate; the likelihood ratio compares with independent residuals at their own best tau.",
            text_of(dance_traj_residual_statement(fj)), fixed = TRUE),
      "the phi statement says its interval is at tau-hat and how independence was fitted",
      "the phi statement does not say what the joint interval and test are conditional on")
  rpj <- text_of(dance_apa_report(values, input, "T"))
  smj <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("(at the estimated time constant); against independent residuals (with the time constant re-estimated)",
            rpj, fixed = TRUE) &&
        grepl("With tau estimated too, phi's interval is computed at tau's estimate", smj, fixed = TRUE),
      "the report and the assumptions list say what the joint interval and test are conditional on",
      "the report or the assumptions list is silent on the joint interval and test")
  # independence, with tau estimated, is the independent model at ITS best tau:
  # the maximum of the tau profile an independent-residual run computes
  run(session, 5, resid = "independent", trend = "exp_sat", tau_mode = "profile")
  fi <- values$harmonic_model$traj
  at_tau_hat <- aj$profile$logLik[aj$profile$phi == 0][1]
  chk(identical(as.integer(fi$tau_profile$rung), as.integer(aj$rung)) &&
        isTRUE(all.equal(aj$lrt$ll0, max(max(fi$tau_profile$profile$logLik), at_tau_hat))) &&
        aj$lrt$ll0 >= at_tau_hat &&
        isTRUE(all.equal(aj$lrt$statistic, 2 * (aj$lrt$ll_max - aj$lrt$ll0))) &&
        grepl("tau re-estimated", aj$lrt$method, fixed = TRUE),
      sprintf(paste("the likelihood ratio for phi compares with independence at its own best tau",
                    "(log-likelihood %.3f, against %.3f at the AR(1)'s tau): chi2(1) = %.2f"),
              aj$lrt$ll0, at_tau_hat, aj$lrt$statistic),
      sprintf("the joint likelihood ratio's null is not the independent model's maximum (%s vs %s)",
              format(aj$lrt$ll0), format(max(fi$tau_profile$profile$logLik))))

  # ====================================================== per cell
  cat("\n-- saturating trend, tau held, A and tau per cell ----------------------\n")
  run(session, 4, resid = "car1", trend = "exp_sat", tau_mode = "fixed", cells = TRUE)
  fc <- values$harmonic_model$traj
  chk(isTRUE(fc$sat$ok) && isTRUE(fc$sat$phi > 0) &&
        inherits(fc$sat$fit$modelStruct$corStruct, "corCAR1"),
      sprintf("the per-cell nonlinear model carries CAR(1) residuals (phi %.3f)", fc$sat$phi %||% NA),
      sprintf("the per-cell model has no AR(1): %s", fc$sat$message %||% ""))
  sm <- text_of(capture.output(.print_harmonic_summary()))
  chk(grepl("The per-cell nonlinear model has the same CAR(1) residuals", sm, fixed = TRUE),
      "the assumptions list says so", "the assumptions list omits the per-cell residual model")

  # ====================================================== restore
  chk("harmonic_traj_resid" %in% RESTORE_INPUTS$select,
      "the residual setting is restored with a session", "the residual setting is not restored")
})

cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
