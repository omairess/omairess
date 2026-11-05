# ==============================================================================
# tests/single_traj_app_test.R -- one trajectory, no between or within factor
#
# With no factor there is nothing to compare, and the questions become: what
# is the trajectory, with what precision, and is there a rhythm and a trend at
# all. Under both approaches: every estimate has its 95% interval (the mixed
# model's amplitudes and acrophases included, which it used to omit), the tests
# against zero are reported, the residual-diagnostics tab lists them with the
# variance components, the overview shows them, the two approaches agree, and
# the exported script reproduces the table, tests included.
# Run with:   Rscript tests/single_traj_app_test.R      (from the DANCE directory)
# ==============================================================================
src <- readLines(if (dir.exists("tests")) "tests/structures_app_test.R" else "DANCE/tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^testServer", src) - 1)]))
txt <- function(x) gsub("\\s+", " ", gsub("<[^>]+>", " ", paste(unlist(x), collapse = " ")))

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  res <- list()
  k <- 0L
  for (ap in c("mixed", "two_stage")) {
    k <- k + 1L
    session$setInputs(harmonic_approach = ap, harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                      harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                      harmonic_group_var = "_none_", harmonic_design = "between", harmonic_between_var = "_none_",
                      harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                      harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                      harmonic_traj_tau = 18, harmonic_traj_tau_mode = "fixed", harmonic_ts_tau_mode = "fixed",
                      harmonic_tau_lo = 2, harmonic_tau_hi = 60, run_harmonic = k)
    oc <- harmonic_outcomes(); res[[ap]] <- oc
    es <- oc$estimates; ts <- oc$tests
    a1 <- es[es$outcome == "amplitude1", ]; p1 <- es[es$outcome == "phase1", ]
    chk(nrow(a1) == 1 && all(is.finite(c(a1$lo, a1$hi))) && nrow(p1) == 1 && is.finite(p1$lo),
        sprintf("%s: the single trajectory has its amplitude and acrophase with intervals", ap),
        sprintf("%s: amplitude or acrophase missing", ap))
    want <- c("circadian", "rhythm1", "rhythm2", "change")
    chk(all(want %in% ts$outcome) && all(is.finite(ts$p[ts$outcome %in% want])),
        sprintf("%s: tests against zero for the rhythm (joint and per harmonic) and the trend", ap),
        sprintf("%s: tests against zero missing: %s", ap, paste(setdiff(want, ts$outcome), collapse = ", ")))
    chk(ts$p[ts$outcome == "rhythm1"] < .001 && ts$p[ts$outcome == "rhythm2"] > .05 && ts$p[ts$outcome == "change"] < .001,
        sprintf("%s: the simulated 24 h rhythm and trend are found, the absent 12 h rhythm is not", ap),
        sprintf("%s: the tests against zero contradict the simulation", ap))
    t5 <- txt(render(paste(ap, "tab 5 panel"), output$harmonic_single_traj))
    chk(grepl("The single trajectory: estimates, intervals and tests", t5) && grepl("Tests against zero", t5) &&
          grepl(dance_outcome_value_text(a1, oc$period, oc$clock_origin, digits = 3), t5, fixed = TRUE) &&
          grepl(if (ap == "mixed") "Variance components" else "Between-participant spread", t5),
        sprintf("%s: the residual-diagnostics tab lists the estimates, the tests and the spread", ap),
        sprintf("%s: the residual-diagnostics panel is incomplete", ap))
    sv <- txt(render(paste(ap, "overview"), output$harmonic_streamlined))
    chk(grepl("Test against zero", sv) && grepl(dance_outcome_stat_text(ts[ts$outcome == "rhythm1", ]), sv, fixed = TRUE),
        sprintf("%s: the overview shows the tests against zero", ap), sprintf("%s: the overview lacks them", ap))
    ctx <- dance_eq_context(values$harmonic_model)
    chk(all(vapply(paste0("zero:", want), function(k) grepl("H<sub>0</sub>", dance_equation(k, ctx)$html), TRUE)),
        sprintf("%s: each test against zero has its equation", ap), sprintf("%s: a test against zero has no equation", ap))
    # the exported script reproduces the table, the tests against zero included
    code <- generate_analysis_code()
    env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
    ok <- tryCatch({ suppressWarnings(suppressMessages(capture.output(eval(parse(text = code), envir = env)))); TRUE },
                   error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
    so <- if (ap == "mixed") env$traj_outcomes else env$ts_outcomes
    chk(ok && !is.null(so) && isTRUE(all.equal(so$estimates$estimate, es$estimate, tolerance = 1e-6)) &&
          identical(so$tests$outcome, ts$outcome) && isTRUE(all.equal(so$tests$statistic, ts$statistic, tolerance = 1e-6)),
        sprintf("%s: the exported script reproduces every estimate and test", ap),
        sprintf("%s: the script's table differs from the app's", ap))
  }
  # the two approaches estimate the same trajectory
  am <- res$mixed$estimates; at <- res$two_stage$estimates
  for (q in c("fitted_mean", "change", "amplitude1")) {
    a <- am$estimate[am$outcome == q][1]; b <- at$estimate[at$outcome == q][1]
    chk(abs(a - b) < 0.05 * max(abs(c(a, b))), sprintf("%s: mixed %.3f and two-stage %.3f agree", q, a, b),
        sprintf("%s: the approaches disagree (%.3f vs %.3f)", q, a, b))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
