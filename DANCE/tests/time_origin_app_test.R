# ==============================================================================
# tests/time_origin_app_test.R -- a time origin before the first observation
#
# The recording starts at 08:00; the origin is set at 06:00 (e.g. wake time).
# With a linear or a saturating trend (tau held) this is a reparameterisation:
# the fitted curve, every window quantity (time-averaged levels, change, rates),
# the amplitudes and the acrophases must not move, and only the intercept and
# the trend coefficient are re-expressed at the origin:
#   linear:   b0' = b0 - 2 b                    b' = b
#   exp_sat:  A'  = A e^(2/tau)                 b0' = b0 + A - A'
# Under both approaches, and the exported script reproduces the run with the
# origin. Run with:   Rscript tests/time_origin_app_test.R   (from DANCE/)
# ==============================================================================
src <- readLines(if (dir.exists("tests")) "tests/structures_app_test.R" else "DANCE/tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^testServer", src) - 1)]))
tau <- 14

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  run <- function(ap, tr, origin) {
    k <<- k + 1L
    session$setInputs(harmonic_approach = ap, harmonic_structure = "combined", harmonic_trend_type = tr,
                      harmonic_period = 24, n_harmonics = 1, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between", harmonic_between_var = "group",
                      harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                      harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                      harmonic_traj_tau = tau, harmonic_traj_tau_mode = "fixed", harmonic_ts_tau_mode = "fixed",
                      harmonic_tau_lo = 2, harmonic_tau_hi = 60,
                      harmonic_origin = if (is.null(origin)) "first" else "clock",
                      harmonic_origin_clock = origin %||% "", run_harmonic = k)
    harmonic_outcomes()
  }
  val <- function(oc, key) { e <- oc$estimates; r <- e[e$outcome == key, ]; stats::setNames(r$estimate, r$cell) }
  for (ap in c("mixed", "two_stage")) for (tr in c("linear", "exp_sat")) {
    lab <- sprintf("%s, %s", ap, tr)
    o1 <- run(ap, tr, NULL)
    o2 <- run(ap, tr, "06:00")
    mod <- values$harmonic_model
    chk(isTRUE(all.equal(mod$origin_lead, 2)) && isTRUE(all.equal(mod$origin_shift, 6)) &&
          isTRUE(all.equal(min(mod$time_vec), 2)),
        sprintf("%s: the origin is 06:00, 2 h before the first observation", lab),
        sprintf("%s: origin lead %s, shift %s", lab, mod$origin_lead, mod$origin_shift))
    same <- c("fitted_mean", "mesor", "change", "avg_rate", "rate0", "amplitude1", "phase1")
    same <- intersect(same, o1$estimates$outcome)
    # an acrophase is stored in model hours: compare it on the clock
    clk <- function(oc, q) if (q == "phase1") (val(oc, q) + oc$clock_origin) %% 24 else val(oc, q)
    dev <- vapply(same, function(q) max(abs(clk(o1, q) - clk(o2, q))), numeric(1))
    # mixed: the reparameterised model is refitted by lme4, whose optimiser stops
    # within its convergence tolerance (differences of about 1e-4 seen)
    tol <- if (ap == "mixed") 1e-3 else 1e-5
    chk(all(dev < tol),
        sprintf("%s: %s are unchanged by the origin", lab, paste(same, collapse = ", ")),
        sprintf("%s: moved by the origin: %s", lab, paste(sprintf("%s %.2g", same, dev)[dev >= tol], collapse = ", ")))
    # the tests over the sampled window are the same tests
    # Mixed: the random-slope variance of the simulated data is near zero (lme4
    # reports a boundary fit while fitting), and near that boundary the
    # optimiser's stopping point -- and with it the Satterthwaite df -- depends
    # on the parameterisation. The F statistics then agree to about 1%, the
    # estimates to 1e-4. Two-stage refits are exact. The deviation is printed.
    t1 <- o1$tests; t2 <- o2$tests
    tk <- intersect(c("mesor", "change", "rhythm1", "circadian", "fitted_mean"), t1$outcome)
    sing <- ap == "mixed"
    stol <- if (sing) 0.02 else 10 * tol
    rel <- vapply(tk, function(q) { a <- t1$statistic[t1$outcome == q][1]; b <- t2$statistic[t2$outcome == q][1]
                                    abs(a - b) / max(1, abs(a)) }, numeric(1))
    chk(all(rel < stol),
        sprintf("%s: the group tests of %s are unchanged%s", lab, paste(tk, collapse = ", "),
                if (sing) sprintf(" (F within %.2f%%; df not compared)", 100 * max(rel)) else ""),
        sprintf("%s: a group test changed with the origin (%s)", lab,
                paste(sprintf("%s %.2g", tk, rel), collapse = ", ")))
    b1 <- val(o1, "intercept"); b2 <- val(o2, "intercept")
    if (tr == "linear") {
      s1 <- val(o1, "avg_rate")
      want <- b1 - 2 * s1
    } else {
      ch1 <- val(o1, "change"); Tw <- diff(range(mod$time_vec))
      A1 <- ch1 / (1 - exp(-Tw / tau))
      want <- b1 + A1 - A1 * exp(2 / tau)
    }
    chk(max(abs(b2 - want)) < 10 * tol && max(abs(b2 - b1)) > 1e-2,
        sprintf("%s: the intercept is re-expressed at the origin (%s)", lab,
                paste(sprintf("%.3f -> %.3f", b1, b2), collapse = "; ")),
        sprintf("%s: intercept %s, expected %s", lab, paste(round(b2, 4), collapse = ", "),
                paste(round(want, 4), collapse = ", ")))
    # THE PLOT drawn from the time origin (trend view): the curve starts at the
    # origin at 0, the annotations start on the curve at the first observation,
    # and the interval before it is marked as extrapolated. Unticked, the curve
    # and its annotations both start at 0 at the first observation.
    if (tr == "exp_sat") for (fo in c(TRUE, FALSE)) {
      session$setInputs(harmonic_traj_component = "trend", harmonic_traj_band = "pointwise",
                        harmonic_traj_annot = c("change", "tangent", "asymptote"),
                        harmonic_traj_from_origin = fo, harmonic_traj_raw = FALSE)
      b <- plotly::plotly_build(harmonic_plot_finish(harmonic_traj_curves_build(), FALSE,
                                                     pre_sample = harmonic_traj_from_origin()))$x
      tr_d <- Filter(function(d) identical(d$mode, "lines") && length(unlist(d$x)) > 20 &&
                       identical(d$fill %||% "none", "none"), b$data)
      xs <- lapply(tr_d, function(d) as.numeric(unlist(d$x))); ys <- lapply(tr_d, function(d) as.numeric(unlist(d$y)))
      x0 <- min(unlist(xs)); first <- vapply(seq_along(xs), function(i) ys[[i]][which.min(xs[[i]])], numeric(1))
      # the curves at the first observation (t = 2), from the solid and the
      # dotted pieces together
      at1 <- vapply(seq_along(xs), function(i) { ok <- is.finite(ys[[i]])
        if (sum(ok) < 2) NA_real_ else stats::approx(xs[[i]][ok], ys[[i]][ok], 2, rule = 2)$y }, numeric(1))
      at1 <- at1[is.finite(at1)]
      ann_l <- Filter(function(d) identical(d$mode, "lines") && length(unlist(d$x)) == 2 &&
                        any(grepl("average trend rate", unlist(d$hovertemplate %||% ""))), b$data)
      a_x <- vapply(ann_l, function(d) as.numeric(unlist(d$x))[1], numeric(1))
      a_y <- vapply(ann_l, function(d) as.numeric(unlist(d$y))[1], numeric(1))
      labs <- vapply(b$layout$annotations %||% list(), function(a) a$text %||% "", character(1))
      if (fo) {
        chk(abs(x0) < 1e-9 && all(abs(first) < 1e-6) && all(at1 > 0.05) &&
              length(a_x) >= 2 && all(abs(a_x - 2) < 1e-9) && all(vapply(a_y, function(y) min(abs(y - at1)) < 0.01, TRUE)) &&
              "before the first sample" %in% labs,
            sprintf("%s: drawn from the origin, the trend is 0 at 06:00, the change line starts on the curve at 08:00, and the pre-sample interval is marked", lab),
            sprintf("%s: the plot from the origin is wrong (x0 %.3g, first %s, rate line from %s at %s, labels %s)", lab, x0,
                    paste(round(first, 3), collapse = ","), paste(a_x, collapse = ","), paste(round(a_y, 3), collapse = ","),
                    paste(labs, collapse = "|")))
      } else {
        chk(abs(x0 - 2) < 1e-9 && all(abs(first) < 1e-6) && length(a_x) >= 2 && all(abs(a_x - 2) < 1e-9) && all(abs(a_y) < 1e-9) &&
              !"before the first sample" %in% labs,
            sprintf("%s: unticked, the curve and the change line start at 0 at the first observation", lab),
            sprintf("%s: the sampled-window plot is wrong (x0 %.3g, rate line from %s at %s)", lab, x0,
                    paste(a_x, collapse = ","), paste(a_y, collapse = ",")))
      }
    }
    # stated: in the assumptions and the summary
    asm <- paste(dance_cosinor_assumptions_text(dance_cosinor_assumptions(mod, if (ap == "mixed") harmonic_traj())),
                 collapse = " ")
    chk(grepl("time origin chosen at 06:00, 2 h before", asm) && grepl("extrapolation", asm),
        sprintf("%s: the assumptions state the origin and that the intercept is extrapolated", lab),
        sprintf("%s: the origin is not stated", lab))
    ctx <- dance_eq_context(mod)
    chk(grepl("time origin", dance_equation("intercept", ctx)$html) &&
          grepl("t</i><sub>1</sub>", dance_equation("mesor", ctx)$html),
        sprintf("%s: the equations show the sampled window since the origin", lab),
        sprintf("%s: the equations still integrate from 0", lab))
    # the exported script reproduces the run with the origin
    code <- generate_analysis_code()
    env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
    ok <- tryCatch({ suppressWarnings(suppressMessages(capture.output(eval(parse(text = code), envir = env)))); TRUE },
                   error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
    so <- if (ap == "mixed") env$traj_outcomes else env$ts_outcomes
    chk(ok && !is.null(so) && isTRUE(all.equal(so$estimates$estimate, o2$estimates$estimate, tolerance = 1e-6)) &&
          grepl("time origin chosen at 06:00", code, fixed = TRUE) &&
          (ap == "mixed" || grepl("t_origin = 0", code, fixed = TRUE)),
        sprintf("%s: the exported script reproduces every estimate with the origin", lab),
        sprintf("%s: the script's table differs from the app's", lab))
  }
  # a time that is not a time of day is refused, not guessed
  m0 <- values$harmonic_model
  run("mixed", "linear", "25:70")
  chk(identical(values$harmonic_model, m0), "an invalid origin is refused and nothing is refitted",
      "an invalid origin was accepted")
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
