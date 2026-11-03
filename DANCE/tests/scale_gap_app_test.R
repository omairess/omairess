# ==============================================================================
# tests/scale_gap_app_test.R -- the log scale and the unsampled intervals
#
# Cortisol-like data sampled 08:00-18:00 and 06:00-12:00 the next day: nothing
# overnight. Group A peaks in the morning; group B's rhythm peaks at midnight,
# inside the gap. Under both approaches and both scales:
#   - on the original scale the two-harmonic curve dips below zero overnight,
#     and the summary says so and flags the minimum as extrapolated;
#   - on the log scale (the user's choice) every fit is on ln(y); the plots read
#     back on the original scale are positive, carry the raw data unchanged and
#     hover the drawn value; the tables give geometric means and peak/trough
#     ratios; a non-positive value is refused unless c makes it positive;
#   - the gap is shaded, the curve dashed across it, group B's acrophase is
#     flagged; and the exported script reproduces the outcome tables, notes
#     included.
# Run with:   Rscript tests/scale_gap_app_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({ library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "shinyWidgets")) if (requireNamespace(p, quietly = TRUE)) library(p, character.only = TRUE) })
app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }
render <- function(label, expr) {
  v <- tryCatch(force(expr), shiny.silent.error = function(e) NULL,
                error = function(e) structure(conditionMessage(e), class = "render_err"))
  if (inherits(v, "render_err")) { cat("FAIL:", label, "->", as.character(v), "\n"); failures <<- failures + 1L; return(NULL) }
  v
}

set.seed(5)
clock <- c(8, 10, 12, 14, 16, 18, 6, 8, 10, 12); tt <- c(8, 10, 12, 14, 16, 18, 30, 32, 34, 36)
n <- 16; grp <- rep(c("A", "B"), each = n / 2)
peak <- ifelse(grp == "A", 8.5, 24)
Y <- t(sapply(seq_len(n), function(i) round(exp(log(1.2) + 0.8 * cos(2 * pi * (tt - peak[i]) / 24) +
                                                  0.25 * cos(4 * pi * (tt - 6) / 24) + rnorm(10, 0, 0.12)), 3)))
d <- data.frame(ID = sprintf("C%02d", 1:n), group = grp, Y)
names(d)[-(1:2)] <- c(sprintf("Cort_%02d", clock[1:6]), sprintf("Cort_%02d_D2", clock[7:10]))
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE)
server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())

traces <- function(js) js$x$data
curve_y <- function(js) unlist(lapply(traces(js), function(tr)
  if (grepl("lines", tr$mode %||% "") && identical(tr$fill %||% "none", "none")) unlist(tr$y)))
marker_y <- function(js) unlist(lapply(traces(js), function(tr) if (identical(tr$mode, "markers")) unlist(tr$y)))
ytitle <- function(js) { yt <- js$x$layout$yaxis$title; if (is.list(yt)) yt$text %||% "" else yt %||% "" }

testServer(server, {
  session$setInputs(datafile = list(name = "cort.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(d)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  for (ap in c("two_stage", "mixed")) for (sc in c("original", "log")) {
    k <- k + 1L; tag <- sprintf("%s / %s", ap, sc)
    cat("\n--", tag, "--\n")
    session$setInputs(harmonic_approach = ap, harmonic_structure = "rhythm", harmonic_trend_type = "linear",
                      harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between", harmonic_between_var = "group",
                      harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                      harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                      harmonic_show_data = TRUE, harmonic_scale = sc, harmonic_log_offset = 0, harmonic_gap_h = NA,
                      harmonic_show_original = TRUE, harmonic_traj_original = TRUE, harmonic_traj_raw = TRUE,
                      harmonic_dv_name = "Cortisol", harmonic_dv_units = "nmol/L", harmonic_dv_min = 0,
                      run_harmonic = k)
    hm <- values$harmonic_model
    ok_fit <- !is.null(hm) && (ap == "two_stage" || isTRUE(hm$traj$ok))
    chk(ok_fit && identical(hm$scale$type, sc), sprintf("%s: fitted on the %s scale", tag, sc),
        sprintf("%s: the fit is missing or on the wrong scale", tag))
    if (!ok_fit) next
    if (sc == "log")
      chk(isTRUE(all.equal(hm$Y, log(values$data), check.attributes = FALSE)),
          sprintf("%s: the model was given ln(y)", tag), sprintf("%s: the model's data are not ln(y)", tag))
    # ---- the unsampled interval ----
    gp <- hm$gaps
    chk(nrow(gp) == 1 && gp$start == 10 && gp$end == 22,
        sprintf("%s: one unsampled interval, 18:00-06:00 (model %s-%s h)", tag, gp$start[1], gp$end[1]),
        sprintf("%s: unsampled intervals %s", tag, paste(gp$start, gp$end, collapse = "; ")))
    sm <- render(paste(tag, "summary"), output$harmonic_summary)
    chk(grepl(if (sc == "log") "Scale: natural log, ln(Cortisol)" else "Scale: original values", sm, fixed = TRUE) &&
          grepl("Unsampled intervals 18:00-06:00", sm, fixed = TRUE) && grepl("extrapolated, not observed", sm, fixed = TRUE),
        sprintf("%s: the summary states the scale and the gap, and flags an extreme inside it", tag),
        sprintf("%s: the summary does not state the scale/gap or flag the extreme", tag))
    if (sc == "original") {
      # the problem itself: the two-harmonic curve goes below zero overnight
      cy0 <- curve_y(jsonlite::fromJSON(render(paste(tag, "fit plot"), output$harmonic_fit_plot), simplifyVector = FALSE))
      chk(min(cy0, na.rm = TRUE) < 0, sprintf("%s: on the original scale the curve dips to %.2f overnight", tag, min(cy0, na.rm = TRUE)),
          sprintf("%s: the original-scale curve does not dip below zero here", tag))
      if (ap == "mixed")
        chk(grepl("leaves the admissible range", sm, fixed = TRUE),
            sprintf("%s: and the summary flags it as outside the admissible range", tag),
            sprintf("%s: the dip below zero is not flagged", tag))
    }
    else chk(!grepl("leaves the admissible range", sm, fixed = TRUE),
             sprintf("%s: on the log scale the curve stays admissible", tag),
             sprintf("%s: a log-scale curve leaves the admissible range", tag))
    for (o in c("harmonic_fit_plot", "harmonic_traj_curves")) {
      pj <- render(paste(tag, o), output[[o]])
      if (is.null(pj)) next
      js <- jsonlite::fromJSON(pj, simplifyVector = FALSE)
      tr <- traces(js)
      sh <- js$x$layout$shapes
      rect <- Filter(function(s) identical(s$type, "rect") && isTRUE(all.equal(as.numeric(s$x0), 10)) && isTRUE(all.equal(as.numeric(s$x1), 22)), sh)
      dash <- Filter(function(t) identical(t$line$dash, "dot"), tr)
      dx <- unlist(lapply(dash, function(t) { x <- unlist(t$x); y <- unlist(t$y); x[!vapply(t$y, is.null, TRUE)] }))
      chk(length(rect) == 1 && length(dash) >= 1 && all(dx >= 9.5 & dx <= 22.5),
          sprintf("%s %s: the gap is shaded and the curve dashed only across it", tag, o),
          sprintf("%s %s: shading %d, dashed traces %d", tag, o, length(rect), length(dash)))
      if (sc == "log") {
        cy <- curve_y(js); my <- marker_y(js)
        chk(min(cy, na.rm = TRUE) > 0 && grepl("read back from the log scale", ytitle(js)),
            sprintf("%s %s: read back, the curve is positive (min %.2f)", tag, o, min(cy, na.rm = TRUE)),
            sprintf("%s %s: the read-back curve reaches %.2f, axis '%s'", tag, o, min(cy, na.rm = TRUE), ytitle(js)))
        if (length(my))
          chk(isTRUE(all.equal(sort(my), sort(as.numeric(values$data[is.finite(values$data)])), tolerance = 1e-8)) ||
                isTRUE(all.equal(range(my), range(values$data, na.rm = TRUE), tolerance = 1e-8)),
              sprintf("%s %s: the raw data are drawn as observed", tag, o),
              sprintf("%s %s: the raw data are not drawn as observed (range %s)", tag, o, paste(round(range(my), 3), collapse = "-")))
        wt <- Filter(function(t) is.character(unlist(t$text)) && length(unlist(t$text)) == length(unlist(t$y)), tr)
        if (length(wt)) {
          t1 <- wt[[1]]; yy <- unlist(t1$y); tx <- unlist(t1$text); j <- which(is.finite(yy))[1]
          num <- as.numeric(sub(".*<br>[^<:]+: (-?[0-9.]+).*", "\\1", tx[j]))
          chk(isTRUE(abs(num - yy[j]) < 0.006), sprintf("%s %s: the hover shows the drawn value", tag, o),
              sprintf("%s %s: hover %s vs drawn %s", tag, o, num, yy[j]))
        }
      }
    }
    if (sc == "log") {
      session$setInputs(harmonic_show_original = FALSE)
      js <- jsonlite::fromJSON(render(paste(tag, "log-scale plot"), output$harmonic_fit_plot), simplifyVector = FALSE)
      chk(ytitle(js) == "ln(Cortisol)" && min(curve_y(js), na.rm = TRUE) < 0,
          sprintf("%s: on request the plot stays on the log scale, labelled ln(Cortisol)", tag),
          sprintf("%s: the log-scale view is not labelled or not on the log scale ('%s')", tag, ytitle(js)))
      session$setInputs(harmonic_show_original = TRUE)
    }
    # ---- the outcome tables ----
    oc <- render(paste(tag, "outcomes"), harmonic_outcomes())
    es <- oc$estimates
    if (sc == "log") {
      m <- es[es$outcome == "fitted_mean", ][1, ]
      gm <- as.numeric(sub(".*geometric mean ([0-9.]+).*", "\\1", m$note))
      a1 <- es[es$outcome == "amplitude1", ][1, ]
      rt <- as.numeric(sub(".*exp\\(2A\\) ([0-9.]+).*", "\\1", a1$note))
      chk(isTRUE(abs(gm - exp(m$estimate)) < 0.001) && isTRUE(abs(rt - exp(2 * a1$estimate)) < 0.001),
          sprintf("%s: the table reads back geometric means (%.3f) and peak/trough ratios (%.3f)", tag, gm, rt),
          sprintf("%s: the original-scale notes are missing or wrong", tag))
    } else chk(!any(grepl("geometric mean", es$note)), sprintf("%s: no log readings on the original scale", tag),
               sprintf("%s: log readings on an original-scale run", tag))
    pb <- es[es$outcome == "phase1" & es$cell == "B", ]
    pa <- es[es$outcome == "phase1" & es$cell == "A", ]
    chk(nrow(pb) == 1 && grepl("unsampled interval", pb$note) && !isTRUE(grepl("unsampled interval", pa$note)),
        sprintf("%s: group B's midnight acrophase is flagged as extrapolated, group A's morning one is not", tag),
        sprintf("%s: acrophase flags wrong (B: %s; A: %s)", tag, pb$note %||% "-", pa$note %||% "-"))
    rep <- render(paste(tag, "APA report"), dance_apa_report(values, input, "T"))
    chk(any(grepl("Unsampled intervals", rep)) && (sc == "original" || any(grepl("geometric mean", rep))),
        sprintf("%s: the report states the gap%s", tag, if (sc == "log") " and reads back the levels" else ""),
        sprintf("%s: the report lacks the gap or the read-back", tag))
    # ---- the exported script reproduces the tables, notes included ----
    code <- render(paste(tag, "R script"), generate_analysis_code())
    if (!is.null(code)) {
      env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
      ok <- tryCatch({ suppressWarnings(suppressMessages(capture.output(eval(parse(text = code), envir = env)))); TRUE },
                     error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
      so <- if (ap == "mixed") env$traj_outcomes else env$ts_outcomes
      ses <- so$estimates
      if (ok && !is.null(ses)) {
        cat(sprintf("   app %d rows, script %d rows\n", nrow(es), nrow(ses)))
        if (nrow(ses) == nrow(es)) { dd <- which(abs(ses$estimate - es$estimate) > 1e-6 | !(ses$note %in% es$note & es$note %in% ses$note) | xor(is.na(ses$note), is.na(es$note)))
          if (length(dd)) print(cbind(es[dd, c("outcome", "cell", "estimate", "note")], script_est = ses$estimate[dd], script_note = ses$note[dd])) }
      } else if (ok) cat("   script has no outcome table; objects:", paste(head(ls(env), 40), collapse = " "), "\n")
      same <- ok && !is.null(ses) && nrow(ses) == nrow(es) &&
        isTRUE(all.equal(ses$estimate, es$estimate, tolerance = 1e-6)) && identical(ses$note, es$note)
      chk(same, sprintf("%s: the exported script reproduces every estimate and note", tag),
          sprintf("%s: the script's outcome table differs from the app's", tag))
      # the acrophases as clock times: the script's models carry the clock origin
      ph <- which(es$kind == "circular")
      ctx <- function(o, i) dance_outcome_value_text(o$estimates[i, ], o$period, o$clock_origin)
      chk(ok && length(ph) && identical(vapply(ph, function(i) ctx(oc, i), ""), vapply(ph, function(i) ctx(so, i), "")),
          sprintf("%s: the script prints the same acrophase clock times (%s)", tag, ctx(oc, ph[1])),
          sprintf("%s: the script's acrophase clock times differ from the app's", tag))
    }
  }
  # ---- a value that has no logarithm ----
  cat("\n-- a non-positive value --\n")
  d0 <- d; d0[3, "Cort_14"] <- 0
  csv0 <- tempfile(fileext = ".csv"); write.csv(d0, csv0, row.names = FALSE)
  session$setInputs(datafile = list(name = "cort0.csv", datapath = csv0), load_data = 2)
  session$setInputs(sel_data_vars = names(d0)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID", apply_selection = 2)
  values$harmonic_model <- NULL
  session$setInputs(harmonic_approach = "two_stage", harmonic_scale = "log", harmonic_log_offset = 0, run_harmonic = 10)
  chk(is.null(values$harmonic_model), "a zero is refused on the log scale with c = 0", "a zero was logged")
  session$setInputs(harmonic_log_offset = 0.5, run_harmonic = 11)
  hm <- values$harmonic_model
  chk(!is.null(hm) && hm$scale$offset == 0.5 && isTRUE(all.equal(hm$Y, log(values$data + 0.5), check.attributes = FALSE)),
      "with c = 0.5 it runs on ln(y + 0.5), and says so",
      "c was not applied")
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
