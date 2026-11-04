# ==============================================================================
# tests/structures_app_test.R -- the three model structures through the app
#
# Rhythmic component only, non-periodic trend only, and both, under each
# approach (two-stage and mixed-effects) and each trend family, run through the
# server and every cosinor output forced to render: the summary, the tables,
# every plot, the comparison tab, the APA report, the exported script and the
# pairwise CSV. An output that errors fails the test, and so does an output that
# offers a component the model does not have.
#
# Run with:   Rscript tests/structures_app_test.R      (from the DANCE directory)
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
txt <- function(x) gsub("\\s+", " ", gsub("<[^>]+>", " ", paste(unlist(x), collapse = " ")))

# 24 participants in two groups, 16 two-hourly samples from 08:00 over 30 h,
# a saturating rise that differs by group plus a 24 h rhythm; a few gaps
set.seed(41)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 14)
tt <- seq(0, 30, by = 2)
n <- 24; grp <- rep(c("ctrl", "sd"), each = n / 2)
Y <- t(vapply(seq_len(n), function(i)
  3 + rnorm(1, 0, 0.6) + (2 + 1.5 * (grp[i] == "sd")) * (1 - exp(-tt / 14)) +
    1.2 * cos(2 * pi * (tt + 8 - 4) / 24) + rnorm(16, 0, 0.5), numeric(16)))
Y[2, 5] <- NA; Y[15, c(3, 11)] <- NA
d <- data.frame(ID = sprintf("P%02d", seq_len(n)), group = grp, round(Y, 2))
names(d)[-(1:2)] <- c(paste0("KSS_", clock[1:12]), paste0("KSS_", clock[13:16], "_SD"))
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE, na = "")

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
# validate() and req() are deliberate: "this does not apply to the model" is a
# correct answer. Anything else that errors is a failure.
render <- function(label, expr) {
  v <- tryCatch(force(expr),
                shiny.silent.error = function(e) structure(conditionMessage(e), class = "render_msg"),
                error = function(e) structure(conditionMessage(e), class = "render_err"))
  if (inherits(v, "render_msg")) return(invisible(structure(as.character(v), class = "render_msg")))
  if (inherits(v, "render_err")) {
    cat("FAIL:", label, "->", as.character(v), "\n"); failures <<- failures + 1L
    return(invisible(NULL))
  }
  invisible(v)
}

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group",
                    subject_id_var = "ID", apply_selection = 1)
  runs <- list(list("rhythm", "none"),
               list("trend", "linear"), list("trend", "log"), list("trend", "exp_sat"),
               list("combined", "linear"), list("combined", "log"), list("combined", "exp_sat"),
               # tau estimated: per participant (two-stage); profiled, and per cell (mixed)
               list("combined", "exp_sat", tau = "estimated"),
               # the model comparison over every structure
               list("trend", "linear", ms = TRUE))
  k <- 0L
  # DANCE_STRUCT_ONLY="two_stage/trend" (a regular expression on the tag)
  # runs a subset while debugging; unset, every combination runs
  only <- Sys.getenv("DANCE_STRUCT_ONLY", "")
  for (ap in c("two_stage", "mixed")) for (r in runs) {
    st <- r[[1]]; tr <- r[[2]]; k <- k + 1L
    tau_est <- identical(r$tau, "estimated"); do_ms <- isTRUE(r$ms)
    tag <- sprintf("%s / %s / %s%s%s", ap, st, tr, if (tau_est) " / tau estimated" else "",
                   if (do_ms) " / model comparison" else "")
    if (nzchar(only) && !grepl(only, gsub(" ", "", tag))) next
    cat("\n--", tag, "--\n")
    session$setInputs(harmonic_approach = ap, harmonic_structure = st,
                      harmonic_trend_type = if (tr == "none") "linear" else tr,
                      harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between",
                      harmonic_between_var = "group", harmonic_df_method = "satterthwaite",
                      harmonic_traj_adjust = "holm", harmonic_traj_what = "intercept",
                      harmonic_traj_resid = "independent", harmonic_tau_lo = 2, harmonic_tau_hi = 60,
                      harmonic_traj_tau = 14,
                      harmonic_traj_tau_mode = if (tau_est) "profile" else "fixed",
                      harmonic_traj_sat_cells = tau_est,
                      harmonic_ts_tau_mode = if (tau_est) "estimate" else "fixed", harmonic_show_data = TRUE,
                      harmonic_model_selection = do_ms,
                      run_harmonic = k)
    hm <- values$harmonic_model
    K_exp <- if (st == "trend") 0L else 2L
    tr_exp <- if (st == "rhythm") "none" else tr
    ok_fit <- !is.null(hm) && (if (ap == "mixed") isTRUE(hm$traj$ok) else length(hm$individual_fits) > 0)
    chk(ok_fit && identical(as.integer(hm$n_harmonics), K_exp) && identical(hm$trend_type, tr_exp),
        sprintf("%s: fitted with %d harmonic(s) and trend '%s'", tag, K_exp, tr_exp),
        sprintf("%s: the fit is missing or has K = %s, trend = %s", tag,
                format(hm$n_harmonics), format(hm$trend_type)))
    if (!ok_fit) next
    sm <- render(paste(tag, "summary"), output$harmonic_summary)
    if (do_ms && !is.null(sm))
      chk(grepl("Candidates: intercept only", sm, fixed = TRUE) && grepl("<-- reported", sm, fixed = TRUE),
          sprintf("%s: the summary lists every candidate structure and marks the reported one", tag),
          sprintf("%s: the model comparison is missing from the summary", tag))
    if (!is.null(sm) && st == "trend") {
      ln <- strsplit(sm, "\n")[[1]]
      bad <- grep("Amplitude|Acrophase|acrophase", ln[!grepl("no rhythm", ln)], value = TRUE)
      chk(!length(bad), sprintf("%s: the summary reports no rhythm", tag),
          sprintf("%s: the summary of a trend-only model mentions amplitude/acrophase: %s", tag,
                  paste(head(trimws(bad), 3), collapse = " | ")))
    }
    for (o in c("harmonic_summary_box", "harmonic_parameters_table", "harmonic_gof_stats",
                "harmonic_fit_plot", "harmonic_residual_plot", "harmonic_qq_plot",
                "harmonic_polar_plot", "harmonic_density_plot", "harmonic_density_note",
                "harmonic_subject_selector", "harmonic_selector_polar",
                "harmonic_traj_header", "harmonic_traj_primary", "harmonic_traj_components",
                "harmonic_traj_curves", "harmonic_traj_band_note", "harmonic_traj_levels_plot",
                "harmonic_traj_levels_note", "harmonic_traj_effect_ui", "harmonic_traj_diff_a_ui",
                "harmonic_traj_diff_b_ui", "harmonic_traj_diff_plot", "harmonic_traj_diff_note",
                "harmonic_traj_pairwise", "harmonic_traj_sat_header", "harmonic_traj_sat_plot",
                "harmonic_traj_sat_tests"))
      render(paste(tag, o), output[[o]])
    # the three questions: only the domains the model has
    oc <- render(paste(tag, "outcomes"), harmonic_outcomes())
    render(paste(tag, "outcomes panel"), output$harmonic_outcomes)
    render(paste(tag, "outcomes CSV"), output$download_harmonic_outcomes)
    if (!is.null(oc) && isTRUE(oc$ok)) {
      doms <- unique(c(oc$estimates$domain, oc$tests$domain))
      chk(("rhythm" %in% doms) == (st != "trend") && ("accumulation" %in% doms) == (st != "rhythm") &&
            "level" %in% doms,
          sprintf("%s: the outcome domains are those of the model (%s)", tag, paste(doms, collapse = ", ")),
          sprintf("%s: outcome domains %s do not match the structure", tag, paste(doms, collapse = ", ")))
      chk(any(oc$tests$domain == "omnibus"), sprintf("%s: the omnibus trajectory test is reported", tag),
          sprintf("%s: no omnibus test", tag))
    }
    # THE STREAMLINED VIEW, the default: the omnibus test as the primary
    # analysis, a row per question, and the same numbers as the outcome table
    sv <- render(paste(tag, "streamlined view"), output$harmonic_streamlined)
    if (!is.null(oc) && isTRUE(oc$ok) && !is.null(sv)) {
      svt <- gsub("\\s+", " ", gsub("<[^>]+>", " ", paste(unlist(sv), collapse = " ")))
      om <- oc$tests[oc$tests$domain == "omnibus", ][1, ]
      lvp <- oc$estimates[oc$estimates$domain == "level" & oc$estimates$primary & is.finite(oc$estimates$estimate), ]
      fm <- lvp[1, ]
      chk(grepl("PRIMARY ANALYSIS", svt) && grepl(dance_outcome_stat_text(om), svt, fixed = TRUE) &&
            grepl(dance_outcome_value_text(fm, oc$period, oc$clock_origin, digits = 2), svt, fixed = TRUE) &&
            grepl("Overall response level", svt) && (st == "rhythm" || grepl("Non-periodic trend", svt)) &&
            (st == "trend" || grepl("RHYTHM DESCRIPTORS", svt)),
          sprintf("%s: the streamlined view shows the omnibus test and the table's own estimates", tag),
          sprintf("%s: the streamlined view disagrees with the outcome table [%s]", tag, paste(c(
            if (!grepl("PRIMARY ANALYSIS", svt)) "primary", if (!grepl(dance_outcome_stat_text(om), svt, fixed = TRUE)) dance_outcome_stat_text(om),
            if (!grepl(dance_outcome_value_text(fm, oc$period, oc$clock_origin, digits = 2), svt, fixed = TRUE))
              paste("fm", dance_outcome_value_text(fm, oc$period, oc$clock_origin, digits = 2)),
            if (!grepl("Non-periodic trend", svt)) "trend row", if (!grepl("RHYTHM DESCRIPTORS", svt)) "rhythm"), collapse = "; ")))
      # the overview and the detailed table carry the same tests: every
      # level/accumulation primary and every rhythm test of the table
      want <- oc$tests[(oc$tests$domain == "rhythm" & oc$tests$outcome %in% c("circadian", sprintf("rhythm%d", 1:5))) |
                         (oc$tests$domain %in% c("level", "accumulation") & oc$tests$primary), , drop = FALSE]
      miss <- want$label[!vapply(seq_len(nrow(want)), function(i)
        grepl(dance_outcome_stat_text(want[i, ]), svt, fixed = TRUE), TRUE)]
      chk(!length(miss), sprintf("%s: every rhythm and primary test of the table is in the overview (%d)", tag, nrow(want)),
          sprintf("%s: the overview lacks %s", tag, paste(miss, collapse = "; ")))
      # every quantity reported has its own equation, written for this model
      ctx <- dance_eq_context(hm)
      keys <- unique(c(oc$estimates$outcome, paste0("test:", stats::na.omit(oc$tests$outcome))))
      bad <- keys[vapply(keys, function(k) {
        q <- dance_equation(k, ctx); grepl("No equation recorded", q$html) || identical(q$title, sub("^test:", "", k)) }, TRUE)]
      chk(!length(bad), sprintf("%s: every reported quantity has its equation (%d keys)", tag, length(keys)),
          sprintf("%s: no equation for %s", tag, paste(bad, collapse = ", ")))
      # the individual tables: every row and column has its equation
      for (o in c("harmonic_parameters_table", "harmonic_individual_eq")) {
        v <- render(paste(tag, o), output[[o]])
        vv <- gsub("&#39;", "'", paste(unlist(v), collapse = " "), fixed = TRUE)
        ks <- regmatches(vv, gregexpr("key: '[^']+'", vv))[[1]]
        ks <- unique(sub("key: '([^']+)'", "\\1", ks))
        bad <- ks[vapply(ks, function(k) { q <- dance_equation(k, ctx); grepl("No equation recorded", q$html) || identical(q$title, k) }, TRUE)]
        chk(length(ks) > 2 && !length(bad), sprintf("%s: %s has an equation for every entry (%d)", tag, o, length(ks)),
            sprintf("%s: %s: %d links, without an equation: %s", tag, o, length(ks), paste(bad, collapse = ", ")))
      }
      ok_modal <- tryCatch({ session$setInputs(dance_eq = list(key = "mesor", nonce = k)); TRUE }, error = function(e) conditionMessage(e))
      chk(isTRUE(ok_modal), sprintf("%s: an (i) click opens its equation", tag), sprintf("%s: the (i) click failed: %s", tag, ok_modal))
    }
    # a row of the table, clicked: its quantity is selected, its annotation drawn,
    # the graph moved to its view; choosing another quantity by hand ends it
    if (!is.null(oc) && isTRUE(oc$ok)) {
      es <- oc$estimates
      q <- es$outcome[es$domain == (if (st == "rhythm") "rhythm" else "accumulation")][1]
      session$setInputs(harmonic_traj_annot = character(0))
      session$setInputs(harmonic_outcome_click = list(outcome = q, nonce = k))
      rq <- harmonic_outcome_request()
      ann <- harmonic_traj_annot_set()
      chk(identical(harmonic_outcome_selected(), q) && !is.null(rq$view) && length(ann) > 0,
          sprintf("%s: clicking '%s' selects it (menu: %s), shows view '%s' and annotation %s", tag, q,
                  format(rq$what %||% "-"), format(rq$view %||% "-"), paste(ann, collapse = "+")),
          sprintf("%s: clicking '%s' did nothing (view %s, annotation %s)", tag, q,
                  format(rq$view %||% "-"), paste(ann, collapse = "+")))
      session$setInputs(harmonic_traj_what = "mesor")
      chk(is.null(harmonic_outcome_selected()) || identical(rq$what, "mesor"),
          sprintf("%s: a hand-picked quantity ends the row's selection", tag),
          sprintf("%s: the clicked row stays selected after another quantity is chosen", tag))
      session$setInputs(harmonic_traj_what = "intercept")
    }
    # the three graphs, every annotation on, and the drawn numbers against the
    # outcome table they are read from
    session$setInputs(harmonic_traj_annot = c("mean", "change", "tangent", "asymptote", "acro"))
    for (cmp in intersect(c("full", "trend", "harmonics"),
                          if (ap == "mixed") dance_traj_components(hm$traj) else c("full", "harmonics", if (st != "rhythm") "trend"))) {
      session$setInputs(harmonic_traj_component = cmp)
      pj <- render(paste(tag, "graph", cmp), output$harmonic_traj_curves)
      if (is.null(pj) || is.null(oc) || !isTRUE(oc$ok)) next
      js <- jsonlite::fromJSON(pj, simplifyVector = FALSE)$x$data
      ht <- vapply(js, function(tr) paste(unlist(tr$hovertemplate), collapse = " "), "")
      es <- oc$estimates
      cl1 <- es$cell[1]
      if (cmp == "full") {
        m <- es$estimate[es$outcome %in% c("fitted_mean") & es$cell == cl1]
        if (!length(m)) m <- es$estimate[es$outcome == "mesor" & es$cell == cl1]
        tr <- js[[which(grepl(sprintf("^%s time-averaged fitted response", cl1), ht))[1]]]
        chk(!is.null(tr) && isTRUE(all.equal(unlist(tr$y)[1], m[1])),
            sprintf("%s: the mean line is the table's time-averaged fitted response", tag),
            sprintf("%s: the mean line does not match the table", tag))
      }
      if (cmp == "trend") {
        ch <- es$estimate[es$outcome %in% c("sat_change", "change") & es$cell == cl1]
        tr <- js[[which(grepl(sprintf("^%s total trend change", cl1), ht))[1]]]
        chk(!is.null(tr) && isTRUE(all.equal(unlist(tr$y)[1], ch[1])),
            sprintf("%s: the end marker is the table's total trend change", tag),
            sprintf("%s: the end marker does not match the table", tag))
      }
    }
    session$setInputs(harmonic_traj_component = "full")
    pw <- render(paste(tag, "pairwise export"), harmonic_pairwise_all())
    if (!is.null(pw)) {
      chk(!any(grepl("^(amplitude|phase)", pw$quantity)) || st != "trend",
          sprintf("%s: the pairwise export has no rhythm row", tag),
          sprintf("%s: a trend-only pairwise export compares amplitude/acrophase", tag))
      chk(!any(pw$quantity %in% c("mesor", "change", "avg_rate", "rate0")) || st != "rhythm",
          sprintf("%s: the pairwise export has no trend row", tag),
          sprintf("%s: a rhythm-only pairwise export compares a trend quantity", tag))
    }
    for (w in (if (ap == "mixed") harmonic_traj_whats(hm$traj) else dance_ts_pair_whats(harmonic_ts()$mod))) {
      session$setInputs(harmonic_traj_what = w)
      render(paste(tag, "pairwise", w), output$harmonic_traj_pairwise)
    }
    session$setInputs(harmonic_traj_what = "intercept")
    render(paste(tag, "APA report"), dance_apa_report(values, input, "T"))
    render(paste(tag, "R script"), generate_analysis_code())
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
