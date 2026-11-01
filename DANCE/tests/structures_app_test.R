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
render <- function(label, expr) {
  v <- tryCatch(force(expr), error = function(e) structure(conditionMessage(e), class = "render_err"))
  if (inherits(v, "render_err") && !grepl("^$", v) &&
      !inherits(attr(v, "condition"), "shiny.silent.error")) {
    # validate()/req() messages are deliberate refusals, reported as such
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
               list("combined", "linear"), list("combined", "log"), list("combined", "exp_sat"))
  k <- 0L
  for (ap in c("two_stage", "mixed")) for (r in runs) {
    st <- r[[1]]; tr <- r[[2]]; k <- k + 1L
    tag <- sprintf("%s / %s / %s", ap, st, tr)
    cat("\n--", tag, "--\n")
    session$setInputs(harmonic_approach = ap, harmonic_structure = st,
                      harmonic_trend_type = if (tr == "none") "linear" else tr,
                      harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between",
                      harmonic_between_var = "group", harmonic_df_method = "satterthwaite",
                      harmonic_traj_adjust = "holm", harmonic_traj_what = "intercept",
                      harmonic_traj_resid = "independent", harmonic_tau_lo = 2, harmonic_tau_hi = 60,
                      harmonic_traj_tau = 14, harmonic_traj_tau_mode = "fixed",
                      harmonic_ts_tau_mode = "estimate", harmonic_show_data = TRUE,
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
    sm <- render(paste(tag, "summary"), paste(capture.output(output$harmonic_summary), collapse = "\n"))
    if (!is.null(sm) && st == "trend")
      chk(!grepl("Amplitude|Acrophase|acrophase", sm), sprintf("%s: the summary reports no rhythm", tag),
          sprintf("%s: the summary of a trend-only model mentions amplitude/acrophase", tag))
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
