# ==============================================================================
# tests/dlmo_app_test.R -- the DLMO module through the app
#
# Two layouts of the same kind of study (a baseline night and a sleep-
# deprivation night, two groups): both nights in one row, and one row per night
# named by a column. Onsets are known by construction. Every DLMO output is
# rendered, the report and the exported script are generated, and the script's
# onsets are checked against the app's.
#
# Run with:   Rscript tests/dlmo_app_test.R      (from the DANCE directory)
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
  v <- tryCatch(force(expr),
                shiny.silent.error = function(e) structure(conditionMessage(e), class = "render_msg"),
                error = function(e) structure(conditionMessage(e), class = "render_err"))
  if (inherits(v, "render_err")) { cat("FAIL:", label, "->", as.character(v), "\n"); failures <<- failures + 1L; return(invisible(NULL)) }
  invisible(v)
}

# 16 participants, hourly saliva samples 18:00-02:00 on two nights; the rise
# (6 pg/mL per hour from 1 pg/mL) starts at 20:00 or 21:00 at baseline and 1 h
# later after sleep deprivation, 2 h later in group "B". The rise starts on a
# sample, so linear interpolation is exact; a rise starting between hourly
# samples is interpolated early, a property of the method, tested in
# tests/testthat/test-dlmo.R.
set.seed(7)
n <- 16; clock <- c(18:23, 0:2); tt <- 18:26
grp <- rep(c("A", "B"), each = n / 2)
on1 <- sample(c(20, 21), n, TRUE); on2 <- on1 + ifelse(grp == "B", 2, 1)
prof <- function(on) round(1 + 6 * pmax(0, tt - on) + rnorm(9, 0, 0.05), 3)
Y1 <- t(sapply(on1, prof)); Y2 <- t(sapply(on2, prof))
bed <- sprintf("%02d:%02d", 23, sample(c(0, 15, 30), n, TRUE))
true_dlmo <- function(on) on + 2 / 6            # 3 pg/mL is reached 1/3 h after the rise starts

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())

run_layout <- function(label, d, sel, night_var, csv_name) {
  csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE, na = "")
  testServer(server, {
    cat("\n--", label, "--\n")
    session$setInputs(datafile = list(name = csv_name, datapath = csv), header = TRUE, sep = ",",
                      data_format = "wide", load_data = 1)
    session$setInputs(sel_data_vars = sel, sel_cov_vars = intersect(c("group", "night", "bed"), names(d)),
                      subject_id_var = "ID", apply_selection = 1)
    chk(!is.null(values$time_clock), sprintf("%s: clock times parsed from the column names", label),
        sprintf("%s: no clock times parsed", label))
    # the hockey-stick without dlmoR is refused, never approximated
    if (!requireNamespace("dlmoR", quietly = TRUE)) {
      session$setInputs(dlmo_time_src = "_shared_", dlmo_methods = c("fixed", "hockey"), dlmo_threshold = 3,
                        run_dlmo = 1)
      chk(isFALSE(values$dlmo_results$ok) && grepl("dlmoR", values$dlmo_results$message),
          sprintf("%s: the hockey-stick is refused without dlmoR", label),
          sprintf("%s: the hockey-stick ran without dlmoR", label))
    }
    # no bedtime first: no phase angle anywhere
    session$setInputs(dlmo_time_src = "_shared_", dlmo_units = "pg/mL", dlmo_night_var = night_var,
                      dlmo_group_var = "group", dlmo_gap_h = 6, dlmo_methods = c("fixed", "variable"),
                      dlmo_threshold = 3, dlmo_n_base = 3, dlmo_k_sd = 2, dlmo_confirm = 1, dlmo_max_gap = 2,
                      dlmo_lloq = NA, dlmo_lloq_rule = "half", dlmo_bed_src = "none", run_dlmo = 2)
    r <- values$dlmo_results
    chk(isTRUE(r$ok), sprintf("%s: the run succeeded", label), sprintf("%s: the run failed: %s", label, r$message))
    if (!isTRUE(r$ok)) return(invisible())
    tab <- r$tab
    chk(nrow(tab) == 2 * 2 * n && setequal(unique(tab$night), if (night_var == "_none_") c("Night 1", "Night 2") else c("baseline", "SD")),
        sprintf("%s: %d profiles x 2 methods, two nights", label, 2 * n),
        sprintf("%s: %d rows, nights %s", label, nrow(tab), paste(unique(tab$night), collapse = ",")))
    fx <- tab[tab$method == "fixed", ]
    n1 <- fx[fx$night %in% c("Night 1", "baseline"), ]; n2 <- fx[fx$night %in% c("Night 2", "SD"), ]
    err <- max(abs(c(n1$dlmo_h[match(sprintf("P%02d", 1:n), n1$id)] - true_dlmo(on1),
                     n2$dlmo_h[match(sprintf("P%02d", 1:n), n2$id)] - true_dlmo(on2))), na.rm = TRUE)
    chk(all(fx$status == "determined") && err < 0.05,
        sprintf("%s: every fixed-threshold onset recovered (max error %.3f h)", label, err),
        sprintf("%s: onsets not recovered (max error %.3f h; statuses %s)", label, err,
                paste(unique(fx$status), collapse = ",")))
    chk(all(is.na(tab$phase_angle_h)) && !r$has_bed, sprintf("%s: no bedtime, no phase angle", label),
        sprintf("%s: a phase angle without a bedtime", label))
    cp <- dance_dlmo_comparisons(r, "fixed")
    ct <- cp$test$tests
    chk(isTRUE(cp$test$ok) && "group:night" %in% ct$effect &&
          abs(ct$estimate[grepl("in A$", ct$effect)] - 1) < 0.2 && abs(ct$estimate[grepl("in B$", ct$effect)] - 2) < 0.2,
        sprintf("%s: night difference within each group recovered (A %.2f h, B %.2f h)", label,
                ct$estimate[grepl("in A$", ct$effect)], ct$estimate[grepl("in B$", ct$effect)]),
        sprintf("%s: the comparison did not recover the night differences", label))
    for (o in c("dlmo_profile_select_ui", "dlmo_profile_plot", "dlmo_profile_note", "dlmo_summary",
                "dlmo_table", "dlmo_compare_ui", "dlmo_desc_table", "dlmo_test_note", "dlmo_test_table",
                "dlmo_agree_table", "dlmo_hockey_status", "dlmo_var_ui"))
      render(paste(label, o), output[[o]])
    session$setInputs(dlmo_profile_id = "_all_")
    render(paste(label, "plot, all participants"), output$dlmo_profile_plot)
    session$setInputs(dlmo_compare_method = "variable")
    render(paste(label, "comparison, variable threshold"), output$dlmo_test_table)
    f <- render(paste(label, "CSV"), output$dlmo_download)
    chk(!is.null(f) && nrow(read.csv(f)) == nrow(tab), sprintf("%s: the CSV holds every onset", label),
        sprintf("%s: the CSV is missing or short", label))

    # now with a bedtime column: the phase angle is bedtime - DLMO
    session$setInputs(dlmo_bed_src = "column", dlmo_bed_var = "bed", run_dlmo = 3)
    r <- values$dlmo_results
    pa <- r$tab[r$tab$method == "fixed", ]
    bh <- dance_dlmo_clock_hours(pa$bedtime_h)
    exp_pa <- (23 + c(0, 0.25, 0.5)[match(bed, sprintf("23:%02d", c(0, 15, 30)))])[match(pa$id, sprintf("P%02d", 1:n))] - pa$dlmo_h
    chk(isTRUE(r$has_bed) && max(abs(pa$phase_angle_h - exp_pa)) < 1e-8,
        sprintf("%s: phase angle = bedtime - DLMO for every profile", label),
        sprintf("%s: the phase angle is wrong", label))
    for (o in c("dlmo_compare_ui", "dlmo_pa_desc_table", "dlmo_pa_test_note", "dlmo_pa_test_table", "dlmo_table"))
      render(paste(label, o, "(bedtime)"), output[[o]])
    rep <- render(paste(label, "APA report"), dance_apa_report(values, input, "T"))
    chk(any(grepl("Dim-light melatonin onset", rep)) && any(grepl("Phase angle", rep)) &&
          any(grepl("Voultsios", rep)),
        sprintf("%s: the report has the DLMO section, with the phase angle and the sources", label),
        sprintf("%s: the report lacks the DLMO section", label))
    # the exported script reproduces the app's onsets
    code <- render(paste(label, "R script"), generate_analysis_code())
    if (!is.null(code)) {
      env <- new.env(parent = globalenv())
      env$dance_input_data <- values$data
      ok <- tryCatch({ eval(parse(text = code), envir = env); TRUE }, error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
      chk(ok && isTRUE(all.equal(env$dlmo_tab$dlmo_h, r$tab$dlmo_h)) &&
            isTRUE(all.equal(env$dlmo_tab$phase_angle_h, r$tab$phase_angle_h)),
          sprintf("%s: the exported script reproduces every onset and phase angle", label),
          sprintf("%s: the exported script does not reproduce the app", label))
    }
  })
}

# Layout 1: both nights in one row
d1 <- data.frame(ID = sprintf("P%02d", 1:n), group = grp, bed = bed, Y1, Y2)
names(d1)[-(1:3)] <- c(sprintf("Mel_%02d", clock), sprintf("Mel_%02d_N2", clock))
run_layout("one row per participant", d1, names(d1)[-(1:3)], "_none_", "mel_wide.csv")

# Layout 2: one row per night, named by a column. A missing baseline sample
# (18:00) changes no onset: the crossing is bracketed by observed samples.
d2 <- data.frame(ID = rep(sprintf("P%02d", 1:n), 2), group = rep(grp, 2), bed = rep(bed, 2),
                 night = rep(c("baseline", "SD"), each = n), rbind(Y1, Y2))
names(d2)[-(1:4)] <- sprintf("Mel_%02d", clock)
d2[3, "Mel_18"] <- NA
run_layout("one row per night", d2, names(d2)[-(1:4)], "night", "mel_long.csv")

# Column names whose first number is not the hour: two columns read as the
# same time, and the run is refused rather than run on wrong times
d3 <- d1; names(d3)[-(1:3)] <- c(sprintf("Mel_%02d", clock), sprintf("Mel_N2_%02d", clock))
csv3 <- tempfile(fileext = ".csv"); write.csv(d3, csv3, row.names = FALSE)
testServer(server, {
  session$setInputs(datafile = list(name = "bad.csv", datapath = csv3), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(d3)[-(1:3)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  session$setInputs(dlmo_time_src = "_shared_", dlmo_methods = "fixed", dlmo_threshold = 3, run_dlmo = 1)
  chk(isFALSE(values$dlmo_results$ok) && grepl("same time", values$dlmo_results$message),
      "times read twice from the column names are refused, and the user told to enter them",
      "a run went ahead on repeated times")
  # entered by hand, the same file runs
  session$setInputs(dlmo_time_src = "_manual_",
                    dlmo_manual_times = paste(c(sprintf("%02d:00", clock), sprintf("%02d:00", clock)), collapse = ","),
                    run_dlmo = 2)
  chk(isTRUE(values$dlmo_results$ok) && length(values$dlmo_results$profiles) == 2 * n,
      "with the times entered by hand, the two nights are found",
      sprintf("manual times failed: %s", values$dlmo_results$message %||% ""))
})

cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
