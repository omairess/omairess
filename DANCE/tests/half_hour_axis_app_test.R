# ==============================================================================
# tests/half_hour_axis_app_test.R -- half-hourly columns and a 20 h gap
#
# A protocol with one evening every 30 min (MT_20h00 ... MT_00h00), then the
# next evening into a night of sleep deprivation (MT_20h00_SD ... MT_07h00_SD,
# MT_16h00_SD). Three things went wrong on such a file:
#   * the manual-times suggestion took the MINUTES: 0, 30, 0, 30, ...
#   * the data plot drew the curves by column index but put the hour labels at
#     real-time positions, so the labels did not sit under the data
#   * (the hover text interpolated clock hours across midnight)
#
# Run with:   Rscript tests/half_hour_axis_app_test.R      (from the DANCE directory)
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

ev1 <- c(20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 0)
ev2 <- c(20, 21, 21.5, 22, 22.5, 23, 23.5, 0, 1, 2, 3, 4, 5, 6, 7, 16)
hm <- function(h) sprintf("%02dh%02d", floor(h), round((h - floor(h)) * 60))
cols <- c(paste0("MT_", hm(ev1)), paste0("MT_", hm(ev2), "_SD"))
clock <- c(ev1, ev2)
elapsed <- c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 24, 25, 25.5, 26, 26.5, 27, 27.5, 28, 29, 30, 31, 32, 33, 34, 35, 44)
set.seed(8)
n <- 16
Y <- t(vapply(1:n, function(i) 5 + 4 * cos(2 * pi * ((20 + elapsed) - 3) / 24) + rnorm(25, 0, 1), numeric(25)))
d <- data.frame(subject = 1:n, group = rep(1:2, n / 2), round(Y, 2))
names(d)[-(1:2)] <- cols
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE)

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
testServer(server, {
  session$setInputs(datafile = list(name = "mt.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = cols, sel_cov_vars = "group", subject_id_var = "subject",
                    apply_selection = 1)
  chk(isTRUE(all.equal(values$time_clock, clock)), "the import reads every half hour",
      paste("time_clock:", paste(values$time_clock, collapse = ",")))
  chk(isTRUE(all.equal(dance_cumulative_hours(values$time_labels), elapsed)),
      "elapsed time: 4 h, a 20 h gap, then 20 h more (44 h)", "the elapsed hours are wrong")

  ui <- paste(as.character(output$harmonic_var_select_ui$html), collapse = " ")
  chk(grepl(paste(clock, collapse = ","), ui, fixed = TRUE) && !grepl("0, 30, 0, 30", ui, fixed = TRUE),
      "the manual-times suggestion is the clock times, not the minutes",
      sprintf("suggestion: %s", regmatches(ui, regexpr("Detected[^<]*", ui))))

  # the data plot: the hour labels sit where the columns are drawn
  pj <- jsonlite::fromJSON(output$data_plot, simplifyVector = FALSE)
  xa <- pj$x$layout$xaxis
  mean_tr <- Filter(function(tr) identical(tr$name, "Mean"), pj$x$data)[[1]]
  x_data <- unlist(mean_tr$x)
  tv <- unlist(xa$tickvals); tt <- unlist(xa$ticktext)
  # "Every hour": a label on each whole-hour column, at that column's x
  keep <- abs(clock - round(clock)) < 1e-6
  lab <- vapply(clock, function(h) sprintf("%dh%02d", floor(h), round((h - floor(h)) * 60)), "")
  chk(isTRUE(all.equal(tv, x_data[keep])) && identical(tt, lab[keep]),
      "every hour label sits on the column it names", "the hour labels are not under their columns")

  # the cosinor on the import's clock times: 44 h, and the APA report carries
  # the pairwise tables of both approaches
  for (ap in c("two_stage", "mixed")) {
    session$setInputs(harmonic_approach = ap, harmonic_period = 24, n_harmonics = 1,
                      harmonic_trend_type = "none", harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between",
                      harmonic_between_var = "group", harmonic_df_method = "satterthwaite",
                      harmonic_traj_adjust = "holm", harmonic_traj_what = "intercept",
                      harmonic_traj_resid = "independent", harmonic_tau_lo = 2,
                      harmonic_tau_hi = 48, harmonic_traj_tau = 18, harmonic_traj_tau_mode = "fixed",
                      run_harmonic = if (ap == "mixed") 2 else 1)
    tv <- values$harmonic_model$time_vec %||% values$harmonic_model$traj$spec$data$t
    chk(isTRUE(abs(diff(range(tv)) - 44) < 1e-8), sprintf("%s: the cosinor's times span 44 h", ap),
        sprintf("%s: time span %s", ap, paste(range(tv), collapse = " to ")))
    rp <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
    pw <- harmonic_pairwise_all()
    chk(grepl("**Pairwise comparisons.**", rp, fixed = TRUE) && !is.null(pw) &&
          grepl("| 1 vs 2 |", rp, fixed = TRUE),
        sprintf("%s: the APA report has the pairwise table (%d rows)", ap, if (is.null(pw)) 0L else nrow(pw)),
        sprintf("%s: no pairwise table in the APA report", ap))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
