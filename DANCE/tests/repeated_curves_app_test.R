# ==============================================================================
# tests/repeated_curves_app_test.R -- several curves per participant, through
# import and the harmonic module
#
# A file shaped like a repeated-cycles study: a participant column, a cycle
# column (each participant measured in 3 to 7 cycles), two-hourly columns from
# 8h to 10h the next day with the 8h and 10h names repeated, and many missing
# values. Checks that
#   * the data status counts curves AND participants;
#   * the design readout reads 'cycle' as within-participant from the
#     participant column (it used row names, so every curve was a participant);
#   * the two-stage run completes with two harmonics when a curve is flat (its
#     NULL confidence regions deleted their slots: "subscript out of bounds");
#   * the mixed run computes the cycle tests and the pairwise table.
#
# Run with:   Rscript tests/repeated_curves_app_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({ library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "shinyWidgets")) if (requireNamespace(p, quietly = TRUE)) library(p, character.only = TRUE) })
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  pickerInput <- function(inputId, label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(), ...)
    selectInput(inputId, label, choices = choices, selected = selected, multiple = multiple)
  updatePickerInput <- function(session, inputId, label = NULL, selected = NULL, choices = NULL, ...)
    updateSelectInput(session, inputId, label = label, choices = choices, selected = selected)
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }
txt <- function(x) gsub("\\s+", " ", gsub("<[^>]+>", " ", paste(unlist(x), collapse = " ")))

# ---- a repeated-cycles file -------------------------------------------------
set.seed(11)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10)
tt <- clock; for (i in 2:length(tt)) if (tt[i] < tt[i - 1]) tt[i:length(tt)] <- tt[i:length(tt)] + 24
ids <- sprintf("P%02d", 1:13)
rows <- do.call(rbind, lapply(ids, function(p) {
  cyc <- sort(sample(1:7, sample(3:7, 1)))
  u <- rnorm(1, 0, 2); a <- rnorm(1, 4, 1)
  do.call(rbind, lapply(cyc, function(k) {
    y <- 8 + u + 0.3 * k + a * cos(2 * pi * (tt - 16) / 24) + rnorm(length(tt), 0, 1.5)
    y[sample(length(y), sample(2:6, 1))] <- NA
    data.frame(cycle = k, participant = p, t(round(y, 1)), check.names = FALSE)
  }))
}))
names(rows)[-(1:2)] <- paste0(clock, "h")          # 8h and 10h appear twice
# a participant at the floor of the scale: 1 at every time point, the night
# missing -- zero amplitude, so no harmonic has a confidence region (the old
# code then deleted the slots and the run failed: "subscript out of bounds")
rows[1, 3:16] <- c(rep(1, 8), NA, NA, NA, NA, 1, 1)
csv <- tempfile(fileext = ".csv"); write.csv(rows, csv, row.names = FALSE, na = "")
n_curves <- nrow(rows)

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
run_harmonic <- function(session, approach, id) {
  session$setInputs(harmonic_approach = approach, harmonic_period = 24, n_harmonics = 2,
                    harmonic_trend_type = "none", harmonic_time_var = "_shared_",
                    harmonic_design = "within", harmonic_within_var = "cycle",
                    harmonic_df_method = "kr", harmonic_traj_adjust = "holm",
                    harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                    harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_traj_tau = 18,
                    harmonic_traj_tau_mode = "fixed", run_harmonic = id)
}
testServer(server, {
  session$setInputs(datafile = list(name = "cycles.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "cycle",
                    subject_id_var = "_none_", apply_selection = 1)
  st <- paste(capture.output(output$data_status), collapse = " ")
  chk(grepl(sprintf("%d curves x 14 time points", n_curves), st, fixed = TRUE) &&
        grepl("Participants: 13 -- repeated measures", st, fixed = TRUE),
      "the data status counts curves and participants", sprintf("the data status: %s", substr(st, 1, 200)))

  run_harmonic(session, "two_stage", 1)
  rd <- txt(output$harmonic_design_readout$html)
  chk(grepl("the data look within", rd, fixed = TRUE) &&
        grepl(sprintf("%d curves from 13 participants", n_curves), rd, fixed = TRUE),
      "the design readout reads cycle as within, from the participant column",
      sprintf("the design readout: %s", substr(rd, 1, 200)))
  chk(isTRUE(values$harmonic_model$ok %||% !is.null(values$harmonic_model$individual_fits)),
      "two-stage runs with two harmonics on sparse curves", "the two-stage run failed")
  for (o in c("harmonic_traj_primary", "harmonic_traj_pairwise")) {
    x <- tryCatch(txt(output[[o]]), error = function(e) paste("Error:", conditionMessage(e)))
    chk(!grepl("Error", x) && grepl("p", x), sprintf("two-stage %s renders", o), sprintf("two-stage %s: %s", o, substr(x, 1, 160)))
  }

  run_harmonic(session, "mixed", 2)
  ff <- values$harmonic_model$traj
  chk(isTRUE(ff$ok) && identical(as.character(ff$spec$roles[["cycle"]] %||% "within"), "within"),
      "the mixed model fits cycle as a within-participant factor", sprintf("mixed fit: %s", ff$message %||% ""))
  x <- txt(output$harmonic_traj_primary)
  chk(grepl("Kenward-Roger", x) && !grepl("Error", x), "the mixed cycle test is computed (Kenward-Roger)",
      sprintf("mixed primary: %s", substr(x, 1, 160)))
  x <- txt(output$harmonic_traj_pairwise)
  chk(grepl("1 vs 2", x, fixed = TRUE) && !grepl("Error", x), "the mixed pairwise table renders",
      sprintf("mixed pairwise: %s", substr(x, 1, 160)))
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
