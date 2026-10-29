# ==============================================================================
# tests/excel_time_names_app_test.R -- time-formatted header cells
#
# A sheet whose header row holds 8:00, 10:00, ..., 0:00, ..., 16:00 (next day)
# formatted as times reaches the app through readxl as fractions of a day,
# "0.33333333333333331...5" and so on. Imported with exactly those names, the
# app must read them as the clock times they display.
#
# Run with:   Rscript tests/excel_time_names_app_test.R      (from the DANCE directory)
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

xl <- c("0.33333333333333331...3", "0.41666666666666669...4", "0.5...5", "0.58333333333333337",
        "0.66666666666666663...7", "0.75", "0.83333333333333337", "0.91666666666666663", "0",
        "8.3333333333333329E-2", "0.16666666666666666", "0.25", "0.33333333333333331...15",
        "0.41666666666666669...16", "0.5...17", "0.66666666666666663...18")
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 16)
set.seed(5)
d <- data.frame(subject = 1:12, group = rep(1:2, 6),
                t(vapply(1:12, function(i) 5 + 2 * cos(2 * pi * (clock - 16) / 24) + rnorm(16, 0, .5), numeric(16))))
names(d)[-(1:2)] <- xl
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE)

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
testServer(server, {
  session$setInputs(datafile = list(name = "times.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  chk(identical(names(values$raw_df)[3:18],
                c("08:00", "10:00", "12:00", "14:00", "16:00", "18:00", "20:00", "22:00", "00:00",
                  "02:00", "04:00", "06:00", "08:00 #2", "10:00 #2", "12:00 #2", "16:00 #2")),
      "the day-fraction names are read as the clock times they display",
      paste("names:", paste(names(values$raw_df), collapse = ", ")))
  session$setInputs(sel_data_vars = names(values$raw_df)[3:18], sel_cov_vars = "group",
                    subject_id_var = "subject", apply_selection = 1)
  chk(isTRUE(all.equal(values$time_clock, clock)), "the clock times are parsed at import",
      paste("time_clock:", paste(values$time_clock, collapse = ",")))
  chk(isTRUE(all.equal(dance_cumulative_hours(values$time_labels), c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32))),
      "the recording runs 32 h, across midnight", "the elapsed hours are wrong")
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
