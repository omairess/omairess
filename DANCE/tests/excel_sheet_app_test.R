# ==============================================================================
# tests/excel_sheet_app_test.R -- choosing the sheet (tab) of an Excel file
#
# The sheet picker and the separator used to sit in conditionalPanels whose
# condition read input.datafile.name in the browser, where a file input has no
# such value: neither was ever shown, and every Excel file was read from its
# first sheet. They are now rendered by the server. Checked here:
#   * the picker lists every sheet and the preview describes the chosen one;
#   * loading reads the chosen sheet, not the first;
#   * the file and sheet are recorded, and named in the report and the script;
#   * no condition in the UI reads input.datafile again.
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
if (!requireNamespace("writexl", quietly = TRUE) || !requireNamespace("readxl", quietly = TRUE)) {
  cat("SKIP: writexl and readxl are needed to build and read the test workbook\n"); quit(status = 0)
}

tt <- sprintf("%02d:00", seq(8, 22, 2))
mk <- function(n, off) { set.seed(n); m <- round(matrix(rnorm(n * length(tt), off, 1), n), 2)
  d <- data.frame(ID = sprintf("P%02d", 1:n), group = rep(c("a", "b"), length.out = n), m); names(d)[-(1:2)] <- tt; d }
xl <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(list(Notes = data.frame(info = "the data are on sheet Ratings"),
                         Ratings = mk(12, 5), Other = mk(6, 50)), xl)

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
testServer(server, {
  session$setInputs(datafile = list(name = "study.xlsx", datapath = xl), header = TRUE, dec = "auto")
  ui <- output$excel_sheet_selector
  chk(all(vapply(c("Notes", "Ratings", "Other"), grepl, logical(1), x = ui$html, fixed = TRUE)),
      "the picker lists every sheet", "the picker does not list the sheets")
  session$setInputs(excel_sheet = "Ratings")
  chk(grepl("12 rows, 10 columns: ID, group, 08:00", output$excel_sheet_preview$html, fixed = TRUE),
      "the preview describes the chosen sheet", "the preview is wrong")
  session$setInputs(data_format = "wide", load_data = 1)
  chk(identical(dim(values$raw_df), c(12L, 10L)) && identical(names(values$raw_df)[3], "08:00"),
      "loading reads the chosen sheet, not the first", sprintf("loaded %s", paste(dim(values$raw_df), collapse = " x ")))
  chk(identical(values$data_source$sheet, "Ratings") && identical(values$data_source$file, "study.xlsx"),
      "the file and sheet are recorded", "the data source was not recorded")
  session$setInputs(excel_sheet = "Other", load_data = 2)
  chk(nrow(values$raw_df) == 6L && mean(as.matrix(values$raw_df[, -(1:2)])) > 40,
      "another sheet can be chosen and loaded", "switching sheets did not load the other sheet")
  session$setInputs(excel_sheet = "Ratings", load_data = 3)
  session$setInputs(sel_data_vars = tt, sel_cov_vars = "group", subject_id_var = "ID", apply_selection = 1)
  rep <- paste(dance_apa_report(values, input), collapse = "\n")
  chk(grepl("It was read from study.xlsx, sheet \"Ratings\".", rep, fixed = TRUE),
      "the report names the file and sheet", "the report does not name the sheet")
  code <- paste(generate_analysis_code(full = TRUE), collapse = "\n")
  chk(grepl("# Data: study.xlsx, sheet \"Ratings\"", code, fixed = TRUE),
      "the script header names the file and sheet", "the script header does not name the sheet")

  # a text file gets the separator instead
  csv <- tempfile(fileext = ".csv"); write.csv(mk(4, 1), csv, row.names = FALSE)
  session$setInputs(datafile = list(name = "x.csv", datapath = csv))
  chk(grepl("name=\"sep\"", output$excel_sheet_selector$html, fixed = TRUE),
      "a text file is offered the separator", "the separator is not offered for a text file")
})

ui_src <- paste(unlist(lapply(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"), readLines)), collapse = "\n")
chk(!grepl("input.datafile", ui_src, fixed = TRUE),
    "no browser-side condition reads the file input", "a conditionalPanel still reads input.datafile")
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
