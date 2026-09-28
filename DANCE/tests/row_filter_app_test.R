# ==============================================================================
# tests/row_filter_app_test.R — the participant/group selection, driven through
# the Data Import tab's own observers
#
# tests/testthat/test-row-filter.R checks the pure rule. This presses the
# buttons: confirm a variable selection, untick a group and a participant,
# apply, re-confirm, try a selection that leaves too little, include everyone
# again, and do the same on the sample data -- checking at every step that the
# analysis frame, its parallel vectors and the cleared results agree.
#
# Run with:   Rscript tests/row_filter_app_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "mgcv", "ggplot2")) if (requireNamespace(p, quietly = TRUE))
    library(p, character.only = TRUE)
})
# The pickers are shinyWidgets widgets. Where that package is not installed a
# multi-select stands in: the server logic under test reads the same input
# values either way.
if (requireNamespace("shinyWidgets", quietly = TRUE)) {
  suppressPackageStartupMessages(library(shinyWidgets))
} else {
  pickerInput <- function(inputId, label = NULL, choices = NULL, selected = NULL,
                          multiple = FALSE, options = list(), ...)
    selectInput(inputId, label, choices = choices, selected = selected, multiple = multiple)
  updatePickerInput <- function(session, inputId, label = NULL, selected = NULL,
                                choices = NULL, ...)
    updateSelectInput(session, inputId, label = label, choices = choices, selected = selected)
}

app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))

failures <- 0L
fail <- function(...) { cat("FAIL:", ..., "\n"); failures <<- failures + 1L }
ok   <- function(...) cat("ok  :", ..., "\n")
chk  <- function(cond, good, bad) if (isTRUE(cond)) ok(good) else fail(bad)

# ---- a raw file as the loader leaves it -------------------------------------
# 30 rows from 24 participants: six of them contribute a second session. Row 7
# is entirely missing and is dropped at import, so file rows and curve
# positions disagree from there on -- which is the case a positional record
# gets wrong. "Visit" is in the file but never selected as a scalar variable.
set.seed(11)
pid  <- c(sprintf("P%02d", 1:24), sprintf("P%02d", c(5, 10, 15, 20, 22, 23)))
num  <- as.integer(sub("P", "", pid))
age  <- c("YOUTH", "ADULT", "ELDERLY")[(num - 1) %% 3 + 1]
sex  <- c("F", "M")[(num %% 2) + 1]
hrs  <- 8:19
curves <- t(vapply(seq_along(pid), function(i)
  50 + 8 * cos(2 * pi * (hrs - 15 - (age[i] == "ELDERLY")) / 24) + rnorm(length(hrs), 0, 1.5),
  numeric(length(hrs))))
curves[7, ] <- NA
raw <- data.frame(ID = pid, AGE = age, Sex = sex, Score = round(rnorm(30, 100, 15), 1),
                  Visit = ifelse(duplicated(pid), "second", "first"),
                  stringsAsFactors = FALSE)
time_cols <- sprintf("%02d:00", hrs)
for (j in seq_along(hrs)) raw[[time_cols[j]]] <- curves[, j]

server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$raw_df <- raw
  values$uploaded_data <- raw
}

testServer(server, {
  cat("\n-- confirm the variable selection -----------------------------------\n")
  session$setInputs(sel_data_vars = time_cols, sel_cov_vars = c("AGE", "Sex"),
                    subject_id_var = "ID", apply_selection = 1)
  session$flushReact()
  fr <- values$import_full
  chk(!is.null(fr) && nrow(fr$data) == 29L,
      "the full frame holds the 29 curves left after the all-missing row",
      sprintf("the full frame was not kept (%s curves)", if (is.null(fr)) "no" else nrow(fr$data)))
  chk(identical(values$row_index, c(1:6, 8:30)),
      "each curve knows its row in the file, across the dropped row",
      paste("row_index is wrong:", paste(head(values$row_index, 10), collapse = ",")))
  chk(is.null(values$row_filter) && nrow(values$data) == 29L,
      "nothing is excluded until the analyst excludes it",
      "rows were excluded without a selection")

  ui <- tryCatch(output$row_filter_ui, error = function(e) e)
  if (inherits(ui, "error")) fail("the selection box does not render:", conditionMessage(ui))
  else {
    html <- paste(as.character(ui$html), collapse = "")
    chk(grepl("filter_levels", html) && grepl("filter_participants", html) &&
          grepl("AGE::ELDERLY", html, fixed = TRUE) && grepl("id:P05", html, fixed = TRUE),
        "the box lists the groups (per variable) and the participants",
        "the box does not list the groups and participants")
    chk(!grepl("Score", html, fixed = TRUE),
        "a variable that was not selected is not offered",
        "an unselected variable is offered as a group")
  }

  cat("\n-- untick a group and a participant ---------------------------------\n")
  vars  <- dance_filter_vars(fr)
  codes <- unlist(lapply(names(vars), function(nm) dance_filter_code(nm, vars[[nm]]$level)),
                  use.names = FALSE)
  keys  <- dance_filter_participants(fr)$key
  # a placeholder result, to see it cleared
  values$harmonic_model <- list(placeholder = TRUE)
  session$setInputs(filter_levels = setdiff(codes, "AGE::ELDERLY"),
                    filter_participants = setdiff(keys, "id:P05"))
  prev <- tryCatch(output$row_filter_preview, error = function(e) conditionMessage(e))
  chk(grepl("would be analysed", prev, fixed = TRUE) && grepl("Applied:", prev, fixed = TRUE),
      "the preview shows the consequence before anything is cleared",
      paste("the preview does not show the pending selection:", prev))
  chk(nrow(values$data) == 29L, "ticking alone changes nothing",
      "ticking changed the analysed data before Apply")

  session$setInputs(apply_row_filter = 1)
  session$flushReact()
  expect_ids <- fr$subject_ids[fr$covariates$AGE != "ELDERLY" & fr$subject_ids != "P05"]
  chk(identical(values$subject_ids, expect_ids),
      sprintf("%d curves remain, the right ones (ELDERLY and both P05 sessions gone)",
              length(expect_ids)),
      "the analysed curves are not the ones the selection names")
  n_now <- nrow(values$data)
  aligned <- all(c(length(values$subject_ids), nrow(values$covariates),
                   nrow(values$group_variables), length(values$group_labels),
                   length(values$row_index)) == n_now)
  chk(aligned, "every parallel vector was cut the same way",
      "a parallel vector is out of step with the curves")
  chk(identical(levels(values$group_labels), c("ADULT", "YOUTH")) &&
        !("ELDERLY" %in% levels(values$group_variables$AGE)),
      "the excluded group is gone from the factors, not left as an empty level",
      paste("group levels after exclusion:", paste(levels(values$group_labels), collapse = ",")))
  raw_rows <- which(raw$AGE != "ELDERLY" & raw$ID != "P05" & seq_len(30) != 7)
  chk(identical(values$row_index, raw_rows) &&
        isTRUE(all.equal(unname(values$data), unname(as.matrix(raw[raw_rows, time_cols])))),
      "the curves are the file's own rows, value for value",
      "the analysed curves are not the file rows they claim to be")
  chk(is.null(values$harmonic_model), "every result was cleared by the change",
      "a result computed on the old selection survived it")
  chk(identical(values$row_filter$excluded_ids, "id:P05") &&
        identical(values$row_filter$excluded_levels, list(AGE = "ELDERLY")),
      "the selection is stored as exclusions, by participant and level",
      "the stored selection is not the keyed exclusion record")
  st <- paste(tryCatch(output$data_status, error = function(e) conditionMessage(e)),
              collapse = "\n")
  chk(grepl("Participant selection", st, fixed = TRUE),
      "the status panel says a selection is in force",
      "the status panel does not mention the selection")

  # a module reading a column straight from the raw file gets the analysed rows
  vis <- tryCatch(dance_rm_column(values, "Visit"), error = function(e) NULL)
  chk(identical(as.character(vis), raw$Visit[raw_rows]),
      "a raw-only column is read at the analysed rows (repeated-measures pickers)",
      "a raw-only column is misaligned with the analysed curves")

  cat("\n-- re-confirm with another covariate ---------------------------------\n")
  session$setInputs(sel_cov_vars = c("AGE", "Sex", "Score"), apply_selection = 2)
  session$flushReact()
  chk(identical(values$subject_ids, expect_ids) && "Score" %in% names(values$covariates),
      "the exclusions are re-applied to the same people after a re-confirm",
      "a re-confirm brought excluded people back or took others out")

  cat("\n-- a selection that leaves too little is refused ---------------------\n")
  before <- values$subject_ids
  session$setInputs(filter_participants = keys[1], apply_row_filter = 2)
  session$flushReact()
  chk(identical(values$subject_ids, before),
      "a selection leaving one curve is refused and nothing changes",
      "a selection leaving one curve was applied")

  cat("\n-- the report and the exported script say what was excluded ----------\n")
  md <- tryCatch(paste(dance_apa_report(values, input, "T"), collapse = "\n"),
                 error = function(e) paste("ERROR", conditionMessage(e)))
  chk(grepl("were excluded by the analyst", md, fixed = TRUE) &&
        grepl("AGE = ELDERLY", md, fixed = TRUE),
      "the report states the exclusions and their criterion",
      paste("the report does not state the exclusions:", substr(md, 1, 200)))
  code <- tryCatch(generate_analysis_code(full = TRUE), error = function(e) paste("ERROR", conditionMessage(e)))
  chk(grepl("analysed_rows <- c(", code, fixed = TRUE) &&
        grepl(sprintf("if (nrow(data_matrix) == %d)", nrow(raw)), code, fixed = TRUE),
      "the exported script carries the analysed rows",
      "the exported script does not carry the selection")
  pc <- tryCatch(parse(text = code), error = function(e) e)
  chk(!inherits(pc, "error"), "the exported script still parses",
      paste("the exported script does not parse:", if (inherits(pc, "error")) conditionMessage(pc)))

  cat("\n-- include everyone again ------------------------------------------------\n")
  session$setInputs(reset_row_filter = 1)
  session$flushReact()
  chk(is.null(values$row_filter) && nrow(values$data) == 29L &&
        identical(levels(values$group_labels), c("ADULT", "ELDERLY", "YOUTH")),
      "'Include everyone' restores all 29 curves and every level",
      "'Include everyone' did not restore the full frame")

  cat("\n-- the same on the sample data -------------------------------------------\n")
  session$setInputs(generate_with_groups = TRUE, n_groups = 3, generate_sample = 1)
  session$flushReact()
  fs <- values$import_full
  chk(!is.null(fs) && nrow(fs$data) == 50L && is.null(values$row_filter) &&
        is.null(values$subject_ids),
      "the sample data gets its own frame, with nothing carried over",
      "the sample data inherited the previous file's selection or identifiers")
  vs <- dance_filter_vars(fs)
  chk(all(c("Group", "Sex") %in% names(vs)) && !any(c("ID", "Age", "Score") %in% names(vs)),
      "its groups are offered (Group, Sex, Outcome), its continuous variables are not",
      paste("sample-data groups offered:", paste(names(vs), collapse = ", ")))
  cs <- unlist(lapply(names(vs), function(nm) dance_filter_code(nm, vs[[nm]]$level)), use.names = FALSE)
  session$setInputs(filter_levels = setdiff(cs, "Group::Group3"),
                    filter_participants = dance_filter_participants(fs)$key,
                    apply_row_filter = 3)
  session$flushReact()
  n_g3 <- sum(fs$group_labels == "Group3")
  chk(nrow(values$data) == 50L - n_g3 && !("Group3" %in% levels(values$group_labels)) &&
        nrow(values$covariates) == 50L - n_g3,
      sprintf("excluding Group3 leaves %d sample curves, aligned", 50L - n_g3),
      "excluding a sample group did not cut the data consistently")
  chk(identical(rownames(values$data), as.character(values$row_index)),
      "sample curves keep their row numbers after an exclusion",
      "sample curves were renumbered after an exclusion")
})

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nParticipant/group selection tests passed.\n")
