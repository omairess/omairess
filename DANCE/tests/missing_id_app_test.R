# ==============================================================================
# tests/missing_id_app_test.R -- a row whose participant identifier is empty
#
# A file with an empty cell in its 'subject' column (row 22 of a real one) showed
# that participant as "Row 22" in the picker and passed NA on as its identifier.
# The mixed models then dropped every observation of that participant (lme4
# drops rows whose grouping factor is NA) while the other tabs kept the curve.
# The import now names such a row by its data row, says so, and every analysis
# sees the same participants.
#
# Run with:   Rscript tests/missing_id_app_test.R      (from the DANCE directory)
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

set.seed(23)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22)
n <- 24; grp <- rep(1:2, length.out = n)
Y <- t(vapply(seq_len(n), function(i)
  5 + rnorm(1, 0, 1) + 0.8 * grp[i] + 2 * cos(2 * pi * (clock - 16) / 24) + rnorm(8, 0, 0.7), numeric(8)))
d <- data.frame(subject = seq_len(n), group = grp, round(Y, 2))
names(d)[-(1:2)] <- paste0("V_", clock)
d$subject[22] <- NA                                   # the empty cell
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE, na = "")

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
testServer(server, {
  session$setInputs(datafile = list(name = "ids.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group",
                    subject_id_var = "subject", apply_selection = 1)
  chk(!anyNA(values$subject_ids) && identical(values$subject_ids[22], "Row 22 (no ID)") &&
        length(unique(values$subject_ids)) == n,
      "the empty identifier becomes 'Row 22 (no ID)', one participant of its own",
      sprintf("ids: %s", paste(values$subject_ids, collapse = ",")))
  st <- paste(capture.output(output$data_status), collapse = " ")
  chk(grepl("1 row(s) have an empty 'subject' cell (data row 22", st, fixed = TRUE),
      "the data status says which row has no identifier", sprintf("status: %s", substr(st, 1, 300)))
  lab <- values$import_full$label[22]
  chk(identical(lab, "Row 22 (no ID)"), "the participant picker shows the same label",
      sprintf("picker label: %s", lab))

  session$setInputs(harmonic_approach = "mixed", harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = "none", harmonic_time_var = "_shared_",
                    harmonic_design = "between", harmonic_between_var = "group",
                    harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                    harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                    harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_traj_tau = 18,
                    harmonic_traj_tau_mode = "fixed", run_harmonic = 1)
  ff <- values$harmonic_model$traj
  if (!isTRUE(ff$ok)) {
    cat("traj:", ff$message %||% "NULL", "| clock:", paste(values$time_clock, collapse = ","), "\n")
    print(tryCatch(output$harmonic_summary, error = function(e) conditionMessage(e)))
  }
  # the participants lme4 actually used: an NA identifier is dropped by the
  # model frame, so counting the labels in the data would not see it
  np <- if (isTRUE(ff$ok)) unname(lme4::ngrps(ff$model)[["subject"]]) else NA
  chk(identical(as.integer(np), as.integer(n)),
      sprintf("the mixed cosinor keeps all %d participants", n),
      sprintf("the mixed cosinor has %s participants", np))
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
