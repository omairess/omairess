# ==============================================================================
# tests/fanova_missing_group_app_test.R -- a curve with no group, through the app
#
# 30 curves, one without a group value. Before: the observed F skipped it, the
# permutations did not, every pointwise p came out missing (an empty p plot, no
# significant region) and the L2 p read 1. Now the curve is left out and the
# readout says so; the test has its p-values; the F plot's red line is labelled
# as the uncorrected criterion and the FDR-significant stretches are shaded;
# the report and the exported script leave the same curve out.
# Run with:   Rscript tests/fanova_missing_group_app_test.R   (from the DANCE directory)
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

set.seed(3)
hrs <- 8:22; n <- 30
grp <- c(rep("1", 12), rep("2", 17), NA)
Y <- t(sapply(seq_len(n), function(i) 1 + 2 * exp(-((hrs - 8) / 2)^2) +
  (if (isTRUE(grp[i] == "1")) 1.2 * exp(-((hrs - 16) / 1.2)^2) else 0) + rnorm(length(hrs), 0, 0.3)))
d <- data.frame(ID = sprintf("S%02d", 1:n), group = grp, round(Y, 3))
names(d)[-(1:2)] <- sprintf("V_%02dh", hrs)
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE, na = "")
server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())

suppressWarnings(testServer(server, {
  session$setInputs(datafile = list(name = "fa.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(d)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  session$setInputs(fanova_design = "between", fanova_group_var = "group", n_permutations = 999,
                    alpha_level = 0.05, fanova_test_type = "both", tick_freq_fanova = "2", run_fanova = 1)
  session$flushReact()
  res <- values$fanova_results
  chk(!is.null(res) && identical(res$n_dropped, 1L) && identical(as.numeric(res$group_sizes), c(12, 17)) &&
        res$df_within == 27,
      "the curve without a group is left out: groups of 12 and 17, df within 27",
      sprintf("dropped %s, sizes %s, df %s", format(res$n_dropped), paste(res$group_sizes, collapse = ","), format(res$df_within)))
  chk(!anyNA(res$p_values_pointwise) && is.finite(res$p_value_L2) && res$p_value_L2 < 0.05,
      sprintf("the test has its p-values (L2 p = %.4f)", res$p_value_L2),
      "pointwise or L2 p-values are missing")
  rd <- paste(output$fanova_global_results, collapse = "\n")
  chk(grepl("Curves left out: 1", rd, fixed = TRUE) && grepl("FDR-adjusted", rd, fixed = TRUE) &&
        grepl("UNCORRECTED", rd, fixed = TRUE),
      "the readout counts the left-out curve and says which criterion 'significant' uses",
      "the readout does not state the left-out curve or the criterion")
  fj <- jsonlite::fromJSON(output$fanova_fstat_plot, simplifyVector = FALSE)$x
  nms <- vapply(fj$data, function(t) t$name %||% "", "")
  shp <- Filter(function(s) identical(s$type, "rect"), fj$layout$shapes %||% list())
  chk(any(grepl("Uncorrected critical F", nms)) && length(shp) == sum(diff(c(FALSE, res$sig_regions)) == 1) &&
        any(res$sig_regions),
      sprintf("the F plot labels its red line as uncorrected and shades the %d FDR-significant stretch(es)", length(shp)),
      sprintf("F plot: line labels %s, %d shaded, %d significant points", paste(nms, collapse = "|"), length(shp), sum(res$sig_regions)))
  pj <- jsonlite::fromJSON(output$fanova_pvalue_plot, simplifyVector = FALSE)$x
  chk(length(pj$layout$annotations %||% list()) == 0 && any(vapply(pj$data, function(t) length(unlist(t$y)) > 2, TRUE)),
      "the p-value plot draws its p-values", "the p-value plot is empty")
  rep <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
  chk(grepl("1 curve(s) with a missing group value were left out", rep, fixed = TRUE),
      "the report says a curve was left out", "the report does not mention the left-out curve")
  code <- generate_analysis_code()
  env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
  ok <- tryCatch({ suppressMessages(capture.output(eval(parse(text = code), envir = env))); TRUE },
                 error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
  fr <- env$fanova_results
  if (ok && !is.null(fr)) cat(sprintf("   script: dropped %s, max |dF| %.3g, L2 p %s; fd coefs app %s script %s\n", format(fr$n_dropped),
      max(abs(fr$F_stat - res$F_stat)), format(fr$p_value_L2), paste(dim(res$fd_used$coefs), collapse = "x"),
      paste(dim(env$fd_obj$coefs), collapse = "x")))
  chk(ok && identical(fr$n_dropped, 1L) && isTRUE(all.equal(fr$F_stat, res$F_stat)) && is.finite(fr$p_value_L2),
      "the exported script leaves the same curve out and gives the same F",
      sprintf("the script differs (dropped %s)", format(fr$n_dropped)))
}))
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
