# ==============================================================================
# tests/no_mesor_app_test.R -- "MESOR" appears nowhere a reader sees
#
# The time-averaged midline level is the app's name for the quantity the
# classical cosinor calls the MESOR; the word is gone from the interface, the
# tables, the plots, the summary, the downloadable files and the report, and
# survives only in the documentation that says how the two relate. Internal
# keys ("mesor") are not shown: a menu's option VALUE is not text a reader sees,
# so tags and attributes are stripped before searching.
# Run with:   Rscript tests/no_mesor_app_test.R      (from the DANCE directory)
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

# visible text only: no tags, no option values, no attributes
visible <- function(x) {
  x <- paste(unlist(x), collapse = "\n")
  x <- gsub("<option[^>]*>", "<option>", x)
  x <- gsub("data-value=\"[^\"]*\"", "", x)
  gsub("<[^>]+>", " ", x)
}
hit <- function(x) {
  m <- regmatches(x, gregexpr("[^\n]{0,40}(?i:mesor)[^\n]{0,40}", x, perl = TRUE))[[1]]
  unique(m)
}

# the UI, as built
ui_src <- new.env()
for (f in sort(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$")))
  dance_source(f, ui_src)
ui_txt <- visible(as.character(ui_src$ui_tab_harmonic))
h <- hit(ui_txt)
chk(!length(h), "the cosinor tab's interface does not say MESOR",
    paste("the interface says:", paste(head(h, 5), collapse = " | ")))

set.seed(41)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 14)
tt <- seq(0, 30, by = 2); n <- 20; grp <- rep(c("ctrl", "sd"), each = n / 2)
Y <- t(vapply(seq_len(n), function(i) 3 + rnorm(1, 0, .6) + (2 + 1.5 * (grp[i] == "sd")) * (1 - exp(-tt / 14)) +
  1.2 * cos(2 * pi * (tt + 4) / 24) + rnorm(16, 0, .5), numeric(16)))
d <- data.frame(ID = sprintf("P%02d", seq_len(n)), group = grp, round(Y, 2))
names(d)[-(1:2)] <- c(paste0("KSS_", clock[1:12]), paste0("KSS_", clock[13:16], "_SD"))
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE)
server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(d)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  for (ap in c("two_stage", "mixed")) for (st in c("rhythm", "combined")) {
    k <- k + 1L; tag <- paste(ap, st)
    session$setInputs(harmonic_approach = ap, harmonic_structure = st, harmonic_trend_type = "linear",
                      harmonic_period = 24, n_harmonics = 1, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between", harmonic_between_var = "group",
                      harmonic_df_method = "satterthwaite", harmonic_traj_adjust = "holm",
                      harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                      harmonic_show_data = TRUE, run_harmonic = k)
    txt <- character(0)
    for (o in grep("^harmonic|^download_harmonic|^export_harmonic", names(outputOptions(output)), value = TRUE)) {
      v <- tryCatch(output[[o]], error = function(e) "")
      if (is.list(v) && !is.null(v$html)) v <- v$html
      if (is.character(v) && length(v) == 1L && file.exists(v)) v <- readLines(v, warn = FALSE)
      txt <- c(txt, visible(v))
    }
    for (w in (if (ap == "mixed") harmonic_traj_whats(values$harmonic_model$traj) else dance_ts_pair_whats(harmonic_ts()$mod))) {
      session$setInputs(harmonic_traj_what = w)
      txt <- c(txt, visible(tryCatch(output$harmonic_traj_pairwise, error = function(e) "")))
    }
    session$setInputs(harmonic_traj_what = "intercept")
    txt <- c(txt, visible(dance_apa_report(values, input, "T")))
    # the exported script's comments (its code -- the app's functions, written
    # out -- keeps the internal keys, which a reader does not see printed)
    code <- strsplit(generate_analysis_code(), "\n")[[1]]
    txt <- c(txt, grep("^\\s*#", code, value = TRUE))
    h <- hit(paste(txt, collapse = "\n"))
    chk(!length(h), sprintf("%s: no MESOR in any output, the report, the downloads or the script's text", tag),
        sprintf("%s: still shown: %s", tag, paste(head(h, 8), collapse = " | ")))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
