# ==============================================================================
# tests/p19_import_cyclic_app_test.R -- P19.1, P19.3 and P19.4 through the app
#
# Drives the real server: a semicolon file with decimal commas and one flagged
# cell is loaded, the curves are selected and smoothed with the cyclic option,
# and the values the analyses would see are checked.
#   * the matrix holds the numbers in the file, not factor codes;
#   * the flagged cell is missing, the rest of its column is kept;
#   * "Point" as the decimal mark refuses columns with no readable number;
#   * cyclic smoothing of 24 hourly columns uses a 24-column cycle, reproduces
#     an exactly 24-hour signal, and the fd object is that smoothed curve.
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

tt <- 0:23
truth <- t(vapply(1:10, function(i)
  round(10 + i / 10 + 5 * cos(2 * pi * (tt - 3) / 24) + 2 * sin(2 * pi * tt / 24), 3), numeric(24)))
txt <- matrix(sub(".", ",", formatC(truth, format = "f", digits = 3), fixed = TRUE), nrow = 10)
txt[4, 7] <- "<LOQ"
d <- data.frame(subject = sprintf("S%02d", 1:10), group = rep(c("a", "b"), 5), txt, stringsAsFactors = FALSE)
names(d)[-(1:2)] <- sprintf("%02d:00", tt)
csv <- tempfile(fileext = ".csv")
write.table(d, csv, sep = ";", row.names = FALSE, quote = FALSE)

server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
testServer(server, {
  cols <- sprintf("%02d:00", tt)

  # ---- decimal mark "Point": nothing is repaired, and nothing is invented ----
  session$setInputs(datafile = list(name = "comma.csv", datapath = csv), header = TRUE, sep = ";",
                    dec = "point", data_format = "wide", load_data = 1)
  chk(is.character(values$raw_df[["00:00"]]), "with 'Point' the comma columns stay text",
      "the columns were converted although 'Point' was chosen")
  session$setInputs(sel_data_vars = cols, sel_cov_vars = "group", subject_id_var = "subject",
                    apply_selection = 1)
  chk(is.null(values$data), "columns with no readable number are refused, not coded",
      sprintf("a matrix was built from text columns; its first cells: %s",
              paste(utils::head(as.vector(values$data), 4), collapse = ", ")))

  # ---- auto-detect ----------------------------------------------------------
  session$setInputs(dec = "auto", load_data = 2)
  chk(is.numeric(values$raw_df[["00:00"]]) && is.character(values$raw_df[["06:00"]]),
      "decimal-comma columns are read as numbers; the column with a flag stays text",
      "the decimal-comma repair did not behave as specified")
  session$setInputs(sel_data_vars = cols, sel_cov_vars = "group", subject_id_var = "subject",
                    apply_selection = 2)
  want <- truth; want[4, 7] <- NA
  chk(!is.null(values$data) && isTRUE(all.equal(unname(values$data), want)),
      "the matrix holds the numbers in the file; the flagged cell is missing, its column is kept",
      sprintf("the matrix is wrong; column 06:00 holds: %s",
              paste(utils::head(values$data[, 7], 5), collapse = ", ")))
})

# ---- a clean decimal-comma file, then cyclic smoothing -----------------------
d2 <- d; d2[4, "06:00"] <- sub(".", ",", formatC(truth[4, 7], format = "f", digits = 3), fixed = TRUE)
csv2 <- tempfile(fileext = ".csv")
write.table(d2, csv2, sep = ";", row.names = FALSE, quote = FALSE)
testServer(server, {
  cols <- sprintf("%02d:00", tt)
  session$setInputs(datafile = list(name = "comma2.csv", datapath = csv2), header = TRUE, sep = ";",
                    dec = "auto", data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = cols, sel_cov_vars = "group", subject_id_var = "subject",
                    apply_selection = 1)
  chk(isTRUE(all.equal(unname(values$data), truth)), "all 240 cells equal the file",
      "the clean decimal-comma file was not read exactly")

  for (rt in c(FALSE, TRUE)) {
    session$setInputs(is_cyclic = TRUE, use_real_time = rt, smooth_method = "manual",
                      n_basis_manual = 13, smooth_factor = 8, constrain_bounds = FALSE,
                      allow_extrapolation = FALSE, apply_smooth = as.integer(rt) + 1L)
    lab <- if (rt) "real clock time" else "column index"
    chk(isTRUE(all.equal(values$smooth_fit_metrics$cyclic_period, 24)),
        sprintf("%s: the cycle is 24 (%s)", lab, values$smooth_fit_metrics$cyclic_period_source),
        sprintf("%s: the cycle is %s", lab, format(values$smooth_fit_metrics$cyclic_period)))
    err <- max(abs(values$smooth_data - truth))
    chk(is.finite(err) && err < 0.01,
        sprintf("%s: an exactly 24-hour signal is reproduced (max error %.2g)", lab, err),
        sprintf("%s: the smooth misses a 24-hour signal by %.3f", lab, err))
    ax <- dance_smoothing_axis(input, values)
    tp01 <- (ax$t_full - ax$t_rng[1]) / diff(ax$t_rng)
    gap <- max(abs(t(fda::eval.fd(tp01, values$fd_obj)) - values$smooth_data))
    chk(is.finite(gap) && gap < 1e-6,
        sprintf("%s: the fd object is the smoothed curve (max gap %.2g)", lab, gap),
        sprintf("%s: the fd object differs from the smoothed data by %.3f", lab, gap))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
