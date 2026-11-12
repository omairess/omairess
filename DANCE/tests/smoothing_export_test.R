# ==============================================================================
# tests/smoothing_export_test.R -- the smoothing step, app against exported script
#
# The script used to retype a simpler smooth (every curve at once on 0-1), so
# a cyclic (Fourier) run, a real-time run or any run with missing values could
# not be reproduced from it. Both now call dance_smooth_curves()
# (server/04_helpers_fd.R). For each setting below, the REAL smoothing tab is
# run, the script is generated, run in a clean R session on the same data, and
# its smoothed curves, lambda and fd object are compared with the app's.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({ library(shiny); library(shinydashboard); library(DT); library(plotly)
  for (p in c("fda", "shinyWidgets")) if (requireNamespace(p, quietly = TRUE)) library(p, character.only = TRUE) })
app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }
tmp <- file.path(tempdir(), "dance_smooth_export"); dir.create(tmp, showWarnings = FALSE)

# 23 clock columns with one hour missing (so real time differs from the index),
# gaps at the ends and inside, one curve with a single observation
tt <- c(0:11, 13:23)
set.seed(5)
Y <- t(vapply(1:14, function(i) 10 + i / 5 + 5 * cos(2 * pi * (tt - 3 - i / 7) / 24) + rnorm(length(tt), 0, .6),
              numeric(length(tt))))
Y[2, 1:3] <- NA; Y[5, 20:23] <- NA; Y[7, c(4, 9, 15)] <- NA; Y[9, ] <- NA; Y[9, 6] <- 12
d <- data.frame(ID = sprintf("S%02d", 1:14), g = rep(c("a", "b"), 7), round(Y, 3))
names(d)[-(1:2)] <- sprintf("%02d:00", tt)
csv <- file.path(tmp, "in.csv"); write.csv(d, csv, row.names = FALSE, na = "")

base <- list(min_bound = 0, max_bound = 1, n_basis = 20, n_basis_manual = 20, smooth_factor = 1,
             allow_extrapolation = FALSE, constrain_bounds = FALSE, use_real_time = FALSE, is_cyclic = FALSE)
cases <- list(
  "B-spline, manual"                  = list(smooth_method = "manual", n_basis_manual = 9, smooth_factor = 1),
  "B-spline, real time, extrapolated" = list(smooth_method = "manual", use_real_time = TRUE, n_basis_manual = 9,
                                             smooth_factor = 2, allow_extrapolation = TRUE),
  "B-spline, GCV, bounds"             = list(smooth_method = "auto", n_basis = 12, constrain_bounds = TRUE,
                                             min_bound = 6, max_bound = 14),
  "Fourier (cyclic), column index"    = list(smooth_method = "manual", is_cyclic = TRUE, n_basis_manual = 9,
                                             smooth_factor = 3),
  "Fourier (cyclic), real time, GCV"  = list(smooth_method = "auto", is_cyclic = TRUE, use_real_time = TRUE,
                                             allow_extrapolation = TRUE),
  "none"                              = list(smooth_method = "none"))

app <- list()
testServer(server, {
  session$setInputs(datafile = list(name = "in.csv", datapath = csv), header = TRUE, dec = "auto",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = sprintf("%02d:00", tt), sel_cov_vars = "g", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  for (nm in names(cases)) {
    k <- k + 1L
    do.call(session$setInputs, c(modifyList(base, cases[[nm]]), list(apply_smooth = k)))
    code <- generate_analysis_code(full = TRUE)
    f <- file.path(tmp, sprintf("script_%d.R", k)); writeLines(code, f)
    app[[nm]] <<- list(file = f, data = values$data, smooth = values$smooth_data,
                       coefs = values$fd_obj$coefs, period01 = values$fd_obj$basis$params,
                       lambda = values$smooth_fit_metrics$lambda,
                       cyclic_period = values$smooth_fit_metrics$cyclic_period, code = code)
  }
})

drv <- file.path(tmp, "driver.R")
writeLines(c('.libPaths(c("~/Rlib", .libPaths()))',
             'a <- commandArgs(trailingOnly = TRUE)',
             'dance_input_data <- readRDS(a[2])',
             'pdf(NULL)',
             'source(a[1], echo = FALSE)',
             'saveRDS(list(smooth = smooth_curves, coefs = fd_obj$coefs, period01 = fd_obj$basis$params,',
             '             lambda = if (exists("lambda")) lambda else NA), a[3])',
             'cat("SCRIPT-RAN-CLEAN\\n")'), drv)

for (nm in names(app)) {
  a <- app[[nm]]
  din <- file.path(tmp, "data.rds"); saveRDS(a$data, din)
  out <- file.path(tmp, "out.rds"); unlink(out)
  res <- system2("Rscript", c("--vanilla", shQuote(drv), shQuote(a$file), shQuote(din), shQuote(out)),
                 stdout = TRUE, stderr = TRUE)
  if (!any(grepl("SCRIPT-RAN-CLEAN", res, fixed = TRUE))) {
    chk(FALSE, "", sprintf("%s: the script did not run:\n%s", nm, paste(tail(res, 15), collapse = "\n")))
    next
  }
  s <- readRDS(out)
  same <- function(x, y) isTRUE(all.equal(unname(as.matrix(x)), unname(as.matrix(y)), tolerance = 1e-10))
  chk(same(a$smooth, s$smooth), sprintf("%-36s smoothed curves equal", nm),
      sprintf("%s: smoothed curves differ (max %.3g)", nm, max(abs(a$smooth - s$smooth), na.rm = TRUE)))
  chk(same(a$coefs, s$coefs) && isTRUE(all.equal(a$period01, s$period01)),
      sprintf("%-36s fd object equal (basis and coefficients)", nm), sprintf("%s: the fd object differs", nm))
  if (!identical(nm, "none"))
    chk(isTRUE(all.equal(a$lambda, s$lambda, tolerance = 1e-12)),
        sprintf("%-36s lambda equal (%.4g)", nm, a$lambda), sprintf("%s: lambda %s vs %s", nm, a$lambda, s$lambda))
  if (is.finite(a$cyclic_period %||% NA))
    chk(grepl(sprintf("one cycle = %s", format(signif(a$cyclic_period, 6))), a$code, fixed = TRUE),
        sprintf("%-36s the script states the cycle (%s)", nm, format(a$cyclic_period)),
        sprintf("%s: the script does not state the cycle", nm))
}
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
