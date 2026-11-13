# ==============================================================================
# tests/fpca_rotation_app_test.R -- VARIMAX-rotated functional PCA, app and script
#
# Runs the fPCA tab with and without the VARIMAX rotation, renders its
# outputs, reads the report, and runs the exported script on the same data:
# the script's components and scores must equal the app's, rotated or not.
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

set.seed(3); tt <- 0:23
Y <- t(sapply(1:30, function(i) 10 + 2 * rnorm(1) * sin(2 * pi * tt / 24) + rnorm(1) * cos(2 * pi * tt / 24) +
                1.5 * rnorm(1) * exp(-((tt - 17) / 2)^2) + rnorm(24, 0, 0.2)))
d <- data.frame(ID = sprintf("S%02d", 1:30), g = rep(c("a", "b"), 15), round(Y, 3))
names(d)[-(1:2)] <- sprintf("%02d:00", tt)
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE)

testServer(server, {
  session$setInputs(datafile = list(name = "f.csv", datapath = csv), header = TRUE, dec = "auto",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = sprintf("%02d:00", tt), sel_cov_vars = "g", subject_id_var = "ID",
                    apply_selection = 1)
  session$setInputs(smooth_method = "manual", n_basis_manual = 12, smooth_factor = 4, is_cyclic = FALSE,
                    use_real_time = FALSE, allow_extrapolation = FALSE, constrain_bounds = FALSE,
                    min_bound = 0, max_bound = 1, n_basis = 12, apply_smooth = 1)
  k <- 0L
  for (rot in c("none", "varimax")) {
    k <- k + 1L
    session$setInputs(pca_type = "fpca", n_components = 3, pca_rotation = rot, tick_freq_results = 2,
                      effect_n_comp = 3, run_analysis = k)
    pr <- values$pca_results
    want <- fda::pca.fd(values$fd_obj, nharm = 3); if (rot == "varimax") want <- fda::varmx.pca.fd(want)
    chk(identical(pr$rotation, rot) && isTRUE(all.equal(pr$scores, want$scores)),
        sprintf("%s: the app's components are fda's", rot), sprintf("%s: the app's components are not fda's", rot))
    summ <- paste(output$pca_summary, collapse = "\n")
    lp <- output$loadings_plot; vp <- output$variance_plot; sp <- output$scores_plot
    if (rot == "varimax")
      chk(grepl("Rotation: VARIMAX", summ) && grepl("RC1", summ) && grepl("[PC1:", summ, fixed = TRUE) &&
            grepl("VARIMAX-rotated", lp) && grepl("not ordered by variance", vp),
          "varimax: the summary and the plots say the components are rotated",
          "varimax: the rotation is not stated in the summary or the plots")
    else chk(!grepl("VARIMAX", summ), "none: no rotation is mentioned", "none: a rotation is mentioned")
    rep <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
    chk(identical(grepl("rotated by VARIMAX", rep), rot == "varimax") &&
          identical(grepl("Unrotated PC (%)", rep, fixed = TRUE), rot == "varimax"),
        sprintf("%s: the report states the rotation as it was", rot), sprintf("%s: the report misstates the rotation", rot))
    code <- generate_analysis_code(full = TRUE)
    env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
    ok <- tryCatch({ pdf(NULL); on.exit(dev.off(), add = TRUE); eval(parse(text = code), envir = env); TRUE },
                   error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
    sp_ <- if (ok) env$pca_results else NULL
    chk(ok && identical(sp_$rotation, rot) && isTRUE(all.equal(sp_$scores, pr$scores, tolerance = 1e-10)) &&
          isTRUE(all.equal(sp_$harmonics$coefs, pr$harmonics$coefs, tolerance = 1e-10)) &&
          isTRUE(all.equal(sp_$varprop, pr$varprop, tolerance = 1e-10)),
        sprintf("%s: the exported script reproduces the components, scores and variance", rot),
        sprintf("%s: the exported script does not reproduce the fPCA", rot))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
