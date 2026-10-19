# ==============================================================================
# tests/warp_report_export_test.R -- the shift model and the per-participant
# warps in the APA report and the exported script
#
# Drives Run (no smoothing, time-warped PCA with a linear shift, the shift
# model), then checks that the report carries the Methods and Results
# sentences, and that the exported script, run in a clean `Rscript --vanilla`
# session on the same data, reproduces the shift model and every participant's
# warp exactly.
#
# Run with:   Rscript tests/warp_report_export_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(shiny); library(shinydashboard); library(DT); library(plotly)
  library(fda) })
if (requireNamespace("shinyWidgets", quietly = TRUE)) suppressPackageStartupMessages(library(shinyWidgets))
app_dir <- if (dir.exists("server")) "." else "DANCE"
dance_source <- function(file, envir) { eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL) }
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }

make_data <- function() {
  set.seed(2); t <- 1:18; n <- 24; grp <- factor(rep(c("A", "B"), each = 12))
  f <- function(x) 10 + 50 * plogis((x - 9) / 1.5)
  s <- rnorm(n, ifelse(grp == "B", 1.5, 0), 1)
  raw <- t(sapply(1:n, function(i) rnorm(1, 0, 5) + f(t - s[i]) + rnorm(18, 0, 5)))
  rownames(raw) <- sprintf("S%02d", 1:n); colnames(raw) <- sprintf("%dh00", t)
  list(raw = raw, grp = grp)
}
dd <- make_data()
server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data <- dd$raw; values$time_labels <- colnames(dd$raw)
  values$group_labels <- dd$grp; values$subject_ids <- rownames(dd$raw)
}
tmp <- tempfile("warp_export_"); dir.create(tmp)
testServer(server, {
  session$setInputs(smooth_method = "none", is_cyclic = FALSE, use_real_time = TRUE, apply_smooth = 1)
  session$setInputs(pca_type = "twpca", n_components = 3, warping_method = "linear_shift",
                    periodic_shift = FALSE, shift_reference = "mean", run_analysis = 1)
  session$setInputs(shift_model_df = 6, run_shift_model = 1)
  sm <- values$shift_model
  chk(isTRUE(sm$ok), sprintf("the shift model ran: SD(s) = %.2f h", sm$sd_shift), "the shift model did not run")
  rp <- paste(dance_apa_report(values, input, "T"), collapse = "\n")
  chk(grepl("**Timing between participants.**", rp, fixed = TRUE) &&
        grepl("shape-invariant (self-modelling) mixed model", rp, fixed = TRUE),
      "Methods describe the shift model", "Methods omit the shift model")
  chk(grepl(sprintf("standard deviation of %s h", dance_apa_num(sm$sd_shift, 2)), rp, fixed = TRUE) &&
        grepl("against a model without shifts", rp, fixed = TRUE),
      "Results report SD(s) and the test against no shifts", "Results omit the shift model")
  chk(grepl("compared between groups with Welch's one-way ANOVA", rp, fixed = TRUE) &&
        grepl("Welch one-way ANOVA, *F*(", rp, fixed = TRUE) && grepl("Kruskal-Wallis, *H*(", rp, fixed = TRUE),
      "the report compares the per-participant shifts between groups", "the report omits the group comparison")
  chk(!grepl("Kneip and Ramsay (2008)", rp, fixed = TRUE),
      "no amplitude/phase decomposition claimed for a non-periodic shift",
      "the report claims a decomposition that is not defined here")
  writeLines(generate_analysis_code(), file.path(tmp, "script.R"))
  saveRDS(list(sd = sm$sd_shift, ci = sm$sd_shift_ci, p = sm$lrt$p, shifts = sm$shifts$shift,
               est = warp_params()$estimate), file.path(tmp, "app.rds"))
})

runner <- file.path(tmp, "run.R")
writeLines(c(
  '.libPaths(c("~/Rlib", .libPaths()))',
  sprintf('dance_input_data <- readRDS("%s")', file.path(tmp, "data.rds")),
  'pdf(NULL)',
  sprintf('source("%s")', file.path(tmp, "script.R")),
  sprintf('saveRDS(list(sd = shift_fit$sd_shift, ci = shift_fit$sd_shift_ci, p = shift_fit$lrt$p, shifts = shift_fit$shifts$shift, est = warp_est$estimate), "%s")',
          file.path(tmp, "script.rds"))), runner)
saveRDS(dd$raw, file.path(tmp, "data.rds"))
st <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", runner), stdout = FALSE, stderr = FALSE)
chk(identical(st, 0L), "the exported script runs in a clean session", "the exported script failed")
if (identical(st, 0L)) {
  a <- readRDS(file.path(tmp, "app.rds")); b <- readRDS(file.path(tmp, "script.rds"))
  chk(isTRUE(all.equal(a$sd, b$sd)) && isTRUE(all.equal(a$ci, b$ci)) && isTRUE(all.equal(a$p, b$p)) &&
        isTRUE(all.equal(a$shifts, b$shifts)),
      "the script reproduces SD(s), its interval, the test and the predicted shifts",
      "the script's shift model differs from the app's")
  chk(isTRUE(all.equal(a$est, b$est)), "the script reproduces every participant's warp",
      "the script's per-participant warps differ")
}
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
