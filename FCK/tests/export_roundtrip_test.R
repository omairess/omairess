# ==============================================================================
# tests/export_roundtrip_test.R — the exported script must RUN, and must produce
# the app's numbers
#
# tests/codegen_test.R checks that the generated script PARSES. Parsing is a low
# bar: a script full of undefined symbols parses perfectly. Before P3.4/P3.5 the
# generated script parsed and could not run -- the cosinor section called eleven
# fck_* helpers it never defined, the repeated-measures kernel carried
# input$rm_global_test and showNotification() in its body, and `covariates`,
# `subject_id` and `rm_factor` were referenced but never assigned. It also
# contained a hand-written FoSR reconstruction that had drifted from the app,
# and an "automatic (REML)" smoothing block that fitted lambda = 0.
#
# This test closes that gap in two stages:
#
#   1. RUNS the exported script end to end in a CLEAN R session (a separate
#      Rscript process, empty environment, nothing but the script and the
#      packages it loads), on the same data the app used.
#   2. Re-runs each estimator inside that clean session under a fixed seed and
#      compares the numbers against the app's, element by element.
#
# Stage 2 is the one that matters. Stage 1 only proves the script does not die.
#
# Run with:   Rscript tests/export_roundtrip_test.R      (from the FCK directory)
# ==============================================================================

.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly); library(fda)
})

app_dir <- if (dir.exists("ui")) "." else "FCK"
fck_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL)
}
ui_files     <- sort(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"))
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
for (f in ui_files) fck_source(f, envir = globalenv())
app_src  <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
ui_names <- unique(regmatches(app_src, gregexpr("ui_tab_[A-Za-z0-9_]+", app_src))[[1]])
ui_obj   <- do.call(tabItems, lapply(ui_names, get))

failures <- 0L
fail <- function(...) { cat("FAIL:", ..., "\n"); failures <<- failures + 1L }
ok   <- function(...) cat("ok  :", ..., "\n")

# Compare two numeric objects elementwise, ignoring names/dim, NA-matched.
same <- function(label, a, b, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) {
    fail(label, ": length", length(a), "vs", length(b)); return(invisible(FALSE))
  }
  na <- is.na(a); nb <- is.na(b)
  if (!identical(na, nb)) {
    fail(label, ": NA pattern differs (", sum(na), "vs", sum(nb), ")")
    return(invisible(FALSE))
  }
  d <- max(abs(a[!na] - b[!nb]))
  if (!is.finite(d) || d > tol) {
    fail(sprintf("%s: max |app - script| = %.3g (tol %.1g)", label, d, tol))
    return(invisible(FALSE))
  }
  ok(sprintf("%-42s max |app - script| = %.2g", label, d))
  invisible(TRUE)
}

# ---------------------------------------------------------------- the data ---
set.seed(20240917)
n_sub <- 12L; n_t <- 24L
hours <- 0:23
truth <- function(m, a, phi) m + a * cos(2 * pi * (hours - phi) / 24)
raw <- t(vapply(seq_len(n_sub), function(i)
  truth(50 + rnorm(1, 0, 4),
        8 + rnorm(1, 0, 1.5),
        16 + rnorm(1, 0, 1.2)) + rnorm(n_t, 0, 1.2),
  numeric(n_t)))
rownames(raw) <- paste0("S", seq_len(n_sub))
groups <- factor(rep(c("A", "B"), each = n_sub / 2))
covs   <- data.frame(Age = round(rnorm(n_sub, 40, 9), 2),
                     Sex = factor(rep(c("F", "M"), n_sub / 2)),
                     stringsAsFactors = FALSE)

# A repeated-measures design on the same generator: 10 subjects seen under two
# conditions, condition B shifted 1.5 h later with a slightly larger amplitude.
n_rm_sub <- 10L
rm_raw <- do.call(rbind, lapply(c("C1", "C2"), function(cond) {
  t(vapply(seq_len(n_rm_sub), function(i)
    truth(50 + rnorm(1, 0, 4),
          8 + rnorm(1, 0, 1.5) + (cond == "C2") * 1.5,
          16 + rnorm(1, 0, 1.2) + (cond == "C2") * 1.5) + rnorm(n_t, 0, 1.2),
    numeric(n_t)))
}))
rm_subject <- rep(paste0("P", seq_len(n_rm_sub)), times = 2)
rm_cond    <- factor(rep(c("C1", "C2"), each = n_rm_sub))
rownames(rm_raw) <- paste0(rm_subject, "_", rm_cond)

tmp <- file.path(tempdir(), "fck_roundtrip"); dir.create(tmp, showWarnings = FALSE)

# ---------------------------------------------------------------- the app ----
server <- function(input, output, session) {
  for (f in server_files) fck_source(f, envir = environment())

  # Straight in the server body, as tests/codegen_test.R does: an observe()
  # here would need a reactive flush that testServer does not reliably give.
  local({
    time_points <- seq(0, 1, length.out = n_t)
    basis  <- create.bspline.basis(rangeval = c(0, 1), nbasis = 12)
    al     <- fck_auto_lambda(raw, time_points, basis)
    lambda <- al$lambda
    fdp    <- fdPar(basis, 2, lambda)
    fd_obj <- smooth.basis(time_points, t(raw), fdp)$fd
    smooth_curves <- t(eval.fd(time_points, fd_obj))

    # ---- the app's numbers -------------------------------------------------
    set.seed(101)
    app_fanova <- perform_functional_anova(fd_obj, groups,
                                           n_permutations = 200, alpha = 0.05)
    app_cos <- lapply(seq_len(n_sub), function(i)
      fit_cosinor(hours, smooth_curves[i, ], period = 24,
                  n_harmonics = 1, trend_type = "none", use_bounds = FALSE))
    set.seed(202)
    app_fosr <- fck_fit_fosr_ols(smooth_curves, covs, c("Age", "Sex"),
                                 use_bootstrap = TRUE, n_boot = 120)

    # ---- state the generator reads ----------------------------------------
    values$data        <- raw
    values$smooth_data <- smooth_curves
    values$time_labels <- sprintf("%02d:00", hours)
    values$group_labels <- groups
    values$covariates  <- covs
    values$smooth_fit_metrics <- list(method = "manual", n_basis = 12,
                                      lambda = lambda, mean_r_squared = NA,
                                      mean_rmse = NA, mean_df = NA,
                                      time_axis = "hours")
    values$fanova_results  <- app_fanova
    values$fanova_results$design <- "between"
    values$harmonic_model <- list(
      period = 24, n_harmonics = 1, trend_type = "none",
      time_vec = hours, using_smoothed = TRUE,
      individual_fits = app_cos,
      fck_settings = list(use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                          amplitude_min = NA, amplitude_max = NA,
                          A_sat_min = NA, A_sat_max = NA,
                          tau_min = NA, tau_max = NA))
    values$reg_model <- c(app_fosr, list(fck_settings = list(
      predictors = c("Age", "Sex"), method = "OLS_nosmooth",
      use_bootstrap = TRUE, n_boot = 120, using_smoothed = TRUE,
      n_subjects = n_sub, n_time = n_t)))

    session$setInputs(smooth_method = "manual", n_basis_manual = 12,
                      smooth_factor = -log10(lambda), n_basis = 12)

    code <- generate_analysis_code(full = TRUE)
    writeLines(code, file.path(tmp, "exported.R"))
    saveRDS(list(fanova = app_fanova, cos = app_cos, fosr = app_fosr,
                 lambda = lambda, smooth_curves = smooth_curves,
                 raw = raw, covs = covs, groups = groups),
            file.path(tmp, "app_results.rds"))
  })

  # ---- and again, for the repeated-measures design ------------------------
  # This is the path P3.5 was about: the deparsed perform_rm_fanova used to
  # carry input$rm_global_test, showNotification() and withProgress() in its
  # body, and the script referenced subject_id and rm_factor without defining
  # them, so it could not run at all.
  local({
    time_points <- seq(0, 1, length.out = n_t)
    basis  <- create.bspline.basis(rangeval = c(0, 1), nbasis = 12)
    lambda <- fck_auto_lambda(rm_raw, time_points, basis)$lambda
    fd_obj <- smooth.basis(time_points, t(rm_raw), fdPar(basis, 2, lambda))$fd

    set.seed(303)
    app_rm <- perform_rm_fanova(fd_obj, rm_subject, rm_cond,
                                n_permutations = 100, alpha = 0.05)

    values$data        <- rm_raw
    values$smooth_data <- t(eval.fd(time_points, fd_obj))
    values$time_labels <- sprintf("%02d:00", hours)
    values$group_labels <- rm_cond
    values$covariates  <- NULL
    values$reg_model   <- NULL
    values$harmonic_model <- NULL
    values$smooth_fit_metrics <- list(method = "manual", n_basis = 12,
                                      lambda = lambda, mean_r_squared = NA,
                                      mean_rmse = NA, mean_df = NA,
                                      time_axis = "hours")
    values$fanova_results <- app_rm
    values$fanova_results$design     <- "within"
    values$fanova_results$subject_id <- as.character(rm_subject)
    values$fanova_results$rm_factor  <- as.character(rm_cond)

    session$setInputs(smooth_method = "manual", n_basis_manual = 12,
                      smooth_factor = -log10(lambda), n_basis = 12)

    writeLines(generate_analysis_code(full = TRUE), file.path(tmp, "exported_rm.R"))
    saveRDS(list(rm = app_rm, lambda = lambda, raw = rm_raw,
                 subject_id = as.character(rm_subject),
                 rm_factor = as.character(rm_cond)),
            file.path(tmp, "app_results_rm.rds"))
  })
}
testServer(shinyApp(ui_obj, server), { NULL })

stopifnot(file.exists(file.path(tmp, "exported.R")))
n_lines <- length(readLines(file.path(tmp, "exported.R")))
cat(sprintf("\nGenerated %d lines of R into %s\n\n", n_lines, tmp))

# --------------------------------------------------- run it, cleanly ---------
driver <- c(
  '.libPaths(c("~/Rlib", .libPaths()))',
  'args <- commandArgs(trailingOnly = TRUE)',
  'tmp  <- args[1]',
  'app  <- readRDS(file.path(tmp, "app_results.rds"))',
  '# Bind the same data the app analysed, then run the script as written.',
  'fck_input_data       <- app$raw',
  'fck_input_covariates <- app$covs',
  'pdf(file.path(tmp, "Rplots.pdf"))',
  'source(file.path(tmp, "exported.R"), echo = FALSE)',
  'dev.off()',
  '# Stage 2: re-run each estimator from the DEFINITIONS THE SCRIPT SUPPLIED,',
  '# under the same seeds the app used, and hand the numbers back.',
  'set.seed(101)',
  'sc_fanova <- perform_functional_anova(fd_obj, group_labels,',
  '                                      n_permutations = 200, alpha = 0.05)',
  'sc_cos <- lapply(seq_len(nrow(smooth_curves)), function(i)',
  '  fit_cosinor(time_vec, smooth_curves[i, ], period = 24,',
  '              n_harmonics = 1, trend_type = "none", use_bounds = FALSE))',
  'set.seed(202)',
  'sc_fosr <- fck_fit_fosr_ols(smooth_curves, covariates, c("Age", "Sex"),',
  '                            use_bootstrap = TRUE, n_boot = 120)',
  'saveRDS(list(fanova = sc_fanova, cos = sc_cos, fosr = sc_fosr,',
  '             lambda = lambda, smooth_curves = smooth_curves),',
  '        file.path(tmp, "script_results.rds"))',
  'cat("SCRIPT-RAN-CLEAN\\n")')
writeLines(driver, file.path(tmp, "driver.R"))

res <- system2("Rscript", c("--vanilla", shQuote(file.path(tmp, "driver.R")), shQuote(tmp)),
               stdout = TRUE, stderr = TRUE)
clean_ok <- any(grepl("SCRIPT-RAN-CLEAN", res, fixed = TRUE))
if (!clean_ok) {
  fail("the exported script did NOT run in a clean R session. Output:")
  cat(paste0("      | ", tail(res, 40), collapse = "\n"), "\n")
  quit(status = 1)
}
ok("the exported script ran end to end in a clean R session")

app <- readRDS(file.path(tmp, "app_results.rds"))
scr <- readRDS(file.path(tmp, "script_results.rds"))

cat("\n-- numeric agreement, app vs exported script --------------------------\n")
same("smoothing lambda",                 app$lambda, scr$lambda)
same("smoothed curves",                  app$smooth_curves, scr$smooth_curves)

same("fANOVA F(t)",                      app$fanova$F_stat,            scr$fanova$F_stat)
same("fANOVA eta-squared(t)",            app$fanova$eta_squared,       scr$fanova$eta_squared)
same("fANOVA L2 statistic",              app$fanova$L2_stat,           scr$fanova$L2_stat)
same("fANOVA pointwise p (permutation)", app$fanova$p_values_pointwise, scr$fanova$p_values_pointwise)
same("fANOVA FDR-adjusted p",            app$fanova$p_values_adjusted, scr$fanova$p_values_adjusted)
same("fANOVA global L2 p",               app$fanova$p_value_L2,        scr$fanova$p_value_L2)

# fit_cosinor returns `amplitudes` / `acrophases` (vectors, one per harmonic).
# The export's reporting block read them as `amplitude` / `acrophase` and got
# NULL -- which is how P3.8 was found. Use the real names here.
getf <- function(fits, nm) vapply(fits, function(f) {
  v <- f[[nm]]
  if (isTRUE(f$success) && length(v)) as.numeric(v)[1] else NA_real_
}, numeric(1))
for (nm in c("mesor", "amplitudes", "acrophases", "acrophases_time",
             "r_squared", "adj_r_squared", "p_value", "f_stat", "aic", "bic",
             "percent_rhythm", "amp_se", "acro_se")) {
  same(paste0("cosinor ", nm, " (", length(app$cos), " subjects)"),
       getf(app$cos, nm), getf(scr$cos, nm))
}
same("cosinor fitted values (all subjects)",
     unlist(lapply(app$cos, `[[`, "fitted")),
     unlist(lapply(scr$cos, `[[`, "fitted")))

same("FoSR beta(t)",              app$fosr$beta.hat,       scr$fosr$beta.hat)
same("FoSR analytical SE(t)",     app$fosr$beta.se,        scr$fosr$beta.se)
same("FoSR raw p(t)",             app$fosr$beta.p.raw,     scr$fosr$beta.p.raw)
same("FoSR FDR-adjusted p(t)",    app$fosr$beta.p,         scr$fosr$beta.p)
same("FoSR R2(t)",                app$fosr$r2_t,           scr$fosr$r2_t)
same("FoSR bootstrap CI lower",   app$fosr$boot_ci_lower,  scr$fosr$boot_ci_lower)
same("FoSR bootstrap CI upper",   app$fosr$boot_ci_upper,  scr$fosr$boot_ci_upper)

# ---------------------------------- the repeated-measures script -------------
cat("\n-- repeated-measures design ------------------------------------------\n")
rm_driver <- c(
  '.libPaths(c("~/Rlib", .libPaths()))',
  'tmp <- commandArgs(trailingOnly = TRUE)[1]',
  'app <- readRDS(file.path(tmp, "app_results_rm.rds"))',
  'fck_input_data <- app$raw',
  'pdf(file.path(tmp, "Rplots_rm.pdf"))',
  'source(file.path(tmp, "exported_rm.R"), echo = FALSE)',
  'dev.off()',
  '# The design vectors must have come from the SCRIPT, not from the app.',
  'stopifnot(identical(as.character(subject_id), app$subject_id))',
  'stopifnot(identical(as.character(rm_factor),  app$rm_factor))',
  'set.seed(303)',
  'sc_rm <- perform_rm_fanova(fd_obj, subject_id, rm_factor,',
  '                           n_permutations = 100, alpha = 0.05)',
  'saveRDS(sc_rm, file.path(tmp, "script_results_rm.rds"))',
  'cat("RM-SCRIPT-RAN-CLEAN\\n")')
writeLines(rm_driver, file.path(tmp, "driver_rm.R"))
res_rm <- system2("Rscript", c("--vanilla", shQuote(file.path(tmp, "driver_rm.R")), shQuote(tmp)),
                  stdout = TRUE, stderr = TRUE)
if (!any(grepl("RM-SCRIPT-RAN-CLEAN", res_rm, fixed = TRUE))) {
  fail("the exported REPEATED-MEASURES script did NOT run in a clean R session. Output:")
  cat(paste0("      | ", tail(res_rm, 40), collapse = "\n"), "\n")
} else {
  ok("the exported repeated-measures script ran end to end in a clean session")
  app_rm <- readRDS(file.path(tmp, "app_results_rm.rds"))$rm
  scr_rm <- readRDS(file.path(tmp, "script_results_rm.rds"))
  same("RM F(t)",                        app_rm$F_stat,             scr_rm$F_stat)
  same("RM partial eta-squared(t)",      app_rm$eta_squared,        scr_rm$eta_squared)
  same("RM SS_condition(t)",             app_rm$SS_condition,       scr_rm$SS_condition)
  same("RM SS_error(t)",                 app_rm$SS_error,           scr_rm$SS_error)
  same("RM pointwise p (permutation)",   app_rm$p_values_pointwise, scr_rm$p_values_pointwise)
  same("RM FDR-adjusted p",              app_rm$p_values_adjusted,  scr_rm$p_values_adjusted)
  same("RM L2 statistic",                app_rm$L2_stat,            scr_rm$L2_stat)
  if (!identical(app_rm$eta_squared_type, "partial"))
    fail("the repeated-measures result is not labelled partial eta-squared")
  else ok("repeated-measures effect size is labelled partial eta-squared (P3.3)")
}

cat("\n")
if (failures > 0L) {
  cat(sprintf("Export round-trip FAILED: %d check(s) did not agree.\n", failures))
  quit(status = 1)
}
cat("Export round-trip passed: the exported script runs clean and reproduces\n")
cat("every checked quantity to 1e-8.\n")
