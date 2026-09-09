# ==============================================================================
# tests/codegen_test.R — the exported R script must be valid R
#
# The Export tab writes a script meant to reproduce the analysis in plain R. A
# script that does not even parse is worse than no script at all: it looks like
# a reproducibility guarantee and is not one. This drives the generator with a
# stub result for each analysis family and checks that what comes out parses,
# and that each family actually contributed a section.
#
# It does not run the emitted script (that needs fda, refund and real data) —
# it checks that the generator produces syntactically valid R for every
# combination of results it might be asked about.
#
# Run with:   Rscript tests/codegen_test.R      (from the DANCE directory)
# ==============================================================================

suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
})

app_dir <- if (dir.exists("ui")) "." else "DANCE"
dance_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL)
}
ui_files     <- sort(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"))
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
for (f in ui_files) dance_source(f, envir = globalenv())

app_src  <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
ui_names <- unique(regmatches(app_src, gregexpr("ui_tab_[A-Za-z0-9_]+", app_src))[[1]])
ui_obj   <- do.call(tabItems, lapply(ui_names, get))

server <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())

  # --- a stub result for every family the generator writes a section for -----
  n_sub <- 6; n_t <- 24
  values$data        <- matrix(rnorm(n_sub * n_t), n_sub, n_t)
  values$smooth_data <- values$data
  values$time_labels <- sprintf("%02d:00", 0:23)
  values$group_labels <- factor(rep(c("A", "B"), each = n_sub / 2))
  values$covariates  <- data.frame(Age = rnorm(n_sub), Sex = factor(rep(c("F", "M"), 3)),
                                   Outcome = rbinom(n_sub, 1, 0.5))
  values$smooth_fit_metrics <- list(method = "auto", n_basis = 12, lambda = 0,
                                    mean_r_squared = 0.9, mean_rmse = 0.1,
                                    mean_df = 8, time_axis = "column index")

  values$harmonic_model <- list(
    period = 24, n_harmonics = 2, trend_type = "none",
    time_vec = 0:23, using_smoothed = TRUE,
    individual_fits = replicate(n_sub, list(success = TRUE), simplify = FALSE),
    dance_settings = list(use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                        amplitude_min = 0, amplitude_max = NA,
                        A_sat_min = NA, A_sat_max = NA,
                        tau_min = 0.5, tau_max = NA))

  values$reg_model <- list(
    method = "OLS (Bootstrap SE)", beta.hat = matrix(0, 3, n_t),
    dance_settings = list(predictors = c("Age", "Sex"), method = "OLS_nosmooth",
                        use_bootstrap = TRUE, n_boot = 200,
                        using_smoothed = TRUE, n_subjects = n_sub, n_time = n_t))

  # P20/R6: Section 8 was never exercised by this test, and it was writing its
  # own approximation of the post-hoc tests -- unpaired Welch t.test() calls
  # over an undefined `n_time_eval`, regardless of the design. A section no test
  # drives is a section that can be wrong for as long as nobody reads it.
  stub_pairwise <- function(design) {
    one <- list(group1 = "A", group2 = "B", n1 = 3, n2 = 3, design = design,
                mean_diff = rep(0, 10), t_stat = rep(0, 10),
                p_values_pointwise = rep(0.5, 10), L2_stat = 1.23,
                p_value_L2 = 0.4, p_value_L2_adjusted = 0.4,
                ci_lower = rep(-1, 10), ci_upper = rep(1, 10),
                cohens_d = rep(0, 10), se_diff = rep(1, 10),
                p_values_adjusted = rep(0.5, 10), sig_regions = rep(FALSE, 10),
                sig_global = FALSE)
    list(results = list(`A vs B` = one), time_points = seq(0, 1, length.out = 10),
         correction_method = "bonferroni", alpha = 0.05, n_permutations = 200,
         groups = c("A", "B"), n_groups = 2, pair_names = "A vs B",
         p_floor = 1 / 201, design = design,
         spec = list(
           design = design, source = "fanova", description = "stub",
           group_labels = if (design == "between") rep(c("A", "B"), each = 3) else NULL,
           subject_id = if (design == "within") sprintf("s%d", rep(1:3, 2)) else NULL,
           rm_factor = if (design == "within") rep(c("A", "B"), each = 3) else NULL,
           used_warped_curves = FALSE, n_permutations = 200,
           correction = "bonferroni", alpha = 0.05, p_floor = 1 / 201,
           family = "1 pairwise comparison", run_at = Sys.time()))
  }

fail <- function(...) { cat("FAIL:", ..., "\n"); quit(status = 1) }

  code <- tryCatch(generate_analysis_code(full = TRUE),
                   error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(code, "err")) fail("generator errored:", code)

  parsed <- tryCatch({ parse(text = code); TRUE },
                     error = function(e) conditionMessage(e))
  if (!isTRUE(parsed)) fail("the exported script is not valid R:", parsed)
  cat(sprintf("ok  : full export parses (%d lines)\n",
              length(strsplit(code, "\n")[[1]])))

  for (marker in c("10. HARMONIC (COSINOR) REGRESSION",
                   "11. FUNCTION-ON-SCALAR REGRESSION")) {
    if (!grepl(marker, code, fixed = TRUE)) fail("missing section:", marker)
  }
  cat("ok  : cosinor and FoSR sections all present\n")

  # the cosinor section must carry the app's real fitting function, not a
  # paraphrase of it — that is the whole point of emitting it via deparse()
  if (!grepl("fit_cosinor <- function", code, fixed = TRUE))
    fail("the cosinor section does not carry the app's own fit_cosinor()")
  cat("ok  : the app's own fit_cosinor() is emitted verbatim\n")

  # --- P20/R6: the post-hoc section emits the app's own kernels -------------
  for (dsn in c("between", "within")) {
    values$pairwise_results <- stub_pairwise(dsn)
    ph <- tryCatch(generate_analysis_code(full = TRUE),
                   error = function(e) structure(conditionMessage(e), class = "err"))
    if (inherits(ph, "err")) fail("generator errored with", dsn, "post-hoc:", ph)
    okp <- tryCatch({ parse(text = ph); TRUE }, error = function(e) conditionMessage(e))
    if (!isTRUE(okp)) fail("the", dsn, "post-hoc export is not valid R:", okp)
    want <- if (dsn == "within") "perform_pairwise_comparisons_rm <- function"
            else "perform_pairwise_comparisons <- function"
    if (!grepl(want, ph, fixed = TRUE))
      fail("the", dsn, "post-hoc section does not carry the app's own kernel")
    # the approximation that used to stand in for it must be gone
    if (grepl("n_time_eval", ph, fixed = TRUE))
      fail("the export still references the undefined n_time_eval")
    if (grepl("tt <- t.test(curves_eval[idx1, t], curves_eval[idx2, t])", ph, fixed = TRUE))
      fail("the export still writes its own unpaired t-test loop")
    if (!grepl("Multiplicity family:", ph, fixed = TRUE))
      fail("the", dsn, "post-hoc section does not name its multiplicity family")
    cat(sprintf("ok  : %s post-hoc section emits the real kernel, not an approximation\n", dsn))
  }
  values$pairwise_results <- NULL

  # ... and with no results at all, it must still produce valid R
  for (nm in c("harmonic_model", "reg_model", "smooth_fit_metrics",
               "group_labels", "covariates", "smooth_data"))
    values[[nm]] <- NULL
  bare <- tryCatch(generate_analysis_code(full = TRUE),
                   error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(bare, "err")) fail("generator errored with no results:", bare)
  ok2 <- tryCatch({ parse(text = bare); TRUE }, error = function(e) conditionMessage(e))
  if (!isTRUE(ok2)) fail("data-only export is not valid R:", ok2)
  cat("ok  : data-only export parses too\n")

  cat("\nCode-generator tests passed.\n")
}

testServer(shinyApp(ui_obj, server), { NULL })
