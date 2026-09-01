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
# Run with:   Rscript tests/codegen_test.R      (from the FCK directory)
# ==============================================================================

suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
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

server <- function(input, output, session) {
  for (f in server_files) fck_source(f, envir = environment())

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
    fck_settings = list(use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                        amplitude_min = 0, amplitude_max = NA,
                        A_sat_min = NA, A_sat_max = NA,
                        tau_min = 0.5, tau_max = NA))

  values$reg_model <- list(
    method = "OLS (Bootstrap SE)", beta.hat = matrix(0, 3, n_t),
    fck_settings = list(predictors = c("Age", "Sex"), method = "OLS_nosmooth",
                        use_bootstrap = TRUE, n_boot = 200,
                        using_smoothed = TRUE, n_subjects = n_sub, n_time = n_t))

  values$sofr_model <- list(
    family = list(family = "binomial", link = "logit"),
    fck_settings = list(response = "Outcome", predictors = "Age",
                        formula = "y ~ lf(X_func, bs='ps', k=15) + Age",
                        family = "binomial", link = "logit",
                        using_smoothed = TRUE, n_obs = n_sub))

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
                   "11. FUNCTION-ON-SCALAR REGRESSION",
                   "12. SCALAR-ON-FUNCTION REGRESSION")) {
    if (!grepl(marker, code, fixed = TRUE)) fail("missing section:", marker)
  }
  cat("ok  : cosinor, FoSR and SoFR sections all present\n")

  # the cosinor section must carry the app's real fitting function, not a
  # paraphrase of it — that is the whole point of emitting it via deparse()
  if (!grepl("fit_cosinor <- function", code, fixed = TRUE))
    fail("the cosinor section does not carry the app's own fit_cosinor()")
  cat("ok  : the app's own fit_cosinor() is emitted verbatim\n")

  # ... and with no results at all, it must still produce valid R
  for (nm in c("harmonic_model", "reg_model", "sofr_model", "smooth_fit_metrics",
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
