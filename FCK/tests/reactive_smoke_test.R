# ==============================================================================
# tests/reactive_smoke_test.R — drive the app's reactives on real-shaped data
#
# WHY THIS EXISTS. Two bugs shipped in the last two rounds and neither was
# visible to any test in this repo:
#
#   * the FoSR GAM branch died with "object 'j' not found" -- a loop whose
#     header was renamed but whose body was not;
#   * the fPCA-ANOVA report died with "$ operator is invalid for atomic
#     vectors" -- R's `$` partial-matching on a list, on the unwarped path,
#     which is the common one.
#
# Both were caught by a person opening the app. Everything the suite had was
# either static (does it parse, does the source contain X) or numerical (does
# the kernel return the right number for the right input). Nothing pressed the
# buttons. This does: it builds a dataset shaped like the real one -- sorted by
# group, unbalanced, four levels -- runs the analyses through the server, and
# FORCES EVERY RELEVANT OUTPUT TO RENDER. An output that errors fails the test.
#
# Run with:   Rscript tests/reactive_smoke_test.R      (from the FCK directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
  library(fda); library(mgcv); library(ggplot2)
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

# Force an output to render and report the error instead of hiding it.
render <- function(label, expr) {
  v <- tryCatch(force(expr), error = function(e) structure(conditionMessage(e), class = "fckerr"))
  if (inherits(v, "fckerr")) { fail(label, "->", as.character(v)); return(invisible(NULL)) }
  ok(label)
  invisible(v)
}

# ---- a dataset shaped like the one that broke the app ----------------------
# Sorted by group, badly unbalanced, four levels: exactly the layout that made
# "the first 200 rows" a single category (P6.5).
set.seed(606)
n_by  <- c(YOUTH = 60, ADULT = 40, MIDDLE_AGE = 18, ELDERLY = 8)
n_sub <- sum(n_by); n_t <- 24; hrs <- 0:23
grp   <- factor(rep(names(n_by), n_by), levels = names(n_by))
shift_by <- c(YOUTH = 0, ADULT = 1.2, MIDDLE_AGE = 2.0, ELDERLY = 3.1)
raw <- t(vapply(seq_len(n_sub), function(i) {
  50 + rnorm(1, 0, 4) +
    (8 + rnorm(1, 0, 1.5)) * cos(2 * pi * (hrs - 16 - shift_by[[as.character(grp[i])]]) / 24) +
    rnorm(n_t, 0, 1.2)
}, numeric(n_t)))
rownames(raw) <- sprintf("S%03d", seq_len(n_sub))
covs <- data.frame(Age = round(rnorm(n_sub, 40, 9), 2),
                   Sex = factor(rep(c("F", "M"), length.out = n_sub)),
                   AGEcategory = grp,
                   stringsAsFactors = FALSE)

server <- function(input, output, session) {
  for (f in server_files) fck_source(f, envir = environment())

  values$data        <- raw
  values$time_labels <- sprintf("%02d:00", hrs)
  values$covariates  <- covs
  values$group_labels <- grp
  values$group_variables <- list(AGEcategory = grp, Sex = covs$Sex)
  values$subject_ids <- rownames(raw)

  tp <- seq(0, 1, length.out = n_t)
  basis  <- create.bspline.basis(rangeval = c(0, 1), nbasis = 12)
  lambda <- fck_auto_lambda(raw, tp, basis)$lambda
  values$fd_obj      <- smooth.basis(tp, t(raw), fdPar(basis, 2, lambda))$fd
  values$smooth_data <- t(eval.fd(tp, values$fd_obj))
  values$smooth_fit_metrics <- list(method = "manual", n_basis = 12, lambda = lambda,
                                    mean_r_squared = NA, mean_rmse = NA,
                                    mean_df = NA, time_axis = "hours")

  session$setInputs(
    n_components = 5, pca_type = "fpca", effect_n_comp = 3, effect_size = 1,
    smooth_method = "manual", n_basis_manual = 12, n_basis = 12,
    smooth_factor = -log10(lambda), is_cyclic = FALSE,
    pca_anova_group_var = "AGEcategory", pca_anova_ncomp = 5,
    pca_anova_posthoc = "holm", pca_anova_across = "holm",
    pca_anova_gate = 0.05, pca_anova_conf = 0.95, pca_anova_omnibus = "auto",
    reg_predictors = c("Age", "AGEcategory"), reg_color_var = "AGEcategory",
    reg_method = "OLS_nosmooth", use_bootstrap = FALSE, n_boot = 0,
    tick_freq_fanova = 4, tick_freq_results = 4, tick_freq_settings = 4,
    tick_freq_preprocess = 4, tick_freq_pairwise = 4, tick_freq_kmeans = 4)

  # ---------------------------------------------------------------- fPCA ----
  cat("\n-- fPCA, 5 components requested ---------------------------------------\n")
  values$pca_results <- pca.fd(values$fd_obj, nharm = 5)
  session$flushReact()   # make the assignment visible to the outputs below

  n_disp <- fck_n_harmonics(values$pca_results)
  if (n_disp != 5L) fail(sprintf("fck_n_harmonics reports %d, expected 5", n_disp))
  else ok("the PCA result reports all 5 components")

  s <- render("pca_summary renders", output$pca_summary)
  if (!is.null(s)) {
    txt <- paste(as.character(s), collapse = "\n")
    for (pc in paste0("PC", 1:5))
      if (!grepl(pc, txt, fixed = TRUE)) fail(sprintf("the summary does not mention %s", pc))
    if (grepl("PC1", txt) && grepl("PC5", txt)) ok("the summary lists PC1 through PC5")
  }
  render("loadings_plot renders", output$loadings_plot)
  render("scores_plot renders", output$scores_plot)
  render("variance_plot renders", output$variance_plot)
  render("scores_table renders", output$scores_table)

  # ------------------------------------------------------- fPCA group ANOVA --
  cat("\n-- fPCA-ANOVA on an UNWARPED run (the path that errored) ---------------\n")
  if (!is.null(values$warping_results)) fail("this run should not be warped")
  session$setInputs(run_pca_anova = 1)
  session$flushReact()
  render("pca_anova_results renders", output$pca_anova_results)

  # ---------------------------------------------------------------- FoSR ----
  cat("\n-- FoSR, pointwise OLS -------------------------------------------------\n")
  session$setInputs(run_fosr = 1)
  session$flushReact()
  if (is.null(values$reg_model)) fail("the OLS fit produced no model")
  else ok(sprintf("OLS fitted: %s", values$reg_model$method))
  render("reg_observed_plot renders", output$reg_observed_plot)
  render("fosr_model_summary renders", output$fosr_model_summary)

  if (!is.null(values$reg_model)) {
    session$setInputs(reg_coeff_select = rownames(values$reg_model$beta.hat)[2])
    render("reg_coeff_plot renders", output$reg_coeff_plot)
    render("reg_pvalue_plot renders", output$reg_pvalue_plot)
    render("reg_r2_plot renders", output$reg_r2_plot)

    # the prediction path that P5.1 broke -- exercise the shared builder
    nd <- as.data.frame(lapply(covs[, c("Age", "AGEcategory")], function(v)
      if (is.numeric(v)) mean(v) else factor(levels(v)[1], levels = levels(v))),
      check.names = FALSE)
    X <- tryCatch(fck_fosr_design(values$reg_model, nd),
                  error = function(e) structure(conditionMessage(e), class = "fckerr"))
    if (inherits(X, "fckerr")) fail("fck_fosr_design ->", as.character(X))
    else if (ncol(X) != nrow(values$reg_model$beta.hat))
      fail("the prediction design has the wrong number of columns")
    else {
      yh <- as.vector(X %*% values$reg_model$beta.hat)
      if (!all(is.finite(yh))) fail("the predicted curve is not finite")
      else ok("a prediction design and curve can be built from the fitted model")
    }
  }

  cat("\n-- FoSR, smoothed OLS (GAM) -------------------------------------------\n")
  session$setInputs(reg_method = "OLS_smooth")
  session$setInputs(run_fosr = 2)
  session$flushReact()
  m <- values$reg_model
  if (is.null(m) || is.null(m$gam_obj)) fail("the GAM fit produced no model")
  else {
    ok(sprintf("GAM fitted: %d coefficient curves", nrow(m$beta.hat)))
    if (!all(is.na(m$beta.se))) fail("the GAM reports non-NA standard errors it did not compute")
    else ok("GAM SE/p are NA, not zero")
    # one contrast curve per non-reference level of the factor
    if (!all(c("AGEcategoryADULT", "AGEcategoryMIDDLE_AGE", "AGEcategoryELDERLY")
             %in% rownames(m$beta.hat)))
      fail("the GAM did not emit a curve per factor level: ",
           paste(rownames(m$beta.hat), collapse = ", "))
    else ok("the GAM emits one coefficient curve per non-reference level")
  }
  render("fosr_model_summary renders (GAM)", output$fosr_model_summary)
  render("reg_pvalue_plot renders (GAM, all-NA p)", output$reg_pvalue_plot)
  render("reg_coeff_plot renders (GAM, all-NA se)", output$reg_coeff_plot)

  # P8.1: the GAM's formula refers to x1..xp, so a prediction frame keyed by the
  # user's column names cannot be found by predict.gam(). The previous round
  # tested that the GAM FITS; it did not test that it PREDICTS, and that is
  # exactly where the bug was. Drive the real helper.
  if (!is.null(m) && !is.null(m$gam_obj)) {
    iv <- list(Age = mean(covs$Age), AGEcategory = "ADULT")
    pd <- tryCatch(build_gam_pred_df(seq(0, 1, length.out = n_t),
                                     m$gam_predictors, m$gam_model_names,
                                     m$gam_long_data, m$gam_factor_levels, iv),
                   error = function(e) structure(conditionMessage(e), class = "fckerr"))
    if (inherits(pd, "fckerr")) {
      fail("build_gam_pred_df ->", as.character(pd))
    } else if (!all(m$gam_model_names %in% names(pd))) {
      fail("the GAM prediction frame is not keyed by the model's own names: ",
           paste(names(pd), collapse = ", "))
    } else {
      ok(sprintf("GAM prediction frame is keyed by the model names (%s)",
                 paste(m$gam_model_names, collapse = ", ")))
      pv <- tryCatch(predict(m$gam_obj, newdata = pd, se.fit = TRUE),
                     error = function(e) structure(conditionMessage(e), class = "fckerr"))
      if (inherits(pv, "fckerr")) fail("predict.gam ->", as.character(pv))
      else if (!all(is.finite(pv$fit)) || !all(is.finite(pv$se.fit)))
        fail("the GAM prediction is not finite")
      else ok("predict.gam returns a finite curve with standard errors")
    }
    # and the mapping is not optional any more
    bad <- tryCatch({
      build_gam_pred_df(seq(0, 1, length.out = n_t), m$gam_predictors, NULL,
                        m$gam_long_data, m$gam_factor_levels, iv); "no error" },
      error = function(e) "refused")
    if (!identical(bad, "refused"))
      fail("build_gam_pred_df accepted a call with no model-name mapping")
    else ok("a call without the name mapping is refused, not silently wrong")
  }

  # ----------------------------------------------- fANOVA and its post-hocs --
  # P7.1: the pairwise permutation box was overridden by the omnibus count, so
  # a value typed there never reached the test. The check that catches that
  # class of bug is: set a control to a DISTINCTIVE value and assert the result
  # carries it. A static "is input$x referenced anywhere" sweep would not have
  # caught it -- the reference existed, in dead code.
  cat("\n-- fANOVA, then pairwise with a distinct permutation count -------------\n")
  session$setInputs(fanova_design = "between", fanova_group_var = "AGEcategory",
                    n_permutations = 200L, alpha_level = 0.05,
                    fanova_data_source = "smoothed", rm_global_test = FALSE)
  session$setInputs(run_fanova = 1); session$flushReact()
  if (is.null(values$fanova_results)) {
    fail("the omnibus fANOVA produced no result")
  } else {
    ok(sprintf("omnibus fANOVA ran with B = %d", values$fanova_results$n_permutations))
    if (!identical(as.integer(values$fanova_results$n_permutations), 200L))
      fail("the omnibus did not use the permutation count it was given")

    session$setInputs(pairwise_permutations = 777L, pairwise_correction = "bonferroni",
                      pairwise_alpha = 0.05, pairwise_confidence_bands = TRUE,
                      posthoc_source = "fanova")
    session$setInputs(run_pairwise = 1); session$flushReact()
    pr <- values$pairwise_results
    if (is.null(pr)) {
      fail("the pairwise comparisons produced no result")
    } else if (!identical(as.integer(pr$n_permutations), 777L)) {
      fail(sprintf("the pairwise tests used B = %s, not the 777 that was asked for",
                   as.character(pr$n_permutations)))
    } else {
      ok("the pairwise tests used the permutation count the control was set to")
      if (!identical(as.integer(pr$omnibus_permutations), 200L))
        fail("the omnibus count was not recorded alongside it")
      else ok("both counts are recorded, so a difference is visible")
    }
    render("pairwise summary renders", output$pairwise_summary)
    render("pairwise global table renders", output$pairwise_global_table)
    render("fanova effect size plot renders", output$fanova_effect_size_plot)
    render("fanova effect summary renders", output$fanova_effect_summary)
  }

  # ------------------------------------------------- the APA report --------
  cat("\n-- APA report ----------------------------------------------------------\n")
  session$setInputs(apa_report_title = "Diurnal profiles by age group")
  md <- tryCatch(fck_apa_report(values, input, "Diurnal profiles by age group"),
                 error = function(e) structure(conditionMessage(e), class = "fckerr"))
  if (inherits(md, "fckerr")) {
    fail("fck_apa_report ->", as.character(md))
  } else {
    ok(sprintf("report generated (%d lines)", length(md)))
    txt <- paste(md, collapse = "\n")
    # it must describe the analyses that ran, and only those
    must <- c("Statistical analysis", "Results", "Reproducibility",
              "Smoothing", "Functional ANOVA", "Post-hoc",
              "Function-on-scalar regression", "principal component",
              "What these numbers do not establish")
    for (k in must)
      if (!grepl(k, txt, fixed = TRUE)) fail("the report does not mention:", k)
    # and it must NOT describe analyses that did not run
    for (k in c("Functional clustering", "Registration diagnostics"))
      if (grepl(k, txt, fixed = TRUE))
        fail("the report describes an analysis that was not run:", k)
    ok("the report covers every analysis that ran, and none that did not")

    # APA formatting: no leading zero on a bounded quantity, p to 3 dp
    if (grepl("\\*p\\* = 0\\.", txt)) fail("a p-value carries a leading zero")
    else ok("p-values follow APA 7 (no leading zero)")
    # APA 7 sec. 6.36 applies to RANGES of a bounded quantity too
    if (grepl("range 0\\.", txt)) fail("a bounded quantity's range carries a leading zero")
    else ok("bounded quantities drop the leading zero in ranges as well")
    if (!grepl("smallest attainable", txt, fixed = TRUE))
      fail("the report does not state the permutation resolution floor")
    else ok("the permutation resolution floor is stated")
    if (!grepl("*B* = 200", txt, fixed = TRUE) &&
        !grepl("*B* = 777", txt, fixed = TRUE))
      fail("the report does not state the permutation counts that were used")
    else ok("the permutation counts used are stated")

    hh <- tryCatch(fck_apa_html(md, "t"), error = function(e)
                   structure(conditionMessage(e), class = "fckerr"))
    if (inherits(hh, "fckerr")) fail("fck_apa_html ->", as.character(hh))
    else {
      h1 <- paste(hh, collapse = "\n")
      if (!grepl("<!DOCTYPE html>", h1, fixed = TRUE)) fail("the HTML has no doctype")
      else if (!grepl("<table>", h1, fixed = TRUE)) fail("no table survived the HTML rendering")
      else if (grepl("| ---:", h1, fixed = TRUE)) fail("a Markdown table separator leaked into the HTML")
      else ok(sprintf("HTML rendering produced %d lines with tables intact", length(hh)))
    }
    writeLines(md, "/tmp/claude-0/fck_apa_report.md")
    writeLines(hh, "/tmp/claude-0/fck_apa_report.html")
    cat("   (written to /tmp/claude-0/fck_apa_report.md)\n")
  }
  render("apa_report_preview renders", output$apa_report_preview)

  cat("\n")
  if (failures) { cat(sprintf("Reactive smoke test FAILED (%d).\n", failures)); quit(status = 1) }
  cat("Reactive smoke test passed.\n")
}

testServer(shinyApp(ui_obj, server), { NULL })
