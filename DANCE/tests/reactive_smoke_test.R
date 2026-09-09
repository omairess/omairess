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
# Run with:   Rscript tests/reactive_smoke_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
  library(fda); library(mgcv); library(ggplot2)
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
  for (f in server_files) dance_source(f, envir = environment())

  values$data        <- raw
  values$time_labels <- sprintf("%02d:00", hrs)
  values$covariates  <- covs
  values$group_labels <- grp
  values$group_variables <- list(AGEcategory = grp, Sex = covs$Sex)
  values$subject_ids <- rownames(raw)

  tp <- seq(0, 1, length.out = n_t)
  basis  <- create.bspline.basis(rangeval = c(0, 1), nbasis = 12)
  lambda <- dance_auto_lambda(raw, tp, basis)$lambda
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

  n_disp <- dance_n_harmonics(values$pca_results)
  if (n_disp != 5L) fail(sprintf("dance_n_harmonics reports %d, expected 5", n_disp))
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
    X <- tryCatch(dance_fosr_design(values$reg_model, nd),
                  error = function(e) structure(conditionMessage(e), class = "fckerr"))
    if (inherits(X, "fckerr")) fail("dance_fosr_design ->", as.character(X))
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

  # ------------------------------------------------- curve registration -----
  #
  # AUDIT (P11.2). This is the gap that let a severe defect ship. Every warping
  # test in the suite called the estimators DIRECTLY -- extracting them from the
  # source and invoking them by hand -- so all of them passed while the app
  # could not run landmark registration at all: dance_landmark_warp() was defined
  # inside a reactive observer where its own caller could not see it, the
  # tryCatch turned the lookup failure into NULL, and the app quietly
  # substituted a LINEAR SHIFT. A user choosing "landmark" got a different
  # method's answer with no indication on screen.
  #
  # A direct call cannot detect that, because a direct call supplies the scope
  # the app does not. Only pressing the button does. So this section presses it,
  # once per method, and requires the method that comes back to be the method
  # that was asked for.
  cat("\n-- registration: all three methods, through the button ------------------\n")
  for (meth in c("linear_shift", "parametric", "landmark")) {
    values$warping_results <- NULL
    values$pca_results     <- NULL
    session$setInputs(pca_type = "twpca", warping_method = meth,
                      periodic_shift = TRUE, shift_reference = "mean",
                      parametric_family = "power", param_range = c(0.5, 2),
                      run_analysis = which(c("linear_shift","parametric","landmark") == meth) + 100)
    session$flushReact()
    wr <- values$warping_results
    if (is.null(wr)) {
      fail(sprintf("registration '%s' produced no result at all", meth))
    } else if (!identical(wr$method, meth)) {
      # exactly the P11.2 symptom: asked for one method, got another
      fail(sprintf("asked for '%s' registration, the app ran '%s'", meth, wr$method))
    } else {
      ok(sprintf("registration '%s' ran and reported itself as '%s'", meth, wr$method))
      if (is.null(wr$warp_functions)) {
        fail(sprintf("'%s' returned no warp functions", meth))
      } else if (!all(is.finite(wr$warp_functions))) {
        fail(sprintf("'%s' returned non-finite warp values", meth))
      } else {
        ok(sprintf("'%s' warp functions are finite", meth))
      }
    }
    render(sprintf("warping_plot renders (%s)", meth), output$warping_plot)
    render(sprintf("warping stats render (%s)", meth), output$warping_model_criteria)
  }
  # the resolution note the estimate now carries (P11.4)
  mc <- render("registration stats render", output$warping_model_criteria)
  if (!is.null(mc)) {
    txt <- paste(as.character(mc), collapse = "\n")
    if (!grepl("Resolution of the phase estimate", txt, fixed = TRUE))
      fail("the registration summary does not report its own resolution")
    else ok("the registration summary reports the resolution of the estimate")
  }
  # P13: generate the report while a registration result is present, so the
  # registration subsection is produced from a real run rather than never being
  # exercised at all -- the gap that let two whole sections ship wrong.
  reg_md <- tryCatch(dance_apa_report(values, input, "Registration check"),
                     error = function(e) structure(conditionMessage(e), class = "fckerr"))
  if (inherits(reg_md, "fckerr")) {
    fail(paste("the report errored with a registration result present:", as.character(reg_md)))
  } else {
    rtxt <- paste(reg_md, collapse = "\n")
    if (!grepl("Registration diagnostics", rtxt, fixed = TRUE))
      fail("the report has a registration result but no registration section")
    else ok("the registration section is generated from a real registration")
    for (bad in c("AIC", "BIC", "variance explained by warping"))
      if (grepl(bad, rtxt, fixed = TRUE))
        fail(sprintf("the registration section mentions '%s', which this module does not compute", bad))
    ok("the registration section reports no criterion it cannot justify")
    if (nzchar(Sys.getenv("DANCE_DUMP_REG"))) writeLines(reg_md, Sys.getenv("DANCE_DUMP_REG"))
  }
  values$warping_results <- NULL

  # ------------------------------------------------- cosinor + groups -------
  #
  # AUDIT (P12.2). The cosinor tab was never driven here, and as a result the
  # publication report's cosinor section had been shipped without anyone
  # generating it from a real fit with a grouping variable: it reported pooled
  # descriptives only, and said nothing about between-group differences even
  # when the app had computed them. Same shape of gap as P11.2 -- the machinery
  # existed and nothing pressed the button.
  cat("\n-- cosinor regression, with a grouping variable -------------------------\n")
  session$setInputs(harmonic_data_source = "raw", harmonic_period = 24,
                    n_harmonics = 1, harmonic_trend_type = "none",
                    harmonic_model_selection = FALSE,
                    harmonic_time_var = "_columns_",
                    harmonic_dv_name = "Activity", harmonic_dv_units = "counts/min",
                    harmonic_group_var = "AGEcategory",
                    run_harmonic = 1)
  session$flushReact()
  hm <- values$harmonic_model
  if (is.null(hm)) {
    fail("the cosinor fit produced no model")
  } else {
    ok(sprintf("cosinor fitted: period %s, %d harmonic(s), %d series",
               format(hm$period), hm$n_harmonics, length(hm$individual_fits)))
    if (is.null(hm$group_var_name)) fail("the grouping variable was not recorded")
    else ok(sprintf("grouped by '%s' with %d group fit(s)",
                    hm$group_var_name, length(hm$group_fits)))
    if (is.null(hm$bingham_summary) || is.null(hm$bingham_summary[[1]]))
      fail("no Bingham joint confidence regions were computed")
    else ok(sprintf("Bingham regions computed for %d series, %d with an identified acrophase",
                    hm$bingham_summary[[1]]$n, hm$bingham_summary[[1]]$n_identified))

    # the pairwise comparison the report needs, on BOTH an amplitude and the
    # acrophase -- the acrophase run is what lets the report CHECK Bingham's
    # caveat instead of reciting it
    for (prm in c("acrophase_time_1", "amplitude_1")) {
      session$setInputs(hp_param = prm, hp_correction = "holm",
                        hp_show_ci = TRUE, hp_show_effect_size = TRUE,
                        hp_run = which(c("acrophase_time_1", "amplitude_1") == prm))
      session$flushReact()
      pr <- values$hp_pairwise_results
      if (is.null(pr) || !nrow(pr)) fail(sprintf("pairwise on '%s' produced nothing", prm))
      else ok(sprintf("pairwise on '%s': %d comparison(s), correction '%s'",
                      prm, nrow(pr), values$hp_pairwise_correction))
    }
    if (is.null(values$hp_acrophase_differs))
      fail("the acrophase verdict was not recorded for the Bingham caveat")
    else ok(sprintf("acrophase verdict recorded: groups differ = %s",
                    values$hp_acrophase_differs))
  }
  render("pairwise cosinor results render", output$hp_results)

  # ---- the diagnostics panel's own text and the tau check (P15) -------------
  cat("\n-- harmonic diagnostics: help text and the free-vs-fixed tau check ------\n")
  for (k in c(1L, 2L)) {
    session$setInputs(n_harmonics = k)
    session$flushReact()
    h <- render(sprintf("model-selection help renders (n_harmonics = %d)", k),
                output$harmonic_model_selection_help)
    if (!is.null(h)) {
      htxt <- paste(as.character(h), collapse = " ")
      want <- if (k == 1L) "4 fits per subject" else "8 fits per subject"
      if (!grepl(want, htxt, fixed = TRUE))
        fail(sprintf("the help text does not say '%s'", want))
      else ok(sprintf("the help text reports the real fit count (%s)", want))
      if (grepl("{1,2,3}", htxt, fixed = TRUE) || grepl("9 fits", htxt, fixed = TRUE))
        fail("the help text still describes the pre-P14 grid")
      else ok("the help text no longer describes the old fixed grid")
      if (!grepl("logarithmic", htxt, fixed = TRUE))
        fail("the help text omits the logarithmic trend, which the set fits")
      else ok("the help text names every trend the set fits")
    }
  }
  session$setInputs(n_harmonics = 2)

  # ---- the nested model-selection diagnostic (P14) --------------------------
  # The grid used to be hardcoded as {none, linear, exp_sat} x 1:3: it omitted
  # the `log` trend the UI offers -- so a user who chose it was shown a table
  # their own model was not in, with Akaike weights normalised over a set that
  # excluded it -- and it ran to three harmonics whatever had been selected.
  # Driven here with the user's own example: saturating exponential, 2 harmonics.
  cat("\n-- nested model selection: exp_sat + 2 harmonics -----------------------\n")
  t0 <- Sys.time()
  session$setInputs(n_harmonics = 2, harmonic_trend_type = "exp_sat",
                    harmonic_model_selection = TRUE, run_harmonic = 2)
  session$flushReact()
  el <- as.numeric(Sys.time() - t0, units = "secs")
  ms <- values$harmonic_model$model_selection
  if (is.null(ms)) {
    fail("the model-selection diagnostic produced no table")
  } else {
    ok(sprintf("model set fitted in %.1f s: %d candidates", el, nrow(ms)))
    want <- c("none + H1", "none + H1-H2", "linear + H1", "linear + H1-H2",
              "log + H1", "log + H1-H2", "exp_sat + H1", "exp_sat + H1-H2")
    missing <- setdiff(want, ms$model)
    extra   <- setdiff(ms$model, want)
    if (length(missing)) fail(paste("model set is missing:", paste(missing, collapse = ", ")))
    else ok("every trend the UI offers appears, crossed with H1 and H1-H2")
    if (length(extra)) fail(paste("model set fits models beyond the selection:",
                                  paste(extra, collapse = ", ")))
    else ok("nothing beyond the 2 harmonics selected was fitted")
    # a higher harmonic must never appear without the lower ones
    if (any(grepl("\\+ H[23]$", ms$model)))
      fail("a higher harmonic is being fitted on its own")
    else ok("no harmonic is fitted without the ones below it")
    if (!isTRUE(all.equal(sum(ms$weight), 1)))
      fail(sprintf("Akaike weights sum to %.6f, not 1", sum(ms$weight)))
    else ok("Akaike weights sum to 1 over the candidate set")
    if (!identical(attr(ms, "selected"), "exp_sat + H1-H2"))
      fail(paste("the reported model is not marked; got", attr(ms, "selected")))
    else ok("the specification actually reported is marked in the table")
    if (nzchar(Sys.getenv("DANCE_DUMP_MS"))) {
      con <- file(Sys.getenv("DANCE_DUMP_MS"), "w")
      writeLines(sprintf("%-18s %12s %10s %8s%s", ms$model,
                         sprintf("%.2f", ms$AICc), sprintf("%.2f", ms$dAICc),
                         sprintf("%.3f", ms$weight),
                         ifelse(ms$model == attr(ms, "selected"), "  <-- reported", "")),
                 con)
      close(con)
    }
  }

  # ------------------------------------------------- the APA report --------
  cat("\n-- APA report ----------------------------------------------------------\n")
  session$setInputs(apa_report_title = "Diurnal profiles by age group")
  md <- tryCatch(dance_apa_report(values, input, "Diurnal profiles by age group"),
                 error = function(e) structure(conditionMessage(e), class = "fckerr"))
  if (inherits(md, "fckerr")) {
    fail("dance_apa_report ->", as.character(md))
  } else {
    ok(sprintf("report generated (%d lines)", length(md)))
  if (nzchar(Sys.getenv("DANCE_DUMP_REPORT"))) writeLines(md, Sys.getenv("DANCE_DUMP_REPORT"))
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

    hh <- tryCatch(dance_apa_html(md, "t"), error = function(e)
                   structure(conditionMessage(e), class = "fckerr"))
    if (inherits(hh, "fckerr")) fail("dance_apa_html ->", as.character(hh))
    else {
      h1 <- paste(hh, collapse = "\n")
      if (!grepl("<!DOCTYPE html>", h1, fixed = TRUE)) fail("the HTML has no doctype")
      else if (!grepl("<table>", h1, fixed = TRUE)) fail("no table survived the HTML rendering")
      else if (grepl("| ---:", h1, fixed = TRUE)) fail("a Markdown table separator leaked into the HTML")
      else ok(sprintf("HTML rendering produced %d lines with tables intact", length(hh)))
    }
    # The dump is a convenience for a human reading the rendered report, not an
    # assertion. It used to write to a hard-coded scratch directory that exists
    # only on the machine the test was written on, so on any other machine the
    # whole APA block died at this line and the checks BELOW it never ran.
    # tempdir() exists everywhere, and the path is printed so it can be found.
    dump_dir <- Sys.getenv("DANCE_TEST_DUMP_DIR", unset = tempdir())
    dumped <- tryCatch({
      dir.create(dump_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(md, file.path(dump_dir, "dance_apa_report.md"))
      writeLines(hh, file.path(dump_dir, "dance_apa_report.html"))
      TRUE
    }, error = function(e) FALSE)
    if (dumped) cat(sprintf("   (written to %s)\n", file.path(dump_dir, "dance_apa_report.md")))
    else cat("   (report dump skipped: no writable directory)\n")
  }
  render("apa_report_preview renders", output$apa_report_preview)

  cat("\n")
  if (failures) { cat(sprintf("Reactive smoke test FAILED (%d).\n", failures)); quit(status = 1) }
  cat("Reactive smoke test passed.\n")
}

testServer(shinyApp(ui_obj, server), { NULL })
