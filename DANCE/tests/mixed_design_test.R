# ==============================================================================
# tests/mixed_design_test.R — the mixed (between x within) tab, driven
#
# A design with one between-subject and one within-subject factor cannot be run
# by either existing kernel, and running the between-subjects one on a long file
# is anticonservative: a subject contributing two conditions appears as two rows
# and the permutation treats them as independent curves.
#
# This file drives the new tab through its actual controls on a dataset shaped
# like the one that prompted it: 27 participants x 2 conditions, a 2-level
# between factor, 17 unevenly spaced time points spanning 25 hours, ~1% missing.
#
# Run with:   Rscript tests/mixed_design_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
need <- c("shiny", "shinydashboard", "DT", "plotly", "fda", "ggplot2", "mgcv", "lme4")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) { cat("SKIP: missing", paste(miss, collapse = ", "), "\n"); quit(status = 0) }
suppressMessages({
  library(shiny); library(shinydashboard); library(DT); library(plotly)
  library(fda); library(ggplot2); library(mgcv); library(lme4)
})

app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/06_helpers_mixed.R"), local = e)

# ---- a dataset with the shape that prompted this ----------------------------
set.seed(42)
n_sub <- 27
hrs <- c(9,10,11,12,14,16,18,20,21,22,23,0,2,4,6,8,10)
cum <- 0; for (i in 2:length(hrs))
  cum[i] <- cum[i-1] + ifelse(hrs[i] >= hrs[i-1], hrs[i]-hrs[i-1], 24-hrs[i-1]+hrs[i])
between <- rep(sample(c("g1","g2"), n_sub, TRUE, c(.4,.6)), 2)
subj    <- rep(sprintf("S%02d", seq_len(n_sub)), 2)
within  <- rep(c("A","B"), each = n_sub)
# a real interaction: the within effect is larger in g2
Y <- t(vapply(seq_along(subj), function(i) {
  amp  <- 20 + 10 * (within[i] == "A") + 8 * (between[i] == "g2" & within[i] == "A")
  base <- 40 + rnorm(1, 0, 6)
  base + amp * sin(2*pi*(cum - 4)/24) + rnorm(length(cum), 0, 5)
}, numeric(length(cum))))
Y[sample(length(Y), 11)] <- NA          # ~1% missing, as in the real file

long <- e$dance_mixed_long(Y, cum, subject = subj, between = between, within = within)
chk(nrow(long) > 0 && all(is.finite(long$y)),
    sprintf("long form built: %d observations, %d cells", nrow(long), nlevels(long$cell)),
    "the long form is empty or carries non-finite values")
chk(nlevels(long$subject) == n_sub,
    sprintf("the participant identifier survives reshaping (%d subjects)", nlevels(long$subject)),
    "subjects were lost in the reshaping")
chk(length(e$dance_mixed_check(long)) == 0,
    "the design validates as mixed",
    paste("design rejected:", paste(e$dance_mixed_check(long), collapse = " ")))

# ---- the checks must REFUSE designs that are not mixed ----------------------
b2 <- long; b2$between <- b2$within        # between varies within subject
chk(any(grepl("not a between-subject factor", e$dance_mixed_check(b2))),
    "a 'between' factor that varies within a subject is refused",
    "a within-varying factor was accepted as between-subject")
w1 <- long[long$within == levels(long$within)[1], ]; w1$within <- droplevels(w1$within)
chk(length(e$dance_mixed_check(w1)) > 0,
    "a design with nothing repeated within subject is refused",
    "a non-repeated design was accepted as mixed")

# ---- mixed functional ANOVA -------------------------------------------------
cat("\n-- mixed functional model ----------------------------------------------\n")
t0 <- Sys.time(); r <- e$dance_mixed_fanova(long); el <- as.numeric(Sys.time()-t0, units="secs")
chk(isTRUE(r$ok), sprintf("fitted in %.1f s, deviance explained %.1f%%", el, 100*r$dev_expl),
    paste("the model failed:", r$message))
if (isTRUE(r$ok)) {
  chk(!is.null(r$s_table) && nrow(r$s_table) >= nlevels(long$cell),
      sprintf("a smooth per cell plus the subject term (%d smooth rows)", nrow(r$s_table)),
      "the smooth table does not carry one term per cell")
  chk(any(grepl("subject", rownames(r$s_table))),
      "the per-subject random smooth is in the model",
      "there is no subject term: the repeated measurement is not modelled")
  chk(is.finite(r$aic_delta),
      sprintf("the interaction is testable: AIC delta %+.1f", r$aic_delta),
      "no additive comparison was available, so the interaction cannot be judged")
  # the data were built WITH an interaction, so it should be preferred
  chk(r$aic_delta > 2,
      sprintf("the planted interaction is detected (delta %+.1f favours the full model)", r$aic_delta),
      sprintf("the planted interaction was NOT detected (delta %+.1f)", r$aic_delta))
  cv <- e$dance_mixed_fanova_curves(r, 40)
  chk(!is.null(cv) && nrow(cv) == 40 * nlevels(long$cell),
      "population cell curves are predictable with the subject term excluded",
      "the fitted cell curves could not be produced")
}

# ---- mixed cosinor ----------------------------------------------------------
cat("\n-- mixed cosinor -------------------------------------------------------\n")
t0 <- Sys.time(); rc <- e$dance_mixed_cosinor(long, period = 24); el <- as.numeric(Sys.time()-t0, units="secs")
chk(isTRUE(rc$ok), sprintf("fitted in %.1f s (random rhythm: %s)", el, rc$random_rhythm),
    paste("the mixed cosinor failed:", rc$message))
if (isTRUE(rc$ok)) {
  chk(nrow(rc$cells) == nlevels(long$cell),
      sprintf("MESOR, amplitude and acrophase for each of the %d cells", nrow(rc$cells)),
      "not every cell got rhythm parameters")
  chk(all(is.finite(rc$cells$amplitude_1)) && all(rc$cells$amplitude_1 >= 0),
      "amplitudes are finite and non-negative",
      "an amplitude is non-finite or negative, which a norm cannot be")
  chk(all(rc$cells$acrophase_1 >= 0 & rc$cells$acrophase_1 < 24),
      "acrophases lie inside one period",
      "an acrophase falls outside the period")
  chk(!is.null(rc$tests) && nrow(rc$tests) == 3,
      "all three likelihood-ratio tests were computed",
      "the likelihood-ratio tests did not run")
  if (!is.null(rc$tests)) {
    chk(all(rc$tests$df %% 2 == 0),
        "each test drops a cosine/sine PAIR, so its df is even",
        "a test has odd df, which means a lone cosine or sine was dropped")
    ix <- grep("interaction", rc$tests$term)
    chk(length(ix) == 1 && rc$tests$p[ix] < 0.05,
        sprintf("the planted rhythm interaction is detected: chi2(%d) = %.2f, p = %.4g",
                rc$tests$df[ix], rc$tests$chisq[ix], rc$tests$p[ix]),
        "the planted rhythm interaction was not detected")
  }
  # amplitude should be larger where it was planted larger (A, and more so in g2)
  a_g2A <- rc$cells$amplitude_1[rc$cells$within == "A" & rc$cells$between == "g2"]
  a_g2B <- rc$cells$amplitude_1[rc$cells$within == "B" & rc$cells$between == "g2"]
  chk(a_g2A > a_g2B,
      sprintf("recovers the planted amplitude ordering in g2 (%.1f > %.1f)", a_g2A, a_g2B),
      sprintf("the planted amplitude ordering is not recovered (%.1f vs %.1f)", a_g2A, a_g2B))
}

# ---- the TAB, driven through its own controls -------------------------------
# The kernels above are pure and were called directly. That is exactly the shape
# of test that let P11.2 ship: a direct call supplies the scope the app does not.
# So the tab is also driven end to end, through setInputs and the real outputs.
cat("\n-- the tab, through its controls ---------------------------------------\n")
dance_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir); invisible(NULL)
}
ui_files     <- sort(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"))
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
for (f in ui_files) dance_source(f, envir = globalenv())

server_fn <- function(input, output, session) {
  for (f in server_files) dance_source(f, envir = environment())
  values$data        <- Y
  values$time_labels <- sprintf("%02d:00", hrs)
  values$covariates  <- data.frame(Condition = within, Group = between,
                                   stringsAsFactors = FALSE)
  values$subject_ids <- subj
}

suppressWarnings(shiny::testServer(server_fn, {
  session$setInputs(mixed_between = "Group", mixed_within = "Condition",
                    mixed_real_time = FALSE, mixed_analysis = "fanova",
                    mixed_k_time = 10, mixed_k_subject = 5,
                    mixed_period = 24, mixed_harmonics = 1)
  session$flushReact()

  d <- dance_mixed_frame()
  chk(!is.null(d) && nrow(d) > 0,
      sprintf("the tab builds the long form from the session (%s rows)",
              if (is.null(d)) "0" else nrow(d)),
      "the tab could not build a long form from the loaded data")

  rr <- tryCatch(force(output$mixed_design_summary), error = function(e) e)
  chk(!inherits(rr, "error"), "the live design summary renders before fitting",
      paste("the design summary errored:", conditionMessage(rr)))

  session$setInputs(run_mixed = 1); session$flushReact()
  chk(!is.null(values$mixed_results) && identical(values$mixed_results$kind, "fanova"),
      "the button runs the mixed functional model",
      "pressing the button produced no functional result")
  for (o in c("mixed_results", "mixed_plot")) {
    v <- tryCatch(force(output[[o]]), error = function(e) e)
    chk(!inherits(v, "error"), paste(o, "renders (functional)"),
        paste0(o, " errored: ", if (inherits(v, "error")) conditionMessage(v) else ""))
  }

  session$setInputs(mixed_analysis = "cosinor", run_mixed = 2); session$flushReact()
  chk(!is.null(values$mixed_results) && identical(values$mixed_results$kind, "cosinor"),
      "the button runs the mixed cosinor",
      "switching to the cosinor produced no result")
  for (o in c("mixed_results", "mixed_plot")) {
    v <- tryCatch(force(output[[o]]), error = function(e) e)
    chk(!inherits(v, "error"), paste(o, "renders (cosinor)"),
        paste0(o, " errored: ", if (inherits(v, "error")) conditionMessage(v) else ""))
  }

  # the readout must not silently look like a one-factor result
  txt <- paste(capture.output(print(force(output$mixed_results))), collapse = " ")
  chk(grepl("Mixed cosinor", txt, fixed = TRUE),
      "the readout names the model it fitted",
      "the readout does not say which model produced these numbers")

  # the publication report must carry it, for both kinds
  for (kind in c("cosinor", "fanova")) {
    session$setInputs(mixed_analysis = kind, run_mixed = 10 + (kind == "fanova"))
    session$flushReact()
    md <- tryCatch(dance_apa_report(values, input, "Mixed design"),
                   error = function(e) e)
    if (inherits(md, "error")) {
      chk(FALSE, "", paste("the report errored with a", kind, "result:", conditionMessage(md)))
      next
    }
    rt <- paste(md, collapse = "\n")
    head <- if (kind == "fanova") "Mixed functional model" else "Mixed cosinor"
    chk(grepl(head, rt, fixed = TRUE),
        paste("the report carries the", kind, "results section"),
        paste("the report has a", kind, "result but no section for it"))
    chk(grepl("varied between participants", rt, fixed = TRUE),
        paste("the report states the design in Methods (", kind, ")"),
        paste("the report does not describe the mixed design (", kind, ")"))
    chk(grepl("the *p* values are approximate", rt, fixed = TRUE),
        paste("the report says the p values are approximate (", kind, ")"),
        paste("the report presents approximate p values as exact (", kind, ")"))
    chk(grepl("What these numbers do not establish", rt, fixed = TRUE),
        paste("the section closes with its caveats (", kind, ")"),
        paste("the", kind, "section has no caveat paragraph"))
  }
}))

# ---- the exported script must RUN and reproduce the fit (P18.2) -------------
# The mixed module shipped with no export branch at all, so the script covered
# every family except the newest one. Parsing is not the bar -- a script full of
# the wrong algorithm parses -- so the emitted text is executed in a clean
# session and its numbers compared against the app's.
cat("\n-- the exported script, executed ---------------------------------------\n")
suppressWarnings(shiny::testServer(server_fn, {
  session$setInputs(mixed_between = "Group", mixed_within = "Condition",
                    mixed_real_time = FALSE, mixed_analysis = "cosinor",
                    mixed_period = 24, mixed_harmonics = 1,
                    mixed_k_time = 10, mixed_k_subject = 5, run_mixed = 1)
  session$flushReact()
  chk(!is.null(values$mixed_results), "a mixed cosinor result exists to export",
      "no mixed result to export")

  code <- tryCatch(generate_analysis_code(full = TRUE), error = function(e) e)
  if (inherits(code, "error")) {
    chk(FALSE, "", paste("the generator errored:", conditionMessage(code)))
  } else {
    chk(grepl("12. MIXED DESIGN", code, fixed = TRUE),
        "the script carries a mixed-design section",
        "the exported script has no mixed-design section")
    for (fn in c("dance_mixed_long", "dance_mixed_check", "dance_mixed_cosinor"))
      chk(grepl(paste0(fn, " <- function"), code, fixed = TRUE),
          paste("the script defines", fn, "from the live function"),
          paste("the script calls", fn, "but never defines it"))
    chk(grepl("fitted to the RAW observations", code, fixed = TRUE),
        "the script says the mixed fit uses raw observations",
        "the script does not say which data the mixed model used")

    f <- tempfile(fileext = ".R"); writeLines(code, f)
    chk(!inherits(tryCatch(parse(f), error = function(e) e), "error"),
        "the exported script parses", "the exported script does not parse")

    # Run just the mixed section in a clean session, against the app's numbers.
    # generate_analysis_code() returns ONE string, so it is split first; the
    # section runs from its own header to the end marker, and emit_kernel()
    # writes the kernel definitions inside that range.
    live  <- values$mixed_results
    lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
    i0 <- grep("12. MIXED DESIGN", lines, fixed = TRUE)[1]
    i1 <- grep("END OF ANALYSIS CODE", lines, fixed = TRUE)[1]
    chk(!is.na(i0) && !is.na(i1) && i1 > i0,
        sprintf("the mixed section is locatable in the script (lines %s-%s)", i0, i1),
        "the mixed section could not be located in the script")
    rf <- tempfile(fileext = ".rds"); saveRDS(Y, rf)
    runner <- tempfile(fileext = ".R")
    writeLines(c(
      '.libPaths(c("~/Rlib", .libPaths()))',
      'suppressMessages(library(lme4))',
      sprintf('raw_data <- as.matrix(readRDS("%s"))', rf),
      lines[i0:(i1 - 1)],
      'cat(sprintf("%.10f", mixed_fit$cells$amplitude_1), sep="\n")'), runner)
    out <- suppressWarnings(system2("Rscript", runner, stdout = TRUE, stderr = FALSE))
    got <- suppressWarnings(as.numeric(out[grepl("^[0-9.]+$", out)]))
    chk(length(got) == nrow(live$cells),
        sprintf("the exported script runs and returns %d cell amplitudes", length(got)),
        sprintf("the exported script returned %d values, expected %d",
                length(got), nrow(live$cells)))
    if (length(got) == nrow(live$cells)) {
      d <- max(abs(got - live$cells$amplitude_1))
      chk(d < 1e-6,
          sprintf("exported amplitudes match the app (max diff %.2e)", d),
          sprintf("exported amplitudes differ from the app by %.6f", d))
    }
  }
}))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (bad == 0) "Mixed design tests PASSED" else "Mixed design tests FAILED", ok_n, bad))
quit(status = if (bad == 0) 0 else 1)
