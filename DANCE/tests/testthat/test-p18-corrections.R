# ==============================================================================
# tests/testthat/test-p18-corrections.R — the ninth review
#
#   P18.1  values$mixed_results was created dynamically by its module, never
#          declared in the state, and never cleared by dance_reset_analyses().
#          Load dataset A, run a mixed model, load dataset B: the old result
#          survives, and server/93_apa_report.R reads it.
#
#          Checking that turned up four MORE slots with the same defect, all
#          older and all mine: hp_pairwise_all, hp_pairwise_param,
#          hp_pairwise_correction, hp_acrophase_param and hp_acrophase_differs.
#          Those are worse, because the report reads hp_pairwise_all as its
#          PRIMARY source for cosinor group comparisons and hp_acrophase_differs
#          as the evidence for the Bingham amplitude caution -- so a stale pair
#          could put dataset A's comparisons, and a caution derived from them,
#          into a publication report about dataset B.
#
# The name-by-name fix is not the interesting part. The test below is
# STRUCTURAL: it reads every `values$x <- ` assignment in the server modules and
# requires each derived slot to be cleared, so the next module that forgets
# fails here instead of shipping.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"

state_env <- function() {
  e <- new.env(parent = globalenv())
  # reactiveValues is not available outside Shiny; the reset only assigns NULL
  # to named elements, which a plain environment models exactly.
  e$reactiveValues <- function(...) list2env(list(...), parent = emptyenv())
  eval(parse(file.path(app_dir, "server/00_state.R"), encoding = "UTF-8"), envir = e)
  e
}

# ---- what the reset actually clears -----------------------------------------
reset_clears <- function() {
  src <- readLines(file.path(app_dir, "server/00_state.R"), warn = FALSE)
  st  <- grep("^dance_reset_analyses <- function", src)
  en  <- st - 1 + which(src[st:length(src)] == "}")[1]
  body <- src[st:en]
  body <- body[!grepl("^\\s*#", body)]
  unique(sub(".*values\\$([A-Za-z0-9_]+)\\s*<-.*", "\\1",
             grep("values\\$[A-Za-z0-9_]+\\s*<-\\s*NULL", body, value = TRUE)))
}

# ---- what the modules write --------------------------------------------------
# Everything a module assigns to values$ is either OWNED by the data step or
# DERIVED from it. Derived things must be cleared; owned things must not be.
owned_by_data_step <- c(
  "data", "raw_df", "raw_data", "covariates", "time_labels", "time_values",
  "group_labels", "group_variables", "subject_ids", "fill_status",
  "landmark_points", "selected_vars", "file_name", "data_source",
  "n_subjects", "n_timepoints", "clock_hours", "cumulative_hours",
  "sample_data_used", "import_notes", "id_column", "excluded_rows",
  # kept deliberately across a re-smooth, cleared with the data
  "smooth_data", "fd_obj", "smooth_fit_metrics", "smoothing_avg_metrics",
  # smoothing diagnostics
  "cv_results", "nbasis_gcv", "reml_profile",
  # the import step's own bookkeeping
  "uploaded_data", "selected_group_vars", "time_numeric", "time_clock",
  # UI-only: which curve the viewer is looking at, not a result
  "selected_curve",
  "current_tab", "last_error", "status_message"
)

module_writes <- function() {
  fs <- list.files(file.path(app_dir, "server"), pattern = "[.]R$", full.names = TRUE)
  fs <- fs[!grepl("00_state[.]R$", fs)]
  out <- character(0)
  for (f in fs) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[!grepl("^\\s*#", ln)]
    m <- regmatches(ln, gregexpr("values\\$[A-Za-z0-9_]+\\s*<<?-", ln))
    out <- c(out, unlist(m))
  }
  unique(sub("values\\$([A-Za-z0-9_]+)\\s*<<?-", "\\1", out))
}

test_that("P18.1: every derived result is cleared when the data changes", {
  cleared <- reset_clears()
  written <- module_writes()
  derived <- setdiff(written, owned_by_data_step)
  missing <- setdiff(derived, cleared)
  # Name the offenders: a bare count tells the next person nothing.
  expect_identical(missing, character(0),
                   info = paste("written by a module, never cleared:",
                                paste(missing, collapse = ", ")))
})

test_that("P18.1: the slots this review found are declared and cleared", {
  cleared <- reset_clears()
  for (nm in c("mixed_results", "hp_pairwise_all", "hp_pairwise_param",
               "hp_pairwise_correction", "hp_acrophase_param",
               "hp_acrophase_differs"))
    expect_true(nm %in% cleared, info = nm)
  # and declared in the state, so they exist before a module invents them
  decl <- paste(readLines(file.path(app_dir, "server/00_state.R"), warn = FALSE),
                collapse = "\n")
  for (nm in c("mixed_results", "hp_pairwise_all", "hp_acrophase_differs"))
    expect_true(grepl(paste0(nm, "\\s*=\\s*NULL"), decl), info = nm)
})

test_that("P18.1: the reset actually nulls them, run for real", {
  e <- state_env()
  v <- new.env(parent = emptyenv())
  slots <- c("mixed_results", "hp_pairwise_all", "hp_acrophase_differs",
             "pca_results", "harmonic_model", "fanova_results")
  for (nm in slots) assign(nm, list(ok = TRUE), envir = v)
  e$dance_reset_analyses(v)
  for (nm in slots) expect_null(get(nm, envir = v), info = nm)
})

test_that("P18.1: a re-smooth keeps the mixed result but clears the rest", {
  # The one deliberate asymmetry: mixed models are fitted to values$data, which
  # a re-smooth does not touch, so they stay valid. Everything computed from the
  # SMOOTHED curves does not.
  e <- state_env()
  v <- new.env(parent = emptyenv())
  for (nm in c("mixed_results", "pca_results", "harmonic_model", "hp_pairwise_all"))
    assign(nm, list(ok = TRUE), envir = v)
  e$dance_reset_analyses(v, keep_smoothing = TRUE)
  expect_false(is.null(get("mixed_results", envir = v)))
  for (nm in c("pca_results", "harmonic_model", "hp_pairwise_all"))
    expect_null(get(nm, envir = v), info = nm)

  # and smoothing itself never rewrites values$data, which is what makes that
  # asymmetry sound rather than convenient
  sm <- readLines(file.path(app_dir, "server/20_smoothing.R"), warn = FALSE)
  sm <- sm[!grepl("^\\s*#", sm)]
  expect_length(grep("values\\$data\\s*<-", sm), 0)
})

# ============================================ P18.2 the export gap ============
test_that("P18.2: the generator emits the mixed kernels", {
  ex <- paste(readLines(file.path(app_dir, "server/90_export.R"), warn = FALSE),
              collapse = "\n")
  expect_true(grepl("12. MIXED DESIGN", ex, fixed = TRUE))
  for (k in c("dance_mixed_long", "dance_mixed_check", "dance_mixed_balance",
              "dance_mixed_fanova", "dance_mixed_fanova_curves", "dance_mixed_cosinor"))
    expect_true(grepl(paste0('emit_kernel("', k, '")'), ex, fixed = TRUE), info = k)
  # the script must say which data the mixed model saw, since the rest of the
  # pipeline runs on smoothed curves
  expect_true(grepl("fitted to the RAW observations", ex, fixed = TRUE))
  # tests/mixed_design_test.R is what proves the emitted script RUNS and matches
})

# ==================================== P18.3 session restore and versions ======
test_that("P18.3: the mixed controls are restored with the result", {
  ss <- paste(readLines(file.path(app_dir, "server/91_session.R"), warn = FALSE),
              collapse = "\n")
  for (nm in c("mixed_between", "mixed_within", "mixed_analysis", "mixed_real_time",
               "mixed_k_time", "mixed_k_subject", "mixed_period", "mixed_harmonics"))
    expect_true(grepl(paste0('"', nm, '"'), ss, fixed = TRUE), info = nm)
})

test_that("P18.3: lme4 is in the version stamp, because it fits the cosinor", {
  ss <- paste(readLines(file.path(app_dir, "server/91_session.R"), warn = FALSE),
              collapse = "\n")
  st <- sub(".*dance_package_versions <- function\\(\\) \\{(.*?)\\}.*", "\\1", ss)
  expect_true(grepl('"lme4"', st, fixed = TRUE))
  # every package a kernel actually calls must be stamped
  expect_true(grepl('"mgcv"', st, fixed = TRUE))
})

# ================================ P18.4 which data the mixed models use =======
test_that("P18.4: the raw-observations choice is stated where it is read", {
  # P19 moved the mixed analyses into the Functional ANOVA and Cosinor tabs, so
  # the statement has to appear where those users read it, not in a tab that no
  # longer exists. A guard that keeps reading the old file goes VACUOUS, which is
  # the failure mode this suite exists to prevent.
  ui <- paste(readLines(file.path(app_dir, "ui/50_fanova.R"), warn = FALSE), collapse = "\n")
  rp <- paste(readLines(file.path(app_dir, "server/93_apa_report.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("fitted to the <b>raw observations</b>", ui, fixed = TRUE))
  expect_true(grepl("RAW observations rather than to the", rp, fixed = TRUE))
  expect_false(file.exists(file.path(app_dir, "ui/55_mixed.R")))

  # and every entry point really reads values$data, not values$smooth_data
  for (f in c("server/51b_fanova_mixed_views.R", "server/73_cosinor_pairwise.R")) {
    src <- paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
    if (!grepl("dance_mixed_long(", src, fixed = TRUE)) next
    expect_true(grepl("dance_mixed_long(values$data", src, fixed = TRUE), info = f)
    expect_false(grepl("dance_mixed_long(values$smooth_data", src, fixed = TRUE), info = f)
  }
})

# ================================== P18.5 ML for the model comparison =========
test_that("P18.5: the interaction comparison is refitted by ML", {
  src <- paste(readLines(file.path(app_dir, "server/06_helpers_mixed.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl('mgcv::bam(form, data = d, method = "ML")', src, fixed = TRUE))
  expect_true(grepl("aic_basis", src, fixed = TRUE))
  # the reported model keeps its fREML fit: that is the better basis for the
  # smoothing parameters, and hence for the curves that get plotted
  expect_true(grepl('method = "fREML"', src, fixed = TRUE))
  # and the readout must not dress a continuous weight of evidence as a verdict.
  # P19 extracted that readout into the shared helper, so both hosts describe a
  # fit the same way; the guard follows it there.
  ro <- paste(readLines(file.path(app_dir, "server/06_helpers_mixed.R"), warn = FALSE),
              collapse = "\n")
  expect_true(grepl("not a test of a null", ro, fixed = TRUE))
  expect_false(grepl("the additive model is worse:", ro, fixed = TRUE))
})

# =================================== P18.6 the documentation caught up ========
test_that("P18.6: the README describes DANCE, not the app it replaced", {
  rd <- paste(readLines(file.path(app_dir, "README.md"), warn = FALSE), collapse = "\n")
  expect_true(grepl("Direct Analysis of multivariate", rd, fixed = TRUE))
  expect_false(grepl("**F**unctional data analysis, **C**ircadian regression", rd, fixed = TRUE))
  # the claim that stopped being true several hundred commits ago
  expect_false(grepl("Everything downstream of smoothing is the original code",
                     rd, fixed = TRUE))
  # The tab list must not skip a number. One entry legitimately covers two tabs
  # in its own heading ("4. ... and 5. ..."), so a skipped number is allowed
  # only when the heading before it names it -- which is exactly the case the
  # reviewer's 10-to-12 jump was NOT.
  ls   <- strsplit(rd, "\n")[[1]]
  head <- grep("^[0-9]+\\. \\*\\*", ls, value = TRUE)
  nums <- as.integer(sub("^([0-9]+)\\..*", "\\1", head))
  expect_identical(nums, sort(nums))
  for (i in seq_along(nums)[-1]) {
    gap <- nums[i] - nums[i - 1]
    if (gap == 1) next
    for (k in (nums[i - 1] + 1):(nums[i] - 1))
      expect_true(grepl(paste0("\\b", k, "\\. "), head[i - 1]),
                  info = paste("tab", k, "is missing from the list"))
  }
  # P19 folded the mixed analyses into the Functional ANOVA and Cosinor entries
  # rather than giving them a tab of their own, so that is where the README has
  # to describe them.
  expect_true(grepl("mixed\n   (between", rd, fixed = TRUE) ||
              grepl("mixed (between", rd, fixed = TRUE) ||
              grepl("**and mixed", rd, fixed = TRUE))
  # P20/R5: the README used to call this "Exact permutation (default)". Only the
  # WITHIN-participant scheme is exact; the two that relabel across groups are
  # exact under exchangeability and asymptotic otherwise, and the README now has
  # to say which is which rather than claim the stronger thing for all three.
  expect_true(grepl("Permutation (default)", rd, fixed = TRUE))
  expect_false(grepl("Exact permutation (default)", rd, fixed = TRUE))
  expect_true(grepl("How exact is \"exact\"?", rd, fixed = TRUE))
  expect_true(grepl("Welch-type studentised", rd, fixed = TRUE))
  expect_true(grepl("Population-mean cosinor", rd, fixed = TRUE))
  expect_true(grepl("| `lme4` |", rd, fixed = TRUE))
})
