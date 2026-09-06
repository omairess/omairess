# ==============================================================================
# tests/testthat/test-p6-runtime-bugs.R
#
# Seven bugs found by RUNNING the app on real data. Two of them are mine, from
# the previous two rounds. All seven are the kind that a parse test, a source
# guard and even a numerical round-trip all miss, because they live in the
# reactive layer or in a display cap.
#
#   P6.1  FoSR GAM: "object 'j' not found"  (regression, P5.8)
#   P6.2  fPCA-ANOVA: "$ operator is invalid for atomic vectors"
#   P6.3  fPCA capped every display at 3 components, whatever was extracted
#   P6.4  the "Components to show" slider drew nonsense tick labels
#   P6.5  the FoSR observed-data plot showed one group out of four
#   P6.6  pairwise p-values were all identical (the Monte Carlo floor)
#   P6.7  three harmonic-regression defaults were the wrong way round
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"

code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  strip <- function(l) {
    ch <- strsplit(l, "", fixed = TRUE)[[1]]
    if (!length(ch)) return(l)
    inq <- ""; out <- character(0)
    for (i in seq_along(ch)) {
      c1 <- ch[i]; esc <- i > 1 && ch[i - 1] == "\\"
      if (nzchar(inq)) { if (c1 == inq && !esc) inq <- "" }
      else if (c1 %in% c("\"", "'")) inq <- c1
      else if (c1 == "#") break
      out <- c(out, c1)
    }
    paste(out, collapse = "")
  }
  paste(vapply(ln, strip, character(1), USE.NAMES = FALSE), collapse = "\n")
}

# ================================================== P6.1 the GAM kernel runs ==
test_that("P6.1: the GAM kernel fits, with a factor and a non-syntactic name", {
  skip_if_not_installed("mgcv")
  suppressPackageStartupMessages(library(mgcv))
  source(file.path(app_dir, "server/07b_helpers_fosr_gam.R"), local = TRUE)

  set.seed(51)
  n <- 40; nt <- 12
  df <- data.frame(rnorm(n, 40, 8), factor(rep(c("A", "B", "C", "D"), each = n / 4)),
                   check.names = FALSE)
  names(df) <- c("Age (years)", "AGEcategory")
  Y <- matrix(rnorm(n * nt), n, nt) + df[[1]] %o% seq(0, 0.05, length.out = nt)

  fit <- fck_fit_fosr_gam(Y, df, c("Age (years)", "AGEcategory"))
  # one row per numeric predictor plus one per NON-REFERENCE factor level
  expect_equal(nrow(fit$beta.hat), 5L)
  expect_setequal(rownames(fit$beta.hat),
                  c("(Intercept)", "Age (years)",
                    "AGEcategoryB", "AGEcategoryC", "AGEcategoryD"))
  expect_true(all(is.finite(fit$beta.hat)))
  # P5.9: unknown inference is NA, never zero
  expect_true(all(is.na(fit$beta.se)))
  expect_identical(fit$inference, "none")

  # the loop index that broke it must not shadow the spline argument either
  src <- code_of("server/07b_helpers_fosr_gam.R")
  expect_false(grepl("for (k in seq_along(preds))", src, fixed = TRUE))
  expect_false(grepl("preds[j]", substr(src, 1, regexpr("for (j in seq_along", src, fixed = TRUE)),
                     fixed = TRUE))
})

# ============================================ P6.2 the partial-matching trap ==
test_that("P6.2: `$` on a list partial-matches, and NULL assignment deletes", {
  # the exact mechanism, pinned so nobody re-derives it from an error message
  res <- list(k = 3)
  res$warped <- FALSE
  res$warping_method <- if (res$warped) "x" else NULL   # DELETES, does not store
  res$warp <- NULL                                      # also a no-op
  expect_false("warping_method" %in% names(res))
  expect_false("warp" %in% names(res))
  # only one name now begins with "warp", so `$` partial-matches it
  expect_identical(res$warp, FALSE)
  expect_false(is.null(res$warp))
  expect_error(res$warp$k, "atomic")
  # exact indexing does not
  expect_null(res[["warp"]])
})

test_that("P6.2: the fPCA-ANOVA report survives an unwarped run", {
  src <- code_of("server/42_fpca_anova.R")
  # the slots are always present, so no dangling prefix exists
  expect_true(grepl('res$warping_method <- if (res$warped) values$warping_results$method else NA_character_',
                    src, fixed = TRUE))
  expect_true(grepl("res$warp <- NA", src, fixed = TRUE))
  expect_true(grepl("res$warp_note <- NA_character_", src, fixed = TRUE))
  # and every read of the warp family is exact
  expect_false(grepl("res$warp$", src, fixed = TRUE))
  expect_true(grepl('is.list(res[["warp"]])', src, fixed = TRUE))
  expect_false(grepl("} else if (!is.null(res$warp)) {", src, fixed = TRUE))

  # the branch decision itself, on an unwarped result
  res <- list(k = 3, warped = FALSE,
              warping_method = NA_character_, warp = NA, warp_note = NA_character_)
  took <- if (!is.null(res[["warp_note"]]) && !all(is.na(res[["warp_note"]]))) "note"
          else if (is.list(res[["warp"]])) "phase" else "skipped"
  expect_identical(took, "skipped")
})

# ============================================ P6.3 the component count ========
test_that("P6.3: length() of an fd object is 3, so it cannot count harmonics", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)

  set.seed(52)
  n <- 40; nt <- 24; tp <- seq(0, 1, length.out = nt)
  Y <- t(vapply(seq_len(n), function(i)
    50 + 8 * cos(2 * pi * (tp - 0.6)) + rnorm(nt, 0, 1), numeric(nt)))
  b  <- create.bspline.basis(c(0, 1), nbasis = 12)
  fd <- smooth.basis(tp, t(Y), fdPar(b, 2, 1e-4))$fd

  for (nh in c(2, 3, 5, 7)) {
    p <- pca.fd(fd, nharm = nh)
    expect_equal(length(p$harmonics), 3L)          # the trap: always 3
    expect_equal(fck_n_harmonics(p), nh)           # the fix
    expect_equal(ncol(p$scores), nh)
  }
  expect_equal(fck_n_harmonics(NULL), 0L)
  expect_equal(fck_n_harmonics(list()), 0L)
})

test_that("P6.3: no display path counts components with length(harmonics)", {
  src <- code_of("server/40_fpca.R")
  expect_false(grepl("length(pca_res$harmonics)", src, fixed = TRUE))
  expect_false(grepl("length(pr$harmonics)", src, fixed = TRUE))
  # and the summary no longer stops at three
  expect_false(grepl("for(i in 1:min(3, length(pca_res$varprop)))", src, fixed = TRUE))
  expect_true(grepl("nk <- min(fck_n_harmonics(pca_res), length(pca_res$varprop))",
                    src, fixed = TRUE))
  expect_gte(length(gregexpr("fck_n_harmonics(", src, fixed = TRUE)[[1]]), 5L)
})

# ============================================ P6.4 the slider ================
test_that("P6.4: the slider update sends a complete specification", {
  src <- code_of("server/40_fpca.R")
  expect_true(grepl('min = 1, max = max(1L, n_avail), step = 1', src, fixed = TRUE))
  expect_false(grepl('max = n_avail, value = min(cur, n_avail))', src, fixed = TRUE))
})

# ============================================ P6.5 the stratified subsample ===
test_that("P6.5: the FoSR observed plot no longer takes the first N rows", {
  src <- code_of("server/70_fosr.R")
  expect_false(grepl("subset_idx <- 1:min(nrow(Y), 200)", src, fixed = TRUE))
  expect_true(grepl("idx_by <- split(seq_len(n_all), f)[lv]", src, fixed = TRUE))
  expect_true(grepl("sample_note", src, fixed = TRUE))
  # and it uses the app palette, not a fourth one
  expect_false(grepl('scale_color_brewer(palette = "Set1")', src, fixed = TRUE))
  expect_true(grepl("fck_group_ramp(nlevels(df_plot$Color))", src, fixed = TRUE))
})

test_that("P6.5: the stratified draw reaches every level on a sorted sample", {
  # the real shape: sorted by group, one level >> the cap
  n_by <- c(YOUTH = 654, ADULT = 410, MIDDLE_AGE = 181, ELDERLY = 59)
  col_all <- factor(rep(names(n_by), n_by), levels = names(n_by))
  n_all <- length(col_all); max_curves <- 300L

  # the OLD rule
  old_idx <- 1:min(n_all, 200)
  expect_equal(nlevels(droplevels(col_all[old_idx])), 1L)   # one group, as reported

  # the NEW rule, transcribed from the module
  f <- as.factor(col_all); lv <- levels(droplevels(f))
  idx_by <- split(seq_len(n_all), f)[lv]
  take <- pmax(5L, round(max_curves * lengths(idx_by) / n_all))
  take <- pmin(take, lengths(idx_by))
  while (sum(take) > max_curves) {
    big <- which.max(take); if (take[big] <= 5L) break; take[big] <- take[big] - 1L
  }
  set.seed(1)
  new_idx <- sort(unlist(Map(function(ii, k) if (length(ii) <= k) ii else sample(ii, k),
                             idx_by, take), use.names = FALSE))
  expect_equal(nlevels(droplevels(col_all[new_idx])), 4L)
  expect_lte(length(new_idx), max_curves)
  expect_true(all(table(droplevels(col_all[new_idx])) >= 5))
})

# ============================================ P6.6 the Monte Carlo floor ======
test_that("P6.6: the summary reports the smallest attainable p", {
  src <- code_of("server/50_fanova.R")
  expect_true(grepl(".floor_p <- if (is.finite(.B) && .B > 0) 1 / (.B + 1) else NA_real_",
                    src, fixed = TRUE))
  expect_true(grepl("sit AT that floor", src, fixed = TRUE))
  ui <- code_of("ui/51_posthoc.R")
  expect_true(grepl("value = 5000, min = 500, max = 50000", ui, fixed = TRUE))
  expect_false(grepl("value = 200, min = 100, max = 5000", ui, fixed = TRUE))
  # the arithmetic the note explains
  expect_equal(1 / (100 + 1), 0.009901, tolerance = 1e-6)   # the reported value
  expect_equal(1 / (5000 + 1), 2e-4, tolerance = 1e-5)
})

# ============================================ P6.7 defaults ==================
test_that("P6.7: the harmonic tab's defaults are the ones you would want", {
  ui <- code_of("ui/72_harmonic.R")
  expect_true(grepl('selected = "raw")', ui, fixed = TRUE))
  expect_true(grepl('selected = "first_observation")', ui, fixed = TRUE))
  expect_true(grepl('checkboxInput("harmonic_show_data", "Show Raw Data Points", FALSE)',
                    ui, fixed = TRUE))
  # the server fallbacks must agree, or a session that never touched the control
  # runs a different model from the one the UI shows
  srv <- code_of("server/72_harmonic.R")
  expect_true(grepl('input$harmonic_data_source %||% "raw"', srv, fixed = TRUE))
  expect_true(grepl('input$harmonic_time_origin %||% "first_observation"', srv, fixed = TRUE))
  expect_false(grepl('input$harmonic_data_source %||% "smoothed"', srv, fixed = TRUE))
  expect_false(grepl('input$harmonic_time_origin %||% "midnight"', srv, fixed = TRUE))
})

# ============================================ P7.1 the dead control ==========
test_that("P7.1: the pairwise permutation control is not overridden", {
  src <- code_of("server/50_fanova.R")
  # the override, and the "Fallback" comment that described the LIVE branch as
  # the dead one
  expect_false(grepl("n_perm_to_use <- if(!is.null(values$fanova_results$n_permutations))",
                     src, fixed = TRUE))
  expect_true(grepl("n_perm_to_use <- suppressWarnings(as.integer(input$pairwise_permutations))",
                    src, fixed = TRUE))
  # both counts travel with the result so a difference is visible
  expect_true(grepl("values$pairwise_results$omnibus_permutations <- omnibus_perm",
                    src, fixed = TRUE))
  expect_true(grepl("the omnibus fANOVA used %d", src, fixed = TRUE))
})

test_that("P7.1: a correction floor above alpha makes significance impossible", {
  # m/(B+1) is the smallest adjusted p Bonferroni can produce. If that exceeds
  # alpha, nothing can be significant however strong the effect -- the app now
  # refuses the run rather than reporting a null result.
  m <- 6; alpha <- 0.05
  expect_gt(m / (100 + 1), alpha)      # B = 100: impossible
  expect_lt(m / (1000 + 1), alpha)     # B = 1000: fine
  expect_equal(ceiling(m / alpha), 120)
  src <- code_of("server/50_fanova.R")
  expect_true(grepl("smallest attainable", src, fixed = TRUE))
  expect_true(grepl("Raise B to at least %d", src, fixed = TRUE))
})

# ============================================ P7.2 two more dead controls =====
test_that("P7.2: every UI control is read somewhere in the server", {
  # The sweep that found max_landmarks and display_options. It would NOT have
  # found P7.1 -- that reference existed, in a branch that never ran -- which is
  # why tests/reactive_smoke_test.R also asserts a distinctive value reaches the
  # result. Both checks are needed; neither subsumes the other.
  ids <- unique(unlist(lapply(
    list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"),
    function(f) {
      s <- paste(readLines(f, warn = FALSE), collapse = "\n")
      c(regmatches(s, gregexpr('(?<=Input\\(")[A-Za-z0-9_]+', s, perl = TRUE))[[1]],
        regmatches(s, gregexpr('(?<=actionButton\\(")[A-Za-z0-9_]+', s, perl = TRUE))[[1]])
    })))
  srv <- paste(unlist(lapply(
    list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"),
    readLines, warn = FALSE)), collapse = "\n")
  unread <- ids[!vapply(ids, function(id)
    grepl(paste0("input\\$", id, "\\b"), srv), logical(1))]
  expect_identical(sort(unread), character(0))
})

test_that("P7.2: max_landmarks and display_options now do something", {
  src <- code_of("server/40_fpca.R")
  expect_true(grepl("cap <- suppressWarnings(as.integer(input$max_landmarks))",
                    src, fixed = TRUE))
  expect_true(grepl("Maximum of %d landmarks reached", src, fixed = TRUE))
  expect_true(grepl(".opt <- input$display_options", src, fixed = TRUE))
  expect_true(grepl('show_curves <- is.null(.opt) || "curves" %in% .opt', src, fixed = TRUE))
  expect_true(grepl("if(show_means && n_display > 0) {", src, fixed = TRUE))
})
