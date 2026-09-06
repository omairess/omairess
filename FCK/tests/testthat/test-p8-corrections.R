# ==============================================================================
# tests/testthat/test-p8-corrections.R
#
# A fifth review. Four findings, all of which held.
#
#   P8.1  GAM PREDICTION still used the user's column names against a model
#         fitted on x1..xp -- the same defect as P5.1, in the branch P5.1 did
#         not touch. My omission.
#   P8.2  periodic-shift warp amplitude and velocity were representation
#         artefacts, and they fed the group ANOVA on warping parameters
#   P8.3  the diagnostics compared an mgcv REML lambda with an fda CV lambda by
#         ratio, and offered the REML one as a smoothing factor
#   P8.4  the CV panel's estimand was mislabelled
#   P8.5  the warped-PCA UI still announced EFDA / variance decomposition /
#         model-selection criteria over output that says none exist
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

# ========================================================= P8.1 GAM predict ==
test_that("P8.1: a GAM fit can be predicted from, with a non-syntactic name", {
  skip_if_not_installed("mgcv")
  suppressPackageStartupMessages(library(mgcv))
  source(file.path(app_dir, "server/07b_helpers_fosr_gam.R"), local = TRUE)
  e <- new.env(); eval(parse(file.path(app_dir, "server/02_helpers_gam.R"),
                             encoding = "UTF-8"), envir = e)
  build_gam_pred_df <- get("build_gam_pred_df", envir = e)

  set.seed(81)
  n <- 40; nt <- 12
  df <- data.frame(rnorm(n, 40, 8), factor(rep(c("A", "B"), each = n / 2)),
                   check.names = FALSE)
  names(df) <- c("Age (years)", "AGEcategory")
  Y <- matrix(rnorm(n * nt), n, nt) + df[[1]] %o% seq(0, 0.05, length.out = nt)
  m <- fck_fit_fosr_gam(Y, df, c("Age (years)", "AGEcategory"))

  # the model refers to the internal names, not the user's
  expect_setequal(setdiff(all.vars(formula(m$gam_obj)), c("Y_val", "time")),
                  m$gam_model_names)

  iv <- list(`Age (years)` = 40, AGEcategory = "B")
  pd <- build_gam_pred_df(seq(0, 1, length.out = nt), m$gam_predictors,
                          m$gam_model_names, m$gam_long_data,
                          m$gam_factor_levels, iv)
  expect_true(all(m$gam_model_names %in% names(pd)))
  expect_false(any(m$gam_predictors %in% names(pd)))
  p <- predict(m$gam_obj, newdata = pd, se.fit = TRUE)
  expect_true(all(is.finite(p$fit)))
  expect_true(all(is.finite(p$se.fit)))

  # the frame the OLD helper built cannot be predicted from -- the defect
  old_pd <- data.frame(time = seq(0, 1, length.out = nt))
  old_pd[["Age (years)"]] <- rep(40, nt)
  old_pd[["AGEcategory"]] <- factor(rep("B", nt), levels = c("A", "B"))
  expect_error(suppressWarnings(predict(m$gam_obj, newdata = old_pd)), "not found")

  # and the mapping is a required argument now
  expect_error(build_gam_pred_df(seq(0, 1, length.out = nt), m$gam_predictors,
                                 NULL, m$gam_long_data, m$gam_factor_levels, iv),
               "model_names is required")
})

test_that("P8.1: every caller passes the mapping", {
  src <- code_of("server/70_fosr.R")
  expect_equal(length(gregexpr("build_gam_pred_df(", src, fixed = TRUE)[[1]]), 3L)
  expect_equal(length(gregexpr("mod$gam_model_names", src, fixed = TRUE)[[1]]), 3L)
  expect_false(grepl("build_gam_pred_df(time_vec, preds, long_data", src, fixed = TRUE))
})

# ================================================== P8.2 periodic warp metric =
test_that("P8.2: RMS distance from the identity is an artefact on a wrapped shift", {
  tp <- seq(0, 1, length.out = 100)
  rms <- function(h) sqrt(mean((h - tp)^2))
  # a periodic shift of ZERO already reads as a 0.1 deformation, because the
  # wrap sends the last grid point back to the first
  expect_equal(rms((tp - 0) %% 1), 0.1, tolerance = 1e-9)
  # ... which is exactly what a REAL 0.1 shift measures non-periodically
  expect_equal(rms(tp - 0.1), 0.1, tolerance = 1e-9)
  # and it is not even monotone in the true displacement
  expect_gt(rms((tp - 0.05) %% 1), rms((tp - 0.10) %% 1) - 0.09)
  # the velocity variance is the same huge number at every shift, zero included
  vv <- function(h) var(c(diff(h) / diff(tp), 0))
  expect_equal(vv((tp - 0) %% 1), vv((tp - 0.1) %% 1), tolerance = 1e-6)
  expect_gt(vv((tp - 0) %% 1), 50)
  # the replacement: shortest circular displacement, exact and monotone
  circ <- function(s, span = 1) abs(((s + span / 2) %% span) - span / 2)
  expect_equal(circ(0), 0)
  expect_equal(circ(0.1), 0.1)
  expect_equal(circ(-0.15), 0.15)
  expect_equal(circ(0.9), 0.1)     # 0.9 forward is 0.1 backward
})

test_that("P8.2: the warping parameters for a shift come from the shift itself", {
  source(file.path(app_dir, "server/09_helpers_pcanova.R"), local = TRUE)
  tp <- seq(0, 1, length.out = 100)
  s  <- c(0, 0.05, 0.10, 0.25, -0.15)
  mk <- function(periodic) {
    h <- vapply(s, function(z) if (periodic) (tp - z) %% 1 else tp - z, numeric(length(tp)))
    list(method = "linear_shift", shifts = s, warp_functions = h, time_points = tp,
         boundary = if (periodic) "periodic wrap" else "constant extrapolation")
  }
  for (per in c(TRUE, FALSE)) {
    wp <- fck_warp_params(mk(per))
    # the artefact column is gone; the signed shift is what is tested
    expect_identical(names(wp), "shift")
    expect_equal(wp$shift, s)
  }
  # an endpoint-preserving warp still measures its distance from the identity
  h <- vapply(c(1, 1.3, 2), function(a) tp^a, numeric(length(tp)))
  wp2 <- fck_warp_params(list(method = "parametric", warp_functions = h,
                              time_points = tp))
  expect_identical(names(wp2), "warp_amplitude")
  expect_equal(wp2$warp_amplitude[1], 0, tolerance = 1e-12)   # identity
  expect_gt(wp2$warp_amplitude[3], wp2$warp_amplitude[2])
})

test_that("P8.2: the module measures shift intensity in the right geometry", {
  src <- code_of("server/40_fpca.R")
  expect_true(grepl('if (identical(geom, "shift")) {\n          warp_amplitude[i]    <- phase_displacement[i]',
                    src, fixed = TRUE))
  expect_true(grepl("warp_velocity_var[i] <- if (isTRUE(periodic)) NA_real_ else",
                    src, fixed = TRUE))
  pc <- code_of("server/09_helpers_pcanova.R")
  expect_true(grepl('periodic <- identical(warping_results$boundary, "periodic wrap")',
                    pc, fixed = TRUE))
  expect_true(grepl("abs(((as.numeric(sh) + span / 2) %% span) - span / 2)", pc, fixed = TRUE))
})

# ============================================ P8.3 / P8.4 the diagnostics ====
test_that("P8.3: no ratio is taken between an mgcv lambda and an fda lambda", {
  src <- code_of("server/30_diagnostics.R")
  expect_false(grepl("ratio_min <- cv_lambda / reml_lambda", src, fixed = TRUE))
  expect_false(grepl("ratio_1se <- cv_lambda_1se / reml_lambda", src, fixed = TRUE))
  expect_false(grepl("REML and CV agree well on optimal smoothing", src, fixed = TRUE))
  expect_false(grepl("CV prefers more smoothing than REML", src, fixed = TRUE))
  # the REML lambda is no longer offered as a smoothing factor
  expect_false(grepl("sf_reml <- -log10", src, fixed = TRUE))
  expect_true(grepl("NOT transferable", src, fixed = TRUE))
  expect_true(grepl("No ratio is reported", src, fixed = TRUE))
})

test_that("P8.4: the CV panel names the estimand it actually optimises", {
  src <- code_of("server/30_diagnostics.R")
  expect_true(grepl("Between-subject population-curve prediction CV", src, fixed = TRUE))
  expect_false(grepl("For prediction tasks: lambda", src, fixed = TRUE))
  # the CV really does smooth the TRAINING-GROUP MEAN, which is why
  expect_true(grepl("train_mean <- colMeans(train_data, na.rm = TRUE)", src, fixed = TRUE))
  expect_true(grepl("fd_train <- smooth.basis(time_points[valid_train],", src, fixed = TRUE))
  # ... and it is on the fda scale, so its smoothing factor IS transferable
  expect_true(grepl("fdParobj <- fdPar(basis, 2, lambda)", src, fixed = TRUE))
})

# ============================================ P8.5 the UI stops contradicting =
test_that("P8.5: the warped-PCA UI no longer promises what the server refuses", {
  ui <- code_of("ui/41_results.R")
  for (bad in c("based on EFDA methodology", "Variance Decomposition",
                "Model Selection Criteria", "Warping Fit Statistics"))
    expect_false(grepl(bad, ui, fixed = TRUE), info = bad)
  expect_true(grepl("Registration Diagnostics", ui, fixed = TRUE))
  expect_true(grepl("Pre/post registration dispersion", ui, fixed = TRUE))
  expect_true(grepl("MAGNITUDE OF THE TRANSFORMATION", ui, fixed = TRUE))
})

test_that("P8.5: the README summary matches its own detailed section", {
  rm <- paste(readLines(file.path(app_dir, "README.md"), warn = FALSE), collapse = "\n")
  expect_false(grepl("automatic REML or manual lambda", rm, fixed = TRUE))
  expect_true(grepl("automatic GCV-selected lambda or manual lambda", rm, fixed = TRUE))
  expect_false(grepl("warping fit statistics and variance", rm, fixed = TRUE))
  # the one surviving mention is the detailed section saying it is NOT one
  expect_true(grepl("It is **not**\nan amplitude/phase variance decomposition", rm, fixed = TRUE))
})
