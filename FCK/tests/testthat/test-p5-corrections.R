# ==============================================================================
# tests/testthat/test-p5-corrections.R
#
# P5: a fourth review. Nine findings, all of which held.
#
#   P5.1  FoSR OLS PREDICTION was broken by my own P4.3 fix -- a regression
#   P5.2  warped-PCA AIC/BIC were not a likelihood; removed
#   P5.3  the "EFDA variance decomposition" was not a decomposition; relabelled
#   P5.4  the Fisher-Rao phase distance was applied to translations, where it
#         is identically 0, and to wrapped periodic shifts, where it is a
#         constant artefact -- 0.142 even at zero shift
#   P5.5  manual landmark registration had no monotonicity requirement
#   P5.6  the two landmark branches stored INVERSE warps of each other
#   P5.7  RM-fANOVA turned an undefined F into F = 0, p = 1
#   P5.8  the GAM export was a reconstruction that had drifted
#   P5.9  GAM SE/p were zeros labelled "placeholders"
#   P5.10 stale "REML" labels on a GCV smoother
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

# AUDIT (P11.1): the registration kernels moved out of server/40_fpca.R into
# server/05_helpers_warp.R. A guard that keeps reading only the old file does
# not fail when that happens -- it goes VACUOUS, because every expect_false()
# passes once the code is simply not in that file any more, and every
# expect_true() fails for a reason that has nothing to do with the property.
# Registration guards therefore read the registration source wherever it lives.
warp_code <- function() paste(code_of("server/40_fpca.R"),
                              code_of("server/05_helpers_warp.R"), sep = "\n")

# =========================================================== P5.1 prediction ==
test_that("P5.1: a fitted FoSR model can build a prediction design for new data", {
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(21)
  n <- 24; nt <- 8
  df <- data.frame(rnorm(n, 40, 8), factor(rep(c("F", "M"), n / 2)),
                   check.names = FALSE)
  names(df) <- c("Age (years)", "Sex")          # deliberately non-syntactic
  Y <- matrix(rnorm(n * nt), n, nt) + df[[1]] %o% seq(0, 0.05, length.out = nt)

  fit <- fck_fit_fosr_ols(Y, df, c("Age (years)", "Sex"))
  # the mapping travels with the fit
  expect_identical(fit$model_names, c("x1", "x2"))
  expect_identical(fit$predictor_names, c("Age (years)", "Sex"))
  expect_true(!is.null(fit$xlevels$x2))

  nd <- data.frame(45, factor("M", levels = c("F", "M")), check.names = FALSE)
  names(nd) <- c("Age (years)", "Sex")
  X <- fck_fosr_design(fit, nd)

  expect_equal(ncol(X), nrow(fit$beta.hat))
  expect_identical(colnames(X), rownames(fit$beta.hat))
  y_hat <- as.vector(X %*% fit$beta.hat)
  expect_length(y_hat, nt)
  expect_true(all(is.finite(y_hat)))

  # ... and it is the RIGHT number: same fit through lm on renamed columns
  d2 <- df; names(d2) <- c("x1", "x2")
  ref <- lm(Y[, 1] ~ x1 + x2, data = d2)
  expect_equal(y_hat[1],
               unname(predict(ref, newdata = data.frame(
                 x1 = 45, x2 = factor("M", levels = c("F", "M"))))),
               tolerance = 1e-10)

  # the failure this replaces: model.matrix on the raw terms cannot find x1
  expect_error(model.matrix(delete.response(fit$terms), data = nd))
})

test_that("P5.1: prediction refuses new data it cannot honour, instead of blanking", {
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(22)
  n <- 20; nt <- 5
  df <- data.frame(a = rnorm(n), g = factor(rep(c("x", "y"), n / 2)))
  Y  <- matrix(rnorm(n * nt), n, nt)
  fit <- fck_fit_fosr_ols(Y, df, c("a", "g"))

  expect_error(fck_fosr_design(fit, data.frame(a = 1)), "missing predictor")
  expect_error(fck_fosr_design(fit, data.frame(a = 1, g = factor("z"))),
               "not seen when the model was fitted")
  # a single-level new frame still yields the full design (xlev is applied)
  ok <- fck_fosr_design(fit, data.frame(a = 1, g = factor("y", levels = c("x", "y"))))
  expect_equal(ncol(ok), 3L)
})

test_that("P5.1: both prediction sites go through the shared builder", {
  src <- code_of("server/70_fosr.R")
  expect_false(grepl("model.matrix(f_clean, data = pred_df)", src, fixed = TRUE))
  expect_false(grepl("model.matrix(f_clean, data = df_min)", src, fixed = TRUE))
  expect_true(grepl("fck_fosr_design(mod, pred_df)", src, fixed = TRUE))
  expect_true(grepl("fck_fosr_design(mod, df_min)", src, fixed = TRUE))
  expect_true(grepl("fck_fosr_design(mod, df_max)", src, fixed = TRUE))
})

# ============================================================ P5.2 AIC/BIC ====
test_that("P5.2: the warping panel no longer reports AIC, BIC or a log-likelihood", {
  src <- warp_code()
  expect_false(grepl("aic <- -2 * log_lik + 2 * k_params", src, fixed = TRUE))
  expect_false(grepl("bic <- -2 * log_lik + k_params * log(n_obs)", src, fixed = TRUE))
  expect_false(grepl("k_params <- 2 * n_subjects", src, fixed = TRUE))
  expect_false(grepl("summary_stats$AIC", src, fixed = TRUE))
  expect_false(grepl("stats$AIC", src, fixed = TRUE))
  # and the instruction to select a method by them is gone
  expect_false(grepl("Lower AIC/BIC values indicate better model fit", src, fixed = TRUE))
})

# ================================================= P5.3 not a decomposition ===
test_that("P5.3: the dispersion statistics are not presented as a decomposition", {
  src <- warp_code()
  expect_false(grepl("variance_explained_by_warping", src, fixed = TRUE))
  expect_false(grepl("total_amp_variance", src, fixed = TRUE))
  expect_false(grepl("total_phase_variance", src, fixed = TRUE))
  expect_true(grepl("dispersion_reduction", src, fixed = TRUE))
  expect_true(grepl("total_dispersion_pre", src, fixed = TRUE))
  expect_true(grepl("total_dispersion_post", src, fixed = TRUE))
  # the arithmetic point: the three old quantities did not sum
  set.seed(31)
  x <- matrix(rnorm(50 * 6), 50, 6); dt <- 1 / 49
  m <- rowMeans(x)
  v_pre <- sum(apply(x, 2, function(c) sum((c - m)^2) * dt))
  parts <- v_pre * c(0.4, 0.4)   # any independently computed pair
  expect_false(isTRUE(all.equal(sum(parts), v_pre)))
})

# ============================================ P5.4 phase metric by geometry ===
test_that("P5.4: the Fisher-Rao phase distance is blind to translation", {
  t <- seq(0, 1, length.out = 100); dt <- diff(t[1:2])
  fr <- function(h) {
    d <- c(diff(h) / diff(t), 0); d[d < 0] <- 0
    acos(min(1, max(-1, sum(sqrt(pmax(0, d))) * dt)))
  }
  # a translation has h' == 1, so the metric is 0 for EVERY shift
  for (s in c(0, 0.05, 0.15, 0.25))
    expect_equal(fr(t - s), 0, tolerance = 1e-12)
  # a wrapped periodic shift gives the same nonzero artefact at every shift,
  # INCLUDING zero -- which is the sharpest form of the defect
  w <- vapply(c(0, 0.05, 0.15, 0.25), function(s) fr((t - s) %% 1), numeric(1))
  expect_gt(w[1], 0.1)
  expect_equal(diff(range(w)), 0, tolerance = 1e-9)
  # on a genuine endpoint-preserving warp it behaves: 0 at the identity,
  # increasing with deformation
  expect_equal(fr(t^1), 0, tolerance = 1e-12)
  expect_gt(fr(t^2), fr(t^1.3))
})

test_that("P5.4: the module computes the phase summary per geometry", {
  src <- warp_code()
  expect_true(grepl('geom <- if (identical(method, "linear_shift")) "shift" else "interval"',
                    src, fixed = TRUE))
  expect_true(grepl('if (identical(geom, "shift")) {', src, fixed = TRUE))
  expect_true(grepl("elastic_phase_dist[i] <- NA_real_", src, fixed = TRUE))
  expect_true(grepl("endpoints_ok <- abs(warp_i[1] - time_points[1]) < 1e-6",
                    src, fixed = TRUE))
  expect_true(grepl("phase_displacement", src, fixed = TRUE))
})

# ================================== P5.5 / P5.6 landmark warps ================
test_that("P5.5: a landmark warp is built and validated in one place", {
  src <- warp_code()
  expect_true(grepl("fck_landmark_warp <- function(ref, own, time_points)", src, fixed = TRUE))
  # neither branch builds its own knots any more
  expect_false(grepl("all_curve_landmarks <- c(0, curve_landmarks, 1)", src, fixed = TRUE))
  expect_false(grepl("kx <- c(0, own, 1); ky <- c(0, ref_lm, 1)", src, fixed = TRUE))
  # one call in each of the two landmark branches (the definition reads
  # "fck_landmark_warp <- function(", so it does not match this pattern)
  expect_equal(length(gregexpr("fck_landmark_warp(", src, fixed = TRUE)[[1]]), 2L)
})

test_that("P5.5: crossed or duplicated landmarks are rejected, not folded", {
  # transcribe the builder and check its contract
  tp <- seq(0, 1, length.out = 101)
  # Pull the real function out of the file by PARSING it and taking the
  # assignment expression, rather than by regex on the text -- the point of the
  # test is that the shipped function has this contract.
  src_env <- new.env(parent = globalenv())
  # P11.1: fck_landmark_warp() moved to the pure kernel file. It is a plain
  # top-level definition there, so it no longer has to be dug out of an AST --
  # which is itself the P11.2 fix: while it was nested inside a reactive
  # observer, its own caller could not see it and every landmark registration
  # silently fell back to a linear shift.
  exprs <- parse(file.path(app_dir, "server/05_helpers_warp.R"), encoding = "UTF-8")
  got <- FALSE
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (length(e) == 3 && as.character(e[[1]]) %in% c("<-", "=") &&
        is.name(e[[2]]) && identical(as.character(e[[2]]), "fck_landmark_warp")) {
      eval(e, envir = src_env); got <<- TRUE; return(invisible(NULL))
    }
    for (i in seq_along(e)) try(walk(e[[i]]), silent = TRUE)
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  expect_true(got)
  f <- get("fck_landmark_warp", envir = src_env)

  # a well-formed pair gives a strictly increasing map fixing both endpoints
  h <- f(c(0.25, 0.5, 0.75), c(0.20, 0.55, 0.70), tp)
  expect_false(is.null(h))
  expect_equal(h[1], 0); expect_equal(h[length(h)], 1)
  expect_true(all(diff(h) > 0))

  # the identity in, the identity out
  hid <- f(c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75), tp)
  expect_equal(max(abs(hid - tp)), 0, tolerance = 1e-12)

  # a known displacement is carried exactly at the knots
  h2 <- f(c(0.5), c(0.6), tp)
  expect_equal(h2[which.min(abs(tp - 0.5))], 0.6, tolerance = 1e-9)

  # crossed landmarks -> rejected
  expect_null(f(c(0.25, 0.5, 0.75), c(0.20, 0.70, 0.55), tp))
  # duplicated -> rejected
  expect_null(f(c(0.25, 0.5, 0.75), c(0.3, 0.3, 0.8), tp))
  # a landmark on the boundary -> rejected (it would collide with the endpoint)
  expect_null(f(c(0.0, 0.5), c(0.1, 0.5), tp))
  # mismatched lengths -> rejected
  expect_null(f(c(0.25, 0.5), c(0.3), tp))
})

test_that("P5.6: every registration method returns h in one direction", {
  src <- warp_code()
  expect_equal(length(gregexpr('warp_direction = "registered -> original"',
                               src, fixed = TRUE)[[1]]), 3L)
  # the inverse application in the automatic branch is gone
  expect_false(grepl("approx(warp_functions[, i], curves[, i],", src, fixed = TRUE))
  # all three apply the warp the same way
  expect_gte(length(gregexpr("approx(time_points, curves[, i],", src, fixed = TRUE)[[1]]), 2L)
})

# ======================================================== P5.7 RM degenerate ==
test_that("P5.7: an undefined repeated-measures F stays NA", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("F_stat[is.na(F_stat)] <- 0", src, fixed = TRUE))
  expect_false(grepl("p_values_pointwise[is.na(p_values_pointwise)] <- 1", src, fixed = TRUE))
  expect_false(grepl("F_stat_perm[t, perm] <- 0", src, fixed = TRUE))
  expect_true(grepl("F_stat_perm[t, perm] <- NA_real_", src, fixed = TRUE))
  expect_true(grepl("if (!is.finite(F_stat[t]) || .np < 1) NA_real_ else", src, fixed = TRUE))
  # both branches now define significance the same way
  expect_equal(length(gregexpr("sig_regions <- !is.na(p_values_adjusted) & p_values_adjusted < alpha",
                               src, fixed = TRUE)[[1]]), 2L)
})

test_that("P5.7: entering 0 for a missing null draw biases the p-value down", {
  # the arithmetic reason the old code was not merely cosmetic
  set.seed(41)
  obs <- 2.5
  perm <- rnorm(999, 1, 1)
  perm[1:200] <- NA_real_                      # draws with no residual variation
  p_na   <- (1 + sum(perm >= obs, na.rm = TRUE)) / (1 + sum(is.finite(perm)))
  perm0  <- perm; perm0[is.na(perm0)] <- 0     # the old behaviour
  p_zero <- (1 + sum(perm0 >= obs)) / (1 + length(perm0))
  expect_lt(p_zero, p_na)
})

# ==================================================== P5.8 / P5.9 GAM branch ==
test_that("P5.8: the GAM estimator is a kernel the export emits verbatim", {
  helper <- code_of("server/07b_helpers_fosr_gam.R")
  expect_true(grepl("fck_fit_fosr_gam <- function", helper, fixed = TRUE))
  for (bad in c("input$", "showNotification(", "withProgress(", "values$", "session$"))
    expect_false(grepl(bad, helper, fixed = TRUE), info = bad)
  # the spline dimension is an argument, not three hard-coded literals
  expect_true(grepl('sprintf("Y_val ~ s(time, bs = \'ps\', k = %d)", k)', helper, fixed = TRUE))
  # the GUI calls it, and the export emits it
  expect_true(grepl("fck_fit_fosr_gam(", code_of("server/70_fosr.R"), fixed = TRUE))
  exp_src <- code_of("server/90_export.R")
  expect_true(grepl('emit_kernel("fck_fit_fosr_gam")', exp_src, fixed = TRUE))
  # the reconstruction, with uploaded names pasted into formula text, is gone
  expect_false(grepl("gam_formula <- as.formula(paste('value ~ s(time) +',", exp_src, fixed = TRUE))
  expect_false(grepl("sprintf('s(time, by = %s)', fosr_predictors)", exp_src, fixed = TRUE))
})

test_that("P5.9: unknown GAM standard errors are NA, not zero", {
  src <- code_of("server/07b_helpers_fosr_gam.R")
  expect_false(grepl("beta.se = beta_hat*0", src, fixed = TRUE))
  expect_false(grepl("beta.p  = beta_hat*0", src, fixed = TRUE))
  expect_true(grepl("na_like <- matrix(NA_real_", src, fixed = TRUE))
  expect_true(grepl('inference = "none"', src, fixed = TRUE))
  expect_true(grepl("inference_note", src, fixed = TRUE))
  # and the consumer no longer does `if (sum(all-NA) > 0)`, which errors
  gui <- code_of("server/70_fosr.R")
  expect_false(grepl("sum(mod$beta.se[idx, ]) > 0", gui, fixed = TRUE))
  expect_true(grepl("any(is.finite(mod$beta.se[idx, ]) & mod$beta.se[idx, ] > 0)",
                    gui, fixed = TRUE))
  expect_error(if (sum(c(NA_real_, NA_real_)) > 0) TRUE)   # why it had to change
})

# ============================================================ P5.10 labels ====
test_that("P5.10: nothing calls the GCV smoother REML any more", {
  for (f in c("server/21_smoothing_views.R", "ui/30_diagnostics.R",
              "server/30_diagnostics.R", "server/20_smoothing.R",
              "server/90_export.R")) {
    src <- code_of(f)
    expect_false(grepl("Automatic (REML)", src, fixed = TRUE), info = f)
    expect_false(grepl("automatic REML optimization", src, fixed = TRUE), info = f)
    expect_false(grepl("Lambda = 0 uses REML optimization", src, fixed = TRUE), info = f)
    expect_false(grepl("For automatic REML: lambda = 0", src, fixed = TRUE), info = f)
  }
  expect_true(grepl("Automatic (lambda by GCV)",
                    code_of("server/21_smoothing_views.R"), fixed = TRUE))
})

test_that("P5.10: no scalar-on-function text survives in the app or its UI", {
  for (f in c("app.R", "ui/10_import.R")) {
    src <- paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
    expect_false(grepl("scalar-on-function", src, fixed = TRUE), info = f)
    expect_false(grepl("seven analys", src, fixed = TRUE), info = f)
  }
})
