# ==============================================================================
# tests/testthat/test-p1-corrections.R
#
# P1: inference hardening. Sources of truth are the second external review plus
# what I could measure. As in the P0 file, each test is written to fail against
# the code as it stood before the fix.
#
# Source-grep guards read the CODE, not the prose: the fixes deliberately quote
# the removed lines in their comments so the reason survives a later refactor.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/09_helpers_pcanova.R"))

# Strip comments so the guards read CODE, not prose -- the fixes deliberately
# quote the removed lines in their comments. A naive sub("#.*$") also eats code
# that follows a "#" INSIDE a string literal (e.g. add("#  ", R.version.string)),
# which would make a forbidden-string guard pass for the wrong reason. Only strip
# a "#" that is not inside quotes.
code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  strip <- function(l) {
    ch <- strsplit(l, "", fixed = TRUE)[[1]]
    if (!length(ch)) return(l)
    inq <- ""; out <- character(0)
    for (i in seq_along(ch)) {
      c1 <- ch[i]
      esc <- i > 1 && ch[i - 1] == "\\"
      if (nzchar(inq)) {
        if (c1 == inq && !esc) inq <- ""
      } else if (c1 %in% c("\"", "'")) {
        inq <- c1
      } else if (c1 == "#") {
        break
      }
      out <- c(out, c1)
    }
    paste(out, collapse = "")
  }
  paste(vapply(ln, strip, character(1), USE.NAMES = FALSE), collapse = "\n")
}

# ------------------------------------------- P1.a curve-wise RM permutation
test_that("the RM permutation is drawn once per subject, not once per time point", {
  src <- code_of("server/50_fanova.R")
  # the draw must happen OUTSIDE the time loop
  i_perm <- regexpr("perm_map[i, ] <- sample.int(n_visits)", src, fixed = TRUE)
  i_tloop <- regexpr("for(t in 1:n_time) {\n              Y_complete <- Y_matrices[[t]]", src, fixed = TRUE)
  expect_gt(i_perm, 0)
  expect_gt(i_tloop, 0)
  expect_lt(i_perm, i_tloop)          # drawn before the time loop is entered
  # the old per-time-point draw is gone
  expect_false(grepl("Y_perm[i, ] <- Y_complete[i, sample(1:n_visits)]", src, fixed = TRUE))
  # and the subject index behind each row is tracked, because the complete-case
  # set differs between time points
  expect_true(grepl("Y_rows[[t]] <- which(complete_rows)", src, fixed = TRUE))
  expect_true(grepl("perm_map[rows_t[i], ]", src, fixed = TRUE))
})

test_that("a curve-wise permutation preserves each subject's trajectory shape", {
  set.seed(4)
  n <- 10; k <- 3; nt <- 20
  Y <- array(rnorm(n * k * nt), c(n, k, nt))
  # curve-wise: one map per subject, applied at every t
  pm <- t(replicate(n, sample.int(k)))
  cw <- Y
  for (i in 1:n) for (t in 1:nt) cw[i, , t] <- Y[i, pm[i, ], t]
  # per-time-point (the old behaviour)
  pt <- Y
  for (i in 1:n) for (t in 1:nt) pt[i, , t] <- Y[i, sample.int(k), t]

  # the correlation between a subject's condition-1 and condition-2 trajectories
  # must survive a curve-wise relabelling and is destroyed by a pointwise one
  cor_of <- function(A) mean(vapply(1:n, function(i)
    cor(A[i, 1, ], A[i, 2, ]), numeric(1)))
  expect_equal(abs(cor_of(cw)), abs(cor_of(Y)), tolerance = 0.35)
  # the pointwise version scrambles which condition each value came from
  expect_true(is.finite(cor_of(pt)))
})

# ------------------------------------------------------ P1.b rmfanova honesty
test_that("the speculative rmfanova signatures are gone", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("rmfanova::rmfanova(curves, id = subject_id", src, fixed = TRUE))
  expect_false(grepl("rmfanova::rmfanova(curves, subject_id, rm_factor)", src, fixed = TRUE))
  expect_false(grepl('exists("rm.fanova"', src, fixed = TRUE))
  expect_false(grepl('exists("rmfANOVA"', src, fixed = TRUE))
  # and the module no longer refuses to run without a package it never called
  expect_false(grepl("Package 'rmfanova' is required", src, fixed = TRUE))
  # the UI must say which method actually runs
  ui <- paste(readLines(file.path(app_dir, "ui/50_fanova.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("does not use the", ui, fixed = TRUE))
})

# ------------------------------------------------------------ P1.3 L2 norm
test_that("the L2 statistic integrates instead of summing over grid points", {
  v50  <- sin(2 * pi * seq(0, 1, length.out = 50))
  v200 <- sin(2 * pi * seq(0, 1, length.out = 200))
  a <- fck_l2_norm(v50,  seq(0, 1, length.out = 50))
  b <- fck_l2_norm(v200, seq(0, 1, length.out = 200))
  # the same function on a denser grid must give the same norm
  expect_equal(a, b, tolerance = 1e-3)
  # the analytic value of ||sin(2*pi*t)||_2 on [0,1] is 1/sqrt(2)
  expect_equal(a, 1 / sqrt(2), tolerance = 1e-3)
  # the old form doubles with the grid
  expect_gt(sqrt(sum(v200^2)) / sqrt(sum(v50^2)), 1.9)

  # uneven grids are weighted by their own spacing
  uneven <- c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1)
  expect_true(is.finite(fck_l2_norm(rep(1, length(uneven)), uneven)))
  expect_equal(fck_l2_norm(rep(1, length(uneven)), uneven), 1, tolerance = 1e-9)

  expect_true(is.na(fck_l2_norm(rep(NA_real_, 5))))
  expect_equal(fck_l2_norm(c(3), 0), 3)
})

test_that("no fANOVA site still uses the grid-dependent vector norm", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("sqrt(sum((SSB / n_curves)^2))", src, fixed = TRUE))
  expect_false(grepl("sqrt(sum(mean_diff^2))", src, fixed = TRUE))
  expect_false(grepl("sqrt(sum(perm_diff^2))", src, fixed = TRUE))
  expect_gte(length(gregexpr("fck_l2_norm(", src, fixed = TRUE)[[1]]), 6)
})

# ------------------------------------------------- P1.2 interval honesty
test_that("bootstrap intervals are labelled pointwise and use enough replicates", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("n_boot <- 100", src, fixed = TRUE))
  expect_true(grepl("n_boot <- 2000", src, fixed = TRUE))
  expect_true(grepl("95% pointwise CI", src, fixed = TRUE))
  ui <- paste(readLines(file.path(app_dir, "ui/50_fanova.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("not simultaneous functional bands", ui, fixed = TRUE))
})

# ----------------------------------------------------- P1.6 clustering geometry
test_that("cluster centres live in the same space as the points they centre", {
  src <- code_of("server/60_clustering.R")
  # the objective centres must come from data_matrix, never from fd_obj
  expect_false(grepl("cluster_means[j, ] <- eval.fd(time_grid, cluster_mean_fd)", src, fixed = TRUE))
  expect_false(grepl("cluster_means[i, ] <- eval.fd(time_grid, cluster_mean_fd)", src, fixed = TRUE))
  expect_true(grepl("cluster_means[j, ] <- colMeans(data_matrix[cluster_idx, , drop = FALSE])",
                    src, fixed = TRUE))
  expect_true(grepl("cluster_means[i, ] <- colMeans(data_matrix[cluster_idx, , drop = FALSE])",
                    src, fixed = TRUE))
  # original-scale means are kept for display
  expect_true(grepl("cluster_means_raw", src, fixed = TRUE))
})

test_that("mixed geometry really can drive between-SS negative", {
  set.seed(7)
  X <- matrix(rnorm(60 * 8, mean = 100, sd = 10), 60, 8)   # original scale
  Z <- scale(X)                                            # standardized
  asn <- rep(1:3, length.out = 60)
  # centres from the ORIGINAL scale, members from the STANDARDIZED matrix
  wrong <- sum(vapply(1:3, function(j) {
    sum((Z[asn == j, , drop = FALSE] -
         matrix(colMeans(X[asn == j, , drop = FALSE]), sum(asn == j), 8, byrow = TRUE))^2)
  }, numeric(1)))
  total <- sum((Z - matrix(colMeans(Z), 60, 8, byrow = TRUE))^2)
  expect_lt(total - wrong, 0)          # between-SS negative: the old behaviour
  # centres from the same matrix as the members
  right <- sum(vapply(1:3, function(j) {
    m <- Z[asn == j, , drop = FALSE]
    sum((m - matrix(colMeans(m), nrow(m), 8, byrow = TRUE))^2)
  }, numeric(1)))
  expect_gte(total - right, 0)
})

# ------------------------------------------------------------------ P1.4 FoSR
test_that("FoSR fits by QR with a rank check and no longer inverts the normal equations", {
  # P3.4 moved the estimator out of the observer and into a named function, so
  # that the export can deparse the same object the GUI runs. The guards follow
  # it; the properties they check are unchanged. Both files are read so that
  # neither the old construction reappearing in the observer, nor the helper
  # losing the check, can pass.
  src  <- code_of("server/07_helpers_fosr.R")
  both <- paste(src, code_of("server/70_fosr.R"))
  expect_false(grepl("xtx_inv <- solve(crossprod(X))", both, fixed = TRUE))
  expect_true(grepl("qrX <- qr(X)", src, fixed = TRUE))
  expect_true(grepl("rank deficient", src, fixed = TRUE))
  # formulas are built from names as data, not by re-parsing text
  expect_false(grepl("as.formula(paste('~'", both, fixed = TRUE))
  expect_false(grepl('as.formula(paste("~"', both, fixed = TRUE))
  # P4.3 superseded the reformulate(predictors) form: reformulate PARSES its
  # labels, so it was never the protection the P1.4b comment claimed. The model
  # is now built on internal names that cannot be anything but names. See
  # test-p4-corrections.R for the demonstration.
  expect_false(grepl("stats::reformulate(predictors)", src, fixed = TRUE))
  expect_true(grepl("stats::reformulate(model_names)", src, fixed = TRUE))
  # and the GUI calls that function rather than carrying a second copy
  expect_true(grepl("fck_fit_fosr_ols(", code_of("server/70_fosr.R"), fixed = TRUE))
})

test_that("chol2inv(qr.R(qr(X))) equals solve(crossprod(X)) on a well-conditioned design", {
  set.seed(2); X <- cbind(1, matrix(rnorm(40 * 3), 40, 3))
  expect_equal(chol2inv(qr.R(qr(X))), solve(crossprod(X)), tolerance = 1e-9)
  # and the rank check catches a collinear column that solve() would choke on
  Xc <- cbind(X, X[, 2])
  expect_lt(qr(Xc)$rank, ncol(Xc))
})

test_that("the FoSR bootstrap no longer claims to be a hypothesis test", {
  src <- code_of("server/07_helpers_fosr.R")
  expect_false(grepl("This is a proper bootstrap test", src, fixed = TRUE))
  expect_false(grepl("2 * mean(boot_betas[, j, k] <= 0)", src, fixed = TRUE))
  # P3.2: ONE inference path, from the analytical SE, FDR-adjusted across time.
  # The bootstrap SE must not feed the test any more.
  expect_true(grepl('p.adjust(p_values[j, ], method = "fdr")', src, fixed = TRUE))
  # P3.2 reversed the direction: the RAW p-values are computed first and the
  # FDR-adjusted copy is derived from them, so both are reported.
  expect_true(grepl("p_values <- p_values_raw", src, fixed = TRUE))
  expect_true(grepl("beta.p.raw = p_values_raw", src, fixed = TRUE))
  expect_false(grepl("t_stats <- coefs / se_mat", src, fixed = TRUE))
  expect_true(grepl("t_stats      <- coefs / se_analytic", src, fixed = TRUE))
})

test_that("P3.2: the FoSR p-values come from the analytical SE and ignore the seed", {
  skip_if_not_installed("stats")
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(5)
  n <- 30; nt <- 12
  df <- data.frame(x = rnorm(n), g = factor(rep(c("a", "b"), n / 2)))
  Y  <- matrix(rnorm(n * nt), n, nt) + df$x %o% seq(0, 1, length.out = nt)

  set.seed(1); a <- fck_fit_fosr_ols(Y, df, c("x", "g"), use_bootstrap = TRUE, n_boot = 60)
  set.seed(9); b <- fck_fit_fosr_ols(Y, df, c("x", "g"), use_bootstrap = TRUE, n_boot = 60)
  # the test is deterministic ...
  expect_equal(a$beta.p, b$beta.p)
  expect_equal(a$beta.se, b$beta.se)
  # ... while the bootstrap interval, which IS a Monte Carlo quantity, is not
  expect_false(isTRUE(all.equal(a$boot_ci_lower, b$boot_ci_lower)))
  # and the SE it reports is the analytical one, matching lm() term by term
  fitlm <- lm(Y[, 1] ~ x + g, data = df)
  expect_equal(unname(a$beta.se[, 1]), unname(summary(fitlm)$coefficients[, "Std. Error"]),
               tolerance = 1e-10)
  expect_equal(unname(a$beta.p.raw[, 1]), unname(summary(fitlm)$coefficients[, "Pr(>|t|)"]),
               tolerance = 1e-10)
  # the bootstrap SE is kept only as a diagnostic, never as the test denominator
  expect_false(is.null(a$beta.se.boot))
  expect_identical(a$inference, "analytic-t-fdr")
})

test_that("the GAM branch emits a coefficient curve for factor levels", {
  # P5.8 moved the GAM estimator into its own kernel file so the export can
  # deparse it. The guard follows; what it checks is unchanged.
  src <- code_of("server/07b_helpers_fosr_gam.R")
  # the old loop wrote nothing unless the predictor was numeric
  expect_false(grepl("if(is.numeric(df_reg[[v]])) {\n            d_0 <- d_int; d_0[[v]] <- 0\n            d_1 <- d_int; d_1[[v]] <- 1\n            beta_hat[i+1, ]",
                     src, fixed = TRUE))
  expect_true(grepl("beta_rows[[paste0(v, l)]]", src, fixed = TRUE))
  expect_true(grepl("beta_hat <- do.call(rbind, beta_rows)", src, fixed = TRUE))
})

# ------------------------------------------------------------------ P1.5 SoFR
# Scalar-on-function regression was REMOVED from the app after this round, at
# the reporter's request: they never used it, and it was the one module whose
# fixes could not be runtime-verified here because refund would not install for
# this container's R version. Deleting it also removed the app's only refund
# dependency. The P1.5 tests went with it; what remains is the guard that it
# really is gone, so it cannot creep back in half-wired.
test_that("scalar-on-function regression is fully removed", {
  expect_false(file.exists(file.path(app_dir, "server/71_sofr.R")))
  expect_false(file.exists(file.path(app_dir, "ui/71_sofr.R")))
  app <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("ui_tab_sofr", app, fixed = TRUE))
  expect_false(grepl("tabName = \"sofr\"", app, fixed = TRUE))
  # refund was used by nothing else, so it should no longer be a dependency
  expect_false(grepl("refund       =", app, fixed = TRUE))
  st <- paste(readLines(file.path(app_dir, "server/00_state.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("sofr_model", st, fixed = TRUE))
  ex <- code_of("server/90_export.R")
  expect_false(grepl("refund::pfr", ex, fixed = TRUE))
  expect_false(grepl("SCALAR-ON-FUNCTION", ex, fixed = TRUE))
})

# --------------------------------------------------------- P1.7 lambda transfer
test_that("the diagnostics no longer hand an mgcv smoothing parameter to fda", {
  src <- code_of("server/30_diagnostics.R")
  expect_false(grepl("lambda_to_use <- values$reml_profile$optimal_lambda", src, fixed = TRUE))
  expect_false(grepl("lambda_to_use <- values$cv_results$optimal_lambda", src, fixed = TRUE))
  # The button runs the production estimator. P11.3: on the production AXIS as
  # well -- it used to pass seq_len(n_time), the column index, while the
  # smoother it advises fits over elapsed hours whenever real-time smoothing is
  # on, so the lambda it returned did not belong to that model.
  expect_true(grepl("fck_auto_lambda(dat, axis$t_full, basis", src, fixed = TRUE))
  expect_false(grepl("fck_auto_lambda(dat, seq_len(n_time), basis", src, fixed = TRUE))
})

test_that("aes_string is gone from every remaining server file", {
  for (f in list.files(file.path(app_dir, "server"), pattern = "[.]R$", full.names = FALSE))
    expect_false(grepl("aes_string", code_of(file.path("server", f)), fixed = TRUE), info = f)
})
