# ==============================================================================
# tests/testthat/test-p4-corrections.R
#
# P4: a third review's remaining findings. Nine held, one was overstated, and
# one of the app's own audit comments turned out to be false.
#
#   P4.1  parametric warp identities (see tests/warp_family_test.R for the maths)
#   P4.2  the shift warp was described as endpoint-anchored; it is a translation
#   P4.3  reformulate() does NOT quote uploaded column names -- the P1.4b note
#         claimed it did, and the GAM branch still pasted names into text
#   P4.4  chol2inv(qr.R(qr(X))) assumed an unpivoted QR (overstated, guarded)
#   P4.5  a design with n == p passed the rank check and then divided by zero df
#   P4.6  FoSR R2 had no guard for a constant response at a time point
#   P4.7  between-subjects F had no guard for zero within-group variation
#   P4.8  the n-basis diagnostic used B-splines even under cyclic smoothing
#   P4.9  app.R pointed at an renv.lock that does not exist
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

# ------------------------------------------- P4.3 reformulate is not a shield
test_that("P4.3: reformulate() does not protect arbitrary column names", {
  # This is the fact the P1.4b comment got wrong. Pinning it here so the claim
  # cannot be reinstated: reformulate PARSES its termlabels.
  expect_error(stats::reformulate("a b"))
  # worse than an error -- a name with parentheses becomes a CALL, silently
  f <- stats::reformulate("Age (years)")
  expect_true(grepl("Age(years)", deparse(f), fixed = TRUE))
  # and an executable label survives into the formula
  f2 <- stats::reformulate('I(1 + 1)')
  expect_true(is.call(f2[[2]]))
})

test_that("P4.3: the FoSR estimator never puts user text in a formula", {
  src <- code_of("server/07_helpers_fosr.R")
  expect_false(grepl("stats::reformulate(predictors)", src, fixed = TRUE))
  expect_true(grepl('model_names <- paste0("x", seq_along(predictors))', src, fixed = TRUE))
  expect_true(grepl("stats::reformulate(model_names)", src, fixed = TRUE))
  # the GAM branch too (moved to its own kernel file at P5.8)
  gam <- code_of("server/07b_helpers_fosr_gam.R")
  expect_true(grepl('mpred <- paste0("x", seq_along(preds))', gam, fixed = TRUE))
  expect_false(grepl('" + s(time, by = ", p_var, ", bs=', gam, fixed = TRUE))
  expect_false(grepl("input$", gam, fixed = TRUE))
})

test_that("P4.3: a hostile column name is fitted as a column, not executed", {
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(11)
  n <- 24; nt <- 8
  nasty <- c("Age (years)", "a b", 'I(stop("executed"))')
  df <- data.frame(rnorm(n), rnorm(n), rnorm(n))
  names(df) <- nasty
  Y <- matrix(rnorm(n * nt), n, nt)

  fit <- fck_fit_fosr_ols(Y, df, nasty, use_bootstrap = FALSE)
  expect_equal(nrow(fit$beta.hat), 4L)                 # intercept + 3
  # the user's own labels come back on the rows, unmangled
  expect_setequal(rownames(fit$beta.hat), c("(Intercept)", nasty))
  # and the fit is the right one: compare against lm() on renamed columns
  d2 <- df; names(d2) <- paste0("x", 1:3)
  ref <- lm(Y[, 1] ~ x1 + x2 + x3, data = d2)
  expect_equal(unname(fit$beta.hat[, 1]), unname(coef(ref)), tolerance = 1e-10)
})

# ------------------------------------------------------- P4.4 / P4.5 QR + df
test_that("P4.4: the QR pivot is checked rather than assumed", {
  src <- code_of("server/07_helpers_fosr.R")
  expect_true(grepl("if (!identical(qrX$pivot, seq_len(ncol(X))))", src, fixed = TRUE))
  expect_true(grepl("ord <- order(qrX$pivot)", src, fixed = TRUE))
})

test_that("P4.4: the unpivoting branch is arithmetically right when it fires", {
  # Force a permuted decomposition and check the undo reproduces (X'X)^-1.
  set.seed(3); X <- cbind(1, matrix(rnorm(40 * 3), 40, 3))
  perm <- c(3, 1, 4, 2)
  Xp <- X[, perm, drop = FALSE]
  q  <- qr(Xp)
  iv <- chol2inv(qr.R(q))
  ord <- order(q$pivot)
  expect_equal(iv[ord, ord, drop = FALSE][order(perm), order(perm)],
               solve(crossprod(X)), tolerance = 1e-9,
               ignore_attr = TRUE)
})

test_that("P4.5: n == p is refused instead of dividing by zero df", {
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(4)
  df <- data.frame(a = rnorm(3), b = rnorm(3))     # n = 3, p = 3 with intercept
  Y  <- matrix(rnorm(3 * 5), 3, 5)
  expect_error(fck_fit_fosr_ols(Y, df, c("a", "b")), "degrees of freedom")
  # n = p + 1 is the smallest design that CAN be fitted, and it works
  df4 <- data.frame(a = rnorm(4), b = rnorm(4))
  Y4  <- matrix(rnorm(4 * 5), 4, 5)
  f   <- fck_fit_fosr_ols(Y4, df4, c("a", "b"))
  expect_equal(f$df_resid, 1L)
  expect_true(all(is.finite(f$beta.se)))
})

# ---------------------------------------------------- P4.6 / P4.7 degeneracy
test_that("P4.6: a constant response gives R2 = NA, not NaN or -Inf", {
  source(file.path(app_dir, "server/07_helpers_fosr.R"), local = TRUE)
  set.seed(6)
  n <- 20; nt <- 6
  df <- data.frame(a = rnorm(n))
  Y  <- matrix(rnorm(n * nt), n, nt)
  Y[, 3] <- 7                                    # constant across subjects
  f <- fck_fit_fosr_ols(Y, df, "a")
  expect_true(is.na(f$r2_t[3]))
  expect_true(all(is.finite(f$r2_t[-3])))
  expect_false(any(is.nan(f$r2_t)))
})

test_that("P4.7: the between-subjects F is NA where there is no within-group variation", {
  src <- code_of("server/50_fanova.R")
  expect_true(grepl("F_stat <- ifelse(SSW > ss_floor,", src, fixed = TRUE))
  expect_true(grepl("if (!is.finite(F_stat[t]) || .np < 1) NA_real_ else", src, fixed = TRUE))
  expect_true(grepl("sig_regions <- !is.na(p_values_adjusted) & p_values_adjusted < alpha",
                    src, fixed = TRUE))
  # the arithmetic the guard replaces
  n_c <- 6
  SSB <- c(0, 4); SSW <- c(0, 0)
  bare <- (SSB / 1) / (SSW / (n_c - 2))
  expect_true(is.nan(bare[1]))          # 0/0
  expect_true(is.infinite(bare[2]))     # 4/0
  guarded <- ifelse(SSW > 1e-12, bare, NA_real_)
  expect_true(all(is.na(guarded)))
})

# ------------------------------------------------------ P4.8 cyclic diagnostic
test_that("P4.8: the n-basis sweep uses the basis the app would actually fit", {
  src <- code_of("server/30_diagnostics.R")
  expect_true(grepl("cyclic <- isTRUE(input$is_cyclic)", src, fixed = TRUE))
  expect_true(grepl("if (cyclic) create.fourier.basis(rangeval = c(1, n_time), nbasis = nb)",
                    src, fixed = TRUE))
  expect_true(grepl("basis <- make_basis(nb)", src, fixed = TRUE))
  # a Fourier basis needs an ODD count (constant + sin/cos pairs); an even grid
  # would have scored the same model twice
  expect_true(grepl("nb_seq <- if (cyclic) seq(3, min(n_time - 1, 41), by = 2)",
                    src, fixed = TRUE))
})

test_that("P4.8: fda rounds an even Fourier nbasis up, so an even grid duplicates", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  b4 <- create.fourier.basis(c(1, 24), nbasis = 4)
  b5 <- create.fourier.basis(c(1, 24), nbasis = 5)
  expect_equal(b4$nbasis, b5$nbasis)
})

# --------------------------------------------------------------- P4.9 renv ---
test_that("P4.9: nothing points at an renv.lock that is not there", {
  expect_false(file.exists(file.path(app_dir, "renv.lock")))
  app <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("see renv.lock in the project root", app, fixed = TRUE))
  expect_true(grepl("this project is NOT environment-pinned", app, fixed = TRUE))
  readme <- paste(readLines(file.path(app_dir, "README.md"), warn = FALSE), collapse = "\n")
  expect_true(grepl("does **not** ship an `renv.lock`", readme, fixed = TRUE))
})
