# ==============================================================================
# tests/testthat/test-p3-corrections.R
#
# P3: the eight items a second reviewer raised after the SoFR removal, plus two
# defects found while acting on them.
#
#   P3.1  the n-basis GCV diagnostic still swept at lambda = 0 in auto mode
#   P3.2  FoSR: t against a BOOTSTRAP SE, referred to t_{n-p}
#   P3.3  repeated-measures effect size was classical, not partial, eta-squared
#   P3.4  the FoSR export was a hand-written reconstruction that had drifted
#   P3.5  deparsed kernels referenced helpers (and input$) the script never had
#   P3.5b the export's auto-smoothing block still emitted lambda = 0 as "REML"
#   P3.7  README still documented SoFR and refund
#   P3.8  the cosinor export read f$amplitude / f$acrophase, which do not exist
#
# The numeric round-trip that backs P3.4/P3.5/P3.8 lives in
# tests/export_roundtrip_test.R -- it needs a subprocess and fda, so it is a
# standalone script rather than a testthat file.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"

# Quote-aware comment stripper: the fixes quote the removed code in their
# comments, so a guard that reads prose finds the string it is meant to prove
# gone. A naive sub("#.*$") would instead eat code after a "#" inside a string.
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

# ------------------------------------------------------- P3.1 n-basis sweep --
test_that("the n-basis GCV sweep no longer scores auto mode at lambda = 0", {
  src <- code_of("server/30_diagnostics.R")
  expect_false(grepl('lam      <- if(input$smooth_method == "auto") 0', src, fixed = TRUE))
  expect_false(grepl('if(input$smooth_method == "auto") 0 else', src, fixed = TRUE))
  # in auto mode it re-selects lambda by GCV at each candidate basis size,
  # using the same search the app runs
  expect_true(grepl("fck_auto_lambda(data_mat, seq_len(n_time), basis", src, fixed = TRUE))
  expect_true(grepl("lam_used[bi]  <- al$lambda", src, fixed = TRUE))
  # and the stored result records which regime produced the curve
  expect_true(grepl("auto     = auto,", src, fixed = TRUE))
})

test_that("P3.1: with lambda re-selected per basis, GCV is flatter than at lambda = 0", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)

  set.seed(4)
  n_t <- 24; hrs <- 0:23
  dat <- t(vapply(1:8, function(i)
    50 + 8 * cos(2 * pi * (hrs - 16) / 24) + rnorm(n_t, 0, 1.5), numeric(n_t)))

  nb_seq <- seq(4, 20, by = 4)
  gcv_at_zero <- vapply(nb_seq, function(nb) {
    b <- create.bspline.basis(c(1, n_t), nbasis = nb)
    mean(vapply(seq_len(nrow(dat)), function(i)
      suppressWarnings(as.numeric(smooth.basis(1:n_t, dat[i, ], fdPar(b, 2, 0))$gcv)),
      numeric(1)), na.rm = TRUE)
  }, numeric(1))
  gcv_auto <- vapply(nb_seq, function(nb) {
    b <- create.bspline.basis(c(1, n_t), nbasis = nb)
    fck_auto_lambda(dat, seq_len(n_t), b, n_grid = 12)$gcv
  }, numeric(1))

  # Selecting lambda can only improve the GCV it is selected on.
  expect_true(all(gcv_auto <= gcv_at_zero + 1e-8))
  # And the point of the fix: the unpenalised sweep is driven by n_basis, the
  # penalised one much less so.
  rng <- function(x) diff(range(x)) / mean(x)
  expect_lt(rng(gcv_auto), rng(gcv_at_zero))
})

# ------------------------------------------------ P3.3 partial eta-squared ----
test_that("the repeated-measures branch reports partial eta-squared", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("eta_squared <- SSB / (SST + 1e-10)", src, fixed = TRUE))
  expect_true(grepl("SS_visit_t[ok_ss] /", src, fixed = TRUE))
  expect_true(grepl('eta_squared_type <- "partial"', src, fixed = TRUE))
  # the between-subjects branch keeps classical eta-squared, and says so
  expect_true(grepl('eta_squared_type = "classical"', src, fixed = TRUE))
  # the sums of squares the test used are kept rather than recomputed
  expect_true(grepl("SS_visit_t[t]    <- SS_visit", src, fixed = TRUE))
  expect_true(grepl("SS_residual_t[t] <- SS_residual", src, fixed = TRUE))
})

test_that("P3.3: partial eta-squared equals df_c*F/(df_c*F + df_r), and classical does not", {
  set.seed(7)
  n_sub <- 14; k <- 3
  # A design with LARGE between-subject spread: this is the case where the two
  # effect sizes diverge most, and where the old number was most misleading.
  subj  <- rnorm(n_sub, 0, 20)
  cond  <- c(0, 2.5, 5)
  Y <- outer(subj, rep(1, k)) + matrix(cond, n_sub, k, byrow = TRUE) +
       matrix(rnorm(n_sub * k, 0, 2), n_sub, k)

  grand <- mean(Y); cm <- colMeans(Y); sm <- rowMeans(Y)
  E  <- sweep(sweep(Y, 1, sm, "-"), 2, cm, "-") + grand
  SS_cond <- sum(n_sub * (cm - grand)^2)
  SS_err  <- sum(E^2)
  df_c <- k - 1; df_r <- (n_sub - 1) * (k - 1)
  Fv <- (SS_cond / df_c) / (SS_err / df_r)

  partial   <- SS_cond / (SS_cond + SS_err)
  classical <- SS_cond / sum((Y - grand)^2)

  # the identity the app's readout claims
  expect_equal(partial, (df_c * Fv) / (df_c * Fv + df_r), tolerance = 1e-12)
  # and it agrees with stats::aov's own decomposition
  d  <- data.frame(y = as.vector(Y),
                   s = factor(rep(seq_len(n_sub), k)),
                   c = factor(rep(seq_len(k), each = n_sub)))
  av <- summary(aov(y ~ c + Error(s/c), data = d))
  tab <- av[["Error: s:c"]][[1]]
  expect_equal(unname(tab[["F value"]][1]), Fv, tolerance = 1e-8)

  # the defect: with this much subject variance the old number is far smaller
  expect_gt(partial, 5 * classical)
})

# --------------------------------------------------- P3.4 / P3.5 the export ---
test_that("the FoSR estimator is one named function, called by the GUI", {
  helper <- code_of("server/07_helpers_fosr.R")
  gui    <- code_of("server/70_fosr.R")
  expect_true(grepl("fck_fit_fosr_ols <- function", helper, fixed = TRUE))
  expect_true(grepl("fck_fit_fosr_ols(", gui, fixed = TRUE))
  # the helper must be free of Shiny, or the deparsed copy cannot run
  for (bad in c("input$", "showNotification(", "withProgress(", "incProgress(",
                "values$", "session$"))
    expect_false(grepl(bad, helper, fixed = TRUE), info = bad)
})

test_that("the repeated-measures kernel is free of Shiny so it can be deparsed", {
  src <- code_of("server/50_fanova.R")
  lines <- strsplit(src, "\n", fixed = TRUE)[[1]]
  # code_of() strips comments, so the delimiters have to be code: the kernel
  # runs from its own definition to the next function definition in the file.
  from <- grep("^  perform_rm_fanova <- function", lines)
  to   <- grep("^  perform_functional_anova <- function", lines)
  expect_length(from, 1L)
  expect_length(to, 1L)
  expect_gt(to, from)
  body <- paste(lines[from:(to - 1)], collapse = "\n")
  for (bad in c("input$", "showNotification(", "withProgress(", "incProgress(",
                "values$", "session$"))
    expect_false(grepl(bad, body, fixed = TRUE), info = bad)
  # the Shiny pieces became arguments with inert defaults
  expect_true(grepl("run_global_test = FALSE", body, fixed = TRUE))
  expect_true(grepl("progress = NULL, notify = NULL", body, fixed = TRUE))
  # and the caller supplies the real ones
  expect_true(grepl("run_global_test = isTRUE(input$rm_global_test)", src, fixed = TRUE))
})

test_that("the export emits each kernel's helper dependency closure", {
  src <- code_of("server/90_export.R")
  expect_true(grepl("codetools::findGlobals(obj, merge = TRUE)", src, fixed = TRUE))
  expect_true(grepl("fck_fn_order <- function", src, fixed = TRUE))
  # post-order: a helper is written before anything that calls it
  expect_true(grepl("c(seen, nm)", src, fixed = TRUE))
  # non-syntactic names get backticks (the app defines `%||%` that way)
  expect_true(grepl('paste0("`", f, "`")', src, fixed = TRUE))
  # and the design vectors the RM section indexes are actually defined
  expect_true(grepl("Repeated-measures design (one entry per curve", src, fixed = TRUE))
})

test_that("P3.5b: the export no longer calls an unpenalised fit 'REML'", {
  src <- code_of("server/90_export.R")
  expect_false(grepl("Smoothing method: Automatic (REML optimization)", src, fixed = TRUE))
  expect_false(grepl("lambda = 0 triggers automatic optimization", src, fixed = TRUE))
  expect_false(grepl("fdPar(basis, Lfdobj = 2, lambda = 0)", src, fixed = TRUE))
  expect_true(grepl("lambda chosen by GCV (fck_auto_lambda)", src, fixed = TRUE))
  expect_true(grepl('emit_kernel("fck_auto_lambda")', src, fixed = TRUE))
})

test_that("P3.8: the cosinor export reads the fields fit_cosinor actually returns", {
  src <- code_of("server/90_export.R")
  expect_false(grepl("amplitude = f$amplitude[1], acrophase = f$acrophase[1]",
                     src, fixed = TRUE))
  expect_true(grepl("g1(f$amplitudes)", src, fixed = TRUE))
  expect_true(grepl("g1(f$acrophases)", src, fixed = TRUE))
})

# ------------------------------------------------------------ P3.7 the docs ---
test_that("no scalar-on-function or refund reference survives outside the record", {
  # PORTING_NOTES.md is the chronological audit log: it DESCRIBES the removal,
  # so it must still mention it. The README describes the app as it is now.
  readme <- paste(readLines(file.path(app_dir, "README.md"), warn = FALSE),
                  collapse = "\n")
  for (bad in c("SoFR", "refund", "Scalar-on-Function", "scalar-on-function",
                "71_sofr.R"))
    expect_false(grepl(bad, readme, fixed = TRUE), info = bad)
  # and the app itself carries no such tab, state or file
  expect_false(file.exists(file.path(app_dir, "server/71_sofr.R")))
  expect_false(file.exists(file.path(app_dir, "ui/71_sofr.R")))
})

test_that("P3.6: the renv bootstrap refuses to write a lockfile it cannot honour", {
  src <- code_of("tools/renv_bootstrap.R")
  # it reads app.R rather than keeping a second list that can go stale
  expect_true(grepl('as.character(ex[[2]]) %in% c("required_packages", "optional_packages")',
                    src, fixed = TRUE))
  # and it stops rather than snapshotting an incomplete library
  expect_true(grepl("Not writing a lockfile", src, fixed = TRUE))
  expect_true(grepl("missing_req <- required[!vapply(required, have, logical(1))]",
                    src, fixed = TRUE))
})
