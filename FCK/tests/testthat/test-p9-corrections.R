# ==============================================================================
# tests/testthat/test-p9-corrections.R
#
# A sixth review. Four findings, all of which held.
#
#   P9.1  the per-subject "Warping Amplitude Scores" table still carried the
#         periodic-shift artefact -- the THIRD copy of a definition P8.2
#         corrected in the other two
#   P9.2  ui/30_diagnostics.R still told the user REML and CV "should agree",
#         which is what P8.3 removed from the server; and internal audit tags
#         had leaked into user-facing text
#   P9.3  the CV built its own basis on the integer column index, so its lambda
#         was not on the production scale under cyclic or real-time smoothing
#   P9.4  an unavailable standard error was reported as 0
#
# The three regression tests the reviewer asked for are the first, the second
# and the last here.
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

# ============================================ P9.1 one warp-amplitude rule ====
test_that("P9.1: a periodic ZERO shift has warping amplitude exactly 0", {
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)
  tp <- seq(0, 1, length.out = 100)
  s  <- c(0, 0.05, 0.10, 0.25, -0.15, 0.9)
  mk <- function(periodic) list(
    method = "linear_shift", shifts = s, time_points = tp,
    warp_functions = vapply(s, function(z) if (periodic) (tp - z) %% 1 else tp - z,
                            numeric(length(tp))),
    boundary = if (periodic) "periodic wrap" else "constant extrapolation")

  a_per <- fck_warp_amplitude(mk(TRUE))
  expect_equal(a_per[1], 0)                       # the defect: this was 0.1000
  expect_equal(a_per[2:5], c(0.05, 0.10, 0.25, 0.15))
  expect_equal(a_per[6], 0.1)                     # 0.9 forward is 0.1 backward
  expect_equal(fck_warp_amplitude(mk(FALSE))[1:5], c(0, 0.05, 0.10, 0.25, 0.15))

  # the arithmetic the old table used, for the record
  rms0 <- sqrt(mean(((tp - 0) %% 1 - tp)^2))
  expect_equal(rms0, 0.1, tolerance = 1e-9)
  expect_gt(rms0, a_per[1])

  # an endpoint-preserving warp is still measured by distance from identity
  h <- vapply(c(1, 1.3, 2), function(a) tp^a, numeric(length(tp)))
  a_par <- fck_warp_amplitude(list(method = "parametric", warp_functions = h,
                                   time_points = tp))
  expect_equal(a_par[1], 0, tolerance = 1e-12)
  expect_gt(a_par[3], a_par[2])
})

test_that("P9.1: the definition exists once and all three sites call it", {
  fd <- code_of("server/04_helpers_fd.R")
  expect_true(grepl("fck_warp_amplitude <- function(warping_results)", fd, fixed = TRUE))
  fp <- warp_code()
  expect_true(grepl("warp_amplitude <- fck_warp_amplitude(warp_results)", fp, fixed = TRUE))
  # the three old bodies are gone
  expect_false(grepl("warp_amplitude[i] <- sqrt(mean((warp_results$warp_functions[,i]",
                     fp, fixed = TRUE))
  expect_false(grepl("warp_amplitude <- abs(warp_results$alpha_values - 1)", fp, fixed = TRUE))
  pc <- code_of("server/09_helpers_pcanova.R")
  expect_true(grepl("fck_warp_amplitude(warping_results)", pc, fixed = TRUE))
  expect_false(grepl("amp <- apply(wf, 2, function(h) sqrt(mean((h - tp)^2", pc, fixed = TRUE))
})

# ============================================ P9.2 the UI stops contradicting =
test_that("P9.2: the diagnostics UI does not tell the user the lambdas should agree", {
  ui <- code_of("ui/30_diagnostics.R")
  for (bad in c("Should agree on optimal range", "If REML and CV agree",
                "If REML and CV disagree", "Compare REML vs CV"))
    expect_false(grepl(bad, ui, fixed = TRUE), info = bad)
  expect_true(grepl("Do NOT compare the two lambdas numerically", ui, fixed = TRUE))
  expect_true(grepl("not transferable", ui, fixed = TRUE))
  # and the server it sits over says the same thing
  expect_true(grepl("No ratio is reported", code_of("server/30_diagnostics.R"), fixed = TRUE))
})

test_that("P9.2: no internal audit tag appears in user-facing text", {
  # comments may carry them; strings the user reads may not
  pat <- "\\(P[0-9]+\\.[0-9]+[a-z]?\\)"
  for (f in c(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"),
              list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))) {
    src <- code_of(sub(paste0("^", app_dir, "/"), "", f))
    lit <- regmatches(src, gregexpr('"(\\\\.|[^"\\\\])*"', src))[[1]]
    bad <- lit[grepl(pat, lit)]
    expect_identical(bad, character(0), info = basename(f))
  }
})

# ============================================ P9.3 the CV basis ===============
test_that("P9.3: the extracted basis builder reproduces production exactly", {
  skip_if_not_installed("fda")
  suppressPackageStartupMessages(library(fda))
  source(file.path(app_dir, "server/04_helpers_fd.R"), local = TRUE)

  n <- 24; ti <- seq_len(n); tr <- range(ti)
  ax <- function(cyc, rt = FALSE, tf = ti)
    list(n_time = n, t_full = tf, t_rng = range(tf), using_real_time = rt, cyclic = cyc)

  same <- function(a, b) {
    expect_equal(a$nbasis, b$nbasis)
    expect_identical(class(a), class(b))
    expect_equal(a$rangeval, b$rangeval)
  }
  same(fck_smoothing_basis(ax(TRUE), 12, "manual"),
       create.fourier.basis(rangeval = tr, nbasis = min(n, 13)))
  same(fck_smoothing_basis(ax(FALSE), 12, "none"),
       create.bspline.basis(rangeval = tr, breaks = ti, norder = 4))
  same(fck_smoothing_basis(ax(FALSE), 12, "manual"),
       create.bspline.basis(rangeval = tr, nbasis = 12))

  # on real elapsed hours a cyclic basis must have a 24-hour period, not the
  # length of the recording
  rt <- seq(8, 8 + 23)
  b <- fck_smoothing_basis(ax(TRUE, TRUE, rt), 12, "manual")
  expect_equal(as.numeric(b$params[1]), 24)

  # and the label names the model the lambda belongs to
  expect_true(grepl("Fourier", fck_basis_label(ax(TRUE), fck_smoothing_basis(ax(TRUE), 12)),
                    fixed = TRUE))
  expect_true(grepl("elapsed hours",
                    fck_basis_label(ax(FALSE, TRUE, rt),
                                    fck_smoothing_basis(ax(FALSE, TRUE, rt), 12)),
                    fixed = TRUE))
})

test_that("P9.3: production and the CV both go through the shared builder", {
  sm <- code_of("server/20_smoothing.R")
  expect_true(grepl("basis   <- fck_smoothing_basis(axis, nb, input$smooth_method)",
                    sm, fixed = TRUE))
  expect_false(grepl("basis   <- create.bspline.basis(rangeval = t_rng, nbasis = nb)",
                     sm, fixed = TRUE))
  dg <- code_of("server/30_diagnostics.R")
  expect_true(grepl("cv_axis     <- fck_smoothing_axis(input, values)", dg, fixed = TRUE))
  expect_true(grepl("basis <- fck_smoothing_basis(cv_axis, nb,", dg, fixed = TRUE))
  expect_false(grepl("basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)",
                     dg, fixed = TRUE))
  # and the CV records which model its lambda belongs to
  expect_true(grepl("basis_label = tryCatch(", dg, fixed = TRUE))
})

# ============================================ P9.4 unknown SE is not zero =====
test_that("P9.4: an unavailable standard error is NA, and no band is drawn", {
  gm <- code_of("server/02_helpers_gam.R")
  expect_false(grepl("se = rep(0, length(pred))", gm, fixed = TRUE))
  expect_true(grepl("se = rep(NA_real_, length(pred))", gm, fixed = TRUE))
  expect_true(grepl("se_available = FALSE", gm, fixed = TRUE))

  fo <- code_of("server/70_fosr.R")
  expect_false(grepl("se_pred <- rep(0, n_t)", fo, fixed = TRUE))
  expect_false(grepl("se_pred <- rep(0, length(y_hat))", fo, fixed = TRUE))
  expect_true(grepl("se_pred <- rep(NA_real_, n_t)", fo, fixed = TRUE))
  # the band is drawn only where an SE exists, and its absence is stated
  expect_false(grepl("if(any(se_pred > 0)) {", fo, fixed = TRUE))
  expect_true(grepl("se_ok <- is.finite(se_pred) & se_pred > 0", fo, fixed = TRUE))
  expect_true(grepl("No confidence band: the standard error could not be computed",
                    fo, fixed = TRUE))

  # why it matters: a zero-width band and a missing band look the same and mean
  # opposite things
  y <- c(1, 2, 3); se0 <- rep(0, 3); seNA <- rep(NA_real_, 3)
  expect_equal(y - 1.96 * se0, y)              # a band exactly on the estimate
  expect_true(all(is.na(y - 1.96 * seNA)))     # nothing to draw
})
