# ==============================================================================
# tests/periodic_shift_test.R — does periodic shift registration find the shift?
#
# WHY THIS EXISTS. P9.1 established that, GIVEN a correct shift, the warping
# amplitude table computes the correct amplitude. It did not establish that the
# estimator PRODUCES a correct shift, and it did not: `which.max()` is 1-based
# while element 1 of the inverse FFT is zero lag, so every periodic estimate was
# displaced by exactly one grid step. Two identical curves were registered with
# s = -0.01 -- 14.4 minutes on a 24-hour day -- and because that shift is
# APPLIED, the per-subject R-squared, RMSE and correlation in the table were all
# computed on a curve that had been moved when it should not have been.
#
# The lesson is the shape of the gap, not the arithmetic: a test that feeds an
# estimator its own correct answer tests the consumer, not the estimator. This
# file drives linear_shift_alignment(periodic = TRUE) end to end on curves whose
# true displacement is known.
#
# Run with:   Rscript tests/periodic_shift_test.R      (from the FCK directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(fda); library(shiny) })

app_dir <- if (dir.exists("server")) "." else "FCK"
failures <- 0L
fail <- function(...) { cat("FAIL:", ..., "\n"); failures <<- failures + 1L }
ok   <- function(...) cat("ok  :", ..., "\n")
near <- function(label, a, b, tol) {
  d <- max(abs(a - b))
  if (!is.finite(d) || d > tol)
    fail(sprintf("%s: max |diff| = %.6g (tol %.6g)", label, d, tol))
  else ok(sprintf("%-56s max |diff| = %.2g", label, d))
}

source(file.path(app_dir, "server/04_helpers_fd.R"))

# --- the estimator, as shipped ------------------------------------------------
# This used to parse server/40_fpca.R and evaluate just the one assignment,
# because linear_shift_alignment() was defined inside the server function and
# could not be sourced. P11.1 moved it to a pure kernel file, so the test now
# loads exactly what the app loads and what the exported script emits -- one
# definition, no extraction, and nothing that could go stale against a copy.
env <- new.env(parent = globalenv())
source(file.path(app_dir, "server/05_helpers_warp.R"), local = env)
if (!exists("linear_shift_alignment", envir = env, inherits = FALSE)) {
  cat("FAIL: server/05_helpers_warp.R does not define linear_shift_alignment\n")
  quit(status = 1)
}
align <- get("linear_shift_alignment", envir = env)

# --- curves with a KNOWN periodic displacement --------------------------------
n_grid <- 100
tp <- seq(0, 1, length.out = n_grid)
profile <- function(t) exp(cos(2 * pi * t))     # smooth, strictly periodic
deltas <- c(0, 0, 0.05, 0.10, -0.05, -0.10, 0.25)
Y <- vapply(deltas, function(d) profile(tp - d), numeric(n_grid))   # time x curve

# an fd object over the same grid, interpolating so the curves are unchanged
basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(n_grid, 40))
fd    <- smooth.basis(tp, Y, fdPar(basis, 2, 1e-12))$fd

res <- align(fd, periodic = TRUE, allow_dilation = FALSE,
             reference = "first", time_points = tp)
if (is.null(res) || is.null(res$shifts)) {
  cat("FAIL: periodic alignment returned nothing\n"); quit(status = 1)
}

cat("-- shift recovery ------------------------------------------------------\n")
cat(sprintf("%-12s %-14s %-14s\n", "true delta", "estimated s", "|s| reported"))
amp <- fck_warp_amplitude(res)
for (i in seq_along(deltas))
  cat(sprintf("%-12.2f %-14.4f %-14.4f\n", deltas[i], res$shifts[i], amp[i]))
cat("\n")

# 1. IDENTICAL curves must get exactly zero. This alone catches the defect.
if (abs(res$shifts[1]) > 1e-12 || abs(res$shifts[2]) > 1e-12)
  fail(sprintf("two identical periodic curves were shifted by %.4f and %.4f",
               res$shifts[1], res$shifts[2])) else ok("identical periodic curves are not shifted at all")

# 2. a known displacement is recovered to the grid resolution
tol <- 1.5 / (n_grid - 1)
near("known displacements recovered to the grid resolution",
     -res$shifts, deltas, tol)

# 3. the table and the estimator agree
near("fck_warp_amplitude() agrees with the estimated shifts",
     amp, abs(((res$shifts + 0.5) %% 1) - 0.5), 1e-12)

# 4. a curve that needed no shift comes back untouched -- the consequence that
#    made this more than a display error
cat("\n-- the applied registration --------------------------------------------\n")
# Compare against what the registration actually operates on -- the EVALUATED
# fd object -- not against the analytic profile. The basis representation is a
# smoothing round trip of about 3e-6 here, which is a property of the test's own
# fd construction, not of the registration.
orig <- eval.fd(tp, fd)[, 1]; reg <- res$registered_curves[, 1]
near("an unshifted curve is returned unchanged (RMSE)", sqrt(mean((reg - orig)^2)), 0, 1e-9)
if (stats::sd(reg) > 0 && stats::sd(orig) > 0)
  near("... and its correlation with itself is 1", cor(reg, orig), 1, 1e-12)

# 5. the warp stays inside the observed domain, for every shift
h <- res$warp_functions
if (any(h < -1e-12 | h > 1 + 1e-12))
  fail("a periodic warp left the observed domain") else ok("every periodic warp stays inside the observed domain")

# 6. and the old estimator really was wrong -- pinned, so the fix cannot be
#    quietly reverted to something that merely looks similar
cat("\n-- the defect, for the record ------------------------------------------\n")
old_est <- function(y, ref) {
  cc <- Re(fft(Conj(fft(ref - mean(ref))) * fft(y - mean(y)), inverse = TRUE)) / n_grid
  mi <- which.max(cc)
  -(if (mi > n_grid / 2) mi - n_grid else mi) / n_grid
}
s_old <- old_est(Y[, 1], Y[, 1])
if (abs(s_old) < 1e-9) fail("the old estimator was not wrong -- re-check the reproduction") else ok(sprintf("old estimator gave s = %.4f on two identical curves (%.1f min on a 24 h day)",
                s_old, abs(s_old) * 24 * 60))

cat("\n")
if (failures) { cat(sprintf("Periodic shift tests FAILED (%d).\n", failures)); quit(status = 1) }
cat("Periodic shift tests passed.\n")
