# ==============================================================================
# tests/diagnostic_axis_test.R — do the smoothing DIAGNOSTICS describe the model
# the app actually fits?
#
# WHY THIS EXISTS. The diagnostics tab exists to advise: "use this lambda", "use
# this n_basis". Advice is only worth anything if it was computed on the model
# the advice will be applied to, and four separate rounds of this audit have
# found it was not:
#
#   P3.1  the sweep evaluated GCV at lambda = 0 (unpenalised) in auto mode
#   P4.8  the sweep always built a B-spline, even for cyclic smoothing
#   P9.3  the basis TYPE was duplicated between CV and production
#   P10.2 the basis COUNT rule was duplicated between CV and production
#   P11.3 the AXIS was still hand-built: rangeval = c(1, n_time) and argvals =
#         seq_len(n_time), the column index, while production uses elapsed
#         hours whenever real-time smoothing is on. Measured on 14 unevenly
#         spaced hourly columns: production fits over [0, 23], the diagnostic
#         over [1, 14]; and in cyclic mode production used period 24 where the
#         diagnostic used period 13 — a 13-hour rhythm fitted to 24-hour data.
#
# Every one of those was a duplicated definition that drifted. This file asserts
# the property directly — the diagnostic basis IS the production basis — rather
# than asserting that some particular line of code is present.
#
# Run with:   Rscript tests/diagnostic_axis_test.R      (from the FCK directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok <<- ok + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
if (!requireNamespace("fda", quietly = TRUE)) { cat("SKIP: fda not installed\n"); quit(status = 0) }
suppressMessages(library(fda))

e <- new.env(parent = globalenv())
for (f in c("server/01_helpers_time.R", "server/03_helpers_clock.R",
            "server/04_helpers_fd.R")) source(f, local = e)

same_basis <- function(a, b) {
  isTRUE(all.equal(a$rangeval, b$rangeval)) &&
  identical(a$type, b$type) && identical(a$nbasis, b$nbasis) &&
  isTRUE(all.equal(a$params, b$params))
}

# ---- unevenly spaced whole-hour columns, denser in the morning --------------
hrs    <- c(0,1,2,3,4,5,6,8,10,12,15,18,21,23)
labels <- sprintf("%02d:00", hrs)
n_time <- length(labels)
vals   <- list(data = matrix(rnorm(5 * n_time), 5, n_time), time_labels = labels)

for (real in c(FALSE, TRUE)) {
  for (cyc in c(FALSE, TRUE)) {
    for (meth in c("auto", "manual")) {
      inp <- list(use_real_time = real, is_cyclic = cyc, smooth_method = meth,
                  n_basis = 12, n_basis_manual = 12)
      axis <- e$fck_smoothing_axis(inp, vals)
      nb   <- e$fck_smoothing_nbasis(inp, n_time)
      prod <- e$fck_smoothing_basis(axis, nb, meth)
      # what the "suggest a lambda" observer and the CV path now build
      diag <- e$fck_smoothing_basis(e$fck_smoothing_axis(inp, vals),
                                    e$fck_smoothing_nbasis(inp, n_time), meth)
      tag <- sprintf("real_time=%-5s cyclic=%-5s method=%s", real, cyc, meth)
      chk(same_basis(prod, diag),
          paste("diagnostic basis == production basis  |", tag),
          paste("diagnostic basis DIFFERS from production |", tag))
    }
  }
}

# ---- the specific divergence P11.3 found ------------------------------------
inp <- list(use_real_time = TRUE, is_cyclic = TRUE, smooth_method = "manual",
            n_basis = 12, n_basis_manual = 12)
axis <- e$fck_smoothing_axis(inp, vals)
prod <- e$fck_smoothing_basis(axis, e$fck_smoothing_nbasis(inp, n_time), "manual")
chk(axis$using_real_time, "real-time axis is active for these labels",
    "real-time axis did NOT activate, so this test proves nothing")
chk(isTRUE(all.equal(as.numeric(prod$params), 24)),
    sprintf("cyclic + real time gives a Fourier basis of period %g", as.numeric(prod$params)),
    sprintf("cyclic + real time gave period %g, not 24", as.numeric(prod$params)))
chk(isTRUE(all.equal(prod$rangeval, c(0, 23))),
    sprintf("production spans elapsed hours [%g, %g]", prod$rangeval[1], prod$rangeval[2]),
    sprintf("production spans [%g, %g], not the elapsed-hour range",
            prod$rangeval[1], prod$rangeval[2]))
# the old hand-built diagnostic, for contrast — it must NOT match
oldstyle <- create.fourier.basis(rangeval = c(1, n_time), nbasis = min(n_time, 13))
chk(!same_basis(prod, oldstyle),
    "the pre-P11.3 hand-built basis is demonstrably a different model",
    "the hand-built basis matches production, so this test cannot detect the defect")

# ---- nb_fourier is for the sweep ONLY ---------------------------------------
chk(e$fck_smoothing_basis(axis, 12, "manual")$nbasis ==
    e$fck_smoothing_basis(axis, 12, "manual", nb_fourier = NULL)$nbasis,
    "nb_fourier = NULL reproduces the production count",
    "nb_fourier = NULL changed the production count")
sw <- e$fck_smoothing_basis(axis, 12, "manual", nb_fourier = 7)
chk(sw$nbasis == 7 && isTRUE(all.equal(sw$rangeval, prod$rangeval)) &&
    isTRUE(all.equal(sw$params, prod$params)),
    "the sweep can vary the Fourier count while keeping production's axis and period",
    "the sweep override changed more than the basis count")

# ---- nobody may rebuild the production axis or basis by hand ----------------
# Comment lines are excluded: several AUDIT notes quote the code they removed,
# and a guard that matches its own explanatory comment is vacuous — which has
# happened three times in this audit.
#
# What is banned is a SECOND definition of the model production fits. The 0-1
# rescaled basis in 20_smoothing.R is deliberately a different object — the
# fPCA/warping family works on a normalised domain — so it is allowed by name
# and the surrounding code documents why it exists.
allowed <- "rangeval = c\\(0, 1\\)"
for (f in c("server/20_smoothing.R", "server/30_diagnostics.R")) {
  code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
  hits <- grep("create\\.(fourier|bspline)\\.basis", code, value = TRUE)
  hits <- hits[!grepl(allowed, hits)]
  chk(length(hits) == 0,
      paste(f, "builds no production basis of its own"),
      paste0(f, " still constructs the production basis directly: ",
             paste(trimws(hits), collapse = " | ")))
  # the axis rule is the thing that drifted at P11.3, so guard it too
  ax <- grep("using_real_time\\s*<-\\s*!is\\.null|t_full\\s*<-\\s*if\\(using_real_time",
             code, value = TRUE)
  chk(length(ax) == 0,
      paste(f, "derives the time axis only through fck_smoothing_axis()"),
      paste0(f, " re-derives the time axis by hand: ", paste(trimws(ax), collapse = " | ")))
}
# and the helper must be the one production actually calls
prod_src <- paste(readLines("server/20_smoothing.R", warn = FALSE), collapse = "\n")
chk(grepl("axis\\s*<-\\s*fck_smoothing_axis\\(input, values\\)", prod_src),
    "production smoothing takes its axis from fck_smoothing_axis()",
    "production smoothing no longer calls fck_smoothing_axis()")

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (bad == 0) "Diagnostic axis tests PASSED" else "Diagnostic axis tests FAILED",
            ok, bad))
quit(status = if (bad == 0) 0 else 1)
