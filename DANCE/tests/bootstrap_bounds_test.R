# ==============================================================================
# tests/bootstrap_bounds_test.R — P21 phase 1: the three reproduced defects
# ==============================================================================
# A2  the "bootstrap" selected rows with %in%, which de-duplicates, so it was a
#     ~63% subsample WITHOUT replacement -- and the finite-population correction
#     that comes with that made the interval too SHORT, not too long.
# A3  the acrophase interval was a linear quantile on a circular quantity.
# A4  ticking a parameter bound sent a model linear in every parameter to nlsLM.
#
# Run with:  Rscript tests/bootstrap_bounds_test.R      (from the DANCE dir)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/08b_helpers_popcosinor.R"), local = e)
source(file.path(app_dir, "server/08c_helpers_circstat.R"), local = e)

cat("-- A2: the resample keeps repeats ------------------------------------\n")
set.seed(1)
i <- e$dance_boot_index(100)
chk(length(i) == 100,
    "one resample of n = 100 returns 100 draws, not the distinct subset",
    sprintf("the resample returned %d draws", length(i)))
chk(length(unique(i)) < 100 && max(table(i)) > 1,
    sprintf("repeats are present (%d distinct, max multiplicity %d)",
            length(unique(i)), max(table(i))),
    "the resample has no repeats, so it is still a subsample")

# the consequence, measured: the old procedure is anti-conservative
set.seed(7); n <- 60; B <- 3000
x <- rnorm(n, 10, 2)
old_way  <- replicate(B, mean(x[(1:n) %in% sample.int(n, n, TRUE)]))
new_way  <- replicate(B, mean(x[e$dance_boot_index(n)]))
se_true  <- sd(x) / sqrt(n)
chk(abs(sd(new_way) - se_true) < 0.02 * se_true + 0.01,
    sprintf("the corrected bootstrap SE %.4f matches the analytic SE %.4f",
            sd(new_way), se_true),
    sprintf("bootstrap SE %.4f does not match analytic %.4f", sd(new_way), se_true))
chk(sd(old_way) < 0.9 * sd(new_way),
    sprintf("the old %%in%% procedure is anti-conservative: SE %.4f vs %.4f (%.0f%% short)",
            sd(old_way), sd(new_way), 100 * (1 - sd(old_way) / sd(new_way))),
    "the old procedure no longer shows its shortfall -- check the fixture")

cat("\n-- A3: the phase interval is circular --------------------------------\n")
set.seed(3)
ang <- rnorm(4000, 0, 1.2) * 2 * pi / 24      # true acrophase 0 h, SD 1.2 h
ci <- e$dance_boot_circ_ci_time(ang, period = 24, harmonic = 1)
lin <- diff(as.numeric(quantile((ang * 24 / (2 * pi)) %% 24, c(.025, .975))))
chk(ci$width < 6,
    sprintf("circular interval width %.2f h on a rhythm peaking at the origin", ci$width),
    sprintf("circular interval width is %.2f h", ci$width))
chk(lin > 20,
    sprintf("a linear quantile would have reported %.2f h -- nearly the whole cycle", lin),
    sprintf("the linear quantile reported %.2f h; the fixture no longer straddles the wrap", lin))
chk(isTRUE(ci$wraps) && ci$lo > ci$hi,
    "the wrapping interval is flagged, and its endpoints run high to low",
    "a wrapping interval was not flagged")
# and it must be origin-invariant: shifting the clock must shift the interval, not resize it
ci2 <- e$dance_boot_circ_ci_time(ang + pi / 3, period = 24, harmonic = 1)
chk(abs(ci2$width - ci$width) < 1e-9,
    "the interval width does not depend on where the time origin sits",
    sprintf("width moved with the origin: %.4f vs %.4f", ci$width, ci2$width))
# on a rhythm that does NOT straddle the wrap the two agree
set.seed(4)
ang3 <- (12 + rnorm(4000, 0, 1.2)) * 2 * pi / 24
ci3 <- e$dance_boot_circ_ci_time(ang3, 24, 1)
lin3 <- diff(as.numeric(quantile((ang3 * 24 / (2 * pi)) %% 24, c(.025, .975))))
chk(abs(ci3$width - lin3) < 0.05,
    sprintf("away from the wrap the circular and linear widths agree (%.2f vs %.2f h)",
            ci3$width, lin3),
    sprintf("the two disagree away from the wrap: %.2f vs %.2f h", ci3$width, lin3))
chk(isFALSE(ci3$wraps), "and that interval is not flagged as wrapping",
    "a non-wrapping interval was flagged")
# harmonic 2 lives on a 12 h effective period
ci4 <- e$dance_boot_circ_ci_time(ang, 24, 2)
chk(abs(ci4$width - ci$width / 2) < 1e-9 && abs(ci4$effective_period - 12) < 1e-12,
    "H2 intervals are returned on the 12 h effective period",
    "the effective-period divisor is wrong for H2")

cat("\n-- A4: a fixed basis is not sent to a nonlinear optimiser ------------\n")
src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
             collapse = "\n")
chk(!grepl('if(use_bounds && trend_type %in% c("linear", "log", "none")) {\n      return(fit_cosinor_nonlinear(',
           src, fixed = TRUE),
    "bounds no longer route a fixed-basis model straight to fit_cosinor_nonlinear()",
    "the unconditional nonlinear route is back")
chk(grepl(".bounds_violated <- function(fit)", src, fixed = TRUE) &&
      grepl('fit_route <- "closed-form least squares (no bound was binding)"', src, fixed = TRUE),
    "the closed form is tried first and the route taken is recorded",
    "the closed-form-first branch is missing")
chk(grepl("bounds_active <- TRUE", src, fixed = TRUE),
    "an actually-binding bound still falls back to the constrained optimiser",
    "the constrained fallback is missing")

cat("\n")
if (bad) { cat(sprintf("Phase 1 tests FAILED  (%d passed, %d failed)\n", ok_n, bad)); quit(status = 1) }
cat(sprintf("Phase 1 tests PASSED  (%d passed, 0 failed)\n", ok_n))
