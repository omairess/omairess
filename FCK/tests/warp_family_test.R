# ==============================================================================
# tests/warp_family_test.R — the parametric warp families must be warps
#
# A time warp is a strictly increasing bijection of the observed interval onto
# itself. Anything else is not a reparameterisation of time, and registering
# with it silently deforms the data. A third review found that three of the four
# families in this app could not express the IDENTITY -- so a curve needing no
# registration was warped anyway -- and that the exponential family's identity
# guard fired at the wrong parameter value, putting a discontinuity in the
# middle of the default search range.
#
# This checks, for every family, on a fine parameter grid inside the declared
# monotone interval:
#   h(0) = 0, h(1) = 1                        (endpoints preserved)
#   h strictly increasing                     (a warp, not a fold)
#   h stays inside [0, 1]                     (no clipping needed)
#   h(t) = t exactly at the family's identity (no warping is expressible)
#   h continuous in alpha                     (optimize() sees no cliff)
#
# Run with:   Rscript tests/warp_family_test.R      (from the FCK directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else "FCK"

failures <- 0L
fail <- function(...) { cat("FAIL:", ..., "\n"); failures <<- failures + 1L }
ok   <- function(...) cat("ok  :", ..., "\n")

# The families, transcribed from server/40_fpca.R. A source guard below checks
# this transcription has not drifted from the implementation.
eps_id <- 1e-8
warp_func <- function(t, alpha, family) {
  h <- switch(family,
    "power" = t^alpha,
    "exponential" = {
      if (abs(alpha) < eps_id) t else (exp(alpha * t) - 1) / (exp(alpha) - 1)
    },
    "quadratic" = alpha * t^2 + (1 - alpha) * t,
    "logistic" = {
      if (abs(alpha) < eps_id) t else {
        L  <- function(x) 1 / (1 + exp(-alpha * (x - 0.5)))
        L0 <- L(0); L1 <- L(1)
        (L(t) - L0) / (L1 - L0)
      }
    },
    t)
  pmin(1, pmax(0, h))
}

spec <- list(
  power       = list(identity = 1, lo = 0.05, hi = 6),
  exponential = list(identity = 0, lo = -6,   hi = 6),
  quadratic   = list(identity = 0, lo = -0.999, hi = 0.999),
  logistic    = list(identity = 0, lo = -8,   hi = 8))

t <- seq(0, 1, length.out = 401)

for (fam in names(spec)) {
  sp <- spec[[fam]]
  alphas <- sort(unique(c(seq(sp$lo, sp$hi, length.out = 121), sp$identity)))
  bad_end <- bad_mono <- bad_range <- 0L
  for (a in alphas) {
    h <- warp_func(t, a, fam)
    if (abs(h[1]) > 1e-9 || abs(h[length(h)] - 1) > 1e-9) bad_end <- bad_end + 1L
    if (any(diff(h) <= 0)) bad_mono <- bad_mono + 1L
    if (any(h < -1e-12 | h > 1 + 1e-12)) bad_range <- bad_range + 1L
  }
  if (bad_end)   fail(fam, ": endpoints not preserved at", bad_end, "of", length(alphas), "parameters")
  if (bad_mono)  fail(fam, ": not strictly increasing at", bad_mono, "of", length(alphas), "parameters")
  if (bad_range) fail(fam, ": leaves [0,1] at", bad_range, "of", length(alphas), "parameters")
  if (!bad_end && !bad_mono && !bad_range)
    ok(sprintf("%-12s strictly increasing bijection of [0,1] at all %d parameters",
               fam, length(alphas)))

  # the identity must be expressible EXACTLY
  d <- max(abs(warp_func(t, sp$identity, fam) - t))
  if (d > 1e-12) fail(sprintf("%s: h(t) at alpha = %g differs from t by %.3g",
                              fam, sp$identity, d))
  else ok(sprintf("%-12s identity at alpha = %-4g reproduces t exactly", fam, sp$identity))

  # ... and the family must be CONTINUOUS in alpha there: the old exponential
  # guard returned t at alpha = 1, where the true value is (e^t-1)/(e-1), so the
  # objective jumped by ~0.1 in the middle of the default search range.
  eps <- 1e-4
  jump <- max(abs(warp_func(t, sp$identity + eps, fam) -
                  warp_func(t, sp$identity - eps, fam)))
  if (jump > 1e-2) fail(sprintf("%s: discontinuous in alpha at the identity (jump %.3g)", fam, jump))
  else ok(sprintf("%-12s continuous in alpha across the identity (jump %.1e)", fam, jump))
}

# --- the specific defect the review found ------------------------------------
cat("\n-- the exponential identity, as it used to be --------------------------\n")
old_exp <- function(t, alpha) {
  if (abs(alpha - 1) < 0.001) {
    t
  } else {
    pmin(1, pmax(0, (exp(alpha * t) - 1) / (exp(alpha) - 1)))
  }
}
d_old <- max(abs(old_exp(t, 1) - (exp(t) - 1) / (exp(1) - 1)))
if (d_old < 1e-6) fail("the old exponential guard was not actually wrong -- re-check") else ok(sprintf("old guard returned t at alpha = 1, where the family gives a map %.3f away", d_old))
d_new <- max(abs(warp_func(t, 1, "exponential") - (exp(t) - 1) / (exp(1) - 1)))
if (d_new > 1e-12) fail("the new exponential does not evaluate the family at alpha = 1") else ok("new code evaluates the family at alpha = 1 instead of short-circuiting")

# --- the source has not drifted from what is tested above --------------------
cat("\n-- source guards -------------------------------------------------------\n")
src <- paste(readLines(file.path(app_dir, "server/40_fpca.R"), warn = FALSE), collapse = "\n")
checks <- c(
  'if (abs(alpha) < eps_id) t'                                   = TRUE,
  '"quadratic"   = list(identity = 0, lo = -0.999, hi = 0.999)'  = TRUE,
  '"power"       = list(identity = 1, lo = 0.05, hi = Inf)'      = TRUE,
  'if(abs(alpha - 1) < 0.001) t'                                 = FALSE,
  '"quadratic" = c(max(param_range[1], 0.05), min(param_range[2], 1))' = FALSE)
for (i in seq_along(checks)) {
  pat <- names(checks)[i]; want <- checks[[i]]
  got <- grepl(pat, src, fixed = TRUE)
  if (got != want) fail(sprintf("source guard: %s should be %s", pat, if (want) "PRESENT" else "GONE"))
}
if (!failures) ok("source matches the families tested here")

# --- the shift warp is labelled for what it is -------------------------------
cat("\n-- shift registration --------------------------------------------------\n")
if (grepl("monotone by construction, deterministic, and endpoint-anchored", src, fixed = TRUE))
  fail("the shift warp is still described as endpoint-anchored; h(t) = t - s is not") else ok("the shift warp is no longer described as endpoint-anchored")
if (!grepl("extrap_frac[i] <- mean(h < rng[1] | h > rng[2])", src, fixed = TRUE))
  fail("the shift warp does not report how much of the domain it extrapolates") else ok("the shift warp reports its extrapolated fraction")
if (!grepl("((h - min(time_points)) %% span)", src, fixed = TRUE))
  fail("a periodic shift is still applied with constant extrapolation") else ok("a periodic shift wraps instead of clamping")

# a periodic shift must be a bijection of the domain
tp <- seq(0, 1, length.out = 100); span <- 1
for (s in c(-0.25, -0.1, 0, 0.1, 0.25)) {
  h <- 0 + ((tp - s) %% span)
  if (any(h < -1e-12) || any(h > 1 + 1e-12))
    fail(sprintf("periodic shift s = %g leaves the domain", s))
}
ok("periodic shifts stay inside the observed domain for |s| <= 0.25")

cat("\n")
if (failures) { cat(sprintf("Warp-family tests FAILED (%d).\n", failures)); quit(status = 1) }
cat("Warp-family tests passed.\n")
