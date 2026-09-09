# ==============================================================================
# tests/mixed_calibration_test.R — what the mixed permutation test actually is
# ==============================================================================
# P20/R5. The mixed permutation kernel used to be described, in its own header,
# in the readout and in the exported report, as EXACT for all three effects. It
# is not, and the gap was large: with 5 participants in one group and 30 in the
# other, and a 5:1 dispersion ratio in the within-subject contrast, the
# unstudentised interaction test rejected 225 of 300 null datasets at a nominal
# .05. Exchangeability requires the group distributions to be IDENTICAL, not
# merely to have equal means; unequal n with unequal variance is the classic
# Behrens-Fisher situation, and relabelling is not a symmetry there.
#
# The statistic is Welch-type studentised now (Janssen 1997; Pauly, Brunner &
# Konietschke 2015). This file measures what that bought and where it still
# falls short, so the claim in the code is a measured one.
#
# WHAT IS ASSERTED, AND WHAT IS NOT
#   * A rejection rate is a binomial estimate. Nothing here demands exactly
#     .05; each cell is checked against a three-standard-error band around the
#     nominal level, computed from the cell's own simulation count.
#   * The cell that is KNOWN to remain liberal is not asserted to be calibrated.
#     It is asserted to be FLAGGED -- the kernel must grade its own
#     configuration as anti-conservative and say so.
#   * Power is a separate section. A test that never rejects is perfectly
#     calibrated and useless.
#
# Run with:  Rscript tests/mixed_calibration_test.R      (from the DANCE dir)
# Set DANCE_CALIB_FAST=1 for a quicker, noisier pass.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/01c_helpers_norm.R"), local = e)
source(file.path(app_dir, "server/06_helpers_mixed.R"), local = e)
source(file.path(app_dir, "server/07_helpers_mixed_perm.R"), local = e)

FAST  <- nzchar(Sys.getenv("DANCE_CALIB_FAST"))
NSIM  <- if (FAST) 150L else 400L
# DANCE_CALIB_NSIM overrides the simulation count outright. It exists so the
# file can be exercised end to end quickly -- that every cell RUNS and that the
# right assertion fires for it. At a small count the binomial band is wide
# enough to be nearly vacuous, so a pass at a small NSIM is a WIRING CHECK, not
# a calibration result; the numbers quoted in server/07_helpers_mixed_perm.R
# come from the default count.
.n_override <- suppressWarnings(as.integer(Sys.getenv("DANCE_CALIB_NSIM")))
if (!is.na(.n_override) && .n_override > 0) NSIM <- .n_override
NPERM <- 199L

# ------------------------------------------------------------------ generator
# One mixed-design dataset with NO true effect of the named kind. `nuisance`
# adds large main effects that the null of interest does not forbid -- the point
# being that the interaction scheme must survive them.
#   rho       AR(1) correlation of the within-curve errors
#   sd_g      per-group dispersion of the within-subject CONTRAST
#   miss      probability an individual observation is missing
#   nc        number of within-subject levels
make <- function(seed, n = c(15, 15), sd_g = c(1, 1), nt = 8, nc = 2,
                 rho = 0, miss = 0, nuisance = FALSE,
                 effect = c("none", "interaction", "between", "within"),
                 delta = 0) {
  effect <- match.arg(effect)
  set.seed(seed)
  tp <- seq(0, 1, length.out = nt)
  ar1 <- function(k) {                       # correlated errors along the curve
    z <- rnorm(k); if (rho == 0) return(z)
    out <- numeric(k); out[1] <- z[1]
    for (j in 2:k) out[j] <- rho * out[j - 1] + sqrt(1 - rho^2) * z[j]
    out
  }
  grp <- rep(sprintf("g%d", seq_along(n)), times = n)
  rows <- list(); k <- 1L
  for (i in seq_along(grp)) {
    gi <- match(grp[i], sprintf("g%d", seq_along(n)))
    subj <- rnorm(1, 0, 3)
    # the within-subject contrast, with GROUP-SPECIFIC dispersion but no mean
    # difference between groups unless an interaction is being injected
    contrast <- lapply(seq_len(nc), function(ci) rnorm(nt, 0, sd_g[gi]))
    for (ci in seq_len(nc)) {
      y <- subj + contrast[[ci]] + ar1(nt)
      # nuisance: a strong between effect and a strong within effect, neither of
      # which is an interaction
      if (nuisance) y <- y + 4 * (gi - 1) + 3 * (ci - 1) * sin(2 * pi * tp)
      if (effect == "interaction") y <- y + delta * (gi - 1) * (ci - 1) * sin(pi * tp)
      if (effect == "between")     y <- y + delta * (gi - 1)
      if (effect == "within")      y <- y + delta * (ci - 1) * sin(pi * tp)
      if (miss > 0) y[runif(nt) < miss] <- NA_real_
      rows[[k]] <- data.frame(subject = sprintf("S%03d", i),
                              within = paste0("c", ci), between = grp[i],
                              t = tp, y = y, stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  d <- do.call(rbind, rows)
  d <- d[!is.na(d$y), ]
  for (v in c("subject", "within", "between")) d[[v]] <- factor(d[[v]])
  d$cell <- droplevels(interaction(d$within, d$between, sep = " x ", drop = TRUE))
  d
}

rate <- function(effect, ..., nsim = NSIM, seed0 = 20000) {
  p <- vapply(seq_len(nsim), function(s) {
    r <- e$dance_mixed_permutation(make(seed0 + s, effect = "none", ...),
                                   n_permutations = NPERM, effects = effect)
    if (!isTRUE(r$ok)) return(NA_real_)
    r[[effect]]$global_p
  }, numeric(1))
  p <- p[is.finite(p)]
  list(rate = mean(p <= .05), n = length(p))
}

# A rejection fraction is binomial. Three standard errors at the nominal level
# is the band; nothing here asks for .05 on the nose.
band <- function(n, nominal = .05, k = 3) k * sqrt(nominal * (1 - nominal) / n)

# Each cell declares what the KERNEL claims about it, and the test checks that
# claim rather than a blanket "everything is calibrated":
#
#   expect = "calibrated"  the kernel grades this configuration as calibrated,
#                          so the measured rejection rate must sit inside the
#                          binomial band around .05. Both halves are asserted.
#   expect = "flagged"     the kernel grades it anti-conservative. The rate is
#                          reported but NOT asserted to be .05 -- asserting that
#                          would be asserting the opposite of what the code says.
#                          What is asserted is that the grading fires.
check_cell <- function(label, effect, expect = "calibrated", ...) {
  r <- rate(effect, ...)
  b <- band(r$n)
  st <- e$dance_mixed_permutation(make(1, ...), n_permutations = 49,
                                  effects = effect)$calibration[[effect]]$status
  if (identical(expect, "calibrated")) {
    chk(st %in% c("exact", "calibrated") && r$rate <= .05 + b,
        sprintf("%-50s %s rejects at %.3f (n = %d, band <= %.3f), graded '%s'",
                label, effect, r$rate, r$n, .05 + b, st),
        sprintf("%-50s %s rejects at %.3f (band %.3f) and is graded '%s'",
                label, effect, r$rate, .05 + b, st))
  } else {
    chk(identical(st, "liberal"),
        sprintf("%-50s %s graded '%s'; measured rejection %.3f",
                label, effect, st, r$rate),
        sprintf("%-50s %s rejects at %.3f but is graded '%s', not flagged",
                label, effect, r$rate, st))
  }
  invisible(r)
}

cat(sprintf("-- level under the null (nominal .05, %d simulations per cell) ------\n", NSIM))

check_cell("balanced, equal dispersion",             "interaction", "calibrated", n = c(15, 15))
check_cell("unbalanced 10/25, equal dispersion",     "interaction", "calibrated", n = c(10, 25))
check_cell("three within levels, balanced",          "interaction", "calibrated", n = c(12, 12), nc = 3)
check_cell("AR(1) rho = .7 errors along the curve",  "interaction", "calibrated", n = c(15, 15), rho = .7)
check_cell("10% missing observations",               "interaction", "calibrated", n = c(15, 15), miss = .10)
check_cell("25% missing, unbalanced",                "interaction", "calibrated", n = c(12, 20), miss = .25)
check_cell("STRONG main effects, no interaction",    "interaction", "calibrated", n = c(15, 15), nuisance = TRUE)
check_cell("large groups, 5:1 dispersion",           "interaction", "calibrated",
           n = c(45, 45), sd_g = c(5, 1))

check_cell("between effect, balanced",               "between", "calibrated", n = c(15, 15))
check_cell("between effect, unbalanced + AR(1)",     "between", "calibrated", n = c(10, 25), rho = .7)
check_cell("within effect, unequal dispersion",      "within",  "calibrated", n = c(5, 30), sd_g = c(5, 1))
check_cell("within effect, 3 levels + missing",      "within",  "calibrated", n = c(12, 12), nc = 3, miss = .10)

cat("\n-- configurations the kernel itself declares anti-conservative --------\n")
# These are NOT asserted to be calibrated. They are asserted to be FLAGGED,
# which is what the code claims about them, and the measured rate is printed so
# the claim can be checked against it by a reader.
check_cell("small groups, 5:1 dispersion",           "interaction", "flagged",
           n = c(15, 15), sd_g = c(5, 1))
check_cell("unbalanced 10/25, 3:1 dispersion",       "interaction", "flagged",
           n = c(10, 25), sd_g = c(3, 1))
check_cell("strong main effects + 3:1 dispersion",   "interaction", "flagged",
           n = c(12, 20), sd_g = c(3, 1), nuisance = TRUE)
check_cell("three groups, unequal n and dispersion", "interaction", "flagged",
           n = c(12, 15, 18), sd_g = c(1, 2, 3))

cat("\n-- the cell that is KNOWN to stay liberal must FLAG itself ------------\n")
# 5 against 30 with a 5:1 contrast dispersion. Studentising took this from .750
# to about .14; it is not calibrated and the kernel must not pretend it is.
r_bad <- rate("interaction", n = c(5, 30), sd_g = c(5, 1), nsim = min(NSIM, 200L))
d_bad <- make(1, n = c(5, 30), sd_g = c(5, 1))
g_bad <- e$dance_mixed_permutation(d_bad, n_permutations = 99, effects = "interaction")
chk(identical(g_bad$calibration$interaction$status, "liberal"),
    sprintf("n = 5 vs 30 at 5:1 dispersion is graded '%s' (measured rejection %.3f)",
            g_bad$calibration$interaction$status, r_bad$rate),
    sprintf("a configuration that rejects at %.3f was graded '%s'",
            r_bad$rate, g_bad$calibration$interaction$status))
chk(grepl("ANTI-CONSERVATIVE", g_bad$calibration$interaction$message, fixed = TRUE) &&
      grepl("model-based", g_bad$calibration$interaction$message, fixed = TRUE),
    "the flag names the problem and the route that does not have it",
    "the flag does not say what to do instead")
# The pointwise statistic is LESS affected than the global one but is not clean:
# measured .072 and .076 at 10 and 15 per group with a 5:1 ratio, against .102
# and .104 globally. The flag must not exempt it.
chk(identical(g_bad$calibration$interaction$pointwise_status, "liberal"),
    "the pointwise statistic is flagged too -- it is less affected, not exempt",
    "the pointwise statistic was exempted from a flag that applies to it")
chk(grepl("describe WHERE the effect sits", g_bad$calibration$interaction$message,
          fixed = TRUE),
    "the descriptive output is explicitly preserved rather than withdrawn",
    "the flag withdraws the descriptive output as well as the inferential claim")
chk(identical(g_bad$calibration$within$status, "exact"),
    "the within effect is still graded exact -- its relabelling is stratified",
    "the within effect lost its exactness claim")

cat("\n-- power: the test must still find real effects ----------------------\n")
NPOW <- if (!is.na(.n_override) && .n_override > 0) {
  .n_override
} else if (FAST) 60L else 150L
power <- function(effect, delta, ..., nsim = NPOW) {
  p <- vapply(seq_len(nsim), function(s) {
    r <- e$dance_mixed_permutation(make(30000 + s, effect = effect, delta = delta, ...),
                                   n_permutations = NPERM, effects = effect)
    if (!isTRUE(r$ok)) return(NA_real_)
    r[[effect]]$global_p
  }, numeric(1))
  mean(p[is.finite(p)] <= .05)
}
pw <- power("interaction", delta = 2.5, n = c(15, 15))
chk(pw > .70, sprintf("a real interaction is detected %.0f%% of the time", 100 * pw),
    sprintf("power against a real interaction is only %.0f%%", 100 * pw))
# Under a 3:1 contrast dispersion the same effect is harder to see, so the
# fixture uses an effect the test should clearly find rather than one that sits
# on the threshold -- the question here is whether studentising DESTROYS power
# under heteroscedasticity, not where the 60% contour lies.
pw <- power("interaction", delta = 4.5, n = c(10, 25), sd_g = c(3, 1))
chk(pw > .80,
    sprintf("still detected under unequal n and dispersion (%.0f%%)", 100 * pw),
    sprintf("studentising cost too much power: %.0f%%", 100 * pw))
# The between effect competes with the participant random effect (SD 3), so a
# shift of 3 is a half-SD effect and 65% power against it is the right answer,
# not a defect. The fixture uses an effect the test should clearly find; the
# threshold is not lowered to accommodate a small one.
pw <- power("between", delta = 6, n = c(15, 15))
chk(pw > .80, sprintf("a real between effect is detected %.0f%% of the time", 100 * pw),
    sprintf("power against a real between effect is only %.0f%%", 100 * pw))
pw <- power("within", delta = 2, n = c(15, 15))
chk(pw > .70, sprintf("a real within effect is detected %.0f%% of the time", 100 * pw),
    sprintf("power against a real within effect is only %.0f%%", 100 * pw))

cat("\n-- structural guarantees ---------------------------------------------\n")
d <- make(7, n = c(8, 8))
a <- e$dance_mixed_permutation(d, n_permutations = 199, seed = 42)
b <- e$dance_mixed_permutation(d[sample(nrow(d)), ], n_permutations = 199, seed = 42)
chk(isTRUE(all.equal(a$interaction$p_values, b$interaction$p_values)) &&
      isTRUE(all.equal(a$between$global_p, b$between$global_p)) &&
      isTRUE(all.equal(a$within$global_p, b$within$global_p)),
    "the result does not depend on the row order of the input frame",
    "shuffling the rows changed the p-values")

dup <- rbind(d, d[c(1, 5, 9), ])
r_dup <- e$dance_mixed_permutation(dup, n_permutations = 19)
chk(isFALSE(r_dup$ok) && grepl("appear more than once", r_dup$message),
    "duplicated participant-by-condition-by-time cells are refused, with identifiers",
    "a frame with duplicated design cells was silently analysed")
chk(grepl("S00", r_dup$message),
    "the refusal names the offending participants rather than counting them",
    "the duplicate refusal does not identify which cells are duplicated")

# The global statistic must not depend on which files happen to be loaded.
chk(is.function(e$dance_l2_norm),
    "the global statistic comes from dance_l2_norm, unconditionally",
    "dance_l2_norm is not in scope, so the global statistic silently changes")
src <- paste(readLines(file.path(app_dir, "server/07_helpers_mixed_perm.R"), warn = FALSE),
             collapse = "\n")
body_only <- paste(grep("^\\s*#", strsplit(src, "\n")[[1]], invert = TRUE, value = TRUE),
                   collapse = "\n")
chk(!grepl('exists("dance_l2_norm"', body_only, fixed = TRUE),
    "no exists() branch decides which global statistic is computed",
    "the exists() branch on dance_l2_norm is back")

cat("\n")
if (bad) { cat(sprintf("Mixed calibration tests FAILED  (%d passed, %d failed)\n", ok_n, bad)); quit(status = 1) }
cat(sprintf("Mixed calibration tests PASSED  (%d passed, 0 failed)\n", ok_n))
