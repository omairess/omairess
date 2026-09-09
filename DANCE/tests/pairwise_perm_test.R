# ==============================================================================
# tests/pairwise_perm_test.R — the two post-hoc permutation kernels
# ==============================================================================
# P20/R1 and P20/R7. Both kernels in server/50_fanova.R had two defects that a
# reader of their output could not see:
#
#   R1  The GLOBAL p used the add-one Monte Carlo estimator (1 + #)/(1 + B),
#       and the POINTWISE p in the same result object used the plain proportion
#       # / B. The plain proportion can be exactly zero, which asserts an
#       impossible event, and at the app's smaller permutation counts it is
#       badly anti-conservative in the tail: with B = 49 it could report p = 0
#       at a point where the observed statistic was simply the largest of 50
#       exchangeable values, where the right answer is 1/50 = .02.
#
#   R7  Both wrote `t_stat[!is.finite(t_stat)] <- 0`, converting the single most
#       extreme possible result into the least. Twelve participants whose paired
#       difference is identical and non-zero have zero dispersion, so the paired
#       t is +Inf: unanimous separation. Coerced to zero it was reported as
#       t = 0, d = 0 and -- since every permuted statistic was compared against
#       zero as well -- p = 1.
#
# The paired sign-flip null is small enough to ENUMERATE at n = 8 (256 flips),
# so the Monte Carlo p is checked against the exact one rather than against
# another approximation.
#
# Run with:  Rscript tests/pairwise_perm_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressMessages(library(fda))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/01c_helpers_norm.R"), local = e)
source(file.path(app_dir, "server/06_helpers_posthoc.R"), local = e)

# The two kernels live inside the server function, so they are sliced out of the
# source rather than sourced. If this stops finding them the test fails loudly
# instead of quietly testing nothing.
src <- readLines(file.path(app_dir, "server/50_fanova.R"), warn = FALSE)
for (nm in c("perform_pairwise_comparisons_rm", "perform_pairwise_comparisons")) {
  i <- grep(paste0("^  ", nm, " <- function"), src)[1]
  stopifnot(!is.na(i))
  j <- i + which(src[(i + 1):length(src)] == "  }")[1]
  eval(parse(text = paste(src[i:j], collapse = "\n")), envir = e)
}
stopifnot(is.function(e$perform_pairwise_comparisons_rm),
          is.function(e$perform_pairwise_comparisons))

# ---------------------------------------------------------------- helpers ----
# Curves are handed to the kernels as an fd object; a dense B-spline basis over
# 0..1 reproduces whatever we put in closely enough for these checks.
as_fd <- function(M, n_time = 100) {          # M is time x curve
  tp <- seq(0, 1, length.out = nrow(M))
  b <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(nrow(M), 40))
  smooth.basis(tp, M, fdPar(b, 2, 1e-12))$fd
}

cat("-- P20/R1: no p-value below the permutation resolution floor ----------\n")
set.seed(4)
B <- 49L
n_per <- 9L; n_time <- 40L
M <- matrix(rnorm(n_time * 2 * n_per, 0, 1), n_time, 2 * n_per)
# a large, obvious difference so many points would report p = 0 under # / B
M[, (n_per + 1):(2 * n_per)] <- M[, (n_per + 1):(2 * n_per)] + 6
g <- factor(rep(c("a", "b"), each = n_per))
res <- e$perform_pairwise_comparisons(as_fd(M), g, n_permutations = B,
                                      correction_method = "none", alpha = .05)
pw <- res$results[[1]]
floor_b <- 1 / (B + 1)
chk(min(pw$p_values_pointwise) >= floor_b - 1e-12,
    sprintf("smallest raw pointwise p is %.4f, at or above the floor 1/(B+1) = %.4f",
            min(pw$p_values_pointwise), floor_b),
    sprintf("a raw pointwise p of %.4f is below the attainable floor %.4f",
            min(pw$p_values_pointwise), floor_b))
chk(!any(pw$p_values_pointwise == 0),
    "no pointwise p is exactly zero",
    "a pointwise p of exactly 0 was reported, which asserts an impossible event")
chk(pw$p_value_L2 >= floor_b - 1e-12 && identical(res$p_floor, floor_b),
    sprintf("the global p respects the same floor and the result records it (%.4f)", res$p_floor),
    "the global p and the recorded floor disagree with 1/(B+1)")
# the two conventions must now be the SAME convention
chk(all(abs(pw$p_values_pointwise * (B + 1) - round(pw$p_values_pointwise * (B + 1))) < 1e-9),
    "every pointwise p is a multiple of 1/(B+1), i.e. the add-one estimator",
    "the pointwise p-values are not on the (1 + k)/(1 + B) grid")

cat("\n-- P20/R7: unanimous separation is the most extreme result, not the least\n")
# Twelve participants, every paired difference identical and non-zero.
n_sub <- 12L
base <- matrix(rnorm(n_time * n_sub), n_time, n_sub)
delta <- 1.7
cond1 <- base
cond2 <- base - delta                          # every difference is exactly +1.7
M2 <- cbind(cond1, cond2)
subj <- factor(rep(sprintf("s%02d", 1:n_sub), 2))
cond <- factor(rep(c("A", "B"), each = n_sub))
res2 <- e$perform_pairwise_comparisons_rm(as_fd(M2), subj, cond, n_permutations = 199,
                                          correction_method = "none", alpha = .05)
p2 <- res2$results[[1]]
# The differences are constructed identical, but the curves reach the kernel
# through smooth.basis()/eval.fd(), and a B-spline round trip leaves a relative
# rounding of order 1e-16 in each evaluated value. So the dispersion is ~1e-16
# rather than exactly 0 and the statistic is ~1e16 rather than literally Inf.
# What the assertion has to capture is the property that failed before: the
# statistic must be LARGER THAN ANY PERMUTED ONE, not zero. The exactly-zero
# denominator is asserted separately, on dance_studentise(), below.
t_perm_max <- max(abs(p2$t_stat)) / 1e6      # any permuted t is many orders below
chk(all(p2$t_stat > t_perm_max) && all(p2$t_stat > 0),
    sprintf("the paired t is ~%.2g at every point -- unanimous separation, not 0",
            min(p2$t_stat)),
    sprintf("the paired t was coerced (range %.3g to %.3g)",
            min(p2$t_stat), max(p2$t_stat)))
chk(all(p2$cohens_d > max(abs(p2$cohens_d)) / 1e6) && all(p2$cohens_d > 0),
    sprintf("Cohen's d is ~%.2g at every point, not 0", min(p2$cohens_d)),
    sprintf("Cohen's d was coerced (range %.3g to %.3g)",
            min(p2$cohens_d), max(p2$cohens_d)))
chk(max(p2$p_values_pointwise) < 0.05,
    sprintf("the pointwise p is significant everywhere (max %.4f), not 1",
            max(p2$p_values_pointwise)),
    sprintf("unanimous separation gave a maximum pointwise p of %.4f",
            max(p2$p_values_pointwise)))
chk(abs(mean(p2$mean_diff) - delta) < 1e-6,
    sprintf("the mean difference is still reported correctly (%.4f)", mean(p2$mean_diff)),
    "the mean difference is wrong")

# 0/0 is genuinely undefined and must still be 0, not Inf: no difference and no
# dispersion is no signal.
chk(identical(e$dance_studentise(0, 0), 0),
    "0/0 stays 0 -- no difference and no dispersion is no signal",
    "0/0 was not treated as an absent effect")
chk(is.infinite(e$dance_studentise(1.7, 0)) && e$dance_studentise(1.7, 0) > 0,
    "a non-zero mean over zero dispersion is +Inf",
    "a degenerate denominator with a real numerator was not carried through")
chk(is.na(e$dance_studentise(NA_real_, 1)),
    "missing input stays NA rather than becoming a zero effect",
    "NA was silently turned into 0")
# And the exact +Inf case must compare correctly against a permutation set:
# only a permuted statistic that is ALSO infinite can match it, which under
# sign-flipping means the two all-same-sign flips out of 2^n.
t_obs_inf <- e$dance_studentise(rep(1.7, 12), 0)[1]
t_perm_mix <- c(rep(2.4, 30), Inf)
chk(sum(abs(t_perm_mix) >= abs(t_obs_inf)) == 1L,
    "an infinite observed statistic is exceeded only by an infinite permuted one",
    "the infinite statistic did not compare correctly against a permutation set")

cat("\n-- the Monte Carlo p against EXACT enumeration of the sign-flip null --\n")
# n = 8 participants: 2^8 = 256 sign patterns, so the paired randomisation
# distribution can be written down in full. The kernel's Monte Carlo p must
# agree with it to Monte Carlo error, at one time point.
n8 <- 8L
set.seed(19)
d8 <- c(1.4, -0.3, 0.9, 2.1, 0.2, 1.1, -0.7, 0.6)   # the paired differences
t_of <- function(v) { s <- sd(v) / sqrt(length(v)); if (s == 0) return(Inf * sign(mean(v))); mean(v) / s }
t_obs <- t_of(d8)
flips <- as.matrix(expand.grid(rep(list(c(-1, 1)), n8)))
t_all <- apply(flips, 1, function(f) t_of(d8 * f))
p_exact <- mean(abs(t_all) >= abs(t_obs))          # includes the identity flip
# the kernel's estimator on the same differences
mc <- function(B, seed) {
  set.seed(seed)
  cnt <- 0L
  for (b in seq_len(B)) {
    f <- sample(c(-1, 1), n8, replace = TRUE)
    if (abs(t_of(d8 * f)) >= abs(t_obs)) cnt <- cnt + 1L
  }
  e$dance_perm_p(cnt, B)
}
p_mc <- vapply(1:12, function(i) mc(4999L, 700 + i), numeric(1))
err <- abs(mean(p_mc) - p_exact)
chk(err < 3 * sqrt(p_exact * (1 - p_exact) / 4999),
    sprintf("Monte Carlo p averages %.5f against the exact %.5f (|diff| %.5f)",
            mean(p_mc), p_exact, err),
    sprintf("Monte Carlo p (%.5f) does not match enumeration (%.5f)", mean(p_mc), p_exact))
chk(all(p_mc >= 1 / 5000),
    "every Monte Carlo replicate respects the floor",
    "a Monte Carlo p fell below 1/(B+1)")
# the add-one estimator is conservative by construction, never below the truth
# by more than Monte Carlo error -- the plain proportion is what was not
chk(e$dance_perm_p(0, 49) == 1/50 && e$dance_perm_p(49, 49) == 1,
    "the estimator spans exactly [1/(B+1), 1]",
    "the estimator's range is wrong")

cat("\n-- both kernels share ONE implementation of the arithmetic ------------\n")
body_src <- paste(grep("^\\s*#", src, invert = TRUE, value = TRUE), collapse = "\n")
chk(!grepl("t_stat[!is.finite(t_stat)] <- 0", body_src, fixed = TRUE) &&
      !grepl("cohens_d[!is.finite(cohens_d)] <- 0", body_src, fixed = TRUE),
    "the non-finite-to-zero coercions are gone from both kernels",
    "a non-finite-to-zero coercion is still present")
chk(!grepl("mean(abs(t_stat_perm[t, ]) >= abs(t_stat[t]), na.rm = TRUE)", body_src, fixed = TRUE),
    "neither kernel computes a pointwise p as a plain proportion",
    "the plain-proportion pointwise p is back")
chk(length(gregexpr("dance_perm_p(", body_src, fixed = TRUE)[[1]]) >= 4,
    "both kernels route every p through the shared dance_perm_p()",
    "the shared p-value helper is not used in both kernels")

cat("\n")
if (bad) { cat(sprintf("Pairwise permutation tests FAILED  (%d passed, %d failed)\n", ok_n, bad)); quit(status = 1) }
cat(sprintf("Pairwise permutation tests PASSED  (%d passed, 0 failed)\n", ok_n))
