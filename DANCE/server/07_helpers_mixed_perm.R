# ==============================================================================
# server/07_helpers_mixed_perm.R — permutation inference for a mixed design
# ==============================================================================
# The one-way kernels get their authority from exchangeability: under the null,
# relabelling is a symmetry of the data, so the permutation distribution is the
# exact null distribution. The mgcv route added at P17 gives up that property --
# its p-values condition on estimated smoothing parameters -- and for a module
# whose whole character is exact inference, that is a real loss.
#
# It turns out not to be necessary. All three effects in a mixed design have a
# relabelling scheme, the interaction included, which is the one usually said not
# to have one.
#
# WHAT "EXACT" DOES AND DOES NOT COVER (corrected at P20/R5). This header used to
# say all three schemes were exact. A relabelling is a symmetry only if the
# things being relabelled are EXCHANGEABLE, which requires their distributions to
# be identical -- not merely to have equal means. The within scheme relabels
# inside a participant and is exact whatever the groups look like. The two schemes
# that move participants ACROSS GROUPS are exact only when the groups are
# exchangeable, and unequal group sizes with unequal dispersion is precisely the
# case where they are not. Measured, on an interaction null with 5 against 30
# participants and a 5:1 dispersion ratio in the contrast: the unstudentised test
# rejected 225 of 300 times at a nominal .05. The statistic is Welch-type
# studentised below, and each run grades its own configuration -- see
# dance_mixed_calibration().
#
# Write Y[s, c, t] for subject s, condition c, time t, and
#     Z[s, c, t] = Y[s, c, t] - mean_c Y[s, c, t]
# for the subject-centred profile.
#
#   WITHIN main effect. Condition labels are exchangeable WITHIN a subject under
#   the null. Permuting them inside each subject leaves every subject's own
#   level untouched, so the scheme is valid whatever the between-subject
#   structure is.
#
#   BETWEEN main effect. Under the null, whole subjects are exchangeable across
#   groups. The subject's entire record moves together -- which is what the
#   one-way kernel already does, applied to the subject means over conditions.
#
#   INTERACTION. This is the one that is usually the problem, and the way out is
#   a reduction: an interaction is a difference BETWEEN GROUPS in the
#   within-subject contrast. Under the null of no interaction the subject-centred
#   profiles Z are identically distributed across groups -- whatever the two main
#   effects are doing, because centring removes the subject's level and the
#   within effect is common to all groups. So group labels can be permuted on the
#   Z profiles -- exactly when those profiles are identically distributed across
#   groups, asymptotically otherwise.
#
# Calibrated by simulation rather than argued. Under a null with NO interaction
# but strong main effects in BOTH factors and equal dispersions, 1500
# simulations at B = 499, the interaction test rejected at 0.053 against a
# nominal .05 (KS against uniform, p = .758), while the two real main effects
# were detected at 1.000 and 0.950; tests/mixed_permutation_test.R re-runs that.
# tests/mixed_calibration_test.R is the wider grid added at P20 -- balanced and
# unbalanced, equal and unequal dispersion, two and three within levels,
# correlated errors, missingness, strong nuisance effects -- and it asserts a
# binomial band rather than an exact .05, and asserts that the configuration
# known to remain liberal FLAGS ITSELF rather than that it is calibrated.
# DANCE_CALIB_NSIM=15 runs it as a wiring check in a fraction of the time; the
# numbers quoted below come from the default simulation count.
#
# Everything here is PURE. It takes arrays, not reactives.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---- reshape the long frame into the subject x condition x time array --------
# Returns NULL when the design is not complete: an exact permutation scheme
# needs every subject to have every condition, because the exchangeability
# argument is about relabelling WITHIN a subject. An unbalanced design is not a
# failure of the data, it is outside what this estimator can claim, and the
# caller is told to use the model-based route instead.
dance_mixed_array <- function(d) {
  subs <- levels(d$subject); cond <- levels(d$within); tt <- sort(unique(d$t))
  ns <- length(subs); nc <- length(cond); nt <- length(tt)
  Y <- array(NA_real_, c(ns, nc, nt), dimnames = list(subs, cond, NULL))
  si <- match(as.character(d$subject), subs)
  ci <- match(as.character(d$within), cond)
  ti <- match(d$t, tt)
  Y[cbind(si, ci, ti)] <- d$y
  grp <- vapply(subs, function(s) as.character(d$between[match(s, as.character(d$subject))]),
                character(1))

  # AUDIT (P19). This used to require !any(is.na(Y)) -- no missing observation
  # anywhere -- which is far stricter than the permutation argument needs and
  # refuses ordinary data (the dataset this was built for is 1.2% missing).
  #
  # What exactness actually requires is that RELABELLING be a symmetry. It is,
  # per time point, provided a participant contributes at a time point only when
  # they have EVERY condition observed there:
  #
  #   relabelling conditions inside a participant permutes that participant's own
  #   cells, so the set of time points at which they are complete is unchanged;
  #
  #   relabelling groups moves whole participants, so each participant's
  #   contribution travels with them intact.
  #
  # So the contributing set at each time point is invariant under both schemes,
  # and the statistics below are computed over it. What is genuinely required is
  # that every participant HAS every condition -- a missing cell is a missing
  # half of the contrast, and no relabelling can invent it.
  cell_present <- apply(Y, c(1, 2), function(v) any(is.finite(v)))
  usable <- apply(Y, c(1, 3), function(v) all(is.finite(v)))   # subjects x time
  list(Y = Y, group = factor(grp), subjects = subs, conditions = cond,
       time = tt, usable = usable,
       complete = all(cell_present),
       n_usable = colSums(usable),
       any_missing = any(!usable))
}

# ---- the three pointwise statistics -----------------------------------------
# Each is a between-group or between-condition sum of squares at one time point.
# They are not scaled by an error term: a permutation test does not need one,
# because the same statistic is recomputed under every relabelling and the
# scaling would be a monotone transformation shared by all of them.
# `usable[s, t]` marks the participant-by-time cells where every condition is
# observed. It is computed once from the ORIGINAL array and passed in, because it
# must not be recomputed from a permuted array: under condition relabelling the
# pattern is the same set of cells, and recomputing it would make the statistic
# depend on the relabelling in a way the null does not.
# ---- studentised statistics (R5) --------------------------------------------
# The unstudentised sums of squares below are exact ONLY under full
# exchangeability, which requires the group distributions to be identical, not
# merely to have equal means. They are not, and the failure is severe: with
# group sizes 5 and 30 and contrast dispersions 5 and 1, no true interaction,
# the unstudentised interaction test rejected 225 of 300 times at a nominal .05.
# Equal variances with unequal n was fine (.060), and unequal variances with
# equal n was fine (.050) -- it is the combination that breaks it, which is the
# classic Behrens-Fisher situation.
#
# The fix is a WELCH-TYPE STUDENTISED statistic, which is the standard route
# (Janssen 1997; Pauly, Brunner & Konietschke 2015): dividing by a
# group-specific variance estimate makes the permutation distribution
# asymptotically correct under heteroscedasticity while remaining exact when the
# groups really are exchangeable.
#
# HONEST STATUS, and it differs by effect:
#
#   within        EXACT. The relabelling is stratified inside each participant,
#                 so it is a randomisation test conditional on that participant;
#                 heterogeneity between groups cannot break a within-participant
#                 symmetry.
#
#   between,      ASYMPTOTICALLY VALID under unequal variances, exact under
#   interaction   exchangeability. Not exact in finite samples when the groups
#                 differ in dispersion. With very small groups the asymptotics
#                 are thin and the test is approximate; the readout says so.
#
# w_g = n_g / s2_g is the Welch weight; the weighted grand mean uses the same
# weights, which is what makes the statistic scale-free under the null.
dance_welch_between <- function(M, group, usable_rows = NULL) {
  gl <- levels(group); ng <- length(gl); nt <- ncol(M)
  W <- matrix(0, ng, nt); Mg <- matrix(NA_real_, ng, nt)

  # FAST PATH. With no missing cell every group contributes the same rows at
  # every time point, so the per-group mean and variance are column operations
  # on one submatrix instead of nt separate var() calls. This function runs once
  # per permutation per effect -- on the calibration grid, ~180,000 times -- so
  # the loop version dominated the runtime. The arithmetic is identical:
  # var() with denominator n - 1, computed as the centred sum of squares.
  complete <- is.null(usable_rows) || all(usable_rows)
  if (complete) {
    for (j in seq_len(ng)) {
      idx <- which(group == gl[j]); n <- length(idx)
      if (n < 2) next
      sub <- M[idx, , drop = FALSE]
      mu <- colMeans(sub)
      s2 <- colSums(sweep(sub, 2, mu, "-")^2) / (n - 1)
      good <- is.finite(s2) & s2 > 0 & is.finite(mu)
      W[j, good] <- n / s2[good]
      Mg[j, good] <- mu[good]
    }
  } else {
    for (j in seq_len(ng)) {
      idx <- which(group == gl[j])
      for (k in seq_len(nt)) {
        r <- idx[usable_rows[idx, k]]
        n <- length(r)
        if (n < 2) next
        v <- M[r, k]; mu <- mean(v); s2 <- stats::var(v)
        # a group with no variance at a time point carries no Welch weight; the
        # alternative (an infinite weight) would let one degenerate group decide
        # the statistic
        if (!is.finite(s2) || s2 <= 0) next
        W[j, k] <- n / s2; Mg[j, k] <- mu
      }
    }
  }

  ok <- W > 0 & is.finite(Mg)
  Wz <- W; Wz[!ok] <- 0
  Mz <- Mg; Mz[!ok] <- 0
  wsum <- colSums(Wz)
  gm <- ifelse(wsum > 0, colSums(Wz * Mz) / wsum, 0)
  num <- colSums(Wz * (Mz - rep(gm, each = ng))^2 * ok)
  # a time point where fewer than two groups carry a weight has no between-group
  # contrast to measure
  num[colSums(ok) < 2] <- 0
  unname(num)
}

# ---- the two quantities that a GROUP relabelling cannot change ---------------
# Permuting group labels does not touch Y, so the subject means and the
# subject-centred profiles are the same for every one of the B draws. They used
# to be recomputed inside dance_mixed_stats() on every draw, which is where most
# of the between/interaction runtime went. Caching them changes no arithmetic:
# the same numbers are formed once instead of B + 1 times. (The WITHIN scheme
# permutes Y itself, so it must not use the cache -- and does not.)
dance_mixed_prep <- function(Y, usable) {
  subj_mean <- apply(Y, c(1, 3), mean)
  list(subj_mean = subj_mean,
       Z = sweep(Y, c(1, 3), subj_mean, "-"),
       complete = all(usable))
}

# `which` restricts the work to the effects the caller will actually read. The
# statistics are unchanged; the ones not asked for are simply not formed.
dance_mixed_stats <- function(Y, group, usable, studentise = TRUE, prep = NULL,
                              which = c("within", "between", "interaction")) {
  nc <- dim(Y)[2]; nt <- dim(Y)[3]
  gl <- levels(group)
  want <- function(nm) nm %in% which

  # FAST PATH. With no missing cell the masks are all TRUE and every mean is a
  # plain column mean, so the loop below is pure overhead -- and this function is
  # called once per permutation, thousands of times. Handling NA correctly is
  # worth doing; paying for it when there is none is not.
  if (all(usable)) {
    if (is.null(prep)) prep <- dance_mixed_prep(Y, usable)
    subj_mean <- prep$subj_mean; Z <- prep$Z
    ns <- dim(Y)[1]
    s_within <- NULL
    if (want("within")) {
      cond_mean <- apply(Z, c(2, 3), mean)
      # WITHIN stays unstudentised: its relabelling is stratified within a
      # participant, so it is exact and needs no variance correction.
      s_within <- ns * colSums(cond_mean^2)
    }
    if (studentise) {
      s_between <- if (want("between")) dance_welch_between(subj_mean, group) else NULL
      s_inter <- NULL
      if (want("interaction")) {
        # the interaction contrast: for 2 conditions this is the difference curve,
        # and in general the centred profile of the FIRST condition carries the
        # same information once the profiles are centred to sum to zero
        s_inter <- numeric(nt)
        for (ci in seq_len(nc - 1L))
          s_inter <- s_inter + dance_welch_between(Z[, ci, , drop = TRUE], group)
      }
      return(list(within = s_within, between = s_between, interaction = s_inter))
    }
    if (is.null(s_within)) cond_mean <- apply(Z, c(2, 3), mean)
    grand <- colMeans(subj_mean)
    s_between <- numeric(nt); s_inter <- numeric(nt)
    for (g in gl) {
      idx <- base::which(group == g); k <- length(idx)
      s_between <- s_between + k * (colMeans(subj_mean[idx, , drop = FALSE]) - grand)^2
      gm <- apply(Z[idx, , , drop = FALSE], c(2, 3), mean)
      s_inter <- s_inter + k * colSums((gm - cond_mean)^2)
    }
    return(list(within = s_within, between = s_between, interaction = s_inter))
  }

  wmean <- function(M, rows) {          # column means over usable rows only
    out <- numeric(nt)
    for (k in seq_len(nt)) {
      r <- rows[usable[rows, k]]
      out[k] <- if (length(r)) mean(M[r, k]) else NA_real_
    }
    out
  }
  if (is.null(prep)) prep <- dance_mixed_prep(Y, usable)
  subj_mean <- prep$subj_mean                        # NA where a condition is absent
  Z <- prep$Z                                        # subject-centred profiles
  n_us <- colSums(usable)
  all_rows <- seq_len(dim(Y)[1])

  # WITHIN: condition means of the centred profiles
  cond_mean <- matrix(NA_real_, nc, nt)
  for (ci in seq_len(nc)) cond_mean[ci, ] <- wmean(Z[, ci, , drop = TRUE], all_rows)
  s_within <- n_us * colSums(cond_mean^2)

  if (studentise) {
    s_between <- dance_welch_between(subj_mean, group, usable)
    s_inter <- numeric(nt)
    for (ci in seq_len(nc - 1L))
      s_inter <- s_inter + dance_welch_between(Z[, ci, , drop = TRUE], group, usable)
    z0 <- function(v) { v[!is.finite(v)] <- 0; v }
    return(list(within = z0(s_within), between = z0(s_between), interaction = z0(s_inter)))
  }

  # BETWEEN: group means of the subject means
  grand <- wmean(subj_mean, all_rows)
  s_between <- numeric(nt)
  for (g in gl) {
    idx <- base::which(group == g)
    ng_t <- colSums(usable[idx, , drop = FALSE])
    s_between <- s_between + ng_t * (wmean(subj_mean, idx) - grand)^2
  }

  # INTERACTION: group differences in the centred profiles, over conditions
  s_inter <- numeric(nt)
  for (g in gl) {
    idx <- base::which(group == g)
    ng_t <- colSums(usable[idx, , drop = FALSE])
    acc <- numeric(nt)
    for (ci in seq_len(nc))
      acc <- acc + (wmean(Z[, ci, , drop = TRUE], idx) - cond_mean[ci, ])^2
    s_inter <- s_inter + ng_t * acc
  }
  z <- function(v) { v[!is.finite(v)] <- 0; v }
  list(within = z(s_within), between = z(s_between), interaction = z(s_inter))
}

# ==============================================================================
# HOW FAR THE STUDENTISED PERMUTATION CAN BE TRUSTED (P20/R5)
# ==============================================================================
# Studentising fixed most of the miscalibration, not all of it, and pretending
# otherwise would repeat the original mistake in a quieter voice. What follows is
# measured, on a two-group interaction null with no true interaction, 8 time
# points, B = 199, and the two statistics this kernel actually reports: the
# global sqrt(integral S(t)^2 dt) and a single pointwise S(t). Rejection rate at
# a nominal .05; Monte Carlo 95% half-width about +/- .014 at 1000 simulations,
# +/- .021 at 500, +/- .035 at 150.
#
#   group n   dispersion   STUDENTISED (what runs)   UNSTUDENTISED   sims
#             ratio        global      pointwise     global
#   -------------------------------------------------------------------
#     5 / 30     5 : 1       .253          -            .698        200/500
#    10 / 10     5 : 1       .102         .072            -          500
#    15 / 15     5 : 1       .104         .076            -          500
#    15 / 15     2 : 1       .053           -             -         1000
#    20 / 20     1 : 1       .052         .059          .052        1000/500
#    20 / 20     2 : 1       .058 / .065  .066            -         500/1000
#    20 / 20     3 : 1       .066 / .073  .064            -         500/1000
#    20 / 20     5 : 1       .082         .058          .064        1000/500
#    30 / 30     5 : 1       .060 / .070  .058            -         500/1000
#    40 / 40     3 : 1       .052           -             -         1000
#    40 / 40     5 : 1       .056         .041          .042        1000/500
#    50 / 50     5 : 1       .041           -             -         1000
#    80 / 80     5 : 1       .047         .046            -         1000
#    10 / 25     3 : 1       .098         .066            -          500
#    12 / 20     3 : 1       .113           -             -          150  (+ strong
#                                                                     main effects)
#
# WHAT THE TABLE SETTLES.
#
#   The excess appears ONLY when the groups' dispersions differ. At a 1:1 ratio
#   both statistics sit at the nominal level at every size tried.
#
#   The excess SHRINKS WITH GROUP SIZE and is gone by about 40 participants per
#   group even at a 5:1 ratio (.056, .041, .047 at n = 40, 50, 80). That is the
#   signature of a procedure that is asymptotically valid and finite-sample
#   liberal, which is exactly what a studentised permutation is (Janssen 1997;
#   Pauly, Brunner & Konietschke 2015).
#
#   STUDENTISING IS THE RIGHT DEFAULT, AND IT IS NOT FREE. The unstudentised
#   statistic is FINE at equal group sizes even at a 5:1 ratio (.064 at 20 per
#   group) and marginally better than the studentised one there (.082). What it
#   cannot survive is unequal n TOGETHER WITH unequal dispersion -- the classic
#   Behrens-Fisher case -- where it rejected .698 of null datasets at 5 against
#   30. Studentising takes that to .253: an enormous improvement and still not a
#   usable test, which is why that configuration is flagged rather than
#   presented. Paying two points of level in the balanced case to avoid a
#   seventy-point failure in the unbalanced one is the trade this kernel makes,
#   deliberately.
#
#   The GLOBAL statistic is the worse of the two in the bad corner -- .102
#   against .072 at 10 per group, .104 against .076 at 15, .098 against .066 at
#   10/25 -- because sqrt(integral S^2) is dominated by the largest pointwise
#   statistic and a tail converges more slowly than a centre. But the POINTWISE
#   statistic is NOT clean there either: .072 and .076 are above nominal, and
#   an earlier draft of this note that called it "calibrated" was reading four
#   large-n cells and over-generalising from them. Both are reported as what
#   they are.
#
# So the test is not offered as though it were calibrated everywhere. It grades
# ITS OWN configuration against that table. Nothing is suppressed -- the
# statistic and the pointwise curves stay, and they are the right description of
# WHERE the effect sits -- but in the flagged regime the p-values are labelled
# anti-conservative and the model-based estimator is named as the route that
# does not have this failure mode, because it fits the dispersions rather than
# permuting across them.
dance_mixed_calibration <- function(min_n, dispersion_ratio, effect) {
  if (identical(effect, "within"))
    return(list(status = "exact", pointwise_status = "exact", message = paste(
      "Exact, pointwise and globally. Condition labels are permuted INSIDE each",
      "participant, so the reference distribution is a randomisation",
      "distribution conditional on that participant. Differences in dispersion",
      "between the between-participant groups cannot disturb a",
      "within-participant symmetry.")))

  r <- dispersion_ratio
  shown <- if (is.finite(r)) sprintf("%.1f", r) else "not estimable"
  thr <- dance_ratio_null_q95(min_n, effect)
  # Read off the table: an excess appears only when the dispersions differ, and
  # it is gone by 40 participants per group even at 5:1. Between those two the
  # rule is deliberately conservative -- it flags cells measured at .060 to .073
  # rather than waving them through on a confidence interval that includes .05.
  # It fires on EVIDENCE of unequal dispersion, not on a point estimate: see
  # dance_ratio_null_q95() for why a bare "ratio > 2" was wrong.
  ok <- !is.finite(r) || r <= thr || min_n >= 40
  if (ok) return(list(
    status = "calibrated", pointwise_status = "calibrated",
    message = paste0(
      "Asymptotically valid under unequal variances (Welch-type studentised ",
      "permutation; Janssen 1997, Pauly, Brunner & Konietschke 2015), and exact ",
      "if the groups are exchangeable. This configuration -- smallest group ",
      min_n, ", observed dispersion ratio ", shown, " -- is inside the range ",
      "where simulation put the rejection rate at the nominal level.")))

  severity <- if (min_n >= 20)
    paste("about .06 to .08 for the global test and .06 to .07 pointwise,",
          "against a nominal .05")
  else if (min_n >= 10)
    paste("about .10 for the global test and .07 pointwise, against a",
          "nominal .05")
  else
    paste("about .25 for the global test at 5 participants against 30 with a",
          "5:1 ratio -- a quarter of null datasets rejected at a nominal .05,",
          "which is not a usable test at all")
  list(
    status = "liberal", pointwise_status = "liberal",
    message = paste0(
      "CAUTION: the p-values for this effect are ANTI-CONSERVATIVE here. The ",
      "smallest group has ", min_n, " participants and the groups' dispersions ",
      "differ by a factor of ", shown, ". A small group whose variance is ",
      "estimated from few participants gets an over-large Welch weight, and ",
      "permuting the labels mixes the dispersions, so the observed ",
      "configuration is systematically more extreme than its own reference ",
      "set. At this size and ratio, simulation put the rejection rate at ",
      severity, "; the global statistic is the worse of the two, because it is ",
      "dominated by the largest pointwise value. The statistic and the ",
      "pointwise curve still describe WHERE the effect sits. For the ",
      "inferential claim, use the model-based estimator, which fits the ",
      "dispersions rather than permuting across them. The excess also ",
      "disappears with sample size -- it was gone by 40 participants per group ",
      "even at a 5:1 ratio.",
      if (min_n < 10)
        paste0(" At this group size the global p-value should not be reported ",
               "as an inferential result at all; take it from the model-based ",
               "estimator.")
      else "")
  )
}

# ------------------------------------------------------------------------------
# WHEN IS AN OBSERVED DISPERSION RATIO EVIDENCE OF ANYTHING? (P20/R5)
# ------------------------------------------------------------------------------
# The grading above needs to know whether the groups' dispersions REALLY differ,
# and all it has is a ratio of sample variances, which is itself noisy. A first
# version of this rule flagged any observed ratio above 2, and that was wrong in
# a way that showed up immediately: on data with genuinely EQUAL dispersions and
# 15 participants per group, the between-participant effect was graded
# anti-conservative, because the ratio of two variances each estimated from 15
# subject means is above 2 about a fifth of the time by chance alone.
#
# So the threshold is the 95th percentile of the ratio UNDER EQUAL TRUE
# DISPERSIONS, measured rather than assumed: 400 simulated datasets per cell, 8
# time points, 2 conditions, equal group sizes.
#
#   n per group      BETWEEN (subject means)      INTERACTION (contrasts)
#                   median   q90    q95           median   q90    q95
#   -----------------------------------------------------------------
#         5          2.07    5.24   6.35           1.25    1.70   1.94
#         8          1.64    3.36   4.02           1.18    1.53   1.68
#        10          1.49    2.81   3.43           1.17    1.46   1.65
#        15          1.37    2.40   2.77           1.14    1.40   1.47
#        20          1.34    2.10   2.31           1.11    1.33   1.38
#        30          1.28    1.77   1.94           1.09    1.23   1.26
#
# The two effects are on completely different scales of noise, which is the
# point: the BETWEEN ratio compares variances of one number per participant, so
# at 5 per group a ratio of 6 is an ordinary null outcome; the INTERACTION ratio
# compares variances of a curve averaged over 8 time points, so its null spread
# is much tighter and a ratio of 2 there really is evidence. One threshold for
# both would have to be either useless for one or wrong for the other.
#
# These quantiles are for 8 time points. More time points make the interaction
# ratio less noisy still, so using this table then errs toward NOT flagging,
# which is the direction that does not manufacture false alarms.
dance_ratio_null_q95 <- function(min_n, effect) {
  n  <- c(5, 8, 10, 15, 20, 30)
  q  <- if (identical(effect, "between")) c(6.35, 4.02, 3.43, 2.77, 2.31, 1.94)
        else                              c(1.94, 1.68, 1.65, 1.47, 1.38, 1.26)
  if (!is.finite(min_n) || min_n <= n[1]) return(q[1])
  if (min_n >= n[length(n)]) return(q[length(q)])
  stats::approx(n, q, xout = min_n)$y
}

# The dispersion ratio the grading above uses: the largest over the smallest
# group variance of the quantity the effect is actually computed from, averaged
# over time points so one noisy point cannot decide it.
dance_mixed_dispersion_ratio <- function(M, group) {
  gl <- levels(group)
  v <- vapply(gl, function(g) {
    idx <- which(group == g)
    if (length(idx) < 2) return(NA_real_)
    mean(apply(M[idx, , drop = FALSE], 2, stats::var), na.rm = TRUE)
  }, numeric(1))
  v <- v[is.finite(v) & v > 0]
  if (length(v) < 2) return(NA_real_)
  max(v) / min(v)
}

# ---- the permutation test ----------------------------------------------------
# `time` is passed so the global statistic can be integrated by trapezoid rather
# than summed over grid points -- summing makes the global statistic scale with
# grid density, which is the defect dance_l2_norm() exists to avoid.
dance_mixed_permutation <- function(d, n_permutations = 999, alpha = 0.05,
                                    correction = "BH", seed = NULL,
                                    effects = c("within", "between", "interaction")) {
  bad <- dance_mixed_check(d)
  if (length(bad)) return(list(ok = FALSE, message = paste(bad, collapse = " ")))
  arr <- dance_mixed_array(d)
  if (!isTRUE(arr$complete))
    return(list(ok = FALSE, message = paste(
      "The exact permutation scheme needs every participant to have every level of",
      "the within-subject factor, because its validity comes from relabelling",
      "inside a participant, and a missing cell is a missing half of the contrast.",
      "At least one participant is missing a condition entirely. Use the",
      "model-based estimator, which handles an unbalanced design, and report it",
      "as approximate.")))
  # P20/R6: a permutation result is a draw from a random procedure, so the
  # state it was drawn under is part of the result. Recording it makes the run
  # reproducible after the fact -- assign the returned rng_state back to
  # .Random.seed and re-run to get the same p-values -- rather than only
  # reproducible if the user happened to set a seed first.
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    set.seed(NULL)
  rng_state <- get(".Random.seed", envir = globalenv())
  rng_kind <- RNGkind()

  Y <- arr$Y; group <- arr$group
  ns <- dim(Y)[1]; nc <- dim(Y)[2]; nt <- dim(Y)[3]
  usable <- arr$usable
  # subject means and centred profiles: invariant under GROUP relabelling, so
  # formed once and reused for all B draws of the between/interaction schemes
  prep0 <- dance_mixed_prep(Y, usable)
  obs <- dance_mixed_stats(Y, group, usable, prep = prep0)
  grp_effects <- intersect(c("between", "interaction"), effects)

  # P20/R5: grade this configuration against the simulated calibration table
  # above, per effect, from the quantity each effect is actually computed from.
  min_n <- min(as.integer(table(group)))
  ratio_between <- dance_mixed_dispersion_ratio(prep0$subj_mean, group)
  ratio_inter <- dance_mixed_dispersion_ratio(prep0$Z[, 1, , drop = TRUE], group)
  calibration <- list(
    within      = dance_mixed_calibration(min_n, NA_real_, "within"),
    between     = dance_mixed_calibration(min_n, ratio_between, "between"),
    interaction = dance_mixed_calibration(min_n, ratio_inter, "interaction"))
  calibration$within$dispersion_ratio <- NA_real_
  calibration$between$dispersion_ratio <- ratio_between
  calibration$interaction$dispersion_ratio <- ratio_inter
  for (nm in names(calibration)) calibration[[nm]]$min_group_n <- min_n

  # AUDIT (P20/R9). This used to be
  #     if (exists("dance_l2_norm", mode = "function")) dance_l2_norm(v, time)
  #     else sum(diff(time) * (head(v,-1) + tail(v,-1)) / 2)
  # -- and the two branches compute DIFFERENT statistics: sqrt(integral v^2 dt)
  # against integral v dt. Which one ran depended on whether a file later in the
  # sourcing order happened to be loaded, so the same data at the same seed gave
  # p = .09 in the app and p = .20 in a standalone script. dance_l2_norm now
  # lives in server/01c_helpers_norm.R, ahead of every caller, so there is one
  # global statistic and no branch to disagree with itself.
  l2 <- function(v) dance_l2_norm(v, arr$time)
  obs_l2 <- lapply(obs, function(v) if (is.null(v)) NA_real_ else l2(v))

  cnt   <- lapply(obs, function(v) numeric(length(v)))
  names(cnt) <- names(obs)
  cnt_g <- stats::setNames(numeric(length(obs)), names(obs))

  for (b in seq_len(n_permutations)) {
    # WITHIN: relabel conditions inside each subject, independently
    if ("within" %in% effects) {
      Yp <- Y
      for (i in seq_len(ns)) Yp[i, , ] <- Y[i, sample.int(nc), , drop = FALSE]
      # no `prep`: the WITHIN scheme permutes Y itself, so the cached means do
      # not apply and must be recomputed from the permuted array
      s <- dance_mixed_stats(Yp, group, usable, which = "within")$within
      cnt$within <- cnt$within + (s >= obs$within)
      cnt_g["within"] <- cnt_g["within"] + (l2(s) >= obs_l2$within)
    }
    # BETWEEN and INTERACTION: relabel groups, the whole subject moving together
    if (length(grp_effects)) {
      gp <- factor(sample(as.character(group)), levels = levels(group))
      s <- dance_mixed_stats(Y, gp, usable, prep = prep0, which = grp_effects)
      if ("between" %in% effects) {
        cnt$between <- cnt$between + (s$between >= obs$between)
        cnt_g["between"] <- cnt_g["between"] + (l2(s$between) >= obs_l2$between)
      }
      if ("interaction" %in% effects) {
        cnt$interaction <- cnt$interaction + (s$interaction >= obs$interaction)
        cnt_g["interaction"] <- cnt_g["interaction"] + (l2(s$interaction) >= obs_l2$interaction)
      }
    }
  }

  # (1 + #{T* >= T}) / (1 + B): a Monte Carlo p can never be exactly zero
  mkp <- function(x) (1 + x) / (1 + n_permutations)
  out <- list(ok = TRUE, effects = effects, n_permutations = n_permutations,
              alpha = alpha, correction = correction, time = arr$time,
              n_subjects = ns, n_conditions = nc,
              any_missing = arr$any_missing, n_usable = arr$n_usable,
              # P19: the shared readout describes every mixed result the same
              # way, so every kernel must return the same design summary. Leaving
              # it out here is what made the readout fail on a permutation
              # result: `if (b$n_partial > 0)` on a NULL is a length-zero error.
              balance = dance_mixed_balance(d),
              group_sizes = as.integer(table(group)),
              group_names = levels(group), condition_names = arr$conditions,
              # P20/R5: the honest status of each effect's p-values, so the
              # readout, the report and the export all say the same thing about
              # what this procedure does and does not establish.
              calibration = calibration,
              min_group_n = min_n,
              seed = seed, rng_state = rng_state, rng_kind = rng_kind,
              # P20/R12: what the pointwise correction is applied across.
              multiplicity_family = sprintf(
                "%s across the %d evaluation points, separately within each effect; the global L2 test is a single test per effect and is not corrected",
                correction, nt),
              p_floor = 1 / (n_permutations + 1))
  for (nm in names(obs)) {
    if (!nm %in% effects) next
    p  <- mkp(cnt[[nm]])
    pa <- stats::p.adjust(p, method = correction)
    out[[nm]] <- list(
      calibration = calibration[[nm]],
      statistic = obs[[nm]], p_values = p, p_adjusted = pa,
      significant = !is.na(pa) & pa < alpha,
      n_significant = sum(!is.na(pa) & pa < alpha),
      global_statistic = obs_l2[[nm]],
      global_p = unname(mkp(cnt_g[[nm]])))
  }
  out
}
