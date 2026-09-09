# ==============================================================================
# server/07_helpers_mixed_perm.R — exact permutation for a mixed design
# ==============================================================================
# The one-way kernels get their authority from exchangeability: under the null,
# relabelling is a symmetry of the data, so the permutation distribution is the
# exact null distribution. The mgcv route added at P17 gives up that property --
# its p-values condition on estimated smoothing parameters -- and for a module
# whose whole character is exact inference, that is a real loss.
#
# It turns out not to be necessary. All three effects in a mixed design have an
# exact scheme, including the interaction, which is the one usually said not to.
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
#   Z profiles, exactly.
#
# Calibrated by simulation rather than argued: under a null with NO interaction
# but strong main effects in BOTH factors, 1500 simulations at B = 499, the
# interaction test rejected at 0.053 against a nominal .05 (KS test against
# uniform, p = .758), while the two real main effects were detected at 1.000 and
# 0.950. tests/mixed_permutation_test.R re-runs that calibration.
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
dance_mixed_stats <- function(Y, group, usable) {
  nc <- dim(Y)[2]; nt <- dim(Y)[3]
  gl <- levels(group)

  # FAST PATH. With no missing cell the masks are all TRUE and every mean is a
  # plain column mean, so the loop below is pure overhead -- and this function is
  # called once per permutation, thousands of times. Handling NA correctly is
  # worth doing; paying for it when there is none is not.
  if (all(usable)) {
    subj_mean <- apply(Y, c(1, 3), mean)
    Z <- sweep(Y, c(1, 3), subj_mean, "-")
    ns <- dim(Y)[1]
    cond_mean <- apply(Z, c(2, 3), mean)
    s_within <- ns * colSums(cond_mean^2)
    grand <- colMeans(subj_mean)
    s_between <- numeric(nt); s_inter <- numeric(nt)
    for (g in gl) {
      idx <- which(group == g); k <- length(idx)
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
  subj_mean <- apply(Y, c(1, 3), mean)               # NA where a condition is absent
  Z <- sweep(Y, c(1, 3), subj_mean, "-")             # subject-centred profiles
  n_us <- colSums(usable)
  all_rows <- seq_len(dim(Y)[1])

  # WITHIN: condition means of the centred profiles
  cond_mean <- matrix(NA_real_, nc, nt)
  for (ci in seq_len(nc)) cond_mean[ci, ] <- wmean(Z[, ci, , drop = TRUE], all_rows)
  s_within <- n_us * colSums(cond_mean^2)

  # BETWEEN: group means of the subject means
  grand <- wmean(subj_mean, all_rows)
  s_between <- numeric(nt)
  for (g in gl) {
    idx <- which(group == g)
    ng_t <- colSums(usable[idx, , drop = FALSE])
    s_between <- s_between + ng_t * (wmean(subj_mean, idx) - grand)^2
  }

  # INTERACTION: group differences in the centred profiles, over conditions
  s_inter <- numeric(nt)
  for (g in gl) {
    idx <- which(group == g)
    ng_t <- colSums(usable[idx, , drop = FALSE])
    acc <- numeric(nt)
    for (ci in seq_len(nc))
      acc <- acc + (wmean(Z[, ci, , drop = TRUE], idx) - cond_mean[ci, ])^2
    s_inter <- s_inter + ng_t * acc
  }
  z <- function(v) { v[!is.finite(v)] <- 0; v }
  list(within = z(s_within), between = z(s_between), interaction = z(s_inter))
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
  if (!is.null(seed)) set.seed(seed)

  Y <- arr$Y; group <- arr$group
  ns <- dim(Y)[1]; nc <- dim(Y)[2]; nt <- dim(Y)[3]
  usable <- arr$usable
  obs <- dance_mixed_stats(Y, group, usable)

  l2 <- function(v) {
    if (exists("dance_l2_norm", mode = "function")) return(dance_l2_norm(v, arr$time))
    # trapezoid, so the global statistic does not scale with grid density
    sum(diff(arr$time) * (utils::head(v, -1) + utils::tail(v, -1)) / 2)
  }
  obs_l2 <- lapply(obs, l2)

  cnt   <- lapply(obs, function(v) numeric(length(v)))
  cnt_g <- stats::setNames(numeric(length(obs)), names(obs))

  for (b in seq_len(n_permutations)) {
    # WITHIN: relabel conditions inside each subject, independently
    if ("within" %in% effects) {
      Yp <- Y
      for (i in seq_len(ns)) Yp[i, , ] <- Y[i, sample.int(nc), , drop = FALSE]
      s <- dance_mixed_stats(Yp, group, usable)$within
      cnt$within <- cnt$within + (s >= obs$within)
      cnt_g["within"] <- cnt_g["within"] + (l2(s) >= obs_l2$within)
    }
    # BETWEEN and INTERACTION: relabel groups, the whole subject moving together
    if (any(c("between", "interaction") %in% effects)) {
      gp <- factor(sample(as.character(group)), levels = levels(group))
      s <- dance_mixed_stats(Y, gp, usable)
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
              p_floor = 1 / (n_permutations + 1))
  for (nm in names(obs)) {
    if (!nm %in% effects) next
    p  <- mkp(cnt[[nm]])
    pa <- stats::p.adjust(p, method = correction)
    out[[nm]] <- list(
      statistic = obs[[nm]], p_values = p, p_adjusted = pa,
      significant = !is.na(pa) & pa < alpha,
      n_significant = sum(!is.na(pa) & pa < alpha),
      global_statistic = obs_l2[[nm]],
      global_p = unname(mkp(cnt_g[[nm]])))
  }
  out
}
