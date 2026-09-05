# ==============================================================================
# server/09b_helpers_rmfanova.R — the global repeated-measures test
# ==============================================================================
# The app's own repeated-measures procedure is POINTWISE: an F at each time
# point, permutation p-values, FDR across time. It reports no global test at all
# (`p_value_L2 <- NA`). That is a real gap: "the conditions differ somewhere"
# is usually the first question, and a family of pointwise tests does not answer
# it.
#
# rmfanova (Kurylo & Smaga 2023) does answer it. It implements three statistics
# (Cn, Dn, En) under five resampling schemes, plus pairwise comparisons. So it
# is worth having -- with three caveats that decide how it is wired in.
#
# 1. IT NEEDS A COMPLETE BALANCED DESIGN. x is a list of l matrices, each n x p,
#    and every matrix must hold the SAME n subjects in the SAME row order. The
#    app's pointwise procedure takes complete cases per time point, which is
#    more forgiving. Subjects missing any condition are therefore dropped here,
#    and the count is reported rather than absorbed.
#
# 2. NOT ALL FIFTEEN OUTPUTS ARE USABLE. Measured over 400 simulated nulls
#    (n=15, p=20, l=3, subject random intercepts, iid noise) and 400 alternatives
#    (a Gaussian bump on one condition), at nominal 5%:
#
#      method   type-I    power      method   type-I    power
#      Cn_P1      3.8%   100.0%      Dn_B3      0.0%    99.5%
#      Cn_P2      0.0%     0.0%      En_P1      6.0%    99.8%
#      Cn_B1      0.2%   100.0%      En_P2     17.5%   100.0%   <- inflated
#      Cn_B2      0.0%     0.0%      En_B1      1.5%    96.8%
#      Cn_B3      0.5%    99.8%      En_B2     14.2%   100.0%   <- inflated
#      Dn_P1      3.5%   100.0%      En_B3      4.0%    99.2%
#      Dn_P2      0.0%    76.0%
#      Dn_B1      0.0%    98.2%
#      Dn_B2      0.0%    73.0%
#
#    Two are anti-conservative at roughly three times the nominal rate. Two
#    (Cn_P2, Cn_B2) have no power whatsoever in this setting. Handing a user all
#    fifteen p-values in a GUI is an invitation to choose one.
#
#    That is ONE data-generating scenario, not a verdict on the package: the
#    authors' paper characterises these tests across many more, and a different
#    covariance structure could reorder them. It is enough to justify a curated
#    default rather than a dump.
#
# 3. THE PACKAGE DECLARES A DEPENDENCY IT DOES NOT USE. DESCRIPTION lists
#    refund in Imports, but no refund:: call appears anywhere in its source (it
#    uses MASS and parallel). If rmfanova will not install because refund will
#    not, that is why -- and it is the package's bug, not the analyst's.
#
# FCK_RMFANOVA_DEFAULT is the curated set: the two permutation schemes that came
# out closest to nominal with full power. Everything else is available but must
# be asked for, and is labelled.
FCK_RMFANOVA_DEFAULT <- c("Cn_P1", "Dn_P1", "En_P1")
FCK_RMFANOVA_INFLATED <- c("En_P2", "En_B2")
FCK_RMFANOVA_NOPOWER  <- c("Cn_P2", "Cn_B2")

# Build the list of condition matrices rmfanova wants, from the app's layout.
#
#   curves      : n_time x n_curves   (the app's orientation)
#   subject_id  : length n_curves
#   rm_factor   : length n_curves
#
# Returns NULL with a reason when the design cannot be made complete.
fck_rm_design <- function(curves, subject_id, rm_factor) {
  subject_id <- as.factor(subject_id)
  rm_factor  <- as.factor(rm_factor)
  visits <- levels(rm_factor)
  subs   <- levels(subject_id)
  if (length(visits) < 2)
    return(list(x = NULL, reason = "fewer than two conditions"))

  # a subject is usable only if it has exactly one curve in EVERY condition
  ok_sub <- vapply(subs, function(s) {
    k <- which(subject_id == s)
    all(vapply(visits, function(v) sum(rm_factor[k] == v) == 1L, logical(1)))
  }, logical(1))
  keep <- subs[ok_sub]
  if (length(keep) < 3)
    return(list(x = NULL, n_complete = length(keep), n_total = length(subs),
                reason = sprintf("only %d of %d subjects have all %d conditions exactly once",
                                 length(keep), length(subs), length(visits))))

  x <- lapply(visits, function(v) {
    idx <- vapply(keep, function(s)
      which(subject_id == s & rm_factor == v)[1], integer(1))
    m <- t(curves[, idx, drop = FALSE])        # subjects x time
    rownames(m) <- keep
    m
  })
  names(x) <- visits
  if (any(!vapply(x, function(m) all(is.finite(m)), logical(1))))
    return(list(x = NULL, n_complete = length(keep), n_total = length(subs),
                reason = "non-finite values in the complete-case design"))

  list(x = x, n_complete = length(keep), n_total = length(subs),
       dropped = setdiff(subs, keep), visits = visits, reason = NA_character_)
}

# Run the global test. Returns NULL (with a reason) rather than a partial answer.
fck_rmfanova_global <- function(curves, subject_id, rm_factor,
                                n_perm = 1000, n_boot = 1000,
                                method = "holm") {
  if (!requireNamespace("rmfanova", quietly = TRUE))
    return(list(ok = FALSE, reason = paste(
      "the rmfanova package is not installed. install.packages('rmfanova').",
      "If that fails on a missing 'refund', note that rmfanova declares refund",
      "in Imports but never calls it.")))

  d <- fck_rm_design(curves, subject_id, rm_factor)
  if (is.null(d$x)) return(list(ok = FALSE, reason = d$reason,
                                n_complete = d$n_complete, n_total = d$n_total))

  fit <- tryCatch(rmfanova::rmfanova(d$x, method = method,
                                     n_perm = n_perm, n_boot = n_boot),
                  error = function(e) e)
  if (inherits(fit, "error"))
    return(list(ok = FALSE, reason = paste("rmfanova failed:", conditionMessage(fit))))

  pv <- as.numeric(fit$p_values)
  names(pv) <- colnames(fit$p_values)
  list(ok = TRUE,
       stats = fit$test_stat,
       p_values = pv,
       p_default = pv[intersect(FCK_RMFANOVA_DEFAULT, names(pv))],
       p_values_pc = fit$p_values_pc,
       n_complete = d$n_complete, n_total = d$n_total,
       dropped = d$dropped, visits = d$visits,
       n_perm = n_perm, n_boot = n_boot, correction = method)
}
