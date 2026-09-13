# ==============================================================================
# server/08c_helpers_circstat.R — circular INFERENCE, out of the reactive scope
# ==============================================================================
# AUDIT (P21, finding A8). These five functions used to be defined at
# server/72_harmonic.R:416-583, INSIDE the server function, alongside the
# reactive outputs. Three consequences, all of which this file removes:
#
#   they could not be reached from a test or from an exported script without
#   slicing the source file by line number -- which is exactly what
#   tests/testthat/test-group-exclusion.R had to do;
#
#   server/73_cosinor_pairwise.R called watson_williams_test() across a file
#   boundary and it resolved only because both files happen to be sourced into
#   one environment and 73 happens to sort after 72; and
#
#   server/07_helpers_circular.R, whose NAME says circular helpers, holds only
#   the polar-density and ring geometry used by the plots, so a reader looking
#   for the circular tests did not find them there.
#
# Nothing here is a behaviour change: the bodies are the ones that were in the
# observer, de-indented and given the dance_ prefix the rest of the app's pure
# helpers use. tests/circular_inference_test.R pins them against the values the
# pre-extraction code produced.
#
# The file sorts after 08b so that dance_bivariate_manova() -- which
# dance_hotelling_t2() delegates to -- is already defined when a caller sources
# these two files alone.
#
# References: Mardia & Jupp (2000), Directional Statistics; Bingham, Arbogast,
# Cornelissen-Guillaume, Lee & Halberg (1982), Chronobiologia 9(4), 397-439.
# ==============================================================================


# Circular mean (returns radians)
dance_circular_mean <- function(angles_rad) {
  x <- mean(cos(angles_rad), na.rm = TRUE)
  y <- mean(sin(angles_rad), na.rm = TRUE)
  atan2(y, x)
}

# Mean resultant length (measure of concentration, 0-1)
dance_mean_resultant_length <- function(angles_rad) {
  x <- mean(cos(angles_rad), na.rm = TRUE)
  y <- mean(sin(angles_rad), na.rm = TRUE)
  sqrt(x^2 + y^2)
}

# Circular standard deviation (in same units as input)
dance_circular_sd <- function(angles_rad) {
  r_bar <- dance_mean_resultant_length(angles_rad)
  if(is.na(r_bar)) {
    return(NA)
  }
  if(r_bar > 0 && r_bar < 1) {
    sqrt(-2 * log(r_bar))
  } else if(r_bar >= 1) {
    0  # All points identical
  } else {
    NA  # Undefined
  }
}

# Circular standard error (approximate, based on von Mises)
dance_circular_se <- function(angles_rad) {
  n <- sum(!is.na(angles_rad))
  r_bar <- dance_mean_resultant_length(angles_rad)
  if(is.na(r_bar)) {
    return(NA)
  }
  if(r_bar > 0 && n > 1) {
    # Approximate SE for circular mean (Mardia & Jupp, 2000)
    1 / sqrt(n * r_bar^2)
  } else {
    NA
  }
}

# Watson-Williams test for comparing two or more groups of circular data
dance_watson_williams_test <- function(angles_list) {
  # angles_list: list of vectors, each containing angles in radians for one group
  k <- length(angles_list)
  if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))
  
  n <- sapply(angles_list, function(x) sum(!is.na(x)))
  N <- sum(n)
  
  if(any(n < 2)) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Each group needs at least 2 observations"))
  
  # Resultant lengths for each group
  R <- sapply(angles_list, function(x) {
    x <- x[!is.na(x)]
    sqrt(sum(cos(x))^2 + sum(sin(x))^2)
  })
  
  # Total resultant length (pooled)
  all_angles <- unlist(angles_list)
  all_angles <- all_angles[!is.na(all_angles)]
  R_total <- sqrt(sum(cos(all_angles))^2 + sum(sin(all_angles))^2)
  
  r_bar_total <- R_total / N
  
  if(r_bar_total < 0.45) {
    return(list(F = NA, df1 = k - 1, df2 = N - k, p = NA, r_bar = r_bar_total,
                message = "Warning: Data too dispersed (r̄ < 0.45). Consider non-parametric test."))
  }
  
  # Concentration parameter estimate
  kappa <- if(r_bar_total < 0.53) {
    2 * r_bar_total + r_bar_total^3 + 5 * r_bar_total^5 / 6
  } else if(r_bar_total < 0.85) {
    -0.4 + 1.39 * r_bar_total + 0.43 / (1 - r_bar_total)
  } else {
    1 / (r_bar_total^3 - 4 * r_bar_total^2 + 3 * r_bar_total)
  }
  
  g <- 1 - 1 / (3 * 8 * kappa^2)
  sum_R <- sum(R)
  F_stat <- g * (N - k) * (sum_R - R_total) / ((k - 1) * (N - sum_R))
  
  df1 <- k - 1
  df2 <- N - k
  p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
  
  list(F = F_stat, df1 = df1, df2 = df2, p = p_value, kappa = kappa, r_bar = r_bar_total, message = NULL)
}

# Hotelling's T² test on (beta_cos, beta_sin) pairs - amplitude-weighted acrophase comparison
# Tests whether the bivariate rhythmic vector (beta_cos, beta_sin) differs between groups.
# Because amplitude = sqrt(beta_cos² + beta_sin²), subjects with stronger rhythms carry
# more weight. For k>2 groups, a one-way MANOVA F approximation is used.
dance_hotelling_t2 <- function(beta_cos_list, beta_sin_list) {
  k <- length(beta_cos_list)
  if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))

  # Build per-group matrices
  mats <- lapply(seq_len(k), function(i) {
    x <- beta_cos_list[[i]]
    y <- beta_sin_list[[i]]
    ok <- complete.cases(x, y)
    cbind(x[ok], y[ok])
  })

  ns <- sapply(mats, nrow)
  if(any(ns < 3)) return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                              message = "Each group needs at least 3 observations"))
  N <- sum(ns)
  p <- 2  # two variables: beta_cos and beta_sin

  means <- lapply(mats, colMeans)

  # ------------------------------------------------------------------------
  # AUDIT (P20/R2). The k > 2 branch used to roll its own Wilks-to-F step and
  # got BOTH the statistic and its second degrees of freedom wrong, each by a
  # factor of two:
  #
  #     F   <- ((1 - sqrt(lambda))/sqrt(lambda)) * (df2 / df1)   with
  #     df1 <- p * (k - 1)          -- correct
  #     df2 <- N - k - p + 1        -- HALF the correct value
  #
  # For p = 2 response variables Rao's transformation is EXACT, not an
  # approximation, and its constant is s = sqrt((p^2 q^2 - 4)/(p^2 + q^2 - 5))
  # = 2 for every k >= 3, giving
  #
  #     df1 = 2(k - 1),  df2 = 2(N - k - 1),
  #     F   = ((1 - sqrt(lambda))/sqrt(lambda)) * (N - k - 1)/(k - 1).
  #
  # Halving df2 halves the multiplier as well, so the reported F was exactly
  # half the true F and was then referred to a distribution with half the
  # denominator degrees of freedom. Both errors push the same way: on a
  # three-group, twenty-per-group fixture the correct answer is F = 3.644338
  # on (4, 112) df, p = .0079, and the old code reported F = 1.822169 on
  # (4, 56) df, p = .1387 -- a real group difference reported as null.
  #
  # Rather than fix the algebra and leave a second implementation of a
  # standard test in the codebase to drift again, the whole thing now goes
  # through stats::manova() + summary(test = "Wilks"), which is R's own
  # reference implementation. It reduces to the exact two-sample Hotelling
  # T-squared when k = 2 (there q = 1, s = 1, df1 = 2, df2 = N - 3), so the
  # special case is no longer needed either.
  # ------------------------------------------------------------------------
  # The arithmetic is dance_bivariate_manova() in
  # server/08b_helpers_popcosinor.R, called rather than copied: the population
  # cosinor needs the identical test as the joint leg of Bingham's procedure,
  # and two implementations of one standard test is how they drift apart.
  Y <- do.call(rbind, mats)
  grp <- factor(rep(seq_len(k), times = ns))
  res <- dance_bivariate_manova(Y[, 1], Y[, 2], grp)
  if (!isTRUE(res$ok))
    return(list(F = NA, df1 = NA, df2 = NA, p = NA, lambda = NA,
                message = res$message %||% "The joint vector test could not be fitted."))

  list(F = as.numeric(res$F), df1 = res$df1, df2 = res$df2,
       p = as.numeric(res$p), lambda = as.numeric(res$lambda),
       n_groups = k, n_total = N, message = NULL)
}

# ==============================================================================
# BOOTSTRAP SUMMARIES THAT RESPECT WHAT IS BEING RESAMPLED (P21, findings A2/A3)
# ==============================================================================
# Two separate defects lived in the Harmonic Regression bootstrap, and both are
# arithmetic rather than judgement, so both are fixed here rather than argued
# about in the readout.
#
# A2 -- IT WAS NOT A BOOTSTRAP. The loop drew indices with replacement and then
# selected rows with `%in%`:
#
#     boot_idx <- sample(1:n_subjects, n_subjects, replace = TRUE)
#     boot_params <- individual_params[individual_params$subject %in% boot_idx, ]
#
# `%in%` is a set test, so a participant drawn three times contributes ONE row.
# What ran was a subsample WITHOUT replacement of the ~63.2% of participants
# that happened to be drawn at least once. That is not a milder bootstrap; it is
# a different estimator with a finite-population correction attached, and the
# correction makes the interval too NARROW. Measured on a 60-participant
# fixture, 4000 replicates: analytic SE of the mean 0.2543, participant
# bootstrap 0.2515, this procedure 0.1950 -- a 95% interval 21% too short.
#
# dance_boot_index() returns the index VECTOR, repeats included, so the caller
# indexes rows rather than filtering them.
#
# A3 -- THE PHASE INTERVAL USED A LINEAR QUANTILE. Acrophase was summarised with
# quantile() on hours. A rhythm peaking near midnight straddles the wrap, and
# there the linear quantile reports almost the whole cycle: measured on 4000
# draws centred at 0 h with SD 1.2 h, the naive interval was [0.07, 23.92] h --
# width 23.85 h -- against a circular interval of [-2.37, 2.36] h, width 4.73 h.
# Near midnight is the ordinary case in sleep research, not an edge case.
#
# dance_boot_circ_ci() takes the quantiles of the SIGNED ANGULAR DEVIATION from
# the bootstrap mean direction and maps them back, which is the standard
# construction and is invariant to where the origin happens to sit.
# ==============================================================================

# One bootstrap resample of the participant index. Repeats are kept: that is the
# entire difference between a bootstrap and a subsample.
dance_boot_index <- function(n) sample.int(n, n, replace = TRUE)

# Circular percentile interval for an angle in RADIANS.
# Returns the centre and the interval in radians, plus the half-width, which is
# the quantity that is meaningful when the interval straddles the origin.
dance_boot_circ_ci <- function(angles_rad, conf = 0.95) {
  a <- angles_rad[is.finite(angles_rad)]
  if (length(a) < 2)
    return(list(centre = NA_real_, lo = NA_real_, hi = NA_real_,
                half_width = NA_real_, n = length(a)))
  mu <- atan2(mean(sin(a)), mean(cos(a)))
  # signed deviation from the mean direction, wrapped to (-pi, pi]
  dev <- ((a - mu + pi) %% (2 * pi)) - pi
  p <- (1 - conf) / 2
  q <- stats::quantile(dev, c(p, 1 - p), names = FALSE)
  list(centre = mu %% (2 * pi),
       lo = (mu + q[1]) %% (2 * pi),
       hi = (mu + q[2]) %% (2 * pi),
       half_width = (q[2] - q[1]) / 2,
       width = q[2] - q[1],
       n = length(a))
}

# The same interval expressed in TIME on a harmonic's effective period P/h.
# `width` is carried separately because lo > hi whenever the interval wraps, and
# a reader subtracting the two endpoints would otherwise get a negative or
# near-full-cycle number.
dance_boot_circ_ci_time <- function(angles_rad, period = 24, harmonic = 1,
                                    conf = 0.95) {
  ci <- dance_boot_circ_ci(angles_rad, conf)
  k <- period / harmonic / (2 * pi)
  list(centre = ci$centre * k, lo = ci$lo * k, hi = ci$hi * k,
       width = ci$width * k, half_width = ci$half_width * k,
       wraps = is.finite(ci$lo) && is.finite(ci$hi) && ci$lo > ci$hi,
       effective_period = period / harmonic, n = ci$n)
}

# ------------------------------------------------------------------------------
# CLUSTER RESAMPLING WITH NEW CLUSTER IDENTITIES
# ------------------------------------------------------------------------------
# dance_boot_index() gives indices with replacement, which is right for a
# statistic computed by AVERAGING ROWS: a participant drawn three times
# contributes three times, and that is the resampling weight doing its job.
#
# It is NOT enough when the resampled data is REFITTED with the participant as a
# grouping factor. If participant 17 is drawn three times and all three copies
# keep the id "17", the model sees ONE participant with tripled observations --
# a cluster with three times the data and one random effect, instead of three
# independent clusters each with one. The random-effect variance is then
# estimated from fewer effective clusters than the bootstrap intended and the
# resulting intervals are wrong, usually too narrow.
#
# So every drawn copy gets a NEW id: 17 drawn three times becomes 17#1, 17#2,
# 17#3. Each copy carries that participant's COMPLETE set of rows -- all their
# conditions, sessions and time points -- because the unit being resampled is
# the participant, not the curve: breaking a participant apart would destroy the
# within-participant pairing the design is built on.
#
# Returns the row indices to take, the new participant ids for those rows, and
# the new curve ids (rebuilt from the new participant id) so a refit groups by
# something that actually distinguishes the copies.
dance_boot_clusters <- function(subject, curve = NULL, index = NULL) {
  subject <- as.character(subject)
  by_subj <- split(seq_along(subject), subject)
  ids <- names(by_subj)
  take <- if (is.null(index)) dance_boot_index(length(ids)) else index
  rows <- new_subj <- new_curve <- vector("list", length(take))
  for (k in seq_along(take)) {
    r <- by_subj[[ids[take[k]]]]
    rows[[k]] <- r
    tag <- sprintf("%s#%d", ids[take[k]], k)
    new_subj[[k]] <- rep(tag, length(r))
    new_curve[[k]] <- if (is.null(curve)) rep(tag, length(r)) else
      paste(tag, as.character(curve)[r], sep = " | ")
  }
  list(rows = unlist(rows, use.names = FALSE),
       subject = factor(unlist(new_subj, use.names = FALSE)),
       curve = factor(unlist(new_curve, use.names = FALSE)),
       n_clusters = length(take),
       n_distinct_drawn = length(unique(take)),
       note = paste(
         "Each drawn copy of a participant carries a NEW cluster id, so a",
         "participant drawn three times behaves as three independent clusters",
         "rather than as one participant with tripled observations. Every copy",
         "keeps that participant's complete set of conditions and sessions."))
}
