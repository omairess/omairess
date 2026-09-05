# ==============================================================================
# server/09_helpers_pcanova.R — per-component group tests for fPCA scores
#
# WHY THIS EXISTS
# ---------------
# fPCA and warped fPCA produce an n x k matrix of component scores. Comparing
# groups on those scores is the natural next question ("do the age bands differ
# in the dominant mode of variation?") and the app had no way to ask it.
#
# The statistics are ordinary one-way ANOVA per column, but four things make
# this NOT an ordinary ANOVA table, and each of them is handled explicitly
# rather than left to the reader:
#
#  1. TWO families of multiplicity, not one. The pairwise comparisons within a
#     component are one family; the omnibus tests across components are
#     another. Correcting the first and ignoring the second is the usual error
#     and it is the more serious of the two, because k omnibus tests are what
#     produce a "significant PC" by chance. Both are offered, separately.
#
#  2. The scores are orthogonal IN SAMPLE but the tests are not independent:
#     every column depends on the same estimated eigenfunctions, fitted to the
#     same curves. So an across-component Bonferroni is conservative and a BH
#     assumes a positive-dependence structure that holds only approximately.
#     Reported, not asserted away.
#
#  3. Eigenvector stability. When two eigenvalues are close, the rotation
#     between their eigenfunctions is barely determined, and "a group
#     difference on PC2" may be a group difference on an arbitrary mixture of
#     PC2 and PC3. fck_pc_separation() quantifies this so a finding on a poorly
#     separated component can be discounted.
#
#  4. Sign indeterminacy. A PC and its negation are the same component, so the
#     SIGN of a group mean difference is meaningful only relative to the
#     plotted loading. Every direction statement is phrased against the loading.
#
# Pure base R + stats, no Shiny, so it is testable without a session.
# ==============================================================================


# ---------------------------------------------------------------- corrections

# The corrections offered. p.adjust covers most; Tukey and Games-Howell are
# different animals (they adjust the null distribution, not the p-value) and
# are implemented separately below.
FCK_PC_CORRECTIONS <- c(
  "None (raw p)"                              = "none",
  "Bonferroni"                                = "bonferroni",
  "Holm (step-down Bonferroni)"               = "holm",
  "Hochberg (step-up)"                        = "hochberg",
  "Hommel"                                    = "hommel",
  "Benjamini-Hochberg (FDR)"                  = "BH",
  "Benjamini-Yekutieli (FDR, any dependence)" = "BY",
  "Tukey HSD (equal variances)"               = "tukey",
  "Games-Howell (unequal variances)"          = "games-howell"
)

fck_pc_adjust <- function(p, method) {
  if (identical(method, "none")) return(p)
  if (method %in% c("tukey", "games-howell")) return(p)   # already adjusted
  stats::p.adjust(p, method = method)
}


# ------------------------------------------------------------ assumption checks

# Brown-Forsythe (Levene on medians): robust to non-normality, which matters
# because PC scores are frequently heavy-tailed even when the curves are not.
fck_brown_forsythe <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  if (nlevels(g) < 2 || length(x) < nlevels(g) + 1) return(NULL)
  med <- tapply(x, g, stats::median)
  z <- abs(x - med[as.character(g)])
  fit <- stats::anova(stats::lm(z ~ g))
  list(F = fit$`F value`[1], df1 = fit$Df[1], df2 = fit$Df[2],
       p = fit$`Pr(>F)`[1])
}

# Normality of the residuals. Shapiro-Wilk is capped at n = 5000 and is absurdly
# powerful at n > 1000 -- it will reject a trivially non-normal sample that
# ANOVA handles perfectly well. So skewness and excess kurtosis are reported
# alongside, and the verdict leans on those.
fck_normality_check <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  if (length(x) < 8) return(NULL)
  r <- x - tapply(x, g, mean)[as.character(g)]
  n <- length(r); s <- sqrt(mean(r^2))
  skew <- if (s > 0) mean(r^3) / s^3 else NA_real_
  kurt <- if (s > 0) mean(r^4) / s^4 - 3 else NA_real_
  sw <- if (n >= 3 && n <= 5000)
    tryCatch(stats::shapiro.test(r), error = function(e) NULL) else NULL
  list(n = n, skew = skew, excess_kurtosis = kurt,
       shapiro_W = if (is.null(sw)) NA_real_ else unname(sw$statistic),
       shapiro_p = if (is.null(sw)) NA_real_ else sw$p.value,
       shapiro_applicable = n <= 5000,
       # the practical rule: ANOVA's F is robust to this much
       severe = isTRUE(abs(skew) > 2) || isTRUE(abs(kurt) > 7))
}


# ------------------------------------------------------------------- omnibus

# One-way test on a single component's scores. Returns Fisher's F, Welch's F
# (which does not assume equal variances) and Kruskal-Wallis, plus eta^2,
# omega^2 and epsilon^2 -- because which one is appropriate depends on checks
# the caller has not yet made, and computing only one of them forces a choice
# before the evidence is in.
fck_pc_omnibus <- function(scores, g, conf = 0.95) {
  ok <- is.finite(scores) & !is.na(g)
  x <- scores[ok]; g <- droplevels(as.factor(g[ok]))
  k <- nlevels(g); n <- length(x)
  if (k < 2 || n < k + 2) return(NULL)

  means <- tapply(x, g, mean); sds <- tapply(x, g, stats::sd)
  ns <- tapply(x, g, length)

  grand <- mean(x)
  ss_b <- sum(ns * (means - grand)^2)
  ss_t <- sum((x - grand)^2)
  ss_w <- ss_t - ss_b
  df1 <- k - 1; df2 <- n - k
  ms_w <- ss_w / df2
  Ff <- if (ms_w > 0) (ss_b / df1) / ms_w else NA_real_
  pf <- if (is.finite(Ff)) stats::pf(Ff, df1, df2, lower.tail = FALSE) else NA_real_

  eta2 <- if (ss_t > 0) ss_b / ss_t else NA_real_
  omega2 <- if (ss_t + ms_w > 0) (ss_b - df1 * ms_w) / (ss_t + ms_w) else NA_real_

  welch <- tryCatch(stats::oneway.test(x ~ g, var.equal = FALSE),
                    error = function(e) NULL)
  kw <- tryCatch(stats::kruskal.test(x ~ g), error = function(e) NULL)
  eps2 <- if (!is.null(kw)) (unname(kw$statistic) - k + 1) / (n - k) else NA_real_

  list(k = k, n = n, levels = levels(g), means = means, sds = sds, ns = ns,
       F = Ff, df1 = df1, df2 = df2, p = pf,
       welch_F = if (is.null(welch)) NA_real_ else unname(welch$statistic),
       welch_df1 = if (is.null(welch)) NA_real_ else unname(welch$parameter[1]),
       welch_df2 = if (is.null(welch)) NA_real_ else unname(welch$parameter[2]),
       welch_p = if (is.null(welch)) NA_real_ else welch$p.value,
       kw_chisq = if (is.null(kw)) NA_real_ else unname(kw$statistic),
       kw_df = if (is.null(kw)) NA_real_ else unname(kw$parameter),
       kw_p = if (is.null(kw)) NA_real_ else kw$p.value,
       eta2 = eta2, omega2 = omega2, epsilon2 = eps2,
       ms_within = ms_w, df_within = df2)
}


# ------------------------------------------------------------------ post-hocs

# All pairwise contrasts for one component. `method` selects both the interval
# and the p-value:
#   tukey         : studentised range, pooled variance   (equal variances)
#   games-howell  : studentised range, Welch df          (unequal variances)
#   everything else: Welch t per pair, p.adjust()ed by that method
#
# Hedges' g accompanies every contrast, because a table of adjusted p-values
# with n = 654 vs 59 tells you almost nothing on its own.
fck_pc_posthoc <- function(scores, g, method = "holm", conf = 0.95,
                           ms_within = NA, df_within = NA) {
  ok <- is.finite(scores) & !is.na(g)
  x <- scores[ok]; g <- droplevels(as.factor(g[ok]))
  lv <- levels(g); k <- length(lv)
  if (k < 2) return(NULL)

  means <- tapply(x, g, mean); vars <- tapply(x, g, stats::var)
  ns <- tapply(x, g, length)
  n <- length(x)
  if (!is.finite(ms_within)) {
    ms_within <- sum((ns - 1) * vars, na.rm = TRUE) / (n - k)
    df_within <- n - k
  }

  rows <- list()
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    a <- lv[i]; b <- lv[j]
    d <- as.numeric(means[a] - means[b])
    n1 <- ns[a]; n2 <- ns[b]; v1 <- vars[a]; v2 <- vars[b]
    if (!is.finite(v1) || !is.finite(v2) || n1 < 2 || n2 < 2) next

    if (identical(method, "tukey")) {
      se <- sqrt(ms_within / 2 * (1 / n1 + 1 / n2))
      q <- abs(d) / se
      p <- stats::ptukey(q, k, df_within, lower.tail = FALSE)
      crit <- stats::qtukey(conf, k, df_within) * se
      stat <- q; dfu <- df_within; label <- "q"
    } else if (identical(method, "games-howell")) {
      se <- sqrt((v1 / n1 + v2 / n2) / 2)
      dfu <- (v1 / n1 + v2 / n2)^2 /
        ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
      q <- abs(d) / se
      p <- stats::ptukey(q, k, dfu, lower.tail = FALSE)
      crit <- stats::qtukey(conf, k, dfu) * se
      stat <- q; label <- "q"
    } else {
      se <- sqrt(v1 / n1 + v2 / n2)
      dfu <- se^4 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
      tstat <- d / se
      p <- 2 * stats::pt(-abs(tstat), dfu)
      crit <- stats::qt(1 - (1 - conf) / 2, dfu) * se
      stat <- tstat; label <- "t"
    }

    sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    J <- 1 - 3 / (4 * (n1 + n2) - 9)
    rows[[length(rows) + 1]] <- data.frame(
      a = a, b = b, diff = d,
      ci_lo = d - crit, ci_hi = d + crit,
      stat = as.numeric(stat), stat_label = label, df = as.numeric(dfu),
      p_raw = as.numeric(p),
      hedges_g = as.numeric(J * d / sp),
      n_a = as.integer(n1), n_b = as.integer(n2),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj <- fck_pc_adjust(out$p_raw, method)
  out$method <- method
  # Tukey / Games-Howell control the family rate by construction, so p_raw for
  # them is already the adjusted quantity. Say so rather than leaving two
  # identical columns unexplained.
  out$already_adjusted <- method %in% c("tukey", "games-howell")
  out
}


# ------------------------------------------------------- eigenvalue separation

# How well is component j separated from its neighbours? The perturbation bound
# for an eigenvector is governed by the RELATIVE eigengap; when that is small
# the eigenfunction is a nearly arbitrary rotation within the near-degenerate
# subspace, and a group difference "on PC j" is not attributable to PC j.
#
# Returns the gap to the nearer neighbour, relative to the eigenvalue itself.
fck_pc_separation <- function(eigenvalues) {
  ev <- eigenvalues[is.finite(eigenvalues)]
  k <- length(ev)
  if (k < 2) return(NULL)
  gap <- numeric(k)
  for (j in seq_len(k)) {
    nb <- c(if (j > 1) ev[j - 1], if (j < k) ev[j + 1])
    gap[j] <- min(abs(ev[j] - nb)) / max(ev[j], .Machine$double.eps)
  }
  data.frame(component = seq_len(k), eigenvalue = ev, rel_gap = gap,
             stability = ifelse(gap >= 0.5, "well separated",
                         ifelse(gap >= 0.2, "moderate",
                                "POORLY separated - do not attribute to this PC alone")),
             stringsAsFactors = FALSE)
}


# --------------------------------------------------------- repeated measures

# The scores are one row per CURVE. If the same participant contributes more
# than one curve, a between-groups ANOVA on those rows treats correlated
# observations as independent and its p-values are anticonservative.
# This detects it so the report can say so; it does not silently "fix" it,
# because the right fix (a mixed model, or averaging within participant)
# is the analyst's choice.
fck_pc_independence_check <- function(subject_ids) {
  if (is.null(subject_ids)) return(list(checked = FALSE))
  sid <- subject_ids[!is.na(subject_ids)]
  n <- length(sid); u <- length(unique(sid))
  list(checked = TRUE, n_rows = n, n_subjects = u,
       n_repeated = n - u,
       ok = n == u,
       max_per_subject = if (n) max(table(sid)) else 0L)
}


# ----------------------------------------------------------- the whole table

# Run the omnibus for every component, then post-hocs for those that pass a
# gate. `across_pc_correction` is the SECOND family: the k omnibus p-values.
fck_pc_anova_all <- function(scores, g, eigenvalues = NULL, varprop = NULL,
                             posthoc_method = "holm",
                             across_pc_correction = "holm",
                             posthoc_gate = 0.05,
                             conf = 0.95,
                             subject_ids = NULL,
                             use_welch = NA) {
  scores <- as.matrix(scores)
  k <- ncol(scores)
  if (k < 1) return(NULL)
  g <- droplevels(as.factor(g))

  omni <- vector("list", k); assum <- vector("list", k); norm <- vector("list", k)
  for (j in seq_len(k)) {
    # NOTE: `l[[j]] <- NULL` DELETES element j and shifts everything after it
    # down. These functions all return NULL for a component that cannot be
    # tested, so assigning with [[<- would silently renumber the components --
    # every result after the first failure would be attributed to the wrong PC.
    # `l[j] <- list(NULL)` stores a NULL instead of removing the slot.
    omni[j]  <- list(fck_pc_omnibus(scores[, j], g, conf))
    assum[j] <- list(fck_brown_forsythe(scores[, j], g))
    norm[j]  <- list(fck_normality_check(scores[, j], g))
  }
  keep <- !vapply(omni, is.null, logical(1))
  if (!any(keep)) return(NULL)

  # which omnibus p to carry forward: Welch when variances differ, unless the
  # caller has overridden it
  p_omni <- vapply(seq_len(k), function(j) {
    o <- omni[[j]]; if (is.null(o)) return(NA_real_)
    hetero <- !is.null(assum[[j]]) && isTRUE(assum[[j]]$p < 0.05)
    usew <- if (is.na(use_welch)) hetero else isTRUE(use_welch)
    if (usew && is.finite(o$welch_p)) o$welch_p else o$p
  }, numeric(1))

  p_omni_adj <- rep(NA_real_, k)
  p_omni_adj[keep] <- fck_pc_adjust(p_omni[keep],
    if (across_pc_correction %in% c("tukey", "games-howell")) "holm" else across_pc_correction)

  ph <- vector("list", k)
  for (j in seq_len(k)) {
    if (is.null(omni[[j]])) next
    # The gate runs on the ACROSS-COMPONENT-ADJUSTED p, never the raw one: a
    # component that only looks significant before that correction must not get
    # a table of pairwise tests lending it credibility.
    if (is.finite(p_omni_adj[j]) && p_omni_adj[j] > posthoc_gate) next
    ph[j] <- list(fck_pc_posthoc(scores[, j], g, posthoc_method, conf,
                                 omni[[j]]$ms_within, omni[[j]]$df_within))
  }

  list(k = k,
       omnibus = omni, variance = assum, normality = norm, posthoc = ph,
       p_omnibus = p_omni, p_omnibus_adj = p_omni_adj,
       across_pc_correction = across_pc_correction,
       posthoc_method = posthoc_method,
       posthoc_gate = posthoc_gate,
       separation = if (!is.null(eigenvalues)) fck_pc_separation(eigenvalues[seq_len(k)]) else NULL,
       varprop = varprop,
       independence = fck_pc_independence_check(subject_ids),
       levels = levels(g))
}


# ---------------------------------------------------- warping parameters

# WHY WARPING PARAMETERS DESERVE THE SAME TEST AS THE SCORES
#
# Registration splits a curve into WHERE it happens (phase) and HOW BIG it is
# (amplitude). The component scores describe the amplitude half, because they
# are computed on the registered curves. The warping functions are the phase
# half, and they are discarded by every analysis that stops at the scores --
# which is exactly the wrong half to throw away when the question is whether
# groups differ in TIMING.
#
# Two parameters, both defined for any warping method:
#   shift           the per-curve time shift, when the method produces one
#   warp_amplitude  RMS deviation of h(t) from the identity, i.e. how far this
#                   curve had to be moved to align it. Method-agnostic, and the
#                   only one available for a nonlinear warp.
#
# They are their OWN multiplicity family, not extra components: "do the groups
# differ in phase" is a different question from "do they differ in the k-th mode
# of amplitude variation", and pooling the two families would let one borrow the
# other's correction.
fck_warp_params <- function(warping_results, time_points = NULL) {
  if (is.null(warping_results)) return(NULL)
  wf <- warping_results$warp_functions
  sh <- warping_results$shifts
  tp <- time_points %||% warping_results$time_points

  amp <- NULL
  if (!is.null(wf) && !is.null(tp) && nrow(wf) == length(tp)) {
    amp <- apply(wf, 2, function(h) sqrt(mean((h - tp)^2, na.rm = TRUE)))
  }
  if (is.null(sh) && is.null(amp)) return(NULL)

  n <- max(length(sh), length(amp))
  # Build the columns first: assigning into a zero-row data.frame throws
  # "replacement has n rows, data has 0" rather than growing it.
  cols <- list()
  if (!is.null(sh) && length(sh) == n) cols$shift <- as.numeric(sh)
  if (!is.null(amp) && length(amp) == n) cols$warp_amplitude <- as.numeric(amp)
  if (!length(cols)) return(NULL)
  out <- as.data.frame(cols, stringsAsFactors = FALSE)

  # A warp that never moved anything is not a parameter worth testing: it is a
  # constant, and an ANOVA on a constant is undefined rather than null.
  keep <- vapply(out, function(v) {
    v <- v[is.finite(v)]
    length(v) > 1 && stats::sd(v) > .Machine$double.eps^0.5
  }, logical(1))
  if (!any(keep)) return(NULL)
  out[, keep, drop = FALSE]
}

# Human labels for the warping parameters, so the report does not print a
# column name at a reader.
FCK_WARP_LABELS <- c(
  shift = "Time shift (how far this curve was moved)",
  warp_amplitude = "Warp amplitude (RMS distance of h(t) from the identity)"
)

# ==============================================================================
# FUNCTIONAL L2 NORM (P1.3)
# ==============================================================================
# The fANOVA modules computed their global statistic as sqrt(sum(v^2)) over the
# evaluation grid. That is a vector norm, not a functional one: its value scales
# with how densely the grid happens to be sampled, so the same data evaluated on
# 50 versus 100 points gives different numbers, and on an unevenly spaced grid it
# silently weights the dense regions more.
#
# The L2 norm of a function is sqrt(integral v(t)^2 dt). Trapezoidal weights are
# exact for the piecewise-linear interpolant the app already draws, and reduce to
# the old constant-weight form (up to the grid spacing) when the grid is even.
fck_l2_norm <- function(v, argvals = NULL) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  if (!any(ok)) return(NA_real_)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = length(v))
  argvals <- as.numeric(argvals)
  if (length(argvals) != length(v)) return(NA_real_)
  v <- v[ok]; argvals <- argvals[ok]
  if (length(v) < 2) return(abs(v[1]))
  o <- order(argvals); v <- v[o]; argvals <- argvals[o]
  d <- diff(argvals)
  w <- c(d[1], (head(d, -1) + tail(d, -1)), d[length(d)]) / 2   # trapezoid
  sqrt(sum(w * v^2))
}
