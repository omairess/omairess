# ==============================================================================
# server/08f_helpers_trajinf.R — LAYERS C and D: inference and derived parameters
# ==============================================================================
# P21 phase 2/3.
#
# LAYER C -- THE OMNIBUS IS A BLOCK TEST (brief §5)
# ------------------------------------------------
# The module this replaces asked "do the amplitudes differ?" and "do the
# acrophases differ?" as separate marginal tests on separate per-participant
# estimates. Amplitude and acrophase are two coordinates of ONE bivariate
# object, and a rhythm can differ in both while neither marginal reaches
# significance. The primary question is therefore a joint test of a coefficient
# BLOCK:
#
#   A. circadian trajectory difference   H0: every Delta a_k = Delta b_k = 0
#   B. overall trajectory difference     H0: trend block and harmonic block = 0
#
# Both are the same machinery over a different set of terms.
#
# WHICH TEST. Kenward-Roger is the better small-sample approximation and is too
# slow to be unconditional, so the rule is: KR below `kr_max_subjects`
# participants, Satterthwaite above, and the output SAYS WHICH RAN. On glmmTMB
# (the residual-AR(1) engine) neither is available and the test is an asymptotic
# Wald chi-square, which is also said rather than quietly substituted.
#
# LAYER D -- DERIVED PARAMETERS COME FROM THE SAME FIT (brief §6)
# --------------------------------------------------------------
# Amplitude and acrophase for a design cell are nonlinear functions of that
# cell's (a_k, b_k) pair. The pair is a LINEAR function of the fixed effects, so
# emmeans gives it together with an exact covariance, and the delta method turns
# that into intervals for A = sqrt(a^2 + b^2) and phi = atan2(b, a):
#
#   dA/da = a/A,  dA/db = b/A          Var(A)   = g' V g
#   dphi/da = -b/A^2,  dphi/db = a/A^2  Var(phi) = h' V h
#
# ==============================================================================
# CALIBRATION STATUS OF THE OMNIBUS: NOT YET ESTABLISHED
# ==============================================================================
# READ THIS BEFORE QUOTING A P-VALUE FROM dance_traj_omnibus().
#
# tests/traj_framework_test.R measured the type-I error of the circadian block
# under a true null -- 40 simulated datasets, 16 participants (8 per group), two
# within-participant conditions, 12 time points, no group or condition effect on
# the rhythm -- and it rejected 11 of 40 times at a nominal .05, a rate of 0.275.
# That is not a calibrated test and the p-values are not reportable as they
# stand.
#
# WHAT IS ESTABLISHED, AND WHAT IS NOT. The first hypothesis -- that the
# random-effects ladder was falling back to a random intercept and testing a
# between-participant effect against a within-participant residual -- is
# REFUTED: on a representative null dataset the ladder reached rung 1, its most
# complex structure, and the denominator degrees of freedom were 69.9 there
# against 70.3 for a forced (1 + c1 + s1 | subject). So it is not a df problem.
#
# What the same dataset does show is that the random structure is doing a great
# deal of inferential work: at rung 1, which includes a participant-specific
# Condition effect, the circadian block gives F = 2.712 on (6, 69.9), p = .020;
# with (1 + c1 + s1 | subject) forced on the identical data it gives F = 1.156
# on (6, 70.3), p = .340. Rung 1 also shows a -0.909 correlation between the
# participant intercept and the participant Condition effect, which is the shape
# of a random-effects covariance the data cannot separate -- and lme4::isSingular
# does not flag it, because no variance is exactly zero. Brief §4 lists
# "non-identifiable random-effect covariance structures" as something to check,
# and dance_traj_is_singular() does not currently check it.
#
# That is a hypothesis, not a finding. Which structure is right, and what the
# type-I error is once the structure is chosen properly, is the work of the
# validation gate between phases 2 and 3 of the migration plan -- the §20
# simulation grid. Until that gate clears, every result from this function
# carries validated = FALSE and a caveat naming the measured rate, so it cannot
# be read as a finished p-value by someone who did not read this comment.
#
# Layers A, B and D are unaffected: the specification, the fitting log and the
# derived amplitude/acrophase estimates with their delta-method intervals all
# pass their own checks, and amplitude recovery and Bingham's undefined-phase
# rule are verified. It is the OMNIBUS TEST specifically that is not yet
# trustworthy.
# ==============================================================================
#
# BINGHAM'S RULE IS ENFORCED, NOT RECITED. When the amplitude interval covers
# zero the acrophase interval is undefined -- the phase of a vector that may be
# the zero vector is not a quantity -- and this file returns NA with a reason
# instead of a number. That is the constraint the classical literature states
# and that most implementations quietly ignore.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# Which terms form the block for each omnibus question
# ------------------------------------------------------------------------------
# A term belongs to the CIRCADIAN block if it involves a harmonic column AND at
# least one design factor; to the TRAJECTORY block if it involves any basis
# column (trend or harmonic) and a design factor. Main effects of the design
# alone are a level difference, not a trajectory difference, and are excluded
# from both -- that is what makes B "does the SHAPE differ" rather than "does
# anything differ".
dance_traj_block_terms <- function(fit, which = c("circadian", "trajectory", "level")) {
  which <- match.arg(which)
  spec <- fit$spec
  tl <- attr(stats::terms(stats::as.formula(spec$fixed_formula)), "term.labels")
  has <- function(term, set) any(strsplit(term, ":", fixed = TRUE)[[1]] %in% set)
  keep <- vapply(tl, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    with_design <- any(parts %in% spec$design_terms)
    switch(which,
      circadian  = has(term, spec$harm_terms) && with_design,
      trajectory = has(term, spec$basis_terms) && with_design,
      # "level" is the design's effect on the fitted value AT t = 0, not on a
      # MESOR: with the basis centred at t0 the design main effects are exactly
      # that, and brief §3 requires the distinction to be stated rather than
      # assumed. A rhythm-adjusted mean is a separate derived quantity.
      level      = with_design && !has(term, spec$basis_terms))
  }, logical(1))
  tl[keep]
}

# ------------------------------------------------------------------------------
# THE BLOCK TEST
# ------------------------------------------------------------------------------
# The caveat every omnibus result carries until the validation gate clears.
DANCE_TRAJ_OMNIBUS_CAVEAT <- paste(
  "NOT YET VALIDATED. Under a true null (40 simulated datasets, 16 participants,",
  "2 within-participant conditions, 12 time points, no effect on the rhythm) this",
  "block test rejected at 0.275 against a nominal .05. The cause is under",
  "investigation -- it is NOT a degrees-of-freedom problem, and the leading",
  "hypothesis is that the random-effects ladder accepts a covariance structure",
  "the data cannot separate. Read the statistic and the fitted trajectories;",
  "do not report this p-value.")

dance_traj_omnibus <- function(fit, which = c("circadian", "trajectory", "level"),
                               kr_max_subjects = 50) {
  which <- match.arg(which)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  terms_in <- dance_traj_block_terms(fit, which)
  if (!length(terms_in))
    return(list(ok = FALSE, which = which,
                message = sprintf("No %s x design terms in this model -- nothing to test.", which)))

  m <- fit$model
  reduced_formula <- stats::as.formula(
    paste(". ~ . -", paste(terms_in, collapse = " - ")))

  # glmmTMB: asymptotic Wald only, and labelled as such
  if (identical(fit$engine, "glmmTMB")) {
    m0 <- tryCatch(stats::update(m, reduced_formula), error = function(e) NULL)
    if (is.null(m0)) return(list(ok = FALSE, which = which,
                                 message = "The reduced model could not be fitted."))
    an <- tryCatch(stats::anova(m0, m), error = function(e) NULL)
    if (is.null(an)) return(list(ok = FALSE, which = which,
                                 message = "The block comparison failed."))
    return(list(ok = TRUE, which = which, terms = terms_in,
                method = "asymptotic likelihood-ratio (glmmTMB)",
                statistic = an$Chisq[2], df1 = an$`Chi Df`[2], df2 = NA_real_,
                p = an$`Pr(>Chisq)`[2], validated = FALSE,
                calibration = DANCE_TRAJ_OMNIBUS_CAVEAT,
                caveat = paste("Kenward-Roger is unavailable on this engine; this test is",
                               "asymptotic and is optimistic in small samples.")))
  }

  n_subj <- fit$spec$n_participants
  use_kr <- is.finite(n_subj) && n_subj <= kr_max_subjects &&
            requireNamespace("pbkrtest", quietly = TRUE)

  # both KR and the LRT need ML fits when the FIXED effects differ
  refit_ml <- function(mod) if (isTRUE(fit$REML))
    tryCatch(stats::update(mod, REML = FALSE), error = function(e) mod) else mod

  if (use_kr) {
    # KR compares REML fits, which is correct for fixed-effect contrasts
    mL <- if (inherits(m, "lmerModLmerTest")) as(m, "lmerMod") else m
    m0 <- tryCatch(stats::update(mL, reduced_formula), error = function(e) NULL)
    kr <- if (is.null(m0)) NULL else
      tryCatch(pbkrtest::KRmodcomp(mL, m0), error = function(e) NULL)
    if (!is.null(kr)) {
      st <- kr$test
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = sprintf("Kenward-Roger F (%d participants <= %d)",
                                   n_subj, kr_max_subjects),
                  statistic = unname(st["Ftest", "stat"]),
                  df1 = unname(st["Ftest", "ndf"]), df2 = unname(st["Ftest", "ddf"]),
                  p = unname(st["Ftest", "p.value"]), caveat = NULL,
                  validated = FALSE, calibration = DANCE_TRAJ_OMNIBUS_CAVEAT))
    }
  }

  # Satterthwaite via lmerTest, over the whole block at once
  if (inherits(m, "lmerModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    X <- stats::model.matrix(m)
    keep <- attr(X, "assign") %in% which(attr(stats::terms(m), "term.labels") %in% terms_in)
    L <- diag(ncol(X))[keep, , drop = FALSE]
    ct <- tryCatch(lmerTest::contest(m, L, joint = TRUE), error = function(e) NULL)
    if (!is.null(ct))
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = sprintf("Satterthwaite F (%d participants > %d)",
                                   n_subj, kr_max_subjects),
                  statistic = ct[["F value"]], df1 = ct[["NumDF"]],
                  df2 = ct[["DenDF"]], p = ct[["Pr(>F)"]], caveat = NULL,
                  validated = FALSE, calibration = DANCE_TRAJ_OMNIBUS_CAVEAT))
  }

  # last resort: a likelihood-ratio test on ML refits
  mm <- refit_ml(m)
  m0 <- tryCatch(stats::update(mm, reduced_formula), error = function(e) NULL)
  if (is.null(m0)) return(list(ok = FALSE, which = which,
                               message = "The reduced model could not be fitted."))
  an <- tryCatch(stats::anova(m0, mm), error = function(e) NULL)
  if (is.null(an)) return(list(ok = FALSE, which = which, message = "The block comparison failed."))
  list(ok = TRUE, which = which, terms = terms_in,
       method = "likelihood-ratio on ML refits",
       statistic = an$Chisq[2], df1 = an$Df[2], df2 = NA_real_, p = an$`Pr(>Chisq)`[2],
       validated = FALSE, calibration = DANCE_TRAJ_OMNIBUS_CAVEAT,
       caveat = "Asymptotic; prefer the F tests above where they are available.")
}

# ------------------------------------------------------------------------------
# LAYER D: the (cos, sin) pair per design cell, with its covariance
# ------------------------------------------------------------------------------
# emmeans is asked for the cell means of each basis column. Each is a linear
# function of the fixed effects, so the covariance across columns within a cell
# is exact and the delta method below is not an approximation of the sampling
# scheme -- only of the nonlinear map.
dance_traj_cell_coefs <- function(fit, harmonic = 1) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  if (!requireNamespace("emmeans", quietly = TRUE))
    return(list(ok = FALSE, message = "emmeans is required for cell-level estimates."))
  spec <- fit$spec
  h <- max(1L, as.integer(harmonic))
  ck <- paste0("c", h); sk <- paste0("s", h)
  if (!all(c(ck, sk) %in% spec$harm_terms))
    return(list(ok = FALSE, message = sprintf("Harmonic %d was not fitted.", h)))
  if (!length(spec$design_terms))
    return(list(ok = FALSE, message = "No design factors: there are no cells to compare."))

  specs <- stats::as.formula(paste("~", paste(spec$design_terms, collapse = " * ")))
  grab <- function(term) {
    e <- tryCatch(emmeans::emtrends(fit$model, specs = specs, var = term),
                  error = function(e) NULL)
    if (is.null(e)) return(NULL)
    list(est = summary(e)[[paste0(term, ".trend")]],
         V = stats::vcov(e), grid = as.data.frame(e)[, spec$design_terms, drop = FALSE])
  }
  A <- grab(ck); B <- grab(sk)
  if (is.null(A) || is.null(B))
    return(list(ok = FALSE, message = "emmeans could not extract the harmonic coefficients."))

  cells <- do.call(paste, c(A$grid, list(sep = " x ")))
  list(ok = TRUE, harmonic = h, period = spec$period,
       effective_period = spec$period / h,
       cells = cells, grid = A$grid,
       a = A$est, b = B$est, Va = A$V, Vb = B$V,
       # cov(a_i, b_i) is not returned by emmeans across two emtrends calls;
       # it is recovered from the joint contrast below
       joint = dance_traj_joint_cov(fit, ck, sk, specs, spec$design_terms))
}

# The 2x2 covariance of (a_k, b_k) WITHIN each cell, which the delta method
# needs and which two separate emtrends calls cannot supply.
dance_traj_joint_cov <- function(fit, ck, sk, specs, design_terms) {
  Vb <- as.matrix(stats::vcov(fit$model))
  bn <- rownames(Vb)
  # the linear map from fixed effects to each cell's coefficient
  L <- function(term) {
    e <- tryCatch(emmeans::emtrends(fit$model, specs = specs, var = term),
                  error = function(e) NULL)
    if (is.null(e)) return(NULL)
    lin <- e@linfct
    colnames(lin) <- bn[seq_len(ncol(lin))]
    lin
  }
  La <- L(ck); Lb <- L(sk)
  if (is.null(La) || is.null(Lb)) return(NULL)
  n <- nrow(La)
  lapply(seq_len(n), function(i) {
    M <- rbind(La[i, ], Lb[i, ])
    M %*% Vb %*% t(M)
  })
}

# ------------------------------------------------------------------------------
# DELTA-METHOD amplitude and acrophase, with Bingham's rule enforced
# ------------------------------------------------------------------------------
dance_traj_amp_phase <- function(coefs, conf = 0.95) {
  if (!isTRUE(coefs$ok)) return(coefs)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  n <- length(coefs$a)
  k <- coefs$effective_period / (2 * pi)      # radians -> time on P/h
  out <- do.call(rbind, lapply(seq_len(n), function(i) {
    a <- coefs$a[i]; b <- coefs$b[i]
    V <- if (!is.null(coefs$joint)) coefs$joint[[i]] else NULL
    A <- sqrt(a^2 + b^2)
    phi <- atan2(b, a) %% (2 * pi)
    se_A <- NA_real_; se_phi <- NA_real_
    if (!is.null(V) && A > 0) {
      gA <- c(a / A, b / A)
      gP <- c(-b / A^2, a / A^2)
      se_A <- sqrt(max(0, as.numeric(t(gA) %*% V %*% gA)))
      se_phi <- sqrt(max(0, as.numeric(t(gP) %*% V %*% gP)))
    }
    A_lo <- A - z * se_A; A_hi <- A + z * se_A
    # BINGHAM: the acrophase interval is undefined when the amplitude interval
    # covers zero. Return NA and a reason, not a number.
    amp_covers_zero <- !is.finite(A_lo) || A_lo <= 0
    data.frame(
      cell = coefs$cells[i], beta_cos = a, beta_sin = b,
      amplitude = A, amplitude_se = se_A,
      amplitude_lo = A_lo, amplitude_hi = A_hi,
      acrophase_rad = phi,
      acrophase_time = phi * k,
      acrophase_se_time = se_phi * k,
      acrophase_lo = if (amp_covers_zero) NA_real_ else (phi - z * se_phi) * k,
      acrophase_hi = if (amp_covers_zero) NA_real_ else (phi + z * se_phi) * k,
      phase_defined = !amp_covers_zero,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  list(ok = TRUE, harmonic = coefs$harmonic,
       effective_period = coefs$effective_period, table = out,
       method = sprintf("delta method on the fitted (cos, sin) pair, %.0f%% intervals",
                        100 * conf),
       undefined_phase = sum(!out$phase_defined),
       note = if (any(!out$phase_defined)) paste(
         "One or more cells have an amplitude interval covering zero. Their",
         "acrophase interval is UNDEFINED (Bingham et al., 1982): the phase of a",
         "vector that may be the zero vector is not a quantity. Those rows report",
         "the point estimate and no interval.") else NULL)
}

# ------------------------------------------------------------------------------
# PHASE CONTRASTS between two cells
# ------------------------------------------------------------------------------
# The wrapped difference of two angles on the harmonic's effective period, with
# its own delta-method SE from the joint covariance of the FOUR coefficients.
# A contrast can be well determined even when neither phase alone is, provided
# both amplitudes are bounded away from zero -- so the guard is on both
# amplitudes, not on the contrast.
dance_traj_phase_contrast <- function(coefs, i, j, conf = 0.95) {
  if (!isTRUE(coefs$ok) || is.null(coefs$joint)) return(NULL)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  k <- coefs$effective_period / (2 * pi)
  a1 <- coefs$a[i]; b1 <- coefs$b[i]; a2 <- coefs$a[j]; b2 <- coefs$b[j]
  A1 <- sqrt(a1^2 + b1^2); A2 <- sqrt(a2^2 + b2^2)
  d <- atan2(b2, a2) - atan2(b1, a1)
  d <- ((d + pi) %% (2 * pi)) - pi                 # wrap to (-pi, pi]
  g1 <- c(b1 / A1^2, -a1 / A1^2)                   # d(-phi1)
  g2 <- c(-b2 / A2^2, a2 / A2^2)                   # d(phi2)
  V <- matrix(0, 4, 4)
  V[1:2, 1:2] <- coefs$joint[[i]]
  V[3:4, 3:4] <- coefs$joint[[j]]
  g <- c(g1, g2)
  se <- sqrt(max(0, as.numeric(t(g) %*% V %*% g)))
  list(cell1 = coefs$cells[i], cell2 = coefs$cells[j],
       diff_time = d * k, se_time = se * k,
       lo = (d - z * se) * k, hi = (d + z * se) * k,
       effective_period = coefs$effective_period,
       amp1 = A1, amp2 = A2,
       # the cross-cell covariance is dropped here (block-diagonal V), which is
       # CONSERVATIVE when the two cells share fixed effects: it cannot make the
       # interval too narrow. Said rather than assumed.
       note = "cross-cell covariance treated as zero, which widens the interval rather than narrowing it")
}
