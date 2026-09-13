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
# CALIBRATION STATUS OF THE OMNIBUS: MEASURED, ON A STATED GRID
# ==============================================================================
# READ THIS BEFORE QUOTING A P-VALUE FROM dance_traj_omnibus().
#
# This block test once rejected 27.5% of true nulls at a nominal .05, and it
# shipped saying so. The cause turned out to be two things, both in how the
# random-effects ladder was walked, and both now fixed.
#
# 1. THE LADDER HAD NO CURVE LEVEL. Its top rung was
#    (1 + basis + within | subject), which gives each participant their own
#    rhythm and their own OFFSET between conditions, but forces one amplitude
#    and one acrophase across that participant's conditions. Real
#    repeated-measures data varies at participant x condition -- at the level of
#    the observed curve. That variance had nowhere to go but the residual, the
#    residual is shared across the fit, and the standard errors of the basis x
#    design terms under test came out too small.
#
# 2. SINGULARITY WAS DESCENDING THE LADDER, TWICE OVER. The walk dropped a rung
#    when lme4::isSingular flagged it -- and then, after that test was removed,
#    kept dropping it anyway, because lme4 files "boundary (singular) fit" in
#    the same message list it uses for convergence failures. Both paths dropped
#    the curve level in 64-100% of null fits, which is to say the rung added in
#    (1) was almost never the rung used. A boundary variance estimate is not a
#    failed fit; dropping the term is what does the damage.
#
# WHAT WAS MEASURED AFTER THE FIX. 12 time points over 22 h, one harmonic, no
# trend, Gaussian noise, nominal .05, rejection rates:
#
#   NULL   mixed 2 x 2, 16 participants, curve-level variance only   0.040  (n=100)
#          mixed 2 x 2, 16 participants, participant + curve         0.020  (n=100)
#          between only, 16 participants, one curve each             0.050  (n=60)
#          between only, 32 participants, one curve each             0.017  (n=60)
#   POWER  phase shift of 1.2 rad in one cell                        1.000  (n=60)
#          phase shift of 0.5 rad in one cell                        0.850  (n=60)
#          phase shift of 0.3 rad in one cell                        0.333  (n=60)
#
# So the test holds its nominal size, and at 16 participants with both variance
# components present it is somewhat CONSERVATIVE (0.020) rather than merely
# nominal -- the usual price of Kenward-Roger on a fit with a component at the
# boundary. Conservative costs power; it does not cost validity. The power row
# is printed alongside because a size assertion on its own can be satisfied by a
# test that never rejects anything.
#
# WHAT WAS NOT MEASURED, and what this calibration therefore does NOT cover:
# more than one harmonic; any trend term (linear, log, or the profiled
# saturating trend, where tau is estimated from the same data and the reported
# df do not account for it); unbalanced or missing cells; designs with three or
# more factors; non-Gaussian or serially correlated residuals; the glmmTMB
# engine, where Kenward-Roger is unavailable and the block test falls back to an
# asymptotic Wald statistic. Outside that grid the p-value is a reasonable
# default, not a validated one, and dance_traj_omnibus() says which case it is
# in through $validated and $calibration.
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
# The calibration statement every omnibus result carries. Two of them, because
# the simulation grid covers one harmonic with no trend and nothing else, and a
# result from outside that grid must not borrow its warrant.
DANCE_TRAJ_OMNIBUS_CALIBRATED <- paste(
  "CALIBRATED on the grid this fit sits in. Under a true null (12 time points,",
  "one harmonic, no trend, Gaussian noise, nominal .05) this block test rejected",
  "at 0.040 and 0.020 in a 2 x 2 mixed design with 16 participants (n = 100 each,",
  "curve-level variance only / participant and curve variance) and at 0.050 and",
  "0.017 in a between-participant design with 16 and 32 participants (n = 60).",
  "Power against a phase shift of 1.2 / 0.5 / 0.3 rad in one cell was 1.000 /",
  "0.850 / 0.333. At 16 participants with both variance components present the",
  "test is conservative rather than exact, which costs power and not validity.")

DANCE_TRAJ_OMNIBUS_UNCALIBRATED <- paste(
  "OUTSIDE THE CALIBRATED GRID. The type-I error of this block test was measured",
  "only for one harmonic with no trend term, on balanced designs with Gaussian",
  "residuals and the lmer engine. This fit is not one of those: %s. The p-value",
  "is a reasonable default, not a validated one. With a profiled saturating",
  "trend in particular, tau is estimated from the same data and the reported",
  "degrees of freedom do not account for it, so the p-value is optimistic by an",
  "unquantified amount.")

# Which of the two applies, and why. `method_kind` matters as much as the model
# does: the grid was measured with Kenward-Roger throughout, so the Satterthwaite,
# Wald and likelihood-ratio branches are outside it by construction even when the
# model itself is an ordinary one-harmonic fit.
dance_traj_calibration <- function(fit, method_kind = "kr") {
  spec <- fit$spec
  out <- character(0)
  if (!identical(method_kind, "kr"))
    out <- c(out, switch(method_kind,
      satterthwaite = "the sample is past the Kenward-Roger cutoff, so the denominator degrees of freedom come from Satterthwaite, which the grid did not cover",
      wald          = "the test is an asymptotic Wald/likelihood-ratio statistic, not a Kenward-Roger F",
      lrt           = "the F approximations were unavailable and the test fell back to an asymptotic likelihood-ratio statistic",
      sprintf("the test used the '%s' method rather than Kenward-Roger", method_kind)))
  if ((spec$n_harmonics %||% 1L) > 1L)
    out <- c(out, sprintf("it fits %d harmonics", spec$n_harmonics))
  if (!identical(spec$trend %||% "none", "none"))
    out <- c(out, sprintf("it carries a '%s' trend term", spec$trend))
  if (!identical(fit$engine %||% "lmer", "lmer"))
    out <- c(out, sprintf("it uses the %s engine, where Kenward-Roger is unavailable",
                          fit$engine))
  if (isTRUE(fit$ar1)) out <- c(out, "it models a residual AR(1) structure")
  if (length(unique(spec$cells$n_obs)) > 1L)
    out <- c(out, "its design cells are unbalanced")
  if (!length(out))
    return(list(validated = TRUE, calibration = DANCE_TRAJ_OMNIBUS_CALIBRATED))
  list(validated = FALSE,
       calibration = sprintf(DANCE_TRAJ_OMNIBUS_UNCALIBRATED,
                             paste(out, collapse = "; ")))
}

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
                p = an$`Pr(>Chisq)`[2],
                validated = dance_traj_calibration(fit, "wald")$validated,
                calibration = dance_traj_calibration(fit, "wald")$calibration,
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
      cal <- dance_traj_calibration(fit, "kr")
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = sprintf("Kenward-Roger F (%d participants <= %d)",
                                   n_subj, kr_max_subjects),
                  statistic = unname(st["Ftest", "stat"]),
                  df1 = unname(st["Ftest", "ndf"]), df2 = unname(st["Ftest", "ddf"]),
                  p = unname(st["Ftest", "p.value"]), caveat = NULL,
                  validated = cal$validated, calibration = cal$calibration))
    }
  }

  # Satterthwaite via lmerTest, over the whole block at once
  if (inherits(m, "lmerModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    X <- stats::model.matrix(m)
    keep <- attr(X, "assign") %in% which(attr(stats::terms(m), "term.labels") %in% terms_in)
    L <- diag(ncol(X))[keep, , drop = FALSE]
    ct <- tryCatch(lmerTest::contest(m, L, joint = TRUE), error = function(e) NULL)
    cals <- dance_traj_calibration(fit, "satterthwaite")
    if (!is.null(ct))
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = sprintf("Satterthwaite F (%d participants > %d)",
                                   n_subj, kr_max_subjects),
                  statistic = ct[["F value"]], df1 = ct[["NumDF"]],
                  df2 = ct[["DenDF"]], p = ct[["Pr(>F)"]], caveat = NULL,
                  validated = cals$validated, calibration = cals$calibration))
  }

  # last resort: a likelihood-ratio test on ML refits
  mm <- refit_ml(m)
  m0 <- tryCatch(stats::update(mm, reduced_formula), error = function(e) NULL)
  if (is.null(m0)) return(list(ok = FALSE, which = which,
                               message = "The reduced model could not be fitted."))
  an <- tryCatch(stats::anova(m0, mm), error = function(e) NULL)
  if (is.null(an)) return(list(ok = FALSE, which = which, message = "The block comparison failed."))
  call <- dance_traj_calibration(fit, "lrt")
  list(ok = TRUE, which = which, terms = terms_in,
       method = "likelihood-ratio on ML refits",
       statistic = an$Chisq[2], df1 = an$Df[2], df2 = NA_real_, p = an$`Pr(>Chisq)`[2],
       validated = call$validated, calibration = call$calibration,
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
  jc <- dance_traj_joint_cov(fit, ck, sk, specs)
  list(ok = TRUE, harmonic = h, period = spec$period,
       effective_period = spec$period / h,
       cells = cells, grid = A$grid,
       a = A$est, b = B$est, Va = A$V, Vb = B$V,
       # cov(a_i, b_i) is not returned by emmeans across two emtrends calls;
       # it is recovered from the linear maps below
       joint = jc$blocks, linfct = jc$linfct, vcov_beta = jc$vcov_beta)
}

# The covariance of the cell coefficients, WITHIN and BETWEEN cells.
# ------------------------------------------------------------------------------
# Every cell coefficient is a linear function of the fixed effects: a_i = l_i' b
# for a row l_i of the emmeans linear map. So the covariance of ANY set of them
# is L V L' exactly, for the corresponding rows -- within a cell and across
# cells alike. The first version of this function returned only the within-cell
# 2 x 2 blocks and left every contrast to assume the cross-cell block was zero.
# That is not a conservative simplification. For a difference,
#
#     Var(t2 - t1) = V11 + V22 - 2 C12
#
# so dropping C12 overstates the variance when C12 > 0 and UNDERSTATES it when
# C12 < 0 -- and C12 is routinely negative between cells of a factorial design,
# because they are built from the same interaction coefficients with opposite
# signs. An interval that is too narrow is the one failure mode a contrast must
# not have, so the full map is kept and every contrast uses the exact block.
dance_traj_joint_cov <- function(fit, ck, sk, specs) {
  Vb <- tryCatch(as.matrix(stats::vcov(fit$model)), error = function(e) NULL)
  if (is.null(Vb)) return(list(blocks = NULL, linfct = NULL, vcov_beta = NULL))
  L <- function(term) {
    e <- tryCatch(emmeans::emtrends(fit$model, specs = specs, var = term),
                  error = function(e) NULL)
    if (is.null(e)) return(NULL)
    e@linfct
  }
  La <- L(ck); Lb <- L(sk)
  if (is.null(La) || is.null(Lb) || ncol(La) != ncol(Vb))
    return(list(blocks = NULL, linfct = NULL, vcov_beta = NULL))
  n <- nrow(La)
  blocks <- lapply(seq_len(n), function(i) {
    M <- rbind(La[i, ], Lb[i, ])
    M %*% Vb %*% t(M)
  })
  list(blocks = blocks, linfct = list(a = La, b = Lb), vcov_beta = Vb)
}

# The exact 4 x 4 covariance of (a_i, b_i, a_j, b_j), cross-cell block included.
# Falls back to the block-diagonal form, and SAYS it has, only when the linear
# maps are unavailable.
dance_traj_pair_cov <- function(coefs, i, j) {
  lf <- coefs$linfct; Vb <- coefs$vcov_beta
  if (!is.null(lf) && !is.null(Vb)) {
    M <- rbind(lf$a[i, ], lf$b[i, ], lf$a[j, ], lf$b[j, ])
    return(list(V = M %*% Vb %*% t(M), exact = TRUE))
  }
  if (is.null(coefs$joint)) return(NULL)
  V <- matrix(0, 4, 4)
  V[1:2, 1:2] <- coefs$joint[[i]]
  V[3:4, 3:4] <- coefs$joint[[j]]
  list(V = V, exact = FALSE)
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
  if (!isTRUE(coefs$ok)) return(NULL)
  pc <- dance_traj_pair_cov(coefs, i, j)
  if (is.null(pc)) return(NULL)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  k <- coefs$effective_period / (2 * pi)
  a1 <- coefs$a[i]; b1 <- coefs$b[i]; a2 <- coefs$a[j]; b2 <- coefs$b[j]
  A1 <- sqrt(a1^2 + b1^2); A2 <- sqrt(a2^2 + b2^2)
  d <- atan2(b2, a2) - atan2(b1, a1)
  d <- ((d + pi) %% (2 * pi)) - pi                 # wrap to (-pi, pi]
  g1 <- c(b1 / A1^2, -a1 / A1^2)                   # d(-phi1)
  g2 <- c(-b2 / A2^2, a2 / A2^2)                   # d(phi2)
  g <- c(g1, g2)
  se <- sqrt(max(0, as.numeric(t(g) %*% pc$V %*% g)))
  # A contrast can be well determined when neither phase alone is, but not when
  # either amplitude may be zero: the map to an angle is undefined there, so the
  # delta-method SE is meaningless rather than merely wide. Bingham's rule
  # applies to the contrast through its two endpoints.
  sA <- function(a, b, V2) {
    A <- sqrt(a^2 + b^2); if (!(A > 0)) return(NA_real_)
    gA <- c(a / A, b / A); sqrt(max(0, as.numeric(t(gA) %*% V2 %*% gA)))
  }
  seA1 <- sA(a1, b1, pc$V[1:2, 1:2]); seA2 <- sA(a2, b2, pc$V[3:4, 3:4])
  defined <- is.finite(seA1) && is.finite(seA2) &&
             (A1 - z * seA1) > 0 && (A2 - z * seA2) > 0
  list(cell1 = coefs$cells[i], cell2 = coefs$cells[j],
       diff_time = d * k, se_time = if (defined) se * k else NA_real_,
       lo = if (defined) (d - z * se) * k else NA_real_,
       hi = if (defined) (d + z * se) * k else NA_real_,
       effective_period = coefs$effective_period,
       amp1 = A1, amp2 = A2, defined = defined,
       exact_cov = isTRUE(pc$exact),
       note = if (!defined)
         paste("At least one cell's amplitude interval covers zero, so its acrophase",
               "is undefined and so is this contrast (Bingham et al., 1982). The",
               "wrapped point difference is reported without an interval.")
       else if (!isTRUE(pc$exact))
         paste("The emmeans linear map was unavailable, so the cross-cell covariance",
               "was treated as zero. That is NOT conservative in general -- it",
               "understates the variance whenever the two cells covary negatively,",
               "which is common in a factorial design. Treat this interval as",
               "indicative only.")
       else NULL)
}
