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
# WHICH TEST, AND WHY THERE IS NO MAGIC NUMBER IN IT. An earlier version of this
# file chose Kenward-Roger below 50 participants and Satterthwaite at or above
# it, and the method string it printed said "(12 participants <= 50)" as though
# 50 were a statistical fact. It is not. There is no theorem that makes KR
# required below some n and Satterthwaite correct above it; KR is the better
# small-sample approximation at every n, and the only thing that changes with n
# is how long it takes to compute.
#
# So the rule is now stated as what it is: PREFER Kenward-Roger whenever the
# engine supports it and it succeeds. Satterthwaite is a cheaper alternative,
# used when KR is unavailable, fails, or is declined -- and the output always
# names which one ran. A caller that needs to cap runtime may pass
# kr_max_subjects, and that is labelled a COMPUTATIONAL POLICY in the method
# string, not a statistical threshold; it is off by default. On glmmTMB (the
# residual-correlation engine) neither is available and the test is an
# asymptotic Wald chi-square, which is also said rather than quietly
# substituted.
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
DANCE_TRAJ_BLOCKS <- c("full", "shape", "circadian", "level", "trend")

# The English of each, so the screen, the report and the exported script cannot
# describe the same test three different ways.
DANCE_TRAJ_BLOCK_LABEL <- c(
  full      = "Full trajectory difference",
  shape     = "Temporal-shape difference",
  circadian = "Circadian (rhythmic) difference",
  level     = "Level difference at the reference time",
  trend     = "Non-circadian trend difference")

DANCE_TRAJ_BLOCK_H0 <- c(
  full      = "H0: no difference in level, in any trend term, or in any harmonic coefficient -- the trajectories coincide.",
  shape     = "H0: no difference in any trend term or harmonic coefficient. A pure vertical shift is EXCLUDED, so this asks whether the curves are the same shape at possibly different heights.",
  circadian = "H0: every difference in a cos/sin coefficient is zero -- same rhythm, whatever the level or the trend does.",
  level     = "H0: no difference in the fitted value AT THE REFERENCE TIME t0. This is the intercept, NOT the MESOR: with a trend in the model the rhythm-adjusted mean over the window is a different quantity.",
  trend     = "H0: no difference in any non-circadian trend coefficient.")

dance_traj_block_terms <- function(fit, which = DANCE_TRAJ_BLOCKS) {
  which <- match.arg(which, DANCE_TRAJ_BLOCKS)
  spec <- fit$spec
  tl <- attr(stats::terms(stats::as.formula(spec$fixed_formula)), "term.labels")
  has <- function(term, set) any(strsplit(term, ":", fixed = TRUE)[[1]] %in% set)
  keep <- vapply(tl, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    with_design <- any(parts %in% spec$design_terms)
    if (!with_design) return(FALSE)
    switch(which,
      # A. THE PRIMARY QUESTION. Everything the design can move: the height at
      # the reference time, the trend, and every harmonic. A shape-only test
      # cannot answer "do these groups differ in their trajectories", because a
      # pure vertical separation -- one group simply sleepier all day -- lives
      # entirely in the level block and would be missed.
      full      = TRUE,
      # B. Shape at possibly different heights: trend and harmonics, the pure
      # vertical shift excluded.
      shape     = has(term, spec$basis_terms),
      # C. The rhythm alone.
      circadian = has(term, spec$harm_terms),
      trend     = has(term, spec$trend_terms),
      # "level" is the design's effect on the fitted value AT t = t0, not on a
      # MESOR: with the basis anchored at t0 the design main effects are exactly
      # that, and brief §3 requires the distinction to be stated rather than
      # assumed. A rhythm-adjusted mean is a separate derived quantity.
      level     = !has(term, spec$basis_terms))
  }, logical(1))
  tl[keep]
}

# The three questions, run together and reported in the order they should be
# read: the full trajectory difference FIRST, because it is the primary answer,
# and the other two as the scientifically distinct follow-ups they are.
dance_traj_omnibus_set <- function(fit, which = c("full", "shape", "circadian"),
                                   ...) {
  res <- lapply(which, function(w) {
    r <- dance_traj_omnibus(fit, w, ...)
    r$block <- w
    r$label <- unname(DANCE_TRAJ_BLOCK_LABEL[w])
    r$hypothesis <- unname(DANCE_TRAJ_BLOCK_H0[w])
    r
  })
  names(res) <- which
  ok <- vapply(res, function(r) isTRUE(r$ok), logical(1))
  tab <- if (any(ok)) do.call(rbind, lapply(res[ok], function(r) data.frame(
    block = r$block, label = r$label, n_terms = length(r$terms),
    method = r$method, statistic = r$statistic, df1 = r$df1, df2 = r$df2,
    p = r$p, validated = isTRUE(r$validated), stringsAsFactors = FALSE))) else NULL
  if (!is.null(tab)) rownames(tab) <- NULL
  list(ok = any(ok), results = res, table = tab, primary = "full",
       note = paste(
         "Read the FULL trajectory difference first: it is the primary answer to",
         "'do these groups differ in their trajectories', and it is the only one of",
         "the three that a purely vertical separation can reach. SHAPE excludes that",
         "vertical shift and asks whether the curves have the same form at possibly",
         "different heights. CIRCADIAN narrows further to the harmonics alone. They",
         "are nested -- circadian is inside shape is inside full -- so they are three",
         "different scientific questions, not three attempts at one, and choosing",
         "whichever came out smallest is a multiplicity error with a story attached."))
}

# ------------------------------------------------------------------------------
# THE BLOCK TEST
# ------------------------------------------------------------------------------
# The calibration statement every omnibus result carries. Two of them, because
# the simulation grid covers one harmonic with no trend and nothing else, and a
# result from outside that grid must not borrow its warrant.
DANCE_TRAJ_OMNIBUS_PROVISIONAL <- paste(
  "PROVISIONALLY CALIBRATED -- pilot evidence, not a general warrant. The",
  "CIRCADIAN block was simulated under a true null with 12 time points, one",
  "harmonic, no trend, Gaussian noise and a nominal .05. Rejection rates with",
  "exact (Clopper-Pearson) intervals: 2 x 2 mixed, 16 participants, 0.040",
  "[0.011, 0.099] (n = 100) with curve-level variance and 0.020 [0.002, 0.070]",
  "(n = 100) with participant and curve variance; between-participant only,",
  "0.050 [0.010, 0.139] (n = 60) at 16 participants and 0.017 [0.000, 0.089]",
  "(n = 60) at 32. Power against a 1.2 / 0.5 / 0.3 rad phase shift in one cell",
  "was 1.000 / 0.850 / 0.333 (n = 60).",
  "THOSE INTERVALS ARE THE POINT: at n = 100 the upper limit is around 0.10, so",
  "the evidence rules out the 0.275 this test used to have and does NOT",
  "establish that the size is .05 rather than .09. The pre-phase-4 gate raises",
  "the critical cells to at least 1,000 simulations. Until it clears, read these",
  "p-values as sound enough to act on and not yet as certified.")

DANCE_TRAJ_OMNIBUS_UNCALIBRATED <- paste(
  "OUTSIDE THE SIMULATED GRID. The type-I error of this block test was measured",
  "only for the CIRCADIAN block, with one harmonic, no trend, balanced designs,",
  "Gaussian residuals, the lmer engine and Kenward-Roger. This result is not one",
  "of those: %s. The p-value is a reasonable default, not a validated one, and",
  "the pilot numbers for the circadian block do not transfer to it.")

# Which of the two applies, and why. Three things decide it: the BLOCK (only the
# circadian one was simulated), the MODEL (one harmonic, no trend, balanced),
# and the METHOD (Kenward-Roger). A calibration established for one must not
# propagate to the others by silence, which is what amendment 4 asks for and
# what the first version of this function got wrong by keying on the model
# alone.
dance_traj_calibration <- function(fit, method_kind = "kr", block = "circadian") {
  spec <- fit$spec
  out <- character(0)
  if (!identical(block, "circadian"))
    out <- c(out, sprintf(paste("it is the '%s' block, and only the circadian block",
                                "has been simulated"), block))
  if (!identical(method_kind, "kr"))
    out <- c(out, switch(method_kind,
      satterthwaite = "the denominator degrees of freedom come from Satterthwaite, which the grid did not cover",
      wald          = "the test is an asymptotic Wald/likelihood-ratio statistic, not a Kenward-Roger F",
      lrt           = "the F approximations were unavailable and the test fell back to an asymptotic likelihood-ratio statistic",
      sprintf("the test used the '%s' method rather than Kenward-Roger", method_kind)))
  if ((spec$n_harmonics %||% 1L) > 1L)
    out <- c(out, sprintf("it fits %d harmonics", spec$n_harmonics))
  if (!identical(spec$trend %||% "none", "none"))
    out <- c(out, sprintf(paste("it carries a '%s' trend term%s"), spec$trend,
                          if (identical(spec$trend, "exp_sat"))
                            " whose tau is profiled from the same data, so the reported degrees of freedom do not account for estimating it"
                          else ""))
  if (!identical(fit$engine %||% "lmer", "lmer"))
    out <- c(out, sprintf("it uses the %s engine, where Kenward-Roger is unavailable",
                          fit$engine))
  if (isTRUE(fit$ar1) || !is.null(fit$residual_cor))
    out <- c(out, "it models a residual correlation structure")
  if (length(unique(spec$cells$n_obs)) > 1L)
    out <- c(out, "its design cells are unbalanced")
  if (!isTRUE(spec$time_regular %||% TRUE))
    out <- c(out, "its observations are irregularly spaced in time")
  if (!length(out))
    return(list(validated = TRUE, provisional = TRUE,
                calibration = DANCE_TRAJ_OMNIBUS_PROVISIONAL))
  list(validated = FALSE, provisional = FALSE,
       calibration = sprintf(DANCE_TRAJ_OMNIBUS_UNCALIBRATED,
                             paste(out, collapse = "; ")))
}

dance_traj_omnibus <- function(fit, which = DANCE_TRAJ_BLOCKS,
                               df_method = c("auto", "kr", "satterthwaite"),
                               kr_max_subjects = NULL) {
  which <- match.arg(which, DANCE_TRAJ_BLOCKS)
  df_method <- match.arg(df_method)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  terms_in <- dance_traj_block_terms(fit, which)
  if (!length(terms_in))
    return(list(ok = FALSE, which = which,
                label = unname(DANCE_TRAJ_BLOCK_LABEL[which]),
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
                validated = dance_traj_calibration(fit, "wald", which)$validated,
                calibration = dance_traj_calibration(fit, "wald", which)$calibration,
                caveat = paste("Kenward-Roger is unavailable on this engine; this test is",
                               "asymptotic and is optimistic in small samples.")))
  }

  n_subj <- fit$spec$n_participants
  kr_available <- requireNamespace("pbkrtest", quietly = TRUE)
  # kr_max_subjects is a runtime cap, not a statistical rule. Off by default.
  capped <- !is.null(kr_max_subjects) && is.finite(n_subj) && n_subj > kr_max_subjects
  use_kr <- switch(df_method,
    kr            = kr_available,
    satterthwaite = FALSE,
    auto          = kr_available && !capped)
  if (identical(df_method, "kr") && !kr_available)
    return(list(ok = FALSE, which = which, message = paste(
      "Kenward-Roger was requested but pbkrtest is not installed. Install it, or",
      "pass df_method = 'satterthwaite' to accept the cheaper approximation.")))

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
      cal <- dance_traj_calibration(fit, "kr", which)
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = "Kenward-Roger F",
                  df_method = "kr",
                  method_note = paste("Kenward-Roger: the preferred small-sample F",
                                      "approximation, used because the engine supports it"),
                  statistic = unname(st["Ftest", "stat"]),
                  df1 = unname(st["Ftest", "ndf"]), df2 = unname(st["Ftest", "ddf"]),
                  p = unname(st["Ftest", "p.value"]), caveat = NULL,
                  validated = cal$validated, provisional = cal$provisional,
                  calibration = cal$calibration))
    }
  }

  # Satterthwaite via lmerTest, over the whole block at once
  if (inherits(m, "lmerModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    X <- stats::model.matrix(m)
    keep <- attr(X, "assign") %in% which(attr(stats::terms(m), "term.labels") %in% terms_in)
    L <- diag(ncol(X))[keep, , drop = FALSE]
    ct <- tryCatch(lmerTest::contest(m, L, joint = TRUE), error = function(e) NULL)
    cals <- dance_traj_calibration(fit, "satterthwaite", which)
    if (!is.null(ct))
      return(list(ok = TRUE, which = which, terms = terms_in,
                  method = "Satterthwaite F",
                  df_method = "satterthwaite",
                  method_note = if (capped) sprintf(paste(
                      "Satterthwaite, because a COMPUTATIONAL cap of %d participants was",
                      "set and this fit has %d. That cap is a runtime policy, not a",
                      "statistical threshold -- Kenward-Roger remains the better",
                      "approximation at this n."), kr_max_subjects, n_subj)
                    else if (identical(df_method, "satterthwaite"))
                      "Satterthwaite, requested explicitly as the cheaper alternative to Kenward-Roger"
                    else paste("Satterthwaite, because Kenward-Roger was unavailable or",
                               "failed on this fit"),
                  statistic = ct[["F value"]], df1 = ct[["NumDF"]],
                  df2 = ct[["DenDF"]], p = ct[["Pr(>F)"]], caveat = NULL,
                  validated = cals$validated, provisional = cals$provisional,
                  calibration = cals$calibration))
  }

  # last resort: a likelihood-ratio test on ML refits
  mm <- refit_ml(m)
  m0 <- tryCatch(stats::update(mm, reduced_formula), error = function(e) NULL)
  if (is.null(m0)) return(list(ok = FALSE, which = which,
                               message = "The reduced model could not be fitted."))
  an <- tryCatch(stats::anova(m0, mm), error = function(e) NULL)
  if (is.null(an)) return(list(ok = FALSE, which = which, message = "The block comparison failed."))
  call <- dance_traj_calibration(fit, "lrt", which)
  list(ok = TRUE, which = which, terms = terms_in,
       method = "likelihood-ratio on ML refits",
       statistic = an$Chisq[2], df1 = an$Df[2], df2 = NA_real_, p = an$`Pr(>Chisq)`[2],
       validated = call$validated, provisional = call$provisional,
       calibration = call$calibration,
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

# ==============================================================================
# JOINT-DISTRIBUTION AMPLITUDE AND ACROPHASE INFERENCE
# ==============================================================================
# WHY THE DELTA METHOD IS NOT ENOUGH AS A DEFAULT
#
# A = sqrt(a^2 + b^2) and phi = atan2(b, a) are smooth in (a, b) everywhere
# except at the origin, where both are undefined and their derivatives blow up.
# The delta method linearises around the estimate, so it is accurate exactly
# when the harmonic vector is far from the origin relative to its own scatter --
# and increasingly wrong as it approaches, which is the regime where a user most
# wants to know. Worse, the delta-method acrophase interval is a SYMMETRIC
# interval on a circle: it can be wider than the circle itself and still be
# printed as though it meant something.
#
# The default here instead draws from the joint distribution of the coefficients
#
#     (a, b) ~ N(betahat, V)
#
# transforms EACH DRAW into (amplitude, phase), and summarises the draws --
# amplitude by ordinary quantiles, phase CIRCULARLY, as a quantile of the signed
# angular deviation from the draws' own mean direction. Nothing is linearised,
# the amplitude interval cannot go below zero, and the phase interval is an arc
# rather than a symmetric interval pretending to be one.
#
# The delta method remains available (method = "delta") and is the right choice
# when the vector is clearly separated from zero and speed matters: it is exact
# in the limit and needs no draws. Both are labelled in the result.
#
# THIS IS NOT A BOOTSTRAP. It resamples nothing; it draws from the estimated
# sampling distribution of the fixed effects, which is the multivariate normal
# the model already asserts. It therefore inherits that assumption and does NOT
# correct for it -- a participant-level bootstrap would, at a much higher cost,
# and dance_traj_boot_curves() below is the hook for one.
dance_traj_amp_phase_joint <- function(coefs, conf = 0.95, n_draw = 20000,
                                       seed = 20240601) {
  if (!isTRUE(coefs$ok)) return(coefs)
  if (is.null(coefs$joint))
    return(list(ok = FALSE, message = paste(
      "The joint covariance of the (cos, sin) pair is unavailable, so the joint",
      "method cannot run. Use method = 'delta' and read its caveats.")))
  # A FIXED seed, so the same fit gives the same interval twice. A simulated
  # interval that moves between runs of the same analysis is not reportable, and
  # the seed is returned so it can be stated in a methods section.
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)

  k <- coefs$effective_period / (2 * pi)
  p <- (1 - conf) / 2
  n <- length(coefs$a)
  rows <- lapply(seq_len(n), function(i) {
    V <- coefs$joint[[i]]
    mu <- c(coefs$a[i], coefs$b[i])
    ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
    R <- ev$vectors %*% diag(sqrt(pmax(0, ev$values)), 2) %*% t(ev$vectors)
    Z <- matrix(stats::rnorm(2 * n_draw), ncol = 2) %*% R
    A_d <- sqrt((mu[1] + Z[, 1])^2 + (mu[2] + Z[, 2])^2)
    P_d <- atan2(mu[2] + Z[, 2], mu[1] + Z[, 1])
    A <- sqrt(mu[1]^2 + mu[2]^2)
    phi <- atan2(mu[2], mu[1]) %% (2 * pi)
    A_q <- stats::quantile(A_d, c(p, 1 - p), names = FALSE)
    # PHASE, CIRCULARLY: the mean direction of the draws, then quantiles of the
    # signed angular deviation from it. A linear quantile of atan2 output would
    # split every distribution straddling the +/-pi cut into two halves and
    # report an interval covering everything except the truth.
    mdir <- atan2(mean(sin(P_d)), mean(cos(P_d)))
    dev <- ((P_d - mdir + pi) %% (2 * pi)) - pi
    d_q <- stats::quantile(dev, c(p, 1 - p), names = FALSE)
    arc <- d_q[2] - d_q[1]
    # BINGHAM'S RULE, STATED PROPERLY FOR THE JOINT METHOD.
    # ------------------------------------------------------------------
    # The first draft of this line asked whether the lower amplitude quantile
    # was above zero. It is essentially always above zero -- an amplitude is a
    # norm, so it is non-negative and hits zero with probability zero -- and the
    # rule therefore never fired. The delta-method version of the same rule
    # works only because its interval is a symmetric linearisation that CAN go
    # negative.
    #
    # The quantity that actually decides identifiability is whether the joint
    # confidence REGION for (a, b) contains the origin. That is a Hotelling
    # statistic on two degrees of freedom, and it is the exact joint-method
    # statement of what Bingham et al. (1982) require: the phase of a vector
    # that may be the zero vector is not a quantity.
    Vi <- tryCatch(solve(V), error = function(e) NULL)
    T2 <- if (is.null(Vi)) NA_real_ else as.numeric(t(mu) %*% Vi %*% mu)
    excludes_origin <- is.finite(T2) && T2 > stats::qchisq(conf, 2)
    # The arc is kept as a second, independent signal: a region that excludes
    # the origin but only just will still produce an arc near the whole circle.
    defined <- excludes_origin && arc < 2 * pi * 0.99
    data.frame(
      cell = coefs$cells[i], beta_cos = mu[1], beta_sin = mu[2],
      amplitude = A, amplitude_lo = A_q[1], amplitude_hi = A_q[2],
      amplitude_se = stats::sd(A_d),
      acrophase_rad = phi, acrophase_time = phi * k,
      acrophase_lo = if (defined) ((mdir + d_q[1]) %% (2 * pi)) * k else NA_real_,
      acrophase_hi = if (defined) ((mdir + d_q[2]) %% (2 * pi)) * k else NA_real_,
      acrophase_arc_time = if (defined) arc * k else NA_real_,
      acrophase_wraps = defined && ((mdir + d_q[1]) %% (2 * pi)) > ((mdir + d_q[2]) %% (2 * pi)),
      amplitude_T2 = T2, excludes_origin = excludes_origin,
      phase_defined = defined, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  list(ok = TRUE, harmonic = coefs$harmonic, method_kind = "joint",
       effective_period = coefs$effective_period, table = out,
       n_draw = n_draw, seed = seed,
       method = sprintf(paste("%s draws from the joint normal distribution of each",
                              "cell's (cos, sin) pair, transformed per draw; amplitude",
                              "by quantile, phase circularly. %.0f%% intervals."),
                        format(n_draw, big.mark = ","), 100 * conf),
       undefined_phase = sum(!out$phase_defined),
       note = paste(
         "The acrophase interval is an ARC, reported low-to-high around the circle:",
         "when acrophase_wraps is TRUE the low bound is numerically greater than the",
         "high one and the interval runs through the period boundary. Its LENGTH is",
         "acrophase_arc_time, which is the number to quote.",
         sprintf(paste("Identifiability is decided by whether the joint %.0f%% confidence",
                       "region for (cos, sin) EXCLUDES THE ORIGIN -- a Hotelling statistic",
                       "on 2 df, reported as amplitude_T2 against a cutoff of %.2f."),
                 100 * conf, stats::qchisq(conf, 2)),
         if (any(!out$phase_defined)) paste(
           "Cells with phase_defined = FALSE have a region that covers the origin, so",
           "their phase is not identified and no interval is given -- the rule",
           "Bingham et al. (1982) state, here as a joint region rather than inferred",
           "from a linearisation of one coordinate.") else NULL))
}

# The dispatcher. `joint` is the default because it is right in the regime the
# delta method is wrong in, and the two agree where the delta method is right.
dance_traj_amp_phase_ci <- function(coefs, conf = 0.95,
                                    method = c("joint", "delta"), ...) {
  method <- match.arg(method)
  r <- if (identical(method, "joint")) dance_traj_amp_phase_joint(coefs, conf, ...)
       else dance_traj_amp_phase(coefs, conf)
  if (isTRUE(r$ok) && identical(method, "delta"))
    r$method_kind <- "delta"
  r
}

# ==============================================================================
# MODEL SELECTION IS NOT CONFIRMATORY INFERENCE
# ==============================================================================
# AIC, AICc, BIC, Delta-AIC and Akaike weights all stay: choosing between two
# harmonics and one, or between a linear and a saturating trend, is a real
# question and these are the right tools for it.
#
# What does NOT survive the choosing is the p-value afterwards. If the SAME data
# picked the number of harmonics, the trend, the interaction structure or the
# random structure, then the p-value from the selected model is computed as
# though that model had been written down in advance. It was not. The reference
# distribution is the one for a fixed model; the actual procedure included a
# search; and the p-value is optimistic by an amount that depends on how wide
# the search was and is not recoverable from the final fit.
#
# This is not a reason to stop selecting. It is a reason to say which mode the
# analysis is in, so a reader is not handed an exploratory number wearing a
# confirmatory p-value's clothes. Two modes, labelled:
#
#   PRE-SPECIFIED   the model was fixed before seeing these data. p-values mean
#                   what they say (within the calibration grid).
#   DATA-SELECTED   something about the model was chosen using these data. The
#                   p-values are CONDITIONAL ON THE SELECTED MODEL and are not
#                   corrected for the selection.
#
# There is no adjustment applied here, because the honest ones (sample splitting,
# selective inference, a bootstrap that repeats the whole search) all change the
# procedure rather than the number, and inventing a correction factor would be
# worse than saying plainly what happened.
dance_traj_selection_state <- function(selected = character(0)) {
  known <- c(harmonics = "the number of harmonics",
             trend = "the trend type",
             tau = "the saturating trend's tau",
             interaction = "the fixed-effect interaction structure",
             random = "the random-effects structure",
             period = "the period")
  selected <- intersect(selected, names(known))
  if (!length(selected))
    return(list(mode = "pre-specified", selected = character(0),
                post_selection = FALSE,
                label = "Pre-specified model",
                note = paste(
                  "The model was fixed before these data were examined, so the p-values",
                  "are ordinary confirmatory p-values -- within the calibration grid the",
                  "omnibus reports.")))
  list(mode = "data-selected", selected = selected, post_selection = TRUE,
       label = "Data-selected model (exploratory)",
       note = paste(
         sprintf("%s %s chosen using these same data.",
                 paste(unname(known[selected]), collapse = ", "),
                 if (length(selected) == 1L) "was" else "were"),
         "Every p-value below is therefore CONDITIONAL ON THE SELECTED MODEL: it is",
         "computed against the reference distribution of a model written down in",
         "advance, which this one was not, so it is optimistic by an amount that",
         "depends on how wide the search was and cannot be recovered from the final",
         "fit. No correction is applied, because the honest remedies -- sample",
         "splitting, selective inference, or a bootstrap that repeats the entire",
         "search -- change the procedure rather than the number.",
         "Report these as exploratory, or re-specify the model and collect new data."))
}

# The ladder walk is itself a selection over random structures whenever it
# descends, so a fit that simplified is data-selected even if the user
# pre-specified everything else. Read off the fit rather than trusted to memory.
dance_traj_selection_from_fit <- function(fit, selected = character(0)) {
  s <- selected
  if (isTRUE(fit$simplified)) s <- union(s, "random")
  if (isTRUE(fit$tau_estimated)) s <- union(s, "tau")
  dance_traj_selection_state(s)
}
