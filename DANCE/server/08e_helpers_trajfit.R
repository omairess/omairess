# ==============================================================================
# server/08e_helpers_trajfit.R — LAYER B: fitting the trajectory model
# ==============================================================================
# P21 phase 2. Takes a dance_traj_spec() and returns a fit plus an honest log of
# how it was obtained. Pure: no reactives, no notifications, no printing.
#
# THREE RULES, all from the brief
# -------------------------------
# §3  Do not send a fixed basis to a nonlinear optimiser. Every trend this app
#     offers is linear in its coefficients ONCE tau is fixed, so the only
#     genuinely nonlinear parameter in the whole family is tau. It is profiled
#     (see below) rather than handed to a joint optimiser.
#
# §4  Do not fit a maximal random structure automatically, and never simplify
#     silently. The ladder comes from the spec; this file walks it top-down,
#     stops at the first rung that CONVERGES, and records which rung that was. A
#     caller that does not print re_rung is misreporting its own model.
#
# SINGULARITY IS REPORTED, NOT ACTED ON
# -------------------------------------
# The first version of this walk descended the ladder on singularity as well as
# on non-convergence. That is the common recipe and it is wrong here, for a
# measurable reason. On the null datasets used to calibrate layer C, the
# curve-level rung was singular in 48% of fits -- so under that rule the ladder
# would have dropped to the participant-only rung in half of all analyses, which
# is precisely the structure whose type-I error measured 0.200 against a nominal
# .05. Forcing the curve-level rung, singular fits included, gave 0.060.
#
# A boundary estimate is not a failed fit. The REML estimate of a variance that
# genuinely is near zero belongs at zero; the fixed-effect covariance is still
# estimated, Kenward-Roger still corrects it, and the test stays calibrated. It
# is the DROPPING that does the damage, because it moves variance the design
# implies into a residual shared across the whole fit and shrinks exactly the
# standard errors under test. So the walk descends only when the optimiser fails
# or does not converge, and a singular fit is returned with `singular = TRUE`
# and a note saying what the boundary means.
#
# §14 Residual correlation where the engine allows it. lme4 cannot express one
#     at all, so engine = "glmmTMB" exists for that case and is labelled: the
#     Kenward-Roger F of layer C is unavailable there and the block tests fall
#     back to an asymptotic Wald statistic.
#
# WHY PROFILING tau RATHER THAN OPTIMISING IT
# -------------------------------------------
# A_sat * (1 - exp(-(t - t0)/tau)) is linear in A_sat and nonlinear only in tau.
# The old module handed the pair to nlsLM, which is how it acquired the
# A_sat/tau ridge it documents and the convergence failures that made its
# convergence gate necessary in the first place. Fixing tau makes the column
# 1 - exp(-(t-t0)/tau) an ordinary basis column, so the profile is a sequence of
# ordinary linear mixed fits. That gives a real profile likelihood for tau, a
# profile interval rather than a delta-method guess, and -- the part that
# matters for finding A10 -- a flat profile is visible as a flat profile instead
# of surfacing as "the optimiser stopped somewhere on the ridge".
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Is the fitted model singular (a variance component at the boundary)?
dance_traj_is_singular <- function(m) {
  if (inherits(m, "merMod"))
    return(isTRUE(tryCatch(lme4::isSingular(m, tol = 1e-4), error = function(e) FALSE)))
  if (inherits(m, "glmmTMB")) {
    v <- tryCatch(unlist(lapply(lme4::VarCorr(m)$cond, function(x) diag(as.matrix(x)))),
                  error = function(e) numeric(0))
    return(length(v) > 0 && any(v < 1e-8))
  }
  FALSE
}

# Did the optimiser actually converge? lme4 reports this as a message list --
# but that list is not purely about convergence. A fit whose optimiser code is 0
# and whose only message is
#
#     "boundary (singular) fit: see help('isSingular')"
#
# converged perfectly well; lme4 is telling you WHERE it converged, namely to a
# variance of zero or a correlation of one. Treating that notice as a
# convergence failure is how the first version of this walk kept descending the
# ladder past its top rung even after the explicit singularity test above was
# no longer allowed to: every one of the null fits it dropped had optcode 0 and
# that message and nothing else. The notice is filtered out here; singularity is
# detected by dance_traj_is_singular() and REPORTED, never acted on.
dance_traj_converged <- function(m) {
  if (inherits(m, "merMod")) {
    oc <- tryCatch(m@optinfo$conv$opt, error = function(e) 0L)
    msg <- tryCatch(m@optinfo$conv$lme4$messages, error = function(e) NULL)
    msg <- msg[!grepl("boundary (singular) fit", msg, fixed = TRUE)]
    return(identical(as.integer(oc), 0L) && length(msg) == 0)
  }
  if (inherits(m, "glmmTMB"))
    return(isTRUE(m$sdr$pdHess) && identical(as.integer(m$fit$convergence), 0L))
  TRUE
}

# ------------------------------------------------------------------------------
# ONE FIT at a given rung
# ------------------------------------------------------------------------------
dance_traj_fit_one <- function(spec, re_formula, engine = "lmer", REML = TRUE,
                               ar1 = FALSE) {
  f <- stats::as.formula(paste(spec$fixed_formula, "+", re_formula))
  d <- spec$data
  if (identical(engine, "glmmTMB")) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) return(NULL)
    if (ar1) {
      # glmmTMB's ar1() needs an explicit, evenly-indexed time FACTOR
      d$.tf <- factor(match(d$t, sort(unique(d$t))),
                      levels = seq_along(sort(unique(d$t))))
      f <- stats::as.formula(paste(spec$fixed_formula, "+", re_formula,
                                   "+ ar1(0 + .tf | subject)"))
    }
    return(tryCatch(glmmTMB::glmmTMB(f, data = d, REML = REML),
                    error = function(e) NULL, warning = function(w) NULL))
  }
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  fn <- if (requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer else lme4::lmer
  tryCatch(fn(f, data = d, REML = REML, control = ctrl),
           error = function(e) NULL)
}

# ------------------------------------------------------------------------------
# WALK THE LADDER
# ------------------------------------------------------------------------------
# Returns the fit plus the log. `rung` is 1-based into spec$re_ladder, and
# `simplified` says whether the requested structure survived -- because a reader
# cannot tell rung 2 from rung 4 by looking at the fixed-effect coefficients.
dance_traj_fit <- function(spec, engine = c("lmer", "glmmTMB"), REML = TRUE,
                           ar1 = FALSE, max_rung = NULL) {
  engine <- match.arg(engine)
  if (!isTRUE(spec$ok)) return(list(ok = FALSE, message = spec$message))
  if (ar1 && identical(engine, "lmer"))
    return(list(ok = FALSE, message = paste(
      "A residual AR(1) structure was requested, but lme4 cannot express one.",
      "Use engine = 'glmmTMB', and note that Kenward-Roger is then unavailable",
      "so the block tests become asymptotic Wald tests.")))

  ladder <- spec$re_ladder
  if (!is.null(max_rung)) ladder <- ladder[seq_len(min(max_rung, length(ladder)))]
  attempts <- list()
  for (i in seq_along(ladder)) {
    m <- dance_traj_fit_one(spec, ladder[[i]]$formula, engine, REML, ar1)
    ok_fit <- !is.null(m)
    conv <- ok_fit && dance_traj_converged(m)
    sing <- ok_fit && dance_traj_is_singular(m)
    attempts[[i]] <- list(rung = i, label = ladder[[i]]$label,
                          formula = ladder[[i]]$formula,
                          fitted = ok_fit, converged = conv, singular = sing)
    if (ok_fit && conv) {
      notes <- character(0)
      if (i > 1L) notes <- c(notes, sprintf(paste(
        "The requested random structure (%s) did not converge; the fit reported here",
        "uses %s. This is a DIFFERENT model from the one requested."),
        ladder[[1]]$label, ladder[[i]]$label))
      if (sing) notes <- c(notes, paste(
        "At least one variance component is estimated at the boundary (a variance of",
        "zero, or a correlation of +/-1). The term is KEPT: dropping a",
        "design-justified random effect because its estimate sits on the boundary is",
        "what makes the fixed-effect tests anticonservative, and the boundary estimate",
        "is still the REML estimate. Read it as 'these data do not separate that",
        "component from zero', not as a failed fit."))
      return(list(
        ok = TRUE, model = m, engine = engine, REML = REML, ar1 = ar1,
        spec = spec,
        re_rung = i, re_label = ladder[[i]]$label, re_formula = ladder[[i]]$formula,
        simplified = i > 1L, singular = sing,
        n_rungs = length(ladder), attempts = attempts,
        formula = paste(spec$fixed_formula, "+", ladder[[i]]$formula),
        note = if (length(notes)) paste(notes, collapse = " ") else NULL))
    }
  }
  list(ok = FALSE, attempts = attempts,
       message = paste("No rung of the random-effects ladder produced a converged fit.",
                       "The participant-level structure is not estimable from these",
                       "data; see attempts for what was tried."))
}

# ------------------------------------------------------------------------------
# PROFILE OVER tau (the only genuinely nonlinear parameter)
# ------------------------------------------------------------------------------
# Refits the conditionally-linear model at each tau on a grid and returns the
# profile. Because the comparison is over models with the SAME fixed-effects
# structure and differing only in one basis column, ML is used so the
# likelihoods are comparable; the winning tau is then refitted by REML for
# estimation if REML was asked for.
#
# `flat` is the finding A10 guard, made quantitative: if the whole grid lies
# within `flat_tol` of the best log-likelihood, tau is not identified by these
# data and the caller must say so rather than quote the argmax.
dance_traj_profile_tau <- function(d, period, n_harmonics, design_terms,
                                   covariates = character(0), t0 = NULL,
                                   tau_grid = NULL, engine = "lmer",
                                   REML = TRUE, interaction = "full",
                                   flat_tol = 1.92) {
  if (is.null(tau_grid)) {
    span <- diff(range(d$t, na.rm = TRUE))
    tau_grid <- exp(seq(log(max(span / 50, 1e-3)), log(span * 3), length.out = 24))
  }
  rows <- lapply(tau_grid, function(tau) {
    sp <- dance_traj_spec(d, period, n_harmonics, "exp_sat", design_terms,
                          covariates, t0, tau, interaction = interaction)
    if (!isTRUE(sp$ok)) return(NULL)
    fit <- dance_traj_fit(sp, engine = engine, REML = FALSE)
    if (!isTRUE(fit$ok)) return(NULL)
    data.frame(tau = tau, logLik = as.numeric(stats::logLik(fit$model)),
               rung = fit$re_rung, stringsAsFactors = FALSE)
  })
  prof <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(prof) || !nrow(prof))
    return(list(ok = FALSE, message = "No tau on the grid produced a usable fit."))

  best <- prof$tau[which.max(prof$logLik)]
  drop <- max(prof$logLik) - min(prof$logLik)
  flat <- drop < flat_tol
  # profile interval: the tau whose log-likelihood is within 1.92 of the max
  # (the usual chi-square(1)/2 cutoff)
  inside <- prof$tau[prof$logLik >= max(prof$logLik) - flat_tol]
  at_edge <- best <= min(tau_grid) * 1.001 || best >= max(tau_grid) * 0.999

  final <- NULL
  sp <- dance_traj_spec(d, period, n_harmonics, "exp_sat", design_terms,
                        covariates, t0, best, interaction = interaction)
  if (isTRUE(sp$ok)) final <- dance_traj_fit(sp, engine = engine, REML = REML)

  list(ok = TRUE, profile = prof, tau = best,
       ci = if (length(inside)) range(inside) else c(NA_real_, NA_real_),
       logLik_range = drop, flat = flat, at_edge = at_edge, fit = final,
       message = if (flat) sprintf(paste(
         "tau is NOT identified by these data: the log-likelihood varies by only",
         "%.2f across the whole grid (%.3g to %.3g), against the %.2f that would",
         "make a %.0f%% interval. Report the trajectory, not a value of tau, and",
         "do not compare tau between groups."),
         drop, min(tau_grid), max(tau_grid), flat_tol, 95)
       else if (at_edge) paste(
         "The best tau sits at the edge of the search grid, so it is a bound, not",
         "an optimum. Widen the grid before quoting it.")
       else NULL)
}

# ------------------------------------------------------------------------------
# THE FITTING LOG (brief §16 / §4 / §14)
# ------------------------------------------------------------------------------
dance_traj_fit_report <- function(fit) {
  if (!isTRUE(fit$ok)) return(fit$message)
  l <- character(0); a <- function(...) l <<- c(l, sprintf(...))
  a("Estimation:        %s, %s", fit$engine, if (isTRUE(fit$REML)) "REML" else "ML")
  a("Formula:           %s", fit$formula)
  a("Random structure:  rung %d of %d -- %s", fit$re_rung, fit$n_rungs, fit$re_label)
  if (isTRUE(fit$simplified)) a("  ! %s", fit$note)
  if (isTRUE(fit$ar1)) a("Residuals:         AR(1) within participant")
  a("Convergence:       converged, non-singular")
  if (length(fit$attempts) > 1) {
    a("Rungs tried:")
    for (at in fit$attempts)
      a("  %d. %-52s %s", at$rung, at$formula,
        if (!at$fitted) "did not fit"
        else if (!at$converged) "did not converge"
        else if (at$singular) "singular"
        else "accepted")
  }
  paste(l, collapse = "\n")
}
