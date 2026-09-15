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


# Refit this fit's model with a different fixed-effect formula and/or likelihood,
# from the spec rather than from the recorded call. Used wherever a reduced model
# is needed: the glmmTMB comparison and the last-resort likelihood-ratio test.
dance_traj_refit <- function(fit, drop_terms = character(0), REML = fit$REML) {
  spec <- fit$spec
  ff <- spec$fixed_formula
  if (length(drop_terms)) ff <- paste(ff, "-", paste(drop_terms, collapse = " - "))
  dance_traj_fit_one(spec, fit$re_formula, fit$engine %||% "lmer", REML,
                     fit$residual_cor, fixed_formula = ff)
}

# ------------------------------------------------------------------------------
# FIVE STATUSES, NOT ONE
# ------------------------------------------------------------------------------
# "It didn't work" covers five different things, and collapsing them is how the
# type-I error of layer C got to 0.275: a singularity notice was read as a
# convergence failure and silently descended the ladder. They are kept apart:
#
#   optimizer_failure  the fit threw, or the optimiser returned a non-zero code
#   converged          the optimiser reached a stationary point
#   boundary/singular  it converged, ON a boundary: a variance at zero or a
#                      correlation at +/-1. A successfully optimised model.
#   rank_deficient     the FIXED-effect design matrix lost columns; the model
#                      fitted is not the model requested
#
# A singular fit is not descended past automatically -- see the note below --
# but nor is it asserted to be the right model. What it is, is RECORDED: which
# dimensions collapsed, what was fitted, whether anything was simplified and
# why. The final choice is design logic plus simulation evidence, which is what
# the validation gate is for; it is not a blanket rule in either direction.
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

# WHICH random-effect dimensions collapsed. rePCA gives the eigen-decomposition
# of each grouping factor's relative covariance; a zero eigenvalue is a
# dimension the data cannot support, and the count of them is far more useful
# than the yes/no of isSingular. Reported so a reader can see that (say) the
# correlation between the participant intercept and the participant cos term
# went to 1, rather than being told only that "the fit is singular".
dance_traj_re_collapse <- function(m) {
  if (!inherits(m, "merMod")) return(NULL)
  pc <- tryCatch(lme4::rePCA(m), error = function(e) NULL)
  if (is.null(pc)) return(NULL)
  out <- lapply(names(pc), function(g) {
    sd <- pc[[g]]$sdev
    tol <- 1e-4 * max(sd, 1)
    list(group = g, n_dim = length(sd), sdev = sd,
         n_collapsed = sum(sd < tol),
         variance_explained = if (sum(sd^2) > 0) cumsum(sd^2) / sum(sd^2) else NA_real_)
  })
  names(out) <- names(pc)
  out
}

# Did the FIXED-effect design matrix lose columns? A rank-deficient fit is not a
# singular fit: the model that was estimated is missing terms the formula asked
# for, so a block test over those terms is testing something else.
dance_traj_rank_deficient <- function(m) {
  b <- tryCatch(if (inherits(m, "glmmTMB")) lme4::fixef(m)$cond else lme4::fixef(m),
                error = function(e) NULL)
  if (is.null(b)) return(list(deficient = FALSE, dropped = character(0)))
  X <- tryCatch(stats::model.matrix(m), error = function(e) NULL)
  dropped <- if (is.null(X)) character(0) else setdiff(colnames(X), names(b))
  dropped <- c(dropped, names(b)[!is.finite(b)])
  list(deficient = length(dropped) > 0, dropped = unique(dropped))
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

# The whole status of one fitted model, as five separate facts.
dance_traj_status <- function(m) {
  if (is.null(m)) return(list(fitted = FALSE, converged = FALSE, singular = NA,
                              boundary_dims = 0L, rank_deficient = NA,
                              dropped_terms = character(0), collapse = NULL,
                              optimizer_failure = TRUE,
                              summary = "the optimiser failed or threw"))
  conv <- dance_traj_converged(m)
  sing <- dance_traj_is_singular(m)
  coll <- dance_traj_re_collapse(m)
  rk <- dance_traj_rank_deficient(m)
  nb <- if (is.null(coll)) 0L else sum(vapply(coll, function(x) x$n_collapsed, integer(1)))
  list(fitted = TRUE, converged = conv, singular = sing,
       boundary_dims = nb,
       rank_deficient = rk$deficient, dropped_terms = rk$dropped,
       collapse = coll, optimizer_failure = !conv,
       summary = paste(c(
         if (conv) "converged" else "the optimiser did NOT converge",
         if (sing) sprintf("on a boundary (%d random-effect dimension%s collapsed)",
                           nb, if (nb == 1L) "" else "s"),
         if (rk$deficient) sprintf("with a RANK-DEFICIENT fixed-effect matrix (dropped: %s)",
                                   paste(rk$dropped, collapse = ", "))),
         collapse = ", "))
}

# ------------------------------------------------------------------------------
# ONE FIT at a given rung
# ------------------------------------------------------------------------------
dance_traj_fit_one <- function(spec, re_formula, engine = "lmer", REML = TRUE,
                               cor_struct = NULL, fixed_formula = NULL) {
  # fixed_formula overrides the spec's, so a REDUCED model can be built the same
  # way the full one was. stats::update() cannot do this job here: the recorded
  # call refers to `d` and `REML` by name, which exist only in this function's
  # frame, so update() re-evaluates it somewhere those are not found and fails
  # with "object 'd' not found". Rebuilding from the spec has no such dependency.
  fixed_formula <- fixed_formula %||% spec$fixed_formula
  f <- stats::as.formula(paste(fixed_formula, "+", re_formula))
  d <- spec$data
  if (identical(engine, "glmmTMB")) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) return(NULL)
    if (!is.null(cor_struct)) {
      # THE CORRELATION RESETS AT THE CURVE, NOT THE PARTICIPANT. Two
      # measurements in different sessions of the same person are not
      # consecutive points of one autocorrelated series, and grouping the
      # structure by subject would say they were -- borrowing a lag-1
      # correlation across a gap of weeks.
      grp <- spec$curve_var %||% "subject"
      if (identical(cor_struct, "ar1")) {
        # glmmTMB's ar1() needs an explicit, evenly-indexed time FACTOR, which
        # is exactly why it is only offered on an even grid: the index is the
        # lag, so an uneven grid would silently call unequal gaps equal.
        d$.tf <- factor(match(d$t, sort(unique(d$t))),
                        levels = seq_along(sort(unique(d$t))))
        f <- stats::as.formula(paste(fixed_formula, "+", re_formula,
                                     sprintf("+ ar1(0 + .tf | %s)", grp)))
      } else {
        # ou() is the continuous-time (Ornstein-Uhlenbeck) analogue: the
        # correlation decays in ACTUAL elapsed time, so an uneven grid is
        # handled as an uneven grid.
        d$.tn <- glmmTMB::numFactor(d$t)
        f <- stats::as.formula(paste(fixed_formula, "+", re_formula,
                                     sprintf("+ ou(0 + .tn | %s)", grp)))
      }
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
# RESIDUAL CORRELATION: DISCRETE AND CONTINUOUS ARE NOT THE SAME MODEL
# --------------------------------------------------------------------
# A discrete AR(1) indexes observations by POSITION: rows 1 and 2 are one lag
# apart whether that is 20 minutes or six hours. On an even grid position and
# elapsed time are the same thing and AR(1) is the natural choice. On an uneven
# one they are not, and fitting AR(1) anyway asserts that a long gap and a short
# gap carry the same correlation -- a claim about the data that nothing in the
# data supports.
#
# So `residual_cor` takes "ar1", "ou" (the continuous-time Ornstein-Uhlenbeck /
# CAR(1) analogue, where correlation decays in actual elapsed time) or "auto",
# which reads spec$time_regular and picks. Asking for "ar1" on an irregular grid
# is REFUSED rather than quietly honoured. Both are grouped by the CURVE, so the
# series resets at each session rather than running across a participant's
# separate visits.
dance_traj_fit <- function(spec, engine = c("lmer", "glmmTMB"), REML = TRUE,
                           residual_cor = c("none", "auto", "ar1", "ou"),
                           ar1 = FALSE, max_rung = NULL) {
  engine <- match.arg(engine)
  residual_cor <- match.arg(residual_cor)
  if (isTRUE(ar1) && identical(residual_cor, "none")) residual_cor <- "ar1"
  if (!isTRUE(spec$ok)) return(list(ok = FALSE, message = spec$message))

  regular <- isTRUE(spec$time_regular %||% TRUE)
  if (identical(residual_cor, "auto"))
    residual_cor <- if (regular) "ar1" else "ou"
  if (identical(residual_cor, "ar1") && !regular)
    return(list(ok = FALSE, message = sprintf(paste(
      "A discrete AR(1) was requested, but these observations are NOT evenly",
      "spaced (gaps run from %.4g to %.4g %s). A discrete AR(1) indexes by",
      "position, so it would treat the longest gap and the shortest as one lag",
      "apart. Use residual_cor = 'ou' for the continuous-time structure, whose",
      "correlation decays in actual elapsed time, or 'auto' to let the spec",
      "choose."), spec$time_gaps[1], spec$time_gaps[2], spec$time_units %||% "units")))
  cor_struct <- if (identical(residual_cor, "none")) NULL else residual_cor
  if (!is.null(cor_struct) && identical(engine, "lmer"))
    return(list(ok = FALSE, message = sprintf(paste(
      "A residual %s structure was requested, but lme4 cannot express one.",
      "Use engine = 'glmmTMB', and note that Kenward-Roger is then unavailable",
      "so the block tests become asymptotic Wald tests."), toupper(cor_struct))))

  ladder <- spec$re_ladder
  if (!is.null(max_rung)) ladder <- ladder[seq_len(min(max_rung, length(ladder)))]
  attempts <- list()
  for (i in seq_along(ladder)) {
    m <- dance_traj_fit_one(spec, ladder[[i]]$formula, engine, REML, cor_struct)
    st <- dance_traj_status(m)
    attempts[[i]] <- c(list(rung = i, label = ladder[[i]]$label,
                            formula = ladder[[i]]$formula), st)
    if (st$fitted && st$converged) {
      notes <- character(0)
      if (i > 1L) notes <- c(notes, sprintf(paste(
        "The requested random structure (%s) did not converge; the fit reported here",
        "uses %s. This is a DIFFERENT model from the one requested."),
        ladder[[1]]$label, ladder[[i]]$label))
      if (isTRUE(st$singular)) notes <- c(notes, sprintf(paste(
        "%d random-effect dimension%s estimated at the boundary (a variance of zero,",
        "or a correlation of +/-1). The term is KEPT here, because dropping a",
        "design-justified random effect because its estimate sits on the boundary is",
        "what makes the fixed-effect tests anticonservative, and a boundary estimate",
        "is still the REML estimate of a variance that is genuinely near zero. That",
        "is NOT a claim that the maximal structure is the scientifically right one --",
        "read $re_collapse for which dimensions went, and see the validation gate for",
        "the evidence behind the rule."),
        st$boundary_dims, if (st$boundary_dims == 1L) " is" else "s are"))
      if (isTRUE(st$rank_deficient)) notes <- c(notes, sprintf(paste(
        "The FIXED-effect design matrix is rank deficient: %s could not be estimated",
        "and were dropped. This is not a singularity -- the model that was fitted is",
        "missing terms the formula asked for, so a block test over those terms is",
        "testing something other than what it names."),
        paste(st$dropped_terms, collapse = ", ")))
      return(list(
        ok = TRUE, model = m, engine = engine, REML = REML,
        ar1 = identical(cor_struct, "ar1"), residual_cor = cor_struct,
        residual_cor_group = if (is.null(cor_struct)) NULL else (spec$curve_var %||% "subject"),
        spec = spec,
        re_rung = i, re_label = ladder[[i]]$label, re_formula = ladder[[i]]$formula,
        simplified = i > 1L,
        # the five statuses, kept apart
        converged = st$converged, singular = st$singular,
        boundary_dims = st$boundary_dims, re_collapse = st$collapse,
        rank_deficient = st$rank_deficient, dropped_terms = st$dropped_terms,
        status = st$summary,
        simplified_reason = if (i > 1L)
          vapply(attempts[seq_len(i - 1L)], function(a) a$summary, character(1)) else NULL,
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
  # A FLAT PROFILE SUPPRESSES THE FIT, it does not merely annotate it. If the
  # likelihood cannot tell one tau from another, then "the fit at the best tau"
  # is the fit at an arbitrary point of a ridge, and handing it back as a
  # fitted model invites exactly the substantive inference the flatness rules
  # out. The profile itself is still returned, so the flatness is visible.
  if (!flat) {
    sp <- dance_traj_spec(d, period, n_harmonics, "exp_sat", design_terms,
                          covariates, t0, best, interaction = interaction)
    if (isTRUE(sp$ok)) final <- dance_traj_fit(sp, engine = engine, REML = REML)
  }

  # TAU WAS ESTIMATED, AND THE FIT BELOW DOES NOT KNOW THAT.
  # ------------------------------------------------------------------
  # `final` is an ordinary linear mixed fit at tau = best. Conditional on tau it
  # is exactly right, and Kenward-Roger on it is exactly the KR of a model with
  # a known basis column. But tau was not known -- it was read off this same
  # data, and every downstream degree of freedom is computed as though it had
  # been handed down from outside. The effect is one-directional: intervals are
  # too narrow and p-values too small, by an amount this function does not
  # quantify.
  #
  # The profile interval below IS a proper interval for tau (a likelihood-ratio
  # interval, not a delta-method one), so the uncertainty is measured -- it is
  # simply not propagated into the fixed-effect inference. Propagating it needs
  # the profiling repeated inside a participant-level bootstrap, which is the
  # honest fix and is not in this phase. Until then dance_traj_calibration()
  # puts any exp_sat fit OUTSIDE the validated grid by name, and this object
  # carries the warning so a caller cannot pick up `fit` without it.
  if (!is.null(final) && isTRUE(final$ok)) {
    final$tau_estimated <- TRUE
    final$tau_profile_ci <- if (length(inside)) range(inside) else c(NA_real_, NA_real_)
    final$tau_warning <- paste(
      "tau was ESTIMATED from these data by profiling, and the inference in this",
      "fit is CONDITIONAL ON THE SELECTED VALUE. The reported degrees of freedom,",
      "standard errors and p-values do not account for having chosen tau, so they",
      "are optimistic: intervals too narrow, p-values too small, by an amount not",
      "quantified here.",
      if (length(inside)) sprintf(
        "The profile-likelihood interval for tau is [%.3g, %.3g] -- read that as the",
        min(inside), max(inside)) else NULL,
      if (length(inside)) "uncertainty the fixed-effect inference is NOT carrying." else NULL,
      "Propagating it requires repeating the profile inside a participant-level",
      "bootstrap.")
  }

  list(ok = TRUE, profile = prof, tau = best,
       ci = if (length(inside)) range(inside) else c(NA_real_, NA_real_),
       logLik_range = drop, flat = flat, at_edge = at_edge, fit = final,
       tau_estimated = TRUE,
       conditional_inference_warning = if (!is.null(final) && isTRUE(final$ok))
         final$tau_warning else NULL,
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
  # BRIEF 8. This line used to be the literal string "converged, non-singular",
  # printed whatever the fit had actually done -- while dance_traj_status() had
  # already separated the five states and dance_traj_fit() had stored every one
  # of them on the object three lines above. In the configuration the validation
  # grid covers, 98% of fits are singular; every one of them was reported here as
  # non-singular. The states are reported separately because they mean different
  # things: a singular fit is a converged fit AT a boundary and its fixed-effect
  # tests stand, a rank-deficient one is missing fixed-effect columns the formula
  # asked for and its block tests do not test what they name.
  a("Convergence:       %s", if (isTRUE(fit$converged)) "converged" else
                             "DID NOT CONVERGE")
  a("Singular:          %s", if (isTRUE(fit$singular))
      sprintf("yes -- %d random-effect dimension%s at a boundary (term kept)",
              fit$boundary_dims %||% 0L,
              if (identical(fit$boundary_dims %||% 0L, 1L)) "" else "s")
    else "no")
  a("Fixed-effect rank: %s", if (isTRUE(fit$rank_deficient))
      sprintf("DEFICIENT -- dropped: %s", paste(fit$dropped_terms, collapse = ", "))
    else "full")
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
