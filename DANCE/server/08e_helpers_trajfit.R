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
                     fit$residual_cor, fixed_formula = ff, phi = fit$phi)
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
# AUTOCORRELATED RESIDUALS: THE CONTINUOUS-TIME AR(1), BY WHITENING
# ------------------------------------------------------------------------------
# Within a curve the residuals are taken to be correlated as
#
#     corr(e(t), e(t')) = phi^|t - t'|          t in hours, 0 <= phi < 1
#
# so phi is the correlation of two residuals one hour apart -- the
# continuous-time AR(1), CAR(1) (Jones 1993; Pinheiro & Bates 2000, sec. 5.3,
# nlme::corCAR1). It handles uneven spacing and missing observations as they
# are: restricted to the times observed, the process is still first-order
# Markov with r_k = phi^(t_k - t_(k-1)) between neighbours.
#
# Its correlation matrix R has a bidiagonal inverse Cholesky factor W (Prais &
# Winsten 1954, for unequal gaps): along each curve in time order
#
#     e*_1 = e_1,    e*_k = (e_k - r_k e_(k-1)) / sqrt(1 - r_k^2)
#
# and W R W' = I exactly. So for a GIVEN phi the model with CAR(1) residuals IS
# the ordinary mixed model fitted to W y, W X and W Z with independent
# residuals (generalized least squares), and its log-likelihood -- ML or REML --
# is the whitened one plus log|det W|. That is why the fit stays an lme4 fit:
# Kenward-Roger (pbkrtest) and Satterthwaite (lmerTest) apply to it unchanged,
# conditional on phi, as they are conditional on tau. Checked against
# nlme::lme(correlation = corCAR1) in tests/testthat/test-traj-car1.R: same
# fixed effects, standard errors, sigma and log-likelihood.
dance_car1_whitener <- function(t, curve, phi) {
  n <- length(t)
  ord <- order(as.character(curve), t)
  cs <- as.character(curve)[ord]; ts <- t[ord]
  first <- c(TRUE, cs[-1L] != cs[-n])
  dt <- c(NA_real_, diff(ts))
  if (any(!first & dt <= 0))
    stop("Two observations of one curve share a time point; the AR(1) residual model needs distinct times.")
  r <- ifelse(first, 0, phi^dt)
  s <- sqrt(1 - r^2)
  nf <- which(!first)
  W <- Matrix::sparseMatrix(i = c(ord, ord[nf]), j = c(ord, ord[nf - 1L]),
                            x = c(1 / s, -r[nf] / s[nf]), dims = c(n, n))
  attr(W, "logdet") <- -sum(log(s))
  W
}

# log|det W| for a spec's data at phi: what turns the whitened log-likelihood
# into the likelihood of the data themselves
dance_car1_logdet <- function(spec, phi) {
  if (!length(phi) || !is.finite(phi) || phi <= 0) return(0)
  attr(dance_car1_whitener(spec$data$t, spec$data[[spec$curve_var %||% "subject"]], phi), "logdet")
}

# The mixed model with CAR(1) residuals at a given phi: lme4's own modular fit
# (lFormula, mkLmerDevfun, optimizeLmer, mkMerMod -- the steps lme4::lmer
# takes) on the whitened design. The object keeps the formula, the coefficient
# names and the model frame of the data, so predict(newdata), ranef(), emmeans
# and the tests read it like any other fit; fitted() and resid() are on the
# whitened scale (the residuals are the normalized ones). lmerTest's
# Satterthwaite object is built from the same whitened deviance function.
dance_traj_lmer_car1 <- function(f, d, phi, curve_var, REML, ctrl, light = FALSE) {
  lmod <- lme4::lFormula(f, data = d, REML = REML, control = ctrl)
  if (nrow(lmod$fr) != nrow(d)) stop("The model frame dropped rows; the whitening needs them all.")
  W <- dance_car1_whitener(d$t, d[[curve_var]], phi)
  att <- attributes(lmod$X)
  X <- as.matrix(W %*% lmod$X); attributes(X) <- att
  lmod$X <- X
  lmod$fr[[1L]] <- as.numeric(W %*% lmod$fr[[1L]])
  lmod$reTrms$Zt <- methods::as(lmod$reTrms$Zt %*% Matrix::t(W), "CsparseMatrix")
  devfun <- do.call(lme4::mkLmerDevfun, c(lmod, list(control = ctrl)))
  opt <- lme4::optimizeLmer(devfun, optimizer = ctrl$optimizer, restart_edge = ctrl$restart_edge,
                            boundary.tol = ctrl$boundary.tol, control = ctrl$optCtrl,
                            calc.derivs = ctrl$calc.derivs, use.last.params = ctrl$use.last.params)
  cc <- lme4::checkConv(attr(opt, "derivs"), opt$par, ctrl = ctrl$checkConv,
                        lbound = environment(devfun)$lower)
  m <- lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
                      mc = quote(dance_traj_lmer_car1()), lme4conv = cc)
  if (!isTRUE(light) && requireNamespace("lmerTest", quietly = TRUE)) {
    # a fresh deviance function on the same whitened design, so lmerTest's
    # Hessian and Jacobian do not move the fitted model's own state
    devfun2 <- do.call(lme4::mkLmerDevfun, c(lmod, list(control = ctrl)))
    as_lt <- get("as_lmerModLT", envir = asNamespace("lmerTest"))
    m <- tryCatch(as_lt(m, devfun2, tol = 1e-8), error = function(e) m)
  }
  m
}

# ------------------------------------------------------------------------------
# ONE FIT at a given rung
# ------------------------------------------------------------------------------
dance_traj_fit_one <- function(spec, re_formula, engine = "lmer", REML = TRUE,
                               cor_struct = NULL, fixed_formula = NULL, light = FALSE,
                               phi = NULL) {
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
  # CAR(1) residuals at a given phi: the whitened fit (phi = 0 IS independence)
  if (length(phi) == 1L && is.finite(phi) && phi > 0)
    return(tryCatch(dance_traj_lmer_car1(f, d, phi, spec$curve_var %||% "subject", REML,
                                         ctrl, light = light),
                    error = function(e) NULL))
  # `light` fits with lme4 alone: lmerTest's conversion computes what the
  # Satterthwaite tests need, which a fit read only for its log-likelihood (a
  # point of the tau profile) never uses -- same likelihood, half the time.
  fn <- if (!isTRUE(light) && requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer
        else lme4::lmer
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
                           ar1 = FALSE, max_rung = NULL, min_rung = 1L, light = FALSE,
                           phi = NULL) {
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
  # min_rung skips the richer structures already known not to converge -- the tau
  # profile fits every grid point on ONE structure, so that the likelihoods it
  # compares differ only in tau. The rungs keep their numbers in the full ladder.
  min_rung <- max(1L, as.integer(min_rung %||% 1L))
  attempts <- list()
  # the richer rungs min_rung skips are recorded as such, so the fit's log and
  # its simplified_reason still account for every rung above the one used
  for (i in seq_along(ladder)[seq_along(ladder) < min_rung])
    attempts[[i]] <- list(rung = i, label = ladder[[i]]$label, formula = ladder[[i]]$formula,
                          fitted = FALSE, converged = FALSE, singular = NA, skipped = TRUE,
                          summary = "not refitted: the structure was carried over from the fit that chose it")
  for (i in seq_along(ladder)[seq_along(ladder) >= min_rung]) {
    m <- dance_traj_fit_one(spec, ladder[[i]]$formula, engine, REML, cor_struct,
                            light = light, phi = phi)
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
        # CAR(1) residuals at this phi (0: independent); the profile that chose
        # phi is attached by the caller as $residual_ar
        phi = if (length(phi) == 1L && is.finite(phi)) phi else 0,
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
                                   flat_tol = 1.92, roles = NULL, progress = NULL,
                                   phi = 0, rung = NULL) {
  # phi: the CAR(1) residual correlation the profile is computed at (0 for
  # independent residuals). It is the same at every grid point, so log|det W|
  # is a constant of the profile and drops out of its comparisons.
  # rung: the random structure to profile on, when the caller already has one
  # (NULL: chosen at the middle of the grid, as below).
  if (is.null(tau_grid)) {
    span <- diff(range(d$t, na.rm = TRUE))
    tau_grid <- exp(seq(log(max(span / 50, 1e-3)), log(span * 3), length.out = 24))
  }
  spec_at <- function(tau) dance_traj_spec(d, period, n_harmonics, "exp_sat", design_terms,
                                           covariates, t0, tau, interaction = interaction,
                                           roles = roles)
  # ONE RANDOM STRUCTURE FOR THE WHOLE PROFILE. Each grid point used to walk the
  # random-effects ladder on its own, so a tau where the maximal structure failed
  # was fitted with a simpler one -- and its log-likelihood, from a different
  # model, was compared with its neighbours'. The structure is now chosen once,
  # at the middle of the grid, and every grid point is fitted on it; a point
  # where it does not converge is left out rather than refitted on another.
  if (is.null(rung)) {
    mid <- tau_grid[ceiling(length(tau_grid) / 2)]
    sp_mid <- spec_at(mid)
    ref <- if (isTRUE(sp_mid$ok)) dance_traj_fit(sp_mid, engine = engine, REML = FALSE,
                                                 light = TRUE, phi = phi) else NULL
    if (is.null(ref) || !isTRUE(ref$ok))
      return(list(ok = FALSE, message = paste(
        "The model did not fit at the middle of the tau grid, so no random structure",
        "could be chosen for the profile.", ref$message %||% sp_mid$message %||% "")))
    rung <- ref$re_rung
  }
  rows <- lapply(seq_along(tau_grid), function(j) {
    tau <- tau_grid[j]
    if (is.function(progress)) progress(j, length(tau_grid))
    sp <- spec_at(tau)
    if (!isTRUE(sp$ok)) return(NULL)
    fit <- dance_traj_fit(sp, engine = engine, REML = FALSE, min_rung = rung, max_rung = rung,
                          light = TRUE, phi = phi)
    if (!isTRUE(fit$ok)) return(NULL)
    data.frame(tau = tau, logLik = as.numeric(stats::logLik(fit$model)),
               rung = fit$re_rung, stringsAsFactors = FALSE)
  })
  prof <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(prof) || !nrow(prof))
    return(list(ok = FALSE, message = "No tau on the grid produced a usable fit."))
  # grid points whose fit did not converge on the chosen structure are left
  # out of the profile -- counted, so the report can say how many
  n_grid_failed <- length(tau_grid) - nrow(prof)

  best <- prof$tau[which.max(prof$logLik)]
  drop <- max(prof$logLik) - min(prof$logLik)
  flat <- drop < flat_tol
  at_edge <- best <= min(tau_grid) * 1.001 || best >= max(tau_grid) * 0.999

  # REFINED BETWEEN THE GRID'S NEIGHBOURS. Adjacent grid points are about 24%
  # apart, so the best grid point alone can miss the maximum by that much. A
  # one-dimensional search on log tau between the best point's neighbours, on
  # the same random structure, finds it to about 1%; the point is added to the
  # profile, so what is reported is still a point of the profile.
  ll_at <- function(tau) {
    sp <- spec_at(tau)
    if (!isTRUE(sp$ok)) return(NA_real_)
    f <- dance_traj_fit(sp, engine = engine, REML = FALSE, min_rung = rung, max_rung = rung,
                        light = TRUE, phi = phi)
    if (isTRUE(f$ok)) as.numeric(stats::logLik(f$model)) else NA_real_
  }
  if (!flat && !at_edge) {
    j <- match(best, prof$tau)
    lo_t <- prof$tau[max(1L, j - 1L)]; hi_t <- prof$tau[min(nrow(prof), j + 1L)]
    if (hi_t > lo_t) {
      if (is.function(progress)) progress(length(tau_grid), length(tau_grid))
      opt <- tryCatch(stats::optimize(function(lt) {
        v <- ll_at(exp(lt)); if (is.finite(v)) -v else 1e100
      }, interval = log(c(lo_t, hi_t)), tol = 0.01), error = function(e) NULL)
      if (!is.null(opt) && is.finite(opt$objective) && -opt$objective > max(prof$logLik)) {
        prof <- rbind(prof, data.frame(tau = exp(opt$minimum), logLik = -opt$objective,
                                       rung = rung, stringsAsFactors = FALSE))
        prof <- prof[order(prof$tau), , drop = FALSE]; rownames(prof) <- NULL
        best <- exp(opt$minimum)
      }
    }
  }
  # profile interval: where the log-likelihood falls 1.92 below its maximum (the
  # chi-square(1) 95% cutoff), interpolated linearly in log tau between the grid
  # points either side of the crossing. A side that never falls that far within
  # the grid is left open -- the interval is then a bound, and says so.
  cut <- max(prof$logLik) - flat_tol
  jb <- which.max(prof$logLik)
  cross <- function(idx) {
    for (j in idx) {
      a <- j; b <- if (idx[1] < jb) j + 1L else j - 1L
      if (prof$logLik[a] < cut && prof$logLik[b] >= cut) {
        w <- (cut - prof$logLik[a]) / (prof$logLik[b] - prof$logLik[a])
        return(exp(log(prof$tau[a]) + w * (log(prof$tau[b]) - log(prof$tau[a]))))
      }
    }
    NA_real_
  }
  ci_lo <- if (jb > 1L) cross(rev(seq_len(jb - 1L))) else NA_real_
  ci_hi <- if (jb < nrow(prof)) cross(seq.int(jb + 1L, nrow(prof))) else NA_real_
  inside <- c(if (is.finite(ci_lo)) ci_lo else min(prof$tau),
              if (is.finite(ci_hi)) ci_hi else max(prof$tau))
  ci_open <- c(lower = !is.finite(ci_lo), upper = !is.finite(ci_hi))

  final <- NULL
  # A FLAT PROFILE SUPPRESSES THE FIT, it does not merely annotate it. If the
  # likelihood cannot tell one tau from another, then "the fit at the best tau"
  # is the fit at an arbitrary point of a ridge, and handing it back as a
  # fitted model invites exactly the substantive inference the flatness rules
  # out. The profile itself is still returned, so the flatness is visible.
  if (!flat) {
    sp <- spec_at(best)
    if (isTRUE(sp$ok)) {
      final <- dance_traj_fit(sp, engine = engine, REML = REML, min_rung = rung, phi = phi)
      # if the profile's structure does not converge by REML at the optimum, the
      # ladder continues from there rather than reporting no fit
    }
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
       rung = rung, grid = tau_grid, ci_open = ci_open, n_grid_failed = n_grid_failed,
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
# PROFILE OVER phi, the CAR(1) residual correlation at a lag of one hour
# ------------------------------------------------------------------------------
# On ONE random structure (`rung`), as the tau profile is, so the likelihoods
# compared differ only in phi. REML when the fixed-effect design is fixed --
# the estimate nlme's corCAR1 gives -- and ML when tau is estimated alongside,
# since REML likelihoods of different trend shapes are not comparable. Every
# likelihood is the data's: the whitened one plus log|det W|.
#
# 20 grid points from 0 to 0.98, evenly spaced in -log(1 - phi) (denser where
# the likelihood turns fastest); the maximum refined between its neighbours;
# the 95% interval where the log-likelihood falls 1.92 below its maximum,
# interpolated on the same scale; and the likelihood-ratio test of phi = 0
# (independent residuals), whose null value lies on the boundary of the
# parameter space, so its p is half the chi-square(1) tail (Self & Liang 1987).
dance_traj_profile_phi <- function(spec, rung = 1L, REML = TRUE, grid = NULL,
                                   flat_tol = 1.92, progress = NULL) {
  if (!isTRUE(spec$ok)) return(list(ok = FALSE, message = spec$message))
  if (is.null(grid)) grid <- 1 - exp(-seq(0, -log(0.02), length.out = 20))
  u_of <- function(phi) -log(1 - phi)
  phi_of <- function(u) 1 - exp(-u)
  ll_at <- function(phi) {
    f <- dance_traj_fit(spec, REML = REML, min_rung = rung, max_rung = rung,
                        light = TRUE, phi = phi)
    if (!isTRUE(f$ok)) return(NA_real_)
    as.numeric(stats::logLik(f$model)) + dance_car1_logdet(spec, phi)
  }
  ll <- vapply(seq_along(grid), function(j) {
    if (is.function(progress)) progress(j, length(grid))
    ll_at(grid[j])
  }, numeric(1))
  if (!is.finite(ll[1L]))
    return(list(ok = FALSE, message = paste(
      "The model did not fit with independent residuals on the random structure chosen,",
      "so the residual autocorrelation could not be profiled.")))
  prof <- data.frame(phi = grid[is.finite(ll)], logLik = ll[is.finite(ll)])
  n_failed <- sum(!is.finite(ll))
  jb <- which.max(prof$logLik); best <- prof$phi[jb]
  at_upper <- jb == nrow(prof)
  if (!at_upper) {
    lo_u <- u_of(prof$phi[max(1L, jb - 1L)]); hi_u <- u_of(prof$phi[jb + 1L])
    opt <- tryCatch(stats::optimize(function(u) {
      v <- ll_at(phi_of(u)); if (is.finite(v)) -v else 1e100
    }, interval = c(lo_u, hi_u), tol = 0.001), error = function(e) NULL)
    if (!is.null(opt) && is.finite(opt$objective) && -opt$objective > max(prof$logLik)) {
      prof <- rbind(prof, data.frame(phi = phi_of(opt$minimum), logLik = -opt$objective))
      prof <- prof[order(prof$phi), , drop = FALSE]; rownames(prof) <- NULL
      best <- phi_of(opt$minimum)
    }
  }
  # the interval: where the profile falls 1.92 below its maximum; phi cannot
  # go below 0, so a profile still within 1.92 at phi = 0 has its interval
  # start there, and one still within 1.92 at 0.98 is open above
  top <- max(prof$logLik); cut <- top - flat_tol; jb <- which.max(prof$logLik)
  cross <- function(idx) {
    for (j in idx) {
      a <- j; b <- if (idx[1] < jb) j + 1L else j - 1L
      if (prof$logLik[a] < cut && prof$logLik[b] >= cut) {
        w <- (cut - prof$logLik[a]) / (prof$logLik[b] - prof$logLik[a])
        return(phi_of(u_of(prof$phi[a]) + w * (u_of(prof$phi[b]) - u_of(prof$phi[a]))))
      }
    }
    NA_real_
  }
  ci_lo <- if (jb > 1L) cross(rev(seq_len(jb - 1L))) else NA_real_
  ci_hi <- if (jb < nrow(prof)) cross(seq.int(jb + 1L, nrow(prof))) else NA_real_
  ll0 <- prof$logLik[prof$phi == 0][1]
  stat <- max(0, 2 * (top - ll0))
  list(ok = TRUE, phi = best,
       ci = c(if (is.finite(ci_lo)) ci_lo else 0, if (is.finite(ci_hi)) ci_hi else max(prof$phi)),
       ci_open = c(lower = FALSE, upper = !is.finite(ci_hi)),
       at_zero = best < 1e-3, at_upper = at_upper,
       lrt = list(statistic = stat, df = 1L,
                  p = if (stat > 1e-8) 0.5 * stats::pchisq(stat, 1, lower.tail = FALSE) else 1,
                  ll_max = top, ll0 = ll0,
                  method = "likelihood ratio against independent residuals (phi = 0), boundary-corrected"),
       profile = prof, grid = grid, rung = rung, REML = REML, n_grid_failed = n_failed)
}

# ------------------------------------------------------------------------------
# WHAT A FIT WITH CAR(1) RESIDUALS IS, ON THE SCALE OF THE DATA
# ------------------------------------------------------------------------------
# The whitened fit's own logLik(), fitted() and resid() are on the whitened
# scale. These put them back where a reader needs them.

# The log-likelihood of the data, with the parameters it counts: the whitened
# log-likelihood plus log|det W|, and phi as one parameter more.
dance_traj_loglik <- function(fit) {
  ll <- tryCatch(stats::logLik(fit$model), error = function(e) NULL)
  if (is.null(ll)) return(NULL)
  car1 <- isTRUE(fit$phi > 0)
  k <- as.integer(attr(ll, "df")) + as.integer(car1)
  v <- as.numeric(ll) + if (car1) dance_car1_logdet(fit$spec, fit$phi) else 0
  n <- tryCatch(stats::nobs(fit$model), error = function(e) nrow(fit$spec$data))
  list(logLik = v, df = k, n = n, REML = isTRUE(fit$REML),
       AIC = -2 * v + 2 * k, BIC = -2 * v + k * log(n))
}

# Fitted values and residuals on the data's scale (conditional on the random
# effects), and the NORMALIZED residuals: for a CAR(1) fit the whitened
# model's residuals -- independent, variance sigma^2, if the model is right --
# and for an independent-residual fit the residuals themselves.
dance_traj_residuals <- function(fit) {
  m <- fit$model; d <- fit$spec$data
  car1 <- isTRUE(fit$phi > 0)
  fv <- if (car1) as.numeric(stats::predict(m, newdata = d, re.form = NULL))
        else as.numeric(stats::fitted(m))
  list(fitted = fv, raw = d$y - fv, normalized = as.numeric(stats::resid(m)), car1 = car1)
}

# How strongly consecutive residuals of one curve go together: the
# correlation of each residual with the curve's previous one, pooled over
# curves. On the residuals of an independent-residual fit it says whether
# that assumption holds; on the normalized residuals of a CAR(1) fit, whether
# the AR(1) took the autocorrelation out.
dance_traj_resid_lag1 <- function(fit, which = c("normalized", "raw")) {
  which <- match.arg(which)
  rs <- dance_traj_residuals(fit)
  d <- fit$spec$data; cv <- as.character(d[[fit$spec$curve_var %||% "subject"]])
  ord <- order(cv, d$t)
  r <- rs[[which]][ord]; cs <- cv[ord]; ts <- d$t[ord]
  i <- which(c(FALSE, cs[-1L] == cs[-length(cs)]))
  if (length(i) < 3L) return(NULL)
  list(r = stats::cor(r[i], r[i - 1L]), n_pairs = length(i),
       median_gap = stats::median(ts[i] - ts[i - 1L]), which = which)
}

# The residual model in one sentence, for every place a fit is reported.
dance_traj_residual_statement <- function(fit) {
  ar <- fit$residual_ar
  if (!isTRUE(fit$phi > 0) && is.null(ar))
    return(paste("Residuals are taken as independent, with one variance: no autocorrelation",
                 "is modelled."))
  f3 <- function(x) format(signif(x, 3))
  if (is.null(ar))
    return(sprintf(paste("Residuals are autocorrelated within each curve as a continuous-time",
                         "AR(1), phi = %s per hour."), f3(fit$phi)))
  sprintf(paste(
    "Residuals are autocorrelated within each curve as a continuous-time AR(1): two",
    "residuals dt hours apart correlate phi^dt, with phi = %s per hour (95%% profile",
    "interval %s to %s%s; %s). Against independent residuals: likelihood ratio %s,",
    "p %s (boundary-corrected). The tests treat phi as known.%s"),
    f3(ar$phi %||% fit$phi), f3(ar$ci[1]), if (isTRUE(ar$ci_open[["upper"]])) "above " else "",
    f3(ar$ci[2]), ar$method %||% "REML", f3(ar$lrt$statistic),
    if (ar$lrt$p < 0.001) "< .001" else sprintf("= %.3f", ar$lrt$p),
    paste0(if (isTRUE(ar$at_zero)) " The estimate is 0: no positive autocorrelation, and the model is the independent-residual one." else "",
           if (isTRUE(ar$at_upper)) " The estimate sits at the top of the search (0.98): a bound, not an estimate." else "",
           if (!identical(ar$method %||% "REML", "REML")) paste(
             " The interval is computed at tau's estimate; the likelihood ratio compares with",
             "independent residuals at their own best tau.") else "",
           if (isTRUE(ar$n_grid_failed > 0L)) sprintf(
             " %d of the %d values of phi on the grid did not converge and are not in the profile.",
             as.integer(ar$n_grid_failed), length(ar$grid)) else "",
           if (identical(ar$settled, FALSE)) " tau and phi had not settled after four rounds." else ""))
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
  if (isTRUE(fit$phi > 0) || !is.null(fit$residual_ar))
    a("Residuals:         continuous-time AR(1) within each curve, phi = %s per hour",
      format(signif(fit$residual_ar$phi %||% fit$phi, 3)))
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
        if (isTRUE(at$skipped)) "not refitted (structure carried over)"
        else if (!at$fitted) "did not fit"
        else if (!at$converged) "did not converge"
        else if (at$singular) "singular"
        else "accepted")
  }
  paste(l, collapse = "\n")
}

# ------------------------------------------------------------------------------
# THE MIXED-EFFECTS COSINOR AS THE APP RUNS IT
# ------------------------------------------------------------------------------
# From a long frame to the fitted model the cosinor tab reports: tau held at a
# value, or profiled once for all cells within a range (a flat profile holds
# the held value, which must lie in the range), and on request A and tau per
# cell. It used to live inside the module's server code, so the exported
# script could not run it; it is one function now, called by the app and
# written verbatim into the script, so the two cannot fit different models.
#
#   tau_mode   "fixed" (hold tau_value) or "estimate" (profile in tau_lo..tau_hi)
#   progress   optional function(message, value in [0, 1])
dance_traj_run <- function(d, period = 24, n_harmonics = 1, trend = "none", roles = NULL,
                           t0 = NULL, tau_mode = c("fixed", "estimate"),
                           tau_value = NA_real_, tau_lo = NA_real_, tau_hi = NA_real_,
                           sat_cells = FALSE, residuals = c("independent", "car1"),
                           progress = NULL) {
  tau_mode <- match.arg(tau_mode)
  residuals <- match.arg(residuals)
  car1 <- identical(residuals, "car1")
  say <- function(msg, value) if (is.function(progress)) progress(msg, value)
  sat <- identical(trend, "exp_sat")
  value_ok <- length(tau_value) == 1L && isTRUE(is.finite(tau_value) && tau_value > 0)
  range_ok <- length(tau_lo) == 1L && length(tau_hi) == 1L &&
    isTRUE(all(is.finite(c(tau_lo, tau_hi))) && tau_lo > 0 && tau_hi > tau_lo)
  mode <- if (!sat) "none" else if (identical(tau_mode, "estimate")) "profile" else "fixed"
  if (sat && identical(mode, "fixed") && !value_ok)
    return(list(ok = FALSE, message = paste(
      "A saturating trend needs a time constant. Set tau under Model",
      "Specification, or choose to estimate it: the mixed model is linear in its",
      "coefficients only for a given tau.")))
  if (sat && (identical(mode, "profile") || isTRUE(sat_cells)) && !range_ok)
    return(list(ok = FALSE, message = paste(
      "Estimating tau needs the range to search: a lower and an upper bound in",
      "hours under Model Specification.")))
  # a flat profile holds tau at the held value, so that value has to be one tau
  # is allowed to take
  if (sat && identical(mode, "profile") && !(value_ok && tau_value >= tau_lo && tau_value <= tau_hi))
    return(list(ok = FALSE, message = paste(
      "The held value of tau is what a flat profile falls back to, so it must be",
      "a positive number inside the search range.")))
  if (is.null(t0)) t0 <- min(d$t, na.rm = TRUE)

  fit_at <- function(tau) {
    sp <- dance_traj_spec(d, period = period, n_harmonics = n_harmonics, trend = trend,
                          t0 = t0, tau = tau, roles = roles)
    if (!isTRUE(sp$ok)) return(list(ok = FALSE, message = sp$message))
    dance_traj_fit(sp)
  }
  if (identical(mode, "profile")) {
    # the profile searches the range only: 24 points spaced evenly in log tau
    # between its ends, so tau cannot come back at a value that was ruled out
    grid <- exp(seq(log(tau_lo), log(tau_hi), length.out = 24))
    pr <- dance_traj_profile_tau(d, period, n_harmonics, design_terms = NULL, t0 = t0,
                                 roles = roles, tau_grid = grid,
                                 progress = function(j, n) say(
                                   sprintf("Profiling the likelihood over tau: grid point %d of %d", j, n),
                                   0.8 * j / (n + 1)))
    if (!isTRUE(pr$ok)) return(list(ok = FALSE, message = pr$message))
    # the profile without its fitted model, which IS the fit below
    prof <- pr[c("profile", "tau", "ci", "ci_open", "flat", "at_edge", "logLik_range",
                 "message", "rung", "grid", "n_grid_failed")]
    if (!is.null(pr$fit) && isTRUE(pr$fit$ok)) {
      ff <- pr$fit
      ff$tau_mode <- "profile"
    } else {
      # A FLAT PROFILE: the data do not identify tau, so the profile returns no
      # fit at its argmax. The model is fitted at the held value instead, and
      # every report says so.
      say("Fitting the mixed-effects cosinor at the held tau (flat profile)...", 0.85)
      ff <- fit_at(tau_value)
      if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
      ff$tau_mode <- "fixed_after_flat"
    }
    ff$tau_profile <- prof
  } else {
    say("Fitting the mixed-effects cosinor...", 0.5)
    ff <- fit_at(if (sat) tau_value else NULL)
    if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
    ff$tau_mode <- mode
  }

  # AUTOCORRELATED RESIDUALS (CAR(1) within each curve). The random structure
  # is the one the fit above reached with independent residuals; phi is
  # profiled on it, and the model refitted by REML at phi-hat.
  #   tau held (or no saturating trend, or a flat tau profile): phi by REML.
  #   tau estimated: tau and phi by ML, each profiled given the other in turn
  #   until phi moves by less than 0.005 (at most four rounds).
  if (car1) {
    tau_fields <- c("tau_mode", "tau_profile", "tau_estimated", "tau_profile_ci", "tau_warning")
    joint <- identical(ff$tau_mode, "profile")
    rung <- if (joint) ff$tau_profile$rung %||% ff$re_rung else ff$re_rung
    prof_indep <- if (joint) ff$tau_profile else NULL   # the tau profile at phi = 0
    phi_used <- 0; rounds <- 0L; settled <- TRUE
    repeat {
      rounds <- rounds + 1L
      say(sprintf("Estimating the residual autocorrelation%s...",
                  if (joint) sprintf(" (round %d)", rounds) else ""), 0.86)
      php <- dance_traj_profile_phi(ff$spec, rung = rung, REML = !joint)
      if (!isTRUE(php$ok)) return(list(ok = FALSE, message = php$message))
      if (!joint) break
      if (abs(php$phi - phi_used) < 0.005) break
      if (rounds >= 4L) { settled <- FALSE; break }
      phi_used <- php$phi
      say(sprintf("Profiling tau at the residual autocorrelation (round %d)...", rounds), 0.88)
      pr <- dance_traj_profile_tau(d, period, n_harmonics, design_terms = NULL, t0 = t0,
                                   roles = roles, tau_grid = grid, phi = phi_used)
      if (!isTRUE(pr$ok)) return(list(ok = FALSE, message = pr$message))
      prof <- pr[c("profile", "tau", "ci", "ci_open", "flat", "at_edge", "logLik_range",
                   "message", "rung", "grid", "n_grid_failed")]
      if (is.null(pr$fit) || !isTRUE(pr$fit$ok)) {
        # flat at this phi: tau is held at the held value, and phi is
        # estimated given it by REML, as for any held tau
        ff <- fit_at(tau_value)
        if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
        ff$tau_mode <- "fixed_after_flat"; ff$tau_profile <- prof
        joint <- FALSE; rung <- ff$re_rung; next
      }
      ff <- pr$fit; ff$tau_mode <- "profile"; ff$tau_profile <- prof
      rung <- pr$rung
    }
    if (joint) {
      # THE LIKELIHOOD RATIO WITH tau ESTIMATED: independence is the model at
      # phi = 0 with tau at ITS best value, not at the tau the AR(1) chose. The
      # first tau profile is exactly that (ML, phi = 0); profiled again on this
      # random structure if the rounds moved to another one.
      p0 <- if (identical(as.integer(prof_indep$rung), as.integer(rung))) prof_indep
            else dance_traj_profile_tau(d, period, n_harmonics, design_terms = NULL, t0 = t0,
                                        roles = roles, tau_grid = grid, phi = 0, rung = rung)
      ll0 <- max(c(php$lrt$ll0, if (isTRUE(p0$ok %||% TRUE)) p0$profile$logLik), na.rm = TRUE)
      stat <- max(0, 2 * (php$lrt$ll_max - ll0))
      php$lrt <- utils::modifyList(php$lrt, list(
        ll0 = ll0, statistic = stat,
        p = if (stat > 1e-8) 0.5 * stats::pchisq(stat, 1, lower.tail = FALSE) else 1,
        method = paste("likelihood ratio against independent residuals (phi = 0) with tau",
                       "re-estimated, boundary-corrected")))
    }
    say("Refitting at the residual autocorrelation...", 0.9)
    ff2 <- dance_traj_fit(ff$spec, REML = TRUE, min_rung = rung, phi = php$phi)
    if (!isTRUE(ff2$ok)) return(list(ok = FALSE, message = ff2$message))
    for (k in tau_fields) ff2[[k]] <- ff[[k]]
    ff <- ff2
    php$method <- if (joint) "ML, jointly with tau" else "REML"
    php$rounds <- rounds; php$settled <- settled
    ff$residual_ar <- php
  }
  # A AND TAU IN EVERY CELL, on request: a nonlinear mixed model (nlme), the
  # only way to compare tau itself, and the initial slope A/tau, between cells
  if (sat && isTRUE(sat_cells)) {
    say("Estimating A and tau in every cell (nonlinear mixed model)...", 0.9)
    ff$sat <- dance_traj_sat_fit(ff, tau_range = c(tau_lo, tau_hi),
                                 progress = function(msg) say(msg, 0.95), phi = ff$phi %||% 0)
  }
  ff$tau_range <- if (sat && range_ok) c(tau_lo, tau_hi) else NULL
  ff
}
