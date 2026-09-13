# ==============================================================================
# server/08d_helpers_traj.R — LAYER A: trajectory model SPECIFICATION
# ==============================================================================
# P21 phase 2. This file fits nothing. It turns "these columns, this period,
# this many harmonics, this trend" into an immutable object that says exactly
# what model is about to be fitted, and it is pure, so every rule in it is
# testable without a Shiny session.
#
# WHY A SPECIFICATION LAYER AT ALL
# --------------------------------
# The module it replaces built its model inside a 1,251-line observer, and the
# design it could express was fixed at write time: one grouping variable in the
# Harmonic Regression tab (finding A1), exactly one between x one within factor
# in dance_mixed_cosinor() (finding A6). Neither could be asked for a three-group
# design, a covariate, or a three-way interaction, because the formula was
# spelled out in the source rather than derived from the data.
#
# Everything here exists to make that derivation explicit:
#
#   dance_traj_long()     wide curves -> long frame with arbitrary factors
#   dance_traj_classify() reads between/within/mixed OFF THE DATA (brief §17)
#   dance_traj_basis()    builds the trend and Fourier columns
#   dance_traj_formula()  composes basis x design into a fixed-effects formula
#   dance_traj_re_ladder() the random-effects hierarchy of brief §4
#   dance_traj_spec()     the immutable object, with $ok and a refusal message
#
# THE ONE IDEA THAT MAKES IT GENERAL
# ----------------------------------
# A trajectory is a set of basis columns: an intercept, zero or more trend
# columns, and a (cos, sin) pair per harmonic. "Do the groups differ in their
# trajectories" is then "does the design interact with that basis", which in a
# formula is
#
#     y ~ (trend + c1 + s1 + ... ) * f1 * f2 * ... + (basis | subject)
#
# and nothing in it is specific to a 2 x 2. Three groups, four conditions, a
# covariate, a three-way interaction: same call, different term list.
#
# CONDITIONAL LINEARITY, which is why the saturating trend is not special.
# A_sat * (1 - exp(-(t - t0)/tau)) is nonlinear in tau and LINEAR IN A_sat.
# Hold tau fixed and the column 1 - exp(-(t - t0)/tau) is just another fixed
# basis column, so the whole model is linear again. That is what lets layer B
# profile over tau instead of throwing the fit at a general nonlinear optimiser,
# which is where the old module's convergence failures came from.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

DANCE_TRAJ_TRENDS <- c("none", "linear", "log", "exp_sat")

# ------------------------------------------------------------------------------
# Syntactic column names, with the human label kept alongside
# ------------------------------------------------------------------------------
# A factor called "CTI flex" or a level called "sleep deprived" cannot appear in
# a formula unquoted. Renaming here, once, is safer than backticking at every
# call site and forgetting one.
dance_traj_safe_name <- function(x) {
  s <- make.names(as.character(x))
  s <- gsub("[^A-Za-z0-9_.]", "_", s)
  sub("^([^A-Za-z.])", "v\\1", s)
}

# ------------------------------------------------------------------------------
# LONG FRAME with an arbitrary number of design factors
# ------------------------------------------------------------------------------
# `curves` is subjects x time. `factors` is a NAMED LIST of per-row vectors --
# one entry per design factor, each as long as nrow(curves). Covariates that are
# numeric are carried through unconverted; anything else becomes a factor.
dance_traj_long <- function(curves, time_points, subject, factors = list()) {
  curves <- as.matrix(curves)
  n <- nrow(curves); p <- ncol(curves)
  stopifnot(length(time_points) == p, length(subject) == n)
  if (length(factors) && is.null(names(factors)))
    stop("dance_traj_long(): `factors` must be a NAMED list.")
  for (nm in names(factors))
    if (length(factors[[nm]]) != n)
      stop(sprintf("dance_traj_long(): factor '%s' has %d values for %d curves.",
                   nm, length(factors[[nm]]), n))

  d <- data.frame(
    .row    = rep(seq_len(n), each = p),
    subject = factor(rep(as.character(subject), each = p)),
    t       = rep(time_points, times = n),
    y       = as.vector(t(curves)),
    stringsAsFactors = FALSE
  )
  labels <- character(0)
  for (nm in names(factors)) {
    col <- dance_traj_safe_name(nm)
    v <- factors[[nm]]
    d[[col]] <- if (is.numeric(v) && length(unique(v[!is.na(v)])) > 12)
      rep(v, each = p) else factor(rep(as.character(v), each = p))
    labels[[col]] <- nm
  }
  d <- d[is.finite(d$y), , drop = FALSE]
  d$subject <- droplevels(d$subject)
  for (col in names(labels))
    if (is.factor(d[[col]])) d[[col]] <- droplevels(d[[col]])
  attr(d, "labels") <- labels
  d
}

# ------------------------------------------------------------------------------
# CLASSIFY each factor from the data (brief §17)
# ------------------------------------------------------------------------------
# The user names their factors; the app works out what kind they are, by
# counting distinct levels within a participant:
#
#   every participant sees exactly one level   -> between-participant
#   some participant sees more than one        -> within-participant (repeated)
#
# A factor that is within for SOME participants and between for others is
# neither, and is reported as "partial" rather than guessed at -- it is usually
# a data error (an inconsistent identifier) and silently treating it as within
# would put unpaired observations into a paired term.
dance_traj_classify <- function(d, factors = NULL) {
  if (is.null(factors))
    factors <- setdiff(names(d), c(".row", "subject", "t", "y"))
  labels <- attr(d, "labels") %||% stats::setNames(factors, factors)
  out <- do.call(rbind, lapply(factors, function(f) {
    v <- d[[f]]
    if (is.numeric(v) && !is.factor(v))
      return(data.frame(factor = f, label = labels[[f]] %||% f, n_levels = NA_integer_,
                        role = "covariate", min_within = NA_integer_,
                        max_within = NA_integer_, note = "numeric covariate",
                        stringsAsFactors = FALSE))
    per <- tapply(as.character(v), d$subject, function(x) length(unique(x[!is.na(x)])))
    per <- per[is.finite(per) & per > 0]
    mn <- if (length(per)) min(per) else NA_integer_
    mx <- if (length(per)) max(per) else NA_integer_
    role <- if (!is.finite(mx)) "unusable"
            else if (mx == 1) "between"
            else if (mn > 1) "within"
            else "partial"
    note <- switch(role,
      between = "one level per participant",
      within  = sprintf("every participant contributes %s levels",
                        if (mn == mx) as.character(mn) else sprintf("%d-%d", mn, mx)),
      partial = sprintf("%d participant(s) contribute one level and %d contribute more -- this is not a clean within factor",
                        sum(per == 1), sum(per > 1)),
      "no usable values")
    data.frame(factor = f, label = labels[[f]] %||% f,
               n_levels = nlevels(droplevels(as.factor(v))),
               role = role, min_within = mn, max_within = mx, note = note,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# The one-word description of the whole design, for the readout and the report.
dance_traj_design_kind <- function(cls) {
  r <- cls$role[cls$role %in% c("between", "within")]
  if (!length(r)) "single-group"
  else if (all(r == "between")) "between-participant"
  else if (all(r == "within")) "within-participant (repeated measures)"
  else "mixed (between x within)"
}

# ------------------------------------------------------------------------------
# BASIS COLUMNS
# ------------------------------------------------------------------------------
# Adds the trend column(s) and the Fourier pairs, and returns the term names so
# the formula builder never has to guess them. `t0` is the reference time: the
# intercept is the fitted level AT t0, and brief §3 requires that to be stated
# rather than assumed.
#
# `tau` is used only by trend = "exp_sat", where the column is
# 1 - exp(-(t - t0)/tau): conditional on tau the model is linear, which is what
# layer B profiles over.
dance_traj_basis <- function(d, period = 24, n_harmonics = 1, trend = "none",
                             t0 = NULL, tau = NULL) {
  trend <- match.arg(trend, DANCE_TRAJ_TRENDS)
  n_harmonics <- max(1L, as.integer(n_harmonics))
  if (is.null(t0)) t0 <- min(d$t, na.rm = TRUE)
  tt <- d$t - t0

  trend_terms <- character(0)
  if (trend == "linear") {
    d$trend_lin <- tt
    trend_terms <- "trend_lin"
  } else if (trend == "log") {
    d$trend_log <- log(tt + 1)
    trend_terms <- "trend_log"
  } else if (trend == "exp_sat") {
    if (is.null(tau) || !is.finite(tau) || tau <= 0)
      stop("dance_traj_basis(): trend = 'exp_sat' needs a positive `tau`; layer B supplies one per profile grid point.")
    d$trend_sat <- 1 - exp(-tt / tau)
    trend_terms <- "trend_sat"
  }

  harm_terms <- character(0)
  for (h in seq_len(n_harmonics)) {
    w <- 2 * pi * h / period
    d[[paste0("c", h)]] <- cos(w * tt)
    d[[paste0("s", h)]] <- sin(w * tt)
    harm_terms <- c(harm_terms, paste0("c", h), paste0("s", h))
  }
  list(data = d, trend_terms = trend_terms, harm_terms = harm_terms,
       basis_terms = c(trend_terms, harm_terms),
       t0 = t0, tau = tau, period = period, n_harmonics = n_harmonics,
       trend = trend)
}

# ------------------------------------------------------------------------------
# FIXED-EFFECTS FORMULA: basis x design
# ------------------------------------------------------------------------------
# `design_terms` are the factors to cross with the basis; `covariates` are
# additive and are NOT crossed with it unless the caller asks, because crossing
# a continuous covariate with a full harmonic block is rarely what anyone means
# and costs 2K parameters per covariate.
dance_traj_formula <- function(basis_terms, design_terms = character(0),
                               covariates = character(0), response = "y",
                               interaction = c("full", "additive")) {
  interaction <- match.arg(interaction)
  bas <- if (length(basis_terms)) paste0("(", paste(basis_terms, collapse = " + "), ")") else NULL
  if (!length(design_terms)) {
    rhs <- if (is.null(bas)) "1" else bas
  } else {
    des <- paste(design_terms, collapse = if (interaction == "full") " * " else " + ")
    if (length(design_terms) > 1 && interaction == "additive") des <- paste0("(", des, ")")
    rhs <- if (is.null(bas)) des else paste(bas, "*", des)
  }
  if (length(covariates)) rhs <- paste(rhs, "+", paste(covariates, collapse = " + "))
  paste(response, "~", rhs)
}

# ------------------------------------------------------------------------------
# RANDOM-EFFECTS LADDER (brief §4)
# ------------------------------------------------------------------------------
# Tried top-down; layer B stops at the first rung that converges AND is not
# singular, and reports which rung that was. Brief §4: "Do not automatically fit
# an unnecessarily maximal random-effects structure" and "Never silently change
# the requested statistical model" -- so the ladder is short, principled and
# reported, rather than an automatic search.
#
# The order drops, in turn: the participant-specific trend, the correlations
# among the participant-specific rhythm terms, the rhythm itself, and finally
# everything but the level.
dance_traj_re_ladder <- function(basis_terms, harm_terms, within_terms = character(0),
                                 group = "subject") {
  rung <- function(terms, label, corr = TRUE) {
    inner <- if (!length(terms)) "1" else paste(c("1", terms), collapse = " + ")
    list(formula = sprintf("(%s %s %s)", inner, if (corr) "|" else "||", group),
         terms = terms, label = label, correlated = corr)
  }
  out <- list()
  if (length(within_terms))
    out[[length(out) + 1L]] <- rung(c(basis_terms, within_terms),
                                    "participant trend, rhythm and within-factor effect")
  out[[length(out) + 1L]] <- rung(basis_terms, "participant trend and rhythm")
  if (!identical(basis_terms, harm_terms))
    out[[length(out) + 1L]] <- rung(harm_terms, "participant rhythm")
  out[[length(out) + 1L]] <- rung(harm_terms, "participant rhythm, uncorrelated", corr = FALSE)
  out[[length(out) + 1L]] <- rung(character(0), "participant level only")
  # de-duplicate structurally identical rungs (e.g. trend = "none")
  seen <- character(0); keep <- logical(length(out))
  for (i in seq_along(out)) {
    keep[i] <- !(out[[i]]$formula %in% seen)
    seen <- c(seen, out[[i]]$formula)
  }
  out[keep]
}

# ------------------------------------------------------------------------------
# THE DESIGN TABLE (brief §16): observations and participants per cell
# ------------------------------------------------------------------------------
dance_traj_cells <- function(d, design_terms) {
  if (!length(design_terms))
    return(data.frame(cell = "(all)", n_participants = nlevels(droplevels(d$subject)),
                      n_obs = nrow(d), stringsAsFactors = FALSE))
  key <- interaction(lapply(design_terms, function(f) d[[f]]), sep = " x ", drop = TRUE)
  sp <- split(seq_len(nrow(d)), key)
  data.frame(
    cell = names(sp),
    n_participants = vapply(sp, function(i) length(unique(as.character(d$subject[i]))), integer(1)),
    n_obs = vapply(sp, length, integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# THE SPEC
# ------------------------------------------------------------------------------
# Everything a reader needs to reproduce the model, and everything layer B needs
# to fit it. Returns $ok = FALSE with a message rather than a half-built object
# when the design cannot support the requested model.
dance_traj_spec <- function(d, period = 24, n_harmonics = 1, trend = "none",
                            design_terms = NULL, covariates = character(0),
                            t0 = NULL, tau = NULL, response = "y",
                            interaction = c("full", "additive"),
                            time_units = "hours") {
  interaction <- match.arg(interaction)
  bad <- function(msg) list(ok = FALSE, message = msg)

  if (!all(c("subject", "t", "y") %in% names(d)))
    return(bad("The long frame needs subject, t and y columns; build it with dance_traj_long()."))
  if (!nrow(d)) return(bad("No usable observations."))

  cls <- dance_traj_classify(d)
  # When the caller does not name the design terms, every classifiable factor in
  # the frame is one. A "partial" factor is therefore NOT quietly left out of
  # that list: it was put in the frame to be modelled, and dropping it silently
  # would answer a different question from the one asked. It is carried into
  # design_terms precisely so the refusal below fires.
  if (is.null(design_terms))
    design_terms <- cls$factor[cls$role %in% c("between", "within", "partial")]
  covariates <- union(covariates, cls$factor[cls$role == "covariate"])
  covariates <- setdiff(covariates, design_terms)

  partial <- cls$factor[cls$role == "partial"]
  if (length(intersect(partial, design_terms)))
    return(bad(sprintf(paste(
      "'%s' is neither a between- nor a within-participant factor: %s.",
      "Fix the participant identifiers, or drop it, before fitting -- guessing",
      "would put unpaired observations into a paired term."),
      paste(intersect(partial, design_terms), collapse = "', '"),
      cls$note[match(intersect(partial, design_terms)[1], cls$factor)])))

  bas <- tryCatch(dance_traj_basis(d, period, n_harmonics, trend, t0, tau),
                  error = function(e) conditionMessage(e))
  if (is.character(bas)) return(bad(bas))

  # ---- identifiability -------------------------------------------------------
  # TWO separate counts, and the second is the one that is easy to get wrong.
  #
  # Observations per cell must exceed the parameters a cell-specific trajectory
  # needs -- obvious, and necessary.
  #
  # DISTINCT TIME POINTS must ALSO exceed them, and that is the binding
  # constraint in this app's data. Every basis column is a function of t alone,
  # so a cell measured at 4 distinct times spans a 4-dimensional space no matter
  # how many participants contribute: asking for 3 harmonics plus a trend there
  # is asking for 8 columns of a rank-4 space. The first draft of this check
  # counted only observations and accepted exactly that model, because 4 time
  # points x 3 participants is 12 observations and 12 > 8. The rank is what
  # matters, and the rank is set by the distinct times.
  cells <- dance_traj_cells(bas$data, design_terms)
  per_cell_par <- 1L + length(bas$basis_terms)
  thin <- cells[cells$n_obs < per_cell_par + 1L, , drop = FALSE]
  if (nrow(thin))
    return(bad(sprintf(paste(
      "%d design cell(s) have fewer observations than the %d parameters a cell-specific",
      "trajectory needs (%s). Reduce the harmonics, simplify the trend, or drop the cell."),
      nrow(thin), per_cell_par,
      paste(sprintf("%s: %d obs", thin$cell, thin$n_obs), collapse = "; "))))

  n_times <- length(unique(bas$data$t))
  if (n_times < per_cell_par)
    return(bad(sprintf(paste(
      "The basis has %d columns (intercept + %d trend + %d harmonic) but the data",
      "carry only %d distinct time points. Every basis column is a function of",
      "time alone, so the design matrix is rank deficient however many",
      "participants contribute. Fit at most %d harmonic(s) on this time grid, or",
      "simplify the trend."),
      per_cell_par, length(bas$trend_terms), length(bas$harm_terms), n_times,
      max(0L, (n_times - 1L - length(bas$trend_terms)) %/% 2L))))

  within_terms <- intersect(cls$factor[cls$role == "within"], design_terms)
  ladder <- dance_traj_re_ladder(bas$basis_terms, bas$harm_terms, within_terms)
  fixed <- dance_traj_formula(bas$basis_terms, design_terms, covariates,
                              response, interaction)

  structure(list(
    ok = TRUE,
    data = bas$data,
    response = response,
    subject = "subject",
    time_var = "t", time_units = time_units,
    t0 = bas$t0,
    period = period, n_harmonics = n_harmonics,
    trend = trend, tau = bas$tau,
    nonlinear = identical(trend, "exp_sat"),
    basis_terms = bas$basis_terms, harm_terms = bas$harm_terms,
    trend_terms = bas$trend_terms,
    design_terms = design_terms, covariates = covariates,
    within_terms = within_terms,
    between_terms = intersect(cls$factor[cls$role == "between"], design_terms),
    interaction = interaction,
    classification = cls,
    design_kind = dance_traj_design_kind(cls),
    cells = cells,
    n_participants = nlevels(droplevels(bas$data$subject)),
    n_obs = nrow(bas$data),
    fixed_formula = fixed,
    re_ladder = ladder,
    formula = paste(fixed, "+", ladder[[1]]$formula)
  ), class = "dance_traj_spec")
}

# The Methods paragraph, straight off the spec (brief §16). One function, so the
# screen, the report and the exported script cannot describe the same fit
# differently.
dance_traj_describe <- function(spec) {
  if (!isTRUE(spec$ok)) return(spec$message)
  l <- character(0)
  a <- function(...) l <<- c(l, sprintf(...))
  a("Design:            %s", spec$design_kind)
  a("Response:          %s", spec$response)
  a("Participants:      %d, contributing %d observations", spec$n_participants, spec$n_obs)
  a("Time variable:     %s (%s); reference time t = 0 at %s",
    spec$time_var, spec$time_units, format(spec$t0))
  a("Period:            %s %s; %d harmonic(s)", format(spec$period), spec$time_units,
    spec$n_harmonics)
  a("Trend basis:       %s%s", spec$trend,
    if (isTRUE(spec$nonlinear)) sprintf(" (tau = %s, profiled)", format(spec$tau)) else "")
  if (length(spec$design_terms)) {
    for (i in seq_len(nrow(spec$classification))) {
      cc <- spec$classification[i, ]
      if (!cc$factor %in% spec$design_terms) next
      a("  %-16s %s, %d levels -- %s", cc$label, cc$role, cc$n_levels, cc$note)
    }
  } else a("  (no design factors: one trajectory for the whole sample)")
  if (length(spec$covariates)) a("Covariates:        %s", paste(spec$covariates, collapse = ", "))
  a("Fixed effects:     %s", spec$fixed_formula)
  a("Random effects:    %s  [rung 1 of %d; layer B reports which was reached]",
    spec$re_ladder[[1]]$formula, length(spec$re_ladder))
  paste(l, collapse = "\n")
}
