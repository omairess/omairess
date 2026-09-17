# ==============================================================================
# server/08g_helpers_trajband.R — LAYER E: trajectories, bands and contrasts
# ==============================================================================
# P21 phase 3. Everything here is a LINEAR function of the fixed effects of the
# fit made in layer B, so everything here has an exact covariance and none of it
# needs a bootstrap.
#
#   dance_traj_cell_grid()   the design cells, as a data frame
#   dance_traj_predict()     one fitted trajectory per cell, with a band
#   dance_traj_diff_curve()  the difference between two cells' trajectories
#   dance_traj_contrasts()   pairwise contrasts of level, amplitude and phase
#   dance_traj_group_fits()  the cell-keyed adapter the existing plots consume
#
# POINTWISE VS SIMULTANEOUS, WHICH IS THE WHOLE POINT OF A BAND
# ------------------------------------------------------------
# A pointwise 95% interval is correct at each t SEPARATELY. Read across a whole
# curve -- "the two trajectories differ, look, the bands come apart around
# 14:00" -- it is not a 95% statement at all, because the reader has made an
# inference about the curve after looking at every point of it. That is the
# error this module exists to stop being made silently, so a band is labelled
# with which of the two it is, and the simultaneous one is the default for a
# difference curve.
#
# The simultaneous band is Scheffe's, widened from z to sqrt(q * F(q, ddf)),
# where q is the rank of the basis block being plotted. It is exact for the
# whole q-dimensional linear family the curve lives in -- every trajectory that
# the basis can express -- and is therefore conservative for the single curve
# actually drawn. That conservatism is stated rather than tuned away: the
# alternative is a simulation-based band (multivariate-t maximum), which is
# tighter, needs a random draw, and is not reproducible without storing a seed.
#
# WHY THE DIFFERENCE CURVE IS NOT TWO CURVES SUBTRACTED
# -----------------------------------------------------
# Cell trajectories are built from shared coefficients, so their errors are
# correlated -- often strongly, and often negatively. Subtracting two curves and
# combining their bands as if independent is wrong in both directions. The
# difference is formed as ONE linear contrast, (x_i - x_j)' beta, and its
# variance is that contrast's own quadratic form. It is common for a difference
# band to exclude zero where the two individual bands overlap; that is not a
# paradox, it is what the covariance was carrying.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# The design cells, in the order every other function here uses
# ------------------------------------------------------------------------------
dance_traj_cell_grid <- function(spec) {
  dt <- spec$design_terms
  if (!length(dt))
    return(data.frame(.cell = "(all)", stringsAsFactors = FALSE))
  lv <- lapply(dt, function(f) levels(spec$data[[f]]))
  names(lv) <- dt
  g <- expand.grid(lv, stringsAsFactors = TRUE, KEEP.OUT.ATTRS = FALSE)
  g$.cell <- do.call(paste, c(lapply(dt, function(f) as.character(g[[f]])),
                              list(sep = " x ")))
  g
}

# ------------------------------------------------------------------------------
# The fixed-effect design row for one cell at a vector of times
# ------------------------------------------------------------------------------
# Built through model.matrix on the SAME terms object as the fit, so a factor
# contrast, a covariate, or an interaction the user added cannot fall out of
# step with it. Covariates are held at their mean, which is what makes the
# result a trajectory "for an average covariate value" rather than for whichever
# row happened to be first -- and that is said in the returned object.
dance_traj_design_rows <- function(fit, cell_row, times) {
  spec <- fit$spec
  d <- spec$data
  nd <- data.frame(t = times)
  # basis columns, recomputed at the requested times on the SAME reference
  tt <- times - spec$t0
  for (nm in spec$trend_terms) nd[[nm]] <- switch(nm,
    trend_lin = tt,
    trend_log = log1p(pmax(0, tt)),
    trend_sat = 1 - exp(-tt / spec$tau),
    stop(sprintf("unknown trend column '%s'", nm)))
  for (h in seq_len(spec$n_harmonics)) {
    w <- 2 * pi * h * tt / spec$period
    nd[[paste0("c", h)]] <- cos(w); nd[[paste0("s", h)]] <- sin(w)
  }
  for (f in spec$design_terms)
    nd[[f]] <- factor(rep(as.character(cell_row[[f]]), length(times)),
                      levels = levels(d[[f]]))
  for (f in spec$covariates)
    nd[[f]] <- if (is.numeric(d[[f]])) rep(mean(d[[f]], na.rm = TRUE), length(times))
               else factor(rep(levels(d[[f]])[1], length(times)), levels = levels(d[[f]]))
  tm <- stats::delete.response(stats::terms(stats::as.formula(spec$fixed_formula), data = d))
  stats::model.matrix(tm, data = nd, contrasts.arg = NULL)
}

# Fixed-effect coefficients and their covariance, engine-agnostic.
dance_traj_beta <- function(fit) {
  m <- fit$model
  b <- if (inherits(m, "glmmTMB")) lme4::fixef(m)$cond else lme4::fixef(m)
  V <- if (inherits(m, "glmmTMB")) as.matrix(stats::vcov(m)$cond)
       else as.matrix(stats::vcov(m))
  list(beta = b, V = V)
}

# The multiplier for a band: z, or Scheffe's sqrt(q F) for a simultaneous one.
dance_traj_band_mult <- function(band, conf, q, ddf) {
  if (identical(band, "pointwise")) return(stats::qnorm(1 - (1 - conf) / 2))
  if (!is.finite(ddf) || ddf <= 0) ddf <- Inf
  sqrt(q * stats::qf(conf, q, ddf))
}

# Residual degrees of freedom for the Scheffe multiplier. lmerTest gives a
# Satterthwaite df per coefficient; the SMALLEST is taken, which widens the band
# rather than narrowing it. With nothing available the band goes to the
# chi-square limit and says so.
dance_traj_ddf <- function(fit) {
  m <- fit$model
  if (inherits(m, "lmerModLmerTest") && requireNamespace("lmerTest", quietly = TRUE)) {
    s <- tryCatch(stats::coef(summary(m)), error = function(e) NULL)
    if (!is.null(s) && "df" %in% colnames(s)) {
      d <- suppressWarnings(min(s[, "df"], na.rm = TRUE))
      if (is.finite(d) && d > 0) return(d)
    }
  }
  n <- tryCatch(stats::nobs(m), error = function(e) NA_integer_)
  p <- tryCatch(length(lme4::fixef(m)), error = function(e) NA_integer_)
  if (is.finite(n) && is.finite(p) && n > p) return(n - p)
  Inf
}

# ------------------------------------------------------------------------------
# THE FOUR COMPONENT VIEWS
# ------------------------------------------------------------------------------
# One fitted model, four ways of looking at it. They are not four models and not
# four fits: each is the SAME coefficient vector with a different subset of the
# design columns allowed to contribute, so every one of them carries the same
# covariance and gets a band by the same quadratic form.
#
#   full            baseline + trend + harmonics -- everything the model fits
#   harmonics       the periodic part alone, centred on zero
#   baseline_harm   baseline + harmonics, with the non-periodic trend removed
#   trend           f(t) - f(0): the non-periodic change, exactly zero at t0
#
# The last one is a CONTRAST, not a sub-model: subtracting the t0 row is what
# makes it zero at the origin, and it is done to the design matrix so that the
# band narrows to zero there too, as it must.
DANCE_TRAJ_COMPONENTS <- c("full", "harmonics", "baseline_harm", "trend")

DANCE_TRAJ_COMPONENT_LABEL <- c(
  full          = "Full fitted trajectory",
  harmonics     = "Harmonics only (zero baseline)",
  baseline_harm = "Baseline + harmonics",
  trend         = "Nonperiodic change from origin")

# "Harmonics only" is every harmonic summed, which is the right default and the
# wrong view for the question "what is H2 doing on its own". With two harmonics
# fitted, the sum is the only thing that was ever drawable, so the 12 h component
# -- the one that decides how far the visible peak sits from the H1 acrophase --
# had no picture anywhere. Each harmonic gets its own option.
#
# The vocabulary depends on the FIT, not on a constant: a one-harmonic model has
# no H2 to offer and a no-trend model has no trend view, and a selector that
# offers either draws a flat zero line and calls it a component.
dance_traj_components <- function(fit) {
  if (!isTRUE(fit$ok)) return("full")
  K <- fit$spec$n_harmonics %||% 1L
  out <- c("full", "harmonics")
  if (K > 1L) out <- c(out, sprintf("harmonic%d", seq_len(K)))
  if (!identical(fit$spec$trend %||% "none", "none"))
    out <- c(out, "baseline_harm", "trend")
  out
}

dance_traj_component_label <- function(component, period = 24) {
  k <- dance_traj_component_harmonic(component)
  if (!is.na(k)) return(sprintf("Harmonic %d only (period %s h)", k,
                                format(round(period / k, 2), trim = TRUE)))
  # %||% would not help here: a name that is not in the vector indexes to NA,
  # not NULL, and an NA label reaches the selector as a blank option
  if (component %in% names(DANCE_TRAJ_COMPONENT_LABEL))
    unname(DANCE_TRAJ_COMPONENT_LABEL[component]) else component
}

# The harmonic number a component names, or NA for the views that are not one.
dance_traj_component_harmonic <- function(component) {
  m <- regmatches(component, regexpr("^harmonic([0-9]+)$", component))
  if (!length(m)) return(NA_integer_)
  as.integer(sub("^harmonic", "", m))
}

# Which fixed-effect columns survive for a given view. Works off the terms
# object, so a covariate or an interaction the user added is classified by the
# same rule as everything else rather than by name-matching.
dance_traj_component_mask <- function(fit, X, component) {
  spec <- fit$spec
  tl <- attr(stats::terms(stats::as.formula(spec$fixed_formula), data = spec$data),
             "term.labels")
  asg <- attr(X, "assign")
  has <- function(term, set) any(strsplit(term, ":", fixed = TRUE)[[1]] %in% set)
  kh <- dance_traj_component_harmonic(component)
  # one harmonic on its own: its own cos and sin columns and nothing else
  one_harm <- if (is.na(kh)) NULL else intersect(paste0(c("c", "s"), kh), spec$harm_terms)
  keep_term <- vapply(tl, function(term) {
    if (!is.na(kh)) return(has(term, one_harm))
    switch(component,
      full          = TRUE,
      harmonics     = has(term, spec$harm_terms),
      baseline_harm = !has(term, spec$trend_terms),
      trend         = has(term, spec$trend_terms))
  }, logical(1))
  # assign == 0 is the intercept: present for full and baseline+harmonics,
  # absent for the views that are explicitly centred on zero -- a single
  # harmonic among them, since it is a deviation and not a level
  keep_int <- component %in% c("full", "baseline_harm")
  ifelse(asg == 0L, keep_int, keep_term[pmax(asg, 1L)])
}

# ------------------------------------------------------------------------------
# ONE FITTED TRAJECTORY PER CELL, WITH A BAND
# ------------------------------------------------------------------------------
dance_traj_predict <- function(fit, times = NULL, conf = 0.95,
                               band = c("pointwise", "simultaneous"),
                               n_time = 200, component = "full") {
  band <- match.arg(band)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  # validated against what THIS fit can offer, rather than match.arg against a
  # constant: the per-harmonic views exist only when the harmonics do
  ok_comp <- dance_traj_components(fit)
  if (!identical(length(component), 1L) || !component %in% ok_comp)
    return(list(ok = FALSE, message = sprintf(
      "'%s' is not a component of this fit. Available: %s.",
      paste(component, collapse = ", "), paste(ok_comp, collapse = ", "))))
  spec <- fit$spec
  if (is.null(times))
    times <- seq(min(spec$data$t), max(spec$data$t), length.out = n_time)
  bv <- dance_traj_beta(fit)
  grid <- dance_traj_cell_grid(spec)
  q <- length(bv$beta)
  ddf <- dance_traj_ddf(fit)
  mult <- dance_traj_band_mult(band, conf, q, ddf)

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    X <- dance_traj_design_rows(fit, grid[i, , drop = FALSE], times)
    mask <- dance_traj_component_mask(fit, X, component)
    X[, !mask] <- 0
    if (identical(component, "trend")) {
      # f(t) - f(0), as a contrast on the design matrix, so the band is zero at
      # the origin rather than merely the point estimate being zero there
      X0 <- dance_traj_design_rows(fit, grid[i, , drop = FALSE], spec$t0)
      X0[, !mask] <- 0
      X <- sweep(X, 2, X0[1, ], "-")
    }
    keep <- intersect(colnames(X), names(bv$beta))
    X <- X[, keep, drop = FALSE]
    Vk <- bv$V[keep, keep, drop = FALSE]
    est <- as.numeric(X %*% bv$beta[keep])
    se <- sqrt(pmax(0, rowSums((X %*% Vk) * X)))
    data.frame(cell = grid$.cell[i], t = times, fit = est, se = se,
               lo = est - mult * se, hi = est + mult * se,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  list(ok = TRUE, table = out, cells = grid$.cell, grid = grid,
       component = component,
       component_label = dance_traj_component_label(component, spec$period),
       times = times, conf = conf, band = band, multiplier = mult, ddf = ddf,
       covariates_held_at = if (length(spec$covariates))
         sprintf("covariate(s) %s held at their sample mean",
                 paste(spec$covariates, collapse = ", ")) else NULL,
       note = if (identical(band, "pointwise")) paste(
         sprintf("POINTWISE %.0f%% intervals: correct at each time separately.", 100 * conf),
         "Reading them across the whole curve -- deciding two cells differ because",
         "the bands separate somewhere -- is not a", sprintf("%.0f%% statement.", 100 * conf),
         "Use band = 'simultaneous' for that, or the difference curve, which tests",
         "it directly.")
       else sprintf(paste("SIMULTANEOUS %.0f%% Scheffe band (multiplier %.2f on %d df,",
                          "against %.2f pointwise). It covers the whole %d-parameter",
                          "linear family at once, so it is conservative for this one curve."),
                    100 * conf, mult, q, stats::qnorm(1 - (1 - conf) / 2), q))
}

# ------------------------------------------------------------------------------
# THE DIFFERENCE BETWEEN TWO CELLS' TRAJECTORIES
# ------------------------------------------------------------------------------
# One linear contrast per time point, so the covariance between the two cells is
# carried exactly rather than assumed away. Simultaneous by default: a
# difference curve is read across its whole length or it is not read at all.
dance_traj_diff_curve <- function(fit, cell1, cell2, times = NULL, conf = 0.95,
                                  band = c("simultaneous", "pointwise"),
                                  n_time = 200) {
  band <- match.arg(band)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  grid <- dance_traj_cell_grid(spec)
  i <- match(cell1, grid$.cell); j <- match(cell2, grid$.cell)
  if (is.na(i) || is.na(j))
    return(list(ok = FALSE, message = sprintf(
      "Unknown cell. This design has: %s.", paste(grid$.cell, collapse = ", "))))
  if (i == j) return(list(ok = FALSE, message = "A cell cannot be contrasted with itself."))
  if (is.null(times))
    times <- seq(min(spec$data$t), max(spec$data$t), length.out = n_time)

  bv <- dance_traj_beta(fit)
  X1 <- dance_traj_design_rows(fit, grid[i, , drop = FALSE], times)
  X2 <- dance_traj_design_rows(fit, grid[j, , drop = FALSE], times)
  keep <- intersect(colnames(X1), names(bv$beta))
  D <- X1[, keep, drop = FALSE] - X2[, keep, drop = FALSE]
  Vk <- bv$V[keep, keep, drop = FALSE]
  est <- as.numeric(D %*% bv$beta[keep])
  se <- sqrt(pmax(0, rowSums((D %*% Vk) * D)))

  # The rank of the contrast family is the rank of D, not the width of beta:
  # a difference curve lives in the span of the columns that actually differ
  # between the two cells, which is what the Scheffe multiplier should cover.
  q <- tryCatch(qr(D)$rank, error = function(e) ncol(D))
  q <- max(1L, as.integer(q))
  ddf <- dance_traj_ddf(fit)
  mult <- dance_traj_band_mult(band, conf, q, ddf)
  tab <- data.frame(t = times, diff = est, se = se,
                    lo = est - mult * se, hi = est + mult * se,
                    excludes_zero = (est - mult * se) > 0 | (est + mult * se) < 0,
                    stringsAsFactors = FALSE)
  list(ok = TRUE, cell1 = cell1, cell2 = cell2, table = tab,
       conf = conf, band = band, multiplier = mult, rank = q, ddf = ddf,
       any_separation = any(tab$excludes_zero),
       note = paste(
         sprintf("%s %.0f%% band on the DIFFERENCE, formed as one linear contrast",
                 if (identical(band, "simultaneous")) "Simultaneous" else "Pointwise",
                 100 * conf),
         "so the covariance between the two cells is carried exactly. A difference",
         "band can exclude zero where the two cells' own bands overlap; that is the",
         "covariance doing its job, not a contradiction.",
         if (identical(band, "simultaneous"))
           sprintf("The Scheffe multiplier is %.2f on a rank-%d contrast family.", mult, q)
         else "Pointwise: correct at each time separately, not across the curve."))
}

# ------------------------------------------------------------------------------
# PAIRWISE CONTRASTS OVER AN ARBITRARY FACTORIAL STRUCTURE
# ------------------------------------------------------------------------------
# Every pair of design cells, for whichever quantity was asked for. "level" is
# the fitted value at a stated time (t0 by default) and is an exactly linear
# contrast; "amplitude" and "phase" are nonlinear maps of the cell's (cos, sin)
# pair and go through the delta method with the exact 4 x 4 covariance.
#
# MULTIPLICITY. m cells give m(m-1)/2 pairs and the family is the whole set, so
# an unadjusted p is wrong by construction once more than two cells exist. The
# default is Holm: it is uniformly more powerful than Bonferroni, needs no
# assumption about the dependence between contrasts, and -- unlike Tukey -- does
# not require the contrasts to be a balanced set of pairwise mean differences,
# which amplitude and phase contrasts are not.
DANCE_TRAJ_PAIR_WHAT <- c("level", "amplitude", "phase", DANCE_TRAJ_BLOCKS)

# WHICH HARMONIC. "Amplitude" and "Acrophase" are properties of ONE harmonic,
# and the panel took harmonic = 1 from a default argument that the UI never set
# -- so with two harmonics fitted it silently compared the 24 h component and
# said only "Amplitude". The 12 h component could not be compared at all.
#
# The harmonic is now part of the name: "amplitude2" is H2's amplitude. Bare
# "amplitude" and "phase" still mean H1, so every existing caller keeps working,
# but nothing the UI offers is bare -- each entry carries its harmonic and that
# harmonic's period, because a table of amplitudes that does not say which
# rhythm it is about is a table of unidentified numbers.
dance_traj_pair_parse <- function(what) {
  m <- regmatches(what, regexec("^(amplitude|phase)([0-9]*)$", what))[[1]]
  if (!length(m)) return(list(base = what, harmonic = NA_integer_))
  h <- if (nzchar(m[3])) as.integer(m[3]) else 1L
  list(base = m[2], harmonic = h)
}

dance_traj_pair_label <- function(what, period = 24, n_harmonics = 1L) {
  pp <- dance_traj_pair_parse(what)
  if (!is.na(pp$harmonic)) {
    nm <- if (identical(pp$base, "amplitude")) "Amplitude" else "Acrophase"
    # with a single harmonic there is nothing to disambiguate, but the period is
    # still worth stating: it is the number the reader needs to interpret it
    return(sprintf("%s — H%d (%s h)", nm, pp$harmonic,
                   format(round(period / pp$harmonic, 2), trim = TRUE)))
  }
  unname(DANCE_TRAJ_PAIR_LABEL[what]) %||% what
}

# Which comparisons this fit can actually answer. The scalar three always; the
# omnibus blocks only where they have terms and are not duplicates of each other
# -- with no trend, "shape" and "circadian" are the same set of columns, and
# offering both prints one test twice under two names, which reads as
# corroboration. Same rule the omnibus table uses, so the two panels cannot
# disagree about what exists.
dance_traj_pair_whats <- function(fit) {
  if (!isTRUE(fit$ok)) return(c("level", "amplitude", "phase"))
  spec <- fit$spec
  K <- max(1L, as.integer(spec$n_harmonics %||% 1L))
  # every fitted harmonic gets its own amplitude and acrophase entry
  out <- c("level",
           sprintf("amplitude%d", seq_len(K)), sprintf("phase%d", seq_len(K)))
  seen <- list()
  for (b in c("full", "shape", "circadian", "trend")) {
    cp <- DANCE_TRAJ_BLOCK_COMPONENTS(spec, b)
    if (!length(cp)) next
    key <- paste(sort(cp), collapse = "|")
    if (key %in% seen) next
    seen[[length(seen) + 1L]] <- key
    out <- c(out, b)
  }
  out
}

DANCE_TRAJ_PAIR_LABEL <- c(
  level     = "Level (at t = 0)",
  amplitude = "Amplitude",
  phase     = "Acrophase",
  full      = "Full trajectory (joint)",
  shape     = "Temporal shape (joint)",
  circadian = "Rhythmic block (joint)",
  trend     = "Non-periodic trend")

dance_traj_contrasts <- function(fit, what = "level",
                                 harmonic = 1, conf = 0.95, at_time = NULL,
                                 adjust = c("holm", "bonferroni", "none"),
                                 method = c("joint", "delta"), n_draw = 20000,
                                 df_method = c("auto", "kr", "satterthwaite")) {
  adjust <- match.arg(adjust); method <- match.arg(method)
  df_method <- match.arg(df_method)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  pp <- dance_traj_pair_parse(what)
  ok_what <- dance_traj_pair_whats(fit)
  # the harmonic named by the quantity wins over the argument default, which is
  # what the UI was silently relying on. Checked FIRST, so asking for a harmonic
  # the model does not fit says so, rather than reporting the name as unknown.
  if (identical(length(what), 1L) && !is.na(pp$harmonic)) {
    if (pp$harmonic > (spec$n_harmonics %||% 1L))
      return(list(ok = FALSE, message = sprintf(
        "This model fits %d harmonic(s), so there is no H%d to compare.",
        spec$n_harmonics %||% 1L, pp$harmonic)))
    harmonic <- pp$harmonic
    what <- pp$base
  } else if (!identical(length(what), 1L) ||
             !(what %in% ok_what || what %in% DANCE_TRAJ_PAIR_WHAT)) {
    return(list(ok = FALSE, message = sprintf(
      "'%s' is not a comparable quantity. Available: %s.",
      paste(what, collapse = ", "), paste(ok_what, collapse = ", "))))
  }
  if (!length(spec$design_terms))
    return(list(ok = FALSE, message = "No design factors: there are no cells to contrast."))
  grid <- dance_traj_cell_grid(spec)
  m <- nrow(grid)
  if (m < 2) return(list(ok = FALSE, message = "Only one design cell."))
  pairs <- utils::combn(m, 2)
  z <- stats::qnorm(1 - (1 - conf) / 2)

  # ---- THE OMNIBUS BLOCKS, PAIR BY PAIR ------------------------------------
  # Same components and same L machinery as the omnibus row above the table, so
  # a significant block test and an empty pairwise panel cannot come from two
  # different notions of what the block is.
  if (what %in% DANCE_TRAJ_BLOCKS) {
    rows <- lapply(seq_len(ncol(pairs)), function(p) {
      i <- pairs[1, p]; j <- pairs[2, p]
      r <- dance_traj_pair_block_test(fit, grid$.cell[i], grid$.cell[j], what,
                                      df_method = df_method, conf = conf)
      if (!isTRUE(r$ok)) return(NULL)
      data.frame(cell1 = grid$.cell[i], cell2 = grid$.cell[j],
                 estimate = r$estimate %||% NA_real_, se = r$se %||% NA_real_,
                 lo = r$lo %||% NA_real_, hi = r$hi %||% NA_real_,
                 statistic = r$statistic, df1 = r$df1, df2 = r$df2 %||% NA_real_,
                 defined = TRUE, p_draw = NA_real_, joint = r$df1 > 1,
                 method = r$method, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows))
      return(list(ok = FALSE, message = sprintf(
        "No %s contrast could be formed. The block may not be identifiable here.", what)))
    tab <- do.call(rbind, rows)
    tab$p_raw <- vapply(seq_len(nrow(tab)), function(k) {
      d2 <- tab$df2[k]
      if (is.finite(d2)) stats::pf(tab$statistic[k], tab$df1[k], d2, lower.tail = FALSE)
      else stats::pchisq(tab$statistic[k], tab$df1[k], lower.tail = FALSE)
    }, numeric(1))
    tab$p_adj <- tab$p_raw
    if (!identical(adjust, "none"))
      tab$p_adj <- stats::p.adjust(tab$p_raw, method = adjust)
    rownames(tab) <- NULL
    joint <- any(tab$joint)
    return(list(ok = TRUE, what = what, harmonic = harmonic, method = tab$method[1],
                block = what, joint = joint, table = tab, conf = conf, adjust = adjust,
                n_cells = m, n_contrasts = nrow(tab),
                label = dance_traj_pair_label(what, spec$period, spec$n_harmonics),
                unit = if (joint)
                  sprintf("joint test over %d component(s); no single effect size -- read the difference curve for the size",
                          max(tab$df1))
                else "difference in the response's own units",
                note = paste(
                  sprintf("%d pairwise %s contrasts over %d design cells, on the SAME components and the same L matrix as the omnibus test.",
                          nrow(tab), unname(DANCE_TRAJ_BLOCK_LABEL[what]), m),
                  if (joint) paste(
                    "Each row is a JOINT test: the cells are compared on every",
                    "component of the block at once, so there is an F and a p but no",
                    "single estimate. For the size and direction of a difference use",
                    "the difference curve, or one of the scalar quantities.")
                  else "One component, so each row is a scalar contrast with an interval.",
                  if (identical(adjust, "none"))
                    "NO multiplicity adjustment: with more than two cells these p-values are wrong as a family."
                  else sprintf("p_adj is %s-adjusted across the family of %d.", adjust, nrow(tab)),
                  sprintf("Degrees of freedom from %s, as chosen under Advanced.",
                          tab$method[1]))))
  }

  if (identical(what, "level")) {
    tt <- at_time %||% spec$t0
    bv <- dance_traj_beta(fit)
    Xs <- lapply(seq_len(m), function(i)
      dance_traj_design_rows(fit, grid[i, , drop = FALSE], tt))
    keep <- intersect(colnames(Xs[[1]]), names(bv$beta))
    Vk <- bv$V[keep, keep, drop = FALSE]
    rows <- lapply(seq_len(ncol(pairs)), function(p) {
      i <- pairs[1, p]; j <- pairs[2, p]
      d <- Xs[[i]][, keep, drop = FALSE] - Xs[[j]][, keep, drop = FALSE]
      est <- as.numeric(d %*% bv$beta[keep])
      se <- sqrt(max(0, as.numeric(d %*% Vk %*% t(d))))
      data.frame(cell1 = grid$.cell[i], cell2 = grid$.cell[j],
                 estimate = est, se = se, lo = est - z * se, hi = est + z * se,
                 statistic = est / se, defined = TRUE, p_draw = NA_real_,
                 stringsAsFactors = FALSE)
    })
    tab <- do.call(rbind, rows)
    unit <- sprintf("fitted value at t = %.3f %s", tt, spec$time_units %||% "")
  } else {
    co <- dance_traj_cell_coefs(fit, harmonic)
    if (!isTRUE(co$ok)) return(co)
    # emmeans orders its grid its own way; align by cell label rather than hope
    ord <- match(grid$.cell, co$cells)
    if (anyNA(ord))
      return(list(ok = FALSE, message = paste(
        "The cell labels from emmeans do not match the design grid, so a contrast",
        "could be attributed to the wrong pair. Refusing rather than guessing.")))
    rows <- lapply(seq_len(ncol(pairs)), function(p) {
      i <- ord[pairs[1, p]]; j <- ord[pairs[2, p]]
      # BRIEF 7. The joint path draws from the 4-dimensional distribution of both
      # cells' (cos, sin) pairs and reads the contrast off the draws -- amplitude
      # by quantiles, phase circularly. Nothing is linearised, so neither the
      # amplitude interval can cross zero nor the phase interval exceed the
      # circle. The delta method stays as the fast option for a rhythm already
      # far from the origin.
      if (identical(method, "joint")) {
        pj <- dance_traj_pair_joint(co, i, j, conf, n_draw = n_draw)
        if (is.null(pj)) return(NULL)
        if (identical(what, "phase"))
          return(data.frame(cell1 = pj$cell1, cell2 = pj$cell2,
                            estimate = pj$diff_time, se = pj$se_time,
                            lo = pj$lo, hi = pj$hi,
                            statistic = if (isTRUE(pj$defined)) pj$diff_time / pj$se_time else NA_real_,
                            defined = isTRUE(pj$defined), p_draw = pj$phase_p,
                            stringsAsFactors = FALSE))
        return(data.frame(cell1 = pj$cell1, cell2 = pj$cell2,
                          estimate = pj$amp_diff, se = pj$amp_se,
                          lo = pj$amp_lo, hi = pj$amp_hi,
                          statistic = pj$amp_diff / pj$amp_se,
                          defined = TRUE, p_draw = pj$amp_p,
                          stringsAsFactors = FALSE))
      }
      if (identical(what, "phase")) {
        pc <- dance_traj_phase_contrast(co, i, j, conf)
        if (is.null(pc)) return(NULL)
        return(data.frame(cell1 = co$cells[i], cell2 = co$cells[j],
                          estimate = pc$diff_time, se = pc$se_time,
                          lo = pc$lo, hi = pc$hi,
                          statistic = if (isTRUE(pc$defined)) pc$diff_time / pc$se_time else NA_real_,
                          defined = isTRUE(pc$defined), p_draw = NA_real_,
                          stringsAsFactors = FALSE))
      }
      cv <- dance_traj_pair_cov(co, i, j)
      if (is.null(cv)) return(NULL)
      a1 <- co$a[i]; b1 <- co$b[i]; a2 <- co$a[j]; b2 <- co$b[j]
      A1 <- sqrt(a1^2 + b1^2); A2 <- sqrt(a2^2 + b2^2)
      g <- c(-a1 / A1, -b1 / A1, a2 / A2, b2 / A2)   # d(A2 - A1)
      se <- sqrt(max(0, as.numeric(t(g) %*% cv$V %*% g)))
      est <- A2 - A1
      data.frame(cell1 = co$cells[i], cell2 = co$cells[j],
                 estimate = est, se = se, lo = est - z * se, hi = est + z * se,
                 statistic = est / se, defined = A1 > 0 && A2 > 0,
                 p_draw = NA_real_, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(list(ok = FALSE, message = "No contrast could be formed."))
    tab <- do.call(rbind, rows)
    # NAME THE HARMONIC. A column of amplitude differences that does not say
    # which rhythm it is about is a column of unidentified numbers, and this
    # panel showed H1 under a bare "Amplitude" for as long as it existed.
    unit <- if (identical(what, "phase"))
      sprintf("H%d acrophase difference in %s, on an effective period of %.4g h",
              harmonic, spec$time_units %||% "time units", co$effective_period)
    else sprintf("H%d amplitude difference in the response's own units (period %.4g h)",
                 harmonic, co$effective_period)
  }

  tab$p_raw <- ifelse(tab$defined, 2 * stats::pnorm(-abs(tab$statistic)), NA_real_)
  tab$p_adj <- tab$p_raw
  if (!identical(adjust, "none"))
    tab$p_adj[tab$defined] <- stats::p.adjust(tab$p_raw[tab$defined], method = adjust)
  rownames(tab) <- NULL
  list(ok = TRUE, what = what, harmonic = harmonic, unit = unit, method = method,
       # the label the panel shows, carried on the result so the table and its
       # heading cannot name different harmonics
       label = dance_traj_pair_label(
         if (what %in% c("amplitude", "phase")) paste0(what, harmonic) else what,
         spec$period, spec$n_harmonics),
       table = tab, conf = conf, adjust = adjust, n_cells = m,
       n_contrasts = nrow(tab),
       note = paste(
         sprintf("%d pairwise contrasts over %d design cells.", nrow(tab), m),
         if (identical(adjust, "none"))
           paste("NO multiplicity adjustment: with more than two cells these p-values",
                 "are wrong as a family and only p_raw is meaningful, one contrast at",
                 "a time, chosen before looking.")
         else sprintf(paste("p_adj is %s-adjusted across the whole family of %d;",
                            "p_raw is kept beside it so the adjustment is visible."),
                      adjust, nrow(tab)),
         if (identical(what, "level"))
           "Intervals are Wald on a normal reference, not the Kenward-Roger F of the omnibus."
         else if (identical(method, "joint"))
           sprintf(paste("Amplitude and phase are nonlinear maps of the cell's (cos, sin)",
                         "pair. These intervals are drawn from the JOINT distribution of",
                         "both cells' pairs (%d draws, exact cross-cell covariance):",
                         "nothing is linearised, the amplitude interval cannot cross",
                         "zero and the phase interval is an arc. p_draw is the",
                         "proportion of draws on the far side of zero."), n_draw)
         else paste("Amplitude and phase are nonlinear maps of the cell's (cos, sin) pair;",
                    "these are DELTA-METHOD intervals with the exact cross-cell covariance,",
                    "accurate when the rhythm is well separated from the origin."),
         if (any(!tab$defined))
           sprintf(paste("%d contrast(s) are UNDEFINED because a cell's amplitude",
                         "interval covers zero (Bingham's rule) and are excluded from",
                         "the adjustment."), sum(!tab$defined)) else NULL))
}

# ------------------------------------------------------------------------------
# THE CELL-KEYED ADAPTER
# ------------------------------------------------------------------------------
# The existing plots and report blocks read a list keyed by group name, each
# entry carrying an intercept, a coefficient vector, amplitudes and acrophases.
# This builds the same shape from the trajectory fit, keyed by DESIGN CELL
# ("young x sleep-deprived") rather than by a single grouping factor, so the
# plotting code works over a 2 x 3 x 2 without knowing it.
#
# WHAT THIS ADAPTER REFUSES TO FAKE
# ---------------------------------
# The legacy structure carries quantities that only exist because the old path
# fitted one model per participant and then averaged the results. A one-stage
# model does not produce a sample of per-participant estimates, so several of
# those fields have no counterpart, and the temptation is to fill them with the
# nearest-looking number. That would be a rename, not a translation:
#
#   sd_amplitudes    was the SD of per-participant amplitudes -- a DISPERSION.
#                    The model's standard error of the cell amplitude is a
#                    PRECISION. They answer different questions and shrink at
#                    different rates with n. Both are supplied, under their own
#                    names, and the dispersion comes from the random-effects
#                    covariance (delta method on the participant (cos, sin)
#                    spread) rather than from the standard error.
#   amp_arithmetic   the arithmetic mean of per-participant amplitudes, which is
#                    biased upward relative to the vector mean because
#                    |E[v]| <= E[|v|]. There is no sample to average here.
#   resultants       a mean resultant length is a statistic OF A SAMPLE of
#                    angles. The model gives a distribution, not a sample.
#   variance_decomp  per-participant R^2 shares, averaged. The one-stage
#                    analogue is a marginal/conditional R^2 on the whole fit,
#                    which is a different decomposition and is not silently
#                    substituted here.
#
# Those fields are returned as NULL/NA with a reason in $unavailable, so a
# caller that needs one gets an empty panel and an explanation instead of a
# plausible wrong number. dance_traj_group_fits_gaps() lists them.
dance_traj_group_fits <- function(fit, conf = 0.95) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  grid <- dance_traj_cell_grid(spec)
  bv <- dance_traj_beta(fit)
  z <- stats::qnorm(1 - (1 - conf) / 2)

  # participant-level dispersion, from the random-effects covariance of the
  # rung that was actually fitted
  vc <- tryCatch(lme4::VarCorr(fit$model), error = function(e) NULL)
  Sub <- NULL
  if (!is.null(vc)) {
    cand <- vc[[spec$subject]] %||% vc[["subject"]] %||% NULL
    if (!is.null(cand)) Sub <- as.matrix(cand)
  }
  sd_from <- function(nm) if (!is.null(Sub) && nm %in% rownames(Sub))
    sqrt(max(0, Sub[nm, nm])) else NA_real_
  amp_sd <- function(a, b, h) {
    ck <- paste0("c", h); sk <- paste0("s", h)
    if (is.null(Sub) || !all(c(ck, sk) %in% rownames(Sub))) return(NA_real_)
    A <- sqrt(a^2 + b^2); if (!(A > 0)) return(NA_real_)
    g <- c(a / A, b / A)
    sqrt(max(0, as.numeric(t(g) %*% Sub[c(ck, sk), c(ck, sk)] %*% g)))
  }

  # cell coefficients and their standard errors, per harmonic
  co <- lapply(seq_len(spec$n_harmonics), function(h) {
    r <- dance_traj_cell_coefs(fit, h)
    if (isTRUE(r$ok)) r else NULL
  })

  n_per_cell <- stats::setNames(rep(NA_integer_, nrow(grid)), grid$.cell)
  if (length(spec$design_terms)) {
    key <- do.call(paste, c(lapply(spec$design_terms,
                                   function(f) as.character(spec$data[[f]])),
                            list(sep = " x ")))
    sp <- split(as.character(spec$data$subject), key)
    got <- vapply(sp, function(s) length(unique(s)), integer(1))
    n_per_cell[names(got)] <- got
  } else n_per_cell[] <- spec$n_participants

  # The coefficient on a trend column in a cell is the derivative of the linear
  # predictor with respect to that column, which emtrends gives exactly; a
  # finite difference on the design row would be a guess. One call per column,
  # not one per column per cell.
  trend_tr <- stats::setNames(lapply(spec$trend_terms, function(nm) {
    em <- dance_emtrends(fit, stats::as.formula(
      paste("~", paste(spec$design_terms, collapse = " * "))), nm)
    if (is.null(em)) return(NULL)
    s <- summary(em)
    data.frame(.cell = do.call(paste, c(lapply(spec$design_terms,
                                               function(f) as.character(s[[f]])),
                                        list(sep = " x "))),
               .est = s[[paste0(nm, ".trend")]], stringsAsFactors = FALSE)
  }), spec$trend_terms)

  out <- list()
  for (i in seq_len(nrow(grid))) {
    cell <- grid$.cell[i]
    X0 <- dance_traj_design_rows(fit, grid[i, , drop = FALSE], spec$t0)
    keep <- intersect(colnames(X0), names(bv$beta))
    # beta_0 for this cell: the fitted value at t = t0 with every basis column
    # at its t0 value. The harmonics are NOT zero there, so this is not the
    # intercept; it is recovered by subtracting the basis contribution.
    b0 <- NA_real_; amps <- acro <- amp_se <- amp_sdv <- rep(NA_real_, spec$n_harmonics)
    cc <- ss <- rep(NA_real_, spec$n_harmonics)
    for (h in seq_len(spec$n_harmonics)) {
      r <- co[[h]]; if (is.null(r)) next
      k <- match(cell, r$cells); if (is.na(k)) next
      cc[h] <- r$a[k]; ss[h] <- r$b[k]
      A <- sqrt(r$a[k]^2 + r$b[k]^2)
      amps[h] <- A
      ph <- atan2(r$b[k], r$a[k]) %% (2 * pi)
      acro[h] <- ph
      V <- if (!is.null(r$joint)) r$joint[[k]] else NULL
      if (!is.null(V) && A > 0) {
        g <- c(r$a[k] / A, r$b[k] / A)
        amp_se[h] <- sqrt(max(0, as.numeric(t(g) %*% V %*% g)))
      }
      amp_sdv[h] <- amp_sd(r$a[k], r$b[k], h)
    }
    # trend coefficients for the cell, read off the same linear predictor
    tr <- stats::setNames(rep(NA_real_, length(spec$trend_terms)), spec$trend_terms)
    for (nm in spec$trend_terms) {
      s <- trend_tr[[nm]]
      if (is.null(s)) next
      kk <- match(cell, s$.cell)
      if (!is.na(kk)) tr[nm] <- s$.est[kk]
    }
    # the constant: fitted value at t0 minus every basis column's contribution
    basis0 <- 0
    for (h in seq_len(spec$n_harmonics)) {
      if (is.na(cc[h])) next
      basis0 <- basis0 + cc[h] * 1 + ss[h] * 0        # at t = t0, cos = 1, sin = 0
    }
    fit0 <- as.numeric(X0[, keep, drop = FALSE] %*% bv$beta[keep])
    b0 <- fit0 - basis0                                # trend columns are 0 at t0

    coefs <- c(b0, unname(tr))
    for (h in seq_len(spec$n_harmonics)) coefs <- c(coefs, cc[h], ss[h])

    out[[cell]] <- list(
      group = cell, is_unassigned = FALSE,
      n = unname(n_per_cell[cell]),
      mean_mesor = b0, intercept = b0,
      rhythm_adjusted_mean = dance_rhythm_adjusted_mean(
        b0, spec$trend, unname(tr), min(spec$data$t), max(spec$data$t), 0),
      trend_coefs = unname(tr),
      trend_params = if (length(tr)) lapply(tr, function(v) list(mean = v)) else NULL,
      mean_coefs = coefs,
      mean_amplitudes = amps, mean_acrophases_rad = acro,
      mean_acrophases_time = acro * (spec$period / seq_len(spec$n_harmonics)) / (2 * pi),
      mean_amplitude = amps[1], mean_acrophase_rad = acro[1],
      mean_acrophase_time = (acro * (spec$period / seq_len(spec$n_harmonics)) / (2 * pi))[1],
      # PRECISION of the cell estimate -- new, and not a rename of anything
      se_amplitudes = amp_se, se_amplitude = amp_se[1],
      amplitude_lo = amps - z * amp_se, amplitude_hi = amps + z * amp_se,
      # DISPERSION across participants, from the random-effects covariance
      sd_amplitudes = amp_sdv, sd_amplitude = amp_sdv[1],
      sd_mesor = sd_from("(Intercept)"),
      dispersion_source = if (all(is.na(amp_sdv)))
        paste("The fitted random structure carries no participant-specific harmonic",
              "terms, so there is no model-based between-participant amplitude",
              "dispersion to report.")
        else paste("Between-participant SDs come from the random-effects covariance of",
                   sprintf("the fitted structure (%s), not from a sample of",
                           fit$re_label %||% "the fitted model"),
                   "per-participant point estimates."),
      # deliberately absent -- see the header
      amp_arithmetic = NULL, resultants = NULL, variance_decomp = NULL,
      stringsAsFactors = FALSE)
  }
  attr(out, "provenance") <- sprintf(
    "one mixed-effects trajectory fit (%s), cells keyed by %s",
    fit$re_label %||% "fitted", paste(spec$design_terms, collapse = " x "))
  attr(out, "n_fitted") <- spec$n_participants
  attr(out, "n_in_groups") <- sum(n_per_cell, na.rm = TRUE)
  attr(out, "unavailable") <- dance_traj_group_fits_gaps()
  attr(out, "conf") <- conf
  out
}

# What the one-stage model cannot supply, and why. A caller renders this next to
# an empty panel instead of inventing a number for it.
dance_traj_group_fits_gaps <- function() {
  list(
    amp_arithmetic = paste(
      "The arithmetic mean of per-participant amplitudes. A one-stage fit has no",
      "sample of per-participant amplitudes to average, and the quantity is",
      "upward-biased relative to the vector mean in any case, since |E[v]| <= E[|v|]."),
    resultants = paste(
      "The mean resultant length and Rayleigh test. These are statistics OF A",
      "SAMPLE of per-participant angles. The model estimates a distribution of",
      "participant phases, not a sample of them, so there is no resultant to take."),
    variance_decomp = paste(
      "The averaged per-participant R^2 shares. The one-stage analogue is a",
      "marginal and conditional R^2 on the whole fit -- a different decomposition,",
      "over different variance components -- and it is not substituted silently."))
}

# ------------------------------------------------------------------------------
# SIMPLE EFFECTS AND INTERACTION CONTRASTS
# ------------------------------------------------------------------------------
# All pairs is the blunt instrument. In a factorial design the questions people
# actually ask are narrower and more powerful:
#
#   SIMPLE EFFECT     the effect of factor A at ONE level of factor B
#                     ("does caffeine shift the rhythm, in the young group?")
#   INTERACTION       whether that simple effect is the SAME at every level of B
#                     ("does caffeine shift it MORE in the young than the old?")
#
# Both come from the one fitted model, never from refitting a subset: refitting
# within a level throws away the pooled residual and the random structure, and
# the resulting standard error is not comparable with anything else on the page.
#
# The interaction contrast is a difference of differences, which is why it is
# formed as a single linear combination rather than by comparing two p-values.
# "Significant in one group and not the other" is not an interaction, and a
# module that offers simple effects without offering the interaction invites
# exactly that mistake.
dance_traj_simple_effects <- function(fit, effect, at, what = "level",
                                      harmonic = 1, conf = 0.95, at_time = NULL,
                                      adjust = c("holm", "bonferroni", "none"),
                                      df_method = c("auto", "kr", "satterthwaite")) {
  adjust <- match.arg(adjust); df_method <- match.arg(df_method)
  # the same vocabulary as the whole-design table: a slice of the design can be
  # asked every question the design can
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  if (!identical(length(what), 1L) ||
      !(what %in% dance_traj_pair_whats(fit) || what %in% DANCE_TRAJ_PAIR_WHAT))
    return(list(ok = FALSE, message = sprintf(
      "'%s' is not a comparable quantity. Available: %s.",
      paste(what, collapse = ", "), paste(dance_traj_pair_whats(fit), collapse = ", "))))
  spec <- fit$spec
  dt <- spec$design_terms
  if (!effect %in% dt)
    return(list(ok = FALSE, message = sprintf(
      "'%s' is not a design factor. This model has: %s.", effect, paste(dt, collapse = ", "))))
  hold <- setdiff(dt, effect)
  if (!all(names(at) %in% hold) || !setequal(names(at), hold))
    return(list(ok = FALSE, message = sprintf(
      "`at` must name every other design factor exactly once: %s.",
      paste(hold, collapse = ", "))))

  grid <- dance_traj_cell_grid(spec)
  keep <- rep(TRUE, nrow(grid))
  for (f in hold) keep <- keep & as.character(grid[[f]]) == as.character(at[[f]])
  if (sum(keep) < 2)
    return(list(ok = FALSE, message = "That slice of the design holds fewer than two cells."))

  all_ct <- dance_traj_contrasts(fit, what, harmonic, conf, at_time, adjust = "none",
                                 df_method = df_method)
  if (!isTRUE(all_ct$ok)) return(all_ct)
  cells <- grid$.cell[keep]
  tab <- all_ct$table[all_ct$table$cell1 %in% cells & all_ct$table$cell2 %in% cells, ,
                      drop = FALSE]
  if (!nrow(tab)) return(list(ok = FALSE, message = "No contrast lies inside that slice."))
  # The family is the SLICE, not the whole design: that is the point of asking
  # for a simple effect, and it is stated so nobody reads these p-values as
  # though they had been adjusted across everything.
  tab$p_adj <- tab$p_raw
  if (!identical(adjust, "none"))
    tab$p_adj[tab$defined] <- stats::p.adjust(tab$p_raw[tab$defined], method = adjust)
  rownames(tab) <- NULL
  list(ok = TRUE, what = what, effect = effect, at = at, unit = all_ct$unit,
       joint = all_ct$joint %||% FALSE, label = all_ct$label,
       table = tab, conf = conf, adjust = adjust, cells = cells,
       note = paste(
         sprintf("The simple effect of %s at %s.", effect,
                 paste(sprintf("%s = %s", names(at), unlist(at)), collapse = ", ")),
         "Estimated from the ONE fitted model, not by refitting this slice: a",
         "within-slice refit would discard the pooled residual and the random",
         "structure, and its standard error would not be comparable with anything",
         "else reported.",
         sprintf("Multiplicity is adjusted across the %d contrast(s) in this slice only,",
                 nrow(tab)),
         "not across the whole design.",
         "A simple effect that is significant in one slice and not another is NOT an",
         "interaction -- use dance_traj_interaction_contrast() for that question."))
}

# The difference of differences: (cell_a1 - cell_a2) - (cell_b1 - cell_b2).
# Level only, because it is the one of the three that is exactly linear; an
# interaction in amplitude or phase is a difference of two nonlinear maps and
# needs the joint draws, which dance_traj_interaction_contrast() refuses to fake.
dance_traj_interaction_contrast <- function(fit, cell_a1, cell_a2, cell_b1, cell_b2,
                                            conf = 0.95, at_time = NULL) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  grid <- dance_traj_cell_grid(spec)
  idx <- match(c(cell_a1, cell_a2, cell_b1, cell_b2), grid$.cell)
  if (anyNA(idx))
    return(list(ok = FALSE, message = sprintf(
      "Unknown cell (%s). This design has: %s.",
      paste(c(cell_a1, cell_a2, cell_b1, cell_b2)[is.na(idx)], collapse = ", "),
      paste(grid$.cell, collapse = ", "))))
  tt <- at_time %||% spec$t0
  bv <- dance_traj_beta(fit)
  X <- lapply(idx, function(i) dance_traj_design_rows(fit, grid[i, , drop = FALSE], tt))
  keep <- intersect(colnames(X[[1]]), names(bv$beta))
  d <- (X[[1]][, keep, drop = FALSE] - X[[2]][, keep, drop = FALSE]) -
       (X[[3]][, keep, drop = FALSE] - X[[4]][, keep, drop = FALSE])
  est <- as.numeric(d %*% bv$beta[keep])
  se <- sqrt(max(0, as.numeric(d %*% bv$V[keep, keep] %*% t(d))))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  list(ok = TRUE, estimate = est, se = se, lo = est - z * se, hi = est + z * se,
       statistic = est / se, p = 2 * stats::pnorm(-abs(est / se)),
       at_time = tt, conf = conf,
       contrast = sprintf("(%s - %s) - (%s - %s)", cell_a1, cell_a2, cell_b1, cell_b2),
       note = paste(
         "A difference of differences, formed as ONE linear combination, so its",
         "standard error carries every covariance among the four cells. This is the",
         "interaction question. It is NOT answered by observing that one difference",
         "is significant and the other is not: two tests either side of a threshold",
         "say nothing about whether they differ from each other.",
         sprintf("Evaluated at t = %.4g %s.", tt, spec$time_units %||% "")))
}

# ------------------------------------------------------------------------------
# PARTICIPANT-LEVEL TRAJECTORIES: CONDITIONAL MODES, AND THEY ARE SHRUNKEN
# ------------------------------------------------------------------------------
# BLUPs are enormously useful to look at and routinely over-read. A conditional
# mode is not an estimate of that participant's parameter in the sense an
# independent per-participant fit would give: it is pulled toward the population
# mean by an amount that depends on how much data that participant contributed
# and how large the random-effect variance is. So:
#
#   - a histogram of BLUP amplitudes is NARROWER than the true spread of
#     participant amplitudes, and its SD is not an estimate of that spread --
#     the random-effect SD is;
#   - a participant with fewer or noisier observations is pulled further, so the
#     ORDER can differ from the order of independent fits;
#   - the shrinkage is a feature, not a defect: these are better predictions of
#     the individual than the unshrunk estimates. They are just not the same
#     quantity, and a plot that does not say so invites the wrong reading.
#
# Every return here carries `shrunken = TRUE` and that warning, and the
# population spread is reported alongside so the two can be compared.
dance_traj_participant_curves <- function(fit, harmonic = 1, conf = 0.95) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  re <- tryCatch(lme4::ranef(fit$model, condVar = FALSE), error = function(e) NULL)
  if (is.null(re)) return(list(ok = FALSE, message = "Conditional modes are unavailable for this fit."))
  grp <- spec$subject %||% "subject"
  if (!grp %in% names(re))
    return(list(ok = FALSE, message = sprintf("The fitted model has no '%s' random effect.", grp)))
  ck <- paste0("c", harmonic); sk <- paste0("s", harmonic)
  R_subj <- as.data.frame(re[[grp]])
  cv <- spec$curve_var %||% "curve"
  R_curve <- if (cv %in% names(re)) as.data.frame(re[[cv]]) else NULL
  has_subj  <- all(c(ck, sk) %in% names(R_subj))
  has_curve <- !is.null(R_curve) && all(c(ck, sk) %in% names(R_curve))
  if (!has_subj && !has_curve)
    return(list(ok = FALSE, message = paste(
      "The fitted random structure carries no participant-specific harmonic terms,",
      sprintf("so there are no per-participant rhythms to report (rung used: %s).",
              fit$re_label %||% "unknown"))))

  # EACH PARTICIPANT DEVIATES FROM THEIR OWN CELL, not from the grand mean.
  # The first version added every conditional mode to mean(co$a), the average
  # over all design cells -- so with groups that genuinely differ, every
  # individual rhythm was pulled toward the overall average and the
  # per-participant amplitudes of a high-amplitude group all read low. A
  # conditional mode is a deviation from the fixed-effect prediction FOR THAT
  # OBSERVATION'S CELL, and that is what it is added to here.
  d <- spec$data
  dts <- spec$design_terms
  if (length(dts)) {
    co <- dance_traj_cell_coefs(fit, harmonic)
    if (!isTRUE(co$ok)) return(list(ok = FALSE, message = co$message))
  } else {
    # no design: one cell, whose rhythm is the fixed-effect pair itself
    bv <- dance_traj_beta(fit)
    if (!all(c(ck, sk) %in% names(bv$beta)))
      return(list(ok = FALSE, message = sprintf("Harmonic %d was not fitted.", harmonic)))
    co <- list(ok = TRUE, cells = "(all)", a = unname(bv$beta[[ck]]), b = unname(bv$beta[[sk]]))
  }
  cell_of <- if (length(dts))
    do.call(paste, c(lapply(dts, function(f) as.character(d[[f]])), list(sep = " x ")))
    else rep(co$cells[1], nrow(d))
  # one row per CURVE (participant x within-cell): that is the unit that carries
  # its own conditional modes, and a participant measured in two conditions has
  # two rhythms, not one averaged over both
  key <- data.frame(subject = as.character(d[[grp]]),
                    curve = if (cv %in% names(d)) as.character(d[[cv]]) else as.character(d[[grp]]),
                    cell = cell_of, stringsAsFactors = FALSE)
  key <- key[!duplicated(key$curve), , drop = FALSE]
  ci <- match(key$cell, co$cells)
  if (anyNA(ci))
    return(list(ok = FALSE, message = "A curve's design cell is not among the fitted cells."))
  a <- co$a[ci]; b <- co$b[ci]
  if (has_subj) {
    si <- match(key$subject, rownames(R_subj))
    a <- a + ifelse(is.na(si), 0, R_subj[[ck]][si])
    b <- b + ifelse(is.na(si), 0, R_subj[[sk]][si])
  }
  if (has_curve) {
    qi <- match(key$curve, rownames(R_curve))
    a <- a + ifelse(is.na(qi), 0, R_curve[[ck]][qi])
    b <- b + ifelse(is.na(qi), 0, R_curve[[sk]][qi])
  }
  k <- (spec$period / harmonic) / (2 * pi)
  vc <- tryCatch(as.matrix(lme4::VarCorr(fit$model)[[grp]]), error = function(e) NULL)
  tab <- data.frame(subject = key$subject, curve = key$curve, cell = key$cell,
                    beta_cos = a, beta_sin = b,
                    amplitude = sqrt(a^2 + b^2),
                    acrophase_rad = atan2(b, a) %% (2 * pi),
                    acrophase_time = (atan2(b, a) %% (2 * pi)) * k,
                    stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  # the population direction, for the between-participant SD of amplitude: the
  # mean of the CELL vectors, which is a summary and is labelled as one
  pop_a <- mean(co$a); pop_b <- mean(co$b)
  list(ok = TRUE, table = tab, harmonic = harmonic, shrunken = TRUE,
       effective_period = spec$period / harmonic,
       re_label = fit$re_label,
       per_curve = has_curve,
       sd_blup_amplitude = stats::sd(tab$amplitude),
       sd_population_amplitude = if (!is.null(vc) && all(c(ck, sk) %in% rownames(vc))) {
         A <- sqrt(pop_a^2 + pop_b^2)
         if (A > 0) { g <- c(pop_a / A, pop_b / A)
           sqrt(max(0, as.numeric(t(g) %*% vc[c(ck, sk), c(ck, sk)] %*% g))) } else NA_real_
       } else NA_real_,
       label = "model-based shrunken participant estimates (conditional modes / BLUPs)",
       note = paste(
         "THESE ARE SHRUNKEN. A conditional mode is pulled toward its CELL's",
         "fixed-effect rhythm by an amount that depends on how much data that",
         "participant contributed and how large the random-effect variance is.",
         "They are better PREDICTIONS of each individual than independent",
         "per-participant fits, but they are not the same quantity.",
         "In particular the spread of these values UNDERSTATES the spread of true",
         "participant parameters: compare sd_blup_amplitude with",
         "sd_population_amplitude, which is the model's estimate of the real",
         "between-participant SD. A histogram or boxplot of this column must say",
         "so on its face, or it will be read as the population distribution."))
}

# ------------------------------------------------------------------------------
# WHERE THE DRAWN CURVE ACTUALLY PEAKS
# ------------------------------------------------------------------------------
# Three different times get called "the peak" and they are not the same number:
#
#   H1 acrophase        the maximum of the FIRST HARMONIC ALONE. With a second
#                       harmonic in the model this is not where the curve peaks:
#                       H2 adds its own maximum, and the sum peaks somewhere
#                       between them, pulled towards whichever is larger.
#   curve peak          the maximum of the fitted trajectory AS DRAWN -- level,
#                       trend and every harmonic together. This is the time a
#                       reader takes off the plot, and it is what this returns.
#   two-stage acrophase what a per-participant single-harmonic cosinor gives,
#                       averaged over the group. A one-harmonic no-trend fit has
#                       to absorb H2 and the trend into its single cosine, so its
#                       acrophase is displaced again, by however much they carry.
#
# Reading any one of them against a plot of another is how an acrophase gets
# called wrong when nothing is. The peak is found by evaluating the fitted curve
# on a fine grid over ONE period from the origin -- a grid, not calculus, because
# the sum of harmonics plus an arbitrary trend has no closed-form maximum.
dance_traj_curve_peaks <- function(fit, component = "full", n_time = 1441L,
                                   within_period = TRUE) {
  if (!isTRUE(fit$ok)) return(NULL)
  spec <- fit$spec
  t_hi <- if (isTRUE(within_period))
    min(spec$t0 + spec$period, max(spec$data$t, na.rm = TRUE))
  else max(spec$data$t, na.rm = TRUE)
  if (!is.finite(t_hi) || t_hi <= spec$t0) return(NULL)
  times <- seq(spec$t0, t_hi, length.out = max(51L, as.integer(n_time)))
  pr <- dance_traj_predict(fit, times = times, component = component)
  if (!isTRUE(pr$ok)) return(NULL)
  # THE CLOCK TIME IS COMPUTED HERE, not by whoever renders the row. The first
  # version returned only peak_t -- model-elapsed hours -- and the table
  # converted it with `peak_t %% period`, forgetting to add the origin back. On
  # a study starting at 08:00 that reported an 08:00 (+1d) peak as 00:00. The
  # offset now travels with the number, so a caller cannot convert it wrongly,
  # and the test checks the CLOCK value rather than the elapsed one.
  # clock = (t - t0) + origin. The ELAPSED part is t - t0, not t: peak_t is an
  # absolute value on whatever axis the model was built on, and only the offset
  # FROM THE BASIS ORIGIN is what clock_origin is the clock time of. In this app
  # the model axis already starts at zero, so t0 = 0 and the two agree -- which
  # is exactly why leaving t0 out would survive every run here and be wrong for
  # anyone who fitted on raw clock times.
  co <- fit$clock_origin %||% 0
  P <- spec$period
  to_clock <- function(t) (t - spec$t0 + co) %% P
  sp <- split(pr$table, pr$table$cell)
  out <- do.call(rbind, lapply(names(sp), function(cl) {
    d <- sp[[cl]]
    i <- which.max(d$fit); j <- which.min(d$fit)
    data.frame(cell = cl,
               peak_t = d$t[i], peak_fit = d$fit[i],
               peak_clock = to_clock(d$t[i]),
               trough_t = d$t[j], trough_fit = d$fit[j],
               trough_clock = to_clock(d$t[j]),
               # on the boundary the maximum inside the window is not a turning
               # point of the curve, and saying so is the difference between a
               # peak and the end of the search
               peak_interior = i > 1L && i < nrow(d),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  list(ok = TRUE, table = out, component = component,
       window = c(spec$t0, t_hi), clock_origin = co,
       window_clock = to_clock(c(spec$t0, t_hi)),
       resolution_min = 60 * (t_hi - spec$t0) / (length(times) - 1))
}

# ==============================================================================
# PER-PARTICIPANT VIEWS OF THE ONE MIXED MODEL
# ==============================================================================
# The two-stage approach has a fit per participant and derives everything else
# from those. The mixed approach has ONE fit and derives each participant from
# it: the cell's fixed effects plus that participant's conditional modes. These
# two helpers give tab 1 a curve per participant and tab 4 a row per
# participant, both SHRUNKEN and both labelled as such -- see
# dance_traj_participant_curves() for what shrinkage does to the spread.

# The conditional (participant-level) prediction on a time grid, one curve per
# CURVE -- a participant measured in two conditions has two -- through the
# model's own predict(), so the random-effects terms enter exactly as fitted.
dance_traj_participant_predict <- function(fit, times = NULL, n_time = 120) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec; d <- spec$data
  if (is.null(times)) times <- seq(min(d$t, na.rm = TRUE), max(d$t, na.rm = TRUE),
                                   length.out = n_time)
  cv <- spec$curve_var %||% "curve"; grp <- spec$subject %||% "subject"
  cv <- if (cv %in% names(d)) cv else grp
  dts <- spec$design_terms
  curves <- unique(as.character(d[[cv]]))
  tt <- times - spec$t0
  rows <- lapply(curves, function(cl) {
    tmpl <- d[match(cl, as.character(d[[cv]])), , drop = FALSE]
    nd <- tmpl[rep(1L, length(times)), , drop = FALSE]
    nd$t <- times
    for (nm in spec$trend_terms) nd[[nm]] <- switch(nm,
      trend_lin = tt, trend_log = log1p(pmax(0, tt)),
      trend_sat = 1 - exp(-tt / spec$tau), stop(sprintf("unknown trend column '%s'", nm)))
    for (h in seq_len(spec$n_harmonics)) {
      w <- 2 * pi * h * tt / spec$period
      nd[[paste0("c", h)]] <- cos(w); nd[[paste0("s", h)]] <- sin(w)
    }
    pred <- tryCatch(as.numeric(stats::predict(fit$model, newdata = nd, re.form = NULL)),
                     error = function(e) rep(NA_real_, length(times)))
    data.frame(curve = cl, subject = as.character(tmpl[[grp]]),
               cell = if (length(dts))
                 paste(vapply(dts, function(f) as.character(tmpl[[f]]), character(1)),
                       collapse = " x ") else "(all)",
               t = times, fit = pred, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  list(ok = TRUE, table = out, times = times, shrunken = TRUE,
       n_curves = length(curves),
       note = paste("Conditional predictions: each curve is its cell's fixed-effect",
                    "trajectory plus that participant's conditional modes (BLUPs),",
                    "so it is SHRUNKEN toward the cell curve by an amount that depends",
                    "on how much data the participant contributed."))
}

# One row per curve: the cell, the shrunken level and trend, and each
# harmonic's shrunken amplitude and acrophase, with the model's own
# between-participant SD beside them so the spread of the rows is not read as
# the spread of the population.
dance_traj_participant_table <- function(fit) {
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec; d <- spec$data
  re <- tryCatch(lme4::ranef(fit$model, condVar = FALSE), error = function(e) NULL)
  if (is.null(re)) return(list(ok = FALSE, message = "Conditional modes are unavailable for this fit."))
  grp <- spec$subject %||% "subject"
  cv <- spec$curve_var %||% "curve"; cv <- if (cv %in% names(d)) cv else grp
  R_subj  <- if (grp %in% names(re)) as.data.frame(re[[grp]]) else NULL
  R_curve <- if (cv %in% names(re) && !identical(cv, grp)) as.data.frame(re[[cv]]) else NULL
  fx <- dance_traj_cell_functionals(fit)
  bv <- dance_traj_beta(fit)
  beta <- bv$beta[!is.na(bv$beta)]
  dot <- function(row) { nm <- intersect(names(row), names(beta)); sum(row[nm] * beta[nm]) }
  dts <- spec$design_terms
  key <- data.frame(subject = as.character(d[[grp]]), curve = as.character(d[[cv]]),
                    cell = if (length(dts))
                      do.call(paste, c(lapply(dts, function(f) as.character(d[[f]])),
                                       list(sep = " x "))) else "(all)",
                    stringsAsFactors = FALSE)
  key <- key[!duplicated(key$curve), , drop = FALSE]
  ci <- match(key$cell, fx$grid$.cell)
  if (anyNA(ci)) return(list(ok = FALSE, message = "A curve's design cell is not among the fitted cells."))
  # a term's conditional mode for this curve, summed over the levels it lives at
  mode_of <- function(term, k) {
    v <- 0
    if (!is.null(R_subj) && term %in% names(R_subj)) {
      si <- match(key$subject[k], rownames(R_subj)); if (!is.na(si)) v <- v + R_subj[[term]][si] }
    if (!is.null(R_curve) && term %in% names(R_curve)) {
      qi <- match(key$curve[k], rownames(R_curve)); if (!is.na(qi)) v <- v + R_curve[[term]][qi] }
    v
  }
  random_terms <- unique(c(names(R_subj), names(R_curve)))
  rows <- lapply(seq_len(nrow(key)), function(k) {
    cell <- fx$cells[[ci[k]]]
    out <- list(subject = key$subject[k], curve = key$curve[k], cell = key$cell[k])
    out$intercept <- dot(cell$intercept) + mode_of("(Intercept)", k)
    out$level_at_t0 <- dot(cell$level0) + mode_of("(Intercept)", k) +
      sum(vapply(spec$harm_terms[grepl("^c", spec$harm_terms)],
                 function(tm) mode_of(tm, k), numeric(1)))
    for (tm in spec$trend_terms) out[[tm]] <- dot(cell$coef[[tm]]) + mode_of(tm, k)
    for (h in seq_len(spec$n_harmonics)) {
      ck <- paste0("c", h); sk <- paste0("s", h)
      a <- dot(cell$coef[[ck]]) + mode_of(ck, k); b <- dot(cell$coef[[sk]]) + mode_of(sk, k)
      kk <- (spec$period / h) / (2 * pi)
      out[[paste0("beta_cos_", h)]] <- a; out[[paste0("beta_sin_", h)]] <- b
      out[[paste0("amplitude_", h)]] <- sqrt(a^2 + b^2)
      out[[paste0("acrophase_rad_", h)]] <- atan2(b, a) %% (2 * pi)
      out[[paste0("acrophase_time_", h)]] <- (atan2(b, a) %% (2 * pi)) * kk
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  # which of these actually vary per participant: a term with no random effect
  # is the cell value for everyone, and the column must say so
  varies <- c(intercept = "(Intercept)" %in% random_terms,
              stats::setNames(spec$trend_terms %in% random_terms, spec$trend_terms),
              stats::setNames(vapply(seq_len(spec$n_harmonics), function(h)
                any(paste0(c("c", "s"), h) %in% random_terms), logical(1)),
                paste0("H", seq_len(spec$n_harmonics))))
  list(ok = TRUE, table = tab, shrunken = TRUE, varies = varies,
       re_label = fit$re_label, random_terms = random_terms,
       note = paste(
         "Shrunken (conditional-mode) estimates from the ONE mixed model: each row",
         "is its cell's fixed effects plus that participant's deviation. A quantity",
         "with no random effect in the fitted structure does not vary per",
         "participant here -- it is the cell value for everyone -- and is marked.",
         "Do not run a second-stage test on these rows; the model's own contrasts",
         "in tab 6 are the group comparison."))
}

# ------------------------------------------------------------------------------
# VARIANCE EXPLAINED BY THE ONE MODEL
# ------------------------------------------------------------------------------
# Nakagawa & Schielzeth (2013) R-squared for a Gaussian mixed model, with the
# random-slope extension of Johnson (2014): the random-effect variance is the
# MEAN over observations of z_i' Sigma z_i, not the sum of the diagonal of Sigma,
# which is wrong as soon as a slope is in the random part -- and every rung of
# this app's ladder has one.
#
#   marginal    = var(X beta) / (var(X beta) + sigma2_random + sigma2_e)
#   conditional = (var(X beta) + sigma2_random) / (same denominator)
#
# Implemented here rather than through the performance / MuMIn packages
# because neither is a declared dependency; the formula is the published one
# and tests/testthat/test-traj-r2.R checks it against a hand computation.
dance_traj_r2 <- function(fit) {
  if (!isTRUE(fit$ok)) return(NULL)
  m <- fit$model
  if (!inherits(m, "merMod")) return(NULL)
  Xb <- tryCatch(as.numeric(stats::predict(m, re.form = NA)), error = function(e) NULL)
  if (is.null(Xb)) return(NULL)
  var_f <- stats::var(Xb)
  Z <- tryCatch(as.matrix(lme4::getME(m, "Z")), error = function(e) NULL)
  vc <- tryCatch(lme4::VarCorr(m), error = function(e) NULL)
  if (is.null(Z) || is.null(vc)) return(NULL)
  # the full random-effect covariance, block-diagonal over grouping levels, in
  # the column order of Z: lme4 orders Z by grouping factor then level then term
  flist <- lme4::getME(m, "flist"); cnms <- lme4::getME(m, "cnms")
  blocks <- list()
  for (g in names(cnms)) {
    Sg <- as.matrix(vc[[g]])[cnms[[g]], cnms[[g]], drop = FALSE]
    nl <- nlevels(flist[[g]])
    blocks[[g]] <- list(S = Sg, nl = nl)
  }
  # z_i' Sigma z_i for every observation, block by block
  var_r <- 0; col <- 0L
  for (g in names(cnms)) {
    b <- blocks[[g]]; q <- ncol(b$S)
    for (l in seq_len(b$nl)) {
      Zl <- Z[, col + seq_len(q), drop = FALSE]
      var_r <- var_r + rowSums((Zl %*% b$S) * Zl)
      col <- col + q
    }
  }
  var_r <- mean(var_r)
  var_e <- stats::sigma(m)^2
  den <- var_f + var_r + var_e
  list(marginal = var_f / den, conditional = (var_f + var_r) / den,
       var_fixed = var_f, var_random = var_r, var_resid = var_e,
       method = "Nakagawa & Schielzeth (2013) with Johnson's (2014) random-slope extension")
}

# ------------------------------------------------------------------------------
# THE FITTED EQUATION OF EVERY DESIGN CELL, from the fixed effects
# ------------------------------------------------------------------------------
# One row per cell: the constant, the trend coefficient(s), and each harmonic's
# (cos, sin) pair with the amplitude and acrophase read off it -- the numbers
# the summary panel and the publication report print under the symbolic
# equation. They are computed from dance_traj_cell_functionals(), i.e. from the
# same beta map dance_traj_predict() draws the curves from, so the printed
# equation and the drawn curve cannot disagree. Acrophase is in model-elapsed
# hours here (t = 0 at the first observation); callers convert to clock time.
dance_traj_cell_equations <- function(fit) {
  if (!isTRUE(fit$ok)) return(NULL)
  spec <- fit$spec
  fx <- dance_traj_cell_functionals(fit)
  bv <- dance_traj_beta(fit); beta <- bv$beta; beta[is.na(beta)] <- 0
  dotb <- function(row) { nm <- intersect(names(row), names(beta)); sum(row[nm] * beta[nm]) }
  cells_n <- spec$cells
  rows <- lapply(seq_len(nrow(fx$grid)), function(i) {
    cl <- fx$cells[[i]]; cname <- fx$grid$.cell[i]
    co <- vapply(fx$basis, function(j) dotb(cl$coef[[j]]), numeric(1)); names(co) <- fx$basis
    ni <- if (!is.null(cells_n)) match(cname, cells_n$cell) else NA_integer_
    out <- list(cell = cname,
                n_participants = if (!is.na(ni)) cells_n$n_participants[ni] else NA_integer_,
                n_obs = if (!is.na(ni)) cells_n$n_obs[ni] else NA_integer_,
                intercept = dotb(cl$intercept), level_at_t0 = dotb(cl$level0))
    for (tm in spec$trend_terms) out[[tm]] <- unname(co[[tm]])
    for (h in seq_len(spec$n_harmonics)) {
      a <- unname(co[[paste0("c", h)]]); b <- unname(co[[paste0("s", h)]])
      out[[paste0("beta_cos_", h)]] <- a; out[[paste0("beta_sin_", h)]] <- b
      out[[paste0("amplitude_", h)]] <- sqrt(a^2 + b^2)
      out[[paste0("acrophase_rad_", h)]] <- atan2(b, a) %% (2 * pi)
      out[[paste0("acrophase_time_", h)]] <- (atan2(b, a) %% (2 * pi)) * (spec$period / h) / (2 * pi)
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  # the trend coefficients in the order dance_format_equation() takes them
  trend_coefs <- lapply(seq_len(nrow(tab)), function(i) switch(spec$trend %||% "none",
    linear = tab$trend_lin[i], log = tab$trend_log[i],
    exp_sat = c(tab$trend_sat[i], spec$tau), NULL))
  list(ok = TRUE, table = tab, trend_coefs = trend_coefs,
       trend = spec$trend %||% "none", period = spec$period,
       n_harmonics = spec$n_harmonics, tau = spec$tau)
}
