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
# ONE FITTED TRAJECTORY PER CELL, WITH A BAND
# ------------------------------------------------------------------------------
dance_traj_predict <- function(fit, times = NULL, conf = 0.95,
                               band = c("pointwise", "simultaneous"),
                               n_time = 200) {
  band <- match.arg(band)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
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
dance_traj_contrasts <- function(fit, what = c("level", "amplitude", "phase"),
                                 harmonic = 1, conf = 0.95, at_time = NULL,
                                 adjust = c("holm", "bonferroni", "none")) {
  what <- match.arg(what); adjust <- match.arg(adjust)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
  spec <- fit$spec
  if (!length(spec$design_terms))
    return(list(ok = FALSE, message = "No design factors: there are no cells to contrast."))
  grid <- dance_traj_cell_grid(spec)
  m <- nrow(grid)
  if (m < 2) return(list(ok = FALSE, message = "Only one design cell."))
  pairs <- utils::combn(m, 2)
  z <- stats::qnorm(1 - (1 - conf) / 2)

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
                 statistic = est / se, defined = TRUE, stringsAsFactors = FALSE)
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
      if (identical(what, "phase")) {
        pc <- dance_traj_phase_contrast(co, i, j, conf)
        if (is.null(pc)) return(NULL)
        return(data.frame(cell1 = co$cells[i], cell2 = co$cells[j],
                          estimate = pc$diff_time, se = pc$se_time,
                          lo = pc$lo, hi = pc$hi,
                          statistic = if (isTRUE(pc$defined)) pc$diff_time / pc$se_time else NA_real_,
                          defined = isTRUE(pc$defined), stringsAsFactors = FALSE))
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
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(list(ok = FALSE, message = "No contrast could be formed."))
    tab <- do.call(rbind, rows)
    unit <- if (identical(what, "phase"))
      sprintf("acrophase difference in %s on an effective period of %.4g",
              spec$time_units %||% "time units", co$effective_period)
    else "amplitude difference in the response's own units"
  }

  tab$p_raw <- ifelse(tab$defined, 2 * stats::pnorm(-abs(tab$statistic)), NA_real_)
  tab$p_adj <- tab$p_raw
  if (!identical(adjust, "none"))
    tab$p_adj[tab$defined] <- stats::p.adjust(tab$p_raw[tab$defined], method = adjust)
  rownames(tab) <- NULL
  list(ok = TRUE, what = what, harmonic = harmonic, unit = unit,
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
         else paste("Amplitude and phase are nonlinear maps of the cell's (cos, sin) pair;",
                    "these are delta-method intervals with the exact cross-cell covariance."),
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
    em <- tryCatch(emmeans::emtrends(
      fit$model,
      specs = stats::as.formula(paste("~", paste(spec$design_terms, collapse = " * "))),
      var = nm), error = function(e) NULL)
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
dance_traj_simple_effects <- function(fit, effect, at, what = c("level", "amplitude", "phase"),
                                      harmonic = 1, conf = 0.95, at_time = NULL,
                                      adjust = c("holm", "bonferroni", "none")) {
  what <- match.arg(what); adjust <- match.arg(adjust)
  if (!isTRUE(fit$ok)) return(list(ok = FALSE, message = fit$message))
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

  all_ct <- dance_traj_contrasts(fit, what, harmonic, conf, at_time, adjust = "none")
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
  R <- as.data.frame(re[[grp]])
  ck <- paste0("c", harmonic); sk <- paste0("s", harmonic)
  if (!all(c(ck, sk) %in% names(R)))
    return(list(ok = FALSE, message = paste(
      "The fitted random structure carries no participant-specific harmonic terms,",
      sprintf("so there are no per-participant rhythms to report (rung used: %s).",
              fit$re_label %||% "unknown"))))
  # the population (fixed-effect) rhythm each participant deviates from
  co <- dance_traj_cell_coefs(fit, harmonic)
  pop_a <- if (isTRUE(co$ok)) mean(co$a) else NA_real_
  pop_b <- if (isTRUE(co$ok)) mean(co$b) else NA_real_
  a <- pop_a + R[[ck]]; b <- pop_b + R[[sk]]
  k <- (spec$period / harmonic) / (2 * pi)
  vc <- tryCatch(as.matrix(lme4::VarCorr(fit$model)[[grp]]), error = function(e) NULL)
  tab <- data.frame(subject = rownames(R), beta_cos = a, beta_sin = b,
                    amplitude = sqrt(a^2 + b^2),
                    acrophase_rad = atan2(b, a) %% (2 * pi),
                    acrophase_time = (atan2(b, a) %% (2 * pi)) * k,
                    stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  list(ok = TRUE, table = tab, harmonic = harmonic, shrunken = TRUE,
       effective_period = spec$period / harmonic,
       re_label = fit$re_label,
       sd_blup_amplitude = stats::sd(tab$amplitude),
       sd_population_amplitude = if (!is.null(vc) && all(c(ck, sk) %in% rownames(vc))) {
         A <- sqrt(pop_a^2 + pop_b^2)
         if (A > 0) { g <- c(pop_a / A, pop_b / A)
           sqrt(max(0, as.numeric(t(g) %*% vc[c(ck, sk), c(ck, sk)] %*% g))) } else NA_real_
       } else NA_real_,
       label = "model-based shrunken participant estimates (conditional modes / BLUPs)",
       note = paste(
         "THESE ARE SHRUNKEN. A conditional mode is pulled toward the population",
         "mean by an amount that depends on how much data that participant",
         "contributed and how large the random-effect variance is. They are better",
         "PREDICTIONS of each individual than independent per-participant fits, but",
         "they are not the same quantity.",
         "In particular the spread of these values UNDERSTATES the spread of true",
         "participant parameters: compare sd_blup_amplitude with",
         "sd_population_amplitude, which is the model's estimate of the real",
         "between-participant SD. A histogram or boxplot of this column must say",
         "so on its face, or it will be read as the population distribution."))
}
