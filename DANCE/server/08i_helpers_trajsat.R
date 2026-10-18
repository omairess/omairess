# ==============================================================================
# server/08i_helpers_trajsat.R -- the saturating trend PER DESIGN CELL:
# the rise A, the time constant tau and the initial slope A/tau
# ==============================================================================
# The trajectory model (08d-08g) is linear in its coefficients, so it holds tau
# at ONE value for every cell -- set by the user, or profiled -- and estimates
# only the rise A per cell. To compare tau itself between cells, and the initial
# slope A/tau that depends on it, both have to be estimated per cell, jointly,
# and that makes the model nonlinear. It is nlme's model (Pinheiro & Bates
# 2000, chapters 6-8):
#
#   y = b0[c] + A[c] * (1 - exp(-t / exp(lt[c])))
#       + sum_h a_h[c] cos(w_h t) + b_h[c] sin(w_h t)  + participant effects + e
#
# with every parameter per design cell c (cell-means coding), t in hours from
# the first observation, and tau on the log scale (lt), which keeps it positive
# and gives it an interval that is asymmetric in hours, as a time constant's
# should be. The same model with ONE lt is fitted alongside, by maximum
# likelihood on the same random structure, and the two are compared by a
# likelihood-ratio test -- the direct answer to "is one tau enough?".
#
# Inference, all of it standard:
#   A        Wald, from the fixed-effect covariance of the fit
#   tau      Wald on log tau, back-transformed; pairwise as a RATIO tau_i/tau_j
#   A/tau    the initial slope (dy/dt of the trend at t = 0), by the delta
#            method from the joint covariance of (A, log tau)
#   curves   A (1 - exp(-t/tau)) per cell, delta-method pointwise band
#
# Degrees of freedom: participants minus between-participant cells -- the
# between-within rule lme applies to terms that do not vary within a
# participant. nlme itself gives every coefficient of a nonlinear model the
# WITHIN-participant df (597 on 48 participants x 14 times), which is far too
# generous for a comparison between groups of participants.
# ==============================================================================

DANCE_TRAJ_SAT_PARS <- c("A_sat", "tau", "init_slope")

dance_traj_sat_label <- function(what, units = NULL) {
  u <- if (is.null(units) || !nzchar(units)) "response units" else units
  switch(what,
    A_sat      = sprintf("Saturating rise A (%s above the intercept)", u),
    tau        = "Time constant tau (hours)",
    init_slope = sprintf("Initial slope A/tau (%s per hour)", u),
    what)
}

# The random structures tried, richest first; the first that converges for the
# per-cell model is used for BOTH models, so the likelihood-ratio test compares
# fixed effects only. Diagonal covariances: a nonlinear model with a full
# unstructured covariance over six random effects rarely converges, and a
# diagonal one is what Pinheiro & Bates start from.
dance_traj_sat_ladder <- function(harm_pars, curve_level = FALSE) {
  rhs <- function(p) stats::as.formula(paste(paste(p, collapse = " + "), "~ 1"))
  subj <- list(
    list(label = "participant level, rise and rhythm (diagonal)", pars = c("b0", "A", harm_pars)),
    list(label = "participant level and rise (diagonal)", pars = c("b0", "A")),
    list(label = "participant level only", pars = "b0"))
  out <- lapply(subj, function(r) {
    if (curve_level)
      list(label = paste0(r$label, ", plus a curve-specific level"),
           random = list(subject = nlme::pdDiag(rhs(r$pars)), curve = nlme::pdDiag(b0 ~ 1)),
           groups = ~ subject/curve)
    else list(label = r$label, random = nlme::pdDiag(rhs(r$pars)), groups = ~ subject)
  })
  if (curve_level)
    out[[length(out) + 1L]] <- list(label = "participant level only",
                                    random = nlme::pdDiag(b0 ~ 1), groups = ~ subject)
  out
}

# The fixed-effect name nlme gives parameter `par` in cell k. "tau" names the
# parameter that carries tau, whichever scale it is on.
dance_traj_sat_coef <- function(sf, par, k) {
  if (identical(par, "tau")) par <- sf$tau_map$par %||% "lt"
  if (sf$n_cells > 1L) paste0(par, ".cell_g", k) else par
}

# HOW tau IS CARRIED. Without a range, on the log scale: tau = exp(lt). With a
# range [lo, hi] the search is CONSTRAINED to it by carrying tau on a logistic
# scale, tau = lo + (hi - lo) * plogis(et) -- the standard way to bound a
# parameter of a model that has no box constraints. Both maps are monotone, so
# "tau equal across cells" is the same hypothesis on either scale; the tests and
# intervals below go through the map's derivative.
dance_traj_sat_tau_map <- function(range = NULL) {
  if (is.null(range) || length(range) != 2L || !all(is.finite(range)) ||
      range[1] <= 0 || range[2] <= range[1])
    return(list(par = "lt", expr = "exp(lt)", range = NULL,
                tau = function(x) exp(x), dtau = function(x) exp(x),
                start = function(tau) log(tau)))
  lo <- range[1]; w <- range[2] - range[1]
  list(par = "et", range = range,
       expr = sprintf("(%s + %s * plogis(et))", format(lo, digits = 15), format(w, digits = 15)),
       tau = function(x) lo + w * stats::plogis(x),
       dtau = function(x) w * stats::dlogis(x),
       # a start inside the range, away from its ends where the scale flattens
       start = function(tau) stats::qlogis(min(max((tau - lo) / w, 0.05), 0.95)))
}

dance_traj_sat_fit <- function(ff, progress = NULL, tau_range = NULL,
                               control = list(maxIter = 200, pnlsMaxIter = 20,
                                              msMaxIter = 200),
                               phi = ff$phi %||% 0) {
  if (!requireNamespace("nlme", quietly = TRUE))
    return(list(ok = FALSE, message = paste(
      "The per-cell saturating model needs the nlme package, which ships with R",
      "but is not installed here.")))
  if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
  sp <- ff$spec
  if (!identical(sp$trend, "exp_sat"))
    return(list(ok = FALSE, message = "The per-cell model applies to the saturating trend only."))
  grid <- dance_traj_cell_grid(sp)
  G <- nrow(grid)
  d <- as.data.frame(sp$data)
  key <- if (length(sp$design_terms))
    do.call(paste, c(lapply(sp$design_terms, function(f) as.character(d[[f]])),
                     list(sep = " x "))) else rep(grid$.cell[1], nrow(d))
  k <- match(key, grid$.cell)
  if (anyNA(k) || length(unique(k)) < G)
    return(list(ok = FALSE, message = "Every design cell needs data for the per-cell model."))
  d$cell_ <- factor(paste0("g", k), levels = paste0("g", seq_len(G)))
  d$tt <- d$t - sp$t0
  K <- sp$n_harmonics
  harm <- as.vector(rbind(paste0("a", seq_len(K)), paste0("b", seq_len(K))))
  map <- dance_traj_sat_tau_map(tau_range)
  form <- stats::as.formula(paste(
    sprintf("y ~ b0 + A * (1 - exp(-tt / %s))", map$expr),
    paste(sprintf("+ a%d * c%d + b%d * s%d", seq_len(K), seq_len(K), seq_len(K), seq_len(K)),
          collapse = " ")))
  cell_rhs <- if (G > 1L) "~ 0 + cell_" else "~ 1"
  covs <- intersect(sp$covariates %||% character(0), names(d))
  fx_list <- function(lt_rhs) c(
    list(stats::as.formula(paste("b0", paste(c(cell_rhs, covs), collapse = " + ")))),
    list(stats::as.formula(paste("A", cell_rhs))),
    list(stats::as.formula(paste(map$par, lt_rhs))),
    lapply(harm, function(p) stats::as.formula(paste(p, cell_rhs))))

  # starting values: the trajectory model's own cell values at its tau
  fx <- dance_traj_cell_functionals(ff)
  bv <- dance_traj_beta(ff)
  have <- names(bv$beta)[!is.na(bv$beta)]
  val <- function(r) { nm <- intersect(names(r), have); sum(r[nm] * bv$beta[nm]) }
  per <- function(f) vapply(fx$cells, f, numeric(1))
  tau0 <- sp$tau
  # covariates enter the level b0 only, as in the trajectory model; they start
  # at zero, one value per model-matrix column after the cell columns
  n_cov <- if (!length(covs)) 0L else
    ncol(stats::model.matrix(stats::as.formula(paste(
      cell_rhs, "+", paste(covs, collapse = " + "))), d)) -
    (if (G > 1L) G else 1L)
  # the order nlme expects: parameter by parameter, cells within each
  st_harm <- unlist(lapply(seq_len(K), function(h) list(
    per(function(c) val(c$coef[[paste0("c", h)]])),
    per(function(c) val(c$coef[[paste0("s", h)]])))), use.names = FALSE)
  st_b0 <- c(per(function(c) val(c$intercept)), rep(0, n_cov))
  st_A <- per(function(c) val(c$coef[["trend_sat"]]))
  start_cell <- c(st_b0, st_A, rep(map$start(tau0), G), st_harm)
  start_common <- c(st_b0, st_A, map$start(tau0), st_harm)

  ctl <- do.call(nlme::nlmeControl, control)
  # AUTOCORRELATED RESIDUALS: when the trajectory model has CAR(1) residuals,
  # this model has them too (nlme's corCAR1, within each curve), ESTIMATED with
  # everything else and started at the trajectory model's phi -- in the per-cell
  # and the one-tau fit alike, so their likelihood ratio compares the same
  # residual model. Not held fixed: nlme 3.1-164 crashes R itself (an invalid
  # pointer in its C code) when a corCAR1 is held fixed inside nlme().
  car1 <- length(phi) == 1L && is.finite(phi) && phi > 0
  cor_for <- function(rung) if (!car1) NULL else
    nlme::corCAR1(value = min(max(phi, 0.01), 0.99),
                  form = if (grepl("curve", paste(deparse(rung$groups), collapse = "")))
                           ~ tt | subject/curve
                         else ~ tt | subject)
  fit_one <- function(fixed, start, rung) tryCatch({
    m <- nlme::nlme(form, data = d, fixed = fixed, random = rung$random,
                    groups = rung$groups, start = start, method = "ML", control = ctl,
                    correlation = cor_for(rung))
    v <- diag(stats::vcov(m))
    if (any(!is.finite(v)) || any(v < 0)) stop("the fixed-effect covariance is not positive definite")
    m
  }, error = function(e) structure(conditionMessage(e), class = "dance_sat_err"))

  ladder <- dance_traj_sat_ladder(harm, isTRUE(sp$has_curve_level))
  attempts <- list()
  for (i in seq_along(ladder)) {
    rung <- ladder[[i]]
    if (is.function(progress)) progress(sprintf("per-cell model, %s", rung$label))
    m <- fit_one(fx_list(cell_rhs), start_cell, rung)
    if (inherits(m, "dance_sat_err")) {
      attempts[[i]] <- list(rung = i, label = rung$label, message = as.character(m))
      next
    }
    if (is.function(progress)) progress("one-tau model, for the likelihood-ratio test")
    m1 <- if (G > 1L) fit_one(fx_list("~ 1"), start_common, rung) else NULL
    lrt <- if (G < 2L) list(ok = FALSE, message = "One cell: there is no tau to compare.")
      else if (inherits(m1, "dance_sat_err"))
        list(ok = FALSE, message = paste("The one-tau model did not converge on the same",
                                         "random structure:", as.character(m1)))
      else {
        stat <- 2 * (as.numeric(stats::logLik(m)) - as.numeric(stats::logLik(m1)))
        list(ok = TRUE, statistic = max(0, stat), df = G - 1L,
             p = stats::pchisq(max(0, stat), G - 1L, lower.tail = FALSE),
             tau_common = map$tau(unname(nlme::fixef(m1)[[map$par]])))
      }
    between <- intersect(sp$between_terms %||% character(0), sp$design_terms)
    n_subj <- length(unique(as.character(d$subject)))
    n_between <- if (length(between))
      nrow(unique(d[, between, drop = FALSE])) else 1L
    out <- list(ok = TRUE, fit = m, fit_common = if (inherits(m1, "dance_sat_err")) NULL else m1,
                cells = grid$.cell, n_cells = G, grid = grid, rung = i, n_rungs = length(ladder),
                re_label = rung$label, attempts = attempts, lrt = lrt,
                df = max(1L, n_subj - n_between), n_participants = n_subj,
                window = range(d$tt), t0 = sp$t0, period = sp$period,
                tau_traj = sp$tau, method = "ML", tau_map = map, tau_range = map$range,
                # the residual autocorrelation this model estimated (0: none modelled)
                phi = if (car1) tryCatch(as.numeric(stats::coef(m$modelStruct$corStruct,
                                                                 unconstrained = FALSE)),
                                         error = function(e) NA_real_) else 0,
                formula = paste(deparse(form, width.cutoff = 500L), collapse = " "))
    out$table <- dance_traj_sat_table(out)
    return(out)
  }
  # WHY IT FAILED, when a range was set: a cell whose data put tau outside the
  # range drives the constrained estimate into the bound, where the logistic
  # scale goes flat and the fit cannot converge. One unconstrained fit tells the
  # two failures apart, and names the cells -- as a diagnosis, not as a result.
  outside <- NULL
  if (!is.null(map$range)) {
    free <- tryCatch(dance_traj_sat_fit(ff, progress = progress, tau_range = NULL,
                                        control = control, phi = phi), error = function(e) NULL)
    if (isTRUE(free$ok)) {
      ft <- free$table
      out_i <- which(ft$tau < map$range[1] | ft$tau > map$range[2])
      if (length(out_i))
        outside <- sprintf(paste(
          "Without the range the model converges, and puts tau at %s -- outside %s to %s h,",
          "so the constrained fit runs into the bound. Widen the range, or hold tau at a value."),
          paste(sprintf("%.3g h in %s", ft$tau[out_i], ft$cell[out_i]), collapse = ", "),
          format(map$range[1]), format(map$range[2]))
    }
  }
  list(ok = FALSE, attempts = attempts, message = paste(
    "The per-cell saturating model did not converge on any random structure tried",
    sprintf("(%s).", paste(vapply(attempts, function(a) a$message, character(1)),
                           collapse = " | ")),
    outside %||% paste(
      "A and tau are then not separately estimable per cell from these data; compare",
      "the non-periodic change instead, which the trajectory model estimates.")))
}

# The fixed-effect covariance as nlme itself reports it for a maximum-likelihood
# fit: vcov() times N / (N - p), nlme's adjustSigma, which turns the ML residual
# variance into a REML-like one. summary() and anova(..., L = ) both apply it,
# so the standard errors and F statistics here are nlme's own.
dance_traj_sat_V <- function(fit) {
  V <- stats::vcov(fit)
  N <- fit$dims$N; p <- length(nlme::fixef(fit))
  if (identical(fit$method, "ML") && N > p) V * N / (N - p) else V
}

# Per cell: A, tau and A/tau with their intervals, and the correlation of A with
# log tau -- the number that says whether the two are separately identified.
dance_traj_sat_table <- function(sf, conf = 0.95) {
  fe <- nlme::fixef(sf$fit); V <- dance_traj_sat_V(sf$fit)
  q <- stats::qt(1 - (1 - conf) / 2, sf$df)
  map <- sf$tau_map %||% dance_traj_sat_tau_map(NULL)
  rows <- lapply(seq_len(sf$n_cells), function(k) {
    na <- dance_traj_sat_coef(sf, "A", k); nl <- dance_traj_sat_coef(sf, "tau", k)
    A <- unname(fe[na]); x <- unname(fe[nl]); tau <- map$tau(x); s <- A / tau
    g <- map$dtau(x)                                 # dtau/dx
    vA <- V[na, na]; vx <- V[nl, nl]; cAx <- V[na, nl]
    # A/tau through the delta method: gradient (1/tau, -s * g / tau) in (A, x)
    se_s <- sqrt(max(0, vA / tau^2 + (s * g / tau)^2 * vx - 2 * (s * g / tau^2) * cAx))
    data.frame(cell = sf$cells[k], n_cells = sf$n_cells,
               A = A, A_se = sqrt(vA), A_lo = A - q * sqrt(vA), A_hi = A + q * sqrt(vA),
               tau = tau, tau_se = g * sqrt(vx),
               # the interval on the carrying scale, mapped back: inside the range
               tau_lo = map$tau(x - q * sqrt(vx)), tau_hi = map$tau(x + q * sqrt(vx)),
               init_slope = s, init_slope_se = se_s,
               init_slope_lo = s - q * se_s, init_slope_hi = s + q * se_s,
               corr_A_logtau = cAx / sqrt(vA * vx),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  tab
}

# The parameter vector over cells and its covariance, on the scale the tests
# use: A as is, tau as log tau, A/tau through the delta method.
dance_traj_sat_vector <- function(sf, what) {
  fe <- nlme::fixef(sf$fit); V <- dance_traj_sat_V(sf$fit)
  G <- sf$n_cells
  map <- sf$tau_map %||% dance_traj_sat_tau_map(NULL)
  nA <- vapply(seq_len(G), function(k) dance_traj_sat_coef(sf, "A", k), character(1))
  nl <- vapply(seq_len(G), function(k) dance_traj_sat_coef(sf, "tau", k), character(1))
  x <- unname(fe[nl]); tau <- map$tau(x); g <- map$dtau(x)
  switch(what,
    A_sat = list(theta = unname(fe[nA]), V = unname(V[nA, nA, drop = FALSE]), scale = "identity"),
    # log tau, whichever scale carried it: d log tau / dx = g / tau
    tau   = {
      J <- diag(g / tau, nrow = G)
      list(theta = log(tau), V = J %*% V[nl, nl, drop = FALSE] %*% t(J), scale = "log")
    },
    init_slope = {
      A <- unname(fe[nA]); s <- A / tau
      J <- matrix(0, G, 2L * G)
      for (k in seq_len(G)) { J[k, k] <- 1 / tau[k]; J[k, G + k] <- -s[k] * g[k] / tau[k] }
      Vs <- J %*% V[c(nA, nl), c(nA, nl)] %*% t(J)
      list(theta = s, V = Vs, scale = "identity")
    },
    stop(sprintf("unknown saturating parameter '%s'", what)))
}

# A Wald F over the cells for one parameter: every cell equal, or one factorial
# effect, through the same marginal contrasts as the trajectory model's tests.
# The denominator df is the participant count minus the between-participant
# cells for every contrast: the participants are the replication units, whether
# a contrast compares groups of people or conditions within the same people.
dance_traj_sat_test <- function(sf, ff, what, effect_terms = character(0)) {
  C <- dance_traj_effect_contrasts(ff, effect_terms)
  if (is.null(C) || !nrow(C)) return(list(ok = FALSE, message = "Nothing to compare."))
  v <- dance_traj_sat_vector(sf, what)
  est <- C %*% v$theta
  CV <- C %*% v$V %*% t(C)
  qr_ <- qr(CV); q <- qr_$rank
  if (q < 1L) return(list(ok = FALSE, message = "The contrast has no variance."))
  W <- tryCatch(as.numeric(t(est) %*% MASS::ginv(CV) %*% est), error = function(e) NA_real_)
  if (!is.finite(W)) return(list(ok = FALSE, message = "The contrast covariance is singular."))
  Fv <- W / q
  list(ok = TRUE, what = what, effect = effect_terms, statistic = Fv, df1 = q, df2 = sf$df,
       p = stats::pf(Fv, q, sf$df, lower.tail = FALSE),
       method = sprintf("Wald F%s", if (identical(what, "init_slope")) ", delta method" else
                          if (identical(what, "tau")) " on log tau" else ""))
}

# Pairwise between cells, optionally inside a subset of cells (a slice of a
# factorial design). tau is compared as a ratio; A and A/tau as differences.
dance_traj_sat_pairwise <- function(sf, what, adjust = c("holm", "bonferroni", "none"),
                                    conf = 0.95, cells = NULL) {
  adjust <- match.arg(adjust)
  v <- dance_traj_sat_vector(sf, what)
  idx <- if (is.null(cells)) seq_len(sf$n_cells) else match(cells, sf$cells)
  idx <- idx[!is.na(idx)]
  if (length(idx) < 2L) return(list(ok = FALSE, message = "Fewer than two cells to compare."))
  pr <- utils::combn(idx, 2)
  q <- stats::qt(1 - (1 - conf) / 2, sf$df)
  rows <- lapply(seq_len(ncol(pr)), function(j) {
    i1 <- pr[1, j]; i2 <- pr[2, j]
    d <- v$theta[i1] - v$theta[i2]
    se <- sqrt(max(0, v$V[i1, i1] + v$V[i2, i2] - 2 * v$V[i1, i2]))
    tval <- d / se
    ratio <- identical(v$scale, "log")
    data.frame(cell1 = sf$cells[i1], cell2 = sf$cells[i2],
               estimate = if (ratio) exp(d) else d, se = se,
               lo = if (ratio) exp(d - q * se) else d - q * se,
               hi = if (ratio) exp(d + q * se) else d + q * se,
               statistic = tval, df2 = sf$df,
               p_raw = 2 * stats::pt(-abs(tval), sf$df), defined = is.finite(tval),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab$p_adj <- if (identical(adjust, "none")) tab$p_raw else stats::p.adjust(tab$p_raw, adjust)
  rownames(tab) <- NULL
  list(ok = TRUE, what = what, table = tab, adjust = adjust, conf = conf,
       ratio = identical(v$scale, "log"), df = sf$df)
}

# The trend part per cell, A (1 - exp(-t/tau)), on elapsed time from the first
# observation, with a pointwise delta-method band.
dance_traj_sat_curves <- function(sf, n_time = 160, conf = 0.95) {
  fe <- nlme::fixef(sf$fit); V <- dance_traj_sat_V(sf$fit)
  tt <- seq(sf$window[1], sf$window[2], length.out = n_time)
  q <- stats::qt(1 - (1 - conf) / 2, sf$df)
  map <- sf$tau_map %||% dance_traj_sat_tau_map(NULL)
  out <- lapply(seq_len(sf$n_cells), function(k) {
    na <- dance_traj_sat_coef(sf, "A", k); nl <- dance_traj_sat_coef(sf, "tau", k)
    x <- unname(fe[nl]); A <- unname(fe[na]); tau <- map$tau(x)
    ex <- exp(-tt / tau)
    fit <- A * (1 - ex)
    # d fit / d x = d fit / d tau * d tau / d x = -A ex t / tau^2 * dtau
    g1 <- 1 - ex; g2 <- -A * ex * tt / tau^2 * map$dtau(x)
    Vk <- V[c(na, nl), c(na, nl)]
    se <- sqrt(pmax(0, g1^2 * Vk[1, 1] + g2^2 * Vk[2, 2] + 2 * g1 * g2 * Vk[1, 2]))
    data.frame(cell = sf$cells[k], tt = tt, fit = fit, se = se,
               lo = fit - q * se, hi = fit + q * se, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(table = tab, cells = sf$cells, conf = conf)
}

# What the table cannot tell you on its own: where A and tau are not separately
# identified, and where tau runs past the end of the recording.
dance_traj_sat_flags <- function(sf, threshold = 0.9) {
  tab <- sf$table
  T <- diff(sf$window)
  out <- character(0)
  hi <- tab$cell[abs(tab$corr_A_logtau) > threshold]
  if (length(hi))
    out <- c(out, sprintf(paste(
      "In %s, A and tau are correlated beyond %.1f: the data pin down the early slope",
      "and the change over the recording far better than A and tau separately, so",
      "read A/tau and the change before A or tau alone."),
      paste(hi, collapse = ", "), threshold))
  rg <- sf$tau_range
  if (!is.null(rg)) {
    edge <- tab$cell[tab$tau < rg[1] + 0.01 * diff(rg) | tab$tau > rg[2] - 0.01 * diff(rg)]
    if (length(edge))
      out <- c(out, sprintf(paste(
        "In %s, tau ran to the end of the search range (%s to %s h): it is a bound, not an",
        "estimate, and its interval is not meaningful. Widen the range, or hold tau at a value."),
        paste(edge, collapse = ", "), format(rg[1]), format(rg[2])))
  }
  long <- tab$cell[tab$tau > T]
  if (length(long))
    out <- c(out, sprintf(paste(
      "In %s, tau exceeds the %.0f h recording: the rise is still far from its",
      "asymptote at the last observation, so A is an extrapolation beyond the data."),
      paste(long, collapse = ", "), T))
  out
}
