# =============================================================================
# rasch_engine.R  --  A WINSTEPS-style Rasch engine in base R
#
# Implements, in plain R with no modelling-package dependencies:
#   * JMLE (UCON) estimation of the dichotomous Rasch model, the Andrich Rating
#     Scale Model (RSM) and the Masters Partial Credit Model (PCM), with
#     item "grouping" (WINSTEPS ISGROUPS=).
#   * Extreme-score handling with the EXTRSCORE= adjustment.
#   * Model / Real standard errors, INFIT & OUTFIT MNSQ and ZSTD
#     (Wilson-Hilferty), point-measure and point-biserial correlations and
#     their model expectations, post-hoc discrimination.
#   * Separation / strata / reliability (WINSTEPS Table 3.1).
#   * Category ("rating scale") structure statistics (Table 3.2 / 3.3):
#     Andrich thresholds, observed & expected average measures, category
#     outfit/infit, category measures, Rasch-Thurstone 50% cumulative
#     thresholds, coherence M->C and C->M.
#   * Score-to-measure table (Table 20), test characteristic and information
#     curves.
#   * DIF: Rasch-Welch iterative-logit t-test + Mantel-Haenszel (dichotomies)
#     / Mantel (polytomies), stratified by measure (Table 30).
#   * PCA of standardized residuals with the full WINSTEPS variance
#     decomposition (Table 23).
#
# ---------------------------------------------------------------------------
# HONEST SCOPE NOTE. This reproduces the *published* WINSTEPS algorithms
# (Linacre's Winsteps Help; Wright & Masters 1982; Wright & Stone 1979) and in
# testing recovers the same numbers to within rounding for well-behaved data.
# It is NOT byte-identical to WINSTEPS 5.11: WINSTEPS uses proprietary
# convergence acceleration, its own extreme-score and missing-data conventions,
# and optional JMLE bias correction (STBIAS=). Treat small last-decimal
# differences as expected, and report the engine you actually used.
#
# References
#   Andrich D (1978) Psychometrika 43:561-573.
#   Masters GN (1982) Psychometrika 47:149-174.
#   Wright BD & Masters GN (1982) Rating Scale Analysis. MESA Press.
#   Wright BD & Stone MH (1979) Best Test Design. MESA Press.
#   Smith RM (1991) Journal of Applied Measurement / Rasch Meas Trans 5:2.
#   Linacre JM (1999) J Outcome Measurement 3:103-122 (category quality).
#   Linacre JM (2002) J Applied Measurement 3:85-106 (optimizing categories).
#   Mantel N & Haenszel W (1959) JNCI 22:719-748.  Mantel N (1963) JASA 58:690.
#   Raiche G (2005) Rasch Meas Trans 19:1012 (critical eigenvalues).
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

# ---------------------------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------------------------

#' Validate / recode a response matrix for Rasch analysis.
#'
#' @param X data.frame or matrix of integer responses (persons x items). NA =
#'   not administered.
#' @param groups character/numeric vector of length ncol(X). Items sharing a
#'   group share one set of Andrich thresholds (WINSTEPS ISGROUPS=).
#'   Use a single value for RSM, unique values for PCM.
#' @param recode "shift" subtracts the observed minimum of each group so
#'   categories start at 0; "collapse" additionally renumbers observed
#'   categories consecutively (structural/incidental zeroes removed).
rasch_prep <- function(X, groups = NULL, recode = c("shift", "collapse", "none")) {
  recode <- match.arg(recode)
  X <- as.matrix(X)
  if (!is.numeric(X)) storage.mode(X) <- "numeric"
  I <- ncol(X); N <- nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("I", seq_len(I))
  if (is.null(rownames(X))) rownames(X) <- paste0("P", seq_len(N))
  if (is.null(groups)) groups <- rep("R1", I)
  groups <- as.character(groups)
  stopifnot(length(groups) == I)

  notes <- character(0)
  if (recode != "none") {
    for (g in unique(groups)) {
      cols <- which(groups == g)
      v <- X[, cols, drop = FALSE]
      obs <- sort(unique(as.vector(v[!is.na(v)])))
      if (!length(obs)) next
      if (recode == "shift") {
        if (obs[1] != 0) {
          X[, cols] <- v - obs[1]
          notes <- c(notes, sprintf("Group '%s': subtracted %g so categories start at 0.", g, obs[1]))
        }
      } else { # collapse
        newv <- match(v, obs) - 1L
        if (!identical(as.vector(v), as.vector(obs[newv + 1L])) || !identical(obs, seq(0, length(obs) - 1))) {
          X[, cols] <- newv
          notes <- c(notes, sprintf("Group '%s': categories %s renumbered to %s.",
                                    g, paste(obs, collapse = ","),
                                    paste(seq_along(obs) - 1, collapse = ",")))
        }
      }
    }
  }
  mx <- vapply(unique(groups), function(g)
    suppressWarnings(max(X[, groups == g, drop = FALSE], na.rm = TRUE)), numeric(1))
  names(mx) <- unique(groups)
  if (any(!is.finite(mx))) stop("At least one item group has no observed responses.")

  # non-integer / gapped category check
  for (g in unique(groups)) {
    obs <- sort(unique(as.vector(X[, groups == g, drop = FALSE])))
    obs <- obs[!is.na(obs)]
    if (any(obs != round(obs)))
      stop(sprintf("Group '%s' contains non-integer responses. Recode first.", g))
    if (!identical(as.numeric(obs), as.numeric(seq(0, max(obs)))))
      notes <- c(notes, sprintf(
        "Group '%s' has unobserved intermediate categories (%s). Andrich thresholds for empty categories are not estimable; consider recode = 'collapse'.",
        g, paste(setdiff(seq(0, max(obs)), obs), collapse = ",")))
  }

  list(X = X, groups = groups, max_cat = mx, notes = notes,
       n_persons = N, n_items = I)
}

# ---------------------------------------------------------------------------
# 2. Model probabilities and moments
# ---------------------------------------------------------------------------

#' Category probabilities for one item group.
#' @return list of (m+1) N x I matrices, element k+1 = P(X = k)
.probs_group <- function(theta, delta, tau) {
  m  <- length(tau)
  ct <- c(0, cumsum(tau))                 # cumulative thresholds, length m+1
  td <- outer(theta, delta, "-")          # N x I  (theta_n - delta_i)
  psi <- lapply(0:m, function(k) k * td - ct[k + 1])
  mx  <- Reduce(pmax, psi)
  ep  <- lapply(psi, function(p) exp(p - mx))
  s   <- Reduce(`+`, ep)
  lapply(ep, function(e) e / s)
}

#' Moments of the score distribution from category probabilities.
.moments <- function(P) {
  m  <- length(P) - 1
  E  <- Reduce(`+`, lapply(0:m, function(k) k * P[[k + 1]]))
  W  <- Reduce(`+`, lapply(0:m, function(k) (k - E)^2 * P[[k + 1]]))
  C4 <- Reduce(`+`, lapply(0:m, function(k) (k - E)^4 * P[[k + 1]]))
  list(E = E, W = W, C = C4)
}

#' Full-matrix expected score, model variance and 4th moment.
.model_matrices <- function(theta, delta, tau_list, groups) {
  N <- length(theta); I <- length(delta)
  E <- W <- C <- matrix(NA_real_, N, I)
  for (g in names(tau_list)) {
    cols <- which(groups == g)
    P <- .probs_group(theta, delta[cols], tau_list[[g]])
    mo <- .moments(P)
    E[, cols] <- mo$E; W[, cols] <- mo$W; C[, cols] <- mo$C
  }
  list(E = E, W = W, C = C)
}

#' Cumulative probabilities P(X >= k), k = 1..m, for one group.
.cumprobs_group <- function(theta, delta, tau) {
  P <- .probs_group(theta, delta, tau)
  m <- length(tau)
  lapply(1:m, function(k) Reduce(`+`, P[(k + 1):(m + 1)]))
}

# ---------------------------------------------------------------------------
# 3. Extreme-score detection
# ---------------------------------------------------------------------------

#' Iteratively flag persons / items with minimum-possible or maximum-possible
#' scores on the responses actually observed (WINSTEPS excludes these from
#' estimation and reports them with an EXTRSCORE= adjusted measure).
.find_extremes <- function(X, groups, max_cat) {
  M   <- !is.na(X)
  cap <- matrix(max_cat[groups][col(X)], nrow(X), ncol(X))
  keep_p <- rep(TRUE, nrow(X)); keep_i <- rep(TRUE, ncol(X))
  repeat {
    Mi <- M; Mi[!keep_p, ] <- FALSE; Mi[, !keep_i] <- FALSE
    rs  <- rowSums(X * Mi, na.rm = TRUE); rmax <- rowSums(cap * Mi)
    cs  <- colSums(X * Mi, na.rm = TRUE); cmax <- colSums(cap * Mi)
    np  <- rowSums(Mi);  ni <- colSums(Mi)
    new_p <- keep_p & np > 0 & rs > 0 & rs < rmax
    new_i <- keep_i & ni > 0 & cs > 0 & cs < cmax
    if (identical(new_p, keep_p) && identical(new_i, keep_i)) break
    keep_p <- new_p; keep_i <- new_i
    if (!any(keep_p) || !any(keep_i)) break
  }
  list(keep_p = keep_p, keep_i = keep_i)
}

# ---------------------------------------------------------------------------
# 4. PROX starting values (Cohen 1979 normal approximation)
# ---------------------------------------------------------------------------

.prox_start <- function(X, groups, max_cat) {
  M   <- !is.na(X)
  cap <- matrix(max_cat[groups][col(X)], nrow(X), ncol(X))
  rs  <- rowSums(X * M, na.rm = TRUE); rmax <- rowSums(cap * M)
  cs  <- colSums(X * M, na.rm = TRUE); cmax <- colSums(cap * M)
  p_r <- pmin(pmax(rs / pmax(rmax, 1e-9), .005), .995)
  p_c <- pmin(pmax(cs / pmax(cmax, 1e-9), .005), .995)
  th  <- log(p_r / (1 - p_r))
  de  <- -log(p_c / (1 - p_c))
  th[!is.finite(th)] <- 0; de[!is.finite(de)] <- 0
  list(theta = th - 0, delta = de - mean(de[is.finite(de)]))
}

# ---------------------------------------------------------------------------
# 5. JMLE (UCON)
# ---------------------------------------------------------------------------

#' Joint maximum likelihood estimation of the (polytomous) Rasch model.
#'
#' @param prep output of rasch_prep()
#' @param maxit maximum JMLE iterations (WINSTEPS MJMLE=)
#' @param conv convergence: max absolute parameter change (WINSTEPS LCONV=)
#' @param rconv convergence: max absolute score residual (WINSTEPS RCONV=)
#' @param extreme_adj score adjustment for extreme persons/items
#'   (WINSTEPS EXTRSCORE=, default 0.3)
#' @param max_step maximum logit change per iteration (numerical safeguard)
#' @param center "items" (UIMEAN=0 on non-extreme items) or "persons"
rasch_jmle <- function(prep, maxit = 400, conv = 1e-4, rconv = 1e-3,
                       extreme_adj = 0.3, max_step = 1.0,
                       center = c("items", "persons"), verbose = FALSE) {
  center <- match.arg(center)
  X <- prep$X; groups <- prep$groups; max_cat <- prep$max_cat
  N <- nrow(X); I <- ncol(X)

  ex <- .find_extremes(X, groups, max_cat)
  keep_p <- ex$keep_p; keep_i <- ex$keep_i
  if (sum(keep_p) < 2 || sum(keep_i) < 2)
    stop("Fewer than two non-extreme persons or items remain; the data cannot be calibrated.")

  Xc <- X[keep_p, keep_i, drop = FALSE]
  gc_ <- groups[keep_i]
  gl  <- unique(gc_)
  mcc <- max_cat[gl]
  Mc  <- !is.na(Xc)
  r_p <- rowSums(Xc * Mc, na.rm = TRUE)      # person raw scores (sufficient)
  s_i <- colSums(Xc * Mc, na.rm = TRUE)      # item raw scores  (sufficient)

  st <- .prox_start(Xc, gc_, mcc)
  theta <- st$theta; delta <- st$delta
  tau_list <- lapply(gl, function(g) {
    m <- mcc[[g]]
    if (m < 1) stop("An item group has maximum category 0.")
    if (m == 1) numeric(0) else rep(0, m)    # dichotomies: no free thresholds
  })
  names(tau_list) <- gl
  # NB: for m == 1 the single "threshold" is absorbed into delta, so tau is empty
  tau_list <- lapply(gl, function(g) if (mcc[[g]] == 1) 0 else rep(0, mcc[[g]]))
  names(tau_list) <- gl

  # observed cumulative counts per group (sufficient statistics for tau)
  cum_obs <- lapply(gl, function(g) {
    cols <- which(gc_ == g); m <- mcc[[g]]
    v <- Xc[, cols, drop = FALSE]
    vapply(1:m, function(k) sum(v >= k, na.rm = TRUE), numeric(1))
  })
  names(cum_obs) <- gl

  clamp <- function(d) pmax(pmin(d, max_step), -max_step)
  converged <- FALSE; it <- 0; hist <- numeric(0)

  for (it in seq_len(maxit)) {
    maxchg <- 0

    ## --- persons ---
    mm <- .model_matrices(theta, delta, tau_list, gc_)
    num <- r_p - rowSums(mm$E * Mc, na.rm = TRUE)
    den <- rowSums(mm$W * Mc, na.rm = TRUE)
    d   <- clamp(num / pmax(den, 1e-8))
    theta <- theta + d
    maxchg <- max(maxchg, max(abs(d)))
    max_res <- max(abs(num))

    ## --- items ---
    mm <- .model_matrices(theta, delta, tau_list, gc_)
    num <- colSums(mm$E * Mc, na.rm = TRUE) - s_i
    den <- colSums(mm$W * Mc, na.rm = TRUE)
    d   <- clamp(num / pmax(den, 1e-8))
    delta <- delta + d
    maxchg  <- max(maxchg, max(abs(d)))
    max_res <- max(max_res, max(abs(num)))

    ## --- Andrich thresholds (per group) ---
    for (g in gl) {
      m <- mcc[[g]]
      if (m < 2) next                          # dichotomies have none
      cols <- which(gc_ == g)
      Pk   <- .cumprobs_group(theta, delta[cols], tau_list[[g]])
      msk  <- Mc[, cols, drop = FALSE]
      Eobs <- vapply(Pk, function(p) sum(p * msk), numeric(1))
      Vr   <- vapply(Pk, function(p) sum(p * (1 - p) * msk), numeric(1))
      dtau <- clamp((Eobs - cum_obs[[g]]) / pmax(Vr, 1e-8))
      tau_new <- tau_list[[g]] + dtau
      cc <- mean(tau_new)
      tau_list[[g]] <- tau_new - cc            # sum(tau) = 0 identification
      delta[cols]   <- delta[cols] + cc        # reparametrisation: fit invariant
      maxchg  <- max(maxchg, max(abs(dtau)))
      max_res <- max(max_res, max(abs(Eobs - cum_obs[[g]])))
    }

    ## --- centering (fit-invariant shift) ---
    a <- if (center == "items") mean(delta) else mean(theta)
    delta <- delta - a; theta <- theta - a

    hist <- c(hist, maxchg)
    if (verbose && it %% 10 == 0)
      message(sprintf("iter %3d  max|change| = %.6f  max|residual| = %.5f", it, maxchg, max_res))
    if (maxchg < conv && max_res < rconv) { converged <- TRUE; break }
  }

  ## --- standard errors for calibrated parameters ---
  mm <- .model_matrices(theta, delta, tau_list, gc_)
  se_theta_c <- 1 / sqrt(pmax(rowSums(mm$W * Mc, na.rm = TRUE), 1e-8))
  se_delta_c <- 1 / sqrt(pmax(colSums(mm$W * Mc, na.rm = TRUE), 1e-8))

  ## --- SE of Andrich thresholds ---
  se_tau <- lapply(gl, function(g) {
    m <- mcc[[g]]
    if (m < 2) return(numeric(0))
    cols <- which(gc_ == g)
    Pk  <- .cumprobs_group(theta, delta[cols], tau_list[[g]])
    msk <- Mc[, cols, drop = FALSE]
    1 / sqrt(vapply(Pk, function(p) max(sum(p * (1 - p) * msk), 1e-8), numeric(1)))
  })
  names(se_tau) <- gl

  ## --- place calibrated values back into full-length vectors ---
  Theta <- rep(NA_real_, N); Theta[keep_p] <- theta
  Delta <- rep(NA_real_, I); Delta[keep_i] <- delta
  SEt   <- rep(NA_real_, N); SEt[keep_p]   <- se_theta_c
  SEd   <- rep(NA_real_, I); SEd[keep_i]   <- se_delta_c

  ## --- extreme persons: anchored estimation with adjusted score ---
  Mfull <- !is.na(X)
  cap   <- matrix(max_cat[groups][col(X)], N, I)
  if (any(!keep_p)) {
    for (n in which(!keep_p)) {
      obs <- Mfull[n, ] & !is.na(Delta)
      if (!any(obs)) next
      rmax <- sum(cap[n, obs]); r <- sum(X[n, obs])
      radj <- min(max(r, extreme_adj), rmax - extreme_adj)
      es <- .solve_measure_person(radj, Delta[obs], tau_list, groups[obs])
      Theta[n] <- es$measure; SEt[n] <- es$se
    }
  }
  ## --- extreme items: anchored estimation with adjusted score ---
  if (any(!keep_i)) {
    for (i in which(!keep_i)) {
      obs <- Mfull[, i] & !is.na(Theta)
      if (!any(obs)) next
      smax <- sum(cap[obs, i]); s <- sum(X[obs, i])
      sadj <- min(max(s, extreme_adj), smax - extreme_adj)
      es <- .solve_measure_item(sadj, Theta[obs], tau_list[[groups[i]]])
      Delta[i] <- es$measure; SEd[i] <- es$se
    }
  }

  ## --- full-data model matrices at final estimates ---
  full <- .model_matrices(Theta, Delta, tau_list, groups)

  structure(list(
    X = X, mask = Mfull, groups = groups, max_cat = max_cat, cap = cap,
    theta = Theta, se_theta = SEt, delta = Delta, se_delta = SEd,
    tau = tau_list, se_tau = se_tau,
    E = full$E, W = full$W, C = full$C,
    keep_p = keep_p, keep_i = keep_i,
    extreme_persons = which(!keep_p), extreme_items = which(!keep_i),
    iterations = it, converged = converged, change_history = hist,
    settings = list(maxit = maxit, conv = conv, rconv = rconv,
                    extreme_adj = extreme_adj, center = center),
    person_id = rownames(X), item_id = colnames(X)
  ), class = "raschfit")
}

#' Measure for a person with given (possibly adjusted) raw score, items anchored.
.solve_measure_person <- function(r, delta, tau_list, groups, tol = 1e-8, maxit = 200) {
  th <- 0
  for (k in seq_len(maxit)) {
    E <- W <- 0
    for (g in unique(groups)) {
      cols <- which(groups == g)
      P <- .probs_group(th, delta[cols], tau_list[[g]])
      mo <- .moments(P)
      E <- E + sum(mo$E); W <- W + sum(mo$W)
    }
    d <- (r - E) / max(W, 1e-8)
    d <- max(min(d, 1), -1)
    th <- th + d
    if (abs(d) < tol) break
  }
  list(measure = th, se = 1 / sqrt(max(W, 1e-8)))
}

#' Measure for an item with given (possibly adjusted) raw score, persons anchored.
.solve_measure_item <- function(s, theta, tau, tol = 1e-8, maxit = 200) {
  de <- 0
  for (k in seq_len(maxit)) {
    P  <- .probs_group(theta, de, tau)
    mo <- .moments(P)
    E  <- sum(mo$E); W <- sum(mo$W)
    d  <- (E - s) / max(W, 1e-8)
    d  <- max(min(d, 1), -1)
    de <- de + d
    if (abs(d) < tol) break
  }
  list(measure = de, se = 1 / sqrt(max(W, 1e-8)))
}

# ---------------------------------------------------------------------------
# 6. Fit statistics
# ---------------------------------------------------------------------------

#' Wilson-Hilferty cube-root transformation of a mean-square to a t/z statistic.
#' ZSTD = (MNSQ^(1/3) - 1) * (3/q) + q/3   (Wright & Masters 1982, p.100-101)
.zstd <- function(mnsq, q2) {
  q <- sqrt(pmax(q2, 1e-12))
  z <- (mnsq^(1/3) - 1) * (3 / q) + q / 3
  z[!is.finite(z)] <- NA_real_
  z
}

#' Item and person fit statistics.
#' @param margin 2 = items (columns), 1 = persons (rows)
.fit_stats <- function(fit, margin = 2L) {
  X <- fit$X; M <- fit$mask; E <- fit$E; W <- fit$W; C <- fit$C
  Y <- (X - E)                      # raw residual
  Z2 <- (Y^2) / W                   # squared standardized residual
  keep <- if (margin == 2L) fit$keep_i else fit$keep_p
  # exclude extreme persons from item statistics and vice versa (WINSTEPS default)
  Msub <- M
  Msub[!fit$keep_p, ] <- FALSE
  Msub[, !fit$keep_i] <- FALSE

  sums <- function(A) if (margin == 2L) colSums(A * Msub, na.rm = TRUE) else rowSums(A * Msub, na.rm = TRUE)
  n    <- if (margin == 2L) colSums(Msub) else rowSums(Msub)

  sumW  <- sums(W)
  outms <- sums(Z2) / pmax(n, 1)
  infms <- sums(Y^2) / pmax(sumW, 1e-8)

  q2_out <- sums(C / W^2 - 1) / pmax(n^2, 1)
  q2_in  <- sums(C - W^2) / pmax(sumW^2, 1e-8)
  out_z  <- .zstd(outms, q2_out)
  in_z   <- .zstd(infms, q2_in)

  outms[n == 0] <- NA; infms[n == 0] <- NA
  out_z[n == 0]  <- NA; in_z[n == 0]  <- NA
  list(n = n, outfit = outms, outfit_z = out_z, infit = infms, infit_z = in_z)
}

#' Point-measure and point-biserial correlations (observed and expected).
.pt_correlations <- function(fit, margin = 2L) {
  X <- fit$X; E <- fit$E
  Msub <- fit$mask
  Msub[!fit$keep_p, ] <- FALSE
  Msub[, !fit$keep_i] <- FALSE
  meas_other <- if (margin == 2L) fit$theta else fit$delta
  tot <- if (margin == 2L) rowSums(X * Msub, na.rm = TRUE) else colSums(X * Msub, na.rm = TRUE)
  K <- if (margin == 2L) ncol(X) else nrow(X)

  W <- fit$W
  ptm <- ptme <- ptb <- ptbe <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    sel <- if (margin == 2L) Msub[, k] else Msub[k, ]
    if (sum(sel) < 3) next
    x  <- if (margin == 2L) X[sel, k] else X[k, sel]
    e  <- if (margin == 2L) E[sel, k] else E[k, sel]
    w  <- if (margin == 2L) W[sel, k] else W[k, sel]
    mo <- meas_other[sel]
    tt <- tot[sel] - x                                # PTBIS "excluding" convention
    if (stats::sd(x) > 0 && stats::sd(mo) > 0) ptm[k] <- stats::cor(x, mo)
    if (stats::sd(x) > 0 && stats::sd(tt) > 0) ptb[k] <- stats::cor(x, tt)
    # Model expectation of the observed correlation:
    #   Cov(X, y) = Cov(E, y)  and  Var(X) = Var(E) + mean(W)   under the model.
    vx <- stats::var(e) * (length(e) - 1) / length(e) + mean(w)
    if (vx > 0 && stats::sd(mo) > 0)
      ptme[k] <- stats::cov(e, mo) / sqrt(vx * stats::var(mo) * (length(mo) - 1) / length(mo))
    if (vx > 0 && stats::sd(tt) > 0)
      ptbe[k] <- stats::cov(e, tt) / sqrt(vx * stats::var(tt) * (length(tt) - 1) / length(tt))
  }
  ptme[!is.finite(ptme)] <- NA; ptbe[!is.finite(ptbe)] <- NA
  list(pt_measure = ptm, pt_measure_exp = pmin(pmax(ptme, -1), 1),
       pt_biserial = ptb, pt_biserial_exp = pmin(pmax(ptbe, -1), 1))
}

#' Post-hoc estimated discrimination: ML fit of a in P(X=k) with
#' psi_k = k*a*(theta - delta) - sum(tau), theta/delta/tau fixed
#' (WINSTEPS Table 14 "DISCR"; 1.0 = model-conforming, >1 over-discriminating).
.discrimination <- function(fit) {
  I <- ncol(fit$X); disc <- rep(NA_real_, I)
  Msub <- fit$mask; Msub[!fit$keep_p, ] <- FALSE
  for (i in seq_len(I)) {
    if (!fit$keep_i[i]) next
    sel <- Msub[, i]
    if (sum(sel) < 5) next
    th <- fit$theta[sel]; x <- fit$X[sel, i]
    tau <- fit$tau[[fit$groups[i]]]; de <- fit$delta[i]
    nll <- function(a) {
      if (!is.finite(a) || a <= 0.05 || a > 10) return(1e10)
      P <- .probs_group(a * (th - de), 0, tau)
      p <- vapply(seq_along(th), function(j) P[[x[j] + 1]][j, 1], numeric(1))
      -sum(log(pmax(p, 1e-12)))
    }
    o <- try(stats::optimize(nll, c(0.06, 6))$minimum, silent = TRUE)
    if (!inherits(o, "try-error")) disc[i] <- o
  }
  disc
}

# ---------------------------------------------------------------------------
# 7. Measure tables (WINSTEPS Tables 10 / 13-17)
# ---------------------------------------------------------------------------

#' Item measure table (WINSTEPS Table 14 layout).
item_table <- function(fit, discrimination = TRUE) {
  X <- fit$X; M <- fit$mask
  Msub <- M; Msub[!fit$keep_p, ] <- FALSE
  cnt  <- colSums(Msub)
  scr  <- colSums(X * Msub, na.rm = TRUE)
  fs   <- .fit_stats(fit, 2L)
  pc   <- .pt_correlations(fit, 2L)
  real_se <- fit$se_delta * sqrt(pmax(fs$infit, 1, na.rm = FALSE))
  real_se[is.na(fs$infit)] <- fit$se_delta[is.na(fs$infit)]
  out <- data.frame(
    Entry      = seq_along(fit$item_id),
    Item       = fit$item_id,
    Total_Score = scr,
    Total_Count = cnt,
    Measure    = fit$delta,
    Model_SE   = fit$se_delta,
    Real_SE    = real_se,
    Infit_MNSQ = fs$infit,  Infit_ZSTD  = fs$infit_z,
    Outfit_MNSQ = fs$outfit, Outfit_ZSTD = fs$outfit_z,
    PtMeasure_Corr = pc$pt_measure,
    PtMeasure_Exp  = pc$pt_measure_exp,
    PtBiserial     = pc$pt_biserial,
    PtBiserial_Exp = pc$pt_biserial_exp,
    Group      = fit$groups,
    Extreme    = ifelse(fit$keep_i, "", "EXTREME"),
    stringsAsFactors = FALSE
  )
  if (discrimination) out$Discrim <- .discrimination(fit)
  out
}

#' Person measure table (WINSTEPS Table 17 layout).
person_table <- function(fit) {
  X <- fit$X
  Msub <- fit$mask; Msub[, !fit$keep_i] <- FALSE
  cnt <- rowSums(Msub)
  scr <- rowSums(X * Msub, na.rm = TRUE)
  fs  <- .fit_stats(fit, 1L)
  pc  <- .pt_correlations(fit, 1L)
  real_se <- fit$se_theta * sqrt(pmax(fs$infit, 1))
  real_se[is.na(fs$infit)] <- fit$se_theta[is.na(fs$infit)]
  data.frame(
    Entry = seq_along(fit$person_id),
    Person = fit$person_id,
    Total_Score = scr, Total_Count = cnt,
    Measure = fit$theta, Model_SE = fit$se_theta, Real_SE = real_se,
    Infit_MNSQ = fs$infit, Infit_ZSTD = fs$infit_z,
    Outfit_MNSQ = fs$outfit, Outfit_ZSTD = fs$outfit_z,
    PtMeasure_Corr = pc$pt_measure, PtMeasure_Exp = pc$pt_measure_exp,
    Extreme = ifelse(fit$keep_p, "", "EXTREME"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 8. Summary statistics, separation and reliability (Table 3.1)
# ---------------------------------------------------------------------------

.summary_block <- function(meas, se_model, se_real, infit, outfit, score, count, label) {
  n <- length(meas)
  sd_pop <- function(v) stats::sd(v) * sqrt((length(v) - 1) / length(v))
  if (n < 2) return(NULL)
  sdm <- sd_pop(meas)
  rmse_m <- sqrt(mean(se_model^2)); rmse_r <- sqrt(mean(se_real^2))
  tsd_m <- sqrt(max(sdm^2 - rmse_m^2, 0)); tsd_r <- sqrt(max(sdm^2 - rmse_r^2, 0))
  sep_m <- if (rmse_m > 0) tsd_m / rmse_m else NA
  sep_r <- if (rmse_r > 0) tsd_r / rmse_r else NA
  rel_m <- if (is.finite(sep_m)) sep_m^2 / (1 + sep_m^2) else NA
  rel_r <- if (is.finite(sep_r)) sep_r^2 / (1 + sep_r^2) else NA
  data.frame(
    Statistic = label, N = n,
    Mean_Score = mean(score), Mean_Count = mean(count),
    Mean_Measure = mean(meas), SD_Measure = sdm,
    Mean_Model_SE = mean(se_model), Mean_Real_SE = mean(se_real),
    Mean_Infit_MNSQ = mean(infit, na.rm = TRUE),
    Mean_Outfit_MNSQ = mean(outfit, na.rm = TRUE),
    Model_RMSE = rmse_m, Real_RMSE = rmse_r,
    True_SD_Model = tsd_m, True_SD_Real = tsd_r,
    Separation_Model = sep_m, Separation_Real = sep_r,
    Reliability_Model = rel_m, Reliability_Real = rel_r,
    Strata_Model = (4 * sep_m + 1) / 3, Strata_Real = (4 * sep_r + 1) / 3,
    stringsAsFactors = FALSE
  )
}

#' WINSTEPS Table 3.1 summary of persons and items.
summary_table <- function(fit) {
  it <- item_table(fit, discrimination = FALSE)
  pt <- person_table(fit)
  rows <- list(
    .summary_block(pt$Measure[fit$keep_p], pt$Model_SE[fit$keep_p], pt$Real_SE[fit$keep_p],
                   pt$Infit_MNSQ[fit$keep_p], pt$Outfit_MNSQ[fit$keep_p],
                   pt$Total_Score[fit$keep_p], pt$Total_Count[fit$keep_p],
                   "PERSON (non-extreme)"),
    .summary_block(pt$Measure, pt$Model_SE, pt$Real_SE, pt$Infit_MNSQ, pt$Outfit_MNSQ,
                   pt$Total_Score, pt$Total_Count, "PERSON (all)"),
    .summary_block(it$Measure[fit$keep_i], it$Model_SE[fit$keep_i], it$Real_SE[fit$keep_i],
                   it$Infit_MNSQ[fit$keep_i], it$Outfit_MNSQ[fit$keep_i],
                   it$Total_Score[fit$keep_i], it$Total_Count[fit$keep_i],
                   "ITEM (non-extreme)"),
    .summary_block(it$Measure, it$Model_SE, it$Real_SE, it$Infit_MNSQ, it$Outfit_MNSQ,
                   it$Total_Score, it$Total_Count, "ITEM (all)")
  )
  do.call(rbind, Filter(Negate(is.null), rows))
}

#' Cronbach alpha (KR-20) on complete raw scores, for comparison (Table 3.1 footer).
cronbach_alpha <- function(fit) {
  X <- fit$X[, fit$keep_i, drop = FALSE]
  cc <- stats::complete.cases(X)
  if (sum(cc) < 3) return(NA_real_)
  Xc <- X[cc, , drop = FALSE]; k <- ncol(Xc)
  v <- apply(Xc, 2, stats::var); tv <- stats::var(rowSums(Xc))
  if (tv <= 0) return(NA_real_)
  (k / (k - 1)) * (1 - sum(v) / tv)
}

# ---------------------------------------------------------------------------
# 9. Category / rating-scale structure (Tables 3.2, 3.3)
# ---------------------------------------------------------------------------

#' Measure (relative to item difficulty) at which the expected score = target.
.score_to_measure <- function(target, tau, lo = -25, hi = 25) {
  f <- function(x) {
    P <- .probs_group(x, 0, tau)
    sum(vapply(seq_along(P), function(k) (k - 1) * P[[k]][1, 1], numeric(1))) - target
  }
  if (f(lo) > 0 || f(hi) < 0) return(NA_real_)
  stats::uniroot(f, c(lo, hi), tol = 1e-8)$root
}

#' Rasch-Thurstone threshold: measure at which P(X >= k) = 0.5.
.thurstone <- function(k, tau, lo = -25, hi = 25) {
  f <- function(x) {
    P <- .probs_group(x, 0, tau)
    sum(vapply((k + 1):(length(tau) + 1), function(j) P[[j]][1, 1], numeric(1))) - 0.5
  }
  if (f(lo) > 0 || f(hi) < 0) return(NA_real_)
  stats::uniroot(f, c(lo, hi), tol = 1e-8)$root
}

#' WINSTEPS Table 3.2: category structure for one item group.
#' @param cat_extreme score adjustment used for the extreme "category measure"
#'   (WINSTEPS uses 0.25 for the parenthesised extreme category measures)
category_table <- function(fit, group = NULL, cat_extreme = 0.25) {
  gl <- if (is.null(group)) names(fit$tau) else group
  res <- list()
  for (g in gl) {
    cols <- which(fit$groups == g)
    m <- fit$max_cat[[g]]
    tau <- fit$tau[[g]]
    # observations used for threshold estimation: non-extreme persons & items
    sel <- fit$mask
    sel[!fit$keep_p, ] <- FALSE
    sel[, !fit$keep_i] <- FALSE
    selg <- sel[, cols, drop = FALSE]
    xg   <- fit$X[, cols, drop = FALSE]
    bd   <- outer(fit$theta, fit$delta[cols], "-")   # Bn - Di
    Eg   <- fit$E[, cols, drop = FALSE]
    Wg   <- fit$W[, cols, drop = FALSE]
    Cg   <- fit$C[, cols, drop = FALSE]
    P    <- .probs_group(fit$theta, fit$delta[cols], tau)

    cnt <- ez <- oa <- ea <- outms <- infms <- numeric(m + 1)
    for (k in 0:(m)) {
      inc <- selg & !is.na(xg) & xg == k
      cnt[k + 1] <- sum(inc)
      oa[k + 1]  <- if (cnt[k + 1] > 0) mean(bd[inc]) else NA_real_
      w  <- P[[k + 1]] * selg
      ea[k + 1]  <- if (sum(w) > 0) sum(bd * w) / sum(w) else NA_real_
      if (cnt[k + 1] > 0) {
        y <- (xg - Eg)[inc]; wv <- Wg[inc]
        outms[k + 1] <- mean(y^2 / wv)
        infms[k + 1] <- sum(y^2) / sum(wv)
      } else { outms[k + 1] <- NA; infms[k + 1] <- NA }
    }
    andrich <- c(NA_real_, tau)
    se_and  <- c(NA_real_, fit$se_tau[[g]])
    catmeas <- vapply(0:m, function(k) {
      tgt <- if (k == 0) cat_extreme else if (k == m) m - cat_extreme else k
      .score_to_measure(tgt, tau)
    }, numeric(1))
    thur <- c(NA_real_, vapply(1:m, function(k) .thurstone(k, tau), numeric(1)))
    # zones: measure interval where the expected score rounds to k
    zlo <- vapply(0:m, function(k) if (k == 0) -Inf else .score_to_measure(k - 0.5, tau), numeric(1))
    zhi <- vapply(0:m, function(k) if (k == m)  Inf else .score_to_measure(k + 0.5, tau), numeric(1))
    # coherence
    mc <- cm <- rep(NA_real_, m + 1)
    for (k in 0:m) {
      inzone <- selg & !is.na(xg) & bd >= zlo[k + 1] & bd < zhi[k + 1]
      incat  <- selg & !is.na(xg) & xg == k
      mc[k + 1] <- if (sum(inzone) > 0) 100 * sum(inzone & incat) / sum(inzone) else NA
      cm[k + 1] <- if (sum(incat)  > 0) 100 * sum(inzone & incat) / sum(incat)  else NA
    }
    tot <- sum(cnt)
    res[[g]] <- data.frame(
      Group = g, Category = 0:m, Score = 0:m, Count = cnt,
      Percent = 100 * cnt / max(tot, 1),
      Obsvd_Avrge = oa, Sample_Expect = ea,
      Infit_MNSQ = infms, Outfit_MNSQ = outms,
      Andrich_Threshold = andrich, Threshold_SE = se_and,
      Category_Measure = catmeas,
      Thurstone_50pct = thur,
      Zone_Lo = zlo, Zone_Hi = zhi,
      Coherence_M_to_C = mc, Coherence_C_to_M = cm,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, res)
}

#' Linacre's (1999, 2002) rating-scale quality guidelines, applied to a
#' category_table() for one group. Returns a data.frame of checks.
category_diagnostics <- function(ct) {
  out <- list()
  for (g in unique(ct$Group)) {
    d <- ct[ct$Group == g, ]
    m <- nrow(d) - 1
    tau <- d$Andrich_Threshold[-1]
    chk <- function(name, ok, detail) data.frame(Group = g, Guideline = name,
                                                 Status = ifelse(ok, "OK", "REVIEW"),
                                                 Detail = detail, stringsAsFactors = FALSE)
    out[[length(out) + 1]] <- chk(
      "At least 10 observations per category",
      all(d$Count >= 10),
      paste0("min count = ", min(d$Count)))
    out[[length(out) + 1]] <- chk(
      "Regular / uniform category distribution",
      all(d$Percent > 0),
      paste0("percentages: ", paste(sprintf("%.1f", d$Percent), collapse = ", ")))
    oa <- d$Obsvd_Avrge
    out[[length(out) + 1]] <- chk(
      "Observed average measures advance with category",
      all(diff(oa[!is.na(oa)]) > 0),
      paste0("observed averages: ", paste(sprintf("%.2f", oa), collapse = ", ")))
    out[[length(out) + 1]] <- chk(
      "Category outfit mean-squares < 2.0",
      all(d$Outfit_MNSQ < 2, na.rm = TRUE),
      paste0("max outfit = ", sprintf("%.2f", max(d$Outfit_MNSQ, na.rm = TRUE))))
    if (m >= 2) {
      out[[length(out) + 1]] <- chk(
        "Andrich thresholds advance (ordered)",
        all(diff(tau) > 0),
        paste0("thresholds: ", paste(sprintf("%.2f", tau), collapse = ", ")))
      dif <- diff(tau)
      out[[length(out) + 1]] <- chk(
        "Threshold advances >= 1.4 logits (distinct categories)",
        all(dif >= 1.4),
        paste0("advances: ", paste(sprintf("%.2f", dif), collapse = ", ")))
      out[[length(out) + 1]] <- chk(
        "Threshold advances <= 5.0 logits (no gaps in the variable)",
        all(dif <= 5),
        paste0("advances: ", paste(sprintf("%.2f", dif), collapse = ", ")))
    }
    out[[length(out) + 1]] <- chk(
      "Coherence M->C >= 40% (Linacre 2002)",
      all(d$Coherence_M_to_C >= 40, na.rm = TRUE),
      paste0("M->C: ", paste(sprintf("%.0f%%", d$Coherence_M_to_C), collapse = ", ")))
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# 10. Curves: category probabilities, ICC, information
# ---------------------------------------------------------------------------

#' Category probability / cumulative probability / expected score curves for a group.
curve_data <- function(fit, group, delta = 0, from = -6, to = 6, n = 401) {
  tau <- fit$tau[[group]]; m <- fit$max_cat[[group]]
  x <- seq(from, to, length.out = n)
  P <- .probs_group(x, delta, tau)
  prob <- do.call(cbind, lapply(P, function(p) p[, 1]))
  colnames(prob) <- paste0("P", 0:m)
  cum <- do.call(cbind, lapply(1:m, function(k) rowSums(prob[, (k + 1):(m + 1), drop = FALSE])))
  colnames(cum) <- paste0("P_ge_", 1:m)
  E <- as.vector(prob %*% (0:m))
  W <- as.vector(prob %*% ((0:m)^2)) - E^2
  list(measure = x, prob = prob, cum = cum, expected = E, information = W)
}

#' Empirical ICC: observed mean score per measure-bin for one item.
empirical_icc <- function(fit, item, nbins = 10) {
  i <- if (is.character(item)) match(item, fit$item_id) else item
  sel <- fit$mask[, i] & fit$keep_p
  if (sum(sel) < nbins) return(NULL)
  th <- fit$theta[sel]; x <- fit$X[sel, i]
  br <- unique(stats::quantile(th, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  if (length(br) < 3) return(NULL)
  b <- cut(th, br, include.lowest = TRUE)
  data.frame(measure = tapply(th, b, mean), observed = tapply(x, b, mean),
             n = as.vector(table(b)))
}

#' Test characteristic curve and test information (Table 20 basis).
test_curves <- function(fit, from = -6, to = 6, n = 241, items = NULL) {
  idx <- if (is.null(items)) which(!is.na(fit$delta)) else items
  x <- seq(from, to, length.out = n)
  E <- W <- numeric(n)
  for (g in unique(fit$groups[idx])) {
    cols <- idx[fit$groups[idx] == g]
    P <- .probs_group(x, fit$delta[cols], fit$tau[[g]])
    mo <- .moments(P)
    E <- E + rowSums(mo$E); W <- W + rowSums(mo$W)
  }
  data.frame(measure = x, expected_score = E, information = W, SE = 1 / sqrt(pmax(W, 1e-9)))
}

#' Score-to-measure table (WINSTEPS Table 20) over the full raw-score range.
score_table <- function(fit, items = NULL, extreme_adj = NULL) {
  idx <- if (is.null(items)) which(!is.na(fit$delta)) else items
  extreme_adj <- extreme_adj %||% fit$settings$extreme_adj
  maxs <- sum(fit$max_cat[fit$groups[idx]])
  sc <- 0:maxs
  meas <- se <- numeric(length(sc))
  for (j in seq_along(sc)) {
    r <- sc[j]
    radj <- min(max(r, extreme_adj), maxs - extreme_adj)
    es <- .solve_measure_person(radj, fit$delta[idx], fit$tau, fit$groups[idx])
    meas[j] <- es$measure; se[j] <- es$se
  }
  data.frame(Score = sc, Measure = meas, SE = se,
             Extreme = ifelse(sc == 0 | sc == maxs, "EXTREME", ""))
}

# ---------------------------------------------------------------------------
# 11. Wright map data
# ---------------------------------------------------------------------------

wright_data <- function(fit, item_labels = NULL, include_extreme = FALSE) {
  pk <- if (include_extreme) rep(TRUE, length(fit$theta)) else fit$keep_p
  ik <- if (include_extreme) rep(TRUE, length(fit$delta)) else fit$keep_i
  list(
    persons = data.frame(measure = fit$theta[pk & !is.na(fit$theta)]),
    items = data.frame(measure = fit$delta[ik & !is.na(fit$delta)],
                       label = (item_labels %||% fit$item_id)[ik & !is.na(fit$delta)],
                       stringsAsFactors = FALSE)
  )
}

#' Item thresholds (Rasch-Thurstone) for a "threshold" Wright map.
threshold_data <- function(fit) {
  out <- list()
  for (i in seq_along(fit$delta)) {
    g <- fit$groups[i]; tau <- fit$tau[[g]]; m <- fit$max_cat[[g]]
    if (is.na(fit$delta[i])) next
    th <- vapply(1:m, function(k) .thurstone(k, tau), numeric(1))
    out[[length(out) + 1]] <- data.frame(
      item = fit$item_id[i], threshold = 1:m,
      measure = fit$delta[i] + th,
      label = paste0(fit$item_id[i], ".", 1:m), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# 12. DIF (WINSTEPS Table 30)
# ---------------------------------------------------------------------------

#' Rasch-Welch iterative-logit DIF plus Mantel-Haenszel / Mantel tests.
#'
#' Person measures and Andrich thresholds are held at their baseline values;
#' a local item difficulty is re-estimated within each person class
#' (exactly the WINSTEPS "iterative-logit" procedure). The DIF contrast is the
#' difference of the local difficulties; t = contrast / sqrt(SE_A^2 + SE_B^2)
#' with Welch-Satterthwaite degrees of freedom.
#'
#' @param class factor of length nrow(X) classifying persons
#' @param mhslice logit width of the measure strata for Mantel-Haenszel
#'   (WINSTEPS MHSLICE=, default 0.1)
dif_analysis <- function(fit, class, mhslice = 0.1, min_n = 5) {
  class <- as.factor(class)
  stopifnot(length(class) == nrow(fit$X))
  lv <- levels(droplevels(class[fit$keep_p]))
  if (length(lv) < 2) stop("The DIF classification needs at least two groups with non-extreme persons.")
  I <- ncol(fit$X)

  ## ---- local (class-specific) item difficulties ----
  loc <- matrix(NA_real_, I, length(lv), dimnames = list(fit$item_id, lv))
  locse <- obsexp <- nobs <- loc
  for (l in lv) {
    sel <- fit$keep_p & !is.na(class) & class == l
    for (i in seq_len(I)) {
      s <- sel & fit$mask[, i]
      nobs[i, l] <- sum(s)
      if (sum(s) < min_n || !fit$keep_i[i]) next
      x <- fit$X[s, i]; th <- fit$theta[s]
      tau <- fit$tau[[fit$groups[i]]]
      cap <- fit$max_cat[[fit$groups[i]]]
      raw <- sum(x); rmax <- cap * sum(s)
      radj <- min(max(raw, fit$settings$extreme_adj), rmax - fit$settings$extreme_adj)
      es <- .solve_measure_item(radj, th, tau)
      loc[i, l] <- es$measure; locse[i, l] <- es$se
      obsexp[i, l] <- mean(x - fit$E[s, i])
    }
  }

  ## ---- pairwise contrasts ----
  pairs <- utils::combn(lv, 2, simplify = FALSE)
  rows <- list()
  for (pp in pairs) {
    A <- pp[1]; B <- pp[2]
    dA <- loc[, A]; dB <- loc[, B]
    sA <- locse[, A]; sB <- locse[, B]
    nA <- nobs[, A]; nB <- nobs[, B]
    contrast <- dA - dB
    jse <- sqrt(sA^2 + sB^2)
    tval <- contrast / jse
    df <- (sA^2 + sB^2)^2 / (sA^4 / pmax(nA - 1, 1) + sB^4 / pmax(nB - 1, 1))
    pval <- 2 * stats::pt(-abs(tval), df = pmax(df, 1))
    mh <- .mantel_dif(fit, class, A, B, mhslice, min_n)
    rows[[length(rows) + 1]] <- data.frame(
      Item = fit$item_id, Entry = seq_len(I),
      Class_A = A, ObsExp_A = obsexp[, A], DIF_Measure_A = dA, DIF_SE_A = sA, N_A = nA,
      Class_B = B, ObsExp_B = obsexp[, B], DIF_Measure_B = dB, DIF_SE_B = sB, N_B = nB,
      DIF_Contrast = contrast, Joint_SE = jse, t = tval, df = df, p = pval,
      MH_ChiSq = mh$chisq, MH_p = mh$p, MH_Size_CUMLOR = mh$lor,
      ETS_Class = .ets_class(contrast, pval),
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  attr(res, "local_measures") <- loc
  attr(res, "local_se") <- locse
  res
}

#' ETS DIF classification applied to the logit contrast.
#' |contrast| < 0.43 -> A (negligible); 0.43-0.64 -> B; >= 0.64 -> C,
#' with the significance requirement (Zwick 2012; Linacre, Winsteps Help).
.ets_class <- function(contrast, p) {
  cl <- rep("A", length(contrast))
  cl[abs(contrast) >= 0.43 & p < 0.05] <- "B"
  cl[abs(contrast) >= 0.64 & p < 0.05] <- "C"
  cl[is.na(contrast)] <- NA
  cl
}

#' Mantel-Haenszel (dichotomies) / Mantel (polytomies) DIF, stratified by measure.
.mantel_dif <- function(fit, class, A, B, mhslice = 0.1, min_n = 5) {
  I <- ncol(fit$X)
  chisq <- p <- lor <- rep(NA_real_, I)
  keep <- fit$keep_p & !is.na(class) & class %in% c(A, B)
  if (!any(keep)) return(list(chisq = chisq, p = p, lor = lor))
  th <- fit$theta
  stratum <- floor(th / mhslice)
  for (i in seq_len(I)) {
    s <- keep & fit$mask[, i]
    if (sum(s) < 2 * min_n) next
    x <- fit$X[s, i]; gcl <- as.character(class[s]); st <- stratum[s]
    m <- fit$max_cat[[fit$groups[i]]]
    num <- den <- 0; numL <- denL <- 0
    for (u in unique(st)) {
      k <- st == u
      if (sum(k) < 2) next
      xf <- x[k]; gf <- gcl[k]
      nA <- sum(gf == A); nB <- sum(gf == B); nt <- nA + nB
      if (nA == 0 || nB == 0) next
      # Mantel (1963) statistic on the focal-group score total
      Y  <- sum(xf[gf == A])
      mu <- nA * mean(xf)
      vr <- (nA * nB / (nt^2 * (nt - 1))) * nt * sum((xf - mean(xf))^2)
      num <- num + (Y - mu); den <- den + vr
      if (m == 1) {                       # classic MH odds-ratio components
        a <- sum(xf[gf == A] == 1); b <- sum(xf[gf == A] == 0)
        cc <- sum(xf[gf == B] == 1); dd <- sum(xf[gf == B] == 0)
        numL <- numL + a * dd / nt; denL <- denL + b * cc / nt
      }
    }
    if (den > 0) {
      chisq[i] <- (abs(num) - 0.5)^2 / den
      p[i] <- stats::pchisq(chisq[i], 1, lower.tail = FALSE)
    }
    if (m == 1 && denL > 0 && numL > 0) lor[i] <- log(numL / denL)
  }
  list(chisq = chisq, p = p, lor = lor)
}

# ---------------------------------------------------------------------------
# 13. PCA of residuals + variance decomposition (WINSTEPS Table 23)
# ---------------------------------------------------------------------------

#' Standardized (or raw / logit) residuals.
residual_matrix <- function(fit, type = c("standardized", "raw", "score")) {
  type <- match.arg(type)
  R <- fit$X - fit$E
  if (type == "standardized") R <- R / sqrt(fit$W)
  R[!fit$mask] <- NA
  R
}

#' Central person / item measure such that the model total equals the observed total.
.central_measure <- function(fit, what = c("person", "item")) {
  what <- match.arg(what)
  M <- fit$mask; M[!fit$keep_p, ] <- FALSE; M[, !fit$keep_i] <- FALSE
  obs_total <- sum(fit$X * M, na.rm = TRUE)
  f <- function(v) {
    E <- matrix(0, nrow(fit$X), ncol(fit$X))
    for (g in unique(fit$groups)) {
      cols <- which(fit$groups == g)
      th <- if (what == "person") rep(v, nrow(fit$X)) else fit$theta
      de <- if (what == "person") fit$delta[cols] else rep(v, length(cols))
      th[is.na(th)] <- 0; de[is.na(de)] <- 0
      P <- .probs_group(th, de, fit$tau[[g]])
      E[, cols] <- .moments(P)$E
    }
    sum(E * M, na.rm = TRUE) - obs_total
  }
  stats::uniroot(f, c(-20, 20), extendInt = "yes", tol = 1e-8)$root
}

#' Full WINSTEPS Table 23 variance decomposition + PCA of residuals.
#'
#' Follows the 24-step recipe in the Winsteps Help for Tables 23.0 / 24.0:
#' raw residuals for the explained/unexplained split, standardized residuals
#' for the PCA, and rescaling of the whole table so that the unexplained
#' variance equals the number of non-extreme items ("eigenvalue units").
pca_residuals <- function(fit, ncontrast = 5, residual_type = "standardized") {
  M <- fit$mask; M[!fit$keep_p, ] <- FALSE; M[, !fit$keep_i] <- FALSE
  ii <- which(fit$keep_i); pp <- which(fit$keep_p)
  if (length(ii) < 3) stop("Too few non-extreme items for a PCA of residuals.")

  Bc <- .central_measure(fit, "person")
  Dc <- .central_measure(fit, "item")
  # central prediction for every observation
  Ecen <- matrix(NA_real_, nrow(fit$X), ncol(fit$X))
  Epers <- Eitem <- Ecen
  for (g in unique(fit$groups)) {
    cols <- which(fit$groups == g); tau <- fit$tau[[g]]
    th <- fit$theta; th[is.na(th)] <- Bc
    de <- fit$delta[cols]; de[is.na(de)] <- Dc
    Ecen[, cols]  <- .moments(.probs_group(rep(Bc, nrow(fit$X)), rep(Dc, length(cols)), tau))$E
    Epers[, cols] <- .moments(.probs_group(th, rep(Dc, length(cols)), tau))$E
    Eitem[, cols] <- .moments(.probs_group(rep(Bc, nrow(fit$X)), de, tau))$E
  }

  ss <- function(A) sum((A^2) * M, na.rm = TRUE)
  total_obs  <- ss(fit$X - Ecen)                       # step 3
  unexp_obs  <- ss(fit$X - fit$E)                      # step 4
  expl_obs   <- total_obs - unexp_obs                  # step 5
  sp <- ss(Epers - Ecen); si <- ss(Eitem - Ecen)       # steps 6-8
  pers_obs <- expl_obs * sp / max(sp + si, 1e-12)
  item_obs <- expl_obs * si / max(sp + si, 1e-12)

  expl_mod  <- ss(fit$E - Ecen)                        # step 17
  unexp_mod <- sum(fit$W * M, na.rm = TRUE)            # step 21
  total_mod <- expl_mod + unexp_mod                    # step 22
  spm <- sp; sim <- si
  pers_mod <- expl_mod * spm / max(spm + sim, 1e-12)
  item_mod <- expl_mod * sim / max(spm + sim, 1e-12)

  ## ---- PCA of the (standardized) residual correlation matrix ----
  R <- residual_matrix(fit, residual_type)[pp, ii, drop = FALSE]
  Rc <- suppressWarnings(stats::cor(R, use = "pairwise.complete.obs"))
  Rc[!is.finite(Rc)] <- 0; diag(Rc) <- 1
  eg <- eigen(Rc, symmetric = TRUE)
  ev <- eg$values
  nI <- length(ii)
  # step 11: rescale so unexplained variance == sum of eigenvalues == nI
  scale <- nI / unexp_obs
  comp <- data.frame(
    Component = c("Total raw variance in observations",
                  "Raw variance explained by measures",
                  "  Raw variance explained by persons",
                  "  Raw variance explained by items",
                  "Raw unexplained variance (total)",
                  paste0("Unexplned variance in ", c("1st", "2nd", "3rd", "4th", "5th")[seq_len(ncontrast)], " contrast")),
    Eigenvalue = c(total_obs * scale, expl_obs * scale, pers_obs * scale,
                   item_obs * scale, nI, ev[seq_len(ncontrast)]),
    stringsAsFactors = FALSE)
  comp$Pct_of_total <- 100 * comp$Eigenvalue / (total_obs * scale)
  comp$Pct_of_unexplained <- 100 * comp$Eigenvalue / nI
  comp$Pct_of_unexplained[1:4] <- NA
  scale_m <- nI / unexp_mod
  comp$Expected_Pct <- c(100,
                         100 * expl_mod / total_mod,
                         100 * pers_mod / total_mod,
                         100 * item_mod / total_mod,
                         100 * unexp_mod / total_mod,
                         rep(NA_real_, ncontrast))

  ## ---- loadings ----
  ld <- eg$vectors[, seq_len(ncontrast), drop = FALSE] %*%
    diag(sqrt(pmax(ev[seq_len(ncontrast)], 0)), ncontrast)
  colnames(ld) <- paste0("Contrast", seq_len(ncontrast))
  itab <- item_table(fit, discrimination = FALSE)[ii, ]
  loadings <- data.frame(
    Entry = itab$Entry, Item = itab$Item, Measure = itab$Measure,
    Infit_MNSQ = itab$Infit_MNSQ, Outfit_MNSQ = itab$Outfit_MNSQ,
    as.data.frame(ld), stringsAsFactors = FALSE)

  ## ---- item clusters on contrast 1 (Winsteps Table 23.1 uses 3 clusters) ----
  l1 <- ld[, 1]
  cluster <- ifelse(l1 >= stats::quantile(l1, 2/3), 1L,
                    ifelse(l1 <= stats::quantile(l1, 1/3), 3L, 2L))
  loadings$Cluster <- cluster

  ## ---- disattenuated correlations between person measures on clusters ----
  clcor <- .cluster_correlations(fit, ii, cluster)

  ## ---- essential unidimensionality (Liu et al. 2023) ----
  keepc <- ev[seq_len(ncontrast)] >= 1.3
  ess <- (expl_obs * scale) / (expl_obs * scale + sum(ev[seq_len(ncontrast)][keepc]))

  list(variance = comp, eigenvalues = ev, loadings = loadings,
       cluster_correlations = clcor,
       essential_unidimensionality = 100 * ess,
       central = c(person = Bc, item = Dc),
       n_items = nI, n_persons = length(pp),
       residual_type = residual_type)
}

#' Person measures on an arbitrary subset of items, items anchored.
#' Uses a (response-pattern x raw-score) lookup so it stays fast on large data.
subset_measures <- function(fit, cols, extreme_adj = NULL) {
  extreme_adj <- extreme_adj %||% fit$settings$extreme_adj
  N <- nrow(fit$X)
  meas <- se <- rep(NA_real_, N)
  Msub <- fit$mask[, cols, drop = FALSE]
  Xsub <- fit$X[, cols, drop = FALSE]
  caps <- fit$max_cat[fit$groups[cols]]
  key  <- apply(Msub, 1, function(v) paste0(which(v), collapse = ","))
  raw  <- rowSums(Xsub * Msub, na.rm = TRUE)
  for (k in unique(key)) {
    if (!nzchar(k)) next
    idx <- as.integer(strsplit(k, ",", fixed = TRUE)[[1]])
    obs <- cols[idx]
    rmax <- sum(caps[idx])
    rows <- which(key == k)
    for (r in unique(raw[rows])) {
      radj <- min(max(r, extreme_adj), rmax - extreme_adj)
      es <- .solve_measure_person(radj, fit$delta[obs], fit$tau, fit$groups[obs])
      sel <- rows[raw[rows] == r]
      meas[sel] <- es$measure; se[sel] <- es$se
    }
  }
  list(measure = meas, se = se)
}

#' Correlate (and disattenuate) person measures computed separately on each
#' item cluster from the first residual contrast (WINSTEPS Table 23.1 / 23.6).
.cluster_correlations <- function(fit, ii, cluster) {
  cl <- sort(unique(cluster))
  meas <- se <- matrix(NA_real_, nrow(fit$X), length(cl))
  rel  <- rep(NA_real_, length(cl))
  for (j in seq_along(cl)) {
    cols <- ii[cluster == cl[j]]
    if (length(cols) < 2) next
    sm <- subset_measures(fit, cols)
    meas[, j] <- ifelse(fit$keep_p, sm$measure, NA)
    se[, j]   <- ifelse(fit$keep_p, sm$se, NA)
    v <- meas[is.finite(meas[, j]), j]; s <- se[is.finite(meas[, j]), j]
    if (length(v) > 2) {
      sdv <- stats::sd(v) * sqrt((length(v) - 1) / length(v))
      mse <- mean(s^2)
      rel[j] <- max(0, min(1, (sdv^2 - mse) / sdv^2))
    }
  }
  out <- list()
  for (a in seq_along(cl)) for (b in seq_along(cl)) if (a < b) {
    ok <- is.finite(meas[, a]) & is.finite(meas[, b])
    if (sum(ok) < 5) next
    r <- stats::cor(meas[ok, a], meas[ok, b])
    den <- sqrt(rel[a] * rel[b])
    out[[length(out) + 1]] <- data.frame(
      Cluster_A = cl[a], Cluster_B = cl[b], n = sum(ok),
      Reliability_A = rel[a], Reliability_B = rel[b],
      Correlation = r,
      Disattenuated = if (is.finite(den) && den > 0) min(r / den, 1) else NA_real_,
      stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  attr(res, "measures") <- meas
  res
}

# ---------------------------------------------------------------------------
# 14. Misfit / unexpected responses (Table 6.6 style)
# ---------------------------------------------------------------------------

unexpected_responses <- function(fit, cut = 2) {
  Z <- residual_matrix(fit, "standardized")
  idx <- which(abs(Z) >= cut & fit$mask, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame())
  d <- data.frame(
    Person = fit$person_id[idx[, 1]],
    Item = fit$item_id[idx[, 2]],
    Observed = fit$X[idx],
    Expected = fit$E[idx],
    Residual = (fit$X - fit$E)[idx],
    Std_Residual = Z[idx],
    Person_Measure = fit$theta[idx[, 1]],
    Item_Measure = fit$delta[idx[, 2]],
    stringsAsFactors = FALSE)
  d[order(-abs(d$Std_Residual)), ]
}

# ---------------------------------------------------------------------------
# 15. Simulation (SIFILE=): for parallel analysis of contrast eigenvalues
# ---------------------------------------------------------------------------

simulate_rasch <- function(fit, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- fit$X
  for (g in unique(fit$groups)) {
    cols <- which(fit$groups == g)
    th <- fit$theta; th[is.na(th)] <- 0
    de <- fit$delta[cols]; de[is.na(de)] <- 0
    P <- .probs_group(th, de, fit$tau[[g]])
    m <- length(P) - 1
    cum <- Reduce(`+`, P, accumulate = TRUE)
    U <- matrix(stats::runif(length(th) * length(cols)), length(th), length(cols))
    S <- matrix(0L, length(th), length(cols))
    for (k in 1:m) S <- S + (U > cum[[k]])
    X[, cols] <- S
  }
  X[!fit$mask] <- NA
  X
}

#' Parallel analysis: expected sizes of the residual contrasts under the model.
parallel_contrasts <- function(prep, fit, nsim = 20, ncontrast = 5, seed = 1) {
  set.seed(seed)
  out <- matrix(NA_real_, nsim, ncontrast)
  for (s in seq_len(nsim)) {
    Xs <- simulate_rasch(fit)
    ps <- prep; ps$X <- Xs
    fs <- try(rasch_jmle(ps, maxit = 150, conv = 1e-3, rconv = 1e-2), silent = TRUE)
    if (inherits(fs, "try-error")) next
    pc <- try(pca_residuals(fs, ncontrast), silent = TRUE)
    if (inherits(pc, "try-error")) next
    out[s, ] <- pc$eigenvalues[seq_len(ncontrast)]
  }
  data.frame(Contrast = seq_len(ncontrast),
             Mean = colMeans(out, na.rm = TRUE),
             SD = apply(out, 2, stats::sd, na.rm = TRUE),
             P95 = apply(out, 2, stats::quantile, probs = .95, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# 16. Demo data (also used by the app's 'Load built-in demo data' button)
# ---------------------------------------------------------------------------

#' Built-in demo data so the app is usable without a file.
demo_data <- function(n = 400, k = 15, m = 3, seed = 2024, dif = TRUE) {
  set.seed(seed)
  tau <- c(-1.4, 0, 1.4)[seq_len(m)]; tau <- tau - mean(tau)
  th  <- rnorm(n, 0.2, 1.4)
  de  <- seq(-1.8, 1.8, length.out = k); de <- de - mean(de)
  sex <- sample(c("Female", "Male"), n, TRUE)
  X <- matrix(0L, n, k)
  for (i in seq_len(k)) {
    shift <- if (dif && i == 4) ifelse(sex == "Male", 0.9, 0) else 0
    ct <- c(0, cumsum(tau))
    psi <- sapply(0:m, function(kk) kk * (th - de[i] - shift) - ct[kk + 1])
    P <- exp(psi - apply(psi, 1, max)); P <- P / rowSums(P)
    X[, i] <- rowSums(runif(n) > t(apply(P, 1, cumsum)))
  }
  colnames(X) <- sprintf("Item%02d", seq_len(k))
  data.frame(PersonID = sprintf("P%03d", seq_len(n)), Sex = sex,
             Age = sample(18:70, n, TRUE), X, stringsAsFactors = FALSE)
}
