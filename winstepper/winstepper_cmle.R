# =============================================================================
# winstepper_cmle.R -- Conditional maximum likelihood (CML) estimation.
#
# JMLE estimates person and item parameters jointly, which makes the item
# estimates inconsistent: with a fixed test length L the item spread is inflated
# by roughly L/(L-1) and no amount of extra data removes it (WINSTEPS offers
# STBIAS= to patch this after the fact; this project does not implement it).
#
# CML removes the person parameters instead of estimating them. The person raw
# score is sufficient for theta under the Rasch model, so conditioning each
# response pattern on its own raw score leaves a likelihood that contains the
# item parameters alone. The resulting estimates are consistent as N grows with
# the test length fixed - the property JMLE lacks.
#
# No Shiny dependency: usable from the command line, like rasch_engine.R.
# Source AFTER rasch_engine.R (it reuses .find_extremes(), .prox_start(),
# .probs_group(), .solve_measure_person(), .solve_measure_item(),
# .model_matrices()).
#
# References
#   Andersen, E.B. (1972). The numerical solution of a set of conditional
#     estimation equations. JRSS-B 34, 42-54.
#   Fischer, G.H. (1981). On the existence and uniqueness of maximum-likelihood
#     estimates in the Rasch model. Psychometrika 46, 59-77.
#   Fischer, G.H. & Molenaar, I.W. (1995). Rasch Models, ch. 3.
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Elementary symmetric functions
# ---------------------------------------------------------------------------

#' Pairwise log(exp(a) + exp(b)), safe when either side is -Inf.
.lse2 <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[!is.finite(m)] <- -Inf              # both terms were zero on the raw scale
  out
}

#' Elementary symmetric functions of the item parameters, in the log domain.
#'
#' gamma_r is the sum, over every response pattern whose total score is r, of
#' prod_i eps_{i,y_i} with eps_ik = exp(-b_ik). It is built by the summation
#' algorithm: start from the polynomial [1] and convolve one item in at a time,
#'
#'     G^(i)_r = sum_k eps_ik * G^(i-1)_{r-k}
#'
#' Everything is kept as logs and combined with log-sum-exp. That is not
#' decoration: with 17 items on a 0-9 scale the raw gamma values span hundreds
#' of orders of magnitude across r = 0..153 and overflow immediately.
#'
#' @param loge list, one numeric vector per item, element k+1 = -b_ik for
#'   k = 0..m_i (so the first element is always 0, because b_i0 = 0).
#' @return numeric vector of length sum(m_i) + 1, element r+1 = log(gamma_r).
esf_log <- function(loge) {
  G <- 0                                  # log gamma for the empty test: [1]
  for (v in loge) {
    m <- length(v) - 1L
    R <- length(G) - 1L
    new <- rep(-Inf, R + m + 1L)
    for (k in 0:m) {
      idx <- (k + 1L):(k + R + 1L)
      new[idx] <- .lse2(new[idx], v[[k + 1L]] + G)
    }
    G <- new
  }
  G
}

#' Brute-force gamma by enumerating every response pattern. O(prod(m_i + 1)) -
#' only for testing esf_log() on tiny problems.
esf_log_brute <- function(loge) {
  grids <- lapply(loge, function(v) seq_along(v) - 1L)
  pat <- as.matrix(expand.grid(grids))
  w <- vapply(seq_len(nrow(pat)), function(j)
    sum(vapply(seq_along(loge), function(i) loge[[i]][[pat[j, i] + 1L]], numeric(1))),
    numeric(1))
  r <- rowSums(pat)
  R <- sum(vapply(loge, length, integer(1)) - 1L)
  vapply(0:R, function(rr) {
    v <- w[r == rr]
    if (!length(v)) -Inf else { mx <- max(v); mx + log(sum(exp(v - mx))) }
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# 2. Parameter packing
#
# The model is the engine's own parameterisation,
#
#     b_ik = k * delta_i + cumsum(tau_g)[k],      b_i0 = 0
#
# so the result drops straight into fit$delta / fit$tau and every downstream
# table and figure works unchanged. Identification is sum(delta) = 0 and
# sum(tau_g) = 0 per group, imposed by construction rather than by projection:
# the last element of each block is minus the sum of the others, which leaves
# the optimiser no flat directions to wander along.
# ---------------------------------------------------------------------------

.cml_pack <- function(delta, tau_list, gl) {
  c(delta[-length(delta)],
    unlist(lapply(gl, function(g) {
      t_ <- tau_list[[g]]
      if (length(t_) < 2L) numeric(0) else t_[-length(t_)]
    }), use.names = FALSE))
}

.cml_unpack <- function(par, I, gl, mcc) {
  d <- par[seq_len(I - 1L)]
  delta <- c(d, -sum(d))
  pos <- I - 1L
  tau_list <- vector("list", length(gl)); names(tau_list) <- gl
  for (g in gl) {
    m <- mcc[[g]]
    if (m < 2L) {
      tau_list[[g]] <- 0                  # dichotomy: absorbed into delta
    } else {
      q <- par[pos + seq_len(m - 1L)]; pos <- pos + m - 1L
      tau_list[[g]] <- c(q, -sum(q))
    }
  }
  list(delta = delta, tau = tau_list)
}

#' -b_ik for every item, given delta and tau. The input to esf_log().
.cml_loge <- function(delta, tau_list, groups, mcc) {
  lapply(seq_along(delta), function(i) {
    g <- groups[i]
    -(0:mcc[[g]] * delta[i] + c(0, cumsum(tau_list[[g]])))
  })
}

# ---------------------------------------------------------------------------
# 3. Conditional log-likelihood and its gradient
# ---------------------------------------------------------------------------

#' Sufficient statistics, grouped by missing-data pattern.
#'
#' The conditioning set is the items a person actually answered, so persons with
#' different missingness patterns have different gamma polynomials and each
#' pattern must contribute its own. Complete data collapses to a single pattern.
.cml_patterns <- function(X, mcc, groups) {
  M <- !is.na(X)
  key <- apply(M, 1L, function(z) paste0(as.integer(z), collapse = ""))
  cap <- mcc[groups]
  out <- lapply(split(seq_len(nrow(X)), key), function(rows) {
    cols <- which(M[rows[1L], ])
    r <- as.integer(round(rowSums(X[rows, cols, drop = FALSE])))
    R <- as.integer(sum(cap[cols]))
    list(cols = cols, n_r = tabulate(r + 1L, nbins = R + 1L), R = R,
         n_persons = length(rows))
  })
  names(out) <- NULL
  out
}

#' Observed category counts N_ik (rows = items, cols = category 0..max).
.cml_counts <- function(X, mcc, groups) {
  mmax <- max(mcc)
  N <- matrix(0, ncol(X), mmax + 1L)
  for (i in seq_len(ncol(X))) {
    v <- X[, i]; v <- as.integer(round(v[!is.na(v)]))
    N[i, ] <- tabulate(v + 1L, nbins = mmax + 1L)
  }
  N
}

#' Negative conditional log-likelihood, and (as an attribute) its gradient.
#'
#' The gradient has a clean form. Writing E[N_ik] for the conditionally expected
#' category count,
#'
#'     d logL / d b_ik = -(N_ik - E[N_ik])
#'
#' so by the chain rule through b_ik = k*delta_i + cumsum(tau)[k],
#'
#'     d logL / d delta_i = -(observed item score      - expected item score)
#'     d logL / d tau_gj  = -(observed count of X >= j - expected count)
#'
#' i.e. CML is solved by matching exactly the sufficient statistics JMLE already
#' uses - the item raw score and the cumulative category counts. That identity
#' is what test 3 in test_cmle.R checks at the solution.
#'
#' E[N_ik] needs G^(-i), the ESF of every item except i. It is obtained by
#' re-running the convolution without item i, not by the "difference algorithm"
#' that divides item i back out: deconvolution is faster but numerically
#' unstable, and at these test lengths the honest version costs nothing.
.cml_negll <- function(par, I, gl, mcc, groups, pats, N_ik, S_i, C_gj,
                       gradient = TRUE) {
  pr <- .cml_unpack(par, I, gl, mcc)
  delta <- pr$delta; tau_list <- pr$tau
  loge <- .cml_loge(delta, tau_list, groups, mcc)

  ## --- the -sum N_ik b_ik term (does not depend on the pattern) ------------
  ll <- 0
  for (i in seq_len(I)) {
    v <- loge[[i]]
    ll <- ll + sum(N_ik[i, seq_along(v)] * v)     # loge = -b, so this is -N*b
  }

  E_ik <- if (gradient) matrix(0, I, max(mcc) + 1L) else NULL

  for (p in pats) {
    cols <- p$cols
    lg <- loge[cols]
    G <- esf_log(lg)
    keep <- p$n_r > 0
    ll <- ll - sum(p$n_r[keep] * G[keep])

    if (gradient) {
      for (jj in seq_along(cols)) {
        i <- cols[jj]
        Gmi <- if (length(lg) == 1L) 0 else esf_log(lg[-jj])
        mi <- length(lg[[jj]]) - 1L
        Rm <- length(Gmi) - 1L
        for (k in 0:mi) {
          rr <- k:(Rm + k)                        # scores this (i,k) can occur at
          nr <- p$n_r[rr + 1L]
          ok <- nr > 0
          if (!any(ok)) next
          E_ik[i, k + 1L] <- E_ik[i, k + 1L] +
            sum(nr[ok] * exp(lg[[jj]][[k + 1L]] + Gmi[ok] - G[rr[ok] + 1L]))
        }
      }
    }
  }

  if (!gradient) return(-ll)

  ## --- residuals in the two sufficient statistics --------------------------
  kk <- 0:(ncol(E_ik) - 1L)
  E_score <- as.vector(E_ik %*% kk)                       # expected item scores
  g_delta <- -(S_i - E_score)

  g_tau <- lapply(gl, function(g) {
    m <- mcc[[g]]
    if (m < 2L) return(numeric(0))
    ii <- which(groups == g)
    Ec <- vapply(1:m, function(j) sum(E_ik[ii, (j + 1L):(m + 1L), drop = FALSE]), numeric(1))
    -(C_gj[[g]] - Ec)
  })
  names(g_tau) <- gl

  ## --- chain rule through the sum-to-zero parameterisation ------------------
  gr <- c(g_delta[-I] - g_delta[I],
          unlist(lapply(gl, function(g) {
            v <- g_tau[[g]]
            if (length(v) < 2L) numeric(0) else v[-length(v)] - v[[length(v)]]
          }), use.names = FALSE))

  out <- -ll
  attr(out, "gradient") <- -gr                            # of the NEGATIVE ll
  attr(out, "E_ik") <- E_ik
  out
}

# ---------------------------------------------------------------------------
# 4. The estimator
# ---------------------------------------------------------------------------

#' Conditional maximum likelihood estimation of the (polytomous) Rasch model.
#'
#' Handles dichotomous, rating-scale (RSM) and partial-credit (PCM) data, and
#' any mix of item groups, exactly as rasch_jmle() does.
#'
#' @param prep output of rasch_prep()
#' @param maxit maximum BFGS iterations
#' @param conv convergence tolerance: the largest absolute score residual
#'   accepted as converged, i.e. observed minus conditionally expected
#'   sufficient statistic, in units of responses. 1e-3 of a response is
#'   negligible against the thousands that go into each statistic.
#' @param extreme_adj score adjustment for extreme persons/items
#'   (WINSTEPS EXTRSCORE=, default 0.3)
#' @param start starting values. "jmle" runs a short, loose JMLE first - the
#'   conditional likelihood is strictly concave for connected data, so this
#'   changes only the speed, never the answer. "prox" uses the normal
#'   approximation and flat thresholds.
#' @param center "items" (sum(delta) = 0, the CML identification) or "persons"
#'   (shift both margins so the person mean is 0 afterwards)
#' @return a `raschfit` object with the same structure rasch_jmle() returns, so
#'   every table and figure downstream works unchanged.
rasch_cmle <- function(prep, maxit = 400, conv = 1e-3, extreme_adj = 0.3,
                       start = c("jmle", "prox"), center = c("items", "persons"),
                       verbose = FALSE) {
  start <- match.arg(start); center <- match.arg(center)
  # A blank numericInput arrives as NA, and every scalar `if` below would then
  # throw "missing value where TRUE/FALSE needed" instead of estimating.
  if (!isTRUE(is.finite(conv)) || conv <= 0) conv <- 1e-3
  if (!isTRUE(is.finite(maxit)) || maxit < 1) maxit <- 400
  X <- prep$X; groups <- prep$groups; max_cat <- prep$max_cat
  N <- nrow(X); I <- ncol(X)

  ex <- .find_extremes(X, groups, max_cat)
  keep_p <- ex$keep_p; keep_i <- ex$keep_i
  if (sum(keep_p) < 2 || sum(keep_i) < 2)
    stop("Fewer than two non-extreme persons or items remain; the data cannot be calibrated.")

  Xc  <- X[keep_p, keep_i, drop = FALSE]
  gc_ <- groups[keep_i]
  gl  <- unique(gc_)
  mcc <- max_cat[gl]
  Ic  <- ncol(Xc)

  ## --- sufficient statistics ------------------------------------------------
  N_ik <- .cml_counts(Xc, mcc, gc_)
  S_i  <- colSums(Xc, na.rm = TRUE)                 # observed item raw scores
  C_gj <- lapply(gl, function(g) {
    v <- Xc[, gc_ == g, drop = FALSE]
    vapply(1:mcc[[g]], function(k) sum(v >= k, na.rm = TRUE), numeric(1))
  })
  names(C_gj) <- gl
  pats <- .cml_patterns(Xc, mcc, gc_)

  ## --- starting values ------------------------------------------------------
  if (start == "jmle") {
    j0 <- try(rasch_jmle(list(X = Xc, groups = gc_, max_cat = mcc),
                         maxit = 60, conv = 1e-2, rconv = 1, extreme_adj = extreme_adj,
                         center = "items"), silent = TRUE)
    if (inherits(j0, "try-error")) start <- "prox"
  }
  if (start == "jmle") {
    d0 <- j0$delta; d0[!is.finite(d0)] <- 0; d0 <- d0 - mean(d0)
    t0 <- lapply(gl, function(g) {
      v <- j0$tau[[g]]
      if (mcc[[g]] < 2L || !all(is.finite(v))) rep(0, max(mcc[[g]], 1L)) else v - mean(v)
    })
    names(t0) <- gl
  } else {
    d0 <- .prox_start(Xc, gc_, mcc)$delta
    d0[!is.finite(d0)] <- 0; d0 <- d0 - mean(d0)
    t0 <- lapply(gl, function(g) rep(0, max(mcc[[g]], 1L)))
    names(t0) <- gl
  }

  ## --- optimise -------------------------------------------------------------
  trace <- numeric(0); tracing <- TRUE
  fn <- function(p) as.numeric(.cml_negll(p, Ic, gl, mcc, gc_, pats, N_ik, S_i, C_gj,
                                          gradient = FALSE))
  gr <- function(p) {
    g <- attr(.cml_negll(p, Ic, gl, mcc, gc_, pats, N_ik, S_i, C_gj, gradient = TRUE),
              "gradient")
    if (tracing) trace <<- c(trace, max(abs(g)))
    g
  }

  # Convergence is judged on the CML equations themselves - observed equals
  # conditionally expected sufficient statistic - not on optim's return code,
  # which only says its line search stopped moving.
  p0 <- .cml_pack(d0, t0, gl)
  op <- stats::optim(p0, fn, gr, method = "BFGS",
                     control = list(maxit = maxit, reltol = 1e-15,
                                    trace = if (verbose) 1L else 0L))
  nit <- as.integer(op$counts[["gradient"]])
  if (!length(trace)) trace <- max(abs(gr(op$par)))
  tracing <- FALSE            # keep the polish and Hessian out of the trace

  ## --- Newton polish, then the observed information -------------------------
  # BFGS stops on relative change in the function value, which typically leaves
  # the score residual near 1e-4 - fine for the estimates, too loose to call the
  # CML equations solved. The Hessian is needed for the standard errors anyway,
  # so spend it twice: one or two Newton steps finish the job, and the second
  # evaluation is at the point actually reported.
  par_hat <- op$par
  np <- length(par_hat)
  H <- try(stats::optimHess(par_hat, fn, gr), silent = TRUE)
  if (!inherits(H, "try-error")) {
    moved <- FALSE
    for (k in seq_len(3L)) {
      g <- gr(par_hat)
      if (isTRUE(max(abs(g)) < conv * 1e-3)) break
      st <- try(solve(H, g), silent = TRUE)
      if (inherits(st, "try-error")) break
      cand <- par_hat - st
      if (!all(is.finite(cand)) || !is.finite(fn(cand)) || fn(cand) > fn(par_hat)) break
      par_hat <- cand; nit <- nit + 1L; moved <- TRUE
    }
    # Only re-evaluate if a step was taken; each Hessian costs 2*np gradients.
    if (moved) H <- try(stats::optimHess(par_hat, fn, gr), silent = TRUE)
  }
  max_res <- max(abs(gr(par_hat)))
  converged <- isTRUE(max_res < conv)
  trace <- c(trace, max_res)
  llc <- -fn(par_hat)
  if (verbose)
    message(sprintf("CMLE: logL = %.6f, max|score residual| = %.3e", llc, max_res))

  pr <- .cml_unpack(par_hat, Ic, gl, mcc)
  delta_c <- pr$delta; tau_list <- pr$tau

  ## --- standard errors from the observed information ------------------------
  V <- if (inherits(H, "try-error")) matrix(NA_real_, np, np) else
    tryCatch(solve(H), error = function(e) matrix(NA_real_, np, np))
  # Map the free parameters onto (delta, tau): the last element of each block is
  # minus the sum of the others, so its variance is not simply a diagonal entry.
  A_d <- matrix(0, Ic, np)
  if (Ic > 1L) {
    A_d[cbind(seq_len(Ic - 1L), seq_len(Ic - 1L))] <- 1
    A_d[Ic, seq_len(Ic - 1L)] <- -1
  }
  se_delta_c <- sqrt(pmax(diag(A_d %*% V %*% t(A_d)), 0))
  pos <- Ic - 1L
  se_tau <- vector("list", length(gl)); names(se_tau) <- gl
  for (g in gl) {
    m <- mcc[[g]]
    if (m < 2L) { se_tau[[g]] <- numeric(0); next }
    A_t <- matrix(0, m, np)
    A_t[cbind(seq_len(m - 1L), pos + seq_len(m - 1L))] <- 1
    A_t[m, pos + seq_len(m - 1L)] <- -1
    se_tau[[g]] <- sqrt(pmax(diag(A_t %*% V %*% t(A_t)), 0))
    pos <- pos + m - 1L
  }

  ## --- back into full-length vectors ---------------------------------------
  Delta <- rep(NA_real_, I); Delta[keep_i] <- delta_c
  SEd   <- rep(NA_real_, I); SEd[keep_i]   <- se_delta_c
  Theta <- rep(NA_real_, N); SEt <- rep(NA_real_, N)

  ## --- person measures: items anchored at the CML estimates ------------------
  # CML estimates items only. Persons are then estimated by ML from their raw
  # score with the items held fixed - the same anchored solve rasch_jmle() uses
  # for its extreme cases, including the EXTRSCORE= adjustment. Persons sharing
  # a missingness pattern and a raw score share a measure, so the solve is
  # cached rather than repeated once per person.
  Mfull <- !is.na(X)
  cap   <- matrix(max_cat[groups][col(X)], N, I)
  cache <- new.env(parent = emptyenv())
  for (n in seq_len(N)) {
    obs <- Mfull[n, ] & !is.na(Delta)
    if (!any(obs)) next
    rmax <- sum(cap[n, obs]); r <- sum(X[n, obs])
    radj <- if (keep_p[n]) r else min(max(r, extreme_adj), rmax - extreme_adj)
    key <- paste0(paste0(as.integer(obs), collapse = ""), "|", radj)
    if (exists(key, envir = cache, inherits = FALSE)) {
      es <- get(key, envir = cache)
    } else {
      es <- .solve_measure_person(radj, Delta[obs], tau_list, groups[obs])
      assign(key, es, envir = cache)
    }
    Theta[n] <- es$measure; SEt[n] <- es$se
  }

  ## --- extreme items: anchored on the persons, as in JMLE -------------------
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

  if (center == "persons") {
    a <- mean(Theta, na.rm = TRUE)
    Theta <- Theta - a; Delta <- Delta - a
  }

  full <- .model_matrices(Theta, Delta, tau_list, groups)

  structure(list(
    X = X, mask = Mfull, groups = groups, max_cat = max_cat, cap = cap,
    theta = Theta, se_theta = SEt, delta = Delta, se_delta = SEd,
    tau = tau_list, se_tau = se_tau,
    E = full$E, W = full$W, C = full$C,
    keep_p = keep_p, keep_i = keep_i,
    extreme_persons = which(!keep_p), extreme_items = which(!keep_i),
    iterations = nit, converged = converged, change_history = trace,
    cml = list(logLik = llc, max_residual = max_res, optim_code = op$convergence,
               n_patterns = length(pats), n_free = np),
    settings = list(method = "CMLE", maxit = maxit, conv = conv,
                    extreme_adj = extreme_adj, center = center, start = start),
    person_id = rownames(X), item_id = colnames(X)
  ), class = "raschfit")
}
