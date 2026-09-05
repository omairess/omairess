# ==========================================================================
# server/07_helpers_fosr.R
#
# The pointwise-OLS function-on-scalar estimator, as a named function.
#
# AUDIT (P3.4). This code used to live inline in the run_fosr observer in
# server/70_fosr.R, and the Export tab wrote a SEPARATE, hand-written
# reconstruction of it into the generated script. The two drifted: after the
# P1.4a/P1.4b corrections the GUI used stats::reformulate() and a QR with a
# rank check, while the export still emitted
#     f       <- as.formula(paste('~', paste(fosr_predictors, collapse = ' + ')))
#     xtx_inv <- solve(crossprod(X))
# -- the exact two constructions those corrections removed -- and emitted no
# standard errors or p-values at all. So "Export analysis code" produced a
# script that silently reintroduced a fixed defect and could not reproduce the
# inference on screen.
#
# The fANOVA and cosinor exports already solved this by deparse()-ing the live
# function object into the script. That pattern needs a function to deparse.
# This is it: one estimator, called by the GUI and written verbatim into the
# export, so the two cannot diverge again.
#
# It depends on nothing but base R and stats, takes no Shiny input, and its
# progress/notify sinks default to no-ops so it runs unchanged in a plain
# session.
# ==========================================================================

fck_fit_fosr_ols <- function(Y, df_reg, predictors,
                             use_bootstrap = FALSE, n_boot = 500,
                             progress = NULL, notify = NULL) {
  if (is.null(progress)) progress <- function(frac, detail = NULL) invisible(NULL)
  if (is.null(notify))   notify   <- function(msg, ...) invisible(NULL)

  Y <- as.matrix(Y)

  # AUDIT (P1.4b): as.formula(paste(...)) built the model formula by pasting
  # UPLOADED COLUMN NAMES into text and re-parsing it. reformulate() takes
  # the names as data and quotes them, so a column called `a + b` or
  # anything else executable is a variable, not code.
  f <- stats::reformulate(predictors)
  X <- model.matrix(f, data = df_reg)

  # AUDIT (P1.4a): solve(crossprod(X)) inverts the normal equations. That
  # squares the condition number and gives no warning on a rank-deficient
  # design -- two collinear predictors, or a factor level that vanished
  # after row filtering, either error out or return nonsense. A QR of X
  # itself is the standard route and reports the rank.
  qrX <- qr(X)
  if (qrX$rank < ncol(X)) {
    dropped <- colnames(X)[-qrX$pivot[seq_len(qrX$rank)]]
    stop(sprintf(
      "The design matrix is rank deficient (rank %d of %d columns). Collinear or empty terms: %s. Remove or recode them and re-run.",
      qrX$rank, ncol(X), paste(dropped, collapse = ", ")))
  }
  xtx_inv <- chol2inv(qr.R(qrX))
  dimnames(xtx_inv) <- list(colnames(X), colnames(X))
  projector <- xtx_inv %*% t(X)
  coefs <- projector %*% Y
  fitted_vals <- X %*% coefs
  residuals <- Y - fitted_vals
  n <- nrow(Y); p <- ncol(X)
  sigma2 <- colSums(residuals^2) / (n - p)
  
  # AUDIT (P3.2): the analytical standard error is available in closed
  # form here for nothing -- xtx_inv is already formed and sigma2 is
  # already computed -- so it is now ALWAYS the basis of the test, and the
  # bootstrap supplies the interval it is actually good for. See the long
  # note in the bootstrap branch below for why the previous arrangement
  # (t statistic against a BOOTSTRAP SE, referred to t_{n-p}) was the
  # wrong pairing.
  se_analytic <- matrix(NA_real_, nrow = nrow(coefs), ncol = ncol(coefs))
  for(j in 1:nrow(coefs)) se_analytic[j, ] <- sqrt(xtx_inv[j, j] * sigma2)
  dimnames(se_analytic) <- dimnames(coefs)

  t_stats      <- coefs / se_analytic
  p_values_raw <- 2 * stats::pt(abs(t_stats), df = n - p, lower.tail = FALSE)
  # P1.4c: FDR across time points -- every coefficient curve is ~100
  # simultaneous tests and none of them were corrected.
  p_values <- p_values_raw
  for(j in seq_len(nrow(p_values)))
    p_values[j, ] <- p.adjust(p_values[j, ], method = "fdr")

  se_mat        <- se_analytic
  se_boot       <- NULL
  boot_ci_lower <- NULL
  boot_ci_upper <- NULL

  if(use_bootstrap) {
    notify("Bootstrapping...")
    B <- n_boot
    boot_betas <- array(NA, dim = c(B, nrow(coefs), ncol(coefs)))
    
    for(b in 1:B) {
      # Residual bootstrap: resample residuals, add to fitted
      resid_idx <- sample(1:n, n, replace = TRUE)
      Y_boot <- fitted_vals + residuals[resid_idx, ]
      boot_betas[b, , ] <- projector %*% Y_boot
      if(b %% 20 == 0) progress(20 / B, "Bootstrapping")
    }
    
    # SE of the bootstrap distribution, and the percentile interval.
    se_boot <- apply(boot_betas, c(2, 3), sd, na.rm = TRUE)
    boot_ci_lower <- apply(boot_betas, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
    boot_ci_upper <- apply(boot_betas, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)
    dimnames(se_boot) <- dimnames(boot_ci_lower) <-
      dimnames(boot_ci_upper) <- dimnames(coefs)

    # AUDIT (P1.4c, revised at P3.2). Two successive wrong pairings here.
    #
    # Originally this block computed
    #     p <- 2 * mean(boot_betas[, j, k] <= 0)
    # under the comment "This is a proper bootstrap test". It is not: the
    # residuals are resampled around the FITTED ALTERNATIVE, so the
    # bootstrap distribution is centred near beta-hat, not at zero. That
    # inverts a percentile interval rather than sampling a null.
    #
    # P1.4c replaced it with beta-hat / SE_bootstrap referred to t_{n-p}.
    # Better, still mismatched. A t reference distribution is what you get
    # when the denominator is the ANALYTICAL sigma-hat, whose scaled
    # square is chi-square on n-p df and independent of beta-hat. The
    # bootstrap SE is neither of those things: it is a Monte Carlo
    # estimate of the same standard error, carrying its own noise of
    # order 1/sqrt(2B) in relative terms -- at B = 200 that is about 5%,
    # enough to move a borderline p-value -- and referring the ratio to
    # t_{n-p} is an assumption about a quantity that no longer has that
    # distribution. It also made a REPORTED P-VALUE DEPEND ON THE SEED.
    #
    # And it bought nothing. This is a RESIDUAL bootstrap: it resamples
    # rows of the residual matrix from the fitted homoscedastic model, so
    # it estimates exactly the quantity xtx_inv[j,j] * sigma2 estimates,
    # with extra variance. (A CASE/pairs bootstrap would be the one that
    # buys heteroscedasticity robustness. That is a different estimator
    # and is not what this control has ever run, so it is not silently
    # substituted here.)
    #
    # So: the test is the analytical one, computed above and unchanged by
    # this branch; the bootstrap contributes the percentile interval, and
    # its SE is kept only as a MISSPECIFICATION DIAGNOSTIC -- if
    # SE_bootstrap and SE_analytic disagree materially, the homoscedastic
    # normal model is in doubt and neither the interval nor the p-value
    # should be trusted at face value. That ratio is reported in the
    # summary.
    se_ratio_max <- suppressWarnings(
      max(se_boot / se_analytic, na.rm = TRUE))
    se_ratio_min <- suppressWarnings(
      min(se_boot / se_analytic, na.rm = TRUE))
  }
  
  rss <- colSums(residuals^2); y_bar <- colMeans(Y)
  tss <- colSums(sweep(Y, 2, y_bar)^2); r2_t <- 1 - (rss/tss)
  
  list(beta.hat = coefs, fitted.values = fitted_vals, resid = residuals,
       beta.se = se_mat, beta.se.boot = se_boot,
       beta.p = p_values, beta.p.raw = p_values_raw,
       terms = terms(f), r2_t = r2_t,
       method = if(isTRUE(use_bootstrap))
         "OLS (analytical SE and p-values; bootstrap percentile CI)"
       else "OLS (analytical SE and p-values)",
       inference = "analytic-t-fdr", df_resid = n - p, n_obs = n,
       xtx_inv = xtx_inv, sigma2 = sigma2,
       boot_ci_lower = boot_ci_lower, boot_ci_upper = boot_ci_upper,
       se_ratio_range = if(isTRUE(use_bootstrap))
         c(se_ratio_min, se_ratio_max) else NULL,
       n_boot = if(use_bootstrap) n_boot else NULL)
}
