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
  if (!length(predictors)) stop("No predictors were selected.")
  missing_cols <- setdiff(predictors, names(df_reg))
  if (length(missing_cols))
    stop("Predictor(s) not found in the covariates: ",
         paste(missing_cols, collapse = ", "))

  # AUDIT (P1.4b, CORRECTED at P4.3). The P1.4b note claimed that
  # "reformulate() takes the names as data and quotes them". THAT IS FALSE, and
  # a reviewer was right to reject it. stats::reformulate() pastes its
  # termlabels together and PARSES the result -- R's own documentation says the
  # labels must be syntactically valid names or already backquoted. Measured:
  #   reformulate("a b")              -> error, unexpected symbol
  #   reformulate("Age (years)")      -> ~Age(years), a FUNCTION CALL
  #   reformulate('I(cat("PWNED"))')  -> a formula holding a live call
  # So the P1.4b change swapped one text-pasting route for another and the
  # comment asserted a safety property the function does not have. Backquoting
  # would fix the first two cases and still breaks on a name containing a
  # backquote.
  #
  # The robust fix is not to put user text in a formula at all. The model is
  # built on an internal frame whose columns are named x1..xp -- names that
  # cannot be anything but names -- and the mapping back to the user's labels is
  # applied to the coefficient rows afterwards, as data. An uploaded column may
  # then be called anything at all, including something executable, and it is
  # still only ever a column.
  model_names <- paste0("x", seq_along(predictors))
  model_df <- df_reg[, predictors, drop = FALSE]
  names(model_df) <- model_names
  # drop unused factor levels: a level emptied by row filtering makes the design
  # rank deficient in a way the user cannot see in their own data
  for (nm in model_names)
    if (is.factor(model_df[[nm]])) model_df[[nm]] <- droplevels(model_df[[nm]])

  f <- stats::reformulate(model_names)
  X <- model.matrix(f, data = model_df)
  # Factor levels as fitted. Prediction on new data MUST use these, or a new
  # frame that happens to contain only one level of a factor silently builds a
  # design with the wrong number of columns.
  xlev <- stats::.getXlevels(terms(f), model_df)

  # Map the safe names on the design columns back to the user's labels. Done by
  # string surgery on the COLUMN NAMES only -- nothing here is ever parsed.
  restore_labels <- function(z) {
    for (k in rev(seq_along(predictors)))
      z <- sub(paste0("^", model_names[k]), predictors[k], z)
    z
  }
  colnames(X) <- restore_labels(colnames(X))

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

  # AUDIT (P4.4): chol2inv(qr.R(qrX)) is (X'X)^-1 only if the decomposition did
  # not permute the columns; under a pivoted QR, R belongs to the PERMUTED
  # design and the inverse would be silently mis-labelled, so the SE of one
  # coefficient would be reported for another.
  #
  # A reviewer flagged this as a live bug. As stated it is overstated for the
  # default LINPACK qr(): dqrdc2 only cycles a column to the back when its
  # reduced norm falls below the rank tolerance, and it decrements the rank when
  # it does -- so rank == ncol(X) implies pivot == 1:p, and the rank check above
  # has already stopped otherwise. Measured on ill-conditioned designs at three
  # conditioning levels: whenever pivoting would have mattered, the rank test
  # fired first. But that is an UNDOCUMENTED invariant of one LAPACK/LINPACK
  # path, not a guarantee, and it costs one comparison to stop relying on it.
  if (!identical(qrX$pivot, seq_len(ncol(X)))) {
    # Undo the permutation rather than trusting it did not happen.
    R  <- qr.R(qrX)
    iv <- chol2inv(R)
    ord <- order(qrX$pivot)
    xtx_inv <- iv[ord, ord, drop = FALSE]
  } else {
    xtx_inv <- chol2inv(qr.R(qrX))
  }
  dimnames(xtx_inv) <- list(colnames(X), colnames(X))
  projector <- xtx_inv %*% t(X)
  coefs <- projector %*% Y
  fitted_vals <- X %*% coefs
  residuals <- Y - fitted_vals
  n <- nrow(Y); p <- ncol(X)

  # AUDIT (P4.5): a design with n == p is full rank, so the check above passes,
  # and then every downstream quantity is meaningless: sigma2 = 0/0, the SE is
  # Inf and pt(df = 0) is NaN. Measured on a 3x3 design: sigma2 Inf, SE Inf,
  # p NaN, with no warning the user would connect to the cause. There is no
  # residual variance to estimate, so the fit cannot support inference at all
  # and saying so is the only honest option.
  if (n <= p)
    stop(sprintf(
      paste("This design has %d subject%s and %d coefficient%s, leaving %d residual",
            "degrees of freedom. Standard errors and p-values need at least one.",
            "Use fewer predictors, or collapse factor levels, and re-run."),
      n, if (n == 1) "" else "s", p, if (p == 1) "" else "s", n - p))

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
  
  # AUDIT (P4.6): r2_t <- 1 - rss/tss with no guard. A time point at which the
  # response is constant across subjects has tss = 0, so R2 came out NaN or
  # -Inf and was then plotted and averaged. The cosinor code already refuses
  # this case (fck_r_squared returns NA when SST is at the noise floor); the
  # same rule applies here. NA means "no variance to explain", not "zero fit".
  rss <- colSums(residuals^2); y_bar <- colMeans(Y)
  tss <- colSums(sweep(Y, 2, y_bar)^2)
  tss_floor <- .Machine$double.eps * max(1, max(abs(Y), na.rm = TRUE))^2 * nrow(Y)
  r2_t <- ifelse(tss > tss_floor, 1 - (rss / tss), NA_real_)
  
  list(beta.hat = coefs, fitted.values = fitted_vals, resid = residuals,
       beta.se = se_mat, beta.se.boot = se_boot,
       beta.p = p_values, beta.p.raw = p_values_raw,
       terms = terms(f), r2_t = r2_t,
       method = if(isTRUE(use_bootstrap))
         "OLS (analytical SE and p-values; bootstrap percentile CI)"
       else "OLS (analytical SE and p-values)",
       inference = "analytic-t-fdr", df_resid = n - p, n_obs = n,
       # AUDIT (P5.1): the terms object above refers to x1..xp, NOT to the
       # user's column names -- that is the whole point of P4.3. Anything that
       # later builds a design matrix for PREDICTION has to rename its new data
       # the same way first. It did not, so from P4.3 until now every FoSR OLS
       # prediction curve failed with "object 'x1' not found" and was swallowed
       # by a tryCatch that returns NULL, showing the user a blank curve rather
       # than an error. This is a regression I introduced in the fix for P4.3.
       #
       # The mapping and the level set now travel WITH the fit, and
       # fck_fosr_design() below is the one place that turns new data into a
       # design matrix. Fit and prediction cannot drift again because they are
       # the same code path.
       model_names = model_names, predictor_names = predictors, xlevels = xlev,
       xtx_inv = xtx_inv, sigma2 = sigma2,
       boot_ci_lower = boot_ci_lower, boot_ci_upper = boot_ci_upper,
       se_ratio_range = if(isTRUE(use_bootstrap))
         c(se_ratio_min, se_ratio_max) else NULL,
       n_boot = if(use_bootstrap) n_boot else NULL)
}


# ==========================================================================
# fck_fosr_design(fit, newdata)
#
# The ONE way to build a design matrix for a fitted pointwise-OLS model. Takes
# new data keyed by the USER's column names, renames it to the internal
# x1..xp the model was fitted on, and applies the factor levels recorded at fit
# time. Used by the fit itself, by the interactive prediction curve, by the
# min/max reference curves, and by the exported script -- so there is no second
# implementation to fall out of step (P5.1).
# ==========================================================================
fck_fosr_design <- function(fit, newdata) {
  if (is.null(fit$model_names) || is.null(fit$predictor_names))
    stop("This fit carries no predictor-name mapping; it predates P5.1 and cannot be used for prediction. Re-run the FoSR fit.")

  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE,
                           check.names = FALSE)
  missing_cols <- setdiff(fit$predictor_names, names(newdata))
  if (length(missing_cols))
    stop("New data is missing predictor(s): ", paste(missing_cols, collapse = ", "))

  md <- newdata[, fit$predictor_names, drop = FALSE]
  names(md) <- fit$model_names

  # Re-apply the levels the model was fitted with, so a one-row prediction
  # frame produces the same columns as the full design.
  for (nm in names(fit$xlevels)) {
    if (!is.null(md[[nm]])) {
      lv <- fit$xlevels[[nm]]
      bad <- setdiff(as.character(md[[nm]]), lv)
      if (length(bad))
        stop("Level(s) not seen when the model was fitted: ",
             paste(unique(bad), collapse = ", "))
      md[[nm]] <- factor(as.character(md[[nm]]), levels = lv)
    }
  }

  X <- stats::model.matrix(stats::delete.response(fit$terms), data = md,
                           xlev = fit$xlevels)

  # The same terms and levels give the same column order as the fit, so the
  # coefficient row names carry straight over. Check rather than assume: a
  # mismatch here would multiply the right numbers in the wrong order.
  if (!is.null(fit$beta.hat) && ncol(X) != nrow(fit$beta.hat))
    stop(sprintf(
      "Design matrix has %d columns but the fit has %d coefficients. The prediction data does not match the fitted model.",
      ncol(X), nrow(fit$beta.hat)))
  if (!is.null(fit$beta.hat)) colnames(X) <- rownames(fit$beta.hat)
  X
}
