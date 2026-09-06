
# ==========================================================================
# fck_fit_fosr_gam(Y, df_reg, predictors, k = 10, notify = NULL)
#
# The smoothed (GAM) function-on-scalar estimator, as a named function.
#
# AUDIT (P5.8). The pointwise-OLS estimator became a kernel at P3.4/P4.3 so
# that the Export tab could deparse the same object the GUI runs. The GAM
# branch was left inline, and the export kept a hand-written reconstruction of
# it that emitted
#     gam_formula <- as.formula(paste('value ~ s(time) + ',
#                     sprintf('s(time, by = %s)', fosr_predictors), ...))
# which is wrong twice over: it pastes UPLOADED COLUMN NAMES back into parsed
# formula text -- the defect P4.3 removed from the live estimator -- and it
# uses a different GAM specification from the app (no basis or k, and the
# app's factor terms are absent). So the exported GAM script could not be
# expected to reproduce the numbers on screen, and could break outright on a
# column called "Age (years)".
#
# Same treatment as the OLS branch: one estimator, called by the GUI and
# written verbatim into the export. It fits on internal names x1..xp and
# restores the user's labels afterwards, takes no Shiny input, and its notify
# sink defaults to a no-op.
# ==========================================================================

fck_fit_fosr_gam <- function(Y, df_reg, predictors, k = 10, notify = NULL) {
  if (is.null(notify)) notify <- function(msg, ...) invisible(NULL)
  Y <- as.matrix(Y)
  if (!length(predictors)) stop("No predictors were selected.")
  missing_cols <- setdiff(predictors, names(df_reg))
  if (length(missing_cols))
    stop("Predictor(s) not found in the covariates: ",
         paste(missing_cols, collapse = ", "))

  n_time      <- ncol(Y)
  time_points <- seq(0, 1, length.out = n_time)

  
  # IMPORTANT: Store original factor levels BEFORE any subsetting
  # This is needed for prediction to work correctly
  orig_factor_levels <- list()
  for(v in predictors) {
    if(is.factor(df_reg[[v]])) {
      orig_factor_levels[[v]] <- levels(df_reg[[v]])
    } else if(!is.numeric(df_reg[[v]])) {
      orig_factor_levels[[v]] <- sort(unique(as.character(df_reg[[v]])))
    }
  }
  
  # AUDIT (P4.3): this branch pasted UPLOADED COLUMN NAMES straight into
  # formula text and re-parsed it with as.formula(). A column called
  # "Age (years)" becomes a function call; anything executable is
  # executed. The pointwise branch had the same problem, and the P1.4b
  # note there wrongly claimed reformulate() solved it -- it does not,
  # it parses text too.
  #
  # Same fix as the pointwise branch: the model is fitted on internal
  # column names x1..xp, which cannot be anything but names, and the
  # user's labels are restored afterwards as data. `preds` keeps the
  # user-facing names for every readout; `mpred` is what the formula sees.
  preds <- predictors
  mpred <- paste0("x", seq_along(preds))

  df_reg$id_temp <- 1:nrow(df_reg)
  long_cov <- df_reg[rep(seq_len(nrow(df_reg)), each = n_time), ]
  long_data <- long_cov
  # the safe aliases, alongside the original columns (which stay for
  # display and for any downstream code that indexes by the user's name)
  for (k in seq_along(preds)) {
    v <- long_data[[preds[j]]]
    if (is.factor(v)) v <- droplevels(v)
    long_data[[mpred[j]]] <- v
  }
  long_data$time <- rep(time_points, times = nrow(df_reg))
  long_data$Y_val <- as.vector(t(Y))

  # Build Formula: Y ~ s(time) + s(time, by = x_k) ...
  # P5.8: the spline dimension is the k ARGUMENT, not a literal. It used to be
  # hard-coded at 10 in three places, and the export emitted s(time, by = x)
  # with no basis or k at all -- a different model from the one fitted.
  gam_formula_str <- sprintf("Y_val ~ s(time, bs = 'ps', k = %d)", k)
  for(j in seq_along(preds)) {
    if(is.numeric(df_reg[[preds[j]]])) {
      gam_formula_str <- paste0(gam_formula_str, " + s(time, by = ", mpred[j], ", bs='ps', k=", k, ")")
    } else {
      # Factor interaction
      gam_formula_str <- paste0(gam_formula_str, " + ", mpred[j],
                                " + s(time, by = ", mpred[j], ", bs='ps', k=", k, ")")
    }
  }

  gam_fit <- mgcv::gam(stats::as.formula(gam_formula_str), data = long_data, method = "REML")
  
  # Reconstruct Beta(t) for Visualization (Approximation)
  pred_grid <- data.frame(time = time_points)
  beta_hat <- matrix(0, nrow = length(preds) + 1, ncol = n_time)
  rownames(beta_hat) <- c("(Intercept)", preds)
  
  # Intercept approx. The prediction frames are built on the MODEL names,
  # because that is what the fitted formula refers to (P4.3).
  d_int <- pred_grid
  for(j in seq_along(preds)) {
    mv <- mpred[j]
    if(is.numeric(df_reg[[preds[j]]])) d_int[[mv]] <- 0
    else {
      # Use proper factor levels from long_data (which is what GAM was trained on)
      orig_col <- long_data[[mv]]
      if(is.factor(orig_col)) {
        d_int[[mv]] <- factor(rep(levels(orig_col)[1], n_time), levels = levels(orig_col))
      } else {
        lvls <- sort(unique(as.character(orig_col)))
        d_int[[mv]] <- factor(rep(lvls[1], n_time), levels = lvls)
      }
    }
  }
  beta_hat[1, ] <- predict(gam_fit, newdata = d_int, type = "response")
  
  # AUDIT (P1.4d): this loop populated beta_hat only when the predictor was
  # NUMERIC. A factor predictor was fitted by the GAM, appeared in the
  # results table, and its coefficient curve stayed a row of zeros -- a
  # null effect displayed as though it had been estimated. A factor's
  # effect is a contrast against the reference level, so there is one
  # curve per non-reference level, not one per variable.
  beta_rows <- list("(Intercept)" = beta_hat[1, ])
  for(j in seq_along(preds)) {
    v <- preds[j]; mv <- mpred[j]   # v names the ROW, mv drives the model
    if(is.numeric(df_reg[[v]])) {
      d_0 <- d_int; d_0[[mv]] <- 0
      d_1 <- d_int; d_1[[mv]] <- 1
      beta_rows[[v]] <- as.vector(predict(gam_fit, newdata = d_1) -
                                  predict(gam_fit, newdata = d_0))
    } else {
      oc <- long_data[[mv]]
      lv <- if (is.factor(oc)) levels(oc) else sort(unique(as.character(oc)))
      if (length(lv) >= 2) {
        d_ref <- d_int
        d_ref[[mv]] <- factor(rep(lv[1], n_time), levels = lv)
        base <- as.vector(predict(gam_fit, newdata = d_ref))
        for (l in lv[-1]) {
          d_l <- d_int
          d_l[[mv]] <- factor(rep(l, n_time), levels = lv)
          beta_rows[[paste0(v, l)]] <-
            as.vector(predict(gam_fit, newdata = d_l)) - base
        }
      }
    }
  }
  beta_hat <- do.call(rbind, beta_rows)
  rownames(beta_hat) <- names(beta_rows)
  
  fitted_vec <- predict(gam_fit, newdata = long_data)
  fitted_vals <- matrix(fitted_vec, nrow = nrow(Y), ncol = n_time, byrow = TRUE)
  residuals <- Y - fitted_vals
  
  # P4.6: same zero-variance guard as the pointwise branch.
  rss <- colSums(residuals^2); y_bar <- colMeans(Y)
  tss <- colSums(sweep(Y, 2, y_bar)^2)
  tss_floor <- .Machine$double.eps * max(1, max(abs(Y), na.rm = TRUE))^2 * nrow(Y)
  r2_t <- ifelse(tss > tss_floor, 1 - (rss / tss), NA_real_)

  # AUDIT (P5.9): these were beta_hat*0, commented "placeholders". Zero is
  # not a placeholder for an unknown standard error -- SE = 0 asserts an
  # estimate known exactly, and p = 0 asserts overwhelming significance
  # at every time point. Nothing downstream currently plots them as such,
  # but a numeric zero in a field named beta.p is one careless read away
  # from a table of p < 0.001 across the whole curve. NA says what is
  # true: this quantity was not computed.
  na_like <- matrix(NA_real_, nrow(beta_hat), ncol(beta_hat),
                    dimnames = dimnames(beta_hat))

  list(beta.hat = beta_hat, fitted.values = fitted_vals, resid = residuals,
       beta.se = na_like, beta.p = na_like,
       inference = "none",
       inference_note = paste(
         "Pointwise inference on the coefficient CURVES is not computed for",
         "the GAM branch. The curves are prediction contrasts, not fitted",
         "coefficients, and no standard error is propagated through them.",
         "Use summary(gam_obj) for the term-level tests mgcv does report,",
         "or the pointwise-OLS method for coefficient-curve inference."),
       terms = terms(gam_fit), gam_model_names = mpred, 
       r2_t = r2_t, method = "Smoothed OLS (GAM)",
       gam_obj = gam_fit,
       gam_predictors = preds,
       gam_long_data = long_data,  # Store long_data for exact factor levels
       gam_factor_levels = orig_factor_levels,  # Store original factor levels
       gam_n_time = n_time, gam_k = k)
}
