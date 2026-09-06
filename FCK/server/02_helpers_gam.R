# ==========================================================================
# server/02_helpers_gam.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 702-774  (GAM prediction helpers used by FoSR)
# ==========================================================================
  # ==============================================================================
  # HELPER FUNCTION: Build GAM prediction dataframe
  # ==============================================================================
  # AUDIT (P8.1). The fitted GAM's formula refers to the INTERNAL names x1..xp
  # (P4.3: uploaded column names must never reach a parser). This helper built
  # its prediction frame with the USER's names, so predict.gam() could not find
  # its variables:
  #     ERROR: object 'x2' not found
  #     (and a warning: not all required variables have been supplied in newdata)
  # -- on every GAM prediction curve and both min/max reference curves.
  #
  # This is the SAME defect as P5.1, which was fixed for the pointwise-OLS
  # branch by making the mapping travel with the fit and routing every design
  # through one builder. The GAM branch was not fixed with it. The lesson is in
  # the fix: `model_names` is now a required argument, so a caller that has not
  # thought about the mapping cannot call this function at all.
  #
  # `predictors` names the CONTROLS the values come from (input_vals is keyed by
  # them); `model_names` names the COLUMNS the model wants. Levels are read from
  # long_data by the model name, which fck_fit_fosr_gam() aliases alongside the
  # originals.
  build_gam_pred_df <- function(time_vec, predictors, model_names, long_data,
                                factor_levels, input_vals) {
    if (missing(model_names) || is.null(model_names))
      stop("build_gam_pred_df(): model_names is required. The fitted GAM refers to ",
           "internal column names, not the user's; pass mod$gam_model_names.")
    if (length(model_names) != length(predictors))
      stop(sprintf("build_gam_pred_df(): %d predictors but %d model names.",
                   length(predictors), length(model_names)))

    # Create dataframe with time column
    pred_df <- data.frame(time = time_vec)
    n_t <- length(time_vec)

    for(j in seq_along(predictors)) {
      var    <- predictors[j]    # the control / user-facing name
      mvar   <- model_names[j]   # the name the fitted formula uses
      val <- input_vals[[var]]
      if(is.null(val)) return(NULL)

      # Get the column from long_data to check type and levels. Prefer the
      # model-named alias; fall back to the user name for a fit made before
      # P4.3 introduced the aliases.
      orig_col <- if (!is.null(long_data[[mvar]])) long_data[[mvar]] else long_data[[var]]

      if(is.numeric(orig_col)) {
        # Numeric variable: replicate scalar value
        pred_df[[mvar]] <- rep(as.numeric(val), n_t)
      } else if(is.factor(orig_col)) {
        # Factor: use EXACT levels from the fitted model's data
        lvls <- levels(orig_col)
        # Validate that val is in levels
        if(!(val %in% lvls)) {
          warning(paste("Value", val, "not in factor levels for", var))
          val <- lvls[1]  # Fall back to first level
        }
        pred_df[[mvar]] <- factor(rep(val, n_t), levels = lvls)
      } else {
        # Character or other: use stored factor levels if available
        if(!is.null(factor_levels) && !is.null(factor_levels[[var]])) {
          lvls <- factor_levels[[var]]
        } else {
          lvls <- sort(unique(as.character(orig_col)))
        }
        if(!(val %in% lvls)) {
          warning(paste("Value", val, "not in levels for", var))
          val <- lvls[1]
        }
        pred_df[[mvar]] <- factor(rep(val, n_t), levels = lvls)
      }
    }
    return(pred_df)
  }
  
  # ==============================================================================
  # HELPER FUNCTION: Safe GAM prediction with SE
  # ==============================================================================
  safe_gam_predict <- function(gam_obj, newdata) {
    result <- tryCatch({
      pred <- predict(gam_obj, newdata = newdata, type = "response", se.fit = TRUE)
      list(
        fit = as.numeric(pred$fit),
        se = as.numeric(pred$se.fit),
        se_available = TRUE,
        success = TRUE
      )
    }, error = function(e) {
      # Fallback: try without SE
      pred <- tryCatch({
        predict(gam_obj, newdata = newdata, type = "response")
      }, error = function(e2) {
        return(NULL)
      })
      if(is.null(pred)) {
        return(list(fit = NULL, se = NULL, success = FALSE))
      }
      # AUDIT (P9.4): this returned se = rep(0, ...). SE = 0 asserts an
      # estimate known EXACTLY -- it is the strongest possible claim, returned
      # here precisely when the standard error could not be computed at all.
      # Downstream that produced a confidence band of zero width, which reads
      # as a perfectly determined prediction. NA says what is true, and the
      # caller omits the band. Same rule as P5.9 for the GAM coefficient
      # curves: unknown is not zero.
      list(
        fit = as.numeric(pred),
        se = rep(NA_real_, length(pred)),
        se_available = FALSE,
        success = TRUE
      )
    })
    return(result)
  }
