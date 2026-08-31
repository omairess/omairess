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
  build_gam_pred_df <- function(time_vec, predictors, long_data, factor_levels, input_vals) {
    # Create dataframe with time column
    pred_df <- data.frame(time = time_vec)
    n_t <- length(time_vec)
    
    for(var in predictors) {
      val <- input_vals[[var]]
      if(is.null(val)) return(NULL)
      
      # Get the column from long_data to check type and levels
      orig_col <- long_data[[var]]
      
      if(is.numeric(orig_col)) {
        # Numeric variable: replicate scalar value
        pred_df[[var]] <- rep(as.numeric(val), n_t)
      } else if(is.factor(orig_col)) {
        # Factor: use EXACT levels from the fitted model's data
        lvls <- levels(orig_col)
        # Validate that val is in levels
        if(!(val %in% lvls)) {
          warning(paste("Value", val, "not in factor levels for", var))
          val <- lvls[1]  # Fall back to first level
        }
        pred_df[[var]] <- factor(rep(val, n_t), levels = lvls)
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
        pred_df[[var]] <- factor(rep(val, n_t), levels = lvls)
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
      list(
        fit = as.numeric(pred),
        se = rep(0, length(pred)),
        success = TRUE
      )
    })
    return(result)
  }
