# ==========================================================================
# server/70_fosr.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 1658-2242  (Function-on-Scalar regression)
# ==========================================================================
  # ==============================================================================
  # FoSR LOGIC (Function-on-Scalar)
  # ==============================================================================
  
  output$fosr_var_select_ui <- renderUI({
    req(values$covariates)
    tagList(
      selectInput("reg_predictors", "Select Scalar Predictors:", choices = colnames(values$covariates), multiple = TRUE),
      helpText("Response is the FUNCTIONAL data.")
    )
  })
  
  observeEvent(input$run_fosr, {
    req(values$data, values$covariates, input$reg_predictors)
    showNotification("Running FoSR...", type = "message", duration = 2)
    tryCatch({
      Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      mode(Y) <- "numeric"
      df_reg <- values$covariates
      
      # Subset Complete Cases
      keep_idx <- complete.cases(df_reg[, input$reg_predictors, drop=FALSE]) & complete.cases(Y)
      if(sum(keep_idx) < nrow(df_reg)) {
        df_reg <- df_reg[keep_idx, , drop=FALSE]
        Y <- Y[keep_idx, , drop=FALSE]
      }
      
      n_time <- ncol(Y)
      time_points <- seq(0, 1, length.out = n_time)
      
      # --- METHOD A: Pointwise OLS ---
      if (input$reg_method == "OLS_nosmooth") {
        # P3.4: the estimator itself lives in server/07_helpers_fosr.R as a
        # named function, so the Export tab can deparse() the SAME object into
        # the generated script instead of maintaining a second, drifting
        # transcription of it. See the note at the top of that file.
        fit <- withProgress(
          message = "Fitting pointwise OLS...", value = 0,
          fck_fit_fosr_ols(
            Y             = Y,
            df_reg        = df_reg,
            predictors    = input$reg_predictors,
            use_bootstrap = isTRUE(input$use_bootstrap),
            n_boot        = if(isTRUE(input$use_bootstrap)) input$n_boot else 0,
            progress      = function(frac, detail = NULL) incProgress(frac, detail = detail),
            notify        = function(msg, ...) showNotification(msg, type = "message")))

      } else {
        # --- METHOD B: Smoothed OLS (GAM) ---
        # P5.8: the estimator lives in server/07b_helpers_fosr_gam.R as a named
        # function, so the export can deparse the SAME object instead of
        # maintaining a second, drifting transcription of it.
        fit <- fck_fit_fosr_gam(
          Y          = Y,
          df_reg     = df_reg,
          predictors = input$reg_predictors,
          k          = 10,
          notify     = function(msg, ...) showNotification(msg, type = "message"))

      }
      
      # MERGED APP: what this fit actually used, for the code export.
      fit$fck_settings <- list(
        predictors   = input$reg_predictors,
        method       = input$reg_method,
        use_bootstrap = isTRUE(input$use_bootstrap),
        n_boot       = if(isTRUE(input$use_bootstrap)) input$n_boot else NULL,
        using_smoothed = !is.null(values$smooth_data),
        n_subjects   = nrow(Y),
        n_time       = ncol(Y))

      values$reg_model <- fit
      updateSelectInput(session, "reg_color_var", choices = input$reg_predictors)
      updateSelectInput(session, "reg_coeff_select", choices = rownames(fit$beta.hat))
      num_preds <- input$reg_predictors[sapply(values$covariates[input$reg_predictors], is.numeric)]
      updateSelectInput(session, "reg_ref_selector", choices = num_preds)
      
      showNotification("FoSR Fitted!", type = "message")
    }, error = function(e) showNotification(paste("Fit Error:", e$message), type="error"))
  })
  
  output$fosr_model_summary <- renderPrint({
    req(values$reg_model)
    mod <- values$reg_model
    cat("Method:", mod$method, "\n")
    cat("Coefficients:", paste(rownames(mod$beta.hat), collapse=", "), "\n")
    if(!is.null(mod$r2_t)) cat("Mean R2:", round(mean(mod$r2_t, na.rm=TRUE), 3), "\n")
    if(identical(mod$inference, "analytic-t-fdr")) {
      cat("\nInference (P3.2): p-values come from the ANALYTICAL OLS standard\n")
      cat(sprintf("error, beta / sqrt(diag((X'X)^-1) * sigma2(t)), referred to t on %s df,\n",
                  if(is.null(mod$df_resid)) "n - p" else format(mod$df_resid)))
      cat("then FDR-adjusted across time points. They do not depend on a seed.\n")
    }
    if(!is.null(mod$n_boot)) {
      cat("\nBootstrap: residual bootstrap, B =", mod$n_boot, "\n")
      cat("It supplies the 95% PERCENTILE INTERVAL only -- it is not the test.\n")
      if(!is.null(mod$se_ratio_range)) {
        cat(sprintf("Diagnostic: SE_bootstrap / SE_analytic ranges %.2f to %.2f.\n",
                    mod$se_ratio_range[1], mod$se_ratio_range[2]))
        if(is.finite(mod$se_ratio_range[2]) &&
           (mod$se_ratio_range[2] > 1.25 || mod$se_ratio_range[1] < 0.8))
          cat("  WARNING: these should agree closely -- the residual bootstrap\n",
              "  resamples from the same homoscedastic model the analytical SE\n",
              "  assumes. A ratio far from 1 means that model is in doubt (or B is\n",
              "  too small). Treat both the interval and the p-values with caution.\n", sep = "")
      }
      cat("The interval and the test can therefore disagree at the margin; the\n")
      cat("interval is the percentile one, not an inversion of the t test.\n")
    }
    if(!is.null(mod$inference_note)) {
      cat("\n"); cat(strwrap(mod$inference_note, width = 76), sep = "\n"); cat("\n")
    }
    if(!is.null(mod$gam_obj)) {
      cat("\nFamily:", mod$gam_obj$family$family, "\n")
      cat("Link function:", mod$gam_obj$family$link, "\n\n")
      print(summary(mod$gam_obj))
    }
  })
  
  # FoSR Visualization
  output$reg_observed_plot <- renderPlotly({
    req(values$data, values$covariates, input$reg_color_var)
    Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
    col_all <- values$covariates[[input$reg_color_var]]

    # AUDIT (P6.5): this took `1:min(nrow(Y), 200)` -- the FIRST 200 rows in
    # file order. Data exported per group arrives sorted by group, so on a
    # sample of 654 YOUTH, 410 ADULT, 181 MIDDLE_AGE and 59 ELDERLY the plot
    # drew 200 YOUTH curves and one legend entry, and looked exactly like a
    # dataset with one category in it. The cap is there to keep the plot
    # drawable; it must not decide WHICH curves you see.
    #
    # Stratified instead: every level is represented, in proportion, with at
    # least a few curves each, and the caption says what was drawn.
    max_curves <- 300L
    n_all <- nrow(Y)
    if (n_all <= max_curves) {
      subset_idx <- seq_len(n_all)
      sample_note <- sprintf("all %d curves", n_all)
    } else {
      f <- as.factor(col_all)
      lv <- levels(droplevels(f))
      idx_by <- split(seq_len(n_all), f)[lv]
      # proportional, floor 5 per level, then trim the largest levels back to
      # the cap so the total is respected
      take <- pmax(5L, round(max_curves * lengths(idx_by) / n_all))
      take <- pmin(take, lengths(idx_by))
      while (sum(take) > max_curves) {
        big <- which.max(take); if (take[big] <= 5L) break; take[big] <- take[big] - 1L
      }
      set.seed(1)   # the same picture on every redraw of the same data
      subset_idx <- sort(unlist(Map(function(ii, k)
        if (length(ii) <= k) ii else sample(ii, k), idx_by, take), use.names = FALSE))
      sample_note <- sprintf("a stratified sample of %d of %d curves (%s)",
                             length(subset_idx), n_all,
                             paste(sprintf("%s %d", lv, take), collapse = ", "))
    }
    Y_subset <- Y[subset_idx, , drop=FALSE]
    col_vals_raw <- col_all[subset_idx]
    
    df_plot <- data.frame(
      Time = rep(seq(0, 1, length.out = ncol(Y)), times = nrow(Y_subset)),
      Value = as.vector(t(Y_subset)),
      ID = as.factor(rep(subset_idx, each = ncol(Y))),
      Color = rep(col_vals_raw, each = ncol(Y))
    )
    
    is_categorical <- length(unique(col_vals_raw)) <= 5 || is.factor(col_vals_raw) || is.character(col_vals_raw)
    if(is_categorical) df_plot$Color <- as.factor(df_plot$Color)
    
    g <- ggplot(df_plot, aes(x = Time, y = Value, group = ID, color = Color)) +
      geom_line(alpha = 0.8, linewidth = 0.5) + theme_minimal() +
      labs(title = paste("Observed Data colored by", input$reg_color_var),
           subtitle = paste("Showing", sample_note))

    # P6.5: the app's own validated palette, so a group is the same colour here
    # as everywhere else. scale_color_brewer("Set1") was a fourth palette in an
    # app that already has one.
    if(is_categorical) {
      g <- g + scale_color_manual(values = fck_group_ramp(nlevels(df_plot$Color)))
    } else {
      g <- g + scale_color_viridis_c(option = "viridis")
    }
    
    if(!is.null(values$time_labels)) {
      tick_idx <- seq(0, 1, length.out = length(values$time_labels))
      g <- g + scale_x_continuous(breaks = tick_idx, labels = values$time_labels) +
        theme(axis.text.x = element_text(angle = -90, hjust = 0))
    }
    ggplotly(g, tooltip = c("group", "x", "y", "colour"))
  })
  
  output$reg_ref_selector <- renderUI({
    req(values$reg_model)
    num_preds <- input$reg_predictors[sapply(values$covariates[input$reg_predictors], is.numeric)]
    selectInput("ref_predictor", "Select Predictor for Min/Max Reference:", choices = num_preds)
  })
  
  output$reg_prediction_controls <- renderUI({
    req(values$reg_model, input$reg_predictors)
    lapply(input$reg_predictors, function(var) {
      vals <- values$covariates[[var]]
      if(is.numeric(vals)) {
        sliderInput(paste0("pred_", var), label = var, 
                    min = min(vals, na.rm=TRUE), max = max(vals, na.rm=TRUE), 
                    value = mean(vals, na.rm=TRUE))
      } else {
        # For factors, use the factor levels if available
        if(is.factor(vals)) {
          choices <- levels(vals)
        } else {
          choices <- sort(unique(as.character(vals)))
        }
        selectInput(paste0("pred_", var), label = var, choices = choices, selected = choices[1])
      }
    })
  })
  
  # =============================================================================
  # COMPLETELY REWRITTEN: reg_fitted_plot - Prediction plot for GAM/OLS
  # =============================================================================
  output$reg_fitted_plot <- renderPlotly({
    req(values$reg_model, input$reg_predictors)
    
    mod <- values$reg_model
    n_t <- if(!is.null(mod$gam_n_time)) mod$gam_n_time else {
      if(!is.null(mod$fitted.values)) ncol(mod$fitted.values) else ncol(values$data)
    }
    time_vec <- seq(0, 1, length.out = n_t)
    
    # Collect input values
    input_vals <- list()
    for(var in input$reg_predictors) {
      val <- input[[paste0("pred_", var)]]
      if(is.null(val)) return(NULL)
      input_vals[[var]] <- val
    }
    
    # ==========================================================================
    # GAM Method
    # ==========================================================================
    if(!is.null(mod$gam_obj)) {
      preds <- mod$gam_predictors
      long_data <- mod$gam_long_data
      factor_levels <- mod$gam_factor_levels
      
      # Build prediction dataframe using helper function
      # P8.1: the fitted GAM's formula uses mod$gam_model_names, not the
      # user's column names. The mapping is a required argument now.
      pred_df <- build_gam_pred_df(time_vec, preds, mod$gam_model_names,
                                   long_data, factor_levels, input_vals)
      if(is.null(pred_df)) return(NULL)
      
      # Get prediction with SE using safe helper
      pred_result <- safe_gam_predict(mod$gam_obj, pred_df)
      if(!pred_result$success || is.null(pred_result$fit)) {
        showNotification("Prediction failed", type = "error")
        return(NULL)
      }
      
      y_hat <- pred_result$fit
      se_pred <- pred_result$se
      
      # Ensure vectors are the correct length
      if(length(y_hat) != n_t) {
        showNotification(paste("Prediction length mismatch:", length(y_hat), "vs", n_t), type = "error")
        return(NULL)
      }
      if(length(se_pred) != n_t) {
        se_pred <- rep(0, n_t)  # Fallback to no CI if SE has wrong length
      }
      
    } else {
      # ==========================================================================
      # OLS Method
      # ==========================================================================
      pred_data <- list()
      for(var in input$reg_predictors) {
        val <- input_vals[[var]]
        orig <- values$covariates[[var]]
        if(is.numeric(orig)) {
          pred_data[[var]] <- as.numeric(val)
        } else {
          if(is.factor(orig)) {
            pred_data[[var]] <- factor(val, levels = levels(orig))
          } else {
            pred_data[[var]] <- factor(val, levels = sort(unique(as.character(orig))))
          }
        }
      }
      pred_df <- as.data.frame(pred_data, check.names = FALSE)

      # AUDIT (P5.1): this used to call
      #     model.matrix(delete.response(mod$terms), data = pred_df)
      # directly. After P4.3 the model's terms refer to the internal names
      # x1..xp while pred_df is keyed by the user's column names, so
      # model.matrix could not find its variables and the tryCatch below
      # returned NULL -- a blank prediction curve, on EVERY FoSR OLS fit, with
      # no error shown. fck_fosr_design() does the rename and re-applies the
      # fitted factor levels; it is the same code the fit itself uses.
      betas <- mod$beta.hat
      X_new <- tryCatch(fck_fosr_design(mod, pred_df),
                        error = function(e) structure(conditionMessage(e), class = "fck_err"))
      if(inherits(X_new, "fck_err")) {
        # and it is REPORTED now. The old code swallowed the error and returned
        # NULL, which draws nothing: a broken prediction looked like an empty
        # plot, which is the failure mode that let P5.1 survive a release.
        showNotification(paste("Could not build the prediction design:",
                               as.character(X_new)), type = "error", duration = 12)
        return(NULL)
      }
      y_hat <- as.vector(X_new %*% betas)
      
      if(!is.null(mod$xtx_inv) && !is.null(mod$sigma2)) {
        scale_factor <- as.numeric(X_new %*% mod$xtx_inv %*% t(X_new))
        se_pred <- sqrt(scale_factor * mod$sigma2)
      } else {
        se_pred <- rep(0, length(y_hat))
      }
    }
    
    # ==========================================================================
    # Build main plot
    # ==========================================================================
    # Compute CI bounds
    ci_lower <- y_hat - 1.96 * se_pred
    ci_upper <- y_hat + 1.96 * se_pred
    
    # Create plot data frame for consistent handling
    plot_df <- data.frame(
      time = time_vec,
      y_hat = y_hat,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
    
    p <- plot_ly(data = plot_df, x = ~time) %>%
      add_trace(y = ~y_hat, type = 'scatter', mode = 'lines',
                line = list(color = 'red', width = 3), name = "Predicted Mean")
    
    # Only add ribbon if SE is non-zero
    if(any(se_pred > 0)) {
      p <- p %>% add_ribbons(ymin = ~ci_lower, ymax = ~ci_upper,
                             name = "95% CI", line = list(color = 'transparent'), 
                             fillcolor = 'rgba(255, 0, 0, 0.2)')
    }
    
    # ==========================================================================
    # Reference Lines
    # ==========================================================================
    ref_pred <- input$ref_predictor
    if(!is.null(ref_pred) && ref_pred %in% input$reg_predictors) {
      orig_vals <- values$covariates[[ref_pred]]
      
      if(is.numeric(orig_vals)) {
        min_val <- min(orig_vals, na.rm = TRUE)
        max_val <- max(orig_vals, na.rm = TRUE)
        
        # Create modified input values for min/max
        input_vals_min <- input_vals
        input_vals_max <- input_vals
        input_vals_min[[ref_pred]] <- min_val
        input_vals_max[[ref_pred]] <- max_val
        
        if(!is.null(mod$gam_obj)) {
          # GAM reference lines
          # P8.1: same mapping for the reference curves.
          pred_df_min <- build_gam_pred_df(time_vec, preds, mod$gam_model_names,
                                           long_data, factor_levels, input_vals_min)
          pred_df_max <- build_gam_pred_df(time_vec, preds, mod$gam_model_names,
                                           long_data, factor_levels, input_vals_max)
          
          if(!is.null(pred_df_min) && !is.null(pred_df_max)) {
            pred_min <- safe_gam_predict(mod$gam_obj, pred_df_min)
            pred_max <- safe_gam_predict(mod$gam_obj, pred_df_max)
            
            if(pred_min$success && pred_max$success) {
              p <- p %>% 
                add_lines(x = time_vec, y = pred_min$fit, name = paste("Min", ref_pred), 
                          line = list(color = 'blue', dash = 'dot', width = 2)) %>%
                add_lines(x = time_vec, y = pred_max$fit, name = paste("Max", ref_pred), 
                          line = list(color = 'purple', dash = 'dot', width = 2))
            }
          }
        } else {
          # OLS reference lines
          pred_data_min <- pred_data_max <- list()
          for(var in input$reg_predictors) {
            orig <- values$covariates[[var]]
            if(var == ref_pred) {
              pred_data_min[[var]] <- min_val
              pred_data_max[[var]] <- max_val
            } else {
              val <- input_vals[[var]]
              if(is.numeric(orig)) {
                pred_data_min[[var]] <- as.numeric(val)
                pred_data_max[[var]] <- as.numeric(val)
              } else {
                if(is.factor(orig)) {
                  pred_data_min[[var]] <- factor(val, levels = levels(orig))
                  pred_data_max[[var]] <- factor(val, levels = levels(orig))
                } else {
                  lvls <- sort(unique(as.character(orig)))
                  pred_data_min[[var]] <- factor(val, levels = lvls)
                  pred_data_max[[var]] <- factor(val, levels = lvls)
                }
              }
            }
          }
          
          df_min <- as.data.frame(pred_data_min, check.names = FALSE)
          df_max <- as.data.frame(pred_data_max, check.names = FALSE)

          # P5.1: same shared builder as the main prediction curve.
          X_min <- tryCatch(fck_fosr_design(mod, df_min), error = function(e) NULL)
          X_max <- tryCatch(fck_fosr_design(mod, df_max), error = function(e) NULL)
          
          if(!is.null(X_min) && !is.null(X_max)) {
            y_min <- as.vector(X_min %*% mod$beta.hat)
            y_max <- as.vector(X_max %*% mod$beta.hat)
            
            p <- p %>% 
              add_lines(x = time_vec, y = y_min, name = paste("Min", ref_pred), 
                        line = list(color = 'blue', dash = 'dot', width = 2)) %>%
              add_lines(x = time_vec, y = y_max, name = paste("Max", ref_pred), 
                        line = list(color = 'purple', dash = 'dot', width = 2))
          }
        }
      }
    }
    
    # Add time labels if available
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, 
                                     ticktext = values$time_labels, tickangle = -90))
    }
    
    p %>% layout(
      title = paste("Predicted Curve (Y ~", paste(input$reg_predictors, collapse = " + "), ")"),
      xaxis = list(title = "Time (normalized)"),
      yaxis = list(title = "Predicted Value")
    )
  })
  
  output$reg_coeff_plot <- renderPlotly({
    req(values$reg_model, input$reg_coeff_select)
    mod <- values$reg_model
    sel <- input$reg_coeff_select
    idx <- which(rownames(mod$beta.hat) == sel)
    if(length(idx) == 0) return(NULL)
    
    beta <- mod$beta.hat[idx, ]
    t <- seq(0, 1, length.out = length(beta))
    
    p <- plot_ly(x = t, y = beta, type = 'scatter', mode = 'lines', name = 'Beta(t)',
                 line = list(color = '#003366', width = 2)) %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 0, 
                   line = list(color = 'black', dash = 'dash', width = 1), showlegend = FALSE)
    
    # Use bootstrap percentile CIs if available, otherwise parametric
    if(!is.null(mod$boot_ci_lower) && !is.null(mod$boot_ci_upper)) {
      ci_lower <- mod$boot_ci_lower[idx, ]
      ci_upper <- mod$boot_ci_upper[idx, ]
      ci_name <- paste0("95% Bootstrap CI (B=", mod$n_boot, ")")
      p <- p %>% add_ribbons(x = t, ymin = ci_lower, ymax = ci_upper,
                             name = ci_name, line = list(color = 'transparent'), 
                             fillcolor = 'rgba(0, 191, 255, 0.3)')
    } else if(any(is.finite(mod$beta.se[idx, ]) & mod$beta.se[idx, ] > 0)) {
      # P5.9: was sum(beta.se[idx, ]) > 0, which is NA when the SE row is all
      # NA -- and `if (NA)` is an error, not a skip. The GAM branch now stores
      # NA rather than zeros, so this has to test for a usable SE explicitly.
      se <- mod$beta.se[idx, ]
      p <- p %>% add_ribbons(x = t, ymin = beta - 1.96*se, ymax = beta + 1.96*se,
                             name = "95% Parametric CI", line = list(color = 'transparent'), 
                             fillcolor = 'rgba(0, 191, 255, 0.3)')
    }
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    
    # P3.2: name the interval AND say what the test is, because they are now
    # deliberately different procedures and can disagree at the margin.
    ci_type <- if(!is.null(mod$boot_ci_lower))
      paste0("95% bootstrap percentile CI, B=", mod$n_boot) else "95% analytical CI"
    p %>% layout(title = list(text = paste0(
      "Coefficient: ", sel,
      "<br><sub>", ci_type,
      if(identical(mod$inference, "analytic-t-fdr"))
        "; p-values from the analytical SE, FDR-adjusted" else "",
      "</sub>")))
  })
  
  output$reg_pvalue_plot <- renderPlotly({
    req(values$reg_model, input$reg_coeff_select)
    mod <- values$reg_model
    sel <- input$reg_coeff_select
    idx <- which(rownames(mod$beta.hat) == sel)
    pvals <- mod$beta.p[idx, ]
    # P5.9: say why the panel is empty instead of drawing nothing. A blank plot
    # is indistinguishable from a broken one -- that is how the P5.1 prediction
    # regression survived a release.
    if(all(is.na(pvals)) || all(pvals == 0, na.rm = TRUE)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>%
        layout(title = list(text = paste0(
          "No coefficient-curve p-values for this fit",
          "<br><sub>",
          if (!is.null(mod$inference_note))
            "GAM branch: the curves are prediction contrasts, and no standard error is propagated through them. Use summary(gam) below, or the pointwise-OLS method."
          else "This model reports no pointwise inference.",
          "</sub>"))))
    }
    
    t <- seq(0, 1, length.out = length(pvals))
    p <- plot_ly(x = t, y = pvals, type = 'scatter', mode = 'lines', name = 'P-value',
                 line = list(color = 'darkgreen', width = 2)) %>%
      add_lines(y = 0.05, line = list(color = 'red', dash = 'dash'), name = "p=0.05") %>%
      layout(yaxis = list(type = "log", title = "P-value (log scale)"), title = "Significance over Time")
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p
  })
  
  output$reg_r2_plot <- renderPlotly({
    req(values$reg_model)
    mod <- values$reg_model
    t <- seq(0, 1, length.out = length(mod$r2_t))
    mean_r2 <- mean(mod$r2_t, na.rm=TRUE)
    
    p <- plot_ly(x = t, y = mod$r2_t, type = 'scatter', mode = 'lines', 
                 fill = 'tozeroy', name = 'R-squared',
                 line = list(color = 'purple', width = 2)) %>%
      add_lines(y = mean_r2, line = list(color = 'black', dash = 'dash'), name = paste("Mean R2:", round(mean_r2, 2))) %>%
      layout(title = "Model Fit (Functional R-squared)", yaxis = list(range = c(0, 1), title = "R-squared"))
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p
  })
  
  output$reg_residual_plot <- renderPlotly({
    req(values$reg_model)
    resid <- values$reg_model$resid
    t <- seq(0, 1, length.out = ncol(resid)); n <- nrow(resid)
    l2 <- rowSums(resid^2); rks <- rank(l2)
    cols <- colorRampPalette(c("red", "orange", "green", "blue", "purple"))(n)[rks]
    
    idx <- if(n > 200) sample(1:n, 200) else 1:n
    p <- plot_ly(type = 'scatter', mode = 'lines', hoverinfo = "text")
    for(i in idx) {
      p <- add_trace(p, x = t, y = resid[i,], line = list(color = cols[i], width = 1), showlegend = FALSE)
    }
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p %>% layout(title = "Residuals (Rainbow Plot)")
  })
