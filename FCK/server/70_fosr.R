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
        f <- as.formula(paste("~", paste(input$reg_predictors, collapse = " + ")))
        X <- model.matrix(f, data = df_reg)
        
        # Point Estimation
        xtx_inv <- solve(crossprod(X))
        projector <- xtx_inv %*% t(X) 
        coefs <- projector %*% Y
        fitted_vals <- X %*% coefs
        residuals <- Y - fitted_vals
        n <- nrow(Y); p <- ncol(X)
        sigma2 <- colSums(residuals^2) / (n - p)
        
        se_mat <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
        p_values <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
        boot_ci_lower <- NULL
        boot_ci_upper <- NULL
        
        if(input$use_bootstrap) {
          showNotification("Bootstrapping...", type = "message")
          B <- input$n_boot
          boot_betas <- array(NA, dim = c(B, nrow(coefs), ncol(coefs)))
          
          withProgress(message = 'Running Bootstrap...', value = 0, {
            for(b in 1:B) {
              # Residual bootstrap: resample residuals, add to fitted
              resid_idx <- sample(1:n, n, replace = TRUE)
              Y_boot <- fitted_vals + residuals[resid_idx, ]
              boot_betas[b, , ] <- projector %*% Y_boot
              if(b %% 20 == 0) incProgress(20/B)
            }
          })
          
          # Calculate SE from bootstrap distribution
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) se_mat[j, k] <- sd(boot_betas[, j, k])
          }
          
          # Percentile-based 95% CI
          boot_ci_lower <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          boot_ci_upper <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) {
              boot_ci_lower[j, k] <- quantile(boot_betas[, j, k], 0.025)
              boot_ci_upper[j, k] <- quantile(boot_betas[, j, k], 0.975)
            }
          }
          
          # Bootstrap p-values: proportion of bootstrap samples on opposite side of 0
          # This is a proper bootstrap test
          p_values <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) {
              if(coefs[j, k] >= 0) {
                p_values[j, k] <- 2 * mean(boot_betas[, j, k] <= 0)
              } else {
                p_values[j, k] <- 2 * mean(boot_betas[, j, k] >= 0)
              }
              # Ensure p-value is in [0, 1]
              p_values[j, k] <- min(1, max(0, p_values[j, k]))
            }
          }
          
        } else {
          # Parametric SE
          for(j in 1:nrow(coefs)) se_mat[j, ] <- sqrt(xtx_inv[j,j] * sigma2)
          t_stats <- coefs / se_mat
          p_values <- 2 * (1 - pt(abs(t_stats), df = n - p))
        }
        
        rss <- colSums(residuals^2); y_bar <- colMeans(Y)
        tss <- colSums(sweep(Y, 2, y_bar)^2); r2_t <- 1 - (rss/tss)
        
        fit <- list(beta.hat = coefs, fitted.values = fitted_vals, resid = residuals, 
                    beta.se = se_mat, beta.p = p_values, terms = terms(f), r2_t = r2_t, 
                    method = ifelse(input$use_bootstrap, "OLS (Bootstrap SE)", "OLS (Parametric SE)"), 
                    xtx_inv = xtx_inv, sigma2 = sigma2,
                    boot_ci_lower = boot_ci_lower, boot_ci_upper = boot_ci_upper,
                    n_boot = if(input$use_bootstrap) input$n_boot else NULL)
        
      } else {
        # --- METHOD B: Smoothed OLS (GAM) ---
        showNotification("Fitting GAM with Splines...", type = "message")
        
        # IMPORTANT: Store original factor levels BEFORE any subsetting
        # This is needed for prediction to work correctly
        orig_factor_levels <- list()
        for(v in input$reg_predictors) {
          if(is.factor(df_reg[[v]])) {
            orig_factor_levels[[v]] <- levels(df_reg[[v]])
          } else if(!is.numeric(df_reg[[v]])) {
            orig_factor_levels[[v]] <- sort(unique(as.character(df_reg[[v]])))
          }
        }
        
        # Reshape to Long format
        df_reg$id_temp <- 1:nrow(df_reg)
        long_cov <- df_reg[rep(seq_len(nrow(df_reg)), each = n_time), ]
        long_data <- long_cov
        long_data$time <- rep(time_points, times = nrow(df_reg))
        long_data$Y_val <- as.vector(t(Y))
        
        # Build Formula: Y ~ s(time) + s(time, by=cov) ...
        gam_formula_str <- "Y_val ~ s(time, bs = 'ps', k = 10)"
        preds <- input$reg_predictors
        
        for(p_var in preds) {
          if(is.numeric(df_reg[[p_var]])) {
            gam_formula_str <- paste0(gam_formula_str, " + s(time, by = ", p_var, ", bs='ps', k=10)")
          } else {
            # Factor interaction
            gam_formula_str <- paste0(gam_formula_str, " + ", p_var, " + s(time, by = ", p_var, ", bs='ps', k=10)")
          }
        }
        
        gam_fit <- gam(as.formula(gam_formula_str), data = long_data, method = "REML")
        
        # Reconstruct Beta(t) for Visualization (Approximation)
        pred_grid <- data.frame(time = time_points)
        beta_hat <- matrix(0, nrow = length(preds) + 1, ncol = n_time)
        rownames(beta_hat) <- c("(Intercept)", preds)
        
        # Intercept approx
        d_int <- pred_grid
        for(v in preds) {
          if(is.numeric(df_reg[[v]])) d_int[[v]] <- 0
          else {
            # Use proper factor levels from long_data (which is what GAM was trained on)
            orig_col <- long_data[[v]]
            if(is.factor(orig_col)) {
              d_int[[v]] <- factor(rep(levels(orig_col)[1], n_time), levels = levels(orig_col))
            } else {
              lvls <- sort(unique(as.character(orig_col)))
              d_int[[v]] <- factor(rep(lvls[1], n_time), levels = lvls)
            }
          }
        }
        beta_hat[1, ] <- predict(gam_fit, newdata = d_int, type = "response")
        
        # Covariate effects approx
        for(i in 1:length(preds)) {
          v <- preds[i]
          if(is.numeric(df_reg[[v]])) {
            d_0 <- d_int; d_0[[v]] <- 0
            d_1 <- d_int; d_1[[v]] <- 1
            beta_hat[i+1, ] <- predict(gam_fit, newdata = d_1) - predict(gam_fit, newdata = d_0)
          }
        }
        
        fitted_vec <- predict(gam_fit, newdata = long_data)
        fitted_vals <- matrix(fitted_vec, nrow = nrow(Y), ncol = n_time, byrow = TRUE)
        residuals <- Y - fitted_vals
        
        rss <- colSums(residuals^2); y_bar <- colMeans(Y)
        tss <- colSums(sweep(Y, 2, y_bar)^2); r2_t <- 1 - (rss/tss)
        
        fit <- list(beta.hat = beta_hat, fitted.values = fitted_vals, resid = residuals, 
                    beta.se = beta_hat*0, beta.p = beta_hat*0, # placeholders
                    terms = terms(as.formula(paste("~", paste(preds, collapse="+")))), 
                    r2_t = r2_t, method = "Smoothed OLS (GAM)",
                    gam_obj = gam_fit,
                    gam_predictors = preds,
                    gam_long_data = long_data,  # Store long_data for exact factor levels
                    gam_factor_levels = orig_factor_levels,  # Store original factor levels
                    gam_n_time = n_time)
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
    if(!is.null(mod$n_boot)) {
      cat("Bootstrap: Yes (B =", mod$n_boot, "), Percentile CIs\n")
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
    subset_idx <- 1:min(nrow(Y), 200) 
    Y_subset <- Y[subset_idx, , drop=FALSE]
    col_vals_raw <- values$covariates[[input$reg_color_var]][subset_idx]
    
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
      labs(title = paste("Observed Data colored by", input$reg_color_var))
    
    if(is_categorical) g <- g + scale_color_brewer(palette = "Set1") else g <- g + scale_color_viridis_c(option = "viridis")
    
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
      pred_df <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals)
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
      pred_df <- as.data.frame(pred_data)
      
      betas <- mod$beta.hat
      f_clean <- delete.response(mod$terms)
      X_new <- tryCatch({ model.matrix(f_clean, data = pred_df) }, error = function(e) NULL)
      if(is.null(X_new)) return(NULL)
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
          pred_df_min <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals_min)
          pred_df_max <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals_max)
          
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
          
          df_min <- as.data.frame(pred_data_min)
          df_max <- as.data.frame(pred_data_max)
          
          f_clean <- delete.response(mod$terms)
          X_min <- tryCatch({ model.matrix(f_clean, data = df_min) }, error = function(e) NULL)
          X_max <- tryCatch({ model.matrix(f_clean, data = df_max) }, error = function(e) NULL)
          
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
    } else if(sum(mod$beta.se[idx, ]) > 0) {
      se <- mod$beta.se[idx, ]
      p <- p %>% add_ribbons(x = t, ymin = beta - 1.96*se, ymax = beta + 1.96*se,
                             name = "95% Parametric CI", line = list(color = 'transparent'), 
                             fillcolor = 'rgba(0, 191, 255, 0.3)')
    }
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    
    ci_type <- if(!is.null(mod$boot_ci_lower)) "Bootstrap" else "Parametric"
    p %>% layout(title = paste("Coefficient:", sel, "(", ci_type, "CI)"))
  })
  
  output$reg_pvalue_plot <- renderPlotly({
    req(values$reg_model, input$reg_coeff_select)
    mod <- values$reg_model
    sel <- input$reg_coeff_select
    idx <- which(rownames(mod$beta.hat) == sel)
    pvals <- mod$beta.p[idx, ]
    if(all(pvals == 0) || all(is.na(pvals))) return(NULL)
    
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
