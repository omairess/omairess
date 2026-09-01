# ==========================================================================
# server/71_sofr.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 2244-2877  (Scalar-on-Function regression)
# ==========================================================================
  # ==============================================================================
  # SoFR LOGIC (Scalar-on-Function)
  # ==============================================================================
  
  # Store detected family for conditional panels
  output$sofr_is_binary <- reactive({
    if(is.null(values$sofr_model)) return(FALSE)
    !is.null(values$sofr_model$family) && values$sofr_model$family$family == "binomial"
  })
  outputOptions(output, "sofr_is_binary", suspendWhenHidden = FALSE)
  
  output$sofr_var_select_ui <- renderUI({
    req(values$covariates)
    tagList(
      selectInput("sofr_response", "Select Scalar Response Variable (y):", choices = colnames(values$covariates), multiple = FALSE),
      selectInput("sofr_scalar_preds", "Additional Scalar Predictors (z):", choices = colnames(values$covariates), multiple = TRUE),
      helpText("Functional Predictor (X) is the imported curve data.")
    )
  })
  
  # Dynamic family info based on response variable
  output$sofr_family_info <- renderUI({
    req(input$sofr_response, values$covariates)
    y <- values$covariates[[input$sofr_response]]
    
    if(is.null(y)) return(NULL)
    
    # Handle factors
    if(is.factor(y)) {
      n_levels <- nlevels(y)
      if(n_levels == 2) {
        return(div(class = "alert alert-info", style = "padding: 8px; margin-top: 10px;", 
                   icon("info-circle"), " ",
                   paste0("Detected: Binary factor (", n_levels, " levels: ", 
                          paste(levels(y), collapse = ", "), "). ",
                          "Will convert to 0/1. Recommended: Binomial family with logit link.")))
      } else {
        return(div(class = "alert alert-warning", style = "padding: 8px; margin-top: 10px;", 
                   icon("exclamation-triangle"), " ",
                   paste0("Factor with ", n_levels, " levels detected. ",
                          "Only 2-level factors (binary) are supported for SoFR.")))
      }
    }
    
    if(!is.numeric(y)) return(NULL)
    
    unique_vals <- length(unique(na.omit(y)))
    is_binary <- unique_vals == 2
    is_count <- all(y == floor(y), na.rm = TRUE) && all(y >= 0, na.rm = TRUE)
    is_proportion <- all(y >= 0 & y <= 1, na.rm = TRUE)
    
    msg <- ""
    if(is_binary) {
      msg <- paste0("Detected: Binary variable (", unique_vals, " unique values). ",
                    "Recommended: Binomial family with logit link.")
    } else if(is_proportion && !is_binary) {
      msg <- "Detected: Values in [0,1]. Could be proportions (use Binomial) or continuous (use Gaussian)."
    } else if(is_count && min(y, na.rm=TRUE) >= 0) {
      msg <- "Detected: Non-negative integers. Consider Poisson for count data."
    } else if(all(y > 0, na.rm = TRUE)) {
      msg <- "Detected: Positive continuous. Gaussian or Gamma may be appropriate."
    } else {
      msg <- "Detected: Continuous variable. Gaussian family recommended."
    }
    
    div(class = "alert alert-info", style = "padding: 8px; margin-top: 10px;", 
        icon("info-circle"), " ", msg)
  })
  
  observeEvent(input$run_sofr, {
    req(values$data, values$covariates, input$sofr_response)
    showNotification("Running SoFR...", type = "message", duration = 2)
    tryCatch({
      X_func <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      df_sofr <- values$covariates
      y <- df_sofr[[input$sofr_response]]
      
      # Handle factor response for binary
      if(is.factor(y)) {
        y <- as.numeric(y) - 1  # Convert to 0/1
        showNotification("Factor response converted to 0/1 numeric.", type = "message")
      }
      
      if(!is.numeric(y)) stop("Response variable must be numeric (or factor for binary).")
      
      preds <- input$sofr_scalar_preds
      if(is.null(preds)) preds <- character(0)  # Ensure it's at least empty character vector
      
      # Handle complete cases - be careful with empty preds
      if(length(preds) > 0) {
        keep_idx <- complete.cases(df_sofr[, c(input$sofr_response, preds), drop=FALSE]) & complete.cases(X_func)
      } else {
        keep_idx <- complete.cases(df_sofr[, input$sofr_response, drop=FALSE]) & complete.cases(X_func)
      }
      
      y_clean <- y[keep_idx]
      X_func_clean <- X_func[keep_idx, ]
      df_clean <- df_sofr[keep_idx, , drop=FALSE]
      
      # Get link function values with defaults (inputs may be NULL if conditionalPanel not shown)
      link_binom <- if(!is.null(input$sofr_link_binomial)) input$sofr_link_binomial else "logit"
      link_gauss <- if(!is.null(input$sofr_link_gaussian)) input$sofr_link_gaussian else "identity"
      family_choice <- if(!is.null(input$sofr_family)) input$sofr_family else "auto"
      
      # Determine the family and validate/convert y accordingly
      is_binary_outcome <- FALSE
      unique_y <- sort(unique(na.omit(y_clean)))
      
      if(family_choice == "binomial" || (family_choice == "auto" && length(unique_y) == 2)) {
        # For binomial: ensure y is strictly 0/1
        if(length(unique_y) == 2) {
          # Convert to 0/1 if needed (e.g., if values are 1/2 or other pairs)
          y_clean <- as.numeric(y_clean == max(unique_y))
          showNotification(paste("Binary response recoded: ", min(unique_y), "->0, ", max(unique_y), "->1"), type = "message")
        } else if(all(y_clean >= 0 & y_clean <= 1)) {
          # Already proportions - use as is
          showNotification("Using response as proportions (0-1 range).", type = "message")
        } else {
          stop("For binomial family, response must be binary (2 unique values) or proportions (0-1 range).")
        }
        pfr_family <- binomial(link = link_binom)
        is_binary_outcome <- (length(unique_y) == 2)
      } else if(family_choice == "gaussian") {
        pfr_family <- gaussian(link = link_gauss)
      } else if(family_choice == "poisson") {
        if(any(y_clean < 0) || any(y_clean != floor(y_clean))) {
          showNotification("Warning: Poisson expects non-negative integer counts.", type = "warning")
        }
        pfr_family <- poisson(link = "log")
      } else if(family_choice == "Gamma") {
        if(any(y_clean <= 0)) {
          stop("Gamma family requires strictly positive response values.")
        }
        pfr_family <- Gamma(link = "log")
      } else {
        # Default: auto-detect -> gaussian
        pfr_family <- gaussian(link = "identity")
      }
      
      # Prepare data for pfr - IMPORTANT: remove original response column
      # to ensure pfr uses our converted y variable
      df_clean_no_response <- df_clean
      df_clean_no_response[[input$sofr_response]] <- NULL  # Remove original response
      
      pfr_data <- as.list(df_clean_no_response)
      pfr_data$X_func <- X_func_clean
      pfr_data$y <- y_clean  # Use our validated/converted y
      
      # Debug info
      cat("SoFR Debug:\n")
      cat("  Family:", pfr_family$family, "\n")
      cat("  y range:", range(y_clean), "\n")
      cat("  y unique values:", paste(head(sort(unique(y_clean)), 10), collapse=", "), "\n")
      
      formula_str <- "y ~ lf(X_func, bs='ps', k=15)"
      if(!is.null(preds) && length(preds) > 0) {
        formula_str <- paste(formula_str, "+", paste(preds, collapse = " + "))
      }
      pfr_formula <- as.formula(formula_str)
      
      # Final validation for binomial
      if(pfr_family$family == "binomial") {
        if(any(pfr_data$y < 0 | pfr_data$y > 1, na.rm = TRUE)) {
          stop(paste("After conversion, y values are still outside [0,1]. Range:", 
                     paste(range(pfr_data$y, na.rm=TRUE), collapse=" to ")))
        }
        cat("  Binomial validation passed: y in [0,1]\n")
      }
      
      # Use do.call to force evaluation of all arguments
      # This avoids pfr's non-standard evaluation issues
      fit <- do.call(pfr, list(
        formula = pfr_formula, 
        data = pfr_data, 
        family = pfr_family
      ))
      
      # Store additional info for diagnostics
      fit$y_original <- y_clean
      fit$is_binary <- is_binary_outcome

      # MERGED APP: what this fit actually used, for the code export.
      fit$fck_settings <- list(
        response     = input$sofr_response,
        predictors   = preds,
        formula      = formula_str,
        family       = pfr_family$family,
        link         = pfr_family$link,
        using_smoothed = !is.null(values$smooth_data),
        n_obs        = length(y_clean))
      
      # Bootstrap for coefficient CIs if requested
      if(isTRUE(input$sofr_use_bootstrap)) {
        B <- if(!is.null(input$sofr_n_boot)) input$sofr_n_boot else 100
        n_obs <- length(y_clean)
        
        showNotification(paste("Running", B, "bootstrap iterations..."), type = "message")
        
        # Extract coefficient from fitted model to get dimensions
        # The functional coefficient is in the smooth term
        coef_orig <- coef(fit)
        
        # Get the functional coefficient values over a grid
        # Use the smooth object to evaluate
        n_grid <- 100
        t_grid <- seq(0, 1, length.out = n_grid)
        
        # Store bootstrap coefficients
        boot_coefs <- matrix(NA, nrow = B, ncol = n_grid)
        boot_success <- 0
        
        withProgress(message = 'Running SoFR Bootstrap...', value = 0, {
          for(b in 1:B) {
            tryCatch({
              # Case resampling
              boot_idx <- sample(1:n_obs, n_obs, replace = TRUE)
              
              # Create bootstrap data
              pfr_data_boot <- list()
              for(nm in names(pfr_data)) {
                if(nm == "X_func") {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx, , drop = FALSE]
                } else if(nm == "y") {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx]
                } else if(length(pfr_data[[nm]]) == n_obs) {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx]
                } else {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]]
                }
              }
              
              # Fit bootstrap model (suppress warnings)
              fit_boot <- suppressWarnings(do.call(pfr, list(
                formula = pfr_formula, 
                data = pfr_data_boot, 
                family = pfr_family
              )))
              
              # Extract functional coefficient using coef.pfr or predict approach
              # Get smooth coefficients
              sm <- fit_boot$smooth[[1]]
              if(!is.null(sm)) {
                # Create prediction data for the smooth
                Xp <- mgcv::PredictMat(sm, data.frame(X_func.tmat = t_grid))
                coef_sm <- coef(fit_boot)[sm$first.para:sm$last.para]
                boot_coefs[b, ] <- as.vector(Xp %*% coef_sm)
                boot_success <- boot_success + 1
              }
              
              if(b %% 10 == 0) incProgress(10/B)
            }, error = function(e) {
              # Skip failed bootstrap iterations
            })
          }
        })
        
        if(boot_success >= B * 0.5) {  # Need at least 50% successful iterations
          # Calculate percentile CIs
          boot_ci_lower <- apply(boot_coefs, 2, quantile, probs = 0.025, na.rm = TRUE)
          boot_ci_upper <- apply(boot_coefs, 2, quantile, probs = 0.975, na.rm = TRUE)
          boot_se <- apply(boot_coefs, 2, sd, na.rm = TRUE)
          
          fit$boot_ci_lower <- boot_ci_lower
          fit$boot_ci_upper <- boot_ci_upper
          fit$boot_se <- boot_se
          fit$boot_t_grid <- t_grid
          fit$n_boot <- boot_success
          
          showNotification(paste("Bootstrap complete:", boot_success, "of", B, "iterations succeeded"), type = "message")
        } else {
          showNotification(paste("Bootstrap had too many failures:", boot_success, "of", B, "succeeded"), type = "warning")
        }
      }
      
      values$sofr_model <- fit
      
      family_msg <- paste0("SoFR Fitted! Family: ", fit$family$family, "(", fit$family$link, ")")
      showNotification(family_msg, type = "message")
      
    }, error = function(e) {
      showNotification(paste("SoFR Error:", e$message), type = "error")
      print(e)  # Also print to console for debugging
    })
  })
  
  output$sofr_inference_summary <- renderPrint({
    req(values$sofr_model)
    fit <- values$sofr_model
    cat("=== Scalar-on-Function Regression Results ===\n\n")
    cat("Family:", fit$family$family, "\n")
    cat("Link function:", fit$family$link, "\n")
    if(!is.null(fit$n_boot)) {
      cat("Bootstrap CIs: Yes (B =", fit$n_boot, ")\n")
    } else {
      cat("Bootstrap CIs: No (using parametric)\n")
    }
    cat("\n")
    summary(fit)
  })
  
  output$sofr_model_diagnostics <- renderUI({
    req(values$sofr_model)
    fit <- values$sofr_model
    
    # Calculate diagnostics
    y_obs <- fit$y_original
    y_pred_link <- fitted(fit)  # On link scale for GLM
    
    if(fit$family$family == "binomial") {
      # For binomial, fitted values are already probabilities
      y_pred_prob <- y_pred_link
      
      # Classification at 0.5 threshold
      y_pred_class <- ifelse(y_pred_prob > 0.5, 1, 0)
      acc <- mean(y_pred_class == y_obs)
      
      # Pseudo R-squared (McFadden)
      null_dev <- fit$null.deviance
      res_dev <- fit$deviance
      pseudo_r2 <- 1 - (res_dev / null_dev)
      
      tagList(
        h4("Model Diagnostics"),
        tags$table(class = "table table-condensed",
                   tags$tr(tags$td("Null Deviance:"), tags$td(round(null_dev, 2))),
                   tags$tr(tags$td("Residual Deviance:"), tags$td(round(res_dev, 2))),
                   tags$tr(tags$td("McFadden Pseudo-R²:"), tags$td(round(pseudo_r2, 3))),
                   tags$tr(tags$td("Accuracy (at p=0.5):"), tags$td(paste0(round(acc * 100, 1), "%")))
        )
      )
    } else {
      # For Gaussian, compute R-squared
      ss_res <- sum((y_obs - y_pred_link)^2)
      ss_tot <- sum((y_obs - mean(y_obs))^2)
      r2 <- 1 - ss_res/ss_tot
      rmse <- sqrt(mean((y_obs - y_pred_link)^2))
      
      tagList(
        h4("Model Diagnostics"),
        tags$table(class = "table table-condensed",
                   tags$tr(tags$td("R-squared:"), tags$td(round(r2, 3))),
                   tags$tr(tags$td("RMSE:"), tags$td(round(rmse, 3))),
                   tags$tr(tags$td("Residual Deviance:"), tags$td(round(fit$deviance, 2)))
        )
      )
    }
  })
  
  output$sofr_coeff_interpretation <- renderUI({
    req(values$sofr_model)
    fit <- values$sofr_model
    
    if(fit$family$family == "binomial" && fit$family$link == "logit") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Logistic):"), " The coefficient β(t) is on the log-odds scale. ",
          "Positive values at time t indicate that higher functional predictor values at that time ",
          "are associated with increased probability of Y=1. ",
          "exp(β(t)) gives the odds ratio for a unit increase in X(t)."
      )
    } else if(fit$family$family == "binomial" && fit$family$link == "probit") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Probit):"), " The coefficient β(t) represents the change in the ",
          "z-score (standard normal quantile) of P(Y=1) for a unit increase in X(t)."
      )
    } else if(fit$family$family == "poisson") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Poisson):"), " The coefficient β(t) is on the log scale. ",
          "exp(β(t)) gives the multiplicative change in the expected count for a unit increase in X(t)."
      )
    } else if(fit$family$link == "log") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Log link):"), " The coefficient β(t) is on the log scale. ",
          "exp(β(t)) gives the multiplicative effect on E(Y)."
      )
    } else {
      div(class = "alert alert-info", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Identity link):"), " The coefficient β(t) represents the ",
          "direct additive effect of X(t) on E(Y). A unit increase in X(t) changes E(Y) by β(t)."
      )
    }
  })
  
  output$sofr_coeff_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    df_coef <- tryCatch({ coef(fit) }, error = function(e) NULL)
    validate(need(!is.null(df_coef), "Could not extract functional coefficients."))
    
    x_col <- if("X_func.arg" %in% colnames(df_coef)) "X_func.arg" else colnames(df_coef)[1]
    
    # Determine y-axis label based on family
    y_label <- if(fit$family$family == "binomial" && fit$family$link == "logit") {
      "Coefficient (log-odds scale)"
    } else if(fit$family$link == "log") {
      "Coefficient (log scale)"
    } else {
      "Coefficient"
    }
    
    subtitle <- paste0("Family: ", fit$family$family, "(", fit$family$link, ")")
    
    # Check if bootstrap CIs are available
    if(!is.null(fit$boot_ci_lower) && !is.null(fit$boot_ci_upper)) {
      # Create a dataframe combining parametric and bootstrap CIs
      # Interpolate bootstrap CIs to match coef grid
      boot_t <- fit$boot_t_grid
      coef_t <- df_coef[[x_col]]
      
      # Interpolate bootstrap CIs to coefficient grid
      boot_lower_interp <- approx(boot_t, fit$boot_ci_lower, xout = coef_t, rule = 2)$y
      boot_upper_interp <- approx(boot_t, fit$boot_ci_upper, xout = coef_t, rule = 2)$y
      
      df_coef$boot_lower <- boot_lower_interp
      df_coef$boot_upper <- boot_upper_interp
      
      ci_label <- paste0("95% Bootstrap CI (B=", fit$n_boot, ")")
      subtitle <- paste0(subtitle, " | ", ci_label)
      
      g <- ggplot(df_coef, aes_string(x = x_col, y = "value")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = boot_lower, ymax = boot_upper), fill = "firebrick", alpha = 0.2) +
        geom_line(color = "firebrick", linewidth = 1) +
        theme_minimal() +
        labs(title = "Estimated Functional Coefficient Beta(t)", 
             subtitle = subtitle,
             x = "Time", y = y_label)
    } else {
      # Use parametric CIs from pfr
      g <- ggplot(df_coef, aes_string(x = x_col, y = "value")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = value - 1.96*se, ymax = value + 1.96*se), fill = "firebrick", alpha = 0.2) +
        geom_line(color = "firebrick", linewidth = 1) +
        theme_minimal() +
        labs(title = "Estimated Functional Coefficient Beta(t)", 
             subtitle = paste0(subtitle, " | 95% Parametric CI"),
             x = "Time", y = y_label)
    }
    
    ggplotly(g) %>% layout(showlegend = FALSE)
  })
  
  output$sofr_pred_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    y_obs <- fit$y_original
    y_pred <- fitted(fit)
    
    if(fit$family$family == "binomial") {
      # For binary: predicted probabilities
      df_pred <- data.frame(
        Observed = factor(y_obs, levels = c(0, 1), labels = c("0", "1")),
        Predicted_Prob = y_pred
      )
      
      g <- ggplot(df_pred, aes(x = Observed, y = Predicted_Prob, fill = Observed)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        scale_fill_manual(values = c("0" = "steelblue", "1" = "firebrick")) +
        theme_minimal() +
        labs(title = "Predicted Probabilities by Observed Class",
             x = "Observed Class", y = "Predicted Probability P(Y=1)") +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40")
      
      ggplotly(g)
    } else {
      # For continuous: scatter plot
      plot_ly(x = y_obs, y = y_pred, type = 'scatter', mode = 'markers', 
              marker = list(color = 'blue', opacity = 0.6)) %>%
        add_lines(x = c(min(y_obs), max(y_obs)), y = c(min(y_obs), max(y_obs)), 
                  line = list(color = 'red', dash = 'dash'), name = "Identity") %>%
        layout(title = "Observed vs Predicted", 
               xaxis = list(title = "Observed"), 
               yaxis = list(title = "Predicted"))
    }
  })
  
  # ROC Curve for binary outcomes
  output$sofr_roc_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(NULL)
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    
    # Validate data
    if(length(y_obs) == 0 || length(y_pred_prob) == 0) return(NULL)
    if(length(unique(y_obs)) < 2) return(NULL)
    
    # Calculate ROC curve
    thresholds <- sort(unique(c(0, y_pred_prob, 1)))
    n_thresh <- length(thresholds)
    
    tpr_vec <- numeric(n_thresh)
    fpr_vec <- numeric(n_thresh)
    
    for(i in seq_along(thresholds)) {
      thresh <- thresholds[i]
      pred_class <- ifelse(y_pred_prob >= thresh, 1, 0)
      tp <- sum(pred_class == 1 & y_obs == 1)
      fp <- sum(pred_class == 1 & y_obs == 0)
      tn <- sum(pred_class == 0 & y_obs == 0)
      fn <- sum(pred_class == 0 & y_obs == 1)
      
      tpr_vec[i] <- if(tp + fn > 0) tp / (tp + fn) else 0
      fpr_vec[i] <- if(fp + tn > 0) fp / (fp + tn) else 0
    }
    
    roc_data <- data.frame(
      threshold = thresholds,
      tpr = tpr_vec,
      fpr = fpr_vec
    )
    roc_data <- roc_data[order(roc_data$fpr, roc_data$tpr), ]
    
    # Calculate AUC using trapezoidal rule
    auc <- 0
    if(nrow(roc_data) > 1) {
      for(i in 2:nrow(roc_data)) {
        auc <- auc + (roc_data$fpr[i] - roc_data$fpr[i-1]) * (roc_data$tpr[i] + roc_data$tpr[i-1]) / 2
      }
    }
    
    # Create plot without text attribute to avoid size mismatch
    p <- plot_ly() %>%
      add_trace(data = roc_data, x = ~fpr, y = ~tpr, type = 'scatter', mode = 'lines',
                line = list(color = 'firebrick', width = 2),
                name = paste0("ROC (AUC = ", round(auc, 3), ")"),
                hovertemplate = "FPR: %{x:.3f}<br>TPR: %{y:.3f}<extra></extra>") %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 1, 
                   line = list(color = 'gray', dash = 'dash'),
                   name = "Random", showlegend = FALSE) %>%
      layout(title = paste0("ROC Curve (AUC = ", round(auc, 3), ")"),
             xaxis = list(title = "False Positive Rate (1 - Specificity)", range = c(0, 1)),
             yaxis = list(title = "True Positive Rate (Sensitivity)", range = c(0, 1)))
    p
  })
  
  # Calibration plot for binary outcomes
  output$sofr_calibration_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(NULL)
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    
    # Validate data
    if(length(y_obs) == 0 || length(y_pred_prob) == 0) return(NULL)
    
    # Create calibration bins
    n_bins <- 10
    df_cal <- data.frame(pred = y_pred_prob, obs = y_obs)
    df_cal$bin <- cut(df_cal$pred, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
    
    cal_summary <- df_cal %>%
      group_by(bin) %>%
      summarise(
        mean_pred = mean(pred),
        mean_obs = mean(obs),
        n = n(),
        se = sqrt(mean_obs * (1 - mean_obs) / n),
        .groups = "drop"
      ) %>%
      filter(n >= 3)  # Only include bins with enough observations
    
    if(nrow(cal_summary) == 0) {
      return(plotly_empty() %>% layout(title = "Not enough data for calibration plot"))
    }
    
    # Convert to regular data frame to avoid tibble issues
    cal_summary <- as.data.frame(cal_summary)
    
    p <- plot_ly() %>%
      add_trace(data = cal_summary, x = ~mean_pred, y = ~mean_obs, 
                type = 'scatter', mode = 'markers',
                marker = list(size = sqrt(cal_summary$n) * 3, color = 'steelblue', opacity = 0.7),
                error_y = list(array = 1.96 * cal_summary$se, color = 'steelblue'),
                hovertemplate = "Pred: %{x:.3f}<br>Obs: %{y:.3f}<br>n=%{text}<extra></extra>",
                text = cal_summary$n,
                name = "Observed") %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 1, 
                   line = list(color = 'red', dash = 'dash'),
                   name = "Perfect calibration") %>%
      layout(title = "Calibration Plot",
             xaxis = list(title = "Mean Predicted Probability", range = c(0, 1)),
             yaxis = list(title = "Observed Proportion", range = c(0, 1)))
    p
  })
  
  # Classification metrics for binary outcomes
  output$sofr_classification_metrics <- renderPrint({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(cat("Classification metrics only available for binary outcomes."))
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    y_pred_class <- ifelse(y_pred_prob > 0.5, 1, 0)
    
    # Confusion matrix
    tp <- sum(y_pred_class == 1 & y_obs == 1)
    fp <- sum(y_pred_class == 1 & y_obs == 0)
    tn <- sum(y_pred_class == 0 & y_obs == 0)
    fn <- sum(y_pred_class == 0 & y_obs == 1)
    
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    ppv <- if(tp + fp > 0) tp / (tp + fp) else NA
    npv <- if(tn + fn > 0) tn / (tn + fn) else NA
    f1 <- if(ppv + sensitivity > 0) 2 * (ppv * sensitivity) / (ppv + sensitivity) else NA
    
    # AUC calculation
    thresholds <- sort(unique(c(0, y_pred_prob, 1)))
    tpr_vec <- fpr_vec <- numeric(length(thresholds))
    for(i in seq_along(thresholds)) {
      pred_class <- ifelse(y_pred_prob >= thresholds[i], 1, 0)
      tpr_vec[i] <- sum(pred_class == 1 & y_obs == 1) / sum(y_obs == 1)
      fpr_vec[i] <- sum(pred_class == 1 & y_obs == 0) / sum(y_obs == 0)
    }
    ord <- order(fpr_vec, tpr_vec)
    auc <- sum(diff(fpr_vec[ord]) * (tpr_vec[ord][-1] + tpr_vec[ord][-length(tpr_vec)]) / 2)
    
    cat("=== Classification Metrics (threshold = 0.5) ===\n\n")
    cat("Confusion Matrix:\n")
    cat("                 Predicted\n")
    cat("                  0      1\n")
    cat(sprintf("Actual 0     %5d  %5d\n", tn, fp))
    cat(sprintf("       1     %5d  %5d\n", fn, tp))
    cat("\n")
    cat(sprintf("Accuracy:    %.3f\n", accuracy))
    cat(sprintf("Sensitivity: %.3f (Recall)\n", sensitivity))
    cat(sprintf("Specificity: %.3f\n", specificity))
    cat(sprintf("PPV:         %.3f (Precision)\n", ppv))
    cat(sprintf("NPV:         %.3f\n", npv))
    cat(sprintf("F1 Score:    %.3f\n", f1))
    cat(sprintf("AUC:         %.3f\n", auc))
  })
