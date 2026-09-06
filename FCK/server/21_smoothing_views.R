# ==========================================================================
# server/21_smoothing_views.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 1406-1494  (smoothing fit statistics printout)
#   CIRCAREG.R lines 913-930  (compact smoothing fit-metrics panel)
#   WaPaa1_3.R lines 2862-3016  (interactive smoothed-curve plot + curve selection)
# ==========================================================================
  output$smoothing_fit_summary <- renderPrint({
    if(!is.null(values$smooth_fit_metrics)) {
      cat("=== Smoothing Quality Metrics ===\n\n")
      
      # Use CIRCAREG-style metrics structure (list with mean_r_squared, mean_rmse, etc.)
      metrics <- values$smooth_fit_metrics
      
      cat("Mean R2:", sprintf("%.3f", metrics$mean_r_squared), 
          sprintf("(SD: %.3f)", metrics$sd_r_squared), "\n")
      if(!is.na(metrics$mean_r_squared)) {
        if(metrics$mean_r_squared >= 0.9) {
          cat("  -> Excellent fit!\n")
        } else if(metrics$mean_r_squared >= 0.7) {
          cat("  -> Good fit\n")
        } else {
          cat("  -> Fair fit (consider adjusting parameters)\n")
        }
      }
      
      cat("\nMean RMSE:", sprintf("%.3f", metrics$mean_rmse),
          sprintf("(SD: %.3f)", metrics$sd_rmse), "\n")
      if(!is.null(metrics$rel_rmse_pct) && !is.na(metrics$rel_rmse_pct)) {
        cat(sprintf("  Relative RMSE: %.1f%% of data range", metrics$rel_rmse_pct))
        rmse_label <- if(metrics$rel_rmse_pct < 5)  " -> Excellent (<5%)" else
                      if(metrics$rel_rmse_pct < 10) " -> Good (5-10%)" else
                      if(metrics$rel_rmse_pct < 20) " -> Acceptable (10-20%)" else
                                                     " -> Poor (>20%)"
        cat(rmse_label, "\n")
      }

      # Number of subjects
      n_subjects <- length(metrics$r_squared)
      cat("\nNumber of subjects smoothed:", n_subjects, "\n")
      
      # Count subjects with valid metrics
      n_valid <- sum(!is.na(metrics$r_squared))
      if(n_valid < n_subjects) {
        cat("Subjects with valid fit metrics:", n_valid, "/", n_subjects, "\n")
      }
      
      # Display smoothing parameters
      cat("\n--- Smoothing Parameters ---\n")
      
      if(!is.null(metrics$method)) {
        # AUDIT (P5.10): this said "Automatic (REML)". The app's automatic
        # smoother selects lambda by GCV on fda::smooth.basis (fck_auto_lambda);
        # it has never used REML, and smooth.basis has no REML option. A user
        # reading this label would write the wrong method in a paper.
        method_label <- switch(metrics$method,
                               "auto" = "Automatic (lambda by GCV)",
                               "manual" = "Manual",
                               "none" = "None (Raw data)")
        cat("Method:", method_label, "\n")
      }

      # MERGED APP: which time axis these numbers came from.
      if(!is.null(metrics$time_axis)) {
        cat("Time axis:", metrics$time_axis, "\n")
      }
      
      if(!is.null(metrics$n_basis)) {
        cat("Number of B-splines:", metrics$n_basis, "\n")
      }
      
      if(!is.null(metrics$lambda)) {
        if(metrics$lambda == 0) {
          # P5.10: lambda = 0 in fda is the UNPENALISED fit. It is not
          # "automatic optimization" and it is not REML.
          cat("Lambda: 0 (UNPENALISED interpolating fit -- no penalty applied)\n")
        } else {
          smooth_factor <- -log10(metrics$lambda)
          cat(sprintf("Lambda: %.3e (smoothing factor = %.2f)\n",
                      metrics$lambda, smooth_factor))
        }
      }

      # EDF-based n-basis recommendation
      if(!is.null(metrics$mean_df) && !is.na(metrics$mean_df)) {
        cat("\n--- Spline Complexity ---\n")
        cat(sprintf("Mean EDF (effective df): %.1f (SD: %.1f, max: %.1f)\n",
                    metrics$mean_df, metrics$sd_df, metrics$max_df))
        recommended_min <- ceiling(metrics$max_df) + 2
        cat("Current n_basis:", metrics$n_basis, "\n")
        cat(sprintf("Recommended minimum: ceil(max EDF) + 2 = %d\n", recommended_min))
        if(metrics$n_basis < recommended_min) {
          cat(sprintf("  ⚠ n_basis may be too low — consider increasing to at least %d\n",
                      recommended_min))
        } else if(metrics$mean_df < metrics$n_basis / 3) {
          cat("  ℹ Many basis functions are inactive (EDF << n_basis).\n")
          cat("    You may reduce n_basis for a more parsimonious model.\n")
        } else {
          cat("  OK — current n_basis provides adequate flexibility.\n")
        }
      }

    } else {
      cat("No smoothing applied yet.\n\n")
      cat("Click 'Apply Smoothing' to see fit quality metrics.")
    }
  })

  output$smooth_fit_display <- renderUI({
    if(is.null(values$smooth_fit_metrics)) return(NULL)
    
    metrics <- values$smooth_fit_metrics
    
    div(
      style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border: 1px solid #dee2e6;",
      h5(icon("chart-line"), " Smoothing Fit Metrics", style = "margin-top: 0;"),
      p(style = "margin-bottom: 5px;",
        strong("Mean R²: "), sprintf("%.3f (SD: %.3f)", metrics$mean_r_squared, metrics$sd_r_squared)),
      p(style = "margin-bottom: 5px;",
        strong("Mean RMSE: "), sprintf("%.3f (SD: %.3f)", metrics$mean_rmse, metrics$sd_rmse)),
      helpText(style = "font-size: 0.85em; margin-top: 8px;",
               "R² = proportion of variance explained by smooth curve (higher = better fit).",
               br(),
               "RMSE = average deviation from original values (lower = better fit).")
    )
  })

  # Data visualization plot - INTERACTIVE with curve selection
  output$data_plot <- renderPlotly({
    if(is.null(values$data)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>%
               layout(title = "No data loaded"))
    }

    tryCatch({
      data_to_plot <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

      n_time <- ncol(data_to_plot)
      n_subj <- nrow(data_to_plot)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      hover_times <- hover_time_labels(time_points)

      # Get currently selected curve (if any)
      selected_idx <- values$selected_curve

      p <- plot_ly(type = 'scatter', mode = 'lines', source = "data_plot_source")

      # Add individual curves (use 0-1 normalized time)
      # MERGED APP: how many curves to draw is the user's choice. WaPaa fixed
      # this at 50 and gave no sign that anything had been left out.
      n_req <- suppressWarnings(as.numeric(input$data_plot_n %||% 50))
      if (!is.finite(n_req) || n_req <= 0) n_req <- n_subj      # "All"
      n_show <- min(n_subj, as.integer(n_req))
      for(i in 1:n_show) {
        # Determine if this curve is selected
        is_selected <- !is.null(selected_idx) && i == selected_idx

        # Set color and width based on selection
        curve_color <- if(is_selected) 'rgba(0, 100, 255, 0.9)' else 'rgba(100, 100, 100, 0.3)'
        curve_width <- if(is_selected) 3 else 1

        # Get group label if available
        group_info <- if(!is.null(values$group_labels) && i <= length(values$group_labels)) {
          paste0(" (Group: ", values$group_labels[i], ")")
        } else {
          ""
        }

        p <- p %>% add_trace(x = time_points,
                             y = data_to_plot[i,],
                             type = 'scatter',
                             mode = 'lines',
                             name = paste0("Subject ", i, group_info),
                             line = list(color = curve_color, width = curve_width),
                             hovertemplate = paste0("Subject ", i, group_info,
                                                    "<br>Time: %{customdata}<br>Value: %{y:.2f}<extra></extra>"),
                             customdata = hover_times,
                             showlegend = FALSE)
      }

      # Add mean curve
      p <- p %>% add_trace(x = time_points,
                           y = colMeans(data_to_plot),
                           type = 'scatter',
                           mode = 'lines',
                           name = "Mean",
                           line = list(color = 'red', width = 3),
                           hovertemplate = "Mean<br>Time: %{customdata}<br>Value: %{y:.2f}<extra></extra>",
                           customdata = hover_times)

      p <- p %>% layout(title = "Functional Data (click curve to select)",
                        yaxis = list(title = "Value"),
                        showlegend = FALSE,
                        hovermode = 'closest',
                        clickmode = 'event')

      # Apply time axis formatting with hour labels
      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_preprocess))
      
      p
    }, error = function(e) {
      cat("Data plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>%
        layout(title = paste("Error:", e$message))
    })
  })

  # ============================================================================
  # CURVE SELECTION HANDLERS
  # ============================================================================

  # Handle click events on the data plot
  observeEvent(event_data("plotly_click", source = "data_plot_source"), {
    click_data <- event_data("plotly_click", source = "data_plot_source")

    if(!is.null(click_data)) {
      # Get the curve number from the click (curveNumber is 0-indexed)
      curve_num <- click_data$curveNumber + 1  # Convert to 1-indexed

      # The last curve is the mean (index 0 in customdata), skip it
      n_subj <- nrow(values$data)
      # MERGED APP: how many curves to draw is the user's choice. WaPaa fixed
      # this at 50 and gave no sign that anything had been left out.
      n_req <- suppressWarnings(as.numeric(input$data_plot_n %||% 50))
      if (!is.finite(n_req) || n_req <= 0) n_req <- n_subj      # "All"
      n_show <- min(n_subj, as.integer(n_req))

      if(curve_num <= n_show) {
        # It's an individual curve, select it
        values$selected_curve <- curve_num
        cat("Selected curve:", curve_num, "\n")
      } else {
        # Clicked on mean curve, do nothing or clear selection
        cat("Clicked on mean curve\n")
      }
    }
  })

  # Clear curve selection button
  observeEvent(input$clear_curve_selection, {
    values$selected_curve <- NULL
  })

  # Display selected curve information
  output$selected_curve_info <- renderText({
    if(is.null(values$selected_curve) || is.null(values$data)) {
      return("No curve selected.")
    }

    tryCatch({
      idx <- values$selected_curve
      data_to_use <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

      if(idx > nrow(data_to_use)) {
        return("Invalid selection.")
      }

      curve_data <- data_to_use[idx, ]

      # Basic statistics
      curve_mean <- mean(curve_data, na.rm = TRUE)
      curve_sd <- sd(curve_data, na.rm = TRUE)
      curve_min <- min(curve_data, na.rm = TRUE)
      curve_max <- max(curve_data, na.rm = TRUE)
      curve_range <- curve_max - curve_min

      # Group info
      group_info <- if(!is.null(values$group_labels) && idx <= length(values$group_labels)) {
        paste0("Group: ", values$group_labels[idx])
      } else {
        "Group: N/A"
      }

      paste0(
        "Subject: ", idx, "\n",
        group_info, "\n",
        "-------------------\n",
        "Mean: ", sprintf("%.3f", curve_mean), "\n",
        "SD: ", sprintf("%.3f", curve_sd), "\n",
        "Min: ", sprintf("%.3f", curve_min), "\n",
        "Max: ", sprintf("%.3f", curve_max), "\n",
        "Range: ", sprintf("%.3f", curve_range)
      )
    }, error = function(e) {
      paste("Error:", e$message)
    })
  })
