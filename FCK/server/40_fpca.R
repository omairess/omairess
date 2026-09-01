# ==========================================================================
# server/40_fpca.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 3018-4423  (group summary, fPCA/warping analysis + outputs)
#   WaPaa1_3.R lines 6229-6522  (warping / alignment / landmark plots)
# ==========================================================================
  # Group summary table
  output$group_summary <- renderDT({
    if(!is.null(values$group_labels)) {
      group_summary <- data.frame(
        Group = names(table(values$group_labels)),
        Count = as.numeric(table(values$group_labels)),
        Percentage = round(100 * as.numeric(table(values$group_labels)) / length(values$group_labels), 1)
      )
      datatable(group_summary, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
    } else {
      return(NULL)
    }
  })
  
  # Group preview plot
  output$group_preview_plot <- renderPlot({
    if(!is.null(values$data) && !is.null(values$group_labels)) {
      n_time <- ncol(values$data)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      hour_labels <- get_hour_labels()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      
      groups <- levels(values$group_labels)
      n_groups <- length(groups)
      
      # Create scalable color palette
      base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
      colors <- colorRampPalette(base_cols)(n_groups)
      
      matplot(time_points, t(values$data), type = "l", 
              col = rgb(0.5, 0.5, 0.5, 0.1), lty = 1,
              xlab = time_label, ylab = "Value", 
              main = "Functional Data by Group")
      
      for(i in 1:n_groups) {
        group_idx <- which(values$group_labels == groups[i])
        if(length(group_idx) > 0) {
          group_mean <- colMeans(values$data[group_idx, , drop = FALSE])
          lines(time_points_plot, group_mean, col = colors[i], lwd = 3)
        }
      }
      
      legend("topright", legend = groups, col = colors, lwd = 3)
    }
  })
  
  # Run PCA analysis - FIXED with better error handling
  observeEvent(input$run_analysis, {
    cat("Running PCA analysis...\n")
    
    if(is.null(values$data)) {
      showNotification("No data loaded", type = "error", duration = 5)
      return()
    }
    
    tryCatch({
      data_to_use <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      
      n <- ncol(data_to_use)
      m <- nrow(data_to_use)
      
      n_components <- min(input$n_components, m-1, 10)
      time_grid <- seq(0, 1, length.out = n)
      
      if(is.null(values$fd_obj)) {
        # MERGED APP: one rule for this, in FCK/server/04_helpers_fd.R. No
        # smoothing is invented here; an interpolating basis is used and the
        # user is told. (WaPaa silently smoothed onto min(20, n-2) bases.)
        if(!fck_ensure_fd_obj(values)) return()
        cat("Created fd_obj with an interpolating basis (no smoothing)\n")
      }
      
      if(input$pca_type == "fpca") {
        # Regular functional PCA
        cat("Running standard functional PCA...\n")
        values$pca_results <- pca.fd(values$fd_obj, nharm = n_components)
        values$warping_results <- NULL
        showNotification("Functional PCA completed!", type = "message", duration = 5)
        
      } else if(input$pca_type == "twpca") {
        # Time-warped PCA with proper error handling
        cat("Running time-warped PCA...\n")
        showNotification("Running time-warped PCA...", type = "message", duration = 2)
        
        # Get warping method settings
        warping_method <- if(!is.null(input$warping_method)) input$warping_method else "linear_shift"
        
        cat("Using warping method:", warping_method, "\n")
        
        # Perform warping based on method
        warp_fd <- NULL
        
        if(warping_method == "linear_shift") {
          periodic <- if(!is.null(input$periodic_shift)) input$periodic_shift else FALSE
          allow_dilation <- if(!is.null(input$allow_dilation)) input$allow_dilation else FALSE
          dilation_range <- if(!is.null(input$dilation_range)) input$dilation_range else c(0.95, 1.05)
          reference <- if(!is.null(input$shift_reference)) input$shift_reference else "mean"
          
          warp_fd <- linear_shift_alignment(values$fd_obj,
                                            periodic = periodic,
                                            allow_dilation = allow_dilation,
                                            dilation_range = dilation_range,
                                            reference = reference,
                                            time_points = time_grid)
          
        } else if(warping_method == "parametric") {
          family <- if(!is.null(input$parametric_family)) input$parametric_family else "power"
          param_range <- if(!is.null(input$param_range)) input$param_range else c(0.5, 2)
          symmetric <- if(!is.null(input$symmetric_warp)) input$symmetric_warp else FALSE
          
          warp_fd <- parametric_alignment(values$fd_obj, 
                                          family = family,
                                          param_range = param_range,
                                          symmetric = symmetric,
                                          time_points = time_grid)
          
        } else if(warping_method == "landmark") {
          # Simple landmark alignment
          landmarks <- c(0.25, 0.5, 0.75)  # Default landmarks
          warp_fd <- landmark_alignment_simple(values$fd_obj, landmarks, time_grid)
        }
        
        # Check if warping was successful
        if(is.null(warp_fd)) {
          cat("Warping failed, using default alignment\n")
          warp_fd <- linear_shift_alignment(values$fd_obj,
                                            periodic = FALSE,
                                            allow_dilation = FALSE,
                                            reference = "mean",
                                            time_points = time_grid)
        }

        # Calculate warping fit statistics
        if(!is.null(warp_fd$registered_curves) && !is.null(warp_fd$warp_functions)) {
          cat("Calculating warping fit statistics...\n")
          original_curves <- eval.fd(time_grid, values$fd_obj)
          warp_fd$fit_statistics <- calculate_warping_fit_statistics(
            original_curves = original_curves,
            registered_curves = warp_fd$registered_curves,
            warp_functions = warp_fd$warp_functions,
            time_points = time_grid
          )
          if(!is.null(warp_fd$fit_statistics)) {
            cat("Warping fit statistics calculated successfully\n")
            cat("  Mean R²:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_r_squared), "\n")
            cat("  Mean RMSE:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_rmse), "\n")
            cat("  Variance explained by warping:",
                sprintf("%.2f%%", warp_fd$fit_statistics$summary$variance_explained_by_warping * 100), "\n")
          }
        }

        # Store warping results
        values$warping_results <- warp_fd
        
        # Run PCA on warped curves with error checking
        if(!is.null(warp_fd$regfd)) {
          cat("Running PCA on warped fd object\n")
          values$pca_results <- pca.fd(warp_fd$regfd, nharm = n_components)
        } else if(!is.null(warp_fd$registered_curves) && ncol(warp_fd$registered_curves) > 0) {
          cat("Creating fd from registered curves\n")
          # Create fd from registered curves
          # Use lambda=0: registered curves are already processed, just need fd representation
          reg_basis <- create.bspline.basis(c(0, 1), nbasis = min(20, n-2, max(4, floor(n/3))))
          reg_fd <- smooth.basis(time_grid, warp_fd$registered_curves, 
                                 fdPar(reg_basis, 2, 0))$fd
          values$pca_results <- pca.fd(reg_fd, nharm = n_components)
        } else {
          cat("Warning: No valid warped data, using original fd_obj\n")
          values$pca_results <- pca.fd(values$fd_obj, nharm = n_components)
        }
        
        showNotification("Time-warped PCA completed!", type = "message", duration = 5)
      }
      
      # Verify PCA results
      if(!is.null(values$pca_results)) {
        cat("PCA complete - PC1 variance:", 
            round(values$pca_results$varprop[1]*100, 2), "%\n")
      }
      
    }, error = function(e) {
      cat("Error in PCA:", e$message, "\n")
      showNotification(paste("Analysis error:", e$message), type = "error", duration = 10)
      values$pca_results <- NULL
    })
  })
  
  # PCA status output - FIXED
  output$pca_status <- renderText({
    if(is.null(values$pca_results)) {
      "No PCA results available. Run analysis first."
    } else {
      tryCatch({
        pca_res <- values$pca_results
        n_comp <- if(!is.null(pca_res$scores)) ncol(pca_res$scores) else 0
        analysis_type <- if(!is.null(values$warping_results)) "Time-warped PCA" else "Standard PCA"
        paste(analysis_type, "complete:", n_comp, "components extracted")
      }, error = function(e) {
        "PCA results available but incomplete"
      })
    }
  })
  
  # PCA summary
  output$pca_summary <- renderPrint({
    if(is.null(values$pca_results)) {
      cat("Run PCA analysis first\n")
    } else {
      tryCatch({
        pca_res <- values$pca_results
        
        cat("PCA Results Summary:\n")
        cat("-------------------\n")
        
        if(!is.null(pca_res$scores)) {
          cat("Number of components:", ncol(pca_res$scores), "\n")
          cat("Number of subjects:", nrow(pca_res$scores), "\n")
        }
        
        if(!is.null(pca_res$varprop)) {
          cat("Variance explained:\n")
          for(i in 1:min(3, length(pca_res$varprop))) {
            cat(sprintf("  PC%d: %.2f%%\n", i, pca_res$varprop[i] * 100))
          }
        }
        
        if(!is.null(values$warping_results)) {
          cat("\nWarping method:", values$warping_results$method, "\n")
          if(!is.null(values$warping_results$shifts)) {
            cat("Mean shift:", round(mean(abs(values$warping_results$shifts)), 4), "\n")
          }
        }
      }, error = function(e) {
        cat("Error displaying summary:", e$message, "\n")
      })
    }
  })
  
  # Component loadings plot - FIXED
  output$loadings_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      # Check validity
      if(is.null(pca_res$meanfd) || is.null(pca_res$harmonics)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      time_points <- seq(0, 1, length.out = 100)
      mean_vals <- eval.fd(time_points, pca_res$meanfd)
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = time_points, 
                  y = as.vector(mean_vals),
                  type = 'scatter',
                  mode = 'lines',
                  name = "Mean", 
                  line = list(color = 'black', width = 2))
      
      colors <- c('red', 'blue', 'green', 'orange', 'purple')
      n_comp <- min(ncol(pca_res$scores), 5, length(pca_res$harmonics))
      
      for(i in 1:n_comp) {
        if(i <= length(pca_res$harmonics) && !is.null(pca_res$harmonics[i])) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          p <- p %>% add_trace(x = time_points,
                               y = as.vector(loading_vals),
                               type = 'scatter',
                               mode = 'lines',
                               name = paste("PC", i),
                               line = list(color = colors[i], width = 2))
        }
      }
      
      p <- p %>% layout(title = "Principal Component Loadings",
                        yaxis = list(title = "Loading"))

      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_results))
      p
      
    }, error = function(e) {
      cat("Loadings plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Variance explained plot - FIXED
  output$variance_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      if(is.null(pca_res$varprop) || is.null(pca_res$scores)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      n_show <- min(length(pca_res$varprop), ncol(pca_res$scores))
      
      if(n_show == 0) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "No components available"))
      }
      
      var_prop <- pca_res$varprop[1:n_show]
      cum_var <- cumsum(var_prop)
      
      p <- plot_ly() %>%
        add_trace(x = 1:n_show,
                  y = var_prop * 100,
                  type = 'bar',
                  name = 'Individual',
                  marker = list(color = 'lightblue'))
      
      p <- p %>%
        add_trace(x = 1:n_show,
                  y = cum_var * 100,
                  type = 'scatter',
                  mode = 'lines+markers',
                  name = 'Cumulative',
                  yaxis = 'y2',
                  line = list(color = 'red'),
                  marker = list(color = 'red'))
      
      p %>% layout(title = "Variance Explained",
                   xaxis = list(title = "Component", dtick = 1),
                   yaxis = list(title = "Variance Explained (%)"),
                   yaxis2 = list(title = "Cumulative Variance (%)",
                                 overlaying = 'y',
                                 side = 'right'))
      
    }, error = function(e) {
      cat("Variance plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Component scores plot - FIXED
  output$scores_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      # Check if PCA results are valid
      if(is.null(pca_res$meanfd) || is.null(pca_res$harmonics) || is.null(pca_res$values)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      time_points <- seq(0, 1, length.out = 100)
      mean_vals <- eval.fd(time_points, pca_res$meanfd)
      
      effect_mult <- if(!is.null(input$effect_size)) input$effect_size else 1
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = time_points, 
                  y = as.vector(mean_vals),
                  type = 'scatter',
                  mode = 'lines',
                  name = "Mean", 
                  line = list(color = 'black', width = 3))
      
      colors <- c('red', 'blue', 'green')
      n_show <- min(3, ncol(pca_res$scores), length(pca_res$values))
      
      for(i in 1:n_show) {
        if(i <= length(pca_res$harmonics) && !is.null(pca_res$harmonics[i])) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          
          # Plus 2 SD
          plus_2sd <- as.vector(mean_vals + effect_mult * 2 * sqrt(pca_res$values[i]) * loading_vals)
          p <- p %>% add_trace(
            x = time_points,
            y = plus_2sd,
            type = 'scatter',
            mode = 'lines',
            name = paste("PC", i, "+2SD"),
            line = list(color = colors[i], dash = 'solid')
          )
          
          # Minus 2 SD
          minus_2sd <- as.vector(mean_vals - effect_mult * 2 * sqrt(pca_res$values[i]) * loading_vals)
          p <- p %>% add_trace(
            x = time_points,
            y = minus_2sd,
            type = 'scatter',
            mode = 'lines',
            name = paste("PC", i, "-2SD"),
            line = list(color = colors[i], dash = 'dash')
          )
        }
      }
      
      p <- p %>% layout(title = "Effect of Component Scores",
                        yaxis = list(title = "Value"))
      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_results))
      p
      
    }, error = function(e) {
      cat("Scores plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Component scores table - FIXED
  output$scores_table <- renderDT({
    if(is.null(values$pca_results)) {
      return(NULL)
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      if(is.null(pca_res$scores)) {
        return(datatable(data.frame(Message = "No scores available"),
                         options = list(dom = 't'),
                         rownames = FALSE))
      }
      
      scores_df <- data.frame(
        Subject = 1:nrow(pca_res$scores)
      )
      
      # Add group information if available
      if(!is.null(values$group_labels) && length(values$group_labels) == nrow(pca_res$scores)) {
        scores_df$Group <- values$group_labels
      }
      
      n_comp <- ncol(pca_res$scores)
      for(i in 1:n_comp) {
        scores_df[paste0("PC", i)] <- round(pca_res$scores[,i], 3)
      }
      
      datatable(scores_df, 
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
      
    }, error = function(e) {
      cat("Scores table error:", e$message, "\n")
      datatable(data.frame(Error = "Unable to display scores"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })
  
  # Warping scores table
  output$warping_scores <- renderDT({
    if(is.null(values$warping_results)) {
      return(NULL)
    }
    
    tryCatch({
      warp_results <- values$warping_results
      n_subjects <- ncol(values$fd_obj$coefs)
      
      # Calculate warping amplitude (deviation from identity)
      if(!is.null(warp_results$warp_functions)) {
        n_time <- nrow(warp_results$warp_functions)
        time_points <- seq(0, 1, length.out = n_time)
        
        warp_amplitude <- numeric(n_subjects)
        for(i in 1:n_subjects) {
          # Calculate deviation from identity line
          warp_amplitude[i] <- sqrt(mean((warp_results$warp_functions[,i] - time_points)^2))
        }
      } else if(!is.null(warp_results$shifts)) {
        # Use shifts if available
        warp_amplitude <- abs(warp_results$shifts)
      } else if(!is.null(warp_results$alpha_values)) {
        # Use alpha values for parametric warping
        warp_amplitude <- abs(warp_results$alpha_values - 1)
      } else {
        warp_amplitude <- rep(0, n_subjects)
      }
      
      warping_df <- data.frame(
        Subject = 1:n_subjects,
        Warping_Amplitude = round(warp_amplitude, 4),
        Method = warp_results$method
      )

      # Add fit statistics if available
      if(!is.null(warp_results$fit_statistics) && !is.null(warp_results$fit_statistics$per_subject)) {
        fit_stats <- warp_results$fit_statistics$per_subject
        warping_df$R_squared <- fit_stats$R_squared
        warping_df$RMSE <- fit_stats$RMSE
        warping_df$Correlation <- fit_stats$Correlation
      }

      # Add group information if available
      if(!is.null(values$group_labels)) {
        warping_df$Group <- values$group_labels
      }

      datatable(warping_df,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE) %>%
        formatStyle("Warping_Amplitude",
                    backgroundColor = styleInterval(c(0.01, 0.05, 0.1),
                                                    c("white", "#ffffcc", "#ffcccc", "#ff9999"))) %>%
        formatStyle("R_squared",
                    backgroundColor = styleInterval(c(0.5, 0.7, 0.9),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("RMSE",
                    backgroundColor = styleInterval(c(0.05, 0.1, 0.2),
                                                    c("#99ff99", "#ccffcc", "#ffffcc", "#ffcccc")))
      
    }, error = function(e) {
      cat("Warping scores table error:", e$message, "\n")
      datatable(data.frame(Message = "No warping scores available"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })

  # ============================================================================
  # WARPING FIT STATISTICS OUTPUTS
  # ============================================================================

  # Summary statistics output (averaged over subjects)
  output$warping_fit_summary <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see fit statistics.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary
      n_valid <- values$warping_results$fit_statistics$n_valid
      n_total <- values$warping_results$fit_statistics$n_subjects

      paste0(
        "Subjects analyzed: ", n_valid, "/", n_total, "\n\n",
        "R² (warped vs original, functional):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_r_squared),
        " (SD: ", sprintf("%.4f", stats$sd_r_squared), ")\n",
        "  [1 - ∫(f-g)²dt / ∫(f-f̄)²dt]\n\n",
        "RMSE (warped vs original):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_rmse),
        " (SD: ", sprintf("%.4f", stats$sd_rmse), ")\n\n",
        "Correlation (orig vs warped):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_correlation),
        " (SD: ", sprintf("%.4f", stats$sd_correlation), ")\n\n",
        "MAE (warped vs original):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_mae),
        " (SD: ", sprintf("%.4f", stats$sd_mae), ")\n\n",
        "Warping Amplitude:\n",
        "  Mean: ", sprintf("%.4f", stats$mean_warp_amplitude),
        " (SD: ", sprintf("%.4f", stats$sd_warp_amplitude), ")"
      )
    }, error = function(e) {
      paste("Error displaying fit statistics:", e$message)
    })
  })

  # Variance decomposition output (EFDA-style)
  output$warping_variance_decomposition <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see variance decomposition.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary

      # Calculate percentages
      total_var <- stats$total_orig_variance
      amp_pct <- if(total_var > 0) stats$total_amp_variance / total_var * 100 else NA
      phase_pct <- if(total_var > 0) stats$total_phase_variance / total_var * 100 else NA
      explained_pct <- stats$variance_explained_by_warping * 100

      paste0(
        "EFDA Variance Decomposition:\n",
        "============================\n\n",
        "Total Original Variance: ", sprintf("%.4f", total_var), "\n\n",
        "After Alignment:\n",
        "  Amplitude Variance: ", sprintf("%.4f", stats$total_amp_variance),
        " (", sprintf("%.1f%%", amp_pct), " of original)\n",
        "  Phase Variance: ", sprintf("%.4f", stats$total_phase_variance),
        " (", sprintf("%.1f%%", phase_pct), " of original)\n\n",
        "Variance Explained by Warping:\n",
        "  ", sprintf("%.2f%%", explained_pct), "\n\n",
        "Elastic Distances (mean):\n",
        "  Full Distance: ", sprintf("%.4f", stats$mean_full_distance), "\n",
        "  Amplitude Distance: ", sprintf("%.4f", stats$mean_elastic_amp_dist), "\n",
        "  Phase Distance: ", sprintf("%.4f", stats$mean_elastic_phase_dist)
      )
    }, error = function(e) {
      paste("Error displaying variance decomposition:", e$message)
    })
  })

  # Model selection criteria output (AIC, BIC)
  output$warping_model_criteria <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see model criteria.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary
      method <- values$warping_results$method

      paste0(
        "Warping Method: ", method, "\n\n",
        "Model Selection Criteria:\n",
        "=========================\n",
        "AIC: ", sprintf("%.2f", stats$AIC), "\n",
        "BIC: ", sprintf("%.2f", stats$BIC), "\n",
        "Log-Likelihood: ", sprintf("%.2f", stats$log_likelihood), "\n",
        "Number of Parameters: ", stats$n_parameters, "\n\n",
        "Note: Lower AIC/BIC values indicate better model fit.\n",
        "Compare these values across different warping methods\n",
        "to select the optimal alignment approach."
      )
    }, error = function(e) {
      paste("Error displaying model criteria:", e$message)
    })
  })

  # Per-subject fit statistics table
  output$warping_fit_per_subject <- renderDT({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return(NULL)
    }

    tryCatch({
      per_subj <- values$warping_results$fit_statistics$per_subject

      # Add group labels if available
      if(!is.null(values$group_labels) && length(values$group_labels) == nrow(per_subj)) {
        per_subj$Group <- values$group_labels
        # Reorder columns to put Group after Subject
        cols <- c("Subject", "Group", setdiff(names(per_subj), c("Subject", "Group")))
        per_subj <- per_subj[, cols]
      }

      datatable(per_subj,
                options = list(
                  pageLength = 10,
                  scrollX = TRUE,
                  columnDefs = list(
                    list(className = 'dt-center', targets = '_all')
                  )
                ),
                rownames = FALSE) %>%
        formatStyle("R_squared",
                    backgroundColor = styleInterval(c(0.5, 0.7, 0.9),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("RMSE",
                    backgroundColor = styleInterval(c(0.05, 0.1, 0.2),
                                                    c("#99ff99", "#ccffcc", "#ffffcc", "#ffcccc"))) %>%
        formatStyle("Correlation",
                    backgroundColor = styleInterval(c(0.7, 0.85, 0.95),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("Warp_Amplitude",
                    backgroundColor = styleInterval(c(0.01, 0.05, 0.1),
                                                    c("white", "#ffffcc", "#ffcccc", "#ff9999")))

    }, error = function(e) {
      cat("Per-subject fit statistics table error:", e$message, "\n")
      datatable(data.frame(Message = "No per-subject statistics available"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })

  # Download handler for fit statistics CSV
  output$download_warping_fit_stats <- downloadHandler(
    filename = function() {
      paste0("warping_fit_statistics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        if(!is.null(values$warping_results) && !is.null(values$warping_results$fit_statistics)) {
          per_subj <- values$warping_results$fit_statistics$per_subject
          summary_stats <- values$warping_results$fit_statistics$summary

          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == nrow(per_subj)) {
            per_subj$Group <- values$group_labels
          }

          # Create summary row
          summary_row <- data.frame(
            Subject = "AVERAGE",
            R_squared = summary_stats$mean_r_squared,
            RMSE = summary_stats$mean_rmse,
            Correlation = summary_stats$mean_correlation,
            MAE = summary_stats$mean_mae,
            Orig_Variance = summary_stats$total_orig_variance,
            Amp_Variance = summary_stats$total_amp_variance,
            Phase_Variance = summary_stats$total_phase_variance,
            Full_Distance = summary_stats$mean_full_distance,
            Elastic_Amp_Dist = summary_stats$mean_elastic_amp_dist,
            Elastic_Phase_Dist = summary_stats$mean_elastic_phase_dist,
            Warp_Amplitude = summary_stats$mean_warp_amplitude,
            Warp_Velocity_Var = NA
          )

          # Add SD row
          sd_row <- data.frame(
            Subject = "SD",
            R_squared = summary_stats$sd_r_squared,
            RMSE = summary_stats$sd_rmse,
            Correlation = summary_stats$sd_correlation,
            MAE = summary_stats$sd_mae,
            Orig_Variance = NA,
            Amp_Variance = NA,
            Phase_Variance = NA,
            Full_Distance = summary_stats$sd_full_distance,
            Elastic_Amp_Dist = NA,
            Elastic_Phase_Dist = NA,
            Warp_Amplitude = summary_stats$sd_warp_amplitude,
            Warp_Velocity_Var = NA
          )

          # Add metadata
          meta_rows <- data.frame(
            Subject = c("---", "AIC", "BIC", "Variance_Explained_%", "Method"),
            R_squared = c(NA, summary_stats$AIC, summary_stats$BIC,
                          summary_stats$variance_explained_by_warping * 100,
                          NA),
            RMSE = NA, Correlation = NA, MAE = NA,
            Orig_Variance = NA, Amp_Variance = NA, Phase_Variance = NA,
            Full_Distance = NA, Elastic_Amp_Dist = NA, Elastic_Phase_Dist = NA,
            Warp_Amplitude = NA, Warp_Velocity_Var = NA
          )
          meta_rows$Subject[5] <- values$warping_results$method

          # Combine all
          output_df <- rbind(per_subj[, names(summary_row)], summary_row, sd_row, meta_rows)

          write.csv(output_df, file, row.names = FALSE)
        } else {
          write.csv(data.frame(Message = "No fit statistics available"), file, row.names = FALSE)
        }
      }, error = function(e) {
        write.csv(data.frame(Error = e$message), file, row.names = FALSE)
      })
    }
  )

  # Warping download buttons - FIXED
  output$download_warping_plot <- downloadHandler(
    filename = function() {
      paste0("warping_functions_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        png(file, width = 800, height = 600)
        
        if(!is.null(values$warping_results) && !is.null(values$warping_results$warp_functions)) {
          warp_functions <- values$warping_results$warp_functions
          time_points <- values$warping_results$time_points
          
          matplot(time_points, warp_functions[,1:min(30, ncol(warp_functions))], 
                  type = "l", col = rgb(0.4, 0.4, 0.8, 0.3), lty = 1,
                  xlab = "Original Time", ylab = "Warped Time",
                  main = "Time Warping Functions")
          lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
          legend("topleft", legend = c("Individual", "Identity"), 
                 col = c("lightblue", "red"), lty = c(1, 2), lwd = c(1, 2))
        } else {
          plot(1, type = "n", xlab = "", ylab = "", 
               main = "No warping functions available")
        }
        
        dev.off()
      }, error = function(e) {
        png(file)
        plot(1, main = "Error generating plot")
        dev.off()
      })
    }
  )
  
  output$download_alignment_plot <- downloadHandler(
    filename = function() {
      paste0("alignment_comparison_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        png(file, width = 800, height = 600)
        
        if(!is.null(values$fd_obj)) {
          time_points <- seq(0, 1, length.out = 100)
          orig_curves <- eval.fd(time_points, values$fd_obj)
          
          par(mfrow = c(1, 2))
          
          # Original curves
          matplot(time_points, orig_curves[,1:min(30, ncol(orig_curves))], 
                  type = "l", col = rgb(1, 0.4, 0.4, 0.3), lty = 1,
                  xlab = "Time", ylab = "Value", main = "Original Curves")
          lines(time_points, rowMeans(orig_curves), col = "darkred", lwd = 3)
          
          # Aligned curves
          if(!is.null(values$warping_results) && !is.null(values$warping_results$registered_curves)) {
            aligned_curves <- values$warping_results$registered_curves
            matplot(time_points, aligned_curves[,1:min(30, ncol(aligned_curves))], 
                    type = "l", col = rgb(0.4, 0.4, 1, 0.3), lty = 1,
                    xlab = "Time", ylab = "Value", main = "Aligned Curves")
            lines(time_points, rowMeans(aligned_curves), col = "darkblue", lwd = 3)
          } else {
            plot(1, type = "n", main = "No aligned curves available")
          }
        } else {
          plot(1, type = "n", main = "No data available")
        }
        
        dev.off()
      }, error = function(e) {
        png(file)
        plot(1, main = "Error generating plot")
        dev.off()
      })
    }
  )

  # ============================================================================
  # WARPING FIT STATISTICS CALCULATION
  # ============================================================================
  # Calculates fit statistics comparing raw vs warped data per subject and averaged
  # Based on EFDA methodology: variance decomposition and elastic distances
  # ============================================================================

  calculate_warping_fit_statistics <- function(original_curves, registered_curves,
                                                warp_functions, time_points) {
    # Calculate comprehensive fit statistics for time warping alignment
    #
    # Args:
    #   original_curves: matrix (n_time x n_subjects) of original data
    #   registered_curves: matrix (n_time x n_subjects) of warped data
    #   warp_functions: matrix (n_time x n_subjects) of warping functions h(t)
    #   time_points: vector of time points
    #
    # Returns:
    #   List with per-subject and averaged fit statistics

    tryCatch({
      n_time <- nrow(original_curves)
      n_subjects <- ncol(original_curves)
      dt <- diff(time_points[1:2])  # Time step for integration

      # Initialize per-subject statistics vectors
      r_squared <- numeric(n_subjects)
      rmse <- numeric(n_subjects)
      correlation <- numeric(n_subjects)
      mae <- numeric(n_subjects)

      # Variance decomposition (EFDA-style)
      orig_variance <- numeric(n_subjects)
      amp_variance <- numeric(n_subjects)
      phase_variance <- numeric(n_subjects)

      # Elastic distances (EFDA-style)
      full_dist <- numeric(n_subjects)
      elastic_amp_dist <- numeric(n_subjects)
      elastic_phase_dist <- numeric(n_subjects)

      # Warping intensity metrics
      warp_amplitude <- numeric(n_subjects)
      warp_velocity_var <- numeric(n_subjects)  # Variance of warping derivative

      # Calculate mean curves
      orig_mean <- rowMeans(original_curves, na.rm = TRUE)
      reg_mean <- rowMeans(registered_curves, na.rm = TRUE)

      # Calculate SRVF (Square Root Velocity Function) for elastic analysis
      # q(t) = sign(f'(t)) * sqrt(|f'(t)|)
      calc_srvf <- function(f, t) {
        n <- length(f)
        df <- diff(f) / diff(t)
        df <- c(df[1], df)  # Pad to same length
        sign(df) * sqrt(abs(df))
      }

      # Calculate SRVF of registered mean
      reg_mean_srvf <- calc_srvf(reg_mean, time_points)

      for(i in 1:n_subjects) {
        orig_i <- original_curves[,i]
        reg_i <- registered_curves[,i]
        warp_i <- warp_functions[,i]

        # Handle any NA values
        valid_idx <- !is.na(orig_i) & !is.na(reg_i)
        if(sum(valid_idx) < 3) {
          r_squared[i] <- NA
          rmse[i] <- NA
          correlation[i] <- NA
          mae[i] <- NA
          next
        }

        # ---- Basic Fit Statistics ----
        # Using FUNCTIONAL R² via integrals (L²-norm based)
        # R² = 1 - ∫(f(t) - g(t))² dt / ∫(f(t) - f̄)² dt
        # where f = original curve, g = warped curve, f̄ = mean of f over domain
        # This measures: "How much of the total squared energy of f is captured by g?"

        # f̄ = mean of original curve over the domain (scalar)
        f_bar <- mean(orig_i[valid_idx])

        # Numerator: ∫(f - g)² dt (integrated squared difference between original and warped)
        ss_res_functional <- sum((orig_i[valid_idx] - reg_i[valid_idx])^2) * dt

        # Denominator: ∫(f - f̄)² dt (integrated squared deviation from mean)
        ss_tot_functional <- sum((orig_i[valid_idx] - f_bar)^2) * dt

        # Functional R² (always between 0 and 1 when g approximates f well)
        # Note: Can still be negative if warped curve is worse than using the mean
        # but we clamp to 0 as a floor since negative R² is not meaningful here
        r_squared[i] <- if(ss_tot_functional > 0) {
          max(0, 1 - ss_res_functional / ss_tot_functional)
        } else {
          NA
        }

        # Correlation between original and warped (information preservation)
        if(var(orig_i[valid_idx]) > 0 && var(reg_i[valid_idx]) > 0) {
          correlation[i] <- cor(orig_i[valid_idx], reg_i[valid_idx])
        } else {
          correlation[i] <- NA
        }

        # RMSE: Distance from warped curve to ORIGINAL curve (not group mean)
        rmse[i] <- sqrt(mean((orig_i[valid_idx] - reg_i[valid_idx])^2))

        # MAE: Mean absolute error from warped curve to ORIGINAL curve
        mae[i] <- mean(abs(orig_i[valid_idx] - reg_i[valid_idx]))

        # ---- Variance Decomposition (EFDA-style) ----

        # Original variance: integrated squared deviation from original mean
        orig_variance[i] <- sum((orig_i - orig_mean)^2) * dt

        # Amplitude variance: integrated squared deviation after alignment
        amp_variance[i] <- sum((reg_i - reg_mean)^2) * dt

        # Phase variance: variance captured by warping
        # Computed from warped mean evaluated at individual warp functions
        warped_mean_i <- approx(time_points, reg_mean, xout = warp_i, rule = 2)$y
        phase_variance[i] <- sum((warped_mean_i - reg_mean)^2) * dt

        # ---- Elastic Distances (EFDA-style) ----

        # Full distance: L2 distance between aligned curve and mean
        full_dist[i] <- sqrt(sum((reg_i - reg_mean)^2) * dt)

        # Elastic amplitude distance: L2 distance in SRVF space
        reg_i_srvf <- calc_srvf(reg_i, time_points)
        elastic_amp_dist[i] <- sqrt(sum((reg_i_srvf - reg_mean_srvf)^2) * dt)

        # Elastic phase distance: arc-cosine of integrated psi
        # psi = sqrt(d(warp)/dt) is the SRVF of the warping function
        warp_deriv <- c(diff(warp_i) / diff(time_points), 0)
        warp_deriv[warp_deriv < 0] <- 0  # Ensure non-negative (monotonic)
        psi_i <- sqrt(pmax(0, warp_deriv))
        psi_integral <- sum(psi_i) * dt
        elastic_phase_dist[i] <- acos(min(1, max(-1, psi_integral)))

        # ---- Warping Intensity Metrics ----

        # Warping amplitude: RMSE deviation from identity h(t) = t
        warp_amplitude[i] <- sqrt(mean((warp_i - time_points)^2))

        # Warping velocity variance: how variable is the warping speed?
        warp_velocity_var[i] <- var(warp_deriv, na.rm = TRUE)
      }

      # ---- Model Selection Criteria (adapted for functional data) ----
      # AIC and BIC based on residual variance after alignment

      n_obs <- n_time * n_subjects  # Total observations
      residual_var <- mean(rmse^2, na.rm = TRUE)  # Mean squared error

      # Number of parameters: depends on warping method
      # For linear shift: 1 parameter (shift) per subject
      # For parametric: 1 parameter (alpha) per subject
      # For landmark: k parameters (landmark positions) per subject
      # Estimate as 2 parameters per subject (location + scale)
      k_params <- 2 * n_subjects

      # Log-likelihood (assuming Gaussian errors)
      log_lik <- -n_obs/2 * (log(2 * pi) + log(residual_var) + 1)

      # AIC = -2*logLik + 2*k
      aic <- -2 * log_lik + 2 * k_params

      # BIC = -2*logLik + k*log(n)
      bic <- -2 * log_lik + k_params * log(n_obs)

      # ---- Variance Explained by Alignment ----
      # Similar to R-squared but at group level
      total_orig_var <- sum(orig_variance, na.rm = TRUE)
      total_amp_var <- sum(amp_variance, na.rm = TRUE)
      total_phase_var <- sum(phase_variance, na.rm = TRUE)

      # Proportion of variance explained by warping
      var_explained_by_warping <- if(total_orig_var > 0) {
        1 - total_amp_var / total_orig_var
      } else {
        NA
      }

      # Return comprehensive statistics
      return(list(
        # Per-subject statistics
        per_subject = data.frame(
          Subject = 1:n_subjects,
          R_squared = round(r_squared, 4),
          RMSE = round(rmse, 4),
          Correlation = round(correlation, 4),
          MAE = round(mae, 4),
          Orig_Variance = round(orig_variance, 4),
          Amp_Variance = round(amp_variance, 4),
          Phase_Variance = round(phase_variance, 4),
          Full_Distance = round(full_dist, 4),
          Elastic_Amp_Dist = round(elastic_amp_dist, 4),
          Elastic_Phase_Dist = round(elastic_phase_dist, 4),
          Warp_Amplitude = round(warp_amplitude, 4),
          Warp_Velocity_Var = round(warp_velocity_var, 6)
        ),

        # Averaged statistics
        summary = list(
          # Basic fit statistics (mean ± SD)
          mean_r_squared = mean(r_squared, na.rm = TRUE),
          sd_r_squared = sd(r_squared, na.rm = TRUE),
          mean_rmse = mean(rmse, na.rm = TRUE),
          sd_rmse = sd(rmse, na.rm = TRUE),
          mean_correlation = mean(correlation, na.rm = TRUE),
          sd_correlation = sd(correlation, na.rm = TRUE),
          mean_mae = mean(mae, na.rm = TRUE),
          sd_mae = sd(mae, na.rm = TRUE),

          # Variance decomposition
          total_orig_variance = total_orig_var,
          total_amp_variance = total_amp_var,
          total_phase_variance = total_phase_var,
          variance_explained_by_warping = var_explained_by_warping,

          # Elastic distances (averaged)
          mean_full_distance = mean(full_dist, na.rm = TRUE),
          sd_full_distance = sd(full_dist, na.rm = TRUE),
          mean_elastic_amp_dist = mean(elastic_amp_dist, na.rm = TRUE),
          mean_elastic_phase_dist = mean(elastic_phase_dist, na.rm = TRUE),

          # Warping intensity
          mean_warp_amplitude = mean(warp_amplitude, na.rm = TRUE),
          sd_warp_amplitude = sd(warp_amplitude, na.rm = TRUE),

          # Model selection criteria
          AIC = aic,
          BIC = bic,
          log_likelihood = log_lik,
          n_parameters = k_params
        ),

        n_subjects = n_subjects,
        n_valid = sum(!is.na(r_squared))
      ))

    }, error = function(e) {
      cat("Error in calculate_warping_fit_statistics:", e$message, "\n")
      return(NULL)
    })
  }

  # Warping functions with better error handling
  linear_shift_alignment <- function(fd_obj, periodic = FALSE, 
                                     allow_dilation = FALSE, 
                                     dilation_range = c(0.95, 1.05),
                                     reference = "mean",
                                     time_points = seq(0, 1, length.out = 100)) {
    
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      # Evaluate curves
      curves <- eval.fd(time_points, fd_obj)
      
      # Validate curves
      if(is.null(curves) || ncol(curves) == 0) {
        stop("No valid curves to align")
      }
      
      # Get reference curve
      if(reference == "mean") {
        ref_curve <- rowMeans(curves)
      } else if(reference == "median") {
        ref_curve <- apply(curves, 1, median)
      } else {
        ref_curve <- curves[,1]
      }
      
      # Initialize warping
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      shifts <- numeric(n_curves)
      
      # Perform alignment for each curve
      for(i in 1:n_curves) {
        # Find best shift using cross-correlation
        if(periodic) {
          # Circular cross-correlation
          fft_curve <- fft(curves[,i] - mean(curves[,i]))
          fft_ref <- fft(ref_curve - mean(ref_curve))
          cross_corr <- Re(fft(Conj(fft_ref) * fft_curve, inverse = TRUE)) / n_time
          max_idx <- which.max(cross_corr)
          shift_idx <- if(max_idx > n_time/2) max_idx - n_time else max_idx
          shifts[i] <- -shift_idx / n_time
        } else {
          # Standard cross-correlation
          ccf_result <- ccf(curves[,i], ref_curve, lag.max = floor(n_time/4), 
                            plot = FALSE, na.action = na.pass)
          if(!is.null(ccf_result$acf) && length(ccf_result$acf) > 0) {
            best_lag <- ccf_result$lag[which.max(ccf_result$acf)]
            shifts[i] <- best_lag / n_time * 0.1  # Scale down shift
          } else {
            shifts[i] <- 0
          }
        }
        
        # Create warping function with some variation
        base_warp <- time_points - shifts[i] * 0.5
        
        # Add slight S-curve for visualization
        distortion <- sin(pi * time_points) * runif(1, -0.03, 0.03)
        warp_functions[,i] <- pmin(1, pmax(0, base_warp + distortion))
        
        # Ensure endpoints are fixed
        warp_functions[1,i] <- 0
        warp_functions[n_time,i] <- 1
        
        # Apply warping
        if(abs(shifts[i]) > 0.001) {
          # Interpolate curve at warped time points
          registered_curves[,i] <- approx(time_points, curves[,i], 
                                          xout = warp_functions[,i], 
                                          rule = 2)$y
        } else {
          registered_curves[,i] <- curves[,i]
        }
      }
      
      # Create fd objects for output
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      # Create warping function fd objects
      warp_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10)
      warp_smooth <- smooth.basis(time_points, warp_functions, warp_basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        shifts = shifts,
        method = "linear_shift",
        time_points = time_points
      ))
      
    }, error = function(e) {
      cat("Error in linear_shift_alignment:", e$message, "\n")
      # Return identity warping as fallback
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      return(list(
        regfd = fd_obj,
        registered_curves = eval.fd(time_points, fd_obj),
        warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
        shifts = rep(0, n_curves),
        method = "identity",
        time_points = time_points
      ))
    })
  }
  
  # Parametric alignment function
  parametric_alignment <- function(fd_obj, family = "power", 
                                   param_range = c(0.5, 2), 
                                   symmetric = FALSE,
                                   time_points = seq(0, 1, length.out = 100)) {
    
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      curves <- eval.fd(time_points, fd_obj)
      mean_curve <- rowMeans(curves)
      
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      alpha_values <- numeric(n_curves)
      
      # Define warping function
      warp_func <- function(t, alpha) {
        switch(family,
               "power" = pmin(1, pmax(0, t^alpha)),
               "exponential" = {
                 if(abs(alpha - 1) < 0.001) t
                 else pmin(1, pmax(0, (exp(alpha * t) - 1) / (exp(alpha) - 1)))
               },
               "quadratic" = pmin(1, pmax(0, alpha * t^2 + (1 - alpha) * t)),
               "logistic" = {
                 L <- function(x) 1 / (1 + exp(-alpha * (x - 0.5)))
                 L0 <- L(0)
                 L1 <- L(1)
                 pmin(1, pmax(0, (L(t) - L0) / (L1 - L0)))
               },
               t
        )
      }
      
      # Optimize warping for each curve
      for(i in 1:n_curves) {
        # Objective function
        objective <- function(alpha) {
          warped_time <- warp_func(time_points, alpha)
          warped_curve <- approx(time_points, curves[,i], xout = warped_time, rule = 2)$y
          sum((warped_curve - mean_curve)^2, na.rm = TRUE)
        }
        
        # Optimize
        result <- optimize(objective, interval = param_range, tol = 1e-4)
        alpha_values[i] <- result$minimum
        
        # Apply warping
        warp_functions[,i] <- warp_func(time_points, alpha_values[i])
        registered_curves[,i] <- approx(time_points, curves[,i], 
                                        xout = warp_functions[,i], rule = 2)$y
      }
      
      # Create fd objects
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        alpha_values = alpha_values,
        family = family,
        method = "parametric",
        time_points = time_points
      ))
      
    }, error = function(e) {
      cat("Error in parametric_alignment:", e$message, "\n")
      # Return identity warping
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      return(list(
        regfd = fd_obj,
        registered_curves = eval.fd(time_points, fd_obj),
        warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
        alpha_values = rep(1, n_curves),
        method = "identity",
        time_points = time_points
      ))
    })
  }
  
  # Simple landmark alignment
  landmark_alignment_simple <- function(fd_obj, landmarks, time_points) {
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      curves <- eval.fd(time_points, fd_obj)
      
      # If landmarks are provided, use them for alignment
      if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
        landmark_times <- values$landmark_points$x
        n_landmarks <- length(landmark_times)
        
        cat("Using", n_landmarks, "landmarks for alignment\n")
        
        # Find corresponding landmark points in each curve
        warp_functions <- matrix(NA, n_time, n_curves)
        registered_curves <- matrix(NA, n_time, n_curves)
        
        for(i in 1:n_curves) {
          # For each curve, find the actual landmarks (peaks/valleys near the specified times)
          curve_landmarks <- numeric(n_landmarks)
          
          for(j in 1:n_landmarks) {
            # Find local extremum near the landmark time
            search_window <- which(abs(time_points - landmark_times[j]) < 0.1)
            if(length(search_window) > 0) {
              local_values <- curves[search_window, i]
              # Find the peak or valley
              if(j %% 2 == 1) {
                # Look for peak for odd landmarks
                curve_landmarks[j] <- time_points[search_window[which.max(local_values)]]
              } else {
                # Look for valley for even landmarks
                curve_landmarks[j] <- time_points[search_window[which.min(local_values)]]
              }
            } else {
              curve_landmarks[j] <- landmark_times[j]
            }
          }
          
          # Create warping function using piecewise linear interpolation
          # Add boundary points
          all_landmark_times <- c(0, landmark_times, 1)
          all_curve_landmarks <- c(0, curve_landmarks, 1)
          
          # Interpolate warping function
          warp_functions[,i] <- approx(all_landmark_times, all_curve_landmarks, 
                                       xout = time_points, rule = 2)$y
          
          # Apply warping
          registered_curves[,i] <- approx(time_points, curves[,i], 
                                          xout = warp_functions[,i], rule = 2)$y
        }
      } else {
        # No landmarks provided, use automatic detection
        cat("No manual landmarks provided, using automatic landmark detection\n")
        
        # Simple automatic landmark detection: find common peaks
        mean_curve <- rowMeans(curves)
        
        # Find peaks in mean curve
        peaks <- which(diff(sign(diff(mean_curve))) == -2) + 1
        valleys <- which(diff(sign(diff(mean_curve))) == 2) + 1
        
        # Select up to 3 most prominent landmarks
        if(length(peaks) > 0 || length(valleys) > 0) {
          all_extrema <- sort(c(peaks, valleys))
          if(length(all_extrema) > 3) {
            # Select based on prominence
            prominence <- abs(mean_curve[all_extrema] - mean(mean_curve))
            all_extrema <- all_extrema[order(prominence, decreasing = TRUE)[1:3]]
          }
          landmark_times <- time_points[all_extrema]
        } else {
          # Default landmarks at quartiles
          landmark_times <- c(0.25, 0.5, 0.75)
        }
        
        # Apply landmark registration
        warp_functions <- matrix(NA, n_time, n_curves)
        registered_curves <- matrix(NA, n_time, n_curves)
        
        for(i in 1:n_curves) {
          # Simple identity warping with slight variation
          distortion <- sin(2*pi*time_points) * runif(1, -0.02, 0.02)
          warp_functions[,i] <- pmin(1, pmax(0, time_points + distortion))
          warp_functions[1,i] <- 0
          warp_functions[n_time,i] <- 1
          
          registered_curves[,i] <- curves[,i]
        }
      }
      
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        method = "landmark",
        time_points = time_points,
        landmarks_used = if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) 
          values$landmark_points$x else NULL
      ))
      
    }, error = function(e) {
      cat("Error in landmark_alignment:", e$message, "\n")
      return(NULL)
    })
  }

  # All remaining outputs (pairwise, warping, landmark, export) - these are already included
  # but I'll ensure warping_plot and alignment_comparison_plot are here
  
  output$warping_plot <- renderPlotly({
    if(is.null(values$warping_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run time-warped PCA first"))
    }
    
    tryCatch({
      warp_results <- values$warping_results
      
      if(is.null(warp_results$warp_functions)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "No warping functions available"))
      }
      
      warp_functions <- warp_results$warp_functions
      n_time <- nrow(warp_functions)
      n_curves <- ncol(warp_functions)
      
      time_points <- if(!is.null(warp_results$time_points)) {
        warp_results$time_points
      } else {
        seq(0, 1, length.out = n_time)
      }
      hover_times <- hover_time_labels(time_points)

      n_display <- if(!is.null(input$n_curves_display)) {
        min(input$n_curves_display, n_curves)
      } else {
        min(30, n_curves)
      }
      
      p <- plot_ly(type = 'scatter', mode = 'lines')
      
      for(i in 1:n_display) {
        p <- p %>% add_trace(
          x = time_points,
          y = warp_functions[,i],
          type = 'scatter',
          mode = 'lines',
          name = paste("Subject", i),
          line = list(color = 'rgba(100, 100, 200, 0.3)', width = 1),
          showlegend = FALSE,
          hovertemplate = paste("Subject", i, "<br>t: %{customdata}<br>h(t): %{y:.3f}<extra></extra>"),
          customdata = hover_times
        )
      }
      
      p <- p %>% add_trace(
        x = c(0, 1), 
        y = c(0, 1),
        type = 'scatter',
        mode = 'lines',
        name = "Identity (no warping)",
        line = list(color = 'red', width = 2, dash = 'dash'),
        hovertemplate = "Identity line<extra></extra>"
      )
      
      if(n_display > 0) {
        mean_warp <- rowMeans(warp_functions[,1:n_display, drop = FALSE])
        p <- p %>% add_trace(
          x = time_points,
          y = mean_warp,
          type = 'scatter',
          mode = 'lines',
          name = "Mean warping",
          line = list(color = 'black', width = 3),
          hovertemplate = "Mean warping<br>t: %{customdata}<br>h(t): %{y:.3f}<extra></extra>",
          customdata = hover_times
        )
      }
      
      p %>% layout(
        title = "Time Warping Functions",
        xaxis = list(title = "Original Time (t)", range = c(0, 1)),
        yaxis = list(title = "Warped Time h(t)", range = c(0, 1)),
        hovermode = 'x'
      )
      
    }, error = function(e) {
      cat("Warping plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  output$alignment_comparison_plot <- renderPlotly({
    if(is.null(values$fd_obj)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "No data available"))
    }
    
    tryCatch({
      n_time <- if(!is.null(values$warping_results) && !is.null(values$warping_results$time_points)) {
        length(values$warping_results$time_points)
      } else {
        100
      }
      
      time_points <- if(!is.null(values$warping_results) && !is.null(values$warping_results$time_points)) {
        values$warping_results$time_points
      } else {
        seq(0, 1, length.out = n_time)
      }
      hover_times <- hover_time_labels(time_points)

      n_display <- if(!is.null(input$n_curves_display)) {
        min(input$n_curves_display, ncol(values$fd_obj$coefs))
      } else {
        min(30, ncol(values$fd_obj$coefs))
      }
      
      orig_curves <- eval.fd(time_points, values$fd_obj)
      n_show <- min(n_display, ncol(orig_curves))
      
      p <- plot_ly(type = 'scatter', mode = 'lines')
      
      for(i in 1:n_show) {
        p <- p %>% add_trace(
          x = time_points,
          y = orig_curves[,i],
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'rgba(255, 100, 100, 0.2)', width = 1),
          showlegend = (i == 1),
          legendgroup = "original",
          name = "Original curves",
          hoverinfo = 'skip'
        )
      }
      
      orig_mean <- rowMeans(orig_curves[,1:n_show, drop = FALSE])
      p <- p %>% add_trace(
        x = time_points,
        y = orig_mean,
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'darkred', width = 3),
        name = "Original mean",
        hovertemplate = "Original mean<br>Time: %{customdata}<br>Value: %{y:.3f}<extra></extra>",
        customdata = hover_times
      )
      
      if(!is.null(values$warping_results) && !is.null(values$warping_results$registered_curves)) {
        aligned_curves <- values$warping_results$registered_curves
        
        if(ncol(aligned_curves) >= n_show && nrow(aligned_curves) == length(time_points)) {
          for(i in 1:n_show) {
            p <- p %>% add_trace(
              x = time_points,
              y = aligned_curves[,i],
              type = 'scatter',
              mode = 'lines',
              line = list(color = 'rgba(100, 100, 255, 0.2)', width = 1),
              showlegend = (i == 1),
              legendgroup = "aligned",
              name = "Aligned curves",
              hoverinfo = 'skip'
            )
          }
          
          aligned_mean <- rowMeans(aligned_curves[,1:n_show, drop = FALSE])
          p <- p %>% add_trace(
            x = time_points,
            y = aligned_mean,
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'darkblue', width = 3),
            name = "Aligned mean",
            hovertemplate = "Aligned mean<br>Time: %{customdata}<br>Value: %{y:.3f}<extra></extra>",
            customdata = hover_times
          )
        }
      }
      
      p <- p %>% layout(
        title = "Original vs Aligned Curves",
        yaxis = list(title = "Value"),
        hovermode = 'x'
      )
      p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_settings))
      p
      
    }, error = function(e) {
      cat("Alignment plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Landmark plot
  output$landmark_plot <- renderPlotly({
    req(values$data)
    
    tryCatch({
      n_time <- ncol(values$data)
      time_points <- seq(0, 1, length.out = n_time)
      
      # Determine what to plot
      # MERGED APP: NULL-safe ordering (see tools/port_fck.py). Neither input
      # is created by any UI, so both are NULL and the original test raised
      # "invalid length zero argument" instead of drawing the mean curve.
      if(is.null(input$selected_subject) || is.null(input$landmark_target) ||
         input$landmark_target == "mean") {
        # Plot mean curve
        if(!is.null(values$fd_obj)) {
          mean_curve <- rowMeans(eval.fd(time_points, values$fd_obj))
        } else {
          mean_curve <- colMeans(values$data)
        }
        
        p <- plot_ly(source = "landmark_source") %>%
          add_trace(x = time_points, y = mean_curve, 
                    type = 'scatter', mode = 'lines',
                    name = 'Mean curve',
                    line = list(color = 'blue', width = 2))
        
        title_text <- "Click to add landmarks on mean curve"
      } else {
        # Plot individual subject
        subj_idx <- as.numeric(input$selected_subject)
        subj_curve <- values$data[subj_idx,]
        
        p <- plot_ly(source = "landmark_source") %>%
          add_trace(x = time_points, y = subj_curve,
                    type = 'scatter', mode = 'lines',
                    name = paste('Subject', subj_idx),
                    line = list(color = 'blue', width = 2))
        
        title_text <- paste("Click to add landmarks for Subject", subj_idx)
      }
      
      # Add existing landmarks if any
      if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
        p <- p %>% add_trace(x = values$landmark_points$x,
                             y = values$landmark_points$y,
                             type = 'scatter', mode = 'markers',
                             marker = list(color = 'red', size = 10),
                             name = 'Landmarks')
      }
      
      p <- p %>% 
        layout(title = title_text,
               yaxis = list(title = "Value"),
               hovermode = 'closest',
               clickmode = 'event+select') %>%
        config(displayModeBar = FALSE)
      
      # Apply time label formatting (uses existing helper function)
      p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_settings))
      p
      
    }, error = function(e) {
      plot_ly() %>% layout(title = paste("Error:", e$message))
    })
  })
  
  # Handle landmark clicks
  observeEvent(event_data("plotly_click", source = "landmark_source"), {
    click <- event_data("plotly_click", source = "landmark_source")
    if(!is.null(click)) {
      # Only add landmark if it's a click on the curve, not on existing landmarks
      if(is.null(click$curveNumber) || click$curveNumber == 0) {
        new_point <- data.frame(x = click$x, y = click$y)
        values$landmark_points <- rbind(values$landmark_points, new_point)
        
        showNotification(paste("Landmark added at t =", round(click$x, 3)), 
                         duration = 2, type = "message")
      }
    }
  })
  
  # Start landmark selection
  observeEvent(input$start_landmark, {
    values$landmark_points <- data.frame(x = numeric(), y = numeric())
    showNotification("Landmarks cleared. Click on the plot to add new landmarks.", 
                     duration = 3, type = "message")
  })
  
  # Clear landmarks
  observeEvent(input$clear_landmarks, {
    values$landmark_points <- data.frame(x = numeric(), y = numeric())
    showNotification("Landmarks cleared", duration = 2)
  })
  
  # Add landmark info display
  output$landmark_info <- renderPrint({
    if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
      cat("Current landmarks:\n")
      for(i in 1:nrow(values$landmark_points)) {
        cat(sprintf("  Landmark %d: t = %.3f\n", i, values$landmark_points$x[i]))
      }
    } else {
      cat("No landmarks defined. Click on the plot to add landmarks.")
    }
  })
