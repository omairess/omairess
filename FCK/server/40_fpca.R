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
      
      # One palette app-wide, keyed by group NAME so filtering never repaints.
      colors <- unname(fck_group_colors(groups)[groups])
      
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
          # ==========================================================================
  # fck_landmark_warp(ref, own, time_points)
  #
  # The one place a landmark warp is built. Returns h with the module-wide
  # contract (P5.6):
  #
  #     h maps REGISTERED time to ORIGINAL time, and registration is always
  #         registered_i <- approx(time_points, original_i, xout = h_i)$y
  #
  # so h(t) answers "which point of the original curve belongs at registered
  # time t". Every method in this file returns h in that direction. Before P5.6
  # the two landmark branches disagreed: the manual branch used this contract,
  # while the automatic branch built h from (own -> ref) and then applied it as
  # approx(h, curve, xout = t), which is the INVERSE map. The registered curves
  # came out plausible either way, because the second form inverts the map
  # while interpolating -- but warp_functions then held two different objects
  # depending on which branch ran, and every downstream statistic that reads
  # warp_functions - time_points (the phase metrics, the warping plots, the
  # per-subject table) was comparing incomparable things across methods.
  #
  # AUDIT (P5.5): the manual branch also had no monotonicity requirement at
  # all. It took whatever the peak/valley search returned, pasted 0 and 1 on
  # the ends and interpolated. On a noisy curve the detected landmarks can
  # cross -- the search alternates peaks and valleys in fixed windows, and
  # nothing stopped landmark 3 landing before landmark 2 -- and a crossed pair
  # makes h non-monotone, which is not a reparameterisation of time. It is now
  # a hard precondition: crossed or duplicated landmarks are REJECTED and the
  # subject falls back to the identity, reported, rather than being registered
  # with a fold in it.
  #
  # Returns NULL when the landmarks do not define a valid warp.
  # ==========================================================================
  fck_landmark_warp <- function(ref, own, time_points) {
    if (length(ref) != length(own) || !length(ref)) return(NULL)
    if (any(!is.finite(ref)) || any(!is.finite(own))) return(NULL)

    t0 <- time_points[1]; t1 <- time_points[length(time_points)]
    # order both by the REFERENCE position, so a reference given out of order
    # does not silently reorder the correspondence
    o   <- order(ref)
    ref <- ref[o]; own <- own[o]

    # strictly increasing, strictly interior, on BOTH sides
    tol <- 1e-6 * max(1, t1 - t0)
    if (any(diff(ref) <= tol) || any(diff(own) <= tol)) return(NULL)
    if (ref[1] <= t0 + tol || ref[length(ref)] >= t1 - tol) return(NULL)
    if (own[1] <= t0 + tol || own[length(own)] >= t1 - tol) return(NULL)

    kx <- c(t0, ref, t1)     # registered time
    ky <- c(t0, own, t1)     # original time
    h  <- stats::approx(kx, ky, xout = time_points, rule = 2)$y

    # the invariants this function exists to guarantee
    if (!isTRUE(all.equal(h[1], t0, tolerance = 1e-8))) return(NULL)
    if (!isTRUE(all.equal(h[length(h)], t1, tolerance = 1e-8))) return(NULL)
    if (any(diff(h) <= 0)) return(NULL)
    h
  }

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
          # P5.4: the statistics need the geometry. A translation and an
          # endpoint-preserving reparameterisation do not share a phase metric.
          warp_fd$fit_statistics <- calculate_warping_fit_statistics(
            original_curves = original_curves,
            registered_curves = warp_fd$registered_curves,
            warp_functions = warp_fd$warp_functions,
            time_points = time_grid,
            method = warp_fd$method,
            shifts = warp_fd$shifts,
            periodic = identical(warp_fd$boundary, "periodic wrap")
          )
          if(!is.null(warp_fd$fit_statistics)) {
            cat("Warping fit statistics calculated successfully\n")
            cat("  Mean R²:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_r_squared), "\n")
            cat("  Mean RMSE:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_rmse), "\n")
            cat("  Between-curve dispersion reduced by:",
                sprintf("%.2f%%", warp_fd$fit_statistics$summary$dispersion_reduction * 100), "\n")
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
  
  # AUDIT: the "Components to show" slider on the effect-of-scores plot went to
  # 10 regardless of how many components the PCA had actually retained, so
  # asking for 6 silently drew 3 and the title read "3 of 3". A display control
  # must not promise more than the analysis produced: the slider's ceiling now
  # follows nharm, and the plot says so when the two disagree.
  observeEvent(values$pca_results, {
    pr <- values$pca_results
    if(is.null(pr) || is.null(pr$scores)) return()
    # P6.3: length() of an fd object is 3, always. Use the coefficient columns.
    n_avail <- min(ncol(pr$scores), length(pr$values), fck_n_harmonics(pr))
    if(!is.finite(n_avail) || n_avail < 1) return()
    cur <- suppressWarnings(as.integer(isolate(input$effect_n_comp) %||% 3))
    if(!is.finite(cur) || cur < 1) cur <- 3
    # P6.4: updating only `max` leaves the client recomputing its tick
    # positions from the OLD range, which is why the slider drew ~19 ticks
    # labelled 1,1,1,2,2,2,2,3,3,3 on a 1-to-3 range. Send the whole
    # specification -- min, max, step -- so the widget is rebuilt consistently.
    updateSliderInput(session, "effect_n_comp",
                      min = 1, max = max(1L, n_avail), step = 1,
                      value = min(cur, n_avail))
  }, ignoreNULL = TRUE)

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
          # P6.3: this printed `1:min(3, length(varprop))` -- a hard-coded three,
          # so a user who asked for five components was shown three and told
          # nothing. Show every component the PCA retained, with the running
          # total, which is what you need to decide how many to keep.
          nk <- min(fck_n_harmonics(pca_res), length(pca_res$varprop))
          if (nk < 1) nk <- length(pca_res$varprop)
          cat("Variance explained:\n")
          cs <- 0
          for(i in seq_len(nk)) {
            cs <- cs + pca_res$varprop[i]
            cat(sprintf("  PC%d: %6.2f%%   (cumulative %6.2f%%)\n",
                        i, pca_res$varprop[i] * 100, cs * 100))
          }
          cat(sprintf("  ---\n  %d component%s retained, %.2f%% of the total variance.\n",
                      nk, if (nk == 1) "" else "s", cs * 100))
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
      
      colors <- fck_component_colors(max(1, min(ncol(pca_res$scores), fck_n_harmonics(pca_res))))
      # P6.3: was min(ncol(scores), 5, length(harmonics)) -- the last term is
      # always 3, so this plot never drew more than three components, and the
      # hard 5 capped it again for anyone who asked for more.
      n_comp <- min(ncol(pca_res$scores), fck_n_harmonics(pca_res))
      
      for(i in 1:n_comp) {
        if(i <= fck_n_harmonics(pca_res)) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          p <- p %>% add_trace(x = time_points,
                               y = as.vector(loading_vals),
                               type = 'scatter',
                               mode = 'lines',
                               name = paste("PC", i),
                               line = list(color = colors[i], width = 2))
        }
      }
      
      # A centred title with plotly's default top margin lands on the topmost
      # trace. Reserve the space rather than shrinking the title.
      p <- p %>% layout(title = list(text = "Principal component loadings",
                                     x = 0.5, xanchor = "center"),
                        yaxis = list(title = "Loading"),
                        margin = list(t = 60))

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
                  line = list(color = FCK_EMPHASIS),
                  marker = list(color = FCK_EMPHASIS))
      
      p %>% layout(title = list(text = "Variance explained", x = 0.5, xanchor = "center"),
                   margin = list(t = 60),
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
      
      # AUDIT: capped at three components with a red/blue/green palette. The cap
      # is now the user's (up to whatever the PCA retained) and the palette is
      # the app's shared one, so a component is the same colour here as in the
      # component-ANOVA figure. Past five, dash carries identity alongside hue --
      # see server/02b_helpers_palette.R.
      n_avail <- min(ncol(pca_res$scores), length(pca_res$values), fck_n_harmonics(pca_res))
      n_req <- suppressWarnings(as.integer(input$effect_n_comp %||% 3))
      if(!is.finite(n_req) || n_req < 1) n_req <- 3
      n_show <- min(n_req, n_avail)
      # solid/dash already carries the SIGN of the deviation (+2SD vs -2SD), so
      # it is not available to carry component identity as well. Hue is the only
      # channel left, and the validated ramp has five all-pairs-separable slots:
      # past that the note under the plot says so rather than pretending eight
      # components are distinguishable.
      colors <- fck_group_ramp(max(n_show, 1))
      
      for(i in 1:n_show) {
        if(i <= fck_n_harmonics(pca_res)) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          
          # Plus 2 SD
          plus_2sd <- as.vector(mean_vals + effect_mult * 2 * sqrt(pca_res$values[i]) * loading_vals)
          p <- p %>% add_trace(
            x = time_points,
            y = plus_2sd,
            type = 'scatter',
            mode = 'lines',
            name = paste("PC", i, "+2SD"),
            legendgroup = paste0("PC", i),
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
            legendgroup = paste0("PC", i),
            line = list(color = colors[i], dash = 'dash')
          )
        }
      }
      
      # If the request exceeds what the PCA retained, name the remedy on the
      # figure rather than silently drawing fewer curves than were asked for.
      short_note <- if(n_req > n_avail)
        sprintf("<br><sub>Only %d component%s extracted \u2014 raise \"Number of components to extract\" in PCA Settings and re-run to see more</sub>",
                n_avail, if(n_avail == 1) " was" else "s were") else ""
      p <- p %>% layout(
        title = list(text = paste0(sprintf("Effect of component scores (%d of %d component%s shown)",
                                           n_show, n_avail, if(n_avail == 1) "" else "s"),
                                   short_note),
                     x = 0.5, xanchor = "center"),
        yaxis = list(title = "Value"),
        legend = list(orientation = "h", y = -0.18),
        margin = list(t = if(nzchar(short_note)) 84 else 60))
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

      # AUDIT (P5.3): this panel was headed "EFDA Variance Decomposition" and
      # split the total into "amplitude" and "phase" percentages. The three
      # quantities were computed independently and do not sum, nothing here
      # establishes orthogonality, and naming a published methodology asserted
      # a property that had never been demonstrated. What the numbers actually
      # are: between-curve dispersion before and after registration, and a
      # separate descriptive measure of how far the template moves under each
      # subject's warp. Reported as that, with the ratio named for what it is.
      pre  <- stats$total_dispersion_pre
      post <- stats$total_dispersion_post
      g    <- stats$dispersion_reduction

      paste0(
        "Pre/post registration dispersion:\n",
        "=================================\n\n",
        "  V_pre  = sum_i integral (x_i(t) - xbar(t))^2 dt        : ",
        sprintf("%.4f", pre), "\n",
        "  V_post = sum_i integral (x_i(h_i(t)) - xbar_R(t))^2 dt : ",
        sprintf("%.4f", post),
        if (is.finite(pre) && pre > 0)
          paste0("  (", sprintf("%.1f%%", 100 * post / pre), " of V_pre)") else "",
        "\n\n",
        "  G = 1 - V_post / V_pre : ",
        if (is.null(g) || is.na(g)) "NA" else sprintf("%+.1f%%", 100 * g), "\n\n",
        "G is the relative reduction in between-curve dispersion after\n",
        "registration. A negative G means registration made the curves MORE\n",
        "dispersed than they were.\n\n",
        "THIS IS NOT AN AMPLITUDE/PHASE VARIANCE DECOMPOSITION. V_pre and\n",
        "V_post are two dispersions of the same curves, before and after; they\n",
        "are not orthogonal components of a total and they do not sum to one.\n",
        "Do not report G as 'variance explained by phase'.\n\n",
        "Template deformation (descriptive, not a component):\n",
        "  sum_i integral (xbar_R(h_i(t)) - xbar_R(t))^2 dt : ",
        sprintf("%.4f", stats$total_template_deformation), "\n",
        "  How far each subject's warp moves the registered template. Large\n",
        "  values mean the warps are doing a lot; they do not add to V_post.\n"
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

      # P5.2: this panel used to print AIC, BIC, a log-likelihood and a
      # parameter count, and told the user to pick a warping method by
      # comparing them. None of those quantities existed. What can honestly be
      # reported is how much between-curve dispersion the registration removed,
      # and how far it moved the time axis to do it.
      g   <- stats$dispersion_reduction
      geo <- stats$geometry
      paste0(
        "Warping Method: ", method,
        if (!is.null(stats$periodic) && isTRUE(stats$periodic)) "  (periodic)" else "", "\n\n",
        "Between-curve dispersion\n",
        "========================\n",
        "  V_pre  = sum_i integral (x_i - xbar)^2 dt      : ",
        sprintf("%.4g", stats$total_dispersion_pre), "\n",
        "  V_post = sum_i integral (x_i(h_i) - xbar_R)^2 dt: ",
        sprintf("%.4g", stats$total_dispersion_post), "\n",
        "  Relative reduction G = 1 - V_post/V_pre        : ",
        if (is.na(g)) "NA" else sprintf("%+.1f%%", 100 * g), "\n\n",
        "G is the proportion of between-curve dispersion removed by\n",
        "registration. It is NOT an amplitude/phase variance decomposition:\n",
        "nothing here establishes that the two are additive or orthogonal.\n",
        "A NEGATIVE G means registration made the curves MORE dispersed.\n\n",
        # P5.12: the amplitude-leakage signature. Least-squares registration on
        # a sharply peaked curve can absorb AMPLITUDE differences with a tiny
        # time shift -- near a peak, a 1% move in time changes the value a lot.
        # Measured on curves that differ ONLY in amplitude, with the logistic
        # family: G = 27.8% from warps averaging 0.0017 of the domain, with the
        # peak heights unchanged. A deviation-from-identity penalty does not
        # help (the offending warps are already near-identity), so the honest
        # response is to make the signature visible rather than to add a knob
        # that does not work.
        if (is.finite(g) && g > 0.15 &&
            is.finite(stats$mean_phase_displacement) &&
            stats$mean_phase_displacement < 0.02) {
          paste0(
            "WARNING -- possible amplitude leakage.\n",
            "  G is ", sprintf("%.0f%%", 100 * g), " but the warps move time by only ",
            sprintf("%.4f", stats$mean_phase_displacement), " of the domain\n",
            "  on average. A large dispersion reduction from a near-identity warp\n",
            "  usually means the criterion is absorbing AMPLITUDE differences, not\n",
            "  aligning phase -- near a peak, a small time move changes the value a\n",
            "  lot. Check the registered curves against the originals before\n",
            "  reporting G as evidence of phase variation.\n\n")
        } else "",
        "Time-axis displacement\n",
        "======================\n",
        if (identical(geo, "shift")) {
          paste0("  Mean |shift| (",
                 if (isTRUE(stats$periodic)) "circular" else "translation",
                 ", fraction of domain): ",
                 sprintf("%.4f", stats$mean_phase_displacement), "\n",
                 "  Fisher-Rao phase distance is not defined for a translation\n",
                 "  and is reported as NA (see P5.4).\n")
        } else {
          paste0("  RMS |h(t) - t|                 : ",
                 sprintf("%.4f", stats$mean_phase_displacement), "\n",
                 "  Fisher-Rao phase distance (mean): ",
                 sprintf("%.4f", stats$mean_elastic_phase_dist),
                 "  on ", stats$n_elastic_phase_valid, " of ",
                 length(values$warping_results$fit_statistics$per_subject$Subject),
                 " curves\n",
                 if (isTRUE(stats$n_invalid_warp > 0))
                   paste0("  ", stats$n_invalid_warp,
                          " warp(s) were not endpoint-preserving monotone maps and\n",
                          "  were excluded from the phase metric rather than clipped.\n")
                 else "")
        },
        "\n",
        "No AIC/BIC is reported. There is no likelihood for the observed\n",
        "curves under a candidate registration in this module, so there is no\n",
        "criterion to compare methods with. Choose a registration method from\n",
        "what you know about the data, not from a fit index."
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
            Dispersion_Pre = summary_stats$total_dispersion_pre,
            Dispersion_Post = summary_stats$total_dispersion_post,
            Template_Deformation = summary_stats$total_template_deformation,
            Full_Distance = summary_stats$mean_full_distance,
            Elastic_Amp_Dist = summary_stats$mean_elastic_amp_dist,
            Elastic_Phase_Dist = summary_stats$mean_elastic_phase_dist,
            Phase_Displacement = summary_stats$mean_phase_displacement,
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
            Dispersion_Pre = NA,
            Dispersion_Post = NA,
            Template_Deformation = NA,
            Full_Distance = summary_stats$sd_full_distance,
            Elastic_Amp_Dist = NA,
            Elastic_Phase_Dist = NA,
            Phase_Displacement = summary_stats$sd_phase_displacement,
            Warp_Amplitude = summary_stats$sd_warp_amplitude,
            Warp_Velocity_Var = NA
          )

          # Add metadata
          meta_rows <- data.frame(
            # P5.2: AIC and BIC removed -- they were not a likelihood.
            Subject = c("---", "Dispersion_Reduction_G_%", "Geometry", "Method"),
            R_squared = c(NA, summary_stats$dispersion_reduction * 100, NA, NA),
            RMSE = NA, Correlation = NA, MAE = NA,
            Dispersion_Pre = NA, Dispersion_Post = NA, Template_Deformation = NA,
            Full_Distance = NA, Elastic_Amp_Dist = NA, Elastic_Phase_Dist = NA,
            Phase_Displacement = NA,
            Warp_Amplitude = NA, Warp_Velocity_Var = NA
          )
          meta_rows$Subject[3] <- paste("Geometry:", summary_stats$geometry %||% "interval")
          meta_rows$Subject[4] <- paste("Method:", values$warping_results$method)

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
                  type = "l", col = fck_group_rgba(FCK_SERIES1, 0.30), lty = 1,
                  xlab = "Original Time", ylab = "Warped Time",
                  main = "Time Warping Functions")
          lines(c(0,1), c(0,1), col = FCK_NEUTRAL, lwd = 2, lty = 2)
          legend("topleft", legend = c("Individual", "Identity"), 
                 col = c(FCK_SERIES1, FCK_NEUTRAL), lty = c(1, 2), lwd = c(1, 2))
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

  # AUDIT (P5.2/P5.3/P5.4): this panel presented three things it had not
  # earned. See the notes at each site below. The function now needs to know
  # WHICH registration produced the warps, because a translation and an
  # endpoint-preserving reparameterisation are different geometries and the
  # same formula is not valid in both.
  calculate_warping_fit_statistics <- function(original_curves, registered_curves,
                                                warp_functions, time_points,
                                                method = NULL, shifts = NULL,
                                                periodic = FALSE) {
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

      # Pre/post registration dispersion (P5.3 -- not a variance decomposition)
      disp_pre  <- numeric(n_subjects)
      disp_post <- numeric(n_subjects)
      template_deformation <- numeric(n_subjects)

      # Elastic distances (SRVF)
      full_dist <- numeric(n_subjects)
      elastic_amp_dist <- numeric(n_subjects)
      elastic_phase_dist <- numeric(n_subjects)
      phase_displacement <- numeric(n_subjects)   # P5.4, in the right geometry
      n_invalid_warp <- 0L

      # P5.4: which geometry are we in? A translation is not an
      # endpoint-preserving reparameterisation and does not admit the same
      # phase metric.
      geom <- if (identical(method, "linear_shift")) "shift" else "interval"

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

        # ---- Between-curve dispersion, before and after registration -------
        #
        # AUDIT (P5.3): these three quantities were computed separately and
        # presented under the heading "Variance Decomposition (EFDA-style)",
        # with a "variance explained by warping" derived from the first two.
        # They are NOT a decomposition: nothing here establishes
        # V_total = V_amplitude + V_phase, the three are not orthogonal, and
        # they do not sum. Calling them a decomposition, and naming a published
        # methodology, claims a property that was never demonstrated.
        #
        # What they honestly are: the integrated squared deviation of each
        # curve from the sample mean, BEFORE registration and AFTER it, plus a
        # third descriptive quantity (how far the template moves when pushed
        # through this subject's warp). The first two are reported as pre- and
        # post-registration dispersion, and their ratio as a relative
        # reduction. The third is kept, renamed, and explicitly NOT presented
        # as a component of anything.
        disp_pre[i]  <- sum((orig_i - orig_mean)^2) * dt
        disp_post[i] <- sum((reg_i  - reg_mean)^2)  * dt

        warped_mean_i <- approx(time_points, reg_mean, xout = warp_i, rule = 2)$y
        template_deformation[i] <- sum((warped_mean_i - reg_mean)^2) * dt

        # ---- Elastic Distances (EFDA-style) ----

        # Full distance: L2 distance between aligned curve and mean
        full_dist[i] <- sqrt(sum((reg_i - reg_mean)^2) * dt)

        # Elastic amplitude distance: L2 distance in SRVF space
        reg_i_srvf <- calc_srvf(reg_i, time_points)
        elastic_amp_dist[i] <- sqrt(sum((reg_i_srvf - reg_mean_srvf)^2) * dt)

        # ---- Phase displacement, in the geometry that applies --------------
        #
        # AUDIT (P5.4): the Fisher-Rao phase distance from the identity,
        #     d(h, id) = arccos( integral sqrt(h\'(t)) dt ),
        # is valid for an ENDPOINT-PRESERVING monotone reparameterisation of
        # [0,1]. It was applied to every warp type this module produces, and
        # for two of them it is meaningless. Measured on a 100-point grid:
        #
        #   non-periodic shift h(t) = t - s:  h' = 1 everywhere, so psi = 1 and
        #     the distance is 0.000000 for EVERY s -- including a quarter of the
        #     domain. The metric is identically blind to translation.
        #   periodic shift (wrapped):  the wrap puts one large negative step in
        #     the difference quotient, which the line above set to zero. The
        #     result is 0.142254 for every s -- INCLUDING s = 0. An unshifted
        #     curve was reported as having the same nonzero phase distance as a
        #     quarter-cycle shift. That number is pure artefact of the
        #     discontinuity; setting a negative derivative to zero does not turn
        #     a rotation into a Fisher-Rao warp.
        #
        # So the phase summary is now computed per geometry, and the elastic
        # distance is reported only where it means something.
        warp_deriv <- c(diff(warp_i) / diff(time_points), 0)

        if (identical(geom, "shift")) {
          # A translation has no Fisher-Rao phase distance from the identity
          # that distinguishes it from the identity. Report the displacement.
          elastic_phase_dist[i] <- NA_real_
          s_i <- if (!is.null(shifts) && length(shifts) >= i) shifts[i] else
                 mean(time_points - warp_i, na.rm = TRUE)
          phase_displacement[i] <- if (isTRUE(periodic)) {
            # circular: wrap to [-span/2, span/2] and report the magnitude
            span <- diff(range(time_points))
            d <- ((s_i + span / 2) %% span) - span / 2
            abs(d)
          } else abs(s_i)
        } else {
          # Endpoint-preserving monotone warp: check that it IS one before
          # computing a metric that assumes it. A warp that fails these is
          # reported as NA rather than silently repaired by clipping.
          endpoints_ok <- abs(warp_i[1] - time_points[1]) < 1e-6 &&
                          abs(warp_i[length(warp_i)] - time_points[length(time_points)]) < 1e-6
          monotone_ok  <- all(diff(warp_i) > -1e-12)
          if (endpoints_ok && monotone_ok) {
            psi_i <- sqrt(pmax(0, warp_deriv))
            psi_integral <- sum(psi_i) * dt
            elastic_phase_dist[i] <- acos(min(1, max(-1, psi_integral)))
          } else {
            elastic_phase_dist[i] <- NA_real_
            n_invalid_warp <- n_invalid_warp + 1L
          }
          phase_displacement[i] <- sqrt(mean((warp_i - time_points)^2))
        }
        warp_deriv[warp_deriv < 0] <- 0   # only for the velocity summary below

        # ---- Warping Intensity Metrics ----

        # Warping amplitude: RMSE deviation from identity h(t) = t
        warp_amplitude[i] <- sqrt(mean((warp_i - time_points)^2))

        # Warping velocity variance: how variable is the warping speed?
        warp_velocity_var[i] <- var(warp_deriv, na.rm = TRUE)
      }

      # ---- AIC / BIC: REMOVED --------------------------------------------
      #
      # AUDIT (P5.2): this block computed
      #     residual_var <- mean(rmse^2)
      #     log_lik <- -n_obs/2 * (log(2*pi) + log(residual_var) + 1)
      #     AIC <- -2*log_lik + 2*k_params;  k_params <- 2 * n_subjects
      # and the tab told the user "Lower AIC/BIC values indicate better model
      # fit. Compare these values across different warping methods to select
      # the optimal alignment approach." Both halves are wrong.
      #
      # First, `rmse` here is the distance between each curve and its OWN
      # registered version -- how much the registration MOVED the curve, not
      # how well it aligned the curve with anything. Deforming the time axis is
      # exactly what registration is for, so a good registration can and should
      # have a large value; treating it as residual error means the criterion
      # rewards doing nothing. Second, there is no likelihood here: no
      # probability model for the observed functional data under a candidate
      # registration was ever written down, so -2 log L + 2k is a formula
      # applied to a number that is not a log-likelihood. Third, k_params was
      # hard-coded at 2 per subject regardless of method, while a shift fits
      # one parameter per subject, a parametric warp fits one, and a landmark
      # warp fits as many as there are landmarks -- so the penalty term did not
      # even distinguish the methods it was being used to compare.
      #
      # An AIC that cannot be computed is not replaced by a worse AIC. The
      # criterion is gone, and the panel reports the descriptive quantities
      # below, which mean what they say. Restoring model selection here needs a
      # probabilistic registration model with an actual likelihood; that is a
      # separate piece of work, not a relabelling.

      # ---- Pre/post registration dispersion (P5.3) ------------------------
      # V_pre  = sum_i integral (x_i(t)      - xbar(t))^2   dt
      # V_post = sum_i integral (x_i(h_i(t)) - xbar_R(t))^2 dt
      # G      = 1 - V_post / V_pre
      # G is the RELATIVE REDUCTION IN BETWEEN-CURVE DISPERSION after
      # registration. It is not "variance explained by warping" in the sense of
      # an additive phase/amplitude split, and it is not called that any more.
      total_disp_pre  <- sum(disp_pre,  na.rm = TRUE)
      total_disp_post <- sum(disp_post, na.rm = TRUE)
      total_template_deformation <- sum(template_deformation, na.rm = TRUE)

      dispersion_reduction <- if(total_disp_pre > 0) {
        1 - total_disp_post / total_disp_pre
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
          Dispersion_Pre = round(disp_pre, 4),
          Dispersion_Post = round(disp_post, 4),
          Template_Deformation = round(template_deformation, 4),
          Full_Distance = round(full_dist, 4),
          Elastic_Amp_Dist = round(elastic_amp_dist, 4),
          Elastic_Phase_Dist = round(elastic_phase_dist, 4),
          Phase_Displacement = round(phase_displacement, 4),
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

          # Pre/post registration dispersion (P5.3 -- NOT a decomposition)
          total_dispersion_pre  = total_disp_pre,
          total_dispersion_post = total_disp_post,
          total_template_deformation = total_template_deformation,
          dispersion_reduction  = dispersion_reduction,

          # Elastic distances (averaged)
          mean_full_distance = mean(full_dist, na.rm = TRUE),
          sd_full_distance = sd(full_dist, na.rm = TRUE),
          mean_elastic_amp_dist = mean(elastic_amp_dist, na.rm = TRUE),
          mean_elastic_phase_dist = mean(elastic_phase_dist, na.rm = TRUE),
          n_elastic_phase_valid = sum(!is.na(elastic_phase_dist)),
          n_invalid_warp = n_invalid_warp,

          # Phase displacement in the geometry that applies (P5.4)
          geometry = geom,
          periodic = isTRUE(periodic),
          mean_phase_displacement = mean(phase_displacement, na.rm = TRUE),
          sd_phase_displacement = sd(phase_displacement, na.rm = TRUE),

          # Warping intensity
          mean_warp_amplitude = mean(warp_amplitude, na.rm = TRUE),
          sd_warp_amplitude = sd(warp_amplitude, na.rm = TRUE)
          # AIC / BIC / log-likelihood / n_parameters removed at P5.2 -- see
          # the note above. They were not a likelihood.
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
      extrap_frac <- numeric(n_curves)   # P4.2: how much of each warp left the domain

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
            # AUDIT (P0.8): was  best_lag / n_time * 0.1  -- the estimated lag
            # was multiplied by 0.1 here and by a further 0.5 below, so the warp
            # carried 5% of the shift it had just measured. Neither constant was
            # justified anywhere. The measured lag is now used as measured.
            shifts[i] <- best_lag / n_time
          } else {
            shifts[i] <- 0
          }
        }
        
        # ---- AUDIT (P0.8): the warp is the shift, and nothing else -----------
        # This block used to read:
        #     base_warp  <- time_points - shifts[i] * 0.5
        #     distortion <- sin(pi * time_points) * runif(1, -0.03, 0.03)
        #     warp_functions[,i] <- pmin(1, pmax(0, base_warp + distortion))
        # under the comment "# Add slight S-curve for visualization". A random
        # draw was being added to an ESTIMATED transformation, on a module with
        # no set.seed, so the same data gave different answers on every run --
        # measured at up to 5.9x the size of the shift being estimated. The
        # clipping then destroyed monotonicity at the ends: a clipped warp maps
        # a whole interval of new time onto one old time, and approx() hands
        # back the same value for all of it.
        #
        # A shift warp is h(t) = t - s. It is monotone by construction and
        # deterministic. Forcing the endpoints by ASSIGNMENT (as the pre-P0.8
        # code did) is what broke monotonicity when |s| was large, so the shift
        # is limited to what the domain can absorb instead. See P4.2 below for
        # what h actually maps onto, and for the boundary rule.
        # AUDIT (P4.2): the comment above used to claim three properties for
        # this warp -- monotone by construction, deterministic, and (in its own
        # words) anchored at the endpoints. The first two are true. THE THIRD IS
        # FALSE, and a reviewer was right to call it out. (The exact old phrase
        # is deliberately not repeated here: tests/warp_family_test.R greps the
        # source for it, and a fix that quotes the sentence it is removing makes
        # its own guard vacuous -- which has happened twice in this audit.)
        # h(t) = t - s maps [0, 1] onto [-s, 1 - s]. For s = 0.1 that is
        # h(0) = -0.1 and h(1) = 0.9: it is not a map of [0,1] to itself and it
        # does not fix the endpoints. It is a TRANSLATION, which is what shift
        # registration is (Ramsay & Silverman's shift model), and translations
        # are not endpoint-preserving diffeomorphisms. The label was wrong, not
        # the arithmetic.
        #
        # What WAS wrong: the boundary rule. A shift estimated by CIRCULAR
        # cross-correlation (the periodic branch above) was then applied with
        # approx(rule = 2), i.e. constant extrapolation -- so a curve shifted by
        # s had its first |s| of the domain filled with a repeat of the endpoint
        # value, flattening exactly the region the circular estimate said should
        # wrap round from the other end. On 24-hour data with a 2.4 h shift that
        # is a tenth of the cycle replaced by a constant. When the design is
        # periodic the warp now wraps; when it is not, it still clamps, and the
        # amount of extrapolation is reported rather than left implicit.
        s_max <- 0.25   # a shift beyond a quarter of the domain is not identified
        shifts[i] <- max(-s_max, min(s_max, shifts[i]))
        h <- time_points - shifts[i]

        if(periodic) {
          # Wrap into the observed domain: the curve is one period of a cycle.
          span <- diff(range(time_points))
          h_use <- min(time_points) +
                   ((h - min(time_points)) %% span)
          extrap_frac[i] <- 0
        } else {
          h_use <- h
          rng <- range(time_points)
          extrap_frac[i] <- mean(h < rng[1] | h > rng[2])
        }
        warp_functions[,i] <- h_use

        # Apply warping
        if(abs(shifts[i]) > 1e-8) {
          # Interpolate curve at warped time points. rule = 2 is a no-op in the
          # periodic branch, where h_use is inside the domain by construction.
          registered_curves[,i] <- approx(time_points, curves[,i],
                                          xout = h_use,
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
        extrap_frac = extrap_frac,
        warp_direction = "registered -> original",   # P5.6, module-wide
        boundary = if(periodic) "periodic wrap" else "constant extrapolation",
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
      
      # AUDIT (P0.8): the quadratic family alpha*t^2 + (1-alpha)*t has
      # derivative 2*alpha*t + (1-alpha), which is negative near t = 0 whenever
      # alpha > 1. The UI's default search range was c(0.5, 2), so the optimiser
      # could and did return non-monotone maps; clipping them to [0,1] produced
      # a flat leading segment that collapses many new times onto one old time
      # (at alpha = 2, 12 of 24 grid points mapped to time 0). A time warp that
      # is not strictly increasing is not a reparameterisation of time.
      #
      # AUDIT (P4.1): that fix was right about monotonicity and WRONG about
      # where each family's identity lies, and a reviewer caught both halves.
      #
      #   exponential  h(t) = (e^(at) - 1)/(e^a - 1).  Its identity is the LIMIT
      #                a -> 0 (L'Hopital: (at + O(a^2))/(a + O(a^2)) -> t). The
      #                code special-cased  abs(alpha - 1) < 0.001  and returned
      #                t there. At a = 1 the function is (e^t - 1)/(e - 1),
      #                which is not t -- it is the most-curved member of the
      #                family the default range could reach. So the guard put a
      #                DISCONTINUITY in the objective in the middle of the
      #                default search interval [0.5, 2], and optimize() can and
      #                does converge onto it: a curve reported as "alpha = 1,
      #                identity" had in fact been left unwarped by accident
      #                while its neighbours were warped by a real map. And the
      #                genuinely singular point, a = 0 (0/0), had no guard at
      #                all -- it was simply outside the range the UI allowed.
      #   quadratic    identity at a = 0. The P0.8 clamp forced a >= 0.05 and
      #                the UI default gave [0.5, 1], so the identity was NOT
      #                REACHABLE: a curve needing no registration was deformed
      #                anyway, by at least a = 0.5. The 0.05 floor was copied
      #                from the power family, where it is correct; here it is
      #                not. The family is monotone on the whole of (-1, 1).
      #   logistic     identity as its steepness -> 0, also excluded, and also
      #                0/0 at exactly 0 (L1 - L0 = 0). Monotone for every a != 0.
      #   power        h(t) = t^a, identity at a = 1, monotone for a > 0. The
      #                only family whose identity the old range contained.
      #
      # Each family now declares its identity and the open interval on which it
      # is a strictly increasing bijection of [0,1]. The user's range is clamped
      # to that interval AND widened, if necessary, to contain the identity --
      # so "no warping needed" is always inside the search space, in every
      # family. Where the identity is a removable singularity it is returned
      # exactly rather than divided out.
      fam_spec <- switch(family,
        "power"       = list(identity = 1, lo = 0.05, hi = Inf),
        "exponential" = list(identity = 0, lo = -Inf, hi = Inf),
        "quadratic"   = list(identity = 0, lo = -0.999, hi = 0.999),
        "logistic"    = list(identity = 0, lo = -Inf, hi = Inf),
        list(identity = 0, lo = -Inf, hi = Inf))

      param_range <- sort(as.numeric(param_range))
      param_range <- c(max(param_range[1], fam_spec$lo),
                       min(param_range[2], fam_spec$hi))
      # always able to say "this curve needs no warping"
      param_range <- c(min(param_range[1], fam_spec$identity),
                       max(param_range[2], fam_spec$identity))
      param_range <- c(max(param_range[1], fam_spec$lo),
                       min(param_range[2], fam_spec$hi))
      if (param_range[1] >= param_range[2])
        param_range <- c(fam_spec$identity - 1e-3, fam_spec$identity + 1e-3)

      # Define warping function. Every branch returns a map with h(0) = 0,
      # h(1) = 1, strictly increasing on the declared interval, so no clipping
      # to [0,1] is needed -- the pmin/pmax that used to be here only ever
      # masked a non-monotone map instead of rejecting it. They are kept as a
      # floating-point tidy-up (widths of order 1e-16), not as a repair.
      eps_id <- 1e-8
      warp_func <- function(t, alpha) {
        h <- switch(family,
               "power" = t^alpha,
               "exponential" = {
                 # identity is the limit a -> 0, where the expression is 0/0
                 if (abs(alpha) < eps_id) t
                 else (exp(alpha * t) - 1) / (exp(alpha) - 1)
               },
               "quadratic" = alpha * t^2 + (1 - alpha) * t,
               "logistic" = {
                 # identity is the limit a -> 0, where L1 - L0 is 0
                 if (abs(alpha) < eps_id) t
                 else {
                   L  <- function(x) 1 / (1 + exp(-alpha * (x - 0.5)))
                   L0 <- L(0); L1 <- L(1)
                   (L(t) - L0) / (L1 - L0)
                 }
               },
               t)
        pmin(1, pmax(0, h))
      }
      
      # AUDIT (P5.11, found by tests/registration_effectiveness_test.R): the
      # search was a bare optimize() over the whole parameter range.
      # optimize() is golden-section plus parabolic interpolation and it
      # assumes the objective is UNIMODAL on the interval. The registration
      # objective is not. On a sharply peaked curve -- which is what a
      # circadian profile is -- the alignment SSE has a deep, narrow well at
      # the correct parameter surrounded by a wide plateau where the curve has
      # been pushed off its own peak. Measured on an already-aligned sample
      # with the power family: SSE is 0.008 at alpha = 1 and about 20 for
      # alpha anywhere in 0.05-0.5 or 1.5-6, and optimize() on [0.05, 6]
      # returned alpha = 6.000, deforming the time axis by 0.58 of the domain
      # on curves that needed no registration at all.
      #
      # This got WORSE at P4.1, not better: widening the ranges to make each
      # family's identity reachable was right, but a wider interval gives
      # optimize() more plateau to get lost on. Registrations run with the
      # narrow pre-P4.1 ranges were partly protected by luck.
      #
      # A coarse grid scan followed by refinement inside the winning bracket
      # is robust to this and costs about 40 extra evaluations per curve --
      # the same pattern fck_auto_lambda() already uses for the GCV search.
      fck_min_1d <- function(f, lo, hi, n_grid = 41) {
        g <- seq(lo, hi, length.out = n_grid)
        v <- vapply(g, function(z) {
          y <- tryCatch(f(z), error = function(e) NA_real_)
          if (is.finite(y)) y else Inf
        }, numeric(1))
        if (all(!is.finite(v))) return(NA_real_)
        k <- which.min(v)
        a <- g[max(1, k - 1)]; b <- g[min(length(g), k + 1)]
        if (a < b) tryCatch(stats::optimize(f, c(a, b), tol = 1e-5)$minimum,
                            error = function(e) g[k])
        else g[k]
      }

      # Optimize warping for each curve
      for(i in 1:n_curves) {
        # Objective function
        objective <- function(alpha) {
          warped_time <- warp_func(time_points, alpha)
          warped_curve <- approx(time_points, curves[,i], xout = warped_time, rule = 2)$y
          sum((warped_curve - mean_curve)^2, na.rm = TRUE)
        }

        a_hat <- fck_min_1d(objective, param_range[1], param_range[2])
        if (!is.finite(a_hat)) a_hat <- fam_spec$identity
        result <- list(minimum = a_hat)
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
        warp_direction = "registered -> original",   # P5.6, module-wide
        param_range_used = param_range,
        identity_alpha = fam_spec$identity,
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
        n_rejected <- 0L; rejected_ids <- integer(0)

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
          
          # P5.5: built and VALIDATED in one place. A subject whose detected
          # landmarks cross does not get a folded warp; it gets the identity,
          # and is counted.
          h <- fck_landmark_warp(landmark_times, curve_landmarks, time_points)
          if (is.null(h)) {
            n_rejected <- n_rejected + 1L
            rejected_ids <- c(rejected_ids, i)
            h <- time_points
          }
          warp_functions[, i] <- h

          # P5.6: the module-wide contract, x_R(t) = x(h(t)).
          registered_curves[, i] <- approx(time_points, curves[, i],
                                           xout = h, rule = 2)$y
        }

        if (n_rejected > 0)
          showNotification(sprintf(
            paste("Landmark registration: %d of %d curve(s) had crossed or duplicated",
                  "landmarks, which cannot define a monotone time warp. Those curves",
                  "were left UNREGISTERED (identity warp) rather than folded:",
                  "subject%s %s. Adjust the landmark positions or use fewer."),
            n_rejected, n_curves, if (n_rejected == 1) "" else "s",
            paste(rejected_ids, collapse = ", ")),
            type = "warning", duration = 15)
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
        
        # ---- AUDIT (P0.8): this branch used to be ---------------------------
        #     distortion <- sin(2*pi*time_points) * runif(1, -0.02, 0.02)
        #     warp_functions[,i]    <- pmin(1, pmax(0, time_points + distortion))
        #     registered_curves[,i] <- curves[,i]
        # under the comment "# Simple landmark alignment". It computed landmark
        # positions, discarded them, returned a warp made entirely of random
        # numbers, and did not warp the curves at all -- while the UI announced
        # "Time-warped PCA completed!".
        #
        # A landmark warp is the monotone piecewise-linear map that carries each
        # curve's own landmarks onto the reference landmarks. Detecting the
        # landmarks is done by the same rule for every curve (the largest local
        # extremum nearest each reference position), so the result is
        # deterministic.
        ref_lm <- sort(unique(pmin(1 - 1e-6, pmax(1e-6, landmark_times))))
        find_lm <- function(v, target) {
          # search a window around the reference position for the strongest
          # turning point; fall back to the reference position itself
          w <- which(abs(time_points - target) <= 0.12)
          if (length(w) < 3) return(target)
          d <- diff(v[w])
          turn <- which(d[-length(d)] * d[-1] < 0)
          if (!length(turn)) return(target)
          best <- turn[which.max(abs(d[turn]))]
          time_points[w[best + 1]]
        }
        n_rejected <- 0L; rejected_ids <- integer(0)
        for(i in 1:n_curves) {
          own <- vapply(ref_lm, function(tt) find_lm(curves[, i], tt), numeric(1))
          # P5.6: knots are (reference -> own), i.e. h maps REGISTERED time to
          # ORIGINAL time -- the same direction as the manual branch. This
          # branch used to build (own -> reference), the inverse map, and then
          # apply it as approx(h, curve, xout = t), which inverts it again. The
          # registered curves were right; warp_functions held the wrong object,
          # and every statistic computed from it was not comparable with the
          # other methods'. P5.5: crossed landmarks are rejected, not clipped.
          h <- fck_landmark_warp(ref_lm, own, time_points)
          if (is.null(h)) {
            n_rejected <- n_rejected + 1L
            rejected_ids <- c(rejected_ids, i)
            h <- time_points
          }
          warp_functions[, i] <- h
          registered_curves[, i] <- approx(time_points, curves[, i],
                                           xout = h, rule = 2)$y
        }

        if (n_rejected > 0)
          showNotification(sprintf(
            paste("Automatic landmark registration: %d of %d curve(s) gave crossed or",
                  "duplicated landmarks and were left UNREGISTERED (identity warp)",
                  "rather than folded: subject%s %s."),
            n_rejected, n_curves, if (n_rejected == 1) "" else "s",
            paste(rejected_ids, collapse = ", ")),
            type = "warning", duration = 15)
      }
      
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        method = "landmark",
        n_rejected = n_rejected,
        rejected_ids = rejected_ids,
        warp_direction = "registered -> original",   # P5.6, module-wide
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
        line = list(color = FCK_NEUTRAL, width = 2, dash = 'dash'),
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
        # Reserve room under the title: the warping functions run the full
        # height of the panel, so a title with the default top margin sits on
        # the topmost curve.
        title = list(text = "Time warping functions", x = 0.5, xanchor = "center"),
        margin = list(t = 60),
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
                    line = list(color = FCK_EMPHASIS, width = 2))
        
        title_text <- "Click to add landmarks on mean curve"
      } else {
        # Plot individual subject
        subj_idx <- as.numeric(input$selected_subject)
        subj_curve <- values$data[subj_idx,]
        
        p <- plot_ly(source = "landmark_source") %>%
          add_trace(x = time_points, y = subj_curve,
                    type = 'scatter', mode = 'lines',
                    name = paste('Subject', subj_idx),
                    line = list(color = FCK_SERIES1, width = 2))
        
        title_text <- paste("Click to add landmarks for Subject", subj_idx)
      }
      
      # Add existing landmarks if any
      if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
        p <- p %>% add_trace(x = values$landmark_points$x,
                             y = values$landmark_points$y,
                             type = 'scatter', mode = 'markers',
                             marker = list(color = FCK_SERIES2, size = 10),
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
