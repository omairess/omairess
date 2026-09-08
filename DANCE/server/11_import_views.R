# ==========================================================================
# server/11_import_views.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 1817-1843  (recommended n_basis when the data change)
#   WaPaa1_3.R lines 1902-1967  (raw data plot with clock-time axis)
# ==========================================================================
  observe({
    req(values$data)
    
    n_time <- ncol(values$data)
    
    # Calculate recommended n_basis based on number of time points
    recommended_nb <- if(n_time <= 10) {
      max(4, n_time - 2)
    } else if(n_time <= 20) {
      max(10, round(n_time * 0.6))
    } else if(n_time <= 50) {
      max(15, round(n_time * 0.5))
    } else {
      max(20, round(n_time * 0.4))
    }
    
    # Ensure it's within bounds
    recommended_nb <- min(recommended_nb, n_time - 2, 100)
    recommended_nb <- max(recommended_nb, 4)
    
    # Update both n_basis inputs
    updateNumericInput(session, "n_basis", value = recommended_nb)
    updateNumericInput(session, "n_basis_manual", value = recommended_nb)
    
    cat(sprintf("Data loaded: %d time points → Recommended n_basis = %d\n", 
                n_time, recommended_nb))
  })

  output$raw_data_plot <- renderPlot({
    if(is.null(values$data)) {
      plot(1, type = "n", xlab = "", ylab = "", main = "No analysis data available")
      return()
    }
    
    tryCatch({
      n_time <- ncol(values$data)
      n_subj <- nrow(values$data)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      hour_labels <- get_hour_labels()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      
      if(!is.null(values$group_labels) && length(unique(values$group_labels)) > 1) {
        # Plot with groups
        groups <- levels(values$group_labels)
        n_groups <- length(groups)
        
        # Create scalable color palette
        base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
        colors <- colorRampPalette(base_cols)(n_groups)
        
        matplot(time_points_plot, t(values$data), type = "l", 
                col = rgb(0.5, 0.5, 0.5, 0.1), lty = 1,
                xlab = time_label, ylab = "Value", 
                main = paste("Raw Functional Data (", n_subj, "subjects,", 
                             n_groups, "groups)"),
                xaxt = if(!is.null(hour_labels)) "n" else "s")
        
        # Add custom x-axis with hour labels if available
        if(!is.null(hour_labels)) {
          n_labels <- length(hour_labels)
          if(n_labels > 15) {
            step <- ceiling(n_labels / 10)
            indices <- seq(1, n_labels, by = step)
            axis(1, at = time_points_plot[indices], labels = sapply(hour_labels[indices], decimal_to_hhmm))
          } else {
            axis(1, at = time_points_plot, labels = sapply(hour_labels, decimal_to_hhmm))
          }
        }
        
        # Add group means
        for(i in 1:n_groups) {
          group_idx <- which(values$group_labels == groups[i])
          if(length(group_idx) > 0) {
            group_mean <- colMeans(values$data[group_idx, , drop = FALSE])
            lines(time_points_plot, group_mean, col = colors[i], lwd = 3)
          }
        }
        
        legend("topright", legend = groups, col = colors, lwd = 3)
      } else {
        # Plot without groups
        matplot(time_points_plot, t(values$data), type = "l", 
                col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1,
                xlab = time_label, ylab = "Value", 
                main = paste("Raw Functional Data (", n_subj, "subjects)"))
        lines(time_points_plot, colMeans(values$data), col = "red", lwd = 3)
        legend("topright", legend = "Mean", col = "red", lwd = 3)
      }
    }, error = function(e) {
      plot(1, type = "n", xlab = "", ylab = "", main = paste("Error:", e$message))
    })
  })
