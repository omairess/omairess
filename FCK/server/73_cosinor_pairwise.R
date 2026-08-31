# ==========================================================================
# server/73_cosinor_pairwise.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 7264-7859  (cosinor pairwise group tests; ids prefixed hp_)
# ==========================================================================
  # ==============================================================================
  # PAIRWISE COMPARISONS MODULE
  # ==============================================================================

  # Run pairwise comparisons
  observeEvent(input$hp_run, {
    req(values$harmonic_model)
    mod <- values$harmonic_model

    # Check if groups are defined
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      showNotification("Pairwise comparisons require 2 or more groups. Please define groups in Harmonic Regression tab.",
                      type = "error", duration = 5)
      return()
    }

    # Get parameter to compare
    param <- input$hp_param
    correction <- input$hp_correction

    # Get individual parameters with group information
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]
    } else {
      showNotification("No group variable selected. Please select groups in Harmonic Regression tab.",
                      type = "error", duration = 5)
      return()
    }

    # Check if parameter exists in data
    if(!param %in% names(params)) {
      showNotification(paste("Parameter", param, "not available in current model."),
                      type = "error", duration = 5)
      return()
    }

    # Get groups
    groups <- levels(params$group)
    n_groups <- length(groups)

    if(n_groups < 2) {
      showNotification("Need at least 2 groups for pairwise comparisons.",
                      type = "error", duration = 5)
      return()
    }

    # Perform pairwise t-tests
    n_comparisons <- choose(n_groups, 2)
    results <- data.frame(
      comparison = character(n_comparisons),
      group1 = character(n_comparisons),
      group2 = character(n_comparisons),
      mean1 = numeric(n_comparisons),
      sd1 = numeric(n_comparisons),
      n1 = integer(n_comparisons),
      mean2 = numeric(n_comparisons),
      sd2 = numeric(n_comparisons),
      n2 = integer(n_comparisons),
      mean_diff = numeric(n_comparisons),
      t_stat = numeric(n_comparisons),
      df = numeric(n_comparisons),
      p_value = numeric(n_comparisons),
      cohens_d = numeric(n_comparisons),
      ci_lower = numeric(n_comparisons),
      ci_upper = numeric(n_comparisons),
      stringsAsFactors = FALSE
    )

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    idx <- 1
    for(i in 1:(n_groups-1)) {
      for(j in (i+1):n_groups) {
        g1 <- groups[i]
        g2 <- groups[j]

        vals1 <- params[[param]][params$group == g1]
        vals2 <- params[[param]][params$group == g2]

        # Remove NAs
        vals1 <- vals1[!is.na(vals1)]
        vals2 <- vals2[!is.na(vals2)]

        if(length(vals1) < 2 || length(vals2) < 2) {
          next
        }

        if(is_circular) {
          # CIRCULAR STATISTICS for acrophase parameters
          # Convert hours to radians (assuming 24-hour period)
          period <- mod$period
          rad1 <- vals1 * 2 * pi / period
          rad2 <- vals2 * 2 * pi / period

          # Circular means (in radians)
          cmean1_rad <- circular_mean(rad1)
          cmean2_rad <- circular_mean(rad2)

          # Convert back to hours for display
          cmean1 <- cmean1_rad * period / (2 * pi)
          cmean2 <- cmean2_rad * period / (2 * pi)

          # Ensure positive (0-24h range)
          if(cmean1 < 0) cmean1 <- cmean1 + period
          if(cmean2 < 0) cmean2 <- cmean2 + period

          # Circular standard deviations (in radians, then convert to hours)
          csd1_rad <- circular_sd(rad1)
          csd2_rad <- circular_sd(rad2)
          csd1 <- csd1_rad * period / (2 * pi)
          csd2 <- csd2_rad * period / (2 * pi)

          # Angular difference (shortest arc)
          ang_diff_rad <- cmean1_rad - cmean2_rad
          # Normalize to [-pi, pi]
          ang_diff_rad <- atan2(sin(ang_diff_rad), cos(ang_diff_rad))
          ang_diff <- ang_diff_rad * period / (2 * pi)

          # Mean resultant lengths (measure of concentration)
          r1 <- mean_resultant_length(rad1)
          r2 <- mean_resultant_length(rad2)

          # Watson-Williams test for two groups
          ww <- watson_williams_test(list(rad1, rad2))

          # Effect size for circular data: difference in mean resultant lengths
          # (alternative: use V statistic, but difference in r is more interpretable)
          effect_size <- r1 - r2

          # For confidence interval on angular difference, use approximate circular CI
          # (simplified: not implemented here, set to NA)
          ci_lower <- NA
          ci_upper <- NA

          # Store results
          results$comparison[idx] <- paste(g1, "vs", g2)
          results$group1[idx] <- as.character(g1)
          results$group2[idx] <- as.character(g2)
          results$mean1[idx] <- cmean1
          results$sd1[idx] <- csd1
          results$n1[idx] <- length(vals1)
          results$mean2[idx] <- cmean2
          results$sd2[idx] <- csd2
          results$n2[idx] <- length(vals2)
          results$mean_diff[idx] <- ang_diff
          results$t_stat[idx] <- ww$F  # F-statistic from Watson-Williams
          results$df[idx] <- ww$df2  # Store df2 in df column
          results$p_value[idx] <- ww$p
          results$cohens_d[idx] <- effect_size  # Actually difference in mean resultant lengths
          results$ci_lower[idx] <- ci_lower
          results$ci_upper[idx] <- ci_upper

        } else {
          # REGULAR STATISTICS for non-circular parameters
          # Perform t-test
          t_result <- t.test(vals1, vals2, var.equal = FALSE)

          # Calculate Cohen's d
          pooled_sd <- sqrt(((length(vals1)-1)*sd(vals1)^2 + (length(vals2)-1)*sd(vals2)^2) /
                           (length(vals1) + length(vals2) - 2))
          cohens_d <- (mean(vals1) - mean(vals2)) / pooled_sd

          results$comparison[idx] <- paste(g1, "vs", g2)
          results$group1[idx] <- as.character(g1)
          results$group2[idx] <- as.character(g2)
          results$mean1[idx] <- mean(vals1)
          results$sd1[idx] <- sd(vals1)
          results$n1[idx] <- length(vals1)
          results$mean2[idx] <- mean(vals2)
          results$sd2[idx] <- sd(vals2)
          results$n2[idx] <- length(vals2)
          results$mean_diff[idx] <- mean(vals1) - mean(vals2)
          results$t_stat[idx] <- t_result$statistic
          results$df[idx] <- t_result$parameter
          results$p_value[idx] <- t_result$p.value
          results$cohens_d[idx] <- cohens_d
          results$ci_lower[idx] <- t_result$conf.int[1]
          results$ci_upper[idx] <- t_result$conf.int[2]
        }

        idx <- idx + 1
      }
    }

    # Remove empty rows
    results <- results[results$comparison != "", ]

    # Apply multiple comparison correction
    if(correction != "none") {
      results$p_adjusted <- p.adjust(results$p_value, method = correction)
    } else {
      results$p_adjusted <- results$p_value
    }

    # Store results
    values$hp_pairwise_results <- results
    values$hp_pairwise_param <- param
    values$hp_pairwise_correction <- correction

    showNotification("Pairwise comparisons completed!", type = "message", duration = 3)
  })

  # Display pairwise results
  output$hp_results <- renderPrint({
    req(values$hp_pairwise_results)
    results <- values$hp_pairwise_results
    param <- values$hp_pairwise_param
    correction <- values$hp_pairwise_correction

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    cat("=== Pairwise Group Comparisons ===\n\n")
    cat("Parameter:", param, "\n")
    if(is_circular) {
      cat("Data type: Circular (using Watson-Williams test)\n")
    } else {
      cat("Data type: Linear (using Welch's t-test)\n")
    }
    cat("Correction method:", correction, "\n")
    cat("Number of comparisons:", nrow(results), "\n\n")

    for(i in 1:nrow(results)) {
      r <- results[i, ]
      cat("---\n")
      cat(sprintf("%s:\n", r$comparison))

      if(is_circular) {
        # Circular statistics display
        cat(sprintf("  Group 1: Circular mean=%.2f h, Circular SD=%.2f h, n=%d\n", r$mean1, r$sd1, r$n1))
        cat(sprintf("  Group 2: Circular mean=%.2f h, Circular SD=%.2f h, n=%d\n", r$mean2, r$sd2, r$n2))
        cat(sprintf("  Angular difference: %.2f h\n", r$mean_diff))

        # No CI for circular data (not implemented)
        # if(input$hp_show_ci && !is.na(r$ci_lower)) {
        #   cat(sprintf("  95%% CI: [%.3f, %.3f]\n", r$ci_lower, r$ci_upper))
        # }

        cat(sprintf("  Watson-Williams F(%d, %d) = %.3f, p = %.4f", 1, r$df, r$t_stat, r$p_value))

      } else {
        # Regular statistics display
        cat(sprintf("  Group 1: M=%.3f, SD=%.3f, n=%d\n", r$mean1, r$sd1, r$n1))
        cat(sprintf("  Group 2: M=%.3f, SD=%.3f, n=%d\n", r$mean2, r$sd2, r$n2))
        cat(sprintf("  Difference: %.3f\n", r$mean_diff))

        if(input$hp_show_ci) {
          cat(sprintf("  95%% CI: [%.3f, %.3f]\n", r$ci_lower, r$ci_upper))
        }

        cat(sprintf("  t(%.1f) = %.3f, p = %.4f", r$df, r$t_stat, r$p_value))
      }

      if(correction != "none") {
        cat(sprintf(", p_adj = %.4f", r$p_adjusted))
      }

      if(r$p_adjusted < 0.001) {
        cat(" ***")
      } else if(r$p_adjusted < 0.01) {
        cat(" **")
      } else if(r$p_adjusted < 0.05) {
        cat(" *")
      }
      cat("\n")

      if(input$hp_show_effect_size) {
        if(is_circular) {
          # For circular data: difference in mean resultant lengths
          cat(sprintf("  Δr̄ (difference in mean resultant length): %.3f", r$cohens_d))
          if(abs(r$cohens_d) < 0.1) {
            cat(" (small)")
          } else if(abs(r$cohens_d) < 0.3) {
            cat(" (medium)")
          } else {
            cat(" (large)")
          }
        } else {
          # For linear data: Cohen's d
          cat(sprintf("  Cohen's d: %.3f", r$cohens_d))
          if(abs(r$cohens_d) < 0.2) {
            cat(" (negligible)")
          } else if(abs(r$cohens_d) < 0.5) {
            cat(" (small)")
          } else if(abs(r$cohens_d) < 0.8) {
            cat(" (medium)")
          } else {
            cat(" (large)")
          }
        }
        cat("\n")
      }
    }

    cat("\n---\n")
    cat("Significance codes: *** p<0.001, ** p<0.01, * p<0.05\n")
    if(is_circular) {
      cat("\nNote: Circular means are in hours (0-24). Angular difference is the shortest arc.\n")
      cat("Effect size Δr̄ measures difference in concentration (mean resultant lengths).\n")
    }
  })

  # Pairwise comparison plot
  output$hp_plot <- renderPlot({
    req(values$hp_pairwise_results, values$harmonic_model)

    param <- values$hp_pairwise_param
    mod <- values$harmonic_model

    # Get individual parameters with group information
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]
    } else {
      return(NULL)
    }

    if(!param %in% names(params)) {
      return(NULL)
    }

    # Get parameter label
    param_labels <- c(
      "mesor" = "MESOR",
      "amplitude_1" = "H1 Amplitude",
      "amplitude_2" = "H2 Amplitude",
      "amplitude_3" = "H3 Amplitude",
      "acrophase_time_1" = "H1 Acrophase (hours)",
      "acrophase_time_2" = "H2 Acrophase (hours)",
      "acrophase_time_3" = "H3 Acrophase (hours)",
      "r_squared" = "R²",
      "A_sat" = "A_sat",
      "tau" = "τ (tau)",
      "percent_S" = "Process S (%)",
      "percent_C" = "Process C (%)"
    )
    param_label <- ifelse(param %in% names(param_labels), param_labels[param], param)

    # Create boxplot with violin overlay
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

    # Prepare data
    groups <- levels(params$group)
    n_groups <- length(groups)
    plot_data <- list()
    for(g in groups) {
      plot_data[[as.character(g)]] <- params[[param]][params$group == g]
    }

    # Boxplot
    boxplot(plot_data,
            main = paste("Group Comparison:", param_label),
            ylab = param_label,
            xlab = "Group",
            col = rainbow(n_groups, alpha = 0.3),
            border = rainbow(n_groups),
            notch = TRUE,
            las = 1,
            cex.axis = 1.2,
            cex.lab = 1.3,
            cex.main = 1.4)

    # Add individual points with jitter
    for(i in 1:n_groups) {
      g <- groups[i]
      vals <- params[[param]][params$group == g]
      vals <- vals[!is.na(vals)]
      points(jitter(rep(i, length(vals)), amount = 0.1), vals,
             col = rainbow(n_groups, alpha = 0.5)[i],
             pch = 19, cex = 0.8)
    }

    # Add group means
    for(i in 1:n_groups) {
      g <- groups[i]
      vals <- params[[param]][params$group == g]
      vals <- vals[!is.na(vals)]
      points(i, mean(vals), pch = 18, cex = 2.5, col = "black")
    }

    # Add significance brackets
    results <- values$hp_pairwise_results
    y_max <- max(params[[param]], na.rm = TRUE)
    y_min <- min(params[[param]], na.rm = TRUE)
    y_range <- y_max - y_min

    sig_results <- results[results$p_adjusted < 0.05, ]
    if(nrow(sig_results) > 0) {
      bracket_y <- y_max + y_range * 0.05
      for(i in 1:min(nrow(sig_results), 5)) {  # Show max 5 brackets
        g1_idx <- which(groups == sig_results$group1[i])
        g2_idx <- which(groups == sig_results$group2[i])

        y_pos <- bracket_y + (i - 1) * y_range * 0.08

        # Draw bracket
        segments(g1_idx, y_pos, g2_idx, y_pos, lwd = 1.5)
        segments(g1_idx, y_pos, g1_idx, y_pos - y_range * 0.02, lwd = 1.5)
        segments(g2_idx, y_pos, g2_idx, y_pos - y_range * 0.02, lwd = 1.5)

        # Add significance stars
        p_val <- sig_results$p_adjusted[i]
        stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else "*"
        text((g1_idx + g2_idx) / 2, y_pos + y_range * 0.02, stars, cex = 1.2)
      }
    }

    # Add legend
    legend("topleft", legend = c("Mean", "Individual"),
           pch = c(18, 19), col = c("black", "gray"),
           pt.cex = c(2.5, 0.8), bty = "n")
  })

  # Dynamic help text for comparison matrix
  output$hp_matrix_help <- renderUI({
    req(values$hp_pairwise_param)
    param <- values$hp_pairwise_param
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    if(is_circular) {
      helpText("Lower triangle: p-values | Upper triangle: Δr̄ (difference in mean resultant length)")
    } else {
      helpText("Lower triangle: p-values | Upper triangle: Cohen's d (effect sizes)")
    }
  })

  # Pairwise comparison matrix
  output$hp_matrix <- renderTable({
    req(values$hp_pairwise_results)
    results <- values$hp_pairwise_results
    param <- values$hp_pairwise_param

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    # Get unique groups
    groups <- unique(c(results$group1, results$group2))
    n_groups <- length(groups)

    # Create matrix
    mat <- matrix("", nrow = n_groups, ncol = n_groups,
                  dimnames = list(groups, groups))

    for(i in 1:nrow(results)) {
      r <- results[i, ]
      g1_idx <- which(groups == r$group1)
      g2_idx <- which(groups == r$group2)

      # Lower triangle: p-values
      p_str <- sprintf("%.4f", r$p_adjusted)
      if(r$p_adjusted < 0.001) p_str <- paste0(p_str, " ***")
      else if(r$p_adjusted < 0.01) p_str <- paste0(p_str, " **")
      else if(r$p_adjusted < 0.05) p_str <- paste0(p_str, " *")
      mat[g2_idx, g1_idx] <- p_str

      # Upper triangle: effect sizes
      d_str <- sprintf("%.3f", r$cohens_d)
      mat[g1_idx, g2_idx] <- d_str
    }

    # Convert to data frame
    mat_df <- as.data.frame(mat)
    mat_df <- cbind(Group = rownames(mat_df), mat_df)
    mat_df
  }, rownames = FALSE, striped = TRUE, bordered = TRUE)

  # Export pairwise results
  output$hp_export_results <- downloadHandler(
    filename = function() paste0("pairwise_comparisons_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$hp_pairwise_results)
      results <- values$hp_pairwise_results
      param <- values$hp_pairwise_param
      correction <- values$hp_pairwise_correction
      is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

      # Create header with metadata
      con <- file(file, "w")
      writeLines(paste0("# Pairwise Comparisons: ", param), con)
      writeLines(paste0("# Date: ", Sys.Date()), con)
      if(is_circular) {
        writeLines("# Data type: Circular (Watson-Williams test)", con)
        writeLines("# Note: mean1/mean2 are circular means; sd1/sd2 are circular SDs", con)
        writeLines("# mean_diff is angular difference (shortest arc); cohens_d is Δr̄", con)
        writeLines("# t_stat is Watson-Williams F-statistic", con)
      } else {
        writeLines("# Data type: Linear (Welch's t-test)", con)
        writeLines("# Note: mean1/mean2 are arithmetic means; cohens_d is Cohen's d", con)
      }
      writeLines(paste0("# Correction method: ", correction), con)
      close(con)

      # Append results
      write.table(results, file, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)
    }
  )

  # Export pairwise plot
  output$hp_export_plot <- downloadHandler(
    filename = function() paste0("pairwise_plot_", Sys.Date(), ".png"),
    content = function(file) {
      req(values$hp_pairwise_results, values$harmonic_model)
      png(file, width = 1200, height = 800, res = 120)

      # Recreate the plot (same code as renderPlot)
      param <- values$hp_pairwise_param
      mod <- values$harmonic_model

      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]
        params <- mod$individual_params
        params$group <- as.factor(group_var[params$subject])
        params <- params[!is.na(params$group), ]

        if(param %in% names(params)) {
          param_labels <- c(
            "mesor" = "MESOR",
            "amplitude_1" = "H1 Amplitude",
            "amplitude_2" = "H2 Amplitude",
            "amplitude_3" = "H3 Amplitude",
            "acrophase_time_1" = "H1 Acrophase (hours)",
            "acrophase_time_2" = "H2 Acrophase (hours)",
            "acrophase_time_3" = "H3 Acrophase (hours)",
            "r_squared" = "R²",
            "A_sat" = "A_sat",
            "tau" = "τ (tau)",
            "percent_S" = "Process S (%)",
            "percent_C" = "Process C (%)"
          )
          param_label <- ifelse(param %in% names(param_labels), param_labels[param], param)

          groups <- levels(params$group)
          n_groups <- length(groups)
          plot_data <- list()
          for(g in groups) {
            plot_data[[as.character(g)]] <- params[[param]][params$group == g]
          }

          par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
          boxplot(plot_data,
                  main = paste("Group Comparison:", param_label),
                  ylab = param_label,
                  xlab = "Group",
                  col = rainbow(n_groups, alpha = 0.3),
                  border = rainbow(n_groups),
                  notch = TRUE,
                  las = 1,
                  cex.axis = 1.2,
                  cex.lab = 1.3,
                  cex.main = 1.4)

          for(i in 1:n_groups) {
            g <- groups[i]
            vals <- params[[param]][params$group == g]
            vals <- vals[!is.na(vals)]
            points(jitter(rep(i, length(vals)), amount = 0.1), vals,
                   col = rainbow(n_groups, alpha = 0.5)[i],
                   pch = 19, cex = 0.8)
            points(i, mean(vals), pch = 18, cex = 2.5, col = "black")
          }

          results <- values$hp_pairwise_results
          y_max <- max(params[[param]], na.rm = TRUE)
          y_min <- min(params[[param]], na.rm = TRUE)
          y_range <- y_max - y_min

          sig_results <- results[results$p_adjusted < 0.05, ]
          if(nrow(sig_results) > 0) {
            bracket_y <- y_max + y_range * 0.05
            for(i in 1:min(nrow(sig_results), 5)) {
              g1_idx <- which(groups == sig_results$group1[i])
              g2_idx <- which(groups == sig_results$group2[i])
              y_pos <- bracket_y + (i - 1) * y_range * 0.08
              segments(g1_idx, y_pos, g2_idx, y_pos, lwd = 1.5)
              segments(g1_idx, y_pos, g1_idx, y_pos - y_range * 0.02, lwd = 1.5)
              segments(g2_idx, y_pos, g2_idx, y_pos - y_range * 0.02, lwd = 1.5)
              p_val <- sig_results$p_adjusted[i]
              stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else "*"
              text((g1_idx + g2_idx) / 2, y_pos + y_range * 0.02, stars, cex = 1.2)
            }
          }

          legend("topleft", legend = c("Mean", "Individual"),
                 pch = c(18, 19), col = c("black", "gray"),
                 pt.cex = c(2.5, 0.8), bty = "n")
        }
      }
      dev.off()
    }
  )
