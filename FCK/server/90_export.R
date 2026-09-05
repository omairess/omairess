# ==========================================================================
# server/90_export.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 8914-9929  (all WaPaa exports + the reproducible-code generator)
#   CIRCAREG.R lines 7192-7198  (FoSR coefficient export -> export_fosr_coefs_csv)
#   CIRCAREG.R lines 7209-7262  (harmonic parameter + summary exports)
# ==========================================================================
  # Export functions
  output$export_scores_csv <- downloadHandler(
    filename = function() paste0("pca_scores_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pca_results) && !is.null(values$pca_results$scores)) {
        scores_df <- data.frame(
          Subject = 1:nrow(values$pca_results$scores),
          values$pca_results$scores
        )
        colnames(scores_df)[-1] <- paste0("PC", 1:ncol(values$pca_results$scores))
        if(!is.null(values$group_labels)) {
          scores_df$Group <- values$group_labels
        }
        write.csv(scores_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No PCA scores available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_loadings_csv <- downloadHandler(
    filename = function() paste0("pca_loadings_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pca_results) && !is.null(values$pca_results$harmonics)) {
        time_points <- seq(0, 1, length.out = 100)
        loadings_mat <- matrix(NA, 100, length(values$pca_results$harmonics))
        for(i in 1:length(values$pca_results$harmonics)) {
          loadings_mat[,i] <- eval.fd(time_points, values$pca_results$harmonics[i])
        }
        loadings_df <- data.frame(Time = time_points, loadings_mat)
        colnames(loadings_df)[-1] <- paste0("PC", 1:ncol(loadings_mat))
        write.csv(loadings_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No loadings available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_warping_csv <- downloadHandler(
    filename = function() paste0("warping_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$warping_results) && !is.null(values$warping_results$warp_functions)) {
        warp_df <- data.frame(
          Time = values$warping_results$time_points,
          values$warping_results$warp_functions
        )
        colnames(warp_df)[-1] <- paste0("Subject", 1:(ncol(warp_df)-1))
        write.csv(warp_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No warping results available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_fanova_results_csv <- downloadHandler(
    filename = function() paste0("fanova_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$fanova_results)) {
        res <- values$fanova_results
        results_df <- data.frame(
          Time = res$time_points,
          F_statistic = res$F_stat,
          p_value_raw = res$p_values_pointwise,
          p_value_adjusted = res$p_values_adjusted,
          significant = res$sig_regions,
          eta_squared = res$eta_squared
        )
        
        for(i in 1:res$n_groups) {
          results_df[paste0("Mean_", res$groups[i])] <- res$group_means[,i]
        }
        
        write.csv(results_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No FANOVA results available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_pairwise_results_csv <- downloadHandler(
    filename = function() paste0("pairwise_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pairwise_results)) {
        results_list <- list()
        
        for(pair_name in values$pairwise_results$pair_names) {
          pair_result <- values$pairwise_results$results[[pair_name]]
          
          results_list[[pair_name]] <- data.frame(
            Time = values$pairwise_results$time_points,
            Mean_Diff = pair_result$mean_diff,
            t_stat = pair_result$t_stat,
            p_raw = pair_result$p_values_pointwise,
            p_adj = pair_result$p_values_adjusted,
            Cohen_d = pair_result$cohens_d,
            CI_lower = pair_result$ci_lower,
            CI_upper = pair_result$ci_upper,
            Significant = pair_result$sig_regions
          )
        }
        
        all_results <- do.call(rbind, results_list)
        all_results$Comparison <- rep(names(results_list), each = length(values$pairwise_results$time_points))
        
        write.csv(all_results, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No pairwise results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster membership
  output$export_cluster_membership_csv <- downloadHandler(
    filename = function() paste0("cluster_membership_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        membership_df <- data.frame(
          Subject_ID = 1:length(values$clustering_results$cluster_assignments),
          Cluster = values$clustering_results$cluster_assignments
        )

        # Add group labels if available
        if(!is.null(values$group_labels)) {
          membership_df$Group <- values$group_labels
        }

        write.csv(membership_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster mean functions
  output$export_cluster_means_csv <- downloadHandler(
    filename = function() paste0("cluster_means_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        means_df <- as.data.frame(values$clustering_results$cluster_means)

        # Use time labels if available
        if(!is.null(values$time_labels) && length(values$time_labels) == ncol(means_df)) {
          colnames(means_df) <- values$time_labels
        } else {
          colnames(means_df) <- paste0("T", 1:ncol(means_df))
        }

        # Add cluster ID column
        means_df <- cbind(Cluster = 1:nrow(means_df), means_df)

        write.csv(means_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster statistics
  output$export_cluster_stats_csv <- downloadHandler(
    filename = function() paste0("cluster_statistics_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        r <- values$clustering_results

        stats_df <- data.frame(
          Cluster = 1:r$k,
          Size = r$cluster_sizes,
          Percentage = round(100 * r$cluster_sizes / sum(r$cluster_sizes), 2),
          Within_SS = r$wcss_per_cluster,
          Avg_Within_SS = r$wcss_per_cluster / r$cluster_sizes
        )

        # Add overall statistics
        overall_stats <- data.frame(
          Metric = c("Total_SS", "Between_SS", "Within_SS", "R_Squared",
                     "Avg_Silhouette", "Calinski_Harabasz"),
          Value = c(r$total_ss, r$between_ss, r$total_wcss, r$r_squared,
                    ifelse(is.na(r$silhouette_width), NA, r$silhouette_width),
                    r$ch_index)
        )

        # Write both tables to separate sheets would require xlsx package
        # For CSV, concatenate with a separator
        write.csv(stats_df, file, row.names = FALSE)
        write("", file, append = TRUE)
        write("Overall Statistics:", file, append = TRUE)
        write.table(overall_stats, file, append = TRUE, sep = ",",
                    row.names = FALSE, col.names = TRUE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export silhouette data
  output$export_silhouette_csv <- downloadHandler(
    filename = function() paste0("silhouette_data_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results) &&
         !is.null(values$clustering_results$silhouette_data)) {
        sil_data <- values$clustering_results$silhouette_data

        sil_df <- data.frame(
          Subject_ID = 1:nrow(sil_data),
          Cluster = sil_data[, 1],
          Neighbor = sil_data[, 2],
          Silhouette_Width = sil_data[, 3]
        )

        # Add group labels if available
        if(!is.null(values$group_labels)) {
          sil_df$Group <- values$group_labels
        }

        write.csv(sil_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No silhouette data available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster-group composition
  output$export_cluster_group_csv <- downloadHandler(
    filename = function() {
      group_var_name <- get_cluster_group_var_name()
      paste0("cluster_", group_var_name, "_composition_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if(!is.null(values$clustering_results) && !is.null(values$group_labels)) {
        clusters <- values$clustering_results$cluster_assignments
        groups <- get_cluster_group_labels()
        group_var_name <- get_cluster_group_var_name()

        # Create contingency table
        cont_table <- table(Cluster = clusters, Group = groups)

        # Calculate percentages for each cluster
        cluster_pcts <- prop.table(cont_table, margin = 1) * 100

        # Create export dataframe
        export_df <- data.frame(Cluster = rownames(cont_table))

        # Add metadata row with group variable name
        attr(export_df, "group_variable") <- group_var_name

        # Add count and percentage columns for each group
        for(g in colnames(cont_table)) {
          export_df[[paste0(g, "_Count")]] <- cont_table[, g]
          export_df[[paste0(g, "_Percent")]] <- round(cluster_pcts[, g], 2)
        }

        # Add total column
        export_df$Total <- rowSums(cont_table)

        # Add group variable name as comment in first row
        write.csv(export_df, file, row.names = FALSE)

        # Append metadata
        cat(paste0("\n# Group Variable: ", group_var_name, "\n"), file = file, append = TRUE)
      } else {
        write.csv(data.frame(Message = "No group composition data available"), file, row.names = FALSE)
      }
    }
  )

  # Export smoothed curves - Wide format
  output$export_smoothed_csv <- downloadHandler(
    filename = function() paste0("smoothed_curves_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        if(!is.null(values$smooth_data)) {
          # Export smooth_data directly (already smoothed curves)
          df_export <- as.data.frame(values$smooth_data)
          
          # Use original time labels directly (no interpolation!)
          if(!is.null(values$time_labels) && length(values$time_labels) == ncol(values$smooth_data)) {
            colnames(df_export) <- values$time_labels
          } else {
            colnames(df_export) <- paste0("T", 1:ncol(df_export))
          }
          
          # Add subject IDs
          df_export <- cbind(Subject = 1:nrow(df_export), df_export)
          
          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == nrow(df_export)) {
            df_export <- cbind(Group = values$group_labels, df_export)
          }
          
          write.csv(df_export, file, row.names = FALSE)
          cat("Successfully exported smoothed curves (wide format) to:", file, "\n")
          cat("Exported", nrow(df_export), "subjects with", ncol(values$smooth_data), "time points\n")
          
        } else {
          write.csv(data.frame(Message = "No smoothed data available. Please apply smoothing first."), 
                    file, row.names = FALSE)
          cat("No smoothed data available for export\n")
        }
      }, error = function(e) {
        cat("Error in export_smoothed_csv:", e$message, "\n")
        write.csv(data.frame(Error = paste("Export failed:", e$message)), 
                  file, row.names = FALSE)
      })
    }
  )
  
  # Export smoothed curves - Long format
  output$export_smoothed_long_csv <- downloadHandler(
    filename = function() paste0("smoothed_curves_long_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        if(!is.null(values$smooth_data)) {
          n_subjects <- nrow(values$smooth_data)
          n_time <- ncol(values$smooth_data)
          
          # Use original time labels directly (no interpolation!)
          if(!is.null(values$time_labels) && length(values$time_labels) == n_time) {
            time_labels_use <- values$time_labels
          } else {
            time_labels_use <- 1:n_time
          }
          
          # Create normalized time grid
          time_normalized <- seq(0, 1, length.out = n_time)
          
          # Convert to long format
          df_long <- data.frame(
            Subject = rep(1:n_subjects, each = n_time),
            Time = rep(time_labels_use, times = n_subjects),
            Time_Normalized = rep(time_normalized, times = n_subjects),
            Value = as.vector(t(values$smooth_data))
          )
          
          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == n_subjects) {
            df_long$Group <- rep(values$group_labels, each = n_time)
            # Reorder columns
            df_long <- df_long[, c("Subject", "Group", "Time", "Time_Normalized", "Value")]
          }
          
          write.csv(df_long, file, row.names = FALSE)
          cat("Successfully exported smoothed curves (long format) to:", file, "\n")
          cat("Exported", n_subjects, "subjects ×", n_time, "time points =", nrow(df_long), "rows\n")
          
        } else {
          write.csv(data.frame(Message = "No smoothed data available. Please apply smoothing first."), 
                    file, row.names = FALSE)
          cat("No smoothed data available for export\n")
        }
      }, error = function(e) {
        cat("Error in export_smoothed_long_csv:", e$message, "\n")
        write.csv(data.frame(Error = paste("Export failed:", e$message)), 
                  file, row.names = FALSE)
      })
    }
  )
  
  output$export_plots <- downloadHandler(
    filename = function() paste0("all_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$data)) {
        plot(1:10, main = "Plots export not fully implemented")
        text(5, 5, "Export functionality to be implemented", cex = 2)
      }
      
      dev.off()
    }
  )
  
  output$export_fanova_plots <- downloadHandler(
    filename = function() paste0("fanova_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$fanova_results)) {
        res <- values$fanova_results
        par(mfrow = c(2, 2))
        
        matplot(res$time_points, res$group_means, type = "l", 
                col = 1:res$n_groups, lwd = 2,
                xlab = "Time", ylab = "Value",
                main = "Group Mean Functions")
        legend("topright", legend = res$groups, col = 1:res$n_groups, lty = 1, lwd = 2)
        
        plot(res$time_points, res$F_stat, type = "l", col = "blue", lwd = 2,
             xlab = "Time", ylab = "F-statistic",
             main = "Pointwise F-statistics")
        
        plot(res$time_points, res$p_values_adjusted, type = "l", col = "darkgreen", lwd = 2,
             xlab = "Time", ylab = "p-value", log = "y",
             main = "Adjusted p-values")
        
        plot(res$time_points, res$eta_squared, type = "l", col = "purple", lwd = 2,
             xlab = "Time", ylab = "η²",
             main = "Effect Size (η²)")
      }
      
      dev.off()
    }
  )
  
  output$export_pairwise_plots <- downloadHandler(
    filename = function() paste0("pairwise_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$pairwise_results)) {
        par(mfrow = c(2, 2))
        
        for(pair_name in values$pairwise_results$pair_names[1:min(4, length(values$pairwise_results$pair_names))]) {
          pair_result <- values$pairwise_results$results[[pair_name]]
          time_points <- values$pairwise_results$time_points
          
          plot(time_points, pair_result$mean_diff, type = "l", col = "blue", lwd = 2,
               xlab = "Time", ylab = "Mean Difference",
               main = pair_name)
          abline(h = 0, lty = 2)
          
          lines(time_points, pair_result$ci_lower, col = "lightblue", lty = 2)
          lines(time_points, pair_result$ci_upper, col = "lightblue", lty = 2)
        }
      }
      
      dev.off()
    }
  )
  
  # ============================================================================
  # ANALYSIS CODE GENERATION
  # ============================================================================
  # Generate reproducible R code for all performed analyses

  generate_analysis_code <- function(full = TRUE) {
    # Generate comprehensive R code for reproducing the analysis
    # Args:
    #   full: If TRUE, generate complete code; if FALSE, generate preview (truncated)

    code_lines <- c()

    # Helper to add lines
    add <- function(...) {
      code_lines <<- c(code_lines, paste0(...))
    }

    add("# =============================================================================")
    add("# FUNCTIONAL DATA ANALYSIS - REPRODUCIBLE CODE")
    add("# Generated from F*CK")
    add("# Date: ", Sys.Date())
    add("# =============================================================================")
    add("#")
    # AUDIT (P2.2): a script that names no versions is not reproducible, it is
    # only re-runnable. The estimators below come from fda, mgcv and refund, all
    # of which change across releases. The environment this analysis actually ran
    # in is stamped here and CHECKED at the top of the script, so a rerun that
    # would silently produce different numbers says so instead.
    add("# THE ENVIRONMENT THIS ANALYSIS RAN IN")
    add("#   ", R.version.string)
    .pk <- c("fda", "mgcv", "refund", "fda.usc", "minpack.lm", "cluster",
             "ggplot2", "plotly", "dplyr")
    .have <- character(0)
    for (.p in .pk) {
      .v <- tryCatch(as.character(utils::packageVersion(.p)), error = function(e) NA_character_)
      if (!is.na(.v)) { add(sprintf("#   %-12s %s", .p, .v)); .have <- c(.have, sprintf('%s = "%s"', .p, .v)) }
    }
    add("# =============================================================================")
    add("")
    add("# Re-running under different package versions can change these numbers.")
    add("# This reports any that differ; it does not stop the script.")
    add(sprintf("fck_recorded_versions <- c(%s)", paste(.have, collapse = ", ")))
    add("fck_check_env <- function(recorded) {")
    add("  now <- vapply(names(recorded), function(p)")
    add("    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_),")
    add("    character(1))")
    add("  diff <- which(!is.na(now) & now != recorded)")
    add("  if (length(diff)) {")
    add("    message('NOTE: package versions differ from the run that produced these results:')")
    add("    for (i in diff) message(sprintf('  %-12s recorded %s, installed %s',")
    add("                                    names(recorded)[i], recorded[i], now[i]))")
    add("  }")
    add("  miss <- names(recorded)[is.na(now)]")
    add("  if (length(miss)) message('NOTE: not installed: ', paste(miss, collapse = ', '))")
    add("  invisible(NULL)")
    add("}")
    add("fck_check_env(fck_recorded_versions)")
    add("")

    # ---- SECTION 1: LIBRARIES ----
    add("# -----------------------------------------------------------------------------")
    add("# 1. LOAD REQUIRED LIBRARIES")
    add("# -----------------------------------------------------------------------------")
    add("library(fda)        # Functional data analysis")
    add("library(mgcv)       # GAM smoothing with REML")
    add("library(ggplot2)    # Plotting")
    add("library(dplyr)      # Data manipulation")

    if(!is.null(values$fanova_results) && !is.null(values$fanova_results$design) &&
       values$fanova_results$design == "within") {
    }

    if(!is.null(values$clustering_results)) {
      add("library(cluster)    # Clustering diagnostics")
      add("library(fda.usc)    # Functional clustering")
    }
    if(!is.null(values$sofr_model)) {
      add("library(refund)     # Scalar-on-function regression")
    }
    if(!is.null(values$harmonic_model) &&
       identical(values$harmonic_model$trend_type, "exp_sat")) {
      add("library(minpack.lm) # Robust non-linear fits (exponential saturation)")
    }
    add("")

    # ---- SECTION 2: DATA LOADING ----
    add("# -----------------------------------------------------------------------------")
    add("# 2. LOAD AND PREPARE DATA")
    add("# -----------------------------------------------------------------------------")

    if(!is.null(values$data)) {
      n_subj <- nrow(values$data)
      n_time <- ncol(values$data)

      add("# Data dimensions: ", n_subj, " subjects x ", n_time, " time points")
      add("")
      add("# Load your data (replace with your actual file path)")
      add("# The data should be a matrix with subjects in rows and time points in columns")
      add("data_matrix <- read.csv('your_data.csv', row.names = 1)")
      add("data_matrix <- as.matrix(data_matrix)")
      add("")
      add("# Number of subjects and time points")
      add("n_subjects <- nrow(data_matrix)  # ", n_subj)
      add("n_time <- ncol(data_matrix)      # ", n_time)
      add("")
      add("# Create normalized time grid [0, 1]")
      add("time_points <- seq(0, 1, length.out = n_time)")
      add("")

      # Time labels if available
      if(!is.null(values$time_labels)) {
        time_labels_str <- paste0("c('", paste(head(values$time_labels, 5), collapse = "', '"),
                                  if(length(values$time_labels) > 5) "', ..." else "'", ")")
        add("# Original time labels: ", time_labels_str)
        add("")
      }

      # Group labels
      if(!is.null(values$group_labels)) {
        groups <- levels(values$group_labels)
        group_counts <- table(values$group_labels)
        add("# Group structure")
        add("group_labels <- factor(c(")
        # Show first few group assignments
        n_show <- min(10, length(values$group_labels))
        group_sample <- paste0("'", as.character(values$group_labels[1:n_show]), "'", collapse = ", ")
        add("  ", group_sample, if(length(values$group_labels) > 10) ", ..." else "")
        add("))")
        add("# Groups: ", paste(names(group_counts), " (n=", group_counts, ")", sep = "", collapse = ", "))
        add("")
      }
    }

    # ---- SECTION 3: SMOOTHING ----
    add("# -----------------------------------------------------------------------------")
    add("# 3. SMOOTHING / FUNCTIONAL DATA OBJECT CREATION")
    add("# -----------------------------------------------------------------------------")

    if(!is.null(values$smooth_fit_metrics)) {
      metrics <- values$smooth_fit_metrics

      # Determine smoothing method used
      smooth_method <- if(!is.null(input$smooth_method)) input$smooth_method else "auto"

      if(smooth_method == "none") {
        add("# Smoothing method: None (raw data)")
        add("# Data used as-is without smoothing")
        add("")
      } else if(smooth_method == "auto") {
        n_basis <- if(!is.null(input$n_basis)) input$n_basis else 20
        add("# Smoothing method: Automatic (REML optimization)")
        add("# Number of B-spline basis functions: ", n_basis)
        if(!is.null(metrics$lambda)) {
          add("# Estimated lambda (smoothing parameter): ", sprintf("%.6e", metrics$lambda))
        }
        add("")
        add("# Create B-spline basis")
        add("n_basis <- ", n_basis)
        add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
        add("")
        add("# Smooth with REML (lambda = 0 triggers automatic optimization)")
        add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = 0)")
        add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
        add("fd_obj <- smooth_result$fd")
        add("")
      } else if(smooth_method == "manual") {
        n_basis <- if(!is.null(input$n_basis_manual)) input$n_basis_manual else 20
        smooth_factor <- if(!is.null(input$smooth_factor)) input$smooth_factor else 1
        lambda <- 10^(-smooth_factor)
        add("# Smoothing method: Manual")
        add("# Number of B-spline basis functions: ", n_basis)
        add("# Smoothing factor: ", smooth_factor, " (lambda = 10^(-", smooth_factor, ") = ", sprintf("%.6e", lambda), ")")
        add("")
        add("# Create B-spline basis")
        add("n_basis <- ", n_basis)
        add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
        add("")
        add("# Smooth with specified lambda")
        add("lambda <- 10^(-", smooth_factor, ")  # = ", sprintf("%.6e", lambda))
        add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = lambda)")
        add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
        add("fd_obj <- smooth_result$fd")
        add("")
      }

      # Fit statistics
      add("# Smoothing fit statistics:")
      add("#   Mean R²: ", sprintf("%.4f", metrics$mean_r_squared),
          " (SD: ", sprintf("%.4f", metrics$sd_r_squared), ")")
      add("#   Mean RMSE: ", sprintf("%.4f", metrics$mean_rmse),
          " (SD: ", sprintf("%.4f", metrics$sd_rmse), ")")
      add("")

      # Bounds constraint
      if(!is.null(input$constrain_bounds) && input$constrain_bounds) {
        add("# Value bounds constraint applied")
        add("min_bound <- ", input$min_bound)
        add("max_bound <- ", input$max_bound)
        add("# Smoothed values were clamped to [min_bound, max_bound]")
        add("")
      }
    } else {
      add("# No smoothing applied yet")
      add("# Default smoothing code:")
      add("n_basis <- 20")
      add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
      add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = 0)  # REML optimization")
      add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
      add("fd_obj <- smooth_result$fd")
      add("")
    }

    add("# Evaluate smoothed curves")
    add("smooth_curves <- eval.fd(time_points, fd_obj)")
    add("smooth_curves <- t(smooth_curves)  # Back to subjects x time format")
    add("")

    # ---- SECTION 4: TIME WARPING ----
    if(!is.null(values$warping_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 4. TIME WARPING / CURVE REGISTRATION")
      add("# -----------------------------------------------------------------------------")

      warp_method <- values$warping_results$method
      add("# Warping method: ", warp_method)
      add("")

      if(warp_method == "linear_shift") {
        periodic <- if(!is.null(input$periodic_shift)) input$periodic_shift else FALSE
        reference <- if(!is.null(input$shift_reference)) input$shift_reference else "mean"
        allow_dilation <- if(!is.null(input$allow_dilation)) input$allow_dilation else FALSE

        add("# Linear shift alignment parameters:")
        add("#   Reference: ", reference)
        add("#   Periodic: ", periodic)
        add("#   Allow dilation: ", allow_dilation)
        add("")

        add("# Linear shift alignment function")
        add("linear_shift_alignment <- function(fd_obj, periodic = ", periodic,
            ", reference = '", reference, "', time_points = seq(0, 1, length.out = 100)) {")
        add("  n_curves <- ncol(fd_obj$coefs)")
        add("  n_time <- length(time_points)")
        add("  curves <- eval.fd(time_points, fd_obj)")
        add("")
        add("  # Get reference curve")
        add("  ref_curve <- switch(reference,")
        add("    'mean' = rowMeans(curves),")
        add("    'median' = apply(curves, 1, median),")
        add("    curves[,1])")
        add("")
        add("  # Initialize outputs")
        add("  registered_curves <- matrix(NA, n_time, n_curves)")
        add("  warp_functions <- matrix(NA, n_time, n_curves)")
        add("  shifts <- numeric(n_curves)")
        add("")
        add("  for(i in 1:n_curves) {")
        add("    # Find optimal shift using cross-correlation")
        add("    ccf_result <- ccf(curves[,i], ref_curve, lag.max = floor(n_time/4), plot = FALSE)")
        add("    best_lag <- ccf_result$lag[which.max(ccf_result$acf)]")
        add("    shifts[i] <- best_lag / n_time * 0.1")
        add("")
        add("    # Create warping function")
        add("    warp_functions[,i] <- pmin(1, pmax(0, time_points - shifts[i] * 0.5))")
        add("    warp_functions[1,i] <- 0; warp_functions[n_time,i] <- 1")
        add("")
        add("    # Apply warping")
        add("    if(abs(shifts[i]) > 0.001) {")
        add("      registered_curves[,i] <- approx(time_points, curves[,i],")
        add("                                       xout = warp_functions[,i], rule = 2)$y")
        add("    } else {")
        add("      registered_curves[,i] <- curves[,i]")
        add("    }")
        add("  }")
        add("")
        add("  return(list(registered_curves = registered_curves,")
        add("              warp_functions = warp_functions, shifts = shifts))")
        add("}")
        add("")
        add("# Apply alignment")
        add("warp_result <- linear_shift_alignment(fd_obj, periodic = ", periodic, ", reference = '", reference, "')")
        add("registered_curves <- warp_result$registered_curves")
        add("")

      } else if(warp_method == "parametric") {
        family <- if(!is.null(input$parametric_family)) input$parametric_family else "power"
        param_range <- if(!is.null(input$param_range)) input$param_range else c(0.5, 2)

        add("# Parametric warping parameters:")
        add("#   Family: ", family)
        add("#   Parameter range: [", param_range[1], ", ", param_range[2], "]")
        add("")

        add("# Parametric alignment function")
        add("parametric_alignment <- function(fd_obj, family = '", family, "',")
        add("                                  param_range = c(", param_range[1], ", ", param_range[2], "),")
        add("                                  time_points = seq(0, 1, length.out = 100)) {")
        add("  n_curves <- ncol(fd_obj$coefs)")
        add("  n_time <- length(time_points)")
        add("  curves <- eval.fd(time_points, fd_obj)")
        add("  mean_curve <- rowMeans(curves)")
        add("")
        add("  # Warping function based on family")
        add("  warp_func <- function(t, alpha) {")
        add("    switch(family,")
        add("      'power' = pmin(1, pmax(0, t^alpha)),")
        add("      'exponential' = if(abs(alpha-1) < 0.001) t else pmin(1, pmax(0, (exp(alpha*t)-1)/(exp(alpha)-1))),")
        add("      'quadratic' = pmin(1, pmax(0, alpha*t^2 + (1-alpha)*t)),")
        add("      t)")
        add("  }")
        add("")
        add("  registered_curves <- matrix(NA, n_time, n_curves)")
        add("  alpha_values <- numeric(n_curves)")
        add("")
        add("  for(i in 1:n_curves) {")
        add("    # Optimize warping parameter")
        add("    objective <- function(alpha) {")
        add("      warped <- approx(time_points, curves[,i], xout = warp_func(time_points, alpha), rule = 2)$y")
        add("      sum((warped - mean_curve)^2, na.rm = TRUE)")
        add("    }")
        add("    result <- optimize(objective, interval = param_range, tol = 1e-4)")
        add("    alpha_values[i] <- result$minimum")
        add("    registered_curves[,i] <- approx(time_points, curves[,i],")
        add("                                     xout = warp_func(time_points, alpha_values[i]), rule = 2)$y")
        add("  }")
        add("")
        add("  return(list(registered_curves = registered_curves, alpha_values = alpha_values))")
        add("}")
        add("")
        add("# Apply alignment")
        add("warp_result <- parametric_alignment(fd_obj, family = '", family, "')")
        add("registered_curves <- warp_result$registered_curves")
        add("")

      } else if(warp_method == "landmark") {
        add("# Landmark-based alignment")
        add("# (Implementation depends on detected landmarks)")
        add("")
      }

      # Warping fit statistics
      if(!is.null(values$warping_results$fit_statistics)) {
        stats <- values$warping_results$fit_statistics$summary
        add("# Warping fit statistics:")
        add("#   Mean R² (orig vs warped): ", sprintf("%.4f", stats$mean_r_squared))
        add("#   Mean RMSE: ", sprintf("%.4f", stats$mean_rmse))
        add("#   Variance explained by warping: ", sprintf("%.2f%%", stats$variance_explained_by_warping * 100))
        add("#   AIC: ", sprintf("%.2f", stats$AIC))
        add("#   BIC: ", sprintf("%.2f", stats$BIC))
        add("")
      }

      add("# Create fd object from registered curves")
      add("reg_basis <- create.bspline.basis(c(0, 1), nbasis = min(20, n_time - 2))")
      add("reg_fd <- smooth.basis(time_points, registered_curves, fdPar(reg_basis, 2, 0))$fd")
      add("")
    }

    # ---- SECTION 5: PCA ----
    if(!is.null(values$pca_results)) {
      add("# -----------------------------------------------------------------------------")
      if(!is.null(values$warping_results)) {
        add("# 5. TIME-WARPED PRINCIPAL COMPONENT ANALYSIS")
      } else {
        add("# 5. FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
      }
      add("# -----------------------------------------------------------------------------")

      n_comp <- ncol(values$pca_results$scores)
      varprop <- values$pca_results$varprop[1:n_comp]

      add("# Number of components: ", n_comp)
      add("# Variance explained:")
      for(i in 1:min(n_comp, 5)) {
        add("#   PC", i, ": ", sprintf("%.2f%%", varprop[i] * 100),
            " (cumulative: ", sprintf("%.2f%%", sum(varprop[1:i]) * 100), ")")
      }
      add("")

      add("# Perform functional PCA")
      add("n_components <- ", n_comp)
      if(!is.null(values$warping_results)) {
        add("pca_results <- pca.fd(reg_fd, nharm = n_components)")
      } else {
        add("pca_results <- pca.fd(fd_obj, nharm = n_components)")
      }
      add("")
      add("# Extract results")
      add("scores <- pca_results$scores          # PC scores (n_subjects x n_components)")
      add("varprop <- pca_results$varprop        # Proportion of variance explained")
      add("harmonics <- pca_results$harmonics    # Functional principal components")
      add("")
      add("# View variance explained")
      add("print(round(varprop * 100, 2))")
      add("")
      add("# Plot PC loadings")
      add("par(mfrow = c(1, min(3, n_components)))")
      add("for(i in 1:min(3, n_components)) {")
      add("  plot(harmonics[i], main = paste0('PC', i, ' (', round(varprop[i]*100, 1), '%)'),")
      add("       xlab = 'Time', ylab = 'Loading')")
      add("}")
      add("")
    }

    # ---- SECTION 6: CLUSTERING ----
    if(!is.null(values$clustering_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 6. FUNCTIONAL K-MEANS CLUSTERING")
      add("# -----------------------------------------------------------------------------")

      k <- values$clustering_results$k
      add("# Number of clusters: ", k)
      add("")

      add("# Prepare data for clustering")
      if(!is.null(values$warping_results)) {
        add("fdata_obj <- fdata(t(registered_curves), argvals = time_points)")
      } else {
        add("fdata_obj <- fdata(t(eval.fd(time_points, fd_obj)), argvals = time_points)")
      }
      add("")
      add("# Perform functional k-means clustering")
      add("set.seed(123)  # For reproducibility")
      add("k_clusters <- ", k)
      add("kmeans_result <- kmeans.fd(fdata_obj, ncl = k_clusters, metric = 'L2')")
      add("")
      add("# Get cluster assignments")
      add("cluster_assignments <- kmeans_result$cluster")
      add("print(table(cluster_assignments))")
      add("")

      if(!is.null(values$clustering_results$r_squared)) {
        add("# Clustering fit statistics:")
        add("#   R² (variance explained): ", sprintf("%.4f", values$clustering_results$r_squared))
        add("#   Within-cluster SS: ", sprintf("%.4f", values$clustering_results$total_wcss))
        add("#   Between-cluster SS: ", sprintf("%.4f", values$clustering_results$between_ss))
      }
      add("")
    }

    # ---- SECTION 7: FUNCTIONAL ANOVA ----
    if(!is.null(values$fanova_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 7. FUNCTIONAL ANOVA")
      add("# -----------------------------------------------------------------------------")

      design <- if(!is.null(values$fanova_results$design)) values$fanova_results$design else "between"
      n_groups <- values$fanova_results$n_groups

      add("# Design type: ", if(design == "within") "Within-subjects (Repeated Measures)" else "Between-subjects")
      add("# Number of groups: ", n_groups)
      add("")

      # AUDIT (P2.1): this section used to WRITE OUT a reimplementation. For the
      # between-subjects design it emitted
      #     aov_result <- summary(aov(curves_eval[, t] ~ group_labels))
      # and took the parametric p-value, while the app runs a permutation test
      # with FDR adjustment; for repeated measures it emitted a commented
      # sketch of an rmfanova call that matches no version of that API and that
      # the app never made. So "Export analysis code" produced a script whose
      # numbers could differ from the ones on screen -- and after the P0.4
      # residual-SS correction it diverged further, because the exported aov()
      # was closer to the app's OLD behaviour than to its current one.
      #
      # The harmonic export already solved this: it deparses the actual fitting
      # function into the script, so the exported code IS the app's estimator.
      # The same pattern applies here. deparse() of the live closure cannot
      # drift from the implementation, because it is the implementation.
      add("# The estimator below is the app's own function, written out verbatim.")
      add("# It is not a reconstruction: deparse() of the live function object")
      add("# guarantees the exported script computes what the app computed.")
      add("")
      if(design == "between") {
        add(paste("perform_functional_anova <-",
                  paste(deparse(perform_functional_anova), collapse = "\n")))
        add("")
        if(!is.null(values$warping_results)) {
          add("fd_to_analyze <- reg_fd")
        } else {
          add("fd_to_analyze <- fd_obj")
        }
        add(sprintf("fanova_results <- perform_functional_anova(fd_to_analyze, group_labels, n_permutations = %d, alpha = %s)",
                    values$fanova_results$n_permutations %||% 5000,
                    format(values$fanova_results$alpha %||% 0.05)))
        add("")
        add("cat('Significant time points (FDR-adjusted):',")
        add("    sum(fanova_results$p_values_adjusted < fanova_results$alpha, na.rm = TRUE),")
        add("    'of', length(fanova_results$p_values_adjusted), '\\n')")
        add("")

      } else {
        add("# Pointwise repeated-measures functional ANOVA with within-subject,")
        add("# curve-wise permutation. This does NOT use the rmfanova package: the")
        add("# app never called that API successfully and no longer pretends to.")
        add("")
        add(paste("perform_rm_fanova <-",
                  paste(deparse(perform_rm_fanova), collapse = "\n")))
        add("")
        add(sprintf("fanova_results <- perform_rm_fanova(fd_obj, subject_id, rm_factor, n_permutations = %d, alpha = %s)",
                    values$fanova_results$n_permutations %||% 5000,
                    format(values$fanova_results$alpha %||% 0.05)))
        add("")
      }

      if(!is.null(values$fanova_results$global_p)) {
        add("# For reference, the values this script should reproduce:")
        add("#   Global F-statistic: ", sprintf("%.4f", values$fanova_results$global_f))
        add("#   Global p-value: ", sprintf("%.6f", values$fanova_results$global_p))
        add("")
      }
    }

    # ---- SECTION 8: PAIRWISE COMPARISONS ----
    if(!is.null(values$pairwise_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 8. PAIRWISE COMPARISONS")
      add("# -----------------------------------------------------------------------------")

      correction <- values$pairwise_results$correction_method
      add("# Multiple testing correction: ", correction)
      add("")

      add("# Perform pairwise functional t-tests")
      add("groups <- levels(group_labels)")
      add("n_groups <- length(groups)")
      add("pairs <- combn(groups, 2, simplify = FALSE)")
      add("")
      add("pairwise_results <- list()")
      add("for(pair in pairs) {")
      add("  g1 <- pair[1]; g2 <- pair[2]")
      add("  idx1 <- which(group_labels == g1)")
      add("  idx2 <- which(group_labels == g2)")
      add("")
      add("  # Pointwise t-tests")
      add("  p_values <- numeric(n_time_eval)")
      add("  for(t in 1:n_time_eval) {")
      add("    tt <- t.test(curves_eval[idx1, t], curves_eval[idx2, t])")
      add("    p_values[t] <- tt$p.value")
      add("  }")
      add("")
      add("  # Apply correction")
      add("  p_adjusted <- p.adjust(p_values, method = '", tolower(correction), "')")
      add("")
      add("  pairwise_results[[paste(g1, 'vs', g2)]] <- list(")
      add("    p_values = p_values,")
      add("    p_adjusted = p_adjusted,")
      add("    sig_times = eval_points[p_adjusted < 0.05]")
      add("  )")
      add("}")
      add("")

      # Summary of significant comparisons
      if(!is.null(values$pairwise_results$results)) {
        n_sig <- sum(sapply(values$pairwise_results$results, function(x)
          length(x$sig_times) > 0))
        add("# Number of significant pairwise comparisons: ", n_sig, "/", length(values$pairwise_results$results))
      }
      add("")
    }

    # ---- SECTION 9: VISUALIZATION ----
    add("# -----------------------------------------------------------------------------")
    add("# 9. VISUALIZATION CODE")
    add("# -----------------------------------------------------------------------------")
    add("")
    add("# Plot all curves with mean")
    add("matplot(time_points, t(smooth_curves), type = 'l', col = 'gray70', lty = 1,")
    add("        xlab = 'Time', ylab = 'Value', main = 'Functional Data')")
    add("lines(time_points, colMeans(smooth_curves), col = 'red', lwd = 3)")
    add("")

    if(!is.null(values$group_labels)) {
      add("# Plot by group")
      add("library(ggplot2)")
      add("plot_data <- data.frame(")
      add("  time = rep(time_points, n_subjects),")
      add("  value = as.vector(t(smooth_curves)),")
      add("  subject = rep(1:n_subjects, each = n_time),")
      add("  group = rep(group_labels, each = n_time)")
      add(")")
      add("")
      add("ggplot(plot_data, aes(x = time, y = value, group = subject, color = group)) +")
      add("  geom_line(alpha = 0.5) +")
      add("  stat_summary(aes(group = group), fun = mean, geom = 'line', linewidth = 2) +")
      add("  theme_minimal() +")
      add("  labs(title = 'Functional Data by Group', x = 'Time', y = 'Value')")
      add("")
    }

    # ---- SECTION 10: HARMONIC (COSINOR) REGRESSION ----
    # MERGED APP: CIRCAREG had no code generator, so these three sections are
    # new. Settings come from the fit itself (fck_settings), never from the
    # live widgets, so the script matches the model that is on screen.
    if(!is.null(values$harmonic_model)) {
      hm <- values$harmonic_model
      st <- hm$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 10. HARMONIC (COSINOR) REGRESSION")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("period        <- %s", hm$period))
      add(sprintf("n_harmonics   <- %s", hm$n_harmonics))
      add(sprintf('trend_type    <- "%s"', hm$trend_type))
      add(sprintf("time_vec      <- c(%s)",
                  paste(signif(as.numeric(hm$time_vec), 10), collapse = ", ")))
      add(sprintf("cosinor_input <- %s   # %s",
                  if(isTRUE(hm$using_smoothed)) "smooth_curves" else "raw_data",
                  if(isTRUE(hm$using_smoothed)) "the smoothed curves"
                  else "the raw data (no smoothing had been applied)"))
      add("")
      if(!is.null(st)) {
        add(sprintf("use_bounds    <- %s", isTRUE(st$use_bounds)))
        for(nm in c("mesor_min", "mesor_max", "amplitude_min", "amplitude_max",
                    "A_sat_min", "A_sat_max", "tau_min", "tau_max")) {
          v <- st[[nm]]
          add(sprintf("%-13s <- %s", nm,
                      if(is.null(v) || !is.finite(v)) "NA" else format(v)))
        }
        add("")
      }
      add("# The app's own fitting function, emitted verbatim (not a")
      add("# re-implementation) so this script reproduces its numbers exactly:")
      add(paste("fit_cosinor <-", paste(deparse(fit_cosinor), collapse = "\n")))
      add("")
      if(identical(hm$trend_type, "exp_sat") && exists("fit_cosinor_nonlinear")) {
        add(paste("fit_cosinor_nonlinear <-",
                  paste(deparse(fit_cosinor_nonlinear), collapse = "\n")))
        add("")
      }
      add("# The same per-subject call the app makes:")
      add("cosinor_fits <- lapply(seq_len(nrow(cosinor_input)), function(i) {")
      add("  fit_cosinor(time_vec, cosinor_input[i, ], period = period,")
      add("              n_harmonics = n_harmonics, trend_type = trend_type,")
      if(!is.null(st)) {
        add("              use_bounds = use_bounds,")
        add("              mesor_min = mesor_min, mesor_max = mesor_max,")
        add("              amplitude_min = amplitude_min, amplitude_max = amplitude_max,")
        add("              A_sat_min = A_sat_min, A_sat_max = A_sat_max,")
        add("              tau_min = tau_min, tau_max = tau_max)")
      } else {
        add("              use_bounds = FALSE)")
      }
      add("})")
      add("")
      add("# MESOR / amplitude / acrophase per subject (successful fits only):")
      add("ok <- vapply(cosinor_fits, function(f) isTRUE(f$success), logical(1))")
      add("cosinor_params <- do.call(rbind, lapply(which(ok), function(i) {")
      add("  f <- cosinor_fits[[i]]")
      add("  data.frame(subject = i, mesor = f$mesor,")
      add("             amplitude = f$amplitude[1], acrophase = f$acrophase[1],")
      add("             r_squared = f$r_squared)")
      add("}))")
      add("print(head(cosinor_params))")
      if(!is.null(hm$individual_fits)) {
        add(sprintf("# The app fitted %d of %d subjects successfully.",
                    sum(vapply(hm$individual_fits,
                               function(f) isTRUE(f$success), logical(1))),
                    length(hm$individual_fits)))
      }
      add("")
      add("# NOTE: acrophase is CIRCULAR. Do not average it arithmetically --")
      add("# average the (cos, sin) components, as the app's polar plot does.")
      add("")
    }

    # ---- SECTION 11: FUNCTION-ON-SCALAR REGRESSION ----
    if(!is.null(values$reg_model) && !is.null(values$reg_model$fck_settings)) {
      st <- values$reg_model$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 11. FUNCTION-ON-SCALAR REGRESSION (FoSR)")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("fosr_predictors <- c(%s)",
                  paste0('"', st$predictors, '"', collapse = ", ")))
      add(sprintf("Y <- %s   # %d subjects x %d time points",
                  if(isTRUE(st$using_smoothed)) "smooth_curves" else "raw_data",
                  st$n_subjects, st$n_time))
      add("keep   <- complete.cases(covariates[, fosr_predictors, drop = FALSE])")
      add("df_reg <- covariates[keep, , drop = FALSE]")
      add("Y      <- Y[keep, , drop = FALSE]")
      add("")
      if(identical(st$method, "OLS_nosmooth")) {
        add("# Pointwise OLS: one regression per time point, in closed form.")
        add("f           <- as.formula(paste('~', paste(fosr_predictors, collapse = ' + ')))")
        add("X           <- model.matrix(f, data = df_reg)")
        add("xtx_inv     <- solve(crossprod(X))")
        add("projector   <- xtx_inv %*% t(X)")
        add("beta_hat    <- projector %*% Y          # the coefficient curves")
        add("fitted_vals <- X %*% beta_hat")
        add("residuals   <- Y - fitted_vals")
        add("n <- nrow(Y); p <- ncol(X)")
        add("sigma2 <- colSums(residuals^2) / (n - p)")
        add("")
        if(isTRUE(st$use_bootstrap)) {
          add(sprintf("# Residual bootstrap for the coefficient bands (B = %d).", st$n_boot))
          add("# The app does not fix a seed, so its bands differ slightly run to run;")
          add("# this script fixes one so the script itself is reproducible.")
          add("set.seed(1)")
          add(sprintf("B <- %d", st$n_boot))
          add("boot_betas <- array(NA, dim = c(B, nrow(beta_hat), ncol(beta_hat)))")
          add("for (b in 1:B) {")
          add("  resid_idx <- sample(seq_len(n), n, replace = TRUE)")
          add("  boot_betas[b, , ] <- projector %*% (fitted_vals + residuals[resid_idx, ])")
          add("}")
          add("boot_ci_lower <- apply(boot_betas, c(2, 3), quantile, 0.025, na.rm = TRUE)")
          add("boot_ci_upper <- apply(boot_betas, c(2, 3), quantile, 0.975, na.rm = TRUE)")
          add("")
        }
      } else {
        add("# Smoothed OLS: one GAM over the (subject, time) long form, with a")
        add("# penalised spline in time and a by-smooth per predictor.")
        add("long_data <- data.frame(")
        add("  value   = as.vector(t(Y)),")
        add("  time    = rep(seq(0, 1, length.out = ncol(Y)), times = nrow(Y)),")
        add("  subject = rep(seq_len(nrow(Y)), each = ncol(Y)))")
        add("for (v in fosr_predictors) long_data[[v]] <- rep(df_reg[[v]], each = ncol(Y))")
        add("gam_formula <- as.formula(paste('value ~ s(time) +',")
        add("  paste(sprintf('s(time, by = %s)', fosr_predictors), collapse = ' + '),")
        add("  '+', paste(fosr_predictors, collapse = ' + ')))")
        add("gam_fit <- mgcv::gam(gam_formula, data = long_data, method = 'REML')")
        add("summary(gam_fit)")
        add("")
        add("# NOTE: the app turns this fit into coefficient curves by prediction")
        add("# contrasts (predict at x = 1 minus predict at x = 0), which is exact")
        add("# for numeric predictors and an approximation for factors.")
        add("")
      }
      add(sprintf("# Method as fitted: %s", values$reg_model$method))
      add("")
    }

    # ---- SECTION 12: SCALAR-ON-FUNCTION REGRESSION ----
    if(!is.null(values$sofr_model) && !is.null(values$sofr_model$fck_settings)) {
      st <- values$sofr_model$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 12. SCALAR-ON-FUNCTION REGRESSION (SoFR)")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("X_func <- %s   # the curves, as the functional predictor",
                  if(isTRUE(st$using_smoothed)) "smooth_curves" else "raw_data"))
      add(sprintf('y      <- covariates[["%s"]]', st$response))
      if(length(st$predictors)) {
        add(sprintf("sofr_scalars <- c(%s)",
                    paste0('"', st$predictors, '"', collapse = ", ")))
      }
      add("")
      add("keep   <- complete.cases(y) & complete.cases(X_func)")
      add("X_func <- X_func[keep, , drop = FALSE]; y <- y[keep]")
      add("")
      add("pfr_data <- list(X_func = X_func, y = y)")
      if(length(st$predictors)) {
        add("for (v in sofr_scalars) pfr_data[[v]] <- covariates[keep, v]")
      }
      add(sprintf("sofr_fit <- refund::pfr(%s,", st$formula))
      add(sprintf("                        data = pfr_data, family = %s(link = '%s'))",
                  st$family, st$link))
      add("summary(sofr_fit)")
      add("")
      add(sprintf("# Fitted on %d observations; %s family, %s link.",
                  st$n_obs, st$family, st$link))
      add("")
    }

    add("# =============================================================================")
    add("# END OF ANALYSIS CODE")
    add("# =============================================================================")

    return(paste(code_lines, collapse = "\n"))
  }

  # Download handler for analysis code
  output$export_code <- downloadHandler(
    filename = function() paste0("fda_analysis_code_", Sys.Date(), ".R"),
    content = function(file) {
      code <- generate_analysis_code(full = TRUE)
      writeLines(code, file)
    }
  )

  # Code preview (truncated version)
  output$code_preview <- renderText({
    tryCatch({
      code <- generate_analysis_code(full = FALSE)

      # For preview, show first ~80 lines
      lines <- strsplit(code, "\n")[[1]]
      if(length(lines) > 80) {
        preview <- c(
          lines[1:80],
          "",
          "# ... (truncated for preview)",
          paste0("# Full code: ", length(lines), " lines"),
          "# Download the R file for complete code"
        )
        paste(preview, collapse = "\n")
      } else {
        code
      }
    }, error = function(e) {
      paste("# Error generating code preview:", e$message)
    })
  })

  output$export_fosr_coefs_csv <- downloadHandler(
    filename = function() paste0("beta_coefficients_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$reg_model)
      write.csv(values$reg_model$beta.hat, file)
    }
  )

  output$export_harmonic_params <- downloadHandler(
    filename = function() paste0("harmonic_individual_params_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$harmonic_model)
      write.csv(values$harmonic_model$individual_params, file, row.names = FALSE)
    }
  )
  
  output$export_harmonic_summary <- downloadHandler(
    filename = function() paste0("harmonic_summary_", Sys.Date(), ".txt"),
    content = function(file) {
      req(values$harmonic_model)
      mod <- values$harmonic_model
      
      sink(file)
      cat("=== Harmonic Regression (Cosinor Analysis) Summary ===\n")
      cat("Generated:", as.character(Sys.time()), "\n\n")
      cat("Period:", mod$period, "\n")
      cat("Number of harmonics:", mod$n_harmonics, "\n")
      cat("Number of subjects:", length(mod$individual_fits), "\n\n")
      
      if(!is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        cat("--- Population Mean Parameters ---\n")
        cat(sprintf("MESOR:     %.4f\n", pop$mean_mesor))
        cat(sprintf("Amplitude: %.4f\n", pop$mean_amplitude))
        # AUDIT: acrophases are exported in CLOCK time, through the same helper
        # the on-screen report uses, with the model-elapsed value alongside so
        # the export is self-explaining.
        cat(sprintf("Acrophase: %s clock  (%.4f rad, %.2f h model-elapsed, origin %s)\n",
                    fck_acrophase_label(hours = pop$mean_acrophase_time,
                                        period = mod$period, harmonic = 1,
                                        clock_origin = fck_clock_origin(mod)),
                    pop$mean_acrophase_rad, pop$mean_acrophase_time,
                    fck_clock_label(fck_clock_origin(mod), mod$period, show_day = FALSE)))
        cat(sprintf("Rayleigh Z: %.4f (p = %.4f)\n", pop$rayleigh_z, pop$rayleigh_p))
        cat("\n")
      }
      
      cat("--- Individual Parameter Statistics ---\n")
      params <- mod$individual_params
      cat(sprintf("MESOR:     Mean=%.4f, SD=%.4f\n", mean(params$mesor), sd(params$mesor)))
      cat(sprintf("Amplitude (H1): Mean=%.4f, SD=%.4f\n", mean(params$amplitude_1), sd(params$amplitude_1)))
      cat(sprintf("Acrophase (H1): arithmetic mean %s clock (%.2f h model-elapsed)\n",
                  fck_acrophase_label(hours = mean(params$acrophase_time_1),
                                      period = mod$period, harmonic = 1,
                                      clock_origin = fck_clock_origin(mod), all = FALSE),
                  mean(params$acrophase_time_1)))
      cat(sprintf("R-squared: Mean=%.4f, Range=[%.4f, %.4f]\n", 
                  mean(params$r_squared), min(params$r_squared), max(params$r_squared)))
      cat(sprintf("Significant rhythms (p<0.05): %d / %d\n",
                  sum(params$p_value < 0.05), nrow(params)))
      
      if(!is.null(mod$boot_results)) {
        cat(sprintf("\n--- Bootstrap 95%% CIs (B=%d) ---\n", mod$boot_results$B))
        cat(sprintf("MESOR:     [%.4f, %.4f]\n", 
                    mod$boot_results$mesor_ci[1], mod$boot_results$mesor_ci[2]))
        cat(sprintf("Amplitude: [%.4f, %.4f]\n", 
                    mod$boot_results$amplitude_ci[1], mod$boot_results$amplitude_ci[2]))
        cat(sprintf("Acrophase: [%s, %s] clock\n",
                    fck_acrophase_label(hours = mod$boot_results$acrophase_ci[1],
                                        period = mod$period, harmonic = 1,
                                        clock_origin = fck_clock_origin(mod), all = FALSE),
                    fck_acrophase_label(hours = mod$boot_results$acrophase_ci[2],
                                        period = mod$period, harmonic = 1,
                                        clock_origin = fck_clock_origin(mod), all = FALSE)))
      }
      sink()
    }
  )
