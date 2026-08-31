# ==========================================================================
# server/60_clustering.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 6978-8912  (functional clustering, optimisation, DCF, outputs)
# ==========================================================================
  # ============================================================================
  # FUNCTIONAL CLUSTERING
  # ============================================================================

  # Cluster optimization analysis (Elbow & Silhouette)
  observeEvent(input$run_cluster_optimization, {
    req(values$data)

    tryCatch({
      # Determine which data to use
      if(input$opt_data_type == "smoothed") {
        if(is.null(values$smooth_data)) {
          showNotification("No smoothed data available. Please apply smoothing first or use raw data.",
                           type = "error", duration = 5)
          return(NULL)
        }
        data_matrix <- values$smooth_data
      } else {
        data_matrix <- values$data
      }

      # Apply standardization if requested (mirrors the clustering module setting)
      standardize_opt <- if(!is.null(input$cluster_standardize)) input$cluster_standardize else "none"
      if(standardize_opt == "within") {
        data_matrix <- t(apply(data_matrix, 1, function(x) {
          s <- sd(x, na.rm = TRUE)
          if(is.na(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
        }))
      } else if(standardize_opt == "between") {
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        col_sds   <- apply(data_matrix, 2, sd, na.rm = TRUE)
        col_sds[col_sds == 0] <- 1
        data_matrix <- sweep(sweep(data_matrix, 2, col_means, "-"), 2, col_sds, "/")
      }

      # Check for functional k-means requirements
      if(input$opt_kmeans_method == "functional") {
        if(input$opt_data_type != "smoothed") {
          showNotification("Functional K-Means requires smoothed data. Please select 'Smoothed Data' or use Standard method.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(is.null(values$fd_obj)) {
          showNotification("No functional data object available. Please apply smoothing first.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(!requireNamespace("fda.usc", quietly = TRUE)) {
          showNotification("fda.usc package not available. Please install it or use Standard K-Means.",
                           type = "error", duration = 5)
          return(NULL)
        }
      }

      # Check for missing values
      if(any(is.na(data_matrix))) {
        showNotification("Data contains missing values. Please handle missing values before clustering.",
                         type = "error", duration = 5)
        return(NULL)
      }

      max_k <- input$max_clusters_test

      # Initialize storage for metrics
      k_values <- 2:max_k
      wcss_values <- numeric(length(k_values))
      silhouette_values <- numeric(length(k_values))
      ch_values <- numeric(length(k_values))

      # Show progress
      progress_msg <- switch(input$opt_kmeans_method,
        "functional"   = 'Testing functional cluster solutions...',
        "hierarchical" = 'Running hierarchical clustering...',
        'Testing standard cluster solutions...'
      )

      # For hierarchical: compute the dendrogram once before the loop (efficient)
      hc_opt <- NULL
      if(input$opt_kmeans_method == "hierarchical") {
        opt_linkage <- if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2"
        dist_matrix_hc_opt <- dist(data_matrix)
        hc_opt <- hclust(dist_matrix_hc_opt, method = opt_linkage)
      }

      withProgress(message = progress_msg, value = 0, {
        for(i in seq_along(k_values)) {
          k <- k_values[i]

          # Update progress
          incProgress(1/length(k_values), detail = paste("Testing k =", k))

          # Run clustering based on selected method
          set.seed(123)

          opt_nstart   <- if(!is.null(input$kmeans_nstart)) input$kmeans_nstart else 10
          opt_iter_max <- if(!is.null(input$kmeans_iter))  input$kmeans_iter  else 100

          if(input$opt_kmeans_method == "functional") {
            # ===== FUNCTIONAL K-MEANS using fda.usc (with restarts) =====
            time_grid <- seq(0, 1, length.out = ncol(data_matrix))
            fdata_obj <- fda.usc::fdata(mdata = data_matrix, argvals = time_grid)

            best_fkm_opt  <- NULL
            best_wcss_opt <- Inf
            for(restart in seq_len(opt_nstart)) {
              set.seed(100 + restart)
              fkm_try <- tryCatch(
                fda.usc::kmeans.fd(fdata_obj, ncl = k, max.iter = opt_iter_max),
                error = function(e) NULL
              )
              if(is.null(fkm_try)) next
              asn <- fkm_try$cluster
              wcss_try <- sum(sapply(seq_len(k), function(j) {
                idx <- which(asn == j)
                if(length(idx) == 0) return(0)
                m <- data_matrix[idx, , drop = FALSE]
                sum(apply(m, 1, function(x) sum((x - colMeans(m))^2)))
              }))
              if(wcss_try < best_wcss_opt) { best_wcss_opt <- wcss_try; best_fkm_opt <- fkm_try }
            }
            fkm_result <- best_fkm_opt

            cluster_assignments <- fkm_result$cluster

            # Calculate cluster means
            cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))

            for(j in 1:k) {
              cluster_idx <- which(cluster_assignments == j)
              if(length(cluster_idx) > 0) {
                cluster_fd <- values$fd_obj[cluster_idx]
                cluster_mean_fd <- mean.fd(cluster_fd)
                cluster_means[j, ] <- eval.fd(time_grid, cluster_mean_fd)
              }
            }

            # Calculate WCSS
            wcss_per_cluster <- numeric(k)
            for(j in 1:k) {
              cluster_members <- data_matrix[cluster_assignments == j, , drop = FALSE]
              if(nrow(cluster_members) > 0) {
                cluster_center <- cluster_means[j, ]
                wcss_per_cluster[j] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
              }
            }

            total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
            within_ss <- sum(wcss_per_cluster)
            between_ss <- total_ss - within_ss

            # Store WCSS
            wcss_values[i] <- within_ss

          } else if(input$opt_kmeans_method == "hierarchical") {
            # ===== HIERARCHICAL CLUSTERING (tree already built above) =====
            cluster_assignments <- cutree(hc_opt, k = k)

            # Cluster means and WCSS
            cluster_means_opt <- matrix(0, nrow = k, ncol = ncol(data_matrix))
            for(j in 1:k) {
              idx <- which(cluster_assignments == j)
              if(length(idx) > 0)
                cluster_means_opt[j, ] <- colMeans(data_matrix[idx, , drop = FALSE])
            }
            wcss_per_cluster_opt <- numeric(k)
            for(j in 1:k) {
              members <- data_matrix[cluster_assignments == j, , drop = FALSE]
              if(nrow(members) > 0)
                wcss_per_cluster_opt[j] <- sum(apply(members, 1, function(x) sum((x - cluster_means_opt[j, ])^2)))
            }
            total_ss_opt <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
            within_ss  <- sum(wcss_per_cluster_opt)
            between_ss <- total_ss_opt - within_ss
            wcss_values[i] <- within_ss

          } else {
            # ===== STANDARD K-MEANS =====
            kmeans_result <- kmeans(data_matrix, centers = k, nstart = opt_nstart, iter.max = opt_iter_max)
            cluster_assignments <- kmeans_result$cluster

            # Store WCSS
            wcss_values[i] <- kmeans_result$tot.withinss

            # Calculate SS for CH index
            between_ss <- kmeans_result$betweenss
            within_ss <- kmeans_result$tot.withinss
          }

          # Calculate silhouette width (common for all methods)
          if(requireNamespace("cluster", quietly = TRUE)) {
            dist_matrix <- dist(data_matrix)
            sil <- cluster::silhouette(cluster_assignments, dist_matrix)
            silhouette_values[i] <- mean(sil[, 3])
          } else {
            silhouette_values[i] <- NA
          }

          # Calculate Calinski-Harabasz index (common for both methods)
          n <- nrow(data_matrix)
          ch_values[i] <- (between_ss / (k - 1)) / (within_ss / (n - k))
        }
      })

      # Store results
      values$cluster_optimization <- list(
        k_values = k_values,
        wcss = wcss_values,
        silhouette = silhouette_values,
        calinski_harabasz = ch_values,
        data_type = input$opt_data_type,
        method = input$opt_kmeans_method,
        hc_obj = hc_opt  # hclust object for dendrogram (NULL for non-hierarchical)
      )

      method_label <- switch(input$opt_kmeans_method,
        "functional"   = "Functional K-Means",
        "hierarchical" = paste0("Hierarchical Clustering (", if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2", ")"),
        "Standard K-Means"
      )
      showNotification(paste0("Cluster optimization analysis completed using ", method_label, "!"),
                       type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Optimization error:", e$message),
                       type = "error", duration = 5)
      values$cluster_optimization <- NULL
    })
  })

  # Optimization status output
  output$optimization_status <- renderText({
    if(is.null(values$cluster_optimization)) {
      "Not run"
    } else {
      paste0("✓ Tested k=2 to k=", max(values$cluster_optimization$k_values))
    }
  })

  # Elbow plot
  output$elbow_plot <- renderPlotly({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Highlight current k if clustering has been run
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Create elbow plot
    p <- plot_ly() %>%
      add_trace(
        x = opt$k_values,
        y = opt$wcss,
        type = 'scatter',
        mode = 'lines+markers',
        line = list(color = '#377EB8', width = 2),
        marker = list(size = 8),
        name = 'WCSS',
        hovertemplate = paste0(
          "k = %{x}",
          "<br>WCSS = %{y:.0f}",
          "<extra></extra>"
        )
      )

    # Add vertical line for current k
    if(current_k %in% opt$k_values) {
      p <- p %>%
        add_trace(
          x = c(current_k, current_k),
          y = c(0, max(opt$wcss)),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'red', width = 2, dash = 'dash'),
          name = paste0('Current k=', current_k),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )
    }

    # Add method label to title
    method_label <- if(!is.null(opt$method) && opt$method == "functional") {
      " - Functional K-Means"
    } else {
      " - Standard K-Means"
    }

    p <- p %>%
      layout(
        title = paste0("Elbow Method (Within-Cluster Sum of Squares)", method_label),
        xaxis = list(title = "Number of Clusters (k)"),
        yaxis = list(title = "Total Within-Cluster SS"),
        hovermode = "closest"
      )

    p
  })

  # Silhouette plot
  output$silhouette_plot <- renderPlotly({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Highlight current k if clustering has been run
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Create silhouette plot
    p <- plot_ly() %>%
      add_trace(
        x = opt$k_values,
        y = opt$silhouette,
        type = 'scatter',
        mode = 'lines+markers',
        line = list(color = '#4DAF4A', width = 2),
        marker = list(size = 8),
        name = 'Silhouette Score',
        hovertemplate = paste0(
          "k = %{x}",
          "<br>Silhouette = %{y:.3f}",
          "<extra></extra>"
        )
      )

    # Add vertical line for current k
    if(current_k %in% opt$k_values) {
      p <- p %>%
        add_trace(
          x = c(current_k, current_k),
          y = c(min(opt$silhouette, na.rm = TRUE), max(opt$silhouette, na.rm = TRUE)),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'red', width = 2, dash = 'dash'),
          name = paste0('Current k=', current_k),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )
    }

    # Add method label to title
    method_label <- if(!is.null(opt$method) && opt$method == "functional") {
      " - Functional K-Means"
    } else {
      " - Standard K-Means"
    }

    p <- p %>%
      layout(
        title = paste0("Silhouette Analysis (Higher is Better)", method_label),
        xaxis = list(title = "Number of Clusters (k)"),
        yaxis = list(title = "Average Silhouette Width"),
        hovermode = "closest"
      )

    p
  })

  # Optimization recommendation
  output$optimization_recommendation <- renderText({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Find optimal k based on different criteria

    # 1. Elbow method - look for the "elbow" (maximum second derivative)
    wcss_diff <- diff(opt$wcss)
    wcss_diff2 <- diff(wcss_diff)
    elbow_k <- opt$k_values[which.max(abs(wcss_diff2)) + 1]

    # 2. Maximum silhouette
    if(!all(is.na(opt$silhouette))) {
      max_sil_k <- opt$k_values[which.max(opt$silhouette)]
      max_sil_score <- max(opt$silhouette, na.rm = TRUE)
    } else {
      max_sil_k <- NA
      max_sil_score <- NA
    }

    # 3. Maximum Calinski-Harabasz
    max_ch_k <- opt$k_values[which.max(opt$calinski_harabasz)]

    # Current k
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Build recommendation text
    output_text <- "=== OPTIMAL CLUSTER NUMBER RECOMMENDATIONS ===\n\n"

    output_text <- paste0(output_text, "Elbow Method:\n")
    output_text <- paste0(output_text, "  Suggested k = ", elbow_k, "\n")
    output_text <- paste0(output_text, "  (Look for the 'elbow' point where WCSS starts to level off)\n\n")

    if(!is.na(max_sil_k)) {
      output_text <- paste0(output_text, "Silhouette Analysis:\n")
      output_text <- paste0(output_text, "  Optimal k = ", max_sil_k,
                            " (Silhouette = ", round(max_sil_score, 3), ")\n")
      if(max_sil_score > 0.7) {
        output_text <- paste0(output_text, "  Strong clustering structure\n\n")
      } else if(max_sil_score > 0.5) {
        output_text <- paste0(output_text, "  Reasonable clustering structure\n\n")
      } else if(max_sil_score > 0.25) {
        output_text <- paste0(output_text, "  Weak clustering structure\n\n")
      } else {
        output_text <- paste0(output_text, "  No substantial clustering structure\n\n")
      }
    }

    output_text <- paste0(output_text, "Calinski-Harabasz Index:\n")
    output_text <- paste0(output_text, "  Optimal k = ", max_ch_k, "\n")
    output_text <- paste0(output_text, "  (Higher values indicate better-defined clusters)\n\n")

    # Current k assessment
    if(current_k %in% opt$k_values) {
      output_text <- paste0(output_text, "Your Current Selection (k=", current_k, "):\n")
      k_idx <- which(opt$k_values == current_k)
      output_text <- paste0(output_text, "  WCSS: ", round(opt$wcss[k_idx], 2), "\n")
      if(!is.na(opt$silhouette[k_idx])) {
        output_text <- paste0(output_text, "  Silhouette: ", round(opt$silhouette[k_idx], 3), "\n")
      }
      output_text <- paste0(output_text, "  Calinski-Harabasz: ", round(opt$calinski_harabasz[k_idx], 2), "\n\n")
    }

    output_text <- paste0(output_text, "INTERPRETATION:\n")
    output_text <- paste0(output_text, "- The elbow plot shows diminishing returns in variance explained\n")
    output_text <- paste0(output_text, "- Silhouette scores range from -1 to 1 (>0.5 is good)\n")
    output_text <- paste0(output_text, "- Consider domain knowledge and interpretability when choosing k\n")

    output_text
  })

  # ============================================================================
  # DCF (DENSITY CORE FINDING) CLUSTERING INTEGRATION
  # ============================================================================

  # Check DCF Python setup
  observeEvent(input$check_dcf_setup, {
    tryCatch({
      # Check if reticulate is available
      if(!requireNamespace("reticulate", quietly = TRUE)) {
        showNotification("reticulate package not installed. Please install it with: install.packages('reticulate')",
                         type = "error", duration = 10)
        return()
      }

      # Check Python availability
      python_available <- reticulate::py_available(initialize = TRUE)
      if(!python_available) {
        showNotification("Python not found. Please install Python and configure reticulate.",
                         type = "error", duration = 10)
        return()
      }

      python_version <- reticulate::py_config()$version
      showNotification(paste("Python found:", python_version), type = "message", duration = 3)

      # Check if DCFcluster is installed
      dcf_available <- reticulate::py_module_available("DCFcluster")

      if(dcf_available) {
        showNotification("DCFcluster Python package is installed and ready!",
                         type = "message", duration = 5)
      } else {
        showNotification(
          HTML("DCFcluster not found. Install with:<br><code>pip install DCFcluster</code><br>or<br><code>pip install git+https://github.com/tobinjo96/DCFcluster.git</code>"),
          type = "warning", duration = 15)
      }

      # Check numpy and scipy
      numpy_ok <- reticulate::py_module_available("numpy")
      scipy_ok <- reticulate::py_module_available("scipy")
      sklearn_ok <- reticulate::py_module_available("sklearn")

      if(!numpy_ok || !scipy_ok || !sklearn_ok) {
        missing <- c()
        if(!numpy_ok) missing <- c(missing, "numpy")
        if(!scipy_ok) missing <- c(missing, "scipy")
        if(!sklearn_ok) missing <- c(missing, "sklearn")
        showNotification(paste("Missing Python packages:", paste(missing, collapse = ", ")),
                         type = "warning", duration = 10)
      }

    }, error = function(e) {
      showNotification(paste("Error checking DCF setup:", e$message),
                       type = "error", duration = 10)
    })
  })

  # DCF clustering function
  run_dcf_clustering <- function(data_matrix, k_param, beta_param) {
    # Run DCF clustering via Python
    #
    # Args:
    #   data_matrix: numeric matrix (subjects x time points)
    #   k_param: neighborhood parameter for density estimation
    #   beta_param: fluctuation parameter for cluster cores
    #
    # Returns:
    #   List with cluster assignments and metadata

    tryCatch({
      # Import DCFcluster
      dcf <- reticulate::import("DCFcluster")

      # Convert R matrix to numpy array
      np <- reticulate::import("numpy")
      X <- np$array(data_matrix)

      # Run DCF clustering
      result <- dcf$DCFcluster$train(X, as.integer(k_param), beta_param)

      # Extract results
      labels <- as.integer(result$labels) + 1L  # Convert 0-indexed to 1-indexed
      peak_values <- as.numeric(result$peak_values)
      core_sets <- result$core_sets

      # Number of clusters found
      n_clusters <- length(unique(labels))

      return(list(
        labels = labels,
        peak_values = peak_values,
        core_sets = core_sets,
        n_clusters = n_clusters,
        k = k_param,
        beta = beta_param,
        success = TRUE
      ))

    }, error = function(e) {
      return(list(
        success = FALSE,
        error = e$message
      ))
    })
  }

  # Perform clustering (K-means or DCF)
  observeEvent(input$run_clustering, {
    req(values$data)

    tryCatch({
      # Determine which data to use
      if(input$cluster_data_type == "smoothed") {
        if(is.null(values$smooth_data)) {
          showNotification("No smoothed data available. Please apply smoothing first or use raw data.",
                           type = "error", duration = 5)
          return(NULL)
        }
        data_matrix <- values$smooth_data
        data_type_label <- "Smoothed Data"
      } else {
        data_matrix <- values$data
        data_type_label <- "Raw Data"
      }

      # Apply standardization if requested
      standardize_opt <- if(!is.null(input$cluster_standardize)) input$cluster_standardize else "none"
      if(standardize_opt == "within") {
        # Within-participant: z-score each participant's own time series (row-wise)
        data_matrix <- t(apply(data_matrix, 1, function(x) {
          s <- sd(x, na.rm = TRUE)
          if(is.na(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
        }))
        data_type_label <- paste0(data_type_label, " [within-participant standardized]")
      } else if(standardize_opt == "between") {
        # Between-participant: z-score each time point across participants (column-wise)
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        col_sds   <- apply(data_matrix, 2, sd, na.rm = TRUE)
        col_sds[col_sds == 0] <- 1  # avoid division by zero for constant time points
        data_matrix <- sweep(sweep(data_matrix, 2, col_means, "-"), 2, col_sds, "/")
        data_type_label <- paste0(data_type_label, " [between-participant standardized]")
      }

      # Get clustering method
      clustering_method <- if(!is.null(input$clustering_method)) input$clustering_method else "standard"

      # Check for functional k-means requirements
      if(clustering_method == "functional") {
        if(input$cluster_data_type != "smoothed") {
          showNotification("Functional K-Means requires smoothed data. Please select 'Smoothed Data' above.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(is.null(values$fd_obj)) {
          showNotification("No functional data object available. Please apply smoothing first.",
                           type = "error", duration = 5)
          return(NULL)
        }
      }

      # Check for DCF requirements
      if(clustering_method == "dcf") {
        if(!requireNamespace("reticulate", quietly = TRUE)) {
          showNotification("reticulate package not available. Please install it for DCF clustering.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(!reticulate::py_module_available("DCFcluster")) {
          showNotification("DCFcluster Python package not found. Install with: pip install git+https://github.com/tobinjo96/DCFcluster.git",
                           type = "error", duration = 10)
          return(NULL)
        }
      }

      # Check for missing values
      if(any(is.na(data_matrix))) {
        showNotification("Data contains missing values. Please handle missing values before clustering.",
                         type = "error", duration = 5)
        return(NULL)
      }

      # Get parameters
      k <- input$n_clusters
      nstart <- input$kmeans_nstart
      iter_max <- input$kmeans_iter

      # ===== PERFORM CLUSTERING BASED ON METHOD =====

      if(clustering_method == "dcf") {
        # ===== DCF (DENSITY CORE FINDING) CLUSTERING =====
        showNotification("Running DCF clustering...", type = "message", duration = 2)

        dcf_k <- if(!is.null(input$dcf_k)) input$dcf_k else 10
        dcf_beta <- if(!is.null(input$dcf_beta)) input$dcf_beta else 0.1

        dcf_result <- run_dcf_clustering(data_matrix, dcf_k, dcf_beta)

        if(!dcf_result$success) {
          showNotification(paste("DCF clustering failed:", dcf_result$error),
                           type = "error", duration = 10)
          return(NULL)
        }

        cluster_assignments <- dcf_result$labels
        k <- dcf_result$n_clusters  # DCF determines number of clusters automatically

        # Calculate cluster means
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))
        for(i in 1:k) {
          cluster_idx <- which(cluster_assignments == i)
          if(length(cluster_idx) > 0) {
            cluster_means[i, ] <- colMeans(data_matrix[cluster_idx, , drop = FALSE])
          }
        }

        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Calculate within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared <- between_ss / total_ss

        # Create compatible result object
        kmeans_result <- list(
          cluster = cluster_assignments,
          centers = cluster_means,
          size = cluster_sizes,
          tot.withinss = within_ss,
          betweenss = between_ss
        )

        method_label <- paste0("DCF (Density Core Finding, k=", dcf_k, ", β=", dcf_beta, ")")

        # Store DCF-specific results
        dcf_extra <- list(
          peak_values = dcf_result$peak_values,
          core_sets = dcf_result$core_sets,
          k_param = dcf_k,
          beta_param = dcf_beta
        )
        hc_result <- NULL

      } else if(clustering_method == "functional") {
        # ===== FUNCTIONAL K-MEANS using fda.usc =====
        if(!requireNamespace("fda.usc", quietly = TRUE)) {
          showNotification("fda.usc package not available. Please install it or use Standard K-Means.",
                           type = "error", duration = 5)
          return(NULL)
        }

        # Convert data matrix to fdata object for fda.usc
        # fdata expects: rows = observations (subjects), cols = time points
        time_grid <- seq(0, 1, length.out = ncol(data_matrix))
        fdata_obj <- fda.usc::fdata(mdata = data_matrix, argvals = time_grid)

        # Run fda.usc::kmeans.fd nstart times; keep the solution with the lowest WCSS
        best_fkm  <- NULL
        best_wcss <- Inf
        for(restart in seq_len(nstart)) {
          set.seed(100 + restart)
          fkm_try <- tryCatch(
            fda.usc::kmeans.fd(fdata_obj, ncl = k, max.iter = iter_max),
            error = function(e) NULL
          )
          if(is.null(fkm_try)) next
          asn <- fkm_try$cluster
          wcss_try <- sum(sapply(seq_len(k), function(i) {
            idx <- which(asn == i)
            if(length(idx) == 0) return(0)
            m <- data_matrix[idx, , drop = FALSE]
            sum(apply(m, 1, function(x) sum((x - colMeans(m))^2)))
          }))
          if(wcss_try < best_wcss) { best_wcss <- wcss_try; best_fkm <- fkm_try }
        }
        if(is.null(best_fkm)) {
          showNotification("Functional K-Means failed across all restarts.", type = "error", duration = 5)
          return(NULL)
        }
        fkm_result <- best_fkm

        # Extract cluster assignments
        cluster_assignments <- fkm_result$cluster

        # Calculate cluster means by evaluating fd object for each cluster
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))

        for(i in 1:k) {
          cluster_idx <- which(cluster_assignments == i)
          if(length(cluster_idx) > 0) {
            cluster_fd <- values$fd_obj[cluster_idx]
            cluster_mean_fd <- mean.fd(cluster_fd)
            cluster_means[i, ] <- eval.fd(time_grid, cluster_mean_fd)
          }
        }

        # Calculate cluster sizes
        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Calculate within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared <- between_ss / total_ss

        # Create a compatible "kmeans_result" object for consistency
        kmeans_result <- list(
          cluster = cluster_assignments,
          centers = cluster_means,
          size = cluster_sizes,
          tot.withinss = within_ss,
          betweenss = between_ss
        )

        method_label <- "Functional K-Means (fda.usc)"
        dcf_extra <- NULL
        hc_result <- NULL

      } else if(clustering_method == "hierarchical") {
        # ===== HIERARCHICAL CLUSTERING =====
        linkage <- if(!is.null(input$hclust_linkage)) input$hclust_linkage else "ward.D2"

        dist_matrix_hc <- dist(data_matrix)
        hc_result <- hclust(dist_matrix_hc, method = linkage)
        cluster_assignments <- cutree(hc_result, k = k)

        # Cluster means
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))
        for(i in 1:k) {
          idx <- which(cluster_assignments == i)
          if(length(idx) > 0)
            cluster_means[i, ] <- colMeans(data_matrix[idx, , drop = FALSE])
        }
        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(members) > 0)
            wcss_per_cluster[i] <- sum(apply(members, 1, function(x) sum((x - cluster_means[i, ])^2)))
        }

        total_ss   <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss  <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared  <- between_ss / total_ss

        kmeans_result <- list(
          cluster     = cluster_assignments,
          centers     = cluster_means,
          size        = cluster_sizes,
          tot.withinss = within_ss,
          betweenss   = between_ss
        )

        method_label <- paste0("Hierarchical Clustering (", linkage, ")")
        dcf_extra <- NULL

      } else {
        # ===== STANDARD K-MEANS =====
        hc_result <- NULL
        set.seed(123)  # For reproducibility
        kmeans_result <- kmeans(data_matrix, centers = k,
                                nstart = nstart, iter.max = iter_max)

        cluster_assignments <- kmeans_result$cluster
        cluster_means <- kmeans_result$centers
        cluster_sizes <- kmeans_result$size

        # Calculate within-cluster sum of squares for each cluster
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        between_ss <- kmeans_result$betweenss
        within_ss <- kmeans_result$tot.withinss
        r_squared <- between_ss / total_ss
        dcf_extra <- NULL

        method_label <- "Standard (Point-wise) K-Means"
        dcf_extra <- NULL
      }

      # ===== COMMON CALCULATIONS FOR ALL METHODS =====

      # Silhouette width (requires cluster package)
      sil_width <- NA
      sil_data <- NULL
      if(requireNamespace("cluster", quietly = TRUE) && k > 1) {
        dist_matrix <- dist(data_matrix)
        sil <- cluster::silhouette(cluster_assignments, dist_matrix)
        sil_width <- mean(sil[, 3])
        # Store full silhouette object for detailed plotting
        sil_data <- sil
      }

      # Calinski-Harabasz index (variance ratio criterion)
      n <- nrow(data_matrix)
      ch_index <- if(k > 1) (between_ss / (k - 1)) / (within_ss / (n - k)) else NA

      # Store results
      values$clustering_results <- list(
        kmeans = kmeans_result,
        cluster_means = cluster_means,
        cluster_assignments = kmeans_result$cluster,
        cluster_sizes = kmeans_result$size,
        wcss_per_cluster = wcss_per_cluster,
        total_wcss = within_ss,
        between_ss = between_ss,
        total_ss = total_ss,
        r_squared = r_squared,
        silhouette_width = sil_width,
        silhouette_data = sil_data,  # Full silhouette object for detailed plot
        ch_index = ch_index,
        k = k,
        data_type = data_type_label,
        data_matrix = data_matrix,
        method = clustering_method,
        method_label = method_label,
        dcf_extra = dcf_extra,  # DCF-specific results (NULL for other methods)
        hc_obj = hc_result,     # hclust object (NULL for non-hierarchical methods)
        linkage_method = if(clustering_method == "hierarchical") {
          if(!is.null(input$hclust_linkage)) input$hclust_linkage else "ward.D2"
        } else NULL
      )

      # Notification message
      notification_msg <- if(clustering_method == "dcf") {
        paste0("DCF clustering completed: ", k, " clusters found automatically")
      } else {
        paste0("Clustering completed successfully with k=", k)
      }
      showNotification(notification_msg, type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Clustering error:", e$message),
                       type = "error", duration = 5)
      values$clustering_results <- NULL
    })
  })

  # Clustering status output
  output$clustering_status <- renderText({
    if(is.null(values$clustering_results)) {
      "No clustering performed yet. Click 'Run Clustering' to begin."
    } else {
      method_info <- if(!is.null(values$clustering_results$method_label)) {
        values$clustering_results$method_label
      } else {
        "Clustering"
      }
      paste0("✓ Clustering completed with k=", values$clustering_results$k,
             " clusters\n   Method: ", method_info,
             "\n   Data: ", values$clustering_results$data_type)
    }
  })

  # Save cluster membership to dataset
  observeEvent(input$save_cluster_membership, {
    req(values$clustering_results)

    tryCatch({
      # Validate variable name
      var_name <- trimws(input$cluster_var_name)

      if(var_name == "" || is.null(var_name)) {
        showNotification("Please enter a variable name.", type = "error", duration = 3)
        return(NULL)
      }

      # Check for invalid characters in variable name
      if(!grepl("^[a-zA-Z][a-zA-Z0-9._]*$", var_name)) {
        showNotification("Variable name must start with a letter and contain only letters, numbers, dots, or underscores.",
                         type = "error", duration = 5)
        return(NULL)
      }

      # Get cluster assignments
      cluster_assignments <- values$clustering_results$cluster_assignments

      # Check if we have uploaded data
      if(is.null(values$uploaded_data)) {
        showNotification("No uploaded data available. Cluster membership can only be saved if you uploaded data with ID variables.",
                         type = "warning", duration = 5)
        return(NULL)
      }

      # Check dimensions match
      if(nrow(values$uploaded_data) != length(cluster_assignments)) {
        showNotification(paste0("Data dimension mismatch. Uploaded data has ", nrow(values$uploaded_data),
                               " rows but clustering has ", length(cluster_assignments), " observations."),
                         type = "error", duration = 5)
        return(NULL)
      }

      # Check if variable already exists
      if(var_name %in% colnames(values$uploaded_data)) {
        showModal(modalDialog(
          title = "Variable Already Exists",
          paste0("The variable '", var_name, "' already exists in your dataset. Do you want to overwrite it?"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_overwrite", "Overwrite", class = "btn-warning")
          )
        ))
        return(NULL)
      }

      # Add cluster membership to uploaded data
      values$uploaded_data[[var_name]] <- as.factor(cluster_assignments)

      showNotification(paste0("✓ Variable '", var_name, "' added to in-memory dataset (",
                             values$clustering_results$k, " clusters). ",
                             "Use 'Download Dataset' button to export to file."),
                       type = "message", duration = 7)

    }, error = function(e) {
      showNotification(paste("Error saving cluster membership:", e$message),
                       type = "error", duration = 5)
    })
  })

  # Handle overwrite confirmation
  observeEvent(input$confirm_overwrite, {
    var_name <- trimws(input$cluster_var_name)
    cluster_assignments <- values$clustering_results$cluster_assignments

    # Overwrite the variable
    values$uploaded_data[[var_name]] <- as.factor(cluster_assignments)

    showNotification(paste0("✓ Variable '", var_name, "' overwritten in in-memory dataset. ",
                           "Use 'Download Dataset' button to export to file."),
                     type = "message", duration = 7)

    removeModal()
  })

  # Cluster save status
  output$cluster_save_status <- renderText({
    req(values$clustering_results)

    var_name <- trimws(input$cluster_var_name)

    if(is.null(values$uploaded_data)) {
      return("No uploaded data")
    }

    if(var_name == "" || is.null(var_name)) {
      return("Enter variable name")
    }

    if(var_name %in% colnames(values$uploaded_data)) {
      return(paste0("✓ '", var_name, "' exists\n(will prompt to overwrite)"))
    } else {
      return("Ready to add")
    }
  })

  # Download dataset with cluster membership
  output$download_data_with_clusters <- downloadHandler(
    filename = function() {
      paste0("data_with_clusters_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$uploaded_data)

      # Create a copy of uploaded data
      data_to_export <- values$uploaded_data

      # If cluster membership was added, it's already in values$uploaded_data
      # If not yet added, add it now
      var_name <- trimws(input$cluster_var_name)
      if(!is.null(values$clustering_results)) {
        if(!(var_name %in% colnames(data_to_export)) && var_name != "") {
          # Add cluster membership if it wasn't added yet
          cluster_assignments <- values$clustering_results$cluster_assignments
          if(nrow(data_to_export) == length(cluster_assignments)) {
            data_to_export[[var_name]] <- as.factor(cluster_assignments)
          }
        }
      }

      # Write to CSV
      write.csv(data_to_export, file, row.names = FALSE)
    }
  )

  # Cluster summary table
  output$cluster_summary_table <- renderDT({
    req(values$clustering_results)

    wcss <- values$clustering_results$wcss_per_cluster
    sizes <- values$clustering_results$cluster_sizes

    summary_df <- data.frame(
      Cluster = 1:values$clustering_results$k,
      Size = sizes,
      Percentage = round(100 * sizes / sum(sizes), 1),
      Within_SS = round(wcss, 2),
      Avg_Within_SS = round(wcss / sizes, 2)
    )

    colnames(summary_df) <- c("Cluster", "Members", "% of Total",
                               "Within-Cluster SS", "Avg. Within-SS")

    datatable(summary_df,
              options = list(
                dom = 't',
                pageLength = 20,
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(summary_df),
                  fontSize = '14px')
  })

  # Clustering fit statistics
  output$clustering_fit_stats <- renderText({
    req(values$clustering_results)

    r <- values$clustering_results

    # Method info
    method_info <- if(!is.null(r$method_label)) r$method_label else "Clustering"
    output_text <- paste0("Method: ", method_info, "\n",
                          "Number of Clusters: ", r$k, "\n\n")

    # Standard statistics
    output_text <- paste0(output_text,
      "Total Sum of Squares: ", round(r$total_ss, 2), "\n",
      "Between-Cluster Sum of Squares: ", round(r$between_ss, 2), "\n",
      "Within-Cluster Sum of Squares: ", round(r$total_wcss, 2), "\n",
      "R² (Variance Explained): ", round(r$r_squared * 100, 2), "%\n",
      "Between-SS / Total-SS Ratio: ", round(r$between_ss / r$total_ss * 100, 2), "%\n"
    )

    if(!is.na(r$silhouette_width)) {
      output_text <- paste0(output_text,
                            "Average Silhouette Width: ", round(r$silhouette_width, 3), "\n",
                            "  (Range: -1 to 1, higher is better)\n")
    }

    if(!is.na(r$ch_index)) {
      output_text <- paste0(output_text,
                            "Calinski-Harabasz Index: ", round(r$ch_index, 2), "\n",
                            "  (Higher values indicate better-defined clusters)\n")
    }

    # DCF-specific information
    if(!is.null(r$dcf_extra)) {
      output_text <- paste0(output_text, "\n--- DCF Parameters ---\n",
                            "Neighborhood k: ", r$dcf_extra$k_param, "\n",
                            "Fluctuation β: ", r$dcf_extra$beta_param, "\n",
                            "Clusters detected automatically: ", r$k, "\n")
    }

    output_text
  })

  # Cluster membership table
  output$cluster_membership_table <- renderDT({
    req(values$clustering_results)

    # Create membership table
    membership_df <- data.frame(
      Subject_ID = 1:length(values$clustering_results$cluster_assignments),
      Cluster = values$clustering_results$cluster_assignments
    )

    # Add ALL group variables if available (not just the primary one)
    if(!is.null(values$group_variables) && ncol(values$group_variables) > 0) {
      # Add all group variables from group_variables data frame
      for(var_name in colnames(values$group_variables)) {
        membership_df[[var_name]] <- values$group_variables[[var_name]]
      }
    } else if(!is.null(values$group_labels)) {
      # Fallback to primary group_labels
      membership_df$Group <- values$group_labels
    }

    # Generate stronger colors with transparency (same as plots)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Add transparency to colors for table background
    table_colors <- paste0(strong_colors[1:values$clustering_results$k], "40")

    datatable(membership_df,
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                scrollY = "400px",
                searching = TRUE,
                ordering = TRUE
              ),
              rownames = FALSE,
              filter = 'top') %>%
      formatStyle('Cluster',
                  backgroundColor = styleEqual(
                    1:values$clustering_results$k,
                    table_colors
                  ))
  })

  # Check if group labels are available
  output$has_group_labels <- reactive({
    !is.null(values$group_labels) && !is.null(values$clustering_results)
  })
  outputOptions(output, "has_group_labels", suspendWhenHidden = FALSE)

  # UI for selecting which group variable to use in cluster composition analysis
  output$cluster_group_var_selector <- renderUI({
    req(values$group_labels)

    # Get available group variables
    available_vars <- NULL
    if(!is.null(values$selected_group_vars) && length(values$selected_group_vars) >= 1) {
      available_vars <- values$selected_group_vars
    } else if(!is.null(values$group_labels)) {
      # Fallback: use "Group" as default name if only group_labels exists
      available_vars <- "Group"
    }

    if(!is.null(available_vars)) {
      # Get current selection or default to first
      current_selection <- if(!is.null(input$cluster_group_var) &&
                              input$cluster_group_var %in% available_vars) {
        input$cluster_group_var
      } else {
        available_vars[1]
      }

      tagList(
        selectInput("cluster_group_var",
                    "Select Group Variable for Composition Analysis:",
                    choices = available_vars,
                    selected = current_selection),
        helpText(paste("Available group variables:", paste(available_vars, collapse = ", "))),
        hr()
      )
    } else {
      NULL
    }
  })

  # Helper to get the selected group variable for clustering
  get_cluster_group_labels <- reactive({
    # First check if a specific variable is selected and available in group_variables
    if(!is.null(input$cluster_group_var) && !is.null(values$group_variables) &&
       input$cluster_group_var %in% colnames(values$group_variables)) {
      labels <- values$group_variables[[input$cluster_group_var]]
      # Ensure it's a factor
      if(!is.factor(labels)) {
        labels <- factor(labels)
      }
      return(labels)
    }
    # Fallback to primary group_labels
    values$group_labels
  })

  # Get the name of currently selected group variable
  get_cluster_group_var_name <- reactive({
    if(!is.null(input$cluster_group_var)) {
      input$cluster_group_var
    } else if(!is.null(values$selected_group_vars) && length(values$selected_group_vars) > 0) {
      values$selected_group_vars[1]
    } else {
      "Group"
    }
  })

  # Group counts table by cluster
  output$cluster_group_counts_table <- renderDT({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Create counts table
    counts_df <- as.data.frame.matrix(cont_table)
    counts_df <- cbind(Cluster = rownames(counts_df), counts_df)
    counts_df$Total <- rowSums(cont_table)

    datatable(counts_df,
              options = list(
                dom = 't',
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(counts_df),
                  fontSize = '14px')
  })

  # Group percentages table by cluster
  output$cluster_group_pct_table <- renderDT({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Calculate percentages for each cluster (row percentages)
    cluster_pcts <- prop.table(cont_table, margin = 1) * 100

    # Create percentages table
    pct_df <- as.data.frame.matrix(round(cluster_pcts, 1))
    pct_df <- cbind(Cluster = rownames(pct_df), pct_df)

    datatable(pct_df,
              options = list(
                dom = 't',
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(pct_df),
                  fontSize = '14px')
  })

  # Chi-square test for group-cluster independence
  output$cluster_group_test <- renderText({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()
    group_var_name <- get_cluster_group_var_name()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Get table dimensions for reporting
    n_clusters <- nrow(cont_table)
    n_groups <- ncol(cont_table)
    table_cells <- length(cont_table)

    # Perform chi-square test (works for any table size)
    chi_test <- tryCatch({
      suppressWarnings(chisq.test(cont_table))
    }, error = function(e) {
      NULL
    })

    # Check if chi-square assumptions are violated (expected frequencies < 5)
    chi_warning <- FALSE
    if(!is.null(chi_test)) {
      expected <- chi_test$expected
      if(any(expected < 5)) {
        chi_warning <- TRUE
      }
    }

    # Try Fisher's exact test with simulation for larger tables
    # Fisher test with simulation can handle larger tables
    fisher_test <- tryCatch({
      # Use simulation-based p-value for tables of any size
      fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    }, error = function(e) {
      # If still fails (very rare), return NULL
      NULL
    })

    # Build output text
    output_text <- "=== CLUSTER-GROUP ASSOCIATION TEST ===\n\n"

    # Show which group variable is being used
    output_text <- paste0(output_text, "Group Variable: ", group_var_name, "\n")
    output_text <- paste0(output_text, "Table dimensions: ", n_clusters, " clusters × ", n_groups, " groups\n\n")

    # Overall percentages (limit display for many groups)
    output_text <- paste0(output_text, "Overall ", group_var_name, " Distribution:\n")
    overall_pcts <- prop.table(table(groups)) * 100
    group_names <- names(overall_pcts)
    if(length(group_names) > 10) {
      # Show only first 10 groups if there are many
      for(g in group_names[1:10]) {
        output_text <- paste0(output_text, "  ", g, ": ", round(overall_pcts[g], 1), "%\n")
      }
      output_text <- paste0(output_text, "  ... and ", length(group_names) - 10, " more groups\n")
    } else {
      for(g in group_names) {
        output_text <- paste0(output_text, "  ", g, ": ", round(overall_pcts[g], 1), "%\n")
      }
    }
    output_text <- paste0(output_text, "\n")

    # Chi-square test results
    if(!is.null(chi_test)) {
      output_text <- paste0(output_text, "Chi-Square Test of Independence:\n")
      output_text <- paste0(output_text, "  χ² = ", round(chi_test$statistic, 3), "\n")
      output_text <- paste0(output_text, "  df = ", chi_test$parameter, "\n")
      output_text <- paste0(output_text, "  p-value = ", format.pval(chi_test$p.value, digits = 4), "\n")

      if(chi_warning) {
        output_text <- paste0(output_text, "  Note: Some expected frequencies < 5; interpret with caution.\n")
      }
      output_text <- paste0(output_text, "\n")

      if(chi_test$p.value < 0.001) {
        output_text <- paste0(output_text, "  *** HIGHLY SIGNIFICANT (p < 0.001) ***\n")
        output_text <- paste0(output_text, "  Cluster membership is strongly associated with group membership.\n\n")
      } else if(chi_test$p.value < 0.01) {
        output_text <- paste0(output_text, "  ** SIGNIFICANT (p < 0.01) **\n")
        output_text <- paste0(output_text, "  Cluster membership is significantly associated with group membership.\n\n")
      } else if(chi_test$p.value < 0.05) {
        output_text <- paste0(output_text, "  * SIGNIFICANT (p < 0.05) *\n")
        output_text <- paste0(output_text, "  Cluster membership is associated with group membership.\n\n")
      } else {
        output_text <- paste0(output_text, "  Not significant (p >= 0.05)\n")
        output_text <- paste0(output_text, "  No significant association between cluster and group membership.\n\n")
      }
    }

    # Fisher's exact test results (if available)
    if(!is.null(fisher_test)) {
      output_text <- paste0(output_text, "Fisher's Exact Test (Monte Carlo simulation, B=10000):\n")
      output_text <- paste0(output_text, "  p-value = ", format.pval(fisher_test$p.value, digits = 4), "\n")
      if(fisher_test$p.value < 0.05) {
        output_text <- paste0(output_text, "  Significant association (p < 0.05)\n\n")
      } else {
        output_text <- paste0(output_text, "  Not significant (p >= 0.05)\n\n")
      }
    } else if(is.null(chi_test)) {
      output_text <- paste0(output_text, "Statistical tests could not be performed.\n")
      output_text <- paste0(output_text, "This may occur with very sparse tables or unusual data distributions.\n\n")
    }

    # Calculate and report effect size (Cramér's V)
    if(!is.null(chi_test)) {
      n <- sum(cont_table)
      k <- min(nrow(cont_table), ncol(cont_table))
      cramers_v <- sqrt(chi_test$statistic / (n * (k - 1)))
      output_text <- paste0(output_text, "Effect Size (Cramér's V): ", round(cramers_v, 3), "\n")
      if(cramers_v < 0.1) {
        output_text <- paste0(output_text, "  (Negligible effect)\n")
      } else if(cramers_v < 0.3) {
        output_text <- paste0(output_text, "  (Small effect)\n")
      } else if(cramers_v < 0.5) {
        output_text <- paste0(output_text, "  (Medium effect)\n")
      } else {
        output_text <- paste0(output_text, "  (Large effect)\n")
      }
    }

    output_text
  })

  # Cluster means plot
  output$cluster_means_plot <- renderPlotly({
    req(values$clustering_results)

    cluster_means <- values$clustering_results$cluster_means
    k <- values$clustering_results$k
    data_matrix <- values$clustering_results$data_matrix
    clusters <- values$clustering_results$cluster_assignments

    # Get time points using helper functions (same as preprocessing module)
    hour_labels <- get_hour_labels()

    if(!is.null(hour_labels)) {
      # Use actual hour values from column names
      time_positions <- calculate_time_positions(hour_labels)
      if(!is.null(time_positions)) {
        time_points <- time_positions
      } else {
        time_points <- seq(0, 1, length.out = ncol(cluster_means))
      }
      tick_text <- sapply(hour_labels, decimal_to_hhmm)

      # Snap ticks to round-hour multiples
      kmeans_step <- as.numeric(input$tick_freq_kmeans)
      if (!is.na(kmeans_step) && kmeans_step > 0) {
        keep <- abs(round(hour_labels / kmeans_step) * kmeans_step - hour_labels) < 1e-6
        if (any(keep)) {
          tick_vals_subset <- time_points[keep]
          tick_text_subset <- tick_text[keep]
        } else {
          tick_vals_subset <- time_points
          tick_text_subset <- tick_text
        }
      } else {
        tick_vals_subset <- time_points
        tick_text_subset <- tick_text
      }
    } else if(!is.null(values$time_labels)) {
      # Use column names as labels
      time_points <- seq(0, 1, length.out = length(values$time_labels))
      tick_vals_subset <- time_points
      tick_text_subset <- values$time_labels
    } else {
      # Fallback to numeric indices
      time_points <- 1:ncol(cluster_means)
      tick_vals_subset <- time_points
      tick_text_subset <- as.character(time_points)
    }
    hover_times <- hover_time_labels(time_points)

    x_label <- get_time_label()

    # Compute a nice y-axis tick interval targeting ~5 ticks
    y_span <- diff(range(data_matrix, na.rm = TRUE))
    y_dtick <- if (y_span > 0) {
      raw_step <- y_span / 5
      pow10 <- 10^floor(log10(raw_step))
      pow10 * ceiling(raw_step / pow10)
    } else NULL

    # Generate stronger colors (not pastel)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Create plotly object
    p <- plot_ly()

    # Add traces for each cluster
    for(i in 1:k) {
      # Get cluster members
      cluster_idx <- which(clusters == i)
      cluster_data <- data_matrix[cluster_idx, , drop = FALSE]

      # Calculate mean, SD for each time point
      cluster_mean <- cluster_means[i, ]
      cluster_sd <- apply(cluster_data, 2, sd)

      # Add lower SD band FIRST (this will be the bottom of the ribbon)
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean - cluster_sd,
          type = 'scatter',
          mode = 'lines',
          line = list(width = 0),
          name = paste0("Cluster ", i, " ±SD"),
          legendgroup = paste0("cluster_sd_", i),
          showlegend = FALSE,
          hoverinfo = 'skip'
        )

      # Add upper SD band SECOND with fill = 'tonexty' to create ribbon
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean + cluster_sd,
          type = 'scatter',
          mode = 'lines',
          line = list(width = 0),
          fillcolor = paste0(strong_colors[i], '40'),  # 40 = 25% opacity in hex
          fill = 'tonexty',
          name = paste0("Cluster ", i, " ±SD"),
          legendgroup = paste0("cluster_sd_", i),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )

      # Add mean curve (smooth interpolation using spline)
      # Use SEPARATE legendgroup so it can be toggled independently
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean,
          type = 'scatter',
          mode = 'lines',
          line = list(color = strong_colors[i], width = 2, shape = 'spline'),
          name = paste0("Cluster ", i, " Mean"),
          legendgroup = paste0("cluster_mean_", i),
          showlegend = TRUE,
          hovertemplate = paste0(
            "Cluster: ", i,
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<br>Size: ", values$clustering_results$cluster_sizes[i], " subjects",
            "<extra></extra>"
          ),
          customdata = hover_times
        )
    }

    # Update layout with proper time axis formatting
    p <- p %>%
      layout(
        title = list(text = "Cluster Mean Functions with ±1 SD Bands",
                     x = 0.5, xanchor = "center"),
        xaxis = list(
          title = x_label,
          tickmode = 'array',
          tickvals = tick_vals_subset,
          ticktext = tick_text_subset,
          tickangle = -90  # Rotate labels for readability
        ),
        yaxis = list(title = "Value", dtick = y_dtick),
        hovermode = "closest",
        legend = list(
          orientation = "v",
          x = 1.02,
          y = 1,
          tracegroupgap = 0
        )
      )

    p
  })

  # Individual curves by cluster
  output$cluster_individuals_plot <- renderPlotly({
    req(values$clustering_results)

    data_matrix <- values$clustering_results$data_matrix
    clusters <- values$clustering_results$cluster_assignments
    k <- values$clustering_results$k

    # Get time points using helper functions (same as preprocessing module)
    hour_labels <- get_hour_labels()

    if(!is.null(hour_labels)) {
      # Use actual hour values from column names
      time_positions <- calculate_time_positions(hour_labels)
      if(!is.null(time_positions)) {
        time_points <- time_positions
      } else {
        time_points <- seq(0, 1, length.out = ncol(data_matrix))
      }
      tick_text <- sapply(hour_labels, decimal_to_hhmm)

      # Snap ticks to round-hour multiples
      kmeans_step <- as.numeric(input$tick_freq_kmeans)
      if (!is.na(kmeans_step) && kmeans_step > 0) {
        keep <- abs(round(hour_labels / kmeans_step) * kmeans_step - hour_labels) < 1e-6
        if (any(keep)) {
          tick_vals_subset <- time_points[keep]
          tick_text_subset <- tick_text[keep]
        } else {
          tick_vals_subset <- time_points
          tick_text_subset <- tick_text
        }
      } else {
        tick_vals_subset <- time_points
        tick_text_subset <- tick_text
      }
    } else if(!is.null(values$time_labels)) {
      # Use column names as labels
      time_points <- seq(0, 1, length.out = length(values$time_labels))
      tick_vals_subset <- time_points
      tick_text_subset <- values$time_labels
    } else {
      # Fallback to numeric indices
      time_points <- 1:ncol(data_matrix)
      tick_vals_subset <- time_points
      tick_text_subset <- as.character(time_points)
    }
    hover_times <- hover_time_labels(time_points)

    x_label <- get_time_label()

    # Compute a nice y-axis tick interval targeting ~5 ticks
    y_span <- diff(range(data_matrix, na.rm = TRUE))
    y_dtick <- if (y_span > 0) {
      raw_step <- y_span / 5
      pow10 <- 10^floor(log10(raw_step))
      pow10 * ceiling(raw_step / pow10)
    } else NULL

    # Generate stronger colors (same as cluster means plot)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Calculate overall mean across all subjects (as reference)
    overall_mean <- colMeans(data_matrix)

    # Create subplots for each cluster
    plot_list <- list()

    for(i in 1:k) {
      cluster_idx <- which(clusters == i)
      cluster_data <- data_matrix[cluster_idx, , drop = FALSE]

      # Create plotly object for this cluster
      p <- plot_ly()

      # Add individual curves (each as a separate trace for interactivity)
      # Group them together so they can be toggled together
      for(j in 1:nrow(cluster_data)) {
        p <- p %>%
          add_trace(
            x = time_points,
            y = cluster_data[j, ],
            type = 'scatter',
            mode = 'lines',
            line = list(color = strong_colors[i], width = 0.5),
            opacity = 0.15,
            name = paste0("Cluster ", i, " Individuals"),
            legendgroup = paste0("individual_", i),
            showlegend = (j == 1),  # Only show in legend once
            hovertemplate = paste0(
              "Subject: ", cluster_idx[j],
              "<br>Time: %{customdata}",
              "<br>Value: %{y:.2f}",
              "<extra></extra>"
            ),
            customdata = hover_times
          )
      }

      # Add cluster mean - bold and opaque
      p <- p %>%
        add_trace(
          x = time_points,
          y = values$clustering_results$cluster_means[i, ],
          type = 'scatter',
          mode = 'lines',
          line = list(color = strong_colors[i], width = 2.5, shape = 'spline'),
          opacity = 1,
          name = paste0("Cluster ", i, " Mean"),
          legendgroup = paste0("mean_", i),
          showlegend = TRUE,
          hovertemplate = paste0(
            "Cluster ", i, " Mean",
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<extra></extra>"
          ),
          customdata = hover_times
        )

      # Add overall mean as reference - dashed black line
      p <- p %>%
        add_trace(
          x = time_points,
          y = overall_mean,
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'grey', width = 2, dash = 'dot', shape = 'spline'),
          opacity = 0.8,
          name = "Overall Mean",
          legendgroup = "overall_mean",
          showlegend = (i == 1),  # Only show in legend once (first subplot)
          hovertemplate = paste0(
            "Overall Mean (all subjects)",
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<extra></extra>"
          ),
          customdata = hover_times
        )

      # Update layout for this subplot with proper time axis formatting
      p <- p %>%
        layout(
          title = list(text = paste0("Cluster ", i, " (n=", length(cluster_idx), ")"),
                       font = list(size = 12)),
          xaxis = list(
            title = x_label,
            tickmode = 'array',
            tickvals = tick_vals_subset,
            ticktext = tick_text_subset,
            tickangle = -90  # Rotate labels for readability
          ),
          yaxis = list(title = "Value"),
          hovermode = "closest"
        )

      plot_list[[i]] <- p
    }

    # Combine plots
    sp <- subplot(plot_list, nrows = ceiling(k/2), shareX = TRUE, shareY = TRUE,
                  titleX = TRUE, titleY = TRUE)

    # Build layout args, applying dtick to every y-axis (yaxis, yaxis2, yaxis3, ...)
    layout_args <- list(
      sp,
      title = list(text = "Individual Curves by Cluster (click legend to toggle)",
                   x = 0.5, xanchor = "center"),
      showlegend = TRUE,
      legend = list(orientation = "v", x = 1.02, y = 1)
    )
    yaxis_names <- c("yaxis", paste0("yaxis", seq_len(k)[-1]))
    for (nm in yaxis_names) {
      layout_args[[nm]] <- list(dtick = y_dtick)
    }
    do.call(layout, layout_args)
  })

  # Detailed silhouette plot (similar to Python scikit-learn visualization)
  output$detailed_silhouette_plot <- renderPlotly({
    req(values$clustering_results)
    req(values$clustering_results$silhouette_data)

    sil_data <- values$clustering_results$silhouette_data
    k <- values$clustering_results$k

    # Extract silhouette coefficients and cluster assignments
    sil_df <- data.frame(
      cluster = sil_data[, 1],
      neighbor = sil_data[, 2],
      sil_width = sil_data[, 3]
    )

    # Sort within each cluster by silhouette width (descending)
    sil_df <- sil_df %>%
      group_by(cluster) %>%
      arrange(cluster, desc(sil_width)) %>%
      mutate(sample_idx = row_number()) %>%
      ungroup()

    # Calculate cumulative positions for y-axis
    cluster_sizes <- table(sil_df$cluster)
    cumsum_sizes <- c(0, cumsum(cluster_sizes))

    # Add y-position for plotting (stacked bars)
    sil_df <- sil_df %>%
      group_by(cluster) %>%
      mutate(y_position = cumsum_sizes[cluster[1]] + row_number()) %>%
      ungroup()

    # Generate stronger colors (same as other plots)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Create plotly figure
    p <- plot_ly()

    # Add bars for each cluster
    for(i in 1:k) {
      cluster_data <- sil_df %>% filter(cluster == i)

      p <- p %>%
        add_trace(
          x = cluster_data$sil_width,
          y = cluster_data$y_position,
          type = 'bar',
          orientation = 'h',
          marker = list(color = strong_colors[i]),
          name = paste0("Cluster ", i),
          hovertemplate = paste0(
            "Cluster: ", i,
            "<br>Sample: ", 1:nrow(cluster_data),
            "<br>Silhouette: %{x:.3f}",
            "<extra></extra>"
          )
        )
    }

    # Calculate average silhouette width
    avg_sil <- mean(sil_df$sil_width)

    # Add vertical line for average silhouette
    p <- p %>%
      add_trace(
        x = c(avg_sil, avg_sil),
        y = c(0, max(sil_df$y_position) + 1),
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'red', width = 2, dash = 'dash'),
        name = paste0('Average (', round(avg_sil, 3), ')'),
        showlegend = TRUE,
        hoverinfo = 'skip'
      )

    # Add cluster separators
    for(i in 1:(k-1)) {
      sep_y <- cumsum_sizes[i + 1] + 0.5
      p <- p %>%
        add_trace(
          x = c(-1, 1),
          y = c(sep_y, sep_y),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'black', width = 1, dash = 'dot'),
          showlegend = FALSE,
          hoverinfo = 'skip'
        )
    }

    # Update layout
    p <- p %>%
      layout(
        title = list(
          text = paste0("Silhouette Plot for k=", k, " Clusters (Avg: ", round(avg_sil, 3), ")"),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          title = "Silhouette Coefficient",
          range = c(min(c(sil_df$sil_width, -0.1)), 1)
        ),
        yaxis = list(
          title = "Sample Index (sorted within cluster)",
          showticklabels = FALSE
        ),
        barmode = 'overlay',
        showlegend = TRUE,
        legend = list(
          orientation = "v",
          x = 1.02,
          y = 1
        ),
        hovermode = "closest"
      )

    # Add annotations for cluster labels
    for(i in 1:k) {
      cluster_data <- sil_df %>% filter(cluster == i)
      mid_y <- mean(cluster_data$y_position)

      p <- p %>%
        add_annotations(
          x = -0.05,
          y = mid_y,
          text = paste0("C", i),
          showarrow = FALSE,
          xref = "x",
          yref = "y",
          font = list(size = 12, color = strong_colors[i])
        )
    }

    p
  })

  # ===== DENDROGRAM OUTPUTS =====

  # Main clustering tab dendrogram
  output$dendrogram_plot <- renderPlot({
    req(values$clustering_results)
    res <- values$clustering_results
    req(res$method == "hierarchical", !is.null(res$hc_obj))
    hc <- res$hc_obj
    k  <- res$k
    # Determine the height at which to cut for k clusters
    cut_height <- hc$height[length(hc$height) - k + 1]
    par(mar = c(4, 4, 3, 1))
    plot(hc, labels = FALSE, hang = -1,
         main = paste0("Dendrogram — ", res$linkage_method, " linkage (cut at k=", k, ")"),
         xlab = "Participants", ylab = "Height", sub = "")
    abline(h = cut_height, col = "red", lty = 2, lwd = 1.5)
    legend("topright", legend = paste0("Cut height = ", round(cut_height, 3)),
           col = "red", lty = 2, lwd = 1.5, bty = "n")
  })

  # Optimization tab dendrogram
  output$opt_dendrogram_plot <- renderPlot({
    req(values$cluster_optimization)
    opt <- values$cluster_optimization
    req(opt$method == "hierarchical", !is.null(opt$hc_obj))
    hc <- opt$hc_obj
    par(mar = c(4, 4, 3, 1))
    plot(hc, labels = FALSE, hang = -1,
         main = paste0("Full Dendrogram — ", if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2", " linkage"),
         xlab = "Participants", ylab = "Height", sub = "")
  })
