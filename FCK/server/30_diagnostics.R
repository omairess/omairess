# ==========================================================================
# server/30_diagnostics.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 2231-2860  (GAM REML, REML profile, CV, n-basis sweep)
# ==========================================================================
  # ============================================================================
  # SMOOTHING DIAGNOSTICS SERVER CODE
  # ============================================================================
  
  # 1. GAM REML Analysis
  observeEvent(input$run_gam_reml, {
    req(values$data)
    
    showNotification("Fitting GAM with REML...", type = "message", duration = 3)
    
    tryCatch({
      n_time <- ncol(values$data)
      time_points <- 1:n_time
      
      # Select data to analyze
      if(input$diag_subject_type == "mean") {
        y_data <- colMeans(values$data, na.rm = TRUE)
        data_label <- "Mean curve (all subjects)"
      } else {
        subject_id <- min(input$diag_subject_id, nrow(values$data))
        y_data <- values$data[subject_id, ]
        data_label <- paste0("Subject ", subject_id)
      }
      
      # Ensure y_data is numeric and handle NAs
      y_data <- as.numeric(y_data)
      if(all(is.na(y_data))) {
        showNotification("All data points are NA. Cannot perform GAM REML.", type = "error")
        return()
      }
      
      # Remove NAs
      valid_idx <- which(!is.na(y_data))
      if(length(valid_idx) == 0) {
        showNotification("No valid data points for GAM", type = "error")
        return()
      }
      
      y_valid <- y_data[valid_idx]
      t_valid <- time_points[valid_idx]
      
      if(length(y_valid) < 4) {
        showNotification("Not enough valid data points for GAM", type = "error")
        return()
      }
      
      # Fit GAM with REML
      df_gam <- data.frame(y = y_valid, t = t_valid)
      gam_fit <- mgcv::gam(y ~ s(t, bs = "cr"), data = df_gam, method = "REML")
      
      # Store results
      values$gam_reml_fit <- gam_fit
      values$gam_data_label <- data_label
      
      showNotification("GAM REML completed successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("GAM REML error:", e$message), type = "error", duration = 5)
      cat("GAM REML error:", e$message, "\n")
    })
  })
  
  output$gam_reml_summary <- renderPrint({
    if(is.null(values$gam_reml_fit)) {
      cat("No GAM REML analysis performed yet.\n\n")
      cat("Click 'Fit GAM (REML)' to run the analysis.")
      return()
    }
    
    gam_fit <- values$gam_reml_fit
    
    cat("=== GAM REML Analysis ===\n")
    cat("Data:", values$gam_data_label, "\n\n")
    
    # Extract key info
    summary_gam <- summary(gam_fit)
    edf <- summary_gam$edf
    
    # Safely extract k (basis dimension) - column name may vary
    if("k'" %in% colnames(summary_gam$s.table)) {
      k <- summary_gam$s.table[1, "k'"]
    } else if("k" %in% colnames(summary_gam$s.table)) {
      k <- summary_gam$s.table[1, "k"]
    } else {
      k <- NA
    }
    
    reml_score <- gam_fit$gcv.ubre  # REML score
    
    cat("Smoothing parameters:\n")
    cat(sprintf("  Effective Degrees of Freedom (EDF): %.2f\n", edf))
    if(!is.na(k)) {
      cat(sprintf("  Basis dimension (k): %.0f\n", k))
    }
    cat(sprintf("  Estimated lambda: %.3e\n", gam_fit$sp))
    cat(sprintf("  REML score: %.2f\n", reml_score))
    
    cat("\nInterpretation:\n")
    if(edf < 2) {
      cat("  -> Very smooth fit (low complexity)\n")
      cat("  -> May be undersmoothing important features\n")
    } else if(edf < 6) {
      cat("  -> Moderate smoothing (balanced complexity)\n")
      cat("  -> Often optimal for functional data\n")
    } else if(!is.na(k) && edf < k - 2) {
      cat("  -> Flexible fit (high complexity)\n")
      cat("  -> Capturing detailed features\n")
    } else if(!is.na(k) && edf >= k - 2) {
      cat("  -> Minimal smoothing (very high complexity)\n")
      cat("  -> May be overfitting, consider increasing lambda\n")
    } else {
      cat("  -> Flexible fit\n")
      cat("  -> Capturing detailed features\n")
    }
    
    cat("\nModel deviance explained:", sprintf("%.1f%%", summary_gam$dev.expl * 100), "\n")
    
    # Add smoothing term details
    cat("\n--- Smooth term details ---\n")
    print(summary_gam$s.table)
  })
  
  # 2. REML Profile Analysis
  observeEvent(input$run_reml_profile, {
    req(values$data)
    
    showNotification("Computing REML profile...", type = "message", duration = 3)
    
    tryCatch({
      n_time <- ncol(values$data)
      time_points <- 1:n_time
      
      # Select data
      if(input$diag_subject_type == "mean") {
        y_data <- colMeans(values$data, na.rm = TRUE)
      } else {
        subject_id <- min(input$diag_subject_id, nrow(values$data))
        y_data <- values$data[subject_id, ]
      }
      
      valid_idx <- !is.na(y_data)
      y_valid <- y_data[valid_idx]
      t_valid <- time_points[valid_idx]
      
      # Lambda range
      lambda_min <- 10^input$lambda_min_exp
      lambda_max <- 10^input$lambda_max_exp
      lambda_seq <- 10^seq(input$lambda_min_exp, input$lambda_max_exp, 
                           length.out = input$n_lambda)
      
      # Compute REML score for each lambda
      reml_scores <- numeric(length(lambda_seq))
      edf_values <- numeric(length(lambda_seq))
      
      withProgress(message = "Computing REML profile...", value = 0, {
        for(i in seq_along(lambda_seq)) {
          df_gam <- data.frame(y = y_valid, t = t_valid)
          
          # Fit with fixed lambda (sp = smoothing parameter)
          gam_fit <- mgcv::gam(y ~ s(t, bs = "cr"), data = df_gam, 
                               method = "REML", sp = lambda_seq[i])
          
          reml_scores[i] <- gam_fit$gcv.ubre
          edf_values[i] <- summary(gam_fit)$edf
          
          incProgress(1 / length(lambda_seq))
        }
      })
      
      # Store results
      values$reml_profile <- list(
        lambda = lambda_seq,
        reml_score = reml_scores,
        edf = edf_values,
        optimal_idx = which.min(reml_scores),
        optimal_lambda = lambda_seq[which.min(reml_scores)],
        optimal_reml = min(reml_scores)
      )
      
      showNotification(
        sprintf("REML profile complete! Optimal lambda: %.2e", 
                values$reml_profile$optimal_lambda),
        type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(paste("REML profile error:", e$message), type = "error", duration = 5)
      cat("REML profile error:", e$message, "\n")
    })
  })
  
  output$reml_profile_plot <- renderPlotly({
    if(is.null(values$reml_profile)) {
      return(plot_ly() %>% 
               layout(title = "No REML profile computed yet",
                      annotations = list(text = "Click 'Compute REML Profile' to run analysis",
                                         showarrow = FALSE)))
    }
    
    profile <- values$reml_profile
    
    # Create plot
    p <- plot_ly() %>%
      add_trace(x = log10(profile$lambda), 
                y = profile$reml_score,
                type = 'scatter', 
                mode = 'lines+markers',
                name = 'REML Score',
                line = list(color = 'blue', width = 2),
                marker = list(size = 6)) %>%
      add_trace(x = log10(profile$optimal_lambda),
                y = profile$optimal_reml,
                type = 'scatter',
                mode = 'markers',
                name = 'Optimal',
                marker = list(size = 12, color = 'red', symbol = 'star')) %>%
      layout(
        title = "REML Profile: Optimal Smoothing Parameter",
        xaxis = list(title = "log10(Lambda)"),
        yaxis = list(title = "REML Score (lower is better)"),
        hovermode = 'closest',
        showlegend = TRUE
      )
    
    # Add annotation for optimal value
    p <- p %>% layout(
      annotations = list(
        x = log10(profile$optimal_lambda),
        y = profile$optimal_reml,
        text = sprintf("Optimal lambda = %.2e<br>EDF = %.2f", 
                       profile$optimal_lambda,
                       profile$edf[profile$optimal_idx]),
        showarrow = TRUE,
        arrowhead = 2,
        ax = 40,
        ay = -40
      )
    )
    
    p
  })
  
  # 3. Cross-Validation Analysis
  observeEvent(input$run_cv_analysis, {
    req(values$data)
    
    showNotification("Running K-fold cross-validation...", type = "message", duration = 3)
    
    tryCatch({
      n_time <- ncol(values$data)
      n_subjects <- nrow(values$data)
      # P9.3: the same time coordinates production smooths against -- real
      # elapsed hours when the user asked for them, the column index otherwise.
      cv_axis     <- fck_smoothing_axis(input, values)
      time_points <- cv_axis$t_full
      k_folds <- input$cv_k_folds
      
      # Lambda range
      lambda_seq <- 10^seq(input$lambda_min_exp, input$lambda_max_exp, 
                           length.out = input$n_lambda)
      
      # Create fold assignments
      if(input$cv_stratified && !is.null(values$group_labels)) {
        # Stratified by groups
        fold_assignments <- rep(NA, n_subjects)
        for(group in unique(values$group_labels)) {
          group_idx <- which(values$group_labels == group)
          n_group <- length(group_idx)
          fold_assignments[group_idx] <- sample(rep(1:k_folds, length.out = n_group))
        }
      } else {
        # Random assignment
        fold_assignments <- sample(rep(1:k_folds, length.out = n_subjects))
      }
      
      # CV error matrix: subjects x lambdas
      cv_errors <- matrix(NA, nrow = n_subjects, ncol = length(lambda_seq))
      
      withProgress(message = "Running cross-validation...", value = 0, {
        total_iter <- n_subjects * length(lambda_seq)
        iter_count <- 0
        
        for(lambda_idx in seq_along(lambda_seq)) {
          lambda <- lambda_seq[lambda_idx]
          
          # AUDIT (P9.3): this used to build its own basis --
          #     create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
          # on the integer column index -- and the report then told the user to
          # type its smoothing factor into Data Preprocessing. That is the
          # production basis only in the default case. Under cyclic smoothing
          # production fits a FOURIER basis; on real clock times it fits over
          # elapsed hours on an uneven grid, where the roughness penalty is per
          # hour. In either case this lambda was a weight on a different penalty
          # over a different basis, so the advice to transfer it was wrong for
          # the same reason the mgcv/fda ratio was (P8.3), one level less
          # obvious. It now builds the SAME object production does.
          # P10.2: the same count rule as production, not a second one. This
          # capped at n_time - 2 where production caps at n_time, so on a short
          # series the CV was fitting a smaller basis than the app -- which made
          # its lambda not quite comparable after all.
          nb <- fck_smoothing_nbasis(input, n_time)
          basis <- fck_smoothing_basis(cv_axis, nb, input$smooth_method %||% "manual")
          fdParobj <- fdPar(basis, 2, lambda)
          
          for(i in 1:n_subjects) {
            # Get subject's fold
            test_fold <- fold_assignments[i]
            train_idx <- which(fold_assignments != test_fold)
            
            # Train on other subjects
            train_data <- values$data[train_idx, , drop = FALSE]
            
            # Get mean of training data (or could use all training curves)
            train_mean <- colMeans(train_data, na.rm = TRUE)
            
            # Smooth training mean
            valid_train <- !is.na(train_mean)
            if(sum(valid_train) >= 4) {
              fd_train <- smooth.basis(time_points[valid_train], 
                                       train_mean[valid_train], 
                                       fdParobj)$fd
              
              # Predict on test subject
              y_test <- values$data[i, ]
              valid_test <- !is.na(y_test)
              
              if(sum(valid_test) >= 1) {
                pred_test <- eval.fd(time_points[valid_test], fd_train)
                
                # Compute prediction error
                cv_errors[i, lambda_idx] <- sqrt(mean((y_test[valid_test] - pred_test)^2))
              }
            }
            
            iter_count <- iter_count + 1
            if(iter_count %% 50 == 0) {
              incProgress(50 / total_iter)
            }
          }
        }
      })
      
      # Compute mean CV error and SE for each lambda
      mean_cv_error <- colMeans(cv_errors, na.rm = TRUE)
      se_cv_error <- apply(cv_errors, 2, sd, na.rm = TRUE) / sqrt(n_subjects)
      
      optimal_idx <- which.min(mean_cv_error)
      optimal_lambda <- lambda_seq[optimal_idx]
      
      # 1-SE rule: most parsimonious model within 1 SE of minimum
      se_threshold <- mean_cv_error[optimal_idx] + se_cv_error[optimal_idx]
      within_1se <- which(mean_cv_error <= se_threshold)
      lambda_1se <- lambda_seq[max(within_1se)]  # Highest lambda (most smooth) within 1 SE
      
      # Store results
      values$cv_results <- list(
        lambda = lambda_seq,
        mean_error = mean_cv_error,
        se_error = se_cv_error,
        optimal_idx = optimal_idx,
        optimal_lambda = optimal_lambda,
        lambda_1se = lambda_1se,
        k_folds = k_folds,
        # P9.3: which model this lambda belongs to, so the report can say it
        # rather than leave the reader to assume.
        # P10.2: the same count rule again, not a third copy of it.
        basis_label = tryCatch(
          fck_basis_label(cv_axis, fck_smoothing_basis(
            cv_axis, fck_smoothing_nbasis(input, n_time),
            input$smooth_method %||% "manual")),
          error = function(e) NA_character_)
      )
      
      showNotification(
        sprintf("CV complete! Optimal lambda: %.2e (1-SE: %.2e)", 
                optimal_lambda, lambda_1se),
        type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(paste("CV error:", e$message), type = "error", duration = 5)
      cat("CV error:", e$message, "\n")
    })
  })
  
  output$cv_curve_plot <- renderPlotly({
    if(is.null(values$cv_results)) {
      return(plot_ly() %>% 
               layout(title = "No CV analysis performed yet",
                      annotations = list(text = "Click 'Run Cross-Validation' to start",
                                         showarrow = FALSE)))
    }
    
    cv <- values$cv_results
    
    # Create plot with error bars
    p <- plot_ly() %>%
      add_trace(x = log10(cv$lambda), 
                y = cv$mean_error,
                type = 'scatter', 
                mode = 'lines',
                name = 'Mean CV Error',
                line = list(color = 'darkgreen', width = 2)) %>%
      add_trace(x = log10(cv$lambda),
                y = cv$mean_error + cv$se_error,
                type = 'scatter',
                mode = 'lines',
                name = 'Mean + SE',
                line = list(color = 'lightgreen', width = 1, dash = 'dash'),
                showlegend = FALSE) %>%
      add_trace(x = log10(cv$lambda),
                y = cv$mean_error - cv$se_error,
                type = 'scatter',
                mode = 'lines',
                name = 'Mean - SE',
                line = list(color = 'lightgreen', width = 1, dash = 'dash'),
                fill = 'tonexty',
                fillcolor = 'rgba(144, 238, 144, 0.2)',
                showlegend = FALSE) %>%
      add_trace(x = log10(cv$optimal_lambda),
                y = cv$mean_error[cv$optimal_idx],
                type = 'scatter',
                mode = 'markers',
                name = 'Optimal (min)',
                marker = list(size = 12, color = 'red', symbol = 'star')) %>%
      add_trace(x = log10(cv$lambda_1se),
                y = cv$mean_error[which(abs(cv$lambda - cv$lambda_1se) < 1e-10)],
                type = 'scatter',
                mode = 'markers',
                name = '1-SE rule',
                marker = list(size = 12, color = 'orange', symbol = 'diamond')) %>%
      layout(
        title = sprintf("Cross-Validation Curve (K=%d folds)", cv$k_folds),
        xaxis = list(title = "log10(Lambda)"),
        yaxis = list(title = "CV Error (RMSE)"),
        hovermode = 'closest',
        showlegend = TRUE
      )
    
    p
  })
  
  # 5. Comparison Summary
  output$smoothing_comparison_summary <- renderPrint({
    if(is.null(values$reml_profile) && is.null(values$cv_results)) {
      cat("No smoothing diagnostic analysis performed yet.\n\n")
      cat("Run REML Profile and/or Cross-Validation to see comparison.")
      return()
    }
    
    cat("=== Smoothing Parameter Comparison ===\n\n")
    
    if(!is.null(values$reml_profile)) {
      cat("REML Analysis:\n")
      cat(sprintf("  Optimal lambda: %.3e\n", values$reml_profile$optimal_lambda))
      cat(sprintf("  Optimal EDF: %.2f\n", 
                  values$reml_profile$edf[values$reml_profile$optimal_idx]))
      cat(sprintf("  REML score: %.2f\n\n", values$reml_profile$optimal_reml))
    }
    
    if(!is.null(values$cv_results)) {
      cat("Cross-Validation Analysis:\n")
      cat(sprintf("  Optimal lambda (min CV): %.3e\n", values$cv_results$optimal_lambda))
      cat(sprintf("  Lambda (1-SE rule): %.3e\n", values$cv_results$lambda_1se))
      cat(sprintf("  Min CV error: %.3f\n", 
                  values$cv_results$mean_error[values$cv_results$optimal_idx]))
      cat(sprintf("  K-folds: %d\n\n", values$cv_results$k_folds))
    }
    
    # AUDIT (P8.3): this used to print
    #     ratio_min <- cv_lambda / reml_lambda
    #     "CV optimal / REML optimal: 1.20"
    #     "-> REML and CV agree well on optimal smoothing"
    # THE TWO LAMBDAS ARE NOT ON THE SAME SCALE, so their ratio means nothing
    # and "agree" was never a finding. The CV lambda is an fda penalty weight
    # (the CV loop calls fdPar(basis, 2, lambda) and smooth.basis), so it IS
    # comparable with the production smoother and IS transferable to the
    # manual smoothing factor. The REML lambda is mgcv's, on mgcv's own
    # penalty for a different basis: it is an advisory diagnostic about how
    # much structure the data supports, and it cannot be carried across.
    # Its "-log10 -> Smoothing Factor" conversion was removed for the same
    # reason: it produced a number the user could type into a control where it
    # means something else.
    if(!is.null(values$reml_profile) && !is.null(values$cv_results)) {
      cat("Comparing the two:\n")
      cat("  No ratio is reported. The mgcv REML optimum and the fda CV optimum\n")
      cat("  are penalty weights on DIFFERENT penalties over DIFFERENT bases;\n")
      cat("  their ratio is not a measure of agreement and a value near 1 would\n")
      cat("  not mean the two methods concur. Read the REML panel for the\n")
      cat("  effective degrees of freedom the data supports, and the CV panel\n")
      cat("  for a lambda on the scale this app actually smooths with.\n\n")
    }

    cat("\nWhat to take from each panel:\n")
    if(!is.null(values$cv_results)) {
      # P8.4: name the estimand. This CV smooths the TRAINING-GROUP MEAN and
      # scores it against a held-out subject, so it answers "how much should the
      # population curve be smoothed to predict a new person" -- NOT "how much
      # should this person's own trajectory be smoothed", which is what the
      # production smoother does. Recommending its lambda for per-subject
      # smoothing was recommending the answer to a different question.
      cat(sprintf("  Between-subject population-curve prediction CV\n"))
      if (!is.null(values$cv_results$basis_label) &&
          !is.na(values$cv_results$basis_label))
        cat(sprintf("    fitted on: %s\n", values$cv_results$basis_label))
      cat(sprintf("    optimal lambda  = %.2e   (smoothing factor %.2f)\n",
                  values$cv_results$optimal_lambda,
                  -log10(values$cv_results$optimal_lambda)))
      cat(sprintf("    1-SE rule       = %.2e   (smoothing factor %.2f)\n",
                  values$cv_results$lambda_1se, -log10(values$cv_results$lambda_1se)))
      cat("    This is the smoothing that best predicts a HELD-OUT SUBJECT from\n")
      cat("    the mean of the others. It is fitted on exactly the\n")
      cat("    basis and time axis the production smoother uses, so the\n")
      cat("    smoothing factor IS on the same scale and can be typed into Data\n")
      cat("    Preprocessing. It still answers a different question from the\n")
      cat("    production smoother, which asks how much to smooth one person's\n")
      cat("    OWN trajectory, not how to predict a new person from the others.\n\n")
    }
    if(!is.null(values$reml_profile)) {
      cat("  mgcv REML diagnostic (advisory)\n")
      cat(sprintf("    optimal EDF     = %.2f\n",
                  values$reml_profile$edf[values$reml_profile$optimal_idx]))
      cat(sprintf("    REML optimum    = %.2e   (mgcv scale -- NOT transferable)\n",
                  values$reml_profile$optimal_lambda))
      cat("    Read the EDF: it says how many effective degrees of freedom the\n")
      cat("    data supports. Do NOT type this lambda into the smoothing factor.\n\n")
    }
    cat("  Production smoothing\n")
    cat("    Automatic mode selects lambda by GCV on the fda smoother\n")
    cat("    (fck_auto_lambda), per run, on this data. That is what the app\n")
    cat("    fits. The 'Use Diagnostic Results' button runs that same search.\n")
  })
  
  # ============================================================================
  # LINK DIAGNOSTICS TO DATA PREPROCESSING
  # ============================================================================
  
  # Check if diagnostic results are available
  output$diagnostics_available <- reactive({
    !is.null(values$reml_profile) || !is.null(values$cv_results)
  })
  outputOptions(output, "diagnostics_available", suspendWhenHidden = FALSE)
  
  # Apply a lambda to the smoothing factor.
  #
  # AUDIT (P1.7): this observer used to take the optimal smoothing parameter
  # from the REML profile (an mgcv::gam fit with s(t, bs="cr")) or from the
  # leave-one-subject-out CV, and hand it to fda::fdPar() as lambda. Those are
  # different quantities. mgcv's sp multiplies a penalty built from a cubic
  # regression-spline basis with mgcv's own scaling; fda's lambda multiplies an
  # integrated squared second-derivative penalty on a B-spline basis. The number
  # transferred cleanly; its meaning did not. The CV had a second problem: it
  # scored lambda for predicting a held-out subject from the population mean,
  # which is not the problem the smoothing module then solves (smoothing each
  # subject's own trajectory).
  #
  # Since P0.2 the app has a GCV search that uses the production smoother, so
  # the transfer is unnecessary as well as wrong. The button now runs that
  # search. The REML and CV panels remain as diagnostics about how much
  # structure the data support -- they no longer set the estimator's lambda.
  observeEvent(input$use_diagnostic_lambda, {
    req(values$data)
    dat    <- values$data
    n_time <- ncol(dat)
    cyclic <- isTRUE(input$is_cyclic)

    # AUDIT (P11.3): this built its own basis over rangeval = c(1, n_time) and
    # scored GCV against argvals = seq_len(n_time) -- the COLUMN INDEX -- while
    # production goes through fck_smoothing_axis(), which uses elapsed hours
    # whenever real-time smoothing is on. Measured on 14 unevenly spaced hourly
    # columns with real time active: production fits over [0, 23] and this fit
    # over [1, 14], and in cyclic mode production uses a Fourier basis of
    # PERIOD 24 while this one used period 13 -- a 13-hour rhythm fitted to
    # 24-hour data. The lambda that came back was then handed to the smoother as
    # if it belonged to it. It also ignored the user's n_basis in cyclic mode by
    # hardcoding min(n_time, 13) rather than deferring to the production rule.
    # One axis, one basis rule, one count rule; all three now come from the
    # shared helpers, which is the same correction as P9.3 and P10.2.
    axis  <- fck_smoothing_axis(input, values)
    nb    <- fck_smoothing_nbasis(input, n_time)
    basis <- fck_smoothing_basis(axis, nb, input$smooth_method)

    al <- tryCatch(fck_auto_lambda(dat, axis$t_full, basis,
                                   min_points_needed = if(cyclic) 3 else 4),
                   error = function(e) NULL)
    if(is.null(al)) {
      showNotification(
        "Could not score GCV for any subject with the current basis (too few observed points). Set a lambda by hand.",
        type = "warning", duration = 10)
      return()
    }

    smooth_factor <- max(0.1, min(10, -log10(al$lambda)))
    updateSliderInput(session, "smooth_factor", value = smooth_factor)
    updateRadioButtons(session, "smooth_method", selected = "manual")

    showNotification(
      sprintf(paste0("GCV search on the actual smoother: lambda = %.3g (smoothing factor %.2f), ",
                     "from %d subjects. Set as a manual value so you can adjust it; ",
                     "'Automatic smoothing (GCV)' performs the same search on every run."),
              al$lambda, smooth_factor, al$n_used),
      type = "message", duration = 12)
  })

  # ===== GCV vs n-BASIS SWEEP =====
  #
  # AUDIT (P3.1): this sweep used to open with
  #     lam <- if(input$smooth_method == "auto") 0 else 10^(-input$smooth_factor)
  # which is the same defect P0.2 removed from server/20_smoothing.R and from
  # the "suggest a lambda" observer above, left behind in the one place that
  # advertises itself as a diagnostic. lambda = 0 in fda is UNPENALISED, not
  # "chosen automatically": under it GCV is minimised by whatever basis
  # interpolates least badly, so in auto mode the curve this tab drew was the
  # GCV of a model the app would never fit, and its "optimal n_basis" was the
  # optimum for an unpenalised fit. With a penalty selected per candidate
  # basis, the curve is normally much flatter -- which is the real finding:
  # once lambda is chosen by GCV, n_basis stops mattering above some floor.
  observeEvent(input$run_nbasis_analysis, {
    req(values$data)
    data_mat <- values$data
    n_time   <- ncol(data_mat)
    auto     <- identical(input$smooth_method, "auto")
    lam_manual <- if(auto) NA_real_ else 10^(-input$smooth_factor)

    # AUDIT (P4.8): this sweep always built a B-spline basis. When the user has
    # selected CYCLIC smoothing the production smoother fits a FOURIER basis, so
    # for periodic data the diagnostic answered a question about a model the app
    # was not fitting -- "which B-spline count is best?" while the analysis ran
    # "Fourier basis at the selected lambda". The "suggest a lambda" observer
    # above already branches on input$is_cyclic; this one now does too, and the
    # sweep runs over the number of Fourier basis functions, which must be odd
    # (a constant plus sin/cos pairs) -- fda silently rounds an even count up,
    # so an even grid would have scored the same model twice.
    #
    # AUDIT (P11.3): P4.8 fixed the basis TYPE and left the AXIS hand-built, so
    # the sweep still ran over rangeval = c(1, n_time) with the column index as
    # argvals. With real-time smoothing on that is a different domain from the
    # one production fits (measured: [1, 14] against [0, 23]) and, for cyclic
    # data, a different PERIOD (13 against 24). A sweep that recommends an
    # n_basis for a model the app will not fit is the same defect P4.8 and P3.1
    # each removed from this observer; the axis now comes from the shared
    # builder, so there is nothing left in this file that constructs a basis.
    cyclic <- isTRUE(input$is_cyclic)
    axis   <- fck_smoothing_axis(input, values)
    # nb_fourier is the one sanctioned override: the sweep must vary the count
    # to have a curve to draw. Everything else -- rangeval, period, basis type --
    # is production's.
    make_basis <- function(nb) {
      if (cyclic) fck_smoothing_basis(axis, nb, input$smooth_method, nb_fourier = nb)
      else        fck_smoothing_basis(axis, nb, "manual")
    }

    nb_seq <- if (cyclic) seq(3, min(n_time - 1, 41), by = 2)
              else        seq(4, min(n_time - 1, 40), by = 2)
    gcv_means <- numeric(length(nb_seq))
    lam_used  <- numeric(length(nb_seq))

    withProgress(message = if(auto)
                   "GCV for each n-basis, with lambda re-selected each time..."
                 else "Computing GCV for each n-basis...", value = 0, {
      for(bi in seq_along(nb_seq)) {
        incProgress(1 / length(nb_seq), detail = paste("n_basis =", nb_seq[bi]))
        nb    <- nb_seq[bi]
        basis <- make_basis(nb)

        if(auto) {
          # The same GCV search the app performs on every run, repeated for
          # this candidate basis. A coarser grid and a smaller subject cap
          # than the single-shot search: this runs once per basis size, and
          # GCV as a function of log-lambda is smooth enough that a 12-point
          # bracket plus the refinement lands in the same place.
          al <- tryCatch(fck_auto_lambda(data_mat, axis$t_full, basis,
                                         min_points_needed = if(cyclic) 3 else 4,
                                         n_grid = 12),
                         error = function(e) NULL)
          if(is.null(al)) { gcv_means[bi] <- NA_real_; lam_used[bi] <- NA_real_; next }
          lam_used[bi]  <- al$lambda
          gcv_means[bi] <- al$gcv
          next
        }

        lam_used[bi] <- lam_manual
        fdP   <- fdPar(basis, 2, lam_manual)
        min_pts <- if(cyclic) 3 else 4
        gcv_sub <- sapply(seq_len(nrow(data_mat)), function(s) {
          y     <- data_mat[s, ]
          valid <- which(!is.na(y))
          if(length(valid) < min_pts) return(NA_real_)
          # P11.3: argvals must be the SAME axis the basis was built over.
          # These were column indices against a basis now spanning elapsed
          # hours, which is not a mismatch fda would report -- it would simply
          # fit the wrong thing.
          sb <- tryCatch(smooth.basis(axis$t_full[valid], y[valid], fdP),
                         error = function(e) NULL)
          if(is.null(sb)) NA_real_ else sb$gcv
        })
        gcv_means[bi] <- mean(gcv_sub, na.rm = TRUE)
      }
    })

    if(all(!is.finite(gcv_means))) {
      showNotification(
        "No n-basis could be scored: every subject failed to smooth. Check for subjects with too few observed points.",
        type = "error", duration = 12)
      return()
    }

    optimal_idx <- which.min(gcv_means)
    values$nbasis_gcv <- list(
      nb       = nb_seq,
      gcv      = gcv_means,
      lambda   = lam_used,
      auto     = auto,
      cyclic   = cyclic,
      optimal  = nb_seq[optimal_idx]
    )
    showNotification(
      if(auto)
        paste0("n-basis analysis complete. Lowest mean GCV at n_basis = ",
               nb_seq[optimal_idx], " (lambda re-selected by GCV at each basis size; ",
               sprintf("%.3g", lam_used[optimal_idx]), " there). ",
               "The curve is usually flat -- with lambda chosen, n_basis matters ",
               "little above the floor.")
      else
        paste0("n-basis analysis complete. Optimal n_basis = ", nb_seq[optimal_idx],
               " (lowest mean GCV at the fixed lambda ",
               sprintf("%.3g", lam_manual), ")."),
      type = "message", duration = 12)
  })

  output$nbasis_gcv_plot <- renderPlotly({
    req(values$nbasis_gcv)
    ng         <- values$nbasis_gcv
    current_nb <- if(!is.null(input$smooth_method) && input$smooth_method == "manual") {
      input$n_basis_manual
    } else {
      input$n_basis
    }

    p <- plot_ly() %>%
      add_trace(
        x = ng$nb, y = ng$gcv,
        type = "scatter", mode = "lines+markers",
        line   = list(color = "#2196F3", width = 2),
        marker = list(color = "#2196F3", size = 6),
        name   = "Mean GCV"
      ) %>%
      layout(
        title  = list(
          text = paste0(
            if(isTRUE(ng$cyclic))
              "GCV Score vs Number of Fourier Basis Functions"
            else
              "GCV Score vs Number of B-spline Basis Functions",
            if(isTRUE(ng$auto))
              "<br><sub>lambda re-selected by GCV at every basis size (auto mode)</sub>"
            else
              paste0("<br><sub>at the fixed lambda ",
                     sprintf("%.3g", ng$lambda[1]), " (manual mode)</sub>")),
          x = 0.5),
        xaxis  = list(title = "n_basis"),
        yaxis  = list(title = "Mean GCV (lower = better)"),
        shapes = list(
          list(type = "line", x0 = ng$optimal, x1 = ng$optimal,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "green", dash = "dash", width = 1.5)),
          list(type = "line", x0 = current_nb, x1 = current_nb,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "orange", dash = "dot", width = 1.5))
        ),
        annotations = list(
          list(x = ng$optimal, y = 1, yref = "paper", xanchor = "left",
               text = paste0("Optimal (", ng$optimal, ")"),
               showarrow = FALSE, font = list(color = "green", size = 11)),
          list(x = current_nb, y = 0.88, yref = "paper", xanchor = "left",
               text = paste0("Current (", current_nb, ")"),
               showarrow = FALSE, font = list(color = "orange", size = 11))
        ),
        legend = list(orientation = "h", y = -0.2)
      )
    p
  })


  # ============================================================================
  # END DIAGNOSTICS LINK
  # ============================================================================
