# ==============================================================================
# server/20_smoothing.R — THE shared smoothing step
#
# Hand-merged union of:
#   WaPaa1_3.R  1971-2229   observeEvent(input$apply_smooth, ...)
#   CIRCAREG.R   932-1097   observeEvent(input$apply_smooth, ...)
#
# These two were already the SAME algorithm — WaPaa's own comments say
# "Applying smoothing using CIRCAREG method" and "CIRCAREG EXACT METHOD" — so
# there is no statistical choice to make between them, only a union of
# features to take:
#
#   from WaPaa    per-subject EDF and GCV from smooth.basis; relative RMSE as
#                 a % of the data range; separate n_basis for auto vs manual;
#                 and the two-stage fd_obj (smoothing done on 1:n_time for
#                 numerical stability, then re-expressed on 0-1 because fPCA,
#                 warping, fANOVA and clustering all assume a 0-1 range)
#   from CIRCAREG the cyclic option: a Fourier basis for 24-hour data, and the
#                 lower minimum-points threshold that goes with it
#
# Everything else is identical in both: per-subject smoothing so that missing
# values are interpolated rather than dropped, linear-interpolation fallback
# when a subject has too few points or the fit fails, optional clamping to a
# range, and per-subject R2/RMSE against the ORIGINAL (unsmoothed) values.
#
# This runs once. values$smooth_data and values$fd_obj are what every analysis
# tab in the app reads, so the curves behind an fPCA, a cosinor fit and a
# cluster solution are guaranteed to be the same curves.
# ==============================================================================

observeEvent(input$apply_smooth, {
  req(values$data)

  cat("Applying shared smoothing step...\n")

  tryCatch({
    n_time     <- ncol(values$data)
    n_subjects <- nrow(values$data)
    cyclic     <- isTRUE(input$is_cyclic)

    # ---- the time axis to smooth against ------------------------------------
    # Both source apps used the column index, which is only correct when the
    # columns are evenly spaced in real time. Opt in to real clock times when
    # they are not (see FCK/server/03_helpers_clock.R).
    real_time <- NULL
    if(isTRUE(input$use_real_time)) {
      real_time <- fck_cumulative_hours(values$time_labels)
      if(is.null(real_time) || length(real_time) != n_time) {
        real_time <- NULL
        showNotification(
          "Real clock times requested, but the column names do not yield hours in [0, 24). Smoothing against the column index instead.",
          type = "warning", duration = 10)
      }
    }
    using_real_time <- !is.null(real_time)
    t_full <- if(using_real_time) real_time else 1:n_time   # smoothing argument
    t_rng  <- range(t_full)

    if(using_real_time) {
      showNotification(
        sprintf("Smoothing against real elapsed time (%.2f to %.2f h). The roughness penalty is now per hour, so the smoothing factor may need re-tuning.",
                t_rng[1], t_rng[2]),
        type = "message", duration = 8)
    }

    # ---- basis size ---------------------------------------------------------
    if(input$smooth_method == "none") {
      nb <- min(20, n_time - 2)
    } else if(input$smooth_method == "manual") {
      nb <- input$n_basis_manual
      nb <- min(nb, n_time)
      if(nb < 4) nb <- 4
    } else {
      nb <- input$n_basis
      nb <- min(nb, n_time)
      if(nb < 4) nb <- 4
    }

    if(!cyclic && nb >= n_time - 2 && input$smooth_method != "none") {
      showNotification(
        sprintf("Number of B-splines (%d) is very high relative to time points (%d). Consider reducing it.",
                nb, n_time),
        type = "warning", duration = 5
      )
    }

    # ---- basis --------------------------------------------------------------
    # Cyclic: CIRCAREG's Fourier basis with min(n_time, 13) functions. fda
    # rounds an even nbasis up to the next odd number; that is the behaviour
    # CIRCAREG shipped, so it is kept rather than silently changed here.
    if(cyclic) {
      nb_used <- min(n_time, 13)
      # On real elapsed hours the period is 24 h, not "however long the
      # recording happened to be" (fda's default period = diff(rangeval)).
      basis   <- if(using_real_time)
        create.fourier.basis(rangeval = t_rng, nbasis = nb_used, period = 24)
      else
        create.fourier.basis(rangeval = t_rng, nbasis = nb_used)
      nb_used <- basis$nbasis
    } else if(input$smooth_method == "none") {
      basis   <- create.bspline.basis(rangeval = t_rng, breaks = t_full, norder = 4)
      nb_used <- nb
    } else {
      basis   <- create.bspline.basis(rangeval = t_rng, nbasis = nb)
      nb_used <- nb
    }

    # The 0-1 basis used for the fd object every downstream analysis reads.
    make_basis_01 <- function(k) {
      if(cyclic) create.fourier.basis(rangeval = c(0, 1), nbasis = k)
      else       create.bspline.basis(rangeval = c(0, 1), nbasis = k)
    }
    # The fPCA / warping / fANOVA / clustering family works on a 0-1 range. With
    # real times that rescaling has to preserve the spacing, otherwise the fd
    # object would put the columns back on an even grid and undo the point of
    # the option. (This matches WaPaa's calculate_time_positions(), which
    # normalises the same cumulative hours for the plot axes.)
    time_points_01 <- if(using_real_time) (t_full - t_rng[1]) / diff(t_rng)
                      else seq(0, 1, length.out = n_time)

    if(input$smooth_method == "none") {
      # ---- raw data, no smoothing -------------------------------------------
      values$smooth_data <- values$data

      # fd_obj still has to exist for the fPCA family; NAs are mean-imputed
      # for that representation only (values$smooth_data keeps the real NAs).
      temp_data <- values$data
      for(i in 1:n_subjects) {
        na_idx <- is.na(temp_data[i, ])
        if(any(na_idx)) {
          temp_data[i, na_idx] <- mean(temp_data[i, !na_idx], na.rm = TRUE)
        }
      }
      basis_01 <- make_basis_01(min(20, n_time - 2))
      values$fd_obj <- smooth.basis(time_points_01, t(temp_data), basis_01)$fd
      values$smooth_fit_metrics <- NULL
      values$smoothing_avg_metrics <- NULL
      showNotification("Using Raw Data (No Smoothing).", type = "message")

    } else {
      # ---- per-subject smoothing --------------------------------------------
      # lambda = 0 in "auto" mode: smooth.basis then uses REML/GCV internally
      # rather than a lambda search here. This is what BOTH source apps did.
      lam      <- if(input$smooth_method == "auto") 0 else 10^(-input$smooth_factor)
      fdParobj <- fdPar(basis, 2, lam)

      smooth_mat      <- matrix(NA, nrow = n_subjects, ncol = n_time)
      n_failed        <- 0
      failed_subjects <- c()
      df_vec          <- rep(NA_real_, n_subjects)   # effective df per subject
      gcv_vec         <- rep(NA_real_, n_subjects)   # GCV score per subject

      # A Fourier basis needs fewer points than a cubic B-spline (CIRCAREG).
      min_points_needed <- if(cyclic) 3 else 4

      withProgress(message = 'Smoothing subjects...', value = 0, {
        for(i in 1:n_subjects) {
          y_i       <- values$data[i, ]
          valid_idx <- !is.na(y_i)
          n_valid   <- sum(valid_idx)

          if(n_valid >= min_points_needed) {
            tryCatch({
              t_valid <- t_full[valid_idx]
              y_valid <- y_i[valid_idx]

              sb_i <- smooth.basis(t_valid, y_valid, fdParobj)

              # Evaluated at ALL time points: this is what interpolates the NAs.
              smooth_mat[i, ] <- as.vector(eval.fd(t_full, sb_i$fd))
              df_vec[i]       <- sb_i$df
              gcv_vec[i]      <- sb_i$gcv
            }, error = function(e) {
              smooth_mat[i, ] <<- approx(t_full[valid_idx], y_i[valid_idx],
                                         xout = t_full, rule = 2)$y
              n_failed <<- n_failed + 1
              failed_subjects <<- c(failed_subjects, i)
            })
          } else if(n_valid > 0) {
            smooth_mat[i, ] <- approx(t_full[valid_idx], y_i[valid_idx],
                                      xout = t_full, rule = 2)$y
          } else {
            n_failed <- n_failed + 1
            failed_subjects <- c(failed_subjects, i)
          }

          if(i %% 10 == 0) incProgress(10 / n_subjects)
        }
      })

      if(input$constrain_bounds) {
        min_val <- input$min_bound
        max_val <- input$max_bound
        if(!is.na(min_val) && !is.na(max_val) && min_val < max_val) {
          smooth_mat <- pmin(pmax(smooth_mat, min_val), max_val)
        }
      }

      values$smooth_data <- smooth_mat

      # ---- fit metrics: smoothed vs ORIGINAL, where the original was present -
      r_squared_vec <- numeric(n_subjects)
      rmse_vec      <- numeric(n_subjects)

      for(i in 1:n_subjects) {
        orig_i    <- values$data[i, ]
        smooth_i  <- smooth_mat[i, ]
        valid_idx <- !is.na(orig_i)

        if(sum(valid_idx) > 1) {
          orig_valid   <- orig_i[valid_idx]
          smooth_valid <- smooth_i[valid_idx]

          rmse_vec[i] <- sqrt(mean((orig_valid - smooth_valid)^2))

          ss_tot <- sum((orig_valid - mean(orig_valid))^2)
          ss_res <- sum((orig_valid - smooth_valid)^2)
          r_squared_vec[i] <- if(ss_tot > 0) 1 - ss_res / ss_tot else NA
        } else {
          r_squared_vec[i] <- NA
          rmse_vec[i]      <- NA
        }
      }

      data_range   <- max(values$data, na.rm = TRUE) - min(values$data, na.rm = TRUE)
      rel_rmse_pct <- if(data_range > 0) {
        (mean(rmse_vec, na.rm = TRUE) / data_range) * 100
      } else NA_real_

      values$smooth_fit_metrics <- list(
        r_squared      = r_squared_vec,
        rmse           = rmse_vec,
        mean_r_squared = mean(r_squared_vec, na.rm = TRUE),
        mean_rmse      = mean(rmse_vec, na.rm = TRUE),
        sd_r_squared   = sd(r_squared_vec, na.rm = TRUE),
        sd_rmse        = sd(rmse_vec, na.rm = TRUE),
        n_basis        = nb_used,
        lambda         = lam,
        method         = input$smooth_method,
        basis_type     = if(cyclic) "fourier" else "bspline",
        time_axis      = if(using_real_time) "real clock time (hours)" else "column index",
        mean_df        = mean(df_vec, na.rm = TRUE),
        sd_df          = sd(df_vec,   na.rm = TRUE),
        max_df         = max(df_vec,  na.rm = TRUE),
        mean_gcv       = mean(gcv_vec, na.rm = TRUE),
        rel_rmse_pct   = rel_rmse_pct,
        data_range     = data_range
      )

      values$smoothing_avg_metrics <- c(
        R_squared   = mean(r_squared_vec, na.rm = TRUE),
        RMSE        = mean(rmse_vec, na.rm = TRUE),
        Correlation = NA,
        MAE         = NA
      )

      # ---- fd object on 0-1 --------------------------------------------------
      # smooth_mat is ALREADY smoothed, so lambda = 0 here: this step only
      # re-expresses the same curves on the 0-1 range that fPCA, warping,
      # fANOVA and clustering assume. Smoothing twice would over-smooth.
      basis_01 <- make_basis_01(nb_used)
      values$fd_obj <- smooth.basis(time_points_01, t(smooth_mat),
                                    fdPar(basis_01, 2, 0))$fd

      cat(sprintf("Smoothed on %s with a %s basis (%d functions); fd object re-expressed on 0-1.\n",
                  if(using_real_time) sprintf("real elapsed hours %.2f-%.2f", t_rng[1], t_rng[2])
                  else sprintf("1:%d", n_time),
                  if(cyclic) "Fourier" else "B-spline", nb_used))

      # ---- report ------------------------------------------------------------
      n_with_na      <- sum(apply(values$data, 1, function(x) any(is.na(x))))
      n_interpolated <- n_with_na - length(failed_subjects)

      if(n_interpolated > 0) {
        showNotification(
          sprintf("Smoothing applied. %d subjects had missing values - successfully interpolated. Mean R2=%.3f, Mean RMSE=%.3f",
                  n_interpolated, values$smooth_fit_metrics$mean_r_squared,
                  values$smooth_fit_metrics$mean_rmse),
          type = "message", duration = 8)
      } else {
        showNotification(
          sprintf("Smoothing applied. Mean R2=%.3f (SD=%.3f), Mean RMSE=%.3f (SD=%.3f)",
                  values$smooth_fit_metrics$mean_r_squared, values$smooth_fit_metrics$sd_r_squared,
                  values$smooth_fit_metrics$mean_rmse, values$smooth_fit_metrics$sd_rmse),
          type = "message", duration = 8)
      }

      if(length(failed_subjects) > 0) {
        showNotification(
          sprintf("Warning: %d subjects could not be smoothed (all NA or too few points): %s",
                  length(failed_subjects), paste(head(failed_subjects, 10), collapse = ", ")),
          type = "warning", duration = 8)
      }

      cat(sprintf("\n=== Smoothing Results ===\n"))
      cat(sprintf("Mean R2: %.3f (SD: %.3f)\n",
                  values$smooth_fit_metrics$mean_r_squared,
                  values$smooth_fit_metrics$sd_r_squared))
      cat(sprintf("Mean RMSE: %.3f (SD: %.3f)\n",
                  values$smooth_fit_metrics$mean_rmse,
                  values$smooth_fit_metrics$sd_rmse))
    }

    # Re-smoothing changes the curves under EVERY analysis. In the two source
    # apps a stale fPCA or cosinor fit could survive a re-smooth; with one
    # smoothing step feeding seven analyses that would be much easier to miss,
    # so results are cleared and have to be re-run against the new curves.
    fck_reset_analyses(values, keep_smoothing = TRUE)

  }, error = function(e) {
    cat("Error in smoothing:", e$message, "\n")
    showNotification(paste("Smoothing error:", e$message), type = "error", duration = 10)

    # Fall back to the raw data so the app stays usable.
    values$smooth_data <- values$data
    n_basis_fb <- min(10, ncol(values$data) - 2)
    basis_fb   <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis_fb)

    fallback_data <- values$data
    if(typeof(fallback_data) != "double" && typeof(fallback_data) != "integer") {
      mode(fallback_data) <- "numeric"
    }

    values$fd_obj <- smooth.basis(seq(0, 1, length.out = ncol(values$data)),
                                  t(fallback_data), basis_fb)$fd
  })
})
