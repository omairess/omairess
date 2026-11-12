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
    # they are not (see DANCE/server/03_helpers_clock.R).
    #
    # AUDIT (P11.3): this block used to derive the axis itself, and the
    # diagnostics derived it again, and the two drifted -- production ended up
    # fitting over elapsed hours while the diagnostic that recommends a lambda
    # for it fitted over the column index. The derivation now happens ONCE, in
    # dance_smoothing_axis(). Production keeps the notification, because it is the
    # place the user made the request; the helper stays silent so it can be
    # called from a diagnostic without raising notices of its own.
    axis            <- dance_smoothing_axis(input, values)
    using_real_time <- axis$using_real_time
    t_full          <- axis$t_full
    t_rng           <- axis$t_rng
    if(isTRUE(input$use_real_time) && !using_real_time) {
      showNotification(
        "Real clock times requested, but the column names do not yield hours in [0, 24). Smoothing against the column index instead.",
        type = "warning", duration = 10)
    }

    # Where every value in the smoothed curves will come from. Computed from the
    # RAW data, before anything is filled, and kept for the whole session so any
    # later plot or export can say which points are measurements.
    values$fill_status <- dance_fill_status(values$data)
    if(using_real_time) {
      showNotification(
        sprintf("Smoothing against real elapsed time (%.2f to %.2f h). The roughness penalty is now per hour, so the smoothing factor may need re-tuning.",
                t_rng[1], t_rng[2]),
        type = "message", duration = 8)
    }

    # ---- basis size ---------------------------------------------------------
    # P10.2: the rule lives in dance_smoothing_nbasis() alongside the basis
    # constructor, so the cross-validation diagnostic caps identically. The
    # rule itself is unchanged: 20 or n_time - 2 for an unpenalised fit, the
    # user's count capped at n_time otherwise, never below 4.
    nb <- dance_smoothing_nbasis(input, n_time)

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
    # P9.3: built by dance_smoothing_basis() in server/04_helpers_fd.R, so the
    # cross-validation diagnostic can construct the SAME object and its lambda
    # is genuinely on this scale. The rules are unchanged: Fourier when cyclic
    # (CIRCAREG's min(n_time, 13), period 24 h on real elapsed time), an
    # interpolating B-spline through every observation when smoothing is off,
    # and an nb-function B-spline otherwise.
    basis   <- dance_smoothing_basis(axis, nb, input$smooth_method)
    if(cyclic && input$smooth_method != "none")
      showNotification(
        sprintf("Cyclic smoothing: one cycle = %s %s (%s).",
                format(signif(basis$params[1], 6)),
                if(using_real_time) "h" else "columns", axis$period_source),
        type = "message", duration = 10)
    nb_used <- if(cyclic) basis$nbasis else nb

    # THE COMPUTATION is dance_smooth_curves() (server/04_helpers_fd.R), the
    # function the exported script calls too: per-curve smoothing on the axis
    # above, gaps filled by the fit, the ends held flat unless extrapolation is
    # allowed, optional bounds, and the fd object re-expressed on 0-1 with the
    # same basis size and (Fourier) the same cycle (P19.3/P19.4). This block
    # only reports what it did.
    bounds <- if(isTRUE(input$constrain_bounds)) c(input$min_bound, input$max_bound) else NULL
    if(input$smooth_method == "none") {
      # ---- raw data, no smoothing -------------------------------------------
      # fd_obj still has to exist for the fPCA family: the data points joined by
      # straight lines, a missing value read off the line between its
      # neighbours for that representation only (values$smooth_data keeps the
      # real NAs).
      sm <- dance_smooth_curves(values$data, axis, basis, "none")
      values$smooth_data <- sm$smooth_mat
      values$fd_obj <- sm$fd_obj
      values$smooth_fit_metrics <- list(
        method        = "none",
        n_basis       = ncol(values$data),
        basis_type    = "bspline order 2 (the data points joined by straight lines)",
        time_axis     = if(using_real_time) "real clock time (hours)" else "column index",
        fill_headline = dance_fill_headline(values$fill_status),
        extrapolation = paste("missing values read off the straight line between the curve's neighbouring",
                              "observations (the nearest observation before the first or after the last),",
                              "for the curve representation only; values$smooth_data keeps the real NAs"),
        settings      = list(method = "none", nb = nb, use_real_time = using_real_time,
                             cyclic = cyclic, allow_extrapolation = FALSE, bounds = NULL)
      )
      values$smoothing_avg_metrics <- NULL
      showNotification(
        "Using raw data (no smoothing): each curve is your data points joined by straight lines.",
        type = "message", duration = 8)

    } else {
      # AUDIT (P0.2): "auto" is a real GCV search with the same estimator that
      # does the smoothing (dance_auto_lambda), not lambda = 0.
      sm <- withProgress(message = 'Smoothing subjects...', value = 0,
        dance_smooth_curves(values$data, axis, basis, input$smooth_method,
                            lambda = 10^(-input$smooth_factor),
                            allow_extrapolation = isTRUE(input$allow_extrapolation),
                            bounds = bounds,
                            progress = function(i, n) if(i %% 10 == 0) incProgress(10 / n)))
      if(identical(input$smooth_method, "auto")) {
        if(isTRUE(sm$auto_failed))
          showNotification(
            "Automatic smoothing could not evaluate GCV for any subject (too few observed points). Falling back to an unpenalised basis fit; set a lambda by hand on the Manual setting if you want a smooth.",
            type = "warning", duration = 12)
        else
          showNotification(
            sprintf("Automatic smoothing (GCV): lambda = %s, chosen by minimising the mean GCV score across %d subjects.",
                    format(sm$lambda, digits = 3), sm$lambda_n_used),
            type = "message", duration = 10)
      }
      if(sm$n_clamped > 0)
        cat(sprintf("Held %d extrapolated values flat at each subject's last observation.\n", sm$n_clamped))
      smooth_mat      <- sm$smooth_mat
      lam             <- sm$lambda
      lambda_n_used   <- sm$lambda_n_used
      df_vec          <- sm$df_vec
      gcv_vec         <- sm$gcv_vec
      failed_subjects <- sm$failed_subjects
      r_squared_vec   <- sm$r_squared
      rmse_vec        <- sm$rmse
      values$smooth_data <- smooth_mat

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
        lambda_n_used  = lambda_n_used,
        n_subjects     = n_subjects,
        method         = input$smooth_method,
        basis_type     = if(cyclic) "fourier" else "bspline",
        cyclic_period  = if(cyclic) basis$params[1] else NA_real_,
        cyclic_period_source = if(cyclic) axis$period_source else NA_character_,
        time_axis      = if(using_real_time) "real clock time (hours)" else "column index",
        mean_df        = mean(df_vec, na.rm = TRUE),
        sd_df          = sd(df_vec,   na.rm = TRUE),
        max_df         = max(df_vec,  na.rm = TRUE),
        mean_gcv       = mean(gcv_vec, na.rm = TRUE),
        rel_rmse_pct   = rel_rmse_pct,
        data_range     = data_range,
        fill_headline  = dance_fill_headline(values$fill_status),
        extrapolation  = if(isTRUE(input$allow_extrapolation))
          "carried past each subject's observed range by the fitted curve"
        else "held flat at each subject's first/last observation",
        # the settings the exported script passes back to dance_smooth_curves()
        settings       = list(method = input$smooth_method, nb = nb, lambda = lam,
                              use_real_time = using_real_time, cyclic = cyclic,
                              allow_extrapolation = isTRUE(input$allow_extrapolation),
                              bounds = bounds)
      )

      values$smoothing_avg_metrics <- c(
        R_squared   = mean(r_squared_vec, na.rm = TRUE),
        RMSE        = mean(rmse_vec, na.rm = TRUE),
        Correlation = NA,
        MAE         = NA
      )

      values$fd_obj <- sm$fd_obj

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

    # Say how much of what just came out of the smoother is not measurement.
    if(!is.null(values$fill_status)) {
      extrap_frac <- mean(values$fill_status == DANCE_FILL_EXTRAPOLATED)
      showNotification(
        paste0(dance_fill_headline(values$fill_status),
               ". See 'Missing data & filled points' below."),
        type = if(extrap_frac > 0.1) "warning" else "message", duration = 12)
    }

    # Re-smoothing changes the curves under EVERY analysis. In the two source
    # apps a stale fPCA or cosinor fit could survive a re-smooth; with one
    # smoothing step feeding seven analyses that would be much easier to miss,
    # so results are cleared and have to be re-run against the new curves.
    dance_reset_analyses(values, keep_smoothing = TRUE)

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
