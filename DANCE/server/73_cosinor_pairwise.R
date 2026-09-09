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

    # ------------------------------------------------------------------------
    # MIXED COSINOR. One model over all observations, with a random rhythm per
    # participant. Shares this tab's group selection; the within factor is its
    # own, because the other two approaches have no notion of one.
    # ------------------------------------------------------------------------
    if (identical(input$hp_approach, "mixed")) {
      bn <- input$hp_mixed_between; wn <- input$hp_mixed_within
      if (is.null(values$subject_ids) || is.null(bn) || is.null(wn) ||
          !nzchar(bn) || !nzchar(wn)) {
        showNotification(paste("A mixed cosinor needs a participant identifier (Data Import tab)",
                               "and both a between- and a within-participant factor."),
                         type = "error", duration = 12); return()
      }
      axis <- dance_smoothing_axis(
        list(use_real_time = isTRUE(input$hp_mixed_real_time), is_cyclic = FALSE), values)
      d <- tryCatch(dance_mixed_long(values$data, axis$t_full,
                                     subject = values$subject_ids,
                                     between = values$covariates[[bn]],
                                     within  = values$covariates[[wn]]),
                    error = function(e) NULL)
      if (is.null(d)) { showNotification("Could not build the long form from those columns.",
                                         type = "error", duration = 10); return() }
      badm <- dance_mixed_check(d)
      if (length(badm)) { showNotification(paste(badm, collapse = " "),
                                           type = "error", duration = 15); return() }
      withProgress(message = "Fitting the mixed cosinor...", value = 0.4, {
        res <- dance_mixed_cosinor(d, period = mod$period,
                                   n_harmonics = mod$n_harmonics)
      })
      if (!isTRUE(res$ok)) { showNotification(res$message, type = "error", duration = 15); return() }
      res$kind <- "cosinor"; res$between_name <- bn; res$within_name <- wn
      res$time_axis <- if (isTRUE(input$hp_mixed_real_time)) "real elapsed time" else "column index"
      res$time_values <- axis$t_full
      values$mixed_results <- res
      showNotification("Mixed cosinor complete.", type = "message", duration = 4)
      return()
    }

    # ------------------------------------------------------------------------
    # POPULATION-MEAN COSINOR (Bingham et al. 1982). An alternative to the
    # two-stage route below, not a replacement: it compares ALL groups at once
    # on the coefficient PAIRS rather than comparing point estimates pairwise.
    # Handled first, and returns, so the two approaches share this tab's group
    # selection and nothing else.
    # ------------------------------------------------------------------------
    if (identical(input$hp_approach, "population")) {
      gv <- input$harmonic_group_var
      if (is.null(gv) || identical(gv, "_none_")) {
        showNotification("Select a grouping variable on the Harmonic Regression tab first.",
                         type = "error", duration = 8); return()
      }
      h <- max(1L, min(as.integer(input$hp_pop_harmonic %||% 1), mod$n_harmonics))
      pr <- mod$individual_params
      bcol <- paste0("beta_cos_", h); scol <- paste0("beta_sin_", h)
      if (!all(c(bcol, scol) %in% names(pr))) {
        showNotification(sprintf("Harmonic %d was not fitted, so its coefficients are not available.", h),
                         type = "error", duration = 8); return()
      }
      gvec <- values$covariates[[gv]][pr$subject]
      res <- dance_pop_cosinor(pr[[bcol]], pr[[scol]], pr$mesor, gvec,
                               period = mod$period, harmonic = h)
      if (!isTRUE(res$ok)) {
        showNotification(res$message, type = "error", duration = 15); return()
      }
      res$group_var <- gv
      values$pop_cosinor <- res
      # Bingham's caution is about the AMPLITUDE comparison, and the two-stage
      # readout already consumes this flag; keep them in step.
      values$hp_acrophase_differs <- res$acrophase_differs
      showNotification("Population-mean cosinor complete.", type = "message", duration = 4)
      return()
    }

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

    # AUDIT (P20/R3). Harmonic h completes h cycles per period, so its acrophase
    # lives on the EFFECTIVE period T/h -- and acrophase_time_h is already stored
    # on that scale (phi_to_hours() divides by h; server/72_harmonic.R:746).
    # This block converted back to radians with 2*pi/period regardless, so an H2
    # acrophase spanning 0-12 h was mapped onto 0-pi instead of a full turn, and
    # H3 onto a third of one. Everything downstream inherited the compression --
    # the circular mean, the circular SD, the shortest-arc difference, the
    # resultant lengths, and the Watson-Williams test run on them -- and the
    # damage goes BOTH ways depending on where the cluster sits, which is why it
    # is not a rescaling that cancels in a contrast:
    #
    #   A cluster that does not straddle the wrap is squeezed into half the
    #   circle, so it looks more concentrated than it is. Measured: two H2 groups
    #   near 3 h, r-bar 0.930 / 0.885 under the old conversion against the true
    #   0.769 / 0.590, Watson-Williams p .115 against .208.
    #
    #   A cluster that DOES straddle the wrap is torn in half, because 11.9 h and
    #   0.1 h -- 0.2 h apart on a 12 h circle -- land at opposite ends of the
    #   mapped half-circle. Measured: r-bar collapsed from 0.991 to 0.084, the
    #   circular mean moved from 0.03 h to 2.54 h, the angular difference was
    #   reported as 2.21 h instead of 0.23 h, and Watson-Williams went from
    #   p = .017 to p = .574. A real 0.2 h difference was reported as null on
    #   angles the procedure had itself scrambled.
    #
    # The divisor is read off the parameter name, which already carries the
    # harmonic, through the one helper both conversions now share.
    effective_period_param <- dance_effective_period(mod$period, param)

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
          # Hours to radians on THIS harmonic's effective period T/h (see above):
          # a full turn is T/h hours, not T.
          period <- effective_period_param
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

    # Store results, WITH the specification they were computed under (P20/R6).
    # A result that carries only its numbers has to be re-described by whatever
    # reads it, from inputs that may since have changed -- and the readout, the
    # plot and the export each did that separately, from `input$...` rather than
    # from the run. The spec below is written once, at the moment of the run,
    # and everything downstream reads it instead of the live controls.
    values$hp_pairwise_results <- results
    values$hp_pairwise_param <- param
    values$hp_pairwise_correction <- correction
    values$hp_pairwise_spec <- list(
      approach = "two_stage",
      parameter = param,
      circular = is_circular,
      period = mod$period,
      harmonic = if (is_circular) dance_param_harmonic(param) else NA_integer_,
      effective_period = if (is_circular) effective_period_param else NA_real_,
      group_var = input$harmonic_group_var,
      groups = groups,
      correction = correction,
      n_comparisons = nrow(results),
      family = sprintf("%d pairwise comparisons of %s across %d groups",
                       nrow(results), param, length(groups)),
      run_at = Sys.time())

    # AUDIT (P12.2): these three slots hold only the MOST RECENT comparison, so
    # a user who compares MESOR, then amplitude, then acrophase -- which is what
    # a chronobiology paper reports -- had the first two silently overwritten and
    # the publication report showed only the last. The panel above is a live view
    # of one comparison and stays as it is; the report needs the set, so every
    # run is also kept, keyed by parameter, and re-running a parameter replaces
    # its own entry rather than appending a duplicate.
    all_pw <- values$hp_pairwise_all
    if (is.null(all_pw)) all_pw <- list()
    all_pw[[param]] <- list(results = results, param = param,
                            correction = correction, run_at = Sys.time())
    values$hp_pairwise_all <- all_pw

    # AUDIT (P12.2). Bingham et al. (1982) note that a difference in AMPLITUDE
    # is not interpretable when the same groups also differ in ACROPHASE: the
    # amplitude is then being estimated about different phases. The publication
    # report states that caution, and a caution that is merely recited is worth
    # little -- so the acrophase verdict is recorded here, when the user runs
    # that comparison, and the report checks it instead of telling the reader to.
    # It is deliberately NOT computed on the fly from whatever is in memory: it
    # is the result of a comparison the user actually ran, or it is unknown.
    if (grepl("acro", param)) {
      values$hp_acrophase_differs <-
        any(is.finite(results$p_adjusted) & results$p_adjusted < 0.05)
      values$hp_acrophase_param <- param
    }

    showNotification("Pairwise comparisons completed!", type = "message", duration = 3)
  })

  # Display pairwise results
  output$hp_mixed_between_ui <- renderUI({
    ch <- dance_mixed_factor_choices(values)
    if (!length(ch)) return(helpText("No categorical covariate with 2-12 levels was found."))
    selectInput("hp_mixed_between", "Between-participant factor:", choices = ch)
  })
  output$hp_mixed_within_ui <- renderUI({
    ch <- dance_mixed_factor_choices(values)
    if (!length(ch)) return(NULL)
    selectInput("hp_mixed_within", "Within-participant factor (repeated):", choices = ch,
                selected = if (length(ch) > 1) ch[2] else ch[1])
  })

  output$hp_results <- renderPrint({
    if (identical(input$hp_approach, "mixed")) {
      if (is.null(values$mixed_results) ||
          !identical(values$mixed_results$kind, "cosinor")) {
        cat("Run the mixed cosinor to see results.\n"); return(invisible(NULL))
      }
      # rendered by the shared mixed readout, so the two entry points cannot
      # describe the same fit differently
      dance_mixed_readout(values$mixed_results); return(invisible(NULL))
    }
    if (identical(input$hp_approach, "population")) {
      res <- values$pop_cosinor
      if (is.null(res)) { cat("Run the population-mean cosinor to see results.\n"); return(invisible(NULL)) }
      f2 <- function(x) if (is.finite(x)) sprintf("%.2f", x) else "--"
      pf <- function(p) if (!is.finite(p)) "--" else if (p < .001) "< .001"
                        else sub("^0", "", sprintf("%.3f", p))
      cat("=== Population-mean cosinor (Bingham et al., 1982) ===\n\n")
      cat(sprintf("Grouping variable: %s   harmonic: %d   period: %s\n",
                  res$group_var, res$harmonic, format(res$period)))
      cat(sprintf("%d groups, %d participants; F tests on (%d, %d) df\n\n",
                  res$n_groups, res$n_subjects, res$df1, res$df2))
      g <- res$groups
      cat(sprintf("  %-14s %5s %10s %11s %12s\n", "group", "n", "MESOR", "amplitude", "acrophase"))
      for (i in seq_len(nrow(g)))
        cat(sprintf("  %-14s %5d %10s %11s %12s\n", g$group[i], g$n[i],
                    f2(g$mesor[i]), f2(g$amplitude[i]), f2(g$acrophase_time[i])))
      cat("\n  Acrophase is in time units on this harmonic's effective period.\n\n")
      t <- res$tests
      for (i in seq_len(nrow(t)))
        cat(sprintf("  %-10s F(%d, %d) = %8s   p %s\n", t$parameter[i],
                    t$df1[i], t$df2[i], f2(t$F[i]), pf(t$p[i])))
      cat("\nHow these differ from the pairwise tests\n")
      cat("---------------------------------------\n")
      cat("  The coefficient PAIRS are averaged as vectors, so amplitude is\n")
      cat("  tested against the variance along the pooled mean phase direction\n")
      cat("  and acrophase against the variance perpendicular to it. All groups\n")
      cat("  are compared at once, so there is no multiplicity correction to\n")
      cat("  apply and no pairwise family to control.\n")
      # The joint test on the (cosine, sine) vector. Amplitude and acrophase are
      # its marginals, and it is the only one of the three with no blind spot --
      # see the note in server/08b_helpers_popcosinor.R.
      jt <- res$joint
      if (isTRUE(jt$ok)) {
        cat(sprintf("\n  Joint test on the rhythmic vector (cosine, sine):\n"))
        cat(sprintf("    Wilks' Lambda = %s, F(%d, %d) = %s, p %s\n",
                    f2(jt$lambda), jt$df1, jt$df2, f2(jt$F), pf(jt$p)))
        cat("    This is the omnibus the amplitude and acrophase rows decompose;\n")
        cat("    it is the row to read when the two marginals disagree with it.\n")
      }

      if (!isTRUE(res$acrophase_test_supported)) {
        cat("\n  CAUTION -- the acrophase test is out of the range where it has power.\n")
        cat(strwrap(res$acrophase_test_note, width = 74, prefix = "  "), sep = "\n")
        cat("\n  The amplitude row is therefore NOT declared interpretable: Bingham's\n")
        cat("  condition is that the groups share a phase, and that has not been\n")
        cat("  established here -- it has only failed to be rejected by a test that\n")
        cat("  cannot see this difference.\n")
      } else if (isTRUE(res$acrophase_differs)) {
        cat("\n  WARNING -- the acrophases DIFFER (p ", pf(t$p[3]), ").\n", sep = "")
        cat("  Bingham et al. note that an amplitude difference cannot be\n")
        cat("  interpreted in that case: the amplitudes are being compared about\n")
        cat("  different phases. Do not read the amplitude row above as a\n")
        cat("  difference in rhythm strength.\n")
      } else {
        cat(sprintf("\n  The group acrophases span %s h on a %s h effective period, which is\n",
                    f2(res$max_angular_sep_time), f2(res$period / res$harmonic)))
        cat("  inside the quarter cycle where the perpendicular-displacement test is\n")
        cat("  monotone, and they do not differ detectably. The amplitude comparison\n")
        cat("  is interpretable on Bingham's own condition.\n")
      }
      cat("\nWhat this does not establish\n----------------------------\n")
      cat("  This is still a two-stage procedure in one respect: it starts from\n")
      cat("  per-participant coefficients and does not propagate each\n")
      cat("  participant's own estimation error. What it adds over the pairwise\n")
      cat("  route is that the pair is treated as one bivariate object and the\n")
      cat("  within-group covariance is pooled, not that the first stage has\n")
      cat("  gone away. The period was fixed, so everything is conditional on it.\n")
      return(invisible(NULL))
    }
    req(values$hp_pairwise_results)
    results <- values$hp_pairwise_results
    param <- values$hp_pairwise_param
    correction <- values$hp_pairwise_correction

    # Read the design off the stored specification, not off the live controls:
    # the controls can have moved since the run.
    spec <- values$hp_pairwise_spec
    is_circular <- if (is.null(spec)) grepl("acrophase_time", param, ignore.case = TRUE)
                   else isTRUE(spec$circular)

    cat("=== Pairwise Group Comparisons ===\n\n")
    cat("Parameter:", param, "\n")
    if(is_circular) {
      cat("Data type: Circular (using Watson-Williams test)\n")
      if (!is.null(spec) && is.finite(spec$effective_period))
        cat(sprintf("Angles are on harmonic %d's effective period, %.4g h (= %g / %d).\n",
                    spec$harmonic, spec$effective_period, spec$period, spec$harmonic))
    } else {
      cat("Data type: Linear (using Welch's t-test)\n")
    }
    cat("Correction method:", correction, "\n")
    # P20/R12: name the family the correction is applied over, so a reader knows
    # what "adjusted" was adjusted across.
    cat("Number of comparisons:", nrow(results), "\n")
    if (!is.null(spec)) cat("Multiplicity family:", spec$family, "\n")
    cat("\n")

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
      "mesor" = "Constant term (b0)",
      "mesor_adj" = "MESOR (rhythm-adjusted mean)",
      "value_at_start" = "Predicted value at the first observation",
      "amplitude_1" = "H1 Amplitude",
      "amplitude_2" = "H2 Amplitude",
      "amplitude_3" = "H3 Amplitude",
      # AUDIT: these columns hold MODEL-elapsed hours. A DIFFERENCE between two
      # groups is origin-free -- the shift cancels -- so the pairwise tests
      # themselves need no conversion. The label says elapsed so nobody reads a
      # group mean off this table as a time of day; clock acrophases are in the
      # results panel.
      "acrophase_time_1" = "H1 Acrophase (elapsed h; differences are origin-free)",
      "acrophase_time_2" = "H2 Acrophase (elapsed h; differences are origin-free)",
      "acrophase_time_3" = "H3 Acrophase (elapsed h; differences are origin-free)",
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
    .pal <- dance_group_colors(groups)
    plot_data <- list()
    for(g in groups) {
      plot_data[[as.character(g)]] <- params[[param]][params$group == g]
    }

    # Boxplot
    boxplot(plot_data,
            main = paste("Group Comparison:", param_label),
            ylab = param_label,
            xlab = "Group",
            # AUDIT: was rainbow(), which is a rainbow ramp used as a
            # categorical palette -- unordered hues that imply an order, and not
            # colourblind-separable. The shared palette is keyed by group NAME,
            # so a group is the same colour here as in every other figure.
            col = dance_group_fill(unname(.pal[groups]), 0.30),
            border = unname(.pal[groups]),
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
             col = dance_group_fill(unname(.pal[[as.character(g)]]), 0.55),
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
            "mesor" = "Constant term (b0)",
      "mesor_adj" = "MESOR (rhythm-adjusted mean)",
      "value_at_start" = "Predicted value at the first observation",
            "amplitude_1" = "H1 Amplitude",
            "amplitude_2" = "H2 Amplitude",
            "amplitude_3" = "H3 Amplitude",
            "acrophase_time_1" = "H1 Acrophase (elapsed h; differences are origin-free)",
            "acrophase_time_2" = "H2 Acrophase (elapsed h; differences are origin-free)",
            "acrophase_time_3" = "H3 Acrophase (elapsed h; differences are origin-free)",
            "r_squared" = "R²",
            "A_sat" = "A_sat",
            "tau" = "τ (tau)",
            "percent_S" = "Process S (%)",
            "percent_C" = "Process C (%)"
          )
          param_label <- ifelse(param %in% names(param_labels), param_labels[param], param)

          groups <- levels(params$group)
          n_groups <- length(groups)
          .pal <- dance_group_colors(groups)
          plot_data <- list()
          for(g in groups) {
            plot_data[[as.character(g)]] <- params[[param]][params$group == g]
          }

          par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
          boxplot(plot_data,
                  main = paste("Group Comparison:", param_label),
                  ylab = param_label,
                  xlab = "Group",
                  col = dance_group_fill(unname(.pal[groups]), 0.30),
                  border = unname(.pal[groups]),
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
                   col = dance_group_fill(unname(.pal[[as.character(g)]]), 0.55),
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
