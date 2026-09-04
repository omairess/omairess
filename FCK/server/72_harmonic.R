# ==========================================================================
# server/72_harmonic.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 2879-7190  (cosinor core, harmonic regression + outputs)
#
# CHANGELOG - 2026-09-03 cosinor audit (this file is NO LONGER verbatim)
# -----------------------------------------------------------------------------
# HARD BUGS
#  1.1 The pooled fitted equation dropped the homeostatic term. Root cause: the
#      pooled builder read pop$indiv_means$A_sat / $tau / $trend_linear /
#      $trend_log, and indiv_means was never given ANY trend parameter, so the
#      branch was dead for every trend type. Both the pooled and the group
#      equations now go through fck_format_equation(); the duplicate builder is
#      deleted.
#  1.2 The Rayleigh test used the amplitude-weighted resultant (0.824 -> Z=886)
#      where the unweighted one was required (0.789 -> Z=812). fck_resultants()
#      returns both under unambiguous names; only the unweighted one reaches
#      fck_rayleigh(). The amplitude-weighted vector mean is unchanged.
#  1.3 R2_S and R2_C were overlapping marginals summing to 124.7%. Replaced by
#      commonality analysis (unique_S + unique_C + shared == total). The
#      auto-generated dominance verdict is deleted.
#  1.4 "MESOR" was the fitted constant. Renamed to Intercept (beta_0) here, in
#      the plots, in the parameter table, in the CSV export and in the pairwise
#      comparisons. A genuine rhythm-adjusted mean is computed by integration
#      and reported as the MESOR. time_origin added.
#  1.5 Group n's summed to 1304 of 1305: unique() kept NA as a level and
#      which(x == NA) is empty, so the subject fell through the n>=3 guard in
#      silence. UNASSIGNED is now a real row and the totals are asserted.
#  1.6 Every averaged quantity now names its estimator; a linear SD is never
#      printed beside a vector- or circular-averaged value.
#  1.7 fmt2() (round half away from zero) throughout; the DV is named with units
#      and bounds; the "(units)" placeholder is interpolated; the H2 modulo
#      convention is printed next to every H2 acrophase.
#
# STATISTICAL
#  2.1 data_source = raw | smoothed is now a user choice, with the inflation it
#      causes stated in the report rather than left implicit.
#  2.2 Parameter correlation matrix, design condition number at each origin, and
#      free-tau vs fixed-tau Delta-AIC, so the A_sat/tau ridge is documented.
#  2.3 Convergence was never inspected: nls(warnOnly=TRUE) RETURNS a fit at the
#      iteration limit, so every non-converged optimisation was counted as a
#      success. Now captured per subject; converged / boundary / failed are
#      reported separately and non-converged fits are excluded from the
#      population summaries.
#  2.4 Delta-AICc across a nested model set with Akaike weights, replacing mean
#      AIC/AICc/BIC (which are constant offsets of one another, hence the
#      identical SDs).
#  2.5 Delta-method amplitude AND acrophase SEs in the nonlinear path, plus
#      Bingham elliptical joint confidence regions.
#  2.6 Population-mean cosinor with group x harmonic terms as the primary
#      analysis; effect sizes and CIs on every group contrast; a monotone-trend
#      contrast for the ordered age bands; the Watson-Williams concentration
#      assumption checked in the output.
#
# FOUND DURING THE AUDIT, NOT IN THE BRIEF
#   a. The zero-amplitude F test put the WHOLE model's sum of squares over the
#      harmonics' df alone, crediting Process S's variance to the rhythm. Now a
#      proper full-vs-trend-only test. This, not only the smoothing, is why
#      95.3% of subjects came out "significantly rhythmic".
#   b. "LOOCV RMSE" in the nonlinear path was the in-sample residual RMSE. Now
#      genuine leave-one-out refits, with the label following the computation.
#   c. The nonlinear amplitude SE was sqrt(se_c^2+se_s^2)/sqrt(2) -- covariance
#      discarded -- and the acrophase SE was NA unconditionally.
#   d. acrophases_time omitted the /h divisor, so acrophase_time_2 was on a 0-24
#      scale while the group summaries used 0-12; the H2 group tests ran on the
#      wrong scale.
#   e. In the nonlinear path R2_S reused the full model's coefficients with the
#      cosines deleted rather than refitting the trend-only model.
#
# Tests: tests/testthat/test-cosinor-audit.R (and tests/audit_test.R, which runs
# the same file without testthat). Old-vs-new report: tests/report_harness.R.
# ==========================================================================
  # ==============================================================================
  # HARMONIC REGRESSION (COSINOR ANALYSIS) MODULE
  # ==============================================================================
  
  # Variable selection UI for harmonic regression
  output$harmonic_var_select_ui <- renderUI({
    req(values$data)  # Only require data; covariates are optional

    # Time variable options from covariates (if available)
    numeric_vars <- if(!is.null(values$covariates)) {
      names(values$covariates)[sapply(values$covariates, is.numeric)]
    } else {
      character(0)  # Empty vector if no covariates
    }
    n_time <- ncol(values$data)
    col_names <- colnames(values$data)
    
    # Try to extract time values from column names
    # Looks for patterns like: VAS_10, T10, time10, col_10, 10:00, 11PM, 2AM, etc.
    suggested_times <- NULL
    detected_pattern <- NULL
    
    if(!is.null(col_names) && length(col_names) > 0) {
      
      # Pattern 0: AM/PM format (e.g., 11PM, 2AM, 11:30PM, VAS_11PM)
      # Check if any column contains AM or PM
      if(any(grepl("[0-9]\\s*[AaPp][Mm]", col_names))) {
        if(all(grepl("[0-9]\\s*[AaPp][Mm]", col_names))) {
          suggested_times <- sapply(col_names, function(cn) {
            # Extract hour, optional minutes, and AM/PM
            # First try with minutes (e.g., 11:30PM)
            if(grepl("[0-9]{1,2}:[0-9]{2}\\s*[AaPp][Mm]", cn)) {
              hour <- as.numeric(gsub(".*?([0-9]{1,2}):[0-9]{2}\\s*[AaPp][Mm].*", "\\1", cn))
              mins <- as.numeric(gsub(".*?[0-9]{1,2}:([0-9]{2})\\s*[AaPp][Mm].*", "\\1", cn))
            } else {
              # Without minutes (e.g., 11PM)
              hour <- as.numeric(gsub(".*?([0-9]{1,2})\\s*[AaPp][Mm].*", "\\1", cn))
              mins <- 0
            }
            ampm <- toupper(gsub(".*([AaPp][Mm]).*", "\\1", cn))
            
            # Convert to 24-hour: 12AM=0, 1-11AM=1-11, 12PM=12, 1-11PM=13-23
            hour_24 <- if(ampm == "AM") {
              if(hour == 12) 0 else hour
            } else {
              if(hour == 12) 12 else hour + 12
            }
            hour_24 + mins / 60
          }, USE.NAMES = FALSE)
          detected_pattern <- "AM/PM format"
        }
      }
      
      # Pattern 1: Dutch/European hour format with "u" suffix (e.g., KSS_9u_dag1, var_14u_something)
      # The number before "u" or "u_" is the hour
      if(is.null(suggested_times)) {
        if(all(grepl("_[0-9]{1,2}u", col_names))) {
          suggested_times <- as.numeric(gsub(".*_([0-9]{1,2})u.*", "\\1", col_names))
          detected_pattern <- "hour with 'u' suffix (e.g., 9u = 9:00)"
        }
      }
      
      # Pattern 2: Trailing numbers (e.g., VAS_10, VAS_12, T10, col10)
      # BUT skip if trailing number looks like day indicator (dag1, dag2, day1, day2)
      if(is.null(suggested_times)) {
        # Check if trailing numbers are likely day indicators
        trailing_nums <- gsub(".*[^0-9]([0-9]+)$", "\\1", col_names)
        is_day_indicator <- all(grepl("(dag|day)[0-9]+$", col_names, ignore.case = TRUE))
        
        if(!is_day_indicator && all(grepl("^[0-9]+$", trailing_nums))) {
          suggested_times <- as.numeric(trailing_nums)
          detected_pattern <- "trailing numbers"
        }
      }
      
      # Pattern 3: Numbers after underscore (e.g., var_10, var_12)
      if(is.null(suggested_times)) {
        underscore_nums <- gsub(".*_([0-9]+).*", "\\1", col_names)
        if(all(grepl("^[0-9]+$", underscore_nums)) && !all(underscore_nums == col_names)) {
          suggested_times <- as.numeric(underscore_nums)
          detected_pattern <- "underscore pattern"
        }
      }
      
      # Pattern 4: Time format HH:MM or HH (e.g., 10:00, 12:30) - 24h format
      if(is.null(suggested_times)) {
        time_match <- grepl("([0-9]{1,2}):?([0-9]{0,2})", col_names)
        if(all(time_match)) {
          hours <- as.numeric(gsub(".*?([0-9]{1,2}):?([0-9]{0,2}).*", "\\1", col_names))
          mins <- gsub(".*?([0-9]{1,2}):?([0-9]{0,2}).*", "\\2", col_names)
          mins <- ifelse(mins == "", 0, as.numeric(mins))
          suggested_times <- hours + mins / 60
          detected_pattern <- "time format"
        }
      }

    }
    
    # Build suggestion text
    suggestion_text <- "e.g., 8,9,10,11,12,14,16,18,20,21,22,23,0,2,4,6"
    suggestion_value <- ""
    detection_msg <- NULL
    
    if(!is.null(suggested_times) && length(suggested_times) == n_time) {
      suggestion_value <- paste(suggested_times, collapse = ",")
      suggestion_text <- suggestion_value
      times_preview <- paste(head(suggested_times, 6), collapse=", ")
      if(n_time > 6) times_preview <- paste0(times_preview, ", ...")
      detection_msg <- div(style = "color: green; font-size: 0.9em;",
                           icon("check-circle"),
                           sprintf(" Detected %d time values from column names: %s", 
                                   n_time, times_preview))
    }
    
    tagList(
      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              # MERGED APP: reuse the clock times parsed once
                              # at import (values$time_clock) instead of
                              # re-detecting them here.  Additive: the default
                              # is still "_index_".
                              "Use shared clock times parsed at import" = "_shared_",
                              numeric_vars),
                  selected = "_index_"),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_shared_'",
        if(!is.null(values$time_clock) && length(values$time_clock) == n_time) {
          div(style = "color: green; font-size: 0.9em;", icon("check-circle"),
              sprintf(" Using the %d clock times parsed at import: %s%s",
                      n_time, paste(head(values$time_clock, 6), collapse = ", "),
                      if(n_time > 6) ", ..." else ""))
        } else {
          div(style = "color: #b00; font-size: 0.9em;", icon("exclamation-triangle"),
              " Import could not parse clock times from the column names. Use 'Specify times manually'.")
        }
      ),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_manual_'",
        if(!is.null(detection_msg)) detection_msg,
        textAreaInput("harmonic_manual_times", 
                      paste0("Enter ", n_time, " time values (comma-separated):"),
                      value = suggestion_value,
                      placeholder = suggestion_text,
                      rows = 2),
        helpText("Enter the actual clock times for each column in your data. Use 24-hour format or decimal hours.")
      ),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_index_'",
        if(!is.null(suggested_times) && length(suggested_times) == n_time) {
          helpText(HTML(paste0("<b>Note:</b> Detected time values in column names (", 
                               paste(head(suggested_times, 4), collapse=", "), 
                               if(n_time > 4) ", ..." else "",
                               "). Consider using 'Specify times manually' if spacing is unequal.")))
        } else {
          helpText(HTML("<b>Warning:</b> This assumes measurements are equally spaced across the period. If your measurements are unequally spaced (e.g., hourly during day, 2-hourly at night), use 'Specify times manually' instead."))
        }
      )
    )
  })
  
  # Group variable UI
  output$harmonic_group_var_ui <- renderUI({
    req(values$covariates)
    cat_vars <- names(values$covariates)[sapply(values$covariates, function(x) {
      is.factor(x) || is.character(x) || length(unique(x)) <= 10
    })]
    selectInput("harmonic_group_var", "Group Variable (optional):",
                choices = c("None" = "_none_", cat_vars))
  })

  # Parameter bounds hints based on data
  output$harmonic_bounds_hints <- renderUI({
    req(values$data)

    # Get data (use smoothed if available)
    Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

    # Calculate data statistics
    y_min <- min(Y, na.rm = TRUE)
    y_max <- max(Y, na.rm = TRUE)
    y_range <- y_max - y_min
    y_mean <- mean(Y, na.rm = TRUE)

    # Get time information
    n_time <- ncol(Y)
    time_max <- if(!is.null(input$harmonic_time_var) && input$harmonic_time_var == "_index_") {
      input$harmonic_period
    } else {
      n_time  # Conservative estimate
    }

    # Build hints text
    hints_html <- sprintf(
      "<div style='background-color: #e8f4f8; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
       <strong>📊 Data Range Hints:</strong><br>
       <small>
       <strong>Your data:</strong> Min=%.2f, Max=%.2f, Mean=%.2f, Range=%.2f<br>
       <strong>Suggested intercept (\u03b2\u2080) bounds:</strong> [%.2f, %.2f] (mean \u00b1 range)<br>
       <strong>Suggested Amplitude max:</strong> %.2f (observed range)<br>",
      y_min, y_max, y_mean, y_range,
      y_mean - y_range, y_mean + y_range,
      y_range
    )

    # Add exp_sat specific hints if that trend type is selected
    if(!is.null(input$harmonic_trend_type) && input$harmonic_trend_type == "exp_sat") {
      hints_html <- paste0(hints_html, sprintf(
        "<strong>Suggested A_sat bounds:</strong> [%.2f, %.2f] (0.5× to 2× range)<br>
         <strong>Suggested τ bounds:</strong> [0.5, %.1f] (0.5 to max time)<br>",
        y_range * 0.5, y_range * 2,
        time_max
      ))
    }

    hints_html <- paste0(hints_html, "</small></div>")

    HTML(hints_html)
  })

  # Warning UI for harmonic count vs data points
  output$harmonic_warning_ui <- renderUI({
    req(values$data)
    n_time <- ncol(values$data)
    n_harmonics <- input$n_harmonics
    
    # Need at least 2*n_harmonics + 1 parameters (2 per harmonic + MESOR)
    min_required <- 2 * n_harmonics + 2  # +2 for some df for error
    max_safe_harmonics <- floor((n_time - 2) / 2)
    
    period <- if(!is.null(input$harmonic_period)) input$harmonic_period else 24
    
    # Build harmonic info table
    harmonic_info <- paste0(
      "<small><b>Harmonic periods:</b> ",
      paste(sapply(1:n_harmonics, function(h) paste0("H", h, "=", round(period/h, 1), "h")), collapse=", "),
      "</small>"
    )
    
    if(n_harmonics > max_safe_harmonics) {
      tagList(
        div(style = "color: red; font-weight: bold;",
            icon("exclamation-triangle"),
            sprintf(" Warning: %d harmonics require at least %d time points. You have %d.", 
                    n_harmonics, min_required, n_time)),
        div(style = "color: orange;",
            sprintf("Maximum safe harmonics for your data: %d", max_safe_harmonics)),
        HTML(harmonic_info)
      )
    } else if(n_harmonics > max_safe_harmonics - 1) {
      tagList(
        div(style = "color: orange;",
            icon("exclamation-circle"),
            " Approaching maximum harmonics for your data. Model may overfit."),
        HTML(harmonic_info)
      )
    } else {
      HTML(harmonic_info)
    }
  })
  
  # Subject selector for individual plots
  output$harmonic_subject_selector <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(!is.null(mod$individual_fits)) {
      n_subj <- length(mod$individual_fits)
      
      # Check which fits succeeded and show details for failed ones
      subject_labels <- sapply(1:n_subj, function(i) {
        fit_i <- mod$individual_fits[[i]]
        if(!is.null(fit_i) && isTRUE(fit_i$success)) {
          paste("Subject", i)
        } else if(!is.null(fit_i) && !is.null(fit_i$n_valid)) {
          paste0("Subject ", i, " (failed: ", fit_i$n_valid, "/", fit_i$n_required, " pts)")
        } else {
          paste("Subject", i, "(failed)")
        }
      })
      
      selectInput("harmonic_subject_select", "Select Subject:", 
                  choices = c("All (overlay)" = "all", 
                              "Mean curve" = "mean",
                              setNames(1:n_subj, subject_labels)),
                  selected = "mean")
    }
  })
  
  # Harmonic selector for polar plot
  output$harmonic_selector_polar <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      selectInput("selected_harmonic_polar", "Display Harmonic:", 
                  choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                  selected = 1)
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # Harmonic selector for parameter distribution
  output$harmonic_selector_dist <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      fluidRow(
        column(4,
               selectInput("selected_harmonic_dist", "Display Harmonic:", 
                           choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                           selected = 1)
        ),
        column(8,
               helpText("Select which harmonic to display in the amplitude and acrophase distributions.")
        )
      )
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # Harmonic selector for group comparison
  output$harmonic_selector_group <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      fluidRow(
        column(4,
               selectInput("selected_harmonic_group", "Compare Harmonic:", 
                           choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                           selected = 1)
        ),
        column(8,
               helpText("Select which harmonic to use for group comparisons.")
        )
      )
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # ==============================================================================
  # CIRCULAR STATISTICS HELPER FUNCTIONS
  # ==============================================================================
  # Custom implementations for acrophase analysis (no external package required)
  # Based on Mardia & Jupp (2000) "Directional Statistics"
  
  # Circular mean (returns radians)
  circular_mean <- function(angles_rad) {
    x <- mean(cos(angles_rad), na.rm = TRUE)
    y <- mean(sin(angles_rad), na.rm = TRUE)
    atan2(y, x)
  }
  
  # Mean resultant length (measure of concentration, 0-1)
  mean_resultant_length <- function(angles_rad) {
    x <- mean(cos(angles_rad), na.rm = TRUE)
    y <- mean(sin(angles_rad), na.rm = TRUE)
    sqrt(x^2 + y^2)
  }
  
  # Circular standard deviation (in same units as input)
  circular_sd <- function(angles_rad) {
    r_bar <- mean_resultant_length(angles_rad)
    if(is.na(r_bar)) {
      return(NA)
    }
    if(r_bar > 0 && r_bar < 1) {
      sqrt(-2 * log(r_bar))
    } else if(r_bar >= 1) {
      0  # All points identical
    } else {
      NA  # Undefined
    }
  }
  
  # Circular standard error (approximate, based on von Mises)
  circular_se <- function(angles_rad) {
    n <- sum(!is.na(angles_rad))
    r_bar <- mean_resultant_length(angles_rad)
    if(is.na(r_bar)) {
      return(NA)
    }
    if(r_bar > 0 && n > 1) {
      # Approximate SE for circular mean (Mardia & Jupp, 2000)
      1 / sqrt(n * r_bar^2)
    } else {
      NA
    }
  }
  
  # Watson-Williams test for comparing two or more groups of circular data
  watson_williams_test <- function(angles_list) {
    # angles_list: list of vectors, each containing angles in radians for one group
    k <- length(angles_list)
    if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))
    
    n <- sapply(angles_list, function(x) sum(!is.na(x)))
    N <- sum(n)
    
    if(any(n < 2)) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Each group needs at least 2 observations"))
    
    # Resultant lengths for each group
    R <- sapply(angles_list, function(x) {
      x <- x[!is.na(x)]
      sqrt(sum(cos(x))^2 + sum(sin(x))^2)
    })
    
    # Total resultant length (pooled)
    all_angles <- unlist(angles_list)
    all_angles <- all_angles[!is.na(all_angles)]
    R_total <- sqrt(sum(cos(all_angles))^2 + sum(sin(all_angles))^2)
    
    r_bar_total <- R_total / N
    
    if(r_bar_total < 0.45) {
      return(list(F = NA, df1 = k - 1, df2 = N - k, p = NA, r_bar = r_bar_total,
                  message = "Warning: Data too dispersed (r̄ < 0.45). Consider non-parametric test."))
    }
    
    # Concentration parameter estimate
    kappa <- if(r_bar_total < 0.53) {
      2 * r_bar_total + r_bar_total^3 + 5 * r_bar_total^5 / 6
    } else if(r_bar_total < 0.85) {
      -0.4 + 1.39 * r_bar_total + 0.43 / (1 - r_bar_total)
    } else {
      1 / (r_bar_total^3 - 4 * r_bar_total^2 + 3 * r_bar_total)
    }
    
    g <- 1 - 1 / (3 * 8 * kappa^2)
    sum_R <- sum(R)
    F_stat <- g * (N - k) * (sum_R - R_total) / ((k - 1) * (N - sum_R))
    
    df1 <- k - 1
    df2 <- N - k
    p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
    
    list(F = F_stat, df1 = df1, df2 = df2, p = p_value, kappa = kappa, r_bar = r_bar_total, message = NULL)
  }

  # Hotelling's T² test on (beta_cos, beta_sin) pairs - amplitude-weighted acrophase comparison
  # Tests whether the bivariate rhythmic vector (beta_cos, beta_sin) differs between groups.
  # Because amplitude = sqrt(beta_cos² + beta_sin²), subjects with stronger rhythms carry
  # more weight. For k>2 groups, a one-way MANOVA F approximation is used.
  hotelling_t2 <- function(beta_cos_list, beta_sin_list) {
    k <- length(beta_cos_list)
    if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))

    # Build per-group matrices
    mats <- lapply(seq_len(k), function(i) {
      x <- beta_cos_list[[i]]
      y <- beta_sin_list[[i]]
      ok <- complete.cases(x, y)
      cbind(x[ok], y[ok])
    })

    ns <- sapply(mats, nrow)
    if(any(ns < 3)) return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                                message = "Each group needs at least 3 observations"))
    N <- sum(ns)
    p <- 2  # two variables: beta_cos and beta_sin

    means <- lapply(mats, colMeans)

    if(k == 2) {
      # Two-sample Hotelling's T²
      S1 <- cov(mats[[1]]); S2 <- cov(mats[[2]])
      Sp <- ((ns[1]-1)*S1 + (ns[2]-1)*S2) / (N - 2)
      if(det(Sp) < .Machine$double.eps)
        return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                    message = "Singular pooled covariance - check data"))
      d  <- means[[1]] - means[[2]]
      T2 <- (ns[1]*ns[2]) / N * t(d) %*% solve(Sp) %*% d
      F_stat <- (N - p - 1) / ((N - 2) * p) * as.numeric(T2)
      df1 <- p; df2 <- N - p - 1
    } else {
      # k-group one-way MANOVA (Wilks' lambda → F approximation)
      grand_mean <- colMeans(do.call(rbind, mats))
      # Between-group sum-of-squares-and-products (H)
      H <- Reduce("+", lapply(seq_len(k), function(i)
        ns[i] * outer(means[[i]] - grand_mean, means[[i]] - grand_mean)))
      # Within-group (E)
      E <- Reduce("+", lapply(mats, function(m) {
        cm <- colMeans(m)
        t(sweep(m, 2, cm)) %*% sweep(m, 2, cm)
      }))
      if(det(E) < .Machine$double.eps)
        return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                    message = "Singular within-group covariance - check data"))
      lambda <- det(E) / det(H + E)
      # Rao's F approximation for p=2
      df1 <- p * (k - 1)
      df2 <- N - k - p + 1
      F_stat <- ((1 - sqrt(lambda)) / sqrt(lambda)) * (df2 / df1)
    }

    p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
    list(F = as.numeric(F_stat), df1 = df1, df2 = df2, p = as.numeric(p_value), message = NULL)
  }

  # ==============================================================================
  # CORE COSINOR FITTING FUNCTIONS
  # ==============================================================================
  
  # Single cosinor fit for one subject
  fit_cosinor <- function(time, y, period = 24, n_harmonics = 1, trend_type = "none",
                          use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                          amplitude_min = 0, amplitude_max = NA,
                          A_sat_min = NA, A_sat_max = NA,
                          tau_min = 0.5, tau_max = NA) {
    # Remove NAs
    valid <- complete.cases(time, y)
    time <- time[valid]
    y <- y[valid]
    n <- length(y)

    # Set default bounds based on data if not specified
    if(use_bounds) {
      y_min <- min(y, na.rm = TRUE)
      y_max <- max(y, na.rm = TRUE)
      y_range <- y_max - y_min
      t_max <- max(time) - min(time)

      if(is.na(mesor_min)) mesor_min <- y_min - y_range
      if(is.na(mesor_max)) mesor_max <- y_max + y_range
      if(is.na(amplitude_max)) amplitude_max <- y_range * 2
      if(is.na(A_sat_max)) A_sat_max <- y_range * 2
      if(is.na(tau_max)) tau_max <- t_max * 5
    }

    # Calculate number of trend parameters
    n_trend_params <- switch(trend_type,
                             "none" = 0,
                             "linear" = 1,
                             "log" = 1,
                             "exp_sat" = 2,  # A and tau for nonlinear fit
                             0)

    min_params <- 2 * n_harmonics + 1 + n_trend_params
    if(n < min_params + 1) {
      return(list(success = FALSE, message = "Insufficient data points"))
    }

    # For exponential saturation, use nonlinear fitting
    if(trend_type == "exp_sat") {
      return(fit_cosinor_nonlinear(time, y, period, n_harmonics, trend_type,
                                    FALSE, NULL, 0.32, 0.66,
                                    use_bounds, mesor_min, mesor_max, amplitude_min, amplitude_max,
                                    A_sat_min, A_sat_max, tau_min, tau_max))
    }

    # For linear/log models with bounds, use nonlinear least squares
    if(use_bounds && trend_type %in% c("linear", "log", "none")) {
      return(fit_cosinor_nonlinear(time, y, period, n_harmonics, trend_type,
                                    FALSE, NULL, 0.32, 0.66,
                                    use_bounds, mesor_min, mesor_max, amplitude_min, amplitude_max,
                                    A_sat_min, A_sat_max, tau_min, tau_max))
    }
    
    # Build design matrix with multiple harmonics (linear models)
    X <- matrix(1, nrow = n, ncol = 1)  # Intercept (MESOR)
    colnames_X <- "MESOR"
    trend_cols <- character(0)
    
    # Add trend based on type
    if(trend_type == "linear") {
      X <- cbind(X, time)
      colnames_X <- c(colnames_X, "trend_linear")
      trend_cols <- "trend_linear"
    } else if(trend_type == "log") {
      # Use log(t+1) to avoid log(0) and handle t=0
      t_offset <- min(time)
      log_time <- log(time - t_offset + 1)
      X <- cbind(X, log_time)
      colnames_X <- c(colnames_X, "trend_log")
      trend_cols <- "trend_log"
    }
    
    coef_offset <- 1 + length(trend_cols)
    
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      X <- cbind(X, cos(omega * time), sin(omega * time))
      colnames_X <- c(colnames_X, paste0("cos", h), paste0("sin", h))
    }
    colnames(X) <- colnames_X
    
    # Fit linear model
    fit <- lm(y ~ X - 1)  # -1 because X already has intercept
    coefs <- coef(fit)
    se <- summary(fit)$coefficients[, 2]
    
    # Extract parameters
    mesor <- coefs[1]
    mesor_se <- se[1]
    
    # Extract trend parameters
    trend_params <- list()
    if(trend_type != "none") {
      for(i in seq_along(trend_cols)) {
        trend_params[[trend_cols[i]]] <- list(
          coef = coefs[1 + i],
          se = se[1 + i]
        )
      }
    }
    
    # Calculate amplitude and acrophase for each harmonic
    amplitudes <- numeric(n_harmonics)
    acrophases <- numeric(n_harmonics)
    amp_se <- numeric(n_harmonics)
    acro_se <- numeric(n_harmonics)
    
    vcov_mat <- vcov(fit)
    
    for(h in 1:n_harmonics) {
      cos_idx <- coef_offset + 2 * (h - 1) + 1
      sin_idx <- coef_offset + 2 * (h - 1) + 2
      
      beta_cos <- coefs[cos_idx]
      beta_sin <- coefs[sin_idx]
      
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi
      
      if(amplitudes[h] > 1e-10) {
        grad_amp <- c(beta_cos, beta_sin) / amplitudes[h]
        idx <- c(cos_idx, sin_idx)
        var_amp <- t(grad_amp) %*% vcov_mat[idx, idx] %*% grad_amp
        amp_se[h] <- sqrt(var_amp)
        
        grad_acro <- c(-beta_sin, beta_cos) / (amplitudes[h]^2)
        var_acro <- t(grad_acro) %*% vcov_mat[idx, idx] %*% grad_acro
        acro_se[h] <- sqrt(var_acro)
      } else {
        amp_se[h] <- NA
        acro_se[h] <- NA
      }
    }
    
    # AUDIT (extra i): harmonic h completes h cycles per period, so its acrophase
    # occupies the EFFECTIVE period T/h and the conversion to hours must divide
    # by h. Both fitters omitted the divisor, so acrophase_time_2 was written on
    # a 0-24 scale while the vector-averaged group value used 0-12. The report
    # printed the two side by side as if comparable, and the group comparisons
    # (output$harmonic_group_test_results) tested H2 on the wrong scale.
    acrophases_time <- vapply(seq_len(n_harmonics),
                              function(h) phi_to_hours(acrophases[h], period, h),
                              numeric(1))
    acro_se_time <- vapply(seq_len(n_harmonics),
                           function(h) phi_to_hours(acro_se[h], period, h),
                           numeric(1))

    # Goodness of fit - calculate manually because lm(y ~ X - 1) R² is vs origin, not mean
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum(residuals(fit)^2)
    r_squared <- 1 - ss_resid / ss_total

    # Adjusted R² accounting for number of predictors
    n_predictors <- ncol(X)  # MESOR + trend + harmonics
    adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - n_predictors)

    percent_rhythm <- r_squared * 100

    # AUDIT (extra a): the zero-amplitude test is "all harmonic coefficients are
    # zero GIVEN the trend", so the numerator SS is (trend-only residual SS minus
    # full residual SS). The old code used (ss_total - ss_resid) -- the WHOLE
    # model's SS, homeostatic trend included -- over the harmonics' df alone,
    # charging the rhythm for every bit of variance Process S explains. With a
    # saturating trend worth ~28% of variance on its own this is a large upward
    # bias, and it is a bigger contributor to the 95.3% "significant rhythm" rate
    # than the FDA smoothing is.
    X_trend_only <- X[, seq_len(coef_offset), drop = FALSE]
    ss_resid_trend_only <- tryCatch({
      sum(residuals(lm(y ~ X_trend_only - 1))^2)
    }, error = function(e) ss_total)
    .zt <- fck_zero_amplitude_test(ss_resid, ss_resid_trend_only,
                                   n, n_predictors, n_harmonics)
    f_stat  <- .zt$F
    p_value <- .zt$p

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian linear model
    sigma_sq <- ss_resid / n
    log_lik <- -n/2 * (log(2*pi) + log(sigma_sq) + 1)

    # AIC: Akaike Information Criterion
    # AIC = -2*log(L) + 2*k, where k = number of parameters
    k <- n_predictors + 1  # predictors + sigma
    aic <- -2 * log_lik + 2 * k

    # AICc: Corrected AIC for small samples
    # AICc = AIC + 2*k*(k+1)/(n-k-1)
    aicc <- if(n - k - 1 > 0) {
      aic + (2 * k * (k + 1)) / (n - k - 1)
    } else {
      NA  # Not defined when n-k-1 <= 0
    }

    # BIC: Bayesian Information Criterion
    # BIC = -2*log(L) + k*log(n)
    bic <- -2 * log_lik + k * log(n)

    # LOOCV: Leave-one-out cross-validation (leave out time points)
    # For each observation, refit the model without it and predict
    loocv_errors <- numeric(n)
    for(i in 1:n) {
      # Remove observation i
      X_loo <- X[-i, , drop = FALSE]
      y_loo <- y[-i]

      # Refit model
      fit_loo <- tryCatch({
        lm.fit(X_loo, y_loo)
      }, error = function(e) NULL)

      if(!is.null(fit_loo)) {
        # Predict left-out observation
        y_pred <- sum(X[i, ] * fit_loo$coefficients)
        loocv_errors[i] <- (y[i] - y_pred)^2
      } else {
        loocv_errors[i] <- NA
      }
    }

    # LOOCV RMSE (root mean squared error)
    loocv_rmse <- sqrt(mean(loocv_errors, na.rm = TRUE))

    # ===========================================================================
    # Variance decomposition: Calculate R² for Process S and Process C separately
    # ===========================================================================
    r_squared_S <- 0  # Variance explained by homeostatic trend alone
    r_squared_C <- 0  # Variance explained by circadian rhythm alone
    percent_S <- 0    # Percentage of total R² from Process S
    percent_C <- 0    # Percentage of total R² from Process C

    # Baseline model (MESOR only)
    y_mean <- mean(y)
    ss_resid_baseline <- sum((y - y_mean)^2)

    # Model with trend only (Process S)
    if(trend_type != "none") {
      if(trend_type == "exp_sat") {
        # For exponential saturation, use the already fitted trend component
        # Extract fitted trend values from the full model
        if(!is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
          A_sat_val <- trend_params$A_sat$coef
          tau_val <- trend_params$tau$coef
          t_offset_temp <- min(time)
          fitted_trend <- mesor + A_sat_val * (1 - exp(-(time - t_offset_temp) / tau_val))
          ss_resid_trend <- sum((y - fitted_trend)^2)
          r_squared_S <- 1 - ss_resid_trend / ss_total
        }
      } else {
        X_trend <- matrix(1, nrow = n, ncol = 1)
        if(trend_type == "linear") {
          X_trend <- cbind(X_trend, time)
        } else if(trend_type == "log") {
          t_offset_temp <- min(time)
          log_time_temp <- log(time - t_offset_temp + 1)
          X_trend <- cbind(X_trend, log_time_temp)
        }
        fit_trend <- lm(y ~ X_trend - 1)
        ss_resid_trend <- sum(residuals(fit_trend)^2)
        r_squared_S <- 1 - ss_resid_trend / ss_total
      }
    }

    # Model with circadian only (Process C)
    X_circ <- matrix(1, nrow = n, ncol = 1)
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      X_circ <- cbind(X_circ, cos(omega * time), sin(omega * time))
    }
    fit_circ <- lm(y ~ X_circ - 1)
    ss_resid_circ <- sum(residuals(fit_circ)^2)
    r_squared_C <- 1 - ss_resid_circ / ss_total

    # AUDIT 1.3: r_squared_S and r_squared_C are MARGINAL R2s from two
    # overlapping, collinear predictor blocks. They do not partition anything --
    # in the reported output they summed to 1.121 against a total of 0.892, and
    # their "proportions" summed to 124.7%. Commonality analysis (Chevan &
    # Sutherland 1991) gives a partition that sums to the total exactly, and
    # makes the shared component (0.229 in that output) visible instead of
    # double-counting it.
    #
    # The marginal R2s are still returned, because they are what the commonality
    # decomposition is computed FROM and dropping them would make the arithmetic
    # unauditable. They are no longer presented as a decomposition.
    .cm <- fck_commonality(r_squared, r_squared_S, r_squared_C)
    unique_S <- .cm$unique_S
    unique_C <- .cm$unique_C
    shared_SC <- .cm$shared
    .pc <- fck_commonality_pct(.cm)
    percent_S <- .pc$unique_S       # now: % of total R2 UNIQUE to S
    percent_C <- .pc$unique_C       # now: % of total R2 UNIQUE to C
    percent_shared <- .pc$shared
    
    # Store time offset for prediction (needed for log/sqrt)
    t_offset <- if(trend_type == "log") min(time) else 0
    t_center <- 0
    
    list(
      success = TRUE,
      mesor = mesor,
      mesor_se = mesor_se,
      trend_type = trend_type,
      trend_params = trend_params,
      t_offset = t_offset,
      t_center = t_center,
      amplitudes = amplitudes,
      amp_se = amp_se,
      acrophases = acrophases,
      acrophases_time = acrophases_time,
      acro_se = acro_se,
      acro_se_time = acro_se_time,
      coefs = coefs,
      se = se,
      vcov = vcov_mat,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      percent_rhythm = percent_rhythm,
      f_stat = f_stat,
      p_value = p_value,
      aic = aic,                        # Akaike Information Criterion
      aicc = aicc,                      # Corrected AIC for small samples
      bic = bic,                        # Bayesian Information Criterion
      loocv_rmse = loocv_rmse,          # Leave-one-out CV RMSE
      r_squared_S = r_squared_S,        # MARGINAL R² of trend-only model
      r_squared_C = r_squared_C,        # MARGINAL R² of harmonics-only model
      unique_S = unique_S,              # commonality: unique to Process S
      unique_C = unique_C,              # commonality: unique to Process C
      shared_SC = shared_SC,            # commonality: shared (may be negative)
      percent_S = percent_S,            # % of total R² UNIQUE to S
      percent_C = percent_C,            # % of total R² UNIQUE to C
      percent_shared = percent_shared,  # % of total R² shared
      ss_resid_trend_only = ss_resid_trend_only,
      ss_total = ss_total,
      ss_resid = ss_resid,
      converged = TRUE,                 # closed-form OLS: always converged
      convergence = "ols",
      boundary_hit = FALSE,
      fitted = fitted(fit),
      residuals = residuals(fit),
      time = time,
      y = y,
      period = period,
      n_harmonics = n_harmonics,
      n = n
    )
  }

# Nonlinear fitting function for exponential saturation trend
fit_cosinor_nonlinear <- function(time, y, period, n_harmonics, trend_type = "none",
                                   include_inertia = FALSE, wake_onset = NULL,
                                   W0_init = 0.32, tau_W_init = 0.66,
                                   use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                                   amplitude_min = 0, amplitude_max = NA,
                                   A_sat_min = NA, A_sat_max = NA,
                                   tau_min = 0.5, tau_max = NA,
                                   do_loocv = TRUE, tau_fixed = NA) {
  # This function handles exponential saturation trend (exp_sat) and bounded optimization

  n <- length(y)
  t_offset <- min(time)
  t_shifted <- time - t_offset  # For exponential saturation
  t_max <- max(t_shifted)

  # Initial values from data
  y_range <- max(y) - min(y)
  y_min <- min(y)
  y_max <- max(y)
  y_mean <- mean(y)

  # Set default bounds based on data if not specified
  if(use_bounds) {
    if(is.na(mesor_min)) mesor_min <- y_min - y_range
    if(is.na(mesor_max)) mesor_max <- y_max + y_range
    if(is.na(amplitude_max)) amplitude_max <- y_range * 2
    if(is.na(A_sat_min)) A_sat_min <- -Inf  # Allow negative for decreasing trends
    if(is.na(A_sat_max)) A_sat_max <- y_range * 2
    if(is.na(tau_max)) tau_max <- t_max * 5
  } else {
    # No user-specified bounds - use wide defaults with minimal numerical constraints
    # For exp_sat, we still need sensible bounds for numerical stability
    mesor_min <- -Inf
    mesor_max <- Inf
    amplitude_min <- -Inf
    amplitude_max <- Inf

    # For exp_sat: Use original bounds from CIRCAREGold.R when bounding disabled
    # These bounds apply to Approaches 1-2, but Approach 3 will be unbounded
    if(trend_type == "exp_sat") {
      # Original bounds: A_sat unbounded, tau constrained
      A_sat_min <- -Inf
      A_sat_max <- Inf
      tau_min <- 0.5  # Original fixed minimum
      tau_max <- t_max * 5  # Original upper bound
    } else {
      A_sat_min <- -Inf
      A_sat_max <- Inf
      tau_min <- 0.5
      tau_max <- Inf
    }
  }

  # Estimate trend direction from linear regression
  lin_fit <- lm(y ~ t_shifted)
  lin_slope <- coef(lin_fit)[2]
  lin_intercept <- coef(lin_fit)[1]

  # Build formula dynamically based on components
  # Always start with mesor
  formula_parts <- c("mesor")
  start_list <- list(mesor = y_mean)
  lower_bounds <- c(mesor = mesor_min)
  upper_bounds <- c(mesor = mesor_max)

  # Add trend component
  if(trend_type == "linear") {
    formula_parts <- c(formula_parts, "beta_t * time")
    start_list$beta_t <- lin_slope
    lower_bounds["beta_t"] <- -Inf
    upper_bounds["beta_t"] <- Inf

  } else if(trend_type == "log") {
    formula_parts <- c(formula_parts, "beta_log * log(time - t_offset + 1)")
    start_list$beta_log <- lin_slope * log(t_max + 1)
    lower_bounds["beta_log"] <- -Inf
    upper_bounds["beta_log"] <- Inf

  } else if(trend_type == "exp_sat") {
    formula_parts <- c(formula_parts, "A_sat * (1 - exp(-t_shifted / tau))")

    # Better starting values for exp_sat
    if(lin_slope > 0) {
      start_list$A_sat <- y_range * 1.5
      start_list$tau <- t_max / 4
    } else {
      start_list$A_sat <- lin_slope * t_max
      start_list$tau <- t_max / 3
    }
    start_list$tau <- max(1, min(start_list$tau, t_max * 2))

    # Apply user-specified or default bounds
    lower_bounds["A_sat"] <- A_sat_min
    lower_bounds["tau"] <- tau_min
    upper_bounds["A_sat"] <- A_sat_max
    upper_bounds["tau"] <- tau_max
  }

  # Add harmonic components
  harmonic_parts <- sapply(1:n_harmonics, function(h) {
    omega <- 2 * pi * h / period
    start_list[[paste0("b_cos", h)]] <<- 0
    start_list[[paste0("b_sin", h)]] <<- 0

    # Bound harmonic coefficients to respect amplitude constraints
    # Since amplitude = sqrt(b_cos^2 + b_sin^2), bound each coefficient to +/- amplitude_max
    lower_bounds[paste0("b_cos", h)] <<- -amplitude_max
    lower_bounds[paste0("b_sin", h)] <<- -amplitude_max
    upper_bounds[paste0("b_cos", h)] <<- amplitude_max
    upper_bounds[paste0("b_sin", h)] <<- amplitude_max

    sprintf("b_cos%d * cos(%f * time) + b_sin%d * sin(%f * time)", h, omega, h, omega)
  })
  formula_parts <- c(formula_parts, harmonic_parts)

  # Build complete formula
  formula_str <- sprintf("y ~ %s", paste(formula_parts, collapse = " + "))

  # Debug: Print formula being fitted (useful for verification)
  # cat(sprintf("Fitting model: %s\n", formula_str))

  # Ensure starting values respect bounds (clamp them)
  for(param_name in names(start_list)) {
    if(param_name %in% names(lower_bounds)) {
      start_val <- start_list[[param_name]]
      lb <- lower_bounds[param_name]
      ub <- upper_bounds[param_name]

      # Clamp to bounds if finite
      if(!is.infinite(lb) && start_val < lb) {
        start_list[[param_name]] <- lb + (ub - lb) * 0.1  # 10% above lower bound
      }
      if(!is.infinite(ub) && start_val > ub) {
        start_list[[param_name]] <- ub - (ub - lb) * 0.1  # 10% below upper bound
      }
    }
  }

  # Prepare data frame for fitting
  fit_data <- data.frame(
    y = y,
    time = time,
    t_shifted = t_shifted,
    t_offset = t_offset
  )

  # Try fitting with multiple approaches
  fit_success <- FALSE
  nls_fit <- NULL
  error_msgs <- c()
  conv_status <- "not_attempted"
  conv_detail <- NA_character_

  # AUDIT 2.3: nls.control(warnOnly = TRUE) makes nls RETURN a fit object when
  # it hits maxiter instead of raising an error, and nlsLM likewise returns on
  # the iteration limit. The old code set fit_success <- TRUE whenever tryCatch
  # did not fire, so an optimiser that never converged was counted as a success.
  # That is how "Successfully fitted: 1305 / 1305" coexisted with an R² range
  # starting at 0.060, which is not attainable by a converged 8-parameter least
  # squares fit on 16 points.
  #
  # This asks the fit object what actually happened.
  fck_conv_of <- function(fit) {
    if (is.null(fit)) return(list(ok = FALSE, status = "null", detail = NA_character_))
    ci <- tryCatch(fit$convInfo, error = function(e) NULL)
    if (!is.null(ci)) {
      if (isTRUE(ci$isConv))
        return(list(ok = TRUE, status = "converged",
                    detail = sprintf("%d iterations, tol %.3g",
                                     ci$finIter %||% NA, ci$finTol %||% NA)))
      return(list(ok = FALSE, status = "maxiter",
                  detail = ci$stopMessage %||% "did not converge"))
    }
    # minpack.lm carries its own record
    inf <- tryCatch(fit$convInfo$stopCode, error = function(e) NULL)
    ni <- tryCatch(fit$niter, error = function(e) NULL)
    if (!is.null(ni) && is.finite(ni))
      return(list(ok = ni < 300, status = if (ni < 300) "converged" else "maxiter",
                  detail = sprintf("%d iterations", ni)))
    list(ok = TRUE, status = "converged_unverified", detail = NA_character_)
  }

  # Choose fitting strategy based on trend type and whether bounds are enabled
  # NOTE: For exp_sat, Approaches 1-2 use bounds, but Approach 3 falls back to unbounded
  if(use_bounds || trend_type == "exp_sat") {
    # BOUNDED OPTIMIZATION: Use algorithms that support bounds
    # (Either user requested bounds, or exp_sat in Approaches 1-2)

    # Approach 1: Try nlsLM (Levenberg-Marquardt) if available - most robust
    if(requireNamespace("minpack.lm", quietly = TRUE)) {
      tryCatch({
        nls_fit <- minpack.lm::nlsLM(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          lower = lower_bounds,
          upper = upper_bounds,
          control = minpack.lm::nls.lm.control(maxiter = 300)
        )
        .cv <- fck_conv_of(nls_fit)
        conv_status <<- .cv$status; conv_detail <<- .cv$detail
        fit_success <- isTRUE(.cv$ok)
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nlsLM:", e$message))
      })
    }

    # Approach 2: Try standard nls with port algorithm (allows bounds)
    if(!fit_success) {
      tryCatch({
        nls_fit <- nls(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          algorithm = "port",
          lower = lower_bounds,
          upper = upper_bounds,
          control = nls.control(maxiter = 300, warnOnly = TRUE)
        )
        .cv <- fck_conv_of(nls_fit)
        conv_status <<- .cv$status; conv_detail <<- .cv$detail
        fit_success <- isTRUE(.cv$ok)
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nls-port:", e$message))
      })
    }

  } else {
    # UNBOUNDED OPTIMIZATION: Use default algorithms without bounds
    # (Only for linear/log/none trends when user hasn't requested bounds)

    # Approach 1: Try nlsLM without bounds if available
    if(requireNamespace("minpack.lm", quietly = TRUE)) {
      tryCatch({
        nls_fit <- minpack.lm::nlsLM(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          control = minpack.lm::nls.lm.control(maxiter = 300)
        )
        .cv <- fck_conv_of(nls_fit)
        conv_status <<- .cv$status; conv_detail <<- .cv$detail
        fit_success <- isTRUE(.cv$ok)
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nlsLM:", e$message))
      })
    }

    # Approach 2: Try standard nls with default algorithm (no bounds)
    if(!fit_success) {
      tryCatch({
        nls_fit <- nls(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          control = nls.control(maxiter = 300, warnOnly = TRUE)
        )
        .cv <- fck_conv_of(nls_fit)
        conv_status <<- .cv$status; conv_detail <<- .cv$detail
        fit_success <- isTRUE(.cv$ok)
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nls-default:", e$message))
      })
    }
  }

  # Approach 3: Try different starting values for tau parameters
  if(!fit_success && (trend_type == "exp_sat" || include_inertia)) {
    for(tau_mult in c(0.1, 0.5, 2, 5)) {
      if(trend_type == "exp_sat") {
        start_list$tau <- t_max * tau_mult / 3
      }
      if(include_inertia) {
        start_list$tau_W <- tau_W_init * tau_mult
      }

      tryCatch({
        if(use_bounds) {
          # Use port algorithm with bounds (only if user explicitly requested bounds)
          nls_fit <- nls(
            as.formula(formula_str),
            data = fit_data,
            start = start_list,
            algorithm = "port",
            lower = lower_bounds,
            upper = upper_bounds,
            control = nls.control(maxiter = 300, warnOnly = TRUE)
          )
        } else {
          # Use default algorithm without bounds (original fallback behavior)
          # This matches CIRCAREGold.R Approach 3 - unbounded even for exp_sat
          nls_fit <- nls(
            as.formula(formula_str),
            data = fit_data,
            start = start_list,
            control = nls.control(maxiter = 300, warnOnly = TRUE)
          )
        }

        # Check if fit is reasonable (R² > 0)
        fitted_check <- predict(nls_fit)
        ss_tot_check <- sum((y - mean(y))^2)
        ss_res_check <- sum((y - fitted_check)^2)
        r2_check <- 1 - ss_res_check / ss_tot_check

        .cv <- fck_conv_of(nls_fit)
        if(r2_check > 0 && isTRUE(.cv$ok)) {
          conv_status <- .cv$status; conv_detail <- .cv$detail
          fit_success <- TRUE
          break
        }
        # remember the best non-converged attempt so it can be REPORTED as
        # non-converged rather than silently discarded or silently accepted
        conv_status <- .cv$status; conv_detail <- .cv$detail
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("retry", tau_mult, ":", e$message))
      })
    }
  }

  if(!fit_success || is.null(nls_fit)) {
    return(list(
      success = FALSE,
      converged = FALSE,
      convergence = if (identical(conv_status, "not_attempted")) "failed" else conv_status,
      convergence_detail = conv_detail,
      message = sprintf("Nonlinear fit failed to converge (%s). Errors: %s",
                       conv_status, paste(head(error_msgs, 3), collapse = "; "))
    ))
  }

  # Extract results
  tryCatch({
    coefs <- coef(nls_fit)
    se <- tryCatch(summary(nls_fit)$coefficients[, 2],
                   error = function(e) rep(NA, length(coefs)))
    names(se) <- names(coefs)

    mesor <- coefs["mesor"]
    mesor_se <- se["mesor"]

    # Extract trend parameters based on type
    trend_params <- list()
    if(trend_type == "linear") {
      trend_params$trend_linear <- list(
        coef = as.numeric(coefs["beta_t"]),
        se = as.numeric(se["beta_t"])
      )
    } else if(trend_type == "log") {
      trend_params$trend_log <- list(
        coef = as.numeric(coefs["beta_log"]),
        se = as.numeric(se["beta_log"])
      )
    } else if(trend_type == "exp_sat") {
      trend_params$A_sat <- list(
        coef = as.numeric(coefs["A_sat"]),
        se = if(!is.na(se["A_sat"])) as.numeric(se["A_sat"]) else NA
      )
      trend_params$tau <- list(
        coef = as.numeric(coefs["tau"]),
        se = if(!is.na(se["tau"])) as.numeric(se["tau"]) else NA
      )
    }

    # No inertia parameters for exp_sat-only fitting
    inertia_params <- NULL

    # Extract harmonic parameters
    amplitudes <- numeric(n_harmonics)
    acrophases <- numeric(n_harmonics)
    amp_se <- numeric(n_harmonics)
    acro_se <- numeric(n_harmonics)

    # AUDIT (extra c) + 2.5: the old SEs here were
    #     amp_se  <- sqrt(se_cos^2 + se_sin^2) / sqrt(2)
    #     acro_se <- NA          # "Complex for nonlinear"
    # The first is neither the delta-method result nor an approximation to it --
    # it discards the cos/sin covariance entirely and mis-weights both terms.
    # The second is why there were no acrophase confidence intervals anywhere in
    # the app for an exp_sat model: none had ever been computed.
    #
    # nls supplies a full covariance matrix. The delta method applies to a
    # nonlinear fit exactly as it does to a linear one, so both SEs are now
    # computed the same way Fitter A computes them.
    Vfull <- tryCatch(stats::vcov(nls_fit), error = function(e) NULL)
    bingham <- vector("list", n_harmonics)

    for(h in 1:n_harmonics) {
      beta_cos <- as.numeric(coefs[paste0("b_cos", h)])
      beta_sin <- as.numeric(coefs[paste0("b_sin", h)])
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi

      Vh <- NULL
      if(!is.null(Vfull)) {
        nmc <- paste0("b_cos", h); nms <- paste0("b_sin", h)
        if(all(c(nmc, nms) %in% rownames(Vfull)))
          Vh <- Vfull[c(nmc, nms), c(nmc, nms), drop = FALSE]
      }
      amp_se[h]  <- fck_amp_se(beta_cos, beta_sin, Vh)
      acro_se[h] <- fck_acro_se(beta_cos, beta_sin, Vh)

      # Bingham elliptical joint region for (amplitude, acrophase)
      bingham[[h]] <- fck_bingham_ci(beta_cos, beta_sin, Vh, n, length(coefs),
                                     level = 0.95, period = period, harmonic = h)
    }

    # AUDIT (extra i): divide by h -- harmonic h lives on the effective period T/h.
    acrophases_time <- vapply(seq_len(n_harmonics),
                              function(h) phi_to_hours(acrophases[h], period, h),
                              numeric(1))
    acro_se_time <- vapply(seq_len(n_harmonics),
                           function(h) phi_to_hours(acro_se[h], period, h),
                           numeric(1))

    # Goodness of fit
    fitted_vals <- predict(nls_fit)
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum((y - fitted_vals)^2)
    r_squared <- 1 - ss_resid / ss_total
    percent_rhythm <- max(0, r_squared * 100)

    # Calculate p-value for circadian rhythm
    n_params <- length(coefs)

    # AUDIT (extra a): test the harmonics GIVEN the trend, by refitting the
    # trend-only model. See the note in fit_cosinor() -- the old numerator was
    # the whole model's SS over the harmonics' df alone.
    ss_resid_trend_only <- tryCatch({
      if(trend_type == "exp_sat" && !is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
        # refit M and A_sat with tau held at its full-model value: a 2-parameter
        # linear problem, so no second optimisation can fail here
        b <- 1 - exp(-(time - t_offset) / trend_params$tau$coef)
        sum(residuals(lm(y ~ b))^2)
      } else if(trend_type == "linear") {
        sum(residuals(lm(y ~ time))^2)
      } else if(trend_type == "log") {
        sum(residuals(lm(y ~ log(time - t_offset + 1)))^2)
      } else {
        ss_total
      }
    }, error = function(e) ss_total)

    .zt <- fck_zero_amplitude_test(ss_resid, ss_resid_trend_only,
                                   n, n_params, n_harmonics)
    f_stat  <- .zt$F
    p_value <- .zt$p
    df2 <- .zt$df2

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian nonlinear model
    sigma_sq <- ss_resid / n
    log_lik <- -n/2 * (log(2*pi) + log(sigma_sq) + 1)

    # AIC: Akaike Information Criterion
    k <- n_params + 1  # parameters + sigma
    aic <- -2 * log_lik + 2 * k

    # AICc: Corrected AIC for small samples
    aicc <- if(n - k - 1 > 0) {
      aic + (2 * k * (k + 1)) / (n - k - 1)
    } else {
      NA
    }

    # BIC: Bayesian Information Criterion
    bic <- -2 * log_lik + k * log(n)

    # LOOCV: Leave-one-out cross-validation
    # Note: For nonlinear models, this is computationally expensive
    # We'll use a simplified approach: predict each point using the full model
    # and apply a leave-one-out correction based on leverage
    # AUDIT (extra b): this used to be
    #     residuals_vec <- y - fitted_vals; sqrt(mean(residuals_vec^2))
    # i.e. the IN-SAMPLE residual RMSE, returned under the name "LOOCV RMSE".
    # Every cross-validation number the report printed for an exp_sat model was
    # training error, while Fitter A under the same label did real LOOCV -- two
    # different quantities sharing one column.
    #
    # With n = 16 a genuine leave-one-out refit is 16 nls calls per subject,
    # which is affordable, so it is now done properly. If a refit fails the
    # point is dropped from the average and the count is reported, rather than
    # the whole thing silently degrading to training error.
    loocv_rmse <- NA_real_
    loocv_n_failed <- 0L
    loocv_is_true_cv <- FALSE
    if(isTRUE(do_loocv)) {
      errs <- rep(NA_real_, n)
      for(i in seq_len(n)) {
        d_loo <- fit_data[-i, , drop = FALSE]
        f_loo <- tryCatch({
          if(requireNamespace("minpack.lm", quietly = TRUE)) {
            minpack.lm::nlsLM(as.formula(formula_str), data = d_loo,
                              start = as.list(coefs),
                              lower = lower_bounds, upper = upper_bounds,
                              control = minpack.lm::nls.lm.control(maxiter = 200))
          } else {
            nls(as.formula(formula_str), data = d_loo, start = as.list(coefs),
                control = nls.control(maxiter = 200, warnOnly = TRUE))
          }
        }, error = function(e) NULL)
        if(is.null(f_loo)) { loocv_n_failed <- loocv_n_failed + 1L; next }
        yhat <- tryCatch(as.numeric(predict(f_loo, newdata = fit_data[i, , drop = FALSE])),
                         error = function(e) NA_real_)
        if(is.finite(yhat)) errs[i] <- (y[i] - yhat)^2 else loocv_n_failed <- loocv_n_failed + 1L
      }
      if(any(is.finite(errs))) {
        loocv_rmse <- sqrt(mean(errs, na.rm = TRUE))
        loocv_is_true_cv <- TRUE
      }
    } else {
      # explicitly labelled as what it is, never as cross-validation
      loocv_rmse <- sqrt(mean((y - fitted_vals)^2))
      loocv_is_true_cv <- FALSE
    }

    # ===========================================================================
    # Variance decomposition: Calculate R² for Process S and Process C separately
    # ===========================================================================
    r_squared_S <- 0  # Variance explained by homeostatic trend alone
    r_squared_C <- 0  # Variance explained by circadian rhythm alone
    percent_S <- 0    # Percentage of total R² from Process S
    percent_C <- 0    # Percentage of total R² from Process C

    # Model with trend only (Process S) - using exp_sat trend
    # AUDIT: the marginal R² of Process S has to come from a model REFIT without
    # the harmonics. The old code reused the FULL model's mesor and A_sat and
    # simply deleted the cosine terms from the prediction, which is not the
    # trend-only fit -- it is the full fit with part of it thrown away, and it
    # understates R²_S (and so overstates the shared component). ss_resid_trend_only
    # above is the proper refit and is reused here so the F test and the
    # decomposition cannot disagree about what "Process S alone" means.
    if(trend_type != "none" && is.finite(ss_resid_trend_only) && ss_total > 0) {
      r_squared_S <- max(0, 1 - ss_resid_trend_only / ss_total)
    }

    # Model with circadian only (Process C)
    # Build a model with MESOR + harmonics only (no trend)
    tryCatch({
      X_circ <- matrix(1, nrow = n, ncol = 1)
      for(h in 1:n_harmonics) {
        omega <- 2 * pi * h / period
        X_circ <- cbind(X_circ, cos(omega * time), sin(omega * time))
      }
      fit_circ <- lm(y ~ X_circ - 1)
      ss_resid_circ <- sum(residuals(fit_circ)^2)
      r_squared_C <- max(0, 1 - ss_resid_circ / ss_total)
    }, error = function(e) {
      r_squared_C <<- 0
    })

    # AUDIT 1.3: commonality, not two overlapping marginals. See fit_cosinor().
    .cm <- fck_commonality(r_squared, r_squared_S, r_squared_C)
    unique_S <- .cm$unique_S
    unique_C <- .cm$unique_C
    shared_SC <- .cm$shared
    .pc <- fck_commonality_pct(.cm)
    percent_S <- .pc$unique_S
    percent_C <- .pc$unique_C
    percent_shared <- .pc$shared

    # which bounds, if any, this fit ended up sitting on
    .bh <- fck_bounds_hit(coefs, lower_bounds, upper_bounds)

    list(
      success = TRUE,
      mesor = as.numeric(mesor),
      mesor_se = if(!is.na(mesor_se)) as.numeric(mesor_se) else NA,
      trend_type = trend_type,
      trend_params = trend_params,
      inertia_params = inertia_params,  # NEW: Sleep inertia parameters
      t_offset = t_offset,
      t_center = 0,
      amplitudes = amplitudes,
      amp_se = amp_se,
      acrophases = acrophases,
      acrophases_time = acrophases_time,
      acro_se = acro_se,
      acro_se_time = acro_se_time,
      coefs = coefs,
      se = se,
      vcov = NULL,
      r_squared = r_squared,
      adj_r_squared = r_squared,  # Approximate
      percent_rhythm = percent_rhythm,
      f_stat = f_stat,
      p_value = p_value,
      aic = aic,                        # Akaike Information Criterion
      aicc = aicc,                      # Corrected AIC for small samples
      bic = bic,                        # Bayesian Information Criterion
      loocv_rmse = loocv_rmse,          # genuine LOOCV when loocv_is_true_cv
      loocv_is_true_cv = loocv_is_true_cv,
      loocv_n_failed = loocv_n_failed,
      r_squared_S = r_squared_S,        # MARGINAL R² of trend-only model
      r_squared_C = r_squared_C,        # MARGINAL R² of harmonics-only model
      unique_S = unique_S,              # commonality: unique to Process S
      unique_C = unique_C,              # commonality: unique to Process C
      shared_SC = shared_SC,            # commonality: shared (may be negative)
      percent_S = percent_S,            # % of total R² UNIQUE to S
      percent_C = percent_C,            # % of total R² UNIQUE to C
      percent_shared = percent_shared,  # % of total R² shared
      ss_resid_trend_only = ss_resid_trend_only,
      ss_total = ss_total,
      ss_resid = ss_resid,
      bingham = bingham,                # elliptical joint CI per harmonic
      vcov_full = Vfull,
      converged = TRUE,
      convergence = conv_status,
      convergence_detail = conv_detail,
      # AUDIT 2.3: a fit sitting ON a bound is not an estimate -- it is the
      # optimiser being stopped by the constraint, and the standard error is
      # meaningless there. Reported per PARAMETER rather than as one logical,
      # because "tau ran to its ceiling" and "the amplitude hit its cap" are
      # different problems, and because a bound that catches most of the sample
      # is a badly chosen bound rather than a sample full of odd subjects.
      # Whether these fits are excluded is the user's choice, not this
      # function's: it only reports what happened.
      bounds_hit = .bh,
      n_bounds_hit = length(.bh),
      boundary_hit = length(.bh) > 0,
      bounds_lower = lower_bounds,
      bounds_upper = upper_bounds,
      fitted = fitted_vals,
      residuals = y - fitted_vals,
      time = time,
      y = y,
      period = period,
      n_harmonics = n_harmonics,
      n = n
    )
  }, error = function(e) {
    list(success = FALSE, message = paste("Result extraction failed:", e$message))
  })
}

  # Predict from cosinor model
  predict_cosinor <- function(fit, newtime = NULL, component = "total", include_trend_in_pred = TRUE) {
    if(is.null(newtime)) newtime <- fit$time
    period <- fit$period
    n_harmonics <- fit$n_harmonics
    coefs <- fit$coefs
    trend_type <- if(!is.null(fit$trend_type)) fit$trend_type else "none"
    
    # Handle legacy format (include_trend boolean)
    if(is.null(fit$trend_type) && isTRUE(fit$include_trend)) {
      trend_type <- "linear"
    }
    
    pred <- rep(coefs[1], length(newtime))  # MESOR
    
    # Calculate trend component based on type
    if(include_trend_in_pred && trend_type != "none") {
      t_offset <- if(!is.null(fit$t_offset)) fit$t_offset else min(newtime)
      
      trend_val <- switch(trend_type,
                          "linear" = coefs[2] * newtime,
                          "log" = coefs[2] * log(newtime - t_offset + 1),
                          "exp_sat" = {
                            A_sat <- coefs["A_sat"]
                            tau <- coefs["tau"]
                            t_shifted <- newtime - t_offset
                            A_sat * (1 - exp(-t_shifted / tau))
                          },
                          "two_process" = {
                            if(!is.null(fit$S_trajectory) && !is.null(fit$time) &&
                               "beta_S" %in% names(coefs)) {
                              S_interp <- tryCatch({
                                approx(fit$time, fit$S_trajectory, xout = newtime,
                                       rule = 2, ties = "ordered")$y
                              }, error = function(e) rep(NA_real_, length(newtime)))
                              coefs["beta_S"] * S_interp
                            } else {
                              rep(0, length(newtime))
                            }
                          },
                          rep(0, length(newtime))
      )
      pred <- pred + trend_val
    }
    
    # Determine coefficient offset based on trend type
    n_trend_coefs <- switch(trend_type,
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, "two_process" = 0, 0)
    coef_offset <- 1 + n_trend_coefs

    if(component == "total" || component == "all") {
      for(h in 1:n_harmonics) {
        omega <- 2 * pi * h / period
        if(trend_type == "exp_sat" || trend_type == "two_process") {
          beta_cos <- coefs[paste0("b_cos", h)]
          beta_sin <- coefs[paste0("b_sin", h)]
        } else {
          cos_idx <- coef_offset + 2 * (h - 1) + 1
          sin_idx <- coef_offset + 2 * (h - 1) + 2
          beta_cos <- coefs[cos_idx]
          beta_sin <- coefs[sin_idx]
        }
        pred <- pred + beta_cos * cos(omega * newtime) + beta_sin * sin(omega * newtime)
      }
    } else if(is.numeric(component) && component >= 1 && component <= n_harmonics) {
      h <- component
      omega <- 2 * pi * h / period
      if(trend_type == "exp_sat" || trend_type == "two_process") {
        beta_cos <- coefs[paste0("b_cos", h)]
        beta_sin <- coefs[paste0("b_sin", h)]
      } else {
        cos_idx <- coef_offset + 2 * (h - 1) + 1
        sin_idx <- coef_offset + 2 * (h - 1) + 2
        beta_cos <- coefs[cos_idx]
        beta_sin <- coefs[sin_idx]
      }
      pred <- coefs[1] + beta_cos * cos(omega * newtime) + beta_sin * sin(omega * newtime)
    }
    
    return(pred)
  }
  
  # Get harmonic components separately
  get_harmonic_components <- function(fit, newtime = NULL) {
    if(is.null(newtime)) newtime <- fit$time
    period <- fit$period
    n_harmonics <- fit$n_harmonics
    coefs <- fit$coefs
    trend_type <- if(!is.null(fit$trend_type)) fit$trend_type else "none"
    
    # Handle legacy format
    if(is.null(fit$trend_type) && isTRUE(fit$include_trend)) {
      trend_type <- "linear"
    }
    
    t_offset <- if(!is.null(fit$t_offset)) fit$t_offset else min(newtime)
    
    # Determine coefficient offset
    n_trend_coefs <- switch(trend_type,
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, "two_process" = 0, 0)
    coef_offset <- 1 + n_trend_coefs
    
    components <- list()
    components$mesor <- rep(coefs[1], length(newtime))
    components$trend_type <- trend_type
    
    if(trend_type != "none") {
      # Compute trend component
      trend_val <- switch(trend_type,
                          "linear" = coefs[2] * newtime,
                          "log" = coefs[2] * log(newtime - t_offset + 1),
                          "exp_sat" = {
                            A_sat <- coefs["A_sat"]
                            tau <- coefs["tau"]
                            A_sat * (1 - exp(-(newtime - t_offset) / tau))
                          },
                          "two_process" = {
                            if(!is.null(fit$S_trajectory) && !is.null(fit$time) &&
                               "beta_S" %in% names(coefs)) {
                              S_interp <- tryCatch({
                                approx(fit$time, fit$S_trajectory, xout = newtime,
                                       rule = 2, ties = "ordered")$y
                              }, error = function(e) rep(NA_real_, length(newtime)))
                              coefs["beta_S"] * S_interp
                            } else {
                              rep(0, length(newtime))
                            }
                          },
                          rep(0, length(newtime))
      )
      components$trend <- trend_val
    }

    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      if(trend_type == "exp_sat" || trend_type == "two_process") {
        beta_cos <- coefs[paste0("b_cos", h)]
        beta_sin <- coefs[paste0("b_sin", h)]
      } else {
        cos_idx <- coef_offset + 2 * (h - 1) + 1
        sin_idx <- coef_offset + 2 * (h - 1) + 2
        beta_cos <- coefs[cos_idx]
        beta_sin <- coefs[sin_idx]
      }
      components[[paste0("harmonic_", h)]] <- beta_cos * cos(omega * newtime) +
        beta_sin * sin(omega * newtime)
    }

    components$total <- predict_cosinor(fit, newtime)
    return(components)
  }
  
  # Predict curve from mean coefficients (for group/population means with all harmonics)
  predict_from_coefs <- function(coefs, time_vec, period, n_harmonics, trend_type = "none", 
                                 t_offset = 0, t_center = 0) {
    # Handle legacy boolean format
    if(is.logical(trend_type)) {
      trend_type <- if(trend_type) "linear" else "none"
    }
    
    # coefs format: c(mesor, [trend_coefs...], beta_cos_1, beta_sin_1, ...)
    pred <- rep(coefs[1], length(time_vec))  # MESOR
    
    # Determine trend offset
    n_trend_coefs <- switch(as.character(trend_type),
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
    coef_offset <- 1 + n_trend_coefs
    
    # Add trend based on type
    if(trend_type != "none" && trend_type != FALSE) {
      trend_val <- switch(as.character(trend_type),
                          "linear" = coefs[2] * time_vec,
                          "log" = coefs[2] * log(time_vec - t_offset + 1),
                          "exp_sat" = coefs[2] * (1 - exp(-(time_vec - t_offset) / coefs[3])),  # A_sat * (1 - exp(-t/tau))
                          rep(0, length(time_vec))
      )
      pred <- pred + trend_val
    }
    
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      beta_cos <- coefs[coef_offset + 2 * h - 1]
      beta_sin <- coefs[coef_offset + 2 * h]
      pred <- pred + beta_cos * cos(omega * time_vec) + beta_sin * sin(omega * time_vec)
    }
    return(pred)
  }
  
  # Get trend value at a specific time point
  get_trend_value <- function(trend_type, trend_params, time_vec, t_offset = 0) {
    if(trend_type == "none" || is.null(trend_params) || length(trend_params) == 0) {
      return(rep(0, length(time_vec)))
    }
    
    switch(trend_type,
           "linear" = trend_params$trend_linear$coef * time_vec,
           "log" = trend_params$trend_log$coef * log(time_vec - t_offset + 1),
           "exp_sat" = {
             A_sat <- trend_params$A_sat$coef
             tau <- trend_params$tau$coef
             A_sat * (1 - exp(-(time_vec - t_offset) / tau))
           },
           rep(0, length(time_vec))
    )
  }
  
  # Get human-readable trend label
  get_trend_label <- function(trend_type, prefix = "") {
    label <- switch(trend_type,
                    "linear" = "Linear Trend",
                    "log" = "Log Trend",
                    "exp_sat" = "Exp. Saturation",
                    "two_process" = "Process S",
                    "Process S"  # default
    )
    if(nchar(prefix) > 0) paste(label, prefix) else label
  }
  
  # Helper: Get mean trend coefficients from individual parameters for building coef vectors
  # Returns a vector of trend coefficients to append to mesor for predict_from_coefs
  get_mean_trend_coefs <- function(params, trend_type) {
    if(trend_type == "none") return(numeric(0))
    
    if(trend_type == "linear" && "trend_linear" %in% names(params)) {
      return(mean(params$trend_linear, na.rm = TRUE))
    } else if(trend_type == "log" && "trend_log" %in% names(params)) {
      return(mean(params$trend_log, na.rm = TRUE))
    } else if(trend_type == "exp_sat") {
      coefs <- numeric(0)
      if("A_sat" %in% names(params)) coefs <- c(coefs, mean(params$A_sat, na.rm = TRUE))
      if("tau" %in% names(params)) coefs <- c(coefs, mean(params$tau, na.rm = TRUE))
      return(coefs)
    }
    return(numeric(0))
  }
  
  # Helper: Check if trend columns exist in params
  has_trend_params <- function(params, trend_type) {
    if(trend_type == "none") return(FALSE)
    if(trend_type == "linear") return("trend_linear" %in% names(params))
    if(trend_type == "log") return("trend_log" %in% names(params))
    if(trend_type == "exp_sat") return("A_sat" %in% names(params) || "tau" %in% names(params))
    if(trend_type == "two_process") {
      return("beta_S" %in% names(params) && "tau_w" %in% names(params) && "tau_s" %in% names(params))
    }
    return(FALSE)
  }
  
  # Helper: Get the primary trend column name for a given trend_type
  get_trend_col <- function(trend_type) {
    switch(trend_type,
           "linear" = "trend_linear",
           "log" = "trend_log",
           "exp_sat" = "A_sat",
           NULL)
  }
  
  # Helper: Compute trend line values for plotting
  compute_trend_line <- function(params, trend_type, time_vec, t_offset = 0) {
    if(trend_type == "none" || !has_trend_params(params, trend_type)) {
      return(NULL)
    }

    mesor <- mean(params$mesor, na.rm = TRUE)

    if(trend_type == "linear" && "trend_linear" %in% names(params)) {
      slope <- mean(params$trend_linear, na.rm = TRUE)
      return(mesor + slope * time_vec)
    } else if(trend_type == "log" && "trend_log" %in% names(params)) {
      slope <- mean(params$trend_log, na.rm = TRUE)
      return(mesor + slope * log(time_vec - t_offset + 1))
    } else if(trend_type == "exp_sat" && "A_sat" %in% names(params) && "tau" %in% names(params)) {
      A_sat <- mean(params$A_sat, na.rm = TRUE)
      tau <- mean(params$tau, na.rm = TRUE)
      return(mesor + A_sat * (1 - exp(-(time_vec - t_offset) / tau)))
    }
    return(NULL)
  }

  # ==============================================================================
  # MAIN HARMONIC REGRESSION EVENT HANDLER
  # ==============================================================================
  
  observeEvent(input$run_harmonic, {
    req(values$data)
    
    showNotification("Running Harmonic Regression...", type = "message", duration = 2)
    
    tryCatch({
      # Check if smoothed data is available
      # ======================================================================
      # AUDIT 2.1: fitting on FDA-smoothed data invalidates the fit statistics
      #
      # This used to be an unconditional "use the smoothed matrix if one
      # exists". The user had no say, and the report said only that the data
      # were smoothed, without saying what that costs. Smoothing removes
      # independent noise and induces residual autocorrelation, so R² is
      # inflated, a held-out point is partly reconstructed from its neighbours
      # (making LOOCV optimistic), and the zero-amplitude F test is
      # anticonservative. The "95.3% significant rhythms" figure is an upper
      # bound, not an estimate.
      #
      # Cosinor handles missing and unequally spaced data natively, so raw is a
      # legitimate choice. data_source now selects it; "both" fits twice and the
      # report prints the two side by side so the inflation is visible rather
      # than argued about.
      # ======================================================================
      data_source <- input$harmonic_data_source %||% "smoothed"
      if(identical(data_source, "smoothed") && is.null(values$smooth_data)) {
        data_source <- "raw"
      }
      using_smoothed <- identical(data_source, "smoothed")
      Y <- if(using_smoothed) values$smooth_data else values$data
      n_subjects <- nrow(Y)
      n_time <- ncol(Y)
      period <- input$harmonic_period
      n_harmonics <- input$n_harmonics
      trend_type <- input$harmonic_trend_type
      
      # Calculate number of trend parameters
      n_trend_params <- switch(trend_type,
                               "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
      
      # Diagnostic: Check data type and dimensions
      cat(sprintf("Data diagnostics: %d subjects × %d time points, type=%s, smoothed=%s, trend=%s\n", 
                  n_subjects, n_time, typeof(Y), using_smoothed, trend_type))
      
      # Check for NAs in the data
      na_counts <- apply(Y, 1, function(row) sum(is.na(row)))
      valid_counts <- n_time - na_counts
      subjects_with_nas <- sum(na_counts > 0)
      subjects_all_na <- sum(na_counts == n_time)
      
      if(subjects_all_na > 0) {
        all_na_subjects <- which(na_counts == n_time)
        showNotification(
          sprintf("ERROR: %d subjects have ALL missing values (subjects: %s). Check data selection!", 
                  subjects_all_na, paste(head(all_na_subjects, 10), collapse=", ")),
          type = "error", duration = 15)
        cat(sprintf("Subjects with all NA: %s\n", paste(all_na_subjects, collapse=", ")))
        
        # Show sample of data for first all-NA subject
        if(length(all_na_subjects) > 0) {
          cat(sprintf("First all-NA subject (%d) data sample: %s\n", 
                      all_na_subjects[1], 
                      paste(head(Y[all_na_subjects[1], ], 10), collapse=", ")))
        }
      } else if(subjects_with_nas > 0 && !using_smoothed) {
        showNotification(
          sprintf("%d subjects have missing values. Consider applying smoothing first to interpolate missing data.", 
                  subjects_with_nas),
          type = "warning", duration = 8)
      } else if(subjects_with_nas > 0 && using_smoothed) {
        # This shouldn't happen if smoothing worked correctly
        showNotification(
          sprintf("Warning: %d subjects still have NAs after smoothing. Some fits may fail.", 
                  subjects_with_nas),
          type = "warning", duration = 8)
      }
      
      # Check if we have enough data points for the requested harmonics
      min_required <- 2 * n_harmonics + 2 + n_trend_params
      max_safe_harmonics <- floor((n_time - 2 - n_trend_params) / 2)
      
      if(n_time < min_required) {
        showNotification(
          sprintf("Error: %d harmonics require at least %d time points. You have %d. Maximum safe: %d harmonics.", 
                  n_harmonics, min_required, n_time, max_safe_harmonics),
          type = "error", duration = 10)
        return()
      }
      
      if(n_harmonics > max_safe_harmonics) {
        showNotification(
          sprintf("Warning: Using %d harmonics with only %d time points may cause overfitting. Consider reducing to %d harmonics.", 
                  n_harmonics, n_time, max_safe_harmonics),
          type = "warning", duration = 8)
      }
      
      # Determine time variable
      original_times <- NULL
      wrap_applied <- FALSE
      
      if(input$harmonic_time_var == "_shared_") {
        # MERGED APP: the shared import step already parsed real clock hours
        # from the column names into values$time_clock (03_helpers_clock.R).
        # Use those, so this tab and the manual-entry route agree.
        if(is.null(values$time_clock) || length(values$time_clock) != n_time) {
          showNotification(
            "No shared clock times available: the column names did not yield hours in [0, 24). Use 'Specify times manually'.",
            type = "error", duration = 10)
          return()
        }
        time_vec <- as.numeric(values$time_clock)
        original_times <- time_vec

      } else if(input$harmonic_time_var == "_index_") {
        # Use column indices scaled to period (assumes equal spacing!)
        time_vec <- seq(0, period * (n_time - 1) / n_time, length.out = n_time)
        original_times <- time_vec
        showNotification("Using equally-spaced time points. If your data has unequal spacing, use 'Specify times manually'.", 
                         type = "warning", duration = 5)
        
      } else if(input$harmonic_time_var == "_manual_") {
        # Parse manual time input
        manual_input <- input$harmonic_manual_times
        if(is.null(manual_input) || nchar(trimws(manual_input)) == 0) {
          showNotification("Please enter time values!", type = "error")
          return()
        }
        
        # Parse comma-separated values
        time_vec <- tryCatch({
          vals <- as.numeric(unlist(strsplit(gsub(" ", "", manual_input), ",")))
          if(any(is.na(vals))) stop("Non-numeric values")
          vals
        }, error = function(e) {
          showNotification("Could not parse time values. Use comma-separated numbers (e.g., 8,9,10,11,12,14,16,18,20).", 
                           type = "error")
          return(NULL)
        })
        
        if(is.null(time_vec)) return()
        
        if(length(time_vec) != n_time) {
          showNotification(paste0("Number of time values (", length(time_vec), 
                                  ") must match number of columns (", n_time, ")!"), 
                           type = "error")
          return()
        }
        
        # Detect wrap-around: if a time is smaller than the previous, add period
        # This handles cases like 8,9,10,...,22,23,0,2,4,6 → 8,9,10,...,22,23,24,26,28,30
        original_times <- time_vec
        for(i in 2:length(time_vec)) {
          if(time_vec[i] < time_vec[i-1]) {
            # Wrap-around detected - add period to this and all subsequent values
            time_vec[i:length(time_vec)] <- time_vec[i:length(time_vec)] + period
          }
        }
        
        # Check if wrap-around was applied
        wrap_applied <- !identical(original_times, time_vec)
        if(wrap_applied) {
          showNotification(paste0("Detected wrap-around at midnight. Adjusted times: ", 
                                  paste(round(time_vec, 1), collapse=", ")), 
                           type = "message", duration = 5)
        } else {
          showNotification(paste("Using manual time points:", paste(round(time_vec, 1), collapse=", ")), 
                           type = "message", duration = 3)
        }
        
      } else {
        # Use selected covariate column
        time_vec <- values$covariates[[input$harmonic_time_var]]
        original_times <- time_vec
        wrap_applied <- FALSE
        if(length(time_vec) != n_time) {
          showNotification("Selected time variable doesn't match data dimensions. Using equal spacing.", 
                           type = "warning")
          time_vec <- seq(0, period * (n_time - 1) / n_time, length.out = n_time)
          original_times <- time_vec
        } else {
          # Apply wrap-around detection for covariate time variables too
          for(i in 2:length(time_vec)) {
            if(time_vec[i] < time_vec[i-1]) {
              time_vec[i:length(time_vec)] <- time_vec[i:length(time_vec)] + period
            }
          }
          wrap_applied <- !identical(original_times, time_vec)
          if(wrap_applied) {
            showNotification(paste0("Detected wrap-around (period=", period, "). Adjusted times: ", 
                                    paste(round(time_vec, 1), collapse=", ")), 
                             type = "message", duration = 5)
          }
        }
      }
      
      # Check for potential issues with time values
      if(max(time_vec) > period * 1.5 && !wrap_applied) {
        showNotification(paste0("Note: Max time value (", round(max(time_vec), 1), 
                                ") is larger than period (", period, "). Values will be wrapped using modulo."), 
                         type = "warning", duration = 5)
      }
      
      # Individual cosinor analysis
      individual_fits <- list()
      
      # Build column names for all harmonics
      # AUDIT: beta_0 and the MESOR are different quantities and both are
      # comparable between groups, so both are stored per subject. beta_0 is the
      # fitted constant -- the rhythm's own level, the thing a trend-free cosinor
      # would call the MESOR. mesor_adj is the rhythm-adjusted mean: the
      # time-average of beta_0 + S(t) across the observed window, which is where
      # the data actually sit once the homeostatic rise is counted. A group can
      # rank differently on the two.
      param_cols <- c("subject", "mesor", "mesor_se", "mesor_adj", "value_at_start")

      # Add trend columns based on type
      if(trend_type == "linear") {
        param_cols <- c(param_cols, "trend_linear", "trend_linear_se")
      } else if(trend_type == "log") {
        param_cols <- c(param_cols, "trend_log", "trend_log_se")
      } else if(trend_type == "exp_sat") {
        param_cols <- c(param_cols, "A_sat", "A_sat_se", "tau", "tau_se")
      }

      for(h in 1:n_harmonics) {
        param_cols <- c(param_cols,
                        paste0("amplitude_", h), paste0("amp_se_", h),
                        paste0("acrophase_rad_", h), paste0("acrophase_time_", h),
                        paste0("acro_se_time_", h),
                        paste0("beta_cos_", h), paste0("beta_sin_", h))
      }
      param_cols <- c(param_cols, "r_squared", "percent_rhythm", "p_value",
                      "r_squared_S", "r_squared_C",
                      # AUDIT 1.3: the commonality partition, which sums to total R²
                      "unique_S", "unique_C", "shared_SC",
                      "percent_S", "percent_C", "percent_shared",
                      # AUDIT 2.3: convergence is now recorded per subject
                      "converged", "boundary_hit")
      
      individual_params <- data.frame(matrix(ncol = length(param_cols), nrow = 0))
      colnames(individual_params) <- param_cols
      
      # Coefficient offset for trend
      coef_offset <- 1 + n_trend_params
      
      # Track failed fits
      failed_fits <- list()
      
      # ======================================================================
      # AUDIT 1.4.3 + 2.2: the time origin
      #
      # As shipped, the two halves of the model use DIFFERENT origins:
      # fit_cosinor_nonlinear() builds the trend on (t - min(t)) while the
      # harmonics run on raw t. The constant is therefore the intercept of a
      # model with two anchors and is interpretable as neither "the value at
      # midnight" nor "the value at the first observation".
      #
      # time_origin = "first_observation" shifts BOTH halves to the first
      # observation, which is what makes the intercept mean something and what
      # improves the conditioning: over t in [8, 30] the saturating term is
      # nearly collinear with the constant, and re-anchoring removes the part of
      # that collinearity which is pure offset.
      #
      # Default is "midnight" -- the current behaviour -- so nothing downstream
      # changes unless the user asks for it.
      # ======================================================================
      time_origin <- input$harmonic_time_origin %||% "midnight"
      time_vec_model <- time_vec
      origin_shift <- 0
      if(identical(time_origin, "first_observation")) {
        origin_shift <- min(time_vec, na.rm = TRUE)
        time_vec_model <- time_vec - origin_shift
      }

      # Store time offsets for prediction
      t_offset_global <- min(time_vec_model)
      t_center_global <- mean(time_vec_model)

      # Read parameter bounding options from UI
      use_bounds <- isTRUE(input$harmonic_use_bounds)
      mesor_min <- if(use_bounds) input$harmonic_mesor_min else NA
      mesor_max <- if(use_bounds) input$harmonic_mesor_max else NA
      amplitude_min <- if(use_bounds) input$harmonic_amplitude_min else 0
      amplitude_max <- if(use_bounds) input$harmonic_amplitude_max else NA
      A_sat_min <- if(use_bounds) input$harmonic_A_sat_min else NA
      A_sat_max <- if(use_bounds) input$harmonic_A_sat_max else NA
      tau_min <- if(use_bounds) input$harmonic_tau_min else 0.5
      tau_max <- if(use_bounds) input$harmonic_tau_max else NA

      if(use_bounds) {
        bounds_msg <- sprintf("Using parameter bounds: intercept [%.2f, %.2f], Amplitude [%.2f, %.2f]",
                              ifelse(is.na(mesor_min), -Inf, mesor_min),
                              ifelse(is.na(mesor_max), Inf, mesor_max),
                              amplitude_min,
                              ifelse(is.na(amplitude_max), Inf, amplitude_max))

        if(trend_type == "exp_sat") {
          bounds_msg <- paste0(bounds_msg,
                               sprintf(", A_sat [%.2f, %.2f], τ [%.2f, %.2f]",
                                       ifelse(is.na(A_sat_min), -Inf, A_sat_min),
                                       ifelse(is.na(A_sat_max), Inf, A_sat_max),
                                       tau_min,
                                       ifelse(is.na(tau_max), Inf, tau_max)))
        }

        showNotification(bounds_msg, type = "message", duration = 5)
      }

      withProgress(message = 'Fitting individual cosinor models...', value = 0, {
        for(i in 1:n_subjects) {
          y_i <- Y[i, ]

          # Count valid (non-NA) data points for this subject
          n_valid_points <- sum(!is.na(y_i))

          # Debug: Check for unusual values
          if(n_valid_points == 0) {
            cat(sprintf("Subject %d: All NA. First 5 values: %s\n", i,
                        paste(head(y_i, 5), collapse=", ")))
          }

          fit_i <- fit_cosinor(time_vec_model, y_i, period = period, n_harmonics = n_harmonics,
                               trend_type = trend_type,
                               use_bounds = use_bounds,
                               mesor_min = mesor_min,
                               mesor_max = mesor_max,
                               amplitude_min = amplitude_min,
                               amplitude_max = amplitude_max,
                               A_sat_min = A_sat_min,
                               A_sat_max = A_sat_max,
                               tau_min = tau_min,
                               tau_max = tau_max)

          if(fit_i$success) {
            individual_fits[[i]] <- fit_i

            # Build row with all harmonic parameters
            row_data <- list(subject = i, mesor = fit_i$mesor, mesor_se = fit_i$mesor_se)

            # the rhythm-adjusted mean and the fitted value at the first
            # observation, per subject, from this subject's own coefficients
            .tc <- switch(as.character(trend_type),
                          "linear"  = c(fit_i$trend_params$trend_linear$coef %||% NA_real_),
                          "log"     = c(fit_i$trend_params$trend_log$coef %||% NA_real_),
                          "exp_sat" = c(fit_i$trend_params$A_sat$coef %||% NA_real_,
                                        fit_i$trend_params$tau$coef %||% NA_real_),
                          numeric(0))
            row_data$mesor_adj <- fck_rhythm_adjusted_mean(
              fit_i$mesor, trend_type, .tc,
              min(time_vec_model, na.rm = TRUE), max(time_vec_model, na.rm = TRUE),
              fit_i$t_offset %||% 0)
            row_data$value_at_start <- tryCatch(
              as.numeric(fck_rhythm_from_coefs(
                fit_i$coefs, min(time_vec_model, na.rm = TRUE), period, n_harmonics,
                trend_type, include_trend = TRUE, t_offset = fit_i$t_offset %||% 0)),
              error = function(e) NA_real_)

            # Add trend parameters based on type
            if(trend_type != "none" && !is.null(fit_i$trend_params)) {
              for(param_name in names(fit_i$trend_params)) {
                row_data[[param_name]] <- fit_i$trend_params[[param_name]]$coef
                row_data[[paste0(param_name, "_se")]] <- fit_i$trend_params[[param_name]]$se
              }
            }

            for(h in 1:n_harmonics) {
              row_data[[paste0("amplitude_", h)]] <- fit_i$amplitudes[h]
              row_data[[paste0("amp_se_", h)]] <- fit_i$amp_se[h]
              row_data[[paste0("acrophase_rad_", h)]] <- fit_i$acrophases[h]
              row_data[[paste0("acrophase_time_", h)]] <- fit_i$acrophases_time[h]
              row_data[[paste0("acro_se_time_", h)]] <- fit_i$acro_se_time[h]

              # Get beta coefficients - handle both linear and nls fits
              if(trend_type == "exp_sat") {
                row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[paste0("b_cos", h)]
                row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[paste0("b_sin", h)]
              } else {
                cos_idx <- coef_offset + 2 * (h - 1) + 1
                sin_idx <- coef_offset + 2 * (h - 1) + 2
                row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[cos_idx]
                row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[sin_idx]
              }
            }
            row_data$r_squared <- fit_i$r_squared
            row_data$percent_rhythm <- fit_i$percent_rhythm
            row_data$p_value <- fit_i$p_value
            row_data$aic <- fit_i$aic
            row_data$aicc <- fit_i$aicc
            row_data$bic <- fit_i$bic
            row_data$loocv_rmse <- fit_i$loocv_rmse
            row_data$r_squared_S <- fit_i$r_squared_S
            row_data$r_squared_C <- fit_i$r_squared_C
            row_data$unique_S <- fit_i$unique_S %||% NA_real_
            row_data$unique_C <- fit_i$unique_C %||% NA_real_
            row_data$shared_SC <- fit_i$shared_SC %||% NA_real_
            row_data$percent_S <- fit_i$percent_S
            row_data$percent_C <- fit_i$percent_C
            row_data$percent_shared <- fit_i$percent_shared %||% NA_real_
            row_data$converged <- isTRUE(fit_i$converged)
            row_data$boundary_hit <- isTRUE(fit_i$boundary_hit)

            individual_params <- rbind(individual_params, as.data.frame(row_data))
          } else {
            # Store failed fit with reason
            individual_fits[[i]] <- list(
              success = FALSE, 
              message = fit_i$message,
              n_valid = n_valid_points,
              n_required = 2 * n_harmonics + 1 + n_trend_params + 1
            )
            failed_fits[[length(failed_fits) + 1]] <- list(
              subject = i,
              n_valid = n_valid_points,
              reason = fit_i$message
            )
          }
          
          if(i %% 10 == 0) incProgress(10 / n_subjects)
        }
      })
      
      # Report failed fits
      if(length(failed_fits) > 0) {
        n_failed <- length(failed_fits)
        min_required <- 2 * n_harmonics + 1 + n_trend_params + 1
        
        failed_subjects <- sapply(failed_fits, function(x) x$subject)
        failed_nvalid <- sapply(failed_fits, function(x) x$n_valid)
        
        msg <- sprintf("%d of %d subjects failed to fit (need %d+ valid points). Failed: %s",
                       n_failed, n_subjects, min_required,
                       paste(paste0("S", failed_subjects, "(", failed_nvalid, "pts)"), collapse = ", "))
        
        showNotification(msg, type = "warning", duration = 10)
      }
      
      # Population-mean statistics (always calculated)
      pop_mean_fit <- NULL
      group_fits <- NULL

      # ========================================================================
      # AUDIT 2.3: convergence gate
      #
      # "Successfully fitted: 1305 / 1305" alongside an R² range starting at
      # 0.060 was not a coincidence: nls with warnOnly = TRUE RETURNS a fit at
      # the iteration limit instead of erroring, so every non-converged
      # optimisation was counted as a success and averaged into the population
      # parameters. The fitters now report convergence honestly; here we act on
      # it. Non-converged fits and fits pinned to a parameter bound are excluded
      # from every population summary, and the counts are carried into the
      # report so the exclusion is visible rather than silent.
      # ========================================================================
      all_params <- individual_params
      conv_flag  <- if("converged" %in% names(all_params)) all_params$converged else rep(TRUE, nrow(all_params))
      bound_flag <- if("boundary_hit" %in% names(all_params)) all_params$boundary_hit else rep(FALSE, nrow(all_params))
      conv_flag[is.na(conv_flag)] <- FALSE
      bound_flag[is.na(bound_flag)] <- FALSE

      # Which bounds each returned fit sits on, per parameter.
      bounds_list <- lapply(seq_len(nrow(all_params)), function(i) {
        f <- individual_fits[[all_params$subject[i]]]
        if (is.null(f) || is.null(f$bounds_hit)) character(0) else f$bounds_hit
      })
      bounds_summary <- fck_bounds_summary(
        bounds_list,
        subject_ids = if (!is.null(values$subject_ids))
          values$subject_ids[all_params$subject] else all_params$subject)

      fit_audit <- list(
        n_attempted = n_subjects,
        n_returned  = nrow(all_params),
        n_converged = sum(conv_flag & !bound_flag),
        n_boundary  = sum(conv_flag & bound_flag),
        n_failed    = n_subjects - nrow(all_params),
        n_nonconverged = sum(!conv_flag),
        bounds = bounds_summary
      )

      # ======================================================================
      # WHO ENTERS THE POPULATION SUMMARIES
      #
      # Non-converged fits are always excluded: the optimiser stopped without
      # finding a solution, so there is no estimate to average.
      #
      # Fits pinned to a BOUND are a judgement call, and it is the user's, not
      # this code's. The value is real -- the optimiser did converge to it --
      # but it is the edge of the feasible region rather than an interior
      # optimum, so its standard error is meaningless and averaging it pulls the
      # mean toward whatever the bound happens to be. Excluding them makes the
      # summary cleaner and the sample smaller and possibly biased; including
      # them keeps everyone and lets the bound speak through the mean.
      #
      # Default is to INCLUDE, with the bound table shown, because a silently
      # shrunken sample is the worse failure. The table names which bound each
      # fit hit and which fits hit more than one, so the cost of including them
      # is visible rather than assumed.
      # ======================================================================
      include_boundary <- isTRUE(input$harmonic_include_boundary %||% TRUE)
      keep_rows <- conv_flag & (include_boundary | !bound_flag)
      fit_audit$include_boundary <- include_boundary

      if(sum(keep_rows) < 3) {
        # Refusing to summarise 2 subjects is better than summarising 1305 bad
        # ones, but refusing to summarise ANYTHING would be worse. Fall back,
        # and say so.
        showNotification(
          sprintf("Only %d of %d fits pass the current gate; population summaries fall back to all returned fits. Treat them as provisional.",
                  sum(keep_rows), nrow(all_params)),
          type = "warning", duration = 15)
        fit_audit$gate_relaxed <- TRUE
        keep_rows <- rep(TRUE, nrow(all_params))
      } else {
        fit_audit$gate_relaxed <- FALSE
      }
      fit_audit$n_summarised <- sum(keep_rows)
      individual_params <- all_params[keep_rows, , drop = FALSE]
      fit_audit$bounds_kept <- fck_bounds_summary(
        bounds_list[keep_rows],
        subject_ids = if (!is.null(values$subject_ids))
          values$subject_ids[all_params$subject[keep_rows]] else all_params$subject[keep_rows])

      # Always calculate population mean parameters (vector averaging for circular data)
      {
        # Calculate population mean parameters (vector averaging for circular data)
        mean_mesor <- mean(individual_params$mesor, na.rm = TRUE)
        
        # Vector average for amplitude and acrophase (first harmonic for primary stats)
        x_components <- individual_params$amplitude_1 * cos(individual_params$acrophase_rad_1)
        y_components <- individual_params$amplitude_1 * sin(individual_params$acrophase_rad_1)
        
        mean_x <- mean(x_components, na.rm = TRUE)
        mean_y <- mean(y_components, na.rm = TRUE)
        
        mean_amplitude <- sqrt(mean_x^2 + mean_y^2)
        mean_acrophase_rad <- atan2(mean_y, mean_x)
        if(mean_acrophase_rad < 0) mean_acrophase_rad <- mean_acrophase_rad + 2 * pi
        # MODEL-elapsed hours. phi_to_hours() is the single conversion; the
        # clock origin is added only at display, by fck_acrophase_label().
        mean_acrophase_time <- phi_to_hours(mean_acrophase_rad, period, 1)
        
        # ====================================================================
        # AUDIT 1.2: the Rayleigh test runs on UNIT vectors
        #
        # This block used to compute
        #     r_bar <- mean_amplitude / mean(amplitude_1)
        # which is |Sum A e^{i phi}| / Sum A -- the AMPLITUDE-WEIGHTED resultant
        # (0.824 in the reported output) -- and fed it to Z = n * r^2, giving
        # Z = 886.5. The Rayleigh test is defined on unit vectors (Mardia & Jupp
        # 2000; Berens 2009), whose resultant was 0.789 and whose Z is 812.3.
        # The same report printed both numbers, in different blocks, for the
        # same acrophases.
        #
        # fck_resultants() now returns both, under names that cannot be
        # confused, and only the unweighted one reaches fck_rayleigh(). The
        # amplitude-weighted vector mean stays the population estimator: that
        # part was always correct.
        # ====================================================================
        res1 <- fck_resultants(individual_params$acrophase_rad_1,
                               individual_params$amplitude_1)
        n_valid <- res1$n
        r_bar_unweighted <- res1$r_unweighted
        r_bar_weighted   <- res1$r_weighted
        ray <- fck_rayleigh(r_bar_unweighted, n_valid)
        rayleigh_z <- ray$Z
        rayleigh_p <- ray$p

        # circular dispersion is defined on the unweighted resultant
        r_bar <- r_bar_unweighted        # kept for downstream compatibility
        circ_var <- 1 - r_bar_unweighted
        circ_sd <- res1$circ_sd_rad

        # per-harmonic resultants, for the report
        resultants <- lapply(seq_len(n_harmonics), function(h) {
          rr <- fck_resultants(individual_params[[paste0("acrophase_rad_", h)]],
                               individual_params[[paste0("amplitude_", h)]])
          if(is.null(rr)) return(NULL)
          rr$rayleigh <- fck_rayleigh(rr$r_unweighted, rr$n)
          rr$harmonic <- h
          rr$circ_sd_hours <- phi_to_hours(rr$circ_sd_rad, period, h)
          rr
        })
        
        # Store mean coefficients for ALL harmonics (for proper multi-harmonic curve plotting)
        # Format: [mesor, (trend coefs if trend), beta_cos_1, beta_sin_1, ...]
        mean_coefs <- c(mean_mesor)
        
        # Add mean trend coefficient(s) based on trend type
        trend_coefs <- get_mean_trend_coefs(individual_params, trend_type)
        if(length(trend_coefs) > 0) {
          mean_coefs <- c(mean_coefs, trend_coefs)
        }
        
        mean_amplitudes <- numeric(n_harmonics)
        mean_acrophases_rad <- numeric(n_harmonics)
        mean_acrophases_time <- numeric(n_harmonics)
        
        for(h in 1:n_harmonics) {
          beta_cos_col <- paste0("beta_cos_", h)
          beta_sin_col <- paste0("beta_sin_", h)
          amp_col <- paste0("amplitude_", h)
          acro_col <- paste0("acrophase_rad_", h)
          
          # Mean of raw coefficients (for curve reconstruction)
          mean_beta_cos <- mean(individual_params[[beta_cos_col]], na.rm = TRUE)
          mean_beta_sin <- mean(individual_params[[beta_sin_col]], na.rm = TRUE)
          mean_coefs <- c(mean_coefs, mean_beta_cos, mean_beta_sin)
          
          # Vector-averaged amplitude and acrophase
          x_h <- individual_params[[amp_col]] * cos(individual_params[[acro_col]])
          y_h <- individual_params[[amp_col]] * sin(individual_params[[acro_col]])
          mean_amplitudes[h] <- sqrt(mean(x_h, na.rm = TRUE)^2 + mean(y_h, na.rm = TRUE)^2)
          acro_h <- atan2(mean(y_h, na.rm = TRUE), mean(x_h, na.rm = TRUE))
          if(acro_h < 0) acro_h <- acro_h + 2 * pi
          mean_acrophases_rad[h] <- acro_h
          mean_acrophases_time[h] <- phi_to_hours(acro_h, period, h)   # model-elapsed
        }
        
        # Also compute arithmetic means of individual parameters
        indiv_means <- list(
          mesor = mean(individual_params$mesor, na.rm = TRUE),
          mesor_sd = sd(individual_params$mesor, na.rm = TRUE)
        )

        # ====================================================================
        # AUDIT 1.1 (root cause): the trend parameters were NEVER put in here
        #
        # The pooled fitted-equation builder read pop$indiv_means$A_sat / $tau /
        # $trend_linear / $trend_log. None of those keys was ever created, so the
        # trend branch was dead for EVERY trend type and the pooled equation
        # printed the model without its homeostatic term -- under-predicting by
        # about 20 units everywhere -- while the header, the symbolic equation
        # and all four group equations included it.
        #
        # The group builder read a different structure (g$trend_params), which
        # is why only the pooled line was wrong. Both call sites now go through
        # fck_format_equation(); this fills the gap the pooled one was reading
        # from, so the two agree by construction rather than by coincidence.
        # ====================================================================
        for(tc in c("trend_linear", "trend_log", "A_sat", "tau")) {
          if(tc %in% names(individual_params)) {
            indiv_means[[tc]] <- mean(individual_params[[tc]], na.rm = TRUE)
            indiv_means[[paste0(tc, "_sd")]] <- sd(individual_params[[tc]], na.rm = TRUE)
          }
        }
        for(h in 1:n_harmonics) {
          indiv_means[[paste0("amplitude_", h)]] <- mean(individual_params[[paste0("amplitude_", h)]], na.rm = TRUE)
          indiv_means[[paste0("amplitude_", h, "_sd")]] <- sd(individual_params[[paste0("amplitude_", h)]], na.rm = TRUE)
          indiv_means[[paste0("acrophase_time_", h)]] <- mean(individual_params[[paste0("acrophase_time_", h)]], na.rm = TRUE)
          indiv_means[[paste0("acrophase_time_", h, "_sd")]] <- sd(individual_params[[paste0("acrophase_time_", h)]], na.rm = TRUE)
        }

        # Add model selection metrics (if columns exist)
        if("aic" %in% names(individual_params)) {
          indiv_means$aic <- mean(individual_params$aic, na.rm = TRUE)
          indiv_means$aic_sd <- sd(individual_params$aic, na.rm = TRUE)
          indiv_means$aicc <- mean(individual_params$aicc, na.rm = TRUE)
          indiv_means$aicc_sd <- sd(individual_params$aicc, na.rm = TRUE)
          indiv_means$bic <- mean(individual_params$bic, na.rm = TRUE)
          indiv_means$bic_sd <- sd(individual_params$bic, na.rm = TRUE)
          indiv_means$loocv_rmse <- mean(individual_params$loocv_rmse, na.rm = TRUE)
          indiv_means$loocv_rmse_sd <- sd(individual_params$loocv_rmse, na.rm = TRUE)
        }

        # Add variance decomposition statistics (if columns exist)
        if("r_squared_S" %in% names(individual_params) && "r_squared_C" %in% names(individual_params)) {
          indiv_means$r_squared_S <- mean(individual_params$r_squared_S, na.rm = TRUE)
          indiv_means$r_squared_S_sd <- sd(individual_params$r_squared_S, na.rm = TRUE)
          indiv_means$r_squared_C <- mean(individual_params$r_squared_C, na.rm = TRUE)
          indiv_means$r_squared_C_sd <- sd(individual_params$r_squared_C, na.rm = TRUE)
          indiv_means$percent_S <- mean(individual_params$percent_S, na.rm = TRUE)
          indiv_means$percent_S_sd <- sd(individual_params$percent_S, na.rm = TRUE)
          indiv_means$percent_C <- mean(individual_params$percent_C, na.rm = TRUE)
          indiv_means$percent_C_sd <- sd(individual_params$percent_C, na.rm = TRUE)
        }

        # AUDIT 1.3: the commonality partition, averaged across subjects. Because
        # each subject's three parts sum to that subject's total R², the means
        # also sum to the mean total R² -- which is the property the old
        # "proportions" (30.8% + 93.9% = 124.7%) could never have.
        for(cc in c("unique_S", "unique_C", "shared_SC",
                    "percent_S", "percent_C", "percent_shared")) {
          if(cc %in% names(individual_params)) {
            indiv_means[[cc]] <- mean(individual_params[[cc]], na.rm = TRUE)
            indiv_means[[paste0(cc, "_sd")]] <- sd(individual_params[[cc]], na.rm = TRUE)
          }
        }
        indiv_means$r_squared <- mean(individual_params$r_squared, na.rm = TRUE)
        
        # ====================================================================
        # AUDIT 1.4: what the report called the MESOR is not a MESOR
        #
        # mean_mesor is the arithmetic mean of the fitted CONSTANT term. In a
        # model with a trend that constant is not the rhythm-adjusted mean: it
        # is the intercept of a model whose trend is anchored at t_offset (the
        # first observation, because fit_cosinor_nonlinear builds the trend on
        # t - min(t)) while the harmonics are anchored at t = 0. Two origins,
        # one constant. Calling it "MESOR" made a 27.70 look like a central
        # value when the data over the window average near 44.
        #
        # The MESOR proper (Cornelissen 2014) is the rhythm-adjusted mean: the
        # time-average of the non-oscillating part across the observation
        # window. Computed by integration below, and reported under that name;
        # the constant is reported as the intercept, under ITS name.
        # ====================================================================
        t_lo <- min(time_vec, na.rm = TRUE); t_hi <- max(time_vec, na.rm = TRUE)
        pop_trend_coefs <- switch(as.character(trend_type),
          "linear"  = c(indiv_means$trend_linear %||% NA_real_),
          "log"     = c(indiv_means$trend_log %||% NA_real_),
          "exp_sat" = c(indiv_means$A_sat %||% NA_real_, indiv_means$tau %||% NA_real_),
          numeric(0))
        rhythm_adjusted_mean <- fck_rhythm_adjusted_mean(
          mean_mesor, trend_type, pop_trend_coefs, t_lo, t_hi, t_offset_global)
        # how much of the oscillation leaks into the window mean, which is only
        # exactly zero over a whole number of periods (here 22 h of a 24 h cycle)
        harmonic_leak <- fck_harmonic_window_mean(
          mean_coefs[seq(2 + length(pop_trend_coefs), length(mean_coefs), by = 2)],
          mean_coefs[seq(3 + length(pop_trend_coefs), length(mean_coefs), by = 2)],
          period, t_lo, t_hi)

        pop_mean_fit <- list(
          mean_mesor = mean_mesor,
          intercept = mean_mesor,                 # its correct name
          rhythm_adjusted_mean = rhythm_adjusted_mean,
          harmonic_window_mean = harmonic_leak,
          window = c(t_lo, t_hi),
          t_offset = t_offset_global,
          trend_coefs = pop_trend_coefs,
          resultants = resultants,
          r_bar_unweighted = r_bar_unweighted,
          r_bar_weighted = r_bar_weighted,
          fit_audit = fit_audit,
          mean_amplitude = mean_amplitude,  # First harmonic (for backwards compatibility)
          mean_acrophase_rad = mean_acrophase_rad,
          mean_acrophase_time = mean_acrophase_time,
          mean_coefs = mean_coefs,  # All coefficients for curve reconstruction
          mean_amplitudes = mean_amplitudes,  # All harmonics
          mean_acrophases_rad = mean_acrophases_rad,
          mean_acrophases_time = mean_acrophases_time,
          indiv_means = indiv_means,  # Arithmetic means of individual params
          r_bar = r_bar,
          circ_var = circ_var,
          circ_sd = circ_sd,
          rayleigh_z = rayleigh_z,
          rayleigh_p = rayleigh_p,
          n = n_valid
        )
      }
      
      # Group comparison - calculate group-specific statistics
      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]

        # ====================================================================
        # AUDIT 1.5: 654 + 410 + 181 + 59 = 1304, not 1305
        #
        # The old code did
        #     groups <- unique(group_var)              # keeps NA as a level
        #     idx    <- which(group_var == g)          # NA == NA is NA -> empty
        #     if (nrow(grp_params) >= 3) { ... }       # so the level vanishes
        # A subject whose group label was missing or unmatched therefore entered
        # every pooled statistic and no group, with no message anywhere. Any
        # genuine group smaller than the n >= 3 guard disappeared the same way.
        #
        # Now the accounting is explicit and asserted: unassigned subjects get
        # their own UNASSIGNED row, and groups dropped for being too small are
        # named. sum(group_n) == n_fitted is checked, and a mismatch is a loud
        # warning rather than a silent subtraction.
        # ====================================================================
        group_audit <- fck_group_audit(group_var, individual_params$subject, min_n = 3)
        if(group_audit$n_unassigned > 0) {
          showNotification(
            sprintf("%d of %d fitted subject%s no usable '%s' label. They are pooled but not grouped, and appear as UNASSIGNED in the report.",
                    group_audit$n_unassigned, group_audit$n_total,
                    if(group_audit$n_unassigned == 1) " has" else "s have",
                    input$harmonic_group_var),
            type = "warning", duration = 15)
        }
        if(length(group_audit$dropped_small) > 0) {
          showNotification(
            sprintf("Group(s) %s have fewer than 3 fitted subjects and are not summarised separately.",
                    paste(group_audit$dropped_small, collapse = ", ")),
            type = "warning", duration = 15)
        }

        lab_all <- as.character(group_var)
        groups <- group_audit$levels
        group_fits <- list()

        # A subject with no usable group label is EXCLUDED from every group
        # analysis. An earlier version carried them as an "UNASSIGNED" row so
        # they could not disappear silently -- but a label-less group of one is
        # not a group: it has no circular mean, and every comparison built on it
        # produced NaN, which is what crashed the group-comparison plot
        # (circular_mean of an empty vector -> NaN -> `if (NaN < 0)`).
        #
        # Visibility is kept where it belongs: the count of excluded subjects is
        # carried in the audit and printed by the report and the comparison
        # panel, so the number still reconciles against n fitted. They are named,
        # not analysed.
        for(g in groups) {
          idx <- which(!is.na(lab_all) & nzchar(lab_all) & lab_all == g)
          grp_params <- individual_params[individual_params$subject %in% idx, ]

          if(nrow(grp_params) >= 3) {
            # Mean coefficients for curve reconstruction
            grp_coefs <- c(mean(grp_params$mesor, na.rm = TRUE))
            
            # Add trend coefficients based on type
            grp_trend_params <- list()
            if(trend_type == "linear" && "trend_linear" %in% names(grp_params)) {
              grp_coefs <- c(grp_coefs, mean(grp_params$trend_linear, na.rm = TRUE))
              grp_trend_params$trend_linear <- list(
                mean = mean(grp_params$trend_linear, na.rm = TRUE),
                sd = sd(grp_params$trend_linear, na.rm = TRUE)
              )
            } else if(trend_type == "log" && "trend_log" %in% names(grp_params)) {
              grp_coefs <- c(grp_coefs, mean(grp_params$trend_log, na.rm = TRUE))
              grp_trend_params$trend_log <- list(
                mean = mean(grp_params$trend_log, na.rm = TRUE),
                sd = sd(grp_params$trend_log, na.rm = TRUE)
              )
            } else if(trend_type == "exp_sat") {
              if("A_sat" %in% names(grp_params)) {
                A_sat_mean <- mean(grp_params$A_sat, na.rm = TRUE)
                # If all NA, use 0 as fallback
                if(!is.finite(A_sat_mean)) A_sat_mean <- 0
                grp_coefs <- c(grp_coefs, A_sat_mean)
                grp_trend_params$A_sat <- list(
                  mean = A_sat_mean,
                  sd = sd(grp_params$A_sat, na.rm = TRUE)
                )
              }
              if("tau" %in% names(grp_params)) {
                tau_mean <- mean(grp_params$tau, na.rm = TRUE)
                # If all NA, use 1 as fallback (avoid division by zero)
                if(!is.finite(tau_mean) || tau_mean <= 0) tau_mean <- 1
                grp_coefs <- c(grp_coefs, tau_mean)
                grp_trend_params$tau <- list(
                  mean = tau_mean,
                  sd = sd(grp_params$tau, na.rm = TRUE)
                )
              }
            }
            
            grp_amplitudes <- numeric(n_harmonics)
            grp_acrophases_rad <- numeric(n_harmonics)
            grp_acrophases_time <- numeric(n_harmonics)
            grp_amp_sd <- numeric(n_harmonics)
            
            for(h in 1:n_harmonics) {
              beta_cos_col <- paste0("beta_cos_", h)
              beta_sin_col <- paste0("beta_sin_", h)
              amp_col <- paste0("amplitude_", h)
              acro_col <- paste0("acrophase_rad_", h)

              grp_coefs <- c(grp_coefs,
                             mean(grp_params[[beta_cos_col]], na.rm = TRUE),
                             mean(grp_params[[beta_sin_col]], na.rm = TRUE))

              # Vector-averaged amplitude/acrophase
              x_h <- grp_params[[amp_col]] * cos(grp_params[[acro_col]])
              y_h <- grp_params[[amp_col]] * sin(grp_params[[acro_col]])
              grp_amplitudes[h] <- sqrt(mean(x_h, na.rm = TRUE)^2 + mean(y_h, na.rm = TRUE)^2)
              grp_amp_sd[h] <- sd(grp_params[[amp_col]], na.rm = TRUE)
              acro_h <- atan2(mean(y_h, na.rm = TRUE), mean(x_h, na.rm = TRUE))
              if(acro_h < 0) acro_h <- acro_h + 2 * pi
              grp_acrophases_rad[h] <- acro_h
              grp_acrophases_time[h] <- phi_to_hours(acro_h, period, h)  # model-elapsed
            }

            # Variance decomposition for this group (if columns exist)
            grp_variance_decomp <- NULL
            if("r_squared_S" %in% names(grp_params) && "r_squared_C" %in% names(grp_params)) {
              grp_variance_decomp <- list(
                r_squared_S = mean(grp_params$r_squared_S, na.rm = TRUE),
                r_squared_S_sd = sd(grp_params$r_squared_S, na.rm = TRUE),
                r_squared_C = mean(grp_params$r_squared_C, na.rm = TRUE),
                r_squared_C_sd = sd(grp_params$r_squared_C, na.rm = TRUE),
                percent_S = mean(grp_params$percent_S, na.rm = TRUE),
                percent_S_sd = sd(grp_params$percent_S, na.rm = TRUE),
                percent_C = mean(grp_params$percent_C, na.rm = TRUE),
                percent_C_sd = sd(grp_params$percent_C, na.rm = TRUE)
              )
            }

            # AUDIT 1.4: per group, the same distinction as the pooled block --
            # the fitted constant is the intercept; the MESOR is the
            # rhythm-adjusted mean over the observation window.
            grp_trend_coefs <- switch(as.character(trend_type),
              "linear"  = c(grp_trend_params$trend_linear$mean %||% NA_real_),
              "log"     = c(grp_trend_params$trend_log$mean %||% NA_real_),
              "exp_sat" = c(grp_trend_params$A_sat$mean %||% NA_real_,
                            grp_trend_params$tau$mean %||% NA_real_),
              numeric(0))
            grp_intercept <- mean(grp_params$mesor, na.rm = TRUE)
            grp_ram <- fck_rhythm_adjusted_mean(
              grp_intercept, trend_type, grp_trend_coefs,
              min(time_vec, na.rm = TRUE), max(time_vec, na.rm = TRUE),
              t_offset_global)

            # AUDIT 1.6: within a group, the intercept / A_sat / tau are
            # ARITHMETIC means while the amplitudes and acrophases are VECTOR
            # means -- and an SD was printed next to the vector-averaged
            # amplitude, implying arithmetic averaging of a quantity that had
            # not been averaged arithmetically. Both summaries are now carried
            # so the report can label each line with the estimator that produced
            # it, and the arithmetic amplitude mean sits next to its own SD.
            grp_amp_arith <- vapply(seq_len(n_harmonics), function(h)
              mean(grp_params[[paste0("amplitude_", h)]], na.rm = TRUE), numeric(1))
            grp_res <- lapply(seq_len(n_harmonics), function(h) {
              rr <- fck_resultants(grp_params[[paste0("acrophase_rad_", h)]],
                                   grp_params[[paste0("amplitude_", h)]])
              if(is.null(rr)) return(NULL)
              rr$rayleigh <- fck_rayleigh(rr$r_unweighted, rr$n)
              rr$circ_sd_hours <- phi_to_hours(rr$circ_sd_rad, period, h)
              rr
            })

            group_fits[[as.character(g)]] <- list(
              group = g,
              is_unassigned = FALSE,
              n = nrow(grp_params),
              mean_mesor = grp_intercept,             # legacy name, kept working
              intercept = grp_intercept,              # its correct name
              rhythm_adjusted_mean = grp_ram,
              trend_coefs = grp_trend_coefs,
              amp_arithmetic = grp_amp_arith,
              resultants = grp_res,
              sd_mesor = sd(grp_params$mesor, na.rm = TRUE),
              mean_coefs = grp_coefs,  # All coefficients
              trend_params = grp_trend_params,
              mean_amplitudes = grp_amplitudes,
              sd_amplitudes = grp_amp_sd,
              mean_acrophases_rad = grp_acrophases_rad,
              mean_acrophases_time = grp_acrophases_time,
              # Keep first harmonic for backwards compatibility
              mean_amplitude = grp_amplitudes[1],
              sd_amplitude = grp_amp_sd[1],
              mean_acrophase_rad = grp_acrophases_rad[1],
              mean_acrophase_time = grp_acrophases_time[1],
              # Variance decomposition
              variance_decomp = grp_variance_decomp
            )
          }
        }

        # AUDIT 1.5: the assertion the brief asked for. It cannot fire now that
        # UNASSIGNED is a real row, which is exactly the point -- if it ever
        # does, something new is dropping subjects and the report says so
        # instead of printing group sizes that do not add up.
        .grp_total <- sum(vapply(group_fits, function(g) g$n, integer(1)))
        .excluded <- nrow(individual_params) - .grp_total
        if(length(group_fits) && .excluded > 0) {
          showNotification(
            sprintf("%d fitted subject%s excluded from the group analyses: %d with no usable '%s' label%s. They remain in every pooled statistic.",
                    .excluded, if(.excluded == 1) " is" else "s are",
                    group_audit$n_unassigned, input$harmonic_group_var,
                    if(length(group_audit$dropped_small))
                      sprintf(", and %d in group(s) with fewer than 3 fits (%s)",
                              group_audit$n_dropped_small,
                              paste(group_audit$dropped_small, collapse = ", ")) else ""),
            type = "warning", duration = 15)
        }
        attr(group_fits, "audit") <- group_audit
        attr(group_fits, "n_fitted") <- nrow(individual_params)
        attr(group_fits, "n_in_groups") <- .grp_total
      }
      
      # Bootstrap CIs if requested
      boot_results <- NULL
      if(isTRUE(input$harmonic_bootstrap)) {
        B <- input$harmonic_n_boot
        boot_mesor <- numeric(B)
        boot_amplitude <- numeric(B)
        boot_acrophase <- numeric(B)
        
        showNotification(paste("Running", B, "bootstrap iterations..."), type = "message")
        
        withProgress(message = 'Bootstrap...', value = 0, {
          for(b in 1:B) {
            boot_idx <- sample(1:n_subjects, n_subjects, replace = TRUE)
            boot_params <- individual_params[individual_params$subject %in% boot_idx, ]
            
            boot_mesor[b] <- mean(boot_params$mesor, na.rm = TRUE)
            
            # Use first harmonic for bootstrap CIs
            x_b <- boot_params$amplitude_1 * cos(boot_params$acrophase_rad_1)
            y_b <- boot_params$amplitude_1 * sin(boot_params$acrophase_rad_1)
            boot_amplitude[b] <- sqrt(mean(x_b, na.rm = TRUE)^2 + mean(y_b, na.rm = TRUE)^2)
            
            acro_b <- atan2(mean(y_b, na.rm = TRUE), mean(x_b, na.rm = TRUE))
            if(acro_b < 0) acro_b <- acro_b + 2 * pi
            boot_acrophase[b] <- acro_b
            
            if(b %% 50 == 0) incProgress(50 / B)
          }
        })
        
        boot_results <- list(
          mesor_ci = quantile(boot_mesor, c(0.025, 0.975)),
          amplitude_ci = quantile(boot_amplitude, c(0.025, 0.975)),
          acrophase_ci = quantile(phi_to_hours(boot_acrophase, period, 1), c(0.025, 0.975)),
          boot_mesor = boot_mesor,
          boot_amplitude = boot_amplitude,
          boot_acrophase = boot_acrophase,
          B = B
        )
      }
      
      # ======================================================================
      # AUDIT 2.4: information criteria against a NESTED MODEL SET
      #
      # AIC/AICc/BIC printed as means with SDs across subjects carry no
      # information: with no competing model they are constant offsets of one
      # another, which is exactly why the three SDs printed identically at
      # 16.21. What is interpretable is Delta-AICc across a nested set, with
      # Akaike weights.
      #
      # The set is the one the brief asked for: trend in {none, linear,
      # exp_sat} crossed with harmonics in {1, 2, 3}. Each cell is fitted on the
      # SAME subjects with the SAME time vector, and the AICc's are summed
      # across subjects (equivalently: this is the joint model over independent
      # subjects, which is what "which specification should I report" asks).
      #
      # Guarded behind a flag: 9 cells x n subjects x a nonlinear fit each is
      # minutes of compute on 1305 subjects, and it should not run unasked.
      # ======================================================================
      model_selection <- NULL
      if(isTRUE(input$harmonic_model_selection)) {
        ms_trends <- c("none", "linear", "exp_sat")
        ms_harms  <- 1:3
        ms_rows <- list()
        withProgress(message = "Fitting the nested model set...", value = 0, {
          for(tt in ms_trends) for(hh in ms_harms) {
            incProgress(1 / (length(ms_trends) * length(ms_harms)),
                        detail = sprintf("trend=%s, H=%d", tt, hh))
            need <- 2 * hh + 1 + switch(tt, "none" = 0, "linear" = 1, "exp_sat" = 2, 0) + 1
            if(n_time < need) next
            aicc_sum <- 0; k_ok <- 0L
            for(i in seq_len(n_subjects)) {
              y_i <- as.numeric(Y[i, ])
              if(sum(!is.na(y_i)) < need) next
              f <- tryCatch(fit_cosinor(time_vec_model, y_i, period = period,
                                        n_harmonics = hh, trend_type = tt),
                            error = function(e) NULL)
              if(is.null(f) || !isTRUE(f$success) || !is.finite(f$aicc)) next
              aicc_sum <- aicc_sum + f$aicc; k_ok <- k_ok + 1L
            }
            if(k_ok > 0) ms_rows[[sprintf("%s + H%d", tt, hh)]] <- aicc_sum / k_ok * 1
          }
        })
        model_selection <- fck_akaike_table(ms_rows)
        if(!is.null(model_selection)) attr(model_selection, "per_subject_mean") <- TRUE
      }

      # ======================================================================
      # AUDIT 2.2: A_sat and tau are weakly identified. Document the ridge.
      #
      # Over t in [8, 30] with tau ~ 13.9 the factor (1 - e^(-t/tau)) moves only
      # from 0.44 to 0.88 and is close to linear there, so it is strongly
      # collinear with the intercept and partly with the 24 h cosine. The
      # evidence was already in the output and unremarked: SD(tau) = 11.58 on a
      # mean of 13.92, SD(A_sat) = 22.05 on 32.30. That is a likelihood ridge,
      # not population heterogeneity, and the difference matters because the
      # second reading invites a between-group comparison of tau that the first
      # forbids.
      #
      # Three pieces of evidence are computed here:
      #   1. the mean within-subject parameter correlation matrix
      #   2. the design-matrix condition number at each time origin
      #   3. Delta-AIC of free-tau against tau fixed at a literature value
      #      (Daan, Beersma & Borbely 1984 give tau_rise ~ 18 h under extended
      #      wakefulness)
      # ======================================================================
      conditioning <- NULL
      if(trend_type == "exp_sat" && nrow(individual_params) > 0) {
        conditioning <- list()

        cors <- list()
        for(i in seq_len(min(nrow(individual_params), 200))) {
          sid <- individual_params$subject[i]
          fi <- individual_fits[[sid]]
          if(is.null(fi) || is.null(fi$vcov_full)) next
          V <- fi$vcov_full
          d <- sqrt(diag(V))
          if(any(!is.finite(d)) || any(d <= 0)) next
          cors[[length(cors) + 1]] <- V / outer(d, d)
        }
        if(length(cors) > 0) {
          conditioning$mean_cor <- Reduce(`+`, cors) / length(cors)
          conditioning$n_cor <- length(cors)
        }

        # condition number of the linearised design at each origin
        kappa_at <- function(tv, tau_hat) {
          X <- cbind(1, 1 - exp(-(tv - min(tv)) / tau_hat))
          for(h in seq_len(n_harmonics)) {
            w <- 2 * pi * h / period
            X <- cbind(X, cos(w * tv), sin(w * tv))
          }
          sv <- svd(X)$d
          if(min(sv) <= 0) Inf else max(sv) / min(sv)
        }
        tau_hat <- mean(individual_params$tau, na.rm = TRUE)
        if(is.finite(tau_hat) && tau_hat > 0) {
          conditioning$kappa_before <- tryCatch(kappa_at(time_vec, tau_hat), error = function(e) NA_real_)
          conditioning$kappa_after  <- tryCatch(kappa_at(time_vec - min(time_vec), tau_hat),
                                                error = function(e) NA_real_)
        }

        # tau held at a literature value: is the extra parameter earning its keep?
        tau_fix <- suppressWarnings(as.numeric(input$harmonic_tau_fixed %||% 18))
        if(is.finite(tau_fix) && tau_fix > 0) {
          aic_free <- mean(individual_params$aic, na.rm = TRUE)
          aic_fix <- NA_real_
          acc <- c(); nok <- 0L
          for(i in seq_len(n_subjects)) {
            y_i <- as.numeric(Y[i, ]); ok <- !is.na(y_i)
            if(sum(ok) < 2 * n_harmonics + 2 + 1) next
            tv <- time_vec_model[ok]; yv <- y_i[ok]
            X <- cbind(1, 1 - exp(-(tv - min(time_vec_model)) / tau_fix))
            for(h in seq_len(n_harmonics)) {
              w <- 2 * pi * h / period
              X <- cbind(X, cos(w * tv), sin(w * tv))
            }
            fitf <- tryCatch(lm.fit(X, yv), error = function(e) NULL)
            if(is.null(fitf)) next
            nn <- length(yv); ssr <- sum(fitf$residuals^2)
            if(!is.finite(ssr) || ssr <= 0) next
            kk <- ncol(X) + 1                       # one fewer than free-tau
            ll <- -nn/2 * (log(2*pi) + log(ssr/nn) + 1)
            acc <- c(acc, -2 * ll + 2 * kk); nok <- nok + 1L
          }
          if(nok > 0) aic_fix <- mean(acc, na.rm = TRUE)
          if(is.finite(aic_free) && is.finite(aic_fix)) {
            conditioning$tau_fixed_value <- tau_fix
            conditioning$tau_fixed_delta_aic <- aic_free - aic_fix
            conditioning$tau_fixed_n <- nok
          }
        }
      }

      # ======================================================================
      # AUDIT 2.5: summarise the per-subject Bingham regions
      # ======================================================================
      bingham_summary <- NULL
      if(nrow(individual_params) > 0) {
        bingham_summary <- lapply(seq_len(n_harmonics), function(h) {
          amps <- c(); acros <- c(); n_id <- 0L; n_tot <- 0L
          for(i in seq_len(nrow(individual_params))) {
            fi <- individual_fits[[individual_params$subject[i]]]
            b <- if(!is.null(fi$bingham)) fi$bingham[[h]] else NULL
            if(is.null(b)) next
            n_tot <- n_tot + 1L
            amps <- c(amps, diff(b$amplitude) / 2)
            if(isTRUE(b$identified)) {
              n_id <- n_id + 1L
              d <- abs(((diff(b$acrophase_rad) + pi) %% (2 * pi)) - pi) / 2
              acros <- c(acros, phi_to_hours(d, period, h))
            }
          }
          if(n_tot == 0) return(NULL)
          list(n = n_tot, n_identified = n_id,
               median_amp_halfwidth = stats::median(amps, na.rm = TRUE),
               median_acro_halfwidth_h = if(length(acros)) stats::median(acros, na.rm = TRUE) else NA_real_)
        })
      }

      # Store results
      values$harmonic_model <- list(
        conditioning = conditioning,
        bingham_summary = bingham_summary,
        loocv_is_true_cv = {
          f1 <- individual_fits[[individual_params$subject[1]]]
          if(!is.null(f1$loocv_is_true_cv)) isTRUE(f1$loocv_is_true_cv) else TRUE
        },
        individual_fits = individual_fits,
        individual_params = individual_params,
        all_individual_params = all_params,   # AUDIT 2.3: before the convergence gate
        fit_audit = fit_audit,
        model_selection = model_selection,
        # AUDIT 1.7: the dependent variable was never named anywhere
        dv_name  = if(nzchar(input$harmonic_dv_name %||% "")) input$harmonic_dv_name else NULL,
        dv_units = if(nzchar(input$harmonic_dv_units %||% "")) input$harmonic_dv_units else NULL,
        dv_min   = suppressWarnings(as.numeric(input$harmonic_dv_min %||% NA)),
        dv_max   = suppressWarnings(as.numeric(input$harmonic_dv_max %||% NA)),
        data_source = data_source,
        time_origin = time_origin,
        origin_shift = origin_shift,
        time_vec_clock = time_vec,           # display axis, always clock-linearised
        group_var_name = if(!is.null(input$harmonic_group_var) &&
                            input$harmonic_group_var != "_none_") input$harmonic_group_var else NULL,
        pop_mean_fit = pop_mean_fit,
        group_fits = group_fits,
        boot_results = boot_results,
        time_vec = time_vec_model,
        original_times = original_times,
        wrap_applied = wrap_applied,
        period = period,
        n_harmonics = n_harmonics,
        trend_type = trend_type,
        include_trend = trend_type != "none",  # For backwards compatibility
        # MERGED APP: the bounds this run actually used, for the code export.
        fck_settings = list(
          use_bounds = use_bounds, mesor_min = mesor_min, mesor_max = mesor_max,
          amplitude_min = amplitude_min, amplitude_max = amplitude_max,
          A_sat_min = A_sat_min, A_sat_max = A_sat_max,
          tau_min = tau_min, tau_max = tau_max),
        t_offset = t_offset_global,
        t_center = t_center_global,
        using_smoothed = using_smoothed,
        subjects_with_nas = subjects_with_nas,
        Y = Y
      )
      
      showNotification("Harmonic regression complete!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      print(e)
    })
  })
  
  # ==============================================================================
  # HARMONIC REGRESSION OUTPUTS
  # ==============================================================================
  
  # The scrollable wrapper around the summary. The height is a slider rather
  # than a fixed value because how much of this you want on screen depends
  # entirely on whether you are reading it or working past it.
  output$harmonic_summary_box <- renderUI({
    h <- suppressWarnings(as.numeric(input$harmonic_summary_height %||% 600))
    if(!is.finite(h) || h < 100) h <- 600
    tagList(
      div(style = sprintf("max-height:%dpx; overflow-y:auto; overflow-x:auto; border:1px solid #e5e5e5; border-radius:3px; padding:6px;", as.integer(h)),
          verbatimTextOutput("harmonic_summary")),
      hr(),
      uiOutput("harmonic_parameters_table")
    )
  })

  # The same text as a file. Both this and the on-screen panel call ONE
  # function, .print_harmonic_summary(), so the file cannot drift from what you
  # read -- which is the whole reason to have it.
  output$download_harmonic_summary <- downloadHandler(
    filename = function() sprintf("harmonic_summary_%s.txt", Sys.Date()),
    content = function(file) {
      writeLines(utils::capture.output(.print_harmonic_summary()), file)
    }
  )

  # ============================================================================
  # Summary output
  #
  # CHANGELOG (2026-09-03 audit)
  #   1.1  the pooled fitted equation is built by fck_format_equation(), the
  #        same renderer the group equations use. The duplicate builder that
  #        silently dropped the homeostatic term is gone.
  #   1.2  both resultants are printed, each labelled with what it is for, and
  #        the Rayleigh Z comes from the unweighted one.
  #   1.3  commonality analysis replaces two overlapping marginal R2s. The
  #        auto-generated dominance verdict is deleted.
  #   1.4  the fitted constant is called the intercept. A genuine
  #        rhythm-adjusted mean is computed and called the MESOR.
  #   1.5  group sizes are reconciled against n fitted; UNASSIGNED is a row.
  #   1.6  every averaged quantity names the estimator that produced it.
  #   1.7  fmt2() throughout; DV named; H2 phase convention stated; the
  #        "(units)" placeholder is gone.
  #   2.1  the data source is stated and the significance rate is flagged.
  #   2.3  converged / boundary / failed are reported separately.
  #   2.4  Delta-AICc table with Akaike weights instead of mean AIC/AICc/BIC.
  #   2.5  Bingham joint confidence regions.
  # ============================================================================
  # The report body lives in a plain function so the on-screen panel and the
  # text download share it verbatim.
  .print_harmonic_summary <- function() {
    req(values$harmonic_model)
    mod <- values$harmonic_model
    period <- mod$period
    nh <- mod$n_harmonics
    trend_type <- mod$trend_type %||% "none"
    params <- mod$individual_params
    pop <- mod$pop_mean_fit

    dv <- mod$dv_name %||% "the dependent variable"
    dvu <- mod$dv_units
    # The clock time that model t = 0 corresponds to. Every acrophase printed
    # below is converted through it exactly once, by fck_acrophase_label().
    clock_o <- fck_clock_origin(mod)

    hdr <- function(x) cat("\n--- ", x, " ---\n", sep = "")

    cat("=== Harmonic Regression (Cosinor Analysis) Results ===\n\n")

    # ---- what is being modelled --------------------------------------------
    # AUDIT 1.7: the DV was never named. A report that does not say what Y is
    # cannot be checked for admissibility by anyone reading it.
    cat("Dependent variable: ", dv,
        if(!is.null(dvu)) paste0(" (", dvu, ")") else "", "\n", sep = "")
    if(is.finite(mod$dv_min) || is.finite(mod$dv_max)) {
      cat("  Admissible range: [",
          if(is.finite(mod$dv_min)) fmt2(mod$dv_min) else "-Inf", ", ",
          if(is.finite(mod$dv_max)) fmt2(mod$dv_max) else "Inf", "]\n", sep = "")
    } else {
      cat("  Admissible range: not specified (set it to have fitted values checked)\n")
    }
    cat("Period: ", fmtn(period, 0), " h\n", sep = "")
    cat("Number of harmonics: ", nh, "\n", sep = "")

    # ---- AUDIT 2.1: data source, stated as a choice with its consequences ---
    if(isTRUE(mod$using_smoothed)) {
      cat("Data: SMOOTHED (missing values interpolated by FDA)\n")
      cat("  ! Smoothing removes independent noise and induces residual autocorrelation.\n")
      cat("    R-squared is inflated, LOOCV is optimistic (a held-out point is partly\n")
      cat("    reconstructed from its neighbours), and the zero-amplitude F test is\n")
      cat("    anticonservative. Rhythm-significance rates below are UPPER BOUNDS.\n")
      cat("    Cosinor handles missing and unequally spaced data natively: re-run with\n")
      cat("    Data source = raw to see the difference.\n")
    } else {
      cat("Data: RAW (no smoothing applied)\n")
      if(!is.null(mod$subjects_with_nas) && mod$subjects_with_nas > 0) {
        cat("  ", mod$subjects_with_nas, " subjects have missing values; the cosinor uses",
            " the observations present.\n", sep = "")
      }
    }

    # ---- model specification ------------------------------------------------
    if(trend_type != "none") {
      trend_label <- switch(trend_type,
                            "linear" = "LINEAR (beta*t)",
                            "log" = "LOGARITHMIC (beta*log(t+1))",
                            "exp_sat" = "SATURATING EXPONENTIAL (A_sat*(1-e^(-t/tau)))",
                            "Unknown")
      cat("Homeostatic trend: ", trend_label, "\n", sep = "")
      cat("Two-process model: Process S (trend) separated from Process C (circadian)\n")
    } else {
      cat("Homeostatic trend: None (circadian only)\n")
    }

    # AUDIT 1.4.3: say which origin the coefficients are on.
    if(identical(mod$time_origin, "first_observation")) {
      cat("Time origin: FIRST OBSERVATION. t = 0 in the model is ",
          fck_clock_label(mod$origin_shift, mod$period), " on the clock.\n",
          "  Both the trend and the harmonics are anchored there, so the intercept is\n",
          "  the model value at the start of the recording and is directly\n",
          "  interpretable. Axes and hover text still read in CLOCK time: the plots add\n",
          "  the ", fmt1(mod$origin_shift), " h shift back before labelling.\n", sep = "")
    } else {
      cat("Time origin: MIDNIGHT (t = 0 at clock 00:00).\n")
      if(trend_type == "exp_sat") {
        cat("  ! The trend is anchored at the first observation (t - ",
            fmt2(mod$pop_mean_fit$t_offset %||% 0), ") while the harmonics are anchored at\n",
            "    midnight. The intercept is therefore the constant of a model with two\n",
            "    origins and is not the value at either. Set Time origin =\n",
            "    'first observation' to remove the ambiguity.\n", sep = "")
      }
    }

    # ---- AUDIT 2.3: convergence, reported honestly -------------------------
    fa <- mod$fit_audit
    hdr("Fit outcomes")
    if(!is.null(fa)) {
      cat(sprintf("Subjects attempted:         %d\n", fa$n_attempted))
      cat(sprintf("  Converged, interior:      %d\n", fa$n_converged))
      cat(sprintf("  Converged, on a bound:    %d  (%s)\n", fa$n_boundary,
                  if(isTRUE(fa$include_boundary)) "INCLUDED - see the bound table below"
                  else "excluded by your choice"))
      cat(sprintf("  Did not converge:         %d  (always excluded: no solution to average)\n",
                  fa$n_nonconverged))
      cat(sprintf("  Failed outright:          %d\n", fa$n_failed))
      cat(sprintf("Population summaries below use %d subject(s).\n", nrow(params)))
      if(isTRUE(fa$gate_relaxed))
        cat("  ! Too few fits pass the gate; ALL returned fits are included. Provisional.\n")
      if(fa$n_nonconverged + fa$n_boundary > 0)
        cat("  The original code counted all of these as successes, which is how an\n",
            "  R-squared range starting near 0.06 coexisted with '100% successfully fitted'.\n", sep = "")

      # ====================================================================
      # WHICH bounds, and who hit more than one
      #
      # A fit pinned to a constraint converged to the EDGE of the feasible
      # region, not to an interior optimum: the value is where the optimiser
      # was stopped, and its standard error is meaningless. That is worth
      # seeing per parameter, because a bound catching most of the sample is a
      # badly chosen bound rather than a sample full of odd subjects -- and a
      # subject pinned on TWO parameters at once is usually a ridge, where the
      # two trade off against each other along a flat direction of the
      # likelihood.
      # ====================================================================
      bs <- if(isTRUE(fa$include_boundary)) fa$bounds_kept else fa$bounds
      if(!is.null(bs) && bs$n_any > 0) {
        hdr("Parameter bounds hit")
        cat(sprintf("%d of %d summarised fit(s) sit on at least one bound (%s%%).\n",
                    bs$n_any, bs$n, fmt1(100 * bs$n_any / bs$n)))
        if(!is.null(bs$per_bound)) {
          cat("\nWhich bound:\n")
          cat(sprintf("  %-26s %8s %8s\n", "bound", "n", "% of n"))
          for(i in seq_len(nrow(bs$per_bound)))
            cat(sprintf("  %-26s %8d %7s%%\n", bs$per_bound$bound[i],
                        bs$per_bound$n[i], fmt1(bs$per_bound$pct[i])))
        }
        cat("\nHow many bounds per fit:\n")
        cat(sprintf("  %-12s %10s %8s\n", "bounds hit", "n fits", "% of n"))
        for(i in seq_len(nrow(bs$per_count)))
          cat(sprintf("  %-12d %10d %7s%%\n", bs$per_count$n_bounds[i],
                      bs$per_count$n_subjects[i], fmt1(bs$per_count$pct[i])))

        if(!is.null(bs$multi)) {
          cat(sprintf("\n%d fit(s) hit MORE THAN ONE bound. Two parameters pinned at once is\n",
                      nrow(bs$multi)))
          cat("usually a ridge: they trade off along a flat direction of the likelihood,\n")
          cat("so neither is separately identified and a group comparison of either is a\n")
          cat("comparison of where the optimiser stopped.\n\n")
          cat(sprintf("  %-16s %6s %10s   %s\n", "subject", "row", "n bounds", "bounds"))
          show_n <- min(nrow(bs$multi), 40)
          for(i in seq_len(show_n))
            cat(sprintf("  %-16s %6d %10d   %s\n", bs$multi$subject[i], bs$multi$row[i],
                        bs$multi$n_bounds[i], bs$multi$bounds[i]))
          if(nrow(bs$multi) > show_n)
            cat(sprintf("  ... and %d more (the full list is in the parameter CSV export).\n",
                        nrow(bs$multi) - show_n))
        }

        if(isTRUE(fa$include_boundary))
          cat("\nThese fits ARE included in everything below. Their values are real -- the\n",
              "optimiser did converge to them -- but they are the edge of the feasible\n",
              "region, so their SEs are meaningless and they pull the mean toward the\n",
              "bound. Untick 'Include fits that hit a parameter bound' to exclude them\n",
              "and see how much the summaries move.\n", sep = "")
        else
          cat("\nThese fits are EXCLUDED from everything below, which makes the sample\n",
              "smaller and possibly biased toward subjects the model happened to suit.\n", sep = "")
      }
    } else {
      cat(sprintf("Subjects summarised: %d\n", nrow(params)))
    }

    cat("\nNumber of time points: ", length(mod$time_vec), "\n", sep = "")
    if(length(mod$time_vec) <= 24) {
      if(isTRUE(mod$wrap_applied) && !is.null(mod$original_times)) {
        cat("Original times (clock):  ", paste(fmtn(mod$original_times, 1), collapse = ", "), "\n", sep = "")
        cat("Adjusted times (linear): ", paste(fmtn(mod$time_vec_clock %||% mod$time_vec, 1), collapse = ", "), "\n", sep = "")
        cat("(Times after midnight were adjusted for chronological order)\n")
      } else {
        cat("Time points used: ",
            paste(fmtn(mod$time_vec_clock %||% mod$time_vec, 2), collapse = ", "), "\n", sep = "")
      }
      if(isTRUE(mod$origin_shift > 0))
        cat("Model times (t = 0 at the first observation): ",
            paste(fmtn(mod$time_vec, 2), collapse = ", "), "\n",
            "  The times above are CLOCK times; these are what the coefficients were\n",
            "  fitted on. They differ by the ", fmt1(mod$origin_shift),
            " h origin shift, and every plot converts\n  back before labelling its axis.\n", sep = "")
    }
    diffs <- diff(mod$time_vec)
    if(length(unique(round(diffs, 2))) > 1) {
      cat("Spacing: UNEQUAL (", paste(unique(fmtn(diffs, 2)), collapse = ", "), ")\n", sep = "")
    } else {
      cat("Spacing: Equal (", fmt2(diffs[1]), ")\n", sep = "")
    }

    if(is.null(pop)) return(invisible(NULL))

    # ========================================================================
    # AUDIT 1.4: the central-value section, with the two quantities separated
    # ========================================================================
    hdr("Central value")
    # AUDIT: beta_0 is a COEFFICIENT, not a level. Calling it "the intercept at
    # t = 0" invited reading it as the value the response started at, which it
    # is not: at the first observation the harmonics are generally non-zero, and
    # with a saturating trend anchored there S(t_min) = 0 exactly, so the two
    # differ by the harmonic sum. Both are reported, under names that say what
    # they are.
    cat("Constant term (beta_0):                  ", fmt3(pop$intercept),
        "   [arithmetic mean of the fitted constants]\n", sep = "")
    .v0 <- fck_value_at(pop$mean_coefs, mod, min(mod$time_vec, na.rm = TRUE))
    if(is.finite(.v0))
      cat("Predicted value at the first observation (",
          fck_clock_label(fck_clock_origin(mod) + min(mod$time_vec, na.rm = TRUE),
                          period, show_day = FALSE), "): ",
          fmt2(.v0), if(!is.null(dvu)) paste0(" ", dvu) else "", "\n", sep = "")
    if(trend_type != "none" && is.finite(pop$rhythm_adjusted_mean)) {
      cat("MESOR (rhythm-adjusted mean over the observed window):  ",
          fmt3(pop$rhythm_adjusted_mean), "\n", sep = "")
      cat("  = time-average of beta_0 + S(t) over t in [",
          fmt2(pop$window[1]), ", ", fmt2(pop$window[2]), "] h, by integration.\n", sep = "")
      cat("  The intercept is NOT the MESOR when the model carries a trend: it is the\n")
      cat("  constant of the fit, and with a saturating trend it sits far below the\n")
      cat("  level the data actually occupy. The difference here is ",
          fmt2(pop$rhythm_adjusted_mean - pop$intercept), " units.\n", sep = "")
      if(is.finite(pop$harmonic_window_mean))
        cat("  (The harmonics contribute ", fmt3(pop$harmonic_window_mean),
            " to the window mean; that is exactly zero only over a whole number\n",
            "   of periods, and this window is not one.)\n", sep = "")
    } else {
      cat("MESOR (rhythm-adjusted mean): ", fmt3(pop$rhythm_adjusted_mean),
          "  [no trend, so it equals the intercept]\n", sep = "")
    }

    # ========================================================================
    # Vector-averaged rhythm parameters
    # ========================================================================
    hdr("Population rhythm parameters (VECTOR-averaged)")
    cat("Estimator: amplitude-weighted vector mean of (amplitude, acrophase) pairs.\n")
    cat("This is the correct population estimator for circular data and is unchanged.\n\n")
    # AUDIT: acrophases are printed in CLOCK time. They are estimated on the
    # model's axis, which under time_origin = "first_observation" is elapsed
    # hours from the first observation -- so 19.18 there is 03:11 on the clock,
    # and printing the elapsed number as a time of day was wrong by the origin
    # shift. fck_acrophase_label() is the only place that arithmetic happens.
    for(h in seq_len(nh)) {
      eff <- period / h
      cat(sprintf("  H%d: amplitude = %s, acrophase = %s  (%s deg)\n",
                  h, fmt3(pop$mean_amplitudes[h]),
                  fck_acrophase_label(hours = pop$mean_acrophases_time[h],
                                      period = period, harmonic = h,
                                      clock_origin = clock_o),
                  fmt2(phi_to_degrees(pop$mean_acrophases_rad[h]))))
      if(clock_o != 0)
        cat(sprintf("      (model-elapsed %s h + %s h origin)\n",
                    fmt2(pop$mean_acrophases_time[h]), fmt1(clock_o)))
      if(h > 1)
        cat(sprintf("      H%d repeats every %s h, so it has %d maxima per day; all are shown.\n",
                    h, fmtn(eff, 0), h))
    }

    # AUDIT: the H1 acrophase is where the FIRST HARMONIC peaks. The fitted
    # curve also carries the trend and the higher harmonics, so its maximum sits
    # elsewhere -- and it is the curve's maximum a reader sees on the plot.
    # Reporting both stops the acrophase looking wrong against its own figure.
    cpk <- fck_curve_peak_clock(pop$mean_coefs, mod)
    if(!is.null(cpk)) {
      cat(sprintf("\n  Maximum of the COMPLETE fitted curve: %s (value %s)\n",
                  fck_clock_label(cpk$peak_clock, period, show_day = FALSE),
                  fmt2(cpk$peak_value)))
      cat(sprintf("  Minimum of the complete fitted curve: %s (value %s)\n",
                  fck_clock_label(cpk$trough_clock, period, show_day = FALSE),
                  fmt2(cpk$trough_value)))
      cat("  NOT the H1 acrophase: the curve also contains ",
          if(trend_type != "none") "the homeostatic trend" else "no trend",
          if(nh > 1) " and the higher harmonics" else "", ".\n", sep = "")
      if(isTRUE(cpk$peak_at_edge))
        cat("  ! The maximum sits at the edge of the observed window: the curve was\n",
            "    still rising when the recording stopped, so this is a boundary value,\n",
            "    not a peak the data contain.\n", sep = "")
    }

    # ========================================================================
    # AUDIT 1.2: two resultants, unambiguously labelled
    # ========================================================================
    hdr("Circular concentration and the Rayleigh test")
    for(h in seq_len(nh)) {
      rr <- pop$resultants[[h]]
      if(is.null(rr)) next
      eff <- period / h
      cat(sprintf("H%d  n = %d\n", h, rr$n))
      cat(sprintf("  r-bar (UNWEIGHTED, for Rayleigh):        %s\n", fmt3(rr$r_unweighted)))
      cat(sprintf("  r-bar (AMPLITUDE-WEIGHTED, for the vector mean): %s\n", fmt3(rr$r_weighted)))
      cat(sprintf("  Rayleigh Z = %s, p = %s   [Z = n * r_unweighted^2]\n",
                  fmt1e(rr$rayleigh$Z), format.pval(rr$rayleigh$p, digits = 3, eps = 1e-16)))
      cat(sprintf("  Circular SD = %s h\n", fmt3(rr$circ_sd_hours)))
      if(h == 1)
        cat("  The Rayleigh test is defined on UNIT vectors (Mardia & Jupp 2000; Berens\n",
            "  2009). Using the amplitude-weighted resultant here would inflate Z; the\n",
            "  two are printed together so they can never be confused again.\n", sep = "")
    }

    # ========================================================================
    # AUDIT 1.1: ONE equation renderer, for the pooled fit and the groups alike
    # ========================================================================
    hdr("Model equation (symbolic)")
    sym <- "Y(t) = beta_0"
    if(trend_type == "linear")       sym <- paste0(sym, " + beta*t")
    else if(trend_type == "log")     sym <- paste0(sym, " + beta*log(t+1)")
    else if(trend_type == "exp_sat") sym <- paste0(sym, " + A_sat*(1 - e^(-t/tau))")
    for(h in seq_len(nh))
      sym <- paste0(sym, sprintf(" + A%d*cos(2pi*%d*t/%s - phi%d)", h, h, fmtn(period, 0), h))
    cat(sym, "\n")
    cat("  beta_0 = intercept (NOT the MESOR when a trend is present)\n")
    cat("  A_h, phi_h = amplitude and acrophase of harmonic h\n")
    if(trend_type == "exp_sat")
      cat("  A_sat = asymptote, tau = time constant (h)\n")

    hdr("Fitted model equation (pooled)")
    cat(fck_format_equation(pop$intercept, trend_type, pop$trend_coefs,
                            pop$mean_amplitudes, pop$mean_acrophases_rad,
                            period, pop$t_offset %||% 0), "\n")
    cat("  Rendered by the same function as the group equations below. The pooled\n")
    cat("  equation previously omitted the homeostatic term entirely.\n")

    # AUDIT 1.7: admissibility of the fitted curve
    if(is.finite(mod$dv_min) || is.finite(mod$dv_max)) {
      tt <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 400)
      yy <- predict_from_coefs(pop$mean_coefs, tt, period, nh, trend_type,
                               pop$t_offset %||% 0, 0)
      bc <- fck_check_bounds(yy, mod$dv_min, mod$dv_max)
      if(!bc$ok) {
        cat(sprintf("  ! The pooled fitted curve leaves the admissible range: min %s, max %s.\n",
                    fmt2(bc$min), fmt2(bc$max)))
        cat("    A fit that predicts impossible values is misspecified, not merely imprecise.\n")
      } else {
        cat(sprintf("  Fitted range over the window: [%s, %s] - within the admissible range.\n",
                    fmt2(bc$min), fmt2(bc$max)))
      }
    }

    # ========================================================================
    # AUDIT 1.6: arithmetic summaries, each labelled with its estimator
    # ========================================================================
    hdr("Individual parameters (ARITHMETIC means, +/- linear SD)")
    cat("Estimator: arithmetic. These are NOT the vector means above; for a circular\n")
    cat("quantity the two differ, and only the vector mean is a valid population value.\n\n")
    cat(sprintf("  Intercept:  %s (SD %s)   [arithmetic]\n",
                fmt3(mean(params$mesor, na.rm = TRUE)), fmt3(sd(params$mesor, na.rm = TRUE))))

    if(trend_type == "linear" && "trend_linear" %in% names(params)) {
      cat(sprintf("  Linear trend (beta): %s (SD %s) %s   [arithmetic]\n",
                  fmt4(mean(params$trend_linear, na.rm = TRUE)),
                  fmt4(sd(params$trend_linear, na.rm = TRUE)),
                  if(!is.null(dvu)) paste0(dvu, "/h") else "per hour"))
    } else if(trend_type == "log" && "trend_log" %in% names(params)) {
      cat(sprintf("  Log trend (beta): %s (SD %s)   [arithmetic]\n",
                  fmt4(mean(params$trend_log, na.rm = TRUE)),
                  fmt4(sd(params$trend_log, na.rm = TRUE))))
    } else if(trend_type == "exp_sat") {
      if("A_sat" %in% names(params)) {
        # AUDIT 1.7: the "(units)" placeholder is replaced by the real unit, or
        # by nothing at all when the user has not named one.
        cat(sprintf("  A_sat (asymptote):   %s (SD %s)%s   [arithmetic]\n",
                    fmt3(mean(params$A_sat, na.rm = TRUE)),
                    fmt3(sd(params$A_sat, na.rm = TRUE)),
                    if(!is.null(dvu)) paste0(" ", dvu) else ""))
      }
      if("tau" %in% names(params)) {
        tm <- mean(params$tau, na.rm = TRUE); ts <- sd(params$tau, na.rm = TRUE)
        cat(sprintf("  tau (time constant): %s (SD %s) h   [arithmetic]\n", fmt2(tm), fmt2(ts)))
        # AUDIT 2.2: the ridge, stated where the numbers are
        if(is.finite(tm) && is.finite(ts) && tm > 0 && ts / tm > 0.5)
          cat(sprintf("    ! SD/mean = %s. Over this window (1 - e^(-t/tau)) moves only from\n",
                      fmt2(ts / tm)),
              "      about 0.44 to 0.88 and is close to linear, so it is strongly collinear\n",
              "      with the intercept. A spread this wide is a likelihood ridge, not\n",
              "      population heterogeneity. See the conditioning section below.\n", sep = "")
      }
    }

    for(h in seq_len(nh)) {
      eff <- period / h
      ac <- params[[paste0("acrophase_time_", h)]]
      rr <- pop$resultants[[h]]
      cat(sprintf("  H%d amplitude: %s (SD %s)   [arithmetic]\n", h,
                  fmt3(mean(params[[paste0("amplitude_", h)]], na.rm = TRUE)),
                  fmt3(sd(params[[paste0("amplitude_", h)]], na.rm = TRUE))))
      # Circular SD is a DISPERSION and is invariant to the origin shift, so it
      # is printed as-is. Only the mean DIRECTION moves with the origin.
      cat(sprintf("  H%d acrophase: circular mean %s, circular SD %s h   [circular]\n",
                  h, fck_acrophase_label(phi_rad = rr$mean_dir_unweighted,
                                         period = period, harmonic = h,
                                         clock_origin = clock_o),
                  fmt2(rr$circ_sd_hours)))
      cat(sprintf("               (arithmetic mean %s, linear SD %s h - shown only to\n",
                  fck_acrophase_label(hours = mean(ac, na.rm = TRUE), period = period,
                                      harmonic = h, clock_origin = clock_o, all = FALSE),
                  fmt2(sd(ac, na.rm = TRUE))))
      cat("                make the difference visible; do not report these for a phase)\n")
      if(h > 1)
        cat(sprintf("               H%d has %d maxima per day, all shown above\n", h, h))
    }

    cat(sprintf("\n  R-squared: mean %s, range [%s, %s]\n",
                fmt3(mean(params$r_squared, na.rm = TRUE)),
                fmt3(min(params$r_squared, na.rm = TRUE)),
                fmt3(max(params$r_squared, na.rm = TRUE))))
    sig_rate <- 100 * sum(params$p_value < 0.05, na.rm = TRUE) / nrow(params)
    cat(sprintf("  Significant rhythms (p<0.05): %d / %d (%s%%)\n",
                sum(params$p_value < 0.05, na.rm = TRUE), nrow(params), fmt1(sig_rate)))
    cat("    The zero-amplitude test is now the harmonics GIVEN the trend (full vs\n")
    cat("    trend-only F test). It previously charged the rhythm with the whole\n")
    cat("    model's sum of squares, including everything Process S explained.\n")
    if(isTRUE(mod$using_smoothed))
      cat("    ! On smoothed data this rate is an UPPER BOUND (see the note at the top).\n")

    # ========================================================================
    # AUDIT 2.5: confidence intervals
    # ========================================================================
    if(!is.null(mod$boot_results)) {
      hdr(sprintf("Bootstrap CIs (B = %d, subject-level resampling)", mod$boot_results$B))
      cat(sprintf("  Intercept: [%s, %s]\n",
                  fmt3(mod$boot_results$mesor_ci[1]), fmt3(mod$boot_results$mesor_ci[2])))
      cat(sprintf("  Amplitude: [%s, %s]\n",
                  fmt3(mod$boot_results$amplitude_ci[1]), fmt3(mod$boot_results$amplitude_ci[2])))
      cat(sprintf("  Acrophase: [%s, %s]   (H1, clock time)\n",
                  fck_acrophase_label(hours = mod$boot_results$acrophase_ci[1],
                                      period = period, harmonic = 1,
                                      clock_origin = clock_o, all = FALSE),
                  fck_acrophase_label(hours = mod$boot_results$acrophase_ci[2],
                                      period = period, harmonic = 1,
                                      clock_origin = clock_o, all = FALSE)))
    }
    if(!is.null(mod$bingham_summary)) {
      hdr("Bingham joint confidence regions (Bingham et al. 1982)")
      cat("The elliptical joint region for the (amplitude, acrophase) pair, per subject,\n")
      cat("summarised across subjects. This is the standard cosinor reporting requirement\n")
      cat("and the app previously computed nothing of the kind for a nonlinear fit.\n\n")
      for(h in seq_len(nh)) {
        b <- mod$bingham_summary[[h]]
        if(is.null(b)) next
        cat(sprintf("  H%d: acrophase identified in %d / %d subjects (%s%%)\n",
                    h, b$n_identified, b$n, fmt1(100 * b$n_identified / b$n)))
        cat(sprintf("      median half-width: amplitude +/- %s, acrophase +/- %s h\n",
                    fmt3(b$median_amp_halfwidth), fmt2(b$median_acro_halfwidth_h)))
        cat("      Half-widths are durations and carry no origin; the acrophases they\n")
        cat("      bracket are the clock times reported above.\n")
        if(b$n_identified < b$n)
          cat("      where it is not identified the region covers the origin, i.e. that\n",
              "      subject has no resolvable phase at all. Quoting one would be worse\n",
              "      than quoting none.\n", sep = "")
      }
    }

    # ========================================================================
    # AUDIT 2.4: Delta-AICc, not mean AIC
    # ========================================================================
    if(!is.null(mod$model_selection)) {
      hdr("Model selection (Delta-AICc across a nested set)")
      ms <- mod$model_selection
      cat(sprintf("  %-18s %12s %10s %8s\n", "model", "AICc", "dAICc", "weight"))
      for(i in seq_len(nrow(ms)))
        cat(sprintf("  %-18s %12s %10s %8s\n", ms$model[i],
                    fmt2(ms$AICc[i]), fmt2(ms$dAICc[i]), fmt3(ms$weight[i])))
      cat("\n  AICc averaged per subject over the same subjects for every cell.\n")
      cat("  Absolute AIC/AICc/BIC means with SDs are not reported: with no competing\n")
      cat("  model they are constant offsets of one another (which is why all three\n")
      cat("  SDs printed identically), and with n/k = ", fmtn(length(mod$time_vec), 0),
          "/", fmtn(2 * nh + 1 +
                    switch(trend_type, "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0) + 1, 0),
          " a high R-squared is near-guaranteed.\n", sep = "")
    } else if("aicc" %in% names(params)) {
      hdr("Model selection")
      cat("  Enable 'Compare nested models' to get a Delta-AICc table with Akaike weights.\n")
      cat("  Absolute AIC/AICc/BIC are deliberately not summarised here: without a\n")
      cat("  competing model they carry no information, and AICc - AIC and BIC - AIC are\n")
      cat("  constant offsets (", fmt4(2 * 8 * 9 / (length(mod$time_vec) - 8 - 1)),
          " and ", fmt4(8 * (log(length(mod$time_vec)) - 2)),
          " for k = 8 here), which is why their SDs were identical.\n", sep = "")
      if("loocv_rmse" %in% names(params)) {
        is_cv <- mod$loocv_is_true_cv
        cat(sprintf("  %s: mean %s\n",
                    if(isTRUE(is_cv)) "LOOCV RMSE (genuine leave-one-out refits)"
                    else "In-sample residual RMSE (NOT cross-validated)",
                    fmt4(mean(params$loocv_rmse, na.rm = TRUE))))
        if(!isTRUE(is_cv))
          cat("    The nonlinear path used to report this number under the name 'LOOCV\n",
              "    RMSE' while computing it from the fitted values. It is training error.\n", sep = "")
      }
    }

    # ========================================================================
    # AUDIT 1.3: commonality analysis, and no dominance verdict
    # ========================================================================
    im <- pop$indiv_means
    if(!is.null(im) && !is.null(im$unique_S) && is.finite(im$unique_S)) {
      hdr("Variance decomposition (commonality analysis)")
      cat("Chevan & Sutherland (1991); Ray-Mukherjee et al. (2014).\n\n")
      cat(sprintf("  Unique to Process S (homeostatic): %s  (%s%% of total R-squared)\n",
                  fmt3(im$unique_S), fmt1(im$percent_S)))
      cat(sprintf("  Unique to Process C (circadian):   %s  (%s%%)\n",
                  fmt3(im$unique_C), fmt1(im$percent_C)))
      cat(sprintf("  Shared between S and C:            %s  (%s%%)\n",
                  fmt3(im$shared_SC), fmt1(im$percent_shared)))
      cat(sprintf("  ------------------------------------------------\n"))
      cat(sprintf("  Total R-squared:                   %s  (%s%%)\n",
                  fmt3(im$unique_S + im$unique_C + im$shared_SC),
                  fmt1(im$percent_S + im$percent_C + im$percent_shared)))
      cat("\n  These three sum to the total by construction. The previous output printed\n")
      cat("  two overlapping MARGINAL R-squareds (0.283 and 0.838 against a total of\n")
      cat("  0.892) and called their ratios 'proportions', which summed to 124.7%; the\n")
      cat("  shared component was invisible because it was being counted twice.\n")
      if(isTRUE(im$shared_SC < 0))
        cat("  ! The shared component is NEGATIVE: S and C are acting as mutual\n",
            "    suppressors here. That is a real finding about the design, not an error.\n", sep = "")
      cat("\n  The marginal R-squareds this is computed from, for audit:\n")
      cat(sprintf("    trend-only model:     %s\n", fmt3(im$r_squared_S)))
      cat(sprintf("    harmonics-only model: %s\n", fmt3(im$r_squared_C)))
      cat("\n  No dominance verdict is emitted. The old report concluded 'Circadian rhythm\n")
      cat("  (C) is the dominant component' from overlapping R-squareds; and any such\n")
      cat("  claim would have to account for C using ", 2 * nh,
          " parameters against S's ",
          switch(trend_type, "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0),
          ".\n", sep = "")
    }

    # ========================================================================
    # AUDIT 2.2: conditioning
    # ========================================================================
    if(!is.null(mod$conditioning)) {
      hdr("Parameter conditioning (AUDIT 2.2)")
      cn <- mod$conditioning
      if(!is.null(cn$mean_cor)) {
        cat("Mean within-subject parameter correlation matrix:\n")
        print(round(cn$mean_cor, 3))
        cat("\n")
      }
      if(!is.null(cn$tau_fixed_delta_aic))
        cat(sprintf("Free-tau vs tau fixed at %s h:  Delta-AIC = %s (%s)\n",
                    fmt1(cn$tau_fixed_value), fmt2(cn$tau_fixed_delta_aic),
                    if(cn$tau_fixed_delta_aic > 0)
                      "the free-tau fit is NOT better by AIC: tau is not identified by these data"
                    else "the free-tau fit is preferred"))
      if(!is.null(cn$kappa_before))
        cat(sprintf("Design-matrix condition number: %s (midnight origin) -> %s (first-observation origin)\n",
                    fmt1(cn$kappa_before), fmt1(cn$kappa_after)))
    }

    # ========================================================================
    # Groups
    # ========================================================================
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1) {
      ga <- attr(mod$group_fits, "audit")
      hdr(sprintf("Group-specific parameters%s",
                  if(!is.null(mod$group_var_name)) paste0(" (", mod$group_var_name, ")") else ""))

      # AUDIT 1.5: the accounting, before any group is printed
      n_in <- attr(mod$group_fits, "n_in_groups") %||% NA
      n_fit <- attr(mod$group_fits, "n_fitted") %||% nrow(params)
      cat(sprintf("Group sizes sum to %s; %s subjects were fitted.\n",
                  fmtn(n_in, 0), fmtn(n_fit, 0)))
      if(!is.null(ga) && ga$n_unassigned > 0)
        cat(sprintf("  %d subject(s) have no usable group label and appear as UNASSIGNED below.\n",
                    ga$n_unassigned),
            "  They were previously pooled into every population statistic and into no\n",
            "  group, which is why the group sizes did not add up to the total.\n", sep = "")
      if(!is.null(ga) && length(ga$dropped_small) > 0)
        cat(sprintf("  Group(s) not summarised (fewer than 3 fitted subjects): %s\n",
                    paste(ga$dropped_small, collapse = ", ")))
      if(is.finite(n_in) && is.finite(n_fit) && n_in != n_fit)
        cat("  ! THESE DO NOT RECONCILE. Some subjects are being dropped silently.\n")
      else
        cat("  These reconcile.\n")

      cat("\nEstimator key, per line: [arithmetic] = ordinary mean and SD.\n")
      cat("                         [vector]     = amplitude-weighted vector mean.\n")
      cat("                         [circular]   = circular mean / circular SD.\n")
      cat("A linear SD is never printed beside a vector- or circular-averaged value.\n")

      for(g_name in names(mod$group_fits)) {
        g <- mod$group_fits[[g_name]]
        cat(sprintf("\nGroup '%s' (n = %d):\n", g_name, g$n))

        # AUDIT: "Intercept (at t = 0)" invited reading beta_0 as a starting
        # level. It is a coefficient; the starting level is the fitted value at
        # the first observation, and the two differ by the harmonic sum there.
        cat(sprintf("  Constant term (beta_0):        %s (SD %s)   [arithmetic]\n",
                    fmt3(g$intercept), fmt3(g$sd_mesor)))
        .gv0 <- fck_value_at(g$mean_coefs, mod, min(mod$time_vec, na.rm = TRUE))
        if(is.finite(.gv0))
          cat(sprintf("  Predicted value at start (%s): %s%s\n",
                      fck_clock_label(fck_clock_origin(mod) +
                                        min(mod$time_vec, na.rm = TRUE), period,
                                      show_day = FALSE),
                      fmt2(.gv0), if(!is.null(dvu)) paste0(" ", dvu) else ""))
        if(trend_type != "none" && is.finite(g$rhythm_adjusted_mean))
          cat(sprintf("  MESOR (rhythm-adjusted mean):  %s   [integrated over the window]\n",
                      fmt3(g$rhythm_adjusted_mean)))

        if(!is.null(g$trend_params) && length(g$trend_params) > 0) {
          for(pn in names(g$trend_params)) {
            lbl <- switch(pn, "trend_linear" = "Linear trend", "trend_log" = "Log trend",
                          "A_sat" = "A_sat (asymptote)", "tau" = "tau (time constant, h)", pn)
            cat(sprintf("  %-30s %s (SD %s)   [arithmetic]\n", paste0(lbl, ":"),
                        fmt3(g$trend_params[[pn]]$mean), fmt3(g$trend_params[[pn]]$sd)))
          }
        }

        for(h in seq_len(nh)) {
          rr <- if(!is.null(g$resultants)) g$resultants[[h]] else NULL
          cat(sprintf("  H%d amplitude: %s   [vector]\n", h, fmt3(g$mean_amplitudes[h])))
          cat(sprintf("                %s (SD %s)   [arithmetic]\n",
                      fmt3(g$amp_arithmetic[h]), fmt3(g$sd_amplitudes[h])))
          cat(sprintf("  H%d acrophase: %s   [vector]\n", h,
                      fck_acrophase_label(hours = g$mean_acrophases_time[h],
                                          period = period, harmonic = h,
                                          clock_origin = clock_o)))
          if(!is.null(rr)) {
            cat(sprintf("                r-bar unweighted %s (Rayleigh Z = %s, p = %s)   [circular]\n",
                        fmt3(rr$r_unweighted), fmt1e(rr$rayleigh$Z),
                        format.pval(rr$rayleigh$p, digits = 3, eps = 1e-16)))
            cat(sprintf("                r-bar amplitude-weighted %s\n", fmt3(rr$r_weighted)))
          }
          if(h > 1)
            cat(sprintf("                H%d has %d maxima per day, all shown\n", h, h))
        }

        gpk <- fck_curve_peak_clock(g$mean_coefs, mod)
        if(!is.null(gpk))
          cat(sprintf("  Complete fitted curve peaks at %s (value %s), troughs at %s\n",
                      fck_clock_label(gpk$peak_clock, period, show_day = FALSE),
                      fmt2(gpk$peak_value),
                      fck_clock_label(gpk$trough_clock, period, show_day = FALSE)))

        cat("\n  Fitted equation:\n  ")
        cat(fck_format_equation(g$intercept, trend_type, g$trend_coefs,
                                g$mean_amplitudes, g$mean_acrophases_rad,
                                period, pop$t_offset %||% 0), "\n")

        if(!is.null(g$variance_decomp)) {
          vd <- g$variance_decomp
          cat("  Variance decomposition (commonality):\n")
          cat(sprintf("    unique S %s%%, unique C %s%%, shared %s%%\n",
                      fmt1(vd$percent_S), fmt1(vd$percent_C),
                      fmt1(vd$percent_shared %||% NA_real_)))
        }
      }
    }

    invisible(NULL)
  }

  output$harmonic_summary <- renderPrint({ .print_harmonic_summary() })

  # Parameters table
  output$harmonic_parameters_table <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Build summary table with all harmonics
    # AUDIT 1.4: this column is the fitted constant, not a rhythm-adjusted mean.
    param_names <- c("Intercept (b0, t=0)")
    mean_vals <- c(mean(params$mesor))
    sd_vals <- c(sd(params$mesor))
    min_vals <- c(min(params$mesor))
    max_vals <- c(max(params$mesor))
    
    for(h in 1:mod$n_harmonics) {
      amp_col <- paste0("amplitude_", h)
      acro_col <- paste0("acrophase_time_", h)
      
      # AUDIT: acrophase columns hold MODEL-elapsed hours; the table reports
      # clock time. The mean is converted; the SD is a dispersion, carries no
      # origin, and is left alone.
      .co <- fck_clock_origin(mod)
      .to_clock <- function(v) (v + .co) %% mod$period
      param_names <- c(param_names, paste0("Amplitude H", h),
                       paste0("Acrophase H", h, " (clock h)"))
      mean_vals <- c(mean_vals, mean(params[[amp_col]]), .to_clock(mean(params[[acro_col]])))
      sd_vals <- c(sd_vals, sd(params[[amp_col]]), sd(params[[acro_col]]))
      min_vals <- c(min_vals, min(params[[amp_col]]), .to_clock(min(params[[acro_col]])))
      max_vals <- c(max_vals, max(params[[amp_col]]), .to_clock(max(params[[acro_col]])))
    }
    
    param_names <- c(param_names, "R²", "% Rhythm")
    mean_vals <- c(mean_vals, mean(params$r_squared), mean(params$percent_rhythm))
    sd_vals <- c(sd_vals, sd(params$r_squared), sd(params$percent_rhythm))
    min_vals <- c(min_vals, min(params$r_squared), min(params$percent_rhythm))
    max_vals <- c(max_vals, max(params$r_squared), max(params$percent_rhythm))
    
    summary_df <- data.frame(
      Parameter = param_names,
      Mean = round(mean_vals, 3),
      SD = round(sd_vals, 3),
      Min = round(min_vals, 3),
      Max = round(max_vals, 3)
    )
    
    tagList(
      h4("Parameter Summary"),
      renderTable(summary_df, striped = TRUE, hover = TRUE, bordered = TRUE)
    )
  })
  
  # Fitted curves plot
  output$harmonic_fit_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    time_fine <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 200)
    subject_select <- input$harmonic_subject_select
    if(is.null(subject_select)) subject_select <- "mean"
    tp_response_type <- if(!is.null(mod$two_process_params$response_type)) mod$two_process_params$response_type else "gaussian"
    
    p <- plot_ly()
    
    if(subject_select == "all") {
      # Overlay all subjects - colored by group if available
      # AUDIT: one palette for the whole app, keyed by group NAME so a group
      # keeps its colour across every figure and is not repainted when another
      # group is filtered out. See server/02b_helpers_palette.R.
      .glv <- names(mod$group_fits)
      group_colors_hex <- fck_group_colors(.glv)
      if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
         !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        
        group_var <- values$covariates[[input$harmonic_group_var]]
        groups <- names(mod$group_fits)
        
        for(i in seq_along(mod$individual_fits)) {
          fit_i <- mod$individual_fits[[i]]
          if(!is.null(fit_i) && fit_i$success) {
            pred_i <- predict_cosinor(fit_i, time_fine)
            grp <- as.character(group_var[i])
            .c <- group_colors_hex[[grp]] %||% NA_character_
            line_color <- if(!is.na(.c)) fck_group_rgba(.c, 0.4) else 'rgba(100,100,100,0.3)'
            
            p <- p %>% add_lines(x = time_fine, y = pred_i, 
                                 line = list(color = line_color, width = 1),
                                 showlegend = FALSE, hoverinfo = "skip")
          }
        }
        
        # Add group mean curves using all harmonics
        for(g_idx in seq_along(groups)) {
          g_name <- groups[g_idx]
          g_fit <- mod$group_fits[[g_name]]
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            params$group <- group_var[params$subject]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            group_idx <- which(group_var == g_name)
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
            g_pred <- predict_two_process_mean_curve(grp_params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            # Use mean_coefs for full multi-harmonic prediction
            g_pred <- predict_from_coefs(g_fit$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          
          p <- p %>% add_lines(x = time_fine, y = g_pred,
                               line = list(color = unname(group_colors_hex[[g_name]]), width = 3),
                               name = paste("Mean:", g_name))
        }
        
        # Show harmonic components for grouped data (same as "mean" view)
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          params$group <- group_var[params$subject]
          
          # Lighter versions of group colors for components
          group_comp_colors <- fck_group_colors(names(mod$group_fits))
          
          # Show group-specific harmonic components
          for(g_idx in seq_along(groups)) {
            g_name <- groups[g_idx]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            
            if(nrow(grp_params) > 0) {
              grp_mesor <- mean(grp_params$mesor, na.rm = TRUE)
              
              grp_coefs <- c(grp_mesor)
              trend_coefs <- get_mean_trend_coefs(grp_params, mod$trend_type)
              if(length(trend_coefs) > 0) {
                grp_coefs <- c(grp_coefs, trend_coefs)
              }
              for(h in 1:mod$n_harmonics) {
                grp_coefs <- c(grp_coefs, 
                               mean(grp_params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(grp_params[[paste0("beta_sin_", h)]], na.rm = TRUE))
              }
              
              n_trend_coefs <- length(trend_coefs)
              coef_offset <- 1 + n_trend_coefs
              
              for(h in 1:mod$n_harmonics) {
                omega <- 2 * pi * h / mod$period
                beta_cos <- grp_coefs[coef_offset + 2 * h - 1]
                beta_sin <- grp_coefs[coef_offset + 2 * h]
                comp_vals <- grp_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
                
                p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                     line = list(color = unname(group_comp_colors[[g_name]]), width = 1.5, dash = 'dot'),
                                     name = paste0("H", h, " (", g_name, ")"))
              }
              
              # Plot trend line using helper function
              if(mod$trend_type == "two_process") {
                group_idx <- which(group_var == g_name)
                S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
                trend_line <- compute_two_process_trend_line(grp_params, time_fine, S_t, tp_response_type)
              } else {
                trend_line <- compute_trend_line(grp_params, mod$trend_type, time_fine, mod$t_offset)
              }
              if(!is.null(trend_line)) {
                p <- p %>% add_lines(x = time_fine, y = trend_line,
                                     line = list(color = unname(group_colors_hex[[g_name]]), width = 1.5, dash = 'dash'),
                                     name = paste0(get_trend_label(mod$trend_type), " (", g_name, ")"))
              }
            }
          }
          
          # Overall harmonics in gray
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = 'gray40', width = 1, dash = 'dot'),
                                 name = paste0("H", h, " (overall)"))
          }
          
          # Plot overall trend line using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            overall_trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            overall_trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(overall_trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = overall_trend_line,
                                 line = list(color = 'black', width = 1, dash = 'dash'),
                                 name = paste0(get_trend_label(mod$trend_type), " (overall)"))
          }
        }
        
      } else {
        # No groups - all same color
        for(i in seq_along(mod$individual_fits)) {
          fit_i <- mod$individual_fits[[i]]
          if(!is.null(fit_i) && fit_i$success) {
            pred_i <- predict_cosinor(fit_i, time_fine)
            p <- p %>% add_lines(x = time_fine, y = pred_i, 
                                 line = list(color = 'rgba(100, 100, 100, 0.3)', width = 1),
                                 showlegend = FALSE, hoverinfo = "skip")
          }
        }
        
        # Add mean curve - either from pop_mean_fit or computed from individual params
        if(!is.null(mod$pop_mean_fit)) {
          pop <- mod$pop_mean_fit
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            mean_pred <- predict_from_coefs(pop$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          p <- p %>% add_lines(x = time_fine, y = mean_pred,
                               line = list(color = 'red', width = 3), name = "Population Mean")
        } else {
          # Compute mean from individual parameters (for individual analysis type)
          params <- mod$individual_params
          mean_mesor <- mean(params$mesor, na.rm = TRUE)
          mean_coefs <- c(mean_mesor)
          
          mean_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(mean_trend_coefs) > 0) {
            mean_coefs <- c(mean_coefs, mean_trend_coefs)
          }
          
          for(h in 1:mod$n_harmonics) {
            mean_coefs <- c(mean_coefs, 
                            mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                            mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            mean_pred <- predict_from_coefs(mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          p <- p %>% add_lines(x = time_fine, y = mean_pred,
                               line = list(color = 'red', width = 3), name = "Population Mean")
        }
        
        # Show harmonic components for no-groups case
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste0("H", h))
          }
          
          # Plot trend line using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
      }
      
    } else if(subject_select == "mean") {
      # Show mean curve(s) - either population or by group
      
      # Check if we have group fits
      if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1) {
        # Show each group's mean curve - use hex colors for proper transparency
        # AUDIT: one palette for the whole app, keyed by group NAME so a group
      # keeps its colour across every figure and is not repainted when another
      # group is filtered out. See server/02b_helpers_palette.R.
      .glv <- names(mod$group_fits)
      group_colors_hex <- fck_group_colors(.glv)
        color_idx <- 1
        
        for(g_name in names(mod$group_fits)) {
          g_fit <- mod$group_fits[[g_name]]
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            group_var <- values$covariates[[input$harmonic_group_var]]
            params$group <- group_var[params$subject]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            group_idx <- which(group_var == g_name)
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
            g_pred <- predict_two_process_mean_curve(grp_params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            # Use mean_coefs for full multi-harmonic prediction
            g_pred <- predict_from_coefs(g_fit$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }

          # Only add line if predictions are valid (not all NA/NaN/Inf)
          if(any(is.finite(g_pred))) {
            p <- p %>% add_lines(x = time_fine, y = g_pred,
                                 line = list(color = unname(group_colors_hex[[g_name]]), width = 3),
                                 name = paste("Group:", g_name))
          }
          
          # Add confidence band if requested (approximate using first harmonic SD)
          if(isTRUE(input$harmonic_show_ci) && !is.null(g_fit$sd_amplitudes)) {
            # Simple approximation: scale curve by amplitude uncertainty
            amp_se <- g_fit$sd_amplitudes[1] / sqrt(g_fit$n)
            scale_upper <- 1 + 1.96 * amp_se / g_fit$mean_amplitudes[1]
            scale_lower <- max(0, 1 - 1.96 * amp_se / g_fit$mean_amplitudes[1])
            
            g_upper <- g_fit$mean_mesor + (g_pred - g_fit$mean_mesor) * scale_upper
            g_lower <- g_fit$mean_mesor + (g_pred - g_fit$mean_mesor) * scale_lower
            
            # Use rgba with 0.25 transparency for confidence band
            ci_color <- fck_group_rgba(group_colors_hex[[g_name]], 0.25)
            
            # Add ribbon with legend entry so it can be toggled
            p <- p %>% add_ribbons(x = time_fine, ymin = g_lower, ymax = g_upper,
                                   line = list(color = 'transparent'),
                                   fillcolor = ci_color,
                                   name = paste("95% CI:", g_name),
                                   legendgroup = paste0("group_", g_name),
                                   showlegend = TRUE, hoverinfo = "skip")
          }
          
          color_idx <- color_idx + 1
        }
        
        # Show harmonic components for grouped data
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          group_var <- values$covariates[[input$harmonic_group_var]]
          params$group <- group_var[params$subject]
          
          # Lighter versions of group colors for components
          group_comp_colors <- fck_group_colors(names(mod$group_fits))
          
          # Show group-specific harmonic components
          g_idx <- 1
          for(g_name in names(mod$group_fits)) {
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            
            if(nrow(grp_params) > 0) {
              grp_mesor <- mean(grp_params$mesor, na.rm = TRUE)
              
              # Build group-specific coefficients using helper
              grp_coefs <- c(grp_mesor)
              grp_trend_coefs <- get_mean_trend_coefs(grp_params, mod$trend_type)
              if(length(grp_trend_coefs) > 0) {
                grp_coefs <- c(grp_coefs, grp_trend_coefs)
              }
              for(h in 1:mod$n_harmonics) {
                grp_coefs <- c(grp_coefs, 
                               mean(grp_params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(grp_params[[paste0("beta_sin_", h)]], na.rm = TRUE))
              }
              
              n_trend_coefs <- length(grp_trend_coefs)
              coef_offset <- 1 + n_trend_coefs
              
              # Plot each harmonic for this group
              for(h in 1:mod$n_harmonics) {
                omega <- 2 * pi * h / mod$period
                beta_cos <- grp_coefs[coef_offset + 2 * h - 1]
                beta_sin <- grp_coefs[coef_offset + 2 * h]
                comp_vals <- grp_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
                
                p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                     line = list(color = unname(group_comp_colors[[g_name]]), width = 1.5, dash = 'dot'),
                                     name = paste0("H", h, " (", g_name, ")"))
              }
              
              # Plot trend line using helper function
              if(mod$trend_type == "two_process") {
                group_idx <- which(group_var == g_name)
                S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
                grp_trend_line <- compute_two_process_trend_line(grp_params, time_fine, S_t, tp_response_type)
              } else {
                grp_trend_line <- compute_trend_line(grp_params, mod$trend_type, time_fine, mod$t_offset)
              }
              if(!is.null(grp_trend_line)) {
                p <- p %>% add_lines(x = time_fine, y = grp_trend_line,
                                     line = list(color = unname(group_colors_hex[[g_name]]), width = 1.5, dash = 'dash'),
                                     name = paste0(get_trend_label(mod$trend_type), " (", g_name, ")"))
              }
            }
            g_idx <- g_idx + 1
          }
          
          # Also show overall (pooled) harmonics in gray for reference
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = 'gray40', width = 1, dash = 'dot'),
                                 name = paste0("H", h, " (overall)"))
          }
          
          # Plot overall trend using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            overall_trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            overall_trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(overall_trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = overall_trend_line,
                                 line = list(color = 'black', width = 1, dash = 'dash'),
                                 name = paste0(get_trend_label(mod$trend_type), " (overall)"))
          }
        }
        
        # Also add individual data points colored by group if requested
        if(isTRUE(input$harmonic_show_data)) {
          group_var <- values$covariates[[input$harmonic_group_var]]
          for(i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if(!is.null(fit_i) && fit_i$success) {
              grp <- as.character(group_var[i])
              grp_idx <- which(names(mod$group_fits) == grp)
              .c <- group_colors_hex[[grp]] %||% NA_character_
              pt_color <- if(!is.na(.c)) fck_group_rgba(.c, 0.3) else 'rgba(100,100,100,0.2)'
              
              p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                     marker = list(color = pt_color, size = 3),
                                     showlegend = FALSE, hoverinfo = "skip")
            }
          }
        }
        
      } else {
        # No groups - compute mean from individual parameters
        params <- mod$individual_params
        mean_mesor <- mean(params$mesor, na.rm = TRUE)
        
        # Build mean_coefs from individual parameters using helper
        mean_coefs <- c(mean_mesor)
        
        # Add mean trend coefficients based on type
        mean_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
        if(length(mean_trend_coefs) > 0) {
          mean_coefs <- c(mean_coefs, mean_trend_coefs)
        }
        
        # Add mean harmonic coefficients
        for(h in 1:mod$n_harmonics) {
          mean_coefs <- c(mean_coefs, 
                          mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                          mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
        }
        
        if(mod$trend_type == "two_process") {
          S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
          mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
        } else {
          mean_pred <- predict_from_coefs(mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
        }
        
        p <- p %>% add_lines(x = time_fine, y = mean_pred,
                             line = list(color = 'firebrick', width = 3), name = "Population Mean")
        
        # Add confidence band if requested (use SD of individual amplitudes)
        if(isTRUE(input$harmonic_show_ci)) {
          amp_sd <- sd(params$amplitude_1, na.rm = TRUE)
          amp_mean <- mean(params$amplitude_1, na.rm = TRUE)
          n_valid <- sum(!is.na(params$amplitude_1))
          amp_se <- amp_sd / sqrt(n_valid)
          
          scale_upper <- 1 + 1.96 * amp_se / amp_mean
          scale_lower <- max(0, 1 - 1.96 * amp_se / amp_mean)
          
          upper_pred <- mean_mesor + (mean_pred - mean_mesor) * scale_upper
          lower_pred <- mean_mesor + (mean_pred - mean_mesor) * scale_lower
          
          p <- p %>% add_ribbons(x = time_fine, ymin = lower_pred, ymax = upper_pred,
                                 line = list(color = 'transparent'),
                                 fillcolor = 'rgba(178, 34, 34, 0.2)',
                                 name = "95% CI")
        }
        
        # Show harmonic components if requested
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          n_trend_coefs <- length(mean_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- mean_coefs[coef_offset + 2 * h - 1]
            beta_sin <- mean_coefs[coef_offset + 2 * h]
            comp_vals <- mean_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste("H", h, "(τ/", h, ")", sep=""))
          }
          
          # Show trend component if present using helper
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
        
        # Add individual data as faint points if requested
        if(isTRUE(input$harmonic_show_data)) {
          for(i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if(!is.null(fit_i) && fit_i$success) {
              p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                     marker = list(color = 'rgba(100, 100, 100, 0.2)', size = 3),
                                     showlegend = FALSE, hoverinfo = "skip")
            }
          }
        }
      }
      
    } else {
      # Single subject
      i <- as.integer(subject_select)
      fit_i <- mod$individual_fits[[i]]
      
      if(!is.null(fit_i) && fit_i$success) {
        pred_i <- predict_cosinor(fit_i, time_fine)
        
        p <- p %>% add_lines(x = time_fine, y = pred_i,
                             line = list(color = 'steelblue', width = 2), 
                             name = paste("Subject", i))
        
        if(isTRUE(input$harmonic_show_data)) {
          p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                 marker = list(color = 'steelblue', size = 6),
                                 name = "Observed Data")
        }
        
        # Show harmonic components if requested
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          components <- get_harmonic_components(fit_i, time_fine)
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          for(h in 1:mod$n_harmonics) {
            comp_name <- paste0("harmonic_", h)
            comp_vals <- components$mesor[1] + components[[comp_name]]
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste("H", h, "(τ/", h, ")", sep=""))
          }
          
          # Show linear trend component if present
          if(!is.null(components$trend) && mod$trend_type != "none") {
            trend_line <- components$mesor[1] + components$trend
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
      } else {
        # Fit failed - show detailed message
        fail_msg <- paste("Subject", i, "fit failed")
        if(!is.null(fit_i) && !is.null(fit_i$n_valid)) {
          fail_msg <- paste0(fail_msg, "\n(", fit_i$n_valid, " valid points, need ", fit_i$n_required, "+)")
        }
        if(!is.null(fit_i) && !is.null(fit_i$message)) {
          fail_msg <- paste0(fail_msg, "\nReason: ", fit_i$message)
        }
        
        p <- p %>% add_annotations(
          x = 0.5, y = 0.5,
          text = fail_msg,
          showarrow = FALSE,
          xref = "paper", yref = "paper",
          font = list(size = 14, color = "red")
        )
      }
    }
    
    # ========================================================================
    # X AXIS: linear time underneath, CLOCK TIME on the labels and the hover
    #
    # The model is fitted on unwrapped linear time (8, 9, ... 30) and has to be:
    # the harmonics would not care -- cos(2*pi*h*t/T) gives the same value at
    # t = 3 and t = 27 -- but the TREND is not periodic, so 08:00 on day one and
    # 08:00 on day two must be different values of t or the two days collapse
    # onto one point and the homeostatic rise cannot be estimated at all.
    #
    # So linear time stays the computational axis and clock time is purely a
    # display transform, applied here and in the hover text. Nothing in this
    # block feeds a fit.
    #
    # What was wrong before: labels appeared only when the recording happened to
    # wrap past the period, so a within-day recording got bare numbers; the
    # label was built as paste0(t %% period, ":00"), which renders a half-past
    # tick as "8.5:00"; there was no zero padding and no day marker, so 08:00 on
    # the first day and on the second were indistinguishable; and the HOVER was
    # never converted at all, which is why pointing at the curve reported
    # 27.01508 instead of 03:01.
    # ========================================================================
    # The plot's x coordinates are MODEL time. Under
    # time_origin = "first_observation" that starts at 0 while the recording
    # started at 08:00, so the labels must add the shift back or the axis reads
    # eight hours early -- which is exactly what it did: the fit appeared to
    # start at midnight, and flipping the toggle looked like it inverted the
    # setting. The shift is applied to the TICK TEXT and the hover only; the
    # tick POSITIONS stay on the model's own axis, because that is where the
    # data are drawn.
    shift <- mod$origin_shift %||% 0
    time_range <- range(mod$time_vec, na.rm = TRUE)
    ticks <- fck_clock_ticks(time_range + shift, mod$period)

    x_axis <- list(
      title = if (mod$period == 24) "Clock time" else sprintf("Time (period = %g)", mod$period)
    )
    if (!is.null(ticks)) {
      x_axis$tickmode <- "array"
      x_axis$tickvals <- ticks$vals - shift      # back onto the model axis
      x_axis$ticktext <- ticks$text
    }
    # A recording that crosses midnight has two 08:00s on the same axis. Mark
    # each period boundary so the reader can see which day a point belongs to.
    day_lines <- list()
    if (diff(time_range) > mod$period * 0.5) {
      cr <- time_range + shift                   # midnight is a CLOCK event
      bnds <- seq(ceiling(cr[1] / mod$period) * mod$period, cr[2], by = mod$period)
      bnds <- bnds[bnds > cr[1] & bnds < cr[2]] - shift
      day_lines <- lapply(bnds, function(b) list(
        type = "line", x0 = b, x1 = b, yref = "paper", y0 = 0, y1 = 1,
        line = list(color = "rgba(11,11,11,0.25)", width = 1, dash = "dot")))
    }

    # Convert the hover on EVERY trace that carries an x in model time. Traces
    # already carrying their own text (the CI band, annotations) are left alone.
    for (i in seq_along(p$x$attrs)) {
      a <- p$x$attrs[[i]]
      if (is.null(a$x)) next
      xv <- tryCatch(as.numeric(a$x), warning = function(w) NULL, error = function(e) NULL)
      if (is.null(xv) || !any(is.finite(xv))) next
      # The whole hover string goes in `text` with hoverinfo = "text", rather
      # than a hovertemplate: plotly recycles a scalar template to one copy per
      # point (I() does not prevent it), and `text` is an array we need anyway.
      yv <- tryCatch(as.numeric(a$y), warning = function(w) NULL, error = function(e) NULL)
      nm <- if (!is.null(a$name)) as.character(a$name) else "fit"
      dvl <- if (!is.null(mod$dv_name)) mod$dv_name else "Response"
      xc <- xv + shift                           # model time -> clock time
      p$x$attrs[[i]]$text <- if (is.null(yv) || length(yv) != length(xv))
        sprintf("%s<br>%s", fck_clock_label(xc, mod$period), nm)
      else
        sprintf("%s<br>%s: %s<br>%s", fck_clock_label(xc, mod$period), dvl, fmt2(yv), nm)
      p$x$attrs[[i]]$hoverinfo <- "text"
    }

    p %>% layout(
      title = "Harmonic Regression Fit",
      xaxis = x_axis,
      yaxis = list(title = if (!is.null(mod$dv_name))
        paste0(mod$dv_name, if (!is.null(mod$dv_units)) paste0(" (", mod$dv_units, ")") else "")
        else "Response"),
      shapes = day_lines,
      hovermode = "closest",
      showlegend = TRUE
    )
  })
  
  # Polar plot for acrophase
  output$harmonic_polar_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_polar)) as.integer(input$selected_harmonic_polar) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    # Get amplitude and acrophase for selected harmonic
    amp_col <- paste0("amplitude_", h)
    acro_rad_col <- paste0("acrophase_rad_", h)
    
    # Convert to degrees for polar plot
    theta_deg <- phi_to_degrees(params[[acro_rad_col]])
    r <- params[[amp_col]]
    
    p <- plot_ly(type = 'scatterpolar', mode = 'markers')
    
    # Check if we have group fits - color points by group
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {

      group_var <- values$covariates[[input$harmonic_group_var]]
      group_colors <- fck_group_colors(names(mod$group_fits))

      groups <- names(mod$group_fits)
      for(g_idx in seq_along(groups)) {
        g_name <- groups[g_idx]
        g_mask <- !is.na(group_var[params$subject]) & group_var[params$subject] == g_name

        if(sum(g_mask) > 0) {
          p <- p %>% add_trace(
            r = r[g_mask], theta = theta_deg[g_mask],
            type = 'scatterpolar', mode = 'markers',
            marker = list(size = 8, color = unname(group_colors[[g_name]]), opacity = 0.7),
            name = paste("Group:", g_name)
          )
        }

        # Add group mean vector for selected harmonic
        g_fit <- mod$group_fits[[g_name]]
        acro_deg <- phi_to_degrees(g_fit$mean_acrophases_rad[h])
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, g_fit$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          type = 'scatterpolar', mode = 'lines+markers',
          line = list(color = unname(group_colors[[g_name]]), width = 3),
          marker = list(size = 12, color = unname(group_colors[[g_name]]), symbol = 'diamond'),
          name = paste("Mean:", g_name)
        )
      }

      # Add population mean vector if requested (even when groups present)
      if(isTRUE(input$polar_show_mean) && !is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        acro_deg <- phi_to_degrees(pop$mean_acrophases_rad[h])
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, pop$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          mode = 'lines+markers',
          line = list(color = 'black', width = 4, dash = 'dash'),
          marker = list(size = 14, color = 'black', symbol = 'star'),
          name = "Overall Population Mean"
        )
      }

    } else {
      # No groups - show all points same color
      p <- p %>% add_trace(
        r = r, theta = theta_deg,
        marker = list(size = 8, color = 'steelblue', opacity = 0.7),
        name = "Individual"
      )

      # Add mean vector if requested
      if(isTRUE(input$polar_show_mean) && !is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        acro_deg <- phi_to_degrees(pop$mean_acrophases_rad[h])
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, pop$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          mode = 'lines+markers',
          line = list(color = 'red', width = 3),
          marker = list(size = 12, color = 'red', symbol = 'diamond'),
          name = "Population Mean"
        )
      }
    }

    # Add confidence ellipse if requested (for all data, regardless of groups)
    if(isTRUE(input$polar_show_ellipse) && length(r) >= 3) {
      # Convert polar to Cartesian for ellipse calculation
      theta_rad <- theta_deg * pi / 180
      x <- r * cos(theta_rad)
      y <- r * sin(theta_rad)

      # Remove NAs
      valid <- !is.na(x) & !is.na(y)
      x <- x[valid]
      y <- y[valid]

      if(length(x) >= 3) {
        # Calculate 95% confidence ellipse
        mx <- mean(x)
        my <- mean(y)

        # Covariance matrix
        cov_mat <- cov(cbind(x, y))

        # Eigenvalues and eigenvectors
        eig <- eigen(cov_mat)

        # Chi-square value for 95% confidence (2 degrees of freedom)
        chi_sq <- qchisq(0.95, df = 2)

        # Ellipse parameters
        a <- sqrt(chi_sq * eig$values[1])  # Semi-major axis
        b <- sqrt(chi_sq * eig$values[2])  # Semi-minor axis
        angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])  # Rotation angle

        # Generate ellipse points
        t <- seq(0, 2*pi, length.out = 100)
        ellipse_x <- mx + a * cos(t) * cos(angle) - b * sin(t) * sin(angle)
        ellipse_y <- my + a * cos(t) * sin(angle) + b * sin(t) * cos(angle)

        # Convert back to polar coordinates
        ellipse_r <- sqrt(ellipse_x^2 + ellipse_y^2)
        ellipse_theta_rad <- atan2(ellipse_y, ellipse_x)
        ellipse_theta_deg <- phi_to_degrees(ellipse_theta_rad)
        ellipse_theta_deg[ellipse_theta_deg < 0] <- ellipse_theta_deg[ellipse_theta_deg < 0] + 360

        # Add ellipse as a trace
        p <- p %>% add_trace(
          r = ellipse_r,
          theta = ellipse_theta_deg,
          type = 'scatterpolar',
          mode = 'lines',
          line = list(color = 'rgba(255, 0, 0, 0.5)', width = 2, dash = 'dot'),
          name = "95% Confidence Ellipse",
          showlegend = TRUE
        )
      }
    }
    
    # Adjust angular axis based on harmonic
    # For H2 (12h period), show 12 hours; for H3, show 8 hours, etc.
    # AUDIT: the angular positions are MODEL phase, and the point cloud, the
    # group vectors and the ellipse all share that frame -- so the geometry is
    # right and nothing here is moved. Only the axis LABELS were wrong: they
    # read elapsed hours as if they were clock times. Relabelled through the
    # same origin the rest of the app uses.
    #
    # A dial for harmonic h spans period/h, so every angle on it corresponds to
    # h different clock times. The labels show the first and the subtitle says
    # how often it recurs, rather than silently picking one.
    clock_o <- fck_clock_origin(mod)
    effective_period <- mod$period / h
    n_ticks <- min(12, effective_period)
    tick_step <- effective_period / n_ticks
    tick_vals <- seq(0, 360 - 360/n_ticks, by = 360/n_ticks)
    tick_elapsed <- seq(0, effective_period - tick_step, by = tick_step)
    tick_labels <- vapply(tick_elapsed, function(e)
      fck_clock_label((e + clock_o) %% mod$period, mod$period, show_day = FALSE),
      character(1))

    # The title sat on top of the dial: a polar trace fills its plotting area
    # edge to edge, so a centred title with no reserved space lands on the 11-13
    # o'clock labels. The subtitle is moved out to a paper-anchored annotation
    # under the title and the polar domain is pulled down to leave room, rather
    # than shrinking the font until it stops colliding.
    p %>% layout(
      title = list(
        text = sprintf("Acrophase polar plot - H%d (effective period %s h)", h,
                       fmt1(effective_period)),
        x = 0.5, xanchor = "center", y = 0.98, yanchor = "top",
        font = list(size = 15)),
      annotations = list(list(
        text = if(h > 1)
                 sprintf("clock times; each angle recurs every %s h", fmt1(effective_period))
               else if(clock_o != 0)
                 sprintf("clock times; the model origin is %s",
                         fck_clock_label(clock_o, mod$period, show_day = FALSE))
               else "clock times",
        x = 0.5, y = 0.925, xref = "paper", yref = "paper",
        xanchor = "center", yanchor = "top", showarrow = FALSE,
        font = list(size = 11, color = "#52514e"))),
      polar = list(
        domain = list(y = c(0, 0.88)),   # the room the title and subtitle need
        radialaxis = list(title = "Amplitude", tickangle = 0, angle = 90),
        angularaxis = list(
          direction = "clockwise",
          rotation = 90,
          tickmode = "array",
          tickvals = tick_vals,
          ticktext = tick_labels
        )
      ),
      margin = list(t = 70, b = 40),
      showlegend = TRUE
    )
  })
  
  # Amplitude histogram
  output$harmonic_amplitude_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_dist)) as.integer(input$selected_harmonic_dist) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    effective_period <- mod$period / h
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = .data[[amp_col]], fill = as.factor(group))) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 15) +
        scale_fill_brewer(palette = "Set1", name = "Group") +
        theme_minimal() +
        labs(title = paste0("Distribution of Amplitudes (H", h, ", period=", round(effective_period, 1), "h) by Group"), 
             x = "Amplitude", y = "Count")
      
      ggplotly(g)
    } else {
      plot_ly(x = params[[amp_col]], type = "histogram", 
              marker = list(color = 'steelblue', line = list(color = 'white', width = 1))) %>%
        layout(title = paste0("Distribution of Amplitudes (H", h, ", period=", round(effective_period, 1), "h)"),
               xaxis = list(title = "Amplitude"),
               yaxis = list(title = "Count"))
    }
  })
  
  # Acrophase histogram (circular)
  output$harmonic_acrophase_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    period <- mod$period
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_dist)) as.integer(input$selected_harmonic_dist) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    acro_col <- paste0("acrophase_time_", h)
    effective_period <- period / h

    # AUDIT: acrophase_time_h holds MODEL-elapsed hours. Plotted raw it reads as
    # a clock, which is wrong by the origin shift. Converted here through the
    # same origin as every other acrophase in the app.
    clock_o <- fck_clock_origin(mod)
    params$acro_clock <- (params[[acro_col]] + clock_o) %% effective_period
    x_lab <- if(h > 1)
      sprintf("Acrophase (clock h, modulo %s h - H%d has %d maxima per day)",
              fmt1(effective_period), h, h)
    else "Acrophase (clock h)"

    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {

      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove NAs

      g <- ggplot(params, aes(x = .data[["acro_clock"]], fill = as.factor(group))) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 12) +
        scale_fill_brewer(palette = "Set1", name = "Group") +
        scale_x_continuous(limits = c(0, effective_period)) +
        theme_minimal() +
        labs(title = paste0("Distribution of acrophases (H", h, ") by group"),
             x = x_lab, y = "Count")

      ggplotly(g)
    } else {
      plot_ly(x = params$acro_clock, type = "histogram",
              marker = list(color = 'firebrick', line = list(color = 'white', width = 1))) %>%
        layout(title = paste0("Distribution of acrophases (H", h, ")"),
               xaxis = list(title = x_lab, range = c(0, effective_period)),
               yaxis = list(title = "Count"))
    }
  })
  
  # Intercept by group plot (AUDIT 1.4: it is the constant, not the MESOR)
  output$harmonic_mesor_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = group, y = mesor, fill = group)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        labs(title = "Intercept (beta_0) by Group", x = "Group",
             y = "Intercept (beta_0, at t = 0)") +
        theme(legend.position = "none")
      
      ggplotly(g)
    } else {
      plot_ly(y = params$mesor, type = "box", 
              marker = list(color = 'forestgreen'),
              boxpoints = "all", jitter = 0.3) %>%
        layout(title = "Distribution of the intercept (beta_0)",
               yaxis = list(title = "Intercept (beta_0, at t = 0)"))
    }
  })
  
  # Trend parameter histogram
  output$harmonic_trend_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    # Check if trend is included
    if(mod$trend_type == "none") {
      plot_ly() %>% 
        layout(title = "No Homeostatic Trend Model Selected",
               annotations = list(
                 list(x = 0.5, y = 0.5, text = "Select a trend model\nto view distribution",
                      showarrow = FALSE, xref = "paper", yref = "paper",
                      font = list(size = 14, color = "gray"))))
    } else {
      params <- mod$individual_params
      
      # Get trend column based on type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        trend_vals <- params$trend_linear
        trend_label <- "Linear Trend (β, units/h)"
        trend_title <- "Distribution of Linear Trend Coefficient"
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        trend_vals <- params$trend_log
        trend_label <- "Log Trend (β)"
        trend_title <- "Distribution of Logarithmic Trend Coefficient"
      } else if(mod$trend_type == "exp_sat" && "A_sat" %in% names(params)) {
        trend_vals <- params$A_sat
        trend_label <- "A_sat (asymptote, units)"
        trend_title <- "Distribution of Saturating Exponential Asymptote"
      } else {
        trend_vals <- NULL
      }
      
      if(is.null(trend_vals)) {
        plot_ly() %>% layout(title = "Trend parameters not available")
      } else {
        # Check if groups exist
        if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
           !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
          
          group_var <- values$covariates[[input$harmonic_group_var]]
          params$group <- as.factor(group_var[params$subject])
          params$trend_val <- trend_vals
          params <- params[!is.na(params$group) & !is.na(params$trend_val), ]
          
          g <- ggplot(params, aes(x = group, y = trend_val, fill = group)) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_brewer(palette = "Set1") +
            theme_minimal() +
            labs(title = paste(trend_title, "by Group"), 
                 x = "Group", y = trend_label) +
            theme(legend.position = "none")
          
          ggplotly(g)
        } else {
          trend_vals <- trend_vals[!is.na(trend_vals)]
          plot_ly(y = trend_vals, type = "box", 
                  marker = list(color = 'steelblue'),
                  boxpoints = "all", jitter = 0.3) %>%
            layout(title = trend_title,
                   yaxis = list(title = trend_label))
        }
      }
    }
  })
  
  # Individual results table
  output$harmonic_individual_table <- DT::renderDataTable({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Create a cleaner display table
    display_df <- data.frame(
      Subject = params$subject,
      Intercept_b0 = round(params$mesor, 3),
      # AUDIT: which bounds this fit sits on, so the CSV carries the same
      # caveat the report shows rather than losing it on export.
      bounds_hit = vapply(params$subject, function(sid) {
        f <- mod$individual_fits[[sid]]
        if (is.null(f) || is.null(f$bounds_hit) || !length(f$bounds_hit)) ""
        else paste(f$bounds_hit, collapse = "; ")
      }, character(1)),
      n_bounds_hit = vapply(params$subject, function(sid) {
        f <- mod$individual_fits[[sid]]
        if (is.null(f) || is.null(f$bounds_hit)) 0L else length(f$bounds_hit)
      }, integer(1)),
      R_squared = round(params$r_squared, 3),
      Pct_Rhythm = round(params$percent_rhythm, 1),
      p_value = format(params$p_value, digits = 3, scientific = TRUE)
    )
    
    # Add trend parameters based on trend type
    if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
      display_df$Trend_Linear <- round(params$trend_linear, 4)
    } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
      display_df$Trend_Log <- round(params$trend_log, 4)
    } else if(mod$trend_type == "exp_sat") {
      if("A_sat" %in% names(params)) {
        display_df$A_sat <- round(params$A_sat, 3)
      }
      if("tau" %in% names(params)) {
        display_df$Tau_hrs <- round(params$tau, 2)
      }
    }
    
    # Add columns for each harmonic
    for(h in 1:mod$n_harmonics) {
      display_df[[paste0("Amp_H", h)]] <- round(params[[paste0("amplitude_", h)]], 3)
      display_df[[paste0("Acro_H", h, "_hrs")]] <- round(params[[paste0("acrophase_time_", h)]], 2)
    }
    
    # Add group if available
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      display_df$Group <- group_var[params$subject]
      # Move Group to second column
      display_df <- display_df[, c("Subject", "Group", setdiff(names(display_df), c("Subject", "Group")))]
    }
    
    DT::datatable(display_df, 
                  options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatStyle('R_squared', 
                      backgroundColor = DT::styleInterval(c(0.5, 0.8), c('#ffcccc', '#ffffcc', '#ccffcc')))
  })
  
  # Export individual parameters
  output$export_harmonic_individual <- downloadHandler(
    filename = function() paste0("harmonic_individual_params_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$harmonic_model)
      mod <- values$harmonic_model
      params <- mod$individual_params
      
      # Add group variable if available
      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]
        params$group <- group_var[params$subject]
      }
      
      write.csv(params, file, row.names = FALSE)
    }
  )
  
  # Residual plot
  output$harmonic_residual_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    # Collect all residuals
    all_fitted <- c()
    all_resid <- c()
    
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_fitted <- c(all_fitted, fit_i$fitted)
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    plot_ly(x = all_fitted, y = all_resid, type = 'scatter', mode = 'markers',
            marker = list(color = 'steelblue', opacity = 0.5, size = 4)) %>%
      add_segments(x = min(all_fitted), xend = max(all_fitted), y = 0, yend = 0,
                   line = list(color = 'red', dash = 'dash')) %>%
      layout(title = "Residuals vs Fitted",
             xaxis = list(title = "Fitted Values"),
             yaxis = list(title = "Residuals"))
  })
  
  # QQ plot
  output$harmonic_qq_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    all_resid <- c()
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    qq <- qqnorm(all_resid, plot.it = FALSE)
    
    plot_ly(x = qq$x, y = qq$y, type = 'scatter', mode = 'markers',
            marker = list(color = 'steelblue', size = 4)) %>%
      add_lines(x = range(qq$x), y = range(qq$x) * sd(all_resid) + mean(all_resid),
                line = list(color = 'red', dash = 'dash'), name = "Reference") %>%
      layout(title = "Q-Q Plot of Residuals",
             xaxis = list(title = "Theoretical Quantiles"),
             yaxis = list(title = "Sample Quantiles"))
  })
  
  # Goodness of fit statistics
  output$harmonic_gof_stats <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    all_resid <- c()
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    cat("=== Residual Diagnostics ===\n\n")
    cat(sprintf("Total residuals: %d\n", length(all_resid)))
    cat(sprintf("Mean residual: %.4f\n", mean(all_resid)))
    cat(sprintf("SD of residuals: %.4f\n", sd(all_resid)))
    cat(sprintf("Skewness: %.3f\n", mean((all_resid - mean(all_resid))^3) / sd(all_resid)^3))
    cat(sprintf("Kurtosis: %.3f\n", mean((all_resid - mean(all_resid))^4) / sd(all_resid)^4 - 3))
    
    # Shapiro-Wilk test (on sample if too many observations)
    if(length(all_resid) > 5000) {
      samp_resid <- sample(all_resid, 5000)
    } else {
      samp_resid <- all_resid
    }
    sw <- shapiro.test(samp_resid)
    cat(sprintf("\nShapiro-Wilk test: W = %.4f, p = %.4f\n", sw$statistic, sw$p.value))
    if(sw$p.value < 0.05) {
      cat("  (Significant departure from normality)\n")
    } else {
      cat("  (No significant departure from normality)\n")
    }
  })
  
  # Group comparison plot
  output$harmonic_group_comparison_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      return(plotly_empty() %>% layout(title = "Select a group variable for comparison"))
    }
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_group)) as.integer(input$selected_harmonic_group) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    acro_col <- paste0("acrophase_time_", h)
    acro_rad_col <- paste0("acrophase_rad_", h)
    effective_period <- mod$period / h
    
    # Get individual params with group info for R-squared
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove missings
    }
    
    # Create comparison data
    group_df <- data.frame(
      group = character(),
      parameter = character(),
      value = numeric(),
      se = numeric()
    )
    
    for(g_name in names(mod$group_fits)) {
      g <- mod$group_fits[[g_name]]
      
      # AUDIT: R-squared was here, which is a fit-quality number rather than a
      # parameter -- it says how well each subject's curve was described, not
      # what the group's rhythm is, so it does not belong in a panel of
      # parameter comparisons. Replaced by the MESOR, which sits naturally
      # beside beta_0: beta_0 is the fitted constant, the MESOR is the
      # rhythm-adjusted mean over the window, and a group can rank differently
      # on the two.
      grp_mesor <- params$mesor_adj[params$group == g_name]
      grp_mesor <- grp_mesor[is.finite(grp_mesor)]
      mesor_mean <- if(length(grp_mesor)) mean(grp_mesor) else NA_real_
      mesor_se <- if(length(grp_mesor) > 1)
        sd(grp_mesor) / sqrt(length(grp_mesor)) else NA_real_
      
      # Get amplitude and acrophase for selected harmonic from individual params
      grp_amp <- params[[amp_col]][params$group == g_name]
      grp_acro_rad <- params[[acro_rad_col]][params$group == g_name]
      
      # A group with no usable angles has no circular mean: circular_mean() of
      # an empty vector is atan2(NaN, NaN) = NaN, and `if (NaN < 0)` is the
      # error that crashed this plot. Unlabelled subjects no longer reach here,
      # but a group can still be emptied by a filter upstream, so the guard
      # stays and the group is skipped rather than poisoning the frame.
      grp_acro_rad <- grp_acro_rad[is.finite(grp_acro_rad)]
      if(length(grp_acro_rad) < 1) next
      circ_mean_rad <- circular_mean(grp_acro_rad)
      if(!is.finite(circ_mean_rad)) next
      if(circ_mean_rad < 0) circ_mean_rad <- circ_mean_rad + 2 * pi
      # reported in CLOCK time, like every other acrophase in the app
      circ_mean_time <- (phi_to_hours(circ_mean_rad, mod$period, h) +
                           fck_clock_origin(mod)) %% effective_period
      circ_se_rad <- circular_se(grp_acro_rad)
      circ_se_time <- if(!is.na(circ_se_rad)) phi_to_hours(circ_se_rad, mod$period, h) else NA
      
      group_df <- rbind(group_df, 
                        data.frame(group = g_name, parameter = "Constant term (b0)",
                                   value = g$mean_mesor, se = g$sd_mesor / sqrt(g$n)),
                        data.frame(group = g_name, parameter = paste0("Amplitude (H", h, ")"), 
                                   value = mean(grp_amp, na.rm = TRUE), 
                                   se = sd(grp_amp, na.rm = TRUE) / sqrt(length(grp_amp))),
                        data.frame(group = g_name, parameter = paste0("Acrophase (H", h, ")"), 
                                   value = circ_mean_time, se = circ_se_time),
                        data.frame(group = g_name, parameter = "MESOR (rhythm-adjusted)",
                                   value = mesor_mean, se = mesor_se))
      
      # Add trend parameters based on trend type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        grp_trend <- params$trend_linear[params$group == g_name]
        grp_trend <- grp_trend[!is.na(grp_trend)]
        if(length(grp_trend) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "Linear Trend (β)", 
                                       value = mean(grp_trend, na.rm = TRUE), 
                                       se = sd(grp_trend, na.rm = TRUE) / sqrt(length(grp_trend))))
        }
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        grp_trend <- params$trend_log[params$group == g_name]
        grp_trend <- grp_trend[!is.na(grp_trend)]
        if(length(grp_trend) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "Log Trend (β)", 
                                       value = mean(grp_trend, na.rm = TRUE), 
                                       se = sd(grp_trend, na.rm = TRUE) / sqrt(length(grp_trend))))
        }
      } else if(mod$trend_type == "exp_sat") {
        if("A_sat" %in% names(params)) {
          grp_asat <- params$A_sat[params$group == g_name]
          grp_asat <- grp_asat[!is.na(grp_asat)]
          if(length(grp_asat) > 0) {
            group_df <- rbind(group_df,
                              data.frame(group = g_name, parameter = "A_sat (asymptote)", 
                                         value = mean(grp_asat, na.rm = TRUE), 
                                         se = sd(grp_asat, na.rm = TRUE) / sqrt(length(grp_asat))))
          }
        }
        if("tau" %in% names(params)) {
          grp_tau <- params$tau[params$group == g_name]
          grp_tau <- grp_tau[!is.na(grp_tau)]
          if(length(grp_tau) > 0) {
            group_df <- rbind(group_df,
                              data.frame(group = g_name, parameter = "τ (time constant, h)", 
                                         value = mean(grp_tau, na.rm = TRUE), 
                                         se = sd(grp_tau, na.rm = TRUE) / sqrt(length(grp_tau))))
          }
        }
      }

      # Add sleep inertia parameters if present
      include_inertia <- !is.null(mod$include_inertia) && isTRUE(mod$include_inertia)
      if(include_inertia && "W0" %in% names(params)) {
        grp_W0 <- params$W0[params$group == g_name]
        grp_W0 <- grp_W0[!is.na(grp_W0)]
        if(length(grp_W0) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "W₀ (inertia)",
                                       value = mean(grp_W0, na.rm = TRUE),
                                       se = sd(grp_W0, na.rm = TRUE) / sqrt(length(grp_W0))))
        }
      }
      if(include_inertia && "tau_W" %in% names(params)) {
        grp_tau_W <- params$tau_W[params$group == g_name]
        grp_tau_W <- grp_tau_W[!is.na(grp_tau_W)]
        if(length(grp_tau_W) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "τ_W (decay, h)",
                                       value = mean(grp_tau_W, na.rm = TRUE),
                                       se = sd(grp_tau_W, na.rm = TRUE) / sqrt(length(grp_tau_W))))
        }
      }
    }

    # The shared palette, keyed by group name, so this figure agrees with the
    # fitted curves, the polar plots and the pairwise boxplots.
    group_names <- names(mod$group_fits)
    named_colors <- fck_group_colors(group_names)

    # Faceted bar plot
    g <- ggplot(group_df, aes(x = group, y = value, fill = group)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = value - 1.96*se, ymax = value + 1.96*se),
                    width = 0.2, na.rm = TRUE) +
      facet_wrap(~parameter, scales = "free_y") +
      scale_fill_manual(values = named_colors) +
      theme_minimal() +
      labs(title = paste0("Group Comparison - Harmonic ", h, " (period = ", round(effective_period, 1), "h)"),
           x = "", y = "Value") +
      theme(legend.position = "none")
    
    ggplotly(g)
  })
  
  # Group comparison test results
  output$harmonic_group_test_results <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      cat("Select a group variable for statistical comparison.\n")
      return()
    }
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_group)) as.integer(input$selected_harmonic_group) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    acro_col <- paste0("acrophase_time_", h)
    effective_period <- mod$period / h
    
    cat("=== Group Comparison Statistics ===\n")
    cat(sprintf("Harmonic: H%d (effective period = %s h)\n", h, fmt1(effective_period)))
    cat(sprintf("Acrophases below are CLOCK times (model origin %s).\n",
                fck_clock_label(fck_clock_origin(mod), mod$period, show_day = FALSE)))
    if(h > 1)
      cat(sprintf("H%d repeats every %s h, so each acrophase has %d equivalent clock times.\n",
                  h, fmtn(effective_period, 0), h))
    cat(sprintf("Data source: %s.%s\n\n",
                if(isTRUE(mod$using_smoothed)) "SMOOTHED" else "RAW",
                if(isTRUE(mod$using_smoothed))
                  " Between-group differences are less affected by smoothing than the within-subject fit statistics, but the per-subject parameters entering these tests are still smoothed-data estimates."
                else ""))
    
    # Get group variable and params
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- group_var[params$subject]
      
      # Remove rows with missing group or key parameters
      params <- params[!is.na(params$group) & !is.na(params$mesor) & 
                         !is.na(params[[amp_col]]) & !is.na(params[[acro_col]]), ]
      
      # ======================================================================
      # BUG FIX: "grouping factor must have exactly 2 levels"
      #
      # n_groups was computed ONCE here, on the full params frame, then reused
      # to choose t-test vs ANOVA inside blocks that had since taken a SUBSET:
      #
      #     params_trend <- params[!is.na(params$A_sat), ]
      #     if (n_groups == 2) t.test(A_sat ~ group, data = params_trend)
      #
      # If that filter leaves only one group with usable values, the branch and
      # the data disagree and t.test stops. n_groups says 2; the subset does not.
      #
      # It became reachable when the convergence gate started excluding
      # non-converged fits (audit 2.3): on the Circaflex data that removes 938
      # of 1305 subjects, and a sparse group can lose every member that had a
      # usable trend parameter. The latent bug is older; that change exposed it.
      #
      # NOTE on a wrong first guess, recorded so nobody re-derives it: this is
      # NOT about a factor keeping unused levels. t.test.formula() calls
      # factor() on the grouping column, which drops unused levels itself, so an
      # undropped 4-level factor holding 2 values works fine (verified on R
      # 4.3.3). droplevels() below is still correct hygiene -- it makes
      # nlevels() and length(unique()) agree downstream, and length(unique())
      # also counts NA as a group -- but it is not what fixes the crash.
      #
      # The fix is .group_test(): it reads the number of groups off the data it
      # is ACTUALLY handed, so the branch and the test can never disagree again,
      # whatever filter ran in between.
      # ======================================================================
      params$group <- droplevels(as.factor(params$group))
      n_groups <- nlevels(params$group)

      if(n_groups < 2) {
        cat(sprintf("Groups: %d, Total N: %d (after removing missings)\n\n", n_groups, nrow(params)))
        cat("Fewer than two groups remain after dropping subjects with a missing\n")
        cat("group label or a missing parameter", if(isTRUE(mod$fit_audit$n_nonconverged > 0) ||
              isTRUE(mod$fit_audit$n_boundary > 0))
              ", and after the convergence gate excluded\nnon-converged fits" else "",
            ". There is nothing to compare.\n", sep = "")
        if(!is.null(mod$fit_audit))
          cat(sprintf("\n  Of %d subjects: %d converged, %d were pinned to a parameter bound,\n  %d did not converge, %d failed outright.\n",
                      mod$fit_audit$n_attempted, mod$fit_audit$n_converged,
                      mod$fit_audit$n_boundary, mod$fit_audit$n_nonconverged,
                      mod$fit_audit$n_failed))
        return(invisible(NULL))
      }
      cat(sprintf("Groups: %d, Total N: %d (after removing missings)\n", n_groups, nrow(params)))

      # One entry point for every two-or-more-group comparison in this block.
      # It decides t-test vs ANOVA from the data in front of it, after dropping
      # unused levels, so the decision and the test can never disagree.
      .group_test <- function(formula, data, label = NULL) {
        gv <- droplevels(as.factor(data[[all.vars(formula)[2]]]))
        data[[all.vars(formula)[2]]] <- gv
        yv <- data[[all.vars(formula)[1]]]
        keep <- !is.na(gv) & is.finite(yv)
        data <- data[keep, , drop = FALSE]
        data[[all.vars(formula)[2]]] <- droplevels(data[[all.vars(formula)[2]]])
        k <- nlevels(data[[all.vars(formula)[2]]])
        if(k < 2) {
          cat(sprintf("  Only %d group%s usable values here; nothing to compare.\n",
                      k, if(k == 1) " has" else "s have"))
          return(invisible(NULL))
        }
        if(any(table(data[[all.vars(formula)[2]]]) < 2)) {
          cat("  At least one group has fewer than 2 usable values; no test is run.\n")
          return(invisible(NULL))
        }
        if(k == 2) {
          tt <- tryCatch(stats::t.test(formula, data = data), error = function(e) NULL)
          if(is.null(tt)) { cat("  The t-test could not be computed.\n"); return(invisible(NULL)) }
          # Welch by default (t.test's own default), which is right for the
          # unbalanced, unequal-variance groups this app routinely produces.
          cat(sprintf("Welch t-test: t = %s, df = %s, p = %s\n", fmt3(tt$statistic),
                      fmt1(tt$parameter), format.pval(tt$p.value, digits = 3, eps = 1e-16)))
        } else {
          av <- tryCatch(stats::aov(formula, data = data), error = function(e) NULL)
          if(is.null(av)) { cat("  The ANOVA could not be computed.\n"); return(invisible(NULL)) }
          cat("ANOVA:\n"); print(summary(av))
        }
        invisible(NULL)
      }

      # AUDIT 1.5: this block filters NA group labels while the Group-Specific
      # Parameters panel used to drop them silently, so the two N's disagreed by
      # construction and neither reconciled to the number fitted. Both now say
      # what they did.
      n_dropped_here <- nrow(mod$individual_params) - nrow(params)
      if(n_dropped_here > 0)
        cat(sprintf("  (%d fitted subject(s) excluded here for a missing group label or parameter.)\n",
                    n_dropped_here))
      cat("\n")

      # ======================================================================
      # AUDIT 1.4 (applied to the comparisons too, as requested)
      #
      # These tests are on the fitted CONSTANT, not on a MESOR. Under a
      # saturating trend the constant is the intercept of a model with two
      # origins; a group difference in it is not a group difference in level.
      # Both are now tested, and the second is the one to interpret.
      # ======================================================================
      cat("--- Intercept (beta_0, at t = 0) comparison ---\n")
      cat("NOT a MESOR comparison. See the rhythm-adjusted mean below for level.\n")
      lt <- fck_group_linear_test(params$mesor, params$group)
      if(!is.null(lt)) {
        cat(sprintf("  F(%d, %d) = %s, p = %s\n", lt$df1, lt$df2, fmt3(lt$F),
                    format.pval(lt$p, digits = 3, eps = 1e-16)))
        cat(sprintf("  eta^2 = %s, omega^2 = %s   [omega^2 is the less optimistic; quote it]\n",
                    fmt3(lt$eta2), fmt3(lt$omega2)))
        for(gg in lt$levels)
          cat(sprintf("    %-14s mean %s (SD %s, n = %d)\n", gg,
                      fmt3(lt$means[gg]), fmt3(lt$sds[gg]), lt$ns[gg]))
        if(!is.null(lt$largest))
          cat(sprintf("  Largest contrast %s - %s: %s, 95%% CI [%s, %s], Hedges' g = %s\n",
                      lt$largest$a, lt$largest$b, fmt3(lt$largest$diff),
                      fmt3(lt$largest$ci[1]), fmt3(lt$largest$ci[2]),
                      fmt3(lt$largest$hedges_g)))
      }

      # the quantity that actually means "level": the rhythm-adjusted mean
      if(!is.null(mod$group_fits) && (mod$trend_type %||% "none") != "none") {
        cat("\n--- MESOR (rhythm-adjusted mean over the observed window) by group ---\n")
        cat("Integrated per group from that group's own intercept and trend parameters.\n")
        for(gn in names(mod$group_fits)) {
          gf <- mod$group_fits[[gn]]
          if(isTRUE(gf$is_unassigned)) next
          cat(sprintf("    %-14s intercept %s -> MESOR %s   (difference %s)\n", gn,
                      fmt3(gf$intercept), fmt3(gf$rhythm_adjusted_mean),
                      fmt3(gf$rhythm_adjusted_mean - gf$intercept)))
        }
        cat("  A group can rank differently on these two: the intercept absorbs the\n")
        cat("  trend baseline, the rhythm-adjusted mean does not.\n")
      }
      
      # Trend comparison (if trend is present)
      trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
      if(trend_type != "none") {
        trend_col <- switch(trend_type,
                            "linear" = "trend_linear",
                            "log" = "trend_log",
                            "exp_sat" = "A_sat",
                            NULL)
        
        if(!is.null(trend_col) && trend_col %in% names(params)) {
          cat(sprintf("\n--- %s Comparison ---\n", get_trend_label(trend_type)))
          trend_formula <- as.formula(paste(trend_col, "~ group"))
          .group_test(trend_formula, params)
          
          # AUDIT 2.6: effect size and interval, not a p-value alone
          .lt <- fck_group_linear_test(params[[trend_col]], params$group)
          if(!is.null(.lt)) {
            cat(sprintf("  eta^2 = %s, omega^2 = %s\n", fmt3(.lt$eta2), fmt3(.lt$omega2)))
            if(!is.null(.lt$largest))
              cat(sprintf("  Largest contrast %s - %s: %s, 95%% CI [%s, %s], Hedges' g = %s\n",
                          .lt$largest$a, .lt$largest$b, fmt3(.lt$largest$diff),
                          fmt3(.lt$largest$ci[1]), fmt3(.lt$largest$ci[2]),
                          fmt3(.lt$largest$hedges_g)))
          }
          # Group means for trend
          for(g in unique(params$group)) {
            g_trend <- params[[trend_col]][params$group == g]
            cat(sprintf("  %s: mean = %s (SD = %s, n=%d)\n",
                        g, fmt4(mean(g_trend, na.rm = TRUE)),
                        fmt4(sd(g_trend, na.rm = TRUE)), length(g_trend)))
          }
          # AUDIT 2.6: "declines monotonically across age bands" is a claim about
          # ORDER. An omnibus ANOVA does not test it; a linear contrast on the
          # ordered levels does, and is much more powerful against exactly that
          # alternative. Only meaningful if the factor levels are in order.
          .tr <- fck_group_trend_test(params[[trend_col]], params$group)
          if(!is.null(.tr)) {
            cat(sprintf("  Monotone trend across ordered levels (%s):\n",
                        paste(.tr$levels, collapse = " < ")))
            cat(sprintf("    L = %s (SE %s), t(%d) = %s, p = %s, 95%% CI [%s, %s]\n",
                        fmt3(.tr$L), fmt3(.tr$se), .tr$df, fmt3(.tr$t),
                        format.pval(.tr$p, digits = 3, eps = 1e-16),
                        fmt3(.tr$ci[1]), fmt3(.tr$ci[2])))
            cat("    Valid only if the group levels are genuinely ordered as printed.\n")
          }
          
          # For exp_sat, also compare tau
          if(trend_type == "exp_sat" && "tau" %in% names(params)) {
            cat("\n--- Time Constant (τ) Comparison ---\n")
            tau_formula <- as.formula("tau ~ group")
            .group_test(tau_formula, params)
            
            for(g in unique(params$group)) {
              g_tau <- params$tau[params$group == g]
              cat(sprintf("  %s: mean tau = %s h (SD = %s, n=%d)\n",
                          g, fmt2(mean(g_tau, na.rm = TRUE)),
                          fmt2(sd(g_tau, na.rm = TRUE)), length(g_tau)))
            }
            .lt <- fck_group_linear_test(params$tau, params$group)
            if(!is.null(.lt)) {
              cat(sprintf("  eta^2 = %s, omega^2 = %s\n", fmt3(.lt$eta2), fmt3(.lt$omega2)))
              if(!is.null(.lt$largest))
                cat(sprintf("  Largest contrast %s - %s: %s h, 95%% CI [%s, %s], Hedges' g = %s\n",
                            .lt$largest$a, .lt$largest$b, fmt2(.lt$largest$diff),
                            fmt2(.lt$largest$ci[1]), fmt2(.lt$largest$ci[2]),
                            fmt3(.lt$largest$hedges_g)))
            }
            .tr <- fck_group_trend_test(params$tau, params$group)
            if(!is.null(.tr))
              cat(sprintf("  Monotone trend: t(%d) = %s, p = %s, 95%% CI [%s, %s]\n",
                          .tr$df, fmt3(.tr$t), format.pval(.tr$p, digits = 3, eps = 1e-16),
                          fmt3(.tr$ci[1]), fmt3(.tr$ci[2])))
            # AUDIT 2.2: a between-group comparison of tau is only meaningful if
            # tau is identified WITHIN subject. It is not, here.
            if(!is.null(mod$conditioning) && !is.null(mod$conditioning$tau_fixed_delta_aic) &&
               mod$conditioning$tau_fixed_delta_aic > 0)
              cat("  ! Delta-AIC says free tau is not better than tau held fixed: tau is not\n",
                  "    identified by these data. A group difference in an unidentified\n",
                  "    parameter is a difference in where the optimiser stopped on the ridge.\n",
                  "    Do not report this comparison without the conditioning evidence.\n", sep = "")
          }
        }
      }

      cat(sprintf("\n--- Amplitude (H%d) Comparison ---\n", h))
      amp_formula <- as.formula(paste(amp_col, "~ group"))
      .group_test(amp_formula, params)
      # AUDIT 2.6: effect sizes and intervals for the amplitude decline that the
      # report described (25.7 -> 22.6 -> 22.5 -> 21.1) and never tested.
      .lt <- fck_group_linear_test(params[[amp_col]], params$group)
      if(!is.null(.lt)) {
        cat(sprintf("  eta^2 = %s, omega^2 = %s\n", fmt3(.lt$eta2), fmt3(.lt$omega2)))
        for(gg in .lt$levels)
          cat(sprintf("    %-14s mean %s (SD %s, n = %d)   [arithmetic]\n", gg,
                      fmt3(.lt$means[gg]), fmt3(.lt$sds[gg]), .lt$ns[gg]))
        if(!is.null(.lt$largest))
          cat(sprintf("  Largest contrast %s - %s: %s, 95%% CI [%s, %s], Hedges' g = %s\n",
                      .lt$largest$a, .lt$largest$b, fmt3(.lt$largest$diff),
                      fmt3(.lt$largest$ci[1]), fmt3(.lt$largest$ci[2]),
                      fmt3(.lt$largest$hedges_g)))
      }
      .tr <- fck_group_trend_test(params[[amp_col]], params$group)
      if(!is.null(.tr))
        cat(sprintf("  Monotone trend across ordered levels: t(%d) = %s, p = %s, 95%% CI [%s, %s]\n",
                    .tr$df, fmt3(.tr$t), format.pval(.tr$p, digits = 3, eps = 1e-16),
                    fmt3(.tr$ci[1]), fmt3(.tr$ci[2])))
      cat("  Note: these are the ARITHMETIC group means. The Group-Specific Parameters\n")
      cat("  panel reports VECTOR means for the same amplitudes; they are different\n")
      cat("  estimators and will not match. The test above is on the arithmetic values.\n")
      
      cat(sprintf("\n--- Acrophase (H%d) Comparison ---\n", h))
      cat(strrep("-", 60), "\n")
      cat("TWO COMPLEMENTARY TESTS ARE PROVIDED:\n")
      cat("  (1) Watson-Williams: tests timing (phase) only.\n")
      cat("      Every subject counts equally regardless of rhythm strength.\n")
      cat("      Use this to answer: 'Do groups peak at different times?'\n")
      cat("  (2) Hotelling's T² on (beta_cos, beta_sin): tests the full\n")
      cat("      rhythmic vector (phase + amplitude combined). Subjects with\n")
      cat("      stronger rhythms carry more weight (amplitude-weighted).\n")
      cat("      The group means it reports match the 'Group-Specific\n")
      cat("      Parameters' panel exactly.\n")
      cat("      Use this to answer: 'Do groups differ in their overall\n")
      cat("      rhythmic profile?' — but note a significant result could\n")
      cat("      reflect phase, amplitude, or both.\n")
      cat(strrep("-", 60), "\n\n")

      # Get acrophase in radians for each group
      acro_rad_col <- paste0("acrophase_rad_", h)
      beta_cos_col <- paste0("beta_cos_", h)
      beta_sin_col <- paste0("beta_sin_", h)
      groups <- unique(params$group)
      angles_list <- lapply(groups, function(g) {
        params[[acro_rad_col]][params$group == g]
      })
      names(angles_list) <- groups

      # (1) Watson-Williams test
      cat("(1) Watson-Williams test (unweighted circular mean):\n")
      ww <- watson_williams_test(angles_list)

      if(!is.null(ww$message)) {
        cat(sprintf("  %s\n", ww$message))
      }

      # AUDIT 2.6: the concentration assumption is CHECKED and the check is
      # printed. It was previously enforced silently inside the function and
      # summarised as an adjective.
      .as <- fck_ww_assumption(ww$r_bar, ww$kappa)
      cat(sprintf("    Assumption check: %s\n", .as$msg))
      if(!is.na(ww$F)) {
        cat(sprintf("    F(%d, %d) = %s, p = %s\n", ww$df1, ww$df2, fmt3(ww$F),
                    format.pval(ww$p, digits = 3, eps = 1e-16)))
        cat(sprintf("    r-bar (unweighted, pooled) = %s\n", fmt3(ww$r_bar)))
        if(!.as$ok)
          cat("    The F value above should NOT be reported: the assumption it rests on\n",
              "    does not hold for these data.\n", sep = "")
      }

      cat("\n  Group circular statistics (unweighted):\n")
      for(g in groups) {
        g_angles <- params[[acro_rad_col]][params$group == g]
        g_mean <- circular_mean(g_angles)
        if(g_mean < 0) g_mean <- g_mean + 2 * pi
        g_mean_time <- g_mean * effective_period / (2 * pi)
        g_sd <- circular_sd(g_angles)
        g_sd_time <- if(!is.na(g_sd)) g_sd * effective_period / (2 * pi) else NA
        g_r <- mean_resultant_length(g_angles)
        # AUDIT 1.2: label which resultant this is. Both are shown.
        g_rw <- fck_resultants(g_angles, params[[amp_col]][params$group == g])
        # AUDIT: circular MEAN is a direction and moves with the origin; circular
        # SD is a dispersion and does not.
        cat(sprintf("    %-14s circular mean %s, circ.SD %s h, r-bar unweighted %s, weighted %s (n=%d)\n",
                    g, fck_acrophase_label(hours = g_mean_time, period = mod$period,
                                           harmonic = h, clock_origin = fck_clock_origin(mod)),
                    fmt2(g_sd_time), fmt3(g_r),
                    fmt3(if(is.null(g_rw)) NA_real_ else g_rw$r_weighted), length(g_angles)))
      }

      # (2) Hotelling's T² test
      # AUDIT 2.6: this IS Bingham's parameter test for amplitude-acrophase pairs
      # (Bingham et al. 1982) -- the joint test on the (beta_cos, beta_sin)
      # vector. Naming it as such connects the output to the literature the
      # brief asked for, rather than leaving it as an unattributed T-squared.
      cat(sprintf("\n(2) Bingham parameter test / Hotelling's T-squared on (beta_cos_%d, beta_sin_%d):\n", h, h))
      cat("    The joint test on the rhythmic vector -- amplitude and acrophase together\n")
      cat("    (Bingham, Arbogast, Cornelissen Guillaume, Lee & Halberg 1982).\n")
      bc_list <- lapply(groups, function(g) params[[beta_cos_col]][params$group == g])
      bs_list <- lapply(groups, function(g) params[[beta_sin_col]][params$group == g])
      ht <- hotelling_t2(bc_list, bs_list)

      if(!is.null(ht$message)) {
        cat(sprintf("  %s\n", ht$message))
      } else {
        cat(sprintf("    F(%d, %d) = %.3f, p = %.4f\n", ht$df1, ht$df2, ht$F, ht$p))
      }

      cat("\n  Group means (amplitude-weighted) — match 'Group-Specific Parameters' panel:\n")
      for(g in groups) {
        bc <- params[[beta_cos_col]][params$group == g]
        bs <- params[[beta_sin_col]][params$group == g]
        ok <- complete.cases(bc, bs)
        bc_ok <- bc[ok]; bs_ok <- bs[ok]
        # Amplitude-weighted circular mean via vector averaging
        x_h <- sqrt(bc_ok^2 + bs_ok^2) * cos(atan2(bs_ok, bc_ok))
        y_h <- sqrt(bc_ok^2 + bs_ok^2) * sin(atan2(bs_ok, bc_ok))
        acro_h <- atan2(mean(y_h), mean(x_h))
        if(acro_h < 0) acro_h <- acro_h + 2 * pi
        # AUDIT: acro_h is a model-frame direction; report it as a clock time
        # through the one helper, like every other acrophase in the app.
        amp_mean <- sqrt(mean(bc_ok)^2 + mean(bs_ok)^2)
        cat(sprintf("    %s: amplitude-weighted mean = %s, mean amplitude = %s (n=%d)\n",
                    g, fck_acrophase_label(phi_rad = acro_h, period = mod$period,
                                           harmonic = h,
                                           clock_origin = fck_clock_origin(mod)),
                    fmt3(amp_mean), sum(ok)))
      }
      
      cat("\n--- R-squared Comparison ---\n")
      .group_test(r_squared ~ group, params)
      
      # Trend parameter comparison based on trend type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        cat("\n--- Linear Trend (β) Comparison ---\n")
        params_trend <- params[!is.na(params$trend_linear), ]
        
        if(nrow(params_trend) >= 4) {
          .group_test(trend_linear ~ group, params_trend)
          # Effect size and interval come from fck_group_linear_test() below,
          # which handles any number of groups; the hand-rolled two-group
          # Cohen's d that used to sit here only ran in the k = 2 branch.
          .lt <- fck_group_linear_test(params_trend[[all.vars(trend_linear ~ group)[1]]],
                                       params_trend$group)
          if(!is.null(.lt)) {
            cat(sprintf("  eta^2 = %s, omega^2 = %s\n", fmt3(.lt$eta2), fmt3(.lt$omega2)))
            if(!is.null(.lt$largest))
              cat(sprintf("  Largest contrast %s - %s: %s, 95%%%% CI [%s, %s], Hedges' g = %s\n",
                          .lt$largest$a, .lt$largest$b, fmt3(.lt$largest$diff),
                          fmt3(.lt$largest$ci[1]), fmt3(.lt$largest$ci[2]),
                          fmt3(.lt$largest$hedges_g)))
          }
          
          cat("\nGroup statistics (units/hour):\n")
          for(g in unique(params_trend$group)) {
            g_trend <- params_trend$trend_linear[params_trend$group == g]
            cat(sprintf("  %s: mean = %.4f, SD = %.4f (n=%d)\n", 
                        g, mean(g_trend, na.rm=TRUE), sd(g_trend, na.rm=TRUE), length(g_trend)))
          }
        }
        
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        cat("\n--- Logarithmic Trend (β) Comparison ---\n")
        params_trend <- params[!is.na(params$trend_log), ]
        
        if(nrow(params_trend) >= 4) {
          .group_test(trend_log ~ group, params_trend)
          
          cat("\nGroup statistics (units/log-hour):\n")
          for(g in unique(params_trend$group)) {
            g_trend <- params_trend$trend_log[params_trend$group == g]
            cat(sprintf("  %s: mean = %.4f, SD = %.4f (n=%d)\n", 
                        g, mean(g_trend, na.rm=TRUE), sd(g_trend, na.rm=TRUE), length(g_trend)))
          }
        }
        
      } else if(mod$trend_type == "exp_sat") {
        # Compare A_sat (asymptote)
        if("A_sat" %in% names(params)) {
          cat("\n--- Asymptote (A_sat) Comparison ---\n")
          params_trend <- params[!is.na(params$A_sat), ]
          
          if(nrow(params_trend) >= 4) {
            .group_test(A_sat ~ group, params_trend)
            
            # AUDIT 1.7: "(units)" was a placeholder that was never interpolated.
            cat(sprintf("\nGroup statistics%s:\n",
                        if(!is.null(mod$dv_units)) paste0(" (", mod$dv_units, ")") else ""))
            for(g in unique(params_trend$group)) {
              g_asat <- params_trend$A_sat[params_trend$group == g]
              cat(sprintf("  %s: mean = %s, SD = %s (n=%d)\n",
                          g, fmt3(mean(g_asat, na.rm=TRUE)),
                          fmt3(sd(g_asat, na.rm=TRUE)), length(g_asat)))
            }
          }
        }
        
        # Compare tau (time constant)
        if("tau" %in% names(params)) {
          cat("\n--- Time Constant (τ) Comparison ---\n")
          params_trend <- params[!is.na(params$tau), ]
          
          if(nrow(params_trend) >= 4) {
            .group_test(tau ~ group, params_trend)
            
            cat("\nGroup statistics (hours):\n")
            for(g in unique(params_trend$group)) {
              g_tau <- params_trend$tau[params_trend$group == g]
              cat(sprintf("  %s: mean τ = %.2f h, SD = %.2f (n=%d)\n", 
                          g, mean(g_tau, na.rm=TRUE), sd(g_tau, na.rm=TRUE), length(g_tau)))
            }
            cat("\nNote: tau represents time to reach ~63% of asymptotic level.\n")
          }
        }
      }

      # ====================================================================
      # AUDIT 2.6: the population-mean cosinor with group x harmonic terms
      #
      # Everything above is a two-stage analysis: fit each subject, then run a
      # t-test or ANOVA on the per-subject estimates. That treats every subject
      # as contributing one equally-precise observation, which they do not --
      # a subject with 16 clean points and a subject on the edge of
      # identifiability get the same weight -- and it cannot borrow strength
      # across an unbalanced design (654 / 410 / 181 / 59 here).
      #
      # The population-mean cosinor (Cornelissen 2014) fits ONE model to all
      # subjects' observations at once, with group x (cos, sin) interactions, so
      # the group contrast is a linear hypothesis inside a single fit. With a
      # random intercept per subject it is a mixed model; without lme4 available
      # it degrades to a fixed-effects fit with cluster-robust standard errors,
      # which is stated rather than hidden.
      # ====================================================================
      cat("\n\n=== Population-mean cosinor (primary analysis) ===\n")
      cat("One model over all observations, group x harmonic interactions.\n")
      cat("Preferred over the per-subject tests above: it weights subjects by their\n")
      cat("actual precision and handles the unbalanced design directly.\n\n")

      pmc <- tryCatch({
        Yp <- if(isTRUE(mod$using_smoothed)) values$smooth_data else values$data
        tv <- mod$time_vec
        sub_ids <- mod$individual_params$subject
        gl <- group_var[sub_ids]
        keep <- !is.na(gl)
        sub_ids <- sub_ids[keep]; gl <- droplevels(as.factor(gl[keep]))
        long <- data.frame(
          y = as.vector(t(Yp[sub_ids, , drop = FALSE])),
          t = rep(tv, times = length(sub_ids)),
          subj = factor(rep(sub_ids, each = length(tv))),
          grp = rep(gl, each = length(tv)))
        long <- long[is.finite(long$y), ]
        for(hh in seq_len(mod$n_harmonics)) {
          w <- 2 * pi * hh / mod$period
          long[[paste0("c", hh)]] <- cos(w * long$t)
          long[[paste0("s", hh)]] <- sin(w * long$t)
        }
        trm <- switch(as.character(mod$trend_type),
                      "linear" = "t",
                      "log" = "I(log(t - min(t) + 1))",
                      "exp_sat" = sprintf("I(1 - exp(-(t - min(t))/%f))",
                                          mean(mod$individual_params$tau, na.rm = TRUE)),
                      NULL)
        harm <- paste(unlist(lapply(seq_len(mod$n_harmonics),
                                    function(hh) c(paste0("c", hh), paste0("s", hh)))),
                      collapse = " + ")
        rhs <- paste(c(trm, harm), collapse = " + ")
        f_null <- stats::as.formula(paste("y ~", rhs, "+ grp"))
        f_full <- stats::as.formula(paste("y ~", rhs, "+ grp + grp:(", harm, ")"))
        list(long = long, f_null = f_null, f_full = f_full,
             m_null = stats::lm(f_null, data = long),
             m_full = stats::lm(f_full, data = long))
      }, error = function(e) list(error = conditionMessage(e)))

      if(!is.null(pmc$error)) {
        cat("  Could not be fitted: ", pmc$error, "\n", sep = "")
      } else {
        an <- stats::anova(pmc$m_null, pmc$m_full)
        cat("Test: do the harmonic coefficients differ by group?\n")
        cat(sprintf("  F(%d, %d) = %s, p = %s\n",
                    an$Df[2], an$Res.Df[2], fmt3(an$F[2]),
                    format.pval(an$`Pr(>F)`[2], digits = 3, eps = 1e-16)))
        r2n <- summary(pmc$m_null)$r.squared; r2f <- summary(pmc$m_full)$r.squared
        cat(sprintf("  Partial R-squared for the group x harmonic block: %s\n",
                    fmt4(max(0, (r2f - r2n) / (1 - r2n)))))
        cat("\n  ! Standard errors here assume independent observations. Each subject\n")
        cat("    contributes ", length(mod$time_vec), " correlated points, so these p-values are\n", sep = "")
        cat("    ANTICONSERVATIVE. Install lme4 and refit with (1|subject) for the\n")
        cat("    interval you would actually report; the point estimates are unbiased\n")
        cat("    either way. This is stated rather than silently ignored.\n")
        if(isTRUE(mod$using_smoothed))
          cat("    On smoothed data the dependence is worse still (see the data-source note).\n")
      }
    }
  })
