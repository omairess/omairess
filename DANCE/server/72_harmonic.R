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
#      equations now go through dance_format_equation(); the duplicate builder is
#      deleted.
#  1.2 The Rayleigh test used the amplitude-weighted resultant (0.824 -> Z=886)
#      where the unweighted one was required (0.789 -> Z=812). dance_resultants()
#      returns both under unambiguous names; only the unweighted one reaches
#      dance_rayleigh(). The amplitude-weighted vector mean is unchanged.
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
  # ==========================================================================
  # STUDY DESIGN (P21 phase 4)
  # ==========================================================================
  # One helper renders every factor selector, so the between picker in the
  # "between" panel and the between picker in the "mixed" panel cannot drift
  # apart. Shiny needs distinct input ids per panel, hence the pairs.
  .harmonic_cat_vars <- reactive({
    req(values$covariates)
    nm <- names(values$covariates)[vapply(values$covariates, function(x)
      is.factor(x) || is.character(x) || length(unique(x[!is.na(x)])) <= 12, logical(1))]
    nm
  })
  .harmonic_factor_picker <- function(id, label, allow_none = TRUE) {
    v <- .harmonic_cat_vars()
    ch <- if (allow_none) c("None" = "_none_", v) else v
    selectInput(id, label, choices = ch)
  }
  output$harmonic_between_var_ui  <- renderUI(.harmonic_factor_picker(
    "harmonic_between_var",  "Group variable (between-subject):"))
  output$harmonic_between_var_ui2 <- renderUI(.harmonic_factor_picker(
    "harmonic_between_var2", "Between-subject factor:"))
  output$harmonic_within_var_ui   <- renderUI(.harmonic_factor_picker(
    "harmonic_within_var",   "Condition variable (within-subject):"))
  output$harmonic_within_var_ui2  <- renderUI(.harmonic_factor_picker(
    "harmonic_within_var2",  "Within-subject factor:"))

  # Which factors the current design actually selects, in one place, so the
  # readout below and the fitter cannot disagree about what was chosen.
  harmonic_design_terms <- reactive({
    d <- input$harmonic_design %||% "between"
    pick <- function(x) if (!is.null(x) && nzchar(x) && x != "_none_") x else NULL
    switch(d,
      between = list(design = d, between = pick(input$harmonic_between_var),  within = NULL),
      within  = list(design = d, between = NULL, within = pick(input$harmonic_within_var)),
      mixed   = list(design = d, between = pick(input$harmonic_between_var2),
                     within  = pick(input$harmonic_within_var2)))
  })

  # THE ONE GROUPING VARIABLE THE REST OF THE MODULE ALREADY UNDERSTOOD.
  # ------------------------------------------------------------------------
  # Eleven outputs -- the fitted-curve plot, the polar plot, all four parameter
  # histograms, the individual table, the export and the legacy comparison --
  # were written against a single `harmonic_group_var` input. Replacing that
  # control with the Study Design panel left every one of them reading an input
  # that no longer existed, so they silently fell back to "no group" and drew a
  # pooled curve with no error anywhere. That is the failure this derivation
  # prevents: the design panel is now the single source, and the legacy
  # consumers read it through here rather than through an input that is gone.
  #
  # A between-participant factor is the grouping variable when there is one;
  # otherwise the within-participant factor, since those displays split on
  # whatever single factor the design offers. A mixed design shows the between
  # factor, which is what a two-stage display can represent -- the crossing
  # itself is tab 6's job.
  harmonic_group_var_eff <- reactive({
    dance_group_var_from_design(input)
  })

  # THE CLASSIFICATION IS SHOWN BEFORE THE FIT, NOT AFTER IT.
  # The app can read between/within off the data by counting levels per
  # participant, and it does -- but a classification the user cannot see is a
  # guess wearing a confident face, and it decides the random-effects structure.
  # This panel states what the data say, whether it matches what was selected,
  # and how many curves that implies.
  output$harmonic_design_readout <- renderUI({
    dt <- harmonic_design_terms()
    chosen <- c(dt$between, dt$within)
    if (!length(chosen)) return(helpText(HTML(
      "<i>No design factor selected \u2014 one trajectory will be fitted for the whole sample.</i>")))
    req(values$covariates, values$data)
    subj <- rownames(values$data) %||% as.character(seq_len(nrow(values$data)))
    rows <- lapply(chosen, function(f) {
      v <- values$covariates[[f]]
      if (is.null(v)) return(sprintf("<li><b>%s</b> \u2014 not found in the data</li>", f))
      per <- tapply(as.character(v), subj, function(x) length(unique(x[!is.na(x)])))
      per <- per[is.finite(per) & per > 0]
      nlv <- length(unique(v[!is.na(v)]))
      role <- if (nlv < 2) "constant" else if (max(per) == 1) "between"
              else if (min(per) > 1) "within" else "partial"
      want <- if (identical(f, dt$between)) "between" else "within"
      ok <- identical(role, want)
      sprintf("<li><b>%s</b>: %d level(s), the data look <b>%s</b> %s</li>",
              f, nlv, role,
              if (ok) "\u2713"
              else sprintf("\u2014 <span style='color:#b8860b'>you selected it as %s</span>", want))
    })
    n_sub <- length(unique(subj))
    n_cells <- prod(vapply(chosen, function(f)
      max(1L, length(unique(values$covariates[[f]][!is.na(values$covariates[[f]])]))), integer(1)))
    helpText(HTML(sprintf(
      "<b>Read from your data:</b><ul style='margin:4px 0 4px 16px;padding:0'>%s</ul>%d participants, %d design cell(s).%s",
      paste(unlist(rows), collapse = ""), n_sub, n_cells,
      if (!is.null(dt$within))
        sprintf(" Each participant contributes up to %d curve(s).",
                max(1L, length(unique(values$covariates[[dt$within]][!is.na(values$covariates[[dt$within]])]))))
      else " One curve per participant.")))
  })

  # Parameter bounds hints based on data
  output$harmonic_bounds_hints <- renderUI({
    req(values$data)

    # The raw observations, which is what the fit will use: a bounds hint read
    # off a smoothed copy describes a different range from the one being fitted
    Y <- values$data

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
    if (harmonic_is_mixed()) {
      # ONE model: the cell curves are the fit; a participant is that cell's
      # curve plus their conditional modes. Chosen by curve id, not row index,
      # because a within-participant design gives a participant several curves.
      ff <- harmonic_traj(); if (!isTRUE(ff$ok)) return(NULL)
      d <- ff$spec$data
      cv <- if ("curve" %in% names(d)) "curve" else "subject"
      curves <- unique(as.character(d[[cv]]))
      return(selectInput("harmonic_subject_select", "Show:",
                  choices = c("Cell curves (the fitted model)" = "mean",
                              "All participants (shrunken) over the cells" = "all",
                              stats::setNames(curves, paste("Participant", curves))),
                  selected = "mean"))
    }
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

    # ========================================================================
    # AUDIT (P21, finding A4)
    #
    # This used to read
    #
    #     if (use_bounds && trend_type %in% c("linear", "log", "none"))
    #       return(fit_cosinor_nonlinear(...))
    #
    # -- i.e. ticking a bound sent a model that is LINEAR IN EVERY PARAMETER to
    # nlsLM, inheriting its convergence failures, iteration-limit returns and
    # boundary pinning for no gain at all. §3 of the redesign brief names this
    # case directly: do not run a fixed-basis model through a nonlinear
    # optimiser.
    #
    # The bounds are inequality CONSTRAINTS, and the constrained and
    # unconstrained least-squares solutions coincide exactly whenever no
    # constraint is active. So: fit the closed form first, and fall back to the
    # constrained optimiser only when the unconstrained solution actually
    # violates a bound. In the common case the user gets the exact lm fit --
    # with exact standard errors -- instead of an optimiser's approximation of
    # it; in the rare case nothing is lost, because that is the branch that was
    # running before.
    #
    # `bounds_active` records which it was, so the readout can distinguish "no
    # bound was binding" from "a bound was applied".
    # ========================================================================
    .bounds_violated <- function(fit) {
      if (!isTRUE(fit$success)) return(TRUE)
      v <- character(0)
      chk <- function(val, lo, hi, nm) {
        if (is.null(val) || !is.finite(val)) return(invisible(NULL))
        if (is.finite(lo) && val < lo) v <<- c(v, nm)
        if (is.finite(hi) && val > hi) v <<- c(v, nm)
      }
      chk(fit$mesor, mesor_min, mesor_max, "mesor")
      for (a in fit$amplitudes) chk(a, amplitude_min, amplitude_max, "amplitude")
      length(v) > 0
    }

    if(use_bounds && trend_type %in% c("linear", "log", "none")) {
      .free <- fit_cosinor(time, y, period, n_harmonics, trend_type,
                           use_bounds = FALSE)
      if(!.bounds_violated(.free)) {
        .free$bounds_active <- FALSE
        .free$bounds_requested <- TRUE
        .free$fit_route <- "closed-form least squares (no bound was binding)"
        return(.free)
      }
      .con <- fit_cosinor_nonlinear(time, y, period, n_harmonics, trend_type,
                                    FALSE, NULL, 0.32, 0.66,
                                    use_bounds, mesor_min, mesor_max, amplitude_min, amplitude_max,
                                    A_sat_min, A_sat_max, tau_min, tau_max)
      .con$bounds_active <- TRUE
      .con$bounds_requested <- TRUE
      .con$fit_route <- "constrained optimiser (the unconstrained fit violated a bound)"
      return(.con)
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
    # AUDIT (P12.2): the joint elliptical confidence region was computed only by
    # the NONLINEAR fitter, so the default analysis -- a linear cosinor with no
    # trend, which is what most users run -- reported an amplitude and an
    # acrophase with no confidence limits at all. Cornelissen (2014) and Bingham
    # et al. (1982) both prescribe limits read off the error ellipse of the
    # (cosine, sine) pair rather than intervals that treat amplitude and phase as
    # independent, and the ellipse is also what tells you when a phase is not
    # identified (the region contains the origin). dance_bingham_ci() is already a
    # shared helper and this fitter already has the covariance matrix; there was
    # no reason for the two fitters to differ, and a reader cannot tell from the
    # output which fitter produced their numbers.
    bingham <- vector("list", n_harmonics)

    for(h in 1:n_harmonics) {
      cos_idx <- coef_offset + 2 * (h - 1) + 1
      sin_idx <- coef_offset + 2 * (h - 1) + 2
      
      beta_cos <- coefs[cos_idx]
      beta_sin <- coefs[sin_idx]
      
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi

      bingham[[h]] <- dance_bingham_ci(as.numeric(beta_cos), as.numeric(beta_sin),
                                     vcov_mat[c(cos_idx, sin_idx),
                                              c(cos_idx, sin_idx), drop = FALSE],
                                     length(y), length(coefs),
                                     level = 0.95, period = period, harmonic = h)

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
    r_squared <- dance_r_squared(ss_resid, ss_total)   # P0.6: SST = 0 -> NA, not NaN

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
    .zt <- dance_zero_amplitude_test(ss_resid, ss_resid_trend_only,
                                   n, n_predictors, n_harmonics)
    f_stat  <- .zt$F
    p_value <- .zt$p

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian linear model
    # P0.6: a perfect fit gave log(0) -> log_lik = Inf -> AIC = -Inf, so a flat
    # subject beat every rival model. sigma^2 is floored relative to the data.
    sigma_sq <- ss_resid / n
    log_lik <- dance_gaussian_loglik(ss_resid, n, y_scale = stats::sd(y))

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
          r_squared_S <- dance_r_squared(ss_resid_trend, ss_total)
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
        r_squared_S <- dance_r_squared(ss_resid_trend, ss_total)
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
    r_squared_C <- dance_r_squared(ss_resid_circ, ss_total)

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
    .cm <- dance_commonality(r_squared, r_squared_S, r_squared_C)
    unique_S <- .cm$unique_S
    unique_C <- .cm$unique_C
    shared_SC <- .cm$shared
    .pc <- dance_commonality_pct(.cm)
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
      bingham = bingham,                # P12.2: joint elliptical CI per harmonic
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
  dance_conv_of <- function(fit) {
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
        .cv <- dance_conv_of(nls_fit)
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
        .cv <- dance_conv_of(nls_fit)
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
        .cv <- dance_conv_of(nls_fit)
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
        .cv <- dance_conv_of(nls_fit)
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

        .cv <- dance_conv_of(nls_fit)
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

  # ---- AUDIT (P0.5): make the amplitude bound mean what the UI says ----------
  # The UI offers "Amplitude bounds" and the code bounded b_cos and b_sin
  # SEPARATELY at +/- amplitude_max, with the comment
  #     # Since amplitude = sqrt(b_cos^2 + b_sin^2), bound each coefficient
  # A box of half-width A permits sqrt(A^2 + A^2) = A*sqrt(2), so a user who
  # asked for a maximum of 10 could be handed an amplitude of 14.14. The
  # constraint region is a DISC and a box cannot express it.
  #
  # Refitting inside the inscribed square (half-width A/sqrt(2)) does satisfy
  # the constraint, but it also forbids legitimate amplitudes up to A at
  # off-axis phases, so applying it unconditionally would distort every fit to
  # protect the few that violate. It is therefore applied ONLY when the
  # converged fit actually breaks the promise: fits inside the disc keep the
  # full box, violating fits are refitted inside the inscribed square, which is
  # the largest box guaranteed to satisfy the stated bound.
  #
  # amplitude_min is NOT enforced here. It is a floor on a fitted magnitude,
  # which is not a constraint an optimiser should be given -- forcing a flat
  # subject up to a minimum amplitude would invent a rhythm. It is reported
  # instead; see amplitude_below_min in the returned object.
  amp_bound_action <- "not requested"
  if (use_bounds && is.finite(amplitude_max) && n_harmonics >= 1) {
    .amp_of <- function(fit) {
      cf <- coef(fit)
      vapply(seq_len(n_harmonics), function(h) {
        a <- cf[paste0("b_cos", h)]; b <- cf[paste0("b_sin", h)]
        if (is.na(a) || is.na(b)) NA_real_ else sqrt(a^2 + b^2)
      }, numeric(1))
    }
    amps0 <- .amp_of(nls_fit)
    if (any(is.finite(amps0) & amps0 > amplitude_max * (1 + 1e-9))) {
      lb2 <- lower_bounds; ub2 <- upper_bounds
      inscribed <- amplitude_max / sqrt(2)
      for (h in seq_len(n_harmonics)) {
        lb2[paste0("b_cos", h)] <- -inscribed; ub2[paste0("b_cos", h)] <- inscribed
        lb2[paste0("b_sin", h)] <- -inscribed; ub2[paste0("b_sin", h)] <- inscribed
      }
      st2 <- as.list(coef(nls_fit))
      for (nm in names(st2)) {
        if (nm %in% names(lb2) && is.finite(lb2[[nm]]) && st2[[nm]] < lb2[[nm]])
          st2[[nm]] <- lb2[[nm]] * 0.9
        if (nm %in% names(ub2) && is.finite(ub2[[nm]]) && st2[[nm]] > ub2[[nm]])
          st2[[nm]] <- ub2[[nm]] * 0.9
      }
      refit <- tryCatch(minpack.lm::nlsLM(as.formula(formula_str), data = fit_data,
                 start = st2, lower = lb2, upper = ub2,
                 control = minpack.lm::nls.lm.control(maxiter = 300)),
                 error = function(e) NULL)
      if (!is.null(refit) && isTRUE(dance_conv_of(refit)$ok)) {
        nls_fit <- refit
        amp_bound_action <- sprintf(
          "amplitude exceeded the stated maximum (%.4g); refitted inside the inscribed box (+/-%.4g per coefficient)",
          max(amps0, na.rm = TRUE), inscribed)
      } else {
        amp_bound_action <- sprintf(
          "amplitude %.4g exceeds the stated maximum %.4g and the constrained refit did not converge -- this fit does NOT satisfy the bound you set",
          max(amps0, na.rm = TRUE), amplitude_max)
      }
    } else {
      amp_bound_action <- "satisfied by the unconstrained-in-disc fit"
    }
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
      amp_se[h]  <- dance_amp_se(beta_cos, beta_sin, Vh)
      acro_se[h] <- dance_acro_se(beta_cos, beta_sin, Vh)

      # Bingham elliptical joint region for (amplitude, acrophase)
      bingham[[h]] <- dance_bingham_ci(beta_cos, beta_sin, Vh, n, length(coefs),
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
    r_squared <- dance_r_squared(ss_resid, ss_total)   # P0.6: SST = 0 -> NA, not NaN
    percent_rhythm <- max(0, r_squared * 100)

    # Calculate p-value for circadian rhythm
    n_params <- length(coefs)

    # AUDIT (extra a): test the harmonics GIVEN the trend, by refitting the
    # trend-only model. See the note in fit_cosinor() -- the old numerator was
    # the whole model's SS over the harmonics' df alone.
    # AUDIT (P0.3): the exp_sat branch used to freeze tau at its full-model value
    # and refit only M and A_sat by OLS. That is not the nested null -- tau is
    # free under it -- so SSE0 was too large, F was biased up, and the test
    # rejected true nulls at 14.6% against a nominal 5%. It now refits tau.
    # See dance_reduced_exp_sat_sse() in server/08_helpers_cosinor.R.
    ss_resid_trend_only <- tryCatch({
      if(trend_type == "exp_sat" && !is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
        dance_reduced_exp_sat_sse(
          y, time, t_offset,
          start = list(mesor = as.numeric(coefs[["mesor"]]),
                       A_sat = trend_params$A_sat$coef,
                       tau   = trend_params$tau$coef),
          lower = as.list(lower_bounds), upper = as.list(upper_bounds))
      } else if(trend_type == "linear") {
        sum(residuals(lm(y ~ time))^2)
      } else if(trend_type == "log") {
        sum(residuals(lm(y ~ log(time - t_offset + 1)))^2)
      } else {
        ss_total
      }
    }, error = function(e) ss_total)

    .zt <- dance_zero_amplitude_test(ss_resid, ss_resid_trend_only,
                                   n, n_params, n_harmonics)
    f_stat  <- .zt$F
    p_value <- .zt$p
    df2 <- .zt$df2

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian nonlinear model
    # P0.6: a perfect fit gave log(0) -> log_lik = Inf -> AIC = -Inf, so a flat
    # subject beat every rival model. sigma^2 is floored relative to the data.
    sigma_sq <- ss_resid / n
    log_lik <- dance_gaussian_loglik(ss_resid, n, y_scale = stats::sd(y))

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
    .cm <- dance_commonality(r_squared, r_squared_S, r_squared_C)
    unique_S <- .cm$unique_S
    unique_C <- .cm$unique_C
    shared_SC <- .cm$shared
    .pc <- dance_commonality_pct(.cm)
    percent_S <- .pc$unique_S
    percent_C <- .pc$unique_C
    percent_shared <- .pc$shared

    # which bounds, if any, this fit ended up sitting on
    .bh <- dance_bounds_hit(coefs, lower_bounds, upper_bounds)

    list(
      success = TRUE,
      # AUDIT (P0.5): what the amplitude bound actually did to this fit, and
      # whether any harmonic came out below the minimum the user asked for.
      # amplitude_min is reported, never imposed: forcing a flat subject up to
      # a floor would invent a rhythm.
      amp_bound_action = amp_bound_action,
      amplitude_below_min = if (use_bounds && is.finite(amplitude_min) &&
                                amplitude_min > 0 && length(amplitudes))
        which(amplitudes < amplitude_min) else integer(0),
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
                          rep(0, length(newtime))
      )
      pred <- pred + trend_val
    }
    
    # Determine coefficient offset based on trend type
    n_trend_coefs <- switch(trend_type,
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, 0)
    coef_offset <- 1 + n_trend_coefs

    if(component == "total" || component == "all") {
      for(h in 1:n_harmonics) {
        omega <- 2 * pi * h / period
        if(trend_type == "exp_sat") {
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
      if(trend_type == "exp_sat") {
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
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, 0)
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
                          rep(0, length(newtime))
      )
      components$trend <- trend_val
    }

    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      if(trend_type == "exp_sat") {
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
                    trend_type
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
  
  # AUDIT (P15.2). The help text under the nested-model checkbox was written by
  # hand and said "trend in {none, linear, saturating} x harmonics in {1,2,3} ...
  # 9 fits per subject". None of that survived P14: the set follows the harmonic
  # count, and the logarithmic trend was added. Rendering it from the same
  # function that builds the set is the only version that cannot go stale.
  output$harmonic_model_selection_help <- renderUI({
    k  <- suppressWarnings(as.integer(input$n_harmonics %||% 1))
    if (!is.finite(k) || k < 1) k <- 1L
    ms <- dance_cosinor_model_set(k)
    helpText(HTML(sprintf(
      paste0("Fits every trend the app offers (none, linear, logarithmic, ",
             "saturating exponential) crossed with the <b>cumulative</b> harmonic ",
             "sets up to the %d you have selected &mdash; %s &mdash; and reports ",
             "&Delta;AICc with Akaike weights. Harmonics are cumulative, so H1&ndash;H2 ",
             "is the model containing harmonics 1 <i>and</i> 2; a higher harmonic is ",
             "never fitted without the ones below it, and nothing beyond your ",
             "selection is fitted at all. <b>Slow:</b> %d fits per subject."),
      k,
      paste(vapply(ms, function(m) m$label, character(1)), collapse = ", "),
      length(ms))))
  })

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
      # THE CHOICE IS GONE, AND IT ALWAYS RUNS ON RAW OBSERVATIONS.
      #
      # Cosinor is a regression on the observations: it handles missing and
      # unequally spaced data natively, so smoothing first buys nothing here and
      # costs a great deal. Interpolating removes independent noise and induces
      # residual autocorrelation, so R-squared is inflated, LOOCV is optimistic
      # (a held-out point is partly rebuilt from its neighbours) and the
      # zero-amplitude F test is anticonservative. The option existed with a
      # warning attached, which makes a reader responsible for not choosing the
      # wrong one; removing it is the honest version of that warning.
      # ======================================================================
      Y <- values$data
      # WHICH ESTIMATOR. Everything downstream reads mod$approach and shows the
      # output of that estimator only. The per-participant fits below run under
      # both approaches: under mixed-effects they are not displayed, but they
      # are the same cheap OLS pass and keep the model object one shape for the
      # export and the report.
      approach <- input$harmonic_approach %||% "mixed"
      n_subjects <- nrow(Y)
      n_time <- ncol(Y)
      period <- input$harmonic_period
      n_harmonics <- input$n_harmonics
      trend_type <- input$harmonic_trend_type
      
      # Calculate number of trend parameters
      n_trend_params <- switch(trend_type,
                               "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
      
      # Diagnostic: Check data type and dimensions
      cat(sprintf("Data diagnostics: %d subjects × %d time points, type=%s, trend=%s\n", 
                  n_subjects, n_time, typeof(Y), trend_type))
      
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
      } else if(subjects_with_nas > 0) {
        # NOT a reason to smooth first: the cosinor uses the observations that
        # are present, which is what makes interpolation unnecessary here.
        showNotification(
          sprintf("%d subjects have missing values. The cosinor uses the observations present for each subject.", 
                  subjects_with_nas),
          type = "message", duration = 8)
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
      # ALWAYS the first observation. The midnight origin was kept "for
      # continuity with earlier runs"; it put the trend and the harmonics on two
      # different anchors and made the intercept the value at neither. There is
      # no analysis that needs it, so it is not offered.
      time_origin <- "first_observation"
      origin_shift <- min(time_vec, na.rm = TRUE)
      time_vec_model <- time_vec - origin_shift

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
            row_data$mesor_adj <- dance_rhythm_adjusted_mean(
              fit_i$mesor, trend_type, .tc,
              min(time_vec_model, na.rm = TRUE), max(time_vec_model, na.rm = TRUE),
              fit_i$t_offset %||% 0)
            row_data$value_at_start <- tryCatch(
              as.numeric(dance_rhythm_from_coefs(
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
      bounds_summary <- dance_bounds_summary(
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
      fit_audit$bounds_kept <- dance_bounds_summary(
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
        # clock origin is added only at display, by dance_acrophase_label().
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
        # dance_resultants() now returns both, under names that cannot be
        # confused, and only the unweighted one reaches dance_rayleigh(). The
        # amplitude-weighted vector mean stays the population estimator: that
        # part was always correct.
        # ====================================================================
        res1 <- dance_resultants(individual_params$acrophase_rad_1,
                               individual_params$amplitude_1)
        n_valid <- res1$n
        r_bar_unweighted <- res1$r_unweighted
        r_bar_weighted   <- res1$r_weighted
        ray <- dance_rayleigh(r_bar_unweighted, n_valid)
        rayleigh_z <- ray$Z
        rayleigh_p <- ray$p

        # circular dispersion is defined on the unweighted resultant
        r_bar <- r_bar_unweighted        # kept for downstream compatibility
        circ_var <- 1 - r_bar_unweighted
        circ_sd <- res1$circ_sd_rad

        # per-harmonic resultants, for the report
        resultants <- lapply(seq_len(n_harmonics), function(h) {
          rr <- dance_resultants(individual_params[[paste0("acrophase_rad_", h)]],
                               individual_params[[paste0("amplitude_", h)]])
          if(is.null(rr)) return(NULL)
          rr$rayleigh <- dance_rayleigh(rr$r_unweighted, rr$n)
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
        # dance_format_equation(); this fills the gap the pooled one was reading
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
        rhythm_adjusted_mean <- dance_rhythm_adjusted_mean(
          mean_mesor, trend_type, pop_trend_coefs, t_lo, t_hi, t_offset_global)
        # how much of the oscillation leaks into the window mean, which is only
        # exactly zero over a whole number of periods (here 22 h of a 24 h cycle)
        harmonic_leak <- dance_harmonic_window_mean(
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
      if(!is.null(harmonic_group_var_eff()) && harmonic_group_var_eff() != "_none_") {
        group_var <- values$covariates[[harmonic_group_var_eff()]]

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
        group_audit <- dance_group_audit(group_var, individual_params$subject, min_n = 3)
        if(group_audit$n_unassigned > 0) {
          showNotification(
            sprintf("%d of %d fitted subject%s no usable '%s' label. They are pooled but not grouped, and appear as UNASSIGNED in the report.",
                    group_audit$n_unassigned, group_audit$n_total,
                    if(group_audit$n_unassigned == 1) " has" else "s have",
                    harmonic_group_var_eff()),
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
            grp_ram <- dance_rhythm_adjusted_mean(
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
              rr <- dance_resultants(grp_params[[paste0("acrophase_rad_", h)]],
                                   grp_params[[paste0("amplitude_", h)]])
              if(is.null(rr)) return(NULL)
              rr$rayleigh <- dance_rayleigh(rr$r_unweighted, rr$n)
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
                    group_audit$n_unassigned, harmonic_group_var_eff(),
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
      if(isTRUE(input$harmonic_bootstrap) && identical(approach, "two_stage")) {
        B <- input$harmonic_n_boot
        boot_mesor <- numeric(B)
        boot_amplitude <- numeric(B)
        boot_acrophase <- numeric(B)
        
        showNotification(paste("Running", B, "bootstrap iterations..."), type = "message")
        
        # ====================================================================
        # AUDIT (P21, findings A2 and A3)
        #
        # A2. This loop used to draw indices with replacement and then SELECT
        # rows with `individual_params$subject %in% boot_idx`. `%in%` is a set
        # test, so a participant drawn three times contributed one row: the
        # procedure was a subsample WITHOUT replacement of the ~63% of
        # participants drawn at least once, and the finite-population
        # correction that comes with it made the interval 21% too SHORT on a
        # 60-participant fixture. Resampling is now by INDEX, so a participant
        # drawn three times appears three times, which is the whole point.
        #
        # The resampling unit is the PARTICIPANT, which is what §11 of the
        # redesign brief requires and what the repeated-measures structure
        # demands: a participant's whole row travels together.
        #
        # A3. The acrophase interval was `quantile()` on hours -- a linear
        # quantile on a circular quantity. It is now the circular percentile
        # interval (quantiles of the signed deviation from the bootstrap mean
        # direction), computed by dance_boot_circ_ci_time(). The linear one is
        # kept alongside ONLY so the readout can show what it would have said
        # when the two disagree, which they do violently near midnight.
        # ====================================================================
        # NOTE ON CLUSTER IDENTITY. Every statistic in this loop is a MEAN OVER
        # ROWS, so a participant drawn three times correctly contributes three
        # times and their id never enters the arithmetic. That is why plain row
        # indices are enough here. Any bootstrap that REFITS a model with the
        # participant as a grouping factor must instead use
        # dance_boot_clusters(), which gives each drawn copy a new id -- three
        # copies sharing one id would be read as one participant with tripled
        # observations, i.e. one random effect where the bootstrap intended
        # three.
        rows_by_subject <- split(seq_len(nrow(individual_params)),
                                 individual_params$subject)
        subj_ids <- names(rows_by_subject)
        n_boot_subj <- length(subj_ids)

        withProgress(message = 'Bootstrap...', value = 0, {
          for(b in 1:B) {
            take <- dance_boot_index(n_boot_subj)
            boot_rows <- unlist(rows_by_subject[subj_ids[take]], use.names = FALSE)
            boot_params <- individual_params[boot_rows, , drop = FALSE]

            boot_mesor[b] <- mean(boot_params$mesor, na.rm = TRUE)

            # Use first harmonic for bootstrap CIs
            x_b <- boot_params$amplitude_1 * cos(boot_params$acrophase_rad_1)
            y_b <- boot_params$amplitude_1 * sin(boot_params$acrophase_rad_1)
            boot_amplitude[b] <- sqrt(mean(x_b, na.rm = TRUE)^2 + mean(y_b, na.rm = TRUE)^2)

            acro_b <- atan2(mean(y_b, na.rm = TRUE), mean(x_b, na.rm = TRUE))
            if(is.finite(acro_b) && acro_b < 0) acro_b <- acro_b + 2 * pi
            boot_acrophase[b] <- acro_b

            if(b %% 50 == 0) incProgress(50 / B)
          }
        })

        .acro_ci <- dance_boot_circ_ci_time(boot_acrophase, period, 1)
        boot_results <- list(
          mesor_ci = quantile(boot_mesor, c(0.025, 0.975), na.rm = TRUE),
          amplitude_ci = quantile(boot_amplitude, c(0.025, 0.975), na.rm = TRUE),
          # circular, and carrying its own width because lo > hi when it wraps
          acrophase_ci = c(`2.5%` = .acro_ci$lo, `97.5%` = .acro_ci$hi),
          acrophase_ci_circular = .acro_ci,
          # what a linear quantile would have reported, for the readout to
          # contrast when the interval wraps
          acrophase_ci_linear = quantile(phi_to_hours(boot_acrophase, period, 1),
                                         c(0.025, 0.975), na.rm = TRUE),
          resample_unit = "participant",
          n_resampled = n_boot_subj,
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
      if(isTRUE(input$harmonic_model_selection) && identical(approach, "two_stage")) {
        # AUDIT (P14). The grid used to be hardcoded here as
        # c("none", "linear", "exp_sat") x 1:3. It omitted the `log` trend the
        # UI offers -- so a user who selected it saw a table their own model was
        # not in, with weights normalised over a set that excluded it -- and it
        # ran to three harmonics whatever had been selected, spending compute on
        # models the user had already ruled out and diluting the Akaike weight
        # of the one they were reporting. The set now comes from
        # dance_cosinor_model_set(), which crosses every trend the UI offers with
        # the CUMULATIVE harmonic sets up to the number selected.
        ms_models <- dance_cosinor_model_set(n_harmonics)
        ms_rows <- list()
        withProgress(message = "Fitting the nested model set...", value = 0, {
          for(m in ms_models) {
            incProgress(1 / length(ms_models),
                        detail = sprintf("%s, %s", dance_trend_label(m$trend),
                                         dance_harmonic_label(m$k)))
            # parameters + 1: a fit needs at least one observation more than it
            # has parameters for the residual variance to be estimable
            need <- dance_model_npar(m$trend, m$k) + 1L
            if(n_time < need) next
            aicc_sum <- 0; k_ok <- 0L
            for(i in seq_len(n_subjects)) {
              y_i <- as.numeric(Y[i, ])
              if(sum(!is.na(y_i)) < need) next
              f <- tryCatch(fit_cosinor(time_vec_model, y_i, period = period,
                                        n_harmonics = m$k, trend_type = m$trend),
                            error = function(e) NULL)
              if(is.null(f) || !isTRUE(f$success) || !is.finite(f$aicc)) next
              aicc_sum <- aicc_sum + f$aicc; k_ok <- k_ok + 1L
            }
            if(k_ok > 0) ms_rows[[m$label]] <- aicc_sum / k_ok
          }
        })
        model_selection <- dance_akaike_table(ms_rows)
        if(!is.null(model_selection)) {
          attr(model_selection, "per_subject_mean") <- TRUE
          # Which row is the specification the user actually ran, so the readout
          # and the report can point at it instead of leaving them to match a
          # label by eye.
          attr(model_selection, "selected") <- dance_model_label(trend_type, n_harmonics)
          attr(model_selection, "n_harmonics_max") <- n_harmonics
        }
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
      if(identical(approach, "two_stage") && trend_type == "exp_sat" && nrow(individual_params) > 0) {
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
        #
        # AUDIT (P15.3). A cleared numericInput yields NA, not NULL, so `%||%`
        # never fired and `is.finite(NA)` was FALSE: the whole free-vs-fixed
        # check was skipped, and because the readout only prints the line when
        # the result exists, it vanished with nothing said. "I left the box
        # empty" and "this comparison was run and found nothing" then looked
        # identical on screen. The skip is now recorded and reported.
        tau_raw <- input$harmonic_tau_fixed
        tau_fix <- suppressWarnings(as.numeric(if (is.null(tau_raw)) NA else tau_raw))
        if(!is.finite(tau_fix) || tau_fix <= 0) {
          conditioning$tau_fixed_skipped <- TRUE
        }
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
      hm_new <- list(
        approach = approach,
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
        data_source = "raw",
        time_origin = time_origin,
        origin_shift = origin_shift,
        time_vec_clock = time_vec,           # display axis, always clock-linearised
        group_var_name = if(!is.null(harmonic_group_var_eff()) &&
                            harmonic_group_var_eff() != "_none_") harmonic_group_var_eff() else NULL,
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
        dance_settings = list(
          use_bounds = use_bounds, mesor_min = mesor_min, mesor_max = mesor_max,
          amplitude_min = amplitude_min, amplitude_max = amplitude_max,
          A_sat_min = A_sat_min, A_sat_max = A_sat_max,
          tau_min = tau_min, tau_max = tau_max),
        t_offset = t_offset_global,
        t_center = t_center_global,
        subjects_with_nas = subjects_with_nas,
        Y = Y
      )
      if (identical(approach, "mixed")) {
        # ONE fitted model for the whole module. Fitted here, once, when the
        # button is pressed, and carried on the model object -- so the fitted
        # curves, the polar dial, the individual table, the diagnostics and the
        # comparison tab are all views of the same object and cannot disagree.
        hm_new$traj <- fit_harmonic_traj(hm_new)
        if (!isTRUE(hm_new$traj$ok))
          showNotification(paste("Mixed-effects fit:", hm_new$traj$message),
                           type = "error", duration = 15)
      }
      values$harmonic_model <- hm_new

      showNotification(if (identical(approach, "mixed")) "Mixed-effects cosinor complete!"
                       else "Two-stage cosinor complete!", type = "message")
      
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
  # ----------------------------------------------------------------------------
  # One text, two approaches. The panel reads the model object the Run button
  # stored and prints, in this order: the choices, the fit outcomes, the time
  # points, the model equation (symbolic, then with the fitted numbers per
  # group or per design cell), the information criteria of the fitted model,
  # the model comparison, and the diagnostics. What the old report also carried
  # -- Rayleigh tests, arithmetic-versus-vector means, the commonality
  # percentages, the correlation matrix -- lives on the tab that owns it (polar
  # plot, individual table, diagnostics) and is not repeated here.
  #
  # The report body is a plain function so the on-screen panel and the text
  # download share it verbatim.
  # ============================================================================
  .print_harmonic_summary <- function() {
    req(values$harmonic_model)
    mod <- values$harmonic_model
    mixed <- identical(mod$approach %||% "two_stage", "mixed")
    ff <- if (mixed) mod$traj else NULL
    period <- mod$period
    nh <- mod$n_harmonics
    trend_type <- mod$trend_type %||% "none"
    params <- mod$individual_params
    pop <- mod$pop_mean_fit
    dv <- mod$dv_name %||% "the dependent variable"
    dvu <- mod$dv_units
    clock_o <- dance_clock_origin(mod)
    dfm <- input$harmonic_df_method %||% "kr"

    hdr <- function(x) cat("\n--- ", x, " ---\n", sep = "")
    kv  <- function(k, v) cat(sprintf("%-20s %s\n", paste0(k, ":"), v))
    acro_clock <- function(rad, h)
      dance_acrophase_label(hours = (rad %% (2 * pi)) * (period / h) / (2 * pi),
                            period = period, harmonic = h, clock_origin = clock_o)
    t_first <- min(mod$time_vec, na.rm = TRUE)
    clock_first <- dance_clock_label(clock_o + t_first, period, show_day = FALSE)

    cat("=== Cosinor (harmonic regression) results ===\n")

    # ---- 1. the choices -------------------------------------------------------
    hdr("Choices")
    kv("Approach", if (mixed)
         "MIXED-EFFECTS: one model over every observation; participants are random effects"
       else "TWO-STAGE: one cosinor per participant; the estimates are then compared")
    kv("Dependent variable", paste0(dv, if (!is.null(dvu)) paste0(" (", dvu, ")") else ""))
    kv("Admissible range", if (is.finite(mod$dv_min) || is.finite(mod$dv_max))
         sprintf("[%s, %s]", if (is.finite(mod$dv_min)) fmt2(mod$dv_min) else "-Inf",
                 if (is.finite(mod$dv_max)) fmt2(mod$dv_max) else "Inf")
       else "not specified (set one to have the fitted curves checked)")
    kv("Period", sprintf("%s h, %d harmonic%s (%s)", fmtn(period, 0), nh,
                         if (nh == 1) "" else "s", dance_harmonic_label(nh)))
    kv("Trend", switch(trend_type,
         none = "none (rhythm only)",
         linear = "linear, beta*t",
         log = "logarithmic, beta*log(t + 1)",
         exp_sat = paste0("saturating exponential, A_sat*(1 - e^(-t/tau))",
                          if (mixed && isTRUE(is.finite(ff$spec$tau)))
                            sprintf("; tau held at %s h", fmt1(ff$spec$tau)) else ""),
         trend_type))
    if (mixed && isTRUE(ff$ok)) {
      cls <- ff$spec$classification
      des <- if (length(ff$spec$design_terms)) paste(vapply(ff$spec$design_terms, function(f) {
        r <- cls[cls$factor == f, ][1, ]
        sprintf("%s (%s, %d levels)", f, r$role, r$n_levels)
      }, character(1)), collapse = "; ") else "none -- one trajectory for the whole sample"
      kv("Design factors", des)
      kv("Denominator df", switch(dfm, kr = "Kenward-Roger", satterthwaite = "Satterthwaite", dfm))
    } else {
      kv("Grouping factor", mod$group_var_name %||% "none")
    }
    kv("Time origin", sprintf("first observation; t = 0 is %s on the clock",
                              dance_clock_label(clock_o, period, show_day = FALSE)))
    kv("Data", paste0("raw observations",
                      if (isTRUE((mod$subjects_with_nas %||% 0) > 0))
                        sprintf(" (%d participant(s) have missing values; the fit uses the observations present)",
                                as.integer(mod$subjects_with_nas)) else ""))

    # ---- 2. fit outcomes ------------------------------------------------------
    hdr("Fit outcomes")
    if (mixed) {
      if (!isTRUE(ff$ok)) {
        cat("The mixed-effects model did not fit: ", ff$message, "\n", sep = "")
        cat("Nothing below this line is available for this run.\n")
        return(invisible(NULL))
      }
      cat(sprintf("Participants: %d   curves: %d   observations: %d   design cells: %d\n",
                  ff$spec$n_participants, ff$spec$n_curves, ff$spec$n_obs, nrow(ff$spec$cells)))
      cat(dance_traj_fit_report(ff), "\n", sep = "")
    } else {
      fa <- mod$fit_audit
      if (!is.null(fa)) {
        cat(sprintf("Participants attempted:     %d\n", fa$n_attempted))
        cat(sprintf("  converged, interior:      %d\n", fa$n_converged))
        cat(sprintf("  converged, on a bound:    %d  (%s)\n", fa$n_boundary,
                    if (isTRUE(fa$include_boundary)) "included" else "excluded by your choice"))
        cat(sprintf("  did not converge:         %d  (always excluded)\n", fa$n_nonconverged))
        cat(sprintf("  failed outright:          %d\n", fa$n_failed))
        cat(sprintf("Summaries below use %d participant(s).\n", nrow(params)))
        if (isTRUE(fa$gate_relaxed))
          cat("  ! Too few fits pass the gate; ALL returned fits are included. Provisional.\n")
        bs <- if (isTRUE(fa$include_boundary)) fa$bounds_kept else fa$bounds
        if (!is.null(bs) && isTRUE(bs$n_any > 0)) {
          cat(sprintf("Parameter bounds hit by %d of %d summarised fit(s)%s:\n", bs$n_any, bs$n,
                      if (isTRUE(fa$include_boundary)) " (INCLUDED: their SEs are meaningless and they pull the mean toward the bound)"
                      else " (excluded)"))
          if (!is.null(bs$per_bound))
            for (i in seq_len(nrow(bs$per_bound)))
              cat(sprintf("  %-26s %5d  (%s%%)\n", bs$per_bound$bound[i],
                          bs$per_bound$n[i], fmt1(bs$per_bound$pct[i])))
          if (!is.null(bs$multi))
            cat(sprintf("  %d fit(s) sit on more than one bound (a likelihood ridge; see the parameter CSV).\n",
                        nrow(bs$multi)))
        }
      } else {
        cat(sprintf("Participants summarised: %d\n", nrow(params)))
      }
    }

    # ---- 3. time points -------------------------------------------------------
    hdr("Time points")
    cat("Number of time points: ", length(mod$time_vec), "\n", sep = "")
    if (length(mod$time_vec) <= 24) {
      if (isTRUE(mod$wrap_applied) && !is.null(mod$original_times)) {
        cat("Original times (clock):  ", paste(fmtn(mod$original_times, 1), collapse = ", "), "\n", sep = "")
        cat("Adjusted times (linear): ", paste(fmtn(mod$time_vec_clock %||% mod$time_vec, 1), collapse = ", "), "\n", sep = "")
        cat("(Times after midnight were adjusted for chronological order.)\n")
      } else {
        cat("Clock times: ", paste(fmtn(mod$time_vec_clock %||% mod$time_vec, 2), collapse = ", "), "\n", sep = "")
      }
      if (isTRUE(clock_o > 0))
        cat("Model times (t = 0 at the first observation): ",
            paste(fmtn(mod$time_vec, 2), collapse = ", "), "\n", sep = "")
    }
    diffs <- diff(mod$time_vec)
    if (length(unique(round(diffs, 2))) > 1)
      cat("Spacing: UNEQUAL (", paste(unique(fmtn(diffs, 2)), collapse = ", "), " h)\n", sep = "")
    else
      cat("Spacing: equal (", fmt2(diffs[1]), " h)\n", sep = "")

    # ---- 4. the model equation, symbolic --------------------------------------
    hdr("Model equation (symbolic)")
    sym <- "Y(t) = beta_0"
    if (trend_type == "linear")       sym <- paste0(sym, " + beta*t")
    else if (trend_type == "log")     sym <- paste0(sym, " + beta*log(t+1)")
    else if (trend_type == "exp_sat") sym <- paste0(sym, " + A_sat*(1 - e^(-t/tau))")
    for (h in seq_len(nh))
      sym <- paste0(sym, sprintf(" + A%d*cos(2pi*%d*t/%s - phi%d)", h, h, fmtn(period, 0), h))
    if (mixed) {
      cat(sym, " + b_i(t) + e\n", sep = "")
      cat("  beta_0, beta, A_h, phi_h: fixed effects, one set per design cell\n")
      cat("  b_i(t): participant i's deviation from its cell's curve (random effects on the same basis)\n")
      cat("  e: residual error\n")
      cat("  Fixed effects:  ", ff$spec$fixed_formula, "\n", sep = "")
      cat("  Random effects: ", ff$re_formula, "\n", sep = "")
    } else {
      cat(sym, " + e,   fitted separately to each participant's series\n", sep = "")
    }
    cat("  beta_0 = the constant term (NOT the MESOR when a trend is present)\n")
    cat("  A_h, phi_h = amplitude and acrophase of harmonic h; t in hours from the first observation\n")
    if (trend_type == "exp_sat") cat("  A_sat = asymptote, tau = time constant (h)\n")

    # ---- 5. the fitted equations ---------------------------------------------
    show_curve <- function(label, b0, trend_coefs, amps, acros_rad, n_txt, coefs_for_peak = NULL, peak = NULL) {
      cat("\n", label, if (!is.null(n_txt)) paste0(" (", n_txt, ")") else "", ":\n", sep = "")
      cat("  ", dance_format_equation(b0, trend_type, trend_coefs, amps, acros_rad,
                                     period, mod$t_offset %||% 0), "\n", sep = "")
      for (h in seq_len(nh))
        cat(sprintf("  H%d: amplitude %s, acrophase %s%s\n", h, fmt3(amps[h]),
                    acro_clock(acros_rad[h], h),
                    if (h > 1) sprintf(" (%d maxima per day, all shown)", h) else ""))
      if (!is.null(peak))
        cat(sprintf("  Complete fitted curve: maximum %s at %s, minimum %s at %s%s\n",
                    fmt2(peak$max_value), dance_clock_label(peak$max_clock, period, show_day = FALSE),
                    fmt2(peak$min_value), dance_clock_label(peak$min_clock, period, show_day = FALSE),
                    if (isTRUE(peak$edge)) " [maximum on the window edge: not a turning point]" else ""))
    }
    check_range <- function(yy, what) {
      if (!(is.finite(mod$dv_min) || is.finite(mod$dv_max))) return(invisible(NULL))
      bc <- dance_check_bounds(yy, mod$dv_min, mod$dv_max)
      if (!bc$ok)
        cat(sprintf("  ! %s leaves the admissible range: fitted min %s, max %s. A model that predicts\n    impossible values is misspecified, not merely imprecise.\n",
                    what, fmt2(bc$min), fmt2(bc$max)))
      else
        cat(sprintf("  Fitted range over the window [%s, %s]: within the admissible range.\n",
                    fmt2(bc$min), fmt2(bc$max)))
    }

    if (mixed) {
      hdr("Fitted equations (fixed effects, per design cell)")
      cat("These are the curves drawn on the fitted-curves tab and on the comparison tab:\n")
      cat("the same fixed effects of the same model, with no participant deviation.\n")
      ce <- dance_traj_cell_equations(ff)
      pk <- tryCatch(dance_traj_curve_peaks(ff, "full"), error = function(e) NULL)
      pr_full <- tryCatch(dance_traj_predict(ff, component = "full", n_time = 200), error = function(e) NULL)
      for (i in seq_len(nrow(ce$table))) {
        r <- ce$table[i, ]; cname <- r$cell
        amps <- vapply(seq_len(nh), function(h) r[[paste0("amplitude_", h)]], numeric(1))
        acr  <- vapply(seq_len(nh), function(h) r[[paste0("acrophase_rad_", h)]], numeric(1))
        n_txt <- if (is.finite(r$n_participants)) sprintf("%d participant(s), %d observations",
                                                         r$n_participants, r$n_obs) else NULL
        peak <- if (!is.null(pk) && cname %in% pk$cell) {
          q <- pk[pk$cell == cname, ][1, ]
          list(max_value = q$peak_fit, max_clock = q$peak_clock,
               min_value = q$trough_fit, min_clock = q$trough_clock, edge = !isTRUE(q$peak_interior))
        } else NULL
        show_curve(if (identical(cname, "(all)")) "Whole sample" else paste0("Cell ", cname),
                   r$intercept, ce$trend_coefs[[i]], amps, acr, n_txt, peak = peak)
        cat(sprintf("  Level at the first observation (%s): %s%s\n", clock_first,
                    fmt2(r$level_at_t0), if (!is.null(dvu)) paste0(" ", dvu) else ""))
        if (!is.null(pr_full) && isTRUE(pr_full$ok))
          check_range(pr_full$table$fit[pr_full$table$cell == cname], "This cell's curve")
      }
      cat("\nAmplitude and acrophase are read off the cell's (cos, sin) fixed effects; their\n")
      cat("intervals and the between-cell tests are on the comparison tab. Participant-level\n")
      cat("(shrunken) estimates are in the individual table.\n")
    } else {
      if (is.null(pop)) return(invisible(NULL))
      hdr("Fitted equations (per-participant estimates, averaged)")
      cat("The constant and the trend are arithmetic means of the per-participant fits;\n")
      cat("amplitude and acrophase are amplitude-weighted VECTOR means of the per-participant\n")
      cat("(amplitude, acrophase) pairs, the population estimator for a circular quantity.\n")
      peak_of <- function(coefs) {
        cpk <- dance_curve_peak_clock(coefs, mod)
        if (is.null(cpk)) return(NULL)
        list(max_value = cpk$peak_value, max_clock = cpk$peak_clock,
             min_value = cpk$trough_value, min_clock = cpk$trough_clock, edge = isTRUE(cpk$peak_at_edge))
      }
      show_curve("Pooled", pop$intercept, pop$trend_coefs, pop$mean_amplitudes,
                 pop$mean_acrophases_rad, sprintf("%d participants", nrow(params)),
                 peak = peak_of(pop$mean_coefs))
      .v0 <- dance_value_at(pop$mean_coefs, mod, t_first)
      if (is.finite(.v0))
        cat(sprintf("  Level at the first observation (%s): %s%s\n", clock_first, fmt2(.v0),
                    if (!is.null(dvu)) paste0(" ", dvu) else ""))
      if (trend_type != "none" && is.finite(pop$rhythm_adjusted_mean))
        cat(sprintf("  MESOR (rhythm-adjusted mean over the window): %s; the constant term is %s\n",
                    fmt3(pop$rhythm_adjusted_mean), fmt3(pop$intercept)))
      tt <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 400)
      check_range(predict_from_coefs(pop$mean_coefs, tt, period, nh, trend_type,
                                     pop$t_offset %||% 0, 0), "The pooled curve")

      if (!is.null(mod$group_fits) && length(mod$group_fits) >= 1) {
        ga <- attr(mod$group_fits, "audit")
        n_in <- attr(mod$group_fits, "n_in_groups") %||% NA
        n_fit <- attr(mod$group_fits, "n_fitted") %||% nrow(params)
        cat(sprintf("\nBy %s: group sizes sum to %s of %s fitted participants%s.\n",
                    mod$group_var_name %||% "group", fmtn(n_in, 0), fmtn(n_fit, 0),
                    if (!is.null(ga) && ga$n_unassigned > 0)
                      sprintf("; %d without a usable label appear as UNASSIGNED", ga$n_unassigned) else ""))
        if (!is.null(ga) && length(ga$dropped_small) > 0)
          cat(sprintf("  Not summarised (fewer than 3 fitted participants): %s\n",
                      paste(ga$dropped_small, collapse = ", ")))
        for (g_name in names(mod$group_fits)) {
          g <- mod$group_fits[[g_name]]
          show_curve(paste0("Group '", g_name, "'"), g$intercept, g$trend_coefs,
                     g$mean_amplitudes, g$mean_acrophases_rad, sprintf("n = %d", g$n),
                     peak = peak_of(g$mean_coefs))
          .gv0 <- dance_value_at(g$mean_coefs, mod, t_first)
          if (is.finite(.gv0))
            cat(sprintf("  Level at the first observation (%s): %s%s\n", clock_first, fmt2(.gv0),
                        if (!is.null(dvu)) paste0(" ", dvu) else ""))
          if (trend_type != "none" && is.finite(g$rhythm_adjusted_mean))
            cat(sprintf("  MESOR (rhythm-adjusted mean): %s\n", fmt3(g$rhythm_adjusted_mean)))
        }
        cat("\nGroup differences are tested on the comparison tab, on the per-participant\n")
        cat("estimates (MANOVA, population-mean cosinor, Watson-Williams).\n")
      }
    }

    # ---- 6. information criteria of the fitted model ---------------------------
    hdr("Information criteria of the fitted model")
    if (mixed) {
      m <- ff$model
      ll <- tryCatch(stats::logLik(m), error = function(e) NULL)
      if (!is.null(ll)) {
        cat(sprintf("  logLik %s (%s, %d parameters)   AIC %s   BIC %s\n",
                    fmt2(as.numeric(ll)), if (isTRUE(ff$REML)) "REML" else "ML",
                    as.integer(attr(ll, "df")), fmt2(stats::AIC(m)), fmt2(stats::BIC(m))))
        if (isTRUE(ff$REML))
          cat("  REML criteria compare models with the SAME fixed effects only (different random\n",
              "  structures); to compare fixed effects, the comparison tab refits by ML.\n", sep = "")
      } else cat("  not available for this fit\n")
    } else if (any(c("aic", "aicc", "bic") %in% names(params))) {
      k_par <- dance_model_npar(trend_type, nh) + 1L
      cat(sprintf("  Specification %s, fitted per participant on n = %s observations, k = %s parameters\n",
                  dance_model_label(trend_type, nh), fmtn(length(mod$time_vec), 0), fmtn(k_par, 0)))
      for (nmi in c("aic", "aicc", "bic")) {
        if (!nmi %in% names(params)) next
        v <- suppressWarnings(as.numeric(params[[nmi]])); v <- v[is.finite(v)]
        if (!length(v)) next
        cat(sprintf("    mean %-4s = %s  (SD %s over %s participants)\n",
                    toupper(nmi), fmt2(mean(v)), fmt2(stats::sd(v)), fmtn(length(v), 0)))
      }
      cat("  Means over participants of a per-participant criterion, not the criterion of one\n")
      cat("  pooled fit. Compare them only against the same quantity from another specification,\n")
      cat("  which is what the Delta-AICc table does.\n")
    } else cat("  not available\n")

    # ---- 7. model comparison -----------------------------------------------------
    hdr("Model comparison")
    if (mixed) {
      cat("Random-effects ladder (the fullest structure that converged is the one used):\n")
      if (length(ff$attempts)) {
        for (at in ff$attempts)
          cat(sprintf("  %d. %-52s %s\n", at$rung, at$formula,
                      if (!isTRUE(at$fitted)) "did not fit"
                      else if (!isTRUE(at$converged)) "did not converge"
                      else if (isTRUE(at$singular)) "singular"
                      else "accepted"))
      } else cat(sprintf("  rung %d of %d: %s\n", ff$re_rung, ff$n_rungs, ff$re_label))
      if (isTRUE(ff$singular))
        cat("  A singular fit is a converged fit at a boundary (a variance at zero); the fixed-effect\n",
            "  tests stand and the term is kept.\n", sep = "")
      cal <- dance_traj_calibration(ff, dfm, "full")
      cat("\nCalibration of the omnibus tests: ", if (isTRUE(cal$validated)) "provisional" else "NOT validated", "\n", sep = "")
      cat("  ", cal$calibration, "\n", sep = "")
      cat("The design-factor tests (whole trajectory, shape, circadian, trend, level) and the\n")
      cat("pairwise comparisons are on the comparison tab; nested fixed-effect comparisons there\n")
      cat("use ML refits of this model.\n")
    } else {
      if (!is.null(mod$model_selection)) {
        ms <- mod$model_selection; sel <- attr(ms, "selected")
        cat("Delta-AICc across the nested set (AICc averaged per participant, same participants in every cell):\n")
        cat(sprintf("  %-20s %12s %10s %8s\n", "model", "AICc", "dAICc", "weight"))
        for (i in seq_len(nrow(ms)))
          cat(sprintf("  %-20s %12s %10s %8s%s\n", ms$model[i], fmt2(ms$AICc[i]),
                      fmt2(ms$dAICc[i]), fmt3(ms$weight[i]),
                      if (!is.null(sel) && identical(ms$model[i], sel)) " <-- reported" else ""))
        cat("  Harmonics are cumulative (H1-H2 contains harmonics 1 AND 2). Akaike weights\n")
        cat("  redistribute over this set only: 'best of these', not 'probably correct'.\n")
      } else {
        cat("Nested-model comparison not run (enable 'Compare nested models' under Advanced).\n")
      }
      cn <- mod$conditioning
      if (!is.null(cn) && !is.null(cn$tau_fixed_delta_aic)) {
        cat(sprintf("Free tau vs tau held at %s h: Delta-AIC = %s (%s)\n",
                    fmt1(cn$tau_fixed_value), fmt2(cn$tau_fixed_delta_aic),
                    if (cn$tau_fixed_delta_aic > 0)
                      "the free-tau fit is NOT better by AIC: tau is not identified by these data"
                    else "the free-tau fit is preferred"))
        cat(sprintf("  AIC(free) - AIC(fixed), averaged over the %s participant(s) scorable in both fits.\n",
                    fmtn(cn$tau_fixed_n %||% NA, 0)))
      } else if (!is.null(cn) && isTRUE(cn$tau_fixed_skipped)) {
        cat("Free-tau vs fixed-tau check: NOT RUN -- enter a value in 'tau held at (h)'.\n")
      }
    }

    # ---- 8. diagnostics ----------------------------------------------------------
    hdr("Diagnostics")
    if (mixed) {
      m <- ff$model
      cat(sprintf("  Residual SD (sigma): %s\n", fmt4(stats::sigma(m))))
      r2 <- dance_traj_r2(ff)
      if (!is.null(r2))
        cat(sprintf("  R-squared: marginal %s (fixed effects), conditional %s (fixed + random)\n    [%s]\n",
                    fmt3(r2$marginal), fmt3(r2$conditional), r2$method))
      cat(sprintf("  Convergence %s; singular %s; fixed-effect rank %s.\n",
                  if (isTRUE(ff$converged)) "ok" else "FAILED",
                  if (isTRUE(ff$singular)) "yes (boundary, kept)" else "no",
                  if (isTRUE(ff$rank_deficient)) "DEFICIENT" else "full"))
      cat("  Residual normality and the residual plots are on the diagnostics tab.\n")
    } else {
      cat(sprintf("  R-squared per participant: mean %s, range [%s, %s]\n",
                  fmt3(mean(params$r_squared, na.rm = TRUE)),
                  fmt3(min(params$r_squared, na.rm = TRUE)),
                  fmt3(max(params$r_squared, na.rm = TRUE))))
      n_sig <- sum(params$p_value < 0.05, na.rm = TRUE)
      cat(sprintf("  Detectable rhythms (zero-amplitude F test, harmonics given the trend, p < .05): %d / %d (%s%%)\n",
                  n_sig, nrow(params), fmt1(100 * n_sig / nrow(params))))
      if (!is.null(mod$bingham_summary)) {
        for (h in seq_len(nh)) {
          b <- mod$bingham_summary[[h]]; if (is.null(b)) next
          cat(sprintf("  H%d acrophase identified (joint region excludes the origin) in %d / %d participants (%s%%)\n",
                      h, b$n_identified, b$n, fmt1(100 * b$n_identified / b$n)))
        }
      }
      if (!is.null(mod$boot_results)) {
        br <- mod$boot_results
        cat(sprintf("  Bootstrap CIs (B = %d, %s resampling): constant [%s, %s]; H1 amplitude [%s, %s]; H1 acrophase [%s, %s]\n",
                    br$B, br$resample_unit %||% "participant",
                    fmt3(br$mesor_ci[1]), fmt3(br$mesor_ci[2]),
                    fmt3(br$amplitude_ci[1]), fmt3(br$amplitude_ci[2]),
                    dance_acrophase_label(hours = br$acrophase_ci[1], period = period, harmonic = 1,
                                          clock_origin = clock_o, all = FALSE),
                    dance_acrophase_label(hours = br$acrophase_ci[2], period = period, harmonic = 1,
                                          clock_origin = clock_o, all = FALSE)))
        if (!is.null(br$acrophase_ci_circular) && isTRUE(br$acrophase_ci_circular$wraps))
          cat("    (the acrophase interval wraps past the origin, so its endpoints run high to low)\n")
      }
      im <- pop$indiv_means
      if (!is.null(im) && !is.null(im$unique_S) && is.finite(im$unique_S))
        cat(sprintf("  Variance shares (commonality): trend %s%%, rhythm %s%%, shared %s%% of the total R-squared\n",
                    fmt1(im$percent_S), fmt1(im$percent_C), fmt1(im$percent_shared)))
      cat("  Residual normality and the residual plots are on the diagnostics tab.\n")
    }

    invisible(NULL)
  }

  output$harmonic_summary <- renderPrint({ .print_harmonic_summary() })

  # The table under the text: a spread of the per-participant quantities. Under
  # the two-stage approach these are the independent per-participant fits;
  # under the mixed-effects approach they are the conditional (shrunken)
  # participant estimates of the one model, and the heading says which.
  output$harmonic_parameters_table <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    .co <- dance_clock_origin(mod)
    .to_clock <- function(v) (v + .co) %% mod$period
    rows <- list()
    push <- function(name, x, clock = FALSE) {
      x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
      if (!length(x)) return(invisible(NULL))
      rows[[length(rows) + 1L]] <<- data.frame(
        Parameter = name,
        Mean = round(if (clock) .to_clock(mean(x)) else mean(x), 3),
        SD = round(stats::sd(x), 3),
        Min = round(if (clock) .to_clock(min(x)) else min(x), 3),
        Max = round(if (clock) .to_clock(max(x)) else max(x), 3),
        stringsAsFactors = FALSE)
    }
    if (harmonic_is_mixed()) {
      ff <- harmonic_traj(); req(isTRUE(ff$ok))
      pt <- dance_traj_participant_table(ff); req(isTRUE(pt$ok))
      tb <- pt$table
      push("Level at the first observation", tb$level_at_t0)
      for (tm in ff$spec$trend_terms) push(paste0("Trend (", tm, ")"), tb[[tm]])
      for (h in seq_len(ff$spec$n_harmonics)) {
        push(paste0("Amplitude H", h), tb[[paste0("amplitude_", h)]])
        push(paste0("Acrophase H", h, " (clock h)"), tb[[paste0("acrophase_time_", h)]], clock = TRUE)
      }
      title <- sprintf("Participant estimates (conditional modes of the one model, %d curves)", nrow(tb))
      note <- pt$note
    } else {
      params <- mod$individual_params
      push("Constant term (beta_0)", params$mesor)
      for (h in seq_len(mod$n_harmonics)) {
        push(paste0("Amplitude H", h), params[[paste0("amplitude_", h)]])
        push(paste0("Acrophase H", h, " (clock h)"), params[[paste0("acrophase_time_", h)]], clock = TRUE)
      }
      push("R\u00b2", params$r_squared)
      push("% Rhythm", params$percent_rhythm)
      title <- sprintf("Per-participant fits (%d participants)", nrow(params))
      note <- paste("Arithmetic spread of the independent per-participant estimates. The",
                    "acrophase mean here is arithmetic; the circular (vector) mean is in",
                    "the summary above.")
    }
    summary_df <- do.call(rbind, rows)
    tagList(
      h4(title),
      renderTable(summary_df, striped = TRUE, hover = TRUE, bordered = TRUE),
      helpText(note)
    )
  })
  
  # Fitted curves plot
  output$harmonic_fit_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    time_fine <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 200)
    subject_select <- input$harmonic_subject_select
    if(is.null(subject_select)) subject_select <- "mean"
    
    p <- plot_ly()

    if (harmonic_is_mixed()) {
      # ======================================================================
      # MIXED-EFFECTS: the cell curves ARE the ones tab 6 draws -- the same
      # dance_traj_predict() call on the same fit -- so the two tabs cannot
      # show different trajectories. Participants are conditional predictions.
      # ======================================================================
      ff <- harmonic_traj(); req(isTRUE(ff$ok))
      pr <- dance_traj_predict(ff, conf = 0.95, band = "pointwise",
                               component = "full", n_time = 200)
      req(isTRUE(pr$ok))
      cols <- dance_group_colors(pr$cells)
      sel <- subject_select
      one <- !(sel %in% c("all", "mean"))
      if (identical(sel, "all") || one) {
        pp <- dance_traj_participant_predict(ff, times = pr$times)
        if (isTRUE(pp$ok)) {
          keep <- if (one) sel else unique(pp$table$curve)
          for (cl in keep) {
            dcl <- pp$table[pp$table$curve == cl, , drop = FALSE]
            if (!nrow(dcl)) next
            col <- cols[[dcl$cell[1]]] %||% "#777777"
            p <- p %>% add_lines(
              x = dcl$t, y = dcl$fit,
              name = if (one) sprintf("Participant %s (shrunken)", dcl$subject[1])
                     else sprintf("%s participant", dcl$cell[1]),
              legendgroup = if (one) "participant" else paste0("pp_", dcl$cell[1]),
              showlegend = one,
              line = list(color = dance_group_rgba(col, if (one) 0.95 else 0.30),
                          width = if (one) 2.2 else 1))
          }
        }
      }
      for (cl in pr$cells) {
        d <- pr$table[pr$table$cell == cl, , drop = FALSE]
        if (isTRUE(input$harmonic_show_ci))
          p <- p %>% add_ribbons(x = d$t, ymin = d$lo, ymax = d$hi, name = cl,
                                 legendgroup = cl, line = list(width = 0),
                                 fillcolor = dance_group_rgba(cols[[cl]], 0.16),
                                 showlegend = FALSE, hoverinfo = "skip")
        p <- p %>% add_lines(x = d$t, y = d$fit, name = cl, legendgroup = cl,
                             line = list(color = cols[[cl]], width = 2.6,
                                         dash = if (one) "dot" else "solid"))
      }
      if (isTRUE(input$harmonic_show_data)) {
        dd <- ff$spec$data
        dts <- ff$spec$design_terms
        key <- if (length(dts))
          do.call(paste, c(lapply(dts, function(f) as.character(dd[[f]])), list(sep = " x ")))
          else rep(pr$cells[1], nrow(dd))
        cvn <- if ("curve" %in% names(dd)) "curve" else "subject"
        for (cl in pr$cells) {
          i <- key == cl
          if (one) i <- i & as.character(dd[[cvn]]) == sel
          if (!any(i)) next
          p <- p %>% add_markers(x = dd$t[i], y = dd$y[i], name = cl, legendgroup = cl,
                                 showlegend = FALSE,
                                 marker = list(color = dance_group_rgba(cols[[cl]], 0.35), size = 4))
        }
      }
    } else {
      # ======================================================================
      # TWO-STAGE: one OLS cosinor per participant. The group (or pooled) line
      # and its band come from dance_ts_group_curves() -- the SAME call the
      # comparison tab makes -- so the two tabs draw the same curve: the
      # group's mean-coefficient curve with a pointwise band from the spread
      # of the participants' own fitted curves. Participants are their own
      # fits. (The band this replaces scaled the whole curve by the SE of the
      # first-harmonic amplitude, which is not an interval for anything.)
      # ======================================================================
      gv <- harmonic_group_var_eff()
      grouped <- !is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
        !is.null(gv) && !identical(gv, "_none_")
      gvals <- if (grouped) values$covariates[[gv]] else NULL
      one <- !(subject_select %in% c("all", "mean"))
      # the components of a mean-coefficient curve, on the layout
      # c(constant, trend coefficients..., cos_1, sin_1, ...) predict_from_coefs() uses
      comp_from_coefs <- function(coefs) {
        nh <- mod$n_harmonics; K <- length(coefs) - 1L - 2L * nh
        mesor <- unname(coefs[1]); tr <- if (K > 0) unname(coefs[1 + seq_len(K)]) else numeric(0)
        toff <- mod$t_offset %||% 0
        trend <- switch(mod$trend_type %||% "none",
          linear  = if (K >= 1) tr[1] * time_fine,
          log     = if (K >= 1) tr[1] * log(time_fine - toff + 1),
          exp_sat = if (K >= 2 && is.finite(tr[2]) && tr[2] > 0)
                      tr[1] * (1 - exp(-(time_fine - toff) / tr[2])),
          NULL)
        harm <- lapply(seq_len(nh), function(h) {
          a <- unname(coefs[1 + K + 2 * h - 1]); b <- unname(coefs[1 + K + 2 * h])
          w <- 2 * pi * h / mod$period
          mesor + a * cos(w * time_fine) + b * sin(w * time_fine)
        })
        list(mesor = mesor, trend = if (is.null(trend)) NULL else mesor + trend, harm = harm)
      }

      if (one) {
        i <- as.integer(subject_select)
        fit_i <- mod$individual_fits[[i]]
        if (!is.null(fit_i) && isTRUE(fit_i$success)) {
          p <- p %>% add_lines(x = time_fine, y = predict_cosinor(fit_i, time_fine),
                               line = list(color = DANCE_SERIES1, width = 2),
                               name = paste("Subject", i))
          if (isTRUE(input$harmonic_show_data))
            p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                   marker = list(color = DANCE_SERIES1, size = 6),
                                   name = "Observed data")
          if (isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
            components <- get_harmonic_components(fit_i, time_fine)
            ccol <- dance_component_colors(max(mod$n_harmonics, 1))
            for (h in seq_len(mod$n_harmonics))
              p <- p %>% add_lines(x = time_fine, y = components$mesor[1] + components[[paste0("harmonic_", h)]],
                                   line = list(color = ccol[h], width = 1.5, dash = "dot"),
                                   name = sprintf("H%d (%s h)", h, fmtn(mod$period / h, 1)))
            if (!is.null(components$trend) && mod$trend_type != "none")
              p <- p %>% add_lines(x = time_fine, y = components$mesor[1] + components$trend,
                                   line = list(color = "black", width = 1.5, dash = "dash"),
                                   name = get_trend_label(mod$trend_type))
          }
        } else {
          fail_msg <- paste("Subject", i, "fit failed")
          if (!is.null(fit_i) && !is.null(fit_i$n_valid))
            fail_msg <- paste0(fail_msg, "\n(", fit_i$n_valid, " valid points, need ", fit_i$n_required, "+)")
          if (!is.null(fit_i) && !is.null(fit_i$message))
            fail_msg <- paste0(fail_msg, "\nReason: ", fit_i$message)
          p <- p %>% add_annotations(x = 0.5, y = 0.5, text = fail_msg, showarrow = FALSE,
                                     xref = "paper", yref = "paper",
                                     font = list(size = 14, color = DANCE_ALERT))
        }
      } else {
        gc <- dance_ts_group_curves(mod, gvals, time_fine, include_trend = TRUE, conf = 0.95)
        cells <- if (!is.null(gc)) gc$cells else character(0)
        cols <- if (grouped) dance_group_colors(cells)
                else stats::setNames(as.list(rep(DANCE_EMPHASIS, length(cells))), cells)
        pcol <- function(cl) cols[[cl]] %||% "#777777"
        cell_of <- function(i) {
          if (!grouped) return("(all)")
          cl <- as.character(gvals[i]); if (is.na(cl) || !nzchar(cl)) NA_character_ else cl
        }
        if (identical(subject_select, "all")) {
          for (i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if (is.null(fit_i) || !isTRUE(fit_i$success)) next
            cl <- cell_of(i); if (is.na(cl)) next
            p <- p %>% add_lines(x = time_fine, y = predict_cosinor(fit_i, time_fine),
                                 name = if (grouped) sprintf("%s participant", cl) else "participant",
                                 legendgroup = paste0("pp_", cl), showlegend = FALSE,
                                 line = list(color = dance_group_rgba(pcol(cl), 0.35), width = 1))
          }
        }
        for (cl in cells) {
          d <- gc$table[gc$table$cell == cl, , drop = FALSE]
          if (!any(is.finite(d$fit))) next
          if (isTRUE(input$harmonic_show_ci))
            p <- p %>% add_ribbons(x = d$t, ymin = d$lo, ymax = d$hi, name = cl, legendgroup = cl,
                                   line = list(width = 0),
                                   fillcolor = dance_group_rgba(pcol(cl), 0.16),
                                   showlegend = FALSE, hoverinfo = "skip")
          p <- p %>% add_lines(x = d$t, y = d$fit, legendgroup = cl,
                               name = if (grouped) paste("Group:", cl) else "Population mean",
                               line = list(color = pcol(cl), width = 2.6))
        }
        if (isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          ccol <- dance_component_colors(max(mod$n_harmonics, 1))
          for (cl in cells) {
            coefs <- if (grouped) mod$group_fits[[cl]]$mean_coefs else mod$pop_mean_fit$mean_coefs
            if (is.null(coefs)) next
            cp <- comp_from_coefs(coefs)
            for (h in seq_len(mod$n_harmonics))
              p <- p %>% add_lines(x = time_fine, y = cp$harm[[h]], legendgroup = cl,
                                   name = if (grouped) sprintf("H%d (%s)", h, cl)
                                          else sprintf("H%d (%s h)", h, fmtn(mod$period / h, 1)),
                                   line = list(color = if (grouped) pcol(cl) else ccol[h],
                                               width = 1.5, dash = "dot"))
            if (!is.null(cp$trend))
              p <- p %>% add_lines(x = time_fine, y = cp$trend, legendgroup = cl,
                                   name = if (grouped) sprintf("%s (%s)", get_trend_label(mod$trend_type), cl)
                                          else get_trend_label(mod$trend_type),
                                   line = list(color = if (grouped) pcol(cl) else "black",
                                               width = 1.5, dash = "dash"))
          }
        }
        if (isTRUE(input$harmonic_show_data)) {
          for (i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if (is.null(fit_i) || !isTRUE(fit_i$success)) next
            cl <- cell_of(i); if (is.na(cl)) next
            p <- p %>% add_markers(x = fit_i$time, y = fit_i$y, legendgroup = cl,
                                   showlegend = FALSE, hoverinfo = "skip",
                                   marker = list(color = dance_group_rgba(pcol(cl), 0.3), size = 3))
          }
        }
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
    ticks <- dance_clock_ticks(time_range + shift, mod$period)

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
        sprintf("%s<br>%s", dance_clock_label(xc, mod$period), nm)
      else
        sprintf("%s<br>%s: %s<br>%s", dance_clock_label(xc, mod$period), dvl, fmt2(yv), nm)
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
  # ---------------------------------------------------------------------------
  # THE MIXED-EFFECTS CELL VECTORS, for the polar dial
  # ---------------------------------------------------------------------------
  # The dial drew group means of the PER-PARTICIPANT cosinor fits -- a two-stage
  # estimator -- while tab 6 reported the mixed-effects cell estimates for the
  # same groups, and nothing on screen said which was which. They are close on
  # balanced data with many participants (the two-stage group vector is the
  # amplitude of the MEAN (cos, sin) pair, which is what the mixed model also
  # estimates), close enough to agree to several figures and hide the fact that
  # they are different estimators. They are not close when a group is small,
  # unbalanced, or has participants with few observations -- exactly where it
  # matters -- because only the mixed model shrinks and weights them.
  #
  # Returns ok = FALSE WITH A REASON rather than silently falling back: the dial
  # says which estimator it drew.
  harmonic_polar_mixed <- function(h) {
    ff <- harmonic_traj()
    if (!isTRUE(ff$ok))
      return(list(ok = FALSE, why = ff$message %||% "the mixed-effects model is unavailable"))
    if (h > (ff$spec$n_harmonics %||% 1L))
      return(list(ok = FALSE, why = sprintf("the model fits %d harmonic(s)",
                                            ff$spec$n_harmonics %||% 1L)))
    co <- dance_traj_cell_coefs(ff, h)
    if (!isTRUE(co$ok)) return(list(ok = FALSE, why = co$message))
    ap <- dance_traj_amp_phase_ci(co, method = "joint", n_draw = 6000)
    if (!isTRUE(ap$ok)) return(list(ok = FALSE, why = ap$message))
    list(ok = TRUE, table = ap$table, co = co, fit = ff)
  }

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

    clock_o0 <- dance_clock_origin(mod)
    hov <- function(td) dance_polar_hover_clock(td, mod$period, h,
                                                clock_origin = clock_o0)
    HT <- paste0("amplitude %{r:.3f}<br>acrophase %{text}",
                 "<br><span style='font-size:10px'>%{theta:.1f}&deg; on the H",
                 h, " dial</span><extra>%{fullData.name}</extra>")

    want_mixed <- identical(mod$approach %||% "two_stage", "mixed")
    mx <- if (want_mixed) harmonic_polar_mixed(h) else list(ok = FALSE, why = NULL)
    use_mixed <- isTRUE(mx$ok)

    p <- plot_ly(type = 'scatterpolar', mode = 'markers')
    
    # a 95% joint region for one (cos, sin) pair, as a dotted ring
    ell_trace <- function(mx_, my_, V, nm, col) {
      V <- (V + t(V)) / 2
      ev <- eigen(V, symmetric = TRUE)
      if (any(!is.finite(ev$values)) || any(ev$values <= 0)) return(NULL)
      chi_sq <- stats::qchisq(0.95, df = 2)
      a <- sqrt(chi_sq * ev$values[1]); b <- sqrt(chi_sq * ev$values[2])
      ang <- atan2(ev$vectors[2, 1], ev$vectors[1, 1])
      tt <- seq(0, 2 * pi, length.out = 100)
      ex <- mx_ + a * cos(tt) * cos(ang) - b * sin(tt) * sin(ang)
      ey <- my_ + a * cos(tt) * sin(ang) + b * sin(tt) * cos(ang)
      th <- phi_to_degrees(atan2(ey, ex)); th[th < 0] <- th[th < 0] + 360
      list(r = sqrt(ex^2 + ey^2), theta = th, name = nm, col = col)
    }

    if (use_mixed) {
      # ======================================================================
      # MIXED-EFFECTS. The bold vectors are the model's cell estimates; the
      # cloud is each participant's SHRUNKEN rhythm (their cell's vector plus
      # their conditional modes); the dotted ring is the cell's 95% joint
      # region for its (cos, sin) pair -- Bingham's region, from the fixed-
      # effect covariance. Nothing on this dial comes from a per-participant fit.
      # ======================================================================
      cells <- mx$table$cell
      group_colors <- dance_group_colors(cells)
      pc <- dance_traj_participant_curves(mx$fit, h)
      for (g_name in cells) {
        col <- unname(group_colors[[g_name]])
        if (isTRUE(pc$ok)) {
          sub <- pc$table[pc$table$cell == g_name, , drop = FALSE]
          if (nrow(sub)) {
            td <- phi_to_degrees(sub$acrophase_rad)
            p <- p %>% add_trace(
              r = sub$amplitude, theta = td, type = 'scatterpolar', mode = 'markers',
              marker = list(size = 7, color = col, opacity = 0.55),
              name = paste("Participants (shrunken):", g_name), legendgroup = g_name,
              text = hov(td), hovertemplate = HT)
          }
        }
        i_mx <- match(g_name, mx$table$cell)
        acro_deg <- phi_to_degrees(mx$table$acrophase_rad[i_mx])
        if (acro_deg < 0) acro_deg <- acro_deg + 360
        p <- p %>% add_trace(
          r = c(0, mx$table$amplitude[i_mx]), theta = c(0, acro_deg),
          type = 'scatterpolar', mode = 'lines+markers',
          line = list(color = col, width = 3),
          marker = list(size = 12, color = col, symbol = 'diamond'),
          name = paste("Cell:", g_name), legendgroup = g_name,
          text = c(hov(0), hov(acro_deg)), hovertemplate = HT)
        if (isTRUE(input$polar_show_ellipse) && !is.null(mx$co$joint)) {
          ic <- match(g_name, mx$co$cells)
          el <- if (!is.na(ic)) ell_trace(mx$co$a[ic], mx$co$b[ic], mx$co$joint[[ic]],
                                          paste("95% region:", g_name), col) else NULL
          if (!is.null(el))
            p <- p %>% add_trace(r = el$r, theta = el$theta, type = 'scatterpolar',
                                 mode = 'lines', line = list(color = col, width = 1.5, dash = 'dot'),
                                 name = el$name, legendgroup = g_name, showlegend = FALSE,
                                 hoverinfo = "skip")
        }
      }
    } else if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
       !is.null(harmonic_group_var_eff()) && harmonic_group_var_eff() != "_none_") {

      group_var <- values$covariates[[harmonic_group_var_eff()]]
      group_colors <- dance_group_colors(names(mod$group_fits))

      groups <- names(mod$group_fits)
      for(g_idx in seq_along(groups)) {
        g_name <- groups[g_idx]
        g_mask <- !is.na(group_var[params$subject]) & group_var[params$subject] == g_name

        if(sum(g_mask) > 0) {
          p <- p %>% add_trace(
            r = r[g_mask], theta = theta_deg[g_mask],
            type = 'scatterpolar', mode = 'markers',
            marker = list(size = 8, color = unname(group_colors[[g_name]]), opacity = 0.7),
            name = paste("Group:", g_name),
            text = hov(theta_deg[g_mask]), hovertemplate = HT
          )
        }

        # Add group mean vector for selected harmonic, from whichever estimator
        # the reader asked for. Both use atan2(beta_sin, beta_cos) on a
        # zero-based model axis, so the angle means the same thing either way
        # and only the ESTIMATE changes -- no rotation, no reconversion.
        g_fit <- mod$group_fits[[g_name]]
        g_amp <- g_fit$mean_amplitudes[h]
        g_rad <- g_fit$mean_acrophases_rad[h]
        if (use_mixed) {
          i_mx <- match(g_name, mx$table$cell)
          if (!is.na(i_mx)) {
            g_amp <- mx$table$amplitude[i_mx]
            g_rad <- mx$table$acrophase_rad[i_mx]
          }
        }
        acro_deg <- phi_to_degrees(g_rad)
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, g_amp),
          theta = c(0, acro_deg),
          type = 'scatterpolar', mode = 'lines+markers',
          line = list(color = unname(group_colors[[g_name]]), width = 3),
          marker = list(size = 12, color = unname(group_colors[[g_name]]), symbol = 'diamond'),
          name = paste(if (use_mixed) "Mixed:" else "Mean:", g_name),
          text = c(hov(0), hov(acro_deg)), hovertemplate = HT
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
          name = "Overall Population Mean",
          text = c(hov(0), hov(acro_deg)), hovertemplate = HT
        )
      }

    } else {
      # No groups - show all points same color
      p <- p %>% add_trace(
        r = r, theta = theta_deg,
        marker = list(size = 8, color = DANCE_SERIES1, opacity = 0.7),
        name = "Individual",
        text = hov(theta_deg), hovertemplate = HT
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
          line = list(color = DANCE_EMPHASIS, width = 3),
          marker = list(size = 12, color = DANCE_EMPHASIS, symbol = 'diamond'),
          name = "Population Mean",
          text = c(hov(0), hov(acro_deg)), hovertemplate = HT
        )
      }
    }

    # Add confidence ellipse if requested (for all data, regardless of groups).
    # Two-stage only: it is the sample ellipse of the per-participant points,
    # and under the mixed approach each cell already carries its own region.
    if(!use_mixed && isTRUE(input$polar_show_ellipse) && length(r) >= 3) {
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
          line = list(color = dance_group_rgba(DANCE_NEUTRAL, 0.75), width = 2, dash = 'dot'),
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
    clock_o <- dance_clock_origin(mod)
    effective_period <- mod$period / h
    n_ticks <- min(12, effective_period)
    tick_step <- effective_period / n_ticks
    tick_vals <- seq(0, 360 - 360/n_ticks, by = 360/n_ticks)
    tick_elapsed <- seq(0, effective_period - tick_step, by = tick_step)
    tick_labels <- vapply(tick_elapsed, function(e)
      dance_clock_label((e + clock_o) %% mod$period, mod$period, show_day = FALSE),
      character(1))

    # The title sat on top of the dial: a polar trace fills its plotting area
    # edge to edge, so a centred title with no reserved space lands on the 11-13
    # o'clock labels. The subtitle is moved out to a paper-anchored annotation
    # under the title and the polar domain is pulled down to leave room, rather
    # than shrinking the font until it stops colliding.
    #
    # AUDIT (second pass): three things were still stacked at 12 o'clock. The
    # radial axis was pinned at angle = 90 -- the top of the dial, once
    # rotation = 90 and a clockwise direction are applied -- so its "Amplitude"
    # title and its tick labels were drawn straight through the topmost angular
    # label and up into the subtitle. The 12 angular ticks sit at 30 deg steps,
    # so the radial axis moves to 45 deg: exactly midway between two of them,
    # where nothing else is ever drawn. The legend was on plotly's default right
    # edge, where a long group name eats into the dial; it goes underneath, as
    # on the density tab, and the bottom margin makes room for it instead of the
    # dial giving up width. uirevision keeps the reader's zoom and legend state
    # across a re-render, as on the density tab.
    polar_sub <- if(h > 1)
        sprintf("clock times; each angle recurs every %s h", fmt1(effective_period))
      else if(clock_o != 0)
        sprintf("clock times; the model origin is %s",
                dance_clock_label(clock_o, mod$period, show_day = FALSE))
      else "clock times"
    # WHICH ESTIMATOR DREW THE BOLD VECTORS. Not a footnote: the point cloud and
    # the group vectors come from different places, and a dial that does not say
    # so is how a figure quietly disagrees with the table beside it.
    polar_sub <- paste0(polar_sub, " &middot; group vectors: ",
      if (use_mixed) "mixed-effects model"
      else if (want_mixed && !is.null(mx$why))
        sprintf("two-stage (mixed-effects unavailable -- %s)", mx$why)
      else "two-stage, mean of participant fits")

    p %>% layout(
      uirevision = "fck-acrophase-polar",
      # The subtitle is a second line of the TITLE, not a free-floating
      # annotation. As an annotation at a paper y it had no idea where the dial
      # ended, and the topmost angular label -- which plotly draws OUTSIDE the
      # polar domain -- was written straight through it. As part of the title
      # plotly stacks and reserves both lines itself.
      title = list(
        text = sprintf("Acrophase polar plot - H%d (effective period %s h)<br><sub>%s</sub>",
                       h, fmt1(effective_period), polar_sub),
        x = 0.5, xanchor = "center", y = 0.97, yanchor = "top",
        font = list(size = 15)),
      polar = list(
        # Top edge clears the two title lines AND the angular labels drawn
        # outside the dial; bottom edge clears the legend under it.
        domain = list(y = c(0.10, 0.83)),
        radialaxis = list(title = list(text = "Amplitude", font = list(size = 11)),
                          tickangle = 0, angle = 45, tickfont = list(size = 10)),
        angularaxis = list(
          direction = "clockwise",
          rotation = 90,
          tickmode = "array",
          tickvals = tick_vals,
          ticktext = tick_labels
        )
      ),
      legend = list(orientation = "h", x = 0.5, xanchor = "center",
                    y = -0.02, yanchor = "top", font = list(size = 11),
                    uirevision = "fck-acrophase-polar-legend"),
      margin = list(t = 88, b = 78),
      showlegend = TRUE
    )
  })
  
  
  
  
  
  # Individual results table
  # The mixed model's per-participant rows, shaped like the two-stage table so
  # the two approaches read alike: one row per curve, the cell it belongs to,
  # level and trend, and each harmonic's amplitude and clock acrophase. A column
  # with no random effect in the fitted structure is the CELL value for every
  # participant and its header says so.
  harmonic_mixed_participant_display <- function() {
    ff <- harmonic_traj(); if (!isTRUE(ff$ok)) return(NULL)
    pt <- dance_traj_participant_table(ff); if (!isTRUE(pt$ok)) return(NULL)
    tb <- pt$table; sp <- ff$spec
    disp <- data.frame(Subject = tb$subject, stringsAsFactors = FALSE)
    if (any(tb$curve != tb$subject)) disp$Curve <- tb$curve
    disp$Cell <- tb$cell
    tag <- function(nm, ok) if (isTRUE(ok)) nm else paste0(nm, " (cell)")
    disp[[tag("Intercept_b0", pt$varies[["intercept"]])]] <- round(tb$intercept, 3)
    disp[["Level_at_t0"]] <- round(tb$level_at_t0, 3)
    for (tm in sp$trend_terms)
      disp[[tag(tm, pt$varies[[tm]])]] <- round(tb[[tm]], 4)
    for (h in seq_len(sp$n_harmonics)) {
      ok_h <- pt$varies[[paste0("H", h)]]
      disp[[tag(paste0("Amp_H", h), ok_h)]] <- round(tb[[paste0("amplitude_", h)]], 3)
      disp[[tag(paste0("Acro_H", h, "_clock"), ok_h)]] <- dance_clock_label(
        dance_acrophase_clock(hours = tb[[paste0("acrophase_time_", h)]],
                              period = sp$period, harmonic = h,
                              clock_origin = ff$clock_origin %||% 0)$hours,
        sp$period, show_day = FALSE)
    }
    list(display = disp, raw = tb, note = pt$note, re_label = pt$re_label)
  }

  output$harmonic_individual_table <- DT::renderDataTable({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if (harmonic_is_mixed()) {
      md <- harmonic_mixed_participant_display(); req(!is.null(md))
      return(DT::datatable(
        md$display, rownames = FALSE,
        caption = tags$caption(style = "caption-side:top;text-align:left;font-size:12px;color:#555",
          sprintf("Shrunken (conditional-mode) estimates from the mixed model; random structure: %s. %s",
                  md$re_label %||% "", md$note)),
        options = list(pageLength = 15, scrollX = TRUE)))
    }
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
    if(!is.null(harmonic_group_var_eff()) && harmonic_group_var_eff() != "_none_") {
      group_var <- values$covariates[[harmonic_group_var_eff()]]
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
      if (harmonic_is_mixed()) {
        md <- harmonic_mixed_participant_display(); req(!is.null(md))
        out <- md$raw
        out$estimator <- "mixed-effects conditional modes (shrunken)"
        write.csv(out, file, row.names = FALSE)
        return(invisible(NULL))
      }
      params <- mod$individual_params
      
      # Add group variable if available
      if(!is.null(harmonic_group_var_eff()) && harmonic_group_var_eff() != "_none_") {
        group_var <- values$covariates[[harmonic_group_var_eff()]]
        params$group <- group_var[params$subject]
      }
      
      write.csv(params, file, row.names = FALSE)
    }
  )
  
  # Residual plot
  # The residuals of whichever model was fitted: the one mixed model's, or the
  # pooled per-participant ones. The two are not comparable -- a mixed-model
  # residual is measured from a participant's OWN conditional curve, so it is
  # smaller than a two-stage residual by construction -- and the panel names
  # which it is showing.
  harmonic_residuals <- function() {
    mod <- values$harmonic_model
    if (harmonic_is_mixed()) {
      ff <- harmonic_traj(); if (!isTRUE(ff$ok)) return(NULL)
      return(list(fitted = as.numeric(stats::fitted(ff$model)),
                  resid = as.numeric(stats::resid(ff$model)),
                  source = "the mixed-effects model (conditional residuals)"))
    }
    all_fitted <- c(); all_resid <- c()
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_fitted <- c(all_fitted, fit_i$fitted)
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    list(fitted = all_fitted, resid = all_resid,
         source = "the per-participant fits, pooled")
  }

  output$harmonic_residual_plot <- renderPlotly({
    req(values$harmonic_model)
    rs <- harmonic_residuals(); req(!is.null(rs))
    all_fitted <- rs$fitted; all_resid <- rs$resid
    
    plot_ly(x = all_fitted, y = all_resid, type = 'scatter', mode = 'markers',
            marker = list(color = DANCE_SERIES1, opacity = 0.5, size = 4)) %>%
      add_segments(x = min(all_fitted), xend = max(all_fitted), y = 0, yend = 0,
                   line = list(color = DANCE_NEUTRAL, dash = 'dash')) %>%
      layout(title = "Residuals vs Fitted",
             xaxis = list(title = "Fitted Values"),
             yaxis = list(title = "Residuals"))
  })
  
  # QQ plot
  output$harmonic_qq_plot <- renderPlotly({
    req(values$harmonic_model)
    rs <- harmonic_residuals(); req(!is.null(rs))
    all_resid <- rs$resid
    
    qq <- qqnorm(all_resid, plot.it = FALSE)
    
    plot_ly(x = qq$x, y = qq$y, type = 'scatter', mode = 'markers',
            marker = list(color = DANCE_SERIES1, size = 4)) %>%
      add_lines(x = range(qq$x), y = range(qq$x) * sd(all_resid) + mean(all_resid),
                line = list(color = DANCE_NEUTRAL, dash = 'dash'), name = "Reference") %>%
      layout(title = "Q-Q Plot of Residuals",
             xaxis = list(title = "Theoretical Quantiles"),
             yaxis = list(title = "Sample Quantiles"))
  })
  
  # Goodness of fit statistics
  output$harmonic_gof_stats <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    rs <- harmonic_residuals(); req(!is.null(rs))
    all_resid <- rs$resid

    cat("=== Residual Diagnostics ===\n")
    cat("Residuals from ", rs$source, ".\n\n", sep = "")
    if (harmonic_is_mixed()) {
      ff <- harmonic_traj()
      m <- ff$model
      cat(sprintf("Observations: %d   participants: %d   curves: %d\n",
                  length(all_resid), ff$spec$n_participants, ff$spec$n_curves))
      cat(sprintf("Random structure: %s (rung %d of %d)\n", ff$re_label, ff$re_rung, ff$n_rungs))
      cat(sprintf("Convergence: %s; singular: %s; fixed-effect rank: %s\n",
                  if (isTRUE(ff$converged)) "converged" else "DID NOT CONVERGE",
                  if (isTRUE(ff$singular)) "yes (boundary, kept)" else "no",
                  if (isTRUE(ff$rank_deficient)) "DEFICIENT" else "full"))
      cat(sprintf("Residual SD (sigma): %.4f\n", stats::sigma(m)))
      r2 <- dance_traj_r2(ff)
      if (!is.null(r2))
        cat(sprintf("R-squared: marginal %.3f (fixed effects), conditional %.3f (fixed + random)\n  [%s]\n",
                    r2$marginal, r2$conditional, r2$method))
      ll <- tryCatch(stats::logLik(m), error = function(e) NULL)
      if (!is.null(ll))
        cat(sprintf("logLik %.2f (%s, df = %d)   AIC %.2f   BIC %.2f\n",
                    as.numeric(ll), if (isTRUE(ff$REML)) "REML" else "ML", attr(ll, "df"),
                    stats::AIC(m), stats::BIC(m)))
      cat("\n")
    }
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

  # ==========================================================================
  # TAB 6: THE MIXED-EFFECTS TRAJECTORY COMPARISON (P21 phase 4)
  # ==========================================================================
  # ONE fit, in one reactive, feeding every panel below it. That is the point of
  # the redesign: the omnibus, the component table, the curves, the difference
  # curve and the pairwise contrasts are all views of the SAME fitted model, so
  # they cannot disagree with each other the way separate procedures can.
  #
  # It is an eventReactive on the Run button, not a plain reactive: refitting a
  # mixed model on every keystroke in the design panel would be unusable.
  # IT DEPENDS ON THE FITTED MODEL, NOT ON THE BUTTON.
  # ------------------------------------------------------------------------
  # The first version keyed off input$run_harmonic and rebuilt the model frame
  # from scratch -- including the time vector, which it read from a
  # `values$time_points` that DOES NOT EXIST anywhere in this app. It was a name
  # I invented, so it was always NULL, and the fallback quietly used COLUMN
  # INDICES as the time variable. Every acrophase, amplitude and contrast in tab
  # 6 was then computed against 1..ncol(Y) while the period said 24, which is
  # not a display problem: it is a different analysis, and it looked plausible
  # because a 24 h period over a 16-point index range just draws a monotone arc.
  #
  # Keying off values$harmonic_model instead makes the module have ONE model.
  # The time vector, period, harmonic count, trend and clock origin all come
  # from the object the rest of the app already built and reports, so tab 6
  # cannot describe a different fit from tab 1.
  # THE FITTER, as a function: called once from the Run button under the
  # mixed-effects approach and stored on the model object. A design with no
  # factor is allowed -- one trajectory for the whole sample is a legitimate
  # mixed model; only the comparison panels need two cells.
  fit_harmonic_traj <- function(mod) {
    if (is.null(mod)) return(list(ok = FALSE, message = "Press Run Harmonic Regression first."))
    dt <- harmonic_design_terms()
    chosen <- c(dt$between, dt$within)
    Y <- values$data
    if (is.null(Y)) return(list(ok = FALSE, message = "No data."))
    # the model-elapsed axis the fit was built on, and the offset back to clock
    tv <- mod$time_vec
    if (is.null(tv) || length(tv) != ncol(Y))
      return(list(ok = FALSE, message = sprintf(
        "The model's time vector has %s values for %d columns -- re-run the fit.",
        if (is.null(tv)) "no" else length(tv), ncol(Y))))
    subj <- values$subject_ids %||% rownames(Y) %||% as.character(seq_len(nrow(Y)))

    fl <- list()
    for (f in chosen) {
      v <- values$covariates[[f]]
      if (is.null(v) || length(v) != nrow(Y))
        return(list(ok = FALSE, message = sprintf(
          "'%s' has %s values for %d curves -- it cannot be used as a design factor.",
          f, if (is.null(v)) "no" else length(v), nrow(Y))))
      fl[[f]] <- v
    }
    # The roles the USER chose are passed as an override, so a deliberate choice
    # wins over the data-derived guess -- and dance_traj_classify records it as
    # an override rather than pretending the data said it.
    roles <- c(stats::setNames(rep("between", length(dt$between)), dt$between),
               stats::setNames(rep("within",  length(dt$within)),  dt$within))

    d <- try(dance_traj_long(Y, tv, subj, fl), silent = TRUE)
    if (inherits(d, "try-error"))
      return(list(ok = FALSE, message = paste("Could not build the model frame:",
                                              conditionMessage(attr(d, "condition")))))
    trend <- mod$trend_type %||% "none"
    tau <- if (identical(trend, "exp_sat"))
      suppressWarnings(as.numeric(input$harmonic_tau_fixed %||% NA)) else NULL
    if (identical(trend, "exp_sat") && !is.finite(tau))
      return(list(ok = FALSE, message = paste(
        "A saturating trend needs a value for tau. Set one under Advanced, or",
        "choose a different trend: tau is nonlinear and cannot be read off the",
        "design matrix the way the other trends can.")))

    withProgress(message = "Fitting the mixed-effects cosinor...", value = 0.3, {
      sp <- dance_traj_spec(d, period = mod$period %||% 24,
                            n_harmonics = mod$n_harmonics %||% 1,
                            trend = trend, tau = tau, roles = roles)
      if (!isTRUE(sp$ok)) return(list(ok = FALSE, message = sp$message))
      incProgress(0.4)
      ff <- dance_traj_fit(sp)
      if (!isTRUE(ff$ok)) return(list(ok = FALSE, message = ff$message))
      incProgress(0.3)
      # the offset that turns this model's elapsed axis back into clock time,
      # carried on the fit so every panel converts the same way
      ff$clock_origin <- dance_clock_origin(mod)
      ff
    })
  }

  # What every mixed-effects panel reads. Under the two-stage approach there is
  # no trajectory model and the panels say so; nothing here fits anything.
  harmonic_traj <- reactive({
    mod <- values$harmonic_model
    if (is.null(mod)) return(list(ok = FALSE, message = "Press Run Harmonic Regression first."))
    if (!identical(mod$approach %||% "two_stage", "mixed"))
      return(list(ok = FALSE, message = paste(
        "The two-stage approach was run. The mixed-effects trajectory model is",
        "not fitted under it -- choose 'Mixed-effects cosinor' under Approach",
        "and run again.")))
    mod$traj %||% list(ok = FALSE, message = "The mixed-effects model was not fitted.")
  })

  # TRUE when the fitted model is the mixed-effects one
  harmonic_is_mixed <- reactive({
    identical(values$harmonic_model$approach %||% "two_stage", "mixed")
  })

  # ==========================================================================
  # TAB 6 UNDER THE TWO-STAGE APPROACH
  # ==========================================================================
  # The same panels, the same questions, answered from the per-participant
  # estimates: MANOVA on the coefficient vectors for "do they differ at all",
  # ANOVA / Bingham / Watson-Williams for the components, Welch and
  # Watson-Williams pairwise. Every renderer below starts by asking which
  # approach was run and hands off to these when it was two-stage, so the tab
  # never mixes estimators.
  harmonic_ts <- reactive({
    mod <- values$harmonic_model
    if (is.null(mod) || harmonic_is_mixed()) return(NULL)
    gv <- harmonic_group_var_eff()
    gvals <- if (!is.null(gv) && !identical(gv, "_none_")) values$covariates[[gv]] else NULL
    df <- dance_ts_frame(mod, gvals)
    if (is.null(df)) return(list(ok = FALSE, message = "No converged per-participant fits to compare."))
    dt <- harmonic_design_terms()
    list(ok = TRUE, mod = mod, df = df, gv = gv, gvals = gvals,
         groups = levels(df$group), n_groups = nlevels(df$group),
         design = dt$design, pairing_ignored = !identical(dt$design %||% "between", "between"))
  })

  ts_stat_td <- function(r) {
    if (!isTRUE(r$ok)) return(tags$td(style = "padding:3px 12px 3px 0;color:#8a5a12;font-size:12px", r$message))
    tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
            sprintf("F(%g, %.1f) = %.2f", r$df1, r$df2, r$F))
  }

  ts_header <- function() {
    ts <- harmonic_ts(); if (is.null(ts)) return(NULL)
    if (!isTRUE(ts$ok)) return(div(class = "alert alert-warning", ts$message))
    mod <- ts$mod; fa <- mod$fit_audit
    div(style = "background:#f7f7f7;border:1px solid #e3e3e3;border-radius:3px;padding:8px 10px;margin-bottom:10px",
        tags$div(style = "font-family:monospace;font-size:11px;color:#555",
                 sprintf("two-stage: y_i(t) = %s, one fit per curve", dance_model_label(mod$trend_type %||% "none", mod$n_harmonics %||% 1L))),
        tags$div(style = "font-size:12px;color:#555;margin-top:4px",
                 sprintf("%d curve(s) compared (%s converged of %s attempted), %s.",
                         nrow(ts$df), fmtn(fa$n_converged %||% nrow(ts$df), 0),
                         fmtn(fa$n_attempted %||% nrow(mod$individual_params), 0),
                         if (ts$n_groups >= 2) sprintf("%d groups by %s", ts$n_groups, ts$gv)
                         else "no grouping factor -- nothing to compare")),
        tags$div(style = "font-size:11px;color:#8a5a12;margin-top:6px",
                 paste("Each participant's cosinor is fitted independently and the POINT",
                       "ESTIMATES are compared; every estimate is treated as exact, so a",
                       "participant whose rhythm is poorly determined counts as much as one",
                       "whose rhythm is precise.",
                       if (isTRUE(ts$pairing_ignored))
                         "This design has a within-participant factor: two-stage cannot respect the pairing, and a participant's curves enter as if independent. Use the mixed-effects approach for that design."
                       else "")))
  }

  ts_primary <- function() {
    ts <- harmonic_ts(); if (is.null(ts)) return(NULL)
    if (!isTRUE(ts$ok) || ts$n_groups < 2)
      return(helpText("Select a design factor with two or more levels to compare groups."))
    r <- dance_ts_manova(ts$df, dance_ts_coef_cols(ts$mod))
    tagList(tags$table(style = "margin:4px 0 6px 0", tags$tbody(tags$tr(
      tags$td(style = "padding:4px 12px 4px 0;font-weight:600", ts$gv),
      ts_stat_td(r),
      if (isTRUE(r$ok)) tags$td(style = "padding:4px 0;font-family:monospace",
                                sprintf("p %s", dance_fmt_p(r$p)))))),
      tags$div(style = "font-size:11px;color:#777",
               if (isTRUE(r$ok)) sprintf(paste("%s on the %d per-participant coefficients (constant, trend, and",
                                              "each harmonic's cosine and sine); Wilks' lambda = %.3f. The",
                                              "two-stage analogue of the full-trajectory test."),
                                        r$method, r$n_par, r$wilks)
               else "The joint test could not be run."))
  }

  ts_components <- function() {
    ts <- harmonic_ts(); if (is.null(ts)) return(NULL)
    if (!isTRUE(ts$ok) || ts$n_groups < 2) return(NULL)
    mod <- ts$mod; period <- mod$period %||% 24
    comps <- dance_ts_components(ts$df, mod)
    rows <- lapply(comps, function(r) tags$tr(
      tags$td(style = "padding:3px 12px 3px 0", r$label),
      tags$td(style = "padding:3px 12px 3px 0;color:#777;font-size:11px", r$method %||% ""),
      ts_stat_td(r),
      tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
              if (isTRUE(r$ok)) dance_fmt_p(r$p) else ""),
      tags$td(style = "padding:3px 0;font-size:11px;color:#8a5a12", r$note %||% "")))
    # derived per group, per harmonic: the vector-mean rhythm, in clock time
    drows <- list()
    for (h in seq_len(mod$n_harmonics %||% 1L)) {
      gr <- dance_ts_group_rhythm(ts$df, mod, h)
      if (is.null(gr)) next
      for (i in seq_len(nrow(gr))) {
        drows[[length(drows) + 1L]] <- tags$tr(
          tags$td(style = "padding:3px 12px 3px 0", sprintf("H%d amplitude — %s", h, gr$cell[i])),
          tags$td(style = "padding:3px 12px 3px 0;color:#777", "vector mean"),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace", sprintf("%.3f", gr$amplitude[i])),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace;color:#777", sprintf("n = %d", gr$n[i])),
          tags$td(""))
        drows[[length(drows) + 1L]] <- tags$tr(
          tags$td(style = "padding:3px 12px 3px 0", sprintf("H%d acrophase — %s", h, gr$cell[i])),
          tags$td(style = "padding:3px 12px 3px 0;color:#777", "circular"),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  dance_clock_label(dance_acrophase_clock(hours = gr$acrophase_time[i], period = period,
                                                          harmonic = h, clock_origin = dance_clock_origin(mod))$hours,
                                    period, show_day = FALSE)),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace;color:#777",
                  if (is.finite(gr$r_bar[i])) sprintf("r-bar %.2f", gr$r_bar[i]) else ""),
          tags$td(""))
      }
    }
    for (gn in ts$groups) {
      gf <- mod$group_fits[[gn]]
      if (is.null(gf) || is.null(gf$mean_coefs)) next
      pk <- dance_curve_peak_clock(gf$mean_coefs, mod)
      if (is.null(pk)) next
      drows[[length(drows) + 1L]] <- tags$tr(
        tags$td(style = "padding:3px 12px 3px 0", sprintf("Fitted curve peak — %s", gn)),
        tags$td(style = "padding:3px 12px 3px 0;color:#777", "plotted"),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                dance_clock_label(pk$peak_clock, period, show_day = FALSE)),
        tags$td(style = sprintf("padding:3px 12px 3px 0;font-family:monospace;color:%s",
                                if (isTRUE(pk$peak_at_edge)) "#c0392b" else "#777"),
                sprintf("value %.2f%s", pk$peak_value,
                        if (isTRUE(pk$peak_at_edge)) " — at the window edge, not a peak" else "")),
        tags$td(""))
    }
    tagList(tags$table(style = "margin:4px 0", tags$tbody(rows, drows)),
            tags$div(style = "font-size:11px;color:#777;margin-top:6px", HTML(paste(
              "Component tests on the per-participant estimates: one-way ANOVA for the",
              "scalars, Bingham et al.'s (1982) population-mean cosinor for each harmonic's",
              "amplitude and acrophase (with the joint MANOVA on the (cos, sin) vector, which",
              "has no half-cycle blind spot), and Watson-Williams on the unweighted phases.",
              "Derived rows are the amplitude-weighted vector means per group, the same",
              "estimator as the polar dial. <b>An H1 acrophase is not the peak of the drawn",
              "curve</b>; the fitted-curve-peak rows give that."))))
  }

  ts_group_curves <- function(component = "full") {
    ts <- harmonic_ts(); if (is.null(ts) || !isTRUE(ts$ok)) return(NULL)
    mod <- ts$mod
    tv <- mod$time_vec
    times <- seq(min(tv, na.rm = TRUE), max(tv, na.rm = TRUE), length.out = 160)
    dance_ts_group_curves(mod, ts$gvals, times,
                          include_trend = !identical(component, "rhythm"))
  }

  ts_curves_plot <- function() {
    ts <- harmonic_ts(); if (is.null(ts)) return(NULL)
    if (!isTRUE(ts$ok)) return(plotly_empty() %>% layout(title = list(text = ts$message, font = list(size = 12))))
    gc <- ts_group_curves(harmonic_traj_component())
    if (is.null(gc)) return(plotly_empty() %>% layout(title = list(text = "Fewer than two converged fits in every group.", font = list(size = 12))))
    mod <- ts$mod
    cols <- dance_group_colors(gc$cells)
    p <- plot_ly()
    for (cl in gc$cells) {
      d <- gc$table[gc$table$cell == cl, , drop = FALSE]
      p <- p %>%
        add_ribbons(x = d$t, ymin = d$lo, ymax = d$hi, name = cl, legendgroup = cl,
                    line = list(width = 0), fillcolor = dance_group_rgba(cols[[cl]], 0.16),
                    showlegend = FALSE, hoverinfo = "skip") %>%
        add_lines(x = d$t, y = d$fit, name = cl, legendgroup = cl,
                  line = list(color = cols[[cl]], width = 2.4))
    }
    if (isTRUE(input$harmonic_traj_raw) && identical(harmonic_traj_component(), "full")) {
      Y <- mod$Y; g <- ts$df$group; sub <- ts$df$subject
      for (cl in gc$cells) {
        rows <- sub[g == cl]
        if (!length(rows)) next
        yy <- as.numeric(t(Y[rows, , drop = FALSE])); xx <- rep(tv, length(rows))
        ok <- is.finite(yy)
        p <- p %>% add_markers(x = xx[ok], y = yy[ok], name = cl, legendgroup = cl, showlegend = FALSE,
                               marker = list(color = dance_group_rgba(cols[[cl]], 0.35), size = 4))
      }
    }
    co <- dance_clock_origin(mod); P <- mod$period %||% 24
    brk <- pretty(range(gc$times), n = 8); brk <- brk[brk >= min(gc$times) & brk <= max(gc$times)]
    p %>% layout(
      xaxis = list(title = if (co != 0) "Clock time" else "Time", tickmode = "array", tickvals = brk,
                   ticktext = dance_clock_label(brk + co, P, show_day = TRUE, with_minutes = FALSE)),
      yaxis = list(title = if (identical(harmonic_traj_component(), "rhythm")) "Rhythm (trend removed)"
                   else if (!is.null(mod$dv_name)) paste0(mod$dv_name, if (!is.null(mod$dv_units)) paste0(" (", mod$dv_units, ")") else "") else "Response"),
      hovermode = "x unified", legend = list(orientation = "h", y = -0.18))
  }

  ts_diff <- reactive({
    ts <- harmonic_ts(); if (is.null(ts) || !isTRUE(ts$ok)) return(list(ok = FALSE, message = ts$message %||% ""))
    a <- input$harmonic_traj_diff_a; b <- input$harmonic_traj_diff_b
    req(a, b)
    gc <- ts_group_curves("full")
    if (is.null(gc) || !all(c(a, b) %in% gc$cells))
      return(list(ok = FALSE, message = "Both groups need at least two converged fits."))
    da <- gc$table[gc$table$cell == a, ]; db <- gc$table[gc$table$cell == b, ]
    z <- stats::qnorm(0.975)
    se <- sqrt(da$se^2 + db$se^2)
    list(ok = TRUE, cell1 = a, cell2 = b,
         table = data.frame(t = da$t, diff = da$fit - db$fit, se = se,
                            lo = da$fit - db$fit - z * se, hi = da$fit - db$fit + z * se),
         band = "pointwise",
         note = paste("Difference of the two groups' mean-coefficient curves with a POINTWISE",
                      "95% band from the two standard errors added in quadrature",
                      "(independent groups). Not simultaneous: it cannot be read across",
                      "the whole curve as one statement."))
  })

  ts_pairwise <- function() {
    ts <- harmonic_ts(); if (is.null(ts)) return(NULL)
    if (!isTRUE(ts$ok) || ts$n_groups < 2) return(NULL)
    what <- input$harmonic_traj_what %||% "level"
    if (!what %in% dance_ts_pair_whats(ts$mod)) what <- "level"
    adj <- input$harmonic_traj_adjust %||% "holm"
    r <- dance_ts_pairwise(ts$df, ts$mod, what, adjust = adj)
    if (!isTRUE(r$ok)) return(tags$div(style = "color:#8a5a12;font-size:12px", r$message))
    tb <- r$table
    hdr <- if (r$circular) c("Pair", "Difference (h)", "Watson-Williams", "p", "p adj")
           else c("Pair", "Estimate", "95% CI", "p", "p adj", "d")
    rows <- lapply(seq_len(nrow(tb)), function(i) {
      z <- tb[i, ]
      tags$tr(
        tags$td(style = "padding:3px 12px 3px 0", sprintf("%s vs %s", z$cell1, z$cell2)),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace", sprintf("%+.3f", z$estimate)),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                if (r$circular) sprintf("F(1, %.0f) = %.2f%s", z$df2, z$statistic,
                                        if (!isTRUE(z$assumption_ok)) " (r-bar < 0.45)" else "")
                else sprintf("[%.3f, %.3f]", z$lo, z$hi)),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace", dance_fmt_p(z$p_raw)),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace;font-weight:600", dance_fmt_p(z$p_adj)),
        if (!r$circular) tags$td(style = "padding:3px 0;font-family:monospace;color:#777",
                                 if (is.finite(z$d)) sprintf("%.2f", z$d) else ""))
    })
    tagList(
      tags$table(style = "margin:4px 0",
        tags$thead(tags$tr(lapply(hdr, function(h)
          tags$th(style = "text-align:left;padding:2px 12px 2px 0;font-size:11px;color:#777;font-weight:500", h)))),
        tags$tbody(rows)),
      tags$div(style = "font-size:11px;color:#777", r$unit),
      tags$div(style = "font-size:11px;color:#777;margin-top:8px", r$note))
  }

  # The model actually fitted, its random structure, and its calibration status.
  # A reader who does not know which random structure was used is reading an
  # unnamed model, and the ladder can descend.
  output$harmonic_traj_header <- renderUI({
    if (!harmonic_is_mixed()) return(ts_header())
    ff <- harmonic_traj()
    if (!isTRUE(ff$ok)) return(div(class = "alert alert-warning", ff$message))
    cal <- dance_traj_calibration(ff, "kr", "full")
    div(style = "background:#f7f7f7;border:1px solid #e3e3e3;border-radius:3px;padding:8px 10px;margin-bottom:10px",
        tags$div(style = "font-family:monospace;font-size:11px;color:#555",
                 sprintf("%s + %s", ff$spec$fixed_formula, ff$re_formula)),
        tags$div(style = "font-size:12px;color:#555;margin-top:4px",
                 sprintf("%s — %d participants, %d curves, %d cells. Random structure: %s (rung %d of %d).%s",
                         ff$spec$design_kind, ff$spec$n_participants, ff$spec$n_curves,
                         nrow(ff$spec$cells), ff$re_label, ff$re_rung, ff$n_rungs,
                         if (isTRUE(ff$singular))
                           sprintf(" %d variance dimension(s) at the boundary, kept.", ff$boundary_dims)
                         else "")),
        # A COST WARNING, not a hidden threshold. Kenward-Roger refits the
        # reduced model for every block, so on a large fit the tests cost
        # several times the fit itself: measured on 1052 participants with a
        # linear trend, 30 s to fit and 129 s to test. Nobody should discover
        # that by waiting. The number is an estimate from the fit's own size,
        # stated as one, and the remedy is named -- the choice stays the
        # user's, and Kenward-Roger remains the default and the recommendation.
        if (identical(input$harmonic_df_method %||% "kr", "kr") &&
            (ff$spec$n_obs %||% 0) > 4000) {
          est <- max(10, round((ff$spec$n_obs / 16000) *
                               (1 + 6 * (!identical(ff$spec$trend, "none"))) * 25))
          tags$div(style = "font-size:11px;color:#8a5a12;margin-top:6px",
                   sprintf(paste("Large fit (%s observations). Kenward-Roger refits the model",
                                 "for each block, so the tests below may take on the order of",
                                 "%d s. Satterthwaite, under Advanced, is near-instant and is",
                                 "the right choice while exploring; switch back before you report."),
                           format(ff$spec$n_obs, big.mark = ","), est))
        },
        if (!isTRUE(cal$validated))
          tags$div(style = "font-size:11px;color:#8a5a12;margin-top:6px", cal$calibration)
        else tags$div(style = "font-size:11px;color:#8a5a12;margin-top:6px",
                      "Provisional calibration — see the validation notes."))
  })

  # EVERY BLOCK TEST TAB 6 NEEDS, COMPUTED ONCE.
  # ------------------------------------------------------------------------
  # Kenward-Roger refits the reduced model AND computes an adjusted covariance,
  # so one KR test costs about as much as the original fit: measured on 1052
  # participants with two harmonics, the fit took 21.2 s and each KR test 18.6 s.
  # The first version of this tab ran about five of them per render -- one per
  # factorial effect, three for the decomposition, and one more purely to read a
  # method name out of the result for a caption. That is why fitting "seemed to
  # take a long time": it was not the model, it was the same expensive test
  # being computed several times over.
  #
  # All of them now come from one reactive, so each distinct block is computed
  # once per fit and the panels read from the cache.
  harmonic_traj_tests <- reactive({
    ff <- harmonic_traj()
    if (!isTRUE(ff$ok)) return(NULL)
    dfm <- input$harmonic_df_method %||% "kr"
    dts <- ff$spec$design_terms
    tl <- attr(stats::terms(stats::as.formula(ff$spec$fixed_formula)), "term.labels")

    effects <- unlist(lapply(seq_along(dts), function(k)
      utils::combn(dts, k, paste, collapse = ":")), use.names = FALSE)
    want <- if (identical(ff$spec$trend, "none"))
      c("full", "circadian", "level") else
      c("full", "shape", "circadian", "trend", "level")

    # BRIEF 1-3. EVERY TEST HERE IS NOW A MODEL-BASED MARGINAL CONTRAST.
    # It used to select treatment-coded term blocks: "the Group effect" was every
    # model term mentioning Group, dropped and refitted. With a second factor in
    # the model that is the Group effect AT THE REFERENCE LEVEL of the other
    # factor, not a marginal main effect -- and since dance_traj_long() rebuilds
    # design factors from as.character(), the reference level is whichever one
    # sorts first alphabetically, which is not a scientific choice. Measured on a
    # 3 x 3 design with a real interaction, renaming the condition levels moved
    # the "Group effect" from F(6, 97.5) = 4.27 to 35.10, p 7e-4 to 3e-22. The
    # marginal test gives the same answer either way.
    #
    # A request is an (effect, block) pair, and two of them can still be the SAME
    # hypothesis -- in a one-factor design the full trajectory effect of Group and
    # the full block are identical -- so requests are deduplicated on the contrast
    # matrix itself rather than on a name. Kenward-Roger refits per test, so this
    # is worth doing.
    reqs <- c(lapply(effects, function(ef) list(effect = strsplit(ef, ":", fixed = TRUE)[[1]],
                                                block = "full", name = ef)),
              lapply(want, function(w) list(effect = character(0), block = w, name = w)))
    Ls <- lapply(reqs, function(r) tryCatch(
      dance_traj_marginal_L(ff, r$effect, r$block), error = function(e) NULL))
    keys <- vapply(Ls, function(L) if (is.null(L)) "" else
      paste(format(round(unclass(L), 10), scientific = FALSE), collapse = "|"), character(1))
    uniq <- unique(keys[nzchar(keys)])

    withProgress(message = sprintf("Testing %d coefficient block%s...",
                                   length(uniq), if (length(uniq) == 1L) "" else "s"),
                 value = 0, {
      computed <- lapply(uniq, function(k) {
        incProgress(1 / length(uniq))
        i <- match(k, keys)
        dance_traj_marginal_test(ff, reqs[[i]]$effect, reqs[[i]]$block, df_method = dfm)
      })
      names(computed) <- uniq
    })
    names(keys) <- vapply(reqs, function(r) r$name, character(1))
    get1 <- function(nm) if (nzchar(keys[[nm]])) computed[[keys[[nm]]]] else NULL

    list(effects = stats::setNames(lapply(effects, get1), effects),
         blocks  = stats::setNames(lapply(want, get1), want),
         n_computed = length(uniq), n_requested = length(keys),
         # THE METHOD IS PER TEST, NOT PER TABLE. dance_traj_block_test() asks
         # for Kenward-Roger and falls through to Satterthwaite when KRmodcomp
         # fails on that particular reduced model -- silently, and one block at
         # a time. This used to read the FIRST successful result and print its
         # method under the whole table, so a table could carry two different
         # approximations and say it carried one. That is how a denominator df
         # ends up looking impossible: the blocks are not comparable because
         # they were not all computed the same way.
         methods_all = unique(vapply(Filter(function(r) isTRUE(r$ok), computed),
                                     function(r) r$method %||% "F test", character(1))),
         method = (Filter(function(r) isTRUE(r$ok), computed)[[1]]$method) %||% "F test")
  })

  # The cell coefficients, cached per harmonic. dance_traj_cell_coefs() runs
  # emmeans twice internally, and the component table and every pairwise panel
  # each asked for them again -- cheap next to Kenward-Roger, but repeated on
  # every selector change, which is what makes a panel feel unresponsive.
  harmonic_traj_coefs <- reactive({
    ff <- harmonic_traj(); if (!isTRUE(ff$ok)) return(NULL)
    stats::setNames(lapply(seq_len(ff$spec$n_harmonics),
                           function(h) dance_traj_cell_coefs(ff, h)),
                    paste0("H", seq_len(ff$spec$n_harmonics)))
  })

  # ---- 1. PRIMARY: the full trajectory difference, per factorial effect ------
  output$harmonic_traj_primary <- renderUI({
    if (!harmonic_is_mixed()) return(ts_primary())
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    dts <- ff$spec$design_terms
    # each factorial effect gets its own full-trajectory test: the main effects
    # and, when there is more than one factor, their interaction
    tt <- harmonic_traj_tests(); req(!is.null(tt))
    rows <- lapply(names(tt$effects), function(ef) {
      r <- tt$effects[[ef]]
      if (is.null(r) || !isTRUE(r$ok)) return(NULL)
      tags$tr(
        tags$td(style = "padding:4px 12px 4px 0;font-weight:600", gsub(":", " × ", ef)),
        tags$td(style = "padding:4px 12px 4px 0;font-family:monospace",
                sprintf("F(%.0f, %.1f) = %.3f", r$df1, r$df2, r$statistic)),
        tags$td(style = "padding:4px 0;font-family:monospace",
                sprintf("p %s", dance_fmt_p(r$p))))
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(helpText("No testable design terms in this model."))
    tagList(tags$table(style = "margin:4px 0 6px 0", tags$tbody(rows)),
            tags$div(style = "font-size:11px;color:#777",
                     sprintf("%s. Each row tests level, trend and all harmonics jointly.",
                             tt$method)))
  })

  # ---- 2. SECONDARY: the component decomposition ----------------------------
  output$harmonic_traj_components <- renderUI({
    if (!harmonic_is_mixed()) return(ts_components())
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    # Without a trend, the SHAPE block and the CIRCADIAN block contain exactly
    # the same terms, and printing both prints one test twice under two names --
    # which reads as corroboration. Real data made this obvious: both rows came
    # back F(12, 1823.5) = 2.028. Only the blocks that differ are shown, and
    # which those are is decided in harmonic_traj_tests() so the two panels
    # cannot disagree -- or compute the same expensive test twice.
    tt <- harmonic_traj_tests(); req(!is.null(tt))
    blocks <- Filter(function(r) isTRUE(r$ok), tt$blocks)
    lab <- c(full = "Full trajectory", shape = "Temporal shape (trend + harmonics)",
             circadian = "Rhythmic block (all harmonics)",
             trend = "Non-periodic trend", level = "Baseline / reference level at t = 0")
    # Marked per row, and only when the table is not uniform: the denominator df
    # of an F test is a property of the approximation that produced it, so two
    # rows computed by different approximations cannot be read against each
    # other -- and a df that looks impossible next to its neighbour is usually
    # this, not an error in the fit.
    mixed <- length(tt$methods_all %||% character(0)) > 1L
    short_m <- function(m) {
      if (grepl("Kenward", m, fixed = TRUE)) return("KR")
      if (grepl("Satterth", m, fixed = TRUE)) return("Satt.")
      if (grepl("asymptotic", m, fixed = TRUE)) return("Wald")
      if (grepl("likelihood", m, fixed = TRUE)) return("LRT")
      m
    }
    rows <- lapply(names(blocks), function(b) {
      r <- blocks[[b]]
      tags$tr(
        tags$td(style = "padding:3px 12px 3px 0", lab[[b]] %||% b),
        # contrasts, not terms: the hypothesis is an L matrix now, and its row
        # count is what the numerator df actually is
        tags$td(style = "padding:3px 12px 3px 0;color:#777",
                sprintf("%d contrasts", r$n_contrasts %||% length(r$terms))),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                sprintf("F(%.0f, %.1f) = %.2f", r$df1, r$df2, r$statistic)),
        tags$td(style = "padding:3px 12px 3px 0;font-family:monospace", dance_fmt_p(r$p)),
        if (mixed) tags$td(style = "padding:3px 0;font-size:11px;color:#c0392b",
                           short_m(r$method %||% "")))
    })
    # derived per-harmonic quantities, from the SAME fit
    drows <- list()
    cf <- harmonic_traj_coefs()
    for (h in seq_len(ff$spec$n_harmonics)) {
      co <- cf[[paste0("H", h)]]
      if (is.null(co) || !isTRUE(co$ok)) next
      ap <- dance_traj_amp_phase_ci(co, method = "joint", n_draw = 6000)
      if (!isTRUE(ap$ok)) next
      for (i in seq_len(nrow(ap$table))) {
        r <- ap$table[i, ]
        drows[[length(drows) + 1L]] <- tags$tr(
          tags$td(style = "padding:3px 12px 3px 0", sprintf("H%d amplitude — %s", h, r$cell)),
          tags$td(style = "padding:3px 12px 3px 0;color:#777", "derived"),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  sprintf("%.3f", r$amplitude)),
          tags$td(style = "padding:3px 0;font-family:monospace",
                  sprintf("[%.3f, %.3f]", r$amplitude_lo, r$amplitude_hi)))
        drows[[length(drows) + 1L]] <- tags$tr(
          tags$td(style = "padding:3px 12px 3px 0", sprintf("H%d acrophase — %s", h, r$cell)),
          tags$td(style = "padding:3px 12px 3px 0;color:#777",
                  if (isTRUE(r$phase_defined)) "circular" else "circular (weak)"),
          # CLOCK TIME, like every other acrophase this app reports. The derived
          # column holds MODEL-ELAPSED hours on the harmonic's own effective
          # period; converting is what the legacy tables already do, and showing
          # one convention here and another there is how an acrophase gets
          # misread by exactly the offset between them.
          tags$td(style = sprintf("padding:3px 12px 3px 0;font-family:monospace;%s",
                                  if (isTRUE(r$phase_defined)) "" else "color:#aaa"),
                  dance_clock_label(
                    # $hours, NOT $first. dance_acrophase_clock() returns
                    # hours / all_hours / elapsed / effective_period / harmonic /
                    # clock_origin, and its COMMENT on the first element reads
                    # "first maximum on the clock" -- which is where $first came
                    # from. R returns NULL for a name that is not there, so this
                    # cell has been rendering EMPTY, silently, since it was
                    # written: no error, no NA, just a blank column where the
                    # acrophase should be.
                    dance_acrophase_clock(hours = r$acrophase_time,
                                          period = ff$spec$period, harmonic = h,
                                          clock_origin = ff$clock_origin %||% 0)$hours,
                    ff$spec$period, show_day = FALSE)),
          # THE INTERVAL, in the same clock as the estimate. An arc width alone
          # says how PRECISE the phase is but not WHERE it is, so it could not be
          # quoted or compared without doing the arithmetic by hand. Both now:
          # the endpoints to report, the width to compare across cells. An
          # interval crossing the period boundary is marked, because
          # "[22:40, 01:15]" read left to right looks like an empty range.
          tags$td(style = "padding:3px 0;font-family:monospace",
                  if (isTRUE(r$phase_defined)) {
                    to_clock <- function(x) dance_clock_label(
                      dance_acrophase_clock(hours = x, period = ff$spec$period,
                                            harmonic = h,
                                            clock_origin = ff$clock_origin %||% 0)$hours,
                      ff$spec$period, show_day = FALSE)
                    sprintf("[%s, %s]%s  arc %.2f h", to_clock(r$acrophase_lo),
                            to_clock(r$acrophase_hi),
                            if (isTRUE(r$acrophase_wraps)) " (wraps)" else "",
                            r$acrophase_arc_time)
                  } else "no interval: the (cos, sin) region includes the origin"))
      }
    }
    # WHERE THE DRAWN CURVE PEAKS, per cell. Asked for directly: the H1 acrophase
    # did not match the peaks visible in the trajectory plot, and it is not
    # supposed to -- H1 is one harmonic, the plotted curve is level + trend + H1 +
    # H2. Showing both, next to each other, is cheaper than explaining the
    # difference every time.
    pk <- tryCatch(dance_traj_curve_peaks(ff, "full"), error = function(e) NULL)
    if (!is.null(pk) && isTRUE(pk$ok)) {
      for (i in seq_len(nrow(pk$table))) {
        pr <- pk$table[i, ]
        drows[[length(drows) + 1L]] <- tags$tr(
          tags$td(style = "padding:3px 12px 3px 0",
                  sprintf("Fitted curve peak — %s", pr$cell)),
          tags$td(style = "padding:3px 12px 3px 0;color:#777", "plotted"),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  dance_clock_label(pr$peak_clock, ff$spec$period, show_day = FALSE)),
          tags$td(style = sprintf("padding:3px 0;font-family:monospace;color:%s",
                                  if (isTRUE(pr$peak_interior)) "#777" else "#c0392b"),
                  if (isTRUE(pr$peak_interior)) sprintf("value %.2f", pr$peak_fit)
                  else sprintf("value %.2f — still rising here, so this is the end of the search window and NOT a peak",
                               pr$peak_fit)))
      }
    }
    tagList(
      tags$table(style = "margin:4px 0", tags$tbody(rows, drows)),
      tags$div(style = "font-size:11px;color:#777;margin-top:6px",
        HTML(paste(
          "Joint rows are block tests on the fitted model; derived rows come from the",
          "<b>same</b> coefficient vector and covariance, with the acrophase interval",
          "reported as an <b>arc</b> and suppressed where the amplitude interval covers",
          "zero. Acrophases are <b>clock times</b>; the arc beside each is its width, in",
          "hours, on that harmonic's own effective period. The baseline row is the fitted",
          "value at the reference time, which is <b>not</b> a MESOR when a trend is",
          "present.",
          "<br><b>An H1 acrophase is not the peak of the drawn curve.</b> H1 is the first",
          "harmonic alone; the plotted trajectory is level + trend + every harmonic, and",
          "H2 has its own maximum that shifts where the sum peaks. The <i>fitted curve",
          "peak</i> rows give the time a reader takes off the plot. A third number again",
          "is the two-stage acrophase in the polar plots: a per-participant",
          "single-harmonic cosinor has to absorb H2 and the trend into one cosine, which",
          "displaces it further. Compare like with like, and judge any gap against the",
          "<b>arc</b> beside the acrophase, not against the eye."))),
      # emmeans turns its df adjustment off above 3000 observations and says so
      # on the CONSOLE, which nobody running a Shiny app is reading. The choice
      # it makes is a statistical one, so it belongs on the screen.
      local({
        nt <- Filter(Negate(is.null), lapply(cf, function(co) co$df_note))
        if (length(nt)) tags$div(style = "font-size:11px;color:#777;margin-top:4px",
                                 HTML(nt[[1]]))
      }),
      if (mixed) tags$div(style = "font-size:11px;color:#c0392b;margin-top:4px",
        HTML(paste0(
          "<b>Not every row used the same approximation</b> (", 
          paste(vapply(tt$methods_all, short_m, character(1)), collapse = ", "),
          "). Kenward-Roger was attempted and failed on at least one block, which falls",
          " back to Satterthwaite for that block alone. The denominator df belongs to the",
          " approximation that produced it, so rows marked differently are not directly",
          " comparable with each other. Select Satterthwaite under Advanced to compute",
          " every row the same way."))))
  })

  # ---- 3. THE FOUR COMPONENT VIEWS, per design cell -------------------------
  # With no trend in the model, "nonperiodic change from origin" is identically
  # zero and "baseline + harmonics" is identical to the full trajectory. Found by
  # running on real data: the selector offered all four regardless, so two of
  # them drew a flat zero line and a duplicate. The options now follow the model.
  harmonic_traj_component <- reactive({
    cc <- input$harmonic_traj_component %||% "full"
    if (!harmonic_is_mixed()) return(if (cc %in% c("full", "rhythm")) cc else "full")
    ff <- harmonic_traj()
    # the fit decides the vocabulary, so a stale selection from a previous model
    # -- H2 after refitting with one harmonic, a trend view with no trend --
    # falls back rather than asking for a component that does not exist
    if (!isTRUE(ff$ok)) return(cc)
    if (cc %in% dance_traj_components(ff)) cc else "full"
  })
  observe({
    if (!harmonic_is_mixed()) {
      mod <- values$harmonic_model; if (is.null(mod)) return()
      ch <- c("Full fitted trajectory" = "full")
      if (!identical(mod$trend_type %||% "none", "none")) ch <- c(ch, "Rhythm only (trend removed)" = "rhythm")
      sel <- isolate(input$harmonic_traj_component) %||% "full"
      updateSelectInput(session, "harmonic_traj_component", choices = ch,
                        selected = if (sel %in% ch) sel else "full")
      return()
    }
    ff <- harmonic_traj()
    if (!isTRUE(ff$ok)) return()
    comps <- dance_traj_components(ff)
    ch <- stats::setNames(comps, vapply(comps, dance_traj_component_label,
                                        character(1), period = ff$spec$period))
    sel <- isolate(input$harmonic_traj_component) %||% "full"
    updateSelectInput(session, "harmonic_traj_component", choices = ch,
                      selected = if (sel %in% comps) sel else "full")
  })

  output$harmonic_traj_curves <- renderPlotly({
    if (!harmonic_is_mixed()) return(ts_curves_plot())
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    pr <- dance_traj_predict(ff, conf = 0.95,
                             band = input$harmonic_traj_band %||% "pointwise",
                             component = harmonic_traj_component(),
                             n_time = 160)
    req(isTRUE(pr$ok))
    cols <- dance_group_colors(pr$cells)
    p <- plot_ly()
    for (cl in pr$cells) {
      d <- pr$table[pr$table$cell == cl, , drop = FALSE]
      # legendgroup, or clicking a name in the legend hides the LINE and leaves
      # its band behind -- reported from a real session, and it makes the plot
      # worse than useless: a shaded band with no curve reads as a different
      # group's uncertainty. Plotly ties visibility to the group, not the name.
      p <- p %>%
        add_ribbons(x = d$t, ymin = d$lo, ymax = d$hi, name = cl, legendgroup = cl,
                    line = list(width = 0), fillcolor = dance_group_rgba(cols[[cl]], 0.16),
                    showlegend = FALSE, hoverinfo = "skip") %>%
        add_lines(x = d$t, y = d$fit, name = cl, legendgroup = cl,
                  line = list(color = cols[[cl]], width = 2.4))
    }
    if (isTRUE(input$harmonic_traj_raw) &&
        identical(harmonic_traj_component(), "full")) {
      # raw points belong only under the FULL trajectory: the other three views
      # are components of the fit, and no observation corresponds to one
      dd <- ff$spec$data
      key <- if (length(ff$spec$design_terms))
        do.call(paste, c(lapply(ff$spec$design_terms, function(f) as.character(dd[[f]])),
                         list(sep = " x "))) else rep("(all)", nrow(dd))
      for (cl in pr$cells) {
        i <- key == cl
        if (!any(i)) next
        p <- p %>% add_markers(x = dd$t[i], y = dd$y[i], name = cl, legendgroup = cl,
                               showlegend = FALSE,
                               marker = list(color = dance_group_rgba(cols[[cl]], 0.35), size = 4))
      }
    }
    cmp <- harmonic_traj_component()
    kh <- dance_traj_component_harmonic(cmp)
    ylab <- if (!is.na(kh))
      sprintf("H%d component (centred on 0)", kh)
    else switch(cmp,
      full = {
        md <- values$harmonic_model
        if (!is.null(md$dv_name))
          paste0(md$dv_name, if (!is.null(md$dv_units)) paste0(" (", md$dv_units, ")") else "")
        else "Response"
      },
      harmonics = "Periodic component (centred on 0)",
      baseline_harm = "Baseline + harmonics",
      trend = "Change from the reference time")
    # CLOCK TIME ON THE AXIS, not elapsed. The model is fitted on an axis
    # anchored at the first observation, but nobody reads a sleepiness
    # trajectory in "hours since the study started", and every other plot in
    # this module already labels the clock. The tick VALUES stay on the model
    # axis so the traces do not move; only the labels are converted.
    co <- ff$clock_origin %||% 0
    P <- ff$spec$period %||% 24
    brk <- pretty(range(pr$times), n = 8)
    brk <- brk[brk >= min(pr$times) & brk <= max(pr$times)]
    p %>% layout(
      xaxis = list(title = if (co != 0) "Clock time" else
                     sprintf("Time (%s)", ff$spec$time_units %||% "h"),
                   tickmode = "array", tickvals = brk,
                   ticktext = dance_clock_label(brk + co, P, show_day = TRUE,
                                                with_minutes = FALSE)),
      yaxis = list(title = ylab),
      hovermode = "x unified",
      legend = list(orientation = "h", y = -0.18))
  })

  output$harmonic_traj_band_note <- renderUI({
    if (!harmonic_is_mixed()) {
      gc <- ts_group_curves(harmonic_traj_component())
      return(if (is.null(gc)) NULL else tags$div(style = "font-size:11px;color:#777;margin-top:4px",
        paste(gc$note, "Simultaneous (Scheffe) bands are not available under two-stage.")))
    }
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    pr <- dance_traj_predict(ff, band = input$harmonic_traj_band %||% "pointwise",
                             component = harmonic_traj_component(), n_time = 2)
    req(isTRUE(pr$ok))
    tags$div(style = "font-size:11px;color:#777;margin-top:4px", pr$note)
  })

  # ---- 4. DIFFERENCE CURVE: one contrast, not two curves subtracted ---------
  output$harmonic_traj_diff_a_ui <- renderUI({
    if (!harmonic_is_mixed()) {
      ts <- harmonic_ts(); req(!is.null(ts), isTRUE(ts$ok))
      return(selectInput("harmonic_traj_diff_a", "Group A:", choices = ts$groups, selected = ts$groups[1]))
    }
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    cl <- dance_traj_cell_grid(ff$spec)$.cell
    selectInput("harmonic_traj_diff_a", "Cell A:", choices = cl, selected = cl[1])
  })
  output$harmonic_traj_diff_b_ui <- renderUI({
    if (!harmonic_is_mixed()) {
      ts <- harmonic_ts(); req(!is.null(ts), isTRUE(ts$ok))
      return(selectInput("harmonic_traj_diff_b", "Group B:", choices = ts$groups,
                         selected = if (length(ts$groups) > 1) ts$groups[2] else ts$groups[1]))
    }
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    cl <- dance_traj_cell_grid(ff$spec)$.cell
    selectInput("harmonic_traj_diff_b", "Cell B:", choices = cl,
                selected = if (length(cl) > 1) cl[2] else cl[1])
  })

  harmonic_traj_diff <- reactive({
    if (!harmonic_is_mixed()) return(ts_diff())
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    a <- input$harmonic_traj_diff_a; b <- input$harmonic_traj_diff_b
    req(a, b)
    dance_traj_diff_curve(ff, a, b, n_time = 160,
                          band = if (isTRUE(input$harmonic_traj_diff_sim))
                            "simultaneous" else "pointwise")
  })

  output$harmonic_traj_diff_plot <- renderPlotly({
    dc <- harmonic_traj_diff()
    if (!isTRUE(dc$ok)) return(plotly_empty() %>%
                                 layout(title = list(text = dc$message, font = list(size = 12))))
    d <- dc$table
    plot_ly() %>%
      # legendgroup: the band has no legend entry of its own, so without it
      # hiding the difference curve leaves its shading behind
      add_ribbons(x = d$t, ymin = d$lo, ymax = d$hi, line = list(width = 0),
                  legendgroup = "diff",
                  fillcolor = "rgba(31,107,74,0.16)", showlegend = FALSE, hoverinfo = "skip") %>%
      add_lines(x = d$t, y = d$diff, line = list(color = "#1F6B4A", width = 2.4),
                legendgroup = "diff",
                name = sprintf("%s − %s", dc$cell1, dc$cell2)) %>%
      add_lines(x = range(d$t), y = c(0, 0),
                line = list(color = "#999", width = 1, dash = "dot"),
                showlegend = FALSE, hoverinfo = "skip") %>%
      layout(xaxis = local({
               ff <- harmonic_traj(); co <- ff$clock_origin %||% 0
               brk <- pretty(range(d$t), n = 8)
               brk <- brk[brk >= min(d$t) & brk <= max(d$t)]
               list(title = if (co != 0) "Clock time" else "Time",
                    tickmode = "array", tickvals = brk,
                    ticktext = dance_clock_label(brk + co, ff$spec$period %||% 24,
                                                 show_day = TRUE, with_minutes = FALSE))
             }),
             yaxis = list(title = sprintf("%s − %s", dc$cell1, dc$cell2)),
             hovermode = "x unified", showlegend = FALSE)
  })

  output$harmonic_traj_diff_note <- renderUI({
    dc <- harmonic_traj_diff()
    if (!isTRUE(dc$ok)) return(NULL)
    d <- dc$table
    sep <- if (any(d$excludes_zero)) {
      r <- range(d$t[d$excludes_zero])
      sprintf("The band excludes zero between %.1f and %.1f h.", r[1], r[2])
    } else "The band covers zero throughout."
    tags$div(style = "font-size:11px;color:#777;margin-top:4px",
             paste(sep, dc$note))
  })

  # ---- 5. PAIRWISE POST-HOC, over an arbitrary factorial --------------------
  output$harmonic_traj_effect_ui <- renderUI({
    if (!harmonic_is_mixed())
      return(selectInput("harmonic_traj_effect", "Effect to inspect:", choices = c("All groups" = "_all_")))
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    dts <- ff$spec$design_terms
    ch <- c("All design cells" = "_all_")
    if (length(dts) > 1) {
      # simple effects: one factor, inside each level of the others
      for (f in dts) ch[sprintf("%s, within each other level", f)] <- f
    }
    selectInput("harmonic_traj_effect", "Effect to inspect:", choices = ch)
  })

  # The pairwise vocabulary follows the fit, like the component selector: the
  # omnibus blocks a model can answer depend on whether it has a trend and more
  # than one harmonic, and a duplicate pair of blocks is one entry, not two.
  observe({
    if (!harmonic_is_mixed()) {
      mod <- values$harmonic_model; if (is.null(mod)) return()
      ws <- dance_ts_pair_whats(mod)
      ch <- stats::setNames(ws, vapply(ws, dance_ts_pair_label, character(1), period = mod$period %||% 24))
      sel <- isolate(input$harmonic_traj_what) %||% "level"
      updateSelectInput(session, "harmonic_traj_what", choices = ch,
                        selected = if (sel %in% ws) sel else "level")
      return()
    }
    ff <- harmonic_traj()
    if (!isTRUE(ff$ok)) return()
    ws <- dance_traj_pair_whats(ff)
    ch <- stats::setNames(ws, vapply(ws, dance_traj_pair_label, character(1),
                                     period = ff$spec$period,
                                     n_harmonics = ff$spec$n_harmonics))
    # "amplitude" from an earlier model means H1 here, so a stale selection maps
    # forward instead of silently resetting to Level
    sel <- isolate(input$harmonic_traj_what) %||% "level"
    if (identical(sel, "amplitude")) sel <- "amplitude1"
    if (identical(sel, "phase")) sel <- "phase1"
    updateSelectInput(session, "harmonic_traj_what", choices = ch,
                      selected = if (sel %in% ws) sel else "level")
  })

  output$harmonic_traj_pairwise <- renderUI({
    if (!harmonic_is_mixed()) return(ts_pairwise())
    ff <- harmonic_traj(); req(isTRUE(ff$ok))
    what <- input$harmonic_traj_what %||% "level"
    adj  <- input$harmonic_traj_adjust %||% "holm"
    eff  <- input$harmonic_traj_effect %||% "_all_"
    dfm  <- input$harmonic_df_method %||% "kr"
    dts  <- ff$spec$design_terms
    if (identical(what, "amplitude")) what <- "amplitude1"
    if (identical(what, "phase")) what <- "phase1"
    if (!what %in% dance_traj_pair_whats(ff)) what <- "level"

    res <- if (identical(eff, "_all_") || !(eff %in% dts)) {
      list(list(title = NULL, r = dance_traj_contrasts(ff, what, adjust = adj,
                                                       df_method = dfm)))
    } else {
      # every combination of the OTHER factors is one slice
      others <- setdiff(dts, eff)
      lv <- lapply(others, function(f) levels(ff$spec$data[[f]]))
      names(lv) <- others
      grid <- expand.grid(lv, stringsAsFactors = FALSE)
      lapply(seq_len(nrow(grid)), function(i) {
        at <- as.list(grid[i, , drop = FALSE])
        list(title = paste(sprintf("%s = %s", names(at), unlist(at)), collapse = ", "),
             r = dance_traj_simple_effects(ff, eff, at = at, what = what, adjust = adj,
                                           df_method = dfm))
      })
    }

    blocks <- lapply(res, function(x) {
      r <- x$r
      if (!isTRUE(r$ok)) return(tags$div(style = "color:#8a5a12;font-size:12px", r$message))
      tb <- r$table
      # A JOINT row compares the cells on several components at once, so it has an
      # F and a p and NO single estimate. Printing a blank Estimate column beside
      # it would read as a missing number rather than an absent quantity, so the
      # header changes with the comparison.
      joint <- isTRUE(r$joint)
      hdr <- if (joint) c("Pair", "Test", "", "p", "p adj")
             else c("Pair", "Estimate", "95% CI", "p", "p adj")
      rows <- lapply(seq_len(nrow(tb)), function(i) {
        z <- tb[i, ]
        tags$tr(
          tags$td(style = "padding:3px 12px 3px 0", sprintf("%s vs %s", z$cell1, z$cell2)),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  if (joint) sprintf("F(%g, %.1f) = %.2f", z$df1, z$df2, z$statistic)
                  else sprintf("%+.3f", z$estimate)),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  if (joint) "" else if (isTRUE(z$defined)) sprintf("[%.3f, %.3f]", z$lo, z$hi) else "—"),
          tags$td(style = "padding:3px 12px 3px 0;font-family:monospace",
                  if (isTRUE(z$defined)) dance_fmt_p(z$p_raw) else "—"),
          tags$td(style = "padding:3px 0;font-family:monospace;font-weight:600",
                  if (isTRUE(z$defined)) dance_fmt_p(z$p_adj) else "undefined"))
      })
      tagList(
        if (!is.null(x$title)) tags$div(style = "font-weight:600;margin-top:10px;font-size:13px", x$title),
        tags$table(style = "margin:4px 0",
          tags$thead(tags$tr(lapply(hdr, function(h)
            tags$th(style = "text-align:left;padding:2px 12px 2px 0;font-size:11px;color:#777;font-weight:500", h)))),
          tags$tbody(rows)),
        tags$div(style = "font-size:11px;color:#777", r$unit))
    })
    tagList(blocks, tags$div(style = "font-size:11px;color:#777;margin-top:8px",
                             res[[1]]$r$note %||% ""))
  })

  
