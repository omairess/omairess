# ==========================================================================
# server/72_harmonic.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 2879-7190  (cosinor core, harmonic regression + outputs)
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
                              # MERGED APP: reuse the clock times detected once
                              # at import (values$time_numeric) instead of
                              # re-detecting them here.  Additive: the default
                              # is still "_index_".
                              "Use shared times detected at import" = "_shared_",
                              numeric_vars),
                  selected = "_index_"),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_shared_'",
        if(!is.null(values$time_numeric) && length(values$time_numeric) == n_time) {
          div(style = "color: green; font-size: 0.9em;", icon("check-circle"),
              sprintf(" Using the %d time values detected at import: %s%s",
                      n_time, paste(head(values$time_numeric, 6), collapse = ", "),
                      if(n_time > 6) ", ..." else ""))
        } else {
          div(style = "color: #b00; font-size: 0.9em;", icon("exclamation-triangle"),
              " Import did not detect usable time values from the column names. Use 'Specify times manually'.")
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
       <strong>Suggested MESOR bounds:</strong> [%.2f, %.2f] (mean ± range)<br>
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
    
    acrophases_time <- acrophases * period / (2 * pi)
    acro_se_time <- acro_se * period / (2 * pi)
    
    # Goodness of fit - calculate manually because lm(y ~ X - 1) R² is vs origin, not mean
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum(residuals(fit)^2)
    r_squared <- 1 - ss_resid / ss_total

    # Adjusted R² accounting for number of predictors
    n_predictors <- ncol(X)  # MESOR + trend + harmonics
    adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - n_predictors)

    percent_rhythm <- r_squared * 100

    df1 <- 2 * n_harmonics
    df2 <- n - 2 * n_harmonics - 1 - n_trend_params
    f_stat <- ((ss_total - ss_resid) / df1) / (ss_resid / df2)
    p_value <- pf(f_stat, df1, df2, lower.tail = FALSE)

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

    # Calculate proportions of total R²
    if(!is.na(r_squared) && r_squared > 0) {
      percent_S <- (r_squared_S / r_squared) * 100
      percent_C <- (r_squared_C / r_squared) * 100
    }
    
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
      r_squared_S = r_squared_S,        # R² from Process S alone
      r_squared_C = r_squared_C,        # R² from Process C alone
      percent_S = percent_S,            # % of total R² from S
      percent_C = percent_C,            # % of total R² from C
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
                                   tau_min = 0.5, tau_max = NA) {
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
        fit_success <- TRUE
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
        fit_success <- TRUE
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
        fit_success <- TRUE
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
        fit_success <- TRUE
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

        if(r2_check > 0) {
          fit_success <- TRUE
          break
        }
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("retry", tau_mult, ":", e$message))
      })
    }
  }

  if(!fit_success || is.null(nls_fit)) {
    return(list(
      success = FALSE,
      message = sprintf("Nonlinear fit failed to converge. Errors: %s",
                       paste(head(error_msgs, 3), collapse = "; "))
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

    for(h in 1:n_harmonics) {
      beta_cos <- coefs[paste0("b_cos", h)]
      beta_sin <- coefs[paste0("b_sin", h)]
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi

      # Approximate SE
      se_cos <- if(!is.na(se[paste0("b_cos", h)])) se[paste0("b_cos", h)] else 0
      se_sin <- if(!is.na(se[paste0("b_sin", h)])) se[paste0("b_sin", h)] else 0
      amp_se[h] <- sqrt(se_cos^2 + se_sin^2) / sqrt(2)
      acro_se[h] <- NA  # Complex for nonlinear
    }

    acrophases_time <- acrophases * period / (2 * pi)
    acro_se_time <- acro_se * period / (2 * pi)

    # Goodness of fit
    fitted_vals <- predict(nls_fit)
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum((y - fitted_vals)^2)
    r_squared <- 1 - ss_resid / ss_total
    percent_rhythm <- max(0, r_squared * 100)

    # Calculate p-value for circadian rhythm
    n_params <- length(coefs)
    df1 <- 2 * n_harmonics
    df2 <- n - n_params
    if(df2 > 0 && ss_resid > 0) {
      f_stat <- ((ss_total - ss_resid) / df1) / (ss_resid / df2)
      p_value <- pf(f_stat, df1, df2, lower.tail = FALSE)
    } else {
      f_stat <- NA
      p_value <- NA
    }

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
    loocv_rmse <- tryCatch({
      # Get leverage values (hat matrix diagonal)
      # For NLS, we approximate using the Jacobian
      if(!is.null(nls_fit)) {
        # Simple LOOCV using prediction errors
        # For nonlinear models, full refit for each point is too expensive
        # Use prediction residuals as approximation
        residuals_vec <- y - fitted_vals
        sqrt(mean(residuals_vec^2))
      } else {
        NA
      }
    }, error = function(e) NA)

    # ===========================================================================
    # Variance decomposition: Calculate R² for Process S and Process C separately
    # ===========================================================================
    r_squared_S <- 0  # Variance explained by homeostatic trend alone
    r_squared_C <- 0  # Variance explained by circadian rhythm alone
    percent_S <- 0    # Percentage of total R² from Process S
    percent_C <- 0    # Percentage of total R² from Process C

    # Model with trend only (Process S) - using exp_sat trend
    if(trend_type == "exp_sat" && !is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
      A_sat_val <- trend_params$A_sat$coef
      tau_val <- trend_params$tau$coef
      fitted_trend <- mesor + A_sat_val * (1 - exp(-t_shifted / tau_val))
      ss_resid_trend <- sum((y - fitted_trend)^2)
      r_squared_S <- max(0, 1 - ss_resid_trend / ss_total)
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

    # Calculate proportions of total R²
    if(!is.na(r_squared) && r_squared > 0) {
      percent_S <- (r_squared_S / r_squared) * 100
      percent_C <- (r_squared_C / r_squared) * 100
    }

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
      loocv_rmse = loocv_rmse,          # Leave-one-out CV RMSE
      r_squared_S = r_squared_S,        # R² from Process S alone
      r_squared_C = r_squared_C,        # R² from Process C alone
      percent_S = percent_S,            # % of total R² from S
      percent_C = percent_C,            # % of total R² from C
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
      using_smoothed <- !is.null(values$smooth_data)
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
        # MERGED APP: the shared import step already extracted numeric clock
        # times from the column names into values$time_numeric.  Use them
        # verbatim so the harmonic tab and the fPCA/fANOVA/clustering tabs
        # all place the same column at the same time.
        if(is.null(values$time_numeric) || length(values$time_numeric) != n_time) {
          showNotification(
            "No shared time values available (import did not detect times from the column names). Use 'Specify times manually'.",
            type = "error", duration = 10)
          return()
        }
        time_vec <- as.numeric(values$time_numeric)
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
      param_cols <- c("subject", "mesor", "mesor_se")

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
                      "r_squared_S", "r_squared_C", "percent_S", "percent_C")
      
      individual_params <- data.frame(matrix(ncol = length(param_cols), nrow = 0))
      colnames(individual_params) <- param_cols
      
      # Coefficient offset for trend
      coef_offset <- 1 + n_trend_params
      
      # Track failed fits
      failed_fits <- list()
      
      # Store time offsets for prediction
      t_offset_global <- min(time_vec)
      t_center_global <- mean(time_vec)

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
        bounds_msg <- sprintf("Using parameter bounds: MESOR [%.2f, %.2f], Amplitude [%.2f, %.2f]",
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

          fit_i <- fit_cosinor(time_vec, y_i, period = period, n_harmonics = n_harmonics,
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
            row_data$percent_S <- fit_i$percent_S
            row_data$percent_C <- fit_i$percent_C

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
        mean_acrophase_time <- mean_acrophase_rad * period / (2 * pi)
        
        # Rayleigh test for uniformity of acrophases (first harmonic)
        n_valid <- sum(!is.na(individual_params$acrophase_rad_1))
        r_bar <- mean_amplitude / mean(individual_params$amplitude_1, na.rm = TRUE)
        rayleigh_z <- n_valid * r_bar^2
        rayleigh_p <- exp(-rayleigh_z)  # Approximation
        
        # SE from circular statistics
        circ_var <- 1 - r_bar
        circ_sd <- if(r_bar > 0 && r_bar < 1) sqrt(-2 * log(r_bar)) else NA
        
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
          mean_acrophases_time[h] <- acro_h * period / (2 * pi) / h  # Adjust for harmonic number
        }
        
        # Also compute arithmetic means of individual parameters
        indiv_means <- list(
          mesor = mean(individual_params$mesor, na.rm = TRUE),
          mesor_sd = sd(individual_params$mesor, na.rm = TRUE)
        )
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
        
        pop_mean_fit <- list(
          mean_mesor = mean_mesor,
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
        groups <- unique(group_var)
        group_fits <- list()
        
        for(g in groups) {
          idx <- which(group_var == g)
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
              grp_acrophases_time[h] <- acro_h * period / (2 * pi) / h
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

            group_fits[[as.character(g)]] <- list(
              group = g,
              n = nrow(grp_params),
              mean_mesor = mean(grp_params$mesor, na.rm = TRUE),
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
          acrophase_ci = quantile(boot_acrophase * period / (2 * pi), c(0.025, 0.975)),
          boot_mesor = boot_mesor,
          boot_amplitude = boot_amplitude,
          boot_acrophase = boot_acrophase,
          B = B
        )
      }
      
      # Store results
      values$harmonic_model <- list(
        individual_fits = individual_fits,
        individual_params = individual_params,
        pop_mean_fit = pop_mean_fit,
        group_fits = group_fits,
        boot_results = boot_results,
        time_vec = time_vec,
        original_times = original_times,
        wrap_applied = wrap_applied,
        period = period,
        n_harmonics = n_harmonics,
        trend_type = trend_type,
        include_trend = trend_type != "none",  # For backwards compatibility
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
  
  # Summary output
  output$harmonic_summary <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    cat("=== Harmonic Regression (Cosinor Analysis) Results ===\n\n")
    cat("Period:", mod$period, "\n")
    cat("Number of harmonics:", mod$n_harmonics, "\n")
    
    # Data source info
    if(isTRUE(mod$using_smoothed)) {
      cat("Data: SMOOTHED (missing values interpolated by FDA)\n")
    } else {
      cat("Data: RAW (no smoothing applied)\n")
      if(!is.null(mod$subjects_with_nas) && mod$subjects_with_nas > 0) {
        cat("  ⚠ ", mod$subjects_with_nas, " subjects have missing values - consider smoothing first!\n")
      }
    }
    
    # Show trend model type and sleep inertia
    trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
    include_inertia <- !is.null(mod$include_inertia) && isTRUE(mod$include_inertia)

    if(trend_type != "none") {
      trend_label <- switch(trend_type,
                            "linear" = "LINEAR (β·t)",
                            "log" = "LOGARITHMIC (β·log(t+1))",
                            "exp_sat" = "SATURATING EXPONENTIAL (A·(1-e^(-t/τ)))",
                            "Unknown")
      cat("Homeostatic trend:", trend_label, "\n")
      cat("Two-process model: Separates Process S (trend) from Process C (circadian)\n")
      cat("  Model equation: Y(t) = M + S(t) + C(t)\n")
    } else {
      cat("Homeostatic trend: None (circadian only)\n")
      cat("  Model equation: Y(t) = M + C(t)\n")
    }
    
    # Count successful and failed fits
    n_total <- length(mod$individual_fits)
    n_success <- sum(sapply(mod$individual_fits, function(f) isTRUE(f$success)))
    n_failed <- n_total - n_success
    
    cat("Number of subjects:", n_total, "\n")
    cat("  - Successfully fitted:", n_success, "\n")
    if(n_failed > 0) {
      cat("  - Failed fits:", n_failed, "(insufficient data points)\n")
      
      # Show which subjects failed
      failed_info <- sapply(seq_along(mod$individual_fits), function(i) {
        f <- mod$individual_fits[[i]]
        if(!isTRUE(f$success)) {
          if(!is.null(f$n_valid)) {
            sprintf("S%d(%d/%d pts)", i, f$n_valid, f$n_required)
          } else {
            sprintf("S%d", i)
          }
        } else {
          NULL
        }
      })
      failed_info <- failed_info[!sapply(failed_info, is.null)]
      if(length(failed_info) <= 10) {
        cat("    Failed subjects:", paste(failed_info, collapse = ", "), "\n")
      } else {
        cat("    Failed subjects:", paste(head(failed_info, 10), collapse = ", "), "...\n")
      }
    }
    
    cat("Number of time points:", length(mod$time_vec), "\n")
    
    # Show time points used
    if(length(mod$time_vec) <= 24) {
      if(isTRUE(mod$wrap_applied) && !is.null(mod$original_times)) {
        cat("Original times (clock):", paste(round(mod$original_times, 1), collapse=", "), "\n")
        cat("Adjusted times (linear):", paste(round(mod$time_vec, 1), collapse=", "), "\n")
        cat("(Times after midnight were adjusted for chronological order)\n")
      } else {
        cat("Time points used:", paste(round(mod$time_vec, 2), collapse=", "), "\n")
      }
    } else {
      cat("Time points: ", round(min(mod$time_vec), 2), "to", round(max(mod$time_vec), 2), 
          "(", length(mod$time_vec), "points)\n")
    }
    
    # Check if spacing is equal
    diffs <- diff(mod$time_vec)
    if(length(unique(round(diffs, 2))) > 1) {
      cat("Spacing: UNEQUAL (", paste(unique(round(diffs, 2)), collapse=", "), ")\n")
    } else {
      cat("Spacing: Equal (", round(diffs[1], 2), ")\n")
    }
    cat("\n")
    
    # Population mean statistics (always calculated)
    if(!is.null(mod$pop_mean_fit)) {
      cat("--- Population Mean Parameters (Vector-Averaged) ---\n")
      pop <- mod$pop_mean_fit
      cat(sprintf("MESOR:     %.3f\n", pop$mean_mesor))
      cat(sprintf("Amplitude (H1): %.3f\n", pop$mean_amplitude))
      cat(sprintf("Acrophase (H1): %.2f (%.2f hours)\n", 
                  pop$mean_acrophase_rad * 180 / pi, pop$mean_acrophase_time))
      
      # Show all harmonics if more than 1
      if(mod$n_harmonics > 1) {
        cat("\n  All Harmonics (vector-averaged):\n")
        for(h in 1:mod$n_harmonics) {
          cat(sprintf("    H%d: Amplitude=%.3f, Acrophase=%.2f hours\n", 
                      h, pop$mean_amplitudes[h], pop$mean_acrophases_time[h]))
        }
      }
      
      cat(sprintf("\nRayleigh test for uniformity (H1):\n"))
      cat(sprintf("  Z = %.3f, p = %.4f\n", pop$rayleigh_z, pop$rayleigh_p))
      cat(sprintf("  Mean resultant length (r̄) = %.3f\n", pop$r_bar))

      # Build and display symbolic equation
      {
        cat("\n--- Model Equation (symbolic) ---\n")
        sym <- "Y(t) = M"
        if(trend_type == "linear")  sym <- paste0(sym, " + \u03b2\u00b7t")
        else if(trend_type == "log")     sym <- paste0(sym, " + \u03b2\u00b7log(t+1)")
        else if(trend_type == "exp_sat") sym <- paste0(sym, " + A_sat\u00b7(1 \u2212 e^(\u2212t/\u03c4))")
        if(mod$n_harmonics == 1) {
          sym <- paste0(sym, sprintf(" + A\u00b7cos(2\u03c0\u00b7t/%.0f \u2212 \u03c6)", mod$period))
        } else {
          for(h_s in 1:mod$n_harmonics) {
            sub_h <- intToUtf8(0x2080 + h_s)
            sym <- paste0(sym, sprintf(" + A%s\u00b7cos(2\u03c0\u00b7%d\u00b7t/%.0f \u2212 \u03c6%s)",
                                       sub_h, h_s, mod$period, sub_h))
          }
        }
        cat(sym, "\n")
        leg <- c("M = MESOR")
        if(mod$n_harmonics == 1) {
          leg <- c(leg, "A = amplitude", "\u03c6 = acrophase (rad)")
        } else {
          leg <- c(leg, "A\u2095 = amplitude of harmonic h", "\u03c6\u2095 = acrophase of harmonic h (rad)")
        }
        leg <- c(leg, sprintf("T = %.0f h (period)", mod$period))
        if(trend_type == "linear")       leg <- c(leg, "\u03b2 = linear trend slope")
        else if(trend_type == "log")     leg <- c(leg, "\u03b2 = log trend slope")
        else if(trend_type == "exp_sat") leg <- c(leg, "A_sat = asymptote", "\u03c4 = time constant (h)")
        cat("  where:", paste(leg, collapse = ", "), "\n")
      }

      # Build and display fitted equation
      cat("\n--- Fitted Model Equation ---\n")
      eq <- sprintf("Y(t) = %.2f", pop$mean_mesor)

      # Add trend component
      trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
      if(trend_type == "linear" && !is.null(pop$indiv_means$trend_linear)) {
        beta <- pop$indiv_means$trend_linear
        if(beta >= 0) {
          eq <- paste0(eq, sprintf(" + %.3f·t", beta))
        } else {
          eq <- paste0(eq, sprintf(" - %.3f·t", abs(beta)))
        }
      } else if(trend_type == "log" && !is.null(pop$indiv_means$trend_log)) {
        beta <- pop$indiv_means$trend_log
        if(beta >= 0) {
          eq <- paste0(eq, sprintf(" + %.3f·log(t+1)", beta))
        } else {
          eq <- paste0(eq, sprintf(" - %.3f·log(t+1)", abs(beta)))
        }
      } else if(trend_type == "exp_sat" && !is.null(pop$indiv_means$A_sat) && !is.null(pop$indiv_means$tau)) {
        A_sat <- pop$indiv_means$A_sat
        tau <- pop$indiv_means$tau
        if(A_sat >= 0) {
          eq <- paste0(eq, sprintf(" + %.2f·(1 - e^(-t/%.1f))", A_sat, tau))
        } else {
          eq <- paste0(eq, sprintf(" - %.2f·(1 - e^(-t/%.1f))", abs(A_sat), tau))
        }
      }

      # Add harmonic components
      for(h in 1:mod$n_harmonics) {
        A <- pop$mean_amplitudes[h]
        phi <- pop$mean_acrophases_rad[h]  # In radians
        omega <- 2 * pi * h / mod$period

        # Convert to cos(ωt - φ) format
        if(A >= 0) {
          eq <- paste0(eq, sprintf(" + %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                   A, h, mod$period, phi))
        } else {
          eq <- paste0(eq, sprintf(" - %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                   abs(A), h, mod$period, phi))
        }
      }

      cat(eq, "\n")
      cat("where t = time in same units as period\n")

      # Show arithmetic means of individual parameters
      if(!is.null(pop$indiv_means)) {
        cat("\n--- Arithmetic Mean of Individual Parameters ---\n")
        cat(sprintf("MESOR:     %.3f (SD=%.3f)\n", pop$indiv_means$mesor, pop$indiv_means$mesor_sd))
        for(h in 1:mod$n_harmonics) {
          amp_key <- paste0("amplitude_", h)
          amp_sd_key <- paste0("amplitude_", h, "_sd")
          acro_key <- paste0("acrophase_time_", h)
          acro_sd_key <- paste0("acrophase_time_", h, "_sd")
          cat(sprintf("H%d Amplitude: %.3f (SD=%.3f)\n", h, pop$indiv_means[[amp_key]], pop$indiv_means[[amp_sd_key]]))
          cat(sprintf("H%d Acrophase: %.2f hours (SD=%.2f)\n", h, pop$indiv_means[[acro_key]], pop$indiv_means[[acro_sd_key]]))
        }
      }
      
      if(!is.null(mod$boot_results)) {
        cat(sprintf("\n95%% Bootstrap CIs (B=%d):\n", mod$boot_results$B))
        cat(sprintf("  MESOR:     [%.3f, %.3f]\n", 
                    mod$boot_results$mesor_ci[1], mod$boot_results$mesor_ci[2]))
        cat(sprintf("  Amplitude: [%.3f, %.3f]\n", 
                    mod$boot_results$amplitude_ci[1], mod$boot_results$amplitude_ci[2]))
        cat(sprintf("  Acrophase: [%.2f, %.2f] hours\n", 
                    mod$boot_results$acrophase_ci[1], mod$boot_results$acrophase_ci[2]))
      }
    }
    
    cat("\n--- Individual Parameter Summary ---\n")
    params <- mod$individual_params
    cat(sprintf("MESOR:     Mean=%.3f, SD=%.3f\n", mean(params$mesor), sd(params$mesor)))
    
    # Show trend parameters based on trend type
    if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
      cat(sprintf("Linear Trend (β): Mean=%.4f, SD=%.4f (units/hour)\n", 
                  mean(params$trend_linear, na.rm=TRUE), sd(params$trend_linear, na.rm=TRUE)))
      mean_trend <- mean(params$trend_linear, na.rm=TRUE)
      if(mean_trend > 0) {
        cat("           (Positive = increasing over time, e.g., increasing sleepiness)\n")
      } else if(mean_trend < 0) {
        cat("           (Negative = decreasing over time)\n")
      }
    } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
      cat(sprintf("Log Trend (β): Mean=%.4f, SD=%.4f (units/log-hour)\n", 
                  mean(params$trend_log, na.rm=TRUE), sd(params$trend_log, na.rm=TRUE)))
      mean_trend <- mean(params$trend_log, na.rm=TRUE)
      if(mean_trend > 0) {
        cat("           (Positive = increasing with diminishing rate)\n")
      } else if(mean_trend < 0) {
        cat("           (Negative = decreasing with diminishing rate)\n")
      }
    } else if(mod$trend_type == "two_process") {
      cat("\n--- Two-Process Homeostatic Parameters ---\n")
      if("beta_S" %in% names(params)) {
        cat(sprintf("β_S (S coefficient): Mean=%.3f, SD=%.3f (units/S-unit)\n",
                    mean(params$beta_S, na.rm=TRUE), sd(params$beta_S, na.rm=TRUE)))
        mean_beta_S <- mean(params$beta_S, na.rm=TRUE)
        if(mean_beta_S > 0) {
          cat("           (Positive = Y increases with sleep pressure)\n")
        } else if(mean_beta_S < 0) {
          cat("           (Negative = Y decreases with sleep pressure)\n")
        }
      }
      if("tau_w" %in% names(params)) {
        cat(sprintf("τ_w (wake time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau_w, na.rm=TRUE), sd(params$tau_w, na.rm=TRUE)))
        cat("           (Time for S to rise to ~63% toward S_max while awake)\n")
      }
      if("tau_s" %in% names(params)) {
        cat(sprintf("τ_s (sleep time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau_s, na.rm=TRUE), sd(params$tau_s, na.rm=TRUE)))
        cat("           (Time for S to decay to ~37% toward S_min while asleep)\n")
      }
    } else if(mod$trend_type == "exp_sat") {
      if("A_sat" %in% names(params)) {
        cat(sprintf("A_sat (asymptote): Mean=%.3f, SD=%.3f (units)\n",
                    mean(params$A_sat, na.rm=TRUE), sd(params$A_sat, na.rm=TRUE)))
      }
      if("tau" %in% names(params)) {
        cat(sprintf("τ (time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau, na.rm=TRUE), sd(params$tau, na.rm=TRUE)))
        cat("           (Time to reach ~63% of asymptote)\n")
      }
    }

    for(h in 1:mod$n_harmonics) {
      amp_col <- paste0("amplitude_", h)
      acro_col <- paste0("acrophase_time_", h)
      acro_rad_col <- paste0("acrophase_rad_", h)
      effective_period <- mod$period / h
      
      cat(sprintf("H%d Amplitude: Mean=%.3f, SD=%.3f\n", h, mean(params[[amp_col]]), sd(params[[amp_col]])))
      
      # Circular statistics for acrophase
      acro_rad <- params[[acro_rad_col]]
      circ_mean <- circular_mean(acro_rad)
      if(circ_mean < 0) circ_mean <- circ_mean + 2 * pi
      circ_mean_time <- circ_mean * effective_period / (2 * pi)
      circ_sd <- circular_sd(acro_rad)
      circ_sd_time <- if(!is.na(circ_sd)) circ_sd * effective_period / (2 * pi) else NA
      r_bar <- mean_resultant_length(acro_rad)
      
      cat(sprintf("H%d Acrophase: Circular mean=%.2f h, Circ.SD=%.2f h, r̄=%.3f\n", 
                  h, circ_mean_time, ifelse(is.na(circ_sd_time), NA, circ_sd_time), r_bar))
      cat(sprintf("            (Arithmetic mean=%.2f h, Linear SD=%.2f h)\n", 
                  mean(params[[acro_col]]), sd(params[[acro_col]])))
    }
    
    cat(sprintf("R-squared: Mean=%.3f, Range=[%.3f, %.3f]\n",
                mean(params$r_squared), min(params$r_squared), max(params$r_squared)))
    cat(sprintf("Significant rhythms (p<0.05): %d / %d (%.1f%%)\n",
                sum(params$p_value < 0.05), nrow(params),
                100 * sum(params$p_value < 0.05) / nrow(params)))

    # Model selection metrics
    if("aic" %in% names(params)) {
      cat("\n--- Model Selection Metrics ---\n")
      cat(sprintf("AIC (Akaike Information Criterion): Mean=%.2f, SD=%.2f\n",
                  mean(params$aic, na.rm = TRUE), sd(params$aic, na.rm = TRUE)))
      cat(sprintf("AICc (Corrected AIC): Mean=%.2f, SD=%.2f\n",
                  mean(params$aicc, na.rm = TRUE), sd(params$aicc, na.rm = TRUE)))
      cat(sprintf("BIC (Bayesian Information Criterion): Mean=%.2f, SD=%.2f\n",
                  mean(params$bic, na.rm = TRUE), sd(params$bic, na.rm = TRUE)))
      cat(sprintf("LOOCV RMSE (Leave-one-out CV): Mean=%.4f, SD=%.4f\n",
                  mean(params$loocv_rmse, na.rm = TRUE), sd(params$loocv_rmse, na.rm = TRUE)))
      cat("Note: Lower AIC/AICc/BIC/LOOCV values indicate better model fit\n")
    }

    # Variance decomposition: Relative importance of Process S and Process C
    if(!is.null(mod$pop_mean_fit) && !is.null(mod$pop_mean_fit$indiv_means) &&
       !is.null(mod$pop_mean_fit$indiv_means$r_squared_S) &&
       !is.null(mod$pop_mean_fit$indiv_means$r_squared_C)) {
      cat("\n--- Variance Decomposition: Relative Importance ---\n")
      indiv <- mod$pop_mean_fit$indiv_means

      cat(sprintf("R² from Process S (homeostatic): Mean=%.3f, SD=%.3f\n",
                  indiv$r_squared_S, indiv$r_squared_S_sd))
      cat(sprintf("R² from Process C (circadian):   Mean=%.3f, SD=%.3f\n",
                  indiv$r_squared_C, indiv$r_squared_C_sd))
      cat(sprintf("\nProportion of total R² explained by:\n"))
      cat(sprintf("  Process S: %.1f%% (SD=%.1f%%)\n",
                  indiv$percent_S, indiv$percent_S_sd))
      cat(sprintf("  Process C: %.1f%% (SD=%.1f%%)\n",
                  indiv$percent_C, indiv$percent_C_sd))

      # Interpretation helper (with NA check)
      if(!is.na(indiv$percent_S) && !is.na(indiv$percent_C)) {
        if(indiv$percent_S > indiv$percent_C) {
          cat("\nInterpretation: Homeostatic process (S) is the dominant component\n")
        } else if(indiv$percent_C > indiv$percent_S) {
          cat("\nInterpretation: Circadian rhythm (C) is the dominant component\n")
        } else {
          cat("\nInterpretation: Both processes contribute equally\n")
        }
      }
    }

    # Show group-specific statistics if groups exist
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 2) {
      cat("\n--- Group-Specific Parameters ---\n")
      # Print symbolic equation once before listing groups
      {
        cat("Model Equation (symbolic):\n  ")
        sym_g <- "Y(t) = M"
        grp_trend <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
        if(grp_trend == "linear")       sym_g <- paste0(sym_g, " + \u03b2\u00b7t")
        else if(grp_trend == "log")     sym_g <- paste0(sym_g, " + \u03b2\u00b7log(t+1)")
        else if(grp_trend == "exp_sat") sym_g <- paste0(sym_g, " + A_sat\u00b7(1 \u2212 e^(\u2212t/\u03c4))")
        if(mod$n_harmonics == 1) {
          sym_g <- paste0(sym_g, sprintf(" + A\u00b7cos(2\u03c0\u00b7t/%.0f \u2212 \u03c6)", mod$period))
        } else {
          for(h_s in 1:mod$n_harmonics) {
            sub_h <- intToUtf8(0x2080 + h_s)
            sym_g <- paste0(sym_g, sprintf(" + A%s\u00b7cos(2\u03c0\u00b7%d\u00b7t/%.0f \u2212 \u03c6%s)",
                                           sub_h, h_s, mod$period, sub_h))
          }
        }
        cat(sym_g, "\n")
        leg_g <- c("M = MESOR")
        if(mod$n_harmonics == 1) {
          leg_g <- c(leg_g, "A = amplitude", "\u03c6 = acrophase (rad)")
        } else {
          leg_g <- c(leg_g, "A\u2095 = amplitude of harmonic h", "\u03c6\u2095 = acrophase of harmonic h (rad)")
        }
        leg_g <- c(leg_g, sprintf("T = %.0f h (period)", mod$period))
        if(grp_trend == "linear")       leg_g <- c(leg_g, "\u03b2 = linear trend slope")
        else if(grp_trend == "log")     leg_g <- c(leg_g, "\u03b2 = log trend slope")
        else if(grp_trend == "exp_sat") leg_g <- c(leg_g, "A_sat = asymptote", "\u03c4 = time constant (h)")
        cat("  where:", paste(leg_g, collapse = ", "), "\n")
      }
      for(g_name in names(mod$group_fits)) {
        g <- mod$group_fits[[g_name]]
        cat(sprintf("\nGroup '%s' (n=%d):\n", g_name, g$n))
        cat(sprintf("  MESOR:     %.3f (SD=%.3f)\n", g$mean_mesor, g$sd_mesor))

        # Show trend parameters if present
        if(!is.null(g$trend_params) && length(g$trend_params) > 0) {
          for(param_name in names(g$trend_params)) {
            param_label <- switch(param_name,
                                 "trend_linear" = "Linear trend",
                                 "trend_log" = "Log trend",
                                 "A_sat" = "A_sat",
                                 "tau" = "τ",
                                 param_name)
            cat(sprintf("  %s: %.3f (SD=%.3f)\n",
                       param_label, g$trend_params[[param_name]]$mean,
                       g$trend_params[[param_name]]$sd))
          }
        }

        # Show sleep inertia parameters if present
        if(!is.null(g$inertia_params)) {
          cat(sprintf("  W₀:        %.3f (SD=%.3f)\n",
                     g$inertia_params$W0$mean, g$inertia_params$W0$sd))
          cat(sprintf("  τ_W:       %.2f h (SD=%.2f)\n",
                     g$inertia_params$tau_W$mean, g$inertia_params$tau_W$sd))
        }

        for(h in 1:mod$n_harmonics) {
          cat(sprintf("  H%d Amplitude: %.3f (SD=%.3f)\n", h, g$mean_amplitudes[h], g$sd_amplitudes[h]))
          cat(sprintf("  H%d Acrophase: %.2f hours\n", h, g$mean_acrophases_time[h]))
        }

        # Build and display fitted equation for this group
        cat("\n  Fitted equation:\n  ")
        g_eq <- sprintf("Y(t) = %.2f", g$mean_mesor)

        # Add trend component
        if(mod$trend_type == "linear" && !is.null(g$trend_params$trend_linear)) {
          beta <- g$trend_params$trend_linear$mean
          if(beta >= 0) {
            g_eq <- paste0(g_eq, sprintf(" + %.3f·t", beta))
          } else {
            g_eq <- paste0(g_eq, sprintf(" - %.3f·t", abs(beta)))
          }
        } else if(mod$trend_type == "log" && !is.null(g$trend_params$trend_log)) {
          beta <- g$trend_params$trend_log$mean
          if(beta >= 0) {
            g_eq <- paste0(g_eq, sprintf(" + %.3f·log(t+1)", beta))
          } else {
            g_eq <- paste0(g_eq, sprintf(" - %.3f·log(t+1)", abs(beta)))
          }
        } else if(mod$trend_type == "exp_sat" && !is.null(g$trend_params$A_sat) && !is.null(g$trend_params$tau)) {
          A_sat <- g$trend_params$A_sat$mean
          tau <- g$trend_params$tau$mean
          if(!is.na(A_sat) && !is.na(tau)) {
            if(A_sat >= 0) {
              g_eq <- paste0(g_eq, sprintf(" + %.2f·(1 - e^(-t/%.1f))", A_sat, tau))
            } else {
              g_eq <- paste0(g_eq, sprintf(" - %.2f·(1 - e^(-t/%.1f))", abs(A_sat), tau))
            }
          }
        }

        # Add harmonic components
        for(h in 1:mod$n_harmonics) {
          A <- g$mean_amplitudes[h]
          phi <- g$mean_acrophases_rad[h]  # In radians

          if(!is.na(A) && !is.na(phi)) {
            if(A >= 0) {
              g_eq <- paste0(g_eq, sprintf(" + %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                       A, h, mod$period, phi))
            } else {
              g_eq <- paste0(g_eq, sprintf(" - %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                       abs(A), h, mod$period, phi))
            }
          }
        }

        cat(g_eq, "\n")

        # Variance decomposition for this group
        if(!is.null(g$variance_decomp)) {
          vd <- g$variance_decomp
          cat(sprintf("  Variance Decomposition:\n"))
          cat(sprintf("    Process S: %.1f%% (SD=%.1f%%)\n", vd$percent_S, vd$percent_S_sd))
          cat(sprintf("    Process C: %.1f%% (SD=%.1f%%)\n", vd$percent_C, vd$percent_C_sd))
        }
      }
    }
  })
  
  # Parameters table
  output$harmonic_parameters_table <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Build summary table with all harmonics
    param_names <- c("MESOR")
    mean_vals <- c(mean(params$mesor))
    sd_vals <- c(sd(params$mesor))
    min_vals <- c(min(params$mesor))
    max_vals <- c(max(params$mesor))
    
    for(h in 1:mod$n_harmonics) {
      amp_col <- paste0("amplitude_", h)
      acro_col <- paste0("acrophase_time_", h)
      
      param_names <- c(param_names, paste0("Amplitude H", h), paste0("Acrophase H", h, " (hours)"))
      mean_vals <- c(mean_vals, mean(params[[amp_col]]), mean(params[[acro_col]]))
      sd_vals <- c(sd_vals, sd(params[[amp_col]]), sd(params[[acro_col]]))
      min_vals <- c(min_vals, min(params[[amp_col]]), min(params[[acro_col]]))
      max_vals <- c(max_vals, max(params[[amp_col]]), max(params[[acro_col]]))
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
      group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')
      group_colors_rgba <- c('rgba(178,34,34,', 'rgba(70,130,180,', 'rgba(34,139,34,', 
                             'rgba(128,0,128,', 'rgba(255,140,0,', 'rgba(139,69,19,')
      
      if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
         !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        
        group_var <- values$covariates[[input$harmonic_group_var]]
        groups <- names(mod$group_fits)
        
        for(i in seq_along(mod$individual_fits)) {
          fit_i <- mod$individual_fits[[i]]
          if(!is.null(fit_i) && fit_i$success) {
            pred_i <- predict_cosinor(fit_i, time_fine)
            grp <- as.character(group_var[i])
            grp_idx <- which(groups == grp)
            line_color <- if(length(grp_idx) > 0) paste0(group_colors_rgba[grp_idx], '0.4)') else 'rgba(100,100,100,0.3)'
            
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
                               line = list(color = group_colors_hex[g_idx], width = 3), 
                               name = paste("Mean:", g_name))
        }
        
        # Show harmonic components for grouped data (same as "mean" view)
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          params$group <- group_var[params$subject]
          
          # Lighter versions of group colors for components
          group_comp_colors <- c('#FF6B6B', '#87CEEB', '#90EE90', '#DDA0DD', '#FFB347', '#D2B48C')
          
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
                                     line = list(color = group_comp_colors[g_idx], width = 1.5, dash = 'dot'),
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
                                     line = list(color = group_colors_hex[g_idx], width = 1.5, dash = 'dash'),
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
        group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')  # firebrick, steelblue, forestgreen, purple, orange, brown
        group_colors_rgba <- c('rgba(178,34,34,', 'rgba(70,130,180,', 'rgba(34,139,34,', 
                               'rgba(128,0,128,', 'rgba(255,140,0,', 'rgba(139,69,19,')
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
                                 line = list(color = group_colors_hex[color_idx], width = 3),
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
            ci_color <- paste0(group_colors_rgba[color_idx], '0.25)')
            
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
          group_comp_colors <- c('#FF6B6B', '#87CEEB', '#90EE90', '#DDA0DD', '#FFB347', '#D2B48C')
          
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
                                     line = list(color = group_comp_colors[g_idx], width = 1.5, dash = 'dot'),
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
                                     line = list(color = group_colors_hex[g_idx], width = 1.5, dash = 'dash'),
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
              pt_color <- if(length(grp_idx) > 0) paste0(group_colors_rgba[grp_idx], '0.3)') else 'rgba(100,100,100,0.2)'
              
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
    
    # Create better x-axis with clock time labels
    time_range <- range(mod$time_vec)
    
    # Generate tick values and labels
    if(max(mod$time_vec) > mod$period) {
      # Times wrap around - create clock-style labels (modulo period)
      tick_interval <- if(mod$period <= 12) 1 else if(mod$period <= 24) 2 else if(mod$period <= 48) 4 else mod$period / 12
      tick_vals <- seq(floor(time_range[1]), ceiling(time_range[2]), by = tick_interval)
      tick_text <- sapply(tick_vals, function(t) {
        t_mod <- t %% mod$period
        if(mod$period == 24) {
          paste0(t_mod, ":00")
        } else {
          round(t_mod, 1)
        }
      })
      x_title <- paste0("Time (period = ", mod$period, ", spanning wrap-around)")
    } else {
      tick_vals <- NULL
      tick_text <- NULL
      x_title <- paste("Time (period =", mod$period, ")")
    }
    
    x_axis <- list(title = x_title)
    if(!is.null(tick_vals)) {
      x_axis$tickmode <- "array"
      x_axis$tickvals <- tick_vals
      x_axis$ticktext <- tick_text
    }
    
    p %>% layout(
      title = "Harmonic Regression Fit",
      xaxis = x_axis,
      yaxis = list(title = "Response"),
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
    theta_deg <- params[[acro_rad_col]] * 180 / pi
    r <- params[[amp_col]]
    
    p <- plot_ly(type = 'scatterpolar', mode = 'markers')
    
    # Check if we have group fits - color points by group
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {

      group_var <- values$covariates[[input$harmonic_group_var]]
      group_colors <- c('firebrick', 'steelblue', 'forestgreen', 'purple', 'orange', 'brown')

      groups <- names(mod$group_fits)
      for(g_idx in seq_along(groups)) {
        g_name <- groups[g_idx]
        g_mask <- !is.na(group_var[params$subject]) & group_var[params$subject] == g_name

        if(sum(g_mask) > 0) {
          p <- p %>% add_trace(
            r = r[g_mask], theta = theta_deg[g_mask],
            type = 'scatterpolar', mode = 'markers',
            marker = list(size = 8, color = group_colors[g_idx], opacity = 0.7),
            name = paste("Group:", g_name)
          )
        }

        # Add group mean vector for selected harmonic
        g_fit <- mod$group_fits[[g_name]]
        acro_deg <- g_fit$mean_acrophases_rad[h] * 180 / pi
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, g_fit$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          type = 'scatterpolar', mode = 'lines+markers',
          line = list(color = group_colors[g_idx], width = 3),
          marker = list(size = 12, color = group_colors[g_idx], symbol = 'diamond'),
          name = paste("Mean:", g_name)
        )
      }

      # Add population mean vector if requested (even when groups present)
      if(isTRUE(input$polar_show_mean) && !is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        acro_deg <- pop$mean_acrophases_rad[h] * 180 / pi
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
        acro_deg <- pop$mean_acrophases_rad[h] * 180 / pi
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
        ellipse_theta_deg <- ellipse_theta_rad * 180 / pi
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
    effective_period <- mod$period / h
    n_ticks <- min(12, effective_period)
    tick_step <- effective_period / n_ticks
    tick_vals <- seq(0, 360 - 360/n_ticks, by = 360/n_ticks)
    tick_labels <- sprintf("%.1fh", seq(0, effective_period - tick_step, by = tick_step))
    
    p %>% layout(
      title = paste("Acrophase Polar Plot - Harmonic", h, "(period =", round(effective_period, 1), "h)"),
      polar = list(
        radialaxis = list(title = "Amplitude"),
        angularaxis = list(
          direction = "clockwise",
          rotation = 90,
          tickmode = "array",
          tickvals = tick_vals,
          ticktext = tick_labels
        )
      ),
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
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = .data[[acro_col]], fill = as.factor(group))) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 12) +
        scale_fill_brewer(palette = "Set1", name = "Group") +
        scale_x_continuous(limits = c(0, effective_period)) +
        theme_minimal() +
        labs(title = paste0("Distribution of Acrophases (H", h, ") by Group"), 
             x = paste("Acrophase (hours, period =", round(effective_period, 1), ")"), y = "Count")
      
      ggplotly(g)
    } else {
      plot_ly(x = params[[acro_col]], type = "histogram",
              marker = list(color = 'firebrick', line = list(color = 'white', width = 1))) %>%
        layout(title = paste0("Distribution of Acrophases (H", h, ")"),
               xaxis = list(title = paste("Acrophase (hours, period =", round(effective_period, 1), ")"),
                            range = c(0, effective_period)),
               yaxis = list(title = "Count"))
    }
  })
  
  # MESOR by group plot
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
        labs(title = "MESOR by Group", x = "Group", y = "MESOR") +
        theme(legend.position = "none")
      
      ggplotly(g)
    } else {
      plot_ly(y = params$mesor, type = "box", 
              marker = list(color = 'forestgreen'),
              boxpoints = "all", jitter = 0.3) %>%
        layout(title = "Distribution of MESOR",
               yaxis = list(title = "MESOR"))
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
      MESOR = round(params$mesor, 3),
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
      
      # Get R-squared stats for this group
      grp_rsq <- params$r_squared[params$group == g_name]
      rsq_mean <- mean(grp_rsq, na.rm = TRUE)
      rsq_se <- sd(grp_rsq, na.rm = TRUE) / sqrt(length(grp_rsq))
      
      # Get amplitude and acrophase for selected harmonic from individual params
      grp_amp <- params[[amp_col]][params$group == g_name]
      grp_acro_rad <- params[[acro_rad_col]][params$group == g_name]
      
      # Compute circular mean and SE for acrophase
      circ_mean_rad <- circular_mean(grp_acro_rad)
      if(circ_mean_rad < 0) circ_mean_rad <- circ_mean_rad + 2 * pi
      circ_mean_time <- circ_mean_rad * effective_period / (2 * pi)
      circ_se_rad <- circular_se(grp_acro_rad)
      circ_se_time <- if(!is.na(circ_se_rad)) circ_se_rad * effective_period / (2 * pi) else NA
      
      group_df <- rbind(group_df, 
                        data.frame(group = g_name, parameter = "MESOR", 
                                   value = g$mean_mesor, se = g$sd_mesor / sqrt(g$n)),
                        data.frame(group = g_name, parameter = paste0("Amplitude (H", h, ")"), 
                                   value = mean(grp_amp, na.rm = TRUE), 
                                   se = sd(grp_amp, na.rm = TRUE) / sqrt(length(grp_amp))),
                        data.frame(group = g_name, parameter = paste0("Acrophase (H", h, ")"), 
                                   value = circ_mean_time, se = circ_se_time),
                        data.frame(group = g_name, parameter = "R-squared", 
                                   value = rsq_mean, se = rsq_se))
      
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

    # Build named color vector matching the line graph palette
    group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')
    group_names <- names(mod$group_fits)
    named_colors <- setNames(group_colors_hex[seq_along(group_names)], group_names)

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
    cat(sprintf("Harmonic: H%d (period = %.1f h)\n\n", h, effective_period))
    
    # Get group variable and params
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- group_var[params$subject]
      
      # Remove rows with missing group or key parameters
      params <- params[!is.na(params$group) & !is.na(params$mesor) & 
                         !is.na(params[[amp_col]]) & !is.na(params[[acro_col]]), ]
      
      n_groups <- length(unique(params$group))
      cat(sprintf("Groups: %d, Total N: %d (after removing missings)\n\n", n_groups, nrow(params)))
      
      cat("--- MESOR Comparison ---\n")
      if(n_groups == 2) {
        tt <- t.test(mesor ~ group, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_mesor <- aov(mesor ~ group, data = params)
        cat("ANOVA:\n")
        print(summary(aov_mesor))
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
          if(n_groups == 2) {
            tt <- t.test(trend_formula, data = params)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
          } else {
            aov_trend <- aov(trend_formula, data = params)
            cat("ANOVA:\n")
            print(summary(aov_trend))
          }
          
          # Group means for trend
          for(g in unique(params$group)) {
            g_trend <- params[[trend_col]][params$group == g]
            cat(sprintf("  %s: mean = %.4f (SD = %.4f, n=%d)\n", 
                        g, mean(g_trend, na.rm = TRUE), sd(g_trend, na.rm = TRUE), length(g_trend)))
          }
          
          # For exp_sat, also compare tau
          if(trend_type == "exp_sat" && "tau" %in% names(params)) {
            cat("\n--- Time Constant (τ) Comparison ---\n")
            tau_formula <- as.formula("tau ~ group")
            if(n_groups == 2) {
              tt <- t.test(tau_formula, data = params)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_tau <- aov(tau_formula, data = params)
              cat("ANOVA:\n")
              print(summary(aov_tau))
            }
            
            for(g in unique(params$group)) {
              g_tau <- params$tau[params$group == g]
              cat(sprintf("  %s: mean τ = %.2f h (SD = %.2f, n=%d)\n", 
                          g, mean(g_tau, na.rm = TRUE), sd(g_tau, na.rm = TRUE), length(g_tau)))
            }
          }
        }
      }

      cat(sprintf("\n--- Amplitude (H%d) Comparison ---\n", h))
      amp_formula <- as.formula(paste(amp_col, "~ group"))
      if(n_groups == 2) {
        tt <- t.test(amp_formula, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_amp <- aov(amp_formula, data = params)
        cat("ANOVA:\n")
        print(summary(aov_amp))
      }
      
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

      if(!is.na(ww$F)) {
        cat(sprintf("    F(%d, %d) = %.3f, p = %.4f\n", ww$df1, ww$df2, ww$F, ww$p))
        cat(sprintf("    Mean resultant length (r̄) = %.3f", ww$r_bar))
        if(ww$r_bar >= 0.7) {
          cat(" (high concentration)\n")
        } else if(ww$r_bar >= 0.45) {
          cat(" (moderate concentration)\n")
        } else {
          cat(" (low concentration - interpret with caution)\n")
        }
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
        cat(sprintf("    %s: mean = %.2f h, circ.SD = %.2f h, r̄ = %.3f (n=%d)\n",
                    g, g_mean_time, ifelse(is.na(g_sd_time), NA, g_sd_time), g_r, length(g_angles)))
      }

      # (2) Hotelling's T² test
      cat(sprintf("\n(2) Hotelling's T² on (beta_cos_%d, beta_sin_%d) (amplitude-weighted):\n", h, h))
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
        acro_time <- acro_h * effective_period / (2 * pi)
        amp_mean <- sqrt(mean(bc_ok)^2 + mean(bs_ok)^2)
        cat(sprintf("    %s: amplitude-weighted mean = %.2f h, mean amplitude = %.3f (n=%d)\n",
                    g, acro_time, amp_mean, sum(ok)))
      }
      
      cat("\n--- R-squared Comparison ---\n")
      if(n_groups == 2) {
        tt <- t.test(r_squared ~ group, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_rsq <- aov(r_squared ~ group, data = params)
        cat("ANOVA:\n")
        print(summary(aov_rsq))
      }
      
      # Trend parameter comparison based on trend type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        cat("\n--- Linear Trend (β) Comparison ---\n")
        params_trend <- params[!is.na(params$trend_linear), ]
        
        if(nrow(params_trend) >= 4) {
          if(n_groups == 2) {
            tt <- t.test(trend_linear ~ group, data = params_trend)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            
            # Effect size (Cohen's d)
            groups <- unique(params_trend$group)
            g1 <- params_trend$trend_linear[params_trend$group == groups[1]]
            g2 <- params_trend$trend_linear[params_trend$group == groups[2]]
            pooled_sd <- sqrt(((length(g1)-1)*sd(g1)^2 + (length(g2)-1)*sd(g2)^2) / (length(g1)+length(g2)-2))
            cohens_d <- (mean(g1) - mean(g2)) / pooled_sd
            cat(sprintf("Cohen's d = %.3f\n", cohens_d))
          } else {
            aov_trend <- aov(trend_linear ~ group, data = params_trend)
            cat("ANOVA:\n")
            print(summary(aov_trend))
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
          if(n_groups == 2) {
            tt <- t.test(trend_log ~ group, data = params_trend)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
          } else {
            aov_trend <- aov(trend_log ~ group, data = params_trend)
            cat("ANOVA:\n")
            print(summary(aov_trend))
          }
          
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
            if(n_groups == 2) {
              tt <- t.test(A_sat ~ group, data = params_trend)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_asat <- aov(A_sat ~ group, data = params_trend)
              cat("ANOVA:\n")
              print(summary(aov_asat))
            }
            
            cat("\nGroup statistics (units):\n")
            for(g in unique(params_trend$group)) {
              g_asat <- params_trend$A_sat[params_trend$group == g]
              cat(sprintf("  %s: mean = %.3f, SD = %.3f (n=%d)\n", 
                          g, mean(g_asat, na.rm=TRUE), sd(g_asat, na.rm=TRUE), length(g_asat)))
            }
          }
        }
        
        # Compare tau (time constant)
        if("tau" %in% names(params)) {
          cat("\n--- Time Constant (τ) Comparison ---\n")
          params_trend <- params[!is.na(params$tau), ]
          
          if(nrow(params_trend) >= 4) {
            if(n_groups == 2) {
              tt <- t.test(tau ~ group, data = params_trend)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_tau <- aov(tau ~ group, data = params_trend)
              cat("ANOVA:\n")
              print(summary(aov_tau))
            }
            
            cat("\nGroup statistics (hours):\n")
            for(g in unique(params_trend$group)) {
              g_tau <- params_trend$tau[params_trend$group == g]
              cat(sprintf("  %s: mean τ = %.2f h, SD = %.2f (n=%d)\n", 
                          g, mean(g_tau, na.rm=TRUE), sd(g_tau, na.rm=TRUE), length(g_tau)))
            }
            cat("\nNote: τ represents time to reach ~63% of asymptotic level.\n")
          }
        }
      }
    }
  })
