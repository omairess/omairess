# ==========================================================================
# server/01_helpers_time.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 1152-1382  (clock-time helpers used by every plot)
# ==========================================================================
  # Helper function to extract hour values from column names
  # Supports multiple formats:
  #   KSS_9u_dag1 -> 9
  #   Base9h -> 9
  #   Base7h30 -> 7.5
  #   R1_5h -> 5
  #   8.0, 8.25, 8.5, 8.75 -> decimal hours (8:00, 8:15, 8:30, 8:45)
  #   X8.0, X8.25 -> decimal hours with X prefix (R's default for numeric column names)
  extract_hour_from_colname <- function(col_name) {
    # Pattern 0: Pure decimal number or with X prefix (e.g., "8.0", "8.25", "X8.0", "X8.25")
    # These represent decimal hours where .25=15min, .5=30min, .75=45min
    # Remove leading X if present (R adds this when column names start with numbers)
    clean_name <- gsub("^X", "", col_name)

    # Check if it's a pure decimal number
    if(grepl("^[0-9]+\\.?[0-9]*$", clean_name)) {
      return(as.numeric(clean_name))
    }

    # Pattern 1: hour with minutes "hMM" (e.g., "7h30"=7.5, "20h15"=20.25, "20h45"=20.75, "20h00"=20)
    if(grepl("[0-9]+h[0-9]+", col_name)) {
      hour_match <- regexpr("[0-9]+h[0-9]+", col_name)
      hour_str <- regmatches(col_name, hour_match)
      parts <- strsplit(hour_str, "h")[[1]]
      hour <- as.numeric(parts[1]) + as.numeric(parts[2]) / 60
      return(hour)
    }

    # Pattern 2: hour with "h" only, no minutes (e.g., "9h" = 9 hours)
    if(grepl("[0-9]+h", col_name)) {
      hour_match <- regexpr("[0-9]+h", col_name)
      hour_str <- regmatches(col_name, hour_match)
      hour <- as.numeric(gsub("h", "", hour_str))
      return(hour)
    }

    # Pattern 3: hour with "u" (e.g., "9u" = 9 hours)
    if(grepl("[0-9]+u", col_name)) {
      hour_match <- regexpr("[0-9]+u", col_name)
      hour_str <- regmatches(col_name, hour_match)
      hour <- as.numeric(gsub("u", "", hour_str))
      return(hour)
    }

    # Fallback: try to extract any number
    num_match <- regexpr("[0-9]+", col_name)
    if(num_match > 0) {
      return(as.numeric(regmatches(col_name, num_match)))
    }

    return(NA)
  }
  
  # Helper function to extract numeric time values from column names
  extract_time_values <- function(col_names) {
    tryCatch({
      # For chronological data (like circadian hours), preserve original order
      # Use sequential numbering to maintain chronology
      # The actual hour labels will be shown separately
      return(1:length(col_names))
    }, error = function(e) {
      return(1:length(col_names))
    })
  }
  
  # Helper function to get plotting time points
  get_plot_time <- function() {
    if(!is.null(values$time_numeric)) {
      return(values$time_numeric)
    } else if(!is.null(values$time_labels)) {
      return(1:length(values$time_labels))
    } else if(!is.null(values$data)) {
      return(1:ncol(values$data))
    } else {
      return(NULL)
    }
  }
  
  # Helper function to get time axis label
  get_time_label <- function() {
    if(!is.null(values$time_labels)) {
      first_label <- values$time_labels[1]
      # Check for circadian/hour data
      if(grepl("[0-9]+h30", first_label)) return("Hour")  # e.g., 7h30
      if(grepl("[0-9]+h", first_label)) return("Hour")    # e.g., 9h, 11h
      if(grepl("[0-9]+u", first_label)) return("Hour")    # e.g., 9u
      if(grepl("hour|hr", tolower(first_label))) return("Hour")
      if(grepl("day|d_", tolower(first_label))) return("Day")
      if(grepl("time|t_|^t[0-9]", tolower(first_label))) return("Time Point")
      if(grepl("month|mon", tolower(first_label))) return("Month")
      if(grepl("year|yr", tolower(first_label))) return("Year")
    }
    return("Measurement Point")
  }
  
  # Helper function to format plotly x-axis with time labels
  # Calculate actual time positions (accounting for day boundaries)
  calculate_time_positions <- function(hour_labels) {
    if(is.null(hour_labels) || length(hour_labels) == 0) {
      return(NULL)
    }
    
    n <- length(hour_labels)
    cumulative_hours <- numeric(n)
    cumulative_hours[1] <- 0  # Start at 0
    
    for(i in 2:n) {
      prev_hour <- hour_labels[i-1]
      curr_hour <- hour_labels[i]
      
      # Calculate time difference
      if(curr_hour >= prev_hour) {
        # Same day progression (e.g., 9 -> 11)
        diff <- curr_hour - prev_hour
      } else {
        # Crossed midnight (e.g., 23 -> 1)
        diff <- (24 - prev_hour) + curr_hour
      }
      
      cumulative_hours[i] <- cumulative_hours[i-1] + diff
    }
    
    # Normalize to 0-1 scale
    if(max(cumulative_hours) > 0) {
      normalized <- cumulative_hours / max(cumulative_hours)
    } else {
      normalized <- cumulative_hours
    }
    
    return(normalized)
  }
  
  format_plotly_time_axis <- function(p, time_grid = NULL, tick_step_hours = NULL) {
    if(!is.null(values$time_labels) && length(values$time_labels) > 0) {
      # Get hour labels if available
      hour_labels <- get_hour_labels()
      
      # Calculate actual time positions (accounting for midnight crossings)
      if(!is.null(hour_labels)) {
        # Calculate positions based on actual time progression
        actual_positions <- calculate_time_positions(hour_labels)
        
        if(!is.null(actual_positions)) {
          time_grid <- actual_positions
        } else if(is.null(time_grid)) {
          time_grid <- seq(0, 1, length.out = length(values$time_labels))
        }
        
        # Use actual hour values as labels
        tick_text <- sapply(hour_labels, decimal_to_hhmm)

        # Snap to round-hour multiples when a step is specified
        step <- as.numeric(tick_step_hours)
        if (!is.null(tick_step_hours) && !is.na(step) && step > 0) {
          keep <- abs(round(hour_labels / step) * step - hour_labels) < 1e-6
          if (any(keep)) {
            tick_vals_subset <- time_grid[keep]
            tick_text_subset <- tick_text[keep]
          } else {
            tick_vals_subset <- time_grid
            tick_text_subset <- tick_text
          }
        } else {
          # "All" (step == 0) or no step: show every point (thin if very many)
          n_labels <- length(tick_text)
          if (n_labels <= 60) {
            tick_vals_subset <- time_grid
            tick_text_subset <- tick_text
          } else {
            thin <- ceiling(n_labels / 50)
            indices <- seq(1, n_labels, by = thin)
            tick_vals_subset <- time_grid[indices]
            tick_text_subset <- tick_text[indices]
          }
        }
      } else {
        # Fall back to evenly-spaced positions and column names
        if(is.null(time_grid)) {
          time_grid <- seq(0, 1, length.out = length(values$time_labels))
        }
        tick_vals_subset <- time_grid
        tick_text_subset <- values$time_labels
      }
      
      # Apply custom axis formatting
      # Always use 90Â° rotation for readability
      p <- p %>% layout(
        xaxis = list(
          tickmode = "array",
          tickvals = tick_vals_subset,
          ticktext = tick_text_subset,
          tickangle = -90,  # Rotate labels vertically (90Â°)
          title = get_time_label()
        )
      )
    }
    return(p)
  }
  
  
  # Helper function to get hour labels for x-axis ticks
  get_hour_labels <- function() {
    if(!is.null(values$time_labels)) {
      # Extract hour values from column names
      hours <- sapply(values$time_labels, extract_hour_from_colname)
      if(!all(is.na(hours))) {
        return(hours)
      }
    }
    return(NULL)
  }

  # Convert decimal hours to "HHhMM" display string (e.g. 20.25 -> "20h15")
  decimal_to_hhmm <- function(h) {
    h_int <- floor(h)
    m_int <- round((h - h_int) * 60)
    if (m_int == 60) { h_int <- h_int + 1; m_int <- 0 }
    sprintf("%dh%02d", h_int, m_int)
  }

  # Map a vector of normalized x positions (0-1) to formatted "HHhMM" hover strings.
  # Uses approx() to interpolate against the original measurement positions,
  # so it works for both raw data grids and smooth 100-point interpolated grids.
  hover_time_labels <- function(x_vals) {
    hl <- get_hour_labels()
    if (is.null(hl)) return(as.character(round(x_vals, 3)))
    tp <- calculate_time_positions(hl)
    if (is.null(tp)) return(as.character(round(x_vals, 3)))
    interp_hours <- approx(tp, hl, xout = x_vals, rule = 2)$y
    sapply(interp_hours, decimal_to_hhmm)
  }
