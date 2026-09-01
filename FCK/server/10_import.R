# ==============================================================================
# server/10_import.R — SHARED data import + variable selection
#
# Hand-merged union of:
#   WaPaa1_3.R  1497-1585 (load_data), 1588-1636 (var_select_container),
#               1639-1743 (apply_selection), 1746-1814 (generate_sample),
#               1385-1403 (data_status), 1846-1899 (data_preview)
#   CIRCAREG.R   777-797  (excel_sheet_selector), 800-829 (load_data),
#                831-863  (generate_sample), 866-878 (var_select_container),
#                880-896  (apply_selection), 899-905 (data_preview)
#
# What the union changes relative to each original, and why:
#   * separator: CIRCAREG made the user pick one, WaPaa sniffed it. Both are
#     kept, with sniffing as the default, so neither app's files break.
#   * Excel: CIRCAREG could pick a sheet, WaPaa always read the first. The
#     sheet picker is kept.
#   * variable selection: WaPaa picked "group variables", CIRCAREG picked
#     "scalar variables". They are the same columns used for two purposes, so
#     there is now ONE picker feeding both values$covariates (original types,
#     for FoSR/SoFR/cosinor predictors) and values$group_variables (factors,
#     for fANOVA/clustering/group comparisons).
#   * time values: WaPaa's extract_time_values() runs here for every dataset,
#     so values$time_numeric is available to all tabs (the cosinor tab can now
#     reuse it instead of re-detecting times from the same column names).
#   * sample data: the two generators produced different datasets (WaPaa: 100
#     points on 0-1 with groups; CIRCAREG: 24 hourly points with covariates and
#     a binary outcome). One dataset has to serve every tab, so the merged
#     generator makes a 24-hour circadian set WITH group structure, scalar
#     covariates and a binary outcome.
# ==============================================================================

# --- Excel sheet picker (from CIRCAREG) --------------------------------------
output$excel_sheet_selector <- renderUI({
  req(input$datafile)
  file_ext <- tools::file_ext(input$datafile$name)

  if(tolower(file_ext) %in% c("xls", "xlsx")) {
    tryCatch({
      if(!requireNamespace("readxl", quietly = TRUE)) {
        return(div(style = "color: red;",
                   "readxl package required for Excel files. Install with: install.packages('readxl')"))
      }
      sheet_names <- readxl::excel_sheets(input$datafile$datapath)
      selectInput("excel_sheet", "Select Sheet:",
                  choices = sheet_names, selected = sheet_names[1])
    }, error = function(e) {
      div(style = "color: red;", paste("Error reading Excel file:", e$message))
    })
  }
})

# --- 1. Load the raw file ----------------------------------------------------
observeEvent(input$load_data, {
  if(is.null(input$datafile)) {
    showNotification("Please select a file first!", type = "warning", duration = 5)
    return()
  }

  ext <- tolower(tools::file_ext(input$datafile$name))

  tryCatch({
    if(ext %in% c("xls", "xlsx")) {
      if(!requireNamespace("readxl", quietly = TRUE)) {
        showNotification("Package 'readxl' is required for Excel files. Install with: install.packages('readxl')",
                         type = "error", duration = 10)
        return()
      }
      # CIRCAREG's sheet choice; WaPaa always took the first sheet.
      sheet_to_read <- if(!is.null(input$excel_sheet)) input$excel_sheet else 1
      raw_data <- suppressWarnings(as.data.frame(readxl::read_excel(
        input$datafile$datapath,
        sheet = sheet_to_read,
        col_names = input$header,
        guess_max = 10000)))
      cat("Read Excel file:", input$datafile$name, "sheet:", sheet_to_read, "\n")

    } else {
      # Separator: explicit (CIRCAREG) or sniffed from the first line (WaPaa).
      sep <- input$sep
      if(is.null(sep) || sep == "auto") {
        sep <- if(ext %in% c("txt", "tsv")) "\t" else ","
        first_line <- readLines(input$datafile$datapath, n = 1)
        if(grepl("\t", first_line)) sep <- "\t"
        else if(grepl(";", first_line)) sep <- ";"
        cat("Auto-detected separator:", if(sep == "\t") "<TAB>" else sep, "\n")
      }
      # check.names = FALSE keeps decimal column names such as 8.25 intact;
      # both apps relied on that for time-of-day column names.
      raw_data <- read.csv(input$datafile$datapath,
                           header = input$header,
                           sep = sep,
                           stringsAsFactors = FALSE,
                           check.names = FALSE,
                           quote = "\"")
    }

    cat("Read raw file dimensions:", nrow(raw_data), "x", ncol(raw_data), "\n")

    # Duplicate column names would break column selection (e.g. "8h" twice).
    orig_names <- colnames(raw_data)
    if(any(duplicated(orig_names))) {
      unique_names <- make.unique(orig_names, sep = "_")
      colnames(raw_data) <- unique_names
      cat("Note:", sum(duplicated(orig_names)),
          "duplicate column name(s) found and made unique\n")
    }

    # "Long" in both apps means subjects in COLUMNS: transpose to subjects in rows.
    if(input$data_format == "long") {
      raw_data <- as.data.frame(t(raw_data))
      cat("Transposed 'Long' format to Wide. New dims:",
          nrow(raw_data), "x", ncol(raw_data), "\n")
    }

    values$raw_df <- raw_data
    values$uploaded_data <- raw_data   # RM-ANOVA pickers read this one

    # A new file voids everything downstream, in every family.
    values$data <- NULL
    values$covariates <- NULL
    values$group_labels <- NULL
    values$group_variables <- NULL
    values$selected_group_vars <- NULL
    values$time_numeric <- NULL
    values$time_clock <- NULL
    fck_reset_analyses(values)

    showNotification("File loaded. Please select variables below.",
                     type = "message", duration = 5)

  }, error = function(e) {
    cat("Error in load_data:", e$message, "\n")
    showNotification(paste("Error loading data:", e$message), type = "error", duration = 10)
  })
})

# --- 2. Variable selection ---------------------------------------------------
output$var_select_container <- renderUI({
  req(values$raw_df)

  cols <- colnames(values$raw_df)
  numeric_cols <- sapply(values$raw_df, is.numeric)
  default_data_cols <- cols[numeric_cols]
  if(length(default_data_cols) == 0 && length(cols) > 1) {
    default_data_cols <- cols[-1]
  }

  tagList(
    h4("Select Variables from Uploaded Data"),

    pickerInput(
      inputId = "sel_data_vars",
      label = "Select Time Series/Function Data Columns (the curves):",
      choices = cols,
      selected = default_data_cols,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `selected-text-format` = "count > 5",
        `preserve-selected-order` = TRUE
      ),
      multiple = TRUE
    ),
    helpText("Column order is preserved as it appears in your data file — that order",
             "is the time order every analysis uses."),

    pickerInput(
      inputId = "sel_cov_vars",
      label = "Select Scalar Variables (grouping factors / predictors / response):",
      choices = cols,
      selected = NULL,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `none-selected-text` = "None selected"
      ),
      multiple = TRUE
    ),
    helpText(HTML(
      "These serve <b>both</b> families: as predictors/response in FoSR, SoFR and
       cosinor regression, and as grouping factors in functional ANOVA, cluster
       composition and cosinor group tests. The first one selected is the primary
       grouping variable; each tab can pick a different one.")),

    actionButton("apply_selection", "Confirm & Process Data",
                 class = "btn-success", icon = icon("check"))
  )
})

observeEvent(input$apply_selection, {
  req(values$raw_df, input$sel_data_vars)

  tryCatch({
    if(length(input$sel_data_vars) < 2) {
      showNotification("Please select at least 2 time points/columns for analysis.",
                       type = "warning")
      return()
    }

    # Keep the columns in FILE order, not selection order (chronology matters).
    all_cols <- colnames(values$raw_df)
    data_cols <- all_cols[all_cols %in% input$sel_data_vars]

    temp_data <- values$raw_df[, data_cols, drop = FALSE]

    temp_data_mat <- tryCatch({
      data.matrix(temp_data)
    }, error = function(e) {
      mat <- matrix(NA, nrow = nrow(temp_data), ncol = ncol(temp_data))
      for(j in 1:ncol(temp_data)) mat[, j] <- as.numeric(as.character(temp_data[, j]))
      mat
    })
    # Force numeric storage: fda chokes on integer/character matrices and it
    # was the cause of WaPaa's -Inf R-squared bug.
    if(typeof(temp_data_mat) != "double" && typeof(temp_data_mat) != "integer") {
      mode(temp_data_mat) <- "numeric"
    }

    # ---- ONE scalar picker, TWO consumers -----------------------------------
    if(!is.null(input$sel_cov_vars) && length(input$sel_cov_vars) > 0) {
      cov_cols <- all_cols[all_cols %in% input$sel_cov_vars]
      # original types: predictors and responses for FoSR / SoFR / cosinor
      values$covariates <- values$raw_df[, cov_cols, drop = FALSE]
      # factor copies: grouping for fANOVA / clustering / cosinor group tests
      values$selected_group_vars <- cov_cols
      values$group_variables <- values$raw_df[, cov_cols, drop = FALSE]
      for(col in colnames(values$group_variables)) {
        values$group_variables[[col]] <- as.factor(values$group_variables[[col]])
      }
      values$group_labels <- values$group_variables[[1]]
      cat("Scalar variables selected:", paste(cov_cols, collapse = ", "), "\n")
      cat("Primary grouping variable:", cov_cols[1], "with",
          length(unique(values$group_labels)), "levels\n")
    } else {
      values$covariates <- NULL
      values$group_labels <- NULL
      values$group_variables <- NULL
      values$selected_group_vars <- NULL
    }

    values$data <- temp_data_mat

    # Time labels + numeric clock times, once, for every tab.
    values$time_labels <- colnames(temp_data)
    # WaPaa's plotting x coordinates (1:n_time — extract_time_values() does not
    # read the column names) ...
    values$time_numeric <- extract_time_values(values$time_labels)
    # ... and, separately, real clock hours when the names actually yield them.
    values$time_clock <- fck_clock_hours(values$time_labels)
    cat("Time labels stored:", paste(head(values$time_labels, 3), collapse = ", "), "...\n")
    if(!is.null(values$time_clock)) {
      cat("Clock times parsed:", paste(head(values$time_clock, 10), collapse = ", "), "...\n")
      if(fck_spacing_is_uneven(values$time_labels))
        cat("NOTE: these time points are NOT evenly spaced.\n")
    } else {
      cat("No clock times could be parsed from the column names.\n")
    }

    # Drop all-NA rows and columns (and keep the group vectors aligned).
    if(!is.null(values$data)) {
      na_rows <- apply(is.na(values$data), 1, all)
      if(any(na_rows)) {
        values$data <- values$data[!na_rows, , drop = FALSE]
        if(!is.null(values$group_labels))
          values$group_labels <- values$group_labels[!na_rows]
        if(!is.null(values$group_variables))
          values$group_variables <- values$group_variables[!na_rows, , drop = FALSE]
        if(!is.null(values$covariates))
          values$covariates <- values$covariates[!na_rows, , drop = FALSE]
      }
      # Rows with barely any measurements are mostly reconstructed by the
      # smoother rather than observed. Dropping them here, where the frame is
      # rebuilt from raw_df each time, keeps every parallel vector aligned and
      # keeps the choice reversible: lower the threshold and press Confirm again.
      min_obs <- input$min_observed_points
      if (!is.null(min_obs) && is.finite(min_obs) && min_obs > 0) {
        n_obs_row <- rowSums(!is.na(values$data))
        too_few <- n_obs_row < min_obs
        if (any(too_few)) {
          values$data <- values$data[!too_few, , drop = FALSE]
          if(!is.null(values$group_labels))
            values$group_labels <- values$group_labels[!too_few]
          if(!is.null(values$group_variables))
            values$group_variables <- values$group_variables[!too_few, , drop = FALSE]
          if(!is.null(values$covariates))
            values$covariates <- values$covariates[!too_few, , drop = FALSE]
          showNotification(
            sprintf("Dropped %d row%s with fewer than %d measured time points (kept %d).",
                    sum(too_few), if(sum(too_few) == 1) "" else "s",
                    as.integer(min_obs), nrow(values$data)),
            type = "warning", duration = 10)
        }
      }

      na_cols <- apply(is.na(values$data), 2, all)
      if(any(na_cols)) {
        values$data <- values$data[, !na_cols, drop = FALSE]
        values$time_labels <- values$time_labels[!na_cols]
        if(!is.null(values$time_numeric) && length(values$time_numeric) == length(na_cols))
          values$time_numeric <- values$time_numeric[!na_cols]
        values$time_clock <- fck_clock_hours(values$time_labels)
      }
    }

    fck_reset_analyses(values)
    showNotification("Data processed successfully!", type = "message")

  }, error = function(e) {
    showNotification(paste("Error processing selection:", e$message), type = "error")
  })
})

# --- 3. Sample data ----------------------------------------------------------
# One dataset has to exercise both families, so this is a 24-hour circadian set
# with group structure (fANOVA / clustering), scalar covariates (FoSR), a binary
# outcome (SoFR logistic) and a real rhythm with group differences in MESOR,
# amplitude and acrophase (cosinor).
observeEvent(input$generate_sample, {
  cat("Generate sample button clicked\n")

  tryCatch({
    set.seed(123)
    n_subjects <- 50
    n_time     <- 24
    hours      <- 0:23                       # clock time, one column per hour

    with_groups <- isTRUE(input$generate_with_groups)
    n_groups    <- if(with_groups) {
      if(is.null(input$n_groups)) 3 else input$n_groups
    } else 1

    grp_idx <- rep(seq_len(n_groups), length.out = n_subjects)
    grp     <- factor(paste0("Group", grp_idx))

    age   <- round(rnorm(n_subjects, 50, 10))
    sex   <- factor(sample(c("Male", "Female"), n_subjects, replace = TRUE),
                    levels = c("Female", "Male"))
    score <- rnorm(n_subjects, 100, 15)

    sample_data <- matrix(NA_real_, nrow = n_subjects, ncol = n_time)
    for(i in 1:n_subjects) {
      g <- grp_idx[i]
      mesor     <- 50 + (g - 1) * 4 + 0.1 * (age[i] - 50) + rnorm(1, 0, 2)
      amplitude <- 10 + (g - 1) * 2.5 + rnorm(1, 0, 1.5)
      acrophase <- 16 + (g - 1) * 1.5 + rnorm(1, 0, 0.8)   # peak hour
      sample_data[i, ] <-
        mesor +
        amplitude * cos(2 * pi * (hours - acrophase) / 24) +
        0.35 * amplitude * cos(2 * pi * 2 * (hours - acrophase) / 24) +
        rnorm(n_time, 0, 1.2)
    }

    # Binary outcome driven by the curve level, for the SoFR logistic branch.
    curve_integral <- rowMeans(sample_data)
    prob_outcome   <- plogis(-0.15 * (curve_integral - mean(curve_integral)) + 0.02 * (age - 50))
    binary_outcome <- rbinom(n_subjects, 1, prob_outcome)

    values$data <- sample_data
    if(typeof(values$data) != "double" && typeof(values$data) != "integer") {
      mode(values$data) <- "numeric"
    }

    values$time_labels  <- sprintf("%02d:00", hours)
    values$time_numeric <- seq_len(n_time)          # plotting x axis, as WaPaa
    values$time_clock   <- as.numeric(hours)        # real clock hours: 0..23

    covs <- data.frame(
      ID      = 1:n_subjects,
      Group   = grp,
      Age     = age,
      Sex     = sex,
      Score   = score,
      Outcome = binary_outcome,
      stringsAsFactors = FALSE
    )
    values$covariates <- covs

    if(with_groups) {
      values$selected_group_vars <- c("Group", "Sex")
      values$group_variables <- data.frame(Group = grp, Sex = sex)
      values$group_labels <- grp
    } else {
      values$selected_group_vars <- "Sex"
      values$group_variables <- data.frame(Sex = sex)
      values$group_labels <- NULL
    }

    values$raw_df <- NULL        # hide the selection UI: nothing to select
    values$uploaded_data <- NULL
    fck_reset_analyses(values)

    showNotification(
      sprintf("Sample data generated: %d subjects x 24 hourly time points%s, with Age/Sex/Score covariates and a binary Outcome.",
              n_subjects, if(with_groups) sprintf(", %d groups", n_groups) else ""),
      type = "message", duration = 8)

  }, error = function(e) {
    cat("Error in generate_sample:", e$message, "\n")
    showNotification(paste("Error generating data:", e$message), type = "error", duration = 10)
  })
})

# --- 4. Status and preview ---------------------------------------------------
output$data_status <- renderPrint({
  if(is.null(values$data)) {
    if(!is.null(values$raw_df)) {
      cat("Raw file loaded. Please select variables to proceed.\n")
    } else {
      cat("No data loaded.\n")
      cat("Click 'Generate Sample Data' or upload a file.\n")
    }
  } else {
    cat("Analysis Data Ready!\n")
    cat("Dimensions:", nrow(values$data), "subjects x", ncol(values$data), "time points\n")
    if(!is.null(values$time_clock)) {
      cat("Clock times parsed from the column names:",
          paste(head(values$time_clock, 8), collapse = ", "),
          if(length(values$time_clock) > 8) "..." else "", "\n")
      if(fck_spacing_is_uneven(values$time_labels)) {
        cat("  These time points are NOT evenly spaced. Smoothing treats every\n")
        cat("  column as one equal step unless you tick 'Space time points by\n")
        cat("  their real clock times' on the smoothing tab.\n")
      }
    } else {
      cat("No clock times could be parsed from the column names.\n")
      cat("  Analyses that need real time (harmonic regression, real-time\n")
      cat("  smoothing) fall back to the column order; harmonic regression can\n")
      cat("  also be given the times manually on its own tab.\n")
    }
    if(!is.null(values$group_labels)) {
      cat("Groups:", length(unique(values$group_labels)), "groups detected\n")
      cat("Group distribution:", table(values$group_labels), "\n")
    } else {
      cat("No grouping variable selected.\n")
    }
    if(!is.null(values$covariates)) {
      cat("Scalar variables:", paste(names(values$covariates), collapse = ", "), "\n")
    } else {
      cat("No scalar variables selected (FoSR/SoFR/cosinor group tests need at least one).\n")
    }
    if(is.null(values$smooth_data)) {
      cat("\nSmoothing: not applied yet — go to 'Data Preprocessing/Smoothing'.\n")
    } else {
      cat("\nSmoothing: applied. Every analysis tab will use the smoothed curves.\n")
    }
  }
})

# WaPaa's preview (row-count control, group column, T1..Tn headers), extended
# to show the scalar variables alongside — CIRCAREG's preview showed those and
# analysts need to see that the right rows line up with the right covariates.
output$data_preview <- renderDT({
  n_rows_display <- if(!is.null(input$data_preview_rows)) {
    as.integer(input$data_preview_rows)
  } else {
    10
  }

  if(is.null(values$data)) {
    if(!is.null(values$raw_df)) {
      n_rows <- if(n_rows_display == -1) nrow(values$raw_df) else min(n_rows_display, nrow(values$raw_df))
      datatable(values$raw_df[1:n_rows, ],
                options = list(pageLength = n_rows, scrollX = TRUE,
                               lengthMenu = c(5, 10, 20, 50, 100)),
                caption = paste("Raw Imported Data (", n_rows, " rows)"))
    } else {
      return(NULL)
    }
  } else {
    tryCatch({
      n_total_rows <- nrow(values$data)
      n_rows <- if(n_rows_display == -1) n_total_rows else min(n_rows_display, n_total_rows)

      n_cols_show <- min(15, ncol(values$data))
      preview_data <- as.data.frame(values$data[1:n_rows, 1:n_cols_show, drop = FALSE])
      colnames(preview_data) <- paste0("T", 1:ncol(preview_data))

      preview_rows <- 1:n_rows
      lead <- data.frame(Subject = preview_rows)
      if(!is.null(values$group_labels)) {
        lead <- cbind(Group = values$group_labels[preview_rows], lead)
      }
      if(!is.null(values$covariates) && nrow(values$covariates) >= n_rows) {
        lead <- cbind(lead, values$covariates[preview_rows, , drop = FALSE])
      }
      preview_data <- cbind(lead, preview_data)

      datatable(preview_data,
                options = list(pageLength = n_rows, scrollX = TRUE,
                               lengthMenu = c(5, 10, 20, 50, 100)),
                rownames = FALSE,
                caption = paste("Processed Analysis Data (", n_rows, "/", n_total_rows,
                                " rows, first ", n_cols_show, " time points)"))
    }, error = function(e) {
      return(NULL)
    })
  }
})
