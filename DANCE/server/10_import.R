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
#     for FoSR/cosinor predictors) and values$group_variables (factors,
#     for fANOVA/clustering/group comparisons).
#   * time values: WaPaa's extract_time_values() runs here for every dataset,
#     so values$time_numeric is available to all tabs (the cosinor tab can now
#     reuse it instead of re-detecting times from the same column names).
#   * participants and groups: neither app could leave anybody out short of
#     editing the file. Section 2b keeps the whole imported frame and applies a
#     keyed selection to it (server/02d_helpers_rowfilter.R), for uploaded and
#     sample data alike.
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

    # A new file voids everything downstream, in every family -- including the
    # participant/group selection, which named people in the OLD file.
    values$data <- NULL
    values$covariates <- NULL
    values$group_labels <- NULL
    values$group_variables <- NULL
    values$selected_group_vars <- NULL
    values$subject_ids <- NULL
    values$time_numeric <- NULL
    values$time_clock <- NULL
    values$import_full <- NULL
    values$row_filter <- NULL
    values$row_index <- NULL
    values$fill_status <- NULL
    dance_reset_analyses(values)

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
      "These serve <b>both</b> families: as predictors/response in FoSR and
       cosinor regression, and as grouping factors in functional ANOVA, cluster
       composition and cosinor group tests. The first one selected is the primary
       grouping variable; each tab can pick a different one.")),

    # Rows are CURVES, not necessarily people. Without an identifier no analysis
    # can tell that two rows are the same participant, and every between-groups
    # test then treats correlated observations as independent.
    selectInput(
      inputId = "subject_id_var",
      label = "Participant identifier (optional):",
      choices = c("Auto-detect / none" = "_none_", cols),
      selected = "_none_"
    ),
    helpText(HTML(
      "Only needed when one participant contributes <b>more than one row</b> —
       repeated days or sessions. The app cannot otherwise know, and a
       between-groups test on repeated rows is anticonservative. Auto-detects a
       column called ID, subject or participant.")),

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
      # original types: predictors and responses for FoSR / cosinor
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
    # which row of the raw frame each curve is, carried through every row drop
    # below: it keys a participant when the file has no identifier, and keeps
    # the modules that read the raw frame aligned (server/02d_helpers_rowfilter.R)
    row_index <- seq_len(nrow(temp_data_mat))

    # ---- participant identifier, if the file carries one --------------------
    # Scores and per-subject parameters are one row per CURVE. When the same
    # participant contributes several curves (repeated days, repeated sessions),
    # every between-groups test that treats rows as independent is
    # anticonservative. Nothing can detect that without an identifier, so one is
    # captured here and carried through the row filters below alongside the
    # group vectors. Auto-detected; the user can override in the scalar picker.
    values$subject_ids <- NULL
    id_pick <- input$subject_id_var
    if (!is.null(id_pick) && nzchar(id_pick) && !identical(id_pick, "_none_") &&
        id_pick %in% names(values$raw_df)) {
      values$subject_ids <- as.character(values$raw_df[[id_pick]])
    } else {
      cand <- names(values$raw_df)[tolower(names(values$raw_df)) %in%
                                     c("id", "subject", "subject_id", "subjectid",
                                       "participant", "participant_id", "pid")]
      if (length(cand)) {
        values$subject_ids <- as.character(values$raw_df[[cand[1]]])
        cat("Participant identifier auto-detected from column:", cand[1], "\n")
      }
    }
    if (!is.null(values$subject_ids) && length(values$subject_ids) != nrow(values$data)) {
      values$subject_ids <- NULL   # cannot vouch for it; better absent than wrong
    }
    if (!is.null(values$subject_ids)) {
      n_dup <- sum(duplicated(values$subject_ids))
      if (n_dup > 0)
        showNotification(
          sprintf("%d of %d rows repeat a participant identifier: %d distinct participants. Analyses that treat rows as independent observations will say so.",
                  n_dup, length(values$subject_ids),
                  length(unique(values$subject_ids))),
          type = "warning", duration = 15)
    }

    # Time labels + numeric clock times, once, for every tab.
    values$time_labels <- colnames(temp_data)
    # WaPaa's plotting x coordinates (1:n_time — extract_time_values() does not
    # read the column names) ...
    values$time_numeric <- extract_time_values(values$time_labels)
    # ... and, separately, real clock hours when the names actually yield them.
    values$time_clock <- dance_clock_hours(values$time_labels)
    cat("Time labels stored:", paste(head(values$time_labels, 3), collapse = ", "), "...\n")
    if(!is.null(values$time_clock)) {
      cat("Clock times parsed:", paste(head(values$time_clock, 10), collapse = ", "), "...\n")
      if(dance_spacing_is_uneven(values$time_labels))
        cat("NOTE: these time points are NOT evenly spaced.\n")
    } else {
      cat("No clock times could be parsed from the column names.\n")
    }

    # Drop all-NA rows and columns (and keep the group vectors aligned).
    if(!is.null(values$data)) {
      na_rows <- apply(is.na(values$data), 1, all)
      if(any(na_rows)) {
        values$data <- values$data[!na_rows, , drop = FALSE]
        row_index <- row_index[!na_rows]
        if(!is.null(values$subject_ids))
          values$subject_ids <- values$subject_ids[!na_rows]
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
          row_index <- row_index[!too_few]
          if(!is.null(values$subject_ids))
            values$subject_ids <- values$subject_ids[!too_few]
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
        values$time_clock <- dance_clock_hours(values$time_labels)
      }
    }

    # ---- which participants and groups (server/02d_helpers_rowfilter.R) -----
    # Everything above built the FULL frame. It is kept whole, so the selection
    # in box 3 can be changed or undone without reading the file again, and the
    # exclusions already made on it are re-applied BY PARTICIPANT AND LEVEL: a
    # re-confirm (another covariate, another time column) must not quietly bring
    # back the people who had been taken out, nor take out somebody else.
    values$import_full <- dance_import_frame(
      values$data, values$subject_ids, values$covariates, values$group_variables,
      values$group_labels, row_index, values$time_labels, values$time_numeric,
      values$time_clock)
    had_selection <- !dance_filter_is_empty(values$row_filter)
    sel <- .import_apply_rows()

    values$fill_status <- NULL
    dance_reset_analyses(values)
    showNotification(
      if (!had_selection) "Data processed successfully!"
      else if (isTRUE(sel$reset))
        paste("Data processed. The participant/group selection made earlier would leave",
              "fewer than two curves on this selection, so it was cleared: every curve is included.")
      else paste0("Data processed. Your participant/group selection was re-applied: ",
                  sel$lines[1],
                  if (sel$n_stale) sprintf(" %d exclusion%s no longer matched anything and %s dropped.",
                                           sel$n_stale, if (sel$n_stale == 1) "" else "s",
                                           if (sel$n_stale == 1) "was" else "were") else ""),
      type = if (isTRUE(sel$reset)) "warning" else "message",
      duration = if (had_selection) 12 else 5)

  }, error = function(e) {
    showNotification(paste("Error processing selection:", e$message), type = "error")
  })
})

# --- 2b. Participants & groups -----------------------------------------------
# The analysed rows are the imported frame (values$import_full) minus the
# exclusions in values$row_filter. This is the ONE place values$data and its
# parallel vectors are rebuilt from that pair; the confirm step, the sample
# generator and the two buttons below all come through here, so they cannot
# disagree about what a selection means.
.import_apply_rows <- function() {
  fr <- values$import_full
  if (is.null(fr)) return(NULL)
  rf <- values$row_filter
  fk <- dance_filter_keep(fr, rf$excluded_ids %||% character(0),
                          rf$excluded_levels %||% list())
  res <- dance_filter_apply(fr, fk$keep)
  reset <- FALSE
  if (!isTRUE(res$ok)) {
    # a stored selection that leaves nothing to analyse on a rebuilt frame is
    # dropped rather than applied; the caller says so
    rf <- NULL
    fk <- dance_filter_keep(fr)
    res <- dance_filter_apply(fr, fk$keep)
    reset <- TRUE
  }
  n_stale <- dance_filter_n_stale(fk)
  rf <- dance_filter_prune(rf, fk)
  values$row_filter      <- rf
  values$data            <- res$data
  values$subject_ids     <- res$subject_ids
  values$covariates      <- res$covariates
  values$group_variables <- res$group_variables
  values$group_labels    <- res$group_labels
  # only a record of real file rows is handed on (see dance_import_frame)
  values$row_index       <- if (isFALSE(fr$file_rows)) NULL else res$row_index
  values$time_labels     <- res$time_labels
  values$time_numeric    <- res$time_numeric
  values$time_clock      <- if (isTRUE(res$cols_dropped)) dance_clock_hours(res$time_labels)
                            else fr$time_clock
  lines <- dance_filter_describe(fr, fk, rf)
  cat("Participant selection:", paste(lines, collapse = "\n  "), "\n")
  list(fk = fk, rf = rf, reset = reset, n_stale = n_stale, lines = lines,
       cols_dropped = isTRUE(res$cols_dropped))
}

output$row_filter_ui <- renderUI({
  fr <- values$import_full
  if (is.null(fr))
    return(helpText(HTML(
      "Confirm a variable selection above (or generate the sample data) first.",
      "Every participant and every group is then listed here, all included; untick",
      "what should not be analysed and press <b>Apply selection</b>.")))
  rf <- values$row_filter
  ex_ids <- rf$excluded_ids %||% character(0)
  ex_lv  <- rf$excluded_levels %||% list()
  vars  <- dance_filter_vars(fr)
  parts <- dance_filter_participants(fr)
  has_ids <- !is.null(fr$subject_ids)

  level_ui <- if (length(vars)) {
    choices <- lapply(names(vars), function(nm) {
      v <- vars[[nm]]
      stats::setNames(dance_filter_code(nm, v$level), sprintf("%s  (n = %d)", v$level, v$n))
    })
    names(choices) <- names(vars)
    selected <- unlist(lapply(names(vars), function(nm) {
      lv <- vars[[nm]]$level
      dance_filter_code(nm, lv[!lv %in% ex_lv[[nm]]])
    }), use.names = FALSE)
    tagList(
      pickerInput("filter_levels", "Groups to include:",
                  choices = choices, selected = selected, multiple = TRUE,
                  options = list(`actions-box` = TRUE, `live-search` = TRUE,
                                 `selected-text-format` = "count > 4",
                                 `count-selected-text` = "{0} of {1} levels included",
                                 `none-selected-text` = "No level included")),
      helpText(HTML(
        "Listed per categorical scalar variable, with the number of curves at each",
        "level. Untick a level to leave out every curve that carries it;",
        "<i>(missing)</i> is a curve with no value for that variable. A curve is",
        "analysed only when all of its levels are ticked.")))
  } else {
    helpText(paste(
      "No categorical scalar variable is selected, so there are no groups to choose",
      "from. Select one under 2. Variable Selection to filter by group."))
  }

  part_ui <- tagList(
    pickerInput("filter_participants",
                if (has_ids) "Participants to include:" else "Curves (rows) to include:",
                choices = stats::setNames(parts$key, dance_filter_participant_text(parts)),
                selected = parts$key[!parts$key %in% ex_ids], multiple = TRUE,
                options = list(`actions-box` = TRUE, `live-search` = TRUE,
                               `selected-text-format` = "count > 4",
                               `count-selected-text` = "{0} of {1} included",
                               `none-selected-text` = "Nobody included")),
    helpText(if (has_ids)
      "Search by identifier. Unticking a participant removes every curve they contribute."
      else paste("The file has no participant identifier (set one under 2. Variable",
                 "Selection), so each curve is listed by its row in the file.")))

  fluidRow(
    column(5, level_ui),
    column(5, part_ui),
    column(2,
           actionButton("apply_row_filter", "Apply selection", class = "btn-success",
                        icon = icon("filter"), width = "100%"),
           br(), br(),
           actionButton("reset_row_filter", "Include everyone", icon = icon("undo"),
                        width = "100%"),
           br(), br(),
           helpText("Applying clears every result, smoothing included: each analysis",
                    "is re-run on the curves that are left."))
  )
})

# What is applied, and -- when the ticks differ from it -- what pressing Apply
# would give, so the consequence is visible before any result is cleared.
output$row_filter_preview <- renderText({
  fr <- values$import_full
  if (is.null(fr)) return("Nothing imported yet.")
  rf <- values$row_filter
  applied <- dance_filter_keep(fr, rf$excluded_ids %||% character(0),
                               rf$excluded_levels %||% list())
  out <- c("Applied:", paste0("  ", dance_filter_describe(fr, applied, rf)))
  if (!is.null(input$filter_participants) || !is.null(input$filter_levels)) {
    pend_rf <- dance_filter_from_ticks(fr, input$filter_levels, input$filter_participants,
                                       has_level_ui = length(dance_filter_vars(fr)) > 0)
    pend <- dance_filter_keep(fr, pend_rf$excluded_ids, pend_rf$excluded_levels)
    if (!identical(pend$keep, applied$keep)) {
      out <- c(out, "", "With the ticks above (press Apply selection to use them):",
               paste0("  ", sub(" analysed", " would be analysed",
                                dance_filter_describe(fr, pend, pend_rf))))
      if (pend$n_kept < 2)
        out <- c(out, "  ! fewer than two curves: this selection cannot be applied")
    }
  }
  paste(out, collapse = "\n")
})

observeEvent(input$apply_row_filter, {
  fr <- values$import_full
  if (is.null(fr)) {
    showNotification("Nothing to select from yet: confirm a variable selection first.",
                     type = "warning", duration = 6)
    return()
  }
  rf <- dance_filter_from_ticks(fr, input$filter_levels, input$filter_participants,
                                has_level_ui = length(dance_filter_vars(fr)) > 0)
  fk <- dance_filter_keep(fr, rf$excluded_ids, rf$excluded_levels)
  chk <- dance_filter_apply(fr, fk$keep)
  if (!isTRUE(chk$ok)) {
    showNotification(paste(chk$message, "Nothing was changed."), type = "error", duration = 10)
    return()
  }
  cur <- values$row_filter
  now <- dance_filter_keep(fr, cur$excluded_ids %||% character(0),
                           cur$excluded_levels %||% list())
  values$row_filter <- if (dance_filter_is_empty(rf)) NULL else rf
  if (identical(now$keep, fk$keep)) {
    # the same curves either way: record the ticks, keep every result
    showNotification("That selection analyses the same curves as before; nothing was cleared.",
                     type = "message", duration = 5)
    return()
  }
  sel <- .import_apply_rows()
  values$fill_status <- NULL
  dance_reset_analyses(values)
  showNotification(
    paste0("Selection applied: ", sel$lines[1],
           if (isTRUE(sel$cols_dropped))
             " Time points none of them was measured at were dropped." else "",
           " Every analysis was cleared; re-run smoothing and the analyses you need."),
    type = "message", duration = 10)
})

observeEvent(input$reset_row_filter, {
  fr <- values$import_full
  if (is.null(fr)) return()
  # put the ticks back as well, even when nothing was applied yet
  vars <- dance_filter_vars(fr)
  if (length(vars))
    updatePickerInput(session, "filter_levels", selected = unlist(lapply(names(vars),
      function(nm) dance_filter_code(nm, vars[[nm]]$level)), use.names = FALSE))
  updatePickerInput(session, "filter_participants",
                    selected = dance_filter_participants(fr)$key)
  if (dance_filter_is_empty(values$row_filter)) {
    showNotification("Every participant and group is already included.",
                     type = "message", duration = 4)
    return()
  }
  values$row_filter <- NULL
  .import_apply_rows()
  values$fill_status <- NULL
  dance_reset_analyses(values)
  showNotification(
    sprintf("Every participant and group is included again (%d curves). Every analysis was cleared; re-run what you need.",
            nrow(values$data)),
    type = "message", duration = 8)
})

# --- 3. Sample data ----------------------------------------------------------
# One dataset has to exercise both families, so this is a 24-hour circadian set
# with group structure (fANOVA / clustering), scalar covariates (FoSR), a binary
# covariate, and a real rhythm with group differences in MESOR, amplitude and
# acrophase (cosinor).
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

    # A binary covariate driven by the curve level. It was generated for the
    # scalar-on-function tab, which has since been removed; it is kept because
    # it is a perfectly good two-level grouping variable for the fANOVA,
    # clustering and cosinor group comparisons.
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
    values$subject_ids <- NULL   # a file's identifiers do not describe this set

    # participants and groups are selectable here too: a new dataset, so a new
    # frame and no exclusions carried over from whatever was loaded before
    values$row_filter <- NULL
    values$import_full <- dance_import_frame(
      values$data, NULL, values$covariates, values$group_variables,
      values$group_labels, seq_len(n_subjects), values$time_labels,
      values$time_numeric, values$time_clock)
    .import_apply_rows()
    values$fill_status <- NULL
    dance_reset_analyses(values)

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
    fr <- values$import_full
    if (!is.null(fr) && !dance_filter_is_empty(values$row_filter)) {
      fk <- dance_filter_keep(fr, values$row_filter$excluded_ids %||% character(0),
                              values$row_filter$excluded_levels %||% list())
      cat("Participant selection (box 3):\n")
      cat(paste0("  ", dance_filter_describe(fr, fk, values$row_filter)), sep = "\n")
    }
    if(!is.null(values$time_clock)) {
      cat("Clock times parsed from the column names:",
          paste(head(values$time_clock, 8), collapse = ", "),
          if(length(values$time_clock) > 8) "..." else "", "\n")
      if(dance_spacing_is_uneven(values$time_labels)) {
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
      cat("No scalar variables selected (FoSR/cosinor group tests need at least one).\n")
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
      # the row of the imported file, so an excluded participant visibly leaves
      # a gap rather than the rows below renumbering into its place
      ri <- values$row_index
      lead <- data.frame(Row = if (!is.null(ri) && length(ri) == n_total_rows)
                                 ri[preview_rows] else preview_rows)
      if(!is.null(values$subject_ids) && length(values$subject_ids) == n_total_rows) {
        lead$ID <- values$subject_ids[preview_rows]
      }
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
