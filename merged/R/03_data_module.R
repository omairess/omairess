# ============================================================================
# 03_data_module.R — shared data ingestion + reshape (house rules 1, 2)
#
# THE DATA CONTRACT for all three analysis tabs:
#   data_bus$raw()   — data.frame exactly as parsed from the uploaded file
#   data_bus$wide()  — canonical ANALYSIS frame: one row per subject, one
#                      column per variable(-per-wave). All tabs consume ONLY
#                      this. If the upload was long, wide() is the pivoted
#                      result; if wide, wide() == raw().
#   data_bus$meta()  — list(filename, n, p, reshaped, id_col, time_col,
#                      sep = "_", missing_prop, missing_policy)
#
# DECISION: canonical format is WIDE, and wide column names use the
#   UNDERSCORE convention var_<time> (tidyr names_sep = "_"), because
#   Dagger's folded-temporal-graph auto-detection greps _t\d+$/_w\d+$/_\d+$
#   and would silently fail on base reshape()'s var.<time> dot naming
#   (PsychoNetrix's old behaviour, dropped).
# DECISION: BOTH raw and pivoted frames get a head+dims preview (rule 2's
#   "inspectable"): PsychoNetrix previously pivoted invisibly — fixed here.
# DECISION: a new upload drops every downstream recorder step (all tabs'
#   analysis/stability/comparison/plot fragments) — an exported script must
#   never mix fragments computed from different datasets.
# ============================================================================

# --- Reshape (rule 2), written in full: this defines the contract -----------
long_to_wide_canonical <- function(df, id_col, time_col, sep = "_") {
  stopifnot(id_col %in% names(df), time_col %in% names(df))
  value_vars <- setdiff(names(df), c(id_col, time_col))
  if (!length(value_vars)) stop("No measurement columns besides id and time.")
  if (length(unique(df[[time_col]])) < 2)
    stop("Time column must have at least 2 distinct values.")
  tidyr::pivot_wider(
    df,
    id_cols     = tidyr::all_of(id_col),
    names_from  = tidyr::all_of(time_col),
    values_from = tidyr::all_of(value_vars),
    names_sep   = sep
  )
}

# Code fragment mirroring long_to_wide_canonical for the exported script.
reshape_code_fragment <- function(id_col, time_col, sep = "_") {
  sprintf(
    paste(
      'dat_wide <- tidyr::pivot_wider(dat_raw,',
      '  id_cols     = "%s",',
      '  names_from  = "%s",',
      '  values_from = setdiff(names(dat_raw), c("%s", "%s")),',
      '  names_sep   = "%s")',
      'dat_wide <- as.data.frame(dat_wide)', sep = "\n"),
    id_col, time_col, id_col, time_col, sep)
}

# --- File loading (rule 1) ---------------------------------------------------
# TODO(stage3-loader): implement read_any_file(path, filename, opts) with
# extension dispatch:
#   .csv/.txt  -> readr::read_delim(delim = opts$delim, locale =
#                 readr::locale(decimal_mark = opts$decimal), guess_max = 10000)
#                 (user-set delimiter + decimal mark; European ; + , support)
#   .tsv       -> readr::read_tsv
#   .xls/.xlsx -> readxl::read_excel(sheet = opts$sheet)
#   .sav       -> haven::read_sav, keep labels, offer haven::as_factor
#   .rds       -> readRDS (must be a data.frame)
#   .RData     -> load() into new.env(); if >1 object, user picks which
# Must return list(df, read_code) where read_code is the plain-R fragment
# (with the explicit read args baked in) for the recorder. Always coerce to
# plain data.frame; convert haven_labelled via haven::as_factor.
read_any_file <- function(path, filename, opts) {
  stop("TODO(stage3-loader): extension-dispatch file reader not yet implemented")
}

# --- Shiny module ------------------------------------------------------------
dataModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(ns("file"), "Upload data",
                     accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx",
                                ".sav", ".rds", ".RData")),
    # TODO(stage3-loader-ui): delimiter + decimal-mark selects, Excel sheet
    # numericInput, .RData object picker (all conditionalPanels on extension).
    shiny::radioButtons(ns("layout"), "Data layout",
                        c("Wide (subjects x variables)" = "wide",
                          "Long (stacked / repeated measures)" = "long")),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'long'", ns("layout")),
      shiny::selectInput(ns("id_col"),   "Subject / ID column",   choices = NULL),
      shiny::selectInput(ns("time_col"), "Time / wave column",    choices = NULL)
    ),
    shiny::radioButtons(ns("missing_policy"), "Missing data policy",
                        c("Listwise (complete cases)" = "listwise",
                          "Pairwise (where supported)" = "pairwise",
                          "Keep as-is (FIML / estimator handles NA)" = "keep")),
    shiny::hr(),
    shiny::h5("Raw data preview (head + dims)"),
    shiny::verbatimTextOutput(ns("raw_dims")),
    DT::dataTableOutput(ns("raw_head")),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'long'", ns("layout")),
      shiny::h5("Pivoted (wide) preview — this is what every analysis uses"),
      shiny::verbatimTextOutput(ns("wide_dims")),
      DT::dataTableOutput(ns("wide_head"))
    )
  )
}

dataModuleServer <- function(id, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    raw <- shiny::reactive({
      shiny::req(input$file)
      res <- read_any_file(input$file$datapath, input$file$name,
                           opts = list())  # TODO(stage3-loader): pass UI opts
      # New dataset voids every downstream fragment, in every tab:
      rec_drop_prefix(rec, "bootnet_")
      rec_drop_prefix(rec, "dag_")
      rec_drop_prefix(rec, "psynet_")
      rec_drop_prefix(rec, "nct_")
      rec_drop_prefix(rec, "shared_reshape")
      rec_upsert(
        rec, "shared_data_load", "data",
        description = sprintf(
          "Loaded '%s' (%d rows x %d columns); %.1f%% missing values; missing-data policy: %s.",
          input$file$name, nrow(res$df), ncol(res$df),
          100 * mean(is.na(res$df)), input$missing_policy),
        code = res$read_code
      )
      res$df
    })

    # Populate id/time pickers when data arrives
    shiny::observeEvent(raw(), {
      shiny::updateSelectInput(session, "id_col",   choices = names(raw()))
      shiny::updateSelectInput(session, "time_col", choices = names(raw()))
    })

    wide <- shiny::reactive({
      df <- raw()
      if (input$layout == "long") {
        shiny::req(input$id_col, input$time_col)
        out <- as.data.frame(
          long_to_wide_canonical(df, input$id_col, input$time_col, sep = "_"))
        rec_upsert(
          rec, "shared_reshape", "reshape",
          description = sprintf(
            "Reshaped long -> wide by id '%s' and time '%s' (columns named var_<time>).",
            input$id_col, input$time_col),
          code = reshape_code_fragment(input$id_col, input$time_col, "_")
        )
        out
      } else {
        rec_upsert(
          rec, "shared_reshape", "reshape",
          description = "Data already wide; used as-is.",
          code = "dat_wide <- dat_raw"
        )
        df
      }
    })

    meta <- shiny::reactive({
      df <- raw()
      list(
        filename = input$file$name,
        n = nrow(wide()), p = ncol(wide()),
        reshaped = identical(input$layout, "long"),
        id_col = if (identical(input$layout, "long")) input$id_col else NULL,
        time_col = if (identical(input$layout, "long")) input$time_col else NULL,
        sep = "_",
        missing_prop = mean(is.na(df)),
        missing_policy = input$missing_policy
      )
    })

    output$raw_dims  <- shiny::renderText({
      df <- raw()
      sprintf("%d rows x %d columns | %.1f%% missing",
              nrow(df), ncol(df), 100 * mean(is.na(df)))
    })
    output$raw_head  <- DT::renderDataTable({
      DT::datatable(utils::head(raw(), 10),
                    options = list(scrollX = TRUE, dom = "t"))
    })
    output$wide_dims <- shiny::renderText({
      df <- wide()
      sprintf("%d rows x %d columns after pivot", nrow(df), ncol(df))
    })
    output$wide_head <- DT::renderDataTable({
      DT::datatable(utils::head(wide(), 10),
                    options = list(scrollX = TRUE, dom = "t"))
    })

    list(raw = raw, wide = wide, meta = meta)
  })
}
