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
# Extension-dispatch reader. Returns list(df, read_code, objects) where
# read_code is the plain-R fragment (explicit args baked in) for the
# recorder, and `objects` is only set for .RData (candidate object names).
# Always coerces to a plain data.frame; haven_labelled columns become factors
# so downstream estimators never see a labelled-vector class they don't expect.
read_any_file <- function(path, filename, opts) {
  ext <- tolower(tools::file_ext(filename))

  delabel <- function(df) {
    if (any(vapply(df, inherits, TRUE, what = "haven_labelled"))) {
      df <- as.data.frame(haven::as_factor(df))
    }
    as.data.frame(df)
  }

  switch(ext,
    csv = ,
    txt = {
      delim <- opts$delim %||% ","
      dec   <- opts$decimal %||% "."
      df <- as.data.frame(readr::read_delim(
        path, delim = delim,
        locale = readr::locale(decimal_mark = dec),
        guess_max = 10000, show_col_types = FALSE))
      list(df = delabel(df), read_code = sprintf(
        paste('dat_raw <- as.data.frame(readr::read_delim("%s", delim = "%s",',
              '  locale = readr::locale(decimal_mark = "%s"), guess_max = 10000))',
              sep = "\n"),
        filename, delim, dec))
    },
    tsv = {
      df <- as.data.frame(readr::read_tsv(path, guess_max = 10000,
                                          show_col_types = FALSE))
      list(df = delabel(df), read_code = sprintf(
        'dat_raw <- as.data.frame(readr::read_tsv("%s", guess_max = 10000))',
        filename))
    },
    xls = ,
    xlsx = {
      sheet <- opts$sheet %||% 1L
      df <- as.data.frame(readxl::read_excel(path, sheet = sheet))
      list(df = delabel(df), read_code = sprintf(
        'dat_raw <- as.data.frame(readxl::read_excel("%s", sheet = %s))',
        filename, sheet))
    },
    sav = {
      df <- haven::read_sav(path)
      list(df = delabel(df), read_code = sprintf(
        paste('dat_raw <- haven::read_sav("%s")',
              '# Labelled columns converted to factors (haven::as_factor):',
              'dat_raw <- as.data.frame(haven::as_factor(dat_raw))', sep = "\n"),
        filename))
    },
    rds = {
      df <- readRDS(path)
      if (!is.data.frame(df)) stop(".rds file must contain a data.frame.")
      list(df = delabel(df), read_code = sprintf(
        'dat_raw <- readRDS("%s")', filename))
    },
    rdata = {
      e <- new.env()
      load(path, envir = e)
      objs <- ls(e)
      df_objs <- objs[vapply(objs, function(o) is.data.frame(get(o, e)), TRUE)]
      if (!length(df_objs)) stop(".RData contains no data.frame objects.")
      chosen <- opts$rdata_object %||% df_objs[1]
      if (!(chosen %in% df_objs))
        stop(sprintf("Object '%s' not found (or not a data.frame) in .RData.", chosen))
      list(df = delabel(get(chosen, e)), objects = df_objs, read_code = sprintf(
        paste('.e <- new.env(); load("%s", envir = .e)',
              'dat_raw <- as.data.frame(get("%s", .e))', sep = "\n"),
        filename, chosen))
    },
    stop(sprintf(
      "Unsupported file type: .%s (supported: .csv .txt .tsv .xls .xlsx .sav .rds .RData)",
      ext))
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Per-module data transforms (request: z-score / centering / npn) --------
# Applied AFTER variable selection and BEFORE the analysis, separately per
# analysis tab. Only numeric columns are transformed. The nonparanormal is
# the rank-based normal-scores transform (Liu et al. 2009), ported from
# PsychoNetrix.R:930-937 (handles NAs; no extra dependency).
apply_house_transform <- function(dat, method = c("none", "center", "zscore", "npn")) {
  method <- match.arg(method)
  if (method == "none") return(dat)
  num <- vapply(dat, is.numeric, TRUE)
  m <- as.matrix(dat[, num, drop = FALSE])
  m_tr <- switch(method,
    center = scale(m, center = TRUE, scale = FALSE),
    zscore = scale(m, center = TRUE, scale = TRUE),
    npn = apply(m, 2, function(x) {
      n     <- sum(!is.na(x))
      delta <- 1 / (4 * n^0.25 * sqrt(pi * log(n)))
      r     <- rank(x, ties.method = "average", na.last = "keep")
      stats::qnorm(pmin(pmax(r / n, delta), 1 - delta))
    })
  )
  dat[, num] <- as.data.frame(m_tr)
  dat
}

TRANSFORM_LABELS <- c(
  "None"                                   = "none",
  "Mean-center"                            = "center",
  "Z-score (standardize)"                  = "zscore",
  "Nonparanormal (rank -> normal scores)"  = "npn")

# Shared UI element: one transform select per analysis tab.
transformUI <- function(ns_id) {
  shiny::selectInput(ns_id, "Data transformation (before analysis)",
                     TRANSFORM_LABELS, selected = "none")
}

# Script fragment reproducing apply_house_transform for the chosen method.
transform_code_fragment <- function(method, dat_var) {
  switch(method,
    none   = sprintf("# no data transformation applied to %s", dat_var),
    center = sprintf(
      "num <- vapply(%s, is.numeric, TRUE)\n%s[, num] <- as.data.frame(scale(as.matrix(%s[, num]), center = TRUE, scale = FALSE))",
      dat_var, dat_var, dat_var),
    zscore = sprintf(
      "num <- vapply(%s, is.numeric, TRUE)\n%s[, num] <- as.data.frame(scale(as.matrix(%s[, num])))",
      dat_var, dat_var, dat_var),
    npn = paste(
      sprintf("num <- vapply(%s, is.numeric, TRUE)", dat_var),
      sprintf("%s[, num] <- as.data.frame(apply(as.matrix(%s[, num]), 2, function(x) {", dat_var, dat_var),
      "  n     <- sum(!is.na(x))",
      "  delta <- 1 / (4 * n^0.25 * sqrt(pi * log(n)))",
      '  r     <- rank(x, ties.method = "average", na.last = "keep")',
      "  qnorm(pmin(pmax(r / n, delta), 1 - delta))",
      "}))  # nonparanormal (Liu et al. 2009)", sep = "\n"))
}

# JS helper: true when the uploaded file's name ends with any of `exts`.
.ext_condition <- function(ns_file, exts) {
  checks <- vapply(exts, function(e)
    sprintf("input['%s'].name.toLowerCase().endsWith('%s')", ns_file, e), "")
  sprintf("input['%s'] != null && (%s)", ns_file, paste(checks, collapse = " || "))
}

# --- Shiny module ------------------------------------------------------------
dataModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(ns("file"), "Upload data",
                     accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx",
                                ".sav", ".rds", ".RData")),
    shiny::conditionalPanel(
      .ext_condition(ns("file"), c(".csv", ".txt")),
      shiny::fluidRow(
        shiny::column(6, shiny::selectInput(ns("delim"), "Delimiter",
          c("Comma (,)" = ",", "Semicolon (;)" = ";",
            "Tab" = "\t", "Space" = " "))),
        shiny::column(6, shiny::selectInput(ns("decimal"), "Decimal mark",
          c("Period (.)" = ".", "Comma (,)" = ",")))
      )
    ),
    shiny::conditionalPanel(
      .ext_condition(ns("file"), c(".xls", ".xlsx")),
      shiny::numericInput(ns("sheet"), "Excel sheet number", value = 1, min = 1)
    ),
    shiny::conditionalPanel(
      .ext_condition(ns("file"), c(".rdata")),
      shiny::selectInput(ns("rdata_object"), "Object to use (.RData contains multiple)",
                         choices = NULL)
    ),
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
    ),
    shiny::hr(),
    shiny::h5("Zero-order correlation heatmap"),
    shinyWidgets::pickerInput(
      ns("heat_vars"), "Variables in the heatmap",
      choices = NULL, multiple = TRUE,
      options = shinyWidgets::pickerOptions(actionsBox = TRUE,
                                            liveSearch = TRUE)),
    shiny::plotOutput(ns("cor_heatmap"), height = "520px")
  )
}

dataModuleServer <- function(id, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    # .RData needs a peek before the "real" read so the object picker can be
    # populated with the candidate data.frame names.
    shiny::observeEvent(input$file, {
      ext <- tolower(tools::file_ext(input$file$name))
      if (identical(ext, "rdata")) {
        e <- new.env()
        load(input$file$datapath, envir = e)
        objs <- ls(e)
        df_objs <- objs[vapply(objs, function(o) is.data.frame(get(o, e)), TRUE)]
        shiny::updateSelectInput(session, "rdata_object", choices = df_objs)
      }
    })

    raw <- shiny::reactive({
      shiny::req(input$file)
      ext <- tolower(tools::file_ext(input$file$name))
      if (identical(ext, "rdata")) shiny::req(input$rdata_object)
      res <- read_any_file(input$file$datapath, input$file$name, opts = list(
        delim        = input$delim,
        decimal      = input$decimal,
        sheet        = input$sheet,
        rdata_object = input$rdata_object
      ))
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

    # --- Zero-order correlation heatmap over user-chosen variables ----------
    shiny::observeEvent(wide(), {
      num <- names(wide())[vapply(wide(), is.numeric, TRUE)]
      shinyWidgets::updatePickerInput(session, "heat_vars",
                                      choices = num, selected = num)
    })

    output$cor_heatmap <- shiny::renderPlot({
      vs <- input$heat_vars
      shiny::req(length(vs) >= 2)
      vs <- intersect(vs, names(wide()))
      cm <- stats::cor(wide()[, vs, drop = FALSE],
                       use = "pairwise.complete.obs")
      p  <- ncol(cm)
      cols <- grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(11, "RdBu")))(200)
      op <- graphics::par(mar = c(8, 8, 2, 2))
      on.exit(graphics::par(op), add = TRUE)
      graphics::image(seq_len(p), seq_len(p), t(cm[p:1, , drop = FALSE]),
                      zlim = c(-1, 1), col = cols, axes = FALSE,
                      xlab = "", ylab = "")
      graphics::axis(1, at = seq_len(p), labels = colnames(cm), las = 2,
                     cex.axis = 0.85)
      graphics::axis(2, at = seq_len(p), labels = rev(colnames(cm)), las = 2,
                     cex.axis = 0.85)
      for (i in seq_len(p)) for (j in seq_len(p))
        graphics::text(j, p - i + 1, sprintf("%.2f", cm[i, j]), cex = 0.75)
      graphics::box()
    })

    list(raw = raw, wide = wide, meta = meta)
  })
}
