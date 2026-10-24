# ============================================================================
# 03_data_module.R — file import + layout + preview (house rules 1, 2)
#
# Several files at once (one per patient in a CBT-I group, or one study file),
# in three layouts:
#   long           one row per night, one column per diary item
#   items_in_rows  diary items in rows, one column per day (paper-style
#                  CBT-I diaries; weeks stacked in a sheet or one per sheet)
#   wide_days      one row per person, day-numbered columns (QS_1, QS_2...)
# Excel, CSV/TXT, JSON (sleep-diary app exports), SPSS, RDS and RData.
# ============================================================================

LAYOUT_CHOICES <- c(
  "One row per night (long)" = "long",
  "Items in rows, days in columns (paper-style diary)" = "items_in_rows",
  "One row per person, day-numbered columns (wide)" = "wide_days"
)

dataModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 4,
      shiny::fileInput(ns("files"), "Diary file(s) \u2014 select several to pool patients",
        multiple = TRUE,
        accept = c(".xlsx", ".xls", ".xlsm", ".csv", ".txt", ".tsv", ".json",
                   ".sav", ".rds", ".RData", ".rda")),
      shiny::actionLink(ns("demo"), "or use the package's demo data (mcompd)"),
      shiny::hr(),
      shiny::radioButtons(ns("layout"), "Layout of each sheet / table",
                          LAYOUT_CHOICES, selected = "long"),
      shiny::uiOutput(ns("layout_hint")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'wide_days'", ns("layout")),
        shiny::radioButtons(ns("day_position"), "Day number position",
          c("After the item (QS_1, bedtijd_2)" = "suffix",
            "Before the item (D1_QS, dag2_bedtijd)" = "prefix"), inline = FALSE),
        shiny::textInput(ns("wide_pattern"),
          "Advanced: regex with two groups (blank = default)", "")
      ),
      shiny::selectInput(ns("sheets"), "Excel sheets",
        c("First sheet" = "first", "All sheets (stacked)" = "all",
          "Sheets whose name matches..." = "regex")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'regex'", ns("sheets")),
        shiny::textInput(ns("sheet_regex"), "Sheet-name pattern", "semaine|week|sem"),
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::checkboxInput(ns("skip_auto"), "Detect title rows", TRUE)),
        shiny::column(6, shiny::conditionalPanel(
          sprintf("!input['%s']", ns("skip_auto")),
          shiny::numericInput(ns("skip"), "Rows to skip", 0, min = 0, step = 1)))
      ),
      shiny::selectInput(ns("delim"), "Text-file delimiter",
        c("Detect" = "auto", "Semicolon ;" = ";", "Comma ," = ",", "Tab" = "\t")),
      shiny::uiOutput(ns("json_ui")),
      shiny::radioButtons(ns("id_from"), "Participant ID comes from",
        c("File name (one file per person)" = "file",
          "Sheet name (one sheet per person)" = "sheet",
          "A column in the data" = "column"), selected = "file"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'column'", ns("id_from")),
        shiny::selectInput(ns("id_column"), "ID column", choices = NULL)
      )
    ),
    shiny::mainPanel(
      width = 8,
      shiny::h4("Import log"),
      shiny::verbatimTextOutput(ns("log")),
      shiny::h4("Preview (one row per participant-night)"),
      shiny::helpText("Check that each row is one night and each column one diary",
                      "item before going on. Column names were made syntactic",
                      "(accents removed) so they can enter model formulas."),
      DT::DTOutput(ns("preview")),
      shiny::h4("Column names: original -> used"),
      DT::DTOutput(ns("labels"))
    )
  )
}

dataModuleServer <- function(id, rec) {
  shiny::moduleServer(id, function(input, output, session) {
    use_demo <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$demo, use_demo(TRUE))
    shiny::observeEvent(input$files, use_demo(FALSE))

    output$layout_hint <- shiny::renderUI({
      txt <- switch(input$layout,
        long = "e.g. a study file with columns Participant, Datum, Naar bed gegaan, ..., or a JSON app export (one record per night).",
        items_in_rows = "e.g. a CBT-I diary where rows are 'Heure de coucher', 'QS', 'QR', 'FJ' and columns are the days. Weeks stacked under each other are detected.",
        wide_days = "e.g. one row per participant with QS_1 ... QS_14, bedtijd_1 ... bedtijd_14.")
      shiny::helpText(txt)
    })

    sheets_arg <- shiny::reactive({
      if (identical(input$sheets, "regex")) input$sheet_regex %||% "first"
      else input$sheets
    })
    skip_arg <- shiny::reactive(if (isTRUE(input$skip_auto)) NA else input$skip)

    # JSON: offer every record table found in the first JSON file
    json_choices <- shiny::reactive({
      f <- input$files
      if (is.null(f)) return(NULL)
      js <- which(tolower(tools::file_ext(f$name)) == "json")
      if (!length(js)) return(NULL)
      doc <- tryCatch(jsonlite::fromJSON(f$datapath[js[1]], simplifyVector = TRUE,
                                         flatten = TRUE), error = function(e) NULL)
      if (is.null(doc)) return(NULL)
      tabs <- json_tables(doc)
      stats::setNames(names(tabs), sprintf("%s (%d records)", names(tabs),
                                           vapply(tabs, nrow, numeric(1))))
    })
    output$json_ui <- shiny::renderUI({
      ch <- json_choices()
      if (is.null(ch)) return(NULL)
      sizes <- as.numeric(gsub("^.*\\((\\d+) records\\)$", "\\1", names(ch)))
      shiny::selectInput(session$ns("json_table"), "JSON: table of diary entries",
                         ch, selected = ch[which.max(sizes)])
    })

    run_import <- function(id_from, id_column = NULL) {
      f <- input$files
      import_diary_files(
        paths = f$datapath, names = f$name, layout = input$layout,
        sheets = sheets_arg(), skip = skip_arg(), delim = input$delim,
        json_table = input$json_table, id_from = id_from, id_column = id_column,
        wide_pattern = if (nzchar(input$wide_pattern %||% "")) input$wide_pattern else NULL,
        day_position = input$day_position %||% "suffix")
    }

    # first pass (IDs from file names) only to list the columns for the ID picker
    base_import <- shiny::reactive({
      shiny::req(input$files, !use_demo())
      tryCatch(run_import("file"), error = function(e) e)
    })
    shiny::observe({
      b <- base_import()
      if (inherits(b, "error")) return()
      cols <- setdiff(names(b), c("ID", "night", ".file", ".sheet", ".row"))
      guess <- grep("^(id|pp|participant|subject|patient|code|naam|nom|name)",
                    cols, ignore.case = TRUE, value = TRUE)
      shiny::updateSelectInput(session, "id_column", choices = cols,
        selected = if (isTRUE(input$id_column %in% cols)) input$id_column
                   else guess[1] %||% cols[1])
    })

    imported <- shiny::reactive({
      if (use_demo()) {
        d <- as.data.frame(multilevelcoda::mcompd)
        d$ID <- as.character(d$ID)
        d$night <- d$Time
        attr(d, "log") <- sprintf("Demo data multilevelcoda::mcompd: %d participants, %d days",
                                  length(unique(d$ID)), nrow(d))
        attr(d, "labels") <- stats::setNames(names(d), names(d))
        return(d)
      }
      shiny::req(input$files)
      if (identical(input$id_from, "column")) shiny::req(input$id_column)
      out <- tryCatch(
        run_import(input$id_from,
                   if (identical(input$id_from, "column")) input$id_column else NULL),
        error = function(e) e)
      shiny::validate(shiny::need(!inherits(out, "error"),
                                  if (inherits(out, "error")) conditionMessage(out) else ""))
      out
    })

    # record the import step with every choice baked in (rule 6)
    shiny::observe({
      d <- imported()
      rec$clear_from("import")
      if (use_demo()) {
        rec$record("import", "Loaded the demo data multilevelcoda::mcompd (one row per person-day)",
          c("raw <- as.data.frame(multilevelcoda::mcompd)",
            "raw$ID <- as.character(raw$ID)",
            "raw$night <- raw$Time"))
        return()
      }
      f <- shiny::isolate(input$files)
      idc <- if (identical(input$id_from, "column")) input$id_column else NULL
      wp  <- if (nzchar(input$wide_pattern %||% "")) input$wide_pattern else NULL
      rec$record("import",
        sprintf("Imported %d file(s) (%s) as '%s' layout; IDs from %s; %d participants, %d rows",
                nrow(f), paste(f$name, collapse = ", "), input$layout, input$id_from,
                length(unique(d$ID)), nrow(d)),
        c('data_dir <- "."   # folder holding the files below',
          sprintf("files <- %s", r_lit(f$name)),
          "raw <- import_diary_files(",
          "  paths = file.path(data_dir, files), names = files,",
          sprintf("  layout = %s, sheets = %s, skip = %s, delim = %s,",
                  r_lit(input$layout), r_lit(sheets_arg()), r_lit(skip_arg()),
                  r_lit(input$delim)),
          sprintf("  json_table = %s, id_from = %s, id_column = %s,",
                  r_lit(input$json_table), r_lit(input$id_from), r_lit(idc)),
          sprintf("  wide_pattern = %s, day_position = %s)",
                  r_lit(wp), r_lit(input$day_position %||% "suffix")),
          "cat(attr(raw, \"log\"), sep = \"\\n\")"))
    }, priority = 100)

    output$log <- shiny::renderText(paste(attr(imported(), "log"), collapse = "\n"))
    output$preview <- DT::renderDT(
      DT::datatable(utils::head(imported(), 500), rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10)))
    output$labels <- DT::renderDT({
      lab <- attr(imported(), "labels")
      DT::datatable(data.frame(original = unname(lab), used = names(lab)),
                    rownames = FALSE, options = list(pageLength = 5, dom = "tp"))
    })

    list(data = imported, is_demo = use_demo)
  })
}
