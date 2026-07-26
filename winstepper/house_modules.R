# =============================================================================
# house_modules.R  —  Reusable Shiny building blocks for the "shiny-house-style"
# skill. Source this from any app: source("house_modules.R").
#
# Implements the ten house rules:
#   0 ensure_packages / pkg_versions   (install+load, record versions; no auto-update)
#   1 mod_data_*                        (load & preview txt/csv/xls/xlsx/sav/rds/RData)
#   2 mod_reshape_*                     (wide<->long)
#   3 mod_varselect_*                   (multi-select with search)
#   4/5/9 mod_network_viz_*             (colour pickers, size sliders, export, pastels)
#   6/8 new_recorder / build_script / build_pipeline_log  (R export + pipeline log)
#   7 mod_exports_*                     (export data, results, script, log)
#
# Verified against package APIs on 2026-07-05. See SKILL.md "Provenance".
# =============================================================================

# ---- Rule 0: packages -------------------------------------------------------

#' Ensure packages are installed and loaded. Does NOT update by default.
#' @param pkgs character vector of package names
#' @param update if TRUE, run update.packages() for pkgs (opt-in only)
ensure_packages <- function(pkgs, update = FALSE) {
  installed <- rownames(installed.packages())
  missing   <- setdiff(pkgs, installed)
  if (length(missing)) {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing)
  }
  if (isTRUE(update)) {
    message("Updating (opt-in): ", paste(pkgs, collapse = ", "))
    update.packages(oldPkgs = pkgs, ask = FALSE)
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(require(p, character.only = TRUE))))
}

#' Named character vector of installed versions (for the reproducible header).
pkg_versions <- function(pkgs) {
  vapply(pkgs, function(p)
    tryCatch(as.character(packageVersion(p)),
             error = function(e) "NOT INSTALLED"),
    character(1))
}

# ---- Rules 6 + 8: the step recorder ----------------------------------------

#' Create a recorder. Each stage upserts a step by id, supplying a plain-language
#' description (pipeline log) and a code fragment (R export). Upsert-by-id keeps
#' the log and script in sync with the app's current state, with no duplicates.
new_recorder <- function() {
  steps <- shiny::reactiveVal(list())
  list(
    record = function(id, description, code) {
      # isolate(): record() must not take a reactive dependency on the value it
      # writes, otherwise calling it from a reactive context self-invalidates.
      cur <- shiny::isolate(steps())
      cur[[id]] <- list(id = id, description = description,
                        code = paste(code, collapse = "\n"),
                        time = Sys.time())
      steps(cur)
    },
    clear = function(id) { cur <- steps(); cur[[id]] <- NULL; steps(cur) },
    reset = function() steps(list()),
    steps = steps
  )
}

#' Build a standalone, reproducible R script from recorded steps.
build_script <- function(steps_list, pkgs) {
  vers <- pkg_versions(pkgs)
  header <- c(
    "# ---------------------------------------------------------------------",
    "# Auto-generated reproducible script (shiny-house-style)",
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "# Package versions at generation time:",
    paste0("#   ", format(names(vers), width = 14), vers),
    "# ---------------------------------------------------------------------",
    "",
    paste0("library(", pkgs, ")"),
    ""
  )
  body <- unlist(lapply(steps_list, function(s)
    c(paste0("## ", s$description), s$code, "")), use.names = FALSE)
  paste(c(header, body), collapse = "\n")
}

#' Build a plain-language, ordered, timestamped pipeline log (Markdown).
build_pipeline_log <- function(steps_list) {
  if (!length(steps_list))
    return("_No steps recorded yet._")
  lines <- vapply(seq_along(steps_list), function(i) {
    s <- steps_list[[i]]
    sprintf("%d. **%s** _(%s)_", i, s$description,
            format(s$time, "%Y-%m-%d %H:%M:%S"))
  }, character(1))
  paste(c("# Analysis pipeline", "", lines), collapse = "\n")
}

# ---- Rule 9: pastel palette -------------------------------------------------

#' Pastel palette of length n (interpolates beyond RColorBrewer maxima).
pastel_palette <- function(n, which = c("Pastel1", "Pastel2")) {
  which <- match.arg(which)
  base_n <- RColorBrewer::brewer.pal.info[which, "maxcolors"]
  base   <- RColorBrewer::brewer.pal(base_n, which)
  if (n <= base_n) base[seq_len(max(n, 3))][seq_len(n)]
  else grDevices::colorRampPalette(base)(n)
}

# ---- Rule 1: data load + preview -------------------------------------------

mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(ns("file"), "Data file",
      accept = c(".txt", ".csv", ".tsv", ".xls", ".xlsx",
                 ".sav", ".rds", ".RData", ".rdata")),
    shiny::conditionalPanel(
      sprintf("output['%s'] == 'delim'", ns("kind")),
      shiny::fluidRow(
        shiny::column(6, shiny::textInput(ns("sep"), "Delimiter", ",")),
        shiny::column(6, shiny::textInput(ns("dec"), "Decimal mark", ".")))),
    shiny::uiOutput(ns("robj_ui")),
    shiny::helpText("Confirm the preview and dimensions before analysing."),
    shiny::verbatimTextOutput(ns("dims")),
    DT::DTOutput(ns("preview"))
  )
}

#' @return reactive data.frame (records a `data_load` step on the recorder)
mod_data_server <- function(id, recorder) {
  shiny::moduleServer(id, function(input, output, session) {
    ext <- shiny::reactive({
      shiny::req(input$file)
      tolower(tools::file_ext(input$file$name))
    })
    output$kind <- shiny::reactive(
      if (isTRUE(ext() %in% c("txt", "csv", "tsv"))) "delim" else "other")
    shiny::outputOptions(output, "kind", suspendWhenHidden = FALSE)

    rdata_objs <- shiny::reactiveVal(NULL)
    output$robj_ui <- shiny::renderUI({
      objs <- rdata_objs()
      if (is.null(objs)) return(NULL)
      shiny::selectInput(session$ns("robj"), "Object in .RData", choices = objs)
    })

    data <- shiny::reactive({
      shiny::req(input$file)
      path <- input$file$datapath
      switch(ext(),
        "csv" = , "tsv" = , "txt" =
          readr::read_delim(path, delim = input$sep %||% ",",
            locale = readr::locale(decimal_mark = input$dec %||% "."),
            guess_max = 10000, show_col_types = FALSE),
        "xls" = , "xlsx" = readxl::read_excel(path),
        "sav" = haven::read_sav(path),
        "rds" = readRDS(path),
        "rdata" = , "RData" = {
          e <- new.env(); load(path, envir = e)
          rdata_objs(ls(e))
          shiny::req(input$robj); get(input$robj, envir = e)
        },
        stop("Unsupported file type: ", ext()))
    })

    output$dims <- shiny::renderText({
      d <- data(); sprintf("%d rows x %d columns", nrow(d), ncol(d)) })
    output$preview <- DT::renderDT(
      utils::head(data(), 50), options = list(scrollX = TRUE))

    shiny::observeEvent(data(), {
      recorder$record("data_load", "Loaded and previewed input data",
        c(sprintf('# Load the same file you used in the app (%s)',
                  input$file$name),
          switch(ext(),
            "csv" = , "tsv" = , "txt" = sprintf(
              'data <- readr::read_delim("%s", delim = "%s",\n  locale = readr::locale(decimal_mark = "%s"), show_col_types = FALSE)',
              input$file$name, input$sep %||% ",", input$dec %||% "."),
            "xls" = , "xlsx" = sprintf('data <- readxl::read_excel("%s")', input$file$name),
            "sav" = sprintf('data <- haven::read_sav("%s")', input$file$name),
            "rds" = sprintf('data <- readRDS("%s")', input$file$name),
            sprintf('load("%s"); data <- %s', input$file$name, input$robj %||% "OBJECT"))))
    })
    data
  })
}

# ---- Rule 2: wide <-> long --------------------------------------------------

mod_reshape_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::radioButtons(ns("dir"), "Reshape",
      c("None" = "none", "Wide -> Long" = "long", "Long -> Wide" = "wide"),
      inline = TRUE),
    shiny::uiOutput(ns("controls")))
}

mod_reshape_server <- function(id, data, recorder) {
  shiny::moduleServer(id, function(input, output, session) {
    output$controls <- shiny::renderUI({
      ns <- session$ns; cols <- names(data())
      if (input$dir == "long")
        shinyWidgets::pickerInput(ns("cols"), "Columns to pivot to long",
          choices = cols, multiple = TRUE,
          options = list(`live-search` = TRUE, `actions-box` = TRUE))
      else if (input$dir == "wide")
        shiny::tagList(
          shiny::selectInput(ns("names_from"), "names_from", cols),
          shiny::selectInput(ns("values_from"), "values_from", cols))
    })

    out <- shiny::reactive({
      d <- data()
      if (input$dir == "long" && length(input$cols))
        tidyr::pivot_longer(d, cols = tidyselect::all_of(input$cols),
                            names_to = "name", values_to = "value")
      else if (input$dir == "wide")
        tidyr::pivot_wider(d, names_from = input$names_from,
                           values_from = input$values_from)
      else d
    })

    shiny::observeEvent(out(), {
      code <- if (input$dir == "long" && length(input$cols))
        sprintf('data <- tidyr::pivot_longer(data, cols = c(%s),\n  names_to = "name", values_to = "value")',
                paste(sprintf('"%s"', input$cols), collapse = ", "))
      else if (input$dir == "wide")
        sprintf('data <- tidyr::pivot_wider(data, names_from = %s, values_from = %s)',
                input$names_from, input$values_from)
      else "# (no reshape applied)"
      recorder$record("reshape",
        sprintf("Reshaped data (%s)", input$dir), code)
    })
    out
  })
}

# ---- Rule 3: multi variable selection --------------------------------------

mod_varselect_ui <- function(id, label = "Variables (node set)") {
  ns <- shiny::NS(id)
  shinyWidgets::pickerInput(ns("vars"), label, choices = NULL, multiple = TRUE,
    options = list(`live-search` = TRUE, `actions-box` = TRUE,
                   `selected-text-format` = "count > 3"))
}

#' @return reactive character vector of selected variable names
mod_varselect_server <- function(id, data, recorder, numeric_only = TRUE) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(data(), {
      d <- data()
      choices <- if (numeric_only) names(d)[vapply(d, is.numeric, logical(1))]
                 else names(d)
      shinyWidgets::updatePickerInput(session, "vars", choices = choices)
    })
    shiny::observeEvent(input$vars, {
      recorder$record("varselect", "Selected analysis variables",
        sprintf('vars <- c(%s)\ndata <- data[, vars, drop = FALSE]',
                paste(sprintf('"%s"', input$vars), collapse = ", ")))
    })
    shiny::reactive(input$vars)
  })
}

# ---- Rules 4 + 5 + 9: network visualisation --------------------------------

mod_network_viz_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(3, colourpicker::colourInput(ns("node"),  "Node fill",  "#FBB4AE")),
      shiny::column(3, colourpicker::colourInput(ns("border"),"Node border","#7A7A7A")),
      shiny::column(3, colourpicker::colourInput(ns("pos"),   "Positive edge","#B3CDE3")),
      shiny::column(3, colourpicker::colourInput(ns("neg"),   "Negative edge","#FDDAEC"))),
    shiny::fluidRow(
      shiny::column(4, shiny::sliderInput(ns("vsize"), "Node size", 2, 20, 8)),
      shiny::column(4, shiny::sliderInput(ns("esize"), "Edge scaling", 1, 25, 12)),
      shiny::column(4, shiny::sliderInput(ns("lcex"),  "Label size", 0.4, 3, 1, 0.1))),
    shiny::plotOutput(ns("plot"), height = "560px"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("w"), "Export width (px)", 1600)),
      shiny::column(3, shiny::numericInput(ns("h"), "Export height (px)", 1600)),
      shiny::column(3, shiny::numericInput(ns("res"), "Export res (dpi)", 300)),
      shiny::column(3, shiny::radioButtons(ns("fmt"), "Format",
                        c("png","pdf","svg"), inline = TRUE))),
    shiny::downloadButton(ns("dl"), "Download figure"),
    shiny::helpText("Clipboard copy isn't reliable from Shiny; use Download.")
  )
}

#' @param net_input reactive supplying a bootnet/qgraph-compatible object
#'   (a weights matrix, or an object with $graph). Records a `plot` step.
mod_network_viz_server <- function(id, net_input, recorder) {
  shiny::moduleServer(id, function(input, output, session) {
    weights <- shiny::reactive({
      x <- net_input()
      if (is.list(x) && !is.null(x$graph)) x$graph else as.matrix(x)
    })

    graph_obj <- shiny::reactive({
      shiny::req(weights())
      qgraph::qgraph(weights(), DoNotPlot = TRUE,
        color = input$node, border.color = input$border,
        posCol = input$pos, negCol = input$neg,
        vsize = input$vsize, esize = input$esize, label.cex = input$lcex,
        layout = "spring")
    })

    output$plot <- shiny::renderPlot({ plot(graph_obj()) })

    output$dl <- shiny::downloadHandler(
      filename = function() paste0("network.", input$fmt),
      content = function(file) {
        switch(input$fmt,
          png = grDevices::png(file, width = input$w, height = input$h, res = input$res),
          pdf = grDevices::pdf(file, width = input$w/300, height = input$h/300),
          svg = grDevices::svg(file, width = input$w/300, height = input$h/300))
        plot(graph_obj()); grDevices::dev.off()
      })

    shiny::observeEvent(graph_obj(), {
      recorder$record("plot", "Rendered network figure (house style)",
        sprintf(paste0(
          'qgraph::qgraph(weights,\n',
          '  color = "%s", border.color = "%s",\n',
          '  posCol = "%s", negCol = "%s",\n',
          '  vsize = %s, esize = %s, label.cex = %s, layout = "spring")'),
          input$node, input$border, input$pos, input$neg,
          input$vsize, input$esize, input$lcex))
    })
    graph_obj
  })
}

# ---- Rules 6 + 7 + 8: exports ----------------------------------------------

mod_exports_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h4("Exports"),
    shiny::downloadButton(ns("script"), "Download reproducible R script"),
    shiny::downloadButton(ns("data"),   "Download analysis data (.rds)"),
    shiny::downloadButton(ns("results"),"Download results (.rds)"),
    shiny::downloadButton(ns("log"),    "Download pipeline log (.md)"),
    shiny::h4("Pipeline"),
    shiny::uiOutput(ns("log_view")),
    shiny::h4("Reproducible R script"),
    shiny::verbatimTextOutput(ns("script_view")))
}

mod_exports_server <- function(id, recorder, pkgs, data, results) {
  shiny::moduleServer(id, function(input, output, session) {
    script <- shiny::reactive(build_script(recorder$steps(), pkgs))
    logmd  <- shiny::reactive(build_pipeline_log(recorder$steps()))

    output$script_view <- shiny::renderText(script())
    output$log_view    <- shiny::renderUI(
      shiny::pre(shiny::HTML(logmd())))

    output$script <- shiny::downloadHandler(
      function() "analysis_script.R", function(f) writeLines(script(), f))
    output$log <- shiny::downloadHandler(
      function() "pipeline_log.md", function(f) writeLines(logmd(), f))
    output$data <- shiny::downloadHandler(
      function() "analysis_data.rds", function(f) saveRDS(data(), f))
    output$results <- shiny::downloadHandler(
      function() "analysis_results.rds",
      function(f) saveRDS(if (is.function(results)) results() else results, f))
  })
}

# small null-coalescing helper
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
