# ============================================================================
# 05_diary_module.R — diary items -> sleep parameters -> composition
#
# Mapping is guessed from Dutch / French / English column names and can be
# overridden. Sleep parameters follow the Consensus Sleep Diary (Carney et
# al., 2012). Zeros are replaced multiplicatively (Martín-Fernández et al.,
# 2003) or the night is dropped; the choice and its count are logged.
# ============================================================================

BOOKKEEPING <- c("ID", "night", ".file", ".sheet", ".row")

diaryModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 4,
      shiny::radioButtons(ns("mode"), "Rows contain",
        c("Sleep-diary items (derive the sleep parameters)" = "diary",
          "Minutes per part already (e.g. accelerometry)" = "custom")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'diary'", ns("mode")),
        shiny::h5("Which column holds which diary item?"),
        shiny::uiOutput(ns("map_ui")),
        shiny::radioButtons(ns("clock_numeric"), "Clock times typed as numbers",
          c("23.30 means 23:30" = "hh.mm", "23.5 means 23:30" = "decimal_hours"),
          inline = FALSE),
        shiny::selectInput(ns("duration_unit"), "Durations (SOL, WASO, naps) are in",
          c("Detect (minutes, '1h30', '0:45', '10-15')" = "auto",
            "Minutes" = "minutes", "Hours" = "hours",
            "Excel time (fraction of a day)" = "day_fraction")),
        shiny::checkboxInput(ns("waso_blank_zero"),
          "Blank WASO = 0 min (no awakenings)", FALSE),
        shiny::selectInput(ns("scheme"), "Composition",
          stats::setNames(c("sleep24", "detail24", "tib4"),
                          COMPOSITION_SCHEMES[c("sleep24", "detail24", "tib4")])),
        shiny::checkboxInput(ns("naps"), "Count naps as sleep (removed from out-of-bed time)", FALSE)
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'custom'", ns("mode")),
        shinyWidgets::pickerInput(ns("custom_parts"), "Compositional parts (minutes)",
          choices = NULL, multiple = TRUE,
          options = list(`live-search` = TRUE, `actions-box` = TRUE)),
        shiny::numericInput(ns("custom_total"), "Total (e.g. 1440 min per day)", 1440, min = 1)
      ),
      shiny::radioButtons(ns("zeros"), "Parts equal to 0",
        c("Multiplicative replacement" = "multiplicative", "Drop the night" = "drop")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'multiplicative'", ns("zeros")),
        shiny::numericInput(ns("dl"), "Detection limit (min); zeros become 0.65 x this",
                            1, min = 0.01, step = 0.5)),
      shiny::hr(),
      shinyWidgets::pickerInput(ns("keep"), "Outcomes and covariates to keep",
        choices = NULL, multiple = TRUE,
        options = list(`live-search` = TRUE, `actions-box` = TRUE)),
      shinyWidgets::pickerInput(ns("lead"), "Take from the NEXT row (+1 night)",
        choices = NULL, multiple = TRUE),
      shinyWidgets::pickerInput(ns("lag"), "Take from the PREVIOUS row (-1 night)",
        choices = NULL, multiple = TRUE),
      shiny::helpText("Align each rating with the night it should be modelled",
        "against. In most morning diaries, sleep quality (QS) and quality of",
        "awakening (QR) describe the night on the same row. If the form of the",
        "day (FJ) is written down the evening AFTER that night but sits on the",
        "next row/column, choose it under NEXT row.")
    ),
    shiny::mainPanel(
      width = 8,
      shiny::tabsetPanel(
        shiny::tabPanel("Sleep parameters",
          shiny::verbatimTextOutput(ns("diary_log")),
          DT::DTOutput(ns("diary_table")),
          shiny::h4("Nights excluded"),
          DT::DTOutput(ns("invalid_table"))),
        shiny::tabPanel("Analysis data",
          shiny::verbatimTextOutput(ns("comp_log")),
          DT::DTOutput(ns("ana_table"))),
        shiny::tabPanel("Describe",
          shiny::h4("Compositional mean per participant (closed geometric mean)"),
          DT::DTOutput(ns("center_table")),
          partColoursUI(ns),
          shiny::sliderInput(ns("base_size"), "Text size", 8, 24, 13),
          shiny::plotOutput(ns("by_id"), height = "420px"),
          plotExportUI(ns, "by_id"),
          shiny::hr(),
          shiny::selectInput(ns("night_id"), "Night-by-night composition for", choices = NULL),
          shiny::plotOutput(ns("nights"), height = "380px"),
          plotExportUI(ns, "nights"),
          shiny::hr(),
          shiny::h4("Ternary (triangle) diagram"),
          shiny::helpText("Three parts at a time. With more than three parts the",
            "chosen three are re-closed (a subcomposition)."),
          shiny::fluidRow(
            shiny::column(6, shiny::selectizeInput(ns("tern_parts"),
              "Parts, in vertex order (bottom-left, bottom-right, top)",
              choices = NULL, multiple = TRUE,
              options = list(maxItems = 3, plugins = list("drag_drop", "remove_button")))),
            shiny::column(6, shiny::checkboxGroupInput(ns("tern_layers"), "Show",
              c("Nights" = "nights", "Participant means" = "persons",
                "Overall centre" = "centre"),
              selected = c("nights", "persons", "centre"), inline = TRUE))),
          shiny::fluidRow(
            shiny::column(4, shiny::radioButtons(ns("tern_ellipse"), "95% region of",
              c("Nights" = "nights", "Participant means" = "persons", "None" = "none"))),
            shiny::column(4,
              shiny::checkboxInput(ns("tern_centre"), "Centre the data", TRUE),
              shiny::helpText("Perturbs by the inverse compositional centre so",
                "small parts (e.g. awake in bed) become visible; gridlines keep",
                "the original percentages (von Eynatten et al., 2002).")),
            shiny::column(4,
              shiny::radioButtons(ns("tern_colour_by"), "Colour",
                c("By participant" = "participant", "Single colour" = "single")),
              colourpicker::colourInput(ns("tern_single"), "Single colour", "#B3CDE3"),
              colourpicker::colourInput(ns("tern_ell_col"), "Region colour", "#E5A0A0"))),
          shiny::fluidRow(
            shiny::column(6, shiny::sliderInput(ns("tern_size"), "Point size", 0.5, 5, 1.6, 0.1)),
            shiny::column(6, shiny::sliderInput(ns("tern_alpha"), "Point opacity", 0.05, 1, 0.45, 0.05))),
          shiny::plotOutput(ns("tern"), height = "620px"),
          plotExportUI(ns, "tern", w = 2400, h = 2100))
      )
    )
  )
}

diaryModuleServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {
    raw <- data_bus$data

    cols <- shiny::reactive(setdiff(names(raw()), BOOKKEEPING))
    num_cols <- shiny::reactive({
      d <- raw()
      setdiff(names(d)[vapply(d, is.numeric, logical(1))], BOOKKEEPING)
    })

    shiny::observeEvent(data_bus$is_demo(), {
      shiny::updateRadioButtons(session, "mode",
        selected = if (data_bus$is_demo()) "custom" else "diary")
    })

    output$map_ui <- shiny::renderUI({
      guess <- guess_diary_map(names(raw()))
      ch <- c("(not recorded)" = "", cols())
      lapply(names(DIARY_FIELDS), function(f)
        shiny::selectInput(session$ns(paste0("map_", f)), DIARY_FIELDS[[f]], ch,
                           selected = guess[[f]] %||% ""))
    })

    the_map <- shiny::reactive({
      m <- lapply(names(DIARY_FIELDS), function(f) input[[paste0("map_", f)]] %||% "")
      names(m) <- names(DIARY_FIELDS)
      m[nzchar(unlist(m))]
    })

    # outcome/covariate candidates: numeric columns not used as diary items
    shiny::observe({
      nc <- num_cols()
      mapped <- unlist(the_map())
      cand <- setdiff(nc, mapped)
      if (data_bus$is_demo()) cand <- setdiff(nc, c("TST", "WAKE", "MVPA", "LPA", "SB", "Time"))
      sel <- grep("^(QS|QR|FJ)$|qualit|kwalit|uitgerust|rested|forme|mood|humeur|stress|fatigue|vermoeid|slaperig|sleepiness|^Stress$|^Age$|^Female$",
                  cand, ignore.case = TRUE, value = TRUE)
      shinyWidgets::updatePickerInput(session, "keep", choices = cand, selected = sel)
      shinyWidgets::updatePickerInput(session, "custom_parts", choices = nc,
        selected = if (data_bus$is_demo()) c("TST", "WAKE", "MVPA", "LPA", "SB") else NULL)
    })
    shiny::observe({
      k <- input$keep
      shinyWidgets::updatePickerInput(session, "lead", choices = k,
                                      selected = intersect(shiny::isolate(input$lead), k))
      shinyWidgets::updatePickerInput(session, "lag", choices = k,
                                      selected = intersect(shiny::isolate(input$lag), k))
    })

    diary <- shiny::reactive({
      d <- raw()
      if (input$mode == "custom") return(d)
      shiny::req(length(the_map()) > 0)
      derive_sleep_diary(d, the_map(), clock_numeric = input$clock_numeric,
                         duration_unit = input$duration_unit,
                         waso_blank_zero = input$waso_blank_zero)
    })

    shiny::observe({
      d <- diary()
      if (input$mode == "custom") {
        rec$record("diary", "Rows already hold minutes per part; no diary derivation",
                   "diary <- raw")
        return()
      }
      m <- the_map()
      rec$record("diary",
        sprintf("Derived TIB, TST, SOL, WASO, EMA and SE per night (Consensus Sleep Diary definitions); %d of %d nights valid",
                sum(d$valid), nrow(d)),
        c(sprintf("diary_map <- list(%s)",
                  paste(sprintf("%s = %s", names(m), vapply(m, r_lit, "")), collapse = ",\n                  ")),
          sprintf("diary <- derive_sleep_diary(raw, diary_map, clock_numeric = %s,\n  duration_unit = %s, waso_blank_zero = %s)",
                  r_lit(input$clock_numeric), r_lit(input$duration_unit),
                  r_lit(input$waso_blank_zero)),
          "print(table(diary$invalid_reason[!diary$valid]))  # nights excluded and why"))
    }, priority = 90)

    shift <- shiny::reactive({
      s <- c(stats::setNames(rep(1L, length(input$lead)), input$lead),
             stats::setNames(rep(-1L, length(input$lag)), input$lag))
      s[!duplicated(names(s))]
    })

    scheme <- shiny::reactive(if (input$mode == "custom") "custom" else input$scheme)

    ana <- shiny::reactive({
      d <- diary()
      if (scheme() == "custom") shiny::req(length(input$custom_parts) >= 2)
      out <- tryCatch(build_composition(
        d, scheme = scheme(), naps_as_sleep = isTRUE(input$naps),
        custom_parts = input$custom_parts, custom_total = input$custom_total,
        zeros = input$zeros, dl = input$dl %||% 1,
        keep = input$keep %||% character(0), shift = shift()),
        error = function(e) e)
      shiny::validate(shiny::need(!inherits(out, "error"),
        if (inherits(out, "error")) conditionMessage(out) else ""))
      shiny::validate(shiny::need(nrow(out) >= 3, "Fewer than 3 usable nights."))
      out
    })
    parts <- shiny::reactive(attr(ana(), "parts"))
    total <- shiny::reactive(attr(ana(), "total"))

    shiny::observe({
      a <- ana()
      rec$record("composition",
        paste(attr(a, "log"), collapse = "; "),
        c(sprintf("ana <- build_composition(diary, scheme = %s, naps_as_sleep = %s,",
                  r_lit(scheme()), r_lit(isTRUE(input$naps))),
          sprintf("  custom_parts = %s, custom_total = %s,",
                  r_lit(if (scheme() == "custom") input$custom_parts else NULL),
                  r_lit(input$custom_total)),
          sprintf("  zeros = %s, dl = %s, keep = %s,",
                  r_lit(input$zeros), r_lit(input$dl %||% 1), r_lit(input$keep %||% character(0))),
          sprintf("  shift = %s)", if (length(shift())) r_lit(shift()) else "integer(0)"),
          "parts <- attr(ana, \"parts\"); total <- attr(ana, \"total\")",
          "cat(attr(ana, \"log\"), sep = \"\\n\")"))
    }, priority = 80)

    # ---- outputs ----
    output$diary_log <- shiny::renderText({
      d <- diary()
      if (input$mode == "custom") return("Parts are taken as recorded.")
      sprintf("%d nights, %d valid. Mean TIB %.0f min, mean TST %.0f min, mean SE %.1f%%.",
              nrow(d), sum(d$valid), mean(d$TIB[d$valid]), mean(d$TST[d$valid]),
              mean(d$SE[d$valid]))
    })
    output$diary_table <- DT::renderDT({
      d <- diary()
      show <- intersect(c("ID", "night", "date", "TIB", "TST", "SOL", "WASO", "EMA",
                          "PRE", "NAP", "SE", "valid"), names(d))
      DT::datatable(d[, show, drop = FALSE], rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10))
    })
    output$invalid_table <- DT::renderDT({
      d <- diary()
      shiny::req("valid" %in% names(d))
      bad <- d[!d$valid, , drop = FALSE]
      DT::datatable(bad[, c("ID", "night", "invalid_reason",
                            intersect(unlist(the_map()), names(bad))), drop = FALSE],
                    rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    })
    output$comp_log <- shiny::renderText(paste(attr(ana(), "log"), collapse = "\n"))
    output$ana_table <- DT::renderDT(
      DT::datatable(ana(), rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 10)) |>
        DT::formatRound(parts(), 1))

    output$center_table <- DT::renderDT({
      a <- ana(); p <- parts()
      ids <- unique(a$ID)
      tab <- t(vapply(ids, function(i) comp_center(a[a$ID == i, p, drop = FALSE], total()),
                      numeric(length(p))))
      tab <- rbind(tab, All = comp_center(a[, p, drop = FALSE], total()))
      df <- data.frame(ID = c(ids, "All participants"), nights = c(as.numeric(table(a$ID)[ids]), nrow(a)),
                       round(tab, 1), check.names = FALSE)
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    render_part_colours(output, session, parts)
    cols_r <- shiny::reactive(read_part_colours(input, parts()))

    by_id_plot <- shiny::reactive(
      plot_composition_by_id(ana(), parts(), total(), cols_r(), input$base_size))
    output$by_id <- shiny::renderPlot(by_id_plot())
    plotExportServer(input, output, "by_id", by_id_plot, "composition_by_participant")

    shiny::observe({
      ids <- unique(ana()$ID)
      shiny::updateSelectInput(session, "night_id", choices = ids,
        selected = if (isTRUE(shiny::isolate(input$night_id) %in% ids)) shiny::isolate(input$night_id) else ids[1])
    })
    nights_plot <- shiny::reactive({
      shiny::req(input$night_id %in% ana()$ID)
      plot_composition_nights(ana(), parts(), input$night_id, cols_r(), input$base_size)
    })
    output$nights <- shiny::renderPlot(nights_plot())
    plotExportServer(input, output, "nights", nights_plot, "composition_nights")

    # ---- ternary diagram ----
    shiny::observe({
      p <- parts()
      cur <- shiny::isolate(input$tern_parts)
      sel <- if (length(cur) == 3 && all(cur %in% p)) cur else utils::head(p, 3)
      shiny::updateSelectizeInput(session, "tern_parts", choices = p, selected = sel)
    })
    tern_args <- shiny::reactive({
      shiny::req(length(input$tern_parts) == 3, all(input$tern_parts %in% parts()))
      list(parts3 = input$tern_parts, layers = input$tern_layers %||% character(0),
           ellipse = input$tern_ellipse, centre_data = isTRUE(input$tern_centre),
           colour_by = input$tern_colour_by, single_colour = input$tern_single,
           ellipse_colour = input$tern_ell_col, point_size = input$tern_size,
           alpha = input$tern_alpha, base_size = input$base_size)
    })
    tern_plot <- shiny::reactive(do.call(ternary_plot, c(list(ana()), tern_args())))
    output$tern <- shiny::renderPlot(tern_plot())
    plotExportServer(input, output, "tern", tern_plot, "ternary")

    shiny::observe({
      a <- tern_args()
      rec$record("ternary",
        sprintf("Ternary diagram of %s%s (%s; 95%% region: %s)",
                paste(a$parts3, collapse = " / "),
                if (length(parts()) > 3) " (re-closed subcomposition)" else "",
                if (a$centre_data) "centred by the inverse compositional centre" else "not centred",
                a$ellipse),
        sprintf(paste0("print(ternary_plot(ana, parts3 = %s, layers = %s,\n",
                       "  ellipse = %s, centre_data = %s, colour_by = %s,\n",
                       "  single_colour = %s, ellipse_colour = %s, point_size = %s,\n",
                       "  alpha = %s, base_size = %s))"),
                r_lit(a$parts3), r_lit(a$layers), r_lit(a$ellipse), r_lit(a$centre_data),
                r_lit(a$colour_by), r_lit(a$single_colour), r_lit(a$ellipse_colour),
                a$point_size, a$alpha, a$base_size))
    }, priority = 65)

    shiny::observe({
      cl <- cols_r()
      rec$record("describe",
        "Described compositions: closed geometric means per participant and figures (pastel colours chosen in the app)",
        c(sprintf("part_colours <- %s", colours_literal(cl)),
          "centres <- t(sapply(split(ana[, parts], ana$ID), comp_center, total = total))",
          "print(round(centres, 1))",
          sprintf("print(plot_composition_by_id(ana, parts, total, part_colours, base_size = %s))",
                  input$base_size),
          if (!is.null(input$night_id))
            sprintf("print(plot_composition_nights(ana, parts, %s, part_colours, base_size = %s))",
                    r_lit(input$night_id), input$base_size)))
    }, priority = 70)

    list(ana = ana, parts = parts, total = total, diary = diary, colours = cols_r)
  })
}
