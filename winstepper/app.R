# =============================================================================
# app.R  --  WINSTEPPER: a fresh Shiny front end for WINSTEPS-style Rasch
#            measurement, built on the shared rasch_engine.R.
#
# This is a redesigned UI (bslib: grouped nav menus, a value-box dashboard,
# cards) that drives the *same* estimation engine, plotting helpers and
# house-style building blocks as R-Winsteps. The engine is deliberately
# untouched so the numbers remain the audited ones.
#
# Built to the shiny-house-style rules:
#   0 packages ensured + versions recorded    1 multi-format data load + preview
#   2 explicit wide<->long reshape            3 searchable multi-select variables
#   4 colour pickers for every figure         5 size sliders + PNG/PDF/SVG export
#   6 auto-generated reproducible R script    7 export of data and results
#   8 plain-language pipeline log             9 pastel palettes by default
#
# Run with:  shiny::runApp("winstepper")   (or  shiny::runApp("app.R")  here)
# Requires rasch_engine.R, winsteps_plots.R and house_modules.R alongside.
# =============================================================================

source("house_modules.R")
source("rasch_engine.R")
source("winsteps_plots.R")
source("winstepper_extras.R")

PKGS <- c("shiny", "bslib", "shinyWidgets", "colourpicker", "DT", "readr",
          "readxl", "haven", "tidyr", "tidyselect", "RColorBrewer")
try(ensure_packages(PKGS), silent = TRUE)
invisible(lapply(PKGS, function(p)
  suppressPackageStartupMessages(try(require(p, character.only = TRUE), silent = TRUE))))

APP_VERSION <- "1.0.0"

# ---------------------------------------------------------------------------
# Small helpers (shared with R-Winsteps; kept identical for parity)
# ---------------------------------------------------------------------------

#' Round only the numeric columns of a data.frame (DT-safe).
round_num <- function(d, digits = 2) {
  d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, digits) else v)
  d
}

#' Apply WINSTEPS user scaling (UIMEAN= / USCALE=) to measure-like columns.
rescale_cols <- function(d, umean = 0, uscale = 1,
                         meas = c("Measure", "Andrich_Threshold", "Category_Measure",
                                  "Thurstone_50pct", "Obsvd_Avrge", "Sample_Expect",
                                  "DIF_Measure_A", "DIF_Measure_B"),
                         sds = c("Model_SE", "Real_SE", "SE", "Threshold_SE",
                                 "DIF_SE_A", "DIF_SE_B", "Joint_SE", "DIF_Contrast")) {
  if (umean == 0 && uscale == 1) return(d)
  for (cn in intersect(meas, names(d))) d[[cn]] <- umean + uscale * d[[cn]]
  for (cn in intersect(sds,  names(d))) d[[cn]] <- uscale * d[[cn]]
  d
}

dt <- function(d, digits = 2, scrollY = "58vh", ...) {
  DT::datatable(round_num(as.data.frame(d), digits), rownames = FALSE,
                extensions = "Buttons", filter = "top",
                options = list(pageLength = 25, scrollX = TRUE, scrollY = scrollY,
                               dom = "Blfrtip", buttons = c("copy", "csv", "excel"),
                               lengthMenu = c(10, 25, 50, 100, 500)), ...)
}

#' Parse a WINSTEPS-style CODES= string into a numeric vector of valid codes.
#' Accepts comma-, space- or semicolon-separated values, e.g. "0,1,2,3" or
#' "0 1 2 3". Returns NULL when blank (meaning: keep every observed value).
parse_codes <- function(txt) {
  if (is.null(txt)) return(NULL)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  parts <- strsplit(txt, "[[:space:],;]+")[[1]]
  parts <- parts[nzchar(parts)]
  v <- suppressWarnings(as.numeric(parts))
  v <- v[!is.na(v)]
  if (!length(v)) NULL else sort(unique(v))
}

#' Apply CODES=: any response not among `codes` becomes NA (not administered).
#' Existing NA stays NA. NULL/empty codes = keep the matrix unchanged.
apply_codes <- function(X, codes) {
  if (!length(codes)) return(X)
  X[!(X %in% codes)] <- NA
  X
}

#' Parse a WINSTEPS-style NEWSCORE= string (positional recode aligned to CODES=).
#' e.g. CODES="0,1,2,3" NEWSCORE="0,1,1,2" collapses categories 1 and 2.
parse_newscore <- function(txt) {
  if (is.null(txt)) return(NULL)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  parts <- strsplit(txt, "[[:space:],;]+")[[1]]
  parts <- parts[nzchar(parts)]
  v <- suppressWarnings(as.numeric(parts))
  if (any(is.na(v))) return(NULL)
  v
}

#' Apply CODES= then (optionally) NEWSCORE= to a matrix / sub-matrix.
#' Returns a list(X, note): note is a human-readable summary of what happened.
#' Order (as in WINSTEPS): invalidate out-of-CODES values, then recode.
apply_codes_recode <- function(X, codes, newscore = NULL, label = "") {
  note <- character(0)
  if (length(codes)) {
    before <- sum(!is.na(X))
    X <- apply_codes(X, codes)
    dropped <- before - sum(!is.na(X))
    note <- c(note, sprintf("%sCODES = %s (%d responses set to missing).",
                            if (nzchar(label)) paste0(label, ": ") else "",
                            paste(codes, collapse = ","), dropped))
    if (length(newscore)) {
      if (length(newscore) != length(codes)) {
        note <- c(note, sprintf("%sNEWSCORE ignored: it has %d values but CODES has %d.",
                                if (nzchar(label)) paste0(label, ": ") else "",
                                length(newscore), length(codes)))
      } else if (!identical(as.numeric(newscore), as.numeric(codes))) {
        X[] <- newscore[match(X, codes)]   # NA -> NA
        note <- c(note, sprintf("%sNEWSCORE = %s (recoded from %s).",
                                if (nzchar(label)) paste0(label, ": ") else "",
                                paste(newscore, collapse = ","),
                                paste(codes, collapse = ",")))
      }
    }
  } else if (length(newscore)) {
    note <- c(note, "NEWSCORE ignored: it requires CODES to be set.")
  }
  list(X = X, note = note)
}

# ---------------------------------------------------------------------------
# Reusable figure module (house rules 5 + honesty note about clipboard)
# ---------------------------------------------------------------------------

mod_fig_ui <- function(id, height = "620px") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("plot"), height = height),
    shiny::fluidRow(
      shiny::column(2, shiny::numericInput(ns("w"), "Export width (px)", 2000, step = 100)),
      shiny::column(2, shiny::numericInput(ns("h"), "Export height (px)", 1400, step = 100)),
      shiny::column(2, shiny::numericInput(ns("res"), "Export res (dpi)", 200, step = 25)),
      shiny::column(3, shiny::radioButtons(ns("fmt"), "Format", c("png", "pdf", "svg"), inline = TRUE)),
      shiny::column(3, shiny::br(), shiny::downloadButton(ns("dl"), "Download figure"))),
    shiny::helpText("One-click clipboard copy is not reliably possible from Shiny; use Download.")
  )
}

mod_fig_server <- function(id, draw, filename = "figure") {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot({ draw() })
    output$dl <- shiny::downloadHandler(
      filename = function() paste0(filename, "_", format(Sys.Date()), ".", input$fmt),
      content = function(file) {
        wi <- input$w / input$res; hi <- input$h / input$res
        switch(input$fmt,
               png = grDevices::png(file, width = input$w, height = input$h, res = input$res),
               pdf = grDevices::pdf(file, width = wi, height = hi),
               svg = grDevices::svg(file, width = wi, height = hi))
        on.exit(grDevices::dev.off())
        draw()
      })
  })
}

# a small card wrapper so every figure sits in a titled bslib card
fig_card <- function(title, id, height = "560px", subtitle = NULL) {
  bslib::card(
    bslib::card_header(title),
    if (!is.null(subtitle)) shiny::div(class = "text-muted small px-1", subtitle),
    mod_fig_ui(id, height),
    full_screen = TRUE
  )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ws_theme <- bslib::bs_theme(
  version = 5,
  bootswatch = "litera",
  primary = "#5B7FA6", secondary = "#B3CDE3", success = "#8FBF9F",
  base_font = bslib::font_google("Inter", local = FALSE),
  heading_font = bslib::font_google("Inter", local = FALSE)
)

## ------------------------------------------------------------------- panels
panel_data <- bslib::nav_panel(
  "Data", icon = shiny::icon("table"),
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360, title = "Load & reshape",
      shiny::strong("1. Load data (rule 1)"),
      mod_data_ui("data"),
      shiny::hr(),
      shiny::actionButton("use_demo", "Load built-in demo data",
                          class = "btn-info w-100", icon = shiny::icon("flask")),
      shiny::helpText("400 persons, 15 items, 4 categories, with 0.9-logit DIF on Item04."),
      shiny::hr(),
      shiny::strong("2. Reshape (rule 2)"),
      mod_reshape_ui("reshape"),
      shiny::helpText("Rasch analysis needs WIDE data: one row per person, one column per item.")
    ),
    bslib::layout_columns(
      col_widths = c(4, 8),
      bslib::value_box(
        title = "Parsed data",
        value = shiny::textOutput("vb_dim"),
        showcase = shiny::icon("database"), theme = "primary"),
      bslib::card(
        bslib::card_header("Data notes"),
        shiny::verbatimTextOutput("data_notes"))
    ),
    bslib::card(
      bslib::card_header("Parsed data preview"),
      DT::DTOutput("data_head"), full_screen = TRUE)
  ))

panel_estimate <- bslib::nav_panel(
  "Estimate", icon = shiny::icon("gears"),
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360, title = "Specify the model",
      shiny::strong("Items (rule 3)"),
      mod_varselect_ui("items", "Item (response) columns"),
      shiny::selectInput("person_id", "Person label column", choices = NULL),
      shiny::hr(),
      shiny::strong("Model"),
      shiny::radioButtons("model", NULL, inline = FALSE,
                          c("Rating Scale Model - all items share one structure" = "rsm",
                            "Partial Credit Model - every item its own structure" = "pcm",
                            "Grouped / mixed - assign items to scales below" = "grp")),
      shiny::conditionalPanel(
        "input.model == 'grp'",
        shiny::numericInput("n_scales", "Number of scales", 2, 1, 26),
        shiny::helpText(shiny::HTML(
          "Define each scale: pick its items and, if needed, its own ",
          "<code>CODES</code> / <code>NEWSCORE</code>. Items left unassigned go ",
          "to the first scale; an item picked in two scales takes the later one.")),
        shiny::uiOutput("scale_editor")),
      shiny::hr(),
      shiny::strong("Valid category codes (CODES=)"),
      shiny::textInput("codes", NULL, "", placeholder = "e.g. 0,1,2,3  (blank = keep all observed)"),
      shiny::helpText(shiny::HTML(
        "As in WINSTEPS <code>CODES=</code>: the response codes that are valid ",
        "categories. Any value <b>not</b> listed is treated as missing / not ",
        "administered. Applied before recoding. For the grouped model these are ",
        "the defaults; a scale's own CODES override them.")),
      shiny::strong("Recode categories (NEWSCORE=)"),
      shiny::textInput("newscore", NULL, "", placeholder = "e.g. 0,1,1,2  (align to CODES; blank = none)"),
      shiny::helpText(shiny::HTML(
        "Positional recode aligned to <code>CODES=</code>, e.g. CODES ",
        "<code>0,1,2,3</code> with NEWSCORE <code>0,1,1,2</code> collapses ",
        "categories 1 and 2.")),
      shiny::verbatimTextOutput("codes_note"),
      shiny::selectInput("recode", "Category recoding",
                         c("Shift so categories start at 0" = "shift",
                           "Collapse to consecutive observed categories" = "collapse",
                           "None (data already 0..m)" = "none")),
      shiny::hr(),
      shiny::strong("Estimation controls"),
      shiny::numericInput("maxit", "Max JMLE iterations (MJMLE=)", 400, 10, 5000),
      shiny::numericInput("lconv", "Convergence: max logit change (LCONV=)", 1e-4),
      shiny::numericInput("rconv", "Convergence: max score residual (RCONV=)", 1e-3),
      shiny::numericInput("extrsc", "Extreme score adjustment (EXTRSCORE=)", 0.3, 0.05, 0.5, 0.05),
      shiny::radioButtons("center", "Centring", c("Items (UIMEAN=0)" = "items",
                                                  "Persons (UPMEAN=0)" = "persons"), inline = TRUE),
      shiny::hr(),
      shiny::strong("User scaling (display only)"),
      shiny::numericInput("umean", "UIMEAN= (origin)", 0),
      shiny::numericInput("uscale", "USCALE= (units per logit)", 1),
      shiny::helpText("Tables are rescaled; graphs stay in logits and say so."),
      shiny::hr(),
      shiny::actionButton("run", "Estimate", class = "btn-primary btn-lg w-100",
                          icon = shiny::icon("play"))
    ),
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      shiny::uiOutput("vb_persons"),
      shiny::uiOutput("vb_items"),
      shiny::uiOutput("vb_iter"),
      shiny::uiOutput("vb_conv")),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("Estimation report"),
                  shiny::verbatimTextOutput("est_report")),
      bslib::card(bslib::card_header("Data notes"),
                  shiny::verbatimTextOutput("prep_notes"))),
    fig_card("Convergence trace", "fig_conv", "300px",
             "Max |parameter change| per JMLE iteration; dashed line = LCONV=.")
  ))

panel_summary <- bslib::nav_panel(
  "Summary", icon = shiny::icon("clipboard-list"),
  bslib::layout_columns(
    col_widths = c(3, 3, 3, 3),
    shiny::uiOutput("vb_prel"),
    shiny::uiOutput("vb_irel"),
    shiny::uiOutput("vb_psep"),
    shiny::uiOutput("vb_isep")),
  bslib::card(
    bslib::card_header("Summary of person and item measures (Table 3.1)"),
    shiny::helpText(shiny::HTML(
      "Separation = TrueSD / RMSE. Reliability = Separation&sup2;/(1+Separation&sup2;). ",
      "Strata = (4&times;Separation+1)/3. 'Real' inflates the model S.E. by ",
      "&radic;max(INFIT MNSQ, 1). Extreme scores are excluded from the ",
      "'non-extreme' rows, as in WINSTEPS.")),
    DT::DTOutput("tbl_summary"), full_screen = TRUE),
  bslib::card(
    bslib::card_header("Raw-score internal consistency (for comparison)"),
    shiny::verbatimTextOutput("alpha")))

panel_wright <- bslib::nav_panel(
  "Wright map", icon = shiny::icon("map"),
  bslib::layout_columns(
    col_widths = c(3, 3, 3, 3),
    shiny::radioButtons("wright_what", "Right-hand column",
                        c("Item measures" = "items", "Item thresholds" = "thresholds")),
    shiny::sliderInput("wright_bins", "Histogram bins", 5, 80, 30),
    shiny::numericInput("wright_maxlab", "Max item labels", 60, 5, 400),
    shiny::checkboxInput("wright_msT", "Show M / S / T markers", TRUE)),
  fig_card("Wright map (Table 1)", "fig_wright", "760px"))

panel_items <- bslib::nav_panel(
  "Items", icon = shiny::icon("list-ol"),
  bslib::card(
    bslib::card_header("Item statistics (Tables 10/14)"),
    bslib::layout_columns(
      col_widths = c(3, 3, 6),
      shiny::selectInput("item_order", "Order",
                         c("Entry number" = "entry", "Measure" = "measure",
                           "Outfit MNSQ" = "outfit", "Infit MNSQ" = "infit",
                           "Point-measure correlation" = "ptmea")),
      shiny::checkboxInput("item_discr", "Estimate discrimination (slower)", TRUE),
      shiny::helpText(shiny::HTML(
        "Productive MNSQ range on a rating scale: <b>0.6-1.4</b> ",
        "(0.7-1.3 for high-stakes dichotomies). MNSQ &gt; 2.0 degrades measurement. ",
        "ZSTD is only interpretable when MNSQ is also deviant."))),
    DT::DTOutput("tbl_items"), full_screen = TRUE),
  bslib::layout_columns(
    col_widths = c(3, 3, 6),
    shiny::selectInput("path_margin", "Pathway margin", c("items", "persons")),
    shiny::selectInput("path_stat", "Fit statistic",
                       c("Infit ZSTD" = "infit_zstd", "Outfit ZSTD" = "outfit_zstd",
                         "Infit MNSQ" = "infit_mnsq", "Outfit MNSQ" = "outfit_mnsq"))),
  fig_card("Pathway (bubble) chart", "fig_path", "560px"))

panel_persons <- bslib::nav_panel(
  "Persons", icon = shiny::icon("users"),
  bslib::card(
    bslib::card_header("Person statistics (Tables 17/18)"),
    shiny::selectInput("person_order", "Order",
                       c("Entry number" = "entry", "Measure" = "measure",
                         "Outfit MNSQ" = "outfit", "Infit MNSQ" = "infit")),
    DT::DTOutput("tbl_persons"), full_screen = TRUE),
  bslib::card(
    bslib::card_header("Include / exclude persons and re-estimate"),
    shiny::helpText(shiny::HTML(
      "Drop misfitting persons and recalibrate. Excluded persons are removed ",
      "before estimation, so item difficulties, thresholds and all downstream ",
      "tables are re-estimated without them. <b>Report this:</b> removing ",
      "misfitting persons improves fit almost by construction &mdash; say how many ",
      "were removed, on what criterion, and what it changed.")),
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      shiny::selectInput("excl_stat", "Flag persons by",
                         c("Outfit MNSQ" = "outfit_mnsq", "Infit MNSQ" = "infit_mnsq",
                           "Outfit ZSTD" = "outfit_zstd", "Infit ZSTD" = "infit_zstd")),
      shiny::numericInput("excl_cut", "Greater than", 2, 0, 100, 0.1),
      shiny::div(shiny::br(),
                 shiny::actionButton("excl_flag", "Add flagged to exclusion list",
                                     class = "btn-warning w-100"))),
    shinyWidgets::pickerInput("excl_persons", "Excluded persons",
                              choices = NULL, multiple = TRUE,
                              options = list(`live-search` = TRUE, `actions-box` = TRUE,
                                             `selected-text-format` = "count > 3")),
    shiny::verbatimTextOutput("excl_note"),
    bslib::layout_columns(
      col_widths = c(6, 6),
      shiny::actionButton("reestimate", "Re-estimate without excluded persons",
                          class = "btn-primary w-100", icon = shiny::icon("rotate")),
      shiny::actionButton("excl_clear", "Clear exclusions", class = "btn-secondary w-100"))),
  bslib::card(
    bslib::card_header("Most unexpected responses (Table 6.6)"),
    shiny::sliderInput("unexp_cut", "|standardized residual| at least", 1.5, 5, 2, 0.1),
    DT::DTOutput("tbl_unexp"), full_screen = TRUE))

panel_ratingscale <- bslib::nav_panel(
  "Rating scale", icon = shiny::icon("sliders"),
  bslib::layout_columns(
    col_widths = c(4, 4, 4),
    shiny::selectInput("cat_group", "Item group / rating scale", NULL),
    shiny::selectInput("cat_item", "Draw curves at item (optional)", NULL),
    shiny::sliderInput("cat_range", "Measure range (logits)", -10, 10, c(-6, 6))),
  bslib::layout_columns(
    col_widths = c(6, 6),
    bslib::card(bslib::card_header("Category structure (Table 3.2)"),
                DT::DTOutput("tbl_cat"), full_screen = TRUE),
    bslib::card(bslib::card_header("Category-quality guidelines (Linacre 1999, 2002)"),
                DT::DTOutput("tbl_catdiag"), full_screen = TRUE)),
  fig_card("Category probability curves", "fig_catprob", "520px"),
  fig_card("Cumulative probabilities and Rasch-Thurstone 50% thresholds", "fig_cumul", "520px"),
  fig_card("Observed vs expected average measures; threshold ordering", "fig_catdiag", "460px"))

panel_score <- bslib::nav_panel(
  "Score table", icon = shiny::icon("arrow-right-arrow-left"),
  bslib::card(
    bslib::card_header("Raw score to measure (Table 20)"),
    shiny::helpText("Conversion for a person answering every non-extreme item."),
    DT::DTOutput("tbl_score"), full_screen = TRUE))

panel_graphs <- bslib::nav_panel(
  "Graphs", icon = shiny::icon("chart-line"),
  bslib::layout_columns(
    col_widths = c(4, 4, 4),
    shiny::selectInput("graph_item", "Item", NULL),
    mod_varselect_ui("info_items", "Items for information plot"),
    shiny::sliderInput("graph_range", "Measure range", -12, 12, c(-6, 6))),
  bslib::card(bslib::card_header("Expected score curve with empirical ICC"),
              shiny::sliderInput("icc_bins", "Empirical bins", 3, 30, 10),
              mod_fig_ui("fig_icc", "520px"), full_screen = TRUE),
  fig_card("Item information functions", "fig_iteminfo", "480px"),
  fig_card("Test characteristic curve and test information", "fig_test", "480px"))

panel_dim <- bslib::nav_panel(
  "Dimensionality", icon = shiny::icon("layer-group"),
  bslib::layout_columns(
    col_widths = c(3, 3, 3, 3),
    shiny::selectInput("pca_resid", "Residual type (PRCOMP=)",
                       c("Standardized (S)" = "standardized", "Raw (R)" = "raw")),
    shiny::numericInput("pca_ncon", "Contrasts to report", 5, 1, 10),
    shiny::selectInput("pca_show", "Loading plot contrast", 1:5),
    shiny::div(shiny::actionButton("run_pca", "Run PCA of residuals", class = "btn-primary w-100"),
               shiny::br(), shiny::br(),
               shiny::actionButton("run_parallel", "Parallel analysis (slow)", class = "btn-warning w-100"))),
  bslib::card(
    bslib::card_header("Variance decomposition, in eigenvalue (item) units (Table 23)"),
    shiny::helpText(shiny::HTML(
      "WINSTEPS Table 23.0 recipe: raw residuals for the explained/unexplained split, ",
      "standardized residuals for the PCA, rescaled so the unexplained total equals the ",
      "number of non-extreme items. If <b>Observed</b> and <b>Expected</b> differ noticeably ",
      "that indicates an estimation problem, not multidimensionality.")),
    DT::DTOutput("tbl_pca"),
    shiny::verbatimTextOutput("pca_note"), full_screen = TRUE),
  fig_card("Scree plot of contrasts", "fig_scree", "460px"),
  fig_card("Contrast loading plot (Table 23.2)", "fig_contrast", "560px"),
  bslib::layout_columns(
    col_widths = c(6, 6),
    bslib::card(bslib::card_header("Contrast loadings and item clusters"),
                DT::DTOutput("tbl_load"), full_screen = TRUE),
    bslib::card(bslib::card_header("Person measures correlated across item clusters (Table 23.6)"),
                DT::DTOutput("tbl_clust"), full_screen = TRUE)),
  bslib::card(bslib::card_header("Parallel analysis of contrast sizes"),
              DT::DTOutput("tbl_parallel"), full_screen = TRUE))

panel_dif <- bslib::nav_panel(
  "DIF", icon = shiny::icon("code-branch"),
  bslib::layout_columns(
    col_widths = c(3, 2, 2, 2, 3),
    shiny::selectInput("dif_var", "Person classification (DIF=)", NULL),
    shiny::numericInput("mhslice", "MHSLICE= (logits)", 0.1, 0.01, 2, 0.05),
    shiny::numericInput("dif_minn", "Min N per cell", 5, 2, 100),
    shiny::div(shiny::br(), shiny::actionButton("run_dif", "Run DIF", class = "btn-primary w-100")),
    shiny::helpText(shiny::HTML(
      "<b>Rasch-Welch</b> anchors persons + thresholds and re-estimates a local item ",
      "difficulty per class. <b>Mantel-Haenszel</b>/<b>Mantel</b> stratify by measure. ",
      "ETS on the logit contrast: A &lt; 0.43, B 0.43-0.64, C &ge; 0.64, each p &lt; .05."))),
  bslib::card(bslib::card_header("DIF table (Table 30)"),
              DT::DTOutput("tbl_dif"), full_screen = TRUE),
  fig_card("DIF measures by class", "fig_dif", "520px"),
  fig_card("DIF contrast with ETS reference lines", "fig_difc", "520px"))

panel_keyform <- bslib::nav_panel(
  "Keyform", icon = shiny::icon("table-cells"),
  bslib::layout_columns(
    col_widths = c(5, 3, 4),
    shiny::selectInput("keyform_kind", "Keyform type",
                       c("Expected score (Table 2.2)" = "expected",
                         "Rasch-Thurstone 50% (Table 2.3)" = "thurstone",
                         "Most probable / modal (Table 2.1)" = "modal")),
    shiny::numericInput("keyform_maxitems", "Max items to show", 80, 5, 400),
    shiny::div(shiny::br(),
               shiny::checkboxInput("keyform_persons",
                                    "Show person distribution underneath", TRUE))),
  fig_card("General keyform (Table 2.2)", "fig_keyform", "760px",
           "Items are rows ordered by measure; each category number sits on the shared logit axis where it keys. The person histogram below shares that same axis."),
  bslib::card(bslib::card_header("Keyform coordinates"),
              DT::DTOutput("tbl_keyform"), full_screen = TRUE))

panel_dgf <- bslib::nav_panel(
  "DGF", icon = shiny::icon("object-group"),
  bslib::layout_columns(
    col_widths = c(4, 3, 5),
    shiny::selectInput("dgf_var", "Person classification (DIF=)", NULL),
    shiny::numericInput("dgf_minn", "Min N per cell", 5, 2, 500),
    shiny::div(shiny::br(), shiny::actionButton("run_dgf", "Run DGF", class = "btn-primary w-100"))),
  shiny::helpText(shiny::HTML(
    "Differential Group Functioning (Table 33): the interaction between the ",
    "<b>item scales</b> defined on the Estimate tab and a <b>person classification</b>. ",
    "Each item-class &times; person-class cell gets one uniform difficulty shift on top of ",
    "the baseline calibration; the DGF contrast is the difference across person classes ",
    "(ETS A/B/C as for DIF). Use a grouped model for meaningful item classes.")),
  bslib::card(bslib::card_header("DGF table (Table 33)"),
              DT::DTOutput("tbl_dgf"), full_screen = TRUE),
  fig_card("DGF: item class x person class", "fig_dgf", "520px"))

panel_style <- bslib::nav_panel(
  "Figure style", icon = shiny::icon("palette"),
  bslib::card(
    bslib::card_header("App appearance"),
    bslib::layout_columns(
      col_widths = c(6, 6),
      shiny::sliderInput("ui_font", "Interface font size (px)", 9, 20, 13, 0.5),
      shiny::sliderInput("table_font", "Table / console font size (rem)", 0.6, 1.4, 0.85, 0.05)),
    shiny::helpText(shiny::HTML(
      "Scales the whole interface (Bootstrap sizes are rem-based, so headings, ",
      "inputs and spacing scale together). The browser default is 16px. Figure ",
      "text is controlled separately by the size sliders below."))),
  bslib::card(
    bslib::card_header("Colours (rule 4) - pastel defaults (rule 9)"),
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      colourpicker::colourInput("col_person", "Persons / observed", "#B3CDE3"),
      colourpicker::colourInput("col_item", "Items / model", "#FBB4AE"),
      colourpicker::colourInput("col_border", "Borders", "#666666"),
      colourpicker::colourInput("col_line", "Reference lines", "#4D4D4D")),
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      colourpicker::colourInput("col_ref", "Grid / neutral", "#CCCCCC"),
      colourpicker::colourInput("col_pos", "Ordered / positive", "#CCEBC5"),
      colourpicker::colourInput("col_neg", "Disordered / negative", "#FDDAEC"),
      shiny::selectInput("palette", "Curve palette",
                         c("Pastel1", "Pastel2", "Set3", "Accent")))),
  bslib::card(
    bslib::card_header("Sizes (rule 5)"),
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      shiny::sliderInput("cex_label", "Label size", 0.5, 2, 0.85, 0.05),
      shiny::sliderInput("cex_point", "Point size", 0.4, 4, 1.2, 0.1),
      shiny::sliderInput("lwd", "Line width", 0.5, 6, 2, 0.25)),
    shiny::checkboxInput("show_grid", "Show grid", TRUE),
    shiny::helpText("Every figure reads these settings; export size/format are set under each figure.")))

panel_exports <- bslib::nav_panel(
  "Exports", icon = shiny::icon("file-export"),
  bslib::card(
    bslib::card_header("Exports & pipeline (rules 6, 7, 8)"),
    shiny::helpText(shiny::HTML(
      "The R script below reproduces the current analysis natively. Keep it in the ",
      "same folder as <code>rasch_engine.R</code> and <code>winsteps_plots.R</code> and ",
      "<code>source()</code> it in a clean session &mdash; you should get the same numbers.")),
    mod_exports_ui("exports")))

panel_about <- bslib::nav_panel(
  "About", icon = shiny::icon("circle-info"),
  bslib::card(shiny::htmlOutput("about")))

ui <- bslib::page_navbar(
  title = shiny::span(shiny::icon("ruler-combined"),
                      sprintf("WINSTEPPER %s", APP_VERSION)),
  id = "nav", theme = ws_theme, fillable = FALSE, collapsible = TRUE,
  window_title = "WINSTEPPER",
  header = shiny::tagList(
    shiny::tags$head(shiny::tags$style(shiny::HTML(
      ".shiny-notification{position:fixed;top:60px;right:20px;width:340px}
       .card{margin-bottom:14px}"))),
    shiny::uiOutput("font_css")),
  panel_data,
  panel_estimate,
  bslib::nav_menu(
    "Results", icon = shiny::icon("chart-column"),
    panel_summary, panel_wright, panel_keyform, panel_items, panel_persons,
    panel_ratingscale, panel_score),
  bslib::nav_menu(
    "Advanced", icon = shiny::icon("flask-vial"),
    panel_graphs, panel_dim, panel_dif, panel_dgf),
  bslib::nav_spacer(),
  bslib::nav_menu(
    "Settings", icon = shiny::icon("gear"), align = "right",
    panel_style, panel_exports, panel_about)
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # Canonical order of the exported script. Observer firing order is NOT
  # guaranteed in Shiny, so every recorded step is re-sorted into this order;
  # otherwise the exported script could plot before it estimates.
  STEP_ORDER <- c("setup", "data_load", "reshape", "prep", "estimate", "style",
                  "tables", "wright", "keyform", "pca", "parallel", "dif", "dgf")

  rec <- new_recorder()
  record_step <- function(id, description, code) {
    rec$record(id, description, code)
    cur <- shiny::isolate(rec$steps())
    known <- intersect(STEP_ORDER, names(cur))
    rec$steps(cur[c(known, setdiff(names(cur), known))])
  }
  # a no-op recorder for helper widgets that must not pollute the exported script
  null_rec <- list(record = function(...) invisible(NULL),
                   clear = function(...) invisible(NULL),
                   reset = function() invisible(NULL),
                   steps = shiny::reactiveVal(list()))
  record_step("setup", "Loaded the Rasch engine, plotting helpers and WINSTEPPER extras",
             c('source("rasch_engine.R")', 'source("winsteps_plots.R")',
               'source("winstepper_extras.R")   # keyform, DGF, corrected category diagnostics'))

  ## ---- data (rules 1 + 2) --------------------------------------------------
  demo <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$use_demo, {
    demo(demo_data())
    record_step("data_load", "Loaded the built-in demo data set",
               "data <- demo_data()   # see app.R for the generator")
    shiny::showNotification("Demo data loaded.", type = "message")
  })

  file_data <- mod_data_server("data", rec)
  raw_data <- shiny::reactive({
    d <- tryCatch(file_data(), error = function(e) NULL)
    if (!is.null(d) && nrow(as.data.frame(d)) > 0) return(as.data.frame(d))
    shiny::req(demo()); demo()
  })
  wide_data <- mod_reshape_server("reshape", raw_data, rec)

  output$data_notes <- shiny::renderText({
    d <- wide_data()
    sprintf("%d rows x %d columns\nNumeric columns: %d",
            nrow(d), ncol(d), sum(vapply(d, is.numeric, logical(1))))
  })
  output$vb_dim <- shiny::renderText({
    d <- wide_data(); sprintf("%d x %d", nrow(d), ncol(d))
  })
  output$data_head <- DT::renderDT(dt(utils::head(wide_data(), 200), 3, "40vh"))

  ## ---- variable selection (rule 3) ----------------------------------------
  item_vars <- mod_varselect_server("items", wide_data, null_rec, numeric_only = TRUE)
  shiny::observeEvent(wide_data(), {
    nm <- names(wide_data())
    shiny::updateSelectInput(session, "person_id", choices = c("(row number)", nm))
    shiny::updateSelectInput(session, "dif_var", choices = nm)
    shiny::updateSelectInput(session, "dgf_var", choices = nm)
  })

  ## ---- grouped-model scale editor (feature: assign items to scales) --------
  output$scale_editor <- shiny::renderUI({
    ns_n <- max(1, as.integer(input$n_scales %||% 1))
    iv <- item_vars()
    shiny::req(length(iv) >= 1)
    lapply(seq_len(ns_n), function(s) {
      # default split: deal items round-robin across scales for a sensible start
      def <- iv[seq_along(iv) %% ns_n == (s %% ns_n)]
      shiny::wellPanel(
        style = "padding:8px;margin-bottom:8px;",
        shiny::strong(sprintf("Scale %d", s)),
        shiny::textInput(paste0("scale_name_", s), "Label", sprintf("S%d", s)),
        shinyWidgets::pickerInput(paste0("scale_items_", s), "Items in this scale",
                                  choices = iv, selected = def, multiple = TRUE,
                                  options = list(`live-search` = TRUE, `actions-box` = TRUE,
                                                 `selected-text-format` = "count > 3")),
        shiny::textInput(paste0("scale_codes_", s), "CODES for this scale (blank = global)", ""),
        shiny::textInput(paste0("scale_newscore_", s), "NEWSCORE for this scale (blank = global)", ""))
    })
  })

  ## ---- preparation --------------------------------------------------------
  # Person labels, computed identically to prep() so the exclusion picker and
  # the estimation agree on what a person is called.
  person_labels <- shiny::reactive({
    d <- wide_data()
    rn <- if (!is.null(input$person_id) && input$person_id %in% names(d))
      as.character(d[[input$person_id]]) else sprintf("P%04d", seq_len(nrow(d)))
    make.unique(rn)
  })

  # Both "Estimate" and "Re-estimate" drive the same trigger. ignoreInit keeps
  # it from firing before the user has chosen anything.
  refit <- shiny::reactiveVal(0)
  shiny::observeEvent(input$run,        refit(refit() + 1))
  shiny::observeEvent(input$reestimate, refit(refit() + 1))

  prep <- shiny::eventReactive(refit(), ignoreInit = TRUE, {
    d <- wide_data(); iv <- item_vars(); model <- input$model
    shiny::validate(shiny::need(length(iv) >= 3, "Select at least three item columns on the Estimate tab."))
    X <- as.matrix(d[, iv, drop = FALSE]); storage.mode(X) <- "numeric"
    rownames(X) <- person_labels(); colnames(X) <- iv

    ## ---- person exclusions (re-estimation without misfitting persons) ----
    excl <- intersect(input$excl_persons %||% character(0), rownames(X))
    keep_rows <- which(!(rownames(X) %in% excl))
    shiny::validate(shiny::need(
      length(keep_rows) >= 10,
      "Fewer than 10 persons would remain after exclusions. Clear some exclusions."))
    X <- X[keep_rows, , drop = FALSE]

    global_codes <- parse_codes(input$codes)
    global_ns    <- parse_newscore(input$newscore)
    codes_notes  <- character(0)

    ## ---- build the item->group (scale) mapping --------------------------
    scale_codes <- scale_ns <- list()
    if (model == "rsm") {
      grp <- rep("R1", length(iv))
    } else if (model == "pcm") {
      grp <- iv
    } else {                                   # grouped / mixed
      ns_n <- max(1, as.integer(input$n_scales %||% 1))
      grp <- rep(NA_character_, length(iv))
      for (s in seq_len(ns_n)) {
        lab <- trimws(input[[paste0("scale_name_", s)]] %||% sprintf("S%d", s))
        if (!nzchar(lab)) lab <- sprintf("S%d", s)
        its <- input[[paste0("scale_items_", s)]]
        idx <- match(its, iv); idx <- idx[!is.na(idx)]
        if (length(idx)) grp[idx] <- lab       # a later scale wins on overlap
        sc <- parse_codes(input[[paste0("scale_codes_", s)]])
        sn <- parse_newscore(input[[paste0("scale_newscore_", s)]])
        scale_codes[[lab]] <- if (length(sc)) sc else global_codes
        scale_ns[[lab]]    <- if (length(sn)) sn else global_ns
      }
      first_lab <- trimws(input[["scale_name_1"]] %||% "S1")
      if (!nzchar(first_lab)) first_lab <- "S1"
      if (any(is.na(grp))) {
        codes_notes <- c(codes_notes, sprintf(
          "%d item(s) were unassigned and placed in scale '%s'.", sum(is.na(grp)), first_lab))
        grp[is.na(grp)] <- first_lab
      }
    }

    ## ---- apply CODES= / NEWSCORE= (per scale, or global) ----------------
    if (model == "grp") {
      for (g in unique(grp)) {
        cols <- which(grp == g)
        res <- apply_codes_recode(X[, cols, drop = FALSE], scale_codes[[g]],
                                  scale_ns[[g]], label = sprintf("Scale '%s'", g))
        X[, cols] <- res$X; codes_notes <- c(codes_notes, res$note)
      }
    } else {
      res <- apply_codes_recode(X, global_codes, global_ns)
      X <- res$X; codes_notes <- c(codes_notes, res$note)
    }

    p <- rasch_prep(X, groups = grp, recode = input$recode)
    p$keep_rows <- keep_rows          # rows of wide_data() that entered the fit
    p$excluded  <- excl
    if (length(excl))
      codes_notes <- c(sprintf("Excluded %d person(s) before estimation: %s%s",
                               length(excl), paste(utils::head(excl, 15), collapse = ", "),
                               if (length(excl) > 15) ", ..." else ""), codes_notes)
    p$notes <- c(codes_notes, p$notes)

    ## ---- reproducible code fragment -------------------------------------
    gen_cr <- function(codes, ns, colsexpr = NULL) {
      if (!length(codes)) return(character(0))
      recode <- length(ns) && length(ns) == length(codes) &&
        !identical(as.numeric(ns), as.numeric(codes))
      if (is.null(colsexpr)) {
        c(sprintf('X[!(X %%in%% c(%s))] <- NA   # CODES=', paste(codes, collapse = ", ")),
          if (recode) sprintf('X[] <- c(%s)[match(X, c(%s))]   # NEWSCORE=',
                              paste(ns, collapse = ", "), paste(codes, collapse = ", ")))
      } else {
        c(sprintf('.cols <- %s; .sub <- X[, .cols, drop = FALSE]', colsexpr),
          sprintf('.sub[!(.sub %%in%% c(%s))] <- NA   # CODES=', paste(codes, collapse = ", ")),
          if (recode) sprintf('.sub[] <- c(%s)[match(.sub, c(%s))]   # NEWSCORE=',
                              paste(ns, collapse = ", "), paste(codes, collapse = ", ")),
          'X[, .cols] <- .sub')
      }
    }
    cr_lines <- if (model == "grp")
      unlist(lapply(unique(grp), function(g)
        gen_cr(scale_codes[[g]], scale_ns[[g]], sprintf('which(groups == "%s")', g))))
    else gen_cr(global_codes, global_ns)

    record_step("prep",
               sprintf("Selected %d items; model = %s; recode = %s; scales = %s; persons excluded = %d",
                       length(iv), toupper(model), input$recode,
                       paste(sort(unique(grp)), collapse = ","), length(excl)),
               c(sprintf('items <- c(%s)', paste(sprintf('"%s"', iv), collapse = ", ")),
                 sprintf('X <- as.matrix(data[, items]); rownames(X) <- %s',
                         if (input$person_id %in% names(d))
                           sprintf('make.unique(as.character(data$`%s`))', input$person_id)
                         else 'sprintf("P%04d", seq_len(nrow(data)))'),
                 if (length(excl))
                   c(sprintf('excluded <- c(%s)   # persons dropped in the app',
                             paste(sprintf('"%s"', excl), collapse = ", ")),
                     'keep <- !(rownames(X) %in% excluded)',
                     'X <- X[keep, , drop = FALSE]')
                 else 'keep <- rep(TRUE, nrow(X))   # no persons excluded',
                 sprintf('groups <- c(%s)', paste(sprintf('"%s"', grp), collapse = ", ")),
                 cr_lines,
                 sprintf('prep <- rasch_prep(X, groups = groups, recode = "%s")', input$recode)))
    p
  })

  ## ---- estimation ---------------------------------------------------------
  fit <- shiny::eventReactive(refit(), ignoreInit = TRUE, {
    p <- prep()
    shiny::withProgress(message = "Running JMLE...", value = 0.3, {
      f <- rasch_jmle(p, maxit = input$maxit, conv = input$lconv, rconv = input$rconv,
                      extreme_adj = input$extrsc, center = input$center)
      shiny::setProgress(1)
      f
    })
  })

  shiny::observeEvent(fit(), {
    f <- fit()
    record_step("estimate",
               sprintf("Estimated the Rasch model by JMLE (%d iterations, converged = %s)",
                       f$iterations, f$converged),
               sprintf('fit <- rasch_jmle(prep, maxit = %d, conv = %g, rconv = %g,\n  extreme_adj = %g, center = "%s")',
                       input$maxit, input$lconv, input$rconv, input$extrsc, input$center))
    shiny::updateSelectInput(session, "cat_group", choices = names(f$tau))
    shiny::updateSelectInput(session, "cat_item", choices = c("(group average)", f$item_id))
    shiny::updateSelectInput(session, "graph_item", choices = f$item_id)
    record_step("tables", "Produced the WINSTEPS-style measure, fit and structure tables",
               c("summary_tab  <- summary_table(fit)      # Table 3.1",
                 "cats         <- category_table(fit)     # Table 3.2",
                 "cat_quality  <- category_diagnostics(cats)",
                 sprintf("items_tab    <- item_table(fit, discrimination = %s)  # Tables 10/14",
                         isTRUE(input$item_discr)),
                 "persons_tab  <- person_table(fit)       # Tables 17/18",
                 "scores       <- score_table(fit)        # Table 20",
                 sprintf("unexpected   <- unexpected_responses(fit, cut = %s)   # Table 6.6",
                         input$unexp_cut %||% 2)))
    shiny::showNotification(sprintf("Estimation finished in %d iterations.", f$iterations),
                            type = if (f$converged) "message" else "warning")
  })

  ## ---- value-box dashboard -------------------------------------------------
  vb <- function(title, value, icon, theme) bslib::value_box(
    title = title, value = value, showcase = shiny::icon(icon), theme = theme)

  output$vb_persons <- shiny::renderUI({
    shiny::req(fit()); f <- fit()
    vb("Persons", sprintf("%d", nrow(f$X)), "users", "secondary")
  })
  output$vb_items <- shiny::renderUI({
    shiny::req(fit()); f <- fit()
    vb("Items", sprintf("%d", ncol(f$X)), "list-ol", "secondary")
  })
  output$vb_iter <- shiny::renderUI({
    shiny::req(fit()); f <- fit()
    vb("JMLE iterations", sprintf("%d", f$iterations), "rotate", "primary")
  })
  output$vb_conv <- shiny::renderUI({
    shiny::req(fit()); f <- fit()
    vb("Converged", if (f$converged) "Yes" else "No", "circle-check",
       if (f$converged) "success" else "warning")
  })

  summ <- shiny::reactive({ shiny::req(fit()); summary_table(fit()) })
  pick <- function(tab, stat, label) {
    v <- tab[tab$Statistic == label, stat]
    if (!length(v) || is.na(v)) NA_real_ else v
  }
  output$vb_prel <- shiny::renderUI({
    s <- summ()
    vb("Person reliability (Real)", sprintf("%.2f", pick(s, "Reliability_Real", "PERSON (non-extreme)")),
       "user-check", "primary")
  })
  output$vb_irel <- shiny::renderUI({
    s <- summ()
    vb("Item reliability (Real)", sprintf("%.2f", pick(s, "Reliability_Real", "ITEM (non-extreme)")),
       "square-check", "primary")
  })
  output$vb_psep <- shiny::renderUI({
    s <- summ()
    vb("Person separation", sprintf("%.2f", pick(s, "Separation_Real", "PERSON (non-extreme)")),
       "arrows-left-right-to-line", "secondary")
  })
  output$vb_isep <- shiny::renderUI({
    s <- summ()
    vb("Item separation", sprintf("%.2f", pick(s, "Separation_Real", "ITEM (non-extreme)")),
       "arrows-left-right-to-line", "secondary")
  })

  output$est_report <- shiny::renderText({
    f <- fit()
    paste(
      sprintf("Persons: %d (%d extreme)   Items: %d (%d extreme)",
              nrow(f$X), length(f$extreme_persons), ncol(f$X), length(f$extreme_items)),
      sprintf("Persons excluded before estimation: %d",
              length(tryCatch(prep()$excluded, error = function(e) character(0)))),
      sprintf("Observations: %d   Missing: %.1f%%", sum(f$mask), 100 * mean(!f$mask)),
      sprintf("Rating-scale groups: %d   Max categories: %s",
              length(f$tau), paste(f$max_cat, collapse = ", ")),
      sprintf("JMLE iterations: %d   Converged: %s   Final max|change|: %.2e",
              f$iterations, f$converged, utils::tail(f$change_history, 1)),
      sprintf("Item measure mean/SD: %.3f / %.3f", mean(f$delta, na.rm = TRUE),
              stats::sd(f$delta, na.rm = TRUE)),
      sprintf("Person measure mean/SD: %.3f / %.3f", mean(f$theta, na.rm = TRUE),
              stats::sd(f$theta, na.rm = TRUE)),
      if (!f$converged)
        "\nWARNING: JMLE did not reach the convergence criteria. Increase MJMLE= or loosen LCONV=/RCONV= and check for sparse or disconnected data." else "",
      sep = "\n")
  })
  output$prep_notes <- shiny::renderText({
    n <- prep()$notes
    if (!length(n)) "No recoding was necessary." else paste("-", n, collapse = "\n")
  })

  # Live preview of the CODES= filter (does not require pressing Estimate).
  output$codes_note <- shiny::renderText({
    codes <- parse_codes(input$codes)
    if (!length(codes)) return("All observed values kept as valid categories.")
    iv <- tryCatch(item_vars(), error = function(e) NULL)
    d  <- tryCatch(wide_data(), error = function(e) NULL)
    msg <- sprintf("Valid codes: %s", paste(codes, collapse = ", "))
    if (!is.null(d) && length(iv) && all(iv %in% names(d))) {
      X <- suppressWarnings(as.matrix(d[, iv, drop = FALSE]))
      storage.mode(X) <- "numeric"
      obs <- !is.na(X)
      outside <- obs & !(X %in% codes)
      msg <- paste0(msg, sprintf(
        "\n%d of %d observed responses (%.1f%%) fall outside CODES= and become missing.",
        sum(outside), sum(obs), 100 * sum(outside) / max(sum(obs), 1)))
      excl <- sort(unique(X[outside]))
      if (length(excl))
        msg <- paste0(msg, "\nExcluded values: ", paste(utils::head(excl, 20), collapse = ", "),
                      if (length(excl) > 20) " ..." else "")
    }
    msg
  })

  ## ---- interface font size -------------------------------------------------
  # Bootstrap 5 sizes are rem-based, so scaling the root font size scales the
  # whole interface proportionally (browser default = 16px).
  output$font_css <- shiny::renderUI({
    px <- input$ui_font %||% 13
    rem <- input$table_font %||% 0.85
    shiny::tags$style(shiny::HTML(sprintf(
      "html{font-size:%.1fpx !important}
       .dataTables_wrapper, table.dataTable{font-size:%.2frem}
       pre, .shiny-text-output{font-size:%.2frem}
       .form-control, .form-select, .btn, label, .help-block{font-size:%.2frem}
       .bslib-value-box .value-box-value{font-size:1.6rem}",
      px, rem, rem, max(rem, 0.75))))
  })

  ## ---- shared figure style (rules 4, 5, 9) --------------------------------
  style <- shiny::reactive(ws_style(
    col_person = input$col_person, col_item = input$col_item,
    col_border = input$col_border, col_line = input$col_line,
    col_ref = input$col_ref, col_pos = input$col_pos, col_neg = input$col_neg,
    palette = input$palette, cex_label = input$cex_label,
    cex_point = input$cex_point, lwd = input$lwd, show_grid = input$show_grid,
    bins = input$wright_bins %||% 30, show_msT = isTRUE(input$wright_msT)))

  # NB: recorder$record() reads the steps reactiveVal, so it must only ever be
  # called from an isolated context (observeEvent handler / eventReactive body).
  # A plain observe() or render*() here would depend on the value it writes and
  # loop forever.
  shiny::observeEvent(style(), {
    s <- style()
    record_step("style", "Set the figure style (colours, sizes, pastel palette)",
               sprintf('style <- ws_style(col_person = "%s", col_item = "%s", col_border = "%s",\n  col_line = "%s", col_ref = "%s", col_pos = "%s", col_neg = "%s",\n  palette = "%s", cex_label = %s, cex_point = %s, lwd = %s, show_grid = %s)',
                       s$col_person, s$col_item, s$col_border, s$col_line, s$col_ref,
                       s$col_pos, s$col_neg, s$palette, s$cex_label, s$cex_point, s$lwd,
                       s$show_grid))
  })

  ## ---- convergence trace ---------------------------------------------------
  mod_fig_server("fig_conv", function() {
    shiny::req(fit()); f <- fit(); s <- style()
    graphics::par(mar = c(4.5, 4.5, 2, 1))
    graphics::plot(seq_along(f$change_history), f$change_history, type = "b", log = "y",
                   pch = 21, bg = s$col_item, col = s$col_border, lwd = s$lwd,
                   xlab = "JMLE iteration", ylab = "max |parameter change| (logits)",
                   main = "Convergence", cex.lab = s$cex_label + .1)
    graphics::abline(h = input$lconv, col = s$col_line, lty = 2, lwd = s$lwd)
    graphics::grid(col = s$col_ref, lty = 3)
  }, "convergence")

  ## ---- Table 3.1 -----------------------------------------------------------
  output$tbl_summary <- DT::renderDT({
    shiny::req(fit())
    d <- rescale_cols(summary_table(fit()), input$umean, input$uscale,
                      meas = c("Mean_Measure"), sds = c("SD_Measure", "Mean_Model_SE",
                                                        "Mean_Real_SE", "Model_RMSE", "Real_RMSE",
                                                        "True_SD_Model", "True_SD_Real"))
    dt(d, 3, "40vh")
  })
  output$alpha <- shiny::renderText({
    a <- cronbach_alpha(fit())
    sprintf("Cronbach's alpha (KR-20) on complete cases = %s\n(Alpha is a raw-score index; the Rasch person reliability above is the measure-based analogue and is usually a little lower.)",
            if (is.na(a)) "not estimable" else sprintf("%.3f", a))
  })

  ## ---- Table 3.2 -----------------------------------------------------------
  cat_tab <- shiny::reactive({ shiny::req(fit()); category_table(fit()) })
  output$tbl_cat <- DT::renderDT(
    dt(rescale_cols(cat_tab(), input$umean, input$uscale), 2, "38vh"))
  output$tbl_catdiag <- DT::renderDT(dt(category_diagnostics(cat_tab()), 2, "30vh"))

  cat_item <- shiny::reactive(
    if (is.null(input$cat_item) || input$cat_item == "(group average)") NULL else input$cat_item)

  mod_fig_server("fig_catprob", function() {
    shiny::req(fit(), input$cat_group)
    plot_category_probs(fit(), input$cat_group, style(), item = cat_item(),
                        from = input$cat_range[1], to = input$cat_range[2])
  }, "category_probability_curves")

  mod_fig_server("fig_cumul", function() {
    shiny::req(fit(), input$cat_group)
    plot_cumulative(fit(), input$cat_group, style(), item = cat_item(),
                    from = input$cat_range[1], to = input$cat_range[2])
  }, "cumulative_probabilities")

  mod_fig_server("fig_catdiag", function() {
    shiny::req(fit(), input$cat_group)
    plot_category_diagnostics(cat_tab(), input$cat_group, style())
  }, "category_diagnostics")

  ## ---- items ---------------------------------------------------------------
  item_tab <- shiny::reactive({
    shiny::req(fit())
    d <- item_table(fit(), discrimination = isTRUE(input$item_discr))
    ord <- switch(input$item_order,
                  entry = order(d$Entry), measure = order(-d$Measure),
                  outfit = order(-d$Outfit_MNSQ), infit = order(-d$Infit_MNSQ),
                  ptmea = order(d$PtMeasure_Corr))
    d[ord, ]
  })
  output$tbl_items <- DT::renderDT(
    dt(rescale_cols(item_tab(), input$umean, input$uscale), 2, "50vh"))

  mod_fig_server("fig_path", function() {
    shiny::req(fit())
    plot_pathway(fit(), style(), margin = input$path_margin, stat = input$path_stat)
  }, "pathway_chart")

  ## ---- persons -------------------------------------------------------------
  person_tab <- shiny::reactive({
    shiny::req(fit())
    d <- person_table(fit())
    ord <- switch(input$person_order,
                  entry = order(d$Entry), measure = order(-d$Measure),
                  outfit = order(-d$Outfit_MNSQ), infit = order(-d$Infit_MNSQ))
    d[ord, ]
  })
  output$tbl_persons <- DT::renderDT(
    dt(rescale_cols(person_tab(), input$umean, input$uscale), 2, "45vh"))
  output$tbl_unexp <- DT::renderDT({
    shiny::req(fit())
    dt(utils::head(unexpected_responses(fit(), input$unexp_cut), 500), 2, "35vh")
  })

  ## ---- person include / exclude + re-estimation ---------------------------
  # Choices come from the full data, not from the current fit, so a person that
  # has been excluded can still be found in the list and put back.
  shiny::observeEvent(person_labels(), {
    shinyWidgets::updatePickerInput(session, "excl_persons",
                                    choices = person_labels(),
                                    selected = intersect(input$excl_persons %||% character(0),
                                                         person_labels()))
  })
  shiny::observeEvent(input$excl_flag, {
    shiny::req(fit())
    d <- person_table(fit())
    v <- switch(input$excl_stat,
                outfit_mnsq = d$Outfit_MNSQ, infit_mnsq = d$Infit_MNSQ,
                outfit_zstd = d$Outfit_ZSTD, infit_zstd = d$Infit_ZSTD)
    hit <- d$Person[is.finite(v) & v > input$excl_cut & fit()$keep_p]
    new <- union(input$excl_persons %||% character(0), hit)
    shinyWidgets::updatePickerInput(session, "excl_persons",
                                    choices = person_labels(), selected = new)
    shiny::showNotification(
      sprintf("%d person(s) flagged (%s > %g); %d now on the exclusion list. Press Re-estimate to apply.",
              length(hit), input$excl_stat, input$excl_cut, length(new)),
      type = "message")
  })
  shiny::observeEvent(input$excl_clear, {
    shinyWidgets::updatePickerInput(session, "excl_persons",
                                    choices = person_labels(), selected = character(0))
    shiny::showNotification("Exclusion list cleared. Press Re-estimate to apply.", type = "message")
  })
  output$excl_note <- shiny::renderText({
    sel <- input$excl_persons %||% character(0)
    n <- length(person_labels())
    applied <- tryCatch(length(prep()$excluded), error = function(e) NA_integer_)
    paste0(
      sprintf("%d of %d persons on the exclusion list (%.1f%%).", length(sel), n,
              100 * length(sel) / max(n, 1)),
      if (!is.na(applied)) sprintf("\nCurrently estimated model excludes %d.", applied) else "",
      if (!is.na(applied) && applied != length(sel))
        "\nThe list has changed since the last estimation - press 'Re-estimate'." else "")
  })

  ## ---- Wright map ----------------------------------------------------------
  mod_fig_server("fig_wright", function() {
    shiny::req(fit())
    plot_wright(fit(), style(), what = input$wright_what,
                max_labels = input$wright_maxlab)
  }, "wright_map")
  shiny::observeEvent(list(fit(), input$wright_what, input$wright_maxlab), {
    record_step("wright", "Drew the Wright (person-item) map",
               sprintf('plot_wright(fit, style, what = "%s", max_labels = %d)',
                       input$wright_what, input$wright_maxlab))
  })

  ## ---- Keyform (Table 2.2) -------------------------------------------------
  mod_fig_server("fig_keyform", function() {
    shiny::req(fit())
    plot_keyform(fit(), style(), kind = input$keyform_kind,
                 max_items = input$keyform_maxitems,
                 show_persons = isTRUE(input$keyform_persons))
  }, "keyform")
  output$tbl_keyform <- DT::renderDT({
    shiny::req(fit())
    dt(rescale_cols(keyform_data(fit(), input$keyform_kind), input$umean, input$uscale,
                    meas = c("Measure", "Item_Measure")), 3, "35vh")
  })
  shiny::observeEvent(list(fit(), input$keyform_kind, input$keyform_maxitems,
                           input$keyform_persons), {
    record_step("keyform", sprintf("Drew the general keyform (Table 2.2, %s)", input$keyform_kind),
               sprintf('plot_keyform(fit, style, kind = "%s", max_items = %d, show_persons = %s)',
                       input$keyform_kind, input$keyform_maxitems,
                       isTRUE(input$keyform_persons)))
  })

  ## ---- graphs --------------------------------------------------------------
  info_items <- mod_varselect_server("info_items", shiny::reactive({
    f <- fit(); as.data.frame(f$X)
  }), null_rec, numeric_only = FALSE)

  mod_fig_server("fig_icc", function() {
    shiny::req(fit(), input$graph_item)
    plot_expected_score(fit(), input$graph_item, style(), nbins = input$icc_bins,
                        from = input$graph_range[1], to = input$graph_range[2])
  }, "expected_score_curve")

  mod_fig_server("fig_iteminfo", function() {
    shiny::req(fit())
    sel <- info_items(); if (!length(sel)) sel <- utils::head(fit()$item_id, 8)
    plot_item_info(fit(), sel, style(), from = input$graph_range[1], to = input$graph_range[2])
  }, "item_information")

  mod_fig_server("fig_test", function() {
    shiny::req(fit())
    plot_test_curves(fit(), style(), "both", input$graph_range[1], input$graph_range[2])
  }, "test_curves")

  ## ---- Table 20 ------------------------------------------------------------
  output$tbl_score <- DT::renderDT({
    shiny::req(fit())
    dt(rescale_cols(score_table(fit()), input$umean, input$uscale), 3, "60vh")
  })

  ## ---- Table 23 ------------------------------------------------------------
  pca <- shiny::eventReactive(input$run_pca, {
    shiny::req(fit())
    shiny::withProgress(message = "PCA of residuals...", value = .4,
      pca_residuals(fit(), ncontrast = input$pca_ncon, residual_type = input$pca_resid))
  })
  shiny::observeEvent(pca(), {
    record_step("pca", sprintf("Ran the PCA of %s residuals with %d contrasts (Table 23)",
                              input$pca_resid, input$pca_ncon),
               sprintf('pca <- pca_residuals(fit, ncontrast = %d, residual_type = "%s")',
                       input$pca_ncon, input$pca_resid))
  })
  output$tbl_pca <- DT::renderDT({
    shiny::validate(shiny::need(!is.null(pca()), "Press 'Run PCA of residuals'."))
    dt(pca()$variance, 2, "35vh")
  })
  output$pca_note <- shiny::renderText({
    p <- pca(); shiny::req(p)
    sprintf(paste0("Non-extreme items: %d   Non-extreme persons: %d\n",
                   "Essential unidimensionality (Rasch / common variance) = %.1f%%\n",
                   "First contrast eigenvalue = %.2f item units.\n",
                   "Guidance: eigenvalues up to ~1.4 are at the random level (Smith & Miao 1994); ",
                   "values as high as 2.0 can occur by chance (Raiche 2005). Compare against the ",
                   "parallel analysis below rather than a fixed cut-off, and always inspect the ",
                   "content of the items at the two poles of the contrast before concluding ",
                   "multidimensionality."),
            p$n_items, p$n_persons, p$essential_unidimensionality, p$eigenvalues[1])
  })
  output$tbl_load <- DT::renderDT(
    dt(rescale_cols(pca()$loadings, input$umean, input$uscale), 3, "40vh"))
  output$tbl_clust <- DT::renderDT(dt(pca()$cluster_correlations, 3, "25vh"))

  mod_fig_server("fig_scree", function() {
    shiny::req(pca())
    plot_scree(pca(), style(), simulated = par_res())
  }, "scree_contrasts")
  mod_fig_server("fig_contrast", function() {
    shiny::req(pca())
    k <- min(as.integer(input$pca_show), input$pca_ncon)
    plot_pca_contrast(pca(), k, style())
  }, "contrast_loadings")

  par_res <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$run_parallel, {
    shiny::withProgress(message = "Simulating Rasch-conforming data...", value = .1, {
      r <- parallel_contrasts(prep(), fit(), nsim = 10, ncontrast = input$pca_ncon)
      par_res(r)
      record_step("parallel", "Ran a parallel analysis: 10 Rasch-conforming simulations of the contrast eigenvalues",
                 sprintf("set.seed(1)\nparallel <- parallel_contrasts(prep, fit, nsim = 10, ncontrast = %d)",
                         input$pca_ncon))
    })
  })
  output$tbl_parallel <- DT::renderDT({
    r <- par_res()
    shiny::validate(shiny::need(!is.null(r), "Press 'Run parallel analysis' to simulate the expected contrast sizes under the Rasch model."))
    dt(r, 3, "25vh")
  })

  ## ---- Table 30 DIF --------------------------------------------------------
  dif_res <- shiny::eventReactive(input$run_dif, {
    d <- wide_data(); shiny::req(input$dif_var)
    # Subset to the rows that actually entered the fit (person exclusions).
    cls <- d[[input$dif_var]][prep()$keep_rows]
    shiny::validate(shiny::need(length(unique(stats::na.omit(cls))) >= 2,
                                "The DIF classification needs at least two groups."))
    shiny::withProgress(message = "DIF analysis...", value = .4,
      dif_analysis(fit(), cls, mhslice = input$mhslice, min_n = input$dif_minn))
  })
  shiny::observeEvent(dif_res(), {
    record_step("dif", sprintf("Ran DIF for '%s' (Rasch-Welch + Mantel-Haenszel, MHSLICE = %g)",
                              input$dif_var, input$mhslice),
               sprintf('dif <- dif_analysis(fit, data$`%s`[keep], mhslice = %g, min_n = %d)',
                       input$dif_var, input$mhslice, input$dif_minn))
  })
  output$tbl_dif <- DT::renderDT({
    d <- dif_res()
    dt(rescale_cols(d, input$umean, input$uscale), 3, "45vh")
  })
  mod_fig_server("fig_dif", function() { shiny::req(dif_res()); plot_dif(dif_res(), style()) }, "dif_measures")
  mod_fig_server("fig_difc", function() { shiny::req(dif_res()); plot_dif_contrast(dif_res(), style()) }, "dif_contrast")

  ## ---- Table 33 DGF --------------------------------------------------------
  dgf_res <- shiny::eventReactive(input$run_dgf, {
    d <- wide_data(); shiny::req(input$dgf_var)
    cls <- d[[input$dgf_var]][prep()$keep_rows]
    shiny::validate(shiny::need(length(unique(stats::na.omit(cls))) >= 2,
                                "The DGF person classification needs at least two groups."))
    f <- fit()
    shiny::validate(shiny::need(length(unique(f$groups[f$keep_i])) >= 2,
                                "DGF needs at least two item classes; estimate a grouped model on the Estimate tab."))
    shiny::withProgress(message = "DGF analysis...", value = .4,
      dgf_analysis(f, cls, item_group = f$groups, min_n = input$dgf_minn))
  })
  shiny::observeEvent(dgf_res(), {
    record_step("dgf", sprintf("Ran DGF for '%s' by item scale (Table 33)", input$dgf_var),
               sprintf('dgf <- dgf_analysis(fit, data$`%s`[keep], item_group = fit$groups, min_n = %d)',
                       input$dgf_var, input$dgf_minn))
  })
  output$tbl_dgf <- DT::renderDT({
    dt(rescale_cols(dgf_res(), input$umean, input$uscale,
                    meas = character(0),
                    sds = c("Shift_A", "Shift_B", "SE_A", "SE_B", "Joint_SE", "DGF_Contrast")), 3, "45vh")
  })
  mod_fig_server("fig_dgf", function() { shiny::req(dgf_res()); plot_dgf(dgf_res(), style()) }, "dgf")

  ## ---- exports (rules 6, 7, 8) --------------------------------------------
  results <- shiny::reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    list(fit = f,
         items = tryCatch(item_table(f), error = function(e) NULL),
         persons = tryCatch(person_table(f), error = function(e) NULL),
         summary = tryCatch(summary_table(f), error = function(e) NULL),
         categories = tryCatch(category_table(f), error = function(e) NULL),
         pca = tryCatch(pca(), error = function(e) NULL),
         dif = tryCatch(dif_res(), error = function(e) NULL),
         dgf = tryCatch(dgf_res(), error = function(e) NULL),
         keyform = tryCatch(keyform_data(f, "expected"), error = function(e) NULL),
         session = utils::sessionInfo())
  })
  mod_exports_server("exports", rec, PKGS, wide_data, results)

  ## ---- about ---------------------------------------------------------------
  output$about <- shiny::renderUI(shiny::HTML(sprintf('
<h3>WINSTEPPER %s</h3>
<p>A Shiny front end for WINSTEPS-style Rasch measurement, implemented in plain R.
This is a redesigned interface built on the same estimation engine
(<code>rasch_engine.R</code>) as R-Winsteps: JMLE (UCON) for the dichotomous Rasch
model, the Andrich Rating Scale Model and the Masters Partial Credit Model, with
item grouping, extreme-score handling, INFIT/OUTFIT statistics, separation and
reliability, category structure statistics, DIF, and the WINSTEPS Table 23 variance
decomposition. The engine is unchanged, so the numbers are the audited ones.</p>

<h4>What has been verified (engine)</h4>
<ul>
<li>JMLE reaches the same solution as a general-purpose maximiser of the joint
likelihood (identical deviance to six decimals on test data).</li>
<li>At convergence, observed and model-expected raw score marginals agree to
&lt; 1e-6 for persons and items, i.e. the sufficient statistics are recovered.</li>
<li>Item and threshold parameters are recovered from simulated data; the mild
spread inflation you will see is the well-known JMLE estimation bias, which
WINSTEPS shares unless STBIAS= correction is enabled.</li>
<li>Fit mean-squares average 1.00 with SD(ZSTD) near 1 on model-conforming data;
the observed and model-expected variance-explained percentages in Table 23 agree
to within a few hundredths of a percent.</li>
</ul>

<h4>Where this differs from WINSTEPS 5.11 - please read</h4>
<ul>
<li><b>Not byte-identical.</b> WINSTEPS uses proprietary convergence acceleration
and its own missing-data and extreme-score conventions. Expect agreement to about
two decimals on well-behaved data, and larger differences on sparse, disconnected
or barely-converged data sets.</li>
<li><b>No JMLE bias correction.</b> There is no STBIAS= equivalent here.</li>
<li><b>Category INFIT MNSQ.</b> Reported as the information-weighted mean squared
residual within the category. Unlike the category OUTFIT it is not centred on 1.0
in the extreme categories; Linacre&#39;s &lt; 2.0 guideline applies to the OUTFIT column.</li>
<li><b>Threshold advance criterion.</b> The category-quality check uses the
category-count-dependent minimum advance
ln((k+1)(m+1&minus;k)/(k(m&minus;k))), not the flat &ldquo;1.4 logits&rdquo; often
quoted &mdash; 1.4 is only the three-category case (2&nbsp;ln&nbsp;2 = 1.386).
For 4 categories the requirement is 1.10, for 5 it is 0.98/0.81/0.98
(Linacre, <i>RMT</i> 2006, 20:1, p. 1052).</li>
<li><b>PCA clusters</b> are tertiles of the first-contrast loadings, an
approximation to the WINSTEPS three-cluster split.</li>
<li><b>Keyform (Table 2.2)</b> is the general expected-score keyform, with
Rasch-Thurstone (2.3) and modal (2.1) variants; the scalogram/keyform data-entry
forms are not reproduced.</li>
<li><b>DGF (Table 33)</b> estimates one uniform difficulty shift per
item-class &times; person-class cell (item classes come from the model scales) and
contrasts them across person classes with the ETS A/B/C rule &mdash; a Rasch-Welch
style summary, not WINSTEPS&#39; exact log-linear DGF.</li>
<li><b>Anchoring</b> (IAFILE=, PAFILE=, SAFILE=), CUTLO/CUTHI, subset detection
and the scalogram table are not implemented.</li>
<li>This is independent software. It is not affiliated with, endorsed by, or
derived from WINSTEPS&reg;, which is John M. Linacre&#39;s software. Cite the engine
you actually ran.</li>
</ul>

<h4>Suggested reporting</h4>
<p>Report the model (RSM/PCM), the estimation method (JMLE), the convergence
criteria, how extreme scores were handled, person and item separation and
reliability, the fit criteria you applied and the number of misfitting items,
the category-quality evidence, the first-contrast eigenvalue <em>with</em> a
parallel analysis, and the DIF method with its effect-size classification.
The Exports tab gives you an R script and a plain-language pipeline log that
together make all of that auditable.</p>

<h4>Key references</h4>
<p style="font-size:12px">
Andrich, D. (1978). A rating formulation for ordered response categories.
<i>Psychometrika, 43</i>, 561-573.<br>
Linacre, J. M. (1999). Investigating rating scale category utility.
<i>Journal of Outcome Measurement, 3</i>, 103-122.<br>
Linacre, J. M. (2002). Optimizing rating scale category effectiveness.
<i>Journal of Applied Measurement, 3</i>, 85-106.<br>
Mantel, N. (1963). Chi-square tests with one degree of freedom.
<i>JASA, 58</i>, 690-700.<br>
Mantel, N., &amp; Haenszel, W. (1959). Statistical aspects of the analysis of data
from retrospective studies. <i>JNCI, 22</i>, 719-748.<br>
Masters, G. N. (1982). A Rasch model for partial credit scoring.
<i>Psychometrika, 47</i>, 149-174.<br>
Raiche, G. (2005). Critical eigenvalue sizes in standardized residual principal
components analysis. <i>Rasch Measurement Transactions, 19</i>, 1012.<br>
Smith, R. M., &amp; Miao, C. Y. (1994). Assessing unidimensionality for Rasch
measurement. In M. Wilson (Ed.), <i>Objective Measurement, Vol. 2</i>.<br>
Wright, B. D., &amp; Masters, G. N. (1982). <i>Rating Scale Analysis</i>. MESA Press.<br>
Wright, B. D., &amp; Linacre, J. M. (1994). Reasonable mean-square fit values.
<i>Rasch Measurement Transactions, 8</i>, 370.
</p>', APP_VERSION)))
}

shiny::shinyApp(ui, server)
