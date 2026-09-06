# ==============================================================================
# ui/10_import.R — SHARED data import + variable selection (tab "import")
#
# Merged by hand from the two source apps' import tabs, which were separate
# implementations of the same step with the same input ids:
#   WaPaa1_3.R lines  84-142   (group-variable picker, sample data with groups,
#                               row-count control on the preview)
#   CIRCAREG.R lines  66-108   (explicit separator choice, Excel sheet picker,
#                               scalar-variable picker)
#
# The union keeps every control from both.  ONE selection step now defines:
#   * the functional/time columns  -> values$data      (all analyses)
#   * the scalar variables         -> values$covariates (FoSR / cosinor)
#                                  -> values$group_variables, values$group_labels
#                                     (fANOVA / clustering / group comparisons)
# so a variable chosen once is available as a predictor AND as a grouping
# factor, instead of being picked twice in two apps that could disagree.
# ==============================================================================

ui_tab_import <- tabItem(
  tabName = "import",

  fluidRow(
    box(
      title = "F*CK — Functional data analysis, Circadian regression, K-means clustering",
      status = "success", solidHeader = TRUE, width = 12, collapsible = TRUE,
      helpText(HTML(
        "One import, one smoothing step, every analysis.<br>",
        "<b>Import &rarr; Preprocessing/Smoothing</b> are shared: do them once and every ",
        "tab below reads the same curves.<br>",
        "<b>F</b> &mdash; functional PCA / time-warped PCA, functional ANOVA and post-hoc tests.<br>",
        "<b>C</b> &mdash; function-on-scalar regression, harmonic (cosinor) ",
        "regression and pairwise tests on circadian parameters.<br>",
        "<b>K</b> &mdash; functional clustering (k-means, hierarchical, DCF) with elbow and ",
        "silhouette diagnostics."))
    )
  ),

  fluidRow(
    box(
      title = "1. Upload Data",
      collapsible = TRUE, collapsed = FALSE,
      status = "primary",
      solidHeader = TRUE,
      width = 4,
      fileInput("datafile", "Choose Data File",
                accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx")),
      helpText("Supported formats: CSV, TXT, TSV, Excel (.xls, .xlsx)"),
      hr(),
      checkboxInput("header", "File has header row", TRUE),

      # CIRCAREG offered an explicit separator, WaPaa sniffed it from the first
      # line. Both are kept: "Auto-detect" is WaPaa's behaviour and the default.
      conditionalPanel(
        condition = "input.datafile && !input.datafile.name.match(/\\.(xls|xlsx)$/i)",
        radioButtons("sep", "Separator:",
                     choices = c("Auto-detect" = "auto", Comma = ",",
                                 Semicolon = ";", Tab = "\t"),
                     selected = "auto", inline = TRUE)
      ),
      conditionalPanel(
        condition = "input.datafile && input.datafile.name.match(/\\.(xls|xlsx)$/i)",
        uiOutput("excel_sheet_selector")
      ),
      hr(),
      radioButtons("data_format", "Data Format:",
                   choices = list("Wide (subjects in rows)" = "wide",
                                  "Long (subjects in columns)" = "long"),
                   selected = "wide"),
      helpText("Wide: Rows=Subjects, Cols=Time/Vars. Long: Cols=Subjects (Transpose)."),
      hr(),
      actionButton("load_data", "Load Raw File", class = "btn-warning"),
      hr(),
      actionButton("generate_sample", "Generate Sample Data", class = "btn-info"),
      checkboxInput("generate_with_groups", "Include group structure", TRUE),
      conditionalPanel(
        condition = "input.generate_with_groups == true",
        numericInput("n_groups", "Number of groups:", value = 3, min = 2, max = 10)
      ),
      helpText("The sample set is a 24-hour circadian dataset (50 subjects, hourly",
               "columns 00:00-23:00) with group structure and scalar variables, so it",
               "exercises every tab: fPCA, fANOVA, clustering, FoSR and cosinor.")
    ),

    box(
      title = "2. Variable Selection",
      collapsible = TRUE, collapsed = FALSE,
      status = "warning",
      solidHeader = TRUE,
      width = 8,
      uiOutput("var_select_container"),
      hr(),
      h4("Data Status:"),
      verbatimTextOutput("data_status")
    )
  ),

  # WaPaa defines group_summary and group_preview_plot in its server but no tab
  # ever placed them, so they were dead code in both the source app and the
  # first cut of this merge. They belong next to the variable selection that
  # produces them.
  fluidRow(
    box(
      title = "Group Structure",
      status = "warning",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      helpText("The primary grouping variable — the first scalar variable you selected.",
               "Each analysis tab can pick a different one."),
      fluidRow(
        column(4, DTOutput("group_summary")),
        column(8, plotOutput("group_preview_plot", height = "280px"))
      )
    )
  ),

  fluidRow(
    box(
      title = "Data Preview",
      collapsible = TRUE, collapsed = FALSE,
      status = "primary",
      width = 12,
      fluidRow(
        column(4,
          selectInput("data_preview_rows", "Rows to display:",
                      choices = c("5" = 5, "10" = 10, "20" = 20, "50" = 50,
                                  "100" = 100, "All" = -1),
                      selected = 10, width = "150px")
        )
      ),
      DTOutput("data_preview"),
      hr(),
      plotOutput("raw_data_plot", height = "300px")
    )
  )
)
