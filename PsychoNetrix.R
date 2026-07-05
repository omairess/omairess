###############################################################################
# Psychonetrics Shiny App — Latent Network Models (LNM / RNM / LRNM / GGM)
# Author: generated with Claude
# Requires: psychonetrics ≥ 0.15, shiny, shinydashboard, shinyWidgets,
#           DT, qgraph, dplyr, readxl, haven, ggplot2, viridis, Matrix
###############################################################################

# ── 0. Packages ──────────────────────────────────────────────────────────────
required_pkgs <- c(

  "shiny", "shinydashboard", "shinyWidgets", "colourpicker", "DT",
  "psychonetrics", "qgraph", "dplyr", "readxl", "haven",
  "ggplot2", "viridis", "Matrix", "corpcor"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── 1. UI ────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Psychonetrics – Latent Network Models"),

  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("1. Data Import",     tabName = "tab_data",    icon = icon("file-import")),
      menuItem("2. Variable Setup",  tabName = "tab_vars",    icon = icon("sliders")),
      menuItem("3. Model",           tabName = "tab_model",   icon = icon("project-diagram")),
      menuItem("4. Results",         tabName = "tab_results", icon = icon("table")),
      menuItem("5. Plots",           tabName = "tab_plots",   icon = icon("network-wired")),
      menuItem("6. Advanced",        tabName = "tab_adv",     icon = icon("cogs")),
      menuItem("About",              tabName = "tab_about",   icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f7f9fc; }
        .box { border-top: 3px solid #3c8dbc; }
        .shiny-output-error { color: #c0392b; }
        pre.r-output { max-height: 500px; overflow-y: auto; background: #fdfdfd;
                        padding: 10px; border: 1px solid #ddd; border-radius: 4px; }
        .nav-tabs-custom > .tab-content { padding: 15px; }
        .lambda-cell { width: 60px; text-align: center; }
      "))
    ),

    tabItems(

      # ── Tab 1: Data Import ────────────────────────────────────────────────
      tabItem(
        tabName = "tab_data",
        fluidRow(
          box(
            title = "Upload Data File", status = "primary", solidHeader = TRUE,
            width = 6,
            fileInput("file_upload", "Choose file",
                      accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls", ".sav", ".por")),
            helpText("Supported: .csv, .txt/.tsv, .xlsx/.xls (Excel), .sav/.por (SPSS)"),

            hr(),
            h4("Import Options"),

            radioButtons("data_format", "Data layout:",
                         choices = c("Wide (subjects × variables)" = "wide",
                                     "Long (stacked / repeated measures)" = "long"),
                         selected = "wide"),

            conditionalPanel(
              condition = "input.data_format == 'long'",
              textInput("id_var",   "Subject / ID variable:", placeholder = "e.g. id"),
              textInput("time_var", "Time / wave variable:",  placeholder = "e.g. wave")
            ),

            hr(),
            h4("CSV / TXT options"),
            fluidRow(
              column(4, selectInput("sep", "Separator",
                                    choices = c("Comma" = ",", "Semicolon" = ";",
                                                "Tab" = "\t", "Space" = " "))),
              column(4, selectInput("dec", "Decimal",
                                    choices = c("Period" = ".", "Comma" = ","))),
              column(4, checkboxInput("header", "First row is header", value = TRUE))
            ),

            conditionalPanel(
              condition = "input.file_upload != null &&
                           (input.file_upload.name.endsWith('.xlsx') ||
                            input.file_upload.name.endsWith('.xls'))",
              numericInput("excel_sheet", "Sheet number:", value = 1, min = 1)
            ),

            hr(),
            h4(icon("star"), " Demo Dataset"),
            actionButton("load_demo", "Load Star Wars (psychonetrics)",
                         icon = icon("jedi"), class = "btn-success btn-sm"),
            helpText("271 respondents · 10 items · 3 factors: Prequel (Q2–Q4),",
                     "Originals (Q5–Q7), Sequel (Q8–Q10). Q1 cross-loads on all three.",
                     "Lambda matrix is pre-configured automatically.")
          ),

          box(
            title = "Data Preview", status = "info", solidHeader = TRUE,
            width = 6,
            verbatimTextOutput("data_summary"),
            hr(),
            DT::dataTableOutput("data_table")
          )
        )
      ),

      # ── Tab 2: Variable Setup ─────────────────────────────────────────────
      tabItem(
        tabName = "tab_vars",
        fluidRow(
          box(
            title = "Select Variables & Specify Types", status = "primary",
            solidHeader = TRUE, width = 12,

            fluidRow(
              column(
                6,
                pickerInput("selected_vars", "Variables to include in the analysis:",
                            choices = NULL, multiple = TRUE,
                            options = pickerOptions(
                              actionsBox = TRUE,
                              liveSearch = TRUE,
                              title = "Select variables…"
                            ))
              ),
              column(
                3,
                selectInput("data_type", "Variable type:",
                            choices = c("Continuous"           = "continuous",
                                        "Ordered categorical"  = "ordered",
                                        "Dichotomous (0/1)"    = "dichotomous"))
              ),
              column(
                3,
                selectInput("estimator_type", "Estimator:",
                            choices = c("ML (auto → FIML if NA)" = "default",
                                        "ML"    = "ML",
                                        "FIML"  = "FIML",
                                        "DWLS (ordinal)"        = "DWLS",
                                        "ULS"   = "ULS",
                                        "WLS"   = "WLS",
                                        "PML (penalised)"       = "PML")),
                helpText(icon("info-circle"),
                         "PML (Penalised ML) applies a penalty to encourage sparse solutions."),
                uiOutput("dwls_ordered_warn")
              ),
              column(
                3,
                selectInput("optimizer", "Optimizer:",
                            choices = c("nlminb (default)" = "nlminb",
                                        "ucminf"           = "ucminf"),
                            selected = "nlminb")
              ),
              column(
                3,
                selectInput("missing_method", "Missing data:",
                            choices = c("Listwise deletion" = "listwise",
                                        "Pairwise deletion" = "pairwise"),
                            selected = "listwise"),
                helpText(icon("info-circle"),
                         "Pairwise uses all available variable pairs (GGM only);
                          listwise uses complete cases.")
              )
            ),

            fluidRow(
              column(
                5,
                selectInput(
                  "data_transform", "Data transformation (applied before fitting):",
                  choices = c(
                    "None"                                  = "none",
                    "Nonparanormal (rank \u2192 normal scores)"  = "npn",
                    "Dichotomize at mean (\u2192 0/1)"           = "dichot_mean",
                    "Dichotomize at median (\u2192 0/1)"         = "dichot_median",
                    "Log (log1p, for positive skew)"        = "log1p",
                    "Square root (for counts)"              = "sqrt",
                    "Z-score (standardize)"                 = "zscore"
                  ),
                  selected = "none"
                )
              ),
              column(
                7,
                conditionalPanel(
                  condition = "input.data_transform != 'none'",
                  br(),
                  helpText(icon("info-circle"),
                           "Transformation applied to all selected continuous variables
                            before model fitting. Original data is not modified.")
                )
              )
            ),

            conditionalPanel(
              condition = "input.data_format == 'long'",
              helpText("For long-format data the app will automatically pivot to wide
                        using the ID and time variables specified on the Data Import tab.")
            ),

            hr(),

            conditionalPanel(
              condition = "input.selected_vars != null && input.selected_vars.length > 0",
              h4("Optional: grouping variable (multi-group models)"),
              selectizeInput("group_var", "Group variable (leave blank for single group):",
                             choices = NULL,
                             options = list(placeholder = "(none)", allowEmptyOption = TRUE))
            )
          )
        ),

        fluidRow(
          box(
            title = "Descriptive Statistics for Selected Variables", status = "info",
            solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = FALSE,
            DT::dataTableOutput("desc_table"),
            hr(),
            plotOutput("cor_heatmap", height = "450px")
          )
        )
      ),

      # ── Tab 3: Model Specification ────────────────────────────────────────
      tabItem(
        tabName = "tab_model",
        fluidRow(
          box(
            title = "Model Family", status = "primary", solidHeader = TRUE, width = 4,
            radioButtons("model_family", "Choose model family:",
                         choices = c(
                           "GGM (Gaussian Graphical Model)"               = "ggm",
                           "Ising (binary / dichotomous)"                  = "ising",
                           "CFA / SEM (lvm)"                              = "lvm",
                           "Latent Network Model (LNM)"                   = "lnm",
                           "Residual Network Model (RNM)"                 = "rnm",
                           "Latent + Residual Network Model (LRNM)"       = "lrnm",
                           "── Panel / Longitudinal ──────────────"       = "sep_panel",
                           "Panel DLVM (latent cross-lagged)"             = "dlvm1",
                           "Panel GVAR (observed cross-lagged)"           = "panelgvar",
                           "RI-CLPM (random intercept)"                   = "ri_clpm"
                         ), selected = "lnm"),

            hr(),

            conditionalPanel(
              condition = "['lvm','lnm','rnm','lrnm','dlvm1','ri_clpm'].includes(input.model_family)",
              h4("Identification"),
              radioButtons("identification", "Identification method:",
                           choices = c("Loading (fix first loading)" = "loadings",
                                       "Variance (fix latent var)"   = "variance"),
                           selected = "variance")
            ),

            hr(),
            h4("Search & Pruning"),
            checkboxInput("do_prune",  "Prune non-significant edges", value = TRUE),
            conditionalPanel(
              condition = "input.do_prune",
              fluidRow(
                column(6, numericInput("prune_alpha", "Prune α:", value = 0.01,
                                       min = 0.001, max = 0.10, step = 0.005)),
                column(6, selectInput("prune_adjust", "Adjust p-values:",
                                      choices = c("None"       = "none",
                                                  "FDR (BH)"   = "fdr",
                                                  "Bonferroni" = "bonferroni",
                                                  "Holm"       = "holm",
                                                  "BY"         = "BY"),
                                      selected = "none"))
              )
            ),
            checkboxInput("do_stepup",  "Stepup search (add edges)", value = FALSE),
            checkboxInput("do_modelsearch", "Full model search (BIC)", value = FALSE),

            hr(),
            actionButton("run_model", "Run Model",
                         icon = icon("play"), class = "btn-primary btn-lg",
                         style = "width:100%;")
          ),

          box(
            title = "Lambda Matrix (Factor Loadings Structure)", status = "warning",
            solidHeader = TRUE, width = 8,
            conditionalPanel(
              condition = "!['ggm','ising'].includes(input.model_family)",
              helpText("Define the factor loading structure. Each column is a latent factor.
                        Enter 1 for a free loading, 0 for fixed-to-zero. Use comma-separated
                        latent names below."),

              textInput("latent_names", "Latent variable names (comma-separated):",
                        placeholder = "e.g. Factor1, Factor2, Factor3"),

              actionButton("gen_lambda", "Generate Lambda Editor", icon = icon("table")),

              hr(),
              uiOutput("lambda_editor_ui"),

              hr(),
              h5("Alternatively: simple structure"),
              helpText("Automatically assigns each variable to one latent in order
                        (requires equal splits)."),
              actionButton("simple_structure", "Auto Simple Structure",
                           icon = icon("magic"))
            ),
            conditionalPanel(
              condition = "input.model_family == 'ggm'",
              helpText("The GGM does not require a lambda matrix. All selected variables
                        are treated as observed nodes in the network."),
              radioButtons("ggm_start", "Starting model:",
                           choices = c("Saturated (full)" = "full",
                                       "Empty (zero)"     = "zero"),
                           selected = "full")
            ),
            conditionalPanel(
              condition = "input.model_family == 'ising'",
              helpText("The Ising model does not require a lambda matrix. Select your
                        binary/dichotomous variables and fit directly.")
            )
          )
        ),

        # ── Panel model configuration (shown only for panel families) ────────
        conditionalPanel(
          condition = "['dlvm1','panelgvar','ri_clpm'].includes(input.model_family)",
          fluidRow(
            box(
              title = "Panel / Longitudinal Setup", status = "danger",
              solidHeader = TRUE, width = 12,
              fluidRow(
                # ── Column 4: Wave structure ────────────────────────────────
                column(4,
                  h4("Wave Structure"),
                  selectInput("wave_detect", "Auto-detect wave pattern:",
                              choices = c("Suffix _t1/_t2/… (e.g. dep_t1)" = "suffix_t",
                                          "Suffix _w1/_w2/… (e.g. dep_w1)" = "suffix_w",
                                          "Suffix _1/_2/…  (e.g. dep_1)"   = "suffix_n",
                                          "Manual assignment"               = "manual"),
                              selected = "suffix_t"),
                  uiOutput("wave_vars_ui"),
                  helpText(icon("info-circle"),
                           "Variables must be in wide format. Each wave should have the
                            same number of indicators in the same order.")
                ),
                # ── Column 4: Temporal structure ────────────────────────────
                column(4,
                  h4("Temporal (Beta) Structure"),
                  conditionalPanel(
                    condition = "input.model_family != 'panelgvar'",
                    checkboxInput("beta_full",
                                  "All cross-lagged paths free (full beta)",
                                  value = TRUE),
                    conditionalPanel(
                      condition = "input.beta_full == false",
                      helpText("Check = free path. Rows = outcome (wave t),
                                Cols = predictor (wave t−1)."),
                      uiOutput("beta_matrix_ui")
                    )
                  ),
                  conditionalPanel(
                    condition = "input.model_family == 'panelgvar'",
                    helpText("Panel GVAR estimates all temporal paths freely (full beta).
                              Use prune/stepup to select significant paths.")
                  ),
                  checkboxInput("panel_stationary",
                                "Constrain beta equal across waves (stationarity)",
                                value = TRUE)
                ),
                # ── Column 4: Within/between structure ─────────────────────
                column(4,
                  h4("Contemporaneous Structure"),
                  radioButtons("within_latent", "Within-person network:",
                               choices = c("GGM (partial correlations)" = "ggm",
                                           "Covariance"                  = "cov"),
                               selected = "ggm", inline = TRUE),
                  conditionalPanel(
                    condition = "input.model_family == 'dlvm1'",
                    radioButtons("between_latent", "Between-person structure:",
                                 choices = c("GGM"        = "ggm",
                                             "Covariance" = "cov",
                                             "None (zero)" = "zero"),
                                 selected = "cov", inline = TRUE)
                  )
                )
              )
            )
          )
        )
      ),

      # ── Tab 4: Results ────────────────────────────────────────────────────
      tabItem(
        tabName = "tab_results",
        fluidRow(
          tabBox(
            title = "Model Output", width = 12,
            tabPanel("Fit Indices",
                     uiOutput("transform_banner"),
                     verbatimTextOutput("fit_output")),
            tabPanel("Parameters",
                     DT::dataTableOutput("param_table")),
            tabPanel("Modification Indices",
                     DT::dataTableOutput("mi_table")),
            tabPanel("Matrices",
                     selectInput("matrix_select", "Select matrix:",
                                 choices = c("omega_zeta", "sigma_zeta", "lambda",
                                             "sigma_epsilon", "omega_epsilon",
                                             "beta", "omega", "sigma", "kappa",
                                             "── Relative Importance ──────────────" = "sep_ri",
                                             "RI Latent Network (raw)"        = "ri_raw",
                                             "RI Latent Network (normalized)" = "ri_norm")),
                     checkboxInput("threshold_matrix", "Threshold (non-significant → 0)",
                                   value = FALSE),
                     verbatimTextOutput("matrix_output")),
            tabPanel("Model Summary",
                     verbatimTextOutput("model_summary")),
            tabPanel("Compare Models",
                     uiOutput("transform_banner"),
                     helpText("After fitting multiple models, compare them here."),
                     verbatimTextOutput("compare_output"),
                     uiOutput("compare_fit_ui")),
            tabPanel("Factor Scores & Loadings",
              fluidRow(
                column(6,
                  h4("Factor Loadings (Lambda)"),
                  uiOutput("fs_group_ui"),
                  DT::dataTableOutput("loadings_table"),
                  br(),
                  downloadButton("download_loadings", "Download Loadings (CSV)",
                                 class = "btn-info btn-sm")
                ),
                column(6,
                  h4("Factor Scores"),
                  helpText("Regression-method scores: F = (X \u2212 \u03bc) \u00b7 \u03a3\u207b\u00b9 \u00b7 \u039b \u00b7 \u03a3\u03b6"),
                  DT::dataTableOutput("factor_scores_table"),
                  br(),
                  downloadButton("download_fscores", "Download Factor Scores (CSV)",
                                 class = "btn-info btn-sm"),
                  br(), br(),
                  downloadButton("download_fscores_merged",
                                 "Download Dataset + Factor Scores (CSV)",
                                 class = "btn-success btn-sm")
                )
              )
            )
          )
        )
      ),

      # ── Tab 5: Plots ──────────────────────────────────────────────────────
      tabItem(
        tabName = "tab_plots",
        fluidRow(
          box(
            title = "Plot Settings", status = "primary", solidHeader = TRUE,
            width = 3,

            selectInput("plot_type", "Plot type:",
                        choices = c(
                          "Factor Structure & Latent Network" = "factor_network",
                          "Factor Structure & Latent RI Network" = "factor_ri_network",
                          "Latent Network (omega_zeta)" = "latent_network",
                          "RI Latent Network" = "ri_latent_network",
                          "Residual Network (omega_epsilon)" = "residual_network",
                          "Full GGM (omega)" = "ggm_network",
                          "EGA \u2013 GGM (community detection)"   = "ega_ggm",
                          "EGA \u2013 Ising (community detection)" = "ega_ising",
                          "Temporal Network (beta)" = "beta_network",
                          "Within-Person Network" = "within_network",
                          "Between-Person Network" = "between_network",
                          "Correlation Heatmap" = "cor_heat",
                          "CI Plot" = "ci_plot",
                          "Path Diagram (semPlot)" = "semplot"
                        )),

            hr(),
            h4("qgraph options"),
            selectInput("qgraph_layout", "Layout:",
                        choices = c("spring", "circle", "groups"),
                        selected = "circle"),
            conditionalPanel(
              condition = "input.qgraph_layout != 'circle'",
              sliderInput("qgraph_repulsion", "Node repulsion:",
                          min = 0.1, max = 5, value = 1, step = 0.1,
                          ticks = FALSE)
            ),
            selectInput("qgraph_theme", "Theme:",
                        choices = c("colorblind", "classic", "gray",
                                    "Hollywood", "Borkulo", "TeamFortress",
                                    "Reddit", "Leuven", "Leuven2"),
                        selected = "colorblind"),
            sliderInput("qgraph_vsize", "Node size:", min = 2, max = 20,
                        value = 5, step = 1),
            sliderInput("qgraph_esize", "Max edge width:", min = 1, max = 20,
                        value = 5, step = 1),
            checkboxInput("qgraph_labels_est", "Show edge estimates", value = TRUE),
            selectInput("qgraph_label_scale", "Node label scaling:",
                        choices = c("Auto (fit to node)"  = "auto",
                                    "Equal (all same size)" = "equal",
                                    "None (fixed size)"     = "none"),
                        selected = "none"),
            sliderInput("qgraph_label_cex", "Label size:", min = 0.3, max = 3,
                        value = 1.5, step = 0.1),
            numericInput("qgraph_cutoff", "Edge cut-off (abs):", value = 0,
                         min = 0, max = 1, step = 0.01),
            sliderInput("edge_label_cex", "Edge label size:", min = 0.3, max = 3,
                        value = 1.0, step = 0.1),

            hr(),
            colourpicker::colourInput("node_color", "Node colour:", value = "#72AFD3",
                        showColour = "both", palette = "square"),
            colourpicker::colourInput("node_border_color", "Node border colour:", value = "#FFFFFF",
                        showColour = "both", palette = "square"),
            sliderInput("darken_min", "Darkness range (0 = very dark, 1 = no darkening):",
                        min = 0.1, max = 1.0, value = 0.75, step = 0.05),
            checkboxInput("disable_node_scaling", "Disable node scaling by severity", value = FALSE),
            uiOutput("latent_color_ui"),

            hr(),
            checkboxInput("show_predictability",
                          "Show node predictability (R\u00b2 rings)",
                          value = FALSE),
            conditionalPanel(
              condition = "input.show_predictability",
              colourpicker::colourInput("pred_ring_color", "Ring color:",
                                        value = "#ADD8E6",
                                        showColour = "both", palette = "square"),
              sliderInput("pred_ring_border",
                          "Ring thickness (0 = full pie, 1 = thin ring):",
                          min = 0, max = 1, value = 0.3, step = 0.05, ticks = FALSE)
            ),

            conditionalPanel(
              condition = "input.plot_type == 'ci_plot'",
              uiOutput("ci_matrix_ui"),
              checkboxInput("ci_split0", "Split zero / non-zero", value = TRUE),
              uiOutput("ci_source_ui")
            ),
            conditionalPanel(
              condition = "input.plot_type == 'ega_ggm' || input.plot_type == 'ega_ising'",
              selectInput("ega_algorithm", "Community detection algorithm:",
                          choices = c(
                            "Walktrap (default)" = "walktrap",
                            "Louvain"            = "louvain",
                            "Fast & Greedy"      = "fast_greedy",
                            "Spinglass"          = "spinglass"
                          ),
                          selected = "walktrap"),
              checkboxInput("ega_pca_layout",
                            "PCA layout \u2014 position nodes by PC1 / PC2",
                            value = FALSE),
              checkboxInput("ega_legend", "Show dimension legend", value = TRUE)
            ),
            conditionalPanel(
              condition = "input.plot_type == 'ri_latent_network' || input.plot_type == 'factor_ri_network'",
              checkboxInput("ri_normalize", "Normalize edges (proportion of R²)", value = FALSE)
            ),

            uiOutput("plot_group_ui"),

            hr(),
            sliderInput("plot_display_height", "Display height (px):",
                        min = 300, max = 1500, value = 650, step = 50),
            numericInput("plot_width",  "Download width (in):",  value = 10, min = 4, max = 24),
            numericInput("plot_height", "Download height (in):", value = 8,  min = 4, max = 20),
            downloadButton("download_plot", "Download Plot (PDF)", class = "btn-success")
          ),

          box(
            title = "Network Plot", status = "info", solidHeader = TRUE,
            width = 9,
            uiOutput("transform_banner"),
            uiOutput("plot_output_area")
          )
        )
      ),

      # ── Tab 6: Advanced ───────────────────────────────────────────────────
      tabItem(
        tabName = "tab_adv",
        fluidRow(
          box(
            title = "Model Comparison", status = "primary", solidHeader = TRUE,
            width = 12,
            fluidRow(
              # ── Column 4: Model 1 ───────────────────────────────────────────
              column(4,
                h4("Model 1 (Current)"),
                textInput("model1_label", "Label:", value = "Model1"),
                downloadButton("save_model1", "Save Model (.rds)", class = "btn-default btn-sm"),
                helpText("Save the current fitted model to reload later.")
              ),
              # ── Column 8: Model 2 ───────────────────────────────────────────
              column(8,
                h4("Model 2"),
                radioButtons("model2_source", "Source:",
                             choices = c("Fit new model" = "fit", "Load saved .rds" = "load"),
                             selected = "fit", inline = TRUE),
                # ── "Fit new" branch ─────────────────────────────────────────
                conditionalPanel(
                  condition = "input.model2_source == 'fit'",
                  radioButtons("model2_family", "Model family:",
                               choices = c("GGM" = "ggm", "CFA" = "lvm",
                                           "LNM" = "lnm", "RNM" = "rnm", "LRNM" = "lrnm"),
                               selected = "lvm", inline = TRUE),
                  conditionalPanel(
                    condition = "input.model2_family != 'ggm'",
                    radioButtons("model2_identification", "Identification:",
                                 choices = c("Loading" = "loadings", "Variance" = "variance"),
                                 selected = "variance", inline = TRUE),
                    helpText(icon("info-circle"),
                             "Use 'Variance' when comparing CFA vs LNM/RNM/LRNM.")
                  ),
                  checkboxInput("model2_override_prune",
                                "Override pruning settings for Model 2", value = FALSE),
                  conditionalPanel(
                    condition = "input.model2_override_prune == true",
                    checkboxInput("model2_do_prune", "Prune", value = FALSE),
                    conditionalPanel(
                      condition = "input.model2_do_prune == true",
                      fluidRow(
                        column(6, numericInput("model2_prune_alpha", "Prune α:",
                                               value = 0.01, min = 0.001, max = 0.2, step = 0.001)),
                        column(6, selectInput("model2_prune_adjust", "Adjust:",
                                              choices = c("None"       = "none",
                                                          "FDR (BH)"   = "fdr",
                                                          "Bonferroni" = "bonferroni",
                                                          "Holm"       = "holm",
                                                          "BY"         = "BY"),
                                              selected = "none"))
                      )
                    ),
                    checkboxInput("model2_do_stepup", "Step-up search", value = FALSE),
                    checkboxInput("model2_do_modelsearch", "Model search", value = FALSE)
                  ),
                  conditionalPanel(
                    condition = "input.group_var != null && input.group_var != ''",
                    checkboxInput("model2_constrain",
                                  "Constrain parameters across groups (invariance test)",
                                  value = FALSE),
                    conditionalPanel(
                      condition = "input.model2_constrain == true",
                      checkboxGroupInput("model2_groupEqual", "Constrain matrices:",
                                         choices = c("Loadings (lambda)"                   = "lambda",
                                                     "Latent covariances (omega_zeta)"      = "omega_zeta",
                                                     "Residual covariances (omega_epsilon)" = "omega_epsilon"),
                                         selected = "lambda")
                    )
                  )
                ),
                # ── "Load saved" branch ──────────────────────────────────────
                conditionalPanel(
                  condition = "input.model2_source == 'load'",
                  fileInput("load_model2_rds", "Load saved model (.rds):", accept = ".rds")
                ),
                # ── Label (always visible, auto-updated by server) ───────────
                textInput("model2_label", "Model 2 label:", value = "Model2")
              )
            ),
            fluidRow(
              column(12,
                actionButton("run_model2", "Fit & Compare", icon = icon("balance-scale"),
                             class = "btn-warning"),
                uiOutput("compare2_estimator_warn"),
                hr(),
                uiOutput("transform_banner"),
                verbatimTextOutput("compare2_output"),
                uiOutput("compare2_fit_ui")
              )
            )
          ),

          box(
            title = "Bootstrapping", status = "warning", solidHeader = TRUE,
            width = 6,
            numericInput("boot_reps", "Number of bootstrap replicates:", value = 100,
                         min = 10, max = 2000, step = 10),
            numericInput("boot_cores", "Cores:", value = 2, min = 1, max = 16),
            actionButton("run_boot", "Run Bootstrap", icon = icon("sync"),
                         class = "btn-info"),
            hr(),
            verbatimTextOutput("boot_output")
          )
        ),

        fluidRow(
          box(
            title = "Export R Code", status = "success", solidHeader = TRUE,
            width = 12,
            helpText("Copy the generated R code to reproduce this analysis in R / Positron."),
            verbatimTextOutput("r_code_output"),
            downloadButton("download_code", "Download R Script", class = "btn-success")
          )
        )
      ),

      # ── Tab About ─────────────────────────────────────────────────────────
      tabItem(
        tabName = "tab_about",
        fluidRow(
          box(
            title = "About this App", status = "info", solidHeader = TRUE, width = 12,
            h3("Psychonetrics Shiny App"),
            p("This application provides a graphical interface for the",
              tags$code("psychonetrics"), "R package (Epskamp, 2020+),
              with a focus on Latent Network Models."),
            h4("References"),
            tags$ul(
              tags$li("Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017).",
                      tags$em("Generalized network psychometrics: Combining network
                               and latent variable models."),
                      "Psychometrika, 82(4), 904-927."),
              tags$li("Epskamp, S. (2020+). psychonetrics: Structural Equation Modeling
                       and Confirmatory Network Analysis. R package.",
                      tags$a(href = "https://psychonetrics.org", "psychonetrics.org"))
            ),
            h4("Supported Models"),
            tags$ul(
              tags$li(tags$b("GGM"), " – Gaussian Graphical Model (pairwise partial correlations)"),
              tags$li(tags$b("CFA / SEM (lvm)"), " – Confirmatory Factor Analysis with covariance latent structure"),
              tags$li(tags$b("LNM"), " – Latent Network Model: latent covariance modeled as GGM"),
              tags$li(tags$b("RNM"), " – Residual Network Model: residual covariance modeled as GGM"),
              tags$li(tags$b("LRNM"), " – Combined Latent + Residual Network Model")
            ),
            h4("How to use"),
            tags$ol(
              tags$li("Upload your data (CSV, TXT, Excel, or SPSS)."),
              tags$li("Select variables, variable type, and estimator."),
              tags$li("Choose model family and define the lambda (loading) matrix."),
              tags$li("Run the model and inspect fit, parameters, and networks."),
              tags$li("Use the Plots tab for publication-ready qgraph figures."),
              tags$li("Use Advanced tab for model comparison and bootstrapping."),
              tags$li("For penalised ML estimation, select 'PML (penalised)' in the Model tab estimator dropdown.")
            )
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# ── 2. SERVER ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive values ──────────────────────────────────────────────────────────
  rv <- reactiveValues(
    raw_data      = NULL,   # original uploaded data
    analysis_data = NULL,   # processed (wide, selected vars)
    model         = NULL,   # fitted psychonetrics object
    model2        = NULL,   # second model for comparison
    model1_label  = "Model1",
    model2_label  = "Model2",
    boot_agg      = NULL,   # aggregated bootstrap
    lambda_mat    = NULL,   # lambda matrix
    latent_names  = NULL,
    vars_by_wave      = NULL,   # named list/matrix for panel models (wave → column names)
    transform_applied = NULL,   # label of data transformation currently in effect
    r_code        = "# Fit a model first to generate reproducible R code."
  )

  # ── Data Import ──────────────────────────────────────────────────────────
  observeEvent(list(input$file_upload, input$sep, input$dec, input$header), {
    req(input$file_upload)
    path <- input$file_upload$datapath
    ext  <- tools::file_ext(input$file_upload$name)

    df <- tryCatch({
      switch(tolower(ext),
        "csv" = read.csv(path, header = input$header, sep = input$sep,
                         dec = input$dec, stringsAsFactors = FALSE),
        "txt" = , "tsv" = read.delim(path, header = input$header, sep = input$sep,
                                      dec = input$dec, stringsAsFactors = FALSE),
        "xlsx" = , "xls" = readxl::read_excel(path, sheet = input$excel_sheet),
        "sav" = haven::read_sav(path),
        "por" = haven::read_por(path),
        stop("Unsupported file type: ", ext)
      )
    }, error = function(e) {
      showNotification(paste("Import error:", e$message), type = "error")
      NULL
    })

    if (!is.null(df)) {
      # Convert haven-labelled to plain R types
      df <- as.data.frame(lapply(df, function(x) {
        if (inherits(x, "haven_labelled")) as.numeric(x) else x
      }))
      rv$raw_data <- df

      # Update variable pickers
      all_vars <- names(df)
      numeric_vars <- names(df)[sapply(df, is.numeric)]

      updatePickerInput(session, "selected_vars",
                        choices = numeric_vars, selected = numeric_vars)
      updateSelectizeInput(session, "group_var",
                           choices = c("(none)" = "", all_vars),
                           selected = "")

      showNotification(paste("Loaded", nrow(df), "rows ×", ncol(df), "columns"),
                       type = "message")
    }
  })

  # ── Demo data: Star Wars ──────────────────────────────────────────────────
  observeEvent(input$load_demo, {
    e <- new.env()
    data("StarWars", package = "psychonetrics", envir = e)
    sw <- as.data.frame(e$StarWars)

    # Keep only the 10 survey items (Q1–Q10); drop demographics (Q11–Q13)
    item_vars <- paste0("Q", 1:10)
    df <- sw[, item_vars, drop = FALSE]

    rv$raw_data <- df

    updatePickerInput(session, "selected_vars",
                      choices = item_vars, selected = item_vars)
    updateSelectizeInput(session, "group_var",
                         choices = c("(none)" = "", item_vars),
                         selected = "")

    # Pre-configure 3-factor lambda:
    #   Q1 cross-loads on all three factors
    #   Q2–Q4  → Prequel
    #   Q5–Q7  → Originals
    #   Q8–Q10 → Sequel
    lats <- c("Prequel", "Originals", "Sequel")
    lam  <- matrix(0, nrow = 10, ncol = 3,
                   dimnames = list(item_vars, lats))
    lam[1, ]   <- 1        # Q1 cross-loads
    lam[2:4,  1] <- 1      # Q2–Q4 → Prequel
    lam[5:7,  2] <- 1      # Q5–Q7 → Originals
    lam[8:10, 3] <- 1      # Q8–Q10 → Sequel

    rv$lambda_mat   <- lam
    rv$latent_names <- lats

    updateTextInput(session, "latent_names",
                    value = paste(lats, collapse = ", "))

    showNotification(
      paste0("Star Wars demo loaded: ", nrow(df),
             " respondents × 10 items. Lambda pre-configured (3 factors)."),
      type = "message", duration = 6
    )
    updateTabItems(session, "tabs", "tab_vars")
  })

  # Data summary
  output$data_summary <- renderPrint({
    req(rv$raw_data)
    cat("Dimensions:", nrow(rv$raw_data), "rows ×", ncol(rv$raw_data), "columns\n")
    cat("Missing values:", sum(is.na(rv$raw_data)), "\n\n")
    str(rv$raw_data, give.attr = FALSE)
  })

  output$data_table <- DT::renderDataTable({
    req(rv$raw_data)
    DT::datatable(head(rv$raw_data, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ── Prepare analysis data ────────────────────────────────────────────────
  analysis_data <- reactive({
    req(rv$raw_data, input$selected_vars)
    df <- rv$raw_data

    # Long → wide pivot if needed
    if (input$data_format == "long") {
      req(input$id_var, input$time_var)
      req(input$id_var %in% names(df), input$time_var %in% names(df))
      # Simple reshape: value vars are selected_vars
      val_vars <- intersect(input$selected_vars, names(df))
      df_long <- df[, c(input$id_var, input$time_var, val_vars), drop = FALSE]
      df <- tryCatch({
        reshape(as.data.frame(df_long),
                idvar = input$id_var, timevar = input$time_var,
                direction = "wide")
      }, error = function(e) {
        showNotification(paste("Pivot error:", e$message), type = "error")
        NULL
      })
      if (is.null(df)) return(NULL)
    }

    # Select only the requested variables (may include group var)
    keep <- intersect(input$selected_vars, names(df))
    if (length(keep) < 2) {
      showNotification("Select at least 2 variables.", type = "warning")
      return(NULL)
    }
    out <- df[, keep, drop = FALSE]

    # Convert to numeric where needed
    out <- as.data.frame(lapply(out, function(x) {
      if (is.factor(x)) as.numeric(as.character(x))
      else if (is.character(x)) suppressWarnings(as.numeric(x))
      else as.numeric(x)
    }))

    # ── Data transformation ──────────────────────────────────────────────────
    transform_choice <- if (!is.null(input$data_transform)) input$data_transform else "none"
    data_type_choice <- if (!is.null(input$data_type)) input$data_type else "continuous"

    # Dichotomization is allowed for continuous + ordered types.
    # npn/log/sqrt/zscore only apply to continuous. Already-binary data is skipped.
    is_dichot_transform <- transform_choice %in% c("dichot_mean", "dichot_median")
    skip_transform <- transform_choice == "none" ||
                      data_type_choice == "dichotomous" ||
                      (!is_dichot_transform && data_type_choice != "continuous")

    if (!skip_transform) {
      out <- tryCatch({
        m <- as.matrix(out)
        m_tr <- switch(
          transform_choice,
          npn = {
            # Truncated normal-score transformation (Liu et al. 2009)
            apply(m, 2, function(x) {
              n     <- sum(!is.na(x))
              delta <- 1 / (4 * n^0.25 * sqrt(pi * log(n)))
              r     <- rank(x, ties.method = "average", na.last = "keep")
              qnorm(pmin(pmax(r / n, delta), 1 - delta))
            })
          },
          dichot_mean = {
            apply(m, 2, function(x) as.numeric(x >= mean(x, na.rm = TRUE)))
          },
          dichot_median = {
            apply(m, 2, function(x) as.numeric(x >= median(x, na.rm = TRUE)))
          },
          log1p  = { log1p(pmax(m, 0)) },
          sqrt   = { sqrt(pmax(m, 0)) },
          zscore = { as.matrix(scale(m)) },
          m
        )
        colnames(m_tr) <- colnames(m)
        as.data.frame(m_tr)
      }, error = function(e) {
        showNotification(
          paste("Transformation failed:", e$message, "— using untransformed data."),
          type = "warning", duration = 6
        )
        transform_choice <<- "none"
        out
      })
    }

    # Record which transformation is in effect
    rv$transform_applied <- if (transform_choice == "none" || skip_transform) {
      if (data_type_choice == "dichotomous" && transform_choice != "none")
        "None (data already binary — transformation skipped)"
      else if (!is_dichot_transform && data_type_choice != "continuous" && transform_choice != "none")
        paste0("None (data type: ", data_type_choice, " — transformation skipped)")
      else
        NULL
    } else {
      switch(transform_choice,
        npn           = "Nonparanormal (rank-based normal scores)",
        dichot_mean   = "Dichotomized at mean (\u2192 0/1)",
        dichot_median = "Dichotomized at median (\u2192 0/1)",
        log1p         = "Log (log1p)",
        sqrt          = "Square root",
        zscore        = "Z-score (standardized)",
        NULL
      )
    }

    rv$analysis_data <- out
    out
  })

  # ── Descriptives ─────────────────────────────────────────────────────────
  output$desc_table <- DT::renderDataTable({
    df <- analysis_data()
    req(df)
    desc <- data.frame(
      Variable = names(df),
      N        = sapply(df, function(x) sum(!is.na(x))),
      Mean     = round(sapply(df, mean, na.rm = TRUE), 3),
      SD       = round(sapply(df, sd,   na.rm = TRUE), 3),
      Min      = sapply(df, min, na.rm = TRUE),
      Max      = sapply(df, max, na.rm = TRUE),
      Missing  = sapply(df, function(x) sum(is.na(x)))
    )
    DT::datatable(desc, rownames = FALSE, options = list(dom = "t", pageLength = 50))
  })

  output$cor_heatmap <- renderPlot({
    df <- analysis_data()
    req(df, ncol(df) >= 2)
    cor_mat <- cor(df, use = "pairwise.complete.obs")
    melted <- as.data.frame(as.table(cor_mat))
    names(melted) <- c("Var1", "Var2", "value")
    ggplot(melted, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), size = 3) +
      scale_fill_viridis(limits = c(-1, 1), option = "D", name = "r") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Correlation Matrix", x = "", y = "")
  })

  # ── Lambda matrix editor ─────────────────────────────────────────────────
  observeEvent(input$gen_lambda, {
    req(input$selected_vars, input$latent_names)
    vars    <- input$selected_vars
    latents <- trimws(strsplit(input$latent_names, ",")[[1]])
    latents <- latents[latents != ""]
    if (length(latents) < 1) {
      showNotification("Enter at least one latent name.", type = "warning")
      return()
    }
    rv$latent_names <- latents

    # Default: zeros
    lam <- matrix(0, nrow = length(vars), ncol = length(latents),
                  dimnames = list(vars, latents))
    rv$lambda_mat <- lam
  })

  output$lambda_editor_ui <- renderUI({
    req(rv$lambda_mat)
    vars    <- rownames(rv$lambda_mat)
    latents <- colnames(rv$lambda_mat)

    header <- tags$tr(
      tags$th("Variable"),
      lapply(latents, function(l) tags$th(l, style = "text-align:center;"))
    )

    rows <- lapply(seq_along(vars), function(i) {
      tags$tr(
        tags$td(tags$b(vars[i])),
        lapply(seq_along(latents), function(j) {
          id <- paste0("lam_", i, "_", j)
          tags$td(
            materialSwitch(inputId = id, label = NULL,
                           value  = as.logical(rv$lambda_mat[i, j]),
                           status = "success"),
            style = "padding: 2px 5px; text-align: center; vertical-align: middle;"
          )
        })
      )
    })

    tags$div(
      style = "overflow-x: auto;",
      tags$table(
        class = "table table-condensed table-bordered",
        tags$thead(header),
        tags$tbody(rows)
      )
    )
  })

  # Read lambda from UI inputs
  get_lambda <- reactive({
    req(rv$lambda_mat)
    vars    <- rownames(rv$lambda_mat)
    latents <- colnames(rv$lambda_mat)
    lam <- matrix(0, nrow = length(vars), ncol = length(latents),
                  dimnames = list(vars, latents))
    for (i in seq_along(vars)) {
      for (j in seq_along(latents)) {
        id  <- paste0("lam_", i, "_", j)
        val <- input[[id]]
        if (!is.null(val)) lam[i, j] <- if (isTRUE(val)) 1L else 0L
      }
    }
    lam
  })

  # Auto simple structure
  observeEvent(input$simple_structure, {
    req(input$selected_vars, input$latent_names)
    vars    <- input$selected_vars
    latents <- trimws(strsplit(input$latent_names, ",")[[1]])
    latents <- latents[latents != ""]
    n_lat   <- length(latents)
    n_var   <- length(vars)

    if (n_lat < 1) {
      showNotification("Enter latent names first.", type = "warning")
      return()
    }

    lam <- matrix(0, n_var, n_lat, dimnames = list(vars, latents))
    items_per <- ceiling(n_var / n_lat)
    for (j in seq_len(n_lat)) {
      rows <- ((j - 1) * items_per + 1):min(j * items_per, n_var)
      lam[rows, j] <- 1
    }
    rv$lambda_mat   <- lam
    rv$latent_names <- latents
    showNotification("Simple structure applied.", type = "message")
  })

  # ── Panel helpers ─────────────────────────────────────────────────────────

  # Detect waves from column names and return named list: wave_id → variable names.
  # Returns NULL if the pattern cannot be applied.
  detect_waves <- function(vars, pattern = "suffix_t") {
    suffix_re <- switch(pattern,
      suffix_t = "_t(\\d+)$",
      suffix_w = "_w(\\d+)$",
      suffix_n = "_(\\d+)$",
      NULL
    )
    if (is.null(suffix_re)) return(NULL)
    m <- regmatches(vars, regexpr(suffix_re, vars, perl = TRUE))
    if (length(m) == 0 || any(nchar(m) == 0)) return(NULL)
    nums  <- as.integer(gsub("[^0-9]", "", m))
    bases <- sub(suffix_re, "", vars, perl = TRUE)
    if (anyDuplicated(bases[nums == nums[1]])) return(NULL)  # ambiguous
    wave_ids <- sort(unique(nums))
    result   <- lapply(wave_ids, function(w) vars[nums == w])
    names(result) <- paste0("wave", wave_ids)
    result
  }

  # Build vars matrix (indicators × waves) from a named list.
  waves_to_matrix <- function(wave_list) {
    lens <- sapply(wave_list, length)
    if (length(unique(lens)) != 1) return(wave_list)  # unequal sizes: return list
    do.call(cbind, wave_list)  # character matrix: rows=indicators, cols=waves
  }

  # Build beta matrix from checkbox inputs (or return "full").
  get_beta <- reactive({
    if (isTRUE(input$beta_full) || is.null(rv$latent_names)) return("full")
    lats <- rv$latent_names
    q    <- length(lats)
    m    <- matrix(0L, q, q, dimnames = list(lats, lats))
    for (i in seq_len(q)) {
      for (j in seq_len(q)) {
        val <- input[[paste0("beta_", i, "_", j)]]
        if (!is.null(val)) m[i, j] <- if (isTRUE(val)) 1L else 0L
      }
    }
    m
  })

  # ── Wave mapping: auto-detect or manual ────────────────────────────────────
  output$wave_vars_ui <- renderUI({
    req(input$selected_vars, length(input$selected_vars) >= 2)
    vars <- input$selected_vars

    if (input$wave_detect == "manual") {
      tagList(
        numericInput("n_waves_manual", "Number of waves:", value = 2, min = 2, max = 10),
        uiOutput("wave_manual_selects")
      )
    } else {
      wave_list <- detect_waves(vars, input$wave_detect)
      if (is.null(wave_list) || length(wave_list) < 2) {
        helpText(icon("exclamation-triangle"),
                 "Could not auto-detect waves from column names.
                  Try a different suffix pattern or use 'Manual assignment'.")
      } else {
        n_ind <- length(wave_list[[1]])
        tagList(
          tags$p(tags$b(sprintf("%d waves detected, %d indicators each:",
                                length(wave_list), n_ind))),
          tags$ul(lapply(names(wave_list), function(nm)
            tags$li(tags$b(nm, ": "), paste(wave_list[[nm]], collapse = ", "))
          ))
        )
      }
    }
  })

  # Manual wave assignment: one pickerInput per wave
  output$wave_manual_selects <- renderUI({
    req(input$selected_vars, input$n_waves_manual)
    vars   <- input$selected_vars
    n_waves <- as.integer(input$n_waves_manual)
    lapply(seq_len(n_waves), function(k) {
      pickerInput(paste0("wave_vars_", k),
                  label    = paste0("Wave ", k, " variables:"),
                  choices  = vars,
                  selected = vars,
                  multiple = TRUE,
                  options  = pickerOptions(actionsBox = TRUE, liveSearch = TRUE,
                                           title = "Select…"))
    })
  })

  # Store rv$vars_by_wave whenever wave inputs change
  observe({
    req(input$model_family %in% c("dlvm1", "panelgvar", "ri_clpm"))
    req(input$selected_vars, length(input$selected_vars) >= 2)

    if (isTRUE(input$wave_detect == "manual")) {
      req(input$n_waves_manual)
      n_waves <- as.integer(input$n_waves_manual)
      wave_list <- setNames(
        lapply(seq_len(n_waves), function(k) {
          v <- input[[paste0("wave_vars_", k)]]
          if (is.null(v)) character(0) else v
        }),
        paste0("wave", seq_len(n_waves))
      )
      if (any(sapply(wave_list, length) == 0)) return()
      rv$vars_by_wave <- waves_to_matrix(wave_list)
    } else {
      wave_list <- detect_waves(input$selected_vars, input$wave_detect)
      if (!is.null(wave_list) && length(wave_list) >= 2) {
        rv$vars_by_wave <- waves_to_matrix(wave_list)
      }
    }
  })

  # ── Beta matrix editor (q × q toggles) ─────────────────────────────────────
  output$beta_matrix_ui <- renderUI({
    req(rv$latent_names, length(rv$latent_names) >= 1)
    lats <- rv$latent_names
    q    <- length(lats)

    header <- tags$tr(
      tags$th("Outcome \\ Predictor"),
      lapply(lats, function(l) tags$th(l, style = "text-align:center;"))
    )
    rows <- lapply(seq_len(q), function(i) {
      tags$tr(
        tags$td(tags$b(lats[i])),
        lapply(seq_len(q), function(j) {
          id <- paste0("beta_", i, "_", j)
          tags$td(
            materialSwitch(inputId = id, label = NULL, value = TRUE, status = "warning"),
            style = "padding: 2px 5px; text-align: center; vertical-align: middle;"
          )
        })
      )
    })

    tags$div(
      style = "overflow-x: auto;",
      tags$table(
        class = "table table-condensed table-bordered",
        tags$thead(header),
        tags$tbody(rows)
      )
    )
  })

  # ── Run Model ────────────────────────────────────────────────────────────
  observeEvent(input$run_model, {
    df <- analysis_data()
    req(df)

    # Warn if missing data is present
    n_total    <- nrow(df)
    n_complete <- sum(complete.cases(df))
    if (n_complete < n_total) {
      n_dropped <- n_total - n_complete
      pct       <- round(100 * n_dropped / n_total)
      is_pairwise <- isTRUE(input$missing_method == "pairwise") &&
                     isTRUE(input$model_family == "ggm") &&
                     !(input$data_type %in% c("ordered", "dichotomous"))
      method_note <- if (is_pairwise) "pairwise deletion (all available pairs used)"
                     else paste0("listwise deletion. nobs = ", n_complete)
      showNotification(
        paste0("Missing data: ", n_dropped, " of ", n_total, " rows (", pct,
               "%) have incomplete cases — using ", method_note, "."),
        type     = "warning",
        duration = 10
      )
    }

    # Determine ordered variables
    ordered_vars <- if (input$data_type %in% c("ordered", "dichotomous")) {
      input$selected_vars
    } else {
      character(0)
    }

    est <- input$estimator_type

    # Group variable
    grp <- NULL
    if (!is.null(input$group_var) && input$group_var != "") {
      gv <- rv$raw_data[[input$group_var]]
      if (!is.null(gv)) {
        df[[input$group_var]] <- gv[seq_len(nrow(df))]
        # Drop rows with NA in the group variable — NA group membership
        # produces NA in the optimizer's parameter vector and crashes it
        # with "missing value where TRUE/FALSE needed"
        n_before <- nrow(df)
        df <- df[!is.na(df[[input$group_var]]), ]
        n_dropped <- n_before - nrow(df)
        if (n_dropped > 0) {
          showNotification(
            paste0("Group variable '", input$group_var, "': dropped ", n_dropped,
                   " row(s) with missing group assignment."),
            type = "warning", duration = 8)
        }
        grp <- input$group_var
      }
    }

    withProgress(message = "Fitting model…", value = 0.2, {

      mod <- tryCatch({
        if (input$model_family == "ggm") {
          # ── GGM ──
          omega_start <- if (input$ggm_start == "full") "full" else "zero"
          use_pairwise <- isTRUE(input$missing_method == "pairwise") &&
                          length(ordered_vars) == 0  # pairwise only for continuous
          if (use_pairwise) {
            vars_ <- input$selected_vars
            if (!is.null(grp) && grp %in% names(df)) {
              # Multi-group: one covariance matrix per group
              grp_levels <- sort(unique(df[[grp]]))
              covs_list  <- lapply(grp_levels, function(gl) {
                sub <- df[df[[grp]] == gl, vars_, drop = FALSE]
                cov(sub, use = "pairwise.complete.obs")
              })
              nobs_vec <- sapply(grp_levels, function(gl) {
                sub <- df[df[[grp]] == gl, vars_, drop = FALSE]
                min(crossprod(!is.na(as.matrix(sub))))
              })
              m <- ggm(covs = covs_list, nobs = as.integer(nobs_vec),
                       vars = vars_, omega = omega_start, estimator = est)
            } else {
              sub      <- df[, vars_, drop = FALSE]
              covs_mat <- cov(sub, use = "pairwise.complete.obs")
              nobs_val <- min(crossprod(!is.na(as.matrix(sub))))
              m <- ggm(covs = covs_mat, nobs = as.integer(nobs_val),
                       vars = vars_, omega = omega_start, estimator = est)
            }
          } else {
            m <- ggm(df, vars = input$selected_vars, omega = omega_start,
                     ordered = ordered_vars, estimator = est,
                     groupvar = grp)
          }
          m <- do_runmodel(m)
          if (input$do_prune)       m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          if (input$do_stepup)      m <- m %>% stepup
          if (input$do_modelsearch) m <- m %>% modelsearch
          m

        } else if (input$model_family == "ising") {
          # ── Ising model (binary data) ─────────────────────────────────────
          m <- Ising(df, vars = input$selected_vars, estimator = est)
          m <- do_runmodel(m)
          if (input$do_prune)       m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          if (input$do_stepup)      m <- m %>% stepup
          if (input$do_modelsearch) m <- m %>% modelsearch
          m

        } else if (input$model_family %in% c("lvm", "lnm", "rnm", "lrnm")) {
          # ── LVM family ──
          lambda <- get_lambda()
          req(lambda)
          latents <- colnames(lambda)
          vars_   <- rownames(lambda)

          # ── Diagnostic: catch name mismatches before psychonetrics sees them ──
          missing_from_df  <- setdiff(vars_, names(df))
          missing_from_lam <- setdiff(names(df), vars_)
          if (length(missing_from_df) > 0) {
            showNotification(
              paste0("Variable name mismatch — in lambda but NOT in data: ",
                     paste(missing_from_df, collapse = ", ")),
              type = "error", duration = 15)
          }
          message(sprintf(
            "[PsychoNetrix] data: %d rows × %d cols | complete cases: %d | vars_ (%d): %s",
            nrow(df), ncol(df), sum(complete.cases(df[, intersect(vars_, names(df))])),
            length(vars_), paste(vars_, collapse = ", ")
          ))

          builder <- switch(input$model_family,
            "lvm"  = function(...) lvm(...),
            "lnm"  = function(...) lnm(...),
            "rnm"  = function(...) rnm(..., omega_epsilon = "full"),
            "lrnm" = function(...) lrnm(..., omega_epsilon = "full")
          )

          m <- builder(
            data           = df,
            lambda         = lambda,
            vars           = vars_,
            latents        = latents,
            identification = input$identification,
            ordered        = ordered_vars,
            estimator      = est,
            groupvar       = grp
          )

          m <- do_runmodel(m)
          if (input$do_prune)       m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          if (input$do_stepup)      m <- m %>% stepup
          if (input$do_modelsearch) m <- m %>% modelsearch
          m

        } else if (input$model_family == "dlvm1") {
          # ── Panel DLVM ──────────────────────────────────────────────────
          req(rv$vars_by_wave)
          req(rv$latent_names)
          lam      <- get_lambda()
          req(lam)
          beta_arg <- get_beta()
          equal_arg <- if (isTRUE(input$panel_stationary)) "beta" else "none"
          m <- dlvm1(
            data           = df,
            vars           = rv$vars_by_wave,
            lambda         = lam,
            latents        = rv$latent_names,
            beta           = beta_arg,
            within_latent  = input$within_latent,
            between_latent = input$between_latent,
            identification = input$identification,
            estimator      = est,
            equal          = equal_arg
          )
          m <- do_runmodel(m)
          if (input$do_prune)       m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          if (input$do_stepup)      m <- m %>% stepup
          if (input$do_modelsearch) m <- m %>% modelsearch
          m

        } else if (input$model_family == "panelgvar") {
          # ── Panel GVAR (observed) ────────────────────────────────────────
          req(rv$vars_by_wave)
          m <- panelgvar(
            data          = df,
            vars          = rv$vars_by_wave,
            within_latent = input$within_latent,
            estimator     = est
          )
          m <- do_runmodel(m)
          if (input$do_prune)       m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          if (input$do_stepup)      m <- m %>% stepup
          if (input$do_modelsearch) m <- m %>% modelsearch
          m

        } else if (input$model_family == "ri_clpm") {
          # ── RI-CLPM ──────────────────────────────────────────────────────
          req(rv$vars_by_wave)
          lam_arg <- if (!is.null(rv$latent_names) && length(rv$latent_names) > 0 &&
                         !is.null(rv$lambda_mat))
                       get_lambda() else NULL
          m_args <- list(
            data      = df,
            vars      = rv$vars_by_wave,
            type      = input$within_latent,
            estimator = est
          )
          if (!is.null(lam_arg)) m_args$lambda <- lam_arg
          m <- do.call(ri_clpm, m_args)
          m <- do_runmodel(m)
          if (input$do_prune) m <- m %>% prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          m
        }
      }, error = function(e) {
        showNotification(paste("Model error:", e$message), type = "error", duration = 10)
        NULL
      })

      setProgress(1)
    })

    if (!is.null(mod)) {
      rv$model        <- mod
      rv$model1_label <- toupper(input$model_family)
      rv$r_code       <- build_r_code()
      updateTextInput(session, "model1_label", value = rv$model1_label)
      showNotification("Model fitted successfully!", type = "message")
      updateTabItems(session, "tabs", "tab_results")
    }
  })

  # ── Data transform banner (shown in every result area) ───────────────────
  output$transform_banner <- renderUI({
    lbl <- rv$transform_applied
    if (is.null(lbl)) {
      div(class = "alert alert-info",
          style = "padding:6px 10px; margin-bottom:8px;",
          icon("database"), " Data used: raw (no transformation applied).")
    } else {
      div(class = "alert alert-success",
          style = "padding:6px 10px; margin-bottom:8px;",
          icon("exchange-alt"), paste(" Data used: transformed \u2014", lbl))
    }
  })

  # ── CI plot: dynamic matrix selector ─────────────────────────────────────
  output$ci_matrix_ui <- renderUI({
    req(rv$model)
    avail <- tryCatch(
      unique(as.character(parameters(rv$model)$matrix)),
      error = function(e) NULL
    )
    if (is.null(avail) || length(avail) == 0) {
      avail <- c("omega_zeta", "omega_epsilon", "lambda", "omega",
                 "sigma_zeta", "sigma_epsilon", "beta",
                 "omega_zeta_within", "omega_zeta_between", "kappa")
    }
    default <- if ("omega_zeta" %in% avail) "omega_zeta" else avail[1]
    selectInput("ci_matrix", "Matrix for CI plot:", choices = avail, selected = default)
  })

  # ── CI plot: model-based vs bootstrap toggle ──────────────────────────────
  output$ci_source_ui <- renderUI({
    req(rv$model)
    has_boot <- !is.null(rv$boot_agg)
    choices  <- c("Model-based (SE)" = "model")
    if (has_boot) choices <- c(choices, "Bootstrap (percentile)" = "boot")
    radioButtons("ci_source", "CI source:", choices = choices,
                 selected = "model", inline = TRUE)
  })

  # ── Fit index interpretation guide ───────────────────────────────────────
  # Thresholds drawn from:
  #   Hu & Bentler (1999) Struct. Equ. Model. 6(1):1-55
  #   Browne & Cudeck (1992) Sociol. Methods Res. 21:230-258
  #   MacCallum, Browne & Sugawara (1996) Psychol. Methods 1(2):130-149
  #   Schermelleh-Engel, Moosbrugger & Müller (2003) MPR-Online 8(2):23-74
  #   Kline (2016) Principles & Practice of SEM, 4th ed. Guilford
  .fit_thresholds <- list(
    CFI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    TLI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    NNFI  = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    NFI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    IFI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    RFI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    RNI   = list(dir = "high", cuts = c(.90, .95, .97),
                 labs = c("Poor (<.90)", "Acceptable (.90-.94)", "Good (.95-.96)", "Excellent (≥.97)")),
    RMSEA = list(dir = "low",  cuts = c(.05, .08, .10),
                 labs = c("Excellent (≤.05)", "Good (.05-.08)", "Acceptable (.08-.10)", "Poor (>.10)")),
    SRMR  = list(dir = "low",  cuts = c(.05, .08, .10),
                 labs = c("Excellent (≤.05)", "Good (.05-.08)", "Acceptable (.08-.10)", "Poor (>.10)")),
    GFI   = list(dir = "high", cuts = c(.85, .90, .95),
                 labs = c("Poor (<.85)", "Acceptable (.85-.89)", "Good (.90-.94)", "Excellent (≥.95)"))
  )

  .print_fit_guide <- function(ft) {
    if (is.null(ft) || !is.data.frame(ft)) return(invisible(NULL))

    get_val <- function(nm) {
      idx <- which(toupper(ft$Measure) == toupper(nm))
      if (length(idx) == 0) return(NA_real_)
      suppressWarnings(as.numeric(ft$Value[idx[1]]))
    }

    classify <- function(val, spec) {
      if (is.na(val)) return(NA_character_)
      cuts <- spec$cuts; labs <- spec$labs
      if (spec$dir == "high") {
        if (val >= cuts[3]) labs[4] else if (val >= cuts[2]) labs[3] else if (val >= cuts[1]) labs[2] else labs[1]
      } else {
        if (val <= cuts[1]) labs[1] else if (val <= cuts[2]) labs[2] else if (val <= cuts[3]) labs[3] else labs[4]
      }
    }

    verdict_str <- function(val, nm) {
      spec <- .fit_thresholds[[toupper(nm)]]
      if (is.null(spec) || is.na(val)) return("")
      label <- classify(val, spec)
      if (is.na(label)) return("")
      tag <- if (grepl("^Excellent", label)) "✔✔"
             else if (grepl("^Good",       label)) "✔ "
             else if (grepl("^Acceptable", label)) "~ "
             else "✖ "
      paste0("  ", tag, " ", label)
    }

    # Print one index row (with optional scaled variant on next line)
    print_idx <- function(nm, label = nm) {
      v <- get_val(nm)
      if (is.na(v)) return(invisible(NULL))
      cat(sprintf("  %-18s %7.3f%s\n", label, v, verdict_str(v, nm)))
      # scaled variant
      v_sc <- get_val(paste0(nm, ".scaled"))
      if (!is.na(v_sc))
        cat(sprintf("  %-18s %7.3f%s\n", paste0(label, " (scaled)"), v_sc, verdict_str(v_sc, nm)))
    }

    SEP <- "  ──────────────────────────────────────────────────────────────\n"

    cat("══════════════════════════════════════════════════════════════════\n")
    cat("  Model Fit Summary\n")
    cat("══════════════════════════════════════════════════════════════════\n")

    # ── Model info ──
    nvar_v <- get_val("nvar"); npar_v <- get_val("npar"); df_v <- get_val("df")
    info_parts <- c(
      if (!is.na(nvar_v)) sprintf("Variables: %g", nvar_v),
      if (!is.na(npar_v)) sprintf("Free parameters: %g", npar_v),
      if (!is.na(df_v))   sprintf("df: %g", df_v)
    )
    if (length(info_parts) > 0) cat("  ", paste(info_parts, collapse = "   "), "\n", sep = "")

    # ── Chi-square ──
    chisq_v <- get_val("chisq"); pv <- get_val("pvalue")
    chisq_sc <- get_val("chisq.scaled"); pv_sc <- get_val("pvalue.scaled")
    if (!is.na(chisq_v) || !is.na(chisq_sc)) {
      cat("\n  Chi-square test of exact fit\n")
      cat(SEP)
      pval_str <- function(p) {
        if (is.na(p)) "" else if (p < .001) "  p < .001" else sprintf("  p = %.3f", p)
      }
      if (!is.na(chisq_v))
        cat(sprintf("  %-18s %7.2f%s\n", sprintf("chi2 (df=%g)", df_v), chisq_v, pval_str(pv)))
      if (!is.na(chisq_sc))
        cat(sprintf("  %-18s %7.2f%s\n", sprintf("chi2.scaled (df=%g)", df_v), chisq_sc, pval_str(pv_sc)))
      cat("  (Significant chi2 is expected for large N; rely on indices below)\n")
    }

    # ── Absolute fit ──
    has_abs <- any(!is.na(c(get_val("rmsea"), get_val("srmr"), get_val("gfi"))))
    if (has_abs) {
      cat("\n  Absolute fit\n")
      cat(SEP)
      # RMSEA with CI
      v_rmsea <- get_val("rmsea")
      if (!is.na(v_rmsea)) {
        lb <- get_val("rmsea.ci.lower"); ub <- get_val("rmsea.ci.upper")
        ci_str <- if (!is.na(lb) && !is.na(ub)) sprintf("  90%% CI [%.3f, %.3f]", lb, ub) else ""
        cat(sprintf("  %-18s %7.3f%s%s\n", "RMSEA", v_rmsea, ci_str, verdict_str(v_rmsea, "RMSEA")))
        v_rmsea_sc <- get_val("rmsea.scaled")
        if (!is.na(v_rmsea_sc))
          cat(sprintf("  %-18s %7.3f%s\n", "RMSEA (scaled)", v_rmsea_sc, verdict_str(v_rmsea_sc, "RMSEA")))
      }
      print_idx("srmr",  "SRMR")
      print_idx("gfi",   "GFI")
    }

    # ── Incremental fit ──
    inc_nms  <- c("cfi","tli","nnfi","nfi","ifi","rfi","rni")
    inc_labs <- c("CFI","TLI / NNFI","NNFI","NFI","IFI","RFI","RNI")
    has_inc  <- any(sapply(inc_nms, function(nm) !is.na(get_val(nm))))
    if (has_inc) {
      cat("\n  Incremental fit\n")
      cat(SEP)
      shown <- character(0)
      for (i in seq_along(inc_nms)) {
        nm <- inc_nms[i]
        if (nm %in% shown) next           # skip NNFI if already shown as TLI alias
        v  <- get_val(nm)
        if (is.na(v)) next
        # TLI and NNFI are often identical — merge label if so
        if (nm == "tli") {
          v_nnfi <- get_val("nnfi")
          lab    <- if (!is.na(v_nnfi) && abs(v - v_nnfi) < 1e-6) { shown <- c(shown, "nnfi"); "TLI / NNFI" }
                    else "TLI"
        } else {
          lab <- inc_labs[i]
        }
        cat(sprintf("  %-18s %7.3f%s\n", lab, v, verdict_str(v, nm)))
        v_sc <- get_val(paste0(nm, ".scaled"))
        if (!is.na(v_sc))
          cat(sprintf("  %-18s %7.3f%s\n", paste0(lab, " (scaled)"), v_sc, verdict_str(v_sc, nm)))
      }
    }

    # ── Information criteria ──
    ic_pairs <- list(c("aic.x","AIC"), c("aic.x2","AIC2"), c("bic","BIC"), c("bic2","BIC2"),
                     c("ebic.25","eBIC(.25)"), c("ebic.5","eBIC(.5)"), c("ebic.75","eBIC(.75)"),
                     c("ebic1","eBIC(1)"))
    ic_vals  <- lapply(ic_pairs, function(p) list(nm = p[1], lab = p[2], v = get_val(p[1])))
    ic_vals  <- Filter(function(x) !is.na(x$v), ic_vals)
    if (length(ic_vals) > 0) {
      cat("\n  Information criteria  (lower = better; use for model comparison)\n")
      cat(SEP)
      for (x in ic_vals) cat(sprintf("  %-18s %10.2f\n", x$lab, x$v))
    }

    # ── Legend & references ──
    cat("\n")
    cat(SEP)
    cat("  Legend: ✔✔ Excellent  ✔ Good  ~ Acceptable  ✖ Poor\n")
    cat(SEP)
    cat("  References:\n")
    cat("  [1] Hu & Bentler (1999) Struct. Equ. Model. 6(1):1-55\n")
    cat("      CFI/TLI ≥.95 good, SRMR ≤.08, RMSEA ≤.06\n")
    cat("  [2] Browne & Cudeck (1992) Sociol. Methods Res. 21:230-258\n")
    cat("      RMSEA: ≤.05 close fit, .05-.08 reasonable, .08-.10 mediocre, >.10 poor\n")
    cat("  [3] MacCallum, Browne & Sugawara (1996) Psychol. Methods 1(2):130-149\n")
    cat("      RMSEA confidence-interval interpretation\n")
    cat("  [4] Schermelleh-Engel et al. (2003) MPR-Online 8(2):23-74\n")
    cat("      CFI/TLI ≥.97 excellent; GFI ≥.95 excellent\n")
    cat("  [5] Kline (2016) Principles & Practice of SEM, 4th ed. Guilford\n")
    cat("      CFI/TLI/NFI ≥.90 acceptable, ≥.95 good\n")
    cat("══════════════════════════════════════════════════════════════════\n")
  }

  # ── Results outputs ──────────────────────────────────────────────────────
  output$fit_output <- renderPrint({
    req(rv$model)
    .print_fit_guide(rv$model %>% fit)
  })

  output$param_table <- DT::renderDataTable({
    req(rv$model)
    params <- tryCatch(
      as.data.frame(rv$model %>% parameters),
      error = function(e) data.frame(Error = e$message)
    )
    DT::datatable(params, options = list(scrollX = TRUE, pageLength = 25))
  })

  output$mi_table <- DT::renderDataTable({
    req(rv$model)
    mi <- tryCatch(
      as.data.frame(rv$model %>% MIs),
      error = function(e) data.frame(Info = "Modification indices not available.")
    )
    DT::datatable(mi, options = list(scrollX = TRUE, pageLength = 25,
                                     order = list(list(ncol(mi) - 1, "desc"))))
  })

  output$matrix_output <- renderPrint({
    req(rv$model)
    mat_name <- input$matrix_select

    if (mat_name %in% c("ri_raw", "ri_norm", "sep_ri")) {
      if (mat_name == "sep_ri") { cat("Select 'RI Latent Network (raw)' or '(normalized)'.\n"); return() }
      sigma_z <- tryCatch(
        collapse_mg(getmatrix(rv$model, "sigma_zeta"), 1L),
        error = function(e) NULL
      )
      if (is.null(sigma_z)) {
        cat("sigma_zeta not available. Fit an LNM / RNM / LRNM model first.\n")
        return()
      }
      R_lat <- tryCatch(cov2cor(sigma_z), error = function(e) NULL)
      if (is.null(R_lat)) { cat("Could not compute latent correlation matrix.\n"); return() }
      lat_names <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(sigma_z))
                     rv$latent_names else colnames(sigma_z)
      if (is.null(lat_names)) lat_names <- paste0("L", seq_len(ncol(sigma_z)))
      dimnames(R_lat) <- list(lat_names, lat_names)
      RI <- johnson_rw_from_cor(R_lat)
      if (mat_name == "ri_norm") {
        col_sums <- colSums(RI)
        col_sums[col_sums < 1e-10] <- 1
        RI <- sweep(RI, 2, col_sums, "/")
        cat("Relative Importance Latent Network — NORMALIZED (each column sums to R² of that outcome)\n\n")
      } else {
        cat("Relative Importance Latent Network — RAW (each column sums to R² of that outcome)\n\n")
      }
      print(round(RI, 4))
      cat("\nColumn R² (total explained variance per outcome factor):\n")
      print(round(colSums(johnson_rw_from_cor(R_lat)), 4))
      return()
    }

    tryCatch({
      getmatrix(rv$model, mat_name, threshold = input$threshold_matrix)
    }, error = function(e) {
      cat("Matrix '", mat_name, "' not available for this model.\n")
      cat("Error:", e$message, "\n")
    })
  })

  output$model_summary <- renderPrint({
    req(rv$model)
    print(rv$model)
    df <- rv$analysis_data
    if (!is.null(df)) {
      n_actual   <- sum(complete.cases(df))
      n_moments  <- tryCatch({
        ft <- fit(rv$model)
        as.numeric(ft[ft$Measure == "nobs", "Value"])
      }, error = function(e) NA)
      cat("\n──────────────────────────────────────────────────\n")
      cat("NOTE: 'nobs' in the fit table above =", n_moments, "\n")
      cat("  This is the number of SAMPLE STATISTICS (means +\n")
      cat("  unique covariance elements), NOT the sample size.\n")
      cat("  Actual N used in estimation:", n_actual, "\n")
      cat("  Formula: p + p(p+1)/2 =", ncol(df), "+",
          ncol(df) * (ncol(df) + 1L) / 2L, "=",
          ncol(df) + ncol(df) * (ncol(df) + 1L) / 2L, "\n")
      cat("──────────────────────────────────────────────────\n")
    }
  })

  # ── Factor Scores & Loadings ─────────────────────────────────────────────
  output$fs_group_ui <- renderUI({
    req(rv$model)
    grp_var <- input$group_var
    if (is.null(grp_var) || grp_var == "") return(NULL)
    grp_levels <- get_grp_levels()
    if (length(grp_levels) == 0) return(NULL)
    named <- setNames(as.character(seq_len(length(grp_levels))),
                      paste0("Group ", seq_len(length(grp_levels)), ": ", grp_levels))
    selectInput("fs_group_sel", "Group:", choices = named, selected = "1")
  })

  loadings_data <- reactive({
    req(rv$model)
    grp <- if (!is.null(input$fs_group_sel)) as.integer(input$fs_group_sel) else 1L
    lam <- tryCatch(collapse_mg(getmatrix(rv$model, "lambda"), grp), error = function(e) NULL)
    req(lam)
    var_names <- rownames(lam)
    if (is.null(var_names) || length(var_names) == 0)
      var_names <- paste0("V", seq_len(nrow(lam)))
    data.frame(Variable = var_names, as.data.frame(round(lam, 4)),
               stringsAsFactors = FALSE, check.names = FALSE)
  })

  fscores_data <- reactive({
    req(rv$model, rv$analysis_data)
    grp     <- if (!is.null(input$fs_group_sel)) as.integer(input$fs_group_sel) else 1L
    mod     <- rv$model
    df      <- rv$analysis_data

    lam     <- tryCatch(collapse_mg(getmatrix(mod, "lambda"),     grp), error = function(e) NULL)
    sigma_z <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp), error = function(e) NULL)
    req(lam)

    vars_    <- rownames(lam)
    if (is.null(vars_) || length(vars_) == 0)
      vars_ <- head(input$selected_vars, nrow(lam))
    data_mat <- as.matrix(df[, intersect(vars_, names(df)), drop = FALSE])
    complete <- complete.cases(data_mat)

    grp_var <- input$group_var
    if (!is.null(grp_var) && grp_var != "" && !is.null(rv$raw_data[[grp_var]])) {
      grp_levels <- get_grp_levels()
      grp_label  <- grp_levels[grp]
      in_grp     <- as.character(rv$raw_data[[grp_var]]) == grp_label
      data_mat   <- data_mat[in_grp & complete, , drop = FALSE]
    } else {
      data_mat   <- data_mat[complete, , drop = FALSE]
    }

    req(nrow(data_mat) > 1)
    mu_hat    <- colMeans(data_mat)
    X_c       <- sweep(data_mat, 2, mu_hat, "-")
    # Use the sample covariance of the group data: always [n_obs × n_obs],
    # avoids dimension mismatches from model-internal sigma representations.
    sigma     <- cov(data_mat)
    sigma_inv <- tryCatch(solve(sigma), error = function(e) MASS::ginv(sigma))
    if (is.null(sigma_z)) sigma_z <- diag(ncol(lam))

    scores    <- X_c %*% sigma_inv %*% lam %*% sigma_z
    lat_names <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(lam))
                   rv$latent_names else colnames(lam)
    colnames(scores) <- lat_names
    as.data.frame(round(scores, 4))
  })

  output$loadings_table <- DT::renderDataTable({
    req(loadings_data())
    DT::datatable(loadings_data(), rownames = FALSE,
                  options = list(dom = "t", pageLength = 100, scrollX = TRUE))
  })

  output$factor_scores_table <- DT::renderDataTable({
    req(fscores_data())
    DT::datatable(fscores_data(), rownames = TRUE,
                  options = list(scrollX = TRUE, pageLength = 15))
  })

  output$download_loadings <- downloadHandler(
    filename = function() paste0("loadings_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(loadings_data(), file, row.names = FALSE)
  )

  output$download_fscores <- downloadHandler(
    filename = function() paste0("factor_scores_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(fscores_data(), file, row.names = TRUE)
  )

  # Reactive: full raw dataset with factor scores appended as prefixed columns.
  # Scores are computed per group (for multi-group models) and written back into
  # the corresponding rows; rows with missing data receive NA for the score cols.
  fscores_merged_data <- reactive({
    req(rv$model, rv$raw_data)
    mod    <- rv$model
    df_raw <- rv$raw_data

    prefix <- switch(input$model_family,
      "lvm"       = "CFA",
      "lnm"       = "LNM",
      "rnm"       = "RNM",
      "lrnm"      = "LRNM",
      "dlvm1"     = "DLVM",
      "ri_clpm"   = "RICLPM",
      "panelgvar" = "PGVAR",
      "FS"
    )

    grp_var    <- input$group_var
    grp_levels <- get_grp_levels()
    n_grp      <- max(1L, length(grp_levels))
    out        <- df_raw

    for (grp in seq_len(n_grp)) {
      lam     <- tryCatch(collapse_mg(getmatrix(mod, "lambda"),     grp), error = function(e) NULL)
      if (is.null(lam)) next
      sigma_z <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp), error = function(e) NULL)

      vars_ <- rownames(lam)
      if (is.null(vars_) || length(vars_) == 0)
        vars_ <- head(input$selected_vars, nrow(lam))

      data_mat <- as.matrix(df_raw[, intersect(vars_, names(df_raw)), drop = FALSE])
      complete <- complete.cases(data_mat)

      in_grp <- if (!is.null(grp_var) && grp_var != "" && length(grp_levels) >= grp) {
        idx <- as.character(df_raw[[grp_var]]) == grp_levels[grp]
        idx[is.na(idx)] <- FALSE
        idx
      } else {
        rep(TRUE, nrow(df_raw))
      }

      row_mask <- in_grp & complete
      if (sum(row_mask) < 2) next

      data_sub  <- data_mat[row_mask, , drop = FALSE]
      mu_hat    <- colMeans(data_sub)
      X_c       <- sweep(data_sub, 2, mu_hat, "-")
      sigma     <- cov(data_sub)
      sigma_inv <- tryCatch(solve(sigma), error = function(e) MASS::ginv(sigma))
      if (is.null(sigma_z)) sigma_z <- diag(ncol(lam))

      scores    <- X_c %*% sigma_inv %*% lam %*% sigma_z
      lat_names <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(lam))
                     rv$latent_names else colnames(lam)
      if (is.null(lat_names) || length(lat_names) == 0)
        lat_names <- paste0("F", seq_len(ncol(lam)))
      colnames(scores) <- paste0(prefix, "_", lat_names)

      for (cn in colnames(scores)) {
        if (!(cn %in% names(out))) out[[cn]] <- NA_real_
      }
      out[row_mask, colnames(scores)] <- round(scores, 4)
    }
    out
  })

  output$download_fscores_merged <- downloadHandler(
    filename = function() {
      paste0("data_", input$model_family, "_fscores_", Sys.Date(), ".csv")
    },
    content = function(file) write.csv(fscores_merged_data(), file, row.names = FALSE)
  )

  # ── Plots ────────────────────────────────────────────────────────────────
  # Run a psychonetrics model with the chosen optimizer via setoptimizer().
  do_runmodel <- function(m) {
    m %>% setoptimizer(input$optimizer) %>% runmodel()
  }

  # ── Node predictability helpers ──────────────────────────────────────────────
  # OLS R² for observed-variable networks (GGM, Ising, Residual, EGA)
  node_predictability_r2 <- function(data_mat) {
    data_mat <- na.omit(as.matrix(data_mat))
    if (nrow(data_mat) < 3 || ncol(data_mat) < 2) return(NULL)
    p <- ncol(data_mat); r2 <- numeric(p); names(r2) <- colnames(data_mat)
    for (i in seq_len(p)) {
      y   <- data_mat[, i]
      fit <- tryCatch(lm.fit(cbind(1, data_mat[, -i, drop = FALSE]), y),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        rss  <- sum(fit$residuals^2, na.rm = TRUE)
        tss  <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
        r2[i] <- if (tss > 1e-10) max(0, min(1, 1 - rss / tss)) else 0
      }
    }
    r2
  }

  # Analytical latent R² from sigma_zeta (LNM, RNM, LRNM, CFA latent network)
  # Uses the latent correlation matrix for scale-invariant, numerically stable R²:
  # R²_i = 1 − 1/K_ii  where K = solve(cov2cor(sigma_zeta))
  latent_predictability_r2 <- function(mod_arg, grp_arg) {
    sg <- tryCatch(collapse_mg(getmatrix(mod_arg, "sigma_zeta"), grp_arg),
                   error = function(e) NULL)
    if (is.null(sg) || nrow(sg) < 2) return(NULL)
    # Convert to correlation matrix (scale-invariant; diag = 1 simplifies formula)
    rg <- tryCatch(cov2cor(sg), error = function(e) NULL)
    if (is.null(rg)) return(NULL)
    kg <- tryCatch(solve(rg),
                   error = function(e) tryCatch(MASS::ginv(rg), error = function(e2) NULL))
    if (is.null(kg)) return(NULL)
    # diag(rg) = 1, so R²_i = 1 - 1/K_ii
    r2 <- pmax(0, pmin(1, 1 - 1 / diag(kg)))
    r2
  }

  # Resolve label-scale inputs into qgraph arguments
  qgraph_label_args <- function() {
    sc <- input$qgraph_label_scale
    list(
      label.scale       = (sc != "none"),
      label.scale.equal = (sc == "equal"),
      label.cex         = input$qgraph_label_cex
    )
  }

  # Resolve colour + edge-label size inputs into qgraph arguments
  qgraph_color_args <- function() {
    rep_val <- if (!is.null(input$qgraph_repulsion)) input$qgraph_repulsion else 1
    list(
      edge.label.cex = input$edge_label_cex,
      repulsion      = rep_val
    )
  }

  # Build pie/pieColor/pieBorder args when predictability rings are enabled.
  # Accepts either:
  #   - plain numeric vector  → rings on all nodes (GGM / EGA / latent_network)
  #   - list with NULL slots  → rings only on non-NULL nodes (factor plots: latents only)
  pie_args <- function(r2_vec) {
    if (is.null(r2_vec) || !isTRUE(input$show_predictability)) return(list())
    if (is.list(r2_vec)) {
      # Mixed list: NULL elements mean "no ring for this node"
      has_any <- any(sapply(r2_vec, function(x) !is.null(x) && length(x) > 0 && x > 0))
      if (!has_any) return(list())
      return(list(
        pie       = r2_vec,
        pieColor  = input$pred_ring_color,
        pieBorder = input$pred_ring_border
      ))
    }
    # Plain numeric vector path
    r2_vec <- unname(r2_vec)
    r2_vec[is.na(r2_vec)] <- 0
    r2_vec <- pmax(0, pmin(1, r2_vec))
    if (all(r2_vec == 0)) return(list())
    list(
      pie       = r2_vec,
      pieColor  = input$pred_ring_color,
      pieBorder = input$pred_ring_border
    )
  }

  # Collapse a multi-group matrix (3D array OR list of matrices) to a single
  # 2D matrix.  getmatrix() in psychonetrics can return either form depending
  # on the model version, so we handle both here.
  collapse_mg <- function(m, grp = 1L) {
    if (is.list(m) && !is.data.frame(m)) return(as.matrix(m[[grp]]))
    if (is.array(m) && length(dim(m)) == 3) return(as.matrix(m[, , grp]))
    as.matrix(m)
  }

  # Returns the authoritative sorted group labels from raw data, trimmed to the
  # number of groups the fitted model actually has (avoids phantom groups from
  # stray values / NA-coded entries in the grouping column).
  get_grp_levels <- function() {
    grp_var <- input$group_var
    if (is.null(grp_var) || grp_var == "") return(character(0))
    req(rv$raw_data)
    raw_levels <- as.character(sort(na.omit(unique(rv$raw_data[[grp_var]]))))
    if (!is.null(rv$model)) {
      lam_raw <- tryCatch(getmatrix(rv$model, "lambda"), error = function(e) NULL)
      if (!is.null(lam_raw)) {
        k <- if (is.array(lam_raw) && length(dim(lam_raw)) == 3) dim(lam_raw)[3]
             else if (is.list(lam_raw) && !is.data.frame(lam_raw)) length(lam_raw)
             else length(raw_levels)
        raw_levels <- raw_levels[seq_len(min(k, length(raw_levels)))]
      }
    }
    raw_levels
  }

  # ── EGA community-detection cache ──────────────────────────────────────────
  # Community detection is expensive and only needs to re-run when the fitted
  # model or algorithm changes — NOT on visual changes like the edge cut-off.
  # We use a plain R environment (not reactiveVal) so writing to the cache
  # inside renderPlot does not trigger another render cycle.
  .ega_cache <- new.env(parent = emptyenv())

  observeEvent(list(rv$model, input$ega_algorithm), {
    rm(list = ls(.ega_cache), envir = .ega_cache)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # When spring layout is used with multiple groups, derive node positions from
  # group 1 and reuse them for all groups so the visual comparison is fair.
  mg_spring_layout <- reactive({
    if (is.null(input$qgraph_layout) || input$qgraph_layout != "spring") return(NULL)
    grp_var <- input$group_var
    if (is.null(grp_var) || grp_var == "") return(NULL)
    mod <- rv$model
    if (is.null(mod)) return(NULL)
    mat_name <- switch(input$plot_type,
      "latent_network"   = "omega_zeta",
      "residual_network" = "omega_epsilon",
      "ggm_network"      = "omega",
      "ega_ggm"          = "omega",
      "ega_ising"        = "omega",
      "beta_network"     = "beta",
      "within_network"   = "omega_zeta_within",
      "between_network"  = "omega_zeta_between",
      NULL
    )
    if (is.null(mat_name)) return(NULL)
    net1 <- tryCatch({
      m <- getmatrix(mod, mat_name)
      if (is.null(m)) stop("null")
      collapse_mg(m, 1L)
    }, error = function(e) NULL)
    if (is.null(net1)) return(NULL)
    rep_val <- if (!is.null(input$qgraph_repulsion)) input$qgraph_repulsion else 1
    g <- tryCatch(
      qgraph(net1, layout = "spring", repulsion = rep_val, DoNotPlot = TRUE),
      error = function(e) NULL
    )
    if (is.null(g)) return(NULL)
    g$layout
  })

  # LMG (Lindeman-Merenda-Gold) relative importance from a latent correlation matrix.
  # R_lat: p×p correlation matrix derived from model-implied sigma_zeta via cov2cor().
  # Returns RI[j,i] = relative importance of predictor j for outcome i (directed: j → i).
  # RI[j,i] >= 0; column i sums to R²_i for that outcome factor.
  # LMG averages marginal R² contributions over all predictor orderings —
  # unlike Johnson's eigendecomposition approach, LMG correctly assigns unequal weights
  # to the two predictors in a 2-predictor case (the common case for 3-factor models).
  johnson_rw_from_cor <- function(R_lat) {
    p  <- ncol(R_lat)
    nm <- colnames(R_lat)
    if (is.null(nm)) nm <- paste0("L", seq_len(p))
    RI <- matrix(0, p, p, dimnames = list(nm, nm))

    # R² from regressing outcome (last position) on a subset of predictors
    r2_sub <- function(sub, rxy, Rxx) {
      if (length(sub) == 0L) return(0)
      rys <- rxy[sub]
      Rxs <- Rxx[sub, sub, drop = FALSE]
      tryCatch(max(0, drop(t(rys) %*% solve(Rxs, rys))),
               error = function(e) max(0, drop(t(rys) %*% (MASS::ginv(Rxs) %*% rys))))
    }

    for (i in seq_len(p)) {
      pred <- seq_len(p)[-i]
      q    <- length(pred)
      if (q < 1L) next
      rxy  <- R_lat[pred, i]
      Rxx  <- R_lat[pred, pred, drop = FALSE]
      if (q == 1L) { RI[pred, i] <- rxy^2; next }

      contrib <- numeric(q)
      for (j in seq_len(q)) {
        others <- seq_len(q)[-j]
        n_oth  <- length(others)
        # Average marginal contribution of j over all 2^(q-1) subsets of others
        for (mask in seq_len(2L^n_oth) - 1L) {
          S  <- if (n_oth == 0L) integer(0) else
                  others[as.logical(intToBits(mask)[seq_len(n_oth)])]
          s  <- length(S)
          wt <- factorial(s) * factorial(q - 1L - s) / factorial(q)
          contrib[j] <- contrib[j] +
            wt * (r2_sub(c(S, j), rxy, Rxx) - r2_sub(S, rxy, Rxx))
        }
      }
      RI[pred, i] <- contrib
    }
    RI
  }

  make_plot <- function(grp = 1L) {
    req(rv$model)
    mod <- rv$model

    # ── Shared helpers ─────────────────────────────────────────────────────────
    darken_hex <- function(hex, amt) {
      v <- col2rgb(hex) / 255
      rgb(v[1] * amt, v[2] * amt, v[3] * amt)
    }

    # Compute mean (uncentered) factor scores per latent factor for the current
    # group. Returns a named numeric vector or NULL on failure.
    mean_fscores_for_plot <- function(lam_mat) {
      tryCatch({
        sigma_z  <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp),
                             error = function(e) NULL)
        if (is.null(sigma_z)) sigma_z <- diag(ncol(lam_mat))
        vars_    <- rownames(lam_mat)
        if (is.null(vars_) || length(vars_) == 0)
          vars_ <- head(input$selected_vars, nrow(lam_mat))
        df_all   <- rv$analysis_data
        data_mat <- as.matrix(df_all[, intersect(vars_, names(df_all)), drop = FALSE])
        complete <- complete.cases(data_mat)
        grp_var  <- input$group_var
        if (!is.null(grp_var) && grp_var != "") {
          gl       <- get_grp_levels()
          in_grp   <- as.character(rv$raw_data[[grp_var]]) == gl[grp]
          data_mat <- data_mat[in_grp & complete, , drop = FALSE]
        } else {
          data_mat <- data_mat[complete, , drop = FALSE]
        }
        if (nrow(data_mat) < 2) stop("too few rows")
        S_inv <- tryCatch(solve(cov(data_mat)), error = function(e) MASS::ginv(cov(data_mat)))
        colMeans(data_mat %*% S_inv %*% lam_mat %*% sigma_z)  # uncentered = severity level
      }, error = function(e) NULL)
    }

    norm01 <- function(x, mn = input$darken_min) {
      if (diff(range(x)) < 1e-6) return(rep(0.75, length(x)))
      mn + (1 - mn) * (x - min(x)) / diff(range(x))
    }

    # Subset a data frame to only the rows belonging to the current group.
    # Returns the full df unchanged when no group variable is active or when
    # the row count doesn't align (safety fallback).
    group_data <- function(df) {
      grp_var <- input$group_var
      if (is.null(grp_var) || grp_var == "" || is.null(rv$raw_data)) return(df)
      gl <- get_grp_levels()
      if (grp > length(gl)) return(df)
      in_grp <- as.character(rv$raw_data[[grp_var]]) == as.character(gl[grp])
      in_grp[is.na(in_grp)] <- FALSE
      if (length(in_grp) != nrow(df)) return(df)
      df[in_grp, , drop = FALSE]
    }

    # For multi-group spring plots: return the group-1 layout so all panels
    # share identical node positions. Falls back to input$qgraph_layout otherwise.
    effective_layout <- function() {
      if (!is.null(input$qgraph_layout) && input$qgraph_layout == "spring" &&
          !is.null(input$group_var) && input$group_var != "") {
        sl <- mg_spring_layout()
        if (!is.null(sl)) return(sl)
      }
      input$qgraph_layout
    }

    plot_fn <- function() {
      type <- input$plot_type

      if (type %in% c("latent_network", "residual_network", "ggm_network")) {
        mat_name <- switch(type,
          "latent_network"   = "omega_zeta",
          "residual_network" = "omega_epsilon",
          "ggm_network"      = "omega"
        )

        net <- tryCatch(
          getmatrix(mod, mat_name),
          error = function(e) NULL
        )
        if (is.null(net)) {
          plot.new()
          text(0.5, 0.5, paste("Matrix", mat_name, "not available"), cex = 1.5)
          return()
        }
        # getmatrix() can return a 3D array OR list for multi-group models
        net <- collapse_mg(net, grp)

        # Apply edge cutoff: zero out edges below threshold (qgraph cut= only fades)
        cutoff <- input$qgraph_cutoff
        if (!is.null(cutoff) && cutoff > 0) net[abs(net) < cutoff] <- 0

        # Determine labels
        labs <- if (type == "latent_network" && !is.null(rv$latent_names) &&
                    length(rv$latent_names) == ncol(net)) {
          rv$latent_names
        } else if (!is.null(colnames(net)) && !all(colnames(net) == "")) {
          colnames(net)
        } else if (type != "latent_network" &&
                   length(input$selected_vars) == ncol(net)) {
          input$selected_vars
        } else {
          paste0("V", seq_len(ncol(net)))
        }

        # For latent_network: size + color encode mean factor score (severity)
        if (type == "latent_network") {
          lam_for_fs <- tryCatch(collapse_mg(getmatrix(mod, "lambda"), grp),
                                 error = function(e) NULL)
          mean_fs    <- if (!is.null(lam_for_fs)) mean_fscores_for_plot(lam_for_fs) else NULL
          importance <- if (!is.null(mean_fs) && !anyNA(mean_fs)) abs(mean_fs)
                        else colSums(abs(net))
          imp_range  <- range(importance)
          default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                           "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
          pal <- sapply(seq_len(ncol(net)), function(i) {
            col <- input[[paste0("lat_color_", i)]]
            if (!is.null(col) && nchar(col) > 0) col
            else if (i <= length(default_pal)) default_pal[i] else "#999999"
          })
          if (isTRUE(input$disable_node_scaling)) {
            node_sizes <- rep(input$qgraph_vsize, ncol(net))
            node_cols  <- pal
          } else {
            node_sizes <- if (diff(imp_range) < 1e-6) rep(input$qgraph_vsize, ncol(net))
                          else input$qgraph_vsize * 1.2 +
                               (importance - imp_range[1]) / diff(imp_range) *
                               input$qgraph_vsize * 1.4
            lat_dark   <- norm01(importance)
            node_cols  <- mapply(darken_hex, pal, lat_dark)
          }
        } else {
          # GGM / Ising / residual: scale nodes by |colMeans| of analysis data
          df_obs   <- group_data(tryCatch(analysis_data(), error = function(e) rv$analysis_data))
          obs_vars <- colnames(net)
          if (is.null(obs_vars) || all(nchar(obs_vars) == 0))
            obs_vars <- if (length(input$selected_vars) == ncol(net)) input$selected_vars else NULL
          means_vec <- if (!is.null(df_obs) && !is.null(obs_vars)) {
            keep <- intersect(obs_vars, names(df_obs))
            if (length(keep) == ncol(net))
              abs(colMeans(df_obs[, keep, drop = FALSE], na.rm = TRUE))
            else NULL
          } else NULL
          if (!is.null(means_vec) && !isTRUE(input$disable_node_scaling)) {
            imp_range  <- range(means_vec, na.rm = TRUE)
            node_sizes <- if (diff(imp_range) < 1e-6) rep(input$qgraph_vsize, ncol(net))
                          else input$qgraph_vsize * 0.6 +
                               (means_vec - imp_range[1]) / diff(imp_range) * input$qgraph_vsize * 1.2
            dark_amt   <- norm01(means_vec)
            node_cols  <- mapply(darken_hex, rep(input$node_color, ncol(net)), dark_amt)
          } else {
            node_sizes <- rep(input$qgraph_vsize, ncol(net))
            node_cols  <- rep(input$node_color,   ncol(net))
          }
        }

        # Node predictability rings
        r2_vec <- if (isTRUE(input$show_predictability)) {
          if (type == "latent_network") {
            latent_predictability_r2(mod, grp)
          } else {
            # psychonetrics matrices have empty dimnames — fall back to selected_vars
            vars_in_net <- colnames(net)
            if (is.null(vars_in_net) || all(nchar(vars_in_net) == 0))
              vars_in_net <- if (length(input$selected_vars) == ncol(net)) input$selected_vars else NULL
            df_pred <- group_data(tryCatch(analysis_data(), error = function(e) rv$analysis_data))
            if (!is.null(df_pred) && !is.null(vars_in_net)) {
              keep <- intersect(vars_in_net, names(df_pred))
              if (length(keep) == length(vars_in_net))
                node_predictability_r2(df_pred[, keep, drop = FALSE])
              else NULL
            } else NULL
          }
        } else NULL

        do.call(qgraph, c(list(
               input       = net,
               labels      = labs,
               layout      = effective_layout(),
               theme       = input$qgraph_theme,
               color       = node_cols,
               vsize       = node_sizes,
               esize        = input$qgraph_esize,
               edge.labels  = input$qgraph_labels_est,
               border.color = input$node_border_color,
               title        = gsub("_", " ", mat_name),
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args(),
               pie_args(r2_vec)))

      } else if (type == "ri_latent_network") {
        sigma_z <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp), error = function(e) NULL)
        if (is.null(sigma_z)) {
          plot.new()
          text(0.5, 0.5, "sigma_zeta not available.\nFit an LNM model first.", cex = 1.3)
          return()
        }
        lat_names <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(sigma_z))
                       rv$latent_names else colnames(sigma_z)
        if (is.null(lat_names)) lat_names <- paste0("L", seq_len(ncol(sigma_z)))
        dimnames(sigma_z) <- list(lat_names, lat_names)

        R_lat <- tryCatch(cov2cor(sigma_z), error = function(e) NULL)
        if (is.null(R_lat)) {
          plot.new()
          text(0.5, 0.5, "Could not compute latent correlation matrix.", cex = 1.3)
          return()
        }

        RI <- johnson_rw_from_cor(R_lat)
        if (isTRUE(input$ri_normalize)) {
          col_sums <- colSums(RI)
          col_sums[col_sums < 1e-10] <- 1
          RI <- sweep(RI, 2, col_sums, "/")
        }

        cutoff <- input$qgraph_cutoff
        if (!is.null(cutoff) && cutoff > 0) RI[RI < cutoff] <- 0

        # Node scaling: same procedure as factor_network
        lam_for_fs <- tryCatch(collapse_mg(getmatrix(mod, "lambda"), grp), error = function(e) NULL)
        mean_fs    <- if (!is.null(lam_for_fs)) mean_fscores_for_plot(lam_for_fs) else NULL
        importance <- if (!is.null(mean_fs) && !anyNA(mean_fs)) abs(mean_fs)
                      else colSums(abs(RI))
        imp_range  <- range(importance)
        default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                         "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
        pal <- sapply(seq_len(ncol(RI)), function(i) {
          col <- input[[paste0("lat_color_", i)]]
          if (!is.null(col) && nchar(col) > 0) col
          else if (i <= length(default_pal)) default_pal[i] else "#999999"
        })
        if (isTRUE(input$disable_node_scaling)) {
          node_sizes <- rep(input$qgraph_vsize, ncol(RI))
          node_cols  <- pal
        } else {
          node_sizes <- if (diff(imp_range) < 1e-6) rep(input$qgraph_vsize, ncol(RI))
                        else input$qgraph_vsize * 1.2 +
                             (importance - imp_range[1]) / diff(imp_range) *
                             input$qgraph_vsize * 1.4
          lat_dark   <- norm01(importance)
          node_cols  <- mapply(darken_hex, pal, lat_dark)
        }

        do.call(qgraph, c(list(
               input        = RI,
               directed     = TRUE,
               labels       = lat_names,
               layout       = effective_layout(),
               theme        = input$qgraph_theme,
               color        = node_cols,
               vsize        = node_sizes,
               esize        = input$qgraph_esize,
               edge.labels  = input$qgraph_labels_est,
               border.color = input$node_border_color,
               title        = "Relative Importance Latent Network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args()))

      } else if (type == "factor_ri_network") {
        lambda      <- tryCatch(getmatrix(mod, "lambda"),     error = function(e) NULL)
        sigma_z_raw <- tryCatch(getmatrix(mod, "sigma_zeta"), error = function(e) NULL)

        if (is.null(lambda)) {
          plot.new()
          text(0.5, 0.5, "Lambda matrix not available.\nFit a factor (LNM/RNM) model first.", cex = 1.3)
          return()
        }
        lambda  <- collapse_mg(lambda, grp)
        sigma_z <- if (!is.null(sigma_z_raw)) collapse_mg(sigma_z_raw, grp) else diag(ncol(lambda))

        lat_names <- colnames(lambda)   # correct order comes from the fitted model
        if (is.null(lat_names) || length(lat_names) == 0 || all(nchar(lat_names) == 0) ||
            all(grepl("^eta_?\\d+$", lat_names, ignore.case = TRUE))) {
          lat_names <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(lambda))
                         rv$latent_names
                       else paste0("L", seq_len(ncol(lambda)))
        }
        dimnames(sigma_z) <- list(lat_names, lat_names)

        R_lat <- tryCatch(cov2cor(sigma_z), error = function(e) NULL)
        if (is.null(R_lat)) {
          plot.new()
          text(0.5, 0.5, "Could not compute latent correlation matrix.", cex = 1.3)
          return()
        }
        RI <- johnson_rw_from_cor(R_lat)
        if (isTRUE(input$ri_normalize)) {
          col_sums <- colSums(RI)
          col_sums[col_sums < 1e-10] <- 1
          RI <- sweep(RI, 2, col_sums, "/")
        }

        cutoff <- input$qgraph_cutoff
        if (!is.null(cutoff) && cutoff > 0) {
          RI[RI < cutoff]            <- 0
          lambda[abs(lambda) < cutoff] <- 0
        }

        n_obs   <- nrow(lambda)
        n_lat   <- ncol(lambda)
        n_total <- n_lat + n_obs

        W <- matrix(0, n_total, n_total)
        W[1:n_lat, 1:n_lat]                   <- RI
        W[1:n_lat, (n_lat + 1):n_total]       <- t(lambda)

        # All edges directed: RI arrows between latents, loading arrows to observed
        dir_mat <- matrix(FALSE, n_total, n_total)
        dir_mat[1:n_lat, 1:n_lat]             <- TRUE
        dir_mat[1:n_lat, (n_lat + 1):n_total] <- TRUE

        obs_labs <- rownames(lambda)   # correct order comes from the fitted model
        if (is.null(obs_labs) || length(obs_labs) == 0 || all(nchar(obs_labs) == 0)) {
          obs_labs <- if (length(input$selected_vars) == n_obs) input$selected_vars
                      else paste0("V", seq_len(n_obs))
        }
        shapes <- c(rep("circle", n_lat), rep("square", n_obs))

        # Node scaling: same as factor_network
        mean_fs    <- mean_fscores_for_plot(lambda)
        importance <- if (!is.null(mean_fs)) abs(mean_fs) else colMeans(abs(lambda))
        imp_range  <- range(importance)
        default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                         "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
        pal <- sapply(seq_len(n_lat), function(i) {
          col <- input[[paste0("lat_color_", i)]]
          if (!is.null(col) && nchar(col) > 0) col
          else if (i <= length(default_pal)) default_pal[i] else "#999999"
        })
        dominant <- apply(abs(lambda), 1, which.max)
        if (isTRUE(input$disable_node_scaling)) {
          lat_sizes <- rep(input$qgraph_vsize * 1.8, n_lat)
          lat_cols  <- pal[seq_len(n_lat)]
          obs_cols  <- pal[dominant]
        } else {
          lat_sizes <- if (diff(imp_range) < 1e-6) rep(input$qgraph_vsize * 1.8, n_lat)
                       else input$qgraph_vsize * 1.2 +
                            (importance - imp_range[1]) / diff(imp_range) * input$qgraph_vsize * 1.4
          lat_dark  <- norm01(importance)
          lat_cols  <- mapply(darken_hex, pal[seq_len(n_lat)], lat_dark)
          load_str  <- apply(abs(lambda), 1, max)
          obs_dark  <- norm01(load_str)
          obs_cols  <- mapply(darken_hex, pal[dominant], obs_dark)
        }
        vsizes    <- c(lat_sizes, rep(input$qgraph_vsize, n_obs))
        node_cols <- c(lat_cols, obs_cols)

        lay <- if (input$qgraph_layout == "circle") {
          angles_lat <- seq(0, 2 * pi, length.out = n_lat + 1)[seq_len(n_lat)]
          lat_pos    <- cbind(cos(angles_lat), sin(angles_lat)) * 0.38
          obs_pos    <- matrix(0, n_obs, 2)
          for (k in seq_len(n_lat)) {
            idx <- which(dominant == k)
            n_k <- length(idx)
            if (n_k == 0) next
            half <- (pi / n_lat) * 0.85
            angs <- if (n_k == 1) angles_lat[k] else
                      seq(angles_lat[k] - half, angles_lat[k] + half, length.out = n_k)
            obs_pos[idx, 1] <- cos(angs)
            obs_pos[idx, 2] <- sin(angs)
          }
          rbind(lat_pos, obs_pos)
        } else {
          effective_layout()
        }

        # Node predictability rings: latents only — no rings on observed (square) nodes
        r2_vec <- if (isTRUE(input$show_predictability)) {
          r2_lat <- latent_predictability_r2(mod, grp)
          if (!is.null(r2_lat) && length(r2_lat) == n_lat)
            c(as.list(r2_lat), vector("list", n_obs))  # NULL slots = no ring on observed squares
          else NULL
        } else NULL

        # Boost near-zero loadings for display so all non-zero cross-loadings are visible.
        lambda_w <- lambda
        max_load <- if (any(lambda_w != 0)) max(abs(lambda_w[lambda_w != 0])) else 0
        if (max_load > 1e-6) {
          min_vis <- max_load * 0.20
          nz      <- lambda_w != 0 & abs(lambda_w) < min_vis
          lambda_w[nz] <- sign(lambda_w[nz]) * min_vis
        }
        W_disp <- W
        W_disp[1:n_lat, (n_lat + 1):n_total] <- t(lambda_w)

        # Edge label matrix: actual fitted values (not boosted), blank where zero
        edge_lab <- matrix("", n_total, n_total)
        edge_lab[1:n_lat, 1:n_lat] <-
          ifelse(RI != 0, as.character(round(RI, 2)), "")
        edge_lab[1:n_lat, (n_lat + 1):n_total] <-
          ifelse(t(lambda) != 0, as.character(round(t(lambda), 2)), "")

        do.call(qgraph, c(list(
               input        = W_disp,
               labels       = c(lat_names, obs_labs),
               shape        = shapes,
               vsize        = vsizes,
               color        = node_cols,
               directed     = dir_mat,
               layout       = lay,
               theme        = input$qgraph_theme,
               edge.labels  = if (isTRUE(input$qgraph_labels_est)) edge_lab else FALSE,
               esize        = input$qgraph_esize,
               border.color = input$node_border_color,
               title        = "Factor structure & latent RI network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args(),
               pie_args(r2_vec)))

      } else if (type == "beta_network") {
        # ── Temporal (cross-lagged) beta network ──────────────────────────
        beta_mat <- tryCatch(getmatrix(mod, "beta"), error = function(e) NULL)
        if (is.null(beta_mat)) {
          plot.new()
          text(0.5, 0.5, "Beta matrix not available.\nFit a panel model first.", cex = 1.3)
          return()
        }
        beta_mat <- collapse_mg(beta_mat, grp)
        lat_labs <- if (!is.null(rv$latent_names) && length(rv$latent_names) == nrow(beta_mat))
                      rv$latent_names else rownames(beta_mat)
        if (is.null(lat_labs)) lat_labs <- paste0("L", seq_len(nrow(beta_mat)))

        do.call(qgraph, c(list(
               input        = beta_mat,
               directed     = TRUE,
               labels       = lat_labs,
               layout       = effective_layout(),
               theme        = input$qgraph_theme,
               edge.labels  = input$qgraph_labels_est,
               esize        = input$qgraph_esize,
               vsize        = input$qgraph_vsize,
               border.color = input$node_border_color,
               title        = "Temporal network (cross-lagged beta)",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args()))

      } else if (type == "within_network") {
        # ── Within-person contemporaneous network ─────────────────────────
        mat_name <- if (input$model_family %in% c("panelgvar", "dlvm1")) "omega_zeta_within"
                    else "omega_zeta"
        net <- tryCatch(collapse_mg(getmatrix(mod, mat_name), grp), error = function(e) NULL)
        if (is.null(net)) {
          # fallback to plain sigma_zeta
          net <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp), error = function(e) NULL)
        }
        if (is.null(net)) {
          plot.new()
          text(0.5, 0.5, "Within-person matrix not available.", cex = 1.3)
          return()
        }
        lat_labs <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(net))
                      rv$latent_names else colnames(net)
        if (is.null(lat_labs)) lat_labs <- paste0("L", seq_len(ncol(net)))

        do.call(qgraph, c(list(
               input        = net,
               directed     = FALSE,
               labels       = lat_labs,
               layout       = effective_layout(),
               theme        = input$qgraph_theme,
               edge.labels  = input$qgraph_labels_est,
               esize        = input$qgraph_esize,
               vsize        = input$qgraph_vsize,
               border.color = input$node_border_color,
               title        = "Within-person contemporaneous network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args()))

      } else if (type == "between_network") {
        # ── Between-person network ────────────────────────────────────────
        mat_name <- if (input$model_family %in% c("panelgvar", "dlvm1")) "omega_zeta_between"
                    else "sigma_zeta"
        net <- tryCatch(collapse_mg(getmatrix(mod, mat_name), grp), error = function(e) NULL)
        if (is.null(net)) {
          # fallback to plain sigma_zeta / omega_zeta
          net <- tryCatch(collapse_mg(getmatrix(mod, "omega_zeta"), grp), error = function(e) NULL)
        }
        if (is.null(net)) {
          plot.new()
          text(0.5, 0.5, "Between-person matrix not available.", cex = 1.3)
          return()
        }
        lat_labs <- if (!is.null(rv$latent_names) && length(rv$latent_names) == ncol(net))
                      rv$latent_names else colnames(net)
        if (is.null(lat_labs)) lat_labs <- paste0("L", seq_len(ncol(net)))

        do.call(qgraph, c(list(
               input        = net,
               directed     = FALSE,
               labels       = lat_labs,
               layout       = effective_layout(),
               theme        = input$qgraph_theme,
               edge.labels  = input$qgraph_labels_est,
               esize        = input$qgraph_esize,
               vsize        = input$qgraph_vsize,
               border.color = input$node_border_color,
               title        = "Between-person network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args()))

      } else if (type %in% c("ega_ggm", "ega_ising")) {
        # ── Exploratory Graph Analysis ──────────────────────────────────────────
        net <- tryCatch(getmatrix(mod, "omega"), error = function(e) NULL)
        if (is.null(net)) {
          plot.new()
          text(0.5, 0.5, "Omega matrix not available.\nFit a GGM or Ising model first.",
               cex = 1.3)
          return()
        }
        net <- collapse_mg(net, grp)

        labs <- if (!is.null(colnames(net)) && !all(nchar(colnames(net)) == 0)) {
          colnames(net)
        } else if (length(input$selected_vars) == ncol(net)) {
          input$selected_vars
        } else {
          paste0("V", seq_len(ncol(net)))
        }

        # Community detection: use cache so it does not re-run on visual-only
        # changes (e.g. edge cut-off, theme, node size). Cache is keyed by group
        # and algorithm; it is cleared by observeEvent whenever rv$model or
        # input$ega_algorithm changes.
        algo      <- isolate(input$ega_algorithm) %||% "walktrap"
        cache_key <- paste0("g", grp, "_", algo)

        if (exists(cache_key, envir = .ega_cache)) {
          comm_result <- get(cache_key, envir = .ega_cache)
        } else {
          # Build igraph from the FULL (unthresholded) network — the cut-off is
          # a display parameter and should not influence community membership.
          abs_net_full <- abs(net)
          diag(abs_net_full) <- 0
          ig <- igraph::graph_from_adjacency_matrix(
            abs_net_full, mode = "undirected", weighted = TRUE, diag = FALSE
          )

          communities <- tryCatch({
            switch(algo,
              walktrap    = igraph::cluster_walktrap(ig),
              louvain     = igraph::cluster_louvain(ig, weights = igraph::E(ig)$weight),
              fast_greedy = igraph::cluster_fast_greedy(ig),
              spinglass   = igraph::cluster_spinglass(ig, weights = igraph::E(ig)$weight),
              igraph::cluster_walktrap(ig)
            )
          }, error = function(e) {
            showNotification(paste("Community detection failed:", e$message),
                             type = "warning", duration = 8)
            NULL
          })

          if (is.null(communities)) {
            plot.new()
            text(0.5, 0.5, "Community detection failed.\nTry a different algorithm.",
                 cex = 1.2)
            return()
          }

          membership <- igraph::membership(communities)
          k          <- max(membership)
          grp_list   <- setNames(
            lapply(seq_len(k), function(d) which(membership == d)),
            paste0("Dim ", seq_len(k))
          )

          comm_result <- list(membership = membership, k = k, grp_list = grp_list)
          assign(cache_key, comm_result, envir = .ega_cache)
        }

        membership <- comm_result$membership
        k          <- comm_result$k
        grp_list   <- comm_result$grp_list

        comm_palette <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                          "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
        grp_colors <- comm_palette[seq_len(k)]

        # Apply edge cut-off for display only (after community detection)
        net_display <- net
        cutoff <- input$qgraph_cutoff
        if (!is.null(cutoff) && cutoff > 0) net_display[abs(net_display) < cutoff] <- 0

        # Node predictability rings (OLS on analysis data)
        r2_vec <- if (isTRUE(input$show_predictability)) {
          vars_in_net <- colnames(net)
          if (is.null(vars_in_net) || all(nchar(vars_in_net) == 0))
            vars_in_net <- if (length(input$selected_vars) == ncol(net)) input$selected_vars else NULL
          df_pred <- group_data(tryCatch(analysis_data(), error = function(e) rv$analysis_data))
          if (!is.null(df_pred) && !is.null(vars_in_net)) {
            keep <- intersect(vars_in_net, names(df_pred))
            if (length(keep) == length(vars_in_net))
              node_predictability_r2(df_pred[, keep, drop = FALSE])
            else NULL
          } else NULL
        } else NULL

        # Node size by colMeans (same fallback for empty dimnames as GGM block)
        ega_vars <- colnames(net)
        if (is.null(ega_vars) || all(nchar(ega_vars) == 0))
          ega_vars <- if (length(input$selected_vars) == ncol(net)) input$selected_vars else NULL
        df_ega <- group_data(tryCatch(analysis_data(), error = function(e) rv$analysis_data))
        ega_sizes <- if (!is.null(df_ega) && !is.null(ega_vars) && !isTRUE(input$disable_node_scaling)) {
          keep <- intersect(ega_vars, names(df_ega))
          if (length(keep) == ncol(net)) {
            m <- abs(colMeans(df_ega[, keep, drop = FALSE], na.rm = TRUE))
            rng <- range(m, na.rm = TRUE)
            if (diff(rng) < 1e-6) rep(input$qgraph_vsize, ncol(net))
            else input$qgraph_vsize * 0.6 +
                 (m - rng[1]) / diff(rng) * input$qgraph_vsize * 1.2
          } else rep(input$qgraph_vsize, ncol(net))
        } else rep(input$qgraph_vsize, ncol(net))

        # PCA layout: position nodes by their loadings on PC1 and PC2 of the
        # correlation matrix (as in qgraph::PCAnet). Falls back to the selected
        # layout if the data/variables are unavailable.
        ega_lay <- effective_layout()
        if (isTRUE(input$ega_pca_layout)) {
          ega_lay <- tryCatch({
            keep_pca <- if (!is.null(df_ega) && !is.null(ega_vars))
                          intersect(ega_vars, names(df_ega)) else NULL
            if (!is.null(keep_pca) && length(keep_pca) == ncol(net)) {
              cor_pca <- cor(df_ega[, keep_pca, drop = FALSE],
                             use = "pairwise.complete.obs")
              eigen(cor_pca)$vectors[, 1:2, drop = FALSE]  # n × 2 layout matrix
            } else effective_layout()
          }, error = function(e) effective_layout())
        }

        pca_suffix <- if (isTRUE(input$ega_pca_layout)) " | PCA layout" else ""
        do.call(qgraph, c(list(
               input        = net_display,
               labels       = labs,
               groups       = grp_list,
               color        = grp_colors,
               layout       = ega_lay,
               theme        = input$qgraph_theme,
               vsize        = ega_sizes,
               esize        = input$qgraph_esize,
               edge.labels  = input$qgraph_labels_est,
               border.color = input$node_border_color,
               legend       = isTRUE(input$ega_legend),
               legend.cex   = 0.7,
               title        = paste0("EGA \u2014 ", k, " dimension(s) [", algo, pca_suffix, "]"),
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args(),
               pie_args(r2_vec)))

      } else if (type == "cor_heat") {
        df <- group_data(analysis_data())
        req(df)
        cor_mat <- cor(df, use = "pairwise.complete.obs")
        do.call(qgraph, c(list(
               input  = cor_mat,
               graph  = "cor",
               layout = effective_layout(),
               theme  = input$qgraph_theme,
               color        = rep(input$node_color, ncol(cor_mat)),
               vsize        = input$qgraph_vsize,
               border.color = input$node_border_color,
               title        = "Correlation Network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args()))

      } else if (type == "factor_network") {
        lambda <- tryCatch(getmatrix(mod, "lambda"), error = function(e) NULL)
        omega_z <- tryCatch(getmatrix(mod, "omega_zeta"), error = function(e) NULL)

        if (is.null(lambda)) {
          plot.new()
          text(0.5, 0.5, "Lambda matrix not available.\nFit a factor (LNM/RNM) model first.", cex = 1.3)
          return()
        }

        # Collapse multi-group matrices (3D array OR list) to 2D
        lambda  <- collapse_mg(lambda, grp)
        omega_z <- if (!is.null(omega_z)) collapse_mg(omega_z, grp) else matrix(0, ncol(lambda), ncol(lambda))

        # Apply edge cutoff to latent network AND loadings before building W
        cutoff <- input$qgraph_cutoff
        if (!is.null(cutoff) && cutoff > 0) {
          omega_z[abs(omega_z) < cutoff] <- 0
          lambda[abs(lambda)   < cutoff] <- 0
        }

        n_obs   <- nrow(lambda)
        n_lat   <- ncol(lambda)
        n_total <- n_lat + n_obs

        # Weight matrix: latents first, then observed
        W <- matrix(0, n_total, n_total)
        W[1:n_lat, 1:n_lat]                     <- omega_z    # symmetric latent network
        W[1:n_lat, (n_lat + 1):n_total]         <- t(lambda)  # loadings lat -> obs only

        # Per-edge directionality:
        #   latent-latent  → FALSE (undirected, plain GGM lines, no arrowheads)
        #   latent-observed → TRUE  (directed arrows)
        dir_mat <- matrix(FALSE, n_total, n_total)
        dir_mat[1:n_lat, (n_lat + 1):n_total] <- TRUE

        # Labels — trust the fitted model's dimnames as the authoritative source
        lat_labs <- colnames(lambda)   # correct order comes from the fitted model
        if (is.null(lat_labs) || length(lat_labs) == 0 || all(nchar(lat_labs) == 0) ||
            all(grepl("^eta_?\\d+$", lat_labs, ignore.case = TRUE))) {
          # colnames are missing or look like internal psychonetrics names → use UI names
          lat_labs <- if (!is.null(rv$latent_names) && length(rv$latent_names) == n_lat)
                        rv$latent_names
                      else paste0("L", seq_len(n_lat))
        }
        obs_labs <- rownames(lambda)   # correct order comes from the fitted model
        if (is.null(obs_labs) || length(obs_labs) == 0 || all(nchar(obs_labs) == 0)) {
          obs_labs <- if (length(input$selected_vars) == n_obs) input$selected_vars
                      else paste0("V", seq_len(n_obs))
        }

        # Node appearance: circles for latents, squares for observed
        shapes  <- c(rep("circle", n_lat), rep("square", n_obs))
        mean_fs    <- mean_fscores_for_plot(lambda)
        importance <- if (!is.null(mean_fs)) abs(mean_fs) else colMeans(abs(lambda))
        imp_range  <- range(importance)
        default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                         "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
        pal <- sapply(seq_len(n_lat), function(i) {
          col <- input[[paste0("lat_color_", i)]]
          if (!is.null(col) && nchar(col) > 0) col
          else if (i <= length(default_pal)) default_pal[i] else "#999999"
        })
        dominant <- apply(abs(lambda), 1, which.max)
        if (isTRUE(input$disable_node_scaling)) {
          lat_sizes <- rep(input$qgraph_vsize * 1.8, n_lat)
          lat_cols  <- pal[seq_len(n_lat)]
          obs_cols  <- pal[dominant]
        } else {
          lat_sizes <- if (diff(imp_range) < 1e-6) rep(input$qgraph_vsize * 1.8, n_lat)
                       else input$qgraph_vsize * 1.2 +
                            (importance - imp_range[1]) / diff(imp_range) * input$qgraph_vsize * 1.4
          lat_dark  <- norm01(importance)
          lat_cols  <- mapply(darken_hex, pal[seq_len(n_lat)], lat_dark)
          load_str  <- apply(abs(lambda), 1, max)
          obs_dark  <- norm01(load_str)
          obs_cols  <- mapply(darken_hex, pal[dominant], obs_dark)
        }
        vsizes    <- c(lat_sizes, rep(input$qgraph_vsize, n_obs))
        node_cols <- c(lat_cols, obs_cols)

        # Layout: "circle" → inner ring (latents) + outer ring (observed grouped by factor)
        lay <- if (input$qgraph_layout == "circle") {
          angles_lat <- seq(0, 2 * pi, length.out = n_lat + 1)[seq_len(n_lat)]
          lat_pos    <- cbind(cos(angles_lat), sin(angles_lat)) * 0.38
          obs_pos    <- matrix(0, n_obs, 2)
          for (k in seq_len(n_lat)) {
            idx <- which(dominant == k)
            n_k <- length(idx)
            if (n_k == 0) next
            half <- (pi / n_lat) * 0.85
            angs <- if (n_k == 1) angles_lat[k] else
                      seq(angles_lat[k] - half, angles_lat[k] + half, length.out = n_k)
            obs_pos[idx, 1] <- cos(angs)
            obs_pos[idx, 2] <- sin(angs)
          }
          rbind(lat_pos, obs_pos)
        } else {
          effective_layout()
        }

        # Node predictability rings: latents only — no rings on observed (square) nodes
        r2_vec <- if (isTRUE(input$show_predictability)) {
          r2_lat <- latent_predictability_r2(mod, grp)
          if (!is.null(r2_lat) && length(r2_lat) == n_lat)
            c(as.list(r2_lat), vector("list", n_obs))  # NULL slots = no ring on observed squares
          else NULL
        } else NULL

        # Boost near-zero loadings for display so all non-zero cross-loadings are visible.
        # Edge labels will show actual (un-boosted) values from the fitted model.
        lambda_w <- lambda
        max_load <- if (any(lambda_w != 0)) max(abs(lambda_w[lambda_w != 0])) else 0
        if (max_load > 1e-6) {
          min_vis <- max_load * 0.20
          nz      <- lambda_w != 0 & abs(lambda_w) < min_vis
          lambda_w[nz] <- sign(lambda_w[nz]) * min_vis
        }
        W_disp <- W
        W_disp[1:n_lat, (n_lat + 1):n_total] <- t(lambda_w)

        # Edge label matrix: actual fitted values (not boosted), blank where zero
        edge_lab <- matrix("", n_total, n_total)
        edge_lab[1:n_lat, 1:n_lat] <-
          ifelse(omega_z != 0, as.character(round(omega_z, 2)), "")
        edge_lab[1:n_lat, (n_lat + 1):n_total] <-
          ifelse(t(lambda) != 0, as.character(round(t(lambda), 2)), "")

        do.call(qgraph, c(list(
               input        = W_disp,
               labels       = c(lat_labs, obs_labs),
               shape        = shapes,
               vsize        = vsizes,
               color        = node_cols,
               directed     = dir_mat,
               layout       = lay,
               theme        = input$qgraph_theme,
               edge.labels  = if (isTRUE(input$qgraph_labels_est)) edge_lab else FALSE,
               esize        = input$qgraph_esize,
               border.color = input$node_border_color,
               title        = "Factor structure & latent network",
               mar          = c(4, 4, 4, 4)),
               qgraph_label_args(),
               qgraph_color_args(),
               pie_args(r2_vec)))

      } else if (type == "ci_plot") {
        is_mg <- !is.null(input$group_var) && input$group_var != ""

        draw_ci <- function(g = NULL) {
          use_boot <- !is.null(input$ci_source) && input$ci_source == "boot" &&
                      !is.null(rv$boot_agg)
          src_mod  <- if (use_boot) rv$boot_agg else mod

          args <- list(src_mod, input$ci_matrix, split0 = input$ci_split0)
          if (!is.null(g)) args$group <- g
          ok <- tryCatch({ do.call(CIplot, args); TRUE }, error = function(e) FALSE)
          if (ok) return(invisible())

          # ── manual fallback ──────────────────────────────────────────────
          params <- tryCatch(as.data.frame(src_mod %>% parameters), error = function(e) NULL)
          if (is.null(params)) {
            plot.new(); text(0.5, 0.5, "CI plot not available for this model."); return()
          }
          mat_col <- grep("^matrix$", names(params), ignore.case = TRUE, value = TRUE)[1]
          grp_col <- grep("^group$",  names(params), ignore.case = TRUE, value = TRUE)[1]
          est_col <- grep("^est",     names(params), ignore.case = TRUE, value = TRUE)[1]
          lo_col  <- grep("lower|lb|ci_l", names(params), ignore.case = TRUE, value = TRUE)[1]
          hi_col  <- grep("upper|ub|ci_u", names(params), ignore.case = TRUE, value = TRUE)[1]

          if (!is.na(mat_col)) params <- params[params[[mat_col]] == input$ci_matrix, , drop = FALSE]
          if (!is.null(g) && !is.na(grp_col))
            params <- params[params[[grp_col]] == g, , drop = FALSE]

          if (nrow(params) == 0 || is.na(est_col) || is.na(lo_col) || is.na(hi_col)) {
            plot.new()
            text(0.5, 0.5, "No CI data available for this matrix / group.")
            return()
          }

          est <- as.numeric(params[[est_col]])
          lo  <- as.numeric(params[[lo_col]])
          hi  <- as.numeric(params[[hi_col]])

          lab_cols <- intersect(c("var1", "var2"), names(params))
          lab <- if (length(lab_cols) == 2)
            paste0(params[[lab_cols[1]]], "\u2013", params[[lab_cols[2]]])
          else seq_len(nrow(params))

          if (input$ci_split0) {
            sig <- lo > 0 | hi < 0
            ord <- c(which(sig), which(!sig))
          } else {
            ord <- order(est)
          }
          est <- est[ord]; lo <- lo[ord]; hi <- hi[ord]; lab <- lab[ord]

          n <- length(est)
          cols <- ifelse(lo[seq_len(n)] > 0 | hi[seq_len(n)] < 0, "#2980B9", "gray60")
          grp_title <- if (!is.null(g)) paste0(" \u2013 Group ", g) else ""
          par(mar = c(4, max(6, max(nchar(lab)) * 0.35), 3, 2))
          plot(est, seq_len(n),
               xlim = range(c(lo, hi), na.rm = TRUE),
               yaxt = "n", ylab = "", xlab = "Estimate",
               pch = 19, col = cols, cex = 0.85,
               main = paste0("CI Plot: ", input$ci_matrix, grp_title))
          abline(v = 0, lty = 2, col = "gray50")
          arrows(lo, seq_len(n), hi, seq_len(n),
                 angle = 90, code = 3, length = 0.04, col = cols)
          axis(2, at = seq_len(n), labels = lab, las = 1, cex.axis = 0.7)
        }

        if (!is_mg) {
          draw_ci()
        } else {
          sel <- input$plot_group_sel
          if (is.null(sel) || sel == "all") {
            req(rv$raw_data)
            grp_levels <- get_grp_levels()
            n_g   <- length(grp_levels)
            old_par <- par(mfrow = c(1, n_g))
            for (gi in seq_len(n_g)) draw_ci(gi)
            par(old_par)
          } else {
            draw_ci(as.integer(sel))
          }
        }

      } else if (type == "semplot") {
        if (requireNamespace("semPlot", quietly = TRUE)) {
          library(semPlot)
          lambdaEst <- collapse_mg(getmatrix(mod, "lambda"), grp)
          psiEst    <- tryCatch(collapse_mg(getmatrix(mod, "sigma_zeta"), grp),
                                error = function(e) collapse_mg(getmatrix(mod, "sigma"), grp))
          thetaEst  <- tryCatch(collapse_mg(getmatrix(mod, "sigma_epsilon"), grp),
                                error = function(e) diag(ncol(lambdaEst)))

          # psychonetrics matrices have empty dimnames — set proper labels
          # Latent names: trust fitted model colnames, else rv$latent_names, else "L1"...
          n_lat_sem <- ncol(lambdaEst)
          n_obs_sem <- nrow(lambdaEst)
          lat_labs_sem <- colnames(lambdaEst)
          if (is.null(lat_labs_sem) || all(nchar(lat_labs_sem) == 0) ||
              all(grepl("^eta_?\\d+$", lat_labs_sem, ignore.case = TRUE))) {
            lat_labs_sem <- if (!is.null(rv$latent_names) && length(rv$latent_names) == n_lat_sem)
                              rv$latent_names
                            else paste0("L", seq_len(n_lat_sem))
          }
          obs_labs_sem <- rownames(lambdaEst)
          if (is.null(obs_labs_sem) || all(nchar(obs_labs_sem) == 0)) {
            obs_labs_sem <- if (length(input$selected_vars) == n_obs_sem)
                              input$selected_vars
                            else paste0("V", seq_len(n_obs_sem))
          }
          dimnames(lambdaEst) <- list(obs_labs_sem, lat_labs_sem)
          dimnames(psiEst)    <- list(lat_labs_sem, lat_labs_sem)
          dimnames(thetaEst)  <- list(obs_labs_sem, obs_labs_sem)

          sem_mod <- lisrelModel(LY = lambdaEst, PS = psiEst, TE = thetaEst)
          semPaths(sem_mod, what = "std", whatLabels = "est",
                   layout = "tree2",
                   theme = "colorblind", sizeLat = 10, sizeMan = 5,
                   edge.label.cex = 0.8, mar = c(6, 1, 6, 1),
                   pastel = TRUE, borders = TRUE)
        } else {
          plot.new()
          text(0.5, 0.5, "Install the 'semPlot' package for path diagrams.", cex = 1.3)
        }
      }
    }

    plot_fn()
  }

  # ── Group selector helper ─────────────────────────────────────────────────
  # Returns the integer group index to pass to make_plot() (or 1L if no group).
  cur_grp <- reactive({
    sel <- input$plot_group_sel
    if (is.null(sel) || sel == "all") 1L else as.integer(sel)
  })

  output$main_plot <- renderPlot({ make_plot(grp = cur_grp()) })

  # Per-latent colour pickers — appear once latent names are defined
  output$latent_color_ui <- renderUI({
    lat_names <- rv$latent_names
    if (is.null(lat_names) || length(lat_names) == 0) return(NULL)
    default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
                     "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
    tagList(
      hr(),
      tags$b("Latent variable colours:"),
      lapply(seq_along(lat_names), function(i) {
        def_col <- if (i <= length(default_pal)) default_pal[i] else "#999999"
        colourpicker::colourInput(
          paste0("lat_color_", i), label = lat_names[i],
          value = def_col, showColour = "both", palette = "square"
        )
      })
    )
  })

  # Render the selectInput only when a group variable is active
  output$plot_group_ui <- renderUI({
    grp_var <- input$group_var
    if (is.null(grp_var) || grp_var == "") return(NULL)
    grp_levels <- get_grp_levels()
    if (length(grp_levels) == 0) return(NULL)
    n <- length(grp_levels)
    named <- setNames(as.character(seq_len(n)),
                      paste0("Group ", seq_len(n), ": ", grp_levels))
    tagList(
      hr(),
      selectInput("plot_group_sel", "Group to display:",
                  choices = c(c("All Groups" = "all"), named),
                  selected = "all")
    )
  })

  # Swap between a single plotOutput and side-by-side outputs
  output$plot_output_area <- renderUI({
    show_all <- !is.null(input$group_var) && input$group_var != "" &&
                !is.null(input$plot_group_sel) && input$plot_group_sel == "all"
    if (show_all) {
      grp_levels <- get_grp_levels()
      n         <- length(grp_levels)
      if (n == 0) return(plotOutput("main_plot", height = paste0(input$plot_display_height, "px")))
      col_width <- max(4, min(6, floor(12 / n)))
      plot_list <- lapply(seq_len(n), function(i) {
        pname <- paste0("grp_plot_", i)
        column(col_width,
               h5(paste("Group:", grp_levels[i])),
               plotOutput(pname, height = paste0(input$plot_display_height, "px")))
      })
      fluidRow(do.call(tagList, plot_list))
    } else {
      plotOutput("main_plot", height = paste0(input$plot_display_height, "px"))
    }
  })

  # Dynamically register one renderPlot per group when "All Groups" is chosen
  observe({
    grp_var <- input$group_var
    sel     <- input$plot_group_sel
    if (is.null(grp_var) || grp_var == "" || is.null(sel) || sel != "all") return()
    grp_levels <- get_grp_levels()
    for (i in seq_along(grp_levels)) {
      local({
        gi <- i
        output[[paste0("grp_plot_", gi)]] <- renderPlot({ make_plot(grp = gi) })
      })
    }
  })

  output$download_plot <- downloadHandler(
    filename = function() paste0("psychonetrics_plot_", Sys.Date(), ".pdf"),
    content  = function(file) {
      pdf(file, width = input$plot_width, height = input$plot_height)
      if (!is.null(input$group_var) && input$group_var != "" &&
          !is.null(input$plot_group_sel) && input$plot_group_sel == "all") {
        grp_levels <- get_grp_levels()
        for (gi in seq_along(grp_levels)) make_plot(grp = gi)
      } else {
        make_plot(grp = cur_grp())
      }
      dev.off()
    }
  )

  # ── Advanced: second model + compare ─────────────────────────────────────
  observeEvent(input$run_model2, {
    req(rv$model)

    # If model2 was loaded from file, it's already set — just refresh output
    if (isTRUE(input$model2_source == "load")) {
      req(rv$model2)
      return()
    }

    df <- analysis_data()
    req(df)
    rv$model2 <- NULL  # clear stale result immediately

    ordered_vars <- if (input$data_type %in% c("ordered", "dichotomous")) {
      input$selected_vars
    } else character(0)

    est <- input$estimator_type
    grp <- if (!is.null(input$group_var) && input$group_var != "") input$group_var else NULL

    # Attach group column to df (same as run_model does)
    if (!is.null(grp)) {
      gv <- rv$raw_data[[grp]]
      if (!is.null(gv)) df[[grp]] <- gv[seq_len(nrow(df))]
    }

    # Determine pruning settings: override or inherit from Model 1 controls
    use_prune        <- if (isTRUE(input$model2_override_prune)) isTRUE(input$model2_do_prune)       else isTRUE(input$do_prune)
    use_prune_alpha  <- if (isTRUE(input$model2_override_prune)) input$model2_prune_alpha            else input$prune_alpha
    use_prune_adjust <- if (isTRUE(input$model2_override_prune)) input$model2_prune_adjust           else input$prune_adjust
    use_stepup       <- if (isTRUE(input$model2_override_prune)) isTRUE(input$model2_do_stepup)      else isTRUE(input$do_stepup)
    use_modelsearch  <- if (isTRUE(input$model2_override_prune)) isTRUE(input$model2_do_modelsearch) else isTRUE(input$do_modelsearch)

    # Group invariance: constrain selected matrices across groups
    groupEqual_arg <- if (!is.null(grp) && isTRUE(input$model2_constrain))
                        input$model2_groupEqual else character(0)

    withProgress(message = "Fitting second model…", value = 0.3, {
      mod2 <- tryCatch({
        if (input$model2_family == "ggm") {
          m <- ggm(df, vars = input$selected_vars, omega = "full", estimator = est,
                   ordered = ordered_vars, groupvar = grp)
          m <- do_runmodel(m)
          if (use_prune)       m <- m %>% prune(alpha = use_prune_alpha, adjust = use_prune_adjust)
          if (use_stepup)      m <- m %>% stepup
          if (use_modelsearch) m <- m %>% modelsearch
          m
        } else {
          lambda  <- get_lambda()
          latents <- colnames(lambda)
          vars_   <- rownames(lambda)

          builder <- switch(input$model2_family,
            "lvm"  = lvm, "lnm" = lnm,
            "rnm"  = function(...) rnm(...,  omega_epsilon = "full"),
            "lrnm" = function(...) lrnm(..., omega_epsilon = "full"))

          ge_args <- if (length(groupEqual_arg) > 0) list(groupEqual = groupEqual_arg) else list()
          m <- do.call(builder,
                       c(list(data = df, lambda = lambda, vars = vars_,
                              latents = latents, identification = input$model2_identification,
                              ordered = ordered_vars, estimator = est, groupvar = grp),
                         ge_args))
          m <- do_runmodel(m)
          if (use_prune)       m <- m %>% prune(alpha = use_prune_alpha, adjust = use_prune_adjust)
          if (use_stepup)      m <- m %>% stepup
          if (use_modelsearch) m <- m %>% modelsearch
          m
        }
      }, error = function(e) {
        showNotification(paste("Model 2 error:", e$message), type = "error")
        NULL
      })
      setProgress(1)
    })

    if (!is.null(mod2)) {
      rv$model2       <- mod2
      rv$model2_label <- toupper(input$model2_family)
      updateTextInput(session, "model2_label", value = rv$model2_label)
      showNotification("Second model fitted!", type = "message")
    }
  })

  .pn_compare <- getFromNamespace("compare", "psychonetrics")

  # ── Save Model 1 as .rds ──────────────────────────────────────────────────
  output$save_model1 <- downloadHandler(
    filename = function() {
      lbl  <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$model1_label))
      if (nchar(lbl) == 0) lbl <- "Model1"
      paste0(lbl, "_", Sys.Date(), ".rds")
    },
    content = function(file) saveRDS(rv$model, file)
  )

  # ── Load Model 2 from .rds ────────────────────────────────────────────────
  observeEvent(input$load_model2_rds, {
    req(input$load_model2_rds)
    loaded <- tryCatch(
      readRDS(input$load_model2_rds$datapath),
      error = function(e) {
        showNotification(paste("Load error:", e$message), type = "error")
        NULL
      }
    )
    if (!is.null(loaded)) {
      rv$model2 <- loaded
      nm <- tools::file_path_sans_ext(input$load_model2_rds$name)
      nm <- sub("_\\d{4}-\\d{2}-\\d{2}$", "", nm)   # strip date suffix
      rv$model2_label <- nm
      updateTextInput(session, "model2_label", value = nm)
      showNotification("Model loaded successfully!", type = "message")
    }
  })

  output$compare2_output <- renderPrint({
    req(rv$model, rv$model2)
    lbl1 <- trimws(input$model1_label); if (nchar(lbl1) == 0) lbl1 <- "Model1"
    lbl2 <- trimws(input$model2_label); if (nchar(lbl2) == 0) lbl2 <- "Model2"
    do.call(.pn_compare, setNames(list(rv$model, rv$model2), c(lbl1, lbl2)))
  })

  output$compare_output <- renderPrint({
    req(rv$model)
    if (!is.null(rv$model2)) {
      lbl1 <- trimws(input$model1_label); if (nchar(lbl1) == 0) lbl1 <- "Model1"
      lbl2 <- trimws(input$model2_label); if (nchar(lbl2) == 0) lbl2 <- "Model2"
      do.call(.pn_compare, setNames(list(rv$model, rv$model2), c(lbl1, lbl2)))
    } else {
      cat("Fit a second model in the Advanced tab to compare.\n\n")
      cat("Current model fit:\n")
      rv$model %>% fit
    }
  })

  # ── DWLS + ordered warning (Variable Setup tab) ──────────────────────────
  output$dwls_ordered_warn <- renderUI({
    est  <- input$estimator_type
    type <- input$data_type
    if (!is.null(est) && !is.null(type) &&
        est == "DWLS" && type %in% c("ordered", "dichotomous")) {
      div(class = "alert alert-warning", style = "margin-top:6px; padding:6px 10px;",
          icon("exclamation-triangle"), " ",
          tags$b("DWLS + ordered data may fail."),
          tags$br(),
          "DWLS requires computing polychoric correlations and their asymptotic weight matrix.
           This fails when: (a) some bivariate response-category cells are empty,
           (b) thresholds are extreme (floor/ceiling effects), (c) N is small relative
           to the number of items, or (d) some polychoric r is near ±1.",
          tags$br(),
          tags$b("Alternative: "), "try ", tags$code("ULS"), " — it fits ordinal models
           without needing the full weight matrix and is much more stable.")
    }
  })

  # ── Estimator warning (Advanced tab) ────────────────────────────────────
  output$compare2_estimator_warn <- renderUI({
    req(rv$model2)
    est <- input$estimator_type
    if (!est %in% c("default", "ML", "FIML")) {
      div(class = "alert alert-warning", style = "margin-top:8px;",
          icon("exclamation-triangle"), " ",
          tags$b("AIC/BIC unavailable:"),
          sprintf(" The '%s' estimator does not produce a log-likelihood. ", est),
          "Switch to ML or FIML in Variable Setup to enable information criteria.")
    }
  })

  # ── Individual fit() blocks below comparison ─────────────────────────────
  .fmt_fit <- function(mod, lbl) {
    cat(paste0("── ", lbl, " (individual fit) ──\n"))
    tryCatch(print(mod %>% fit), error = function(e) cat("Fit not available.\n"))
  }

  output$compare2_fit_ui <- renderUI({
    req(rv$model, rv$model2)
    tagList(
      hr(),
      h5("Individual model fit:"),
      verbatimTextOutput("compare2_fit1"),
      verbatimTextOutput("compare2_fit2")
    )
  })

  output$compare2_fit1 <- renderPrint({
    req(rv$model, rv$model2)
    lbl1 <- trimws(input$model1_label); if (nchar(lbl1) == 0) lbl1 <- "Model1"
    .fmt_fit(rv$model, lbl1)
  })

  output$compare2_fit2 <- renderPrint({
    req(rv$model, rv$model2)
    lbl2 <- trimws(input$model2_label); if (nchar(lbl2) == 0) lbl2 <- "Model2"
    .fmt_fit(rv$model2, lbl2)
  })

  output$compare_fit_ui <- renderUI({
    req(rv$model, rv$model2)
    tagList(
      hr(),
      h5("Individual model fit:"),
      verbatimTextOutput("compare_fit1"),
      verbatimTextOutput("compare_fit2")
    )
  })

  output$compare_fit1 <- renderPrint({
    req(rv$model, rv$model2)
    lbl1 <- trimws(input$model1_label); if (nchar(lbl1) == 0) lbl1 <- "Model1"
    .fmt_fit(rv$model, lbl1)
  })

  output$compare_fit2 <- renderPrint({
    req(rv$model, rv$model2)
    lbl2 <- trimws(input$model2_label); if (nchar(lbl2) == 0) lbl2 <- "Model2"
    .fmt_fit(rv$model2, lbl2)
  })

  # ── Bootstrapping ────────────────────────────────────────────────────────
  observeEvent(input$run_boot, {
    req(rv$model)
    showNotification("Bootstrapping started – this may take a while…", type = "message")

    withProgress(message = "Bootstrapping…", value = 0.1, {
      boot_res <- tryCatch({
        # Re-build model call with storedata = TRUE for bootstrap
        # Use loop_psychonetrics for parallel bootstrapping
        # We need to reconstruct the model with bootstrap = "nonparametric"
        df <- analysis_data()
        est <- input$estimator_type
        ordered_vars <- if (input$data_type %in% c("ordered", "dichotomous")) {
          input$selected_vars
        } else character(0)

        if (input$model_family == "ggm") {
          bootstraps <- loop_psychonetrics({
            ggm(df, vars = input$selected_vars, estimator = est,
                ordered = ordered_vars, bootstrap = "nonparametric") %>%
              runmodel %>%
              prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          }, reps = input$boot_reps, nCores = input$boot_cores)

          sample_mod <- rv$model
          aggregate_bootstraps(sample = sample_mod, bootstraps = bootstraps)

        } else {
          lambda  <- get_lambda()
          latents <- colnames(lambda)
          vars_   <- rownames(lambda)
          fam     <- input$model_family
          ident   <- input$identification

          builder <- switch(fam,
            "lvm" = lvm, "lnm" = lnm, "rnm" = rnm, "lrnm" = lrnm)

          bootstraps <- loop_psychonetrics({
            builder(data = df, lambda = lambda, vars = vars_,
                    latents = latents, identification = ident,
                    ordered = ordered_vars, estimator = est,
                    bootstrap = "nonparametric") %>%
              runmodel %>%
              prune(alpha = input$prune_alpha, adjust = input$prune_adjust)
          }, reps = input$boot_reps, nCores = input$boot_cores)

          aggregate_bootstraps(sample = rv$model, bootstraps = bootstraps)
        }
      }, error = function(e) {
        showNotification(paste("Bootstrap error:", e$message), type = "error")
        NULL
      })
      setProgress(1)
    })

    if (!is.null(boot_res)) {
      rv$boot_agg <- boot_res
      showNotification("Bootstrap complete!", type = "message")
    }
  })

  output$boot_output <- renderPrint({
    if (!is.null(rv$boot_agg)) {
      parameters(rv$boot_agg)
    } else {
      cat("No bootstrap results yet. Click 'Run Bootstrap' above.\n")
    }
  })

  # ── R Code generation ───────────────────────────────────────────────────
  # Called from the fit observers to snapshot the exact settings used.
  build_r_code <- function() {
    vars <- input$selected_vars
    if (is.null(vars) || length(vars) < 2) return("# Select variables first")

    est     <- input$estimator_type
    opt     <- input$optimizer
    fam     <- input$model_family
    grp     <- if (!is.null(input$group_var) && input$group_var != "") input$group_var else NULL
    has_grp <- !is.null(grp)

    # Helpers used repeatedly below
    q <- function(...) paste0('"', ..., '"')
    vec_str <- function(x) paste0('c("', paste(x, collapse = '", "'), '")')
    indent  <- function(s, n = 2) gsub("(?m)^", strrep(" ", n), s, perl = TRUE)

    # ── Collect context ────────────────────────────────────────────────────────
    ordered_str <- if (input$data_type %in% c("ordered", "dichotomous"))
                     paste0('ordered = ', vec_str(vars))
                   else 'ordered = character(0)'

    grp_str <- if (has_grp) paste0(', groupvar = "', grp, '"') else ""

    tr_choice    <- if (!is.null(input$data_transform)) input$data_transform else "none"
    tr_dt        <- if (!is.null(input$data_type)) input$data_type else "continuous"
    tr_is_dichot <- tr_choice %in% c("dichot_mean", "dichot_median")
    tr_applicable <- tr_choice != "none" && tr_dt != "dichotomous" &&
                     (tr_is_dichot || tr_dt == "continuous")

    lam_for_code <- if (fam %in% c("lvm","lnm","rnm","lrnm","dlvm1","ri_clpm"))
                      tryCatch(get_lambda(), error = function(e) NULL) else NULL
    latents      <- if (!is.null(lam_for_code)) colnames(lam_for_code)
                    else if (!is.null(rv$latent_names)) rv$latent_names
                    else character(0)
    has_lam      <- !is.null(lam_for_code) && nrow(lam_for_code) > 0

    has_model2 <- !is.null(rv$model2)
    has_boot   <- !is.null(rv$boot_agg)

    # ── Group levels (for plot loops) ─────────────────────────────────────────
    grp_levels_vec <- if (has_grp) get_grp_levels() else character(0)
    n_grp          <- max(1L, length(grp_levels_vec))

    # ── Transformation snippet ─────────────────────────────────────────────────
    transform_code <- if (tr_applicable) {
      tr_line <- switch(tr_choice,
        npn = paste0(
          'df[vars] <- as.data.frame(apply(as.matrix(df[vars]), 2, function(x) {\n',
          '  n <- sum(!is.na(x))\n',
          '  delta <- 1/(4*n^0.25*sqrt(pi*log(n)))\n',
          '  qnorm(pmin(pmax(rank(x, ties.method="average", na.last="keep")/n,\n',
          '             delta), 1-delta))\n',
          '}))'),
        dichot_mean   = 'df[vars] <- as.data.frame(lapply(df[vars], function(x) as.numeric(x >= mean(x, na.rm=TRUE))))',
        dichot_median = 'df[vars] <- as.data.frame(lapply(df[vars], function(x) as.numeric(x >= median(x, na.rm=TRUE))))',
        log1p  = 'df[vars] <- as.data.frame(log1p(pmax(as.matrix(df[vars]), 0)))',
        sqrt   = 'df[vars] <- as.data.frame(sqrt(pmax(as.matrix(df[vars]), 0)))',
        zscore = 'df[vars] <- as.data.frame(scale(df[vars]))',
        NULL
      )
      tr_label <- switch(tr_choice,
        npn="Nonparanormal (rank-based normal scores)", dichot_mean="Dichotomize at mean",
        dichot_median="Dichotomize at median", log1p="Log (log1p)",
        sqrt="Square root", zscore="Z-score", tr_choice)
      if (!is.null(tr_line)) paste0('# Transformation: ', tr_label, '\n', tr_line, '\n\n') else ""
    } else ""

    # ── Lambda block (shared by LVM family & panel models) ────────────────────
    make_lam_code <- function(lam, lat_nms, row_names = "vars") {
      lc <- "Lambda <- matrix(c(\n"
      for (i in seq_len(nrow(lam)))
        lc <- paste0(lc, "  ", paste(lam[i,], collapse=", "),
                     if (i < nrow(lam)) ",\n" else "\n")
      paste0(lc, '), nrow=', nrow(lam), ', ncol=', ncol(lam), ', byrow=TRUE,\n',
             '  dimnames=list(', row_names, ', ', vec_str(lat_nms), '))\n\n')
    }

    lam_code <- if (has_lam) make_lam_code(lam_for_code, latents) else ""

    # ── Panel wave structure ──────────────────────────────────────────────────
    wave_str <- if (!is.null(rv$vars_by_wave)) {
      m <- rv$vars_by_wave
      if (is.matrix(m)) {
        rows <- apply(m, 1, function(r) paste0('  c("', paste(r, collapse='", "'), '")'))
        paste0("vars_panel <- rbind(\n", paste(rows, collapse=",\n"), "\n)\n")
      } else {
        items <- lapply(names(m), function(nm)
          paste0('  ', nm, ' = c("', paste(m[[nm]], collapse='", "'), '")'))
        paste0("vars_panel <- list(\n", paste(items, collapse=",\n"), "\n)\n")
      }
    } else "# vars_panel <- list(wave1=c(...), wave2=c(...))\n"

    # ── Helper functions to embed in generated code ────────────────────────────
    # collapse_mg: needed whenever multi-group model is used
    collapse_mg_fn <- if (has_grp) paste0(
      'collapse_mg <- function(m, grp=1L) {\n',
      '  if (is.list(m) && !is.data.frame(m)) return(as.matrix(m[[grp]]))\n',
      '  if (is.array(m) && length(dim(m))==3) return(as.matrix(m[,,grp]))\n',
      '  as.matrix(m)\n',
      '}\n\n') else ""

    # johnson_rw: needed for RI network plots
    needs_ri   <- fam %in% c("lnm","lrnm","lvm")
    rw_fn <- if (needs_ri) paste0(
      'johnson_rw_from_cor <- function(R) {\n',
      '  p <- ncol(R); nm <- colnames(R)\n',
      '  if (is.null(nm)) nm <- paste0("L", seq_len(p))\n',
      '  RI <- matrix(0, p, p, dimnames=list(nm,nm))\n',
      '  r2s <- function(sub,rxy,Rxx) {\n',
      '    if (!length(sub)) return(0)\n',
      '    tryCatch(max(0,drop(t(rxy[sub])%*%solve(Rxx[sub,sub,drop=FALSE],rxy[sub]))),\n',
      '             error=function(e) 0)\n',
      '  }\n',
      '  for (i in seq_len(p)) {\n',
      '    pred <- seq_len(p)[-i]; q <- length(pred)\n',
      '    if (!q) next\n',
      '    rxy <- R[pred,i]; Rxx <- R[pred,pred,drop=FALSE]\n',
      '    if (q==1L) { RI[pred,i] <- rxy^2; next }\n',
      '    cc <- numeric(q)\n',
      '    for (j in seq_len(q)) {\n',
      '      oth <- seq_len(q)[-j]\n',
      '      for (mask in seq_len(2^length(oth))-1L) {\n',
      '        S  <- if (!length(oth)) integer(0) else oth[as.logical(intToBits(mask)[seq_len(length(oth))])]\n',
      '        wt <- factorial(length(S))*factorial(q-1L-length(S))/factorial(q)\n',
      '        cc[j] <- cc[j]+wt*(r2s(c(S,j),rxy,Rxx)-r2s(S,rxy,Rxx))\n',
      '      }\n',
      '    }\n',
      '    RI[pred,i] <- cc\n',
      '  }\n',
      '  RI\n',
      '}\n\n') else ""

    # ── qgraph option variables ────────────────────────────────────────────────
    qg_opts_code <- paste0(
      '# ── qgraph visual options (edit as needed) ──────────────────────────\n',
      'qg_layout  <- "', if (!is.null(input$qgraph_layout))  input$qgraph_layout  else "circle", '"\n',
      'qg_theme   <- "', if (!is.null(input$qgraph_theme))   input$qgraph_theme   else "colorblind", '"\n',
      'qg_vsize   <-  ', if (!is.null(input$qgraph_vsize))   input$qgraph_vsize   else 5, '\n',
      'qg_esize   <-  ', if (!is.null(input$qgraph_esize))   input$qgraph_esize   else 5, '\n',
      'qg_cutoff  <-  ', if (!is.null(input$qgraph_cutoff))  input$qgraph_cutoff  else 0, '\n',
      'qg_elabels <- ', if (isTRUE(input$qgraph_labels_est)) "TRUE" else "FALSE", '\n\n'
    )

    # ── Model-fitting code ─────────────────────────────────────────────────────
    fit_code <- if (fam == "ggm") {
      paste0(
        '# ── Gaussian Graphical Model ─────────────────────────────────────────\n',
        'mod <- ggm(df, vars=vars, omega="', input$ggm_start, '",\n',
        '           estimator="', est, '", ', ordered_str, grp_str, ')\n',
        'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')

    } else if (fam == "ising") {
      paste0(
        '# ── Ising Model ──────────────────────────────────────────────────────\n',
        'mod <- Ising(df, vars=vars, estimator="', est, '")\n',
        'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')

    } else if (fam %in% c("lvm","lnm","rnm","lrnm")) {
      omega_eps <- if (fam %in% c("rnm","lrnm")) ',\n         omega_epsilon="full"' else ""
      paste0(
        '# ── ', toupper(fam), ' ────────────────────────────────────────────────────────\n',
        'latents <- ', vec_str(latents), '\n\n',
        lam_code,
        'mod <- ', fam, '(df, lambda=Lambda, vars=vars, latents=latents,\n',
        '         identification="', input$identification, '",\n',
        '         estimator="', est, '", ', ordered_str, grp_str, omega_eps, ')\n',
        'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')

    } else if (fam == "dlvm1") {
      paste0(
        '# ── Panel DLVM ───────────────────────────────────────────────────────\n',
        'latents <- ', vec_str(latents), '\n\n',
        lam_code,
        wave_str, '\n',
        'mod <- dlvm1(df, vars=vars_panel, lambda=Lambda, latents=latents,\n',
        '             beta="full", within_latent="', input$within_latent, '",\n',
        '             between_latent="', input$between_latent, '",\n',
        '             identification="', input$identification, '",\n',
        '             estimator="', est, '")\n',
        'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')

    } else if (fam == "panelgvar") {
      paste0(
        '# ── Panel GVAR ───────────────────────────────────────────────────────\n',
        wave_str, '\n',
        'mod <- panelgvar(df, vars=vars_panel,\n',
        '                 within_latent="', input$within_latent, '",\n',
        '                 estimator="', est, '")\n',
        'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')

    } else if (fam == "ri_clpm") {
      lam_blk <- if (has_lam)
        paste0('latents <- ', vec_str(latents), '\n\n', lam_code,
               'mod <- ri_clpm(df, vars=vars_panel, lambda=Lambda,\n',
               '               type="', input$within_latent, '", estimator="', est, '")\n')
      else
        paste0('mod <- ri_clpm(df, vars=vars_panel,\n',
               '               type="', input$within_latent, '", estimator="', est, '")\n')
      paste0('# ── RI-CLPM ──────────────────────────────────────────────────────────\n',
             wave_str, '\n', lam_blk,
             'mod <- mod %>% setoptimizer("', opt, '") %>% runmodel()\n')
    } else ""

    # ── Pruning / search ───────────────────────────────────────────────────────
    prune_code <- ""
    if (input$do_prune) {
      adj <- if (!is.null(input$prune_adjust) && input$prune_adjust != "none")
               paste0(', adjust="', input$prune_adjust, '"') else ""
      prune_code <- paste0('mod <- mod %>% prune(alpha=', input$prune_alpha, adj, ')\n')
    }
    if (isTRUE(input$do_stepup))      prune_code <- paste0(prune_code, 'mod <- mod %>% stepup\n')
    if (isTRUE(input$do_modelsearch)) prune_code <- paste0(prune_code, 'mod <- mod %>% modelsearch\n')

    # ── Results: fit / parameters / MIs ───────────────────────────────────────
    results_code <- paste0(
      '\n# ── Fit indices ──────────────────────────────────────────────────────\n',
      'mod %>% fit\n\n',
      '# ── Parameters ──────────────────────────────────────────────────────\n',
      'mod %>% parameters\n\n',
      '# ── Modification indices ─────────────────────────────────────────────\n',
      'mod %>% MIs\n\n')

    # ── Matrix extraction ─────────────────────────────────────────────────────
    mg_wrap <- function(body, indent_n = 2) {
      if (!has_grp) return(body)
      paste0(
        'for (gi in seq_len(', n_grp, ')) {  # group: ', paste(grp_levels_vec, collapse=", "), '\n',
        paste0(strsplit(trimws(body, "right"), "\n")[[1]] %>%
                 { paste0(strrep(" ", indent_n), ., collapse="\n") }),
        '\n}\n')
    }

    extract_code <- paste0('# ── Matrix extraction ───────────────────────────────────────────────\n')
    if (has_grp) extract_code <- paste0(extract_code,
      '# collapse_mg() selects the 2D slice for group gi from a 3D array / list\n')

    if (fam %in% c("lnm","lrnm","lvm")) {
      inner <- if (has_grp)
        'omega_lat <- collapse_mg(getmatrix(mod, "omega_zeta"), gi)\n'
      else
        'omega_lat <- as.matrix(getmatrix(mod, "omega_zeta"))\n'
      extract_code <- paste0(extract_code, mg_wrap(inner))
    }
    if (fam %in% c("rnm","lrnm")) {
      inner <- if (has_grp)
        'omega_res <- collapse_mg(getmatrix(mod, "omega_epsilon"), gi)\n'
      else
        'omega_res <- as.matrix(getmatrix(mod, "omega_epsilon"))\n'
      extract_code <- paste0(extract_code, mg_wrap(inner))
    }
    if (fam %in% c("ggm","ising")) {
      inner <- if (has_grp)
        'omega <- collapse_mg(getmatrix(mod, "omega"), gi)\n'
      else
        'omega <- as.matrix(getmatrix(mod, "omega"))\n'
      extract_code <- paste0(extract_code, mg_wrap(inner))
    }
    if (fam %in% c("lvm","lnm","rnm","lrnm")) {
      inner <- if (has_grp)
        'lambda_est <- collapse_mg(getmatrix(mod, "lambda"), gi)\n'
      else
        'lambda_est <- as.matrix(getmatrix(mod, "lambda"))\n'
      extract_code <- paste0(extract_code, mg_wrap(inner))
    }
    if (fam %in% c("dlvm1","panelgvar","ri_clpm")) {
      inner <- paste0(
        if (has_grp) 'beta_mat   <- collapse_mg(getmatrix(mod, "beta"), gi)\n'
        else         'beta_mat   <- as.matrix(getmatrix(mod, "beta"))\n',
        if (has_grp)
          'within_net <- tryCatch(collapse_mg(getmatrix(mod,"omega_zeta_within"),gi),\n                       error=function(e) collapse_mg(getmatrix(mod,"omega_zeta"),gi))\n'
        else
          'within_net <- tryCatch(as.matrix(getmatrix(mod,"omega_zeta_within")),\n                       error=function(e) as.matrix(getmatrix(mod,"omega_zeta")))\n')
      extract_code <- paste0(extract_code, mg_wrap(inner))
    }
    extract_code <- paste0(extract_code, "\n")

    # ── Factor scores (LVM family) ────────────────────────────────────────────
    fscores_code <- ""
    if (fam %in% c("lvm","lnm","rnm","lrnm")) {
      fs_inner <- paste0(
        if (has_grp) {
          paste0(
            'lam_g   <- collapse_mg(getmatrix(mod,"lambda"), gi)\n',
            'sig_z_g <- tryCatch(collapse_mg(getmatrix(mod,"sigma_zeta"),gi),\n',
            '                    error=function(e) diag(ncol(lam_g)))\n',
            'in_g    <- as.character(df[["', grp, '"]]) == grp_levels[gi]\n',
            'data_g  <- df[!is.na(in_g) & in_g, rownames(lam_g), drop=FALSE]\n',
            'data_g  <- data_g[complete.cases(data_g),]\n')
        } else {
          paste0(
            'lam_g   <- as.matrix(getmatrix(mod,"lambda"))\n',
            'sig_z_g <- tryCatch(as.matrix(getmatrix(mod,"sigma_zeta")),\n',
            '                    error=function(e) diag(ncol(lam_g)))\n',
            'data_g  <- df[complete.cases(df[,vars]),vars,drop=FALSE]\n')
        },
        'S_inv      <- tryCatch(solve(cov(data_g)), error=function(e) corpcor::pseudoinverse(cov(data_g)))\n',
        'fscores_g  <- as.matrix(data_g) %*% S_inv %*% lam_g %*% sig_z_g\n',
        'colnames(fscores_g) <- colnames(lam_g)\n',
        if (has_grp) 'cat("Factor scores group", gi, "— first 6 rows:\\n"); print(head(fscores_g))\n'
        else         'print(head(fscores_g))\n'
      )
      fscores_code <- paste0(
        '# ── Factor scores (regression method) ───────────────────────────────\n',
        if (has_grp) paste0(
          'grp_levels <- sort(unique(na.omit(df[["', grp, '"]]))) # ', paste(grp_levels_vec, collapse=", "), '\n',
          'for (gi in seq_along(grp_levels)) {\n',
          indent(fs_inner, 2),
          '}\n\n')
        else paste0(fs_inner, "\n")
      )
    }

    # ── Plot helpers (always embedded when needed) ────────────────────────────
    plot_helpers <- paste0(
      '# ── Plot helpers ────────────────────────────────────────────────────\n',
      'darken_hex <- function(hex, amt) {\n',
      '  v <- col2rgb(hex)/255; rgb(v[1]*amt,v[2]*amt,v[3]*amt)\n}\n',
      'norm01 <- function(x, mn=0.75) {\n',
      '  if (diff(range(x))<1e-6) return(rep(0.75,length(x)))\n',
      '  mn+(1-mn)*(x-min(x))/diff(range(x))\n}\n',
      'default_pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",\n',
      '                 "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")\n\n')

    # ── Build individual plot blocks ─────────────────────────────────────────
    # Helper: make a simple network qgraph call string
    simple_qg <- function(mat_expr, labs_expr, directed="FALSE",
                          title_expr='"Network"', extra="") {
      paste0(
        'net <- ', mat_expr, '\n',
        'if (qg_cutoff>0) net[abs(net)<qg_cutoff] <- 0\n',
        'qgraph(net, labels=', labs_expr, ', layout=qg_layout,\n',
        '       theme=qg_theme, directed=', directed, ', vsize=qg_vsize,\n',
        '       esize=qg_esize, edge.labels=qg_elabels,\n',
        '       title=', title_expr, ', mar=c(4,4,4,4))\n',
        extra)
    }

    grp_loop <- function(body) {
      if (!has_grp) return(body)
      paste0(
        'for (gi in seq_len(', n_grp, ')) {\n',
        '  grp_lbl <- ', if (length(grp_levels_vec)>0) paste0('c(', paste0('"', grp_levels_vec, '"', collapse=", "), ')[gi]') else '"1"', '\n',
        indent(body, 2),
        '}\n')
    }

    plots_code <- '\n# ── Plots ────────────────────────────────────────────────────────────\n'

    # ── Latent network ────────────────────────────────────────────────────────
    if (fam %in% c("lnm","lrnm","lvm")) {
      mat_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"omega_zeta"), gi)'
               else         'as.matrix(getmatrix(mod,"omega_zeta"))'
      title_e <- if (has_grp) 'paste0("Latent Network — Group ", grp_lbl)'
                 else         '"Latent Network"'
      plots_code <- paste0(plots_code,
        '# Latent network (omega_zeta)\n',
        grp_loop(simple_qg(mat_e, 'latents', title_expr=title_e)), '\n')
    }

    # ── Residual network ──────────────────────────────────────────────────────
    if (fam %in% c("rnm","lrnm")) {
      mat_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"omega_epsilon"), gi)'
               else         'as.matrix(getmatrix(mod,"omega_epsilon"))'
      title_e <- if (has_grp) 'paste0("Residual Network — Group ", grp_lbl)'
                 else         '"Residual Network"'
      plots_code <- paste0(plots_code,
        '# Residual network (omega_epsilon)\n',
        grp_loop(simple_qg(mat_e, 'vars', title_expr=title_e)), '\n')
    }

    # ── GGM / Ising network ───────────────────────────────────────────────────
    if (fam %in% c("ggm","ising")) {
      mat_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"omega"), gi)'
               else         'as.matrix(getmatrix(mod,"omega"))'
      ttl   <- if (fam=="ggm") "GGM Network" else "Ising Network"
      title_e <- if (has_grp) paste0('paste0("', ttl, ' — Group ", grp_lbl)')
                 else         paste0('"', ttl, '"')
      plots_code <- paste0(plots_code,
        '# ', ttl, '\n',
        grp_loop(simple_qg(mat_e, 'vars', title_expr=title_e)), '\n')
    }

    # ── Factor structure + latent network ─────────────────────────────────────
    if (fam %in% c("lvm","lnm","rnm","lrnm")) {
      lam_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"lambda"), gi)'
               else         'as.matrix(getmatrix(mod,"lambda"))'
      oz_e  <- if (has_grp) 'collapse_mg(getmatrix(mod,"omega_zeta"), gi)'
               else         'as.matrix(getmatrix(mod,"omega_zeta"))'
      title_e <- if (has_grp) 'paste0("Factor Structure & Latent Network — Group ", grp_lbl)'
                 else         '"Factor Structure & Latent Network"'
      factor_body <- paste0(
        'lam   <- ', lam_e, '; if(qg_cutoff>0) lam[abs(lam)<qg_cutoff]<-0\n',
        'oz    <- tryCatch(', oz_e, ', error=function(e) matrix(0,ncol(lam),ncol(lam)))\n',
        'if(qg_cutoff>0) oz[abs(oz)<qg_cutoff]<-0\n',
        'n_lat <- ncol(lam); n_obs <- nrow(lam); n_tot <- n_lat+n_obs\n',
        'W     <- matrix(0,n_tot,n_tot)\n',
        'max_l <- if(any(lam!=0)) max(abs(lam[lam!=0])) else 0\n',
        'lam_w <- lam; if(max_l>1e-6) { nz <- lam_w!=0 & abs(lam_w)<max_l*0.2;\n',
        '  lam_w[nz] <- sign(lam_w[nz])*max_l*0.2 }\n',
        'W[1:n_lat,1:n_lat]           <- oz\n',
        'W[1:n_lat,(n_lat+1):n_tot]   <- t(lam_w)\n',
        'dir_mat <- matrix(FALSE,n_tot,n_tot)\n',
        'dir_mat[1:n_lat,(n_lat+1):n_tot] <- TRUE\n',
        'ang_l <- seq(0,2*pi,length.out=n_lat+1)[seq_len(n_lat)]\n',
        'dom   <- apply(abs(lam),1,which.max)\n',
        'obs_p <- matrix(0,n_obs,2)\n',
        'for(k in seq_len(n_lat)) {\n',
        '  idx <- which(dom==k); nk <- length(idx); if(!nk) next\n',
        '  half <- (pi/n_lat)*0.85\n',
        '  angs <- if(nk==1) ang_l[k] else seq(ang_l[k]-half,ang_l[k]+half,length.out=nk)\n',
        '  obs_p[idx,1]<-cos(angs); obs_p[idx,2]<-sin(angs)\n}\n',
        'lay <- rbind(cbind(cos(ang_l),sin(ang_l))*0.38, obs_p)\n',
        'edge_lab <- matrix("",n_tot,n_tot)\n',
        'edge_lab[1:n_lat,1:n_lat] <- ifelse(oz!=0,as.character(round(oz,2)),"")\n',
        'edge_lab[1:n_lat,(n_lat+1):n_tot] <- ifelse(t(lam)!=0,as.character(round(t(lam),2)),"")\n',
        'lat_nms <- if(!is.null(colnames(lam))&&!all(nchar(colnames(lam))==0)) colnames(lam) else latents\n',
        'obs_nms <- if(!is.null(rownames(lam))&&!all(nchar(rownames(lam))==0)) rownames(lam) else vars\n',
        'qgraph(W, labels=c(lat_nms,obs_nms),\n',
        '       shape=c(rep("circle",n_lat),rep("square",n_obs)),\n',
        '       color=c(default_pal[seq_len(n_lat)],rep("#CCCCCC",n_obs)),\n',
        '       vsize=c(rep(qg_vsize*1.8,n_lat),rep(qg_vsize,n_obs)),\n',
        '       directed=dir_mat, layout=lay, theme=qg_theme,\n',
        '       edge.labels=if(qg_elabels) edge_lab else FALSE,\n',
        '       esize=qg_esize, title=', title_e, ', mar=c(4,4,4,4))\n')
      plots_code <- paste0(plots_code,
        '# Factor structure + latent network\n',
        grp_loop(factor_body), '\n')
    }

    # ── RI Latent Network ─────────────────────────────────────────────────────
    if (needs_ri) {
      sz_e  <- if (has_grp) 'collapse_mg(getmatrix(mod,"sigma_zeta"), gi)'
               else         'as.matrix(getmatrix(mod,"sigma_zeta"))'
      title_e <- if (has_grp) 'paste0("RI Latent Network — Group ", grp_lbl)'
                 else         '"RI Latent Network"'
      ri_body <- paste0(
        'sig_z <- ', sz_e, '\n',
        'dimnames(sig_z) <- list(latents,latents)\n',
        'R_lat <- tryCatch(cov2cor(sig_z), error=function(e) NULL)\n',
        'if (!is.null(R_lat)) {\n',
        '  RI <- johnson_rw_from_cor(R_lat)\n',
        '  if(qg_cutoff>0) RI[RI<qg_cutoff]<-0\n',
        '  qgraph(RI, directed=TRUE, labels=latents, layout=qg_layout,\n',
        '         theme=qg_theme, vsize=qg_vsize, esize=qg_esize,\n',
        '         edge.labels=qg_elabels, title=', title_e, ', mar=c(4,4,4,4))\n',
        '}\n')
      plots_code <- paste0(plots_code,
        '# RI latent network\n',
        grp_loop(ri_body), '\n')
    }

    # ── Panel network plots ────────────────────────────────────────────────────
    if (fam %in% c("dlvm1","panelgvar","ri_clpm")) {
      lat_expr <- if (length(latents)>0) paste0('latents') else 'NULL'
      beta_title  <- if (has_grp) 'paste0("Temporal Network — Group ",grp_lbl)' else '"Temporal Network"'
      within_title <- if (has_grp) 'paste0("Within-person Network — Group ",grp_lbl)' else '"Within-person Network"'
      between_title <- if (has_grp) 'paste0("Between-person Network — Group ",grp_lbl)' else '"Between-person Network"'

      beta_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"beta"),gi)'
                else         'as.matrix(getmatrix(mod,"beta"))'
      wn_mat <- if (fam %in% c("dlvm1","panelgvar")) '"omega_zeta_within"' else '"omega_zeta"'
      wn_fb  <- if (fam %in% c("dlvm1","panelgvar")) '"omega_zeta"' else '"sigma_zeta"'
      wn_e   <- if (has_grp)
        paste0('tryCatch(collapse_mg(getmatrix(mod,',wn_mat,'),gi),\n           error=function(e) collapse_mg(getmatrix(mod,',wn_fb,'),gi))')
      else
        paste0('tryCatch(as.matrix(getmatrix(mod,',wn_mat,')),\n           error=function(e) as.matrix(getmatrix(mod,',wn_fb,')))')
      bn_mat <- if (fam %in% c("dlvm1","panelgvar")) '"omega_zeta_between"' else '"sigma_zeta"'
      bn_e   <- if (has_grp)
        paste0('tryCatch(collapse_mg(getmatrix(mod,',bn_mat,'),gi), error=function(e) NULL)')
      else
        paste0('tryCatch(as.matrix(getmatrix(mod,',bn_mat,')), error=function(e) NULL)')

      panel_body <- paste0(
        'qgraph(', beta_e, ', directed=TRUE, labels=', lat_expr, ',\n',
        '       layout=qg_layout, theme=qg_theme, vsize=qg_vsize, esize=qg_esize,\n',
        '       edge.labels=qg_elabels, title=', beta_title, ', mar=c(4,4,4,4))\n',
        'wn <- ', wn_e, '\n',
        'qgraph(wn, directed=FALSE, labels=', lat_expr, ',\n',
        '       layout=qg_layout, theme=qg_theme, vsize=qg_vsize, esize=qg_esize,\n',
        '       edge.labels=qg_elabels, title=', within_title, ', mar=c(4,4,4,4))\n',
        'bn <- ', bn_e, '\n',
        'if (!is.null(bn)) {\n',
        '  qgraph(bn, directed=FALSE, labels=', lat_expr, ',\n',
        '         layout=qg_layout, theme=qg_theme, vsize=qg_vsize, esize=qg_esize,\n',
        '         edge.labels=qg_elabels, title=', between_title, ', mar=c(4,4,4,4))\n',
        '}\n')
      plots_code <- paste0(plots_code,
        '# Temporal / within / between networks\n',
        grp_loop(panel_body), '\n')
    }

    # ── EGA (GGM / Ising only) ────────────────────────────────────────────────
    if (fam %in% c("ggm","ising")) {
      om_e <- if (has_grp) 'collapse_mg(getmatrix(mod,"omega"),gi)'
              else         'as.matrix(getmatrix(mod,"omega"))'
      title_e <- if (has_grp) 'paste0("EGA — Group ",grp_lbl)'
                 else         '"EGA"'
      ega_body <- paste0(
        'net_ega <- ', om_e, '\n',
        'abs_net <- abs(net_ega); diag(abs_net)<-0\n',
        'ig  <- igraph::graph_from_adjacency_matrix(abs_net,"undirected",weighted=TRUE,diag=FALSE)\n',
        'com <- igraph::cluster_walktrap(ig)  # change algorithm as needed\n',
        'mem <- igraph::membership(com)\n',
        'k   <- max(mem)\n',
        'grp_list <- setNames(lapply(seq_len(k), function(d) which(mem==d)), paste0("Dim ",seq_len(k)))\n',
        'if(qg_cutoff>0) net_ega[abs(net_ega)<qg_cutoff]<-0\n',
        'qgraph(net_ega, labels=vars, groups=grp_list, color=default_pal[seq_len(k)],\n',
        '       layout=qg_layout, theme=qg_theme, vsize=qg_vsize, esize=qg_esize,\n',
        '       edge.labels=qg_elabels, legend=TRUE, legend.cex=0.7,\n',
        '       title=', title_e, ', mar=c(4,4,4,4))\n')
      plots_code <- paste0(plots_code,
        '# EGA — community detection\n',
        grp_loop(ega_body), '\n')
    }

    # ── semPlot path diagram ──────────────────────────────────────────────────
    if (fam %in% c("lvm","lnm","rnm","lrnm")) {
      lam_s <- if (has_grp) 'collapse_mg(getmatrix(mod,"lambda"),gi)'
               else         'as.matrix(getmatrix(mod,"lambda"))'
      psi_s <- if (has_grp) 'tryCatch(collapse_mg(getmatrix(mod,"sigma_zeta"),gi),error=function(e)diag(ncol(lam_s)))'
               else         'tryCatch(as.matrix(getmatrix(mod,"sigma_zeta")),error=function(e)diag(ncol(lam_s)))'
      tht_s <- if (has_grp) 'tryCatch(collapse_mg(getmatrix(mod,"sigma_epsilon"),gi),error=function(e)diag(nrow(lam_s)))'
               else         'tryCatch(as.matrix(getmatrix(mod,"sigma_epsilon")),error=function(e)diag(nrow(lam_s)))'
      title_e <- if (has_grp) 'paste0("Path Diagram — Group ",grp_lbl)' else '"Path Diagram"'
      sem_body <- paste0(
        'if (requireNamespace("semPlot", quietly=TRUE)) {\n',
        '  library(semPlot)\n',
        '  lam_s <- ', lam_s, '\n',
        '  psi_s <- ', psi_s, '\n',
        '  tht_s <- ', tht_s, '\n',
        '  lat_nm_s <- if(!is.null(colnames(lam_s))&&!all(nchar(colnames(lam_s))==0)&&\n',
        '                 !all(grepl("^eta_?\\\\d+$",colnames(lam_s),ignore.case=TRUE)))\n',
        '               colnames(lam_s) else latents\n',
        '  obs_nm_s <- if(!is.null(rownames(lam_s))&&!all(nchar(rownames(lam_s))==0)) rownames(lam_s) else vars\n',
        '  dimnames(lam_s)<-list(obs_nm_s,lat_nm_s)\n',
        '  dimnames(psi_s)<-list(lat_nm_s,lat_nm_s)\n',
        '  dimnames(tht_s)<-list(obs_nm_s,obs_nm_s)\n',
        '  semPaths(lisrelModel(LY=lam_s,PS=psi_s,TE=tht_s),\n',
        '           what="std", whatLabels="est", layout="tree2",\n',
        '           theme="colorblind", sizeLat=10, sizeMan=5,\n',
        '           edge.label.cex=0.8, mar=c(6,1,6,1),\n',
        '           title=', title_e, ')\n',
        '}\n')
      plots_code <- paste0(plots_code,
        '# semPlot path diagram\n',
        grp_loop(sem_body), '\n')
    }

    # ── Model comparison ──────────────────────────────────────────────────────
    compare_code <- ""
    if (has_model2) {
      lbl1 <- trimws(input$model1_label); if (nchar(lbl1)==0) lbl1 <- toupper(fam)
      lbl2 <- trimws(input$model2_label); if (nchar(lbl2)==0) lbl2 <- "Model2"
      # Reconstruct model 2 fitting code
      fam2  <- input$model2_family
      ident2 <- input$model2_identification
      ge    <- if (!is.null(grp) && isTRUE(input$model2_constrain) && length(input$model2_groupEqual)>0)
                 paste0(',\n           groupEqual=', vec_str(input$model2_groupEqual))
               else ""
      mod2_fit <- if (fam2 == "ggm") {
        paste0('mod2 <- ggm(df, vars=vars, omega="full", estimator="', est, '",\n',
               '            ', ordered_str, grp_str, ')\n',
               'mod2 <- mod2 %>% setoptimizer("', opt, '") %>% runmodel()\n')
      } else {
        omega_eps2 <- if (fam2 %in% c("rnm","lrnm")) ',\n         omega_epsilon="full"' else ""
        paste0('mod2 <- ', fam2, '(df, lambda=Lambda, vars=vars, latents=latents,\n',
               '         identification="', ident2, '",\n',
               '         estimator="', est, '", ', ordered_str, grp_str, omega_eps2, ge, ')\n',
               'mod2 <- mod2 %>% setoptimizer("', opt, '") %>% runmodel()\n')
      }
      compare_code <- paste0(
        '\n# ── Model comparison ─────────────────────────────────────────────────\n',
        '# Model 2: ', toupper(fam2), if (nchar(ge)>0) paste0(' (constrained: ', paste(input$model2_groupEqual,collapse=", "), ')') else "", '\n',
        mod2_fit,
        '.pn_compare <- getFromNamespace("compare","psychonetrics")\n',
        '.pn_compare(', q(lbl1), '=mod, ', q(lbl2), '=mod2)\n',
        'mod  %>% fit   # Model 1\n',
        'mod2 %>% fit   # Model 2\n')
    }

    # ── Bootstrap ─────────────────────────────────────────────────────────────
    boot_code <- ""
    if (has_boot) {
      boot_inner <- if (fam == "ggm") {
        paste0('ggm(df, vars=vars, estimator="', est, '",\n',
               '        ordered=', ordered_str, ', bootstrap="nonparametric") %>%\n',
               '      runmodel %>%\n',
               '      prune(alpha=', input$prune_alpha, ')')
      } else {
        paste0(fam, '(df, lambda=Lambda, vars=vars, latents=latents,\n',
               '          identification="', input$identification, '",\n',
               '          estimator="', est, '", bootstrap="nonparametric") %>%\n',
               '      runmodel %>%\n',
               '      prune(alpha=', input$prune_alpha, ')')
      }
      boot_code <- paste0(
        '\n# ── Bootstrapping ────────────────────────────────────────────────────\n',
        'bootstraps <- loop_psychonetrics({\n',
        '  ', boot_inner, '\n',
        '}, reps=200, nCores=2)  # adjust reps / nCores\n',
        'boot_agg <- aggregate_bootstraps(sample=mod, bootstraps=bootstraps)\n',
        'CIplot(boot_agg, "', if (fam %in% c("lnm","lrnm")) "omega_zeta"
                               else if (fam %in% c("rnm")) "omega_epsilon"
                               else "omega", '", split0=TRUE)\n')
    }

    # ── Assemble final script ─────────────────────────────────────────────────
    lib_code <- paste0(
      '# ── Psychonetrics Analysis ───────────────────────────────────────────\n',
      '# Generated by PsychoNetrix  |  ', Sys.Date(), '\n',
      '# Model: ', toupper(fam),
      if (has_grp) paste0('  |  Multi-group: "', grp, '" (', paste(grp_levels_vec, collapse=", "), ')') else "",
      '\n\n',
      'library(psychonetrics)\n',
      'library(qgraph)\n',
      'library(dplyr)\n',
      if (needs_ri || fam %in% c("ggm","ising")) 'library(igraph)\n' else "",
      '\n',
      '# Load your data:\n',
      '# df <- read.csv("your_data.csv")\n\n',
      'vars <- ', vec_str(vars), '\n\n',
      if (has_grp) paste0('# Group variable: rows are split by df[["', grp, '"]]\n') else "",
      transform_code,
      collapse_mg_fn,
      rw_fn)

    paste0(lib_code,
           qg_opts_code,
           plot_helpers,
           fit_code,
           prune_code,
           results_code,
           extract_code,
           fscores_code,
           plots_code,
           compare_code,
           boot_code)
  }

  generated_code <- reactive({ rv$r_code })

  output$r_code_output <- renderPrint({
    cat(generated_code())
  })

  output$download_code <- downloadHandler(
    filename = function() paste0("psychonetrics_analysis_", Sys.Date(), ".R"),
    content  = function(file) writeLines(generated_code(), file)
  )

} # end server


# ── 3. Run ───────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
