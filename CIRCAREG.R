# ==============================================================================
# Functional Regression Suite (Shiny App) - VERSION 2 FULLY FIXED
# ==============================================================================
# Features:
# 1. Data Import (Wide/Long) & Sample Generator
# 2. Smoothing (B-Splines, Fourier)
# 3. Function-on-Scalar Regression (FoSR)
#    - Pointwise OLS with Matrix Algebra & Residual Bootstrap
#    - Smoothed OLS using mgcv::gam (Penalized Splines)
# 4. Scalar-on-Function Regression (SoFR) using refund::pfr
# 5. Interactive Visualizations (Plotly)
# ==============================================================================

# --- Install & Load Packages ---
packages <- c("shiny", "shinydashboard", "shinyWidgets", "fda", "mgcv", 
              "plotly", "DT", "dplyr", "tidyr", "ggplot2", "refund", "viridis")

# Optional but recommended for robust exponential saturation fitting
optional_packages <- c("minpack.lm")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(packages, install_if_missing))

# Try to load optional packages (don't fail if not available)
for(pkg in optional_packages) {
  try(library(pkg, character.only = TRUE), silent = TRUE)
}

# ==============================================================================
# UI DEFINITION
# ==============================================================================
ui <- dashboardPage(
  dashboardHeader(title = "Functional Regression Suite"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Import", tabName = "import", icon = icon("upload")),
      menuItem("Smoothing", tabName = "preprocess", icon = icon("wave-square")),
      menuItem("Smoothing Diagnostics", tabName = "smooth_diag", icon = icon("chart-area"),
               badgeLabel = "New", badgeColor = "green"),
      menuItem("Function-on-Scalar (FoSR)", tabName = "fosr", icon = icon("chart-line")),
      menuItem("Scalar-on-Function (SoFR)", tabName = "sofr", icon = icon("chart-bar")),
      menuItem("Harmonic Regression", tabName = "harmonic", icon = icon("sync")),
      menuItem("Pairwise Comparisons", tabName = "pairwise", icon = icon("not-equal"),
               badgeLabel = "New", badgeColor = "green"),
      menuItem("Export", tabName = "export", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f4f4f4; }
        .reg-control-panel { background-color: #fff; padding: 15px; border: 1px solid #ddd; border-radius: 5px; margin-bottom: 15px; }
        .box-header { font-weight: bold; }
      "))
    ),
    
    tabItems(
      # --- TAB 1: Data Import ---
      tabItem(
        tabName = "import",
        fluidRow(
          box(
            title = "1. Upload Data", status = "primary", solidHeader = TRUE, width = 4,
            fileInput("datafile", "Choose File", accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx")),
            hr(),
            checkboxInput("header", "File has header row", TRUE),
            conditionalPanel(
              condition = "input.datafile && !input.datafile.name.match(/\\.(xls|xlsx)$/i)",
              radioButtons("sep", "Separator:", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ",", inline = TRUE)
            ),
            conditionalPanel(
              condition = "input.datafile && input.datafile.name.match(/\\.(xls|xlsx)$/i)",
              uiOutput("excel_sheet_selector")
            ),
            hr(),
            radioButtons("data_format", "Data Format:",
                         choices = list("Wide (subjects in rows)" = "wide", "Long (subjects in columns)" = "long"),
                         selected = "wide"),
            helpText("Wide: Rows=Subjects, Cols=Time Points."),
            hr(),
            actionButton("load_data", "Load Raw File", class = "btn-warning"),
            hr(),
            actionButton("generate_sample", "Generate Sample Data", class = "btn-info")
          ),
          box(
            title = "2. Variable Selection", status = "warning", solidHeader = TRUE, width = 8,
            uiOutput("var_select_container"),
            hr(),
            h4("Data Status:"),
            verbatimTextOutput("data_status")
          )
        ),
        fluidRow(
          box(
            title = "Data Preview", status = "primary", width = 12,
            DTOutput("data_preview"),
            hr(),
            plotOutput("raw_data_plot", height = "250px")
          )
        )
      ),
      
      # --- TAB 2: Smoothing ---
      tabItem(
        tabName = "preprocess",
        fluidRow(
          box(
            title = "Smoothing Options", status = "warning", solidHeader = TRUE, width = 4,
            radioButtons("smooth_method", "Smoothing Method:",
                         choices = list("Raw data (No Smoothing)" = "none",
                                        "Automatic smoothing (REML)" = "auto",
                                        "Manual smoothing" = "manual"),
                         selected = "none"),
            checkboxInput("is_cyclic", "Is data cyclic? (e.g. 24h clock)", FALSE),
            checkboxInput("constrain_bounds", "Constrain smoothed values to specific range", FALSE),
            conditionalPanel(
              condition = "input.constrain_bounds == true",
              numericInput("min_bound", "Minimum value:", value = 0, step = 1),
              numericInput("max_bound", "Maximum value:", value = 100, step = 1),
              helpText("Smoothed values will be clamped to [min, max] range.")
            ),
            conditionalPanel(
              condition = "input.smooth_method == 'manual'",
              sliderInput("smooth_factor", "Smoothing Factor (Lambda):", min = 0.1, max = 10, value = 1, step = 0.1),
              conditionalPanel(
                condition = "output.diagnostics_available",
                hr(),
                actionButton("use_diagnostic_lambda", "📊 Use Diagnostic Results",
                             class = "btn-info btn-sm"),
                helpText("Click to automatically set smoothing factor based on REML/CV analysis")
              )
            ),
            conditionalPanel(
              condition = "input.smooth_method == 'auto' || input.smooth_method == 'manual'",
              numericInput("n_basis", "Number of B-spline basis functions:", value = 12, min = 4, max = 100),
              helpText("Higher basis count captures more detail.")
            ),
            actionButton("apply_smooth", "Apply / Update Data", class = "btn-warning"),
            br(), br(),
            uiOutput("smooth_fit_display")
          ),
          box(
            title = "Functional Data Visualization", status = "info", solidHeader = TRUE, width = 8,
            plotlyOutput("data_plot", height = "500px")
          )
        )
      ),

      # --- TAB 3: Smoothing Diagnostics ---
      tabItem(
        tabName = "smooth_diag",
        fluidRow(
          box(
            title = "Smoothing Parameter Exploration",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            helpText("Explore optimal smoothing parameters using REML, cross-validation, and effective degrees of freedom (EDF).")
          )
        ),

        fluidRow(
          box(
            title = "1. GAM REML Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            helpText("Fit a REML GAM to a representative subject/curve and inspect effective degrees of freedom (EDF)."),
            selectInput("diag_subject_type", "Select data to analyze:",
                        choices = list("Mean curve across all subjects" = "mean",
                                       "Specific subject" = "subject"),
                        selected = "mean"),
            conditionalPanel(
              condition = "input.diag_subject_type == 'subject'",
              numericInput("diag_subject_id", "Subject ID:", value = 1, min = 1, step = 1)
            ),
            actionButton("run_gam_reml", "Fit GAM (REML)", class = "btn-primary", icon = icon("play")),
            hr(),
            h4("GAM REML Results:"),
            verbatimTextOutput("gam_reml_summary"),
            helpText(HTML("<b>Interpretation:</b><br>
                           - <b>EDF (Effective Degrees of Freedom):</b> Measures model complexity. Higher EDF = more flexible fit.<br>
                           - <b>EDF ≈ 1-2:</b> Very smooth, may be underfit<br>
                           - <b>EDF ≈ 3-6:</b> Moderate smoothing, often optimal<br>
                           - <b>EDF ≈ k-1:</b> Minimal smoothing, may overfit<br>
                           - <b>REML score:</b> Lower is better (balance fit vs smoothness)"))
          ),

          box(
            title = "2. Smoothing Parameter Selection",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            helpText("Configure lambda range for REML profile and cross-validation analysis."),
            sliderInput("lambda_min_exp", "Log10(Lambda Min):",
                        min = -8, max = 0, value = -6, step = 0.5),
            sliderInput("lambda_max_exp", "Log10(Lambda Max):",
                        min = -2, max = 4, value = 2, step = 0.5),
            numericInput("n_lambda", "Number of lambda values:",
                         value = 20, min = 5, max = 50, step = 1),
            hr(),
            numericInput("cv_k_folds", "K-folds for cross-validation:",
                         value = 5, min = 3, max = 10, step = 1),
            hr(),
            actionButton("run_reml_profile", "Compute REML Profile",
                         class = "btn-warning", icon = icon("chart-line")),
            br(), br(),
            actionButton("run_cv_analysis", "Run Cross-Validation",
                         class = "btn-success", icon = icon("random"))
          )
        ),

        fluidRow(
          box(
            title = "3. REML Profile",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("reml_profile_plot", height = "400px"),
            helpText(HTML("<b>Interpretation:</b><br>
                           - <b>REML score:</b> Restricted Maximum Likelihood criterion<br>
                           - <b>Optimal lambda:</b> Minimum of the REML curve<br>
                           - <b>Flat region:</b> Multiple good lambda values (robust)<br>
                           - <b>Sharp minimum:</b> Sensitive to lambda choice<br>
                           - <b>Compare with auto REML (lambda=0):</b> Should be near minimum"))
          ),

          box(
            title = "4. Cross-Validation Curve",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("cv_curve_plot", height = "400px"),
            helpText(HTML("<b>Interpretation:</b><br>
                           - <b>CV Error:</b> Average prediction error on held-out data<br>
                           - <b>Optimal lambda:</b> Minimum CV error<br>
                           - <b>1-SE rule:</b> Most parsimonious model within 1 SE of minimum<br>
                           - <b>U-shape:</b> Underfitting (left) vs overfitting (right)<br>
                           - <b>Compare REML vs CV:</b> Should agree on optimal range"))
          )
        ),

        fluidRow(
          box(
            title = "5. Comparison Summary",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("smoothing_comparison_summary"),
            helpText(HTML("<b>Decision Guide:</b><br>
                           - <b>If REML and CV agree:</b> Use their recommended lambda<br>
                           - <b>If REML and CV disagree:</b> Prefer CV for prediction, REML for smoothness<br>
                           - <b>For automatic smoothing:</b> Lambda = 0 uses REML optimization<br>
                           - <b>For manual smoothing:</b> See 'Smoothing Factor' values above, or use the button below<br>
                           - <b>Trade-off:</b> Smooth curves (high lambda) vs capturing variability (low lambda)<br>
                           <br>
                           <b>📊 How to apply these results:</b><br>
                           1. Go to 'Smoothing' tab<br>
                           2. Select 'Manual smoothing'<br>
                           3. Click '📊 Use Diagnostic Results' button<br>
                           4. Or manually set the Smoothing Factor value shown above<br>
                           5. Click 'Apply / Update Data' to process your data"))
          )
        )
      ),

      # --- TAB 4: Function-on-Scalar Regression (FoSR) ---
      tabItem(
        tabName = "fosr",
        fluidRow(
          box(
            title = "FoSR Settings (Response: Curve)", status = "primary", solidHeader = TRUE, width = 4,
            uiOutput("fosr_var_select_ui"),
            selectInput("reg_method", "Method:", 
                        choices = c("Pointwise OLS (Raw)" = "OLS_nosmooth",
                                    "Smoothed OLS (GAM/Penalized Splines)" = "GAM_smooth"),
                        selected = "OLS_nosmooth"),
            helpText("Smoothed OLS uses mgcv::gam with penalized splines."),
            
            # Bootstrap Options (Visible only for OLS)
            conditionalPanel(
              condition = "input.reg_method == 'OLS_nosmooth'",
              checkboxInput("use_bootstrap", "Compute Bootstrap CIs (Scalable)", FALSE),
              conditionalPanel(
                condition = "input.use_bootstrap == true",
                numericInput("n_boot", "Number of Simulations (B):", value = 200, min = 10, max = 1000)
              )
            ),
            hr(),
            actionButton("run_fosr", "Run FoSR Model", class = "btn-success", icon = icon("play"))
          ),
          box(
            title = "Model Summary", status = "info", solidHeader = TRUE, width = 8,
            verbatimTextOutput("fosr_model_summary")
          )
        ),
        fluidRow(
          tabBox(
            title = "Diagnostics & Visualization", id = "fosr_tabs", width = 12,
            tabPanel("1. Observed Data", icon = icon("eye"),
                     fluidRow(
                       column(4, selectInput("reg_color_var", "Color curves by:", choices = NULL)),
                       column(8, plotlyOutput("reg_observed_plot", height = "600px"))
                     )
            ),
            tabPanel("2. Interactive Predictions", icon = icon("sliders-h"),
                     fluidRow(
                       column(4, 
                              div(class = "reg-control-panel",
                                  h4("Predictor Values"),
                                  helpText("Adjust sliders to see the predicted curve."),
                                  uiOutput("reg_prediction_controls"),
                                  hr(),
                                  h4("Reference Lines"),
                                  uiOutput("reg_ref_selector")
                              )
                       ),
                       column(8, plotlyOutput("reg_fitted_plot", height = "600px"))
                     )
            ),
            tabPanel("3. Coefficient Functions", icon = icon("wave-square"),
                     fluidRow(
                       column(4, selectInput("reg_coeff_select", "Select Coefficient:", choices = NULL)),
                       column(8, 
                              plotlyOutput("reg_coeff_plot", height = "400px"),
                              hr(),
                              h4("Significance (P-values)"),
                              plotlyOutput("reg_pvalue_plot", height = "300px")
                       )
                     )
            ),
            tabPanel("4. Model Fit (R2)", icon = icon("chart-line"), plotlyOutput("reg_r2_plot", height = "600px")),
            tabPanel("5. Residuals", icon = icon("chart-area"), plotlyOutput("reg_residual_plot", height = "600px"))
          )
        )
      ),
      
      # --- TAB 4: Scalar-on-Function Regression (SoFR) ---
      tabItem(
        tabName = "sofr",
        fluidRow(
          box(
            title = "SoFR Settings (Response: Scalar)", status = "warning", solidHeader = TRUE, width = 4,
            uiOutput("sofr_var_select_ui"),
            selectInput("sofr_method", "Method:", choices = c("Penalized Functional Regression (pfr)" = "pfr")),
            hr(),
            selectInput("sofr_family", "Response Distribution (Family):",
                        choices = c("Auto-detect" = "auto",
                                    "Gaussian (continuous)" = "gaussian",
                                    "Binomial (binary/proportion)" = "binomial",
                                    "Poisson (counts)" = "poisson",
                                    "Gamma (positive continuous)" = "Gamma"),
                        selected = "auto"),
            conditionalPanel(
              condition = "input.sofr_family == 'binomial'",
              selectInput("sofr_link_binomial", "Link Function:",
                          choices = c("logit" = "logit", "probit" = "probit", "cloglog" = "cloglog"),
                          selected = "logit")
            ),
            conditionalPanel(
              condition = "input.sofr_family == 'gaussian'",
              selectInput("sofr_link_gaussian", "Link Function:",
                          choices = c("identity" = "identity", "log" = "log"),
                          selected = "identity")
            ),
            uiOutput("sofr_family_info"),
            hr(),
            checkboxInput("sofr_use_bootstrap", "Compute Bootstrap CIs for β(t)", FALSE),
            conditionalPanel(
              condition = "input.sofr_use_bootstrap == true",
              numericInput("sofr_n_boot", "Number of Bootstrap Samples (B):", value = 100, min = 20, max = 500),
              helpText("Note: Bootstrap for SoFR can be slow. Start with B=50-100.")
            ),
            hr(),
            helpText("Model: g(E(y)) = alpha + Integral(X(t)*Beta(t)) + Z*Gamma"),
            helpText("where g() is the link function."),
            hr(),
            actionButton("run_sofr", "Run SoFR Model", class = "btn-warning", icon = icon("play"))
          ),
          box(
            title = "Inference & Summary", status = "info", solidHeader = TRUE, width = 8,
            verbatimTextOutput("sofr_inference_summary"),
            hr(),
            uiOutput("sofr_model_diagnostics")
          )
        ),
        fluidRow(
          tabBox(
            title = "SoFR Visualizations", id = "sofr_tabs", width = 12,
            tabPanel("1. Estimated Coefficient Beta(t)", icon = icon("wave-square"),
                     plotlyOutput("sofr_coeff_plot", height = "500px"),
                     uiOutput("sofr_coeff_interpretation")
            ),
            tabPanel("2. Observed vs Predicted", icon = icon("chart-line"),
                     plotlyOutput("sofr_pred_plot", height = "500px")
            ),
            tabPanel("3. Classification Diagnostics", icon = icon("bullseye"),
                     conditionalPanel(
                       condition = "output.sofr_is_binary",
                       fluidRow(
                         column(6, plotlyOutput("sofr_roc_plot", height = "400px")),
                         column(6, plotlyOutput("sofr_calibration_plot", height = "400px"))
                       ),
                       hr(),
                       verbatimTextOutput("sofr_classification_metrics")
                     ),
                     conditionalPanel(
                       condition = "!output.sofr_is_binary",
                       helpText("Classification diagnostics are only available for binary outcomes.")
                     )
            )
          )
        )
      ),
      
      # --- TAB 5: Harmonic Regression (Cosinor Analysis) ---
      tabItem(
        tabName = "harmonic",
        fluidRow(
          box(
            title = "Harmonic Regression Settings", status = "success", solidHeader = TRUE, width = 4,
            h4("Data Selection"),
            uiOutput("harmonic_var_select_ui"),
            hr(),
            h4("Model Specification"),
            numericInput("harmonic_period", "Fundamental Period (τ):", value = 24, min = 1, max = 168, step = 1),
            helpText("Period in the same units as your time variable (e.g., 24 for circadian)."),
            sliderInput("n_harmonics", "Number of Harmonics:", min = 1, max = 8, value = 1, step = 1),
            helpText("1 = fundamental (τ), 2 = adds τ/2, 3 = adds τ/3, etc."),
            uiOutput("harmonic_warning_ui"),
            selectInput("harmonic_trend_type", "Homeostatic Trend Model:",
                        choices = c("None (circadian only)" = "none",
                                    "Linear" = "linear",
                                    "Logarithmic" = "log",
                                    "Saturating exponential" = "exp_sat"),
                        selected = "none"),
            conditionalPanel(
              condition = "input.harmonic_trend_type != 'none'",
              helpText(HTML("<b>Two-Process Model:</b> Separates homeostatic sleep pressure (trend) from circadian modulation (sinusoidal)."))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'linear'",
              helpText(HTML("Y(t) = M + β·t + harmonics<br><small>Constant rate of change</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'log'",
              helpText(HTML("Y(t) = M + β·log(t+1) + harmonics<br><small>Rapid initial change, slowing over time</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'exp_sat'",
              helpText(HTML("Y(t) = M + A·(1-e<sup>-t/τ</sup>) + harmonics<br>
                            <small>Asymptotic saturation (classic Process S)</small><br>
                            <small><em>⚠️ Uses nonlinear least squares. If fit is poor, try Log trend instead.</em></small>"))
            ),
            hr(),
            h4("Group Analysis"),
            uiOutput("harmonic_group_var_ui"),
            helpText("Optional: Select a group variable to compare rhythms between groups."),
            hr(),
            h4("Bootstrap Options"),
            checkboxInput("harmonic_bootstrap", "Compute Bootstrap CIs", FALSE),
            conditionalPanel(
              condition = "input.harmonic_bootstrap == true",
              numericInput("harmonic_n_boot", "Bootstrap Iterations (B):", value = 500, min = 100, max = 2000)
            ),
            hr(),
            h4("Parameter Constraints"),
            checkboxInput("harmonic_use_bounds", "Enable Parameter Bounding", FALSE),
            conditionalPanel(
              condition = "input.harmonic_use_bounds == true",
              helpText("Constrain parameters to stay within plausible ranges based on your data scale."),
              uiOutput("harmonic_bounds_hints"),
              hr(),
              h5("Common Parameters"),
              numericInput("harmonic_mesor_min", "MESOR Min:", value = NA, step = 0.1),
              numericInput("harmonic_mesor_max", "MESOR Max:", value = NA, step = 0.1),
              numericInput("harmonic_amplitude_min", "Amplitude Min:", value = 0, step = 0.1),
              numericInput("harmonic_amplitude_max", "Amplitude Max:", value = NA, step = 0.1),
              conditionalPanel(
                condition = "input.harmonic_trend_type == 'exp_sat'",
                hr(),
                h5("Exponential Saturation Parameters"),
                numericInput("harmonic_A_sat_min", "A_sat (Asymptote) Min:", value = NA, step = 0.1),
                numericInput("harmonic_A_sat_max", "A_sat (Asymptote) Max:", value = NA, step = 0.1),
                numericInput("harmonic_tau_min", "τ (Time Constant) Min:", value = 0.5, step = 0.1),
                numericInput("harmonic_tau_max", "τ (Time Constant) Max:", value = NA, step = 1),
                helpText(HTML("<small>τ controls how fast saturation is reached. Smaller τ = faster saturation.</small>"))
              ),
              helpText(HTML("<small>Leave blank (NA) for no bound. Defaults based on data range shown in hints above.</small>"))
            ),
            hr(),
            actionButton("run_harmonic", "Run Harmonic Regression", class = "btn-success", icon = icon("play"))
          ),
          box(
            title = "Model Summary", status = "info", solidHeader = TRUE, width = 8,
            verbatimTextOutput("harmonic_summary"),
            hr(),
            uiOutput("harmonic_parameters_table")
          )
        ),
        fluidRow(
          tabBox(
            title = "Harmonic Regression Results", id = "harmonic_tabs", width = 12,
            tabPanel("1. Fitted Curves", icon = icon("chart-line"),
                     fluidRow(
                       column(8, plotlyOutput("harmonic_fit_plot", height = "500px")),
                       column(4, 
                              h4("Display Options"),
                              uiOutput("harmonic_subject_selector"),
                              checkboxInput("harmonic_show_ci", "Show Confidence Bands", TRUE),
                              checkboxInput("harmonic_show_data", "Show Raw Data Points", TRUE),
                              checkboxInput("harmonic_show_components", "Show Harmonic Components", FALSE)
                       )
                     )
            ),
            tabPanel("2. Polar Plot (Acrophase)", icon = icon("compass"),
                     fluidRow(
                       column(8, plotlyOutput("harmonic_polar_plot", height = "500px")),
                       column(4,
                              h4("Polar Plot Settings"),
                              uiOutput("harmonic_selector_polar"),
                              helpText("Acrophase displayed in polar coordinates. Radius = Amplitude, Angle = Acrophase."),
                              checkboxInput("polar_show_mean", "Show Population Mean Vector", TRUE),
                              checkboxInput("polar_show_ellipse", "Show Confidence Ellipse", TRUE)
                       )
                     )
            ),
            tabPanel("3. Parameter Distribution", icon = icon("chart-bar"),
                     fluidRow(
                       column(12, 
                              uiOutput("harmonic_selector_dist"),
                              hr()
                       )
                     ),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_amplitude_hist", height = "400px")),
                       column(6, plotlyOutput("harmonic_acrophase_hist", height = "400px"))
                     ),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_mesor_plot", height = "300px")),
                       column(6, plotlyOutput("harmonic_trend_hist", height = "300px"))
                     )
            ),
            tabPanel("4. Individual Results", icon = icon("table"),
                     fluidRow(
                       column(12,
                              h4("Individual Subject Parameters"),
                              helpText("R-squared and other parameters for each subject. Use the export button to download."),
                              DT::dataTableOutput("harmonic_individual_table"),
                              br(),
                              downloadButton("export_harmonic_individual", "Export Individual Parameters (CSV)", class = "btn-success")
                       )
                     )
            ),
            tabPanel("5. Residual Diagnostics", icon = icon("microscope"),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_residual_plot", height = "400px")),
                       column(6, plotlyOutput("harmonic_qq_plot", height = "400px"))
                     ),
                     fluidRow(
                       column(12, verbatimTextOutput("harmonic_gof_stats"))
                     )
            ),
            tabPanel("6. Group Comparison", icon = icon("users"),
                     fluidRow(
                       column(12, 
                              uiOutput("harmonic_selector_group"),
                              hr()
                       )
                     ),
                     fluidRow(
                       column(12, plotlyOutput("harmonic_group_comparison_plot", height = "500px"))
                     ),
                     hr(),
                     verbatimTextOutput("harmonic_group_test_results")
            )
          )
        )
      ),

      # --- TAB 6: Pairwise Comparisons ---
      tabItem(
        tabName = "pairwise",
        fluidRow(
          box(
            title = "Pairwise Group Comparisons", status = "primary", solidHeader = TRUE, width = 12,
            helpText("Compare circadian parameters between groups with multiple comparison corrections."),
            helpText("Note: This analysis requires 2 or more groups defined in the Harmonic Regression tab.")
          )
        ),
        fluidRow(
          box(
            title = "Settings", status = "info", solidHeader = TRUE, width = 4,
            selectInput("pairwise_param", "Parameter to Compare:",
                       choices = c("MESOR" = "mesor",
                                   "H1 Amplitude" = "amplitude_1",
                                   "H1 Acrophase (hours)" = "acrophase_time_1",
                                   "H2 Amplitude" = "amplitude_2",
                                   "H2 Acrophase (hours)" = "acrophase_time_2",
                                   "H3 Amplitude" = "amplitude_3",
                                   "H3 Acrophase (hours)" = "acrophase_time_3",
                                   "R-squared" = "r_squared",
                                   "A_sat" = "A_sat",
                                   "τ (tau)" = "tau",
                                   "Process S (%)" = "percent_S",
                                   "Process C (%)" = "percent_C"),
                       selected = "amplitude_1"),
            selectInput("pairwise_correction", "Multiple Comparison Correction:",
                       choices = c("None" = "none",
                                   "Bonferroni" = "bonferroni",
                                   "Holm" = "holm",
                                   "Hochberg" = "hochberg",
                                   "Hommel" = "hommel",
                                   "Benjamini-Hochberg (FDR)" = "BH",
                                   "Benjamini-Yekutieli (FDR)" = "BY"),
                       selected = "holm"),
            hr(),
            checkboxInput("pairwise_show_effect_size", "Show Effect Sizes (Cohen's d)", value = TRUE),
            checkboxInput("pairwise_show_ci", "Show 95% Confidence Intervals", value = TRUE),
            hr(),
            actionButton("run_pairwise", "Run Pairwise Comparisons", class = "btn-warning", icon = icon("play"))
          ),
          box(
            title = "Pairwise Test Results", status = "success", solidHeader = TRUE, width = 8,
            verbatimTextOutput("pairwise_results")
          )
        ),
        fluidRow(
          box(
            title = "Visualization", status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("pairwise_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Comparison Matrix", status = "info", solidHeader = TRUE, width = 6,
            uiOutput("pairwise_matrix_help"),
            tableOutput("pairwise_matrix")
          ),
          box(
            title = "Download", status = "warning", solidHeader = TRUE, width = 6,
            downloadButton("export_pairwise_results", "Download Pairwise Results (CSV)", class = "btn-primary"),
            br(), br(),
            downloadButton("export_pairwise_plot", "Download Plot (PNG)", class = "btn-success")
          )
        )
      ),

      # --- TAB 7: Export ---
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "Export Results", status = "info", solidHeader = TRUE, width = 12,
            h4("FoSR Exports"),
            downloadButton("export_scores_csv", "Download FoSR Coefficients (CSV)", class = "btn-primary"),
            br(), br(),
            downloadButton("export_r_code", "Download Reproduction R Code", class = "btn-warning"),
            hr(),
            h4("Harmonic Regression Exports"),
            downloadButton("export_harmonic_params", "Download Individual Parameters (CSV)", class = "btn-success"),
            br(), br(),
            downloadButton("export_harmonic_summary", "Download Summary Report (TXT)", class = "btn-info")
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {
  
  # --- Reactive State ---
  values <- reactiveValues(
    raw_df = NULL,
    data = NULL,        # Functional Matrix (Y for FoSR, X for SoFR)
    covariates = NULL,  # Scalar DataFrame
    time_labels = NULL, # Original time column names
    smooth_data = NULL, # Smoothed Matrix
    smooth_fit_metrics = NULL, # Smoothing fit metrics (R², RMSE)
    fd_obj = NULL,      # fda object
    reg_model = NULL,   # FoSR Model Result
    sofr_model = NULL,  # SoFR Model Result
    harmonic_model = NULL,  # Harmonic Regression Result
    # Smoothing diagnostics
    gam_reml_result = NULL,
    reml_profile_result = NULL,
    cv_result = NULL,
    diagnostic_lambda = NULL  # Stores recommended lambda from diagnostics
  )
  
  # ==============================================================================
  # HELPER FUNCTION: Build GAM prediction dataframe
  # ==============================================================================
  build_gam_pred_df <- function(time_vec, predictors, long_data, factor_levels, input_vals) {
    # Create dataframe with time column
    pred_df <- data.frame(time = time_vec)
    n_t <- length(time_vec)
    
    for(var in predictors) {
      val <- input_vals[[var]]
      if(is.null(val)) return(NULL)
      
      # Get the column from long_data to check type and levels
      orig_col <- long_data[[var]]
      
      if(is.numeric(orig_col)) {
        # Numeric variable: replicate scalar value
        pred_df[[var]] <- rep(as.numeric(val), n_t)
      } else if(is.factor(orig_col)) {
        # Factor: use EXACT levels from the fitted model's data
        lvls <- levels(orig_col)
        # Validate that val is in levels
        if(!(val %in% lvls)) {
          warning(paste("Value", val, "not in factor levels for", var))
          val <- lvls[1]  # Fall back to first level
        }
        pred_df[[var]] <- factor(rep(val, n_t), levels = lvls)
      } else {
        # Character or other: use stored factor levels if available
        if(!is.null(factor_levels) && !is.null(factor_levels[[var]])) {
          lvls <- factor_levels[[var]]
        } else {
          lvls <- sort(unique(as.character(orig_col)))
        }
        if(!(val %in% lvls)) {
          warning(paste("Value", val, "not in levels for", var))
          val <- lvls[1]
        }
        pred_df[[var]] <- factor(rep(val, n_t), levels = lvls)
      }
    }
    return(pred_df)
  }
  
  # ==============================================================================
  # HELPER FUNCTION: Safe GAM prediction with SE
  # ==============================================================================
  safe_gam_predict <- function(gam_obj, newdata) {
    result <- tryCatch({
      pred <- predict(gam_obj, newdata = newdata, type = "response", se.fit = TRUE)
      list(
        fit = as.numeric(pred$fit),
        se = as.numeric(pred$se.fit),
        success = TRUE
      )
    }, error = function(e) {
      # Fallback: try without SE
      pred <- tryCatch({
        predict(gam_obj, newdata = newdata, type = "response")
      }, error = function(e2) {
        return(NULL)
      })
      if(is.null(pred)) {
        return(list(fit = NULL, se = NULL, success = FALSE))
      }
      list(
        fit = as.numeric(pred),
        se = rep(0, length(pred)),
        success = TRUE
      )
    })
    return(result)
  }
  
  # --- Excel Sheet Selector ---
  output$excel_sheet_selector <- renderUI({
    req(input$datafile)
    file_ext <- tools::file_ext(input$datafile$name)

    if(tolower(file_ext) %in% c("xls", "xlsx")) {
      tryCatch({
        # Check if readxl is available
        if(!requireNamespace("readxl", quietly = TRUE)) {
          return(div(style = "color: red;",
                     "readxl package required for Excel files. Install with: install.packages('readxl')"))
        }

        sheet_names <- readxl::excel_sheets(input$datafile$datapath)
        selectInput("excel_sheet", "Select Sheet:",
                    choices = sheet_names,
                    selected = sheet_names[1])
      }, error = function(e) {
        div(style = "color: red;", paste("Error reading Excel file:", e$message))
      })
    }
  })

  # --- 1. Data Loading ---
  observeEvent(input$load_data, {
    if(is.null(input$datafile)) return()
    tryCatch({
      file_ext <- tools::file_ext(input$datafile$name)

      # Read Excel files
      if(tolower(file_ext) %in% c("xls", "xlsx")) {
        if(!requireNamespace("readxl", quietly = TRUE)) {
          showNotification("readxl package required. Install with: install.packages('readxl')", type = "error")
          return()
        }

        # Use selected sheet or default to first sheet
        sheet_to_read <- if(!is.null(input$excel_sheet)) input$excel_sheet else 1
        raw <- suppressWarnings(as.data.frame(readxl::read_excel(input$datafile$datapath,
                                                 sheet = sheet_to_read,
                                                 col_names = input$header,
                                                 guess_max = 10000)))
      } else {
        # Read CSV/text files
        raw <- read.csv(input$datafile$datapath, header = input$header, sep = input$sep, check.names = FALSE)
      }

      if(any(duplicated(colnames(raw)))) colnames(raw) <- make.unique(colnames(raw))
      if(input$data_format == "long") raw <- as.data.frame(t(raw))
      values$raw_df <- raw
      values$data <- NULL; values$covariates <- NULL
      showNotification(paste("File loaded successfully:", input$datafile$name), type = "message")
    }, error = function(e) showNotification(paste("Load Error:", e$message), type = "error"))
  })
  
  observeEvent(input$generate_sample, {
    set.seed(123)
    n_subjects <- 50; n_time <- 24
    t <- seq(0, 1, length.out = n_time)
    age <- round(rnorm(n_subjects, 50, 10))
    sex <- factor(sample(c("Male", "Female"), n_subjects, replace = TRUE), levels = c("Female", "Male"))
    score <- rnorm(n_subjects, 100, 15) 
    
    mat <- matrix(NA, n_subjects, n_time)
    for(i in 1:n_subjects) {
      y <- 10 * sin(2*pi*t) + (age[i]-50)*0.1 
      if(sex[i] == "Female") y <- y * 1.2
      mat[i,] <- y + rnorm(n_time, 0, 0.5) + 5
    }
    
    # Generate binary outcome based on functional integral (for SoFR testing)
    # Higher curve values -> higher probability of outcome = 1
    curve_integral <- rowMeans(mat)  # Simple summary of functional data
    prob_outcome <- plogis(-2 + 0.3 * curve_integral + 0.02 * age)  # Logistic model
    binary_outcome <- rbinom(n_subjects, 1, prob_outcome)
    
    values$data <- mat
    values$covariates <- data.frame(
      ID = 1:n_subjects, 
      Age = age, 
      Sex = sex, 
      Score = score,
      Outcome = binary_outcome  # Binary outcome for SoFR logistic testing
    )
    values$time_labels <- sprintf("%02d:00", 0:23)
    values$raw_df <- NULL
    showNotification("Sample data generated! (includes binary 'Outcome' for SoFR testing)", type = "message")
  })
  
  # --- 2. Variable Selection ---
  output$var_select_container <- renderUI({
    req(values$raw_df)
    cols <- colnames(values$raw_df)
    num_cols <- cols[sapply(values$raw_df, is.numeric)]
    tagList(
      h4("Define Data Structure"),
      pickerInput("sel_func_vars", "Select Functional Data Columns (Time Points):", 
                  choices = cols, selected = num_cols, multiple = TRUE, options = list(`actions-box` = TRUE)),
      pickerInput("sel_cov_vars", "Select Scalar Variables (Predictors/Response):", 
                  choices = cols, multiple = TRUE, options = list(`actions-box` = TRUE)),
      actionButton("apply_selection", "Confirm Selection", class = "btn-success")
    )
  })
  
  observeEvent(input$apply_selection, {
    req(values$raw_df, input$sel_func_vars)
    tryCatch({
      values$time_labels <- input$sel_func_vars
      temp_data <- values$raw_df[, input$sel_func_vars, drop=FALSE]
      mat_data <- data.matrix(temp_data) 
      if(typeof(mat_data) != "double" && typeof(mat_data) != "integer") mode(mat_data) <- "numeric"
      values$data <- mat_data
      
      if(!is.null(input$sel_cov_vars)) values$covariates <- values$raw_df[, input$sel_cov_vars, drop = FALSE]
      else values$covariates <- NULL
      
      values$smooth_data <- NULL 
      values$smooth_fit_metrics <- NULL
      showNotification("Data structure applied.", type = "message")
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })
  
  # --- 3. Preprocessing / Smoothing ---
  output$data_preview <- renderDT({
    if(is.null(values$data)) return(NULL)
    df_p <- as.data.frame(values$data[, 1:min(5, ncol(values$data))])
    colnames(df_p) <- paste0("T", 1:ncol(df_p))
    if(!is.null(values$covariates)) df_p <- cbind(values$covariates, df_p)
    datatable(df_p, options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$raw_data_plot <- renderPlot({
    if(is.null(values$data)) return(NULL)
    matplot(t(values$data), type='l', lty=1, col=rgb(0,0,0,0.1), main="Raw Data", ylab="Value")
  })
  
  # Display smoothing fit metrics
  output$smooth_fit_display <- renderUI({
    if(is.null(values$smooth_fit_metrics)) return(NULL)
    
    metrics <- values$smooth_fit_metrics
    
    div(
      style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border: 1px solid #dee2e6;",
      h5(icon("chart-line"), " Smoothing Fit Metrics", style = "margin-top: 0;"),
      p(style = "margin-bottom: 5px;",
        strong("Mean R²: "), sprintf("%.3f (SD: %.3f)", metrics$mean_r_squared, metrics$sd_r_squared)),
      p(style = "margin-bottom: 5px;",
        strong("Mean RMSE: "), sprintf("%.3f (SD: %.3f)", metrics$mean_rmse, metrics$sd_rmse)),
      helpText(style = "font-size: 0.85em; margin-top: 8px;",
               "R² = proportion of variance explained by smooth curve (higher = better fit).",
               br(),
               "RMSE = average deviation from original values (lower = better fit).")
    )
  })
  
  observeEvent(input$apply_smooth, {
    req(values$data)
    tryCatch({
      n_time <- ncol(values$data)
      n_subjects <- nrow(values$data)
      t_full <- 1:n_time 
      
      # Create basis
      if(input$is_cyclic) {
        basis <- create.fourier.basis(rangeval = c(1, n_time), nbasis = min(n_time, 13)) 
      } else {
        if(input$smooth_method == "none") {
          basis <- create.bspline.basis(rangeval = c(1, n_time), breaks = t_full, norder=4)
        } else {
          nb <- input$n_basis
          nb <- min(nb, n_time)
          if(nb < 4) nb <- 4
          basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
        }
      }
      
      # Check for NAs in data
      has_any_na <- any(is.na(values$data))
      
      if(input$smooth_method == "none") {
        values$smooth_data <- values$data 
        # For fd_obj, need to handle NAs - use mean imputation temporarily
        temp_data <- values$data
        for(i in 1:n_subjects) {
          na_idx <- is.na(temp_data[i, ])
          if(any(na_idx)) {
            temp_data[i, na_idx] <- mean(temp_data[i, !na_idx], na.rm = TRUE)
          }
        }
        values$fd_obj <- smooth.basis(t_full, t(temp_data), basis)$fd
        values$smooth_fit_metrics <- NULL  # No metrics for raw data
        showNotification("Using Raw Data (No Smoothing).", type="message")
        
      } else {
        # Smooth with NA handling - process each subject individually
        lam <- if(input$smooth_method=="auto") 0 else 10^(-input$smooth_factor)
        fdParobj <- fdPar(basis, 2, lam)
        
        # Initialize smoothed data matrix
        smooth_mat <- matrix(NA, nrow = n_subjects, ncol = n_time)
        fd_coefs_list <- list()
        n_failed <- 0
        failed_subjects <- c()
        
        withProgress(message = 'Smoothing subjects...', value = 0, {
          for(i in 1:n_subjects) {
            y_i <- values$data[i, ]
            valid_idx <- !is.na(y_i)
            n_valid <- sum(valid_idx)
            
            # Need at least 4 points for B-spline smoothing
            min_points_needed <- if(input$is_cyclic) 3 else 4
            
            if(n_valid >= min_points_needed) {
              tryCatch({
                t_valid <- t_full[valid_idx]
                y_valid <- y_i[valid_idx]
                
                # Smooth using only valid points
                fd_i <- smooth.basis(t_valid, y_valid, fdParobj)$fd
                
                # Evaluate at ALL time points (this interpolates over NAs)
                smooth_mat[i, ] <- as.vector(eval.fd(t_full, fd_i))
                fd_coefs_list[[i]] <- fd_i$coefs
              }, error = function(e) {
                # If smoothing fails, use linear interpolation as fallback
                smooth_mat[i, ] <<- approx(t_full[valid_idx], y_i[valid_idx], 
                                           xout = t_full, rule = 2)$y
                n_failed <<- n_failed + 1
                failed_subjects <<- c(failed_subjects, i)
              })
            } else if(n_valid > 0) {
              # Too few points for basis smoothing - use linear interpolation
              smooth_mat[i, ] <- approx(t_full[valid_idx], y_i[valid_idx], 
                                        xout = t_full, rule = 2)$y
            } else {
              # All NA - cannot smooth
              n_failed <- n_failed + 1
              failed_subjects <- c(failed_subjects, i)
            }
            
            if(i %% 10 == 0) incProgress(10 / n_subjects)
          }
        })

        # Apply boundary constraints if requested
        if(input$constrain_bounds) {
          min_val <- input$min_bound
          max_val <- input$max_bound
          if(!is.na(min_val) && !is.na(max_val) && min_val < max_val) {
            smooth_mat <- pmin(pmax(smooth_mat, min_val), max_val)
          }
        }

        values$smooth_data <- smooth_mat
        
        # Create fd object from smoothed data for compatibility
        values$fd_obj <- smooth.basis(t_full, t(smooth_mat), fdParobj)$fd
        
        # Compute fit metrics (comparing smoothed to original where original was not NA)
        r_squared_vec <- numeric(n_subjects)
        rmse_vec <- numeric(n_subjects)
        
        for(i in 1:n_subjects) {
          orig_i <- values$data[i, ]
          smooth_i <- smooth_mat[i, ]
          valid_idx <- !is.na(orig_i)
          
          if(sum(valid_idx) > 1) {
            orig_valid <- orig_i[valid_idx]
            smooth_valid <- smooth_i[valid_idx]
            
            # RMSE
            rmse_vec[i] <- sqrt(mean((orig_valid - smooth_valid)^2))
            
            # R-squared
            ss_tot <- sum((orig_valid - mean(orig_valid))^2)
            ss_res <- sum((orig_valid - smooth_valid)^2)
            r_squared_vec[i] <- if(ss_tot > 0) 1 - ss_res / ss_tot else NA
          } else {
            r_squared_vec[i] <- NA
            rmse_vec[i] <- NA
          }
        }
        
        # Store fit metrics
        values$smooth_fit_metrics <- list(
          r_squared = r_squared_vec,
          rmse = rmse_vec,
          mean_r_squared = mean(r_squared_vec, na.rm = TRUE),
          mean_rmse = mean(rmse_vec, na.rm = TRUE),
          sd_r_squared = sd(r_squared_vec, na.rm = TRUE),
          sd_rmse = sd(rmse_vec, na.rm = TRUE)
        )
        
        # Report results
        n_with_na <- sum(apply(values$data, 1, function(x) any(is.na(x))))
        n_interpolated <- n_with_na - length(failed_subjects)
        
        if(n_interpolated > 0) {
          showNotification(
            sprintf("Smoothing applied. %d subjects had missing values - successfully interpolated. Mean R²=%.3f, Mean RMSE=%.3f", 
                    n_interpolated, values$smooth_fit_metrics$mean_r_squared, values$smooth_fit_metrics$mean_rmse),
            type = "message", duration = 8)
        } else {
          showNotification(
            sprintf("Smoothing applied. Mean R²=%.3f (SD=%.3f), Mean RMSE=%.3f (SD=%.3f)", 
                    values$smooth_fit_metrics$mean_r_squared, values$smooth_fit_metrics$sd_r_squared,
                    values$smooth_fit_metrics$mean_rmse, values$smooth_fit_metrics$sd_rmse),
            type = "message", duration = 8)
        }
        
        if(length(failed_subjects) > 0) {
          showNotification(
            sprintf("Warning: %d subjects could not be smoothed (all NA or too few points): %s", 
                    length(failed_subjects), paste(head(failed_subjects, 10), collapse = ", ")),
            type = "warning", duration = 8)
        }
      }
    }, error = function(e) showNotification(paste("Smoothing Error:", e$message), type="error"))
  })

  # ============================================================================
  # SMOOTHING DIAGNOSTICS SERVER CODE
  # ============================================================================

  # 1. GAM REML Analysis
  observeEvent(input$run_gam_reml, {
    req(values$data)

    showNotification("Fitting GAM with REML...", type = "message", duration = 3)

    tryCatch({
      n_time <- ncol(values$data)
      time_points <- 1:n_time

      # Select data to analyze
      if(input$diag_subject_type == "mean") {
        y_data <- colMeans(values$data, na.rm = TRUE)
        data_label <- "Mean curve (all subjects)"
      } else {
        subject_id <- min(input$diag_subject_id, nrow(values$data))
        y_data <- values$data[subject_id, ]
        data_label <- paste0("Subject ", subject_id)
      }

      # Ensure y_data is numeric and handle NAs
      y_data <- as.numeric(y_data)
      if(all(is.na(y_data))) {
        showNotification("All data points are NA. Cannot perform GAM REML.", type = "error")
        return()
      }

      # Remove NAs
      valid_idx <- which(!is.na(y_data))
      if(length(valid_idx) == 0) {
        showNotification("No valid data points for GAM", type = "error")
        return()
      }

      y_valid <- y_data[valid_idx]
      t_valid <- time_points[valid_idx]

      if(length(y_valid) < 4) {
        showNotification("Not enough valid data points for GAM", type = "error")
        return()
      }

      # Fit GAM with REML
      df_gam <- data.frame(y = y_valid, t = t_valid)
      gam_fit <- mgcv::gam(y ~ s(t, bs = "cr"), data = df_gam, method = "REML")

      # Store results
      values$gam_reml_result <- gam_fit
      values$gam_data_label <- data_label

      showNotification("GAM REML completed successfully!", type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("GAM REML error:", e$message), type = "error", duration = 5)
      cat("GAM REML error:", e$message, "\n")
    })
  })

  output$gam_reml_summary <- renderPrint({
    if(is.null(values$gam_reml_result)) {
      cat("No GAM REML analysis performed yet.\n\n")
      cat("Click 'Fit GAM (REML)' to run the analysis.")
      return()
    }

    gam_fit <- values$gam_reml_result

    cat("=== GAM REML Analysis ===\n")
    cat("Data:", values$gam_data_label, "\n\n")

    # Extract key info
    summary_gam <- summary(gam_fit)
    edf <- summary_gam$edf

    # Safely extract k (basis dimension)
    if("k'" %in% colnames(summary_gam$s.table)) {
      k <- summary_gam$s.table[1, "k'"]
    } else if("k" %in% colnames(summary_gam$s.table)) {
      k <- summary_gam$s.table[1, "k"]
    } else {
      k <- NA
    }

    reml_score <- gam_fit$gcv.ubre  # REML score

    cat("Smoothing parameters:\n")
    cat(sprintf("  Effective Degrees of Freedom (EDF): %.2f\n", edf))
    if(!is.na(k)) {
      cat(sprintf("  Basis dimension (k): %.0f\n", k))
    }
    cat(sprintf("  Estimated lambda: %.3e\n", gam_fit$sp))
    cat(sprintf("  REML score: %.2f\n", reml_score))

    cat("\nInterpretation:\n")
    if(edf < 2) {
      cat("  -> Very smooth fit (low complexity)\n")
      cat("  -> May be undersmoothing important features\n")
    } else if(edf < 6) {
      cat("  -> Moderate smoothing (balanced complexity)\n")
      cat("  -> Often optimal for functional data\n")
    } else if(!is.na(k) && edf < k - 2) {
      cat("  -> Flexible fit (high complexity)\n")
      cat("  -> Capturing detailed features\n")
    } else if(!is.na(k) && edf >= k - 2) {
      cat("  -> Minimal smoothing (very high complexity)\n")
      cat("  -> May be overfitting, consider increasing lambda\n")
    } else {
      cat("  -> Flexible fit\n")
      cat("  -> Capturing detailed features\n")
    }

    cat("\nModel deviance explained:", sprintf("%.1f%%", summary_gam$dev.expl * 100), "\n")

    # Add smoothing term details
    cat("\n--- Smooth term details ---\n")
    print(summary_gam$s.table)
  })

  # 2. REML Profile Analysis
  observeEvent(input$run_reml_profile, {
    req(values$data)

    showNotification("Computing REML profile...", type = "message", duration = 3)

    tryCatch({
      n_time <- ncol(values$data)
      time_points <- 1:n_time

      # Select data
      if(input$diag_subject_type == "mean") {
        y_data <- colMeans(values$data, na.rm = TRUE)
      } else {
        subject_id <- min(input$diag_subject_id, nrow(values$data))
        y_data <- values$data[subject_id, ]
      }

      valid_idx <- !is.na(y_data)
      y_valid <- y_data[valid_idx]
      t_valid <- time_points[valid_idx]

      # Lambda range
      lambda_min <- 10^input$lambda_min_exp
      lambda_max <- 10^input$lambda_max_exp
      lambda_seq <- 10^seq(input$lambda_min_exp, input$lambda_max_exp,
                           length.out = input$n_lambda)

      # Compute REML score for each lambda
      reml_scores <- numeric(length(lambda_seq))
      edf_values <- numeric(length(lambda_seq))

      withProgress(message = "Computing REML profile...", value = 0, {
        for(i in seq_along(lambda_seq)) {
          df_gam <- data.frame(y = y_valid, t = t_valid)

          # Fit with fixed lambda (sp = smoothing parameter)
          gam_fit <- mgcv::gam(y ~ s(t, bs = "cr"), data = df_gam,
                               method = "REML", sp = lambda_seq[i])

          reml_scores[i] <- gam_fit$gcv.ubre
          edf_values[i] <- summary(gam_fit)$edf

          incProgress(1 / length(lambda_seq))
        }
      })

      # Store results
      values$reml_profile_result <- list(
        lambda = lambda_seq,
        reml_score = reml_scores,
        edf = edf_values,
        optimal_idx = which.min(reml_scores),
        optimal_lambda = lambda_seq[which.min(reml_scores)],
        optimal_reml = min(reml_scores)
      )

      showNotification(
        sprintf("REML profile complete! Optimal lambda: %.2e",
                values$reml_profile_result$optimal_lambda),
        type = "message", duration = 5)

    }, error = function(e) {
      showNotification(paste("REML profile error:", e$message), type = "error", duration = 5)
      cat("REML profile error:", e$message, "\n")
    })
  })

  output$reml_profile_plot <- renderPlotly({
    if(is.null(values$reml_profile_result)) {
      return(plot_ly() %>%
               layout(title = "No REML profile computed yet",
                      annotations = list(text = "Click 'Compute REML Profile' to run analysis",
                                         showarrow = FALSE)))
    }

    profile <- values$reml_profile_result

    # Create plot
    p <- plot_ly() %>%
      add_trace(x = log10(profile$lambda),
                y = profile$reml_score,
                type = 'scatter',
                mode = 'lines+markers',
                name = 'REML Score',
                line = list(color = 'blue', width = 2),
                marker = list(size = 6)) %>%
      add_trace(x = log10(profile$optimal_lambda),
                y = profile$optimal_reml,
                type = 'scatter',
                mode = 'markers',
                name = 'Optimal',
                marker = list(size = 12, color = 'red', symbol = 'star')) %>%
      layout(
        title = "REML Profile: Optimal Smoothing Parameter",
        xaxis = list(title = "log10(Lambda)"),
        yaxis = list(title = "REML Score (lower is better)"),
        hovermode = 'closest',
        showlegend = TRUE
      )

    # Add annotation for optimal value
    p <- p %>% layout(
      annotations = list(
        x = log10(profile$optimal_lambda),
        y = profile$optimal_reml,
        text = sprintf("Optimal lambda = %.2e<br>EDF = %.2f",
                       profile$optimal_lambda,
                       profile$edf[profile$optimal_idx]),
        showarrow = TRUE,
        arrowhead = 2,
        ax = 40,
        ay = -40
      )
    )

    p
  })

  # 3. Cross-Validation Analysis
  observeEvent(input$run_cv_analysis, {
    req(values$data)

    showNotification("Running K-fold cross-validation...", type = "message", duration = 3)

    tryCatch({
      n_time <- ncol(values$data)
      n_subjects <- nrow(values$data)
      time_points <- 1:n_time
      k_folds <- input$cv_k_folds

      # Lambda range
      lambda_seq <- 10^seq(input$lambda_min_exp, input$lambda_max_exp,
                           length.out = input$n_lambda)

      # Create fold assignments (random)
      fold_assignments <- sample(rep(1:k_folds, length.out = n_subjects))

      # CV error matrix: subjects x lambdas
      cv_errors <- matrix(NA, nrow = n_subjects, ncol = length(lambda_seq))

      withProgress(message = "Running cross-validation...", value = 0, {
        total_iter <- n_subjects * length(lambda_seq)
        iter_count <- 0

        for(lambda_idx in seq_along(lambda_seq)) {
          lambda <- lambda_seq[lambda_idx]

          # Create basis for this lambda
          nb <- min(input$n_basis, n_time - 2)
          basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
          fdParobj <- fdPar(basis, 2, lambda)

          for(i in 1:n_subjects) {
            # Get subject's fold
            test_fold <- fold_assignments[i]
            train_idx <- which(fold_assignments != test_fold)

            # Train on other subjects
            train_data <- values$data[train_idx, , drop = FALSE]

            # Get mean of training data
            train_mean <- colMeans(train_data, na.rm = TRUE)

            # Smooth training mean
            valid_train <- !is.na(train_mean)
            if(sum(valid_train) >= 4) {
              fd_train <- smooth.basis(time_points[valid_train],
                                       train_mean[valid_train],
                                       fdParobj)$fd

              # Predict on test subject
              y_test <- values$data[i, ]
              valid_test <- !is.na(y_test)

              if(sum(valid_test) >= 1) {
                pred_test <- eval.fd(time_points[valid_test], fd_train)

                # Compute prediction error
                cv_errors[i, lambda_idx] <- sqrt(mean((y_test[valid_test] - pred_test)^2))
              }
            }

            iter_count <- iter_count + 1
            if(iter_count %% 50 == 0) {
              incProgress(50 / total_iter)
            }
          }
        }
      })

      # Compute mean CV error and SE for each lambda
      mean_cv_error <- colMeans(cv_errors, na.rm = TRUE)
      se_cv_error <- apply(cv_errors, 2, sd, na.rm = TRUE) / sqrt(n_subjects)

      optimal_idx <- which.min(mean_cv_error)
      optimal_lambda <- lambda_seq[optimal_idx]

      # 1-SE rule: most parsimonious model within 1 SE of minimum
      se_threshold <- mean_cv_error[optimal_idx] + se_cv_error[optimal_idx]
      within_1se <- which(mean_cv_error <= se_threshold)
      lambda_1se <- lambda_seq[max(within_1se)]  # Highest lambda (most smooth) within 1 SE

      # Store results
      values$cv_result <- list(
        lambda = lambda_seq,
        mean_error = mean_cv_error,
        se_error = se_cv_error,
        optimal_idx = optimal_idx,
        optimal_lambda = optimal_lambda,
        lambda_1se = lambda_1se,
        k_folds = k_folds
      )

      showNotification(
        sprintf("CV complete! Optimal lambda: %.2e (1-SE: %.2e)",
                optimal_lambda, lambda_1se),
        type = "message", duration = 5)

    }, error = function(e) {
      showNotification(paste("CV error:", e$message), type = "error", duration = 5)
      cat("CV error:", e$message, "\n")
    })
  })

  output$cv_curve_plot <- renderPlotly({
    if(is.null(values$cv_result)) {
      return(plot_ly() %>%
               layout(title = "No CV analysis performed yet",
                      annotations = list(text = "Click 'Run Cross-Validation' to start",
                                         showarrow = FALSE)))
    }

    cv <- values$cv_result

    # Create plot with error bars
    p <- plot_ly() %>%
      add_trace(x = log10(cv$lambda),
                y = cv$mean_error,
                type = 'scatter',
                mode = 'lines',
                name = 'Mean CV Error',
                line = list(color = 'darkgreen', width = 2)) %>%
      add_trace(x = log10(cv$lambda),
                y = cv$mean_error + cv$se_error,
                type = 'scatter',
                mode = 'lines',
                name = 'Mean + SE',
                line = list(color = 'lightgreen', width = 1, dash = 'dash'),
                showlegend = FALSE) %>%
      add_trace(x = log10(cv$lambda),
                y = cv$mean_error - cv$se_error,
                type = 'scatter',
                mode = 'lines',
                name = 'Mean - SE',
                line = list(color = 'lightgreen', width = 1, dash = 'dash'),
                fill = 'tonexty',
                fillcolor = 'rgba(144, 238, 144, 0.2)',
                showlegend = FALSE) %>%
      add_trace(x = log10(cv$optimal_lambda),
                y = cv$mean_error[cv$optimal_idx],
                type = 'scatter',
                mode = 'markers',
                name = 'Optimal (min)',
                marker = list(size = 12, color = 'red', symbol = 'star')) %>%
      add_trace(x = log10(cv$lambda_1se),
                y = cv$mean_error[which(abs(cv$lambda - cv$lambda_1se) < 1e-10)],
                type = 'scatter',
                mode = 'markers',
                name = '1-SE rule',
                marker = list(size = 12, color = 'orange', symbol = 'diamond')) %>%
      layout(
        title = sprintf("Cross-Validation Curve (K=%d folds)", cv$k_folds),
        xaxis = list(title = "log10(Lambda)"),
        yaxis = list(title = "CV Error (RMSE)"),
        hovermode = 'closest',
        showlegend = TRUE
      )

    p
  })

  # 4. Comparison Summary
  output$smoothing_comparison_summary <- renderPrint({
    if(is.null(values$reml_profile_result) && is.null(values$cv_result)) {
      cat("No smoothing diagnostic analysis performed yet.\n\n")
      cat("Run REML Profile and/or Cross-Validation to see comparison.")
      return()
    }

    cat("=== Smoothing Parameter Comparison ===\n\n")

    if(!is.null(values$reml_profile_result)) {
      cat("REML Analysis:\n")
      cat(sprintf("  Optimal lambda: %.3e\n", values$reml_profile_result$optimal_lambda))
      cat(sprintf("  Optimal EDF: %.2f\n",
                  values$reml_profile_result$edf[values$reml_profile_result$optimal_idx]))
      cat(sprintf("  REML score: %.2f\n\n", values$reml_profile_result$optimal_reml))
    }

    if(!is.null(values$cv_result)) {
      cat("Cross-Validation Analysis:\n")
      cat(sprintf("  Optimal lambda (min CV): %.3e\n", values$cv_result$optimal_lambda))
      cat(sprintf("  Lambda (1-SE rule): %.3e\n", values$cv_result$lambda_1se))
      cat(sprintf("  Min CV error: %.3f\n",
                  values$cv_result$mean_error[values$cv_result$optimal_idx]))
      cat(sprintf("  K-folds: %d\n\n", values$cv_result$k_folds))
    }

    if(!is.null(values$reml_profile_result) && !is.null(values$cv_result)) {
      reml_lambda <- values$reml_profile_result$optimal_lambda
      cv_lambda <- values$cv_result$optimal_lambda
      cv_lambda_1se <- values$cv_result$lambda_1se

      ratio_min <- cv_lambda / reml_lambda
      ratio_1se <- cv_lambda_1se / reml_lambda

      cat("Comparison:\n")
      cat(sprintf("  CV optimal / REML optimal: %.2f\n", ratio_min))
      cat(sprintf("  CV 1-SE / REML optimal: %.2f\n\n", ratio_1se))

      if(ratio_min > 0.5 && ratio_min < 2) {
        cat("  -> REML and CV agree well on optimal smoothing\n")
      } else if(ratio_min > 2) {
        cat("  -> CV prefers more smoothing than REML\n")
        cat("  -> Consider using CV lambda for better prediction\n")
      } else {
        cat("  -> CV prefers less smoothing than REML\n")
        cat("  -> Consider validating with held-out data\n")
      }
    }

    cat("\nRecommendation:\n")
    if(!is.null(values$cv_result)) {
      cat(sprintf("  For prediction tasks: lambda = %.2e (CV optimal)\n",
                  values$cv_result$optimal_lambda))
      cat(sprintf("  For smooth visualization: lambda = %.2e (1-SE rule)\n",
                  values$cv_result$lambda_1se))
    }
    if(!is.null(values$reml_profile_result)) {
      cat(sprintf("  For automatic REML: lambda = 0 (or %.2e from profile)\n",
                  values$reml_profile_result$optimal_lambda))
    }

    # Add smoothing factor conversion for manual smoothing
    cat("\n--- For Manual Smoothing ---\n")
    cat("To use these lambdas, set 'Smoothing Factor' to:\n")

    if(!is.null(values$cv_result)) {
      sf_cv <- -log10(values$cv_result$optimal_lambda)
      sf_1se <- -log10(values$cv_result$lambda_1se)
      cat(sprintf("  CV optimal: %.2f (lambda = %.2e)\n", sf_cv, values$cv_result$optimal_lambda))
      cat(sprintf("  1-SE rule: %.2f (lambda = %.2e)\n", sf_1se, values$cv_result$lambda_1se))
    }
    if(!is.null(values$reml_profile_result)) {
      sf_reml <- -log10(values$reml_profile_result$optimal_lambda)
      cat(sprintf("  REML optimal: %.2f (lambda = %.2e)\n", sf_reml, values$reml_profile_result$optimal_lambda))
    }

    cat("\nOr use the '📊 Use Diagnostic Results' button in Smoothing tab!")
  })

  # ============================================================================
  # LINK DIAGNOSTICS TO SMOOTHING TAB
  # ============================================================================

  # Check if diagnostic results are available
  output$diagnostics_available <- reactive({
    !is.null(values$reml_profile_result) || !is.null(values$cv_result)
  })
  outputOptions(output, "diagnostics_available", suspendWhenHidden = FALSE)

  # Apply diagnostic lambda to smoothing factor
  observeEvent(input$use_diagnostic_lambda, {
    if(is.null(values$reml_profile_result) && is.null(values$cv_result)) {
      showNotification("No diagnostic results available. Please run REML Profile or Cross-Validation first.",
                       type = "warning", duration = 5)
      return()
    }

    # Determine which lambda to use
    # Priority: CV optimal > REML optimal
    lambda_to_use <- NULL
    lambda_source <- NULL

    if(!is.null(values$cv_result)) {
      lambda_to_use <- values$cv_result$optimal_lambda
      lambda_source <- "CV optimal"
    } else if(!is.null(values$reml_profile_result)) {
      lambda_to_use <- values$reml_profile_result$optimal_lambda
      lambda_source <- "REML profile optimal"
    }

    # Convert lambda to smoothing factor
    # lambda = 10^(-smooth_factor) → smooth_factor = -log10(lambda)
    smooth_factor <- -log10(lambda_to_use)

    # Clamp to slider range [0.1, 10]
    smooth_factor <- max(0.1, min(10, smooth_factor))

    # Update slider
    updateSliderInput(session, "smooth_factor", value = smooth_factor)

    # Also switch to manual mode if not already
    updateRadioButtons(session, "smooth_method", selected = "manual")

    showNotification(
      sprintf("Smoothing factor set to %.2f (lambda = %.2e) from %s",
              smooth_factor, lambda_to_use, lambda_source),
      type = "message", duration = 8)
  })

  # ============================================================================
  # END SMOOTHING DIAGNOSTICS
  # ============================================================================

  output$data_plot <- renderPlotly({
    req(values$data)
    Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
    t <- 1:ncol(Y)
    p <- plot_ly()
    n_show <- min(50, nrow(Y))
    for(i in 1:n_show) {
      p <- add_trace(p, x = t, y = Y[i,], type = 'scatter', mode = 'lines',
                     line = list(color='rgba(0,0,0,0.2)', width=1), 
                     showlegend = FALSE, hoverinfo = "none")
    }
    p <- add_trace(p, x = t, y = colMeans(Y, na.rm=TRUE), 
                   type = 'scatter', mode = 'lines',
                   line = list(color = 'red', width = 3), name = "Mean") %>%
      layout(title = "Functional Data (First 50)", xaxis = list(title = "Time"), yaxis = list(title = "Value"))
    if(!is.null(values$time_labels)) {
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = 1:length(values$time_labels), ticktext = values$time_labels, tickangle = -90))
    }
    p
  })
  
  # ==============================================================================
  # FoSR LOGIC (Function-on-Scalar)
  # ==============================================================================
  
  output$fosr_var_select_ui <- renderUI({
    req(values$covariates)
    tagList(
      selectInput("reg_predictors", "Select Scalar Predictors:", choices = colnames(values$covariates), multiple = TRUE),
      helpText("Response is the FUNCTIONAL data.")
    )
  })
  
  observeEvent(input$run_fosr, {
    req(values$data, values$covariates, input$reg_predictors)
    showNotification("Running FoSR...", type = "message", duration = 2)
    tryCatch({
      Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      mode(Y) <- "numeric"
      df_reg <- values$covariates
      
      # Subset Complete Cases
      keep_idx <- complete.cases(df_reg[, input$reg_predictors, drop=FALSE]) & complete.cases(Y)
      if(sum(keep_idx) < nrow(df_reg)) {
        df_reg <- df_reg[keep_idx, , drop=FALSE]
        Y <- Y[keep_idx, , drop=FALSE]
      }
      
      n_time <- ncol(Y)
      time_points <- seq(0, 1, length.out = n_time)
      
      # --- METHOD A: Pointwise OLS ---
      if (input$reg_method == "OLS_nosmooth") {
        f <- as.formula(paste("~", paste(input$reg_predictors, collapse = " + ")))
        X <- model.matrix(f, data = df_reg)
        
        # Point Estimation
        xtx_inv <- solve(crossprod(X))
        projector <- xtx_inv %*% t(X) 
        coefs <- projector %*% Y
        fitted_vals <- X %*% coefs
        residuals <- Y - fitted_vals
        n <- nrow(Y); p <- ncol(X)
        sigma2 <- colSums(residuals^2) / (n - p)
        
        se_mat <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
        p_values <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
        boot_ci_lower <- NULL
        boot_ci_upper <- NULL
        
        if(input$use_bootstrap) {
          showNotification("Bootstrapping...", type = "message")
          B <- input$n_boot
          boot_betas <- array(NA, dim = c(B, nrow(coefs), ncol(coefs)))
          
          withProgress(message = 'Running Bootstrap...', value = 0, {
            for(b in 1:B) {
              # Residual bootstrap: resample residuals, add to fitted
              resid_idx <- sample(1:n, n, replace = TRUE)
              Y_boot <- fitted_vals + residuals[resid_idx, ]
              boot_betas[b, , ] <- projector %*% Y_boot
              if(b %% 20 == 0) incProgress(20/B)
            }
          })
          
          # Calculate SE from bootstrap distribution
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) se_mat[j, k] <- sd(boot_betas[, j, k])
          }
          
          # Percentile-based 95% CI
          boot_ci_lower <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          boot_ci_upper <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) {
              boot_ci_lower[j, k] <- quantile(boot_betas[, j, k], 0.025)
              boot_ci_upper[j, k] <- quantile(boot_betas[, j, k], 0.975)
            }
          }
          
          # Bootstrap p-values: proportion of bootstrap samples on opposite side of 0
          # This is a proper bootstrap test
          p_values <- matrix(NA, nrow=nrow(coefs), ncol=ncol(coefs))
          for(j in 1:nrow(coefs)) {
            for(k in 1:ncol(coefs)) {
              if(coefs[j, k] >= 0) {
                p_values[j, k] <- 2 * mean(boot_betas[, j, k] <= 0)
              } else {
                p_values[j, k] <- 2 * mean(boot_betas[, j, k] >= 0)
              }
              # Ensure p-value is in [0, 1]
              p_values[j, k] <- min(1, max(0, p_values[j, k]))
            }
          }
          
        } else {
          # Parametric SE
          for(j in 1:nrow(coefs)) se_mat[j, ] <- sqrt(xtx_inv[j,j] * sigma2)
          t_stats <- coefs / se_mat
          p_values <- 2 * (1 - pt(abs(t_stats), df = n - p))
        }
        
        rss <- colSums(residuals^2); y_bar <- colMeans(Y)
        tss <- colSums(sweep(Y, 2, y_bar)^2); r2_t <- 1 - (rss/tss)
        
        fit <- list(beta.hat = coefs, fitted.values = fitted_vals, resid = residuals, 
                    beta.se = se_mat, beta.p = p_values, terms = terms(f), r2_t = r2_t, 
                    method = ifelse(input$use_bootstrap, "OLS (Bootstrap SE)", "OLS (Parametric SE)"), 
                    xtx_inv = xtx_inv, sigma2 = sigma2,
                    boot_ci_lower = boot_ci_lower, boot_ci_upper = boot_ci_upper,
                    n_boot = if(input$use_bootstrap) input$n_boot else NULL)
        
      } else {
        # --- METHOD B: Smoothed OLS (GAM) ---
        showNotification("Fitting GAM with Splines...", type = "message")
        
        # IMPORTANT: Store original factor levels BEFORE any subsetting
        # This is needed for prediction to work correctly
        orig_factor_levels <- list()
        for(v in input$reg_predictors) {
          if(is.factor(df_reg[[v]])) {
            orig_factor_levels[[v]] <- levels(df_reg[[v]])
          } else if(!is.numeric(df_reg[[v]])) {
            orig_factor_levels[[v]] <- sort(unique(as.character(df_reg[[v]])))
          }
        }
        
        # Reshape to Long format
        df_reg$id_temp <- 1:nrow(df_reg)
        long_cov <- df_reg[rep(seq_len(nrow(df_reg)), each = n_time), ]
        long_data <- long_cov
        long_data$time <- rep(time_points, times = nrow(df_reg))
        long_data$Y_val <- as.vector(t(Y))
        
        # Build Formula: Y ~ s(time) + s(time, by=cov) ...
        gam_formula_str <- "Y_val ~ s(time, bs = 'ps', k = 10)"
        preds <- input$reg_predictors
        
        for(p_var in preds) {
          if(is.numeric(df_reg[[p_var]])) {
            gam_formula_str <- paste0(gam_formula_str, " + s(time, by = ", p_var, ", bs='ps', k=10)")
          } else {
            # Factor interaction
            gam_formula_str <- paste0(gam_formula_str, " + ", p_var, " + s(time, by = ", p_var, ", bs='ps', k=10)")
          }
        }
        
        gam_fit <- gam(as.formula(gam_formula_str), data = long_data, method = "REML")
        
        # Reconstruct Beta(t) for Visualization (Approximation)
        pred_grid <- data.frame(time = time_points)
        beta_hat <- matrix(0, nrow = length(preds) + 1, ncol = n_time)
        rownames(beta_hat) <- c("(Intercept)", preds)
        
        # Intercept approx
        d_int <- pred_grid
        for(v in preds) {
          if(is.numeric(df_reg[[v]])) d_int[[v]] <- 0
          else {
            # Use proper factor levels from long_data (which is what GAM was trained on)
            orig_col <- long_data[[v]]
            if(is.factor(orig_col)) {
              d_int[[v]] <- factor(rep(levels(orig_col)[1], n_time), levels = levels(orig_col))
            } else {
              lvls <- sort(unique(as.character(orig_col)))
              d_int[[v]] <- factor(rep(lvls[1], n_time), levels = lvls)
            }
          }
        }
        beta_hat[1, ] <- predict(gam_fit, newdata = d_int, type = "response")
        
        # Covariate effects approx
        for(i in 1:length(preds)) {
          v <- preds[i]
          if(is.numeric(df_reg[[v]])) {
            d_0 <- d_int; d_0[[v]] <- 0
            d_1 <- d_int; d_1[[v]] <- 1
            beta_hat[i+1, ] <- predict(gam_fit, newdata = d_1) - predict(gam_fit, newdata = d_0)
          }
        }
        
        fitted_vec <- predict(gam_fit, newdata = long_data)
        fitted_vals <- matrix(fitted_vec, nrow = nrow(Y), ncol = n_time, byrow = TRUE)
        residuals <- Y - fitted_vals
        
        rss <- colSums(residuals^2); y_bar <- colMeans(Y)
        tss <- colSums(sweep(Y, 2, y_bar)^2); r2_t <- 1 - (rss/tss)
        
        fit <- list(beta.hat = beta_hat, fitted.values = fitted_vals, resid = residuals, 
                    beta.se = beta_hat*0, beta.p = beta_hat*0, # placeholders
                    terms = terms(as.formula(paste("~", paste(preds, collapse="+")))), 
                    r2_t = r2_t, method = "Smoothed OLS (GAM)",
                    gam_obj = gam_fit,
                    gam_predictors = preds,
                    gam_long_data = long_data,  # Store long_data for exact factor levels
                    gam_factor_levels = orig_factor_levels,  # Store original factor levels
                    gam_n_time = n_time)
      }
      
      values$reg_model <- fit
      updateSelectInput(session, "reg_color_var", choices = input$reg_predictors)
      updateSelectInput(session, "reg_coeff_select", choices = rownames(fit$beta.hat))
      num_preds <- input$reg_predictors[sapply(values$covariates[input$reg_predictors], is.numeric)]
      updateSelectInput(session, "reg_ref_selector", choices = num_preds)
      
      showNotification("FoSR Fitted!", type = "message")
    }, error = function(e) showNotification(paste("Fit Error:", e$message), type="error"))
  })
  
  output$fosr_model_summary <- renderPrint({
    req(values$reg_model)
    mod <- values$reg_model
    cat("Method:", mod$method, "\n")
    cat("Coefficients:", paste(rownames(mod$beta.hat), collapse=", "), "\n")
    if(!is.null(mod$r2_t)) cat("Mean R2:", round(mean(mod$r2_t, na.rm=TRUE), 3), "\n")
    if(!is.null(mod$n_boot)) {
      cat("Bootstrap: Yes (B =", mod$n_boot, "), Percentile CIs\n")
    }
    if(!is.null(mod$gam_obj)) {
      cat("\nFamily:", mod$gam_obj$family$family, "\n")
      cat("Link function:", mod$gam_obj$family$link, "\n\n")
      print(summary(mod$gam_obj))
    }
  })
  
  # FoSR Visualization
  output$reg_observed_plot <- renderPlotly({
    req(values$data, values$covariates, input$reg_color_var)
    Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
    subset_idx <- 1:min(nrow(Y), 200) 
    Y_subset <- Y[subset_idx, , drop=FALSE]
    col_vals_raw <- values$covariates[[input$reg_color_var]][subset_idx]
    
    df_plot <- data.frame(
      Time = rep(seq(0, 1, length.out = ncol(Y)), times = nrow(Y_subset)),
      Value = as.vector(t(Y_subset)),
      ID = as.factor(rep(subset_idx, each = ncol(Y))),
      Color = rep(col_vals_raw, each = ncol(Y))
    )
    
    is_categorical <- length(unique(col_vals_raw)) <= 5 || is.factor(col_vals_raw) || is.character(col_vals_raw)
    if(is_categorical) df_plot$Color <- as.factor(df_plot$Color)
    
    g <- ggplot(df_plot, aes(x = Time, y = Value, group = ID, color = Color)) +
      geom_line(alpha = 0.8, linewidth = 0.5) + theme_minimal() +
      labs(title = paste("Observed Data colored by", input$reg_color_var))
    
    if(is_categorical) g <- g + scale_color_brewer(palette = "Set1") else g <- g + scale_color_viridis_c(option = "viridis")
    
    if(!is.null(values$time_labels)) {
      tick_idx <- seq(0, 1, length.out = length(values$time_labels))
      g <- g + scale_x_continuous(breaks = tick_idx, labels = values$time_labels) +
        theme(axis.text.x = element_text(angle = -90, hjust = 0))
    }
    ggplotly(g, tooltip = c("group", "x", "y", "colour"))
  })
  
  output$reg_ref_selector <- renderUI({
    req(values$reg_model)
    num_preds <- input$reg_predictors[sapply(values$covariates[input$reg_predictors], is.numeric)]
    selectInput("ref_predictor", "Select Predictor for Min/Max Reference:", choices = num_preds)
  })
  
  output$reg_prediction_controls <- renderUI({
    req(values$reg_model, input$reg_predictors)
    lapply(input$reg_predictors, function(var) {
      vals <- values$covariates[[var]]
      if(is.numeric(vals)) {
        sliderInput(paste0("pred_", var), label = var, 
                    min = min(vals, na.rm=TRUE), max = max(vals, na.rm=TRUE), 
                    value = mean(vals, na.rm=TRUE))
      } else {
        # For factors, use the factor levels if available
        if(is.factor(vals)) {
          choices <- levels(vals)
        } else {
          choices <- sort(unique(as.character(vals)))
        }
        selectInput(paste0("pred_", var), label = var, choices = choices, selected = choices[1])
      }
    })
  })
  
  # =============================================================================
  # COMPLETELY REWRITTEN: reg_fitted_plot - Prediction plot for GAM/OLS
  # =============================================================================
  output$reg_fitted_plot <- renderPlotly({
    req(values$reg_model, input$reg_predictors)
    
    mod <- values$reg_model
    n_t <- if(!is.null(mod$gam_n_time)) mod$gam_n_time else {
      if(!is.null(mod$fitted.values)) ncol(mod$fitted.values) else ncol(values$data)
    }
    time_vec <- seq(0, 1, length.out = n_t)
    
    # Collect input values
    input_vals <- list()
    for(var in input$reg_predictors) {
      val <- input[[paste0("pred_", var)]]
      if(is.null(val)) return(NULL)
      input_vals[[var]] <- val
    }
    
    # ==========================================================================
    # GAM Method
    # ==========================================================================
    if(!is.null(mod$gam_obj)) {
      preds <- mod$gam_predictors
      long_data <- mod$gam_long_data
      factor_levels <- mod$gam_factor_levels
      
      # Build prediction dataframe using helper function
      pred_df <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals)
      if(is.null(pred_df)) return(NULL)
      
      # Get prediction with SE using safe helper
      pred_result <- safe_gam_predict(mod$gam_obj, pred_df)
      if(!pred_result$success || is.null(pred_result$fit)) {
        showNotification("Prediction failed", type = "error")
        return(NULL)
      }
      
      y_hat <- pred_result$fit
      se_pred <- pred_result$se
      
      # Ensure vectors are the correct length
      if(length(y_hat) != n_t) {
        showNotification(paste("Prediction length mismatch:", length(y_hat), "vs", n_t), type = "error")
        return(NULL)
      }
      if(length(se_pred) != n_t) {
        se_pred <- rep(0, n_t)  # Fallback to no CI if SE has wrong length
      }
      
    } else {
      # ==========================================================================
      # OLS Method
      # ==========================================================================
      pred_data <- list()
      for(var in input$reg_predictors) {
        val <- input_vals[[var]]
        orig <- values$covariates[[var]]
        if(is.numeric(orig)) {
          pred_data[[var]] <- as.numeric(val)
        } else {
          if(is.factor(orig)) {
            pred_data[[var]] <- factor(val, levels = levels(orig))
          } else {
            pred_data[[var]] <- factor(val, levels = sort(unique(as.character(orig))))
          }
        }
      }
      pred_df <- as.data.frame(pred_data)
      
      betas <- mod$beta.hat
      f_clean <- delete.response(mod$terms)
      X_new <- tryCatch({ model.matrix(f_clean, data = pred_df) }, error = function(e) NULL)
      if(is.null(X_new)) return(NULL)
      y_hat <- as.vector(X_new %*% betas)
      
      if(!is.null(mod$xtx_inv) && !is.null(mod$sigma2)) {
        scale_factor <- as.numeric(X_new %*% mod$xtx_inv %*% t(X_new))
        se_pred <- sqrt(scale_factor * mod$sigma2)
      } else {
        se_pred <- rep(0, length(y_hat))
      }
    }
    
    # ==========================================================================
    # Build main plot
    # ==========================================================================
    # Compute CI bounds
    ci_lower <- y_hat - 1.96 * se_pred
    ci_upper <- y_hat + 1.96 * se_pred
    
    # Create plot data frame for consistent handling
    plot_df <- data.frame(
      time = time_vec,
      y_hat = y_hat,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
    
    p <- plot_ly(data = plot_df, x = ~time) %>%
      add_trace(y = ~y_hat, type = 'scatter', mode = 'lines',
                line = list(color = 'red', width = 3), name = "Predicted Mean")
    
    # Only add ribbon if SE is non-zero
    if(any(se_pred > 0)) {
      p <- p %>% add_ribbons(ymin = ~ci_lower, ymax = ~ci_upper,
                             name = "95% CI", line = list(color = 'transparent'), 
                             fillcolor = 'rgba(255, 0, 0, 0.2)')
    }
    
    # ==========================================================================
    # Reference Lines
    # ==========================================================================
    ref_pred <- input$ref_predictor
    if(!is.null(ref_pred) && ref_pred %in% input$reg_predictors) {
      orig_vals <- values$covariates[[ref_pred]]
      
      if(is.numeric(orig_vals)) {
        min_val <- min(orig_vals, na.rm = TRUE)
        max_val <- max(orig_vals, na.rm = TRUE)
        
        # Create modified input values for min/max
        input_vals_min <- input_vals
        input_vals_max <- input_vals
        input_vals_min[[ref_pred]] <- min_val
        input_vals_max[[ref_pred]] <- max_val
        
        if(!is.null(mod$gam_obj)) {
          # GAM reference lines
          pred_df_min <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals_min)
          pred_df_max <- build_gam_pred_df(time_vec, preds, long_data, factor_levels, input_vals_max)
          
          if(!is.null(pred_df_min) && !is.null(pred_df_max)) {
            pred_min <- safe_gam_predict(mod$gam_obj, pred_df_min)
            pred_max <- safe_gam_predict(mod$gam_obj, pred_df_max)
            
            if(pred_min$success && pred_max$success) {
              p <- p %>% 
                add_lines(x = time_vec, y = pred_min$fit, name = paste("Min", ref_pred), 
                          line = list(color = 'blue', dash = 'dot', width = 2)) %>%
                add_lines(x = time_vec, y = pred_max$fit, name = paste("Max", ref_pred), 
                          line = list(color = 'purple', dash = 'dot', width = 2))
            }
          }
        } else {
          # OLS reference lines
          pred_data_min <- pred_data_max <- list()
          for(var in input$reg_predictors) {
            orig <- values$covariates[[var]]
            if(var == ref_pred) {
              pred_data_min[[var]] <- min_val
              pred_data_max[[var]] <- max_val
            } else {
              val <- input_vals[[var]]
              if(is.numeric(orig)) {
                pred_data_min[[var]] <- as.numeric(val)
                pred_data_max[[var]] <- as.numeric(val)
              } else {
                if(is.factor(orig)) {
                  pred_data_min[[var]] <- factor(val, levels = levels(orig))
                  pred_data_max[[var]] <- factor(val, levels = levels(orig))
                } else {
                  lvls <- sort(unique(as.character(orig)))
                  pred_data_min[[var]] <- factor(val, levels = lvls)
                  pred_data_max[[var]] <- factor(val, levels = lvls)
                }
              }
            }
          }
          
          df_min <- as.data.frame(pred_data_min)
          df_max <- as.data.frame(pred_data_max)
          
          f_clean <- delete.response(mod$terms)
          X_min <- tryCatch({ model.matrix(f_clean, data = df_min) }, error = function(e) NULL)
          X_max <- tryCatch({ model.matrix(f_clean, data = df_max) }, error = function(e) NULL)
          
          if(!is.null(X_min) && !is.null(X_max)) {
            y_min <- as.vector(X_min %*% mod$beta.hat)
            y_max <- as.vector(X_max %*% mod$beta.hat)
            
            p <- p %>% 
              add_lines(x = time_vec, y = y_min, name = paste("Min", ref_pred), 
                        line = list(color = 'blue', dash = 'dot', width = 2)) %>%
              add_lines(x = time_vec, y = y_max, name = paste("Max", ref_pred), 
                        line = list(color = 'purple', dash = 'dot', width = 2))
          }
        }
      }
    }
    
    # Add time labels if available
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, 
                                     ticktext = values$time_labels, tickangle = -90))
    }
    
    p %>% layout(
      title = paste("Predicted Curve (Y ~", paste(input$reg_predictors, collapse = " + "), ")"),
      xaxis = list(title = "Time (normalized)"),
      yaxis = list(title = "Predicted Value")
    )
  })
  
  output$reg_coeff_plot <- renderPlotly({
    req(values$reg_model, input$reg_coeff_select)
    mod <- values$reg_model
    sel <- input$reg_coeff_select
    idx <- which(rownames(mod$beta.hat) == sel)
    if(length(idx) == 0) return(NULL)
    
    beta <- mod$beta.hat[idx, ]
    t <- seq(0, 1, length.out = length(beta))
    
    p <- plot_ly(x = t, y = beta, type = 'scatter', mode = 'lines', name = 'Beta(t)',
                 line = list(color = '#003366', width = 2)) %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 0, 
                   line = list(color = 'black', dash = 'dash', width = 1), showlegend = FALSE)
    
    # Use bootstrap percentile CIs if available, otherwise parametric
    if(!is.null(mod$boot_ci_lower) && !is.null(mod$boot_ci_upper)) {
      ci_lower <- mod$boot_ci_lower[idx, ]
      ci_upper <- mod$boot_ci_upper[idx, ]
      ci_name <- paste0("95% Bootstrap CI (B=", mod$n_boot, ")")
      p <- p %>% add_ribbons(x = t, ymin = ci_lower, ymax = ci_upper,
                             name = ci_name, line = list(color = 'transparent'), 
                             fillcolor = 'rgba(0, 191, 255, 0.3)')
    } else if(sum(mod$beta.se[idx, ]) > 0) {
      se <- mod$beta.se[idx, ]
      p <- p %>% add_ribbons(x = t, ymin = beta - 1.96*se, ymax = beta + 1.96*se,
                             name = "95% Parametric CI", line = list(color = 'transparent'), 
                             fillcolor = 'rgba(0, 191, 255, 0.3)')
    }
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    
    ci_type <- if(!is.null(mod$boot_ci_lower)) "Bootstrap" else "Parametric"
    p %>% layout(title = paste("Coefficient:", sel, "(", ci_type, "CI)"))
  })
  
  output$reg_pvalue_plot <- renderPlotly({
    req(values$reg_model, input$reg_coeff_select)
    mod <- values$reg_model
    sel <- input$reg_coeff_select
    idx <- which(rownames(mod$beta.hat) == sel)
    pvals <- mod$beta.p[idx, ]
    if(all(pvals == 0) || all(is.na(pvals))) return(NULL)
    
    t <- seq(0, 1, length.out = length(pvals))
    p <- plot_ly(x = t, y = pvals, type = 'scatter', mode = 'lines', name = 'P-value',
                 line = list(color = 'darkgreen', width = 2)) %>%
      add_lines(y = 0.05, line = list(color = 'red', dash = 'dash'), name = "p=0.05") %>%
      layout(yaxis = list(type = "log", title = "P-value (log scale)"), title = "Significance over Time")
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p
  })
  
  output$reg_r2_plot <- renderPlotly({
    req(values$reg_model)
    mod <- values$reg_model
    t <- seq(0, 1, length.out = length(mod$r2_t))
    mean_r2 <- mean(mod$r2_t, na.rm=TRUE)
    
    p <- plot_ly(x = t, y = mod$r2_t, type = 'scatter', mode = 'lines', 
                 fill = 'tozeroy', name = 'R-squared',
                 line = list(color = 'purple', width = 2)) %>%
      add_lines(y = mean_r2, line = list(color = 'black', dash = 'dash'), name = paste("Mean R2:", round(mean_r2, 2))) %>%
      layout(title = "Model Fit (Functional R-squared)", yaxis = list(range = c(0, 1), title = "R-squared"))
    
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p
  })
  
  output$reg_residual_plot <- renderPlotly({
    req(values$reg_model)
    resid <- values$reg_model$resid
    t <- seq(0, 1, length.out = ncol(resid)); n <- nrow(resid)
    l2 <- rowSums(resid^2); rks <- rank(l2)
    cols <- colorRampPalette(c("red", "orange", "green", "blue", "purple"))(n)[rks]
    
    idx <- if(n > 200) sample(1:n, 200) else 1:n
    p <- plot_ly(type = 'scatter', mode = 'lines', hoverinfo = "text")
    for(i in idx) {
      p <- add_trace(p, x = t, y = resid[i,], line = list(color = cols[i], width = 1), showlegend = FALSE)
    }
    if(!is.null(values$time_labels)) {
      tick_pos <- seq(0, 1, length.out = length(values$time_labels))
      p <- p %>% layout(xaxis = list(tickmode = "array", tickvals = tick_pos, ticktext = values$time_labels, tickangle = -90))
    }
    p %>% layout(title = "Residuals (Rainbow Plot)")
  })
  
  # ==============================================================================
  # SoFR LOGIC (Scalar-on-Function)
  # ==============================================================================
  
  # Store detected family for conditional panels
  output$sofr_is_binary <- reactive({
    if(is.null(values$sofr_model)) return(FALSE)
    !is.null(values$sofr_model$family) && values$sofr_model$family$family == "binomial"
  })
  outputOptions(output, "sofr_is_binary", suspendWhenHidden = FALSE)
  
  output$sofr_var_select_ui <- renderUI({
    req(values$covariates)
    tagList(
      selectInput("sofr_response", "Select Scalar Response Variable (y):", choices = colnames(values$covariates), multiple = FALSE),
      selectInput("sofr_scalar_preds", "Additional Scalar Predictors (z):", choices = colnames(values$covariates), multiple = TRUE),
      helpText("Functional Predictor (X) is the imported curve data.")
    )
  })
  
  # Dynamic family info based on response variable
  output$sofr_family_info <- renderUI({
    req(input$sofr_response, values$covariates)
    y <- values$covariates[[input$sofr_response]]
    
    if(is.null(y)) return(NULL)
    
    # Handle factors
    if(is.factor(y)) {
      n_levels <- nlevels(y)
      if(n_levels == 2) {
        return(div(class = "alert alert-info", style = "padding: 8px; margin-top: 10px;", 
                   icon("info-circle"), " ",
                   paste0("Detected: Binary factor (", n_levels, " levels: ", 
                          paste(levels(y), collapse = ", "), "). ",
                          "Will convert to 0/1. Recommended: Binomial family with logit link.")))
      } else {
        return(div(class = "alert alert-warning", style = "padding: 8px; margin-top: 10px;", 
                   icon("exclamation-triangle"), " ",
                   paste0("Factor with ", n_levels, " levels detected. ",
                          "Only 2-level factors (binary) are supported for SoFR.")))
      }
    }
    
    if(!is.numeric(y)) return(NULL)
    
    unique_vals <- length(unique(na.omit(y)))
    is_binary <- unique_vals == 2
    is_count <- all(y == floor(y), na.rm = TRUE) && all(y >= 0, na.rm = TRUE)
    is_proportion <- all(y >= 0 & y <= 1, na.rm = TRUE)
    
    msg <- ""
    if(is_binary) {
      msg <- paste0("Detected: Binary variable (", unique_vals, " unique values). ",
                    "Recommended: Binomial family with logit link.")
    } else if(is_proportion && !is_binary) {
      msg <- "Detected: Values in [0,1]. Could be proportions (use Binomial) or continuous (use Gaussian)."
    } else if(is_count && min(y, na.rm=TRUE) >= 0) {
      msg <- "Detected: Non-negative integers. Consider Poisson for count data."
    } else if(all(y > 0, na.rm = TRUE)) {
      msg <- "Detected: Positive continuous. Gaussian or Gamma may be appropriate."
    } else {
      msg <- "Detected: Continuous variable. Gaussian family recommended."
    }
    
    div(class = "alert alert-info", style = "padding: 8px; margin-top: 10px;", 
        icon("info-circle"), " ", msg)
  })
  
  observeEvent(input$run_sofr, {
    req(values$data, values$covariates, input$sofr_response)
    showNotification("Running SoFR...", type = "message", duration = 2)
    tryCatch({
      X_func <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      df_sofr <- values$covariates
      y <- df_sofr[[input$sofr_response]]
      
      # Handle factor response for binary
      if(is.factor(y)) {
        y <- as.numeric(y) - 1  # Convert to 0/1
        showNotification("Factor response converted to 0/1 numeric.", type = "message")
      }
      
      if(!is.numeric(y)) stop("Response variable must be numeric (or factor for binary).")
      
      preds <- input$sofr_scalar_preds
      if(is.null(preds)) preds <- character(0)  # Ensure it's at least empty character vector
      
      # Handle complete cases - be careful with empty preds
      if(length(preds) > 0) {
        keep_idx <- complete.cases(df_sofr[, c(input$sofr_response, preds), drop=FALSE]) & complete.cases(X_func)
      } else {
        keep_idx <- complete.cases(df_sofr[, input$sofr_response, drop=FALSE]) & complete.cases(X_func)
      }
      
      y_clean <- y[keep_idx]
      X_func_clean <- X_func[keep_idx, ]
      df_clean <- df_sofr[keep_idx, , drop=FALSE]
      
      # Get link function values with defaults (inputs may be NULL if conditionalPanel not shown)
      link_binom <- if(!is.null(input$sofr_link_binomial)) input$sofr_link_binomial else "logit"
      link_gauss <- if(!is.null(input$sofr_link_gaussian)) input$sofr_link_gaussian else "identity"
      family_choice <- if(!is.null(input$sofr_family)) input$sofr_family else "auto"
      
      # Determine the family and validate/convert y accordingly
      is_binary_outcome <- FALSE
      unique_y <- sort(unique(na.omit(y_clean)))
      
      if(family_choice == "binomial" || (family_choice == "auto" && length(unique_y) == 2)) {
        # For binomial: ensure y is strictly 0/1
        if(length(unique_y) == 2) {
          # Convert to 0/1 if needed (e.g., if values are 1/2 or other pairs)
          y_clean <- as.numeric(y_clean == max(unique_y))
          showNotification(paste("Binary response recoded: ", min(unique_y), "->0, ", max(unique_y), "->1"), type = "message")
        } else if(all(y_clean >= 0 & y_clean <= 1)) {
          # Already proportions - use as is
          showNotification("Using response as proportions (0-1 range).", type = "message")
        } else {
          stop("For binomial family, response must be binary (2 unique values) or proportions (0-1 range).")
        }
        pfr_family <- binomial(link = link_binom)
        is_binary_outcome <- (length(unique_y) == 2)
      } else if(family_choice == "gaussian") {
        pfr_family <- gaussian(link = link_gauss)
      } else if(family_choice == "poisson") {
        if(any(y_clean < 0) || any(y_clean != floor(y_clean))) {
          showNotification("Warning: Poisson expects non-negative integer counts.", type = "warning")
        }
        pfr_family <- poisson(link = "log")
      } else if(family_choice == "Gamma") {
        if(any(y_clean <= 0)) {
          stop("Gamma family requires strictly positive response values.")
        }
        pfr_family <- Gamma(link = "log")
      } else {
        # Default: auto-detect -> gaussian
        pfr_family <- gaussian(link = "identity")
      }
      
      # Prepare data for pfr - IMPORTANT: remove original response column
      # to ensure pfr uses our converted y variable
      df_clean_no_response <- df_clean
      df_clean_no_response[[input$sofr_response]] <- NULL  # Remove original response
      
      pfr_data <- as.list(df_clean_no_response)
      pfr_data$X_func <- X_func_clean
      pfr_data$y <- y_clean  # Use our validated/converted y
      
      # Debug info
      cat("SoFR Debug:\n")
      cat("  Family:", pfr_family$family, "\n")
      cat("  y range:", range(y_clean), "\n")
      cat("  y unique values:", paste(head(sort(unique(y_clean)), 10), collapse=", "), "\n")
      
      formula_str <- "y ~ lf(X_func, bs='ps', k=15)"
      if(!is.null(preds) && length(preds) > 0) {
        formula_str <- paste(formula_str, "+", paste(preds, collapse = " + "))
      }
      pfr_formula <- as.formula(formula_str)
      
      # Final validation for binomial
      if(pfr_family$family == "binomial") {
        if(any(pfr_data$y < 0 | pfr_data$y > 1, na.rm = TRUE)) {
          stop(paste("After conversion, y values are still outside [0,1]. Range:", 
                     paste(range(pfr_data$y, na.rm=TRUE), collapse=" to ")))
        }
        cat("  Binomial validation passed: y in [0,1]\n")
      }
      
      # Use do.call to force evaluation of all arguments
      # This avoids pfr's non-standard evaluation issues
      fit <- do.call(pfr, list(
        formula = pfr_formula, 
        data = pfr_data, 
        family = pfr_family
      ))
      
      # Store additional info for diagnostics
      fit$y_original <- y_clean
      fit$is_binary <- is_binary_outcome
      
      # Bootstrap for coefficient CIs if requested
      if(isTRUE(input$sofr_use_bootstrap)) {
        B <- if(!is.null(input$sofr_n_boot)) input$sofr_n_boot else 100
        n_obs <- length(y_clean)
        
        showNotification(paste("Running", B, "bootstrap iterations..."), type = "message")
        
        # Extract coefficient from fitted model to get dimensions
        # The functional coefficient is in the smooth term
        coef_orig <- coef(fit)
        
        # Get the functional coefficient values over a grid
        # Use the smooth object to evaluate
        n_grid <- 100
        t_grid <- seq(0, 1, length.out = n_grid)
        
        # Store bootstrap coefficients
        boot_coefs <- matrix(NA, nrow = B, ncol = n_grid)
        boot_success <- 0
        
        withProgress(message = 'Running SoFR Bootstrap...', value = 0, {
          for(b in 1:B) {
            tryCatch({
              # Case resampling
              boot_idx <- sample(1:n_obs, n_obs, replace = TRUE)
              
              # Create bootstrap data
              pfr_data_boot <- list()
              for(nm in names(pfr_data)) {
                if(nm == "X_func") {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx, , drop = FALSE]
                } else if(nm == "y") {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx]
                } else if(length(pfr_data[[nm]]) == n_obs) {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]][boot_idx]
                } else {
                  pfr_data_boot[[nm]] <- pfr_data[[nm]]
                }
              }
              
              # Fit bootstrap model (suppress warnings)
              fit_boot <- suppressWarnings(do.call(pfr, list(
                formula = pfr_formula, 
                data = pfr_data_boot, 
                family = pfr_family
              )))
              
              # Extract functional coefficient using coef.pfr or predict approach
              # Get smooth coefficients
              sm <- fit_boot$smooth[[1]]
              if(!is.null(sm)) {
                # Create prediction data for the smooth
                Xp <- mgcv::PredictMat(sm, data.frame(X_func.tmat = t_grid))
                coef_sm <- coef(fit_boot)[sm$first.para:sm$last.para]
                boot_coefs[b, ] <- as.vector(Xp %*% coef_sm)
                boot_success <- boot_success + 1
              }
              
              if(b %% 10 == 0) incProgress(10/B)
            }, error = function(e) {
              # Skip failed bootstrap iterations
            })
          }
        })
        
        if(boot_success >= B * 0.5) {  # Need at least 50% successful iterations
          # Calculate percentile CIs
          boot_ci_lower <- apply(boot_coefs, 2, quantile, probs = 0.025, na.rm = TRUE)
          boot_ci_upper <- apply(boot_coefs, 2, quantile, probs = 0.975, na.rm = TRUE)
          boot_se <- apply(boot_coefs, 2, sd, na.rm = TRUE)
          
          fit$boot_ci_lower <- boot_ci_lower
          fit$boot_ci_upper <- boot_ci_upper
          fit$boot_se <- boot_se
          fit$boot_t_grid <- t_grid
          fit$n_boot <- boot_success
          
          showNotification(paste("Bootstrap complete:", boot_success, "of", B, "iterations succeeded"), type = "message")
        } else {
          showNotification(paste("Bootstrap had too many failures:", boot_success, "of", B, "succeeded"), type = "warning")
        }
      }
      
      values$sofr_model <- fit
      
      family_msg <- paste0("SoFR Fitted! Family: ", fit$family$family, "(", fit$family$link, ")")
      showNotification(family_msg, type = "message")
      
    }, error = function(e) {
      showNotification(paste("SoFR Error:", e$message), type = "error")
      print(e)  # Also print to console for debugging
    })
  })
  
  output$sofr_inference_summary <- renderPrint({
    req(values$sofr_model)
    fit <- values$sofr_model
    cat("=== Scalar-on-Function Regression Results ===\n\n")
    cat("Family:", fit$family$family, "\n")
    cat("Link function:", fit$family$link, "\n")
    if(!is.null(fit$n_boot)) {
      cat("Bootstrap CIs: Yes (B =", fit$n_boot, ")\n")
    } else {
      cat("Bootstrap CIs: No (using parametric)\n")
    }
    cat("\n")
    summary(fit)
  })
  
  output$sofr_model_diagnostics <- renderUI({
    req(values$sofr_model)
    fit <- values$sofr_model
    
    # Calculate diagnostics
    y_obs <- fit$y_original
    y_pred_link <- fitted(fit)  # On link scale for GLM
    
    if(fit$family$family == "binomial") {
      # For binomial, fitted values are already probabilities
      y_pred_prob <- y_pred_link
      
      # Classification at 0.5 threshold
      y_pred_class <- ifelse(y_pred_prob > 0.5, 1, 0)
      acc <- mean(y_pred_class == y_obs)
      
      # Pseudo R-squared (McFadden)
      null_dev <- fit$null.deviance
      res_dev <- fit$deviance
      pseudo_r2 <- 1 - (res_dev / null_dev)
      
      tagList(
        h4("Model Diagnostics"),
        tags$table(class = "table table-condensed",
                   tags$tr(tags$td("Null Deviance:"), tags$td(round(null_dev, 2))),
                   tags$tr(tags$td("Residual Deviance:"), tags$td(round(res_dev, 2))),
                   tags$tr(tags$td("McFadden Pseudo-R²:"), tags$td(round(pseudo_r2, 3))),
                   tags$tr(tags$td("Accuracy (at p=0.5):"), tags$td(paste0(round(acc * 100, 1), "%")))
        )
      )
    } else {
      # For Gaussian, compute R-squared
      ss_res <- sum((y_obs - y_pred_link)^2)
      ss_tot <- sum((y_obs - mean(y_obs))^2)
      r2 <- 1 - ss_res/ss_tot
      rmse <- sqrt(mean((y_obs - y_pred_link)^2))
      
      tagList(
        h4("Model Diagnostics"),
        tags$table(class = "table table-condensed",
                   tags$tr(tags$td("R-squared:"), tags$td(round(r2, 3))),
                   tags$tr(tags$td("RMSE:"), tags$td(round(rmse, 3))),
                   tags$tr(tags$td("Residual Deviance:"), tags$td(round(fit$deviance, 2)))
        )
      )
    }
  })
  
  output$sofr_coeff_interpretation <- renderUI({
    req(values$sofr_model)
    fit <- values$sofr_model
    
    if(fit$family$family == "binomial" && fit$family$link == "logit") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Logistic):"), " The coefficient β(t) is on the log-odds scale. ",
          "Positive values at time t indicate that higher functional predictor values at that time ",
          "are associated with increased probability of Y=1. ",
          "exp(β(t)) gives the odds ratio for a unit increase in X(t)."
      )
    } else if(fit$family$family == "binomial" && fit$family$link == "probit") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Probit):"), " The coefficient β(t) represents the change in the ",
          "z-score (standard normal quantile) of P(Y=1) for a unit increase in X(t)."
      )
    } else if(fit$family$family == "poisson") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Poisson):"), " The coefficient β(t) is on the log scale. ",
          "exp(β(t)) gives the multiplicative change in the expected count for a unit increase in X(t)."
      )
    } else if(fit$family$link == "log") {
      div(class = "alert alert-warning", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Log link):"), " The coefficient β(t) is on the log scale. ",
          "exp(β(t)) gives the multiplicative effect on E(Y)."
      )
    } else {
      div(class = "alert alert-info", style = "margin-top: 15px;",
          icon("info-circle"), " ",
          strong("Interpretation (Identity link):"), " The coefficient β(t) represents the ",
          "direct additive effect of X(t) on E(Y). A unit increase in X(t) changes E(Y) by β(t)."
      )
    }
  })
  
  output$sofr_coeff_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    df_coef <- tryCatch({ coef(fit) }, error = function(e) NULL)
    validate(need(!is.null(df_coef), "Could not extract functional coefficients."))
    
    x_col <- if("X_func.arg" %in% colnames(df_coef)) "X_func.arg" else colnames(df_coef)[1]
    
    # Determine y-axis label based on family
    y_label <- if(fit$family$family == "binomial" && fit$family$link == "logit") {
      "Coefficient (log-odds scale)"
    } else if(fit$family$link == "log") {
      "Coefficient (log scale)"
    } else {
      "Coefficient"
    }
    
    subtitle <- paste0("Family: ", fit$family$family, "(", fit$family$link, ")")
    
    # Check if bootstrap CIs are available
    if(!is.null(fit$boot_ci_lower) && !is.null(fit$boot_ci_upper)) {
      # Create a dataframe combining parametric and bootstrap CIs
      # Interpolate bootstrap CIs to match coef grid
      boot_t <- fit$boot_t_grid
      coef_t <- df_coef[[x_col]]
      
      # Interpolate bootstrap CIs to coefficient grid
      boot_lower_interp <- approx(boot_t, fit$boot_ci_lower, xout = coef_t, rule = 2)$y
      boot_upper_interp <- approx(boot_t, fit$boot_ci_upper, xout = coef_t, rule = 2)$y
      
      df_coef$boot_lower <- boot_lower_interp
      df_coef$boot_upper <- boot_upper_interp
      
      ci_label <- paste0("95% Bootstrap CI (B=", fit$n_boot, ")")
      subtitle <- paste0(subtitle, " | ", ci_label)
      
      g <- ggplot(df_coef, aes_string(x = x_col, y = "value")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = boot_lower, ymax = boot_upper), fill = "firebrick", alpha = 0.2) +
        geom_line(color = "firebrick", linewidth = 1) +
        theme_minimal() +
        labs(title = "Estimated Functional Coefficient Beta(t)", 
             subtitle = subtitle,
             x = "Time", y = y_label)
    } else {
      # Use parametric CIs from pfr
      g <- ggplot(df_coef, aes_string(x = x_col, y = "value")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = value - 1.96*se, ymax = value + 1.96*se), fill = "firebrick", alpha = 0.2) +
        geom_line(color = "firebrick", linewidth = 1) +
        theme_minimal() +
        labs(title = "Estimated Functional Coefficient Beta(t)", 
             subtitle = paste0(subtitle, " | 95% Parametric CI"),
             x = "Time", y = y_label)
    }
    
    ggplotly(g) %>% layout(showlegend = FALSE)
  })
  
  output$sofr_pred_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    y_obs <- fit$y_original
    y_pred <- fitted(fit)
    
    if(fit$family$family == "binomial") {
      # For binary: predicted probabilities
      df_pred <- data.frame(
        Observed = factor(y_obs, levels = c(0, 1), labels = c("0", "1")),
        Predicted_Prob = y_pred
      )
      
      g <- ggplot(df_pred, aes(x = Observed, y = Predicted_Prob, fill = Observed)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        scale_fill_manual(values = c("0" = "steelblue", "1" = "firebrick")) +
        theme_minimal() +
        labs(title = "Predicted Probabilities by Observed Class",
             x = "Observed Class", y = "Predicted Probability P(Y=1)") +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40")
      
      ggplotly(g)
    } else {
      # For continuous: scatter plot
      plot_ly(x = y_obs, y = y_pred, type = 'scatter', mode = 'markers', 
              marker = list(color = 'blue', opacity = 0.6)) %>%
        add_lines(x = c(min(y_obs), max(y_obs)), y = c(min(y_obs), max(y_obs)), 
                  line = list(color = 'red', dash = 'dash'), name = "Identity") %>%
        layout(title = "Observed vs Predicted", 
               xaxis = list(title = "Observed"), 
               yaxis = list(title = "Predicted"))
    }
  })
  
  # ROC Curve for binary outcomes
  output$sofr_roc_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(NULL)
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    
    # Validate data
    if(length(y_obs) == 0 || length(y_pred_prob) == 0) return(NULL)
    if(length(unique(y_obs)) < 2) return(NULL)
    
    # Calculate ROC curve
    thresholds <- sort(unique(c(0, y_pred_prob, 1)))
    n_thresh <- length(thresholds)
    
    tpr_vec <- numeric(n_thresh)
    fpr_vec <- numeric(n_thresh)
    
    for(i in seq_along(thresholds)) {
      thresh <- thresholds[i]
      pred_class <- ifelse(y_pred_prob >= thresh, 1, 0)
      tp <- sum(pred_class == 1 & y_obs == 1)
      fp <- sum(pred_class == 1 & y_obs == 0)
      tn <- sum(pred_class == 0 & y_obs == 0)
      fn <- sum(pred_class == 0 & y_obs == 1)
      
      tpr_vec[i] <- if(tp + fn > 0) tp / (tp + fn) else 0
      fpr_vec[i] <- if(fp + tn > 0) fp / (fp + tn) else 0
    }
    
    roc_data <- data.frame(
      threshold = thresholds,
      tpr = tpr_vec,
      fpr = fpr_vec
    )
    roc_data <- roc_data[order(roc_data$fpr, roc_data$tpr), ]
    
    # Calculate AUC using trapezoidal rule
    auc <- 0
    if(nrow(roc_data) > 1) {
      for(i in 2:nrow(roc_data)) {
        auc <- auc + (roc_data$fpr[i] - roc_data$fpr[i-1]) * (roc_data$tpr[i] + roc_data$tpr[i-1]) / 2
      }
    }
    
    # Create plot without text attribute to avoid size mismatch
    p <- plot_ly() %>%
      add_trace(data = roc_data, x = ~fpr, y = ~tpr, type = 'scatter', mode = 'lines',
                line = list(color = 'firebrick', width = 2),
                name = paste0("ROC (AUC = ", round(auc, 3), ")"),
                hovertemplate = "FPR: %{x:.3f}<br>TPR: %{y:.3f}<extra></extra>") %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 1, 
                   line = list(color = 'gray', dash = 'dash'),
                   name = "Random", showlegend = FALSE) %>%
      layout(title = paste0("ROC Curve (AUC = ", round(auc, 3), ")"),
             xaxis = list(title = "False Positive Rate (1 - Specificity)", range = c(0, 1)),
             yaxis = list(title = "True Positive Rate (Sensitivity)", range = c(0, 1)))
    p
  })
  
  # Calibration plot for binary outcomes
  output$sofr_calibration_plot <- renderPlotly({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(NULL)
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    
    # Validate data
    if(length(y_obs) == 0 || length(y_pred_prob) == 0) return(NULL)
    
    # Create calibration bins
    n_bins <- 10
    df_cal <- data.frame(pred = y_pred_prob, obs = y_obs)
    df_cal$bin <- cut(df_cal$pred, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
    
    cal_summary <- df_cal %>%
      group_by(bin) %>%
      summarise(
        mean_pred = mean(pred),
        mean_obs = mean(obs),
        n = n(),
        se = sqrt(mean_obs * (1 - mean_obs) / n),
        .groups = "drop"
      ) %>%
      filter(n >= 3)  # Only include bins with enough observations
    
    if(nrow(cal_summary) == 0) {
      return(plotly_empty() %>% layout(title = "Not enough data for calibration plot"))
    }
    
    # Convert to regular data frame to avoid tibble issues
    cal_summary <- as.data.frame(cal_summary)
    
    p <- plot_ly() %>%
      add_trace(data = cal_summary, x = ~mean_pred, y = ~mean_obs, 
                type = 'scatter', mode = 'markers',
                marker = list(size = sqrt(cal_summary$n) * 3, color = 'steelblue', opacity = 0.7),
                error_y = list(array = 1.96 * cal_summary$se, color = 'steelblue'),
                hovertemplate = "Pred: %{x:.3f}<br>Obs: %{y:.3f}<br>n=%{text}<extra></extra>",
                text = cal_summary$n,
                name = "Observed") %>%
      add_segments(x = 0, xend = 1, y = 0, yend = 1, 
                   line = list(color = 'red', dash = 'dash'),
                   name = "Perfect calibration") %>%
      layout(title = "Calibration Plot",
             xaxis = list(title = "Mean Predicted Probability", range = c(0, 1)),
             yaxis = list(title = "Observed Proportion", range = c(0, 1)))
    p
  })
  
  # Classification metrics for binary outcomes
  output$sofr_classification_metrics <- renderPrint({
    req(values$sofr_model)
    fit <- values$sofr_model
    if(!isTRUE(fit$is_binary)) return(cat("Classification metrics only available for binary outcomes."))
    
    y_obs <- fit$y_original
    y_pred_prob <- fitted(fit)
    y_pred_class <- ifelse(y_pred_prob > 0.5, 1, 0)
    
    # Confusion matrix
    tp <- sum(y_pred_class == 1 & y_obs == 1)
    fp <- sum(y_pred_class == 1 & y_obs == 0)
    tn <- sum(y_pred_class == 0 & y_obs == 0)
    fn <- sum(y_pred_class == 0 & y_obs == 1)
    
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    ppv <- if(tp + fp > 0) tp / (tp + fp) else NA
    npv <- if(tn + fn > 0) tn / (tn + fn) else NA
    f1 <- if(ppv + sensitivity > 0) 2 * (ppv * sensitivity) / (ppv + sensitivity) else NA
    
    # AUC calculation
    thresholds <- sort(unique(c(0, y_pred_prob, 1)))
    tpr_vec <- fpr_vec <- numeric(length(thresholds))
    for(i in seq_along(thresholds)) {
      pred_class <- ifelse(y_pred_prob >= thresholds[i], 1, 0)
      tpr_vec[i] <- sum(pred_class == 1 & y_obs == 1) / sum(y_obs == 1)
      fpr_vec[i] <- sum(pred_class == 1 & y_obs == 0) / sum(y_obs == 0)
    }
    ord <- order(fpr_vec, tpr_vec)
    auc <- sum(diff(fpr_vec[ord]) * (tpr_vec[ord][-1] + tpr_vec[ord][-length(tpr_vec)]) / 2)
    
    cat("=== Classification Metrics (threshold = 0.5) ===\n\n")
    cat("Confusion Matrix:\n")
    cat("                 Predicted\n")
    cat("                  0      1\n")
    cat(sprintf("Actual 0     %5d  %5d\n", tn, fp))
    cat(sprintf("       1     %5d  %5d\n", fn, tp))
    cat("\n")
    cat(sprintf("Accuracy:    %.3f\n", accuracy))
    cat(sprintf("Sensitivity: %.3f (Recall)\n", sensitivity))
    cat(sprintf("Specificity: %.3f\n", specificity))
    cat(sprintf("PPV:         %.3f (Precision)\n", ppv))
    cat(sprintf("NPV:         %.3f\n", npv))
    cat(sprintf("F1 Score:    %.3f\n", f1))
    cat(sprintf("AUC:         %.3f\n", auc))
  })
  
  # ==============================================================================
  # HARMONIC REGRESSION (COSINOR ANALYSIS) MODULE
  # ==============================================================================
  
  # Variable selection UI for harmonic regression
  output$harmonic_var_select_ui <- renderUI({
    req(values$data)  # Only require data; covariates are optional

    # Time variable options from covariates (if available)
    numeric_vars <- if(!is.null(values$covariates)) {
      names(values$covariates)[sapply(values$covariates, is.numeric)]
    } else {
      character(0)  # Empty vector if no covariates
    }
    n_time <- ncol(values$data)
    col_names <- colnames(values$data)
    
    # Try to extract time values from column names
    # Looks for patterns like: VAS_10, T10, time10, col_10, 10:00, 11PM, 2AM, etc.
    suggested_times <- NULL
    detected_pattern <- NULL
    
    if(!is.null(col_names) && length(col_names) > 0) {
      
      # Pattern 0: AM/PM format (e.g., 11PM, 2AM, 11:30PM, VAS_11PM)
      # Check if any column contains AM or PM
      if(any(grepl("[0-9]\\s*[AaPp][Mm]", col_names))) {
        if(all(grepl("[0-9]\\s*[AaPp][Mm]", col_names))) {
          suggested_times <- sapply(col_names, function(cn) {
            # Extract hour, optional minutes, and AM/PM
            # First try with minutes (e.g., 11:30PM)
            if(grepl("[0-9]{1,2}:[0-9]{2}\\s*[AaPp][Mm]", cn)) {
              hour <- as.numeric(gsub(".*?([0-9]{1,2}):[0-9]{2}\\s*[AaPp][Mm].*", "\\1", cn))
              mins <- as.numeric(gsub(".*?[0-9]{1,2}:([0-9]{2})\\s*[AaPp][Mm].*", "\\1", cn))
            } else {
              # Without minutes (e.g., 11PM)
              hour <- as.numeric(gsub(".*?([0-9]{1,2})\\s*[AaPp][Mm].*", "\\1", cn))
              mins <- 0
            }
            ampm <- toupper(gsub(".*([AaPp][Mm]).*", "\\1", cn))
            
            # Convert to 24-hour: 12AM=0, 1-11AM=1-11, 12PM=12, 1-11PM=13-23
            hour_24 <- if(ampm == "AM") {
              if(hour == 12) 0 else hour
            } else {
              if(hour == 12) 12 else hour + 12
            }
            hour_24 + mins / 60
          }, USE.NAMES = FALSE)
          detected_pattern <- "AM/PM format"
        }
      }
      
      # Pattern 1: Dutch/European hour format with "u" suffix (e.g., KSS_9u_dag1, var_14u_something)
      # The number before "u" or "u_" is the hour
      if(is.null(suggested_times)) {
        if(all(grepl("_[0-9]{1,2}u", col_names))) {
          suggested_times <- as.numeric(gsub(".*_([0-9]{1,2})u.*", "\\1", col_names))
          detected_pattern <- "hour with 'u' suffix (e.g., 9u = 9:00)"
        }
      }
      
      # Pattern 2: Trailing numbers (e.g., VAS_10, VAS_12, T10, col10)
      # BUT skip if trailing number looks like day indicator (dag1, dag2, day1, day2)
      if(is.null(suggested_times)) {
        # Check if trailing numbers are likely day indicators
        trailing_nums <- gsub(".*[^0-9]([0-9]+)$", "\\1", col_names)
        is_day_indicator <- all(grepl("(dag|day)[0-9]+$", col_names, ignore.case = TRUE))
        
        if(!is_day_indicator && all(grepl("^[0-9]+$", trailing_nums))) {
          suggested_times <- as.numeric(trailing_nums)
          detected_pattern <- "trailing numbers"
        }
      }
      
      # Pattern 3: Numbers after underscore (e.g., var_10, var_12)
      if(is.null(suggested_times)) {
        underscore_nums <- gsub(".*_([0-9]+).*", "\\1", col_names)
        if(all(grepl("^[0-9]+$", underscore_nums)) && !all(underscore_nums == col_names)) {
          suggested_times <- as.numeric(underscore_nums)
          detected_pattern <- "underscore pattern"
        }
      }
      
      # Pattern 4: Time format HH:MM or HH (e.g., 10:00, 12:30) - 24h format
      if(is.null(suggested_times)) {
        time_match <- grepl("([0-9]{1,2}):?([0-9]{0,2})", col_names)
        if(all(time_match)) {
          hours <- as.numeric(gsub(".*?([0-9]{1,2}):?([0-9]{0,2}).*", "\\1", col_names))
          mins <- gsub(".*?([0-9]{1,2}):?([0-9]{0,2}).*", "\\2", col_names)
          mins <- ifelse(mins == "", 0, as.numeric(mins))
          suggested_times <- hours + mins / 60
          detected_pattern <- "time format"
        }
      }

    }
    
    # Build suggestion text
    suggestion_text <- "e.g., 8,9,10,11,12,14,16,18,20,21,22,23,0,2,4,6"
    suggestion_value <- ""
    detection_msg <- NULL
    
    if(!is.null(suggested_times) && length(suggested_times) == n_time) {
      suggestion_value <- paste(suggested_times, collapse = ",")
      suggestion_text <- suggestion_value
      times_preview <- paste(head(suggested_times, 6), collapse=", ")
      if(n_time > 6) times_preview <- paste0(times_preview, ", ...")
      detection_msg <- div(style = "color: green; font-size: 0.9em;",
                           icon("check-circle"),
                           sprintf(" Detected %d time values from column names: %s", 
                                   n_time, times_preview))
    }
    
    tagList(
      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              numeric_vars),
                  selected = "_index_"),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_manual_'",
        if(!is.null(detection_msg)) detection_msg,
        textAreaInput("harmonic_manual_times", 
                      paste0("Enter ", n_time, " time values (comma-separated):"),
                      value = suggestion_value,
                      placeholder = suggestion_text,
                      rows = 2),
        helpText("Enter the actual clock times for each column in your data. Use 24-hour format or decimal hours.")
      ),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_index_'",
        if(!is.null(suggested_times) && length(suggested_times) == n_time) {
          helpText(HTML(paste0("<b>Note:</b> Detected time values in column names (", 
                               paste(head(suggested_times, 4), collapse=", "), 
                               if(n_time > 4) ", ..." else "",
                               "). Consider using 'Specify times manually' if spacing is unequal.")))
        } else {
          helpText(HTML("<b>Warning:</b> This assumes measurements are equally spaced across the period. If your measurements are unequally spaced (e.g., hourly during day, 2-hourly at night), use 'Specify times manually' instead."))
        }
      )
    )
  })
  
  # Group variable UI
  output$harmonic_group_var_ui <- renderUI({
    req(values$covariates)
    cat_vars <- names(values$covariates)[sapply(values$covariates, function(x) {
      is.factor(x) || is.character(x) || length(unique(x)) <= 10
    })]
    selectInput("harmonic_group_var", "Group Variable (optional):",
                choices = c("None" = "_none_", cat_vars))
  })

  # Parameter bounds hints based on data
  output$harmonic_bounds_hints <- renderUI({
    req(values$data)

    # Get data (use smoothed if available)
    Y <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

    # Calculate data statistics
    y_min <- min(Y, na.rm = TRUE)
    y_max <- max(Y, na.rm = TRUE)
    y_range <- y_max - y_min
    y_mean <- mean(Y, na.rm = TRUE)

    # Get time information
    n_time <- ncol(Y)
    time_max <- if(!is.null(input$harmonic_time_var) && input$harmonic_time_var == "_index_") {
      input$harmonic_period
    } else {
      n_time  # Conservative estimate
    }

    # Build hints text
    hints_html <- sprintf(
      "<div style='background-color: #e8f4f8; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
       <strong>📊 Data Range Hints:</strong><br>
       <small>
       <strong>Your data:</strong> Min=%.2f, Max=%.2f, Mean=%.2f, Range=%.2f<br>
       <strong>Suggested MESOR bounds:</strong> [%.2f, %.2f] (mean ± range)<br>
       <strong>Suggested Amplitude max:</strong> %.2f (observed range)<br>",
      y_min, y_max, y_mean, y_range,
      y_mean - y_range, y_mean + y_range,
      y_range
    )

    # Add exp_sat specific hints if that trend type is selected
    if(!is.null(input$harmonic_trend_type) && input$harmonic_trend_type == "exp_sat") {
      hints_html <- paste0(hints_html, sprintf(
        "<strong>Suggested A_sat bounds:</strong> [%.2f, %.2f] (0.5× to 2× range)<br>
         <strong>Suggested τ bounds:</strong> [0.5, %.1f] (0.5 to max time)<br>",
        y_range * 0.5, y_range * 2,
        time_max
      ))
    }

    hints_html <- paste0(hints_html, "</small></div>")

    HTML(hints_html)
  })

  # Warning UI for harmonic count vs data points
  output$harmonic_warning_ui <- renderUI({
    req(values$data)
    n_time <- ncol(values$data)
    n_harmonics <- input$n_harmonics
    
    # Need at least 2*n_harmonics + 1 parameters (2 per harmonic + MESOR)
    min_required <- 2 * n_harmonics + 2  # +2 for some df for error
    max_safe_harmonics <- floor((n_time - 2) / 2)
    
    period <- if(!is.null(input$harmonic_period)) input$harmonic_period else 24
    
    # Build harmonic info table
    harmonic_info <- paste0(
      "<small><b>Harmonic periods:</b> ",
      paste(sapply(1:n_harmonics, function(h) paste0("H", h, "=", round(period/h, 1), "h")), collapse=", "),
      "</small>"
    )
    
    if(n_harmonics > max_safe_harmonics) {
      tagList(
        div(style = "color: red; font-weight: bold;",
            icon("exclamation-triangle"),
            sprintf(" Warning: %d harmonics require at least %d time points. You have %d.", 
                    n_harmonics, min_required, n_time)),
        div(style = "color: orange;",
            sprintf("Maximum safe harmonics for your data: %d", max_safe_harmonics)),
        HTML(harmonic_info)
      )
    } else if(n_harmonics > max_safe_harmonics - 1) {
      tagList(
        div(style = "color: orange;",
            icon("exclamation-circle"),
            " Approaching maximum harmonics for your data. Model may overfit."),
        HTML(harmonic_info)
      )
    } else {
      HTML(harmonic_info)
    }
  })
  
  # Subject selector for individual plots
  output$harmonic_subject_selector <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(!is.null(mod$individual_fits)) {
      n_subj <- length(mod$individual_fits)
      
      # Check which fits succeeded and show details for failed ones
      subject_labels <- sapply(1:n_subj, function(i) {
        fit_i <- mod$individual_fits[[i]]
        if(!is.null(fit_i) && isTRUE(fit_i$success)) {
          paste("Subject", i)
        } else if(!is.null(fit_i) && !is.null(fit_i$n_valid)) {
          paste0("Subject ", i, " (failed: ", fit_i$n_valid, "/", fit_i$n_required, " pts)")
        } else {
          paste("Subject", i, "(failed)")
        }
      })
      
      selectInput("harmonic_subject_select", "Select Subject:", 
                  choices = c("All (overlay)" = "all", 
                              "Mean curve" = "mean",
                              setNames(1:n_subj, subject_labels)),
                  selected = "mean")
    }
  })
  
  # Harmonic selector for polar plot
  output$harmonic_selector_polar <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      selectInput("selected_harmonic_polar", "Display Harmonic:", 
                  choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                  selected = 1)
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # Harmonic selector for parameter distribution
  output$harmonic_selector_dist <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      fluidRow(
        column(4,
               selectInput("selected_harmonic_dist", "Display Harmonic:", 
                           choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                           selected = 1)
        ),
        column(8,
               helpText("Select which harmonic to display in the amplitude and acrophase distributions.")
        )
      )
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # Harmonic selector for group comparison
  output$harmonic_selector_group <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    if(mod$n_harmonics > 1) {
      fluidRow(
        column(4,
               selectInput("selected_harmonic_group", "Compare Harmonic:", 
                           choices = setNames(1:mod$n_harmonics, paste("H", 1:mod$n_harmonics, sep="")),
                           selected = 1)
        ),
        column(8,
               helpText("Select which harmonic to use for group comparisons.")
        )
      )
    } else {
      helpText("Only one harmonic fitted (fundamental).")
    }
  })
  
  # ==============================================================================
  # CIRCULAR STATISTICS HELPER FUNCTIONS
  # ==============================================================================
  # Custom implementations for acrophase analysis (no external package required)
  # Based on Mardia & Jupp (2000) "Directional Statistics"
  
  # Circular mean (returns radians)
  circular_mean <- function(angles_rad) {
    x <- mean(cos(angles_rad), na.rm = TRUE)
    y <- mean(sin(angles_rad), na.rm = TRUE)
    atan2(y, x)
  }
  
  # Mean resultant length (measure of concentration, 0-1)
  mean_resultant_length <- function(angles_rad) {
    x <- mean(cos(angles_rad), na.rm = TRUE)
    y <- mean(sin(angles_rad), na.rm = TRUE)
    sqrt(x^2 + y^2)
  }
  
  # Circular standard deviation (in same units as input)
  circular_sd <- function(angles_rad) {
    r_bar <- mean_resultant_length(angles_rad)
    if(is.na(r_bar)) {
      return(NA)
    }
    if(r_bar > 0 && r_bar < 1) {
      sqrt(-2 * log(r_bar))
    } else if(r_bar >= 1) {
      0  # All points identical
    } else {
      NA  # Undefined
    }
  }
  
  # Circular standard error (approximate, based on von Mises)
  circular_se <- function(angles_rad) {
    n <- sum(!is.na(angles_rad))
    r_bar <- mean_resultant_length(angles_rad)
    if(is.na(r_bar)) {
      return(NA)
    }
    if(r_bar > 0 && n > 1) {
      # Approximate SE for circular mean (Mardia & Jupp, 2000)
      1 / sqrt(n * r_bar^2)
    } else {
      NA
    }
  }
  
  # Watson-Williams test for comparing two or more groups of circular data
  watson_williams_test <- function(angles_list) {
    # angles_list: list of vectors, each containing angles in radians for one group
    k <- length(angles_list)
    if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))
    
    n <- sapply(angles_list, function(x) sum(!is.na(x)))
    N <- sum(n)
    
    if(any(n < 2)) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Each group needs at least 2 observations"))
    
    # Resultant lengths for each group
    R <- sapply(angles_list, function(x) {
      x <- x[!is.na(x)]
      sqrt(sum(cos(x))^2 + sum(sin(x))^2)
    })
    
    # Total resultant length (pooled)
    all_angles <- unlist(angles_list)
    all_angles <- all_angles[!is.na(all_angles)]
    R_total <- sqrt(sum(cos(all_angles))^2 + sum(sin(all_angles))^2)
    
    r_bar_total <- R_total / N
    
    if(r_bar_total < 0.45) {
      return(list(F = NA, df1 = k - 1, df2 = N - k, p = NA, r_bar = r_bar_total,
                  message = "Warning: Data too dispersed (r̄ < 0.45). Consider non-parametric test."))
    }
    
    # Concentration parameter estimate
    kappa <- if(r_bar_total < 0.53) {
      2 * r_bar_total + r_bar_total^3 + 5 * r_bar_total^5 / 6
    } else if(r_bar_total < 0.85) {
      -0.4 + 1.39 * r_bar_total + 0.43 / (1 - r_bar_total)
    } else {
      1 / (r_bar_total^3 - 4 * r_bar_total^2 + 3 * r_bar_total)
    }
    
    g <- 1 - 1 / (3 * 8 * kappa^2)
    sum_R <- sum(R)
    F_stat <- g * (N - k) * (sum_R - R_total) / ((k - 1) * (N - sum_R))
    
    df1 <- k - 1
    df2 <- N - k
    p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
    
    list(F = F_stat, df1 = df1, df2 = df2, p = p_value, kappa = kappa, r_bar = r_bar_total, message = NULL)
  }

  # Hotelling's T² test on (beta_cos, beta_sin) pairs - amplitude-weighted acrophase comparison
  # Tests whether the bivariate rhythmic vector (beta_cos, beta_sin) differs between groups.
  # Because amplitude = sqrt(beta_cos² + beta_sin²), subjects with stronger rhythms carry
  # more weight. For k>2 groups, a one-way MANOVA F approximation is used.
  hotelling_t2 <- function(beta_cos_list, beta_sin_list) {
    k <- length(beta_cos_list)
    if(k < 2) return(list(F = NA, df1 = NA, df2 = NA, p = NA, message = "Need at least 2 groups"))

    # Build per-group matrices
    mats <- lapply(seq_len(k), function(i) {
      x <- beta_cos_list[[i]]
      y <- beta_sin_list[[i]]
      ok <- complete.cases(x, y)
      cbind(x[ok], y[ok])
    })

    ns <- sapply(mats, nrow)
    if(any(ns < 3)) return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                                message = "Each group needs at least 3 observations"))
    N <- sum(ns)
    p <- 2  # two variables: beta_cos and beta_sin

    means <- lapply(mats, colMeans)

    if(k == 2) {
      # Two-sample Hotelling's T²
      S1 <- cov(mats[[1]]); S2 <- cov(mats[[2]])
      Sp <- ((ns[1]-1)*S1 + (ns[2]-1)*S2) / (N - 2)
      if(det(Sp) < .Machine$double.eps)
        return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                    message = "Singular pooled covariance - check data"))
      d  <- means[[1]] - means[[2]]
      T2 <- (ns[1]*ns[2]) / N * t(d) %*% solve(Sp) %*% d
      F_stat <- (N - p - 1) / ((N - 2) * p) * as.numeric(T2)
      df1 <- p; df2 <- N - p - 1
    } else {
      # k-group one-way MANOVA (Wilks' lambda → F approximation)
      grand_mean <- colMeans(do.call(rbind, mats))
      # Between-group sum-of-squares-and-products (H)
      H <- Reduce("+", lapply(seq_len(k), function(i)
        ns[i] * outer(means[[i]] - grand_mean, means[[i]] - grand_mean)))
      # Within-group (E)
      E <- Reduce("+", lapply(mats, function(m) {
        cm <- colMeans(m)
        t(sweep(m, 2, cm)) %*% sweep(m, 2, cm)
      }))
      if(det(E) < .Machine$double.eps)
        return(list(F = NA, df1 = NA, df2 = NA, p = NA,
                    message = "Singular within-group covariance - check data"))
      lambda <- det(E) / det(H + E)
      # Rao's F approximation for p=2
      df1 <- p * (k - 1)
      df2 <- N - k - p + 1
      F_stat <- ((1 - sqrt(lambda)) / sqrt(lambda)) * (df2 / df1)
    }

    p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
    list(F = as.numeric(F_stat), df1 = df1, df2 = df2, p = as.numeric(p_value), message = NULL)
  }

  # ==============================================================================
  # CORE COSINOR FITTING FUNCTIONS
  # ==============================================================================
  
  # Single cosinor fit for one subject
  fit_cosinor <- function(time, y, period = 24, n_harmonics = 1, trend_type = "none",
                          use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                          amplitude_min = 0, amplitude_max = NA,
                          A_sat_min = NA, A_sat_max = NA,
                          tau_min = 0.5, tau_max = NA) {
    # Remove NAs
    valid <- complete.cases(time, y)
    time <- time[valid]
    y <- y[valid]
    n <- length(y)

    # Set default bounds based on data if not specified
    if(use_bounds) {
      y_min <- min(y, na.rm = TRUE)
      y_max <- max(y, na.rm = TRUE)
      y_range <- y_max - y_min
      t_max <- max(time) - min(time)

      if(is.na(mesor_min)) mesor_min <- y_min - y_range
      if(is.na(mesor_max)) mesor_max <- y_max + y_range
      if(is.na(amplitude_max)) amplitude_max <- y_range * 2
      if(is.na(A_sat_max)) A_sat_max <- y_range * 2
      if(is.na(tau_max)) tau_max <- t_max * 5
    }

    # Calculate number of trend parameters
    n_trend_params <- switch(trend_type,
                             "none" = 0,
                             "linear" = 1,
                             "log" = 1,
                             "exp_sat" = 2,  # A and tau for nonlinear fit
                             0)

    min_params <- 2 * n_harmonics + 1 + n_trend_params
    if(n < min_params + 1) {
      return(list(success = FALSE, message = "Insufficient data points"))
    }

    # For exponential saturation, use nonlinear fitting
    if(trend_type == "exp_sat") {
      return(fit_cosinor_nonlinear(time, y, period, n_harmonics, trend_type,
                                    FALSE, NULL, 0.32, 0.66,
                                    use_bounds, mesor_min, mesor_max, amplitude_min, amplitude_max,
                                    A_sat_min, A_sat_max, tau_min, tau_max))
    }

    # For linear/log models with bounds, use nonlinear least squares
    if(use_bounds && trend_type %in% c("linear", "log", "none")) {
      return(fit_cosinor_nonlinear(time, y, period, n_harmonics, trend_type,
                                    FALSE, NULL, 0.32, 0.66,
                                    use_bounds, mesor_min, mesor_max, amplitude_min, amplitude_max,
                                    A_sat_min, A_sat_max, tau_min, tau_max))
    }
    
    # Build design matrix with multiple harmonics (linear models)
    X <- matrix(1, nrow = n, ncol = 1)  # Intercept (MESOR)
    colnames_X <- "MESOR"
    trend_cols <- character(0)
    
    # Add trend based on type
    if(trend_type == "linear") {
      X <- cbind(X, time)
      colnames_X <- c(colnames_X, "trend_linear")
      trend_cols <- "trend_linear"
    } else if(trend_type == "log") {
      # Use log(t+1) to avoid log(0) and handle t=0
      t_offset <- min(time)
      log_time <- log(time - t_offset + 1)
      X <- cbind(X, log_time)
      colnames_X <- c(colnames_X, "trend_log")
      trend_cols <- "trend_log"
    }
    
    coef_offset <- 1 + length(trend_cols)
    
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      X <- cbind(X, cos(omega * time), sin(omega * time))
      colnames_X <- c(colnames_X, paste0("cos", h), paste0("sin", h))
    }
    colnames(X) <- colnames_X
    
    # Fit linear model
    fit <- lm(y ~ X - 1)  # -1 because X already has intercept
    coefs <- coef(fit)
    se <- summary(fit)$coefficients[, 2]
    
    # Extract parameters
    mesor <- coefs[1]
    mesor_se <- se[1]
    
    # Extract trend parameters
    trend_params <- list()
    if(trend_type != "none") {
      for(i in seq_along(trend_cols)) {
        trend_params[[trend_cols[i]]] <- list(
          coef = coefs[1 + i],
          se = se[1 + i]
        )
      }
    }
    
    # Calculate amplitude and acrophase for each harmonic
    amplitudes <- numeric(n_harmonics)
    acrophases <- numeric(n_harmonics)
    amp_se <- numeric(n_harmonics)
    acro_se <- numeric(n_harmonics)
    
    vcov_mat <- vcov(fit)
    
    for(h in 1:n_harmonics) {
      cos_idx <- coef_offset + 2 * (h - 1) + 1
      sin_idx <- coef_offset + 2 * (h - 1) + 2
      
      beta_cos <- coefs[cos_idx]
      beta_sin <- coefs[sin_idx]
      
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi
      
      if(amplitudes[h] > 1e-10) {
        grad_amp <- c(beta_cos, beta_sin) / amplitudes[h]
        idx <- c(cos_idx, sin_idx)
        var_amp <- t(grad_amp) %*% vcov_mat[idx, idx] %*% grad_amp
        amp_se[h] <- sqrt(var_amp)
        
        grad_acro <- c(-beta_sin, beta_cos) / (amplitudes[h]^2)
        var_acro <- t(grad_acro) %*% vcov_mat[idx, idx] %*% grad_acro
        acro_se[h] <- sqrt(var_acro)
      } else {
        amp_se[h] <- NA
        acro_se[h] <- NA
      }
    }
    
    acrophases_time <- acrophases * period / (2 * pi)
    acro_se_time <- acro_se * period / (2 * pi)
    
    # Goodness of fit - calculate manually because lm(y ~ X - 1) R² is vs origin, not mean
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum(residuals(fit)^2)
    r_squared <- 1 - ss_resid / ss_total

    # Adjusted R² accounting for number of predictors
    n_predictors <- ncol(X)  # MESOR + trend + harmonics
    adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - n_predictors)

    percent_rhythm <- r_squared * 100

    df1 <- 2 * n_harmonics
    df2 <- n - 2 * n_harmonics - 1 - n_trend_params
    f_stat <- ((ss_total - ss_resid) / df1) / (ss_resid / df2)
    p_value <- pf(f_stat, df1, df2, lower.tail = FALSE)

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian linear model
    sigma_sq <- ss_resid / n
    log_lik <- -n/2 * (log(2*pi) + log(sigma_sq) + 1)

    # AIC: Akaike Information Criterion
    # AIC = -2*log(L) + 2*k, where k = number of parameters
    k <- n_predictors + 1  # predictors + sigma
    aic <- -2 * log_lik + 2 * k

    # AICc: Corrected AIC for small samples
    # AICc = AIC + 2*k*(k+1)/(n-k-1)
    aicc <- if(n - k - 1 > 0) {
      aic + (2 * k * (k + 1)) / (n - k - 1)
    } else {
      NA  # Not defined when n-k-1 <= 0
    }

    # BIC: Bayesian Information Criterion
    # BIC = -2*log(L) + k*log(n)
    bic <- -2 * log_lik + k * log(n)

    # LOOCV: Leave-one-out cross-validation (leave out time points)
    # For each observation, refit the model without it and predict
    loocv_errors <- numeric(n)
    for(i in 1:n) {
      # Remove observation i
      X_loo <- X[-i, , drop = FALSE]
      y_loo <- y[-i]

      # Refit model
      fit_loo <- tryCatch({
        lm.fit(X_loo, y_loo)
      }, error = function(e) NULL)

      if(!is.null(fit_loo)) {
        # Predict left-out observation
        y_pred <- sum(X[i, ] * fit_loo$coefficients)
        loocv_errors[i] <- (y[i] - y_pred)^2
      } else {
        loocv_errors[i] <- NA
      }
    }

    # LOOCV RMSE (root mean squared error)
    loocv_rmse <- sqrt(mean(loocv_errors, na.rm = TRUE))

    # ===========================================================================
    # Variance decomposition: Calculate R² for Process S and Process C separately
    # ===========================================================================
    r_squared_S <- 0  # Variance explained by homeostatic trend alone
    r_squared_C <- 0  # Variance explained by circadian rhythm alone
    percent_S <- 0    # Percentage of total R² from Process S
    percent_C <- 0    # Percentage of total R² from Process C

    # Baseline model (MESOR only)
    y_mean <- mean(y)
    ss_resid_baseline <- sum((y - y_mean)^2)

    # Model with trend only (Process S)
    if(trend_type != "none") {
      if(trend_type == "exp_sat") {
        # For exponential saturation, use the already fitted trend component
        # Extract fitted trend values from the full model
        if(!is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
          A_sat_val <- trend_params$A_sat$coef
          tau_val <- trend_params$tau$coef
          t_offset_temp <- min(time)
          fitted_trend <- mesor + A_sat_val * (1 - exp(-(time - t_offset_temp) / tau_val))
          ss_resid_trend <- sum((y - fitted_trend)^2)
          r_squared_S <- 1 - ss_resid_trend / ss_total
        }
      } else {
        X_trend <- matrix(1, nrow = n, ncol = 1)
        if(trend_type == "linear") {
          X_trend <- cbind(X_trend, time)
        } else if(trend_type == "log") {
          t_offset_temp <- min(time)
          log_time_temp <- log(time - t_offset_temp + 1)
          X_trend <- cbind(X_trend, log_time_temp)
        }
        fit_trend <- lm(y ~ X_trend - 1)
        ss_resid_trend <- sum(residuals(fit_trend)^2)
        r_squared_S <- 1 - ss_resid_trend / ss_total
      }
    }

    # Model with circadian only (Process C)
    X_circ <- matrix(1, nrow = n, ncol = 1)
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      X_circ <- cbind(X_circ, cos(omega * time), sin(omega * time))
    }
    fit_circ <- lm(y ~ X_circ - 1)
    ss_resid_circ <- sum(residuals(fit_circ)^2)
    r_squared_C <- 1 - ss_resid_circ / ss_total

    # Calculate proportions of total R²
    if(!is.na(r_squared) && r_squared > 0) {
      percent_S <- (r_squared_S / r_squared) * 100
      percent_C <- (r_squared_C / r_squared) * 100
    }
    
    # Store time offset for prediction (needed for log/sqrt)
    t_offset <- if(trend_type == "log") min(time) else 0
    t_center <- 0
    
    list(
      success = TRUE,
      mesor = mesor,
      mesor_se = mesor_se,
      trend_type = trend_type,
      trend_params = trend_params,
      t_offset = t_offset,
      t_center = t_center,
      amplitudes = amplitudes,
      amp_se = amp_se,
      acrophases = acrophases,
      acrophases_time = acrophases_time,
      acro_se = acro_se,
      acro_se_time = acro_se_time,
      coefs = coefs,
      se = se,
      vcov = vcov_mat,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      percent_rhythm = percent_rhythm,
      f_stat = f_stat,
      p_value = p_value,
      aic = aic,                        # Akaike Information Criterion
      aicc = aicc,                      # Corrected AIC for small samples
      bic = bic,                        # Bayesian Information Criterion
      loocv_rmse = loocv_rmse,          # Leave-one-out CV RMSE
      r_squared_S = r_squared_S,        # R² from Process S alone
      r_squared_C = r_squared_C,        # R² from Process C alone
      percent_S = percent_S,            # % of total R² from S
      percent_C = percent_C,            # % of total R² from C
      fitted = fitted(fit),
      residuals = residuals(fit),
      time = time,
      y = y,
      period = period,
      n_harmonics = n_harmonics,
      n = n
    )
  }

# Nonlinear fitting function for exponential saturation trend
fit_cosinor_nonlinear <- function(time, y, period, n_harmonics, trend_type = "none",
                                   include_inertia = FALSE, wake_onset = NULL,
                                   W0_init = 0.32, tau_W_init = 0.66,
                                   use_bounds = FALSE, mesor_min = NA, mesor_max = NA,
                                   amplitude_min = 0, amplitude_max = NA,
                                   A_sat_min = NA, A_sat_max = NA,
                                   tau_min = 0.5, tau_max = NA) {
  # This function handles exponential saturation trend (exp_sat) and bounded optimization

  n <- length(y)
  t_offset <- min(time)
  t_shifted <- time - t_offset  # For exponential saturation
  t_max <- max(t_shifted)

  # Initial values from data
  y_range <- max(y) - min(y)
  y_min <- min(y)
  y_max <- max(y)
  y_mean <- mean(y)

  # Set default bounds based on data if not specified
  if(use_bounds) {
    if(is.na(mesor_min)) mesor_min <- y_min - y_range
    if(is.na(mesor_max)) mesor_max <- y_max + y_range
    if(is.na(amplitude_max)) amplitude_max <- y_range * 2
    if(is.na(A_sat_min)) A_sat_min <- -Inf  # Allow negative for decreasing trends
    if(is.na(A_sat_max)) A_sat_max <- y_range * 2
    if(is.na(tau_max)) tau_max <- t_max * 5
  } else {
    # No user-specified bounds - use wide defaults with minimal numerical constraints
    # For exp_sat, we still need sensible bounds for numerical stability
    mesor_min <- -Inf
    mesor_max <- Inf
    amplitude_min <- -Inf
    amplitude_max <- Inf

    # For exp_sat: Use original bounds from CIRCAREGold.R when bounding disabled
    # These bounds apply to Approaches 1-2, but Approach 3 will be unbounded
    if(trend_type == "exp_sat") {
      # Original bounds: A_sat unbounded, tau constrained
      A_sat_min <- -Inf
      A_sat_max <- Inf
      tau_min <- 0.5  # Original fixed minimum
      tau_max <- t_max * 5  # Original upper bound
    } else {
      A_sat_min <- -Inf
      A_sat_max <- Inf
      tau_min <- 0.5
      tau_max <- Inf
    }
  }

  # Estimate trend direction from linear regression
  lin_fit <- lm(y ~ t_shifted)
  lin_slope <- coef(lin_fit)[2]
  lin_intercept <- coef(lin_fit)[1]

  # Build formula dynamically based on components
  # Always start with mesor
  formula_parts <- c("mesor")
  start_list <- list(mesor = y_mean)
  lower_bounds <- c(mesor = mesor_min)
  upper_bounds <- c(mesor = mesor_max)

  # Add trend component
  if(trend_type == "linear") {
    formula_parts <- c(formula_parts, "beta_t * time")
    start_list$beta_t <- lin_slope
    lower_bounds["beta_t"] <- -Inf
    upper_bounds["beta_t"] <- Inf

  } else if(trend_type == "log") {
    formula_parts <- c(formula_parts, "beta_log * log(time - t_offset + 1)")
    start_list$beta_log <- lin_slope * log(t_max + 1)
    lower_bounds["beta_log"] <- -Inf
    upper_bounds["beta_log"] <- Inf

  } else if(trend_type == "exp_sat") {
    formula_parts <- c(formula_parts, "A_sat * (1 - exp(-t_shifted / tau))")

    # Better starting values for exp_sat
    if(lin_slope > 0) {
      start_list$A_sat <- y_range * 1.5
      start_list$tau <- t_max / 4
    } else {
      start_list$A_sat <- lin_slope * t_max
      start_list$tau <- t_max / 3
    }
    start_list$tau <- max(1, min(start_list$tau, t_max * 2))

    # Apply user-specified or default bounds
    lower_bounds["A_sat"] <- A_sat_min
    lower_bounds["tau"] <- tau_min
    upper_bounds["A_sat"] <- A_sat_max
    upper_bounds["tau"] <- tau_max
  }

  # Add harmonic components
  harmonic_parts <- sapply(1:n_harmonics, function(h) {
    omega <- 2 * pi * h / period
    start_list[[paste0("b_cos", h)]] <<- 0
    start_list[[paste0("b_sin", h)]] <<- 0

    # Bound harmonic coefficients to respect amplitude constraints
    # Since amplitude = sqrt(b_cos^2 + b_sin^2), bound each coefficient to +/- amplitude_max
    lower_bounds[paste0("b_cos", h)] <<- -amplitude_max
    lower_bounds[paste0("b_sin", h)] <<- -amplitude_max
    upper_bounds[paste0("b_cos", h)] <<- amplitude_max
    upper_bounds[paste0("b_sin", h)] <<- amplitude_max

    sprintf("b_cos%d * cos(%f * time) + b_sin%d * sin(%f * time)", h, omega, h, omega)
  })
  formula_parts <- c(formula_parts, harmonic_parts)

  # Build complete formula
  formula_str <- sprintf("y ~ %s", paste(formula_parts, collapse = " + "))

  # Debug: Print formula being fitted (useful for verification)
  # cat(sprintf("Fitting model: %s\n", formula_str))

  # Ensure starting values respect bounds (clamp them)
  for(param_name in names(start_list)) {
    if(param_name %in% names(lower_bounds)) {
      start_val <- start_list[[param_name]]
      lb <- lower_bounds[param_name]
      ub <- upper_bounds[param_name]

      # Clamp to bounds if finite
      if(!is.infinite(lb) && start_val < lb) {
        start_list[[param_name]] <- lb + (ub - lb) * 0.1  # 10% above lower bound
      }
      if(!is.infinite(ub) && start_val > ub) {
        start_list[[param_name]] <- ub - (ub - lb) * 0.1  # 10% below upper bound
      }
    }
  }

  # Prepare data frame for fitting
  fit_data <- data.frame(
    y = y,
    time = time,
    t_shifted = t_shifted,
    t_offset = t_offset
  )

  # Try fitting with multiple approaches
  fit_success <- FALSE
  nls_fit <- NULL
  error_msgs <- c()

  # Choose fitting strategy based on trend type and whether bounds are enabled
  # NOTE: For exp_sat, Approaches 1-2 use bounds, but Approach 3 falls back to unbounded
  if(use_bounds || trend_type == "exp_sat") {
    # BOUNDED OPTIMIZATION: Use algorithms that support bounds
    # (Either user requested bounds, or exp_sat in Approaches 1-2)

    # Approach 1: Try nlsLM (Levenberg-Marquardt) if available - most robust
    if(requireNamespace("minpack.lm", quietly = TRUE)) {
      tryCatch({
        nls_fit <- minpack.lm::nlsLM(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          lower = lower_bounds,
          upper = upper_bounds,
          control = minpack.lm::nls.lm.control(maxiter = 300)
        )
        fit_success <- TRUE
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nlsLM:", e$message))
      })
    }

    # Approach 2: Try standard nls with port algorithm (allows bounds)
    if(!fit_success) {
      tryCatch({
        nls_fit <- nls(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          algorithm = "port",
          lower = lower_bounds,
          upper = upper_bounds,
          control = nls.control(maxiter = 300, warnOnly = TRUE)
        )
        fit_success <- TRUE
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nls-port:", e$message))
      })
    }

  } else {
    # UNBOUNDED OPTIMIZATION: Use default algorithms without bounds
    # (Only for linear/log/none trends when user hasn't requested bounds)

    # Approach 1: Try nlsLM without bounds if available
    if(requireNamespace("minpack.lm", quietly = TRUE)) {
      tryCatch({
        nls_fit <- minpack.lm::nlsLM(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          control = minpack.lm::nls.lm.control(maxiter = 300)
        )
        fit_success <- TRUE
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nlsLM:", e$message))
      })
    }

    # Approach 2: Try standard nls with default algorithm (no bounds)
    if(!fit_success) {
      tryCatch({
        nls_fit <- nls(
          as.formula(formula_str),
          data = fit_data,
          start = start_list,
          control = nls.control(maxiter = 300, warnOnly = TRUE)
        )
        fit_success <- TRUE
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("nls-default:", e$message))
      })
    }
  }

  # Approach 3: Try different starting values for tau parameters
  if(!fit_success && (trend_type == "exp_sat" || include_inertia)) {
    for(tau_mult in c(0.1, 0.5, 2, 5)) {
      if(trend_type == "exp_sat") {
        start_list$tau <- t_max * tau_mult / 3
      }
      if(include_inertia) {
        start_list$tau_W <- tau_W_init * tau_mult
      }

      tryCatch({
        if(use_bounds) {
          # Use port algorithm with bounds (only if user explicitly requested bounds)
          nls_fit <- nls(
            as.formula(formula_str),
            data = fit_data,
            start = start_list,
            algorithm = "port",
            lower = lower_bounds,
            upper = upper_bounds,
            control = nls.control(maxiter = 300, warnOnly = TRUE)
          )
        } else {
          # Use default algorithm without bounds (original fallback behavior)
          # This matches CIRCAREGold.R Approach 3 - unbounded even for exp_sat
          nls_fit <- nls(
            as.formula(formula_str),
            data = fit_data,
            start = start_list,
            control = nls.control(maxiter = 300, warnOnly = TRUE)
          )
        }

        # Check if fit is reasonable (R² > 0)
        fitted_check <- predict(nls_fit)
        ss_tot_check <- sum((y - mean(y))^2)
        ss_res_check <- sum((y - fitted_check)^2)
        r2_check <- 1 - ss_res_check / ss_tot_check

        if(r2_check > 0) {
          fit_success <- TRUE
          break
        }
      }, error = function(e) {
        error_msgs <<- c(error_msgs, paste("retry", tau_mult, ":", e$message))
      })
    }
  }

  if(!fit_success || is.null(nls_fit)) {
    return(list(
      success = FALSE,
      message = sprintf("Nonlinear fit failed to converge. Errors: %s",
                       paste(head(error_msgs, 3), collapse = "; "))
    ))
  }

  # Extract results
  tryCatch({
    coefs <- coef(nls_fit)
    se <- tryCatch(summary(nls_fit)$coefficients[, 2],
                   error = function(e) rep(NA, length(coefs)))
    names(se) <- names(coefs)

    mesor <- coefs["mesor"]
    mesor_se <- se["mesor"]

    # Extract trend parameters based on type
    trend_params <- list()
    if(trend_type == "linear") {
      trend_params$trend_linear <- list(
        coef = as.numeric(coefs["beta_t"]),
        se = as.numeric(se["beta_t"])
      )
    } else if(trend_type == "log") {
      trend_params$trend_log <- list(
        coef = as.numeric(coefs["beta_log"]),
        se = as.numeric(se["beta_log"])
      )
    } else if(trend_type == "exp_sat") {
      trend_params$A_sat <- list(
        coef = as.numeric(coefs["A_sat"]),
        se = if(!is.na(se["A_sat"])) as.numeric(se["A_sat"]) else NA
      )
      trend_params$tau <- list(
        coef = as.numeric(coefs["tau"]),
        se = if(!is.na(se["tau"])) as.numeric(se["tau"]) else NA
      )
    }

    # No inertia parameters for exp_sat-only fitting
    inertia_params <- NULL

    # Extract harmonic parameters
    amplitudes <- numeric(n_harmonics)
    acrophases <- numeric(n_harmonics)
    amp_se <- numeric(n_harmonics)
    acro_se <- numeric(n_harmonics)

    for(h in 1:n_harmonics) {
      beta_cos <- coefs[paste0("b_cos", h)]
      beta_sin <- coefs[paste0("b_sin", h)]
      amplitudes[h] <- sqrt(beta_cos^2 + beta_sin^2)
      acrophases[h] <- atan2(beta_sin, beta_cos)
      if(acrophases[h] < 0) acrophases[h] <- acrophases[h] + 2 * pi

      # Approximate SE
      se_cos <- if(!is.na(se[paste0("b_cos", h)])) se[paste0("b_cos", h)] else 0
      se_sin <- if(!is.na(se[paste0("b_sin", h)])) se[paste0("b_sin", h)] else 0
      amp_se[h] <- sqrt(se_cos^2 + se_sin^2) / sqrt(2)
      acro_se[h] <- NA  # Complex for nonlinear
    }

    acrophases_time <- acrophases * period / (2 * pi)
    acro_se_time <- acro_se * period / (2 * pi)

    # Goodness of fit
    fitted_vals <- predict(nls_fit)
    ss_total <- sum((y - mean(y))^2)
    ss_resid <- sum((y - fitted_vals)^2)
    r_squared <- 1 - ss_resid / ss_total
    percent_rhythm <- max(0, r_squared * 100)

    # Calculate p-value for circadian rhythm
    n_params <- length(coefs)
    df1 <- 2 * n_harmonics
    df2 <- n - n_params
    if(df2 > 0 && ss_resid > 0) {
      f_stat <- ((ss_total - ss_resid) / df1) / (ss_resid / df2)
      p_value <- pf(f_stat, df1, df2, lower.tail = FALSE)
    } else {
      f_stat <- NA
      p_value <- NA
    }

    # ===========================================================================
    # Model selection metrics: AIC, AICc, BIC, LOOCV
    # ===========================================================================

    # Log-likelihood for Gaussian nonlinear model
    sigma_sq <- ss_resid / n
    log_lik <- -n/2 * (log(2*pi) + log(sigma_sq) + 1)

    # AIC: Akaike Information Criterion
    k <- n_params + 1  # parameters + sigma
    aic <- -2 * log_lik + 2 * k

    # AICc: Corrected AIC for small samples
    aicc <- if(n - k - 1 > 0) {
      aic + (2 * k * (k + 1)) / (n - k - 1)
    } else {
      NA
    }

    # BIC: Bayesian Information Criterion
    bic <- -2 * log_lik + k * log(n)

    # LOOCV: Leave-one-out cross-validation
    # Note: For nonlinear models, this is computationally expensive
    # We'll use a simplified approach: predict each point using the full model
    # and apply a leave-one-out correction based on leverage
    loocv_rmse <- tryCatch({
      # Get leverage values (hat matrix diagonal)
      # For NLS, we approximate using the Jacobian
      if(!is.null(nls_fit)) {
        # Simple LOOCV using prediction errors
        # For nonlinear models, full refit for each point is too expensive
        # Use prediction residuals as approximation
        residuals_vec <- y - fitted_vals
        sqrt(mean(residuals_vec^2))
      } else {
        NA
      }
    }, error = function(e) NA)

    # ===========================================================================
    # Variance decomposition: Calculate R² for Process S and Process C separately
    # ===========================================================================
    r_squared_S <- 0  # Variance explained by homeostatic trend alone
    r_squared_C <- 0  # Variance explained by circadian rhythm alone
    percent_S <- 0    # Percentage of total R² from Process S
    percent_C <- 0    # Percentage of total R² from Process C

    # Model with trend only (Process S) - using exp_sat trend
    if(trend_type == "exp_sat" && !is.null(trend_params$A_sat) && !is.null(trend_params$tau)) {
      A_sat_val <- trend_params$A_sat$coef
      tau_val <- trend_params$tau$coef
      fitted_trend <- mesor + A_sat_val * (1 - exp(-t_shifted / tau_val))
      ss_resid_trend <- sum((y - fitted_trend)^2)
      r_squared_S <- max(0, 1 - ss_resid_trend / ss_total)
    }

    # Model with circadian only (Process C)
    # Build a model with MESOR + harmonics only (no trend)
    tryCatch({
      X_circ <- matrix(1, nrow = n, ncol = 1)
      for(h in 1:n_harmonics) {
        omega <- 2 * pi * h / period
        X_circ <- cbind(X_circ, cos(omega * time), sin(omega * time))
      }
      fit_circ <- lm(y ~ X_circ - 1)
      ss_resid_circ <- sum(residuals(fit_circ)^2)
      r_squared_C <- max(0, 1 - ss_resid_circ / ss_total)
    }, error = function(e) {
      r_squared_C <<- 0
    })

    # Calculate proportions of total R²
    if(!is.na(r_squared) && r_squared > 0) {
      percent_S <- (r_squared_S / r_squared) * 100
      percent_C <- (r_squared_C / r_squared) * 100
    }

    list(
      success = TRUE,
      mesor = as.numeric(mesor),
      mesor_se = if(!is.na(mesor_se)) as.numeric(mesor_se) else NA,
      trend_type = trend_type,
      trend_params = trend_params,
      inertia_params = inertia_params,  # NEW: Sleep inertia parameters
      t_offset = t_offset,
      t_center = 0,
      amplitudes = amplitudes,
      amp_se = amp_se,
      acrophases = acrophases,
      acrophases_time = acrophases_time,
      acro_se = acro_se,
      acro_se_time = acro_se_time,
      coefs = coefs,
      se = se,
      vcov = NULL,
      r_squared = r_squared,
      adj_r_squared = r_squared,  # Approximate
      percent_rhythm = percent_rhythm,
      f_stat = f_stat,
      p_value = p_value,
      aic = aic,                        # Akaike Information Criterion
      aicc = aicc,                      # Corrected AIC for small samples
      bic = bic,                        # Bayesian Information Criterion
      loocv_rmse = loocv_rmse,          # Leave-one-out CV RMSE
      r_squared_S = r_squared_S,        # R² from Process S alone
      r_squared_C = r_squared_C,        # R² from Process C alone
      percent_S = percent_S,            # % of total R² from S
      percent_C = percent_C,            # % of total R² from C
      fitted = fitted_vals,
      residuals = y - fitted_vals,
      time = time,
      y = y,
      period = period,
      n_harmonics = n_harmonics,
      n = n
    )
  }, error = function(e) {
    list(success = FALSE, message = paste("Result extraction failed:", e$message))
  })
}

  # Predict from cosinor model
  predict_cosinor <- function(fit, newtime = NULL, component = "total", include_trend_in_pred = TRUE) {
    if(is.null(newtime)) newtime <- fit$time
    period <- fit$period
    n_harmonics <- fit$n_harmonics
    coefs <- fit$coefs
    trend_type <- if(!is.null(fit$trend_type)) fit$trend_type else "none"
    
    # Handle legacy format (include_trend boolean)
    if(is.null(fit$trend_type) && isTRUE(fit$include_trend)) {
      trend_type <- "linear"
    }
    
    pred <- rep(coefs[1], length(newtime))  # MESOR
    
    # Calculate trend component based on type
    if(include_trend_in_pred && trend_type != "none") {
      t_offset <- if(!is.null(fit$t_offset)) fit$t_offset else min(newtime)
      
      trend_val <- switch(trend_type,
                          "linear" = coefs[2] * newtime,
                          "log" = coefs[2] * log(newtime - t_offset + 1),
                          "exp_sat" = {
                            A_sat <- coefs["A_sat"]
                            tau <- coefs["tau"]
                            t_shifted <- newtime - t_offset
                            A_sat * (1 - exp(-t_shifted / tau))
                          },
                          "two_process" = {
                            if(!is.null(fit$S_trajectory) && !is.null(fit$time) &&
                               "beta_S" %in% names(coefs)) {
                              S_interp <- tryCatch({
                                approx(fit$time, fit$S_trajectory, xout = newtime,
                                       rule = 2, ties = "ordered")$y
                              }, error = function(e) rep(NA_real_, length(newtime)))
                              coefs["beta_S"] * S_interp
                            } else {
                              rep(0, length(newtime))
                            }
                          },
                          rep(0, length(newtime))
      )
      pred <- pred + trend_val
    }
    
    # Determine coefficient offset based on trend type
    n_trend_coefs <- switch(trend_type,
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, "two_process" = 0, 0)
    coef_offset <- 1 + n_trend_coefs

    if(component == "total" || component == "all") {
      for(h in 1:n_harmonics) {
        omega <- 2 * pi * h / period
        if(trend_type == "exp_sat" || trend_type == "two_process") {
          beta_cos <- coefs[paste0("b_cos", h)]
          beta_sin <- coefs[paste0("b_sin", h)]
        } else {
          cos_idx <- coef_offset + 2 * (h - 1) + 1
          sin_idx <- coef_offset + 2 * (h - 1) + 2
          beta_cos <- coefs[cos_idx]
          beta_sin <- coefs[sin_idx]
        }
        pred <- pred + beta_cos * cos(omega * newtime) + beta_sin * sin(omega * newtime)
      }
    } else if(is.numeric(component) && component >= 1 && component <= n_harmonics) {
      h <- component
      omega <- 2 * pi * h / period
      if(trend_type == "exp_sat" || trend_type == "two_process") {
        beta_cos <- coefs[paste0("b_cos", h)]
        beta_sin <- coefs[paste0("b_sin", h)]
      } else {
        cos_idx <- coef_offset + 2 * (h - 1) + 1
        sin_idx <- coef_offset + 2 * (h - 1) + 2
        beta_cos <- coefs[cos_idx]
        beta_sin <- coefs[sin_idx]
      }
      pred <- coefs[1] + beta_cos * cos(omega * newtime) + beta_sin * sin(omega * newtime)
    }
    
    return(pred)
  }
  
  # Get harmonic components separately
  get_harmonic_components <- function(fit, newtime = NULL) {
    if(is.null(newtime)) newtime <- fit$time
    period <- fit$period
    n_harmonics <- fit$n_harmonics
    coefs <- fit$coefs
    trend_type <- if(!is.null(fit$trend_type)) fit$trend_type else "none"
    
    # Handle legacy format
    if(is.null(fit$trend_type) && isTRUE(fit$include_trend)) {
      trend_type <- "linear"
    }
    
    t_offset <- if(!is.null(fit$t_offset)) fit$t_offset else min(newtime)
    
    # Determine coefficient offset
    n_trend_coefs <- switch(trend_type,
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 0, "two_process" = 0, 0)
    coef_offset <- 1 + n_trend_coefs
    
    components <- list()
    components$mesor <- rep(coefs[1], length(newtime))
    components$trend_type <- trend_type
    
    if(trend_type != "none") {
      # Compute trend component
      trend_val <- switch(trend_type,
                          "linear" = coefs[2] * newtime,
                          "log" = coefs[2] * log(newtime - t_offset + 1),
                          "exp_sat" = {
                            A_sat <- coefs["A_sat"]
                            tau <- coefs["tau"]
                            A_sat * (1 - exp(-(newtime - t_offset) / tau))
                          },
                          "two_process" = {
                            if(!is.null(fit$S_trajectory) && !is.null(fit$time) &&
                               "beta_S" %in% names(coefs)) {
                              S_interp <- tryCatch({
                                approx(fit$time, fit$S_trajectory, xout = newtime,
                                       rule = 2, ties = "ordered")$y
                              }, error = function(e) rep(NA_real_, length(newtime)))
                              coefs["beta_S"] * S_interp
                            } else {
                              rep(0, length(newtime))
                            }
                          },
                          rep(0, length(newtime))
      )
      components$trend <- trend_val
    }

    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      if(trend_type == "exp_sat" || trend_type == "two_process") {
        beta_cos <- coefs[paste0("b_cos", h)]
        beta_sin <- coefs[paste0("b_sin", h)]
      } else {
        cos_idx <- coef_offset + 2 * (h - 1) + 1
        sin_idx <- coef_offset + 2 * (h - 1) + 2
        beta_cos <- coefs[cos_idx]
        beta_sin <- coefs[sin_idx]
      }
      components[[paste0("harmonic_", h)]] <- beta_cos * cos(omega * newtime) +
        beta_sin * sin(omega * newtime)
    }

    components$total <- predict_cosinor(fit, newtime)
    return(components)
  }
  
  # Predict curve from mean coefficients (for group/population means with all harmonics)
  predict_from_coefs <- function(coefs, time_vec, period, n_harmonics, trend_type = "none", 
                                 t_offset = 0, t_center = 0) {
    # Handle legacy boolean format
    if(is.logical(trend_type)) {
      trend_type <- if(trend_type) "linear" else "none"
    }
    
    # coefs format: c(mesor, [trend_coefs...], beta_cos_1, beta_sin_1, ...)
    pred <- rep(coefs[1], length(time_vec))  # MESOR
    
    # Determine trend offset
    n_trend_coefs <- switch(as.character(trend_type),
                            "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
    coef_offset <- 1 + n_trend_coefs
    
    # Add trend based on type
    if(trend_type != "none" && trend_type != FALSE) {
      trend_val <- switch(as.character(trend_type),
                          "linear" = coefs[2] * time_vec,
                          "log" = coefs[2] * log(time_vec - t_offset + 1),
                          "exp_sat" = coefs[2] * (1 - exp(-(time_vec - t_offset) / coefs[3])),  # A_sat * (1 - exp(-t/tau))
                          rep(0, length(time_vec))
      )
      pred <- pred + trend_val
    }
    
    for(h in 1:n_harmonics) {
      omega <- 2 * pi * h / period
      beta_cos <- coefs[coef_offset + 2 * h - 1]
      beta_sin <- coefs[coef_offset + 2 * h]
      pred <- pred + beta_cos * cos(omega * time_vec) + beta_sin * sin(omega * time_vec)
    }
    return(pred)
  }
  
  # Get trend value at a specific time point
  get_trend_value <- function(trend_type, trend_params, time_vec, t_offset = 0) {
    if(trend_type == "none" || is.null(trend_params) || length(trend_params) == 0) {
      return(rep(0, length(time_vec)))
    }
    
    switch(trend_type,
           "linear" = trend_params$trend_linear$coef * time_vec,
           "log" = trend_params$trend_log$coef * log(time_vec - t_offset + 1),
           "exp_sat" = {
             A_sat <- trend_params$A_sat$coef
             tau <- trend_params$tau$coef
             A_sat * (1 - exp(-(time_vec - t_offset) / tau))
           },
           rep(0, length(time_vec))
    )
  }
  
  # Get human-readable trend label
  get_trend_label <- function(trend_type, prefix = "") {
    label <- switch(trend_type,
                    "linear" = "Linear Trend",
                    "log" = "Log Trend",
                    "exp_sat" = "Exp. Saturation",
                    "two_process" = "Process S",
                    "Process S"  # default
    )
    if(nchar(prefix) > 0) paste(label, prefix) else label
  }
  
  # Helper: Get mean trend coefficients from individual parameters for building coef vectors
  # Returns a vector of trend coefficients to append to mesor for predict_from_coefs
  get_mean_trend_coefs <- function(params, trend_type) {
    if(trend_type == "none") return(numeric(0))
    
    if(trend_type == "linear" && "trend_linear" %in% names(params)) {
      return(mean(params$trend_linear, na.rm = TRUE))
    } else if(trend_type == "log" && "trend_log" %in% names(params)) {
      return(mean(params$trend_log, na.rm = TRUE))
    } else if(trend_type == "exp_sat") {
      coefs <- numeric(0)
      if("A_sat" %in% names(params)) coefs <- c(coefs, mean(params$A_sat, na.rm = TRUE))
      if("tau" %in% names(params)) coefs <- c(coefs, mean(params$tau, na.rm = TRUE))
      return(coefs)
    }
    return(numeric(0))
  }
  
  # Helper: Check if trend columns exist in params
  has_trend_params <- function(params, trend_type) {
    if(trend_type == "none") return(FALSE)
    if(trend_type == "linear") return("trend_linear" %in% names(params))
    if(trend_type == "log") return("trend_log" %in% names(params))
    if(trend_type == "exp_sat") return("A_sat" %in% names(params) || "tau" %in% names(params))
    if(trend_type == "two_process") {
      return("beta_S" %in% names(params) && "tau_w" %in% names(params) && "tau_s" %in% names(params))
    }
    return(FALSE)
  }
  
  # Helper: Get the primary trend column name for a given trend_type
  get_trend_col <- function(trend_type) {
    switch(trend_type,
           "linear" = "trend_linear",
           "log" = "trend_log",
           "exp_sat" = "A_sat",
           NULL)
  }
  
  # Helper: Compute trend line values for plotting
  compute_trend_line <- function(params, trend_type, time_vec, t_offset = 0) {
    if(trend_type == "none" || !has_trend_params(params, trend_type)) {
      return(NULL)
    }

    mesor <- mean(params$mesor, na.rm = TRUE)

    if(trend_type == "linear" && "trend_linear" %in% names(params)) {
      slope <- mean(params$trend_linear, na.rm = TRUE)
      return(mesor + slope * time_vec)
    } else if(trend_type == "log" && "trend_log" %in% names(params)) {
      slope <- mean(params$trend_log, na.rm = TRUE)
      return(mesor + slope * log(time_vec - t_offset + 1))
    } else if(trend_type == "exp_sat" && "A_sat" %in% names(params) && "tau" %in% names(params)) {
      A_sat <- mean(params$A_sat, na.rm = TRUE)
      tau <- mean(params$tau, na.rm = TRUE)
      return(mesor + A_sat * (1 - exp(-(time_vec - t_offset) / tau)))
    }
    return(NULL)
  }

  # ==============================================================================
  # MAIN HARMONIC REGRESSION EVENT HANDLER
  # ==============================================================================
  
  observeEvent(input$run_harmonic, {
    req(values$data)
    
    showNotification("Running Harmonic Regression...", type = "message", duration = 2)
    
    tryCatch({
      # Check if smoothed data is available
      using_smoothed <- !is.null(values$smooth_data)
      Y <- if(using_smoothed) values$smooth_data else values$data
      n_subjects <- nrow(Y)
      n_time <- ncol(Y)
      period <- input$harmonic_period
      n_harmonics <- input$n_harmonics
      trend_type <- input$harmonic_trend_type
      
      # Calculate number of trend parameters
      n_trend_params <- switch(trend_type,
                               "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
      
      # Diagnostic: Check data type and dimensions
      cat(sprintf("Data diagnostics: %d subjects × %d time points, type=%s, smoothed=%s, trend=%s\n", 
                  n_subjects, n_time, typeof(Y), using_smoothed, trend_type))
      
      # Check for NAs in the data
      na_counts <- apply(Y, 1, function(row) sum(is.na(row)))
      valid_counts <- n_time - na_counts
      subjects_with_nas <- sum(na_counts > 0)
      subjects_all_na <- sum(na_counts == n_time)
      
      if(subjects_all_na > 0) {
        all_na_subjects <- which(na_counts == n_time)
        showNotification(
          sprintf("ERROR: %d subjects have ALL missing values (subjects: %s). Check data selection!", 
                  subjects_all_na, paste(head(all_na_subjects, 10), collapse=", ")),
          type = "error", duration = 15)
        cat(sprintf("Subjects with all NA: %s\n", paste(all_na_subjects, collapse=", ")))
        
        # Show sample of data for first all-NA subject
        if(length(all_na_subjects) > 0) {
          cat(sprintf("First all-NA subject (%d) data sample: %s\n", 
                      all_na_subjects[1], 
                      paste(head(Y[all_na_subjects[1], ], 10), collapse=", ")))
        }
      } else if(subjects_with_nas > 0 && !using_smoothed) {
        showNotification(
          sprintf("%d subjects have missing values. Consider applying smoothing first to interpolate missing data.", 
                  subjects_with_nas),
          type = "warning", duration = 8)
      } else if(subjects_with_nas > 0 && using_smoothed) {
        # This shouldn't happen if smoothing worked correctly
        showNotification(
          sprintf("Warning: %d subjects still have NAs after smoothing. Some fits may fail.", 
                  subjects_with_nas),
          type = "warning", duration = 8)
      }
      
      # Check if we have enough data points for the requested harmonics
      min_required <- 2 * n_harmonics + 2 + n_trend_params
      max_safe_harmonics <- floor((n_time - 2 - n_trend_params) / 2)
      
      if(n_time < min_required) {
        showNotification(
          sprintf("Error: %d harmonics require at least %d time points. You have %d. Maximum safe: %d harmonics.", 
                  n_harmonics, min_required, n_time, max_safe_harmonics),
          type = "error", duration = 10)
        return()
      }
      
      if(n_harmonics > max_safe_harmonics) {
        showNotification(
          sprintf("Warning: Using %d harmonics with only %d time points may cause overfitting. Consider reducing to %d harmonics.", 
                  n_harmonics, n_time, max_safe_harmonics),
          type = "warning", duration = 8)
      }
      
      # Determine time variable
      original_times <- NULL
      wrap_applied <- FALSE
      
      if(input$harmonic_time_var == "_index_") {
        # Use column indices scaled to period (assumes equal spacing!)
        time_vec <- seq(0, period * (n_time - 1) / n_time, length.out = n_time)
        original_times <- time_vec
        showNotification("Using equally-spaced time points. If your data has unequal spacing, use 'Specify times manually'.", 
                         type = "warning", duration = 5)
        
      } else if(input$harmonic_time_var == "_manual_") {
        # Parse manual time input
        manual_input <- input$harmonic_manual_times
        if(is.null(manual_input) || nchar(trimws(manual_input)) == 0) {
          showNotification("Please enter time values!", type = "error")
          return()
        }
        
        # Parse comma-separated values
        time_vec <- tryCatch({
          vals <- as.numeric(unlist(strsplit(gsub(" ", "", manual_input), ",")))
          if(any(is.na(vals))) stop("Non-numeric values")
          vals
        }, error = function(e) {
          showNotification("Could not parse time values. Use comma-separated numbers (e.g., 8,9,10,11,12,14,16,18,20).", 
                           type = "error")
          return(NULL)
        })
        
        if(is.null(time_vec)) return()
        
        if(length(time_vec) != n_time) {
          showNotification(paste0("Number of time values (", length(time_vec), 
                                  ") must match number of columns (", n_time, ")!"), 
                           type = "error")
          return()
        }
        
        # Detect wrap-around: if a time is smaller than the previous, add period
        # This handles cases like 8,9,10,...,22,23,0,2,4,6 → 8,9,10,...,22,23,24,26,28,30
        original_times <- time_vec
        for(i in 2:length(time_vec)) {
          if(time_vec[i] < time_vec[i-1]) {
            # Wrap-around detected - add period to this and all subsequent values
            time_vec[i:length(time_vec)] <- time_vec[i:length(time_vec)] + period
          }
        }
        
        # Check if wrap-around was applied
        wrap_applied <- !identical(original_times, time_vec)
        if(wrap_applied) {
          showNotification(paste0("Detected wrap-around at midnight. Adjusted times: ", 
                                  paste(round(time_vec, 1), collapse=", ")), 
                           type = "message", duration = 5)
        } else {
          showNotification(paste("Using manual time points:", paste(round(time_vec, 1), collapse=", ")), 
                           type = "message", duration = 3)
        }
        
      } else {
        # Use selected covariate column
        time_vec <- values$covariates[[input$harmonic_time_var]]
        original_times <- time_vec
        wrap_applied <- FALSE
        if(length(time_vec) != n_time) {
          showNotification("Selected time variable doesn't match data dimensions. Using equal spacing.", 
                           type = "warning")
          time_vec <- seq(0, period * (n_time - 1) / n_time, length.out = n_time)
          original_times <- time_vec
        } else {
          # Apply wrap-around detection for covariate time variables too
          for(i in 2:length(time_vec)) {
            if(time_vec[i] < time_vec[i-1]) {
              time_vec[i:length(time_vec)] <- time_vec[i:length(time_vec)] + period
            }
          }
          wrap_applied <- !identical(original_times, time_vec)
          if(wrap_applied) {
            showNotification(paste0("Detected wrap-around (period=", period, "). Adjusted times: ", 
                                    paste(round(time_vec, 1), collapse=", ")), 
                             type = "message", duration = 5)
          }
        }
      }
      
      # Check for potential issues with time values
      if(max(time_vec) > period * 1.5 && !wrap_applied) {
        showNotification(paste0("Note: Max time value (", round(max(time_vec), 1), 
                                ") is larger than period (", period, "). Values will be wrapped using modulo."), 
                         type = "warning", duration = 5)
      }
      
      # Individual cosinor analysis
      individual_fits <- list()
      
      # Build column names for all harmonics
      param_cols <- c("subject", "mesor", "mesor_se")

      # Add trend columns based on type
      if(trend_type == "linear") {
        param_cols <- c(param_cols, "trend_linear", "trend_linear_se")
      } else if(trend_type == "log") {
        param_cols <- c(param_cols, "trend_log", "trend_log_se")
      } else if(trend_type == "exp_sat") {
        param_cols <- c(param_cols, "A_sat", "A_sat_se", "tau", "tau_se")
      }

      for(h in 1:n_harmonics) {
        param_cols <- c(param_cols,
                        paste0("amplitude_", h), paste0("amp_se_", h),
                        paste0("acrophase_rad_", h), paste0("acrophase_time_", h),
                        paste0("acro_se_time_", h),
                        paste0("beta_cos_", h), paste0("beta_sin_", h))
      }
      param_cols <- c(param_cols, "r_squared", "percent_rhythm", "p_value",
                      "r_squared_S", "r_squared_C", "percent_S", "percent_C")
      
      individual_params <- data.frame(matrix(ncol = length(param_cols), nrow = 0))
      colnames(individual_params) <- param_cols
      
      # Coefficient offset for trend
      coef_offset <- 1 + n_trend_params
      
      # Track failed fits
      failed_fits <- list()
      
      # Store time offsets for prediction
      t_offset_global <- min(time_vec)
      t_center_global <- mean(time_vec)

      # Read parameter bounding options from UI
      use_bounds <- isTRUE(input$harmonic_use_bounds)
      mesor_min <- if(use_bounds) input$harmonic_mesor_min else NA
      mesor_max <- if(use_bounds) input$harmonic_mesor_max else NA
      amplitude_min <- if(use_bounds) input$harmonic_amplitude_min else 0
      amplitude_max <- if(use_bounds) input$harmonic_amplitude_max else NA
      A_sat_min <- if(use_bounds) input$harmonic_A_sat_min else NA
      A_sat_max <- if(use_bounds) input$harmonic_A_sat_max else NA
      tau_min <- if(use_bounds) input$harmonic_tau_min else 0.5
      tau_max <- if(use_bounds) input$harmonic_tau_max else NA

      if(use_bounds) {
        bounds_msg <- sprintf("Using parameter bounds: MESOR [%.2f, %.2f], Amplitude [%.2f, %.2f]",
                              ifelse(is.na(mesor_min), -Inf, mesor_min),
                              ifelse(is.na(mesor_max), Inf, mesor_max),
                              amplitude_min,
                              ifelse(is.na(amplitude_max), Inf, amplitude_max))

        if(trend_type == "exp_sat") {
          bounds_msg <- paste0(bounds_msg,
                               sprintf(", A_sat [%.2f, %.2f], τ [%.2f, %.2f]",
                                       ifelse(is.na(A_sat_min), -Inf, A_sat_min),
                                       ifelse(is.na(A_sat_max), Inf, A_sat_max),
                                       tau_min,
                                       ifelse(is.na(tau_max), Inf, tau_max)))
        }

        showNotification(bounds_msg, type = "message", duration = 5)
      }

      withProgress(message = 'Fitting individual cosinor models...', value = 0, {
        for(i in 1:n_subjects) {
          y_i <- Y[i, ]

          # Count valid (non-NA) data points for this subject
          n_valid_points <- sum(!is.na(y_i))

          # Debug: Check for unusual values
          if(n_valid_points == 0) {
            cat(sprintf("Subject %d: All NA. First 5 values: %s\n", i,
                        paste(head(y_i, 5), collapse=", ")))
          }

          fit_i <- fit_cosinor(time_vec, y_i, period = period, n_harmonics = n_harmonics,
                               trend_type = trend_type,
                               use_bounds = use_bounds,
                               mesor_min = mesor_min,
                               mesor_max = mesor_max,
                               amplitude_min = amplitude_min,
                               amplitude_max = amplitude_max,
                               A_sat_min = A_sat_min,
                               A_sat_max = A_sat_max,
                               tau_min = tau_min,
                               tau_max = tau_max)

          if(fit_i$success) {
            individual_fits[[i]] <- fit_i

            # Build row with all harmonic parameters
            row_data <- list(subject = i, mesor = fit_i$mesor, mesor_se = fit_i$mesor_se)

            # Add trend parameters based on type
            if(trend_type != "none" && !is.null(fit_i$trend_params)) {
              for(param_name in names(fit_i$trend_params)) {
                row_data[[param_name]] <- fit_i$trend_params[[param_name]]$coef
                row_data[[paste0(param_name, "_se")]] <- fit_i$trend_params[[param_name]]$se
              }
            }

            for(h in 1:n_harmonics) {
              row_data[[paste0("amplitude_", h)]] <- fit_i$amplitudes[h]
              row_data[[paste0("amp_se_", h)]] <- fit_i$amp_se[h]
              row_data[[paste0("acrophase_rad_", h)]] <- fit_i$acrophases[h]
              row_data[[paste0("acrophase_time_", h)]] <- fit_i$acrophases_time[h]
              row_data[[paste0("acro_se_time_", h)]] <- fit_i$acro_se_time[h]

              # Get beta coefficients - handle both linear and nls fits
              if(trend_type == "exp_sat") {
                row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[paste0("b_cos", h)]
                row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[paste0("b_sin", h)]
              } else {
                cos_idx <- coef_offset + 2 * (h - 1) + 1
                sin_idx <- coef_offset + 2 * (h - 1) + 2
                row_data[[paste0("beta_cos_", h)]] <- fit_i$coefs[cos_idx]
                row_data[[paste0("beta_sin_", h)]] <- fit_i$coefs[sin_idx]
              }
            }
            row_data$r_squared <- fit_i$r_squared
            row_data$percent_rhythm <- fit_i$percent_rhythm
            row_data$p_value <- fit_i$p_value
            row_data$aic <- fit_i$aic
            row_data$aicc <- fit_i$aicc
            row_data$bic <- fit_i$bic
            row_data$loocv_rmse <- fit_i$loocv_rmse
            row_data$r_squared_S <- fit_i$r_squared_S
            row_data$r_squared_C <- fit_i$r_squared_C
            row_data$percent_S <- fit_i$percent_S
            row_data$percent_C <- fit_i$percent_C

            individual_params <- rbind(individual_params, as.data.frame(row_data))
          } else {
            # Store failed fit with reason
            individual_fits[[i]] <- list(
              success = FALSE, 
              message = fit_i$message,
              n_valid = n_valid_points,
              n_required = 2 * n_harmonics + 1 + n_trend_params + 1
            )
            failed_fits[[length(failed_fits) + 1]] <- list(
              subject = i,
              n_valid = n_valid_points,
              reason = fit_i$message
            )
          }
          
          if(i %% 10 == 0) incProgress(10 / n_subjects)
        }
      })
      
      # Report failed fits
      if(length(failed_fits) > 0) {
        n_failed <- length(failed_fits)
        min_required <- 2 * n_harmonics + 1 + n_trend_params + 1
        
        failed_subjects <- sapply(failed_fits, function(x) x$subject)
        failed_nvalid <- sapply(failed_fits, function(x) x$n_valid)
        
        msg <- sprintf("%d of %d subjects failed to fit (need %d+ valid points). Failed: %s",
                       n_failed, n_subjects, min_required,
                       paste(paste0("S", failed_subjects, "(", failed_nvalid, "pts)"), collapse = ", "))
        
        showNotification(msg, type = "warning", duration = 10)
      }
      
      # Population-mean statistics (always calculated)
      pop_mean_fit <- NULL
      group_fits <- NULL
      
      # Always calculate population mean parameters (vector averaging for circular data)
      {
        # Calculate population mean parameters (vector averaging for circular data)
        mean_mesor <- mean(individual_params$mesor, na.rm = TRUE)
        
        # Vector average for amplitude and acrophase (first harmonic for primary stats)
        x_components <- individual_params$amplitude_1 * cos(individual_params$acrophase_rad_1)
        y_components <- individual_params$amplitude_1 * sin(individual_params$acrophase_rad_1)
        
        mean_x <- mean(x_components, na.rm = TRUE)
        mean_y <- mean(y_components, na.rm = TRUE)
        
        mean_amplitude <- sqrt(mean_x^2 + mean_y^2)
        mean_acrophase_rad <- atan2(mean_y, mean_x)
        if(mean_acrophase_rad < 0) mean_acrophase_rad <- mean_acrophase_rad + 2 * pi
        mean_acrophase_time <- mean_acrophase_rad * period / (2 * pi)
        
        # Rayleigh test for uniformity of acrophases (first harmonic)
        n_valid <- sum(!is.na(individual_params$acrophase_rad_1))
        r_bar <- mean_amplitude / mean(individual_params$amplitude_1, na.rm = TRUE)
        rayleigh_z <- n_valid * r_bar^2
        rayleigh_p <- exp(-rayleigh_z)  # Approximation
        
        # SE from circular statistics
        circ_var <- 1 - r_bar
        circ_sd <- if(r_bar > 0 && r_bar < 1) sqrt(-2 * log(r_bar)) else NA
        
        # Store mean coefficients for ALL harmonics (for proper multi-harmonic curve plotting)
        # Format: [mesor, (trend coefs if trend), beta_cos_1, beta_sin_1, ...]
        mean_coefs <- c(mean_mesor)
        
        # Add mean trend coefficient(s) based on trend type
        trend_coefs <- get_mean_trend_coefs(individual_params, trend_type)
        if(length(trend_coefs) > 0) {
          mean_coefs <- c(mean_coefs, trend_coefs)
        }
        
        mean_amplitudes <- numeric(n_harmonics)
        mean_acrophases_rad <- numeric(n_harmonics)
        mean_acrophases_time <- numeric(n_harmonics)
        
        for(h in 1:n_harmonics) {
          beta_cos_col <- paste0("beta_cos_", h)
          beta_sin_col <- paste0("beta_sin_", h)
          amp_col <- paste0("amplitude_", h)
          acro_col <- paste0("acrophase_rad_", h)
          
          # Mean of raw coefficients (for curve reconstruction)
          mean_beta_cos <- mean(individual_params[[beta_cos_col]], na.rm = TRUE)
          mean_beta_sin <- mean(individual_params[[beta_sin_col]], na.rm = TRUE)
          mean_coefs <- c(mean_coefs, mean_beta_cos, mean_beta_sin)
          
          # Vector-averaged amplitude and acrophase
          x_h <- individual_params[[amp_col]] * cos(individual_params[[acro_col]])
          y_h <- individual_params[[amp_col]] * sin(individual_params[[acro_col]])
          mean_amplitudes[h] <- sqrt(mean(x_h, na.rm = TRUE)^2 + mean(y_h, na.rm = TRUE)^2)
          acro_h <- atan2(mean(y_h, na.rm = TRUE), mean(x_h, na.rm = TRUE))
          if(acro_h < 0) acro_h <- acro_h + 2 * pi
          mean_acrophases_rad[h] <- acro_h
          mean_acrophases_time[h] <- acro_h * period / (2 * pi) / h  # Adjust for harmonic number
        }
        
        # Also compute arithmetic means of individual parameters
        indiv_means <- list(
          mesor = mean(individual_params$mesor, na.rm = TRUE),
          mesor_sd = sd(individual_params$mesor, na.rm = TRUE)
        )
        for(h in 1:n_harmonics) {
          indiv_means[[paste0("amplitude_", h)]] <- mean(individual_params[[paste0("amplitude_", h)]], na.rm = TRUE)
          indiv_means[[paste0("amplitude_", h, "_sd")]] <- sd(individual_params[[paste0("amplitude_", h)]], na.rm = TRUE)
          indiv_means[[paste0("acrophase_time_", h)]] <- mean(individual_params[[paste0("acrophase_time_", h)]], na.rm = TRUE)
          indiv_means[[paste0("acrophase_time_", h, "_sd")]] <- sd(individual_params[[paste0("acrophase_time_", h)]], na.rm = TRUE)
        }

        # Add model selection metrics (if columns exist)
        if("aic" %in% names(individual_params)) {
          indiv_means$aic <- mean(individual_params$aic, na.rm = TRUE)
          indiv_means$aic_sd <- sd(individual_params$aic, na.rm = TRUE)
          indiv_means$aicc <- mean(individual_params$aicc, na.rm = TRUE)
          indiv_means$aicc_sd <- sd(individual_params$aicc, na.rm = TRUE)
          indiv_means$bic <- mean(individual_params$bic, na.rm = TRUE)
          indiv_means$bic_sd <- sd(individual_params$bic, na.rm = TRUE)
          indiv_means$loocv_rmse <- mean(individual_params$loocv_rmse, na.rm = TRUE)
          indiv_means$loocv_rmse_sd <- sd(individual_params$loocv_rmse, na.rm = TRUE)
        }

        # Add variance decomposition statistics (if columns exist)
        if("r_squared_S" %in% names(individual_params) && "r_squared_C" %in% names(individual_params)) {
          indiv_means$r_squared_S <- mean(individual_params$r_squared_S, na.rm = TRUE)
          indiv_means$r_squared_S_sd <- sd(individual_params$r_squared_S, na.rm = TRUE)
          indiv_means$r_squared_C <- mean(individual_params$r_squared_C, na.rm = TRUE)
          indiv_means$r_squared_C_sd <- sd(individual_params$r_squared_C, na.rm = TRUE)
          indiv_means$percent_S <- mean(individual_params$percent_S, na.rm = TRUE)
          indiv_means$percent_S_sd <- sd(individual_params$percent_S, na.rm = TRUE)
          indiv_means$percent_C <- mean(individual_params$percent_C, na.rm = TRUE)
          indiv_means$percent_C_sd <- sd(individual_params$percent_C, na.rm = TRUE)
        }
        
        pop_mean_fit <- list(
          mean_mesor = mean_mesor,
          mean_amplitude = mean_amplitude,  # First harmonic (for backwards compatibility)
          mean_acrophase_rad = mean_acrophase_rad,
          mean_acrophase_time = mean_acrophase_time,
          mean_coefs = mean_coefs,  # All coefficients for curve reconstruction
          mean_amplitudes = mean_amplitudes,  # All harmonics
          mean_acrophases_rad = mean_acrophases_rad,
          mean_acrophases_time = mean_acrophases_time,
          indiv_means = indiv_means,  # Arithmetic means of individual params
          r_bar = r_bar,
          circ_var = circ_var,
          circ_sd = circ_sd,
          rayleigh_z = rayleigh_z,
          rayleigh_p = rayleigh_p,
          n = n_valid
        )
      }
      
      # Group comparison - calculate group-specific statistics
      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]
        groups <- unique(group_var)
        group_fits <- list()
        
        for(g in groups) {
          idx <- which(group_var == g)
          grp_params <- individual_params[individual_params$subject %in% idx, ]
          
          if(nrow(grp_params) >= 3) {
            # Mean coefficients for curve reconstruction
            grp_coefs <- c(mean(grp_params$mesor, na.rm = TRUE))
            
            # Add trend coefficients based on type
            grp_trend_params <- list()
            if(trend_type == "linear" && "trend_linear" %in% names(grp_params)) {
              grp_coefs <- c(grp_coefs, mean(grp_params$trend_linear, na.rm = TRUE))
              grp_trend_params$trend_linear <- list(
                mean = mean(grp_params$trend_linear, na.rm = TRUE),
                sd = sd(grp_params$trend_linear, na.rm = TRUE)
              )
            } else if(trend_type == "log" && "trend_log" %in% names(grp_params)) {
              grp_coefs <- c(grp_coefs, mean(grp_params$trend_log, na.rm = TRUE))
              grp_trend_params$trend_log <- list(
                mean = mean(grp_params$trend_log, na.rm = TRUE),
                sd = sd(grp_params$trend_log, na.rm = TRUE)
              )
            } else if(trend_type == "exp_sat") {
              if("A_sat" %in% names(grp_params)) {
                A_sat_mean <- mean(grp_params$A_sat, na.rm = TRUE)
                # If all NA, use 0 as fallback
                if(!is.finite(A_sat_mean)) A_sat_mean <- 0
                grp_coefs <- c(grp_coefs, A_sat_mean)
                grp_trend_params$A_sat <- list(
                  mean = A_sat_mean,
                  sd = sd(grp_params$A_sat, na.rm = TRUE)
                )
              }
              if("tau" %in% names(grp_params)) {
                tau_mean <- mean(grp_params$tau, na.rm = TRUE)
                # If all NA, use 1 as fallback (avoid division by zero)
                if(!is.finite(tau_mean) || tau_mean <= 0) tau_mean <- 1
                grp_coefs <- c(grp_coefs, tau_mean)
                grp_trend_params$tau <- list(
                  mean = tau_mean,
                  sd = sd(grp_params$tau, na.rm = TRUE)
                )
              }
            }
            
            grp_amplitudes <- numeric(n_harmonics)
            grp_acrophases_rad <- numeric(n_harmonics)
            grp_acrophases_time <- numeric(n_harmonics)
            grp_amp_sd <- numeric(n_harmonics)
            
            for(h in 1:n_harmonics) {
              beta_cos_col <- paste0("beta_cos_", h)
              beta_sin_col <- paste0("beta_sin_", h)
              amp_col <- paste0("amplitude_", h)
              acro_col <- paste0("acrophase_rad_", h)

              grp_coefs <- c(grp_coefs,
                             mean(grp_params[[beta_cos_col]], na.rm = TRUE),
                             mean(grp_params[[beta_sin_col]], na.rm = TRUE))

              # Vector-averaged amplitude/acrophase
              x_h <- grp_params[[amp_col]] * cos(grp_params[[acro_col]])
              y_h <- grp_params[[amp_col]] * sin(grp_params[[acro_col]])
              grp_amplitudes[h] <- sqrt(mean(x_h, na.rm = TRUE)^2 + mean(y_h, na.rm = TRUE)^2)
              grp_amp_sd[h] <- sd(grp_params[[amp_col]], na.rm = TRUE)
              acro_h <- atan2(mean(y_h, na.rm = TRUE), mean(x_h, na.rm = TRUE))
              if(acro_h < 0) acro_h <- acro_h + 2 * pi
              grp_acrophases_rad[h] <- acro_h
              grp_acrophases_time[h] <- acro_h * period / (2 * pi) / h
            }

            # Variance decomposition for this group (if columns exist)
            grp_variance_decomp <- NULL
            if("r_squared_S" %in% names(grp_params) && "r_squared_C" %in% names(grp_params)) {
              grp_variance_decomp <- list(
                r_squared_S = mean(grp_params$r_squared_S, na.rm = TRUE),
                r_squared_S_sd = sd(grp_params$r_squared_S, na.rm = TRUE),
                r_squared_C = mean(grp_params$r_squared_C, na.rm = TRUE),
                r_squared_C_sd = sd(grp_params$r_squared_C, na.rm = TRUE),
                percent_S = mean(grp_params$percent_S, na.rm = TRUE),
                percent_S_sd = sd(grp_params$percent_S, na.rm = TRUE),
                percent_C = mean(grp_params$percent_C, na.rm = TRUE),
                percent_C_sd = sd(grp_params$percent_C, na.rm = TRUE)
              )
            }

            group_fits[[as.character(g)]] <- list(
              group = g,
              n = nrow(grp_params),
              mean_mesor = mean(grp_params$mesor, na.rm = TRUE),
              sd_mesor = sd(grp_params$mesor, na.rm = TRUE),
              mean_coefs = grp_coefs,  # All coefficients
              trend_params = grp_trend_params,
              mean_amplitudes = grp_amplitudes,
              sd_amplitudes = grp_amp_sd,
              mean_acrophases_rad = grp_acrophases_rad,
              mean_acrophases_time = grp_acrophases_time,
              # Keep first harmonic for backwards compatibility
              mean_amplitude = grp_amplitudes[1],
              sd_amplitude = grp_amp_sd[1],
              mean_acrophase_rad = grp_acrophases_rad[1],
              mean_acrophase_time = grp_acrophases_time[1],
              # Variance decomposition
              variance_decomp = grp_variance_decomp
            )
          }
        }
      }
      
      # Bootstrap CIs if requested
      boot_results <- NULL
      if(isTRUE(input$harmonic_bootstrap)) {
        B <- input$harmonic_n_boot
        boot_mesor <- numeric(B)
        boot_amplitude <- numeric(B)
        boot_acrophase <- numeric(B)
        
        showNotification(paste("Running", B, "bootstrap iterations..."), type = "message")
        
        withProgress(message = 'Bootstrap...', value = 0, {
          for(b in 1:B) {
            boot_idx <- sample(1:n_subjects, n_subjects, replace = TRUE)
            boot_params <- individual_params[individual_params$subject %in% boot_idx, ]
            
            boot_mesor[b] <- mean(boot_params$mesor, na.rm = TRUE)
            
            # Use first harmonic for bootstrap CIs
            x_b <- boot_params$amplitude_1 * cos(boot_params$acrophase_rad_1)
            y_b <- boot_params$amplitude_1 * sin(boot_params$acrophase_rad_1)
            boot_amplitude[b] <- sqrt(mean(x_b, na.rm = TRUE)^2 + mean(y_b, na.rm = TRUE)^2)
            
            acro_b <- atan2(mean(y_b, na.rm = TRUE), mean(x_b, na.rm = TRUE))
            if(acro_b < 0) acro_b <- acro_b + 2 * pi
            boot_acrophase[b] <- acro_b
            
            if(b %% 50 == 0) incProgress(50 / B)
          }
        })
        
        boot_results <- list(
          mesor_ci = quantile(boot_mesor, c(0.025, 0.975)),
          amplitude_ci = quantile(boot_amplitude, c(0.025, 0.975)),
          acrophase_ci = quantile(boot_acrophase * period / (2 * pi), c(0.025, 0.975)),
          boot_mesor = boot_mesor,
          boot_amplitude = boot_amplitude,
          boot_acrophase = boot_acrophase,
          B = B
        )
      }
      
      # Store results
      values$harmonic_model <- list(
        individual_fits = individual_fits,
        individual_params = individual_params,
        pop_mean_fit = pop_mean_fit,
        group_fits = group_fits,
        boot_results = boot_results,
        time_vec = time_vec,
        original_times = original_times,
        wrap_applied = wrap_applied,
        period = period,
        n_harmonics = n_harmonics,
        trend_type = trend_type,
        include_trend = trend_type != "none",  # For backwards compatibility
        t_offset = t_offset_global,
        t_center = t_center_global,
        using_smoothed = using_smoothed,
        subjects_with_nas = subjects_with_nas,
        Y = Y
      )
      
      showNotification("Harmonic regression complete!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      print(e)
    })
  })
  
  # ==============================================================================
  # HARMONIC REGRESSION OUTPUTS
  # ==============================================================================
  
  # Summary output
  output$harmonic_summary <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    cat("=== Harmonic Regression (Cosinor Analysis) Results ===\n\n")
    cat("Period:", mod$period, "\n")
    cat("Number of harmonics:", mod$n_harmonics, "\n")
    
    # Data source info
    if(isTRUE(mod$using_smoothed)) {
      cat("Data: SMOOTHED (missing values interpolated by FDA)\n")
    } else {
      cat("Data: RAW (no smoothing applied)\n")
      if(!is.null(mod$subjects_with_nas) && mod$subjects_with_nas > 0) {
        cat("  ⚠ ", mod$subjects_with_nas, " subjects have missing values - consider smoothing first!\n")
      }
    }
    
    # Show trend model type and sleep inertia
    trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
    include_inertia <- !is.null(mod$include_inertia) && isTRUE(mod$include_inertia)

    if(trend_type != "none") {
      trend_label <- switch(trend_type,
                            "linear" = "LINEAR (β·t)",
                            "log" = "LOGARITHMIC (β·log(t+1))",
                            "exp_sat" = "SATURATING EXPONENTIAL (A·(1-e^(-t/τ)))",
                            "Unknown")
      cat("Homeostatic trend:", trend_label, "\n")
      cat("Two-process model: Separates Process S (trend) from Process C (circadian)\n")
      cat("  Model equation: Y(t) = M + S(t) + C(t)\n")
    } else {
      cat("Homeostatic trend: None (circadian only)\n")
      cat("  Model equation: Y(t) = M + C(t)\n")
    }
    
    # Count successful and failed fits
    n_total <- length(mod$individual_fits)
    n_success <- sum(sapply(mod$individual_fits, function(f) isTRUE(f$success)))
    n_failed <- n_total - n_success
    
    cat("Number of subjects:", n_total, "\n")
    cat("  - Successfully fitted:", n_success, "\n")
    if(n_failed > 0) {
      cat("  - Failed fits:", n_failed, "(insufficient data points)\n")
      
      # Show which subjects failed
      failed_info <- sapply(seq_along(mod$individual_fits), function(i) {
        f <- mod$individual_fits[[i]]
        if(!isTRUE(f$success)) {
          if(!is.null(f$n_valid)) {
            sprintf("S%d(%d/%d pts)", i, f$n_valid, f$n_required)
          } else {
            sprintf("S%d", i)
          }
        } else {
          NULL
        }
      })
      failed_info <- failed_info[!sapply(failed_info, is.null)]
      if(length(failed_info) <= 10) {
        cat("    Failed subjects:", paste(failed_info, collapse = ", "), "\n")
      } else {
        cat("    Failed subjects:", paste(head(failed_info, 10), collapse = ", "), "...\n")
      }
    }
    
    cat("Number of time points:", length(mod$time_vec), "\n")
    
    # Show time points used
    if(length(mod$time_vec) <= 24) {
      if(isTRUE(mod$wrap_applied) && !is.null(mod$original_times)) {
        cat("Original times (clock):", paste(round(mod$original_times, 1), collapse=", "), "\n")
        cat("Adjusted times (linear):", paste(round(mod$time_vec, 1), collapse=", "), "\n")
        cat("(Times after midnight were adjusted for chronological order)\n")
      } else {
        cat("Time points used:", paste(round(mod$time_vec, 2), collapse=", "), "\n")
      }
    } else {
      cat("Time points: ", round(min(mod$time_vec), 2), "to", round(max(mod$time_vec), 2), 
          "(", length(mod$time_vec), "points)\n")
    }
    
    # Check if spacing is equal
    diffs <- diff(mod$time_vec)
    if(length(unique(round(diffs, 2))) > 1) {
      cat("Spacing: UNEQUAL (", paste(unique(round(diffs, 2)), collapse=", "), ")\n")
    } else {
      cat("Spacing: Equal (", round(diffs[1], 2), ")\n")
    }
    cat("\n")
    
    # Population mean statistics (always calculated)
    if(!is.null(mod$pop_mean_fit)) {
      cat("--- Population Mean Parameters (Vector-Averaged) ---\n")
      pop <- mod$pop_mean_fit
      cat(sprintf("MESOR:     %.3f\n", pop$mean_mesor))
      cat(sprintf("Amplitude (H1): %.3f\n", pop$mean_amplitude))
      cat(sprintf("Acrophase (H1): %.2f (%.2f hours)\n", 
                  pop$mean_acrophase_rad * 180 / pi, pop$mean_acrophase_time))
      
      # Show all harmonics if more than 1
      if(mod$n_harmonics > 1) {
        cat("\n  All Harmonics (vector-averaged):\n")
        for(h in 1:mod$n_harmonics) {
          cat(sprintf("    H%d: Amplitude=%.3f, Acrophase=%.2f hours\n", 
                      h, pop$mean_amplitudes[h], pop$mean_acrophases_time[h]))
        }
      }
      
      cat(sprintf("\nRayleigh test for uniformity (H1):\n"))
      cat(sprintf("  Z = %.3f, p = %.4f\n", pop$rayleigh_z, pop$rayleigh_p))
      cat(sprintf("  Mean resultant length (r̄) = %.3f\n", pop$r_bar))

      # Build and display symbolic equation
      {
        cat("\n--- Model Equation (symbolic) ---\n")
        sym <- "Y(t) = M"
        if(trend_type == "linear")  sym <- paste0(sym, " + \u03b2\u00b7t")
        else if(trend_type == "log")     sym <- paste0(sym, " + \u03b2\u00b7log(t+1)")
        else if(trend_type == "exp_sat") sym <- paste0(sym, " + A_sat\u00b7(1 \u2212 e^(\u2212t/\u03c4))")
        if(mod$n_harmonics == 1) {
          sym <- paste0(sym, sprintf(" + A\u00b7cos(2\u03c0\u00b7t/%.0f \u2212 \u03c6)", mod$period))
        } else {
          for(h_s in 1:mod$n_harmonics) {
            sub_h <- intToUtf8(0x2080 + h_s)
            sym <- paste0(sym, sprintf(" + A%s\u00b7cos(2\u03c0\u00b7%d\u00b7t/%.0f \u2212 \u03c6%s)",
                                       sub_h, h_s, mod$period, sub_h))
          }
        }
        cat(sym, "\n")
        leg <- c("M = MESOR")
        if(mod$n_harmonics == 1) {
          leg <- c(leg, "A = amplitude", "\u03c6 = acrophase (rad)")
        } else {
          leg <- c(leg, "A\u2095 = amplitude of harmonic h", "\u03c6\u2095 = acrophase of harmonic h (rad)")
        }
        leg <- c(leg, sprintf("T = %.0f h (period)", mod$period))
        if(trend_type == "linear")       leg <- c(leg, "\u03b2 = linear trend slope")
        else if(trend_type == "log")     leg <- c(leg, "\u03b2 = log trend slope")
        else if(trend_type == "exp_sat") leg <- c(leg, "A_sat = asymptote", "\u03c4 = time constant (h)")
        cat("  where:", paste(leg, collapse = ", "), "\n")
      }

      # Build and display fitted equation
      cat("\n--- Fitted Model Equation ---\n")
      eq <- sprintf("Y(t) = %.2f", pop$mean_mesor)

      # Add trend component
      trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
      if(trend_type == "linear" && !is.null(pop$indiv_means$trend_linear)) {
        beta <- pop$indiv_means$trend_linear
        if(beta >= 0) {
          eq <- paste0(eq, sprintf(" + %.3f·t", beta))
        } else {
          eq <- paste0(eq, sprintf(" - %.3f·t", abs(beta)))
        }
      } else if(trend_type == "log" && !is.null(pop$indiv_means$trend_log)) {
        beta <- pop$indiv_means$trend_log
        if(beta >= 0) {
          eq <- paste0(eq, sprintf(" + %.3f·log(t+1)", beta))
        } else {
          eq <- paste0(eq, sprintf(" - %.3f·log(t+1)", abs(beta)))
        }
      } else if(trend_type == "exp_sat" && !is.null(pop$indiv_means$A_sat) && !is.null(pop$indiv_means$tau)) {
        A_sat <- pop$indiv_means$A_sat
        tau <- pop$indiv_means$tau
        if(A_sat >= 0) {
          eq <- paste0(eq, sprintf(" + %.2f·(1 - e^(-t/%.1f))", A_sat, tau))
        } else {
          eq <- paste0(eq, sprintf(" - %.2f·(1 - e^(-t/%.1f))", abs(A_sat), tau))
        }
      }

      # Add harmonic components
      for(h in 1:mod$n_harmonics) {
        A <- pop$mean_amplitudes[h]
        phi <- pop$mean_acrophases_rad[h]  # In radians
        omega <- 2 * pi * h / mod$period

        # Convert to cos(ωt - φ) format
        if(A >= 0) {
          eq <- paste0(eq, sprintf(" + %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                   A, h, mod$period, phi))
        } else {
          eq <- paste0(eq, sprintf(" - %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                   abs(A), h, mod$period, phi))
        }
      }

      cat(eq, "\n")
      cat("where t = time in same units as period\n")

      # Show arithmetic means of individual parameters
      if(!is.null(pop$indiv_means)) {
        cat("\n--- Arithmetic Mean of Individual Parameters ---\n")
        cat(sprintf("MESOR:     %.3f (SD=%.3f)\n", pop$indiv_means$mesor, pop$indiv_means$mesor_sd))
        for(h in 1:mod$n_harmonics) {
          amp_key <- paste0("amplitude_", h)
          amp_sd_key <- paste0("amplitude_", h, "_sd")
          acro_key <- paste0("acrophase_time_", h)
          acro_sd_key <- paste0("acrophase_time_", h, "_sd")
          cat(sprintf("H%d Amplitude: %.3f (SD=%.3f)\n", h, pop$indiv_means[[amp_key]], pop$indiv_means[[amp_sd_key]]))
          cat(sprintf("H%d Acrophase: %.2f hours (SD=%.2f)\n", h, pop$indiv_means[[acro_key]], pop$indiv_means[[acro_sd_key]]))
        }
      }
      
      if(!is.null(mod$boot_results)) {
        cat(sprintf("\n95%% Bootstrap CIs (B=%d):\n", mod$boot_results$B))
        cat(sprintf("  MESOR:     [%.3f, %.3f]\n", 
                    mod$boot_results$mesor_ci[1], mod$boot_results$mesor_ci[2]))
        cat(sprintf("  Amplitude: [%.3f, %.3f]\n", 
                    mod$boot_results$amplitude_ci[1], mod$boot_results$amplitude_ci[2]))
        cat(sprintf("  Acrophase: [%.2f, %.2f] hours\n", 
                    mod$boot_results$acrophase_ci[1], mod$boot_results$acrophase_ci[2]))
      }
    }
    
    cat("\n--- Individual Parameter Summary ---\n")
    params <- mod$individual_params
    cat(sprintf("MESOR:     Mean=%.3f, SD=%.3f\n", mean(params$mesor), sd(params$mesor)))
    
    # Show trend parameters based on trend type
    if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
      cat(sprintf("Linear Trend (β): Mean=%.4f, SD=%.4f (units/hour)\n", 
                  mean(params$trend_linear, na.rm=TRUE), sd(params$trend_linear, na.rm=TRUE)))
      mean_trend <- mean(params$trend_linear, na.rm=TRUE)
      if(mean_trend > 0) {
        cat("           (Positive = increasing over time, e.g., increasing sleepiness)\n")
      } else if(mean_trend < 0) {
        cat("           (Negative = decreasing over time)\n")
      }
    } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
      cat(sprintf("Log Trend (β): Mean=%.4f, SD=%.4f (units/log-hour)\n", 
                  mean(params$trend_log, na.rm=TRUE), sd(params$trend_log, na.rm=TRUE)))
      mean_trend <- mean(params$trend_log, na.rm=TRUE)
      if(mean_trend > 0) {
        cat("           (Positive = increasing with diminishing rate)\n")
      } else if(mean_trend < 0) {
        cat("           (Negative = decreasing with diminishing rate)\n")
      }
    } else if(mod$trend_type == "two_process") {
      cat("\n--- Two-Process Homeostatic Parameters ---\n")
      if("beta_S" %in% names(params)) {
        cat(sprintf("β_S (S coefficient): Mean=%.3f, SD=%.3f (units/S-unit)\n",
                    mean(params$beta_S, na.rm=TRUE), sd(params$beta_S, na.rm=TRUE)))
        mean_beta_S <- mean(params$beta_S, na.rm=TRUE)
        if(mean_beta_S > 0) {
          cat("           (Positive = Y increases with sleep pressure)\n")
        } else if(mean_beta_S < 0) {
          cat("           (Negative = Y decreases with sleep pressure)\n")
        }
      }
      if("tau_w" %in% names(params)) {
        cat(sprintf("τ_w (wake time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau_w, na.rm=TRUE), sd(params$tau_w, na.rm=TRUE)))
        cat("           (Time for S to rise to ~63% toward S_max while awake)\n")
      }
      if("tau_s" %in% names(params)) {
        cat(sprintf("τ_s (sleep time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau_s, na.rm=TRUE), sd(params$tau_s, na.rm=TRUE)))
        cat("           (Time for S to decay to ~37% toward S_min while asleep)\n")
      }
    } else if(mod$trend_type == "exp_sat") {
      if("A_sat" %in% names(params)) {
        cat(sprintf("A_sat (asymptote): Mean=%.3f, SD=%.3f (units)\n",
                    mean(params$A_sat, na.rm=TRUE), sd(params$A_sat, na.rm=TRUE)))
      }
      if("tau" %in% names(params)) {
        cat(sprintf("τ (time constant): Mean=%.2f, SD=%.2f (hours)\n",
                    mean(params$tau, na.rm=TRUE), sd(params$tau, na.rm=TRUE)))
        cat("           (Time to reach ~63% of asymptote)\n")
      }
    }

    for(h in 1:mod$n_harmonics) {
      amp_col <- paste0("amplitude_", h)
      acro_col <- paste0("acrophase_time_", h)
      acro_rad_col <- paste0("acrophase_rad_", h)
      effective_period <- mod$period / h
      
      cat(sprintf("H%d Amplitude: Mean=%.3f, SD=%.3f\n", h, mean(params[[amp_col]]), sd(params[[amp_col]])))
      
      # Circular statistics for acrophase
      acro_rad <- params[[acro_rad_col]]
      circ_mean <- circular_mean(acro_rad)
      if(circ_mean < 0) circ_mean <- circ_mean + 2 * pi
      circ_mean_time <- circ_mean * effective_period / (2 * pi)
      circ_sd <- circular_sd(acro_rad)
      circ_sd_time <- if(!is.na(circ_sd)) circ_sd * effective_period / (2 * pi) else NA
      r_bar <- mean_resultant_length(acro_rad)
      
      cat(sprintf("H%d Acrophase: Circular mean=%.2f h, Circ.SD=%.2f h, r̄=%.3f\n", 
                  h, circ_mean_time, ifelse(is.na(circ_sd_time), NA, circ_sd_time), r_bar))
      cat(sprintf("            (Arithmetic mean=%.2f h, Linear SD=%.2f h)\n", 
                  mean(params[[acro_col]]), sd(params[[acro_col]])))
    }
    
    cat(sprintf("R-squared: Mean=%.3f, Range=[%.3f, %.3f]\n",
                mean(params$r_squared), min(params$r_squared), max(params$r_squared)))
    cat(sprintf("Significant rhythms (p<0.05): %d / %d (%.1f%%)\n",
                sum(params$p_value < 0.05), nrow(params),
                100 * sum(params$p_value < 0.05) / nrow(params)))

    # Model selection metrics
    if("aic" %in% names(params)) {
      cat("\n--- Model Selection Metrics ---\n")
      cat(sprintf("AIC (Akaike Information Criterion): Mean=%.2f, SD=%.2f\n",
                  mean(params$aic, na.rm = TRUE), sd(params$aic, na.rm = TRUE)))
      cat(sprintf("AICc (Corrected AIC): Mean=%.2f, SD=%.2f\n",
                  mean(params$aicc, na.rm = TRUE), sd(params$aicc, na.rm = TRUE)))
      cat(sprintf("BIC (Bayesian Information Criterion): Mean=%.2f, SD=%.2f\n",
                  mean(params$bic, na.rm = TRUE), sd(params$bic, na.rm = TRUE)))
      cat(sprintf("LOOCV RMSE (Leave-one-out CV): Mean=%.4f, SD=%.4f\n",
                  mean(params$loocv_rmse, na.rm = TRUE), sd(params$loocv_rmse, na.rm = TRUE)))
      cat("Note: Lower AIC/AICc/BIC/LOOCV values indicate better model fit\n")
    }

    # Variance decomposition: Relative importance of Process S and Process C
    if(!is.null(mod$pop_mean_fit) && !is.null(mod$pop_mean_fit$indiv_means) &&
       !is.null(mod$pop_mean_fit$indiv_means$r_squared_S) &&
       !is.null(mod$pop_mean_fit$indiv_means$r_squared_C)) {
      cat("\n--- Variance Decomposition: Relative Importance ---\n")
      indiv <- mod$pop_mean_fit$indiv_means

      cat(sprintf("R² from Process S (homeostatic): Mean=%.3f, SD=%.3f\n",
                  indiv$r_squared_S, indiv$r_squared_S_sd))
      cat(sprintf("R² from Process C (circadian):   Mean=%.3f, SD=%.3f\n",
                  indiv$r_squared_C, indiv$r_squared_C_sd))
      cat(sprintf("\nProportion of total R² explained by:\n"))
      cat(sprintf("  Process S: %.1f%% (SD=%.1f%%)\n",
                  indiv$percent_S, indiv$percent_S_sd))
      cat(sprintf("  Process C: %.1f%% (SD=%.1f%%)\n",
                  indiv$percent_C, indiv$percent_C_sd))

      # Interpretation helper (with NA check)
      if(!is.na(indiv$percent_S) && !is.na(indiv$percent_C)) {
        if(indiv$percent_S > indiv$percent_C) {
          cat("\nInterpretation: Homeostatic process (S) is the dominant component\n")
        } else if(indiv$percent_C > indiv$percent_S) {
          cat("\nInterpretation: Circadian rhythm (C) is the dominant component\n")
        } else {
          cat("\nInterpretation: Both processes contribute equally\n")
        }
      }
    }

    # Show group-specific statistics if groups exist
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 2) {
      cat("\n--- Group-Specific Parameters ---\n")
      # Print symbolic equation once before listing groups
      {
        cat("Model Equation (symbolic):\n  ")
        sym_g <- "Y(t) = M"
        grp_trend <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
        if(grp_trend == "linear")       sym_g <- paste0(sym_g, " + \u03b2\u00b7t")
        else if(grp_trend == "log")     sym_g <- paste0(sym_g, " + \u03b2\u00b7log(t+1)")
        else if(grp_trend == "exp_sat") sym_g <- paste0(sym_g, " + A_sat\u00b7(1 \u2212 e^(\u2212t/\u03c4))")
        if(mod$n_harmonics == 1) {
          sym_g <- paste0(sym_g, sprintf(" + A\u00b7cos(2\u03c0\u00b7t/%.0f \u2212 \u03c6)", mod$period))
        } else {
          for(h_s in 1:mod$n_harmonics) {
            sub_h <- intToUtf8(0x2080 + h_s)
            sym_g <- paste0(sym_g, sprintf(" + A%s\u00b7cos(2\u03c0\u00b7%d\u00b7t/%.0f \u2212 \u03c6%s)",
                                           sub_h, h_s, mod$period, sub_h))
          }
        }
        cat(sym_g, "\n")
        leg_g <- c("M = MESOR")
        if(mod$n_harmonics == 1) {
          leg_g <- c(leg_g, "A = amplitude", "\u03c6 = acrophase (rad)")
        } else {
          leg_g <- c(leg_g, "A\u2095 = amplitude of harmonic h", "\u03c6\u2095 = acrophase of harmonic h (rad)")
        }
        leg_g <- c(leg_g, sprintf("T = %.0f h (period)", mod$period))
        if(grp_trend == "linear")       leg_g <- c(leg_g, "\u03b2 = linear trend slope")
        else if(grp_trend == "log")     leg_g <- c(leg_g, "\u03b2 = log trend slope")
        else if(grp_trend == "exp_sat") leg_g <- c(leg_g, "A_sat = asymptote", "\u03c4 = time constant (h)")
        cat("  where:", paste(leg_g, collapse = ", "), "\n")
      }
      for(g_name in names(mod$group_fits)) {
        g <- mod$group_fits[[g_name]]
        cat(sprintf("\nGroup '%s' (n=%d):\n", g_name, g$n))
        cat(sprintf("  MESOR:     %.3f (SD=%.3f)\n", g$mean_mesor, g$sd_mesor))

        # Show trend parameters if present
        if(!is.null(g$trend_params) && length(g$trend_params) > 0) {
          for(param_name in names(g$trend_params)) {
            param_label <- switch(param_name,
                                 "trend_linear" = "Linear trend",
                                 "trend_log" = "Log trend",
                                 "A_sat" = "A_sat",
                                 "tau" = "τ",
                                 param_name)
            cat(sprintf("  %s: %.3f (SD=%.3f)\n",
                       param_label, g$trend_params[[param_name]]$mean,
                       g$trend_params[[param_name]]$sd))
          }
        }

        # Show sleep inertia parameters if present
        if(!is.null(g$inertia_params)) {
          cat(sprintf("  W₀:        %.3f (SD=%.3f)\n",
                     g$inertia_params$W0$mean, g$inertia_params$W0$sd))
          cat(sprintf("  τ_W:       %.2f h (SD=%.2f)\n",
                     g$inertia_params$tau_W$mean, g$inertia_params$tau_W$sd))
        }

        for(h in 1:mod$n_harmonics) {
          cat(sprintf("  H%d Amplitude: %.3f (SD=%.3f)\n", h, g$mean_amplitudes[h], g$sd_amplitudes[h]))
          cat(sprintf("  H%d Acrophase: %.2f hours\n", h, g$mean_acrophases_time[h]))
        }

        # Build and display fitted equation for this group
        cat("\n  Fitted equation:\n  ")
        g_eq <- sprintf("Y(t) = %.2f", g$mean_mesor)

        # Add trend component
        if(mod$trend_type == "linear" && !is.null(g$trend_params$trend_linear)) {
          beta <- g$trend_params$trend_linear$mean
          if(beta >= 0) {
            g_eq <- paste0(g_eq, sprintf(" + %.3f·t", beta))
          } else {
            g_eq <- paste0(g_eq, sprintf(" - %.3f·t", abs(beta)))
          }
        } else if(mod$trend_type == "log" && !is.null(g$trend_params$trend_log)) {
          beta <- g$trend_params$trend_log$mean
          if(beta >= 0) {
            g_eq <- paste0(g_eq, sprintf(" + %.3f·log(t+1)", beta))
          } else {
            g_eq <- paste0(g_eq, sprintf(" - %.3f·log(t+1)", abs(beta)))
          }
        } else if(mod$trend_type == "exp_sat" && !is.null(g$trend_params$A_sat) && !is.null(g$trend_params$tau)) {
          A_sat <- g$trend_params$A_sat$mean
          tau <- g$trend_params$tau$mean
          if(!is.na(A_sat) && !is.na(tau)) {
            if(A_sat >= 0) {
              g_eq <- paste0(g_eq, sprintf(" + %.2f·(1 - e^(-t/%.1f))", A_sat, tau))
            } else {
              g_eq <- paste0(g_eq, sprintf(" - %.2f·(1 - e^(-t/%.1f))", abs(A_sat), tau))
            }
          }
        }

        # Add harmonic components
        for(h in 1:mod$n_harmonics) {
          A <- g$mean_amplitudes[h]
          phi <- g$mean_acrophases_rad[h]  # In radians

          if(!is.na(A) && !is.na(phi)) {
            if(A >= 0) {
              g_eq <- paste0(g_eq, sprintf(" + %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                       A, h, mod$period, phi))
            } else {
              g_eq <- paste0(g_eq, sprintf(" - %.2f·cos(2π·%d·t/%.0f - %.2f)",
                                       abs(A), h, mod$period, phi))
            }
          }
        }

        cat(g_eq, "\n")

        # Variance decomposition for this group
        if(!is.null(g$variance_decomp)) {
          vd <- g$variance_decomp
          cat(sprintf("  Variance Decomposition:\n"))
          cat(sprintf("    Process S: %.1f%% (SD=%.1f%%)\n", vd$percent_S, vd$percent_S_sd))
          cat(sprintf("    Process C: %.1f%% (SD=%.1f%%)\n", vd$percent_C, vd$percent_C_sd))
        }
      }
    }
  })
  
  # Parameters table
  output$harmonic_parameters_table <- renderUI({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Build summary table with all harmonics
    param_names <- c("MESOR")
    mean_vals <- c(mean(params$mesor))
    sd_vals <- c(sd(params$mesor))
    min_vals <- c(min(params$mesor))
    max_vals <- c(max(params$mesor))
    
    for(h in 1:mod$n_harmonics) {
      amp_col <- paste0("amplitude_", h)
      acro_col <- paste0("acrophase_time_", h)
      
      param_names <- c(param_names, paste0("Amplitude H", h), paste0("Acrophase H", h, " (hours)"))
      mean_vals <- c(mean_vals, mean(params[[amp_col]]), mean(params[[acro_col]]))
      sd_vals <- c(sd_vals, sd(params[[amp_col]]), sd(params[[acro_col]]))
      min_vals <- c(min_vals, min(params[[amp_col]]), min(params[[acro_col]]))
      max_vals <- c(max_vals, max(params[[amp_col]]), max(params[[acro_col]]))
    }
    
    param_names <- c(param_names, "R²", "% Rhythm")
    mean_vals <- c(mean_vals, mean(params$r_squared), mean(params$percent_rhythm))
    sd_vals <- c(sd_vals, sd(params$r_squared), sd(params$percent_rhythm))
    min_vals <- c(min_vals, min(params$r_squared), min(params$percent_rhythm))
    max_vals <- c(max_vals, max(params$r_squared), max(params$percent_rhythm))
    
    summary_df <- data.frame(
      Parameter = param_names,
      Mean = round(mean_vals, 3),
      SD = round(sd_vals, 3),
      Min = round(min_vals, 3),
      Max = round(max_vals, 3)
    )
    
    tagList(
      h4("Parameter Summary"),
      renderTable(summary_df, striped = TRUE, hover = TRUE, bordered = TRUE)
    )
  })
  
  # Fitted curves plot
  output$harmonic_fit_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    time_fine <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 200)
    subject_select <- input$harmonic_subject_select
    if(is.null(subject_select)) subject_select <- "mean"
    tp_response_type <- if(!is.null(mod$two_process_params$response_type)) mod$two_process_params$response_type else "gaussian"
    
    p <- plot_ly()
    
    if(subject_select == "all") {
      # Overlay all subjects - colored by group if available
      group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')
      group_colors_rgba <- c('rgba(178,34,34,', 'rgba(70,130,180,', 'rgba(34,139,34,', 
                             'rgba(128,0,128,', 'rgba(255,140,0,', 'rgba(139,69,19,')
      
      if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
         !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        
        group_var <- values$covariates[[input$harmonic_group_var]]
        groups <- names(mod$group_fits)
        
        for(i in seq_along(mod$individual_fits)) {
          fit_i <- mod$individual_fits[[i]]
          if(!is.null(fit_i) && fit_i$success) {
            pred_i <- predict_cosinor(fit_i, time_fine)
            grp <- as.character(group_var[i])
            grp_idx <- which(groups == grp)
            line_color <- if(length(grp_idx) > 0) paste0(group_colors_rgba[grp_idx], '0.4)') else 'rgba(100,100,100,0.3)'
            
            p <- p %>% add_lines(x = time_fine, y = pred_i, 
                                 line = list(color = line_color, width = 1),
                                 showlegend = FALSE, hoverinfo = "skip")
          }
        }
        
        # Add group mean curves using all harmonics
        for(g_idx in seq_along(groups)) {
          g_name <- groups[g_idx]
          g_fit <- mod$group_fits[[g_name]]
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            params$group <- group_var[params$subject]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            group_idx <- which(group_var == g_name)
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
            g_pred <- predict_two_process_mean_curve(grp_params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            # Use mean_coefs for full multi-harmonic prediction
            g_pred <- predict_from_coefs(g_fit$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          
          p <- p %>% add_lines(x = time_fine, y = g_pred,
                               line = list(color = group_colors_hex[g_idx], width = 3), 
                               name = paste("Mean:", g_name))
        }
        
        # Show harmonic components for grouped data (same as "mean" view)
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          params$group <- group_var[params$subject]
          
          # Lighter versions of group colors for components
          group_comp_colors <- c('#FF6B6B', '#87CEEB', '#90EE90', '#DDA0DD', '#FFB347', '#D2B48C')
          
          # Show group-specific harmonic components
          for(g_idx in seq_along(groups)) {
            g_name <- groups[g_idx]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            
            if(nrow(grp_params) > 0) {
              grp_mesor <- mean(grp_params$mesor, na.rm = TRUE)
              
              grp_coefs <- c(grp_mesor)
              trend_coefs <- get_mean_trend_coefs(grp_params, mod$trend_type)
              if(length(trend_coefs) > 0) {
                grp_coefs <- c(grp_coefs, trend_coefs)
              }
              for(h in 1:mod$n_harmonics) {
                grp_coefs <- c(grp_coefs, 
                               mean(grp_params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(grp_params[[paste0("beta_sin_", h)]], na.rm = TRUE))
              }
              
              n_trend_coefs <- length(trend_coefs)
              coef_offset <- 1 + n_trend_coefs
              
              for(h in 1:mod$n_harmonics) {
                omega <- 2 * pi * h / mod$period
                beta_cos <- grp_coefs[coef_offset + 2 * h - 1]
                beta_sin <- grp_coefs[coef_offset + 2 * h]
                comp_vals <- grp_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
                
                p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                     line = list(color = group_comp_colors[g_idx], width = 1.5, dash = 'dot'),
                                     name = paste0("H", h, " (", g_name, ")"))
              }
              
              # Plot trend line using helper function
              if(mod$trend_type == "two_process") {
                group_idx <- which(group_var == g_name)
                S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
                trend_line <- compute_two_process_trend_line(grp_params, time_fine, S_t, tp_response_type)
              } else {
                trend_line <- compute_trend_line(grp_params, mod$trend_type, time_fine, mod$t_offset)
              }
              if(!is.null(trend_line)) {
                p <- p %>% add_lines(x = time_fine, y = trend_line,
                                     line = list(color = group_colors_hex[g_idx], width = 1.5, dash = 'dash'),
                                     name = paste0(get_trend_label(mod$trend_type), " (", g_name, ")"))
              }
            }
          }
          
          # Overall harmonics in gray
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = 'gray40', width = 1, dash = 'dot'),
                                 name = paste0("H", h, " (overall)"))
          }
          
          # Plot overall trend line using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            overall_trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            overall_trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(overall_trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = overall_trend_line,
                                 line = list(color = 'black', width = 1, dash = 'dash'),
                                 name = paste0(get_trend_label(mod$trend_type), " (overall)"))
          }
        }
        
      } else {
        # No groups - all same color
        for(i in seq_along(mod$individual_fits)) {
          fit_i <- mod$individual_fits[[i]]
          if(!is.null(fit_i) && fit_i$success) {
            pred_i <- predict_cosinor(fit_i, time_fine)
            p <- p %>% add_lines(x = time_fine, y = pred_i, 
                                 line = list(color = 'rgba(100, 100, 100, 0.3)', width = 1),
                                 showlegend = FALSE, hoverinfo = "skip")
          }
        }
        
        # Add mean curve - either from pop_mean_fit or computed from individual params
        if(!is.null(mod$pop_mean_fit)) {
          pop <- mod$pop_mean_fit
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            mean_pred <- predict_from_coefs(pop$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          p <- p %>% add_lines(x = time_fine, y = mean_pred,
                               line = list(color = 'red', width = 3), name = "Population Mean")
        } else {
          # Compute mean from individual parameters (for individual analysis type)
          params <- mod$individual_params
          mean_mesor <- mean(params$mesor, na.rm = TRUE)
          mean_coefs <- c(mean_mesor)
          
          mean_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(mean_trend_coefs) > 0) {
            mean_coefs <- c(mean_coefs, mean_trend_coefs)
          }
          
          for(h in 1:mod$n_harmonics) {
            mean_coefs <- c(mean_coefs, 
                            mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                            mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            mean_pred <- predict_from_coefs(mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }
          p <- p %>% add_lines(x = time_fine, y = mean_pred,
                               line = list(color = 'red', width = 3), name = "Population Mean")
        }
        
        # Show harmonic components for no-groups case
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste0("H", h))
          }
          
          # Plot trend line using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
      }
      
    } else if(subject_select == "mean") {
      # Show mean curve(s) - either population or by group
      
      # Check if we have group fits
      if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1) {
        # Show each group's mean curve - use hex colors for proper transparency
        group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')  # firebrick, steelblue, forestgreen, purple, orange, brown
        group_colors_rgba <- c('rgba(178,34,34,', 'rgba(70,130,180,', 'rgba(34,139,34,', 
                               'rgba(128,0,128,', 'rgba(255,140,0,', 'rgba(139,69,19,')
        color_idx <- 1
        
        for(g_name in names(mod$group_fits)) {
          g_fit <- mod$group_fits[[g_name]]
          if(mod$trend_type == "two_process") {
            params <- mod$individual_params
            group_var <- values$covariates[[input$harmonic_group_var]]
            params$group <- group_var[params$subject]
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            group_idx <- which(group_var == g_name)
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
            g_pred <- predict_two_process_mean_curve(grp_params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
          } else {
            # Use mean_coefs for full multi-harmonic prediction
            g_pred <- predict_from_coefs(g_fit$mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
          }

          # Only add line if predictions are valid (not all NA/NaN/Inf)
          if(any(is.finite(g_pred))) {
            p <- p %>% add_lines(x = time_fine, y = g_pred,
                                 line = list(color = group_colors_hex[color_idx], width = 3),
                                 name = paste("Group:", g_name))
          }
          
          # Add confidence band if requested (approximate using first harmonic SD)
          if(isTRUE(input$harmonic_show_ci) && !is.null(g_fit$sd_amplitudes)) {
            # Simple approximation: scale curve by amplitude uncertainty
            amp_se <- g_fit$sd_amplitudes[1] / sqrt(g_fit$n)
            scale_upper <- 1 + 1.96 * amp_se / g_fit$mean_amplitudes[1]
            scale_lower <- max(0, 1 - 1.96 * amp_se / g_fit$mean_amplitudes[1])
            
            g_upper <- g_fit$mean_mesor + (g_pred - g_fit$mean_mesor) * scale_upper
            g_lower <- g_fit$mean_mesor + (g_pred - g_fit$mean_mesor) * scale_lower
            
            # Use rgba with 0.25 transparency for confidence band
            ci_color <- paste0(group_colors_rgba[color_idx], '0.25)')
            
            # Add ribbon with legend entry so it can be toggled
            p <- p %>% add_ribbons(x = time_fine, ymin = g_lower, ymax = g_upper,
                                   line = list(color = 'transparent'),
                                   fillcolor = ci_color,
                                   name = paste("95% CI:", g_name),
                                   legendgroup = paste0("group_", g_name),
                                   showlegend = TRUE, hoverinfo = "skip")
          }
          
          color_idx <- color_idx + 1
        }
        
        # Show harmonic components for grouped data
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          params <- mod$individual_params
          group_var <- values$covariates[[input$harmonic_group_var]]
          params$group <- group_var[params$subject]
          
          # Lighter versions of group colors for components
          group_comp_colors <- c('#FF6B6B', '#87CEEB', '#90EE90', '#DDA0DD', '#FFB347', '#D2B48C')
          
          # Show group-specific harmonic components
          g_idx <- 1
          for(g_name in names(mod$group_fits)) {
            grp_params <- params[params$group == g_name & !is.na(params$group), ]
            
            if(nrow(grp_params) > 0) {
              grp_mesor <- mean(grp_params$mesor, na.rm = TRUE)
              
              # Build group-specific coefficients using helper
              grp_coefs <- c(grp_mesor)
              grp_trend_coefs <- get_mean_trend_coefs(grp_params, mod$trend_type)
              if(length(grp_trend_coefs) > 0) {
                grp_coefs <- c(grp_coefs, grp_trend_coefs)
              }
              for(h in 1:mod$n_harmonics) {
                grp_coefs <- c(grp_coefs, 
                               mean(grp_params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(grp_params[[paste0("beta_sin_", h)]], na.rm = TRUE))
              }
              
              n_trend_coefs <- length(grp_trend_coefs)
              coef_offset <- 1 + n_trend_coefs
              
              # Plot each harmonic for this group
              for(h in 1:mod$n_harmonics) {
                omega <- 2 * pi * h / mod$period
                beta_cos <- grp_coefs[coef_offset + 2 * h - 1]
                beta_sin <- grp_coefs[coef_offset + 2 * h]
                comp_vals <- grp_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
                
                p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                     line = list(color = group_comp_colors[g_idx], width = 1.5, dash = 'dot'),
                                     name = paste0("H", h, " (", g_name, ")"))
              }
              
              # Plot trend line using helper function
              if(mod$trend_type == "two_process") {
                group_idx <- which(group_var == g_name)
                S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine, group_idx)
                grp_trend_line <- compute_two_process_trend_line(grp_params, time_fine, S_t, tp_response_type)
              } else {
                grp_trend_line <- compute_trend_line(grp_params, mod$trend_type, time_fine, mod$t_offset)
              }
              if(!is.null(grp_trend_line)) {
                p <- p %>% add_lines(x = time_fine, y = grp_trend_line,
                                     line = list(color = group_colors_hex[g_idx], width = 1.5, dash = 'dash'),
                                     name = paste0(get_trend_label(mod$trend_type), " (", g_name, ")"))
              }
            }
            g_idx <- g_idx + 1
          }
          
          # Also show overall (pooled) harmonics in gray for reference
          overall_mesor <- mean(params$mesor, na.rm = TRUE)
          overall_coefs <- c(overall_mesor)
          overall_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
          if(length(overall_trend_coefs) > 0) {
            overall_coefs <- c(overall_coefs, overall_trend_coefs)
          }
          for(h in 1:mod$n_harmonics) {
            overall_coefs <- c(overall_coefs, 
                               mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                               mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
          }
          
          n_trend_coefs <- length(overall_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- overall_coefs[coef_offset + 2 * h - 1]
            beta_sin <- overall_coefs[coef_offset + 2 * h]
            comp_vals <- overall_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = 'gray40', width = 1, dash = 'dot'),
                                 name = paste0("H", h, " (overall)"))
          }
          
          # Plot overall trend using helper function
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            overall_trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            overall_trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(overall_trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = overall_trend_line,
                                 line = list(color = 'black', width = 1, dash = 'dash'),
                                 name = paste0(get_trend_label(mod$trend_type), " (overall)"))
          }
        }
        
        # Also add individual data points colored by group if requested
        if(isTRUE(input$harmonic_show_data)) {
          group_var <- values$covariates[[input$harmonic_group_var]]
          for(i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if(!is.null(fit_i) && fit_i$success) {
              grp <- as.character(group_var[i])
              grp_idx <- which(names(mod$group_fits) == grp)
              pt_color <- if(length(grp_idx) > 0) paste0(group_colors_rgba[grp_idx], '0.3)') else 'rgba(100,100,100,0.2)'
              
              p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                     marker = list(color = pt_color, size = 3),
                                     showlegend = FALSE, hoverinfo = "skip")
            }
          }
        }
        
      } else {
        # No groups - compute mean from individual parameters
        params <- mod$individual_params
        mean_mesor <- mean(params$mesor, na.rm = TRUE)
        
        # Build mean_coefs from individual parameters using helper
        mean_coefs <- c(mean_mesor)
        
        # Add mean trend coefficients based on type
        mean_trend_coefs <- get_mean_trend_coefs(params, mod$trend_type)
        if(length(mean_trend_coefs) > 0) {
          mean_coefs <- c(mean_coefs, mean_trend_coefs)
        }
        
        # Add mean harmonic coefficients
        for(h in 1:mod$n_harmonics) {
          mean_coefs <- c(mean_coefs, 
                          mean(params[[paste0("beta_cos_", h)]], na.rm = TRUE),
                          mean(params[[paste0("beta_sin_", h)]], na.rm = TRUE))
        }
        
        if(mod$trend_type == "two_process") {
          S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
          mean_pred <- predict_two_process_mean_curve(params, time_fine, mod$period, mod$n_harmonics, S_t, tp_response_type)
        } else {
          mean_pred <- predict_from_coefs(mean_coefs, time_fine, mod$period, mod$n_harmonics, mod$trend_type, mod$t_offset, mod$t_center)
        }
        
        p <- p %>% add_lines(x = time_fine, y = mean_pred,
                             line = list(color = 'firebrick', width = 3), name = "Population Mean")
        
        # Add confidence band if requested (use SD of individual amplitudes)
        if(isTRUE(input$harmonic_show_ci)) {
          amp_sd <- sd(params$amplitude_1, na.rm = TRUE)
          amp_mean <- mean(params$amplitude_1, na.rm = TRUE)
          n_valid <- sum(!is.na(params$amplitude_1))
          amp_se <- amp_sd / sqrt(n_valid)
          
          scale_upper <- 1 + 1.96 * amp_se / amp_mean
          scale_lower <- max(0, 1 - 1.96 * amp_se / amp_mean)
          
          upper_pred <- mean_mesor + (mean_pred - mean_mesor) * scale_upper
          lower_pred <- mean_mesor + (mean_pred - mean_mesor) * scale_lower
          
          p <- p %>% add_ribbons(x = time_fine, ymin = lower_pred, ymax = upper_pred,
                                 line = list(color = 'transparent'),
                                 fillcolor = 'rgba(178, 34, 34, 0.2)',
                                 name = "95% CI")
        }
        
        # Show harmonic components if requested
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          n_trend_coefs <- length(mean_trend_coefs)
          coef_offset <- 1 + n_trend_coefs
          
          for(h in 1:mod$n_harmonics) {
            omega <- 2 * pi * h / mod$period
            beta_cos <- mean_coefs[coef_offset + 2 * h - 1]
            beta_sin <- mean_coefs[coef_offset + 2 * h]
            comp_vals <- mean_mesor + beta_cos * cos(omega * time_fine) + beta_sin * sin(omega * time_fine)
            
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste("H", h, "(τ/", h, ")", sep=""))
          }
          
          # Show trend component if present using helper
          if(mod$trend_type == "two_process") {
            S_t <- compute_mean_S_from_fits(mod$individual_fits, time_fine)
            trend_line <- compute_two_process_trend_line(params, time_fine, S_t, tp_response_type)
          } else {
            trend_line <- compute_trend_line(params, mod$trend_type, time_fine, mod$t_offset)
          }
          if(!is.null(trend_line)) {
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
        
        # Add individual data as faint points if requested
        if(isTRUE(input$harmonic_show_data)) {
          for(i in seq_along(mod$individual_fits)) {
            fit_i <- mod$individual_fits[[i]]
            if(!is.null(fit_i) && fit_i$success) {
              p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                     marker = list(color = 'rgba(100, 100, 100, 0.2)', size = 3),
                                     showlegend = FALSE, hoverinfo = "skip")
            }
          }
        }
      }
      
    } else {
      # Single subject
      i <- as.integer(subject_select)
      fit_i <- mod$individual_fits[[i]]
      
      if(!is.null(fit_i) && fit_i$success) {
        pred_i <- predict_cosinor(fit_i, time_fine)
        
        p <- p %>% add_lines(x = time_fine, y = pred_i,
                             line = list(color = 'steelblue', width = 2), 
                             name = paste("Subject", i))
        
        if(isTRUE(input$harmonic_show_data)) {
          p <- p %>% add_markers(x = fit_i$time, y = fit_i$y,
                                 marker = list(color = 'steelblue', size = 6),
                                 name = "Observed Data")
        }
        
        # Show harmonic components if requested
        if(isTRUE(input$harmonic_show_components) && mod$n_harmonics >= 1) {
          components <- get_harmonic_components(fit_i, time_fine)
          colors <- c("green", "orange", "purple", "brown", "pink", "cyan", "magenta", "olive")
          for(h in 1:mod$n_harmonics) {
            comp_name <- paste0("harmonic_", h)
            comp_vals <- components$mesor[1] + components[[comp_name]]
            p <- p %>% add_lines(x = time_fine, y = comp_vals,
                                 line = list(color = colors[h], width = 1.5, dash = 'dot'),
                                 name = paste("H", h, "(τ/", h, ")", sep=""))
          }
          
          # Show linear trend component if present
          if(!is.null(components$trend) && mod$trend_type != "none") {
            trend_line <- components$mesor[1] + components$trend
            p <- p %>% add_lines(x = time_fine, y = trend_line,
                                 line = list(color = 'black', width = 1.5, dash = 'dash'),
                                 name = get_trend_label(mod$trend_type))
          }
        }
      } else {
        # Fit failed - show detailed message
        fail_msg <- paste("Subject", i, "fit failed")
        if(!is.null(fit_i) && !is.null(fit_i$n_valid)) {
          fail_msg <- paste0(fail_msg, "\n(", fit_i$n_valid, " valid points, need ", fit_i$n_required, "+)")
        }
        if(!is.null(fit_i) && !is.null(fit_i$message)) {
          fail_msg <- paste0(fail_msg, "\nReason: ", fit_i$message)
        }
        
        p <- p %>% add_annotations(
          x = 0.5, y = 0.5,
          text = fail_msg,
          showarrow = FALSE,
          xref = "paper", yref = "paper",
          font = list(size = 14, color = "red")
        )
      }
    }
    
    # Create better x-axis with clock time labels
    time_range <- range(mod$time_vec)
    
    # Generate tick values and labels
    if(max(mod$time_vec) > mod$period) {
      # Times wrap around - create clock-style labels (modulo period)
      tick_interval <- if(mod$period <= 12) 1 else if(mod$period <= 24) 2 else if(mod$period <= 48) 4 else mod$period / 12
      tick_vals <- seq(floor(time_range[1]), ceiling(time_range[2]), by = tick_interval)
      tick_text <- sapply(tick_vals, function(t) {
        t_mod <- t %% mod$period
        if(mod$period == 24) {
          paste0(t_mod, ":00")
        } else {
          round(t_mod, 1)
        }
      })
      x_title <- paste0("Time (period = ", mod$period, ", spanning wrap-around)")
    } else {
      tick_vals <- NULL
      tick_text <- NULL
      x_title <- paste("Time (period =", mod$period, ")")
    }
    
    x_axis <- list(title = x_title)
    if(!is.null(tick_vals)) {
      x_axis$tickmode <- "array"
      x_axis$tickvals <- tick_vals
      x_axis$ticktext <- tick_text
    }
    
    p %>% layout(
      title = "Harmonic Regression Fit",
      xaxis = x_axis,
      yaxis = list(title = "Response"),
      showlegend = TRUE
    )
  })
  
  # Polar plot for acrophase
  output$harmonic_polar_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_polar)) as.integer(input$selected_harmonic_polar) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    # Get amplitude and acrophase for selected harmonic
    amp_col <- paste0("amplitude_", h)
    acro_rad_col <- paste0("acrophase_rad_", h)
    
    # Convert to degrees for polar plot
    theta_deg <- params[[acro_rad_col]] * 180 / pi
    r <- params[[amp_col]]
    
    p <- plot_ly(type = 'scatterpolar', mode = 'markers')
    
    # Check if we have group fits - color points by group
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 &&
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {

      group_var <- values$covariates[[input$harmonic_group_var]]
      group_colors <- c('firebrick', 'steelblue', 'forestgreen', 'purple', 'orange', 'brown')

      groups <- names(mod$group_fits)
      for(g_idx in seq_along(groups)) {
        g_name <- groups[g_idx]
        g_mask <- !is.na(group_var[params$subject]) & group_var[params$subject] == g_name

        if(sum(g_mask) > 0) {
          p <- p %>% add_trace(
            r = r[g_mask], theta = theta_deg[g_mask],
            type = 'scatterpolar', mode = 'markers',
            marker = list(size = 8, color = group_colors[g_idx], opacity = 0.7),
            name = paste("Group:", g_name)
          )
        }

        # Add group mean vector for selected harmonic
        g_fit <- mod$group_fits[[g_name]]
        acro_deg <- g_fit$mean_acrophases_rad[h] * 180 / pi
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, g_fit$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          type = 'scatterpolar', mode = 'lines+markers',
          line = list(color = group_colors[g_idx], width = 3),
          marker = list(size = 12, color = group_colors[g_idx], symbol = 'diamond'),
          name = paste("Mean:", g_name)
        )
      }

      # Add population mean vector if requested (even when groups present)
      if(isTRUE(input$polar_show_mean) && !is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        acro_deg <- pop$mean_acrophases_rad[h] * 180 / pi
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, pop$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          mode = 'lines+markers',
          line = list(color = 'black', width = 4, dash = 'dash'),
          marker = list(size = 14, color = 'black', symbol = 'star'),
          name = "Overall Population Mean"
        )
      }

    } else {
      # No groups - show all points same color
      p <- p %>% add_trace(
        r = r, theta = theta_deg,
        marker = list(size = 8, color = 'steelblue', opacity = 0.7),
        name = "Individual"
      )

      # Add mean vector if requested
      if(isTRUE(input$polar_show_mean) && !is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        acro_deg <- pop$mean_acrophases_rad[h] * 180 / pi
        if(acro_deg < 0) acro_deg <- acro_deg + 360

        p <- p %>% add_trace(
          r = c(0, pop$mean_amplitudes[h]),
          theta = c(0, acro_deg),
          mode = 'lines+markers',
          line = list(color = 'red', width = 3),
          marker = list(size = 12, color = 'red', symbol = 'diamond'),
          name = "Population Mean"
        )
      }
    }

    # Add confidence ellipse if requested (for all data, regardless of groups)
    if(isTRUE(input$polar_show_ellipse) && length(r) >= 3) {
      # Convert polar to Cartesian for ellipse calculation
      theta_rad <- theta_deg * pi / 180
      x <- r * cos(theta_rad)
      y <- r * sin(theta_rad)

      # Remove NAs
      valid <- !is.na(x) & !is.na(y)
      x <- x[valid]
      y <- y[valid]

      if(length(x) >= 3) {
        # Calculate 95% confidence ellipse
        mx <- mean(x)
        my <- mean(y)

        # Covariance matrix
        cov_mat <- cov(cbind(x, y))

        # Eigenvalues and eigenvectors
        eig <- eigen(cov_mat)

        # Chi-square value for 95% confidence (2 degrees of freedom)
        chi_sq <- qchisq(0.95, df = 2)

        # Ellipse parameters
        a <- sqrt(chi_sq * eig$values[1])  # Semi-major axis
        b <- sqrt(chi_sq * eig$values[2])  # Semi-minor axis
        angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])  # Rotation angle

        # Generate ellipse points
        t <- seq(0, 2*pi, length.out = 100)
        ellipse_x <- mx + a * cos(t) * cos(angle) - b * sin(t) * sin(angle)
        ellipse_y <- my + a * cos(t) * sin(angle) + b * sin(t) * cos(angle)

        # Convert back to polar coordinates
        ellipse_r <- sqrt(ellipse_x^2 + ellipse_y^2)
        ellipse_theta_rad <- atan2(ellipse_y, ellipse_x)
        ellipse_theta_deg <- ellipse_theta_rad * 180 / pi
        ellipse_theta_deg[ellipse_theta_deg < 0] <- ellipse_theta_deg[ellipse_theta_deg < 0] + 360

        # Add ellipse as a trace
        p <- p %>% add_trace(
          r = ellipse_r,
          theta = ellipse_theta_deg,
          type = 'scatterpolar',
          mode = 'lines',
          line = list(color = 'rgba(255, 0, 0, 0.5)', width = 2, dash = 'dot'),
          name = "95% Confidence Ellipse",
          showlegend = TRUE
        )
      }
    }
    
    # Adjust angular axis based on harmonic
    # For H2 (12h period), show 12 hours; for H3, show 8 hours, etc.
    effective_period <- mod$period / h
    n_ticks <- min(12, effective_period)
    tick_step <- effective_period / n_ticks
    tick_vals <- seq(0, 360 - 360/n_ticks, by = 360/n_ticks)
    tick_labels <- sprintf("%.1fh", seq(0, effective_period - tick_step, by = tick_step))
    
    p %>% layout(
      title = paste("Acrophase Polar Plot - Harmonic", h, "(period =", round(effective_period, 1), "h)"),
      polar = list(
        radialaxis = list(title = "Amplitude"),
        angularaxis = list(
          direction = "clockwise",
          rotation = 90,
          tickmode = "array",
          tickvals = tick_vals,
          ticktext = tick_labels
        )
      ),
      showlegend = TRUE
    )
  })
  
  # Amplitude histogram
  output$harmonic_amplitude_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_dist)) as.integer(input$selected_harmonic_dist) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    effective_period <- mod$period / h
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = .data[[amp_col]], fill = as.factor(group))) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 15) +
        scale_fill_brewer(palette = "Set1", name = "Group") +
        theme_minimal() +
        labs(title = paste0("Distribution of Amplitudes (H", h, ", period=", round(effective_period, 1), "h) by Group"), 
             x = "Amplitude", y = "Count")
      
      ggplotly(g)
    } else {
      plot_ly(x = params[[amp_col]], type = "histogram", 
              marker = list(color = 'steelblue', line = list(color = 'white', width = 1))) %>%
        layout(title = paste0("Distribution of Amplitudes (H", h, ", period=", round(effective_period, 1), "h)"),
               xaxis = list(title = "Amplitude"),
               yaxis = list(title = "Count"))
    }
  })
  
  # Acrophase histogram (circular)
  output$harmonic_acrophase_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    period <- mod$period
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_dist)) as.integer(input$selected_harmonic_dist) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    acro_col <- paste0("acrophase_time_", h)
    effective_period <- period / h
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = .data[[acro_col]], fill = as.factor(group))) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 12) +
        scale_fill_brewer(palette = "Set1", name = "Group") +
        scale_x_continuous(limits = c(0, effective_period)) +
        theme_minimal() +
        labs(title = paste0("Distribution of Acrophases (H", h, ") by Group"), 
             x = paste("Acrophase (hours, period =", round(effective_period, 1), ")"), y = "Count")
      
      ggplotly(g)
    } else {
      plot_ly(x = params[[acro_col]], type = "histogram",
              marker = list(color = 'firebrick', line = list(color = 'white', width = 1))) %>%
        layout(title = paste0("Distribution of Acrophases (H", h, ")"),
               xaxis = list(title = paste("Acrophase (hours, period =", round(effective_period, 1), ")"),
                            range = c(0, effective_period)),
               yaxis = list(title = "Count"))
    }
  })
  
  # MESOR by group plot
  output$harmonic_mesor_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Check if we have groups
    if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
       !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      
      group_var <- values$covariates[[input$harmonic_group_var]]
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]  # Remove NAs
      
      g <- ggplot(params, aes(x = group, y = mesor, fill = group)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        labs(title = "MESOR by Group", x = "Group", y = "MESOR") +
        theme(legend.position = "none")
      
      ggplotly(g)
    } else {
      plot_ly(y = params$mesor, type = "box", 
              marker = list(color = 'forestgreen'),
              boxpoints = "all", jitter = 0.3) %>%
        layout(title = "Distribution of MESOR",
               yaxis = list(title = "MESOR"))
    }
  })
  
  # Trend parameter histogram
  output$harmonic_trend_hist <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    # Check if trend is included
    if(mod$trend_type == "none") {
      plot_ly() %>% 
        layout(title = "No Homeostatic Trend Model Selected",
               annotations = list(
                 list(x = 0.5, y = 0.5, text = "Select a trend model\nto view distribution",
                      showarrow = FALSE, xref = "paper", yref = "paper",
                      font = list(size = 14, color = "gray"))))
    } else {
      params <- mod$individual_params
      
      # Get trend column based on type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        trend_vals <- params$trend_linear
        trend_label <- "Linear Trend (β, units/h)"
        trend_title <- "Distribution of Linear Trend Coefficient"
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        trend_vals <- params$trend_log
        trend_label <- "Log Trend (β)"
        trend_title <- "Distribution of Logarithmic Trend Coefficient"
      } else if(mod$trend_type == "exp_sat" && "A_sat" %in% names(params)) {
        trend_vals <- params$A_sat
        trend_label <- "A_sat (asymptote, units)"
        trend_title <- "Distribution of Saturating Exponential Asymptote"
      } else {
        trend_vals <- NULL
      }
      
      if(is.null(trend_vals)) {
        plot_ly() %>% layout(title = "Trend parameters not available")
      } else {
        # Check if groups exist
        if(!is.null(mod$group_fits) && length(mod$group_fits) >= 1 && 
           !is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
          
          group_var <- values$covariates[[input$harmonic_group_var]]
          params$group <- as.factor(group_var[params$subject])
          params$trend_val <- trend_vals
          params <- params[!is.na(params$group) & !is.na(params$trend_val), ]
          
          g <- ggplot(params, aes(x = group, y = trend_val, fill = group)) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_brewer(palette = "Set1") +
            theme_minimal() +
            labs(title = paste(trend_title, "by Group"), 
                 x = "Group", y = trend_label) +
            theme(legend.position = "none")
          
          ggplotly(g)
        } else {
          trend_vals <- trend_vals[!is.na(trend_vals)]
          plot_ly(y = trend_vals, type = "box", 
                  marker = list(color = 'steelblue'),
                  boxpoints = "all", jitter = 0.3) %>%
            layout(title = trend_title,
                   yaxis = list(title = trend_label))
        }
      }
    }
  })
  
  # Individual results table
  output$harmonic_individual_table <- DT::renderDataTable({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    params <- mod$individual_params
    
    # Create a cleaner display table
    display_df <- data.frame(
      Subject = params$subject,
      MESOR = round(params$mesor, 3),
      R_squared = round(params$r_squared, 3),
      Pct_Rhythm = round(params$percent_rhythm, 1),
      p_value = format(params$p_value, digits = 3, scientific = TRUE)
    )
    
    # Add trend parameters based on trend type
    if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
      display_df$Trend_Linear <- round(params$trend_linear, 4)
    } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
      display_df$Trend_Log <- round(params$trend_log, 4)
    } else if(mod$trend_type == "exp_sat") {
      if("A_sat" %in% names(params)) {
        display_df$A_sat <- round(params$A_sat, 3)
      }
      if("tau" %in% names(params)) {
        display_df$Tau_hrs <- round(params$tau, 2)
      }
    }
    
    # Add columns for each harmonic
    for(h in 1:mod$n_harmonics) {
      display_df[[paste0("Amp_H", h)]] <- round(params[[paste0("amplitude_", h)]], 3)
      display_df[[paste0("Acro_H", h, "_hrs")]] <- round(params[[paste0("acrophase_time_", h)]], 2)
    }
    
    # Add group if available
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      display_df$Group <- group_var[params$subject]
      # Move Group to second column
      display_df <- display_df[, c("Subject", "Group", setdiff(names(display_df), c("Subject", "Group")))]
    }
    
    DT::datatable(display_df, 
                  options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatStyle('R_squared', 
                      backgroundColor = DT::styleInterval(c(0.5, 0.8), c('#ffcccc', '#ffffcc', '#ccffcc')))
  })
  
  # Export individual parameters
  output$export_harmonic_individual <- downloadHandler(
    filename = function() paste0("harmonic_individual_params_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$harmonic_model)
      mod <- values$harmonic_model
      params <- mod$individual_params
      
      # Add group variable if available
      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]
        params$group <- group_var[params$subject]
      }
      
      write.csv(params, file, row.names = FALSE)
    }
  )
  
  # Residual plot
  output$harmonic_residual_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    # Collect all residuals
    all_fitted <- c()
    all_resid <- c()
    
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_fitted <- c(all_fitted, fit_i$fitted)
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    plot_ly(x = all_fitted, y = all_resid, type = 'scatter', mode = 'markers',
            marker = list(color = 'steelblue', opacity = 0.5, size = 4)) %>%
      add_segments(x = min(all_fitted), xend = max(all_fitted), y = 0, yend = 0,
                   line = list(color = 'red', dash = 'dash')) %>%
      layout(title = "Residuals vs Fitted",
             xaxis = list(title = "Fitted Values"),
             yaxis = list(title = "Residuals"))
  })
  
  # QQ plot
  output$harmonic_qq_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    all_resid <- c()
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    qq <- qqnorm(all_resid, plot.it = FALSE)
    
    plot_ly(x = qq$x, y = qq$y, type = 'scatter', mode = 'markers',
            marker = list(color = 'steelblue', size = 4)) %>%
      add_lines(x = range(qq$x), y = range(qq$x) * sd(all_resid) + mean(all_resid),
                line = list(color = 'red', dash = 'dash'), name = "Reference") %>%
      layout(title = "Q-Q Plot of Residuals",
             xaxis = list(title = "Theoretical Quantiles"),
             yaxis = list(title = "Sample Quantiles"))
  })
  
  # Goodness of fit statistics
  output$harmonic_gof_stats <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    all_resid <- c()
    for(fit_i in mod$individual_fits) {
      if(!is.null(fit_i) && fit_i$success) {
        all_resid <- c(all_resid, fit_i$residuals)
      }
    }
    
    cat("=== Residual Diagnostics ===\n\n")
    cat(sprintf("Total residuals: %d\n", length(all_resid)))
    cat(sprintf("Mean residual: %.4f\n", mean(all_resid)))
    cat(sprintf("SD of residuals: %.4f\n", sd(all_resid)))
    cat(sprintf("Skewness: %.3f\n", mean((all_resid - mean(all_resid))^3) / sd(all_resid)^3))
    cat(sprintf("Kurtosis: %.3f\n", mean((all_resid - mean(all_resid))^4) / sd(all_resid)^4 - 3))
    
    # Shapiro-Wilk test (on sample if too many observations)
    if(length(all_resid) > 5000) {
      samp_resid <- sample(all_resid, 5000)
    } else {
      samp_resid <- all_resid
    }
    sw <- shapiro.test(samp_resid)
    cat(sprintf("\nShapiro-Wilk test: W = %.4f, p = %.4f\n", sw$statistic, sw$p.value))
    if(sw$p.value < 0.05) {
      cat("  (Significant departure from normality)\n")
    } else {
      cat("  (No significant departure from normality)\n")
    }
  })
  
  # Group comparison plot
  output$harmonic_group_comparison_plot <- renderPlotly({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      return(plotly_empty() %>% layout(title = "Select a group variable for comparison"))
    }
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_group)) as.integer(input$selected_harmonic_group) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    acro_col <- paste0("acrophase_time_", h)
    acro_rad_col <- paste0("acrophase_rad_", h)
    effective_period <- mod$period / h
    
    # Get individual params with group info for R-squared
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- group_var[params$subject]
      params <- params[!is.na(params$group), ]  # Remove missings
    }
    
    # Create comparison data
    group_df <- data.frame(
      group = character(),
      parameter = character(),
      value = numeric(),
      se = numeric()
    )
    
    for(g_name in names(mod$group_fits)) {
      g <- mod$group_fits[[g_name]]
      
      # Get R-squared stats for this group
      grp_rsq <- params$r_squared[params$group == g_name]
      rsq_mean <- mean(grp_rsq, na.rm = TRUE)
      rsq_se <- sd(grp_rsq, na.rm = TRUE) / sqrt(length(grp_rsq))
      
      # Get amplitude and acrophase for selected harmonic from individual params
      grp_amp <- params[[amp_col]][params$group == g_name]
      grp_acro_rad <- params[[acro_rad_col]][params$group == g_name]
      
      # Compute circular mean and SE for acrophase
      circ_mean_rad <- circular_mean(grp_acro_rad)
      if(circ_mean_rad < 0) circ_mean_rad <- circ_mean_rad + 2 * pi
      circ_mean_time <- circ_mean_rad * effective_period / (2 * pi)
      circ_se_rad <- circular_se(grp_acro_rad)
      circ_se_time <- if(!is.na(circ_se_rad)) circ_se_rad * effective_period / (2 * pi) else NA
      
      group_df <- rbind(group_df, 
                        data.frame(group = g_name, parameter = "MESOR", 
                                   value = g$mean_mesor, se = g$sd_mesor / sqrt(g$n)),
                        data.frame(group = g_name, parameter = paste0("Amplitude (H", h, ")"), 
                                   value = mean(grp_amp, na.rm = TRUE), 
                                   se = sd(grp_amp, na.rm = TRUE) / sqrt(length(grp_amp))),
                        data.frame(group = g_name, parameter = paste0("Acrophase (H", h, ")"), 
                                   value = circ_mean_time, se = circ_se_time),
                        data.frame(group = g_name, parameter = "R-squared", 
                                   value = rsq_mean, se = rsq_se))
      
      # Add trend parameters based on trend type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        grp_trend <- params$trend_linear[params$group == g_name]
        grp_trend <- grp_trend[!is.na(grp_trend)]
        if(length(grp_trend) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "Linear Trend (β)", 
                                       value = mean(grp_trend, na.rm = TRUE), 
                                       se = sd(grp_trend, na.rm = TRUE) / sqrt(length(grp_trend))))
        }
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        grp_trend <- params$trend_log[params$group == g_name]
        grp_trend <- grp_trend[!is.na(grp_trend)]
        if(length(grp_trend) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "Log Trend (β)", 
                                       value = mean(grp_trend, na.rm = TRUE), 
                                       se = sd(grp_trend, na.rm = TRUE) / sqrt(length(grp_trend))))
        }
      } else if(mod$trend_type == "exp_sat") {
        if("A_sat" %in% names(params)) {
          grp_asat <- params$A_sat[params$group == g_name]
          grp_asat <- grp_asat[!is.na(grp_asat)]
          if(length(grp_asat) > 0) {
            group_df <- rbind(group_df,
                              data.frame(group = g_name, parameter = "A_sat (asymptote)", 
                                         value = mean(grp_asat, na.rm = TRUE), 
                                         se = sd(grp_asat, na.rm = TRUE) / sqrt(length(grp_asat))))
          }
        }
        if("tau" %in% names(params)) {
          grp_tau <- params$tau[params$group == g_name]
          grp_tau <- grp_tau[!is.na(grp_tau)]
          if(length(grp_tau) > 0) {
            group_df <- rbind(group_df,
                              data.frame(group = g_name, parameter = "τ (time constant, h)", 
                                         value = mean(grp_tau, na.rm = TRUE), 
                                         se = sd(grp_tau, na.rm = TRUE) / sqrt(length(grp_tau))))
          }
        }
      }

      # Add sleep inertia parameters if present
      include_inertia <- !is.null(mod$include_inertia) && isTRUE(mod$include_inertia)
      if(include_inertia && "W0" %in% names(params)) {
        grp_W0 <- params$W0[params$group == g_name]
        grp_W0 <- grp_W0[!is.na(grp_W0)]
        if(length(grp_W0) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "W₀ (inertia)",
                                       value = mean(grp_W0, na.rm = TRUE),
                                       se = sd(grp_W0, na.rm = TRUE) / sqrt(length(grp_W0))))
        }
      }
      if(include_inertia && "tau_W" %in% names(params)) {
        grp_tau_W <- params$tau_W[params$group == g_name]
        grp_tau_W <- grp_tau_W[!is.na(grp_tau_W)]
        if(length(grp_tau_W) > 0) {
          group_df <- rbind(group_df,
                            data.frame(group = g_name, parameter = "τ_W (decay, h)",
                                       value = mean(grp_tau_W, na.rm = TRUE),
                                       se = sd(grp_tau_W, na.rm = TRUE) / sqrt(length(grp_tau_W))))
        }
      }
    }

    # Build named color vector matching the line graph palette
    group_colors_hex <- c('#B22222', '#4682B4', '#228B22', '#800080', '#FF8C00', '#8B4513')
    group_names <- names(mod$group_fits)
    named_colors <- setNames(group_colors_hex[seq_along(group_names)], group_names)

    # Faceted bar plot
    g <- ggplot(group_df, aes(x = group, y = value, fill = group)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = value - 1.96*se, ymax = value + 1.96*se),
                    width = 0.2, na.rm = TRUE) +
      facet_wrap(~parameter, scales = "free_y") +
      scale_fill_manual(values = named_colors) +
      theme_minimal() +
      labs(title = paste0("Group Comparison - Harmonic ", h, " (period = ", round(effective_period, 1), "h)"),
           x = "", y = "Value") +
      theme(legend.position = "none")
    
    ggplotly(g)
  })
  
  # Group comparison test results
  output$harmonic_group_test_results <- renderPrint({
    req(values$harmonic_model)
    mod <- values$harmonic_model
    
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      cat("Select a group variable for statistical comparison.\n")
      return()
    }
    
    # Get selected harmonic (default to 1)
    h <- if(!is.null(input$selected_harmonic_group)) as.integer(input$selected_harmonic_group) else 1
    h <- min(h, mod$n_harmonics)  # Safety check
    
    amp_col <- paste0("amplitude_", h)
    acro_col <- paste0("acrophase_time_", h)
    effective_period <- mod$period / h
    
    cat("=== Group Comparison Statistics ===\n")
    cat(sprintf("Harmonic: H%d (period = %.1f h)\n\n", h, effective_period))
    
    # Get group variable and params
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- group_var[params$subject]
      
      # Remove rows with missing group or key parameters
      params <- params[!is.na(params$group) & !is.na(params$mesor) & 
                         !is.na(params[[amp_col]]) & !is.na(params[[acro_col]]), ]
      
      n_groups <- length(unique(params$group))
      cat(sprintf("Groups: %d, Total N: %d (after removing missings)\n\n", n_groups, nrow(params)))
      
      cat("--- MESOR Comparison ---\n")
      if(n_groups == 2) {
        tt <- t.test(mesor ~ group, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_mesor <- aov(mesor ~ group, data = params)
        cat("ANOVA:\n")
        print(summary(aov_mesor))
      }
      
      # Trend comparison (if trend is present)
      trend_type <- if(!is.null(mod$trend_type)) mod$trend_type else "none"
      if(trend_type != "none") {
        trend_col <- switch(trend_type,
                            "linear" = "trend_linear",
                            "log" = "trend_log",
                            "exp_sat" = "A_sat",
                            NULL)
        
        if(!is.null(trend_col) && trend_col %in% names(params)) {
          cat(sprintf("\n--- %s Comparison ---\n", get_trend_label(trend_type)))
          trend_formula <- as.formula(paste(trend_col, "~ group"))
          if(n_groups == 2) {
            tt <- t.test(trend_formula, data = params)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
          } else {
            aov_trend <- aov(trend_formula, data = params)
            cat("ANOVA:\n")
            print(summary(aov_trend))
          }
          
          # Group means for trend
          for(g in unique(params$group)) {
            g_trend <- params[[trend_col]][params$group == g]
            cat(sprintf("  %s: mean = %.4f (SD = %.4f, n=%d)\n", 
                        g, mean(g_trend, na.rm = TRUE), sd(g_trend, na.rm = TRUE), length(g_trend)))
          }
          
          # For exp_sat, also compare tau
          if(trend_type == "exp_sat" && "tau" %in% names(params)) {
            cat("\n--- Time Constant (τ) Comparison ---\n")
            tau_formula <- as.formula("tau ~ group")
            if(n_groups == 2) {
              tt <- t.test(tau_formula, data = params)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_tau <- aov(tau_formula, data = params)
              cat("ANOVA:\n")
              print(summary(aov_tau))
            }
            
            for(g in unique(params$group)) {
              g_tau <- params$tau[params$group == g]
              cat(sprintf("  %s: mean τ = %.2f h (SD = %.2f, n=%d)\n", 
                          g, mean(g_tau, na.rm = TRUE), sd(g_tau, na.rm = TRUE), length(g_tau)))
            }
          }
        }
      }

      cat(sprintf("\n--- Amplitude (H%d) Comparison ---\n", h))
      amp_formula <- as.formula(paste(amp_col, "~ group"))
      if(n_groups == 2) {
        tt <- t.test(amp_formula, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_amp <- aov(amp_formula, data = params)
        cat("ANOVA:\n")
        print(summary(aov_amp))
      }
      
      cat(sprintf("\n--- Acrophase (H%d) Comparison ---\n", h))
      cat(strrep("-", 60), "\n")
      cat("TWO COMPLEMENTARY TESTS ARE PROVIDED:\n")
      cat("  (1) Watson-Williams: tests timing (phase) only.\n")
      cat("      Every subject counts equally regardless of rhythm strength.\n")
      cat("      Use this to answer: 'Do groups peak at different times?'\n")
      cat("  (2) Hotelling's T² on (beta_cos, beta_sin): tests the full\n")
      cat("      rhythmic vector (phase + amplitude combined). Subjects with\n")
      cat("      stronger rhythms carry more weight (amplitude-weighted).\n")
      cat("      The group means it reports match the 'Group-Specific\n")
      cat("      Parameters' panel exactly.\n")
      cat("      Use this to answer: 'Do groups differ in their overall\n")
      cat("      rhythmic profile?' — but note a significant result could\n")
      cat("      reflect phase, amplitude, or both.\n")
      cat(strrep("-", 60), "\n\n")

      # Get acrophase in radians for each group
      acro_rad_col <- paste0("acrophase_rad_", h)
      beta_cos_col <- paste0("beta_cos_", h)
      beta_sin_col <- paste0("beta_sin_", h)
      groups <- unique(params$group)
      angles_list <- lapply(groups, function(g) {
        params[[acro_rad_col]][params$group == g]
      })
      names(angles_list) <- groups

      # (1) Watson-Williams test
      cat("(1) Watson-Williams test (unweighted circular mean):\n")
      ww <- watson_williams_test(angles_list)

      if(!is.null(ww$message)) {
        cat(sprintf("  %s\n", ww$message))
      }

      if(!is.na(ww$F)) {
        cat(sprintf("    F(%d, %d) = %.3f, p = %.4f\n", ww$df1, ww$df2, ww$F, ww$p))
        cat(sprintf("    Mean resultant length (r̄) = %.3f", ww$r_bar))
        if(ww$r_bar >= 0.7) {
          cat(" (high concentration)\n")
        } else if(ww$r_bar >= 0.45) {
          cat(" (moderate concentration)\n")
        } else {
          cat(" (low concentration - interpret with caution)\n")
        }
      }

      cat("\n  Group circular statistics (unweighted):\n")
      for(g in groups) {
        g_angles <- params[[acro_rad_col]][params$group == g]
        g_mean <- circular_mean(g_angles)
        if(g_mean < 0) g_mean <- g_mean + 2 * pi
        g_mean_time <- g_mean * effective_period / (2 * pi)
        g_sd <- circular_sd(g_angles)
        g_sd_time <- if(!is.na(g_sd)) g_sd * effective_period / (2 * pi) else NA
        g_r <- mean_resultant_length(g_angles)
        cat(sprintf("    %s: mean = %.2f h, circ.SD = %.2f h, r̄ = %.3f (n=%d)\n",
                    g, g_mean_time, ifelse(is.na(g_sd_time), NA, g_sd_time), g_r, length(g_angles)))
      }

      # (2) Hotelling's T² test
      cat(sprintf("\n(2) Hotelling's T² on (beta_cos_%d, beta_sin_%d) (amplitude-weighted):\n", h, h))
      bc_list <- lapply(groups, function(g) params[[beta_cos_col]][params$group == g])
      bs_list <- lapply(groups, function(g) params[[beta_sin_col]][params$group == g])
      ht <- hotelling_t2(bc_list, bs_list)

      if(!is.null(ht$message)) {
        cat(sprintf("  %s\n", ht$message))
      } else {
        cat(sprintf("    F(%d, %d) = %.3f, p = %.4f\n", ht$df1, ht$df2, ht$F, ht$p))
      }

      cat("\n  Group means (amplitude-weighted) — match 'Group-Specific Parameters' panel:\n")
      for(g in groups) {
        bc <- params[[beta_cos_col]][params$group == g]
        bs <- params[[beta_sin_col]][params$group == g]
        ok <- complete.cases(bc, bs)
        bc_ok <- bc[ok]; bs_ok <- bs[ok]
        # Amplitude-weighted circular mean via vector averaging
        x_h <- sqrt(bc_ok^2 + bs_ok^2) * cos(atan2(bs_ok, bc_ok))
        y_h <- sqrt(bc_ok^2 + bs_ok^2) * sin(atan2(bs_ok, bc_ok))
        acro_h <- atan2(mean(y_h), mean(x_h))
        if(acro_h < 0) acro_h <- acro_h + 2 * pi
        acro_time <- acro_h * effective_period / (2 * pi)
        amp_mean <- sqrt(mean(bc_ok)^2 + mean(bs_ok)^2)
        cat(sprintf("    %s: amplitude-weighted mean = %.2f h, mean amplitude = %.3f (n=%d)\n",
                    g, acro_time, amp_mean, sum(ok)))
      }
      
      cat("\n--- R-squared Comparison ---\n")
      if(n_groups == 2) {
        tt <- t.test(r_squared ~ group, data = params)
        cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
      } else {
        aov_rsq <- aov(r_squared ~ group, data = params)
        cat("ANOVA:\n")
        print(summary(aov_rsq))
      }
      
      # Trend parameter comparison based on trend type
      if(mod$trend_type == "linear" && "trend_linear" %in% names(params)) {
        cat("\n--- Linear Trend (β) Comparison ---\n")
        params_trend <- params[!is.na(params$trend_linear), ]
        
        if(nrow(params_trend) >= 4) {
          if(n_groups == 2) {
            tt <- t.test(trend_linear ~ group, data = params_trend)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            
            # Effect size (Cohen's d)
            groups <- unique(params_trend$group)
            g1 <- params_trend$trend_linear[params_trend$group == groups[1]]
            g2 <- params_trend$trend_linear[params_trend$group == groups[2]]
            pooled_sd <- sqrt(((length(g1)-1)*sd(g1)^2 + (length(g2)-1)*sd(g2)^2) / (length(g1)+length(g2)-2))
            cohens_d <- (mean(g1) - mean(g2)) / pooled_sd
            cat(sprintf("Cohen's d = %.3f\n", cohens_d))
          } else {
            aov_trend <- aov(trend_linear ~ group, data = params_trend)
            cat("ANOVA:\n")
            print(summary(aov_trend))
          }
          
          cat("\nGroup statistics (units/hour):\n")
          for(g in unique(params_trend$group)) {
            g_trend <- params_trend$trend_linear[params_trend$group == g]
            cat(sprintf("  %s: mean = %.4f, SD = %.4f (n=%d)\n", 
                        g, mean(g_trend, na.rm=TRUE), sd(g_trend, na.rm=TRUE), length(g_trend)))
          }
        }
        
      } else if(mod$trend_type == "log" && "trend_log" %in% names(params)) {
        cat("\n--- Logarithmic Trend (β) Comparison ---\n")
        params_trend <- params[!is.na(params$trend_log), ]
        
        if(nrow(params_trend) >= 4) {
          if(n_groups == 2) {
            tt <- t.test(trend_log ~ group, data = params_trend)
            cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
          } else {
            aov_trend <- aov(trend_log ~ group, data = params_trend)
            cat("ANOVA:\n")
            print(summary(aov_trend))
          }
          
          cat("\nGroup statistics (units/log-hour):\n")
          for(g in unique(params_trend$group)) {
            g_trend <- params_trend$trend_log[params_trend$group == g]
            cat(sprintf("  %s: mean = %.4f, SD = %.4f (n=%d)\n", 
                        g, mean(g_trend, na.rm=TRUE), sd(g_trend, na.rm=TRUE), length(g_trend)))
          }
        }
        
      } else if(mod$trend_type == "exp_sat") {
        # Compare A_sat (asymptote)
        if("A_sat" %in% names(params)) {
          cat("\n--- Asymptote (A_sat) Comparison ---\n")
          params_trend <- params[!is.na(params$A_sat), ]
          
          if(nrow(params_trend) >= 4) {
            if(n_groups == 2) {
              tt <- t.test(A_sat ~ group, data = params_trend)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_asat <- aov(A_sat ~ group, data = params_trend)
              cat("ANOVA:\n")
              print(summary(aov_asat))
            }
            
            cat("\nGroup statistics (units):\n")
            for(g in unique(params_trend$group)) {
              g_asat <- params_trend$A_sat[params_trend$group == g]
              cat(sprintf("  %s: mean = %.3f, SD = %.3f (n=%d)\n", 
                          g, mean(g_asat, na.rm=TRUE), sd(g_asat, na.rm=TRUE), length(g_asat)))
            }
          }
        }
        
        # Compare tau (time constant)
        if("tau" %in% names(params)) {
          cat("\n--- Time Constant (τ) Comparison ---\n")
          params_trend <- params[!is.na(params$tau), ]
          
          if(nrow(params_trend) >= 4) {
            if(n_groups == 2) {
              tt <- t.test(tau ~ group, data = params_trend)
              cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.4f\n", tt$statistic, tt$parameter, tt$p.value))
            } else {
              aov_tau <- aov(tau ~ group, data = params_trend)
              cat("ANOVA:\n")
              print(summary(aov_tau))
            }
            
            cat("\nGroup statistics (hours):\n")
            for(g in unique(params_trend$group)) {
              g_tau <- params_trend$tau[params_trend$group == g]
              cat(sprintf("  %s: mean τ = %.2f h, SD = %.2f (n=%d)\n", 
                          g, mean(g_tau, na.rm=TRUE), sd(g_tau, na.rm=TRUE), length(g_tau)))
            }
            cat("\nNote: τ represents time to reach ~63% of asymptotic level.\n")
          }
        }
      }
    }
  })
  
  # --- Exports ---
  output$export_scores_csv <- downloadHandler(
    filename = function() paste0("beta_coefficients_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$reg_model)
      write.csv(values$reg_model$beta.hat, file)
    }
  )
  
  output$export_r_code <- downloadHandler(
    filename = function() paste0("reproduce_fosr_analysis.R"),
    content = function(file) {
      code <- "# R Code Generation logic here (omitted for brevity in single file app)"
      writeLines(code, file)
    }
  )
  
  # Harmonic regression exports
  output$export_harmonic_params <- downloadHandler(
    filename = function() paste0("harmonic_individual_params_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$harmonic_model)
      write.csv(values$harmonic_model$individual_params, file, row.names = FALSE)
    }
  )
  
  output$export_harmonic_summary <- downloadHandler(
    filename = function() paste0("harmonic_summary_", Sys.Date(), ".txt"),
    content = function(file) {
      req(values$harmonic_model)
      mod <- values$harmonic_model
      
      sink(file)
      cat("=== Harmonic Regression (Cosinor Analysis) Summary ===\n")
      cat("Generated:", as.character(Sys.time()), "\n\n")
      cat("Period:", mod$period, "\n")
      cat("Number of harmonics:", mod$n_harmonics, "\n")
      cat("Number of subjects:", length(mod$individual_fits), "\n\n")
      
      if(!is.null(mod$pop_mean_fit)) {
        pop <- mod$pop_mean_fit
        cat("--- Population Mean Parameters ---\n")
        cat(sprintf("MESOR:     %.4f\n", pop$mean_mesor))
        cat(sprintf("Amplitude: %.4f\n", pop$mean_amplitude))
        cat(sprintf("Acrophase: %.4f radians (%.2f hours)\n", 
                    pop$mean_acrophase_rad, pop$mean_acrophase_time))
        cat(sprintf("Rayleigh Z: %.4f (p = %.4f)\n", pop$rayleigh_z, pop$rayleigh_p))
        cat("\n")
      }
      
      cat("--- Individual Parameter Statistics ---\n")
      params <- mod$individual_params
      cat(sprintf("MESOR:     Mean=%.4f, SD=%.4f\n", mean(params$mesor), sd(params$mesor)))
      cat(sprintf("Amplitude (H1): Mean=%.4f, SD=%.4f\n", mean(params$amplitude_1), sd(params$amplitude_1)))
      cat(sprintf("Acrophase (H1): Mean=%.2f hours\n", mean(params$acrophase_time_1)))
      cat(sprintf("R-squared: Mean=%.4f, Range=[%.4f, %.4f]\n", 
                  mean(params$r_squared), min(params$r_squared), max(params$r_squared)))
      cat(sprintf("Significant rhythms (p<0.05): %d / %d\n",
                  sum(params$p_value < 0.05), nrow(params)))
      
      if(!is.null(mod$boot_results)) {
        cat(sprintf("\n--- Bootstrap 95%% CIs (B=%d) ---\n", mod$boot_results$B))
        cat(sprintf("MESOR:     [%.4f, %.4f]\n", 
                    mod$boot_results$mesor_ci[1], mod$boot_results$mesor_ci[2]))
        cat(sprintf("Amplitude: [%.4f, %.4f]\n", 
                    mod$boot_results$amplitude_ci[1], mod$boot_results$amplitude_ci[2]))
        cat(sprintf("Acrophase: [%.2f, %.2f] hours\n", 
                    mod$boot_results$acrophase_ci[1], mod$boot_results$acrophase_ci[2]))
      }
      sink()
    }
  )

  # ==============================================================================
  # PAIRWISE COMPARISONS MODULE
  # ==============================================================================

  # Run pairwise comparisons
  observeEvent(input$run_pairwise, {
    req(values$harmonic_model)
    mod <- values$harmonic_model

    # Check if groups are defined
    if(is.null(mod$group_fits) || length(mod$group_fits) < 2) {
      showNotification("Pairwise comparisons require 2 or more groups. Please define groups in Harmonic Regression tab.",
                      type = "error", duration = 5)
      return()
    }

    # Get parameter to compare
    param <- input$pairwise_param
    correction <- input$pairwise_correction

    # Get individual parameters with group information
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]
    } else {
      showNotification("No group variable selected. Please select groups in Harmonic Regression tab.",
                      type = "error", duration = 5)
      return()
    }

    # Check if parameter exists in data
    if(!param %in% names(params)) {
      showNotification(paste("Parameter", param, "not available in current model."),
                      type = "error", duration = 5)
      return()
    }

    # Get groups
    groups <- levels(params$group)
    n_groups <- length(groups)

    if(n_groups < 2) {
      showNotification("Need at least 2 groups for pairwise comparisons.",
                      type = "error", duration = 5)
      return()
    }

    # Perform pairwise t-tests
    n_comparisons <- choose(n_groups, 2)
    results <- data.frame(
      comparison = character(n_comparisons),
      group1 = character(n_comparisons),
      group2 = character(n_comparisons),
      mean1 = numeric(n_comparisons),
      sd1 = numeric(n_comparisons),
      n1 = integer(n_comparisons),
      mean2 = numeric(n_comparisons),
      sd2 = numeric(n_comparisons),
      n2 = integer(n_comparisons),
      mean_diff = numeric(n_comparisons),
      t_stat = numeric(n_comparisons),
      df = numeric(n_comparisons),
      p_value = numeric(n_comparisons),
      cohens_d = numeric(n_comparisons),
      ci_lower = numeric(n_comparisons),
      ci_upper = numeric(n_comparisons),
      stringsAsFactors = FALSE
    )

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    idx <- 1
    for(i in 1:(n_groups-1)) {
      for(j in (i+1):n_groups) {
        g1 <- groups[i]
        g2 <- groups[j]

        vals1 <- params[[param]][params$group == g1]
        vals2 <- params[[param]][params$group == g2]

        # Remove NAs
        vals1 <- vals1[!is.na(vals1)]
        vals2 <- vals2[!is.na(vals2)]

        if(length(vals1) < 2 || length(vals2) < 2) {
          next
        }

        if(is_circular) {
          # CIRCULAR STATISTICS for acrophase parameters
          # Convert hours to radians (assuming 24-hour period)
          period <- mod$period
          rad1 <- vals1 * 2 * pi / period
          rad2 <- vals2 * 2 * pi / period

          # Circular means (in radians)
          cmean1_rad <- circular_mean(rad1)
          cmean2_rad <- circular_mean(rad2)

          # Convert back to hours for display
          cmean1 <- cmean1_rad * period / (2 * pi)
          cmean2 <- cmean2_rad * period / (2 * pi)

          # Ensure positive (0-24h range)
          if(cmean1 < 0) cmean1 <- cmean1 + period
          if(cmean2 < 0) cmean2 <- cmean2 + period

          # Circular standard deviations (in radians, then convert to hours)
          csd1_rad <- circular_sd(rad1)
          csd2_rad <- circular_sd(rad2)
          csd1 <- csd1_rad * period / (2 * pi)
          csd2 <- csd2_rad * period / (2 * pi)

          # Angular difference (shortest arc)
          ang_diff_rad <- cmean1_rad - cmean2_rad
          # Normalize to [-pi, pi]
          ang_diff_rad <- atan2(sin(ang_diff_rad), cos(ang_diff_rad))
          ang_diff <- ang_diff_rad * period / (2 * pi)

          # Mean resultant lengths (measure of concentration)
          r1 <- mean_resultant_length(rad1)
          r2 <- mean_resultant_length(rad2)

          # Watson-Williams test for two groups
          ww <- watson_williams_test(list(rad1, rad2))

          # Effect size for circular data: difference in mean resultant lengths
          # (alternative: use V statistic, but difference in r is more interpretable)
          effect_size <- r1 - r2

          # For confidence interval on angular difference, use approximate circular CI
          # (simplified: not implemented here, set to NA)
          ci_lower <- NA
          ci_upper <- NA

          # Store results
          results$comparison[idx] <- paste(g1, "vs", g2)
          results$group1[idx] <- as.character(g1)
          results$group2[idx] <- as.character(g2)
          results$mean1[idx] <- cmean1
          results$sd1[idx] <- csd1
          results$n1[idx] <- length(vals1)
          results$mean2[idx] <- cmean2
          results$sd2[idx] <- csd2
          results$n2[idx] <- length(vals2)
          results$mean_diff[idx] <- ang_diff
          results$t_stat[idx] <- ww$F  # F-statistic from Watson-Williams
          results$df[idx] <- ww$df2  # Store df2 in df column
          results$p_value[idx] <- ww$p
          results$cohens_d[idx] <- effect_size  # Actually difference in mean resultant lengths
          results$ci_lower[idx] <- ci_lower
          results$ci_upper[idx] <- ci_upper

        } else {
          # REGULAR STATISTICS for non-circular parameters
          # Perform t-test
          t_result <- t.test(vals1, vals2, var.equal = FALSE)

          # Calculate Cohen's d
          pooled_sd <- sqrt(((length(vals1)-1)*sd(vals1)^2 + (length(vals2)-1)*sd(vals2)^2) /
                           (length(vals1) + length(vals2) - 2))
          cohens_d <- (mean(vals1) - mean(vals2)) / pooled_sd

          results$comparison[idx] <- paste(g1, "vs", g2)
          results$group1[idx] <- as.character(g1)
          results$group2[idx] <- as.character(g2)
          results$mean1[idx] <- mean(vals1)
          results$sd1[idx] <- sd(vals1)
          results$n1[idx] <- length(vals1)
          results$mean2[idx] <- mean(vals2)
          results$sd2[idx] <- sd(vals2)
          results$n2[idx] <- length(vals2)
          results$mean_diff[idx] <- mean(vals1) - mean(vals2)
          results$t_stat[idx] <- t_result$statistic
          results$df[idx] <- t_result$parameter
          results$p_value[idx] <- t_result$p.value
          results$cohens_d[idx] <- cohens_d
          results$ci_lower[idx] <- t_result$conf.int[1]
          results$ci_upper[idx] <- t_result$conf.int[2]
        }

        idx <- idx + 1
      }
    }

    # Remove empty rows
    results <- results[results$comparison != "", ]

    # Apply multiple comparison correction
    if(correction != "none") {
      results$p_adjusted <- p.adjust(results$p_value, method = correction)
    } else {
      results$p_adjusted <- results$p_value
    }

    # Store results
    values$pairwise_results <- results
    values$pairwise_param <- param
    values$pairwise_correction <- correction

    showNotification("Pairwise comparisons completed!", type = "message", duration = 3)
  })

  # Display pairwise results
  output$pairwise_results <- renderPrint({
    req(values$pairwise_results)
    results <- values$pairwise_results
    param <- values$pairwise_param
    correction <- values$pairwise_correction

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    cat("=== Pairwise Group Comparisons ===\n\n")
    cat("Parameter:", param, "\n")
    if(is_circular) {
      cat("Data type: Circular (using Watson-Williams test)\n")
    } else {
      cat("Data type: Linear (using Welch's t-test)\n")
    }
    cat("Correction method:", correction, "\n")
    cat("Number of comparisons:", nrow(results), "\n\n")

    for(i in 1:nrow(results)) {
      r <- results[i, ]
      cat("---\n")
      cat(sprintf("%s:\n", r$comparison))

      if(is_circular) {
        # Circular statistics display
        cat(sprintf("  Group 1: Circular mean=%.2f h, Circular SD=%.2f h, n=%d\n", r$mean1, r$sd1, r$n1))
        cat(sprintf("  Group 2: Circular mean=%.2f h, Circular SD=%.2f h, n=%d\n", r$mean2, r$sd2, r$n2))
        cat(sprintf("  Angular difference: %.2f h\n", r$mean_diff))

        # No CI for circular data (not implemented)
        # if(input$pairwise_show_ci && !is.na(r$ci_lower)) {
        #   cat(sprintf("  95%% CI: [%.3f, %.3f]\n", r$ci_lower, r$ci_upper))
        # }

        cat(sprintf("  Watson-Williams F(%d, %d) = %.3f, p = %.4f", 1, r$df, r$t_stat, r$p_value))

      } else {
        # Regular statistics display
        cat(sprintf("  Group 1: M=%.3f, SD=%.3f, n=%d\n", r$mean1, r$sd1, r$n1))
        cat(sprintf("  Group 2: M=%.3f, SD=%.3f, n=%d\n", r$mean2, r$sd2, r$n2))
        cat(sprintf("  Difference: %.3f\n", r$mean_diff))

        if(input$pairwise_show_ci) {
          cat(sprintf("  95%% CI: [%.3f, %.3f]\n", r$ci_lower, r$ci_upper))
        }

        cat(sprintf("  t(%.1f) = %.3f, p = %.4f", r$df, r$t_stat, r$p_value))
      }

      if(correction != "none") {
        cat(sprintf(", p_adj = %.4f", r$p_adjusted))
      }

      if(r$p_adjusted < 0.001) {
        cat(" ***")
      } else if(r$p_adjusted < 0.01) {
        cat(" **")
      } else if(r$p_adjusted < 0.05) {
        cat(" *")
      }
      cat("\n")

      if(input$pairwise_show_effect_size) {
        if(is_circular) {
          # For circular data: difference in mean resultant lengths
          cat(sprintf("  Δr̄ (difference in mean resultant length): %.3f", r$cohens_d))
          if(abs(r$cohens_d) < 0.1) {
            cat(" (small)")
          } else if(abs(r$cohens_d) < 0.3) {
            cat(" (medium)")
          } else {
            cat(" (large)")
          }
        } else {
          # For linear data: Cohen's d
          cat(sprintf("  Cohen's d: %.3f", r$cohens_d))
          if(abs(r$cohens_d) < 0.2) {
            cat(" (negligible)")
          } else if(abs(r$cohens_d) < 0.5) {
            cat(" (small)")
          } else if(abs(r$cohens_d) < 0.8) {
            cat(" (medium)")
          } else {
            cat(" (large)")
          }
        }
        cat("\n")
      }
    }

    cat("\n---\n")
    cat("Significance codes: *** p<0.001, ** p<0.01, * p<0.05\n")
    if(is_circular) {
      cat("\nNote: Circular means are in hours (0-24). Angular difference is the shortest arc.\n")
      cat("Effect size Δr̄ measures difference in concentration (mean resultant lengths).\n")
    }
  })

  # Pairwise comparison plot
  output$pairwise_plot <- renderPlot({
    req(values$pairwise_results, values$harmonic_model)

    param <- values$pairwise_param
    mod <- values$harmonic_model

    # Get individual parameters with group information
    if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
      group_var <- values$covariates[[input$harmonic_group_var]]
      params <- mod$individual_params
      params$group <- as.factor(group_var[params$subject])
      params <- params[!is.na(params$group), ]
    } else {
      return(NULL)
    }

    if(!param %in% names(params)) {
      return(NULL)
    }

    # Get parameter label
    param_labels <- c(
      "mesor" = "MESOR",
      "amplitude_1" = "H1 Amplitude",
      "amplitude_2" = "H2 Amplitude",
      "amplitude_3" = "H3 Amplitude",
      "acrophase_time_1" = "H1 Acrophase (hours)",
      "acrophase_time_2" = "H2 Acrophase (hours)",
      "acrophase_time_3" = "H3 Acrophase (hours)",
      "r_squared" = "R²",
      "A_sat" = "A_sat",
      "tau" = "τ (tau)",
      "percent_S" = "Process S (%)",
      "percent_C" = "Process C (%)"
    )
    param_label <- ifelse(param %in% names(param_labels), param_labels[param], param)

    # Create boxplot with violin overlay
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

    # Prepare data
    groups <- levels(params$group)
    n_groups <- length(groups)
    plot_data <- list()
    for(g in groups) {
      plot_data[[as.character(g)]] <- params[[param]][params$group == g]
    }

    # Boxplot
    boxplot(plot_data,
            main = paste("Group Comparison:", param_label),
            ylab = param_label,
            xlab = "Group",
            col = rainbow(n_groups, alpha = 0.3),
            border = rainbow(n_groups),
            notch = TRUE,
            las = 1,
            cex.axis = 1.2,
            cex.lab = 1.3,
            cex.main = 1.4)

    # Add individual points with jitter
    for(i in 1:n_groups) {
      g <- groups[i]
      vals <- params[[param]][params$group == g]
      vals <- vals[!is.na(vals)]
      points(jitter(rep(i, length(vals)), amount = 0.1), vals,
             col = rainbow(n_groups, alpha = 0.5)[i],
             pch = 19, cex = 0.8)
    }

    # Add group means
    for(i in 1:n_groups) {
      g <- groups[i]
      vals <- params[[param]][params$group == g]
      vals <- vals[!is.na(vals)]
      points(i, mean(vals), pch = 18, cex = 2.5, col = "black")
    }

    # Add significance brackets
    results <- values$pairwise_results
    y_max <- max(params[[param]], na.rm = TRUE)
    y_min <- min(params[[param]], na.rm = TRUE)
    y_range <- y_max - y_min

    sig_results <- results[results$p_adjusted < 0.05, ]
    if(nrow(sig_results) > 0) {
      bracket_y <- y_max + y_range * 0.05
      for(i in 1:min(nrow(sig_results), 5)) {  # Show max 5 brackets
        g1_idx <- which(groups == sig_results$group1[i])
        g2_idx <- which(groups == sig_results$group2[i])

        y_pos <- bracket_y + (i - 1) * y_range * 0.08

        # Draw bracket
        segments(g1_idx, y_pos, g2_idx, y_pos, lwd = 1.5)
        segments(g1_idx, y_pos, g1_idx, y_pos - y_range * 0.02, lwd = 1.5)
        segments(g2_idx, y_pos, g2_idx, y_pos - y_range * 0.02, lwd = 1.5)

        # Add significance stars
        p_val <- sig_results$p_adjusted[i]
        stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else "*"
        text((g1_idx + g2_idx) / 2, y_pos + y_range * 0.02, stars, cex = 1.2)
      }
    }

    # Add legend
    legend("topleft", legend = c("Mean", "Individual"),
           pch = c(18, 19), col = c("black", "gray"),
           pt.cex = c(2.5, 0.8), bty = "n")
  })

  # Dynamic help text for comparison matrix
  output$pairwise_matrix_help <- renderUI({
    req(values$pairwise_param)
    param <- values$pairwise_param
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    if(is_circular) {
      helpText("Lower triangle: p-values | Upper triangle: Δr̄ (difference in mean resultant length)")
    } else {
      helpText("Lower triangle: p-values | Upper triangle: Cohen's d (effect sizes)")
    }
  })

  # Pairwise comparison matrix
  output$pairwise_matrix <- renderTable({
    req(values$pairwise_results)
    results <- values$pairwise_results
    param <- values$pairwise_param

    # Detect if parameter is acrophase (circular data)
    is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

    # Get unique groups
    groups <- unique(c(results$group1, results$group2))
    n_groups <- length(groups)

    # Create matrix
    mat <- matrix("", nrow = n_groups, ncol = n_groups,
                  dimnames = list(groups, groups))

    for(i in 1:nrow(results)) {
      r <- results[i, ]
      g1_idx <- which(groups == r$group1)
      g2_idx <- which(groups == r$group2)

      # Lower triangle: p-values
      p_str <- sprintf("%.4f", r$p_adjusted)
      if(r$p_adjusted < 0.001) p_str <- paste0(p_str, " ***")
      else if(r$p_adjusted < 0.01) p_str <- paste0(p_str, " **")
      else if(r$p_adjusted < 0.05) p_str <- paste0(p_str, " *")
      mat[g2_idx, g1_idx] <- p_str

      # Upper triangle: effect sizes
      d_str <- sprintf("%.3f", r$cohens_d)
      mat[g1_idx, g2_idx] <- d_str
    }

    # Convert to data frame
    mat_df <- as.data.frame(mat)
    mat_df <- cbind(Group = rownames(mat_df), mat_df)
    mat_df
  }, rownames = FALSE, striped = TRUE, bordered = TRUE)

  # Export pairwise results
  output$export_pairwise_results <- downloadHandler(
    filename = function() paste0("pairwise_comparisons_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$pairwise_results)
      results <- values$pairwise_results
      param <- values$pairwise_param
      correction <- values$pairwise_correction
      is_circular <- grepl("acrophase_time", param, ignore.case = TRUE)

      # Create header with metadata
      con <- file(file, "w")
      writeLines(paste0("# Pairwise Comparisons: ", param), con)
      writeLines(paste0("# Date: ", Sys.Date()), con)
      if(is_circular) {
        writeLines("# Data type: Circular (Watson-Williams test)", con)
        writeLines("# Note: mean1/mean2 are circular means; sd1/sd2 are circular SDs", con)
        writeLines("# mean_diff is angular difference (shortest arc); cohens_d is Δr̄", con)
        writeLines("# t_stat is Watson-Williams F-statistic", con)
      } else {
        writeLines("# Data type: Linear (Welch's t-test)", con)
        writeLines("# Note: mean1/mean2 are arithmetic means; cohens_d is Cohen's d", con)
      }
      writeLines(paste0("# Correction method: ", correction), con)
      close(con)

      # Append results
      write.table(results, file, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)
    }
  )

  # Export pairwise plot
  output$export_pairwise_plot <- downloadHandler(
    filename = function() paste0("pairwise_plot_", Sys.Date(), ".png"),
    content = function(file) {
      req(values$pairwise_results, values$harmonic_model)
      png(file, width = 1200, height = 800, res = 120)

      # Recreate the plot (same code as renderPlot)
      param <- values$pairwise_param
      mod <- values$harmonic_model

      if(!is.null(input$harmonic_group_var) && input$harmonic_group_var != "_none_") {
        group_var <- values$covariates[[input$harmonic_group_var]]
        params <- mod$individual_params
        params$group <- as.factor(group_var[params$subject])
        params <- params[!is.na(params$group), ]

        if(param %in% names(params)) {
          param_labels <- c(
            "mesor" = "MESOR",
            "amplitude_1" = "H1 Amplitude",
            "amplitude_2" = "H2 Amplitude",
            "amplitude_3" = "H3 Amplitude",
            "acrophase_time_1" = "H1 Acrophase (hours)",
            "acrophase_time_2" = "H2 Acrophase (hours)",
            "acrophase_time_3" = "H3 Acrophase (hours)",
            "r_squared" = "R²",
            "A_sat" = "A_sat",
            "tau" = "τ (tau)",
            "percent_S" = "Process S (%)",
            "percent_C" = "Process C (%)"
          )
          param_label <- ifelse(param %in% names(param_labels), param_labels[param], param)

          groups <- levels(params$group)
          n_groups <- length(groups)
          plot_data <- list()
          for(g in groups) {
            plot_data[[as.character(g)]] <- params[[param]][params$group == g]
          }

          par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
          boxplot(plot_data,
                  main = paste("Group Comparison:", param_label),
                  ylab = param_label,
                  xlab = "Group",
                  col = rainbow(n_groups, alpha = 0.3),
                  border = rainbow(n_groups),
                  notch = TRUE,
                  las = 1,
                  cex.axis = 1.2,
                  cex.lab = 1.3,
                  cex.main = 1.4)

          for(i in 1:n_groups) {
            g <- groups[i]
            vals <- params[[param]][params$group == g]
            vals <- vals[!is.na(vals)]
            points(jitter(rep(i, length(vals)), amount = 0.1), vals,
                   col = rainbow(n_groups, alpha = 0.5)[i],
                   pch = 19, cex = 0.8)
            points(i, mean(vals), pch = 18, cex = 2.5, col = "black")
          }

          results <- values$pairwise_results
          y_max <- max(params[[param]], na.rm = TRUE)
          y_min <- min(params[[param]], na.rm = TRUE)
          y_range <- y_max - y_min

          sig_results <- results[results$p_adjusted < 0.05, ]
          if(nrow(sig_results) > 0) {
            bracket_y <- y_max + y_range * 0.05
            for(i in 1:min(nrow(sig_results), 5)) {
              g1_idx <- which(groups == sig_results$group1[i])
              g2_idx <- which(groups == sig_results$group2[i])
              y_pos <- bracket_y + (i - 1) * y_range * 0.08
              segments(g1_idx, y_pos, g2_idx, y_pos, lwd = 1.5)
              segments(g1_idx, y_pos, g1_idx, y_pos - y_range * 0.02, lwd = 1.5)
              segments(g2_idx, y_pos, g2_idx, y_pos - y_range * 0.02, lwd = 1.5)
              p_val <- sig_results$p_adjusted[i]
              stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else "*"
              text((g1_idx + g2_idx) / 2, y_pos + y_range * 0.02, stars, cex = 1.2)
            }
          }

          legend("topleft", legend = c("Mean", "Individual"),
                 pch = c(18, 19), col = c("black", "gray"),
                 pt.cex = c(2.5, 0.8), bty = "n")
        }
      }
      dev.off()
    }
  )
}

shinyApp(ui, server)
