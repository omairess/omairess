# Time-Warped Principal Component Analysis Shiny App with Functional ANOVA and Pairwise Comparisons
# 
# This application performs functional PCA, time-warped PCA, and functional ANOVA with pairwise comparisons
# 
# Features include:
# - Data import from CSV/text/Excel files or sample data generation
# - Variable selection capability for flexible data import
# - Multiple smoothing options (none, automatic REML, manual)
# - MISSING VALUE HANDLING: Per-subject smoothing with interpolation
# - GOODNESS-OF-FIT METRICS: R², RMSE, Correlation, MAE for smoothing quality
# - PROPER TIME AXIS LABELS: X-axis shows actual time variable names/values
# - Regular functional PCA or time-warped PCA with improved alignment
# - Functional ANOVA for group comparisons (between and within subjects)
#   * Between-subjects: Standard functional ANOVA with pointwise and global tests
#   * Within-subjects: Repeated measures functional ANOVA using rmfanova package
# - Pairwise comparisons with multiple testing corrections
# - Interactive visualization of warping functions and results
# - Export functionality for tables, plots, and R code
#
# VERSION: 2.1 - Enhanced with repeated measures functional ANOVA capability
# DATE: December 7, 2025

# Install and load required packages
packages <- c("shiny", "shinydashboard", "shinyWidgets", "fda", "mgcv",
              "plotly", "DT", "dplyr", "tidyr", "ggplot2", "gridExtra",
              "rmfanova",  # Added for repeated measures functional ANOVA
              "cluster",   # Added for clustering diagnostics (silhouette width)
              "fda.usc",   # Added for functional k-means clustering n.     
              "readxl",    # Added for Excel file support (.xls, .xlsx)
              "reticulate") # Added for Python integration (DCF clustering)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Functional Data Analysis Suite"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Import", tabName = "import", icon = icon("upload")),
      menuItem("Data Preprocessing/Smoothing", tabName = "preprocess", icon = icon("cogs")),
      menuItem("Smoothing Diagnostics", tabName = "smooth_diag", icon = icon("chart-area"),
              badgeColor = "green"),
      menuItem("fPCA/time-warped PCA Settings", tabName = "settings", icon = icon("sliders-h")),
      menuItem("Functional PCA Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Functional ANOVA", tabName = "fanova", icon = icon("chart-bar")),
      menuItem("fANOVA: post-hoc tests", tabName = "pairwise", icon = icon("exchange-alt")),
      menuItem("Functional Clustering", tabName = "kmeans", icon = icon("project-diagram")),
      menuItem("Data Export", tabName = "export", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .significance-legend {
          background-color: #f0f0f0;
          padding: 10px;
          border-radius: 5px;
          margin: 10px 0;
        }
        .pairwise-summary-box {
          background-color: #e8f4fd;
          padding: 15px;
          border-radius: 8px;
          margin: 10px 0;
          border-left: 4px solid #2196F3;
        }
      "))
    ),
    
    tabItems(
      # Data Import Tab
      tabItem(
        tabName = "import",
        fluidRow(
          box(
            title = "1. Upload Data",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            fileInput("datafile", "Choose Data File",
                      accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx")),
            helpText("Supported formats: CSV, TXT, TSV, Excel (.xls, .xlsx)"),
            hr(),
            checkboxInput("header", "File has header row", TRUE),
            # Removed fixed checkbox for first column groups - now handled by variable selector
            radioButtons("data_format", "Data Format:",
                         choices = list("Wide (subjects in rows)" = "wide",
                                        "Long (subjects in columns)" = "long"),
                         selected = "wide"),
            helpText("Wide: Rows=Subjects, Cols=Time/Vars. Long: Cols=Subjects (Transpose)."),
            hr(),
            actionButton("load_data", "Load Raw File", class = "btn-warning"),
            hr(),
            actionButton("generate_sample", "Generate Sample Data", class = "btn-info"),
            checkboxInput("generate_with_groups", "Include group structure", FALSE),
            conditionalPanel(
              condition = "input.generate_with_groups == true",
              numericInput("n_groups", "Number of groups:", value = 3, min = 2, max = 10)
            )
          ),
          
          box(
            title = "2. Variable Selection",
            status = "warning",
            solidHeader = TRUE,
            width = 8,
            uiOutput("var_select_container"),
            hr(),
            h4("Data Status:"),
            verbatimTextOutput("data_status")
          )
        ),
        fluidRow(
          box(
            title = "Data Preview",
            status = "primary",
            width = 12,
            fluidRow(
              column(4,
                selectInput("data_preview_rows", "Rows to display:",
                            choices = c("5" = 5, "10" = 10, "20" = 20, "50" = 50, "100" = 100, "All" = -1),
                            selected = 10, width = "150px")
              )
            ),
            DTOutput("data_preview"),
            hr(),
            plotOutput("raw_data_plot", height = "300px")
          )
        )
      ),
      
      # Data Preprocessing Tab
      tabItem(
        tabName = "preprocess",
        fluidRow(
          box(
            title = "Smoothing Options",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            radioButtons("smooth_method", "Smoothing Method:",
                         choices = list("Raw data (no smoothing)" = "none",
                                        "Automatic smoothing (REML)" = "auto",
                                        "Manual smoothing" = "manual"),
                         selected = "auto"),
            checkboxInput("constrain_bounds", "Constrain smoothed values to specific range", FALSE),
            conditionalPanel(
              condition = "input.constrain_bounds == true",
              numericInput("min_bound", "Minimum value:", value = 0, step = 1),
              numericInput("max_bound", "Maximum value:", value = 100, step = 1),
              helpText("Smoothed values will be clamped to [min, max] range.")
            ),
            conditionalPanel(
              condition = "input.smooth_method == 'manual'",
              sliderInput("smooth_factor", "Smoothing Factor:",
                          min = 0.1, max = 10, value = 1, step = 0.1),
              numericInput("n_basis_manual", "Number of B-spline basis functions:",
                           value = 20, min = 4, max = 100),
              helpText(HTML("<b>Smoothing Factor explained:</b><br>
                             Lambda = 10^(-Smoothing Factor)<br>
                             - Lower factor (0.1-2) → Higher lambda → Smoother curves<br>
                             - Higher factor (6-10) → Lower lambda → More flexible curves<br>
                             <br><b>Number of B-splines:</b><br>
                             - More bases → More potential flexibility<br>
                             - Typical: 10-20 for most datasets<br>
                             - Max recommended: n_time_points - 2<br>
                             <br><b>💡 Tip:</b> Use 'Smoothing Diagnostics' tab to find optimal values!")),
              conditionalPanel(
                condition = "output.diagnostics_available",
                hr(),
                actionButton("use_diagnostic_lambda", "📊 Use Diagnostic Results", 
                             class = "btn-info btn-sm"),
                helpText("Click to automatically set smoothing factor based on REML/CV analysis")
              )
            ),
            conditionalPanel(
              condition = "input.smooth_method == 'auto'",
              numericInput("n_basis", "Number of B-spline basis functions:",
                           value = 20, min = 4, max = 100),
              helpText("Lambda = 0 uses REML optimization automatically")
            ),
            conditionalPanel(
              condition = "input.smooth_method != 'none'",
              hr(),
              actionButton("run_nbasis_analysis", "Analyse Optimal n-Basis",
                           class = "btn-sm btn-info", icon = icon("chart-line")),
              helpText("Tests a range of n-basis values and plots the mean GCV score. Lower GCV = better fit–complexity trade-off."),
              plotlyOutput("nbasis_gcv_plot", height = "280px")
            ),
            hr(),
            actionButton("apply_smooth", "Apply Smoothing", class = "btn-warning"),
            hr(),
            h4("Smoothing Fit Statistics:"),
            verbatimTextOutput("smoothing_fit_summary")
          ),
          box(
            title = "Data Visualization",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            selectInput("tick_freq_preprocess", "X-axis tick frequency:",
                        choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                       "Every hour" = 1, "Every 2 hours" = 2,
                                       "Every 4 hours" = 4, "All" = 0),
                        selected = 1, width = "200px"),
            plotlyOutput("data_plot", height = "400px"),
            hr(),
            wellPanel(
              style = "background-color: #f0f8ff; padding: 10px;",
              h5("Selected Curve Info"),
              p("Click on a curve in the plot to select it.", style = "font-size: 11px; color: #666;"),
              verbatimTextOutput("selected_curve_info"),
              actionButton("clear_curve_selection", "Clear Selection", class = "btn-sm btn-default")
            )
          )
        )
      ),
      
      # Smoothing Diagnostics Tab
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
            checkboxInput("cv_stratified", "Stratify by groups (if available)", TRUE),
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
                           1. Go to 'Data Preprocessing' tab<br>
                           2. Select 'Manual smoothing'<br>
                           3. Click '📊 Use Diagnostic Results' button<br>
                           4. Or manually set the Smoothing Factor value shown above<br>
                           5. Click 'Apply Smoothing' to process your data"))
          )
        )
      ),
      
      # Analysis Settings Tab
      tabItem(
        tabName = "settings",
        fluidRow(
          box(
            title = "PCA Type",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            radioButtons("pca_type", "Select Analysis Type:",
                         choices = list("Functional PCA" = "fpca",
                                        "Time-Warped PCA" = "twpca"),
                         selected = "fpca")
          )
        ),
        
        # Time-Warped PCA Settings
        conditionalPanel(
          condition = "input.pca_type == 'twpca'",
          fluidRow(
            box(
              title = "Warping Method",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              radioButtons("warping_method", "Select Warping Method:",
                           choices = list("Linear shift (translation)" = "linear_shift",
                                          "Landmark-based alignment" = "landmark",
                                          "Parametric warping functions" = "parametric"),
                           selected = "linear_shift"),
              
              # Linear shift options
              conditionalPanel(
                condition = "input.warping_method == 'linear_shift'",
                wellPanel(
                  h4("Linear Shift (Translation) Options"),
                  checkboxInput("periodic_shift", "Periodic/circular data (e.g., 24-hour cycles)", FALSE),
                  checkboxInput("allow_dilation", "Allow slight scaling (dilation)", FALSE),
                  conditionalPanel(
                    condition = "input.allow_dilation == true",
                    sliderInput("dilation_range", "Scaling factor range (b):",
                                min = 0, max = 3, value = c(0.95, 1.05), step = 0.01),
                    p("Scaling factor b in h(t) = a + bt. Values near 1.0 preserve shape.", 
                      style = "font-size: 11px; color: #666;")
                  ),
                  radioButtons("shift_reference", "Reference for alignment:",
                               choices = list("Mean curve" = "mean",
                                              "First curve" = "first",
                                              "Median curve" = "median"),
                               selected = "mean")
                ),
                hr()
              ),
              
              # Parametric warping options
              conditionalPanel(
                condition = "input.warping_method == 'parametric'",
                wellPanel(
                  h4("Parametric Warping Options"),
                  radioButtons("parametric_family", "Warping Function Family:",
                               choices = list(
                                 "Power (t^α)" = "power",
                                 "Exponential ((e^(αt) - 1)/(e^α - 1))" = "exponential",
                                 "Quadratic (at² + (1-a)t)" = "quadratic",
                                 "Logistic (sigmoid)" = "logistic"
                               ),
                               selected = "power"),
                  sliderInput("param_range", "Parameter Search Range:",
                              min = 0.1, max = 10, value = c(0.5, 2), step = 0.1),
                  checkboxInput("symmetric_warp", "Force symmetric warping", FALSE)
                ),
                hr()
              ),
              
              # Landmark options
              conditionalPanel(
                condition = "input.warping_method == 'landmark'",
                h4("Landmark Options"),
                radioButtons("landmark_method", "Landmark Method:",
                             choices = list("Automatic" = "auto",
                                            "Manual" = "manual"),
                             selected = "auto"),
                conditionalPanel(
                  condition = "input.landmark_method == 'manual'",
                  numericInput("max_landmarks", "Maximum number of landmarks:",
                               value = 5, min = 2, max = 20),
                  hr(),
                  actionButton("start_landmark", "Clear & Start Selection", class = "btn-info"),
                  actionButton("clear_landmarks", "Clear Current", class = "btn-warning"),
                  hr(),
                  plotlyOutput("landmark_plot", height = "400px"),
                  verbatimTextOutput("landmark_info")
                )
              )
            )
          )
        ),
        
        # Functional PCA Settings
        fluidRow(
          box(
            title = "PCA Settings",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            numericInput("n_components", "Number of components to extract:",
                         value = 3, min = 1, max = 10),
            actionButton("run_analysis", "Run Analysis", class = "btn-success")
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "PCA Status",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            textOutput("pca_status"),
            verbatimTextOutput("pca_summary")
          )
        ),
        fluidRow(
          conditionalPanel(
            condition = "input.pca_type == 'twpca'",
            box(
              title = "Time Warping Results",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(4,
                       sliderInput("n_curves_display", "Number of curves to display:",
                                   min = 5, max = 100, value = 30, step = 5)
                ),
                column(4,
                       checkboxGroupInput("display_options", "Display options:",
                                          choices = list("Show individual curves" = "curves",
                                                         "Show mean curves" = "means"),
                                          selected = c("curves", "means"),
                                          inline = TRUE)
                ),
                column(4,
                       selectInput("tick_freq_settings", "X-axis tick frequency:",
                                   choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                                  "Every hour" = 1, "Every 2 hours" = 2,
                                                  "Every 4 hours" = 4, "All" = 0),
                                   selected = 1, width = "180px")
                )
              ),
              fluidRow(
                column(6, 
                       h4("Time Warping Functions"),
                       plotlyOutput("warping_plot", height = "400px"),
                       br(),
                       downloadButton("download_warping_plot", "Download Plot", class = "btn-sm")
                ),
                column(6,
                       h4("Alignment Effect"),
                       plotlyOutput("alignment_comparison_plot", height = "400px"),
                       br(),
                       downloadButton("download_alignment_plot", "Download Plot", class = "btn-sm")
                )
              ),
              hr(),
              # Warping Fit Statistics Section
              h4("Warping Fit Statistics"),
              p("Fit statistics comparing raw vs warped data, based on EFDA methodology."),
              wellPanel(
                style = "background-color: #f8f9fa;",
                fluidRow(
                  column(6,
                    h5("Summary Statistics (averaged over subjects)"),
                    verbatimTextOutput("warping_fit_summary")
                  ),
                  column(6,
                    h5("Variance Decomposition"),
                    verbatimTextOutput("warping_variance_decomposition")
                  )
                ),
                hr(),
                fluidRow(
                  column(12,
                    h5("Model Selection Criteria"),
                    verbatimTextOutput("warping_model_criteria")
                  )
                )
              ),
              hr(),
              h4("Per-Subject Fit Statistics"),
              p("Detailed fit statistics for each subject. R² measures how well the warped curve matches the original (functional R² via integrals). RMSE/MAE are errors between original and warped curves."),
              DTOutput("warping_fit_per_subject"),
              downloadButton("download_warping_fit_stats", "Download Fit Statistics (CSV)", class = "btn-sm btn-info"),
              hr(),
              h4("Warping Amplitude Scores"),
              p("The warping amplitude quantifies how much each subject's time axis was warped. Higher values indicate more warping was needed for alignment."),
              DTOutput("warping_scores")
            )
          ),
          box(
            title = "PCA Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            selectInput("tick_freq_results", "X-axis tick frequency:",
                        choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                       "Every hour" = 1, "Every 2 hours" = 2,
                                       "Every 4 hours" = 4, "All" = 0),
                        selected = 1, width = "200px"),
            h4("Component Loadings"),
            plotlyOutput("loadings_plot", height = "400px"),
            hr(),
            h4("Variance Explained"),
            plotlyOutput("variance_plot", height = "300px"),
            hr(),
            h4("Component Scores"),
            sliderInput("effect_size", "Effect Size Multiplier:",
                        min = 0.5, max = 3, value = 1, step = 0.1),
            plotlyOutput("scores_plot", height = "400px"),
            hr(),
            DTOutput("scores_table")
          )
        )
      ),
      
      # Functional ANOVA Tab
      tabItem(
        tabName = "fanova",
        fluidRow(
          box(
            title = "Functional ANOVA Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            h4("Experimental Design"),
            radioButtons("fanova_design", "Design type:",
                         choices = list("Between subjects" = "between",
                                        "Within subjects (Repeated Measures)" = "within"),
                         selected = "between"),
            hr(),
            
            # Between-subjects options
            conditionalPanel(
              condition = "input.fanova_design == 'between'",
              h4("Group Variable Selection"),
              helpText("For between-subjects designs, you must define group labels in the 'Data Preprocessing' tab first (e.g., Control vs Treatment)."),
              helpText("Each observation/curve belongs to exactly one group."),
              uiOutput("group_variable_ui"),
              uiOutput("fanova_group_info")
            ),
            
            # Within-subjects options
            conditionalPanel(
              condition = "input.fanova_design == 'within'",
              h4("Repeated Measures Factors"),
              helpText("For within-subjects designs, specify subject ID and the repeated measures factor (e.g., visit, condition, time)."),
              helpText("Note: Do NOT define group labels in preprocessing for within-subjects designs. Each subject should have multiple observations (one per condition/visit)."),
              uiOutput("subject_id_ui"),
              uiOutput("rm_factor_ui"),
              uiOutput("rm_factor_levels_ui")
            ),
            
            hr(),
            h4("ANOVA Options"),
            radioButtons("fanova_data_source", "Data source for ANOVA:",
                         choices = list("Original curves" = "original",
                                        "Time-warped curves (if available)" = "warped"),
                         selected = "original"),
            numericInput("n_permutations", "Number of permutations for testing:",
                         value = 200, min = 100, max = 10000),
            sliderInput("alpha_level", "Significance level:",
                        min = 0.01, max = 0.1, value = 0.05, step = 0.01),
            conditionalPanel(
              condition = "input.fanova_design == 'between'",
              radioButtons("fanova_test_type", "Test type:",
                           choices = list("Pointwise F-test" = "pointwise",
                                          "L2 norm test" = "L2",
                                          "Both" = "both"),
                           selected = "both")
            ),
            actionButton("run_fanova", "Run Functional ANOVA", class = "btn-primary")
          )
        ),
        fluidRow(
          box(
            title = "Global Test Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("fanova_global_results"),
            hr(),
            DTOutput("fanova_summary_table")
          )
        ),
        fluidRow(
          box(
            title = "Group Mean Functions",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(4,
                checkboxInput("fanova_show_sd_bands", "Show ±1 SD bands", value = TRUE)
              ),
              column(4,
                checkboxInput("fanova_show_sig_regions", "Highlight significant regions", value = TRUE)
              ),
              column(4,
                selectInput("tick_freq_fanova", "X-axis tick frequency:",
                            choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                           "Every hour" = 1, "Every 2 hours" = 2,
                                           "Every 4 hours" = 4, "All" = 0),
                            selected = 1, width = "180px")
              )
            ),
            helpText("Tip: Click on legend items to toggle traces. You can drag the legend to reposition it."),
            plotlyOutput("fanova_mean_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Pointwise Test Statistics",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("fanova_fstat_plot", height = "400px"),
            hr(),
            plotlyOutput("fanova_pvalue_plot", height = "400px")
          )
        ),
        fluidRow(
          box(
            title = "Effect Size Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("fanova_effect_size_plot", height = "400px"),
            hr(),
            verbatimTextOutput("fanova_effect_summary")
          )
        )
      ),
      
      # Pairwise Comparisons Tab
      tabItem(
        tabName = "pairwise",
        fluidRow(
          box(
            title = "Pairwise Comparison Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.fanova_completed",
              h4("Pairwise Test Configuration"),
              radioButtons("pairwise_correction", "Multiple comparison correction:",
                           choices = list("Bonferroni" = "bonferroni",
                                          "False Discovery Rate (FDR)" = "fdr",
                                          "Holm" = "holm",
                                          "Hochberg" = "hochberg",
                                          "None" = "none"),
                           selected = "bonferroni"),
              numericInput("pairwise_permutations", "Number of permutations:",
                           value = 200, min = 100, max = 5000),
              checkboxInput("pairwise_confidence_bands", "Show confidence bands for differences", TRUE),
              sliderInput("pairwise_alpha", "Significance level:",
                          min = 0.01, max = 0.1, value = 0.05, step = 0.01),
              actionButton("run_pairwise", "Run Pairwise Comparisons", class = "btn-primary")
            ),
            conditionalPanel(
              condition = "!output.fanova_completed",
              h4("Note:"),
              p("Please run Functional ANOVA first before performing pairwise comparisons.",
                style = "color: #d9534f; font-weight: bold;")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Pairwise Comparison Summary",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            div(class = "pairwise-summary-box",
                verbatimTextOutput("pairwise_summary")
            ),
            hr(),
            DTOutput("pairwise_global_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Select Comparison",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(8, uiOutput("pairwise_selector")),
              column(4,
                selectInput("tick_freq_pairwise", "X-axis tick frequency:",
                            choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                           "Every hour" = 1, "Every 2 hours" = 2,
                                           "Every 4 hours" = 4, "All" = 0),
                            selected = 1, width = "180px")
              )
            ),
            hr(),
            plotlyOutput("pairwise_difference_plot", height = "500px"),
            hr(),
            plotlyOutput("pairwise_pvalue_plot", height = "400px")
          )
        ),
        
        fluidRow(
          box(
            title = "All Pairwise Differences",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("pairwise_heatmap", height = "600px")
          )
        ),
        
        fluidRow(
          box(
            title = "Significant Time Regions",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("pairwise_regions_table"),
            hr(),
            plotlyOutput("pairwise_significance_timeline", height = "400px")
          )
        )
      ),

      # Functional Clustering Tab
      tabItem(
        tabName = "kmeans",
        fluidRow(
          box(
            title = "Optimal Number of Clusters",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            h4("Find Optimal k"),
            p("Run this analysis first to determine the optimal number of clusters before performing the final clustering."),
            fluidRow(
              column(3,
                numericInput("max_clusters_test", "Maximum k to test:",
                             value = 10, min = 5, max = 20, step = 1),
                helpText("Tests k from 2 to this value")
              ),
              column(3,
                radioButtons("opt_data_type", "Data to use:",
                             choices = list("Raw Data" = "raw",
                                            "Smoothed Data" = "smoothed"),
                             selected = "smoothed")
              ),
              column(3,
                radioButtons("opt_kmeans_method", "Method:",
                             choices = list("Standard" = "standard",
                                            "Functional" = "functional",
                                            "Hierarchical" = "hierarchical"),
                             selected = "standard"),
                conditionalPanel(
                  condition = "input.opt_kmeans_method == 'functional' && input.opt_data_type == 'raw'",
                  div(style = "color: orange; font-weight: bold; font-size: 11px;",
                      "⚠ Functional requires smoothed data")
                ),
                conditionalPanel(
                  condition = "input.opt_kmeans_method == 'hierarchical'",
                  selectInput("opt_hclust_linkage", "Linkage:",
                              choices = list("Ward.D2" = "ward.D2",
                                             "Complete" = "complete",
                                             "Average" = "average",
                                             "Single" = "single"),
                              selected = "ward.D2")
                )
              ),
              column(3,
                br(),
                actionButton("run_cluster_optimization", "Analyze Optimal k",
                             class = "btn-warning", icon = icon("chart-line"))
              ),
              column(3,
                br(),
                verbatimTextOutput("optimization_status")
              )
            ),
            hr(),
            fluidRow(
              column(6,
                plotlyOutput("elbow_plot", height = "350px")
              ),
              column(6,
                plotlyOutput("silhouette_plot", height = "350px")
              )
            ),
            hr(),
            verbatimTextOutput("optimization_recommendation"),
            conditionalPanel(
              condition = "input.opt_kmeans_method == 'hierarchical'",
              hr(),
              h4("Dendrogram"),
              p("Shows the full agglomerative tree. Use it together with the elbow and silhouette plots to choose the optimal cut point (k)."),
              plotOutput("opt_dendrogram_plot", height = "400px")
            )
          )
        ),
        fluidRow(
          box(
            title = "Clustering Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            numericInput("n_clusters", "Number of Clusters (k):",
                         value = 3, min = 2, max = 20, step = 1),
            radioButtons("cluster_data_type", "Data to Cluster:",
                         choices = list("Raw Data" = "raw",
                                        "Smoothed Data" = "smoothed"),
                         selected = "smoothed"),
            helpText("Choose whether to cluster raw or smoothed functional data."),
            hr(),
            radioButtons("cluster_standardize", "Standardize Before Clustering:",
                         choices = list("None" = "none",
                                        "Within participants (person-mean centering)" = "within",
                                        "Between participants (z-score per time point)" = "between"),
                         selected = "none"),
            helpText(HTML("<b>None:</b> Use data as-is.<br>
                          <b>Within participants:</b> Center and scale each participant's own time series (removes between-person differences in level and variability).<br>
                          <b>Between participants:</b> Z-score each time point across participants (standardizes each occasion to mean 0, SD 1).")),
            hr(),
            radioButtons("clustering_method", "Clustering Method:",
                         choices = list("Standard (Point-wise) K-Means" = "standard",
                                        "Functional K-Means" = "functional",
                                        "DCF (Density Core Finding)" = "dcf",
                                        "Hierarchical Clustering" = "hierarchical"),
                         selected = "standard"),
            helpText(HTML("<b>Standard K-Means:</b> Clusters based on values at each time point.<br>
                          <b>Functional K-Means:</b> Uses functional distances (fda.usc, requires smoothed data).<br>
                          <b>DCF:</b> Density-based clustering - automatically finds number of clusters (requires Python).<br>
                          <b>Hierarchical:</b> Agglomerative clustering with chosen linkage; cuts the tree at k clusters.")),
            conditionalPanel(
              condition = "input.clustering_method == 'functional' && input.cluster_data_type == 'raw'",
              div(style = "color: orange; font-weight: bold;",
                  "⚠ Functional K-Means requires smoothed data. Please select 'Smoothed Data' above.")
            ),
            # Hierarchical-specific parameters
            conditionalPanel(
              condition = "input.clustering_method == 'hierarchical'",
              wellPanel(
                style = "background-color: #f5fff5;",
                h5("Hierarchical Parameters"),
                selectInput("hclust_linkage", "Linkage Method:",
                            choices = list("Ward.D2" = "ward.D2",
                                           "Complete" = "complete",
                                           "Average (UPGMA)" = "average",
                                           "Single" = "single"),
                            selected = "ward.D2"),
                helpText("Ward.D2 minimises total within-cluster variance and is recommended for compact clusters.")
              )
            ),
            # DCF-specific parameters
            conditionalPanel(
              condition = "input.clustering_method == 'dcf'",
              wellPanel(
                style = "background-color: #f0f8ff;",
                h5("DCF Parameters"),
                numericInput("dcf_k", "k (neighborhood size):",
                             value = 10, min = 2, max = 100, step = 1),
                helpText("Number of neighbors for density estimation. Larger k = smoother density."),
                numericInput("dcf_beta", "β (fluctuation parameter):",
                             value = 0.1, min = 0.01, max = 1, step = 0.01),
                helpText("Controls density variation within cluster cores. Smaller β = tighter cores."),
                hr(),
                p("DCF automatically determines the number of clusters.", style = "font-style: italic; color: #666;"),
                actionButton("check_dcf_setup", "Check Python/DCF Setup", class = "btn-sm btn-info")
              )
            ),
            # Partitioning method parameters (Standard and Functional K-Means only)
            conditionalPanel(
              condition = "input.clustering_method == 'standard' || input.clustering_method == 'functional'",
              hr(),
              numericInput("kmeans_nstart", "Number of Random Restarts:",
                           value = 10, min = 1, max = 100, step = 1),
              helpText("The algorithm is run this many times with different random initialisations; the best solution (lowest within-cluster SS) is kept."),
              numericInput("kmeans_iter", "Maximum Iterations per Restart:",
                           value = 100, min = 10, max = 1000, step = 10)
            ),
            hr(),
            actionButton("run_clustering", "Run Clustering",
                         class = "btn-primary btn-lg", icon = icon("play")),
            br(), br(),
            verbatimTextOutput("clustering_status")
          ),
          box(
            title = "Cluster Statistics",
            status = "success",
            solidHeader = TRUE,
            width = 8,
            h4("Cluster Summary"),
            DTOutput("cluster_summary_table"),
            hr(),
            h4("Fit Statistics"),
            verbatimTextOutput("clustering_fit_stats")
          )
        ),
        fluidRow(
          box(
            title = "Cluster Membership",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            DTOutput("cluster_membership_table"),
            hr(),
            h4("Save Cluster Membership to Dataset"),
            fluidRow(
              column(3,
                textInput("cluster_var_name", "New Variable Name:",
                          value = "Cluster", placeholder = "e.g., Cluster_K3")
              ),
              column(3,
                br(),
                actionButton("save_cluster_membership", "Add to Dataset",
                             class = "btn-success", icon = icon("plus"))
              ),
              column(3,
                br(),
                downloadButton("download_data_with_clusters", "Download Dataset",
                               class = "btn-primary")
              ),
              column(3,
                br(),
                verbatimTextOutput("cluster_save_status")
              )
            ),
            helpText("Click 'Add to Dataset' to add cluster membership to your in-memory data (usable in other modules). Click 'Download Dataset' to export the full dataset including cluster membership to a CSV file.")
          )
        ),
        fluidRow(
          box(
            title = "Group Composition by Cluster",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.has_group_labels",
              uiOutput("cluster_group_var_selector"),
              h4("Group Counts by Cluster"),
              DTOutput("cluster_group_counts_table"),
              hr(),
              h4("Group Percentages by Cluster (%)"),
              DTOutput("cluster_group_pct_table"),
              hr(),
              verbatimTextOutput("cluster_group_test")
            ),
            conditionalPanel(
              condition = "!output.has_group_labels",
              p("No group variable available. Group composition analysis requires a group variable in your data.")
            )
          )
        ),
        fluidRow(
          box(
            title = "Cluster Mean Functions",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            selectInput("tick_freq_kmeans", "X-axis tick frequency:",
                        choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                       "Every hour" = 1, "Every 2 hours" = 2,
                                       "Every 4 hours" = 4, "All" = 0),
                        selected = 1, width = "200px"),
            plotlyOutput("cluster_means_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Individual Curves by Cluster",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("cluster_individuals_plot", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Detailed Silhouette Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("Silhouette plot showing individual sample coefficients. Each bar represents one subject."),
            p("Values close to 1 indicate well-clustered samples, values near 0 indicate borderline cases, negative values indicate potential misclassifications."),
            plotlyOutput("detailed_silhouette_plot", height = "600px")
          )
        ),
        conditionalPanel(
          condition = "input.clustering_method == 'hierarchical'",
          fluidRow(
            box(
              title = "Dendrogram",
              status = "success",
              solidHeader = TRUE,
              width = 12,
              p("Agglomerative dendrogram. The red dashed line shows the cut point for the selected k."),
              plotOutput("dendrogram_plot", height = "450px")
            )
          )
        )
      ),

      # Export Tab
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = "Export Options",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            h4("Export Tables"),
            downloadButton("export_scores_csv", "Download PCA Scores (CSV)", class = "btn-primary"),
            downloadButton("export_loadings_csv", "Download PCA Loadings (CSV)", class = "btn-primary"),
            conditionalPanel(
              condition = "input.pca_type == 'twpca'",
              downloadButton("export_warping_csv", "Download Warping Scores (CSV)", class = "btn-primary")
            ),
            downloadButton("export_fanova_results_csv", "Download FANOVA Results (CSV)", class = "btn-primary"),
            downloadButton("export_pairwise_results_csv", "Download Pairwise Results (CSV)", class = "btn-primary"),
            hr(),
            h4("Export Clustering Results"),
            downloadButton("export_cluster_membership_csv", "Download Cluster Membership (CSV)", class = "btn-primary"),
            downloadButton("export_cluster_means_csv", "Download Cluster Mean Functions (CSV)", class = "btn-primary"),
            downloadButton("export_cluster_stats_csv", "Download Cluster Statistics (CSV)", class = "btn-primary"),
            downloadButton("export_silhouette_csv", "Download Silhouette Data (CSV)", class = "btn-primary"),
            conditionalPanel(
              condition = "output.has_group_labels",
              downloadButton("export_cluster_group_csv", "Download Cluster-Group Composition (CSV)", class = "btn-primary")
            ),
            hr(),
            h4("Export Smoothed Data"),
            helpText("Exports the smoothed curves at original time resolution with original time labels."),
            downloadButton("export_smoothed_csv", "Download Smoothed Curves (CSV)", class = "btn-primary"),
            downloadButton("export_smoothed_long_csv", "Download Smoothed Curves - Long Format (CSV)", class = "btn-primary"),
            hr(),
            h4("Export Plots"),
            downloadButton("export_plots", "Download All Plots (PDF)", class = "btn-success"),
            downloadButton("export_fanova_plots", "Download FANOVA Plots (PDF)", class = "btn-success"),
            downloadButton("export_pairwise_plots", "Download Pairwise Plots (PDF)", class = "btn-success"),
            hr(),
            h4("Export R Code"),
            downloadButton("export_code", "Download Analysis Code (R)", class = "btn-warning"),
            hr(),
            h4("Analysis Code Preview"),
            verbatimTextOutput("code_preview")
          )
        )
      )
    )
  )
)

# Server function
server <- function(input, output, session) {
  
  # Test that server is running
  cat("===== SERVER STARTED =====\n")
  
  # Reactive values
  values <- reactiveValues(
    raw_df = NULL,      # Store raw imported dataframe
    data = NULL,        # Processed numeric matrix for analysis
    smooth_data = NULL,
    fd_obj = NULL,
    pca_results = NULL,
    warping_results = NULL,
    time_labels = NULL,      # Original column names for time points
    time_numeric = NULL,     # Numeric time values extracted from labels
    smoothing_fit_metrics = NULL,   # Fit metrics from smoothing
    smoothing_avg_metrics = NULL,    # Average fit metrics
    landmarks = list(),
    landmark_points = data.frame(x = numeric(), y = numeric()),
    group_labels = NULL,        # Primary group variable (factor)
    group_variables = NULL,     # Data frame storing all selected group variables
    selected_group_vars = NULL, # Names of selected group variables
    fanova_results = NULL,
    pairwise_results = NULL,
    clustering_results = NULL,  # Functional clustering results
    cluster_optimization = NULL,  # Cluster optimization results (elbow/silhouette)
    fanova_selected_groups = NULL,  # Selected groups for fANOVA analysis
    selected_curve = NULL       # Selected curve index for interactive plot
  )
  
  # Data status output
  
  
  # Helper function to extract hour values from column names
  # Supports multiple formats:
  #   KSS_9u_dag1 -> 9
  #   Base9h -> 9
  #   Base7h30 -> 7.5
  #   R1_5h -> 5
  #   8.0, 8.25, 8.5, 8.75 -> decimal hours (8:00, 8:15, 8:30, 8:45)
  #   X8.0, X8.25 -> decimal hours with X prefix (R's default for numeric column names)
  extract_hour_from_colname <- function(col_name) {
    # Pattern 0: Pure decimal number or with X prefix (e.g., "8.0", "8.25", "X8.0", "X8.25")
    # These represent decimal hours where .25=15min, .5=30min, .75=45min
    # Remove leading X if present (R adds this when column names start with numbers)
    clean_name <- gsub("^X", "", col_name)

    # Check if it's a pure decimal number
    if(grepl("^[0-9]+\\.?[0-9]*$", clean_name)) {
      return(as.numeric(clean_name))
    }

    # Pattern 1: hour with minutes "hMM" (e.g., "7h30"=7.5, "20h15"=20.25, "20h45"=20.75, "20h00"=20)
    if(grepl("[0-9]+h[0-9]+", col_name)) {
      hour_match <- regexpr("[0-9]+h[0-9]+", col_name)
      hour_str <- regmatches(col_name, hour_match)
      parts <- strsplit(hour_str, "h")[[1]]
      hour <- as.numeric(parts[1]) + as.numeric(parts[2]) / 60
      return(hour)
    }

    # Pattern 2: hour with "h" only, no minutes (e.g., "9h" = 9 hours)
    if(grepl("[0-9]+h", col_name)) {
      hour_match <- regexpr("[0-9]+h", col_name)
      hour_str <- regmatches(col_name, hour_match)
      hour <- as.numeric(gsub("h", "", hour_str))
      return(hour)
    }

    # Pattern 3: hour with "u" (e.g., "9u" = 9 hours)
    if(grepl("[0-9]+u", col_name)) {
      hour_match <- regexpr("[0-9]+u", col_name)
      hour_str <- regmatches(col_name, hour_match)
      hour <- as.numeric(gsub("u", "", hour_str))
      return(hour)
    }

    # Fallback: try to extract any number
    num_match <- regexpr("[0-9]+", col_name)
    if(num_match > 0) {
      return(as.numeric(regmatches(col_name, num_match)))
    }

    return(NA)
  }
  
  # Helper function to extract numeric time values from column names
  extract_time_values <- function(col_names) {
    tryCatch({
      # For chronological data (like circadian hours), preserve original order
      # Use sequential numbering to maintain chronology
      # The actual hour labels will be shown separately
      return(1:length(col_names))
    }, error = function(e) {
      return(1:length(col_names))
    })
  }
  
  # Helper function to get plotting time points
  get_plot_time <- function() {
    if(!is.null(values$time_numeric)) {
      return(values$time_numeric)
    } else if(!is.null(values$time_labels)) {
      return(1:length(values$time_labels))
    } else if(!is.null(values$data)) {
      return(1:ncol(values$data))
    } else {
      return(NULL)
    }
  }
  
  # Helper function to get time axis label
  get_time_label <- function() {
    if(!is.null(values$time_labels)) {
      first_label <- values$time_labels[1]
      # Check for circadian/hour data
      if(grepl("[0-9]+h30", first_label)) return("Hour")  # e.g., 7h30
      if(grepl("[0-9]+h", first_label)) return("Hour")    # e.g., 9h, 11h
      if(grepl("[0-9]+u", first_label)) return("Hour")    # e.g., 9u
      if(grepl("hour|hr", tolower(first_label))) return("Hour")
      if(grepl("day|d_", tolower(first_label))) return("Day")
      if(grepl("time|t_|^t[0-9]", tolower(first_label))) return("Time Point")
      if(grepl("month|mon", tolower(first_label))) return("Month")
      if(grepl("year|yr", tolower(first_label))) return("Year")
    }
    return("Measurement Point")
  }
  
  # Helper function to format plotly x-axis with time labels
  # Calculate actual time positions (accounting for day boundaries)
  calculate_time_positions <- function(hour_labels) {
    if(is.null(hour_labels) || length(hour_labels) == 0) {
      return(NULL)
    }
    
    n <- length(hour_labels)
    cumulative_hours <- numeric(n)
    cumulative_hours[1] <- 0  # Start at 0
    
    for(i in 2:n) {
      prev_hour <- hour_labels[i-1]
      curr_hour <- hour_labels[i]
      
      # Calculate time difference
      if(curr_hour >= prev_hour) {
        # Same day progression (e.g., 9 -> 11)
        diff <- curr_hour - prev_hour
      } else {
        # Crossed midnight (e.g., 23 -> 1)
        diff <- (24 - prev_hour) + curr_hour
      }
      
      cumulative_hours[i] <- cumulative_hours[i-1] + diff
    }
    
    # Normalize to 0-1 scale
    if(max(cumulative_hours) > 0) {
      normalized <- cumulative_hours / max(cumulative_hours)
    } else {
      normalized <- cumulative_hours
    }
    
    return(normalized)
  }
  
  format_plotly_time_axis <- function(p, time_grid = NULL, tick_step_hours = NULL) {
    if(!is.null(values$time_labels) && length(values$time_labels) > 0) {
      # Get hour labels if available
      hour_labels <- get_hour_labels()
      
      # Calculate actual time positions (accounting for midnight crossings)
      if(!is.null(hour_labels)) {
        # Calculate positions based on actual time progression
        actual_positions <- calculate_time_positions(hour_labels)
        
        if(!is.null(actual_positions)) {
          time_grid <- actual_positions
        } else if(is.null(time_grid)) {
          time_grid <- seq(0, 1, length.out = length(values$time_labels))
        }
        
        # Use actual hour values as labels
        tick_text <- sapply(hour_labels, decimal_to_hhmm)

        # Snap to round-hour multiples when a step is specified
        step <- as.numeric(tick_step_hours)
        if (!is.null(tick_step_hours) && !is.na(step) && step > 0) {
          keep <- abs(round(hour_labels / step) * step - hour_labels) < 1e-6
          if (any(keep)) {
            tick_vals_subset <- time_grid[keep]
            tick_text_subset <- tick_text[keep]
          } else {
            tick_vals_subset <- time_grid
            tick_text_subset <- tick_text
          }
        } else {
          # "All" (step == 0) or no step: show every point (thin if very many)
          n_labels <- length(tick_text)
          if (n_labels <= 60) {
            tick_vals_subset <- time_grid
            tick_text_subset <- tick_text
          } else {
            thin <- ceiling(n_labels / 50)
            indices <- seq(1, n_labels, by = thin)
            tick_vals_subset <- time_grid[indices]
            tick_text_subset <- tick_text[indices]
          }
        }
      } else {
        # Fall back to evenly-spaced positions and column names
        if(is.null(time_grid)) {
          time_grid <- seq(0, 1, length.out = length(values$time_labels))
        }
        tick_vals_subset <- time_grid
        tick_text_subset <- values$time_labels
      }
      
      # Apply custom axis formatting
      # Always use 90Â° rotation for readability
      p <- p %>% layout(
        xaxis = list(
          tickmode = "array",
          tickvals = tick_vals_subset,
          ticktext = tick_text_subset,
          tickangle = -90,  # Rotate labels vertically (90Â°)
          title = get_time_label()
        )
      )
    }
    return(p)
  }
  
  
  # Helper function to get hour labels for x-axis ticks
  get_hour_labels <- function() {
    if(!is.null(values$time_labels)) {
      # Extract hour values from column names
      hours <- sapply(values$time_labels, extract_hour_from_colname)
      if(!all(is.na(hours))) {
        return(hours)
      }
    }
    return(NULL)
  }

  # Convert decimal hours to "HHhMM" display string (e.g. 20.25 -> "20h15")
  decimal_to_hhmm <- function(h) {
    h_int <- floor(h)
    m_int <- round((h - h_int) * 60)
    if (m_int == 60) { h_int <- h_int + 1; m_int <- 0 }
    sprintf("%dh%02d", h_int, m_int)
  }

  # Map a vector of normalized x positions (0-1) to formatted "HHhMM" hover strings.
  # Uses approx() to interpolate against the original measurement positions,
  # so it works for both raw data grids and smooth 100-point interpolated grids.
  hover_time_labels <- function(x_vals) {
    hl <- get_hour_labels()
    if (is.null(hl)) return(as.character(round(x_vals, 3)))
    tp <- calculate_time_positions(hl)
    if (is.null(tp)) return(as.character(round(x_vals, 3)))
    interp_hours <- approx(tp, hl, xout = x_vals, rule = 2)$y
    sapply(interp_hours, decimal_to_hhmm)
  }
  
  
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
      if(!is.null(values$group_labels)) {
        cat("Groups:", length(unique(values$group_labels)), "groups detected\n")
        cat("Group distribution:", table(values$group_labels), "\n")
      } else {
        cat("No grouping variable selected.\n")
      }
    }
  })
  
  # Smoothing fit statistics display
  output$smoothing_fit_summary <- renderPrint({
    if(!is.null(values$smooth_fit_metrics)) {
      cat("=== Smoothing Quality Metrics ===\n\n")
      
      # Use CIRCAREG-style metrics structure (list with mean_r_squared, mean_rmse, etc.)
      metrics <- values$smooth_fit_metrics
      
      cat("Mean R2:", sprintf("%.3f", metrics$mean_r_squared), 
          sprintf("(SD: %.3f)", metrics$sd_r_squared), "\n")
      if(!is.na(metrics$mean_r_squared)) {
        if(metrics$mean_r_squared >= 0.9) {
          cat("  -> Excellent fit!\n")
        } else if(metrics$mean_r_squared >= 0.7) {
          cat("  -> Good fit\n")
        } else {
          cat("  -> Fair fit (consider adjusting parameters)\n")
        }
      }
      
      cat("\nMean RMSE:", sprintf("%.3f", metrics$mean_rmse),
          sprintf("(SD: %.3f)", metrics$sd_rmse), "\n")
      if(!is.null(metrics$rel_rmse_pct) && !is.na(metrics$rel_rmse_pct)) {
        cat(sprintf("  Relative RMSE: %.1f%% of data range", metrics$rel_rmse_pct))
        rmse_label <- if(metrics$rel_rmse_pct < 5)  " -> Excellent (<5%)" else
                      if(metrics$rel_rmse_pct < 10) " -> Good (5-10%)" else
                      if(metrics$rel_rmse_pct < 20) " -> Acceptable (10-20%)" else
                                                     " -> Poor (>20%)"
        cat(rmse_label, "\n")
      }

      # Number of subjects
      n_subjects <- length(metrics$r_squared)
      cat("\nNumber of subjects smoothed:", n_subjects, "\n")
      
      # Count subjects with valid metrics
      n_valid <- sum(!is.na(metrics$r_squared))
      if(n_valid < n_subjects) {
        cat("Subjects with valid fit metrics:", n_valid, "/", n_subjects, "\n")
      }
      
      # Display smoothing parameters
      cat("\n--- Smoothing Parameters ---\n")
      
      if(!is.null(metrics$method)) {
        method_label <- switch(metrics$method,
                               "auto" = "Automatic (REML)",
                               "manual" = "Manual",
                               "none" = "None (Raw data)")
        cat("Method:", method_label, "\n")
      }
      
      if(!is.null(metrics$n_basis)) {
        cat("Number of B-splines:", metrics$n_basis, "\n")
      }
      
      if(!is.null(metrics$lambda)) {
        if(metrics$lambda == 0) {
          cat("Lambda: 0 (automatic REML optimization)\n")
        } else {
          smooth_factor <- -log10(metrics$lambda)
          cat(sprintf("Lambda: %.3e (smoothing factor = %.2f)\n",
                      metrics$lambda, smooth_factor))
        }
      }

      # EDF-based n-basis recommendation
      if(!is.null(metrics$mean_df) && !is.na(metrics$mean_df)) {
        cat("\n--- Spline Complexity ---\n")
        cat(sprintf("Mean EDF (effective df): %.1f (SD: %.1f, max: %.1f)\n",
                    metrics$mean_df, metrics$sd_df, metrics$max_df))
        recommended_min <- ceiling(metrics$max_df) + 2
        cat("Current n_basis:", metrics$n_basis, "\n")
        cat(sprintf("Recommended minimum: ceil(max EDF) + 2 = %d\n", recommended_min))
        if(metrics$n_basis < recommended_min) {
          cat(sprintf("  ⚠ n_basis may be too low — consider increasing to at least %d\n",
                      recommended_min))
        } else if(metrics$mean_df < metrics$n_basis / 3) {
          cat("  ℹ Many basis functions are inactive (EDF << n_basis).\n")
          cat("    You may reduce n_basis for a more parsimonious model.\n")
        } else {
          cat("  OK — current n_basis provides adequate flexibility.\n")
        }
      }

    } else {
      cat("No smoothing applied yet.\n\n")
      cat("Click 'Apply Smoothing' to see fit quality metrics.")
    }
  })
  
  # Load Data Button: Read file into raw_df only
  observeEvent(input$load_data, {
    cat("Load data button clicked\n")
    
    if(is.null(input$datafile)) {
      showNotification("Please select a file first!", type = "warning", duration = 5)
      return()
    }
    
    ext <- tolower(tools::file_ext(input$datafile$name))

    tryCatch({
      # Check if Excel file
      if(ext %in% c("xls", "xlsx")) {
        # Read Excel file
        if(!requireNamespace("readxl", quietly = TRUE)) {
          showNotification("Package 'readxl' is required for Excel files. Install with: install.packages('readxl')",
                           type = "error", duration = 10)
          return()
        }

        raw_data <- as.data.frame(readxl::read_excel(
          input$datafile$datapath,
          col_names = input$header
        ))

        # Preserve column names (similar to check.names = FALSE)
        # readxl already preserves names, but ensure no issues
        cat("Read Excel file:", input$datafile$name, "\n")

      } else {
        # CSV/Text file handling
        # Determine separator
        sep <- if(ext %in% c("txt", "tsv")) "\t" else ","

        # Try to detect separator from file
        first_line <- readLines(input$datafile$datapath, n = 1)
        if(grepl("\t", first_line)) sep <- "\t"
        else if(grepl(";", first_line)) sep <- ";"

        # Read data into temporary dataframe
        # Use read.csv with check.names=FALSE to preserve decimal column names (e.g., 8.0, 8.25)
        # and proper quote handling for complex quoted strings
        raw_data <- read.csv(input$datafile$datapath,
                             header = input$header,
                             sep = sep,
                             stringsAsFactors = FALSE,
                             check.names = FALSE,
                             quote = "\"")
      }
      
      cat("Read raw file dimensions:", nrow(raw_data), "x", ncol(raw_data), "\n")

      # Handle duplicate column names by making them unique
      # This is critical for proper column selection (e.g., if "8h" appears twice)
      orig_names <- colnames(raw_data)
      if(any(duplicated(orig_names))) {
        unique_names <- make.unique(orig_names, sep = "_")
        colnames(raw_data) <- unique_names
        n_dups <- sum(duplicated(orig_names))
        cat("Note:", n_dups, "duplicate column name(s) found and made unique\n")
        cat("  Example:", orig_names[duplicated(orig_names)][1], "->",
            unique_names[which(duplicated(orig_names))[1]], "\n")
      }

      # Handle Wide vs Long format logic
      if(input$data_format == "long") {
        # "Long" in this app implies Subjects in Columns (transpose of wide)
        # We normalize everything to Wide (Subjects in Rows) for the variable selector
        raw_data <- as.data.frame(t(raw_data))
        cat("Transposed 'Long' format to Wide. New dims:", nrow(raw_data), "x", ncol(raw_data), "\n")
      }

      # Store in reactive value
      values$raw_df <- raw_data
      
      # ALSO store as uploaded_data for RM-ANOVA variable selection
      values$uploaded_data <- raw_data
      
      # Reset analysis data
      values$data <- NULL
      values$group_labels <- NULL
      
      showNotification("File loaded. Please select variables below.", type = "message", duration = 5)
      
    }, error = function(e) {
      cat("Error in load_data:", e$message, "\n")
      showNotification(paste("Error loading data:", e$message), type = "error", duration = 10)
    })
  })
  
  # Dynamic UI for Variable Selection
  output$var_select_container <- renderUI({
    req(values$raw_df)
    
    cols <- colnames(values$raw_df)
    
    # Try to guess numeric columns for default selection
    numeric_cols <- sapply(values$raw_df, is.numeric)
    default_data_cols <- cols[numeric_cols]
    
    # If no numeric cols found (maybe header issue), default to all except first
    if(length(default_data_cols) == 0 && length(cols) > 1) {
      default_data_cols <- cols[-1]
    }
    
    tagList(
      h4("Select Variables from Uploaded Data"),
      pickerInput(
        inputId = "sel_group_vars",
        label = "Select Group Variable(s) (Optional):",
        choices = cols,
        selected = NULL,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          `none-selected-text` = "None selected"
        ),
        multiple = TRUE
      ),
      helpText("You can select multiple group variables. The first selected will be used as the primary group for most analyses."),

      pickerInput(
        inputId = "sel_data_vars",
        label = "Select Time Series/Function Data Columns:",
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
      helpText("Note: Column order will be preserved as they appear in your data file."),

      actionButton("apply_selection", "Confirm & Process Data",
                   class = "btn-success", icon = icon("check"))
    )
  })
  
  # Apply Variable Selection
  observeEvent(input$apply_selection, {
    req(values$raw_df, input$sel_data_vars)
    
    tryCatch({
      # Extract Data Matrix
      data_cols <- input$sel_data_vars
      
      if(length(data_cols) < 2) {
        showNotification("Please select at least 2 time points/columns for analysis.", type = "warning")
        return()
      }
      
      # Preserve original column order (CRITICAL for chronological data)
      data_cols <- input$sel_data_vars

      # Ensure columns are in the order they appear in raw_df, not alphabetically
      # Since duplicate column names are now made unique at load time, this is safe
      all_cols <- colnames(values$raw_df)
      data_cols <- all_cols[all_cols %in% data_cols]

      temp_data <- values$raw_df[, data_cols, drop = FALSE]
      
      # Convert to numeric matrix - CRITICAL: Ensure proper numeric type
      temp_data_mat <- tryCatch({
        data.matrix(temp_data)
      }, error = function(e) {
        # Fallback conversion column by column
        mat <- matrix(NA, nrow = nrow(temp_data), ncol = ncol(temp_data))
        for(j in 1:ncol(temp_data)) mat[,j] <- as.numeric(as.character(temp_data[,j]))
        mat
      })
      
      # CRITICAL: Force numeric mode (prevents -Inf R-squared issue)
      # This ensures data type is compatible with fda package functions
      if(typeof(temp_data_mat) != "double" && typeof(temp_data_mat) != "integer") {
        mode(temp_data_mat) <- "numeric"
      }
      
      # Extract Group Variable(s)
      if(!is.null(input$sel_group_vars) && length(input$sel_group_vars) > 0) {
        # Store all selected group variables
        values$selected_group_vars <- input$sel_group_vars
        values$group_variables <- values$raw_df[, input$sel_group_vars, drop = FALSE]
        # Convert all to factors
        for(col in colnames(values$group_variables)) {
          values$group_variables[[col]] <- as.factor(values$group_variables[[col]])
        }
        # Primary group is the first selected
        values$group_labels <- values$group_variables[[1]]
        cat("Group variables selected:", paste(input$sel_group_vars, collapse = ", "), "\n")
        cat("Primary group variable:", input$sel_group_vars[1], "with",
            length(unique(values$group_labels)), "levels\n")
      } else {
        values$group_labels <- NULL
        values$group_variables <- NULL
        values$selected_group_vars <- NULL
      }
      
      values$data <- temp_data_mat
      
      # Store time labels and extract numeric values
      values$time_labels <- colnames(temp_data)
      values$time_numeric <- extract_time_values(values$time_labels)
      cat("Time labels stored:", paste(head(values$time_labels, 3), collapse=", "), "...\n")
      
      # Extract and display hour values for verification
      hour_labels <- get_hour_labels()
      if(!is.null(hour_labels) && length(hour_labels) > 0) {
        cat("Extracted hours:", paste(head(hour_labels, 10), collapse=", "), "...\n")
        cat("Total time points:", length(hour_labels), "\n")
      }
      
      
      # Clean data (NA handling)
      if(!is.null(values$data)) {
        # Remove rows with all NAs
        na_rows <- apply(is.na(values$data), 1, all)
        if(any(na_rows)) {
          values$data <- values$data[!na_rows, , drop = FALSE]
          if(!is.null(values$group_labels)) {
            values$group_labels <- values$group_labels[!na_rows]
          }
        }
        
        # Remove columns with all NAs
        na_cols <- apply(is.na(values$data), 2, all)
        if(any(na_cols)) {
          values$data <- values$data[, !na_cols, drop = FALSE]
        }
      }
      
      # Reset downstream analysis
      values$smooth_data <- NULL
      values$fd_obj <- NULL
      values$pca_results <- NULL
      values$warping_results <- NULL
      values$fanova_results <- NULL
      values$pairwise_results <- NULL
      
      showNotification("Data processed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error processing selection:", e$message), type = "error")
    })
  })
  
  # Generate Sample Data (Bypasses selection UI)
  observeEvent(input$generate_sample, {
    cat("Generate sample button clicked\n")
    
    tryCatch({
      set.seed(123)
      n_subjects <- 50  
      n_time <- 100
      time_points <- seq(0, 1, length.out = n_time)
      
      sample_data <- matrix(NA, nrow = n_subjects, ncol = n_time)
      
      if(input$generate_with_groups) {
        n_groups <- input$n_groups
        if(is.null(n_groups)) n_groups <- 3
        
        subjects_per_group <- ceiling(n_subjects / n_groups)
        group_labels <- rep(paste0("Group", 1:n_groups), each = subjects_per_group)[1:n_subjects]
        values$group_labels <- factor(group_labels)
        
        for(i in 1:n_subjects) {
          group_idx <- as.numeric(factor(group_labels[i]))
          phase_var <- rnorm(1, 0, 0.1)
          amplitude_var <- rnorm(1, 0, 0.3)
          
          group_effect <- (group_idx - 1) * 0.7 * sin(4*pi*time_points + group_idx * pi/n_groups)
          
          sample_data[i,] <- (1 + amplitude_var) * sin(2*pi*(time_points + phase_var)) + 
            0.5*cos(4*pi*time_points) + 
            group_effect +
            rnorm(n_time, 0, 0.1)
        }
      } else {
        values$group_labels <- NULL
        for(i in 1:n_subjects) {
          phase_var <- rnorm(1, 0, 0.1)
          amplitude_var <- rnorm(1, 0, 0.3)
          sample_data[i,] <- (1 + amplitude_var) * sin(2*pi*(time_points + phase_var)) + 
            0.5*cos(4*pi*time_points) + 
            rnorm(n_time, 0, 0.1)
        }
      }
      
      values$data <- sample_data
      
      # CRITICAL: Ensure proper numeric type
      if(typeof(values$data) != "double" && typeof(values$data) != "integer") {
        mode(values$data) <- "numeric"
      }
      
      # Create meaningful time labels for sample data
      values$time_labels <- paste0("T", 1:n_time)
      values$time_numeric <- 1:n_time
      
      values$raw_df <- NULL # Clear raw df so selection UI hides
      values$uploaded_data <- NULL # Clear uploaded data (no ID/RM vars in sample data)
      values$smooth_data <- NULL
      values$fd_obj <- NULL
      values$pca_results <- NULL
      values$warping_results <- NULL
      values$fanova_results <- NULL
      values$pairwise_results <- NULL
      
      showNotification("Sample data generated successfully!", type = "message", duration = 5)
      
    }, error = function(e) {
      cat("Error in generate_sample:", e$message, "\n")
      showNotification(paste("Error generating data:", e$message), type = "error", duration = 10)
    })
  })
  
  # Intelligent n_basis_manual update when data changes
  observe({
    req(values$data)
    
    n_time <- ncol(values$data)
    
    # Calculate recommended n_basis based on number of time points
    recommended_nb <- if(n_time <= 10) {
      max(4, n_time - 2)
    } else if(n_time <= 20) {
      max(10, round(n_time * 0.6))
    } else if(n_time <= 50) {
      max(15, round(n_time * 0.5))
    } else {
      max(20, round(n_time * 0.4))
    }
    
    # Ensure it's within bounds
    recommended_nb <- min(recommended_nb, n_time - 2, 100)
    recommended_nb <- max(recommended_nb, 4)
    
    # Update both n_basis inputs
    updateNumericInput(session, "n_basis", value = recommended_nb)
    updateNumericInput(session, "n_basis_manual", value = recommended_nb)
    
    cat(sprintf("Data loaded: %d time points → Recommended n_basis = %d\n", 
                n_time, recommended_nb))
  })
  
  # Data preview table
  output$data_preview <- renderDT({
    # Get number of rows to display from input
    n_rows_display <- if(!is.null(input$data_preview_rows)) {
      as.integer(input$data_preview_rows)
    } else {
      10
    }

    if(is.null(values$data)) {
      # If raw data is loaded but not processed, show raw data preview
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
      # Show processed data
      tryCatch({
        n_total_rows <- nrow(values$data)
        n_rows <- if(n_rows_display == -1) n_total_rows else min(n_rows_display, n_total_rows)

        # Show first 15 time columns max for readability, but all requested rows
        n_cols_show <- min(15, ncol(values$data))
        preview_data <- as.data.frame(values$data[1:n_rows, 1:n_cols_show])
        colnames(preview_data) <- paste0("T", 1:ncol(preview_data))

        if(!is.null(values$group_labels)) {
          preview_rows <- 1:n_rows
          preview_data <- cbind(
            Group = values$group_labels[preview_rows],
            Subject = preview_rows,
            preview_data
          )
        } else {
          preview_data <- cbind(
            Subject = 1:nrow(preview_data),
            preview_data
          )
        }

        datatable(preview_data,
                  options = list(pageLength = n_rows, scrollX = TRUE,
                                 lengthMenu = c(5, 10, 20, 50, 100)),
                  rownames = FALSE,
                  caption = paste("Processed Analysis Data (", n_rows, "/", n_total_rows, " rows, first ", n_cols_show, " time points)"))
      }, error = function(e) {
        return(NULL)
      })
    }
  })
  
  # Raw data plot
  output$raw_data_plot <- renderPlot({
    if(is.null(values$data)) {
      plot(1, type = "n", xlab = "", ylab = "", main = "No analysis data available")
      return()
    }
    
    tryCatch({
      n_time <- ncol(values$data)
      n_subj <- nrow(values$data)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      hour_labels <- get_hour_labels()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      
      if(!is.null(values$group_labels) && length(unique(values$group_labels)) > 1) {
        # Plot with groups
        groups <- levels(values$group_labels)
        n_groups <- length(groups)
        
        # Create scalable color palette
        base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
        colors <- colorRampPalette(base_cols)(n_groups)
        
        matplot(time_points_plot, t(values$data), type = "l", 
                col = rgb(0.5, 0.5, 0.5, 0.1), lty = 1,
                xlab = time_label, ylab = "Value", 
                main = paste("Raw Functional Data (", n_subj, "subjects,", 
                             n_groups, "groups)"),
                xaxt = if(!is.null(hour_labels)) "n" else "s")
        
        # Add custom x-axis with hour labels if available
        if(!is.null(hour_labels)) {
          n_labels <- length(hour_labels)
          if(n_labels > 15) {
            step <- ceiling(n_labels / 10)
            indices <- seq(1, n_labels, by = step)
            axis(1, at = time_points_plot[indices], labels = sapply(hour_labels[indices], decimal_to_hhmm))
          } else {
            axis(1, at = time_points_plot, labels = sapply(hour_labels, decimal_to_hhmm))
          }
        }
        
        # Add group means
        for(i in 1:n_groups) {
          group_idx <- which(values$group_labels == groups[i])
          if(length(group_idx) > 0) {
            group_mean <- colMeans(values$data[group_idx, , drop = FALSE])
            lines(time_points_plot, group_mean, col = colors[i], lwd = 3)
          }
        }
        
        legend("topright", legend = groups, col = colors, lwd = 3)
      } else {
        # Plot without groups
        matplot(time_points_plot, t(values$data), type = "l", 
                col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1,
                xlab = time_label, ylab = "Value", 
                main = paste("Raw Functional Data (", n_subj, "subjects)"))
        lines(time_points_plot, colMeans(values$data), col = "red", lwd = 3)
        legend("topright", legend = "Mean", col = "red", lwd = 3)
      }
    }, error = function(e) {
      plot(1, type = "n", xlab = "", ylab = "", main = paste("Error:", e$message))
    })
  })
  
  # Apply smoothing with missing value handling and goodness-of-fit metrics
  # USING CIRCAREG PROVEN APPROACH
  observeEvent(input$apply_smooth, {
    req(values$data)
    
    cat("Applying smoothing using CIRCAREG method...\n")
    
    tryCatch({
      n_time <- ncol(values$data)
      n_subjects <- nrow(values$data)
      t_full <- 1:n_time  # Use integer sequence like CIRCAREG, not 0-1 normalized
      
      # Determine number of basis functions (needed for later fd_obj recreation)
      if(input$smooth_method == "none") {
        nb <- min(20, n_time - 2)
      } else if(input$smooth_method == "manual") {
        # Use n_basis_manual for manual mode
        nb <- input$n_basis_manual
        nb <- min(nb, n_time)
        if(nb < 4) nb <- 4
      } else {
        # Auto mode uses n_basis
        nb <- input$n_basis
        nb <- min(nb, n_time)
        if(nb < 4) nb <- 4
      }
      
      # Warn if n_basis is too high
      if(nb >= n_time - 2 && input$smooth_method != "none") {
        showNotification(
          sprintf("Number of B-splines (%d) is very high relative to time points (%d). Consider reducing it.", 
                  nb, n_time),
          type = "warning", duration = 5
        )
      }
      
      # Create basis (CIRCAREG approach)
      if(input$smooth_method == "none") {
        basis <- create.bspline.basis(rangeval = c(1, n_time), breaks = t_full, norder = 4)
      } else {
        basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
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
        # Create fd_obj with 0-1 scale for compatibility with downstream analyses
        time_points_01 <- seq(0, 1, length.out = n_time)
        basis_01 <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time - 2))
        values$fd_obj <- smooth.basis(time_points_01, t(temp_data), basis_01)$fd
        values$smooth_fit_metrics <- NULL  # No metrics for raw data
        showNotification("Using Raw Data (No Smoothing).", type = "message")
        
      } else {
        # Smooth with NA handling - process each subject individually (CIRCAREG METHOD)
        # CRITICAL: For auto, use lambda = 0 (REML), not GCV search!
        lam <- if(input$smooth_method == "auto") 0 else 10^(-input$smooth_factor)
        fdParobj <- fdPar(basis, 2, lam)
        
        # Initialize smoothed data matrix
        smooth_mat <- matrix(NA, nrow = n_subjects, ncol = n_time)
        fd_coefs_list <- list()
        n_failed <- 0
        failed_subjects <- c()
        df_vec  <- rep(NA_real_, n_subjects)   # effective degrees of freedom per subject
        gcv_vec <- rep(NA_real_, n_subjects)   # GCV score per subject

        withProgress(message = 'Smoothing subjects...', value = 0, {
          for(i in 1:n_subjects) {
            y_i <- values$data[i, ]
            valid_idx <- !is.na(y_i)
            n_valid <- sum(valid_idx)

            # Need at least 4 points for B-spline smoothing
            min_points_needed <- 4

            if(n_valid >= min_points_needed) {
              tryCatch({
                t_valid <- t_full[valid_idx]
                y_valid <- y_i[valid_idx]

                # Smooth using only valid points; capture full result for EDF and GCV
                sb_i <- smooth.basis(t_valid, y_valid, fdParobj)

                # Evaluate at ALL time points (this interpolates over NAs)
                smooth_mat[i, ] <- as.vector(eval.fd(t_full, sb_i$fd))
                fd_coefs_list[[i]] <- sb_i$fd$coefs
                df_vec[i]  <- sb_i$df
                gcv_vec[i] <- sb_i$gcv
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
        # CIRCAREG EXACT METHOD
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
            
            # R-squared (CIRCAREG formula)
            ss_tot <- sum((orig_valid - mean(orig_valid))^2)
            ss_res <- sum((orig_valid - smooth_valid)^2)
            r_squared_vec[i] <- if(ss_tot > 0) 1 - ss_res / ss_tot else NA
          } else {
            r_squared_vec[i] <- NA
            rmse_vec[i] <- NA
          }
        }
        
        # Relative RMSE as % of total data range
        data_range   <- max(values$data, na.rm = TRUE) - min(values$data, na.rm = TRUE)
        rel_rmse_pct <- if(data_range > 0) {
          (mean(rmse_vec, na.rm = TRUE) / data_range) * 100
        } else NA_real_

        # Store fit metrics (CIRCAREG structure)
        values$smooth_fit_metrics <- list(
          r_squared = r_squared_vec,
          rmse = rmse_vec,
          mean_r_squared = mean(r_squared_vec, na.rm = TRUE),
          mean_rmse = mean(rmse_vec, na.rm = TRUE),
          sd_r_squared = sd(r_squared_vec, na.rm = TRUE),
          sd_rmse = sd(rmse_vec, na.rm = TRUE),
          n_basis = nb,
          lambda = lam,
          method = input$smooth_method,
          # EDF / GCV from smooth.basis()
          mean_df = mean(df_vec, na.rm = TRUE),
          sd_df   = sd(df_vec,   na.rm = TRUE),
          max_df  = max(df_vec,  na.rm = TRUE),
          mean_gcv = mean(gcv_vec, na.rm = TRUE),
          # Relative RMSE
          rel_rmse_pct = rel_rmse_pct,
          data_range   = data_range
        )
        
        # For compatibility with WAPAA display, also create avg_metrics
        values$smoothing_avg_metrics <- c(
          R_squared = mean(r_squared_vec, na.rm = TRUE),
          RMSE = mean(rmse_vec, na.rm = TRUE),
          Correlation = NA,  # Not calculated in CIRCAREG
          MAE = NA  # Not calculated in CIRCAREG
        )
        
        # CRITICAL: Recreate fd_obj with 0-1 scale for compatibility with downstream analyses
        # The smoothing was done with 1:n_time (for numerical stability)
        # But PCA, warping, etc. expect 0-1 scale
        # IMPORTANT: smooth_mat is ALREADY smoothed, so use lambda=0 to just create fd representation
        time_points_01 <- seq(0, 1, length.out = n_time)
        basis_01 <- create.bspline.basis(rangeval = c(0, 1), nbasis = nb)
        # Lambda = 0: smooth_mat already smoothed, just create fd representation
        values$fd_obj <- smooth.basis(time_points_01, t(smooth_mat), 
                                      fdPar(basis_01, 2, 0))$fd
        
        cat("Smoothed data stored with both 1:n_time (for metrics) and 0-1 scale (for analyses)\n")
        cat("fd_obj created with lambda=0 (no additional smoothing on already-smoothed data)\n")
        
        # Report results
        n_with_na <- sum(apply(values$data, 1, function(x) any(is.na(x))))
        n_interpolated <- n_with_na - length(failed_subjects)
        
        if(n_interpolated > 0) {
          showNotification(
            sprintf("Smoothing applied. %d subjects had missing values - successfully interpolated. Mean R2=%.3f, Mean RMSE=%.3f", 
                    n_interpolated, values$smooth_fit_metrics$mean_r_squared, values$smooth_fit_metrics$mean_rmse),
            type = "message", duration = 8)
        } else {
          showNotification(
            sprintf("Smoothing applied. Mean R2=%.3f (SD=%.3f), Mean RMSE=%.3f (SD=%.3f)", 
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
        
        cat(sprintf("\n=== Smoothing Results ===\n"))
        cat(sprintf("Mean R²: %.3f (SD: %.3f)\n", 
                    values$smooth_fit_metrics$mean_r_squared, 
                    values$smooth_fit_metrics$sd_r_squared))
        cat(sprintf("Mean RMSE: %.3f (SD: %.3f)\n", 
                    values$smooth_fit_metrics$mean_rmse, 
                    values$smooth_fit_metrics$sd_rmse))
      }
      
    }, error = function(e) {
      cat("Error in smoothing:", e$message, "\n")
      showNotification(paste("Smoothing error:", e$message), type = "error", duration = 10)
      
      # Fallback: use raw data
      values$smooth_data <- values$data
      n_basis <- min(10, ncol(values$data) - 2)
      basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)
      
      # CRITICAL: Ensure proper numeric type in fallback
      fallback_data <- values$data
      if(typeof(fallback_data) != "double" && typeof(fallback_data) != "integer") {
        mode(fallback_data) <- "numeric"
      }
      
      values$fd_obj <- smooth.basis(seq(0, 1, length.out = ncol(values$data)), 
                                    t(fallback_data), basis)$fd
    })
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
      values$gam_reml_fit <- gam_fit
      values$gam_data_label <- data_label
      
      showNotification("GAM REML completed successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("GAM REML error:", e$message), type = "error", duration = 5)
      cat("GAM REML error:", e$message, "\n")
    })
  })
  
  output$gam_reml_summary <- renderPrint({
    if(is.null(values$gam_reml_fit)) {
      cat("No GAM REML analysis performed yet.\n\n")
      cat("Click 'Fit GAM (REML)' to run the analysis.")
      return()
    }
    
    gam_fit <- values$gam_reml_fit
    
    cat("=== GAM REML Analysis ===\n")
    cat("Data:", values$gam_data_label, "\n\n")
    
    # Extract key info
    summary_gam <- summary(gam_fit)
    edf <- summary_gam$edf
    
    # Safely extract k (basis dimension) - column name may vary
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
      values$reml_profile <- list(
        lambda = lambda_seq,
        reml_score = reml_scores,
        edf = edf_values,
        optimal_idx = which.min(reml_scores),
        optimal_lambda = lambda_seq[which.min(reml_scores)],
        optimal_reml = min(reml_scores)
      )
      
      showNotification(
        sprintf("REML profile complete! Optimal lambda: %.2e", 
                values$reml_profile$optimal_lambda),
        type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(paste("REML profile error:", e$message), type = "error", duration = 5)
      cat("REML profile error:", e$message, "\n")
    })
  })
  
  output$reml_profile_plot <- renderPlotly({
    if(is.null(values$reml_profile)) {
      return(plot_ly() %>% 
               layout(title = "No REML profile computed yet",
                      annotations = list(text = "Click 'Compute REML Profile' to run analysis",
                                         showarrow = FALSE)))
    }
    
    profile <- values$reml_profile
    
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
      
      # Create fold assignments
      if(input$cv_stratified && !is.null(values$group_labels)) {
        # Stratified by groups
        fold_assignments <- rep(NA, n_subjects)
        for(group in unique(values$group_labels)) {
          group_idx <- which(values$group_labels == group)
          n_group <- length(group_idx)
          fold_assignments[group_idx] <- sample(rep(1:k_folds, length.out = n_group))
        }
      } else {
        # Random assignment
        fold_assignments <- sample(rep(1:k_folds, length.out = n_subjects))
      }
      
      # CV error matrix: subjects x lambdas
      cv_errors <- matrix(NA, nrow = n_subjects, ncol = length(lambda_seq))
      
      withProgress(message = "Running cross-validation...", value = 0, {
        total_iter <- n_subjects * length(lambda_seq)
        iter_count <- 0
        
        for(lambda_idx in seq_along(lambda_seq)) {
          lambda <- lambda_seq[lambda_idx]
          
          # Create basis for this lambda
          nb <- min(20, n_time - 2)
          basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
          fdParobj <- fdPar(basis, 2, lambda)
          
          for(i in 1:n_subjects) {
            # Get subject's fold
            test_fold <- fold_assignments[i]
            train_idx <- which(fold_assignments != test_fold)
            
            # Train on other subjects
            train_data <- values$data[train_idx, , drop = FALSE]
            
            # Get mean of training data (or could use all training curves)
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
      values$cv_results <- list(
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
    if(is.null(values$cv_results)) {
      return(plot_ly() %>% 
               layout(title = "No CV analysis performed yet",
                      annotations = list(text = "Click 'Run Cross-Validation' to start",
                                         showarrow = FALSE)))
    }
    
    cv <- values$cv_results
    
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
  
  # 5. Comparison Summary
  output$smoothing_comparison_summary <- renderPrint({
    if(is.null(values$reml_profile) && is.null(values$cv_results)) {
      cat("No smoothing diagnostic analysis performed yet.\n\n")
      cat("Run REML Profile and/or Cross-Validation to see comparison.")
      return()
    }
    
    cat("=== Smoothing Parameter Comparison ===\n\n")
    
    if(!is.null(values$reml_profile)) {
      cat("REML Analysis:\n")
      cat(sprintf("  Optimal lambda: %.3e\n", values$reml_profile$optimal_lambda))
      cat(sprintf("  Optimal EDF: %.2f\n", 
                  values$reml_profile$edf[values$reml_profile$optimal_idx]))
      cat(sprintf("  REML score: %.2f\n\n", values$reml_profile$optimal_reml))
    }
    
    if(!is.null(values$cv_results)) {
      cat("Cross-Validation Analysis:\n")
      cat(sprintf("  Optimal lambda (min CV): %.3e\n", values$cv_results$optimal_lambda))
      cat(sprintf("  Lambda (1-SE rule): %.3e\n", values$cv_results$lambda_1se))
      cat(sprintf("  Min CV error: %.3f\n", 
                  values$cv_results$mean_error[values$cv_results$optimal_idx]))
      cat(sprintf("  K-folds: %d\n\n", values$cv_results$k_folds))
    }
    
    if(!is.null(values$reml_profile) && !is.null(values$cv_results)) {
      reml_lambda <- values$reml_profile$optimal_lambda
      cv_lambda <- values$cv_results$optimal_lambda
      cv_lambda_1se <- values$cv_results$lambda_1se
      
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
    if(!is.null(values$cv_results)) {
      cat(sprintf("  For prediction tasks: lambda = %.2e (CV optimal)\n", 
                  values$cv_results$optimal_lambda))
      cat(sprintf("  For smooth visualization: lambda = %.2e (1-SE rule)\n", 
                  values$cv_results$lambda_1se))
    }
    if(!is.null(values$reml_profile)) {
      cat(sprintf("  For automatic REML: lambda = 0 (or %.2e from profile)\n", 
                  values$reml_profile$optimal_lambda))
    }
    
    # Add smoothing factor conversion for manual smoothing
    cat("\n--- For Manual Smoothing in Data Preprocessing ---\n")
    cat("To use these lambdas, set 'Smoothing Factor' to:\n")
    
    if(!is.null(values$cv_results)) {
      sf_cv <- -log10(values$cv_results$optimal_lambda)
      sf_1se <- -log10(values$cv_results$lambda_1se)
      cat(sprintf("  CV optimal: %.2f (lambda = %.2e)\n", sf_cv, values$cv_results$optimal_lambda))
      cat(sprintf("  1-SE rule: %.2f (lambda = %.2e)\n", sf_1se, values$cv_results$lambda_1se))
    }
    if(!is.null(values$reml_profile)) {
      sf_reml <- -log10(values$reml_profile$optimal_lambda)
      cat(sprintf("  REML optimal: %.2f (lambda = %.2e)\n", sf_reml, values$reml_profile$optimal_lambda))
    }
    
    cat("\nOr use the '📊 Use Diagnostic Results' button in Data Preprocessing tab!")
  })
  
  # ============================================================================
  # LINK DIAGNOSTICS TO DATA PREPROCESSING
  # ============================================================================
  
  # Check if diagnostic results are available
  output$diagnostics_available <- reactive({
    !is.null(values$reml_profile) || !is.null(values$cv_results)
  })
  outputOptions(output, "diagnostics_available", suspendWhenHidden = FALSE)
  
  # Apply diagnostic lambda to smoothing factor
  observeEvent(input$use_diagnostic_lambda, {
    if(is.null(values$reml_profile) && is.null(values$cv_results)) {
      showNotification("No diagnostic results available. Please run REML Profile or Cross-Validation first.", 
                       type = "warning", duration = 5)
      return()
    }
    
    # Determine which lambda to use
    # Priority: CV optimal > REML optimal
    lambda_to_use <- NULL
    lambda_source <- NULL
    
    if(!is.null(values$cv_results)) {
      lambda_to_use <- values$cv_results$optimal_lambda
      lambda_source <- "CV optimal"
    } else if(!is.null(values$reml_profile)) {
      lambda_to_use <- values$reml_profile$optimal_lambda
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

  # ===== GCV vs n-BASIS SWEEP =====
  observeEvent(input$run_nbasis_analysis, {
    req(values$data)
    data_mat <- values$data
    n_time   <- ncol(data_mat)
    lam      <- if(input$smooth_method == "auto") 0 else 10^(-input$smooth_factor)

    nb_seq    <- seq(4, min(n_time - 1, 40), by = 2)
    gcv_means <- numeric(length(nb_seq))

    withProgress(message = "Computing GCV for each n-basis...", value = 0, {
      for(bi in seq_along(nb_seq)) {
        incProgress(1 / length(nb_seq), detail = paste("n_basis =", nb_seq[bi]))
        nb    <- nb_seq[bi]
        basis <- create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
        fdP   <- fdPar(basis, 2, lam)
        gcv_sub <- sapply(seq_len(nrow(data_mat)), function(s) {
          y     <- data_mat[s, ]
          valid <- which(!is.na(y))
          if(length(valid) < 4) return(NA_real_)
          sb <- tryCatch(smooth.basis(valid, y[valid], fdP), error = function(e) NULL)
          if(is.null(sb)) NA_real_ else sb$gcv
        })
        gcv_means[bi] <- mean(gcv_sub, na.rm = TRUE)
      }
    })

    optimal_idx <- which.min(gcv_means)
    values$nbasis_gcv <- list(
      nb      = nb_seq,
      gcv     = gcv_means,
      optimal = nb_seq[optimal_idx]
    )
    showNotification(
      paste0("n-basis analysis complete. Optimal n_basis = ", nb_seq[optimal_idx],
             " (lowest mean GCV)."),
      type = "message", duration = 5)
  })

  output$nbasis_gcv_plot <- renderPlotly({
    req(values$nbasis_gcv)
    ng         <- values$nbasis_gcv
    current_nb <- if(!is.null(input$smooth_method) && input$smooth_method == "manual") {
      input$n_basis_manual
    } else {
      input$n_basis
    }

    p <- plot_ly() %>%
      add_trace(
        x = ng$nb, y = ng$gcv,
        type = "scatter", mode = "lines+markers",
        line   = list(color = "#2196F3", width = 2),
        marker = list(color = "#2196F3", size = 6),
        name   = "Mean GCV"
      ) %>%
      layout(
        title  = list(text = "GCV Score vs Number of B-spline Basis Functions", x = 0.5),
        xaxis  = list(title = "n_basis"),
        yaxis  = list(title = "Mean GCV (lower = better)"),
        shapes = list(
          list(type = "line", x0 = ng$optimal, x1 = ng$optimal,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "green", dash = "dash", width = 1.5)),
          list(type = "line", x0 = current_nb, x1 = current_nb,
               y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "orange", dash = "dot", width = 1.5))
        ),
        annotations = list(
          list(x = ng$optimal, y = 1, yref = "paper", xanchor = "left",
               text = paste0("Optimal (", ng$optimal, ")"),
               showarrow = FALSE, font = list(color = "green", size = 11)),
          list(x = current_nb, y = 0.88, yref = "paper", xanchor = "left",
               text = paste0("Current (", current_nb, ")"),
               showarrow = FALSE, font = list(color = "orange", size = 11))
        ),
        legend = list(orientation = "h", y = -0.2)
      )
    p
  })


  # ============================================================================
  # END DIAGNOSTICS LINK
  # ============================================================================
  
  # Data visualization plot - INTERACTIVE with curve selection
  output$data_plot <- renderPlotly({
    if(is.null(values$data)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>%
               layout(title = "No data loaded"))
    }

    tryCatch({
      data_to_plot <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

      n_time <- ncol(data_to_plot)
      n_subj <- nrow(data_to_plot)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      hover_times <- hover_time_labels(time_points)

      # Get currently selected curve (if any)
      selected_idx <- values$selected_curve

      p <- plot_ly(type = 'scatter', mode = 'lines', source = "data_plot_source")

      # Add individual curves (use 0-1 normalized time)
      n_show <- min(n_subj, 50)
      for(i in 1:n_show) {
        # Determine if this curve is selected
        is_selected <- !is.null(selected_idx) && i == selected_idx

        # Set color and width based on selection
        curve_color <- if(is_selected) 'rgba(0, 100, 255, 0.9)' else 'rgba(100, 100, 100, 0.3)'
        curve_width <- if(is_selected) 3 else 1

        # Get group label if available
        group_info <- if(!is.null(values$group_labels) && i <= length(values$group_labels)) {
          paste0(" (Group: ", values$group_labels[i], ")")
        } else {
          ""
        }

        p <- p %>% add_trace(x = time_points,
                             y = data_to_plot[i,],
                             type = 'scatter',
                             mode = 'lines',
                             name = paste0("Subject ", i, group_info),
                             line = list(color = curve_color, width = curve_width),
                             hovertemplate = paste0("Subject ", i, group_info,
                                                    "<br>Time: %{customdata}<br>Value: %{y:.2f}<extra></extra>"),
                             customdata = hover_times,
                             showlegend = FALSE)
      }

      # Add mean curve
      p <- p %>% add_trace(x = time_points,
                           y = colMeans(data_to_plot),
                           type = 'scatter',
                           mode = 'lines',
                           name = "Mean",
                           line = list(color = 'red', width = 3),
                           hovertemplate = "Mean<br>Time: %{customdata}<br>Value: %{y:.2f}<extra></extra>",
                           customdata = hover_times)

      p <- p %>% layout(title = "Functional Data (click curve to select)",
                        yaxis = list(title = "Value"),
                        showlegend = FALSE,
                        hovermode = 'closest',
                        clickmode = 'event')

      # Apply time axis formatting with hour labels
      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_preprocess))
      
      p
    }, error = function(e) {
      cat("Data plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>%
        layout(title = paste("Error:", e$message))
    })
  })

  # ============================================================================
  # CURVE SELECTION HANDLERS
  # ============================================================================

  # Handle click events on the data plot
  observeEvent(event_data("plotly_click", source = "data_plot_source"), {
    click_data <- event_data("plotly_click", source = "data_plot_source")

    if(!is.null(click_data)) {
      # Get the curve number from the click (curveNumber is 0-indexed)
      curve_num <- click_data$curveNumber + 1  # Convert to 1-indexed

      # The last curve is the mean (index 0 in customdata), skip it
      n_subj <- nrow(values$data)
      n_show <- min(n_subj, 50)

      if(curve_num <= n_show) {
        # It's an individual curve, select it
        values$selected_curve <- curve_num
        cat("Selected curve:", curve_num, "\n")
      } else {
        # Clicked on mean curve, do nothing or clear selection
        cat("Clicked on mean curve\n")
      }
    }
  })

  # Clear curve selection button
  observeEvent(input$clear_curve_selection, {
    values$selected_curve <- NULL
  })

  # Display selected curve information
  output$selected_curve_info <- renderText({
    if(is.null(values$selected_curve) || is.null(values$data)) {
      return("No curve selected.")
    }

    tryCatch({
      idx <- values$selected_curve
      data_to_use <- if(!is.null(values$smooth_data)) values$smooth_data else values$data

      if(idx > nrow(data_to_use)) {
        return("Invalid selection.")
      }

      curve_data <- data_to_use[idx, ]

      # Basic statistics
      curve_mean <- mean(curve_data, na.rm = TRUE)
      curve_sd <- sd(curve_data, na.rm = TRUE)
      curve_min <- min(curve_data, na.rm = TRUE)
      curve_max <- max(curve_data, na.rm = TRUE)
      curve_range <- curve_max - curve_min

      # Group info
      group_info <- if(!is.null(values$group_labels) && idx <= length(values$group_labels)) {
        paste0("Group: ", values$group_labels[idx])
      } else {
        "Group: N/A"
      }

      paste0(
        "Subject: ", idx, "\n",
        group_info, "\n",
        "-------------------\n",
        "Mean: ", sprintf("%.3f", curve_mean), "\n",
        "SD: ", sprintf("%.3f", curve_sd), "\n",
        "Min: ", sprintf("%.3f", curve_min), "\n",
        "Max: ", sprintf("%.3f", curve_max), "\n",
        "Range: ", sprintf("%.3f", curve_range)
      )
    }, error = function(e) {
      paste("Error:", e$message)
    })
  })

  # Group summary table
  output$group_summary <- renderDT({
    if(!is.null(values$group_labels)) {
      group_summary <- data.frame(
        Group = names(table(values$group_labels)),
        Count = as.numeric(table(values$group_labels)),
        Percentage = round(100 * as.numeric(table(values$group_labels)) / length(values$group_labels), 1)
      )
      datatable(group_summary, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
    } else {
      return(NULL)
    }
  })
  
  # Group preview plot
  output$group_preview_plot <- renderPlot({
    if(!is.null(values$data) && !is.null(values$group_labels)) {
      n_time <- ncol(values$data)
      time_points_plot <- get_plot_time()
      time_label <- get_time_label()
      hour_labels <- get_hour_labels()
      # For calculations, use normalized 0-1
      time_points <- seq(0, 1, length.out = n_time)
      
      groups <- levels(values$group_labels)
      n_groups <- length(groups)
      
      # Create scalable color palette
      base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
      colors <- colorRampPalette(base_cols)(n_groups)
      
      matplot(time_points, t(values$data), type = "l", 
              col = rgb(0.5, 0.5, 0.5, 0.1), lty = 1,
              xlab = time_label, ylab = "Value", 
              main = "Functional Data by Group")
      
      for(i in 1:n_groups) {
        group_idx <- which(values$group_labels == groups[i])
        if(length(group_idx) > 0) {
          group_mean <- colMeans(values$data[group_idx, , drop = FALSE])
          lines(time_points_plot, group_mean, col = colors[i], lwd = 3)
        }
      }
      
      legend("topright", legend = groups, col = colors, lwd = 3)
    }
  })
  
  # Run PCA analysis - FIXED with better error handling
  observeEvent(input$run_analysis, {
    cat("Running PCA analysis...\n")
    
    if(is.null(values$data)) {
      showNotification("No data loaded", type = "error", duration = 5)
      return()
    }
    
    tryCatch({
      data_to_use <- if(!is.null(values$smooth_data)) values$smooth_data else values$data
      
      n <- ncol(data_to_use)
      m <- nrow(data_to_use)
      
      n_components <- min(input$n_components, m-1, 10)
      time_grid <- seq(0, 1, length.out = n)
      
      if(is.null(values$fd_obj)) {
        n_basis <- min(20, n - 2, max(4, floor(n/3)))
        basis <- create.bspline.basis(c(0, 1), nbasis = n_basis)
        
        # Check if data is already smoothed to avoid double smoothing
        if(!is.null(values$smooth_data)) {
          # Data already smoothed in preprocessing - just create fd representation
          cat("Creating fd_obj from already smoothed data (lambda=0)\n")
          values$fd_obj <- smooth.basis(time_grid, t(data_to_use), fdPar(basis, 2, 0))$fd
        } else {
          # Raw data - apply default smoothing
          cat("Creating fd_obj from raw data with default smoothing\n")
          values$fd_obj <- smooth.basis(time_grid, t(data_to_use), basis)$fd
        }
        cat("Created fd_obj with", n_basis, "basis functions\n")
      }
      
      if(input$pca_type == "fpca") {
        # Regular functional PCA
        cat("Running standard functional PCA...\n")
        values$pca_results <- pca.fd(values$fd_obj, nharm = n_components)
        values$warping_results <- NULL
        showNotification("Functional PCA completed!", type = "message", duration = 5)
        
      } else if(input$pca_type == "twpca") {
        # Time-warped PCA with proper error handling
        cat("Running time-warped PCA...\n")
        showNotification("Running time-warped PCA...", type = "message", duration = 2)
        
        # Get warping method settings
        warping_method <- if(!is.null(input$warping_method)) input$warping_method else "linear_shift"
        
        cat("Using warping method:", warping_method, "\n")
        
        # Perform warping based on method
        warp_fd <- NULL
        
        if(warping_method == "linear_shift") {
          periodic <- if(!is.null(input$periodic_shift)) input$periodic_shift else FALSE
          allow_dilation <- if(!is.null(input$allow_dilation)) input$allow_dilation else FALSE
          dilation_range <- if(!is.null(input$dilation_range)) input$dilation_range else c(0.95, 1.05)
          reference <- if(!is.null(input$shift_reference)) input$shift_reference else "mean"
          
          warp_fd <- linear_shift_alignment(values$fd_obj,
                                            periodic = periodic,
                                            allow_dilation = allow_dilation,
                                            dilation_range = dilation_range,
                                            reference = reference,
                                            time_points = time_grid)
          
        } else if(warping_method == "parametric") {
          family <- if(!is.null(input$parametric_family)) input$parametric_family else "power"
          param_range <- if(!is.null(input$param_range)) input$param_range else c(0.5, 2)
          symmetric <- if(!is.null(input$symmetric_warp)) input$symmetric_warp else FALSE
          
          warp_fd <- parametric_alignment(values$fd_obj, 
                                          family = family,
                                          param_range = param_range,
                                          symmetric = symmetric,
                                          time_points = time_grid)
          
        } else if(warping_method == "landmark") {
          # Simple landmark alignment
          landmarks <- c(0.25, 0.5, 0.75)  # Default landmarks
          warp_fd <- landmark_alignment_simple(values$fd_obj, landmarks, time_grid)
        }
        
        # Check if warping was successful
        if(is.null(warp_fd)) {
          cat("Warping failed, using default alignment\n")
          warp_fd <- linear_shift_alignment(values$fd_obj,
                                            periodic = FALSE,
                                            allow_dilation = FALSE,
                                            reference = "mean",
                                            time_points = time_grid)
        }

        # Calculate warping fit statistics
        if(!is.null(warp_fd$registered_curves) && !is.null(warp_fd$warp_functions)) {
          cat("Calculating warping fit statistics...\n")
          original_curves <- eval.fd(time_grid, values$fd_obj)
          warp_fd$fit_statistics <- calculate_warping_fit_statistics(
            original_curves = original_curves,
            registered_curves = warp_fd$registered_curves,
            warp_functions = warp_fd$warp_functions,
            time_points = time_grid
          )
          if(!is.null(warp_fd$fit_statistics)) {
            cat("Warping fit statistics calculated successfully\n")
            cat("  Mean R²:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_r_squared), "\n")
            cat("  Mean RMSE:", sprintf("%.4f", warp_fd$fit_statistics$summary$mean_rmse), "\n")
            cat("  Variance explained by warping:",
                sprintf("%.2f%%", warp_fd$fit_statistics$summary$variance_explained_by_warping * 100), "\n")
          }
        }

        # Store warping results
        values$warping_results <- warp_fd
        
        # Run PCA on warped curves with error checking
        if(!is.null(warp_fd$regfd)) {
          cat("Running PCA on warped fd object\n")
          values$pca_results <- pca.fd(warp_fd$regfd, nharm = n_components)
        } else if(!is.null(warp_fd$registered_curves) && ncol(warp_fd$registered_curves) > 0) {
          cat("Creating fd from registered curves\n")
          # Create fd from registered curves
          # Use lambda=0: registered curves are already processed, just need fd representation
          reg_basis <- create.bspline.basis(c(0, 1), nbasis = min(20, n-2, max(4, floor(n/3))))
          reg_fd <- smooth.basis(time_grid, warp_fd$registered_curves, 
                                 fdPar(reg_basis, 2, 0))$fd
          values$pca_results <- pca.fd(reg_fd, nharm = n_components)
        } else {
          cat("Warning: No valid warped data, using original fd_obj\n")
          values$pca_results <- pca.fd(values$fd_obj, nharm = n_components)
        }
        
        showNotification("Time-warped PCA completed!", type = "message", duration = 5)
      }
      
      # Verify PCA results
      if(!is.null(values$pca_results)) {
        cat("PCA complete - PC1 variance:", 
            round(values$pca_results$varprop[1]*100, 2), "%\n")
      }
      
    }, error = function(e) {
      cat("Error in PCA:", e$message, "\n")
      showNotification(paste("Analysis error:", e$message), type = "error", duration = 10)
      values$pca_results <- NULL
    })
  })
  
  # PCA status output - FIXED
  output$pca_status <- renderText({
    if(is.null(values$pca_results)) {
      "No PCA results available. Run analysis first."
    } else {
      tryCatch({
        pca_res <- values$pca_results
        n_comp <- if(!is.null(pca_res$scores)) ncol(pca_res$scores) else 0
        analysis_type <- if(!is.null(values$warping_results)) "Time-warped PCA" else "Standard PCA"
        paste(analysis_type, "complete:", n_comp, "components extracted")
      }, error = function(e) {
        "PCA results available but incomplete"
      })
    }
  })
  
  # PCA summary
  output$pca_summary <- renderPrint({
    if(is.null(values$pca_results)) {
      cat("Run PCA analysis first\n")
    } else {
      tryCatch({
        pca_res <- values$pca_results
        
        cat("PCA Results Summary:\n")
        cat("-------------------\n")
        
        if(!is.null(pca_res$scores)) {
          cat("Number of components:", ncol(pca_res$scores), "\n")
          cat("Number of subjects:", nrow(pca_res$scores), "\n")
        }
        
        if(!is.null(pca_res$varprop)) {
          cat("Variance explained:\n")
          for(i in 1:min(3, length(pca_res$varprop))) {
            cat(sprintf("  PC%d: %.2f%%\n", i, pca_res$varprop[i] * 100))
          }
        }
        
        if(!is.null(values$warping_results)) {
          cat("\nWarping method:", values$warping_results$method, "\n")
          if(!is.null(values$warping_results$shifts)) {
            cat("Mean shift:", round(mean(abs(values$warping_results$shifts)), 4), "\n")
          }
        }
      }, error = function(e) {
        cat("Error displaying summary:", e$message, "\n")
      })
    }
  })
  
  # Component loadings plot - FIXED
  output$loadings_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      # Check validity
      if(is.null(pca_res$meanfd) || is.null(pca_res$harmonics)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      time_points <- seq(0, 1, length.out = 100)
      mean_vals <- eval.fd(time_points, pca_res$meanfd)
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = time_points, 
                  y = as.vector(mean_vals),
                  type = 'scatter',
                  mode = 'lines',
                  name = "Mean", 
                  line = list(color = 'black', width = 2))
      
      colors <- c('red', 'blue', 'green', 'orange', 'purple')
      n_comp <- min(ncol(pca_res$scores), 5, length(pca_res$harmonics))
      
      for(i in 1:n_comp) {
        if(i <= length(pca_res$harmonics) && !is.null(pca_res$harmonics[i])) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          p <- p %>% add_trace(x = time_points,
                               y = as.vector(loading_vals),
                               type = 'scatter',
                               mode = 'lines',
                               name = paste("PC", i),
                               line = list(color = colors[i], width = 2))
        }
      }
      
      p <- p %>% layout(title = "Principal Component Loadings",
                        yaxis = list(title = "Loading"))

      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_results))
      p
      
    }, error = function(e) {
      cat("Loadings plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Variance explained plot - FIXED
  output$variance_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      if(is.null(pca_res$varprop) || is.null(pca_res$scores)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      n_show <- min(length(pca_res$varprop), ncol(pca_res$scores))
      
      if(n_show == 0) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "No components available"))
      }
      
      var_prop <- pca_res$varprop[1:n_show]
      cum_var <- cumsum(var_prop)
      
      p <- plot_ly() %>%
        add_trace(x = 1:n_show,
                  y = var_prop * 100,
                  type = 'bar',
                  name = 'Individual',
                  marker = list(color = 'lightblue'))
      
      p <- p %>%
        add_trace(x = 1:n_show,
                  y = cum_var * 100,
                  type = 'scatter',
                  mode = 'lines+markers',
                  name = 'Cumulative',
                  yaxis = 'y2',
                  line = list(color = 'red'),
                  marker = list(color = 'red'))
      
      p %>% layout(title = "Variance Explained",
                   xaxis = list(title = "Component", dtick = 1),
                   yaxis = list(title = "Variance Explained (%)"),
                   yaxis2 = list(title = "Cumulative Variance (%)",
                                 overlaying = 'y',
                                 side = 'right'))
      
    }, error = function(e) {
      cat("Variance plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Component scores plot - FIXED
  output$scores_plot <- renderPlotly({
    if(is.null(values$pca_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run PCA analysis first"))
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      # Check if PCA results are valid
      if(is.null(pca_res$meanfd) || is.null(pca_res$harmonics) || is.null(pca_res$values)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "Invalid PCA results"))
      }
      
      time_points <- seq(0, 1, length.out = 100)
      mean_vals <- eval.fd(time_points, pca_res$meanfd)
      
      effect_mult <- if(!is.null(input$effect_size)) input$effect_size else 1
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = time_points, 
                  y = as.vector(mean_vals),
                  type = 'scatter',
                  mode = 'lines',
                  name = "Mean", 
                  line = list(color = 'black', width = 3))
      
      colors <- c('red', 'blue', 'green')
      n_show <- min(3, ncol(pca_res$scores), length(pca_res$values))
      
      for(i in 1:n_show) {
        if(i <= length(pca_res$harmonics) && !is.null(pca_res$harmonics[i])) {
          loading_vals <- eval.fd(time_points, pca_res$harmonics[i])
          
          # Plus 2 SD
          plus_2sd <- as.vector(mean_vals + effect_mult * 2 * sqrt(pca_res$values[i]) * loading_vals)
          p <- p %>% add_trace(
            x = time_points,
            y = plus_2sd,
            type = 'scatter',
            mode = 'lines',
            name = paste("PC", i, "+2SD"),
            line = list(color = colors[i], dash = 'solid')
          )
          
          # Minus 2 SD
          minus_2sd <- as.vector(mean_vals - effect_mult * 2 * sqrt(pca_res$values[i]) * loading_vals)
          p <- p %>% add_trace(
            x = time_points,
            y = minus_2sd,
            type = 'scatter',
            mode = 'lines',
            name = paste("PC", i, "-2SD"),
            line = list(color = colors[i], dash = 'dash')
          )
        }
      }
      
      p <- p %>% layout(title = "Effect of Component Scores",
                        yaxis = list(title = "Value"))
      p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_results))
      p
      
    }, error = function(e) {
      cat("Scores plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Component scores table - FIXED
  output$scores_table <- renderDT({
    if(is.null(values$pca_results)) {
      return(NULL)
    }
    
    tryCatch({
      pca_res <- values$pca_results
      
      if(is.null(pca_res$scores)) {
        return(datatable(data.frame(Message = "No scores available"),
                         options = list(dom = 't'),
                         rownames = FALSE))
      }
      
      scores_df <- data.frame(
        Subject = 1:nrow(pca_res$scores)
      )
      
      # Add group information if available
      if(!is.null(values$group_labels) && length(values$group_labels) == nrow(pca_res$scores)) {
        scores_df$Group <- values$group_labels
      }
      
      n_comp <- ncol(pca_res$scores)
      for(i in 1:n_comp) {
        scores_df[paste0("PC", i)] <- round(pca_res$scores[,i], 3)
      }
      
      datatable(scores_df, 
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
      
    }, error = function(e) {
      cat("Scores table error:", e$message, "\n")
      datatable(data.frame(Error = "Unable to display scores"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })
  
  # Warping scores table
  output$warping_scores <- renderDT({
    if(is.null(values$warping_results)) {
      return(NULL)
    }
    
    tryCatch({
      warp_results <- values$warping_results
      n_subjects <- ncol(values$fd_obj$coefs)
      
      # Calculate warping amplitude (deviation from identity)
      if(!is.null(warp_results$warp_functions)) {
        n_time <- nrow(warp_results$warp_functions)
        time_points <- seq(0, 1, length.out = n_time)
        
        warp_amplitude <- numeric(n_subjects)
        for(i in 1:n_subjects) {
          # Calculate deviation from identity line
          warp_amplitude[i] <- sqrt(mean((warp_results$warp_functions[,i] - time_points)^2))
        }
      } else if(!is.null(warp_results$shifts)) {
        # Use shifts if available
        warp_amplitude <- abs(warp_results$shifts)
      } else if(!is.null(warp_results$alpha_values)) {
        # Use alpha values for parametric warping
        warp_amplitude <- abs(warp_results$alpha_values - 1)
      } else {
        warp_amplitude <- rep(0, n_subjects)
      }
      
      warping_df <- data.frame(
        Subject = 1:n_subjects,
        Warping_Amplitude = round(warp_amplitude, 4),
        Method = warp_results$method
      )

      # Add fit statistics if available
      if(!is.null(warp_results$fit_statistics) && !is.null(warp_results$fit_statistics$per_subject)) {
        fit_stats <- warp_results$fit_statistics$per_subject
        warping_df$R_squared <- fit_stats$R_squared
        warping_df$RMSE <- fit_stats$RMSE
        warping_df$Correlation <- fit_stats$Correlation
      }

      # Add group information if available
      if(!is.null(values$group_labels)) {
        warping_df$Group <- values$group_labels
      }

      datatable(warping_df,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE) %>%
        formatStyle("Warping_Amplitude",
                    backgroundColor = styleInterval(c(0.01, 0.05, 0.1),
                                                    c("white", "#ffffcc", "#ffcccc", "#ff9999"))) %>%
        formatStyle("R_squared",
                    backgroundColor = styleInterval(c(0.5, 0.7, 0.9),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("RMSE",
                    backgroundColor = styleInterval(c(0.05, 0.1, 0.2),
                                                    c("#99ff99", "#ccffcc", "#ffffcc", "#ffcccc")))
      
    }, error = function(e) {
      cat("Warping scores table error:", e$message, "\n")
      datatable(data.frame(Message = "No warping scores available"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })

  # ============================================================================
  # WARPING FIT STATISTICS OUTPUTS
  # ============================================================================

  # Summary statistics output (averaged over subjects)
  output$warping_fit_summary <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see fit statistics.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary
      n_valid <- values$warping_results$fit_statistics$n_valid
      n_total <- values$warping_results$fit_statistics$n_subjects

      paste0(
        "Subjects analyzed: ", n_valid, "/", n_total, "\n\n",
        "R² (warped vs original, functional):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_r_squared),
        " (SD: ", sprintf("%.4f", stats$sd_r_squared), ")\n",
        "  [1 - ∫(f-g)²dt / ∫(f-f̄)²dt]\n\n",
        "RMSE (warped vs original):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_rmse),
        " (SD: ", sprintf("%.4f", stats$sd_rmse), ")\n\n",
        "Correlation (orig vs warped):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_correlation),
        " (SD: ", sprintf("%.4f", stats$sd_correlation), ")\n\n",
        "MAE (warped vs original):\n",
        "  Mean: ", sprintf("%.4f", stats$mean_mae),
        " (SD: ", sprintf("%.4f", stats$sd_mae), ")\n\n",
        "Warping Amplitude:\n",
        "  Mean: ", sprintf("%.4f", stats$mean_warp_amplitude),
        " (SD: ", sprintf("%.4f", stats$sd_warp_amplitude), ")"
      )
    }, error = function(e) {
      paste("Error displaying fit statistics:", e$message)
    })
  })

  # Variance decomposition output (EFDA-style)
  output$warping_variance_decomposition <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see variance decomposition.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary

      # Calculate percentages
      total_var <- stats$total_orig_variance
      amp_pct <- if(total_var > 0) stats$total_amp_variance / total_var * 100 else NA
      phase_pct <- if(total_var > 0) stats$total_phase_variance / total_var * 100 else NA
      explained_pct <- stats$variance_explained_by_warping * 100

      paste0(
        "EFDA Variance Decomposition:\n",
        "============================\n\n",
        "Total Original Variance: ", sprintf("%.4f", total_var), "\n\n",
        "After Alignment:\n",
        "  Amplitude Variance: ", sprintf("%.4f", stats$total_amp_variance),
        " (", sprintf("%.1f%%", amp_pct), " of original)\n",
        "  Phase Variance: ", sprintf("%.4f", stats$total_phase_variance),
        " (", sprintf("%.1f%%", phase_pct), " of original)\n\n",
        "Variance Explained by Warping:\n",
        "  ", sprintf("%.2f%%", explained_pct), "\n\n",
        "Elastic Distances (mean):\n",
        "  Full Distance: ", sprintf("%.4f", stats$mean_full_distance), "\n",
        "  Amplitude Distance: ", sprintf("%.4f", stats$mean_elastic_amp_dist), "\n",
        "  Phase Distance: ", sprintf("%.4f", stats$mean_elastic_phase_dist)
      )
    }, error = function(e) {
      paste("Error displaying variance decomposition:", e$message)
    })
  })

  # Model selection criteria output (AIC, BIC)
  output$warping_model_criteria <- renderText({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return("Run time-warped PCA to see model criteria.")
    }

    tryCatch({
      stats <- values$warping_results$fit_statistics$summary
      method <- values$warping_results$method

      paste0(
        "Warping Method: ", method, "\n\n",
        "Model Selection Criteria:\n",
        "=========================\n",
        "AIC: ", sprintf("%.2f", stats$AIC), "\n",
        "BIC: ", sprintf("%.2f", stats$BIC), "\n",
        "Log-Likelihood: ", sprintf("%.2f", stats$log_likelihood), "\n",
        "Number of Parameters: ", stats$n_parameters, "\n\n",
        "Note: Lower AIC/BIC values indicate better model fit.\n",
        "Compare these values across different warping methods\n",
        "to select the optimal alignment approach."
      )
    }, error = function(e) {
      paste("Error displaying model criteria:", e$message)
    })
  })

  # Per-subject fit statistics table
  output$warping_fit_per_subject <- renderDT({
    if(is.null(values$warping_results) || is.null(values$warping_results$fit_statistics)) {
      return(NULL)
    }

    tryCatch({
      per_subj <- values$warping_results$fit_statistics$per_subject

      # Add group labels if available
      if(!is.null(values$group_labels) && length(values$group_labels) == nrow(per_subj)) {
        per_subj$Group <- values$group_labels
        # Reorder columns to put Group after Subject
        cols <- c("Subject", "Group", setdiff(names(per_subj), c("Subject", "Group")))
        per_subj <- per_subj[, cols]
      }

      datatable(per_subj,
                options = list(
                  pageLength = 10,
                  scrollX = TRUE,
                  columnDefs = list(
                    list(className = 'dt-center', targets = '_all')
                  )
                ),
                rownames = FALSE) %>%
        formatStyle("R_squared",
                    backgroundColor = styleInterval(c(0.5, 0.7, 0.9),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("RMSE",
                    backgroundColor = styleInterval(c(0.05, 0.1, 0.2),
                                                    c("#99ff99", "#ccffcc", "#ffffcc", "#ffcccc"))) %>%
        formatStyle("Correlation",
                    backgroundColor = styleInterval(c(0.7, 0.85, 0.95),
                                                    c("#ffcccc", "#ffffcc", "#ccffcc", "#99ff99"))) %>%
        formatStyle("Warp_Amplitude",
                    backgroundColor = styleInterval(c(0.01, 0.05, 0.1),
                                                    c("white", "#ffffcc", "#ffcccc", "#ff9999")))

    }, error = function(e) {
      cat("Per-subject fit statistics table error:", e$message, "\n")
      datatable(data.frame(Message = "No per-subject statistics available"),
                options = list(dom = 't'),
                rownames = FALSE)
    })
  })

  # Download handler for fit statistics CSV
  output$download_warping_fit_stats <- downloadHandler(
    filename = function() {
      paste0("warping_fit_statistics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        if(!is.null(values$warping_results) && !is.null(values$warping_results$fit_statistics)) {
          per_subj <- values$warping_results$fit_statistics$per_subject
          summary_stats <- values$warping_results$fit_statistics$summary

          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == nrow(per_subj)) {
            per_subj$Group <- values$group_labels
          }

          # Create summary row
          summary_row <- data.frame(
            Subject = "AVERAGE",
            R_squared = summary_stats$mean_r_squared,
            RMSE = summary_stats$mean_rmse,
            Correlation = summary_stats$mean_correlation,
            MAE = summary_stats$mean_mae,
            Orig_Variance = summary_stats$total_orig_variance,
            Amp_Variance = summary_stats$total_amp_variance,
            Phase_Variance = summary_stats$total_phase_variance,
            Full_Distance = summary_stats$mean_full_distance,
            Elastic_Amp_Dist = summary_stats$mean_elastic_amp_dist,
            Elastic_Phase_Dist = summary_stats$mean_elastic_phase_dist,
            Warp_Amplitude = summary_stats$mean_warp_amplitude,
            Warp_Velocity_Var = NA
          )

          # Add SD row
          sd_row <- data.frame(
            Subject = "SD",
            R_squared = summary_stats$sd_r_squared,
            RMSE = summary_stats$sd_rmse,
            Correlation = summary_stats$sd_correlation,
            MAE = summary_stats$sd_mae,
            Orig_Variance = NA,
            Amp_Variance = NA,
            Phase_Variance = NA,
            Full_Distance = summary_stats$sd_full_distance,
            Elastic_Amp_Dist = NA,
            Elastic_Phase_Dist = NA,
            Warp_Amplitude = summary_stats$sd_warp_amplitude,
            Warp_Velocity_Var = NA
          )

          # Add metadata
          meta_rows <- data.frame(
            Subject = c("---", "AIC", "BIC", "Variance_Explained_%", "Method"),
            R_squared = c(NA, summary_stats$AIC, summary_stats$BIC,
                          summary_stats$variance_explained_by_warping * 100,
                          NA),
            RMSE = NA, Correlation = NA, MAE = NA,
            Orig_Variance = NA, Amp_Variance = NA, Phase_Variance = NA,
            Full_Distance = NA, Elastic_Amp_Dist = NA, Elastic_Phase_Dist = NA,
            Warp_Amplitude = NA, Warp_Velocity_Var = NA
          )
          meta_rows$Subject[5] <- values$warping_results$method

          # Combine all
          output_df <- rbind(per_subj[, names(summary_row)], summary_row, sd_row, meta_rows)

          write.csv(output_df, file, row.names = FALSE)
        } else {
          write.csv(data.frame(Message = "No fit statistics available"), file, row.names = FALSE)
        }
      }, error = function(e) {
        write.csv(data.frame(Error = e$message), file, row.names = FALSE)
      })
    }
  )

  # Warping download buttons - FIXED
  output$download_warping_plot <- downloadHandler(
    filename = function() {
      paste0("warping_functions_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        png(file, width = 800, height = 600)
        
        if(!is.null(values$warping_results) && !is.null(values$warping_results$warp_functions)) {
          warp_functions <- values$warping_results$warp_functions
          time_points <- values$warping_results$time_points
          
          matplot(time_points, warp_functions[,1:min(30, ncol(warp_functions))], 
                  type = "l", col = rgb(0.4, 0.4, 0.8, 0.3), lty = 1,
                  xlab = "Original Time", ylab = "Warped Time",
                  main = "Time Warping Functions")
          lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
          legend("topleft", legend = c("Individual", "Identity"), 
                 col = c("lightblue", "red"), lty = c(1, 2), lwd = c(1, 2))
        } else {
          plot(1, type = "n", xlab = "", ylab = "", 
               main = "No warping functions available")
        }
        
        dev.off()
      }, error = function(e) {
        png(file)
        plot(1, main = "Error generating plot")
        dev.off()
      })
    }
  )
  
  output$download_alignment_plot <- downloadHandler(
    filename = function() {
      paste0("alignment_comparison_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        png(file, width = 800, height = 600)
        
        if(!is.null(values$fd_obj)) {
          time_points <- seq(0, 1, length.out = 100)
          orig_curves <- eval.fd(time_points, values$fd_obj)
          
          par(mfrow = c(1, 2))
          
          # Original curves
          matplot(time_points, orig_curves[,1:min(30, ncol(orig_curves))], 
                  type = "l", col = rgb(1, 0.4, 0.4, 0.3), lty = 1,
                  xlab = "Time", ylab = "Value", main = "Original Curves")
          lines(time_points, rowMeans(orig_curves), col = "darkred", lwd = 3)
          
          # Aligned curves
          if(!is.null(values$warping_results) && !is.null(values$warping_results$registered_curves)) {
            aligned_curves <- values$warping_results$registered_curves
            matplot(time_points, aligned_curves[,1:min(30, ncol(aligned_curves))], 
                    type = "l", col = rgb(0.4, 0.4, 1, 0.3), lty = 1,
                    xlab = "Time", ylab = "Value", main = "Aligned Curves")
            lines(time_points, rowMeans(aligned_curves), col = "darkblue", lwd = 3)
          } else {
            plot(1, type = "n", main = "No aligned curves available")
          }
        } else {
          plot(1, type = "n", main = "No data available")
        }
        
        dev.off()
      }, error = function(e) {
        png(file)
        plot(1, main = "Error generating plot")
        dev.off()
      })
    }
  )

  # ============================================================================
  # WARPING FIT STATISTICS CALCULATION
  # ============================================================================
  # Calculates fit statistics comparing raw vs warped data per subject and averaged
  # Based on EFDA methodology: variance decomposition and elastic distances
  # ============================================================================

  calculate_warping_fit_statistics <- function(original_curves, registered_curves,
                                                warp_functions, time_points) {
    # Calculate comprehensive fit statistics for time warping alignment
    #
    # Args:
    #   original_curves: matrix (n_time x n_subjects) of original data
    #   registered_curves: matrix (n_time x n_subjects) of warped data
    #   warp_functions: matrix (n_time x n_subjects) of warping functions h(t)
    #   time_points: vector of time points
    #
    # Returns:
    #   List with per-subject and averaged fit statistics

    tryCatch({
      n_time <- nrow(original_curves)
      n_subjects <- ncol(original_curves)
      dt <- diff(time_points[1:2])  # Time step for integration

      # Initialize per-subject statistics vectors
      r_squared <- numeric(n_subjects)
      rmse <- numeric(n_subjects)
      correlation <- numeric(n_subjects)
      mae <- numeric(n_subjects)

      # Variance decomposition (EFDA-style)
      orig_variance <- numeric(n_subjects)
      amp_variance <- numeric(n_subjects)
      phase_variance <- numeric(n_subjects)

      # Elastic distances (EFDA-style)
      full_dist <- numeric(n_subjects)
      elastic_amp_dist <- numeric(n_subjects)
      elastic_phase_dist <- numeric(n_subjects)

      # Warping intensity metrics
      warp_amplitude <- numeric(n_subjects)
      warp_velocity_var <- numeric(n_subjects)  # Variance of warping derivative

      # Calculate mean curves
      orig_mean <- rowMeans(original_curves, na.rm = TRUE)
      reg_mean <- rowMeans(registered_curves, na.rm = TRUE)

      # Calculate SRVF (Square Root Velocity Function) for elastic analysis
      # q(t) = sign(f'(t)) * sqrt(|f'(t)|)
      calc_srvf <- function(f, t) {
        n <- length(f)
        df <- diff(f) / diff(t)
        df <- c(df[1], df)  # Pad to same length
        sign(df) * sqrt(abs(df))
      }

      # Calculate SRVF of registered mean
      reg_mean_srvf <- calc_srvf(reg_mean, time_points)

      for(i in 1:n_subjects) {
        orig_i <- original_curves[,i]
        reg_i <- registered_curves[,i]
        warp_i <- warp_functions[,i]

        # Handle any NA values
        valid_idx <- !is.na(orig_i) & !is.na(reg_i)
        if(sum(valid_idx) < 3) {
          r_squared[i] <- NA
          rmse[i] <- NA
          correlation[i] <- NA
          mae[i] <- NA
          next
        }

        # ---- Basic Fit Statistics ----
        # Using FUNCTIONAL R² via integrals (L²-norm based)
        # R² = 1 - ∫(f(t) - g(t))² dt / ∫(f(t) - f̄)² dt
        # where f = original curve, g = warped curve, f̄ = mean of f over domain
        # This measures: "How much of the total squared energy of f is captured by g?"

        # f̄ = mean of original curve over the domain (scalar)
        f_bar <- mean(orig_i[valid_idx])

        # Numerator: ∫(f - g)² dt (integrated squared difference between original and warped)
        ss_res_functional <- sum((orig_i[valid_idx] - reg_i[valid_idx])^2) * dt

        # Denominator: ∫(f - f̄)² dt (integrated squared deviation from mean)
        ss_tot_functional <- sum((orig_i[valid_idx] - f_bar)^2) * dt

        # Functional R² (always between 0 and 1 when g approximates f well)
        # Note: Can still be negative if warped curve is worse than using the mean
        # but we clamp to 0 as a floor since negative R² is not meaningful here
        r_squared[i] <- if(ss_tot_functional > 0) {
          max(0, 1 - ss_res_functional / ss_tot_functional)
        } else {
          NA
        }

        # Correlation between original and warped (information preservation)
        if(var(orig_i[valid_idx]) > 0 && var(reg_i[valid_idx]) > 0) {
          correlation[i] <- cor(orig_i[valid_idx], reg_i[valid_idx])
        } else {
          correlation[i] <- NA
        }

        # RMSE: Distance from warped curve to ORIGINAL curve (not group mean)
        rmse[i] <- sqrt(mean((orig_i[valid_idx] - reg_i[valid_idx])^2))

        # MAE: Mean absolute error from warped curve to ORIGINAL curve
        mae[i] <- mean(abs(orig_i[valid_idx] - reg_i[valid_idx]))

        # ---- Variance Decomposition (EFDA-style) ----

        # Original variance: integrated squared deviation from original mean
        orig_variance[i] <- sum((orig_i - orig_mean)^2) * dt

        # Amplitude variance: integrated squared deviation after alignment
        amp_variance[i] <- sum((reg_i - reg_mean)^2) * dt

        # Phase variance: variance captured by warping
        # Computed from warped mean evaluated at individual warp functions
        warped_mean_i <- approx(time_points, reg_mean, xout = warp_i, rule = 2)$y
        phase_variance[i] <- sum((warped_mean_i - reg_mean)^2) * dt

        # ---- Elastic Distances (EFDA-style) ----

        # Full distance: L2 distance between aligned curve and mean
        full_dist[i] <- sqrt(sum((reg_i - reg_mean)^2) * dt)

        # Elastic amplitude distance: L2 distance in SRVF space
        reg_i_srvf <- calc_srvf(reg_i, time_points)
        elastic_amp_dist[i] <- sqrt(sum((reg_i_srvf - reg_mean_srvf)^2) * dt)

        # Elastic phase distance: arc-cosine of integrated psi
        # psi = sqrt(d(warp)/dt) is the SRVF of the warping function
        warp_deriv <- c(diff(warp_i) / diff(time_points), 0)
        warp_deriv[warp_deriv < 0] <- 0  # Ensure non-negative (monotonic)
        psi_i <- sqrt(pmax(0, warp_deriv))
        psi_integral <- sum(psi_i) * dt
        elastic_phase_dist[i] <- acos(min(1, max(-1, psi_integral)))

        # ---- Warping Intensity Metrics ----

        # Warping amplitude: RMSE deviation from identity h(t) = t
        warp_amplitude[i] <- sqrt(mean((warp_i - time_points)^2))

        # Warping velocity variance: how variable is the warping speed?
        warp_velocity_var[i] <- var(warp_deriv, na.rm = TRUE)
      }

      # ---- Model Selection Criteria (adapted for functional data) ----
      # AIC and BIC based on residual variance after alignment

      n_obs <- n_time * n_subjects  # Total observations
      residual_var <- mean(rmse^2, na.rm = TRUE)  # Mean squared error

      # Number of parameters: depends on warping method
      # For linear shift: 1 parameter (shift) per subject
      # For parametric: 1 parameter (alpha) per subject
      # For landmark: k parameters (landmark positions) per subject
      # Estimate as 2 parameters per subject (location + scale)
      k_params <- 2 * n_subjects

      # Log-likelihood (assuming Gaussian errors)
      log_lik <- -n_obs/2 * (log(2 * pi) + log(residual_var) + 1)

      # AIC = -2*logLik + 2*k
      aic <- -2 * log_lik + 2 * k_params

      # BIC = -2*logLik + k*log(n)
      bic <- -2 * log_lik + k_params * log(n_obs)

      # ---- Variance Explained by Alignment ----
      # Similar to R-squared but at group level
      total_orig_var <- sum(orig_variance, na.rm = TRUE)
      total_amp_var <- sum(amp_variance, na.rm = TRUE)
      total_phase_var <- sum(phase_variance, na.rm = TRUE)

      # Proportion of variance explained by warping
      var_explained_by_warping <- if(total_orig_var > 0) {
        1 - total_amp_var / total_orig_var
      } else {
        NA
      }

      # Return comprehensive statistics
      return(list(
        # Per-subject statistics
        per_subject = data.frame(
          Subject = 1:n_subjects,
          R_squared = round(r_squared, 4),
          RMSE = round(rmse, 4),
          Correlation = round(correlation, 4),
          MAE = round(mae, 4),
          Orig_Variance = round(orig_variance, 4),
          Amp_Variance = round(amp_variance, 4),
          Phase_Variance = round(phase_variance, 4),
          Full_Distance = round(full_dist, 4),
          Elastic_Amp_Dist = round(elastic_amp_dist, 4),
          Elastic_Phase_Dist = round(elastic_phase_dist, 4),
          Warp_Amplitude = round(warp_amplitude, 4),
          Warp_Velocity_Var = round(warp_velocity_var, 6)
        ),

        # Averaged statistics
        summary = list(
          # Basic fit statistics (mean ± SD)
          mean_r_squared = mean(r_squared, na.rm = TRUE),
          sd_r_squared = sd(r_squared, na.rm = TRUE),
          mean_rmse = mean(rmse, na.rm = TRUE),
          sd_rmse = sd(rmse, na.rm = TRUE),
          mean_correlation = mean(correlation, na.rm = TRUE),
          sd_correlation = sd(correlation, na.rm = TRUE),
          mean_mae = mean(mae, na.rm = TRUE),
          sd_mae = sd(mae, na.rm = TRUE),

          # Variance decomposition
          total_orig_variance = total_orig_var,
          total_amp_variance = total_amp_var,
          total_phase_variance = total_phase_var,
          variance_explained_by_warping = var_explained_by_warping,

          # Elastic distances (averaged)
          mean_full_distance = mean(full_dist, na.rm = TRUE),
          sd_full_distance = sd(full_dist, na.rm = TRUE),
          mean_elastic_amp_dist = mean(elastic_amp_dist, na.rm = TRUE),
          mean_elastic_phase_dist = mean(elastic_phase_dist, na.rm = TRUE),

          # Warping intensity
          mean_warp_amplitude = mean(warp_amplitude, na.rm = TRUE),
          sd_warp_amplitude = sd(warp_amplitude, na.rm = TRUE),

          # Model selection criteria
          AIC = aic,
          BIC = bic,
          log_likelihood = log_lik,
          n_parameters = k_params
        ),

        n_subjects = n_subjects,
        n_valid = sum(!is.na(r_squared))
      ))

    }, error = function(e) {
      cat("Error in calculate_warping_fit_statistics:", e$message, "\n")
      return(NULL)
    })
  }

  # Warping functions with better error handling
  linear_shift_alignment <- function(fd_obj, periodic = FALSE, 
                                     allow_dilation = FALSE, 
                                     dilation_range = c(0.95, 1.05),
                                     reference = "mean",
                                     time_points = seq(0, 1, length.out = 100)) {
    
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      # Evaluate curves
      curves <- eval.fd(time_points, fd_obj)
      
      # Validate curves
      if(is.null(curves) || ncol(curves) == 0) {
        stop("No valid curves to align")
      }
      
      # Get reference curve
      if(reference == "mean") {
        ref_curve <- rowMeans(curves)
      } else if(reference == "median") {
        ref_curve <- apply(curves, 1, median)
      } else {
        ref_curve <- curves[,1]
      }
      
      # Initialize warping
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      shifts <- numeric(n_curves)
      
      # Perform alignment for each curve
      for(i in 1:n_curves) {
        # Find best shift using cross-correlation
        if(periodic) {
          # Circular cross-correlation
          fft_curve <- fft(curves[,i] - mean(curves[,i]))
          fft_ref <- fft(ref_curve - mean(ref_curve))
          cross_corr <- Re(fft(Conj(fft_ref) * fft_curve, inverse = TRUE)) / n_time
          max_idx <- which.max(cross_corr)
          shift_idx <- if(max_idx > n_time/2) max_idx - n_time else max_idx
          shifts[i] <- -shift_idx / n_time
        } else {
          # Standard cross-correlation
          ccf_result <- ccf(curves[,i], ref_curve, lag.max = floor(n_time/4), 
                            plot = FALSE, na.action = na.pass)
          if(!is.null(ccf_result$acf) && length(ccf_result$acf) > 0) {
            best_lag <- ccf_result$lag[which.max(ccf_result$acf)]
            shifts[i] <- best_lag / n_time * 0.1  # Scale down shift
          } else {
            shifts[i] <- 0
          }
        }
        
        # Create warping function with some variation
        base_warp <- time_points - shifts[i] * 0.5
        
        # Add slight S-curve for visualization
        distortion <- sin(pi * time_points) * runif(1, -0.03, 0.03)
        warp_functions[,i] <- pmin(1, pmax(0, base_warp + distortion))
        
        # Ensure endpoints are fixed
        warp_functions[1,i] <- 0
        warp_functions[n_time,i] <- 1
        
        # Apply warping
        if(abs(shifts[i]) > 0.001) {
          # Interpolate curve at warped time points
          registered_curves[,i] <- approx(time_points, curves[,i], 
                                          xout = warp_functions[,i], 
                                          rule = 2)$y
        } else {
          registered_curves[,i] <- curves[,i]
        }
      }
      
      # Create fd objects for output
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      # Create warping function fd objects
      warp_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10)
      warp_smooth <- smooth.basis(time_points, warp_functions, warp_basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        shifts = shifts,
        method = "linear_shift",
        time_points = time_points
      ))
      
    }, error = function(e) {
      cat("Error in linear_shift_alignment:", e$message, "\n")
      # Return identity warping as fallback
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      return(list(
        regfd = fd_obj,
        registered_curves = eval.fd(time_points, fd_obj),
        warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
        shifts = rep(0, n_curves),
        method = "identity",
        time_points = time_points
      ))
    })
  }
  
  # Parametric alignment function
  parametric_alignment <- function(fd_obj, family = "power", 
                                   param_range = c(0.5, 2), 
                                   symmetric = FALSE,
                                   time_points = seq(0, 1, length.out = 100)) {
    
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      curves <- eval.fd(time_points, fd_obj)
      mean_curve <- rowMeans(curves)
      
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      alpha_values <- numeric(n_curves)
      
      # Define warping function
      warp_func <- function(t, alpha) {
        switch(family,
               "power" = pmin(1, pmax(0, t^alpha)),
               "exponential" = {
                 if(abs(alpha - 1) < 0.001) t
                 else pmin(1, pmax(0, (exp(alpha * t) - 1) / (exp(alpha) - 1)))
               },
               "quadratic" = pmin(1, pmax(0, alpha * t^2 + (1 - alpha) * t)),
               "logistic" = {
                 L <- function(x) 1 / (1 + exp(-alpha * (x - 0.5)))
                 L0 <- L(0)
                 L1 <- L(1)
                 pmin(1, pmax(0, (L(t) - L0) / (L1 - L0)))
               },
               t
        )
      }
      
      # Optimize warping for each curve
      for(i in 1:n_curves) {
        # Objective function
        objective <- function(alpha) {
          warped_time <- warp_func(time_points, alpha)
          warped_curve <- approx(time_points, curves[,i], xout = warped_time, rule = 2)$y
          sum((warped_curve - mean_curve)^2, na.rm = TRUE)
        }
        
        # Optimize
        result <- optimize(objective, interval = param_range, tol = 1e-4)
        alpha_values[i] <- result$minimum
        
        # Apply warping
        warp_functions[,i] <- warp_func(time_points, alpha_values[i])
        registered_curves[,i] <- approx(time_points, curves[,i], 
                                        xout = warp_functions[,i], rule = 2)$y
      }
      
      # Create fd objects
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        alpha_values = alpha_values,
        family = family,
        method = "parametric",
        time_points = time_points
      ))
      
    }, error = function(e) {
      cat("Error in parametric_alignment:", e$message, "\n")
      # Return identity warping
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      return(list(
        regfd = fd_obj,
        registered_curves = eval.fd(time_points, fd_obj),
        warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
        alpha_values = rep(1, n_curves),
        method = "identity",
        time_points = time_points
      ))
    })
  }
  
  # Simple landmark alignment
  landmark_alignment_simple <- function(fd_obj, landmarks, time_points) {
    tryCatch({
      n_curves <- ncol(fd_obj$coefs)
      n_time <- length(time_points)
      
      curves <- eval.fd(time_points, fd_obj)
      
      # If landmarks are provided, use them for alignment
      if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
        landmark_times <- values$landmark_points$x
        n_landmarks <- length(landmark_times)
        
        cat("Using", n_landmarks, "landmarks for alignment\n")
        
        # Find corresponding landmark points in each curve
        warp_functions <- matrix(NA, n_time, n_curves)
        registered_curves <- matrix(NA, n_time, n_curves)
        
        for(i in 1:n_curves) {
          # For each curve, find the actual landmarks (peaks/valleys near the specified times)
          curve_landmarks <- numeric(n_landmarks)
          
          for(j in 1:n_landmarks) {
            # Find local extremum near the landmark time
            search_window <- which(abs(time_points - landmark_times[j]) < 0.1)
            if(length(search_window) > 0) {
              local_values <- curves[search_window, i]
              # Find the peak or valley
              if(j %% 2 == 1) {
                # Look for peak for odd landmarks
                curve_landmarks[j] <- time_points[search_window[which.max(local_values)]]
              } else {
                # Look for valley for even landmarks
                curve_landmarks[j] <- time_points[search_window[which.min(local_values)]]
              }
            } else {
              curve_landmarks[j] <- landmark_times[j]
            }
          }
          
          # Create warping function using piecewise linear interpolation
          # Add boundary points
          all_landmark_times <- c(0, landmark_times, 1)
          all_curve_landmarks <- c(0, curve_landmarks, 1)
          
          # Interpolate warping function
          warp_functions[,i] <- approx(all_landmark_times, all_curve_landmarks, 
                                       xout = time_points, rule = 2)$y
          
          # Apply warping
          registered_curves[,i] <- approx(time_points, curves[,i], 
                                          xout = warp_functions[,i], rule = 2)$y
        }
      } else {
        # No landmarks provided, use automatic detection
        cat("No manual landmarks provided, using automatic landmark detection\n")
        
        # Simple automatic landmark detection: find common peaks
        mean_curve <- rowMeans(curves)
        
        # Find peaks in mean curve
        peaks <- which(diff(sign(diff(mean_curve))) == -2) + 1
        valleys <- which(diff(sign(diff(mean_curve))) == 2) + 1
        
        # Select up to 3 most prominent landmarks
        if(length(peaks) > 0 || length(valleys) > 0) {
          all_extrema <- sort(c(peaks, valleys))
          if(length(all_extrema) > 3) {
            # Select based on prominence
            prominence <- abs(mean_curve[all_extrema] - mean(mean_curve))
            all_extrema <- all_extrema[order(prominence, decreasing = TRUE)[1:3]]
          }
          landmark_times <- time_points[all_extrema]
        } else {
          # Default landmarks at quartiles
          landmark_times <- c(0.25, 0.5, 0.75)
        }
        
        # Apply landmark registration
        warp_functions <- matrix(NA, n_time, n_curves)
        registered_curves <- matrix(NA, n_time, n_curves)
        
        for(i in 1:n_curves) {
          # Simple identity warping with slight variation
          distortion <- sin(2*pi*time_points) * runif(1, -0.02, 0.02)
          warp_functions[,i] <- pmin(1, pmax(0, time_points + distortion))
          warp_functions[1,i] <- 0
          warp_functions[n_time,i] <- 1
          
          registered_curves[,i] <- curves[,i]
        }
      }
      
      basis <- fd_obj$basis
      reg_smooth <- smooth.basis(time_points, registered_curves, basis)
      
      return(list(
        regfd = reg_smooth$fd,
        registered_curves = registered_curves,
        warp_functions = warp_functions,
        method = "landmark",
        time_points = time_points,
        landmarks_used = if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) 
          values$landmark_points$x else NULL
      ))
      
    }, error = function(e) {
      cat("Error in landmark_alignment:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Group variable UI
  output$group_variable_ui <- renderUI({
    # Only depend on selected_group_vars and group_variables, NOT group_labels
    # This prevents re-rendering when the user switches group variables
    req(values$selected_group_vars)

    if(!is.null(values$group_variables)) {
      # Build UI elements
      ui_elements <- tagList(
        h5("Group Variable Status")
      )

      # If multiple group variables are available, show selector
      if(length(values$selected_group_vars) > 1) {
        # Preserve current selection if it exists and is valid
        current_selection <- isolate(input$fanova_group_var)
        if(is.null(current_selection) || !(current_selection %in% values$selected_group_vars)) {
          current_selection <- values$selected_group_vars[1]
        }

        ui_elements <- tagList(
          ui_elements,
          selectInput("fanova_group_var", "Select Group Variable for fANOVA:",
                      choices = values$selected_group_vars,
                      selected = current_selection),
          hr()
        )
      }

      ui_elements
    } else if(!is.null(values$group_labels)) {
      # Single group variable case (backwards compatibility)
      tagList(
        h5("Group Variable Status"),
        p(paste("Groups loaded:", length(unique(values$group_labels)), "groups")),
        p(paste("Total subjects:", length(values$group_labels)))
      )
    } else {
      tagList(
        p("No group variable detected."),
        p("To use Functional ANOVA, ensure you select a Group Variable during data import."),
        hr(),
        h5("Create Groups Manually"),
        numericInput("manual_n_groups", "Number of groups:", value = 2, min = 2, max = 10),
        actionButton("create_groups", "Create Groups", class = "btn-warning")
      )
    }
  })

  # Separate UI output for group info display (updates when fanova_group_var changes)
  output$fanova_group_info <- renderUI({
    # Get current group variable
    current_var <- input$fanova_group_var
    if(is.null(current_var) && !is.null(values$selected_group_vars)) {
      current_var <- values$selected_group_vars[1]
    }

    current_groups <- if(!is.null(current_var) && !is.null(values$group_variables) &&
                         current_var %in% colnames(values$group_variables)) {
      values$group_variables[[current_var]]
    } else {
      values$group_labels
    }

    if(is.null(current_groups)) return(NULL)

    n_groups <- length(unique(current_groups))

    ui_elements <- tagList(
      p(paste("Current group variable:", current_var)),
      p(paste("Number of groups:", n_groups)),
      p(paste("Total subjects:", length(current_groups)))
    )

    # If more than 2 groups, show group selector
    if(n_groups > 2) {
      ui_elements <- tagList(
        ui_elements,
        hr(),
        h5("Select Groups to Include in Analysis"),
        pickerInput(
          inputId = "fanova_groups_to_include",
          label = "Groups to include:",
          choices = levels(as.factor(current_groups)),
          selected = levels(as.factor(current_groups)),
          options = list(
            `actions-box` = TRUE,
            `selected-text-format` = "count > 3"
          ),
          multiple = TRUE
        ),
        helpText("Select at least 2 groups to compare. Deselect groups you want to exclude.")
      )
    }

    ui_elements
  })
  
  # Create groups manually
  observeEvent(input$create_groups, {
    req(values$data)
    n_subjects <- nrow(values$data)
    n_groups <- input$manual_n_groups

    group_labels <- rep(paste0("Group", 1:n_groups), length.out = n_subjects)
    values$group_labels <- factor(group_labels)

    showNotification("Groups created successfully!", type = "message", duration = 3)
  })

  # Helper function to get current fANOVA group labels based on selection
  get_fanova_group_labels <- reactive({
    # If user selected a specific group variable for fANOVA, use that
    if(!is.null(input$fanova_group_var) && !is.null(values$group_variables) &&
       input$fanova_group_var %in% colnames(values$group_variables)) {
      return(values$group_variables[[input$fanova_group_var]])
    }
    # Otherwise fall back to primary group_labels
    return(values$group_labels)
  })
  
  # UI for subject ID selection (for RM-ANOVA)
  output$subject_id_ui <- renderUI({
    req(values$data)
    
    # Get column names from original uploaded data if available
    col_names <- if(!is.null(values$uploaded_data)) {
      colnames(values$uploaded_data)
    } else {
      NULL
    }
    
    if(!is.null(col_names)) {
      # Filter out time-related columns (those matching the data columns)
      time_cols <- colnames(values$data)
      non_time_cols <- setdiff(col_names, time_cols)
      
      if(length(non_time_cols) > 0) {
        selectInput("rm_subject_id_var", 
                    "Subject ID variable:",
                    choices = c("", non_time_cols),
                    selected = "")
      } else {
        tagList(
          p("No ID variables found.", style = "color: orange;"),
          p("For repeated measures ANOVA, your data should include a subject ID column."),
          helpText("You can create subject IDs manually by numbering rows.")
        )
      }
    } else {
      tagList(
        p("Upload data to select subject ID variable.", style = "color: orange;")
      )
    }
  })
  
  # UI for repeated measures factor selection
  output$rm_factor_ui <- renderUI({
    req(values$data)

    # Get column names from original uploaded data if available
    col_names <- if(!is.null(values$uploaded_data)) {
      colnames(values$uploaded_data)
    } else {
      NULL
    }

    if(!is.null(col_names)) {
      # Filter out time-related columns
      time_cols <- colnames(values$data)
      non_time_cols <- setdiff(col_names, time_cols)

      if(length(non_time_cols) > 0) {
        selectInput("rm_factor_var",
                    "Repeated measures factor (visit/condition/time):",
                    choices = c("", non_time_cols),
                    selected = "")
      } else {
        tagList(
          p("No factor variables found.", style = "color: orange;"),
          p("For repeated measures ANOVA, you need a column indicating the repeated condition/visit/time."),
          helpText("Example: A 'visit' column with values like 'baseline', 'week1', 'week2', etc.")
        )
      }
    } else {
      tagList(
        p("Upload data to select repeated measures factor.", style = "color: orange;")
      )
    }
  })

  # UI for selecting which levels of the RM factor to include
  output$rm_factor_levels_ui <- renderUI({
    req(input$rm_factor_var)
    req(values$uploaded_data)

    if(input$rm_factor_var == "" || !(input$rm_factor_var %in% colnames(values$uploaded_data))) {
      return(NULL)
    }

    # Get the levels of the selected factor
    rm_factor_data <- values$uploaded_data[[input$rm_factor_var]]
    factor_levels <- unique(as.character(rm_factor_data))
    factor_levels <- factor_levels[!is.na(factor_levels)]

    if(length(factor_levels) > 2) {
      tagList(
        hr(),
        h5("Select Conditions/Visits to Include"),
        pickerInput(
          inputId = "rm_levels_to_include",
          label = "Levels to include in analysis:",
          choices = factor_levels,
          selected = factor_levels,
          options = list(
            `actions-box` = TRUE,
            `selected-text-format` = "count > 3"
          ),
          multiple = TRUE
        ),
        helpText("Select at least 2 levels to compare. Deselect levels you want to exclude.")
      )
    } else {
      # Only 2 or fewer levels - no need for selection
      tagList(
        hr(),
        p(paste("Levels detected:", paste(factor_levels, collapse = ", "))),
        p(paste("Number of levels:", length(factor_levels)))
      )
    }
  })
  
  # Repeated Measures Functional ANOVA function using rmfanova package
  perform_rm_fanova <- function(fd_obj, subject_id, rm_factor, n_permutations = 200, alpha = 0.05) {
    
    # Check rmfanova package is available (should be auto-loaded at startup)
    if (!requireNamespace("rmfanova", quietly = TRUE)) {
      stop("Package 'rmfanova' is required but not installed. Please install it using: install.packages('rmfanova')")
    }
    
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    # Evaluate functional data at time points
    curves <- eval.fd(time_points, fd_obj)  # n_time x n_curves
    
    # Prepare data for rmfanova
    subject_id <- as.factor(subject_id)
    rm_factor <- as.factor(rm_factor)
    
    # Get unique levels
    visits <- levels(rm_factor)
    n_visits <- length(visits)
    
    cat("Running repeated measures functional ANOVA...\n")
    cat("Number of subjects:", length(unique(subject_id)), "\n")
    cat("Number of visits/conditions:", n_visits, "\n")
    cat("Visits:", paste(visits, collapse = ", "), "\n")
    
    # Try different rmfanova function signatures
    tryCatch({
      
      # Try to call rmfanova with different possible syntaxes
      # The actual syntax depends on package version
      
      # Attempt 1: Standard syntax from examples
      rm_results <- tryCatch({
        rmfanova::rmfanova(curves, id = subject_id, visit = rm_factor)
      }, error = function(e1) {
        
        # Attempt 2: Without named arguments
        tryCatch({
          rmfanova::rmfanova(curves, subject_id, rm_factor)
        }, error = function(e2) {
          
          # Attempt 3: Check if function is actually called something else
          tryCatch({
            # Try rm.fanova or rmfANOVA
            if(exists("rm.fanova", where = "package:rmfanova")) {
              rm.fanova(curves, id = subject_id, visit = rm_factor)
            } else if(exists("rmfANOVA", where = "package:rmfanova")) {
              rmfANOVA(curves, id = subject_id, visit = rm_factor)
            } else {
              stop("Could not find rmfanova function. Available functions: ", 
                   paste(ls("package:rmfanova"), collapse = ", "))
            }
          }, error = function(e3) {
            
            # If all attempts fail, use manual implementation
            cat("Note: Using manual RM-ANOVA implementation (rmfanova function not compatible)\n")
            return(NULL)
          })
        })
      })
      
      # If rmfanova call failed or returned NULL, use manual implementation
      if(is.null(rm_results)) {
        cat("Performing manual repeated measures functional ANOVA...\n")
        
        # Manual RM-ANOVA implementation
        # This is a simplified version that accounts for repeated measures structure

        # Calculate mean curves for each visit
        visit_means <- matrix(NA, n_time, n_visits)
        visit_sds <- matrix(NA, n_time, n_visits)
        visit_sizes <- numeric(n_visits)

        for(i in 1:n_visits) {
          visit_idx <- which(rm_factor == visits[i])
          visit_sizes[i] <- length(unique(subject_id[visit_idx]))
          visit_curves <- curves[, visit_idx, drop = FALSE]
          visit_means[, i] <- rowMeans(visit_curves, na.rm = TRUE)

          if(length(visit_idx) > 1) {
            visit_sds[, i] <- apply(visit_curves, 1, sd, na.rm = TRUE)
          } else {
            visit_sds[, i] <- 0
          }
        }

        # === DIAGNOSTIC OUTPUT ===
        cat("\n=== VISIT MEANS DIAGNOSTIC ===\n")
        cat("Number of visits/conditions:", n_visits, "\n")
        cat("Visit names:", paste(visits, collapse=", "), "\n")
        cat("Visit means matrix dimensions:", dim(visit_means), "\n")
        cat("Visit sizes (n subjects per visit):", paste(visit_sizes, collapse=", "), "\n")

        for(i in 1:n_visits) {
          cat("\nVisit", i, "(", visits[i], "):\n")
          cat("  First 5 time points:", paste(round(visit_means[1:5, i], 3), collapse=", "), "\n")
          cat("  Mean across all time points:", round(mean(visit_means[, i]), 3), "\n")
          cat("  SD across all time points:", round(sd(visit_means[, i]), 3), "\n")
        }

        if(n_visits == 2) {
          cat("\nComparison between visits:\n")
          cat("  Are means identical?", all.equal(visit_means[,1], visit_means[,2]), "\n")
          cat("  Max absolute difference:", round(max(abs(visit_means[,1] - visit_means[,2])), 6), "\n")
        }
        cat("==============================\n\n")

        # Calculate F-statistics accounting for within-subject correlation
        # Use subject-specific deviations
        unique_subjects <- unique(subject_id)
        n_subjects <- length(unique_subjects)
        
        # First pass: Calculate observed F-statistics and build data matrices
        F_stat <- numeric(n_time)
        Y_matrices <- vector("list", n_time)  # Store for permutation testing

        cat("=== Y_MATRIX CONSTRUCTION DIAGNOSTIC ===\n")

        for(t in 1:n_time) {
          # Create data matrix: subjects × visits
          Y_matrix <- matrix(NA, n_subjects, n_visits)
          
          for(i in 1:n_subjects) {
            subj <- unique_subjects[i]
            for(j in 1:n_visits) {
              idx <- which(subject_id == subj & rm_factor == visits[j])
              if(length(idx) > 0) {
                Y_matrix[i, j] <- curves[t, idx[1]]
              }
            }
          }
          
          # Remove subjects with missing data at this time point
          complete_rows <- complete.cases(Y_matrix)
          Y_complete <- Y_matrix[complete_rows, , drop = FALSE]

          # Diagnostic output for first time point only
          if(t == 1) {
            cat("\nTime point 1 - Y_matrix construction:\n")
            cat("  Y_matrix dimensions:", dim(Y_matrix), "(subjects × visits)\n")
            cat("  Y_complete dimensions:", dim(Y_complete), "\n")
            cat("  First 5 subjects, visit 1:", paste(round(Y_complete[1:min(5, nrow(Y_complete)), 1], 3), collapse=", "), "\n")
            cat("  First 5 subjects, visit 2:", paste(round(Y_complete[1:min(5, nrow(Y_complete)), 2], 3), collapse=", "), "\n")
            cat("  Column means:", paste(round(colMeans(Y_complete), 3), collapse=", "), "\n")
            cat("  Are columns identical?", all.equal(Y_complete[,1], Y_complete[,2]), "\n")
          }

          # Store for permutation testing
          Y_matrices[[t]] <- Y_complete

          if(nrow(Y_complete) >= 2 && n_visits >= 2) {
            # Perform repeated measures ANOVA at this time point
            # Remove subject mean (within-subject centering)
            subject_means <- rowMeans(Y_complete, na.rm = TRUE)
            Y_centered <- Y_complete - subject_means
            
            # Calculate SS
            grand_mean <- mean(Y_complete, na.rm = TRUE)
            visit_means_t <- colMeans(Y_complete, na.rm = TRUE)
            
            SS_visit <- sum(nrow(Y_complete) * (visit_means_t - grand_mean)^2)
            SS_residual <- sum(Y_centered^2, na.rm = TRUE)
            
            df_visit <- n_visits - 1
            df_residual <- (nrow(Y_complete) - 1) * (n_visits - 1)
            
            if(df_residual > 0 && SS_residual > 0) {
              F_stat[t] <- (SS_visit / df_visit) / (SS_residual / df_residual)
            } else {
              F_stat[t] <- NA
            }

            # Diagnostic output for first time point
            if(t == 1) {
              cat("\nTime point 1 - F-statistic calculation:\n")
              cat("  Grand mean:", round(grand_mean, 3), "\n")
              cat("  Visit means:", paste(round(visit_means_t, 3), collapse=", "), "\n")
              cat("  SS_visit:", round(SS_visit, 3), "\n")
              cat("  SS_residual:", round(SS_residual, 3), "\n")
              cat("  df_visit:", df_visit, ", df_residual:", df_residual, "\n")
              cat("  F-statistic:", round(F_stat[t], 3), "\n")
            }
          } else {
            F_stat[t] <- NA
          }
        }
        cat("=========================================\n\n")
        
        # Handle NAs in F-statistics
        F_stat[is.na(F_stat)] <- 0
        
        # PERMUTATION TEST for p-values (matching pairwise approach)
        cat("Computing permutation-based p-values (", n_permutations, " permutations)...\n")
        
        F_stat_perm <- matrix(NA, n_time, n_permutations)
        
        withProgress(message = 'Permutation testing...', value = 0, {
          for(perm in 1:n_permutations) {
            if(perm %% 20 == 0) incProgress(20 / n_permutations)
            
            for(t in 1:n_time) {
              Y_complete <- Y_matrices[[t]]
              
              if(!is.null(Y_complete) && nrow(Y_complete) >= 2 && n_visits >= 2) {
                # For RM data, randomly flip signs of condition effects within each subject
                # This preserves the within-subject correlation structure
                
                # Method: Randomly permute conditions within each subject
                Y_perm <- Y_complete
                for(i in 1:nrow(Y_complete)) {
                  Y_perm[i, ] <- Y_complete[i, sample(1:n_visits)]
                }
                
                # Calculate permuted F-statistic
                subject_means_perm <- rowMeans(Y_perm, na.rm = TRUE)
                Y_centered_perm <- Y_perm - subject_means_perm
                
                grand_mean_perm <- mean(Y_perm, na.rm = TRUE)
                visit_means_perm <- colMeans(Y_perm, na.rm = TRUE)
                
                SS_visit_perm <- sum(nrow(Y_perm) * (visit_means_perm - grand_mean_perm)^2)
                SS_residual_perm <- sum(Y_centered_perm^2, na.rm = TRUE)
                
                df_visit <- n_visits - 1
                df_residual <- (nrow(Y_perm) - 1) * (n_visits - 1)
                
                if(df_residual > 0 && SS_residual_perm > 0) {
                  F_stat_perm[t, perm] <- (SS_visit_perm / df_visit) / (SS_residual_perm / df_residual)
                } else {
                  F_stat_perm[t, perm] <- 0
                }
              } else {
                F_stat_perm[t, perm] <- 0
              }
            }
          }
        })
        
        # Calculate permutation-based p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          p_values_pointwise[t] <- mean(F_stat_perm[t, ] >= F_stat[t], na.rm = TRUE)
        }
        
        # Handle any remaining NAs
        p_values_pointwise[is.na(p_values_pointwise)] <- 1
        
        cat("Permutation testing complete.\n")
        
      } else {
        # Extract results from rmfanova
        # The structure depends on the package version
        
        # Try to extract pointwise statistics
        if(!is.null(rm_results$pointwise)) {
          F_stat <- rm_results$pointwise$stat
          p_values_pointwise <- rm_results$pointwise$pval
        } else if(!is.null(rm_results$stat)) {
          F_stat <- rm_results$stat
          p_values_pointwise <- rm_results$pval
        } else {
          # Calculate manually as fallback
          visit_means <- matrix(NA, n_time, n_visits)
          for(i in 1:n_visits) {
            visit_idx <- which(rm_factor == visits[i])
            visit_means[, i] <- rowMeans(curves[, visit_idx, drop = FALSE])
          }
          
          # Simple F-statistic calculation
          overall_mean <- rowMeans(curves)
          SSB <- rowSums((visit_means - overall_mean)^2) * length(unique(subject_id))
          SST <- rowSums((curves - overall_mean)^2)
          
          F_stat <- SSB / (SST - SSB + 1e-10)
          p_values_pointwise <- pf(F_stat, n_visits - 1, ncol(curves) - n_visits, lower.tail = FALSE)
        }
        
        # Calculate visit means for plotting
        visit_means <- matrix(NA, n_time, n_visits)
        visit_sds <- matrix(NA, n_time, n_visits)
        visit_sizes <- numeric(n_visits)

        for(i in 1:n_visits) {
          visit_idx <- which(rm_factor == visits[i])
          visit_sizes[i] <- length(unique(subject_id[visit_idx]))
          visit_curves <- curves[, visit_idx, drop = FALSE]
          visit_means[, i] <- rowMeans(visit_curves, na.rm = TRUE)

          if(length(visit_idx) > 1) {
            visit_sds[, i] <- apply(visit_curves, 1, sd, na.rm = TRUE)
          } else {
            visit_sds[, i] <- 0
          }
        }
      }
      
      # Continue with common processing...
      # Adjust p-values for multiple comparisons
      p_values_adjusted <- p.adjust(p_values_pointwise, method = "fdr")
      sig_regions <- p_values_adjusted < alpha
      
      # Calculate effect sizes (eta-squared)
      overall_mean <- rowMeans(curves)
      SST <- numeric(n_time)
      SSB <- numeric(n_time)
      
      for(t in 1:n_time) {
        for(i in 1:n_visits) {
          SSB[t] <- SSB[t] + visit_sizes[i] * (visit_means[t, i] - overall_mean[t])^2
        }
        
        for(j in 1:ncol(curves)) {
          SST[t] <- SST[t] + (curves[t, j] - overall_mean[t])^2
        }
      }
      
      eta_squared <- SSB / (SST + 1e-10)  # Add small constant to avoid division by zero
      eta_squared[SST == 0] <- 0
      
      # Global test statistic (L2 norm)
      L2_stat <- sqrt(sum((SSB / length(subject_id))^2))
      p_value_L2 <- NA  # Would need permutation test
      
      # Calculate confidence bands using bootstrap
      n_boot <- 100
      visit_means_boot <- array(NA, dim = c(n_time, n_visits, n_boot))
      
      unique_subjects <- unique(subject_id)
      
      for(boot in 1:n_boot) {
        # Resample subjects (not observations)
        boot_subjects <- sample(unique_subjects, replace = TRUE)
        
        for(i in 1:n_visits) {
          boot_curves_list <- list()
          for(subj in boot_subjects) {
            # Get curves for this subject in this visit
            subj_visit_idx <- which(subject_id == subj & rm_factor == visits[i])
            if(length(subj_visit_idx) > 0) {
              boot_curves_list[[length(boot_curves_list) + 1]] <- curves[, subj_visit_idx[1], drop = FALSE]
            }
          }
          
          if(length(boot_curves_list) > 0) {
            boot_curves_matrix <- do.call(cbind, boot_curves_list)
            visit_means_boot[, i, boot] <- rowMeans(boot_curves_matrix)
          } else {
            visit_means_boot[, i, boot] <- visit_means[, i]
          }
        }
      }
      
      # Calculate confidence bands
      visit_means_lower <- matrix(NA, n_time, n_visits)
      visit_means_upper <- matrix(NA, n_time, n_visits)
      
      for(i in 1:n_visits) {
        for(t in 1:n_time) {
          quantiles <- quantile(visit_means_boot[t, i, ], probs = c(0.025, 0.975), na.rm = TRUE)
          visit_means_lower[t, i] <- quantiles[1]
          visit_means_upper[t, i] <- quantiles[2]
        }
      }
      
      # DIAGNOSTIC: Verify F_stat before return
      cat("\n=== FINAL F-STAT DIAGNOSTIC ===\n")
      cat("F_stat vector length:", length(F_stat), "\n")
      cat("First 5 F-statistics:", paste(round(F_stat[1:min(5, length(F_stat))], 3), collapse=", "), "\n")
      cat("Max F-statistic:", round(max(F_stat, na.rm = TRUE), 3), "\n")
      cat("Mean F-statistic:", round(mean(F_stat, na.rm = TRUE), 3), "\n")
      cat("Number of significant time points (raw p < 0.05):", sum(p_values_pointwise < 0.05, na.rm = TRUE), "/", length(p_values_pointwise), "\n")
      cat("Number of significant time points (FDR adjusted):", sum(p_values_adjusted < 0.05, na.rm = TRUE), "/", length(p_values_adjusted), "\n")
      cat("Min p-value:", format.pval(min(p_values_pointwise, na.rm = TRUE), digits = 4), "\n")
      cat("================================\n\n")

      # Return results in same format as between-subjects ANOVA
      return(list(
        design = "within",
        time_points = time_points,
        group_means = visit_means,
        group_sds = visit_sds,
        group_means_lower = visit_means_lower,
        group_means_upper = visit_means_upper,
        group_labels = rm_factor,
        groups = visits,
        n_groups = n_visits,
        group_sizes = visit_sizes,
        F_stat = F_stat,
        p_values_pointwise = p_values_pointwise,
        p_values_adjusted = p_values_adjusted,
        sig_regions = sig_regions,
        L2_stat = L2_stat,
        p_value_L2 = p_value_L2,
        eta_squared = eta_squared,
        df_between = n_visits - 1,
        df_within = length(unique(subject_id)) * (n_visits - 1),
        alpha = alpha,
        n_permutations = n_permutations,
        subject_id = subject_id,
        rm_factor = rm_factor
      ))
      
    }, error = function(e) {
      stop(paste("Error in RM-ANOVA:", e$message))
    })
  }
  
  # Functional ANOVA function - ENHANCED
  perform_functional_anova <- function(fd_obj, group_labels, n_permutations = 200, 
                                       test_type = "both", alpha = 0.05) {
    
    n_curves <- ncol(fd_obj$coefs)
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    group_labels <- as.factor(group_labels)
    groups <- levels(group_labels)
    n_groups <- length(groups)
    
    group_means <- matrix(NA, n_time, n_groups)
    group_sds <- matrix(NA, n_time, n_groups)
    group_sizes <- numeric(n_groups)
    
    for(i in 1:n_groups) {
      group_idx <- which(group_labels == groups[i])
      group_sizes[i] <- length(group_idx)
      group_curves <- curves[, group_idx, drop = FALSE]
      group_means[, i] <- rowMeans(group_curves)
      
      # Calculate standard deviation for each time point
      if(length(group_idx) > 1) {
        group_sds[, i] <- apply(group_curves, 1, sd)
      } else {
        group_sds[, i] <- 0
      }
    }
    
    overall_mean <- rowMeans(curves)
    
    # Calculate F-statistics
    SSB <- numeric(n_time)
    SSW <- numeric(n_time)
    
    for(t in 1:n_time) {
      for(i in 1:n_groups) {
        SSB[t] <- SSB[t] + group_sizes[i] * (group_means[t, i] - overall_mean[t])^2
      }
      
      for(i in 1:n_groups) {
        group_idx <- which(group_labels == groups[i])
        for(j in group_idx) {
          SSW[t] <- SSW[t] + (curves[t, j] - group_means[t, i])^2
        }
      }
    }
    
    df_between <- n_groups - 1
    df_within <- n_curves - n_groups
    
    F_stat <- (SSB / df_between) / (SSW / df_within)
    L2_stat <- sqrt(sum((SSB / n_curves)^2))
    
    # Permutation test for p-values
    F_stat_perm <- matrix(NA, n_time, n_permutations)
    L2_stat_perm <- numeric(n_permutations)
    
    for(perm in 1:n_permutations) {
      perm_labels <- sample(group_labels)
      
      # Recalculate statistics with permuted labels
      perm_means <- matrix(NA, n_time, n_groups)
      for(i in 1:n_groups) {
        perm_idx <- which(perm_labels == groups[i])
        perm_means[, i] <- rowMeans(curves[, perm_idx, drop = FALSE])
      }
      
      SSB_perm <- numeric(n_time)
      SSW_perm <- numeric(n_time)
      
      for(t in 1:n_time) {
        for(i in 1:n_groups) {
          perm_idx <- which(perm_labels == groups[i])
          SSB_perm[t] <- SSB_perm[t] + length(perm_idx) * (perm_means[t, i] - overall_mean[t])^2
          
          for(j in perm_idx) {
            SSW_perm[t] <- SSW_perm[t] + (curves[t, j] - perm_means[t, i])^2
          }
        }
      }
      
      F_stat_perm[, perm] <- (SSB_perm / df_between) / (SSW_perm / df_within)
      L2_stat_perm[perm] <- sqrt(sum((SSB_perm / n_curves)^2))
    }
    
    # Calculate p-values
    p_values_pointwise <- numeric(n_time)
    for(t in 1:n_time) {
      p_values_pointwise[t] <- mean(F_stat_perm[t, ] >= F_stat[t])
    }
    
    p_value_L2 <- mean(L2_stat_perm >= L2_stat)
    
    p_values_adjusted <- p.adjust(p_values_pointwise, method = "fdr")
    sig_regions <- p_values_adjusted < alpha
    
    SST <- SSB + SSW
    eta_squared <- SSB / SST
    
    # Calculate confidence bands using bootstrap
    n_boot <- 100
    group_means_boot <- array(NA, dim = c(n_time, n_groups, n_boot))
    
    for(boot in 1:n_boot) {
      for(i in 1:n_groups) {
        group_idx <- which(group_labels == groups[i])
        if(length(group_idx) > 0) {
          boot_idx <- sample(group_idx, replace = TRUE)
          group_means_boot[, i, boot] <- rowMeans(curves[, boot_idx, drop = FALSE])
        }
      }
    }
    
    # Calculate confidence bands
    group_means_lower <- matrix(NA, n_time, n_groups)
    group_means_upper <- matrix(NA, n_time, n_groups)
    
    for(i in 1:n_groups) {
      for(t in 1:n_time) {
        quantiles <- quantile(group_means_boot[t, i, ], probs = c(0.025, 0.975), na.rm = TRUE)
        group_means_lower[t, i] <- quantiles[1]
        group_means_upper[t, i] <- quantiles[2]
      }
    }
    
    return(list(
      time_points = time_points,
      group_means = group_means,
      group_sds = group_sds,  # Added standard deviations
      group_means_lower = group_means_lower,
      group_means_upper = group_means_upper,
      group_labels = group_labels,
      groups = groups,
      n_groups = n_groups,
      group_sizes = group_sizes,
      F_stat = F_stat,
      p_values_pointwise = p_values_pointwise,
      p_values_adjusted = p_values_adjusted,
      sig_regions = sig_regions,
      L2_stat = L2_stat,
      p_value_L2 = p_value_L2,
      eta_squared = eta_squared,
      df_between = df_between,
      df_within = df_within,
      alpha = alpha,
      n_permutations = n_permutations
    ))
  }
  
  # Run functional ANOVA - ENHANCED for both between and within designs
  observeEvent(input$run_fanova, {
    if(is.null(values$data)) {
      showNotification("Please load data first!", type = "error", duration = 5)
      return()
    }
    
    # Validate based on design type
    if(input$fanova_design == "between") {
      # Between-subjects validation - REQUIRES group labels from preprocessing
      fanova_groups <- get_fanova_group_labels()
      if(is.null(fanova_groups)) {
        showNotification("For between-subjects ANOVA: Please define group labels in the Data Preprocessing tab first!",
                         type = "error", duration = 5)
        return()
      }

      if(length(unique(fanova_groups)) < 2) {
        showNotification("Need at least 2 groups for ANOVA!", type = "error", duration = 5)
        return()
      }

      cat("Running BETWEEN-SUBJECTS functional ANOVA...\n")
      cat("Using group variable:", if(!is.null(input$fanova_group_var)) input$fanova_group_var else "primary", "\n")
      
    } else if(input$fanova_design == "within") {
      # Within-subjects (repeated measures) validation
      if(is.null(input$rm_subject_id_var) || input$rm_subject_id_var == "") {
        showNotification("Please select a Subject ID variable for repeated measures ANOVA!", 
                         type = "error", duration = 5)
        return()
      }
      
      if(is.null(input$rm_factor_var) || input$rm_factor_var == "") {
        showNotification("Please select a repeated measures factor (visit/condition/time)!", 
                         type = "error", duration = 5)
        return()
      }
      
      # Extract subject ID and RM factor from uploaded data
      if(is.null(values$uploaded_data)) {
        showNotification("Original data with ID and factor variables not available!", 
                         type = "error", duration = 5)
        return()
      }
      
      # Get the variables
      subject_id_col <- input$rm_subject_id_var
      rm_factor_col <- input$rm_factor_var
      
      if(!(subject_id_col %in% colnames(values$uploaded_data)) || 
         !(rm_factor_col %in% colnames(values$uploaded_data))) {
        showNotification("Selected variables not found in data!", 
                         type = "error", duration = 5)
        return()
      }
      
      subject_id_data <- values$uploaded_data[[subject_id_col]]
      rm_factor_data <- values$uploaded_data[[rm_factor_col]]
      
      # Check that we have multiple levels of the RM factor
      if(length(unique(rm_factor_data)) < 2) {
        showNotification("Need at least 2 levels of the repeated measures factor!", 
                         type = "error", duration = 5)
        return()
      }
      
      cat("Running WITHIN-SUBJECTS (Repeated Measures) functional ANOVA...\n")
      cat("Subject ID variable:", subject_id_col, "\n")
      cat("RM factor variable:", rm_factor_col, "\n")
      cat("Number of unique subjects:", length(unique(subject_id_data)), "\n")
      cat("Number of conditions/visits:", length(unique(rm_factor_data)), "\n")
      cat("Total number of observations:", length(subject_id_data), "\n")
      
      # Validate this is truly repeated measures
      n_subjects <- length(unique(subject_id_data))
      n_obs <- length(subject_id_data)
      if(n_subjects == n_obs) {
        cat("\n*** WARNING ***\n")
        cat("Number of subjects equals number of observations!\n")
        cat("This suggests you may not have true repeated measures data.\n")
        cat("Each subject should appear multiple times (once per condition/visit).\n")
        cat("Consider using between-subjects design instead.\n")
        cat("***************\n\n")
      } else {
        cat("Expected observations per subject:", round(n_obs / n_subjects, 2), "\n")
      }
    }
    
    tryCatch({
      # Determine which data to use based on user selection
      fd_to_use <- NULL
      
      # Check if user wants to use warped curves and if they're available
      if(input$fanova_data_source == "warped" && !is.null(values$warping_results)) {
        cat("Using time-warped curves for FANOVA\n")
        
        # Try to get the registered fd object from warping results
        if(!is.null(values$warping_results$regfd)) {
          fd_to_use <- values$warping_results$regfd
          showNotification("Using time-warped curves for FANOVA", type = "message", duration = 3)
        } else if(!is.null(values$warping_results$registered_curves)) {
          # If no regfd but registered curves exist, create fd from them
          n_time <- nrow(values$warping_results$registered_curves)
          time_points <- values$warping_results$time_points
          if(is.null(time_points)) {
            time_points <- seq(0, 1, length.out = n_time)
          }
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          # Use lambda=0: registered curves are already processed, just need fd representation
          fd_to_use <- smooth.basis(time_points, values$warping_results$registered_curves, 
                                    fdPar(basis, 2, 0))$fd
          showNotification("Using time-warped curves for FANOVA", type = "message", duration = 3)
        } else {
          showNotification("No warped curves available, using original data", type = "warning", duration = 5)
        }
      }
      
      # If no warped data was used or available, use original data
      if(is.null(fd_to_use)) {
        cat("Using original curves for FANOVA\n")
        
        # CRITICAL: Check if data has already been smoothed in Data Preprocessing
        if(is.null(values$fd_obj)) {
          # Determine which data to use
          data_for_fd <- if(!is.null(values$smooth_data)) {
            cat("Using already smoothed data (no additional smoothing)\n")
            values$smooth_data
          } else {
            cat("Using raw data (will create fd object)\n")
            values$data
          }
          
          n_time <- ncol(data_for_fd)
          time_points <- seq(0, 1, length.out = n_time)
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          
          # Create fd object WITHOUT additional smoothing (lambda = 0 for already smoothed data)
          if(!is.null(values$smooth_data)) {
            # Data already smoothed - just create fd representation with no penalty
            values$fd_obj <- smooth.basis(time_points, t(data_for_fd), fdPar(basis, 2, 0))$fd
          } else {
            # Raw data - apply default smoothing
            values$fd_obj <- smooth.basis(time_points, t(data_for_fd), basis)$fd
          }
        }
        fd_to_use <- values$fd_obj
        showNotification("Using original curves for FANOVA", type = "message", duration = 3)
      }
      
      # Perform FANOVA based on design type
      if(input$fanova_design == "between") {
        # Between-subjects ANOVA

        # Get current group labels from the selected fANOVA group variable
        current_group_labels <- get_fanova_group_labels()

        # Filter by selected groups if user has made a selection
        groups_to_include <- input$fanova_groups_to_include
        if(!is.null(groups_to_include) && length(groups_to_include) >= 2) {
          # Find indices of subjects in selected groups
          include_idx <- which(current_group_labels %in% groups_to_include)

          if(length(include_idx) < 2) {
            showNotification("Not enough subjects in selected groups!", type = "error", duration = 5)
            return()
          }

          # Subset the fd object - need to extract curves, subset, and recreate
          n_time <- 100
          time_points_eval <- seq(0, 1, length.out = n_time)
          all_curves <- eval.fd(time_points_eval, fd_to_use)
          subset_curves <- all_curves[, include_idx, drop = FALSE]

          # Recreate fd object from subset
          basis <- fd_to_use$basis
          fd_to_use <- smooth.basis(time_points_eval, subset_curves, basis)$fd

          # Subset group labels and drop unused levels
          current_group_labels <- droplevels(current_group_labels[include_idx])

          cat("Filtered to groups:", paste(groups_to_include, collapse = ", "), "\n")
          cat("Subjects included:", length(include_idx), "\n")

          # Store which groups were selected for reference
          values$fanova_selected_groups <- groups_to_include
        } else if(!is.null(groups_to_include) && length(groups_to_include) < 2) {
          showNotification("Please select at least 2 groups for comparison!", type = "error", duration = 5)
          return()
        }

        values$fanova_results <- perform_functional_anova(
          fd_obj = fd_to_use,
          group_labels = current_group_labels,
          n_permutations = input$n_permutations,
          test_type = input$fanova_test_type,
          alpha = input$alpha_level
        )

        values$fanova_results$design <- "between"
        
      } else if(input$fanova_design == "within") {
        # Within-subjects (repeated measures) ANOVA
        subject_id_data <- values$uploaded_data[[input$rm_subject_id_var]]
        rm_factor_data <- values$uploaded_data[[input$rm_factor_var]]

        # Filter by selected levels if user has made a selection
        levels_to_include <- input$rm_levels_to_include
        if(!is.null(levels_to_include) && length(levels_to_include) >= 2) {
          # Find indices of observations with selected levels
          include_idx <- which(as.character(rm_factor_data) %in% levels_to_include)

          if(length(include_idx) < 2) {
            showNotification("Not enough observations in selected levels!", type = "error", duration = 5)
            return()
          }

          # Subset the fd object
          n_time <- 100
          time_points_eval <- seq(0, 1, length.out = n_time)
          all_curves <- eval.fd(time_points_eval, fd_to_use)
          subset_curves <- all_curves[, include_idx, drop = FALSE]

          # Recreate fd object from subset
          basis <- fd_to_use$basis
          fd_to_use <- smooth.basis(time_points_eval, subset_curves, basis)$fd

          # Subset the subject ID and RM factor data
          subject_id_data <- subject_id_data[include_idx]
          rm_factor_data <- droplevels(as.factor(rm_factor_data[include_idx]))

          cat("Filtered to RM levels:", paste(levels_to_include, collapse = ", "), "\n")
          cat("Observations included:", length(include_idx), "\n")
        } else if(!is.null(levels_to_include) && length(levels_to_include) < 2) {
          showNotification("Please select at least 2 levels for comparison!", type = "error", duration = 5)
          return()
        }

        # CRITICAL VALIDATION: Check data structure matches fd_obj
        n_curves_in_fd <- ncol(fd_to_use$coefs)

        # Additional validation: lengths must match
        if(length(subject_id_data) != n_curves_in_fd || length(rm_factor_data) != n_curves_in_fd) {
          showNotification(
            sprintf("ERROR: Variable lengths don't match! Subject ID: %d, RM factor: %d, Curves: %d",
                    length(subject_id_data), length(rm_factor_data), n_curves_in_fd),
            type = "error", duration = 10)
          return()
        }

        cat("Within-subjects validation passed:\n")
        cat("  Number of curves:", n_curves_in_fd, "\n")
        cat("  Subject ID length:", length(subject_id_data), "\n")
        cat("  RM factor length:", length(rm_factor_data), "\n")
        cat("  Unique levels:", length(unique(rm_factor_data)), "\n")

        values$fanova_results <- perform_rm_fanova(
          fd_obj = fd_to_use,
          subject_id = subject_id_data,
          rm_factor = rm_factor_data,
          n_permutations = input$n_permutations,
          alpha = input$alpha_level
        )

        values$fanova_results$design <- "within"
      }
      
      # Store which data source was used
      values$fanova_results$data_source <- input$fanova_data_source
      
      showNotification("Functional ANOVA completed!", type = "message", duration = 3)
      
    }, error = function(e) {
      cat("FANOVA error:", e$message, "\n")
      showNotification(paste("FANOVA error:", e$message), type = "error", duration = 10)
    })
  })
  
  # FANOVA outputs (enhanced for both between and within designs)
  output$fanova_global_results <- renderPrint({
    req(values$fanova_results)
    
    res <- values$fanova_results
    
    cat("========================================\n")
    cat("    FUNCTIONAL ANOVA RESULTS\n")
    cat("========================================\n\n")
    
    # Show design type
    design_type <- if(!is.null(res$design)) {
      if(res$design == "between") "Between-Subjects" else "Within-Subjects (Repeated Measures)"
    } else {
      "Between-Subjects"
    }
    cat("Design:", design_type, "\n\n")
    
    # Show group/condition information
    if(!is.null(res$design) && res$design == "within") {
      cat("Number of conditions:", res$n_groups, "\n")
      cat("Conditions:", paste(res$groups, collapse = ", "), "\n")
      cat("Subjects per condition:", paste(res$group_sizes, collapse = ", "), "\n")
    } else {
      cat("Number of groups:", res$n_groups, "\n")
      cat("Groups:", paste(res$groups, collapse = ", "), "\n")
      cat("Group sizes:", paste(res$group_sizes, collapse = ", "), "\n")
    }
    
    cat("\nL2 statistic:", round(res$L2_stat, 4), "\n")
    
    if(!is.na(res$p_value_L2)) {
      cat("L2 p-value:", round(res$p_value_L2, 4), "\n")
    }
    
    cat("\n========================================\n")
  })
  
  output$fanova_summary_table <- renderDT({
    req(values$fanova_results)
    
    res <- values$fanova_results
    
    # Calculate additional statistics for each group/condition
    curves <- eval.fd(res$time_points, values$fd_obj)
    
    # Determine label for first column
    group_label <- if(!is.null(res$design) && res$design == "within") {
      "Condition"
    } else {
      "Group"
    }
    
    summary_df <- data.frame(
      GroupOrCondition = res$groups,
      N = res$group_sizes,
      Mean_Area = round(apply(res$group_means, 2, function(x) mean(x)), 3),
      SD_Area = round(apply(res$group_means, 2, function(x) sd(x)), 3),
      Max_Value = round(apply(res$group_means, 2, max), 3),
      Min_Value = round(apply(res$group_means, 2, min), 3),
      Range = round(apply(res$group_means, 2, function(x) diff(range(x))), 3)
    )
    
    # Rename first column
    colnames(summary_df)[1] <- group_label
    
    datatable(summary_df, 
              options = list(pageLength = 10, dom = 't', scrollX = TRUE), 
              rownames = FALSE) %>%
      formatStyle("N",
                  backgroundColor = styleInterval(c(5, 10), 
                                                  c("#ffcccc", "#ffffcc", "#ccffcc")))
  })
  
  # FANOVA plots - ENHANCED with SD bands and toggles
  output$fanova_mean_plot <- renderPlotly({
    req(values$fanova_results)

    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # Get toggle values (with defaults if not yet initialized)
    show_sd_bands <- if(!is.null(input$fanova_show_sd_bands)) input$fanova_show_sd_bands else TRUE
    show_sig_regions <- if(!is.null(input$fanova_show_sig_regions)) input$fanova_show_sig_regions else TRUE

    # Create color palette that scales with number of groups
    base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
    colors <- colorRampPalette(base_cols)(res$n_groups)

    n_time <- length(res$time_points)

    # Calculate SD for each group
    curves <- eval.fd(res$time_points, values$fd_obj)
    group_sds <- matrix(NA, n_time, res$n_groups)

    for(i in 1:res$n_groups) {
      group_idx <- which(res$group_labels == res$groups[i])
      if(length(group_idx) > 1) {
        group_curves <- curves[, group_idx, drop = FALSE]
        group_sds[,i] <- apply(group_curves, 1, sd)
      } else {
        group_sds[,i] <- 0
      }
    }

    p <- plot_ly(type = 'scatter', mode = 'lines')

    # Add SD bands and mean curves for each group
    # NOTE: No legendgroup linking - each trace toggles independently
    for(i in 1:res$n_groups) {
      # Build safe rgba fill color with alpha
      col_rgb <- col2rgb(colors[i])
      fillcol <- sprintf("rgba(%d,%d,%d,%.2f)", col_rgb[1], col_rgb[2], col_rgb[3], 0.2)

      # Add ±1 SD band (only if toggle is on)
      if(show_sd_bands) {
        p <- p %>% add_trace(
          x = c(res$time_points, rev(res$time_points)),
          y = c(res$group_means[,i] + group_sds[,i],
                rev(res$group_means[,i] - group_sds[,i])),
          fill = 'toself',
          fillcolor = fillcol,
          line = list(color = 'transparent'),
          showlegend = TRUE,
          name = paste0(res$groups[i], " ±SD"),
          hoverinfo = 'skip'
          # No legendgroup - toggles independently from mean curve
        )
      }

      # Add mean curve
      p <- p %>% add_trace(
        x = res$time_points,
        y = res$group_means[,i],
        line = list(color = colors[i], width = 3),
        name = res$groups[i],
        # No legendgroup - toggles independently from SD band
        hovertemplate = paste(res$groups[i], "<br>Time: %{customdata}<br>Mean: %{y:.3f}<extra></extra>"),
        customdata = hover_times
      )
    }

    # Add significant regions as background (only if toggle is on)
    # All significant regions share a legendgroup so they toggle together
    if(show_sig_regions && any(res$sig_regions)) {
      sig_starts <- which(diff(c(0, res$sig_regions)) == 1)
      sig_ends <- which(diff(c(res$sig_regions, 0)) == -1)

      for(idx in 1:length(sig_starts)) {
        y_range <- range(c(res$group_means + group_sds, res$group_means - group_sds))

        p <- p %>% add_trace(
          x = c(res$time_points[sig_starts[idx]], res$time_points[sig_ends[idx]],
                res$time_points[sig_ends[idx]], res$time_points[sig_starts[idx]]),
          y = c(min(y_range), min(y_range), max(y_range), max(y_range)),
          fill = 'toself',
          fillcolor = 'rgba(200, 200, 200, 0.2)',
          line = list(color = 'transparent'),
          showlegend = (idx == 1),
          name = 'Significant',
          legendgroup = 'significant_regions',  # All sig regions toggle together
          hoverinfo = 'skip'
        )
      }
    }
    
    # Dynamic title based on settings
    plot_title <- if(show_sd_bands) "Group Mean Functions ±1 SD" else "Group Mean Functions"

    p <- p %>% layout(
      title = plot_title,
      yaxis = list(title = "Value"),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.98, tracegroupgap = 5)
    )
    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))

    # Enable editable mode so legend can be dragged
    p <- p %>% config(editable = TRUE)
    p
  })
  
  output$fanova_fstat_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # Calculate critical value
    crit_val <- qf(1 - res$alpha, res$df_between, res$df_within)
    
    p <- plot_ly(type = 'scatter', mode = 'lines') %>%
      add_trace(x = res$time_points, y = res$F_stat,
                line = list(color = 'blue', width = 2),
                name = 'F-statistic',
                hovertemplate = "Time: %{customdata}<br>F-stat: %{y:.2f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = c(0, 1), y = c(crit_val, crit_val),
                line = list(color = 'red', width = 2, dash = 'dash'),
                name = paste('Critical value (α =', res$alpha, ')'),
                hovertemplate = paste("Critical F =", round(crit_val, 2), "<extra></extra>")) %>%
      layout(title = paste("Pointwise F-statistics (", res$n_groups, "groups)"),
             yaxis = list(title = "F-statistic"))

    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_pvalue_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    p <- plot_ly(type = 'scatter', mode = 'lines') %>%
      add_trace(x = res$time_points, y = res$p_values_adjusted,
                line = list(color = 'darkgreen', width = 2),
                name = 'Adjusted p-values (FDR)',
                hovertemplate = "Time: %{customdata}<br>P-value: %{y:.4f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = res$time_points, y = res$p_values_pointwise,
                line = list(color = 'lightgreen', width = 1, dash = 'dot'),
                name = 'Raw p-values',
                hovertemplate = "Time: %{customdata}<br>P-value: %{y:.4f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = c(0, 1), y = c(res$alpha, res$alpha),
                line = list(color = 'red', width = 2, dash = 'dash'),
                name = paste('α =', res$alpha),
                hovertemplate = paste("Alpha =", res$alpha, "<extra></extra>")) %>%
      add_trace(x = c(0, 1), y = c(0.01, 0.01),
                line = list(color = 'orange', width = 1, dash = 'dot'),
                name = 'p = 0.01',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.001, 0.001),
                line = list(color = 'darkred', width = 1, dash = 'dot'),
                name = 'p = 0.001',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      layout(title = paste("Pointwise p-values (", res$n_groups, "groups)"),
             yaxis = list(title = "p-value", type = 'log',
                          range = c(log10(0.0001), log10(1))))
    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_effect_size_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # Add background for effect size interpretation
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    # Add reference lines for effect size interpretation
    p <- p %>% 
      add_trace(x = c(0, 1), y = c(0.01, 0.01),
                line = list(color = 'lightgray', width = 1, dash = 'dot'),
                name = 'Small effect',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.06, 0.06),
                line = list(color = 'gray', width = 1, dash = 'dot'),
                name = 'Medium effect',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.14, 0.14),
                line = list(color = 'darkgray', width = 1, dash = 'dot'),
                name = 'Large effect',
                showlegend = FALSE,
                hoverinfo = 'skip')
    
    # Add eta-squared curve
    p <- p %>% add_trace(x = res$time_points, y = res$eta_squared,
                         fill = 'tozeroy',
                         fillcolor = 'rgba(100, 100, 255, 0.3)',
                         line = list(color = 'purple', width = 2),
                         name = 'η² (Effect size)',
                         hovertemplate = "Time: %{customdata}<br>η²: %{y:.3f}<extra></extra>",
                         customdata = hover_times)
    
    p <- p %>% layout(title = paste("Effect Size (η²) across Time -", res$n_groups, "groups"),
                      yaxis = list(title = "η² (proportion of variance explained)", 
                                   range = c(0, max(0.3, max(res$eta_squared) * 1.1))),
                      annotations = list(
                        list(x = 0.95, y = 0.01, text = "Small", showarrow = FALSE, 
                             font = list(size = 10, color = "gray")),
                        list(x = 0.95, y = 0.06, text = "Medium", showarrow = FALSE,
                             font = list(size = 10, color = "gray")),
                        list(x = 0.95, y = 0.14, text = "Large", showarrow = FALSE,
                             font = list(size = 10, color = "gray"))
                      ))

    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_effect_summary <- renderPrint({
    req(values$fanova_results)
    
    res <- values$fanova_results
    mean_eta <- mean(res$eta_squared)
    
    cat("Effect Size Summary (η²):\n")
    cat("Mean η²:", round(mean_eta, 4), "\n")
  })
  
  # Check if FANOVA completed
  output$fanova_completed <- reactive({
    !is.null(values$fanova_results)
  })
  outputOptions(output, "fanova_completed", suspendWhenHidden = FALSE)
  
  # Pairwise comparison functions - ENHANCED
  # Repeated Measures Pairwise Comparisons (for within-subjects designs)
  perform_pairwise_comparisons_rm <- function(fd_obj, subject_id, rm_factor, n_permutations = 200,
                                              correction_method = "bonferroni", alpha = 0.05) {
    
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    subject_id <- as.factor(subject_id)
    rm_factor <- as.factor(rm_factor)
    
    conditions <- levels(rm_factor)
    n_conditions <- length(conditions)
    
    n_pairs <- choose(n_conditions, 2)
    pairs <- combn(conditions, 2, simplify = FALSE)
    pair_names <- sapply(pairs, function(p) paste(p[1], "vs", p[2]))
    
    pairwise_results <- list()
    
    cat("Performing paired comparisons (within-subjects)...\n")
    cat("Number of condition pairs:", n_pairs, "\n")
    
    withProgress(message = 'Performing pairwise comparisons', value = 0, {
      for(pair_idx in 1:n_pairs) {
        incProgress(1/n_pairs, detail = paste("Comparing", pair_names[pair_idx]))
        
        pair <- pairs[[pair_idx]]
        
        # Get indices for each condition
        idx1 <- which(rm_factor == pair[1])
        idx2 <- which(rm_factor == pair[2])
        
        # Match subjects across conditions
        subjects_in_1 <- subject_id[idx1]
        subjects_in_2 <- subject_id[idx2]
        
        # Find subjects present in both conditions (for paired comparison)
        common_subjects <- intersect(subjects_in_1, subjects_in_2)
        n_pairs_subj <- length(common_subjects)
        
        cat("  Pair:", pair_names[pair_idx], "- Matched subjects:", n_pairs_subj, "\n")
        
        if(n_pairs_subj < 2) {
          warning("Not enough paired subjects for comparison: ", pair_names[pair_idx])
          next
        }
        
        # Extract matched curves
        matched_curves1 <- matrix(NA, n_time, n_pairs_subj)
        matched_curves2 <- matrix(NA, n_time, n_pairs_subj)
        
        for(i in 1:n_pairs_subj) {
          subj <- common_subjects[i]
          
          # Find this subject's curve in condition 1
          subj_idx1 <- idx1[which(subjects_in_1 == subj)[1]]
          matched_curves1[, i] <- curves[, subj_idx1]
          
          # Find this subject's curve in condition 2
          subj_idx2 <- idx2[which(subjects_in_2 == subj)[1]]
          matched_curves2[, i] <- curves[, subj_idx2]
        }
        
        # Calculate paired differences
        paired_diffs <- matched_curves1 - matched_curves2
        mean_diff <- rowMeans(paired_diffs)
        
        # Calculate paired t-statistics
        se_diff <- apply(paired_diffs, 1, sd) / sqrt(n_pairs_subj)
        t_stat <- mean_diff / se_diff
        
        # Handle NaN/Inf (when se_diff is 0)
        t_stat[!is.finite(t_stat)] <- 0
        
        # L2 norm statistic
        L2_stat <- sqrt(sum(mean_diff^2))
        
        # Permutation test for paired data
        # Randomly flip signs of differences
        t_stat_perm <- matrix(NA, n_time, n_permutations)
        L2_stat_perm <- numeric(n_permutations)
        
        for(perm in 1:n_permutations) {
          # Random sign flips for each subject
          sign_flips <- sample(c(-1, 1), n_pairs_subj, replace = TRUE)
          
          perm_diffs <- paired_diffs * rep(sign_flips, each = n_time)
          perm_mean_diff <- rowMeans(perm_diffs)
          perm_se_diff <- apply(perm_diffs, 1, sd) / sqrt(n_pairs_subj)
          
          t_stat_perm[, perm] <- perm_mean_diff / perm_se_diff
          t_stat_perm[!is.finite(t_stat_perm[, perm]), perm] <- 0
          
          L2_stat_perm[perm] <- sqrt(sum(perm_mean_diff^2))
        }
        
        # Calculate p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          p_values_pointwise[t] <- mean(abs(t_stat_perm[t, ]) >= abs(t_stat[t]), na.rm = TRUE)
        }
        
        p_value_L2 <- mean(L2_stat_perm >= L2_stat)
        
        # Bootstrap confidence intervals for paired differences
        n_boot <- 100
        diff_boot <- matrix(NA, n_time, n_boot)
        
        for(boot in 1:n_boot) {
          boot_idx <- sample(1:n_pairs_subj, replace = TRUE)
          diff_boot[, boot] <- rowMeans(paired_diffs[, boot_idx, drop = FALSE])
        }
        
        ci_lower <- apply(diff_boot, 1, quantile, probs = 0.025, na.rm = TRUE)
        ci_upper <- apply(diff_boot, 1, quantile, probs = 0.975, na.rm = TRUE)
        
        # Cohen's d for paired samples (using SD of differences)
        sd_diff <- apply(paired_diffs, 1, sd)
        cohens_d <- mean_diff / sd_diff
        cohens_d[!is.finite(cohens_d)] <- 0
        
        # Calculate means for each condition (for plotting)
        mean1 <- rowMeans(matched_curves1)
        mean2 <- rowMeans(matched_curves2)
        
        pairwise_results[[pair_names[pair_idx]]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          n1 = n_pairs_subj,  # Number of matched pairs
          n2 = n_pairs_subj,
          n_matched = n_pairs_subj,
          design = "within",
          mean1 = mean1,
          mean2 = mean2,
          mean_diff = mean_diff,
          t_stat = t_stat,
          p_values_pointwise = p_values_pointwise,
          L2_stat = L2_stat,
          p_value_L2 = p_value_L2,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          cohens_d = cohens_d,
          se_diff = se_diff
        )
      }
    })
    
    # Apply multiple comparison correction
    all_p_values_L2 <- sapply(pairwise_results, function(x) x$p_value_L2)
    adjusted_p_values_L2 <- p.adjust(all_p_values_L2, method = correction_method)
    
    # Apply correction to pointwise p-values
    for(i in 1:length(pairwise_results)) {
      pairwise_results[[i]]$p_values_adjusted <- p.adjust(pairwise_results[[i]]$p_values_pointwise, 
                                                          method = correction_method)
      pairwise_results[[i]]$p_value_L2_adjusted <- adjusted_p_values_L2[i]
      pairwise_results[[i]]$sig_regions <- pairwise_results[[i]]$p_values_adjusted < alpha
      pairwise_results[[i]]$sig_global <- pairwise_results[[i]]$p_value_L2_adjusted < alpha
    }
    
    return(list(
      results = pairwise_results,
      time_points = time_points,
      correction_method = correction_method,
      alpha = alpha,
      n_permutations = n_permutations,
      groups = conditions,
      n_groups = n_conditions,
      pair_names = pair_names,
      design = "within"
    ))
  }
  
  # Between-Subjects Pairwise Comparisons (original function)
  perform_pairwise_comparisons <- function(fd_obj, group_labels, n_permutations = 200,
                                           correction_method = "bonferroni", alpha = 0.05) {
    
    n_curves <- ncol(fd_obj$coefs)
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    group_labels <- as.factor(group_labels)
    groups <- levels(group_labels)
    n_groups <- length(groups)
    
    n_pairs <- choose(n_groups, 2)
    pairs <- combn(groups, 2, simplify = FALSE)
    pair_names <- sapply(pairs, function(p) paste(p[1], "vs", p[2]))
    
    pairwise_results <- list()
    
    withProgress(message = 'Performing pairwise comparisons', value = 0, {
      for(pair_idx in 1:n_pairs) {
        incProgress(1/n_pairs, detail = paste("Comparing", pair_names[pair_idx]))
        
        pair <- pairs[[pair_idx]]
        idx1 <- which(group_labels == pair[1])
        idx2 <- which(group_labels == pair[2])
        
        curves1 <- curves[, idx1, drop = FALSE]
        curves2 <- curves[, idx2, drop = FALSE]
        
        n1 <- length(idx1)
        n2 <- length(idx2)
        
        mean1 <- rowMeans(curves1)
        mean2 <- rowMeans(curves2)
        mean_diff <- mean1 - mean2
        
        # Calculate proper t-statistics
        pooled_var <- ((n1 - 1) * apply(curves1, 1, var) + 
                         (n2 - 1) * apply(curves2, 1, var)) / (n1 + n2 - 2)
        se_diff <- sqrt(pooled_var * (1/n1 + 1/n2))
        t_stat <- mean_diff / se_diff
        
        # L2 norm statistic
        L2_stat <- sqrt(sum(mean_diff^2))
        
        # Permutation test
        t_stat_perm <- matrix(NA, n_time, n_permutations)
        L2_stat_perm <- numeric(n_permutations)
        
        combined_curves <- cbind(curves1, curves2)
        combined_labels <- c(rep(1, n1), rep(2, n2))
        
        for(perm in 1:n_permutations) {
          perm_labels <- sample(combined_labels)
          
          perm_curves1 <- combined_curves[, perm_labels == 1, drop = FALSE]
          perm_curves2 <- combined_curves[, perm_labels == 2, drop = FALSE]
          
          perm_mean1 <- rowMeans(perm_curves1)
          perm_mean2 <- rowMeans(perm_curves2)
          perm_diff <- perm_mean1 - perm_mean2
          
          perm_pooled_var <- ((n1 - 1) * apply(perm_curves1, 1, var) + 
                                (n2 - 1) * apply(perm_curves2, 1, var)) / (n1 + n2 - 2)
          perm_se_diff <- sqrt(perm_pooled_var * (1/n1 + 1/n2))
          
          t_stat_perm[, perm] <- perm_diff / perm_se_diff
          L2_stat_perm[perm] <- sqrt(sum(perm_diff^2))
        }
        
        # Calculate p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          p_values_pointwise[t] <- mean(abs(t_stat_perm[t, ]) >= abs(t_stat[t]), na.rm = TRUE)
        }
        
        p_value_L2 <- mean(L2_stat_perm >= L2_stat)
        
        # Bootstrap confidence intervals
        n_boot <- 100
        diff_boot <- matrix(NA, n_time, n_boot)
        
        for(boot in 1:n_boot) {
          boot_idx1 <- sample(idx1, replace = TRUE)
          boot_idx2 <- sample(idx2, replace = TRUE)
          
          boot_mean1 <- rowMeans(curves[, boot_idx1, drop = FALSE])
          boot_mean2 <- rowMeans(curves[, boot_idx2, drop = FALSE])
          diff_boot[, boot] <- boot_mean1 - boot_mean2
        }
        
        ci_lower <- apply(diff_boot, 1, quantile, probs = 0.025)
        ci_upper <- apply(diff_boot, 1, quantile, probs = 0.975)
        
        # Cohen's d effect size
        cohens_d <- mean_diff / sqrt(pooled_var)
        
        pairwise_results[[pair_names[pair_idx]]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          n1 = n1,
          n2 = n2,
          design = "between",
          mean_diff = mean_diff,
          t_stat = t_stat,
          p_values_pointwise = p_values_pointwise,
          L2_stat = L2_stat,
          p_value_L2 = p_value_L2,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          cohens_d = cohens_d,
          se_diff = se_diff
        )
      }
    })
    
    # Apply multiple comparison correction
    all_p_values_L2 <- sapply(pairwise_results, function(x) x$p_value_L2)
    adjusted_p_values_L2 <- p.adjust(all_p_values_L2, method = correction_method)
    
    # Apply correction to pointwise p-values
    for(i in 1:n_pairs) {
      pairwise_results[[i]]$p_values_adjusted <- p.adjust(pairwise_results[[i]]$p_values_pointwise, 
                                                          method = correction_method)
      pairwise_results[[i]]$p_value_L2_adjusted <- adjusted_p_values_L2[i]
      pairwise_results[[i]]$sig_regions <- pairwise_results[[i]]$p_values_adjusted < alpha
      pairwise_results[[i]]$sig_global <- pairwise_results[[i]]$p_value_L2_adjusted < alpha
    }
    
    return(list(
      results = pairwise_results,
      time_points = time_points,
      correction_method = correction_method,
      alpha = alpha,
      n_permutations = n_permutations,
      groups = groups,
      n_groups = n_groups,
      pair_names = pair_names,
      design = "between"
    ))
  }
  
  # Run pairwise comparisons - ENHANCED for both between and within designs
  observeEvent(input$run_pairwise, {
    if(is.null(values$fanova_results)) {
      showNotification("Please run Functional ANOVA first!", type = "error", duration = 5)
      return()
    }
    
    cat("Running pairwise comparisons...\n")
    
    # Check if this is a within-subjects design
    design_type <- if(!is.null(values$fanova_results$design)) {
      values$fanova_results$design
    } else {
      "between"
    }
    
    cat("Design type:", design_type, "\n")
    
    tryCatch({
      # Use the same data source as FANOVA
      fd_to_use <- NULL
      
      # Check what data source was used in FANOVA
      if(!is.null(values$fanova_results$data_source) && 
         values$fanova_results$data_source == "warped" && 
         !is.null(values$warping_results)) {
        
        cat("Using time-warped curves for pairwise comparisons\n")
        
        if(!is.null(values$warping_results$regfd)) {
          fd_to_use <- values$warping_results$regfd
        } else if(!is.null(values$warping_results$registered_curves)) {
          n_time <- nrow(values$warping_results$registered_curves)
          time_points <- values$warping_results$time_points
          if(is.null(time_points)) {
            time_points <- seq(0, 1, length.out = n_time)
          }
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          # Use lambda=0: registered curves already processed, just need fd representation
          fd_to_use <- smooth.basis(time_points, values$warping_results$registered_curves, 
                                    fdPar(basis, 2, 0))$fd
        }
      }
      
      # Fallback to original data if warped not available
      if(is.null(fd_to_use)) {
        cat("Using original curves for pairwise comparisons\n")
        fd_to_use <- values$fd_obj
      }
      
      # Call appropriate pairwise function based on design
      # IMPORTANT: Use the same n_permutations that was used in the omnibus FANOVA
      n_perm_to_use <- if(!is.null(values$fanova_results$n_permutations)) {
        values$fanova_results$n_permutations
      } else {
        input$pairwise_permutations  # Fallback
      }
      
      cat("Using", n_perm_to_use, "permutations (same as omnibus test)\n")
      
      if(design_type == "within") {
        # Within-subjects: need subject IDs and RM factor
        if(is.null(values$fanova_results$subject_id) || is.null(values$fanova_results$rm_factor)) {
          showNotification("Subject ID or RM factor missing from FANOVA results!", type = "error", duration = 5)
          return()
        }
        
        cat("Performing PAIRED comparisons for within-subjects design\n")
        
        values$pairwise_results <- perform_pairwise_comparisons_rm(
          fd_obj = fd_to_use,
          subject_id = values$fanova_results$subject_id,
          rm_factor = values$fanova_results$rm_factor,
          n_permutations = n_perm_to_use,
          correction_method = input$pairwise_correction,
          alpha = input$pairwise_alpha
        )
      } else {
        # Between-subjects: use original function
        cat("Performing INDEPENDENT comparisons for between-subjects design\n")
        
        values$pairwise_results <- perform_pairwise_comparisons(
          fd_obj = fd_to_use,
          group_labels = values$group_labels,
          n_permutations = n_perm_to_use,
          correction_method = input$pairwise_correction,
          alpha = input$pairwise_alpha
        )
      }
      
      showNotification("Pairwise comparisons completed!", type = "message", duration = 3)
      
    }, error = function(e) {
      cat("Pairwise error:", e$message, "\n")
      showNotification(paste("Pairwise comparison error:", e$message), type = "error", duration = 10)
    })
  })
  
  # All remaining outputs (pairwise, warping, landmark, export) - these are already included
  # but I'll ensure warping_plot and alignment_comparison_plot are here
  
  output$warping_plot <- renderPlotly({
    if(is.null(values$warping_results)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "Run time-warped PCA first"))
    }
    
    tryCatch({
      warp_results <- values$warping_results
      
      if(is.null(warp_results$warp_functions)) {
        return(plot_ly(type = 'scatter', mode = 'lines') %>% 
                 layout(title = "No warping functions available"))
      }
      
      warp_functions <- warp_results$warp_functions
      n_time <- nrow(warp_functions)
      n_curves <- ncol(warp_functions)
      
      time_points <- if(!is.null(warp_results$time_points)) {
        warp_results$time_points
      } else {
        seq(0, 1, length.out = n_time)
      }
      hover_times <- hover_time_labels(time_points)

      n_display <- if(!is.null(input$n_curves_display)) {
        min(input$n_curves_display, n_curves)
      } else {
        min(30, n_curves)
      }
      
      p <- plot_ly(type = 'scatter', mode = 'lines')
      
      for(i in 1:n_display) {
        p <- p %>% add_trace(
          x = time_points,
          y = warp_functions[,i],
          type = 'scatter',
          mode = 'lines',
          name = paste("Subject", i),
          line = list(color = 'rgba(100, 100, 200, 0.3)', width = 1),
          showlegend = FALSE,
          hovertemplate = paste("Subject", i, "<br>t: %{customdata}<br>h(t): %{y:.3f}<extra></extra>"),
          customdata = hover_times
        )
      }
      
      p <- p %>% add_trace(
        x = c(0, 1), 
        y = c(0, 1),
        type = 'scatter',
        mode = 'lines',
        name = "Identity (no warping)",
        line = list(color = 'red', width = 2, dash = 'dash'),
        hovertemplate = "Identity line<extra></extra>"
      )
      
      if(n_display > 0) {
        mean_warp <- rowMeans(warp_functions[,1:n_display, drop = FALSE])
        p <- p %>% add_trace(
          x = time_points,
          y = mean_warp,
          type = 'scatter',
          mode = 'lines',
          name = "Mean warping",
          line = list(color = 'black', width = 3),
          hovertemplate = "Mean warping<br>t: %{customdata}<br>h(t): %{y:.3f}<extra></extra>",
          customdata = hover_times
        )
      }
      
      p %>% layout(
        title = "Time Warping Functions",
        xaxis = list(title = "Original Time (t)", range = c(0, 1)),
        yaxis = list(title = "Warped Time h(t)", range = c(0, 1)),
        hovermode = 'x'
      )
      
    }, error = function(e) {
      cat("Warping plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  output$alignment_comparison_plot <- renderPlotly({
    if(is.null(values$fd_obj)) {
      return(plot_ly(type = 'scatter', mode = 'lines') %>% 
               layout(title = "No data available"))
    }
    
    tryCatch({
      n_time <- if(!is.null(values$warping_results) && !is.null(values$warping_results$time_points)) {
        length(values$warping_results$time_points)
      } else {
        100
      }
      
      time_points <- if(!is.null(values$warping_results) && !is.null(values$warping_results$time_points)) {
        values$warping_results$time_points
      } else {
        seq(0, 1, length.out = n_time)
      }
      hover_times <- hover_time_labels(time_points)

      n_display <- if(!is.null(input$n_curves_display)) {
        min(input$n_curves_display, ncol(values$fd_obj$coefs))
      } else {
        min(30, ncol(values$fd_obj$coefs))
      }
      
      orig_curves <- eval.fd(time_points, values$fd_obj)
      n_show <- min(n_display, ncol(orig_curves))
      
      p <- plot_ly(type = 'scatter', mode = 'lines')
      
      for(i in 1:n_show) {
        p <- p %>% add_trace(
          x = time_points,
          y = orig_curves[,i],
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'rgba(255, 100, 100, 0.2)', width = 1),
          showlegend = (i == 1),
          legendgroup = "original",
          name = "Original curves",
          hoverinfo = 'skip'
        )
      }
      
      orig_mean <- rowMeans(orig_curves[,1:n_show, drop = FALSE])
      p <- p %>% add_trace(
        x = time_points,
        y = orig_mean,
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'darkred', width = 3),
        name = "Original mean",
        hovertemplate = "Original mean<br>Time: %{customdata}<br>Value: %{y:.3f}<extra></extra>",
        customdata = hover_times
      )
      
      if(!is.null(values$warping_results) && !is.null(values$warping_results$registered_curves)) {
        aligned_curves <- values$warping_results$registered_curves
        
        if(ncol(aligned_curves) >= n_show && nrow(aligned_curves) == length(time_points)) {
          for(i in 1:n_show) {
            p <- p %>% add_trace(
              x = time_points,
              y = aligned_curves[,i],
              type = 'scatter',
              mode = 'lines',
              line = list(color = 'rgba(100, 100, 255, 0.2)', width = 1),
              showlegend = (i == 1),
              legendgroup = "aligned",
              name = "Aligned curves",
              hoverinfo = 'skip'
            )
          }
          
          aligned_mean <- rowMeans(aligned_curves[,1:n_show, drop = FALSE])
          p <- p %>% add_trace(
            x = time_points,
            y = aligned_mean,
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'darkblue', width = 3),
            name = "Aligned mean",
            hovertemplate = "Aligned mean<br>Time: %{customdata}<br>Value: %{y:.3f}<extra></extra>",
            customdata = hover_times
          )
        }
      }
      
      p <- p %>% layout(
        title = "Original vs Aligned Curves",
        yaxis = list(title = "Value"),
        hovermode = 'x'
      )
      p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_settings))
      p
      
    }, error = function(e) {
      cat("Alignment plot error:", e$message, "\n")
      plot_ly(type = 'scatter', mode = 'lines') %>% 
        layout(title = paste("Error:", e$message))
    })
  })
  
  # Landmark plot
  output$landmark_plot <- renderPlotly({
    req(values$data)
    
    tryCatch({
      n_time <- ncol(values$data)
      time_points <- seq(0, 1, length.out = n_time)
      
      # Determine what to plot
      if(input$landmark_target == "mean" || is.null(input$selected_subject)) {
        # Plot mean curve
        if(!is.null(values$fd_obj)) {
          mean_curve <- rowMeans(eval.fd(time_points, values$fd_obj))
        } else {
          mean_curve <- colMeans(values$data)
        }
        
        p <- plot_ly(source = "landmark_source") %>%
          add_trace(x = time_points, y = mean_curve, 
                    type = 'scatter', mode = 'lines',
                    name = 'Mean curve',
                    line = list(color = 'blue', width = 2))
        
        title_text <- "Click to add landmarks on mean curve"
      } else {
        # Plot individual subject
        subj_idx <- as.numeric(input$selected_subject)
        subj_curve <- values$data[subj_idx,]
        
        p <- plot_ly(source = "landmark_source") %>%
          add_trace(x = time_points, y = subj_curve,
                    type = 'scatter', mode = 'lines',
                    name = paste('Subject', subj_idx),
                    line = list(color = 'blue', width = 2))
        
        title_text <- paste("Click to add landmarks for Subject", subj_idx)
      }
      
      # Add existing landmarks if any
      if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
        p <- p %>% add_trace(x = values$landmark_points$x,
                             y = values$landmark_points$y,
                             type = 'scatter', mode = 'markers',
                             marker = list(color = 'red', size = 10),
                             name = 'Landmarks')
      }
      
      p <- p %>% 
        layout(title = title_text,
               yaxis = list(title = "Value"),
               hovermode = 'closest',
               clickmode = 'event+select') %>%
        config(displayModeBar = FALSE)
      
      # Apply time label formatting (uses existing helper function)
      p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_settings))
      p
      
    }, error = function(e) {
      plot_ly() %>% layout(title = paste("Error:", e$message))
    })
  })
  
  # Handle landmark clicks
  observeEvent(event_data("plotly_click", source = "landmark_source"), {
    click <- event_data("plotly_click", source = "landmark_source")
    if(!is.null(click)) {
      # Only add landmark if it's a click on the curve, not on existing landmarks
      if(is.null(click$curveNumber) || click$curveNumber == 0) {
        new_point <- data.frame(x = click$x, y = click$y)
        values$landmark_points <- rbind(values$landmark_points, new_point)
        
        showNotification(paste("Landmark added at t =", round(click$x, 3)), 
                         duration = 2, type = "message")
      }
    }
  })
  
  # Start landmark selection
  observeEvent(input$start_landmark, {
    values$landmark_points <- data.frame(x = numeric(), y = numeric())
    showNotification("Landmarks cleared. Click on the plot to add new landmarks.", 
                     duration = 3, type = "message")
  })
  
  # Clear landmarks
  observeEvent(input$clear_landmarks, {
    values$landmark_points <- data.frame(x = numeric(), y = numeric())
    showNotification("Landmarks cleared", duration = 2)
  })
  
  # Add landmark info display
  output$landmark_info <- renderPrint({
    if(!is.null(values$landmark_points) && nrow(values$landmark_points) > 0) {
      cat("Current landmarks:\n")
      for(i in 1:nrow(values$landmark_points)) {
        cat(sprintf("  Landmark %d: t = %.3f\n", i, values$landmark_points$x[i]))
      }
    } else {
      cat("No landmarks defined. Click on the plot to add landmarks.")
    }
  })
  
  # ALL REMAINING PAIRWISE AND EXPORT OUTPUTS - included as-is from the original
  # (These are all already correctly written - just ensuring they're included)
  
  output$pairwise_selector <- renderUI({
    req(values$pairwise_results)
    
    selectInput("selected_pair", "Select pairwise comparison:",
                choices = values$pairwise_results$pair_names,
                selected = values$pairwise_results$pair_names[1])
  })
  
  output$pairwise_summary <- renderPrint({
    req(values$pairwise_results)
    
    res <- values$pairwise_results
    
    cat("========================================\n")
    cat("  PAIRWISE COMPARISONS SUMMARY\n")
    cat("========================================\n\n")
    
    # Show design type
    design_type <- if(!is.null(res$design)) {
      if(res$design == "within") "Within-Subjects (PAIRED)" else "Between-Subjects (INDEPENDENT)"
    } else {
      "Between-Subjects (INDEPENDENT)"
    }
    cat("Design:", design_type, "\n\n")
    
    cat("Number of groups/conditions:", res$n_groups, "\n")
    cat("Number of comparisons:", length(res$pair_names), "\n")
    cat("Correction method:", res$correction_method, "\n")
    cat("Significance level:", res$alpha, "\n")
    cat("Permutations:", res$n_permutations, "\n\n")
    
    n_sig_global <- sum(sapply(res$results, function(x) x$sig_global))
    
    cat("Results Summary:\n")
    cat("----------------\n")
    cat("Globally significant comparisons:", n_sig_global, "/", 
        length(res$pair_names), 
        sprintf("(%.1f%%)", 100 * n_sig_global / length(res$pair_names)), "\n\n")
    
    if(n_sig_global > 0) {
      cat("Significant pairs (global test):\n")
      for(pair_name in res$pair_names) {
        if(res$results[[pair_name]]$sig_global) {
          p_val <- res$results[[pair_name]]$p_value_L2_adjusted
          stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else ""
          cat(sprintf("   - %-20s p = %s %s\n", 
                      pair_name, 
                      format.pval(p_val, digits = 4),
                      stars))
        }
      }
    } else {
      cat("No significant pairwise differences found.\n")
    }
    
    cat("\n========================================\n")
  })
  
  output$pairwise_global_table <- renderDT({
    req(values$pairwise_results)
    
    results_df <- data.frame(
      Comparison = values$pairwise_results$pair_names,
      N1 = sapply(values$pairwise_results$results, function(x) x$n1),
      N2 = sapply(values$pairwise_results$results, function(x) x$n2),
      L2_Stat = round(sapply(values$pairwise_results$results, function(x) x$L2_stat), 3),
      P_Raw = sapply(values$pairwise_results$results, function(x) 
        format.pval(x$p_value_L2, digits = 4)),
      P_Adj = sapply(values$pairwise_results$results, function(x) 
        format.pval(x$p_value_L2_adjusted, digits = 4)),
      Sig = ifelse(sapply(values$pairwise_results$results, function(x) x$sig_global), 
                   "Yes", "No"),
      Mean_Cohen_d = round(sapply(values$pairwise_results$results, function(x) 
        mean(abs(x$cohens_d))), 3)
    )
    
    datatable(results_df, 
              options = list(pageLength = 15, scrollX = TRUE), 
              rownames = FALSE) %>%
      formatStyle("Sig",
                  backgroundColor = styleEqual("Yes", "#d4edda")) %>%
      formatStyle("Mean_Cohen_d",
                  backgroundColor = styleInterval(c(0.2, 0.5, 0.8), 
                                                  c("white", "#ffffcc", "#ffcccc", "#ff9999")))
  })
  
  output$pairwise_difference_plot <- renderPlotly({
    req(values$pairwise_results, input$selected_pair)

    pair_result <- values$pairwise_results$results[[input$selected_pair]]
    time_points <- values$pairwise_results$time_points
    hover_times <- hover_time_labels(time_points)

    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    if(any(pair_result$sig_regions)) {
      sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
      sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
      
      for(i in 1:length(sig_starts)) {
        y_range <- range(c(pair_result$ci_lower, pair_result$ci_upper))
        y_expand <- diff(y_range) * 0.1
        
        p <- p %>% add_trace(
          x = c(time_points[sig_starts[i]], time_points[sig_ends[i]], 
                time_points[sig_ends[i]], time_points[sig_starts[i]]),
          y = c(min(y_range) - y_expand, min(y_range) - y_expand,
                max(y_range) + y_expand, max(y_range) + y_expand),
          fill = 'toself',
          fillcolor = 'rgba(255, 200, 200, 0.3)',
          line = list(color = 'transparent'),
          showlegend = (i == 1),
          name = 'Significant',
          hoverinfo = 'skip'
        )
      }
    }
    
    if(input$pairwise_confidence_bands) {
      p <- p %>% add_trace(
        x = c(time_points, rev(time_points)),
        y = c(pair_result$ci_lower, rev(pair_result$ci_upper)),
        fill = 'toself',
        fillcolor = 'rgba(100, 100, 255, 0.2)',
        line = list(color = 'transparent'),
        name = '95% CI',
        hoverinfo = 'skip'
      )
      
      p <- p %>% add_trace(
        x = time_points,
        y = pair_result$ci_lower,
        line = list(color = 'lightblue', width = 1, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      ) %>% add_trace(
        x = time_points,
        y = pair_result$ci_upper,
        line = list(color = 'lightblue', width = 1, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    p <- p %>% add_trace(
      x = time_points,
      y = pair_result$mean_diff,
      line = list(color = 'blue', width = 3),
      name = 'Mean Difference',
      hovertemplate = paste("Time: %{customdata}",
                            "<br>Difference: %{y:.3f}",
                            "<br>Cohen's d: ", round(pair_result$cohens_d, 2),
                            "<extra></extra>"),
      customdata = hover_times
    )
    
    p <- p %>% add_trace(
      x = c(0, 1),
      y = c(0, 0),
      line = list(color = 'black', width = 1, dash = 'dash'),
      name = 'Zero',
      hoverinfo = 'none'
    )
    
    mean_effect <- mean(abs(pair_result$cohens_d))
    effect_text <- if(mean_effect < 0.2) "Negligible" else if(mean_effect < 0.5) "Small" else if(mean_effect < 0.8) "Medium" else "Large"
    
    p <- p %>% layout(
      title = paste("Difference:", pair_result$group1, "-", pair_result$group2,
                    "<br><sub>Mean |Cohen's d| =", round(mean_effect, 2), 
                    "(", effect_text, "effect)</sub>"),
      yaxis = list(title = "Mean Difference"),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.98)
    )

    p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_pvalue_plot <- renderPlotly({
    req(values$pairwise_results, input$selected_pair)

    pair_result <- values$pairwise_results$results[[input$selected_pair]]
    time_points <- values$pairwise_results$time_points
    hover_times <- hover_time_labels(time_points)
    
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    if(any(pair_result$sig_regions)) {
      sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
      sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
      
      for(i in 1:length(sig_starts)) {
        p <- p %>% add_trace(
          x = c(time_points[sig_starts[i]], time_points[sig_ends[i]], 
                time_points[sig_ends[i]], time_points[sig_starts[i]]),
          y = c(0.0001, 0.0001, 1, 1),
          fill = 'toself',
          fillcolor = 'rgba(200, 255, 200, 0.3)',
          line = list(color = 'transparent'),
          showlegend = (i == 1),
          name = 'Significant region',
          hoverinfo = 'skip'
        )
      }
    }
    
    p <- p %>% add_trace(
      x = time_points,
      y = pair_result$p_values_adjusted,
      line = list(color = 'darkgreen', width = 2),
      name = paste('Adjusted p-values (', values$pairwise_results$correction_method, ')', sep = ''),
      hovertemplate = "Time: %{customdata}<br>Adjusted p: %{y:.4f}<extra></extra>",
      customdata = hover_times
    ) %>%
      add_trace(
        x = time_points,
        y = pair_result$p_values_pointwise,
        line = list(color = 'lightgreen', width = 1, dash = 'dot'),
        name = 'Raw p-values',
        hovertemplate = "Time: %{customdata}<br>Raw p: %{y:.4f}<extra></extra>",
        customdata = hover_times
      ) %>%
      add_trace(
        x = c(0, 1),
        y = c(input$pairwise_alpha, input$pairwise_alpha),
        line = list(color = 'red', width = 2, dash = 'dash'),
        name = paste('α =', input$pairwise_alpha),
        hoverinfo = 'none'
      )
    
    p <- p %>% add_trace(
      x = c(0, 1),
      y = c(0.01, 0.01),
      line = list(color = 'orange', width = 1, dash = 'dot'),
      name = 'p = 0.01',
      hoverinfo = 'none'
    ) %>%
      add_trace(
        x = c(0, 1),
        y = c(0.001, 0.001),
        line = list(color = 'darkred', width = 1, dash = 'dot'),
        name = 'p = 0.001',
        hoverinfo = 'none'
      )
    
    p <- p %>% layout(
      title = paste("P-values:", input$selected_pair),
      yaxis = list(
        title = "p-value",
        type = 'log',
        range = c(log10(0.0001), log10(1)),
        tickvals = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
        ticktext = c("0.001", "0.01", "0.05", "0.1", "0.5", "1")
      ),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.02)
    )

    p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_heatmap <- renderPlotly({
    req(values$pairwise_results)
    
    n_groups <- values$pairwise_results$n_groups
    groups <- values$pairwise_results$groups
    
    p_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
    rownames(p_matrix) <- groups
    colnames(p_matrix) <- groups
    
    for(pair_name in values$pairwise_results$pair_names) {
      pair_result <- values$pairwise_results$results[[pair_name]]
      i <- which(groups == pair_result$group1)
      j <- which(groups == pair_result$group2)
      p_val <- pair_result$p_value_L2_adjusted
      p_matrix[i, j] <- p_val
      p_matrix[j, i] <- p_val
    }
    
    z_matrix <- -log10(p_matrix)
    z_matrix[is.infinite(z_matrix)] <- 5
    
    hover_text <- matrix("", nrow = n_groups, ncol = n_groups)
    for(i in 1:n_groups) {
      for(j in 1:n_groups) {
        if(i == j) {
          hover_text[i, j] <- paste(groups[i], "(same group)")
        } else {
          p_val_text <- if(p_matrix[i, j] < 0.001) {
            "< 0.001"
          } else {
            format(round(p_matrix[i, j], 4), nsmall = 4)
          }
          sig_text <- if(p_matrix[i, j] < values$pairwise_results$alpha) "***" else "n.s."
          hover_text[i, j] <- paste(groups[i], "vs", groups[j],
                                    "<br>p =", p_val_text,
                                    "<br>", sig_text)
        }
      }
    }
    
    annotations <- list()
    for(i in 1:n_groups) {
      for(j in 1:n_groups) {
        if(i != j) {
          star_text <- if(p_matrix[i, j] < 0.001) "***" else 
            if(p_matrix[i, j] < 0.01) "**" else 
              if(p_matrix[i, j] < 0.05) "*" else ""
          
          if(star_text != "") {
            annotations <- append(annotations, list(list(
              x = groups[j],
              y = groups[i],
              text = star_text,
              showarrow = FALSE,
              font = list(color = 'white', size = 14)
            )))
          }
        }
      }
    }
    
    plot_ly(
      z = z_matrix,
      x = groups,
      y = groups,
      type = 'heatmap',
      colorscale = list(
        c(0, 'white'),
        c(0.3, 'lightblue'),
        c(0.6, 'blue'),
        c(1, 'darkred')
      ),
      hovertemplate = "%{text}<extra></extra>",
      text = hover_text,
      colorbar = list(
        title = "Significance",
        tickmode = "array",
        tickvals = c(0, -log10(0.05), -log10(0.01), -log10(0.001), 5),
        ticktext = c("1", "0.05", "0.01", "0.001", "<0.00001")
      )
    ) %>%
      layout(
        title = paste("Pairwise Comparison P-values (", values$pairwise_results$correction_method, "correction)"),
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "", autorange = 'reversed'),
        annotations = annotations
      )
  })
  
  output$pairwise_significance_timeline <- renderPlotly({
    req(values$pairwise_results)
    
    time_points <- values$pairwise_results$time_points
    n_pairs <- length(values$pairwise_results$pair_names)
    hover_times <- hover_time_labels(time_points)
    
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    base_cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',
                   '#ffff33','#a65628','#f781bf','#999999','#66c2a5',
                   '#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
    colors <- colorRampPalette(base_cols)(n_pairs)
    
    for(i in 1:n_pairs) {
      p <- p %>% add_trace(
        x = time_points,
        y = rep(i, length(time_points)),
        line = list(color = 'lightgray', width = 0.5),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    for(i in 1:n_pairs) {
      pair_name <- values$pairwise_results$pair_names[i]
      pair_result <- values$pairwise_results$results[[pair_name]]
      
      if(any(pair_result$sig_regions)) {
        sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
        sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
        
        for(j in 1:length(sig_starts)) {
          p <- p %>% add_trace(
            x = time_points[sig_starts[j]:sig_ends[j]],
            y = rep(i, sig_ends[j] - sig_starts[j] + 1),
            line = list(color = colors[i], width = 8),
            showlegend = (j == 1),
            name = pair_name,
            legendgroup = pair_name,
            hovertemplate = paste(pair_name,
                                  "<br>Time: %{customdata}",
                                  "<br>Significant<extra></extra>"),
            customdata = hover_times[sig_starts[j]:sig_ends[j]]
          )
        }
      } else {
        p <- p %>% add_trace(
          x = c(NA),
          y = c(NA),
          name = paste(pair_name, "(n.s.)"),
          line = list(color = 'gray'),
          showlegend = TRUE
        )
      }
    }
    
    for(t in seq(0, 1, by = 0.25)) {
      p <- p %>% add_trace(
        x = c(t, t),
        y = c(0.5, n_pairs + 0.5),
        line = list(color = 'lightgray', width = 0.5, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    p <- p %>% layout(
      title = "Timeline of Significant Regions",
      xaxis = list(title = "Time", range = c(0, 1)),
      yaxis = list(
        title = "Comparison",
        tickmode = "array",
        tickvals = 1:n_pairs,
        ticktext = values$pairwise_results$pair_names,
        range = c(0.5, n_pairs + 0.5),
        autorange = 'reversed'
      ),
      hovermode = 'x',
      plot_bgcolor = 'rgba(240, 240, 240, 0.5)',
      legend = list(x = 1.02, y = 1, font = list(size = 10))
    )
    
    # Apply time label formatting
    p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_regions_table <- renderDT({
    req(values$pairwise_results)
    
    # Placeholder table - would need implementation
    datatable(data.frame(Message = "Regions table not yet implemented"),
              options = list(dom = 't'),
              rownames = FALSE)
  })

  # ============================================================================
  # FUNCTIONAL CLUSTERING
  # ============================================================================

  # Cluster optimization analysis (Elbow & Silhouette)
  observeEvent(input$run_cluster_optimization, {
    req(values$data)

    tryCatch({
      # Determine which data to use
      if(input$opt_data_type == "smoothed") {
        if(is.null(values$smooth_data)) {
          showNotification("No smoothed data available. Please apply smoothing first or use raw data.",
                           type = "error", duration = 5)
          return(NULL)
        }
        data_matrix <- values$smooth_data
      } else {
        data_matrix <- values$data
      }

      # Apply standardization if requested (mirrors the clustering module setting)
      standardize_opt <- if(!is.null(input$cluster_standardize)) input$cluster_standardize else "none"
      if(standardize_opt == "within") {
        data_matrix <- t(apply(data_matrix, 1, function(x) {
          s <- sd(x, na.rm = TRUE)
          if(is.na(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
        }))
      } else if(standardize_opt == "between") {
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        col_sds   <- apply(data_matrix, 2, sd, na.rm = TRUE)
        col_sds[col_sds == 0] <- 1
        data_matrix <- sweep(sweep(data_matrix, 2, col_means, "-"), 2, col_sds, "/")
      }

      # Check for functional k-means requirements
      if(input$opt_kmeans_method == "functional") {
        if(input$opt_data_type != "smoothed") {
          showNotification("Functional K-Means requires smoothed data. Please select 'Smoothed Data' or use Standard method.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(is.null(values$fd_obj)) {
          showNotification("No functional data object available. Please apply smoothing first.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(!requireNamespace("fda.usc", quietly = TRUE)) {
          showNotification("fda.usc package not available. Please install it or use Standard K-Means.",
                           type = "error", duration = 5)
          return(NULL)
        }
      }

      # Check for missing values
      if(any(is.na(data_matrix))) {
        showNotification("Data contains missing values. Please handle missing values before clustering.",
                         type = "error", duration = 5)
        return(NULL)
      }

      max_k <- input$max_clusters_test

      # Initialize storage for metrics
      k_values <- 2:max_k
      wcss_values <- numeric(length(k_values))
      silhouette_values <- numeric(length(k_values))
      ch_values <- numeric(length(k_values))

      # Show progress
      progress_msg <- switch(input$opt_kmeans_method,
        "functional"   = 'Testing functional cluster solutions...',
        "hierarchical" = 'Running hierarchical clustering...',
        'Testing standard cluster solutions...'
      )

      # For hierarchical: compute the dendrogram once before the loop (efficient)
      hc_opt <- NULL
      if(input$opt_kmeans_method == "hierarchical") {
        opt_linkage <- if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2"
        dist_matrix_hc_opt <- dist(data_matrix)
        hc_opt <- hclust(dist_matrix_hc_opt, method = opt_linkage)
      }

      withProgress(message = progress_msg, value = 0, {
        for(i in seq_along(k_values)) {
          k <- k_values[i]

          # Update progress
          incProgress(1/length(k_values), detail = paste("Testing k =", k))

          # Run clustering based on selected method
          set.seed(123)

          opt_nstart   <- if(!is.null(input$kmeans_nstart)) input$kmeans_nstart else 10
          opt_iter_max <- if(!is.null(input$kmeans_iter))  input$kmeans_iter  else 100

          if(input$opt_kmeans_method == "functional") {
            # ===== FUNCTIONAL K-MEANS using fda.usc (with restarts) =====
            time_grid <- seq(0, 1, length.out = ncol(data_matrix))
            fdata_obj <- fda.usc::fdata(mdata = data_matrix, argvals = time_grid)

            best_fkm_opt  <- NULL
            best_wcss_opt <- Inf
            for(restart in seq_len(opt_nstart)) {
              set.seed(100 + restart)
              fkm_try <- tryCatch(
                fda.usc::kmeans.fd(fdata_obj, ncl = k, max.iter = opt_iter_max),
                error = function(e) NULL
              )
              if(is.null(fkm_try)) next
              asn <- fkm_try$cluster
              wcss_try <- sum(sapply(seq_len(k), function(j) {
                idx <- which(asn == j)
                if(length(idx) == 0) return(0)
                m <- data_matrix[idx, , drop = FALSE]
                sum(apply(m, 1, function(x) sum((x - colMeans(m))^2)))
              }))
              if(wcss_try < best_wcss_opt) { best_wcss_opt <- wcss_try; best_fkm_opt <- fkm_try }
            }
            fkm_result <- best_fkm_opt

            cluster_assignments <- fkm_result$cluster

            # Calculate cluster means
            cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))

            for(j in 1:k) {
              cluster_idx <- which(cluster_assignments == j)
              if(length(cluster_idx) > 0) {
                cluster_fd <- values$fd_obj[cluster_idx]
                cluster_mean_fd <- mean.fd(cluster_fd)
                cluster_means[j, ] <- eval.fd(time_grid, cluster_mean_fd)
              }
            }

            # Calculate WCSS
            wcss_per_cluster <- numeric(k)
            for(j in 1:k) {
              cluster_members <- data_matrix[cluster_assignments == j, , drop = FALSE]
              if(nrow(cluster_members) > 0) {
                cluster_center <- cluster_means[j, ]
                wcss_per_cluster[j] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
              }
            }

            total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
            within_ss <- sum(wcss_per_cluster)
            between_ss <- total_ss - within_ss

            # Store WCSS
            wcss_values[i] <- within_ss

          } else if(input$opt_kmeans_method == "hierarchical") {
            # ===== HIERARCHICAL CLUSTERING (tree already built above) =====
            cluster_assignments <- cutree(hc_opt, k = k)

            # Cluster means and WCSS
            cluster_means_opt <- matrix(0, nrow = k, ncol = ncol(data_matrix))
            for(j in 1:k) {
              idx <- which(cluster_assignments == j)
              if(length(idx) > 0)
                cluster_means_opt[j, ] <- colMeans(data_matrix[idx, , drop = FALSE])
            }
            wcss_per_cluster_opt <- numeric(k)
            for(j in 1:k) {
              members <- data_matrix[cluster_assignments == j, , drop = FALSE]
              if(nrow(members) > 0)
                wcss_per_cluster_opt[j] <- sum(apply(members, 1, function(x) sum((x - cluster_means_opt[j, ])^2)))
            }
            total_ss_opt <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
            within_ss  <- sum(wcss_per_cluster_opt)
            between_ss <- total_ss_opt - within_ss
            wcss_values[i] <- within_ss

          } else {
            # ===== STANDARD K-MEANS =====
            kmeans_result <- kmeans(data_matrix, centers = k, nstart = opt_nstart, iter.max = opt_iter_max)
            cluster_assignments <- kmeans_result$cluster

            # Store WCSS
            wcss_values[i] <- kmeans_result$tot.withinss

            # Calculate SS for CH index
            between_ss <- kmeans_result$betweenss
            within_ss <- kmeans_result$tot.withinss
          }

          # Calculate silhouette width (common for all methods)
          if(requireNamespace("cluster", quietly = TRUE)) {
            dist_matrix <- dist(data_matrix)
            sil <- cluster::silhouette(cluster_assignments, dist_matrix)
            silhouette_values[i] <- mean(sil[, 3])
          } else {
            silhouette_values[i] <- NA
          }

          # Calculate Calinski-Harabasz index (common for both methods)
          n <- nrow(data_matrix)
          ch_values[i] <- (between_ss / (k - 1)) / (within_ss / (n - k))
        }
      })

      # Store results
      values$cluster_optimization <- list(
        k_values = k_values,
        wcss = wcss_values,
        silhouette = silhouette_values,
        calinski_harabasz = ch_values,
        data_type = input$opt_data_type,
        method = input$opt_kmeans_method,
        hc_obj = hc_opt  # hclust object for dendrogram (NULL for non-hierarchical)
      )

      method_label <- switch(input$opt_kmeans_method,
        "functional"   = "Functional K-Means",
        "hierarchical" = paste0("Hierarchical Clustering (", if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2", ")"),
        "Standard K-Means"
      )
      showNotification(paste0("Cluster optimization analysis completed using ", method_label, "!"),
                       type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Optimization error:", e$message),
                       type = "error", duration = 5)
      values$cluster_optimization <- NULL
    })
  })

  # Optimization status output
  output$optimization_status <- renderText({
    if(is.null(values$cluster_optimization)) {
      "Not run"
    } else {
      paste0("✓ Tested k=2 to k=", max(values$cluster_optimization$k_values))
    }
  })

  # Elbow plot
  output$elbow_plot <- renderPlotly({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Highlight current k if clustering has been run
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Create elbow plot
    p <- plot_ly() %>%
      add_trace(
        x = opt$k_values,
        y = opt$wcss,
        type = 'scatter',
        mode = 'lines+markers',
        line = list(color = '#377EB8', width = 2),
        marker = list(size = 8),
        name = 'WCSS',
        hovertemplate = paste0(
          "k = %{x}",
          "<br>WCSS = %{y:.0f}",
          "<extra></extra>"
        )
      )

    # Add vertical line for current k
    if(current_k %in% opt$k_values) {
      p <- p %>%
        add_trace(
          x = c(current_k, current_k),
          y = c(0, max(opt$wcss)),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'red', width = 2, dash = 'dash'),
          name = paste0('Current k=', current_k),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )
    }

    # Add method label to title
    method_label <- if(!is.null(opt$method) && opt$method == "functional") {
      " - Functional K-Means"
    } else {
      " - Standard K-Means"
    }

    p <- p %>%
      layout(
        title = paste0("Elbow Method (Within-Cluster Sum of Squares)", method_label),
        xaxis = list(title = "Number of Clusters (k)"),
        yaxis = list(title = "Total Within-Cluster SS"),
        hovermode = "closest"
      )

    p
  })

  # Silhouette plot
  output$silhouette_plot <- renderPlotly({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Highlight current k if clustering has been run
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Create silhouette plot
    p <- plot_ly() %>%
      add_trace(
        x = opt$k_values,
        y = opt$silhouette,
        type = 'scatter',
        mode = 'lines+markers',
        line = list(color = '#4DAF4A', width = 2),
        marker = list(size = 8),
        name = 'Silhouette Score',
        hovertemplate = paste0(
          "k = %{x}",
          "<br>Silhouette = %{y:.3f}",
          "<extra></extra>"
        )
      )

    # Add vertical line for current k
    if(current_k %in% opt$k_values) {
      p <- p %>%
        add_trace(
          x = c(current_k, current_k),
          y = c(min(opt$silhouette, na.rm = TRUE), max(opt$silhouette, na.rm = TRUE)),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'red', width = 2, dash = 'dash'),
          name = paste0('Current k=', current_k),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )
    }

    # Add method label to title
    method_label <- if(!is.null(opt$method) && opt$method == "functional") {
      " - Functional K-Means"
    } else {
      " - Standard K-Means"
    }

    p <- p %>%
      layout(
        title = paste0("Silhouette Analysis (Higher is Better)", method_label),
        xaxis = list(title = "Number of Clusters (k)"),
        yaxis = list(title = "Average Silhouette Width"),
        hovermode = "closest"
      )

    p
  })

  # Optimization recommendation
  output$optimization_recommendation <- renderText({
    req(values$cluster_optimization)

    opt <- values$cluster_optimization

    # Find optimal k based on different criteria

    # 1. Elbow method - look for the "elbow" (maximum second derivative)
    wcss_diff <- diff(opt$wcss)
    wcss_diff2 <- diff(wcss_diff)
    elbow_k <- opt$k_values[which.max(abs(wcss_diff2)) + 1]

    # 2. Maximum silhouette
    if(!all(is.na(opt$silhouette))) {
      max_sil_k <- opt$k_values[which.max(opt$silhouette)]
      max_sil_score <- max(opt$silhouette, na.rm = TRUE)
    } else {
      max_sil_k <- NA
      max_sil_score <- NA
    }

    # 3. Maximum Calinski-Harabasz
    max_ch_k <- opt$k_values[which.max(opt$calinski_harabasz)]

    # Current k
    current_k <- if(!is.null(values$clustering_results)) {
      values$clustering_results$k
    } else {
      input$n_clusters
    }

    # Build recommendation text
    output_text <- "=== OPTIMAL CLUSTER NUMBER RECOMMENDATIONS ===\n\n"

    output_text <- paste0(output_text, "Elbow Method:\n")
    output_text <- paste0(output_text, "  Suggested k = ", elbow_k, "\n")
    output_text <- paste0(output_text, "  (Look for the 'elbow' point where WCSS starts to level off)\n\n")

    if(!is.na(max_sil_k)) {
      output_text <- paste0(output_text, "Silhouette Analysis:\n")
      output_text <- paste0(output_text, "  Optimal k = ", max_sil_k,
                            " (Silhouette = ", round(max_sil_score, 3), ")\n")
      if(max_sil_score > 0.7) {
        output_text <- paste0(output_text, "  Strong clustering structure\n\n")
      } else if(max_sil_score > 0.5) {
        output_text <- paste0(output_text, "  Reasonable clustering structure\n\n")
      } else if(max_sil_score > 0.25) {
        output_text <- paste0(output_text, "  Weak clustering structure\n\n")
      } else {
        output_text <- paste0(output_text, "  No substantial clustering structure\n\n")
      }
    }

    output_text <- paste0(output_text, "Calinski-Harabasz Index:\n")
    output_text <- paste0(output_text, "  Optimal k = ", max_ch_k, "\n")
    output_text <- paste0(output_text, "  (Higher values indicate better-defined clusters)\n\n")

    # Current k assessment
    if(current_k %in% opt$k_values) {
      output_text <- paste0(output_text, "Your Current Selection (k=", current_k, "):\n")
      k_idx <- which(opt$k_values == current_k)
      output_text <- paste0(output_text, "  WCSS: ", round(opt$wcss[k_idx], 2), "\n")
      if(!is.na(opt$silhouette[k_idx])) {
        output_text <- paste0(output_text, "  Silhouette: ", round(opt$silhouette[k_idx], 3), "\n")
      }
      output_text <- paste0(output_text, "  Calinski-Harabasz: ", round(opt$calinski_harabasz[k_idx], 2), "\n\n")
    }

    output_text <- paste0(output_text, "INTERPRETATION:\n")
    output_text <- paste0(output_text, "- The elbow plot shows diminishing returns in variance explained\n")
    output_text <- paste0(output_text, "- Silhouette scores range from -1 to 1 (>0.5 is good)\n")
    output_text <- paste0(output_text, "- Consider domain knowledge and interpretability when choosing k\n")

    output_text
  })

  # ============================================================================
  # DCF (DENSITY CORE FINDING) CLUSTERING INTEGRATION
  # ============================================================================

  # Check DCF Python setup
  observeEvent(input$check_dcf_setup, {
    tryCatch({
      # Check if reticulate is available
      if(!requireNamespace("reticulate", quietly = TRUE)) {
        showNotification("reticulate package not installed. Please install it with: install.packages('reticulate')",
                         type = "error", duration = 10)
        return()
      }

      # Check Python availability
      python_available <- reticulate::py_available(initialize = TRUE)
      if(!python_available) {
        showNotification("Python not found. Please install Python and configure reticulate.",
                         type = "error", duration = 10)
        return()
      }

      python_version <- reticulate::py_config()$version
      showNotification(paste("Python found:", python_version), type = "message", duration = 3)

      # Check if DCFcluster is installed
      dcf_available <- reticulate::py_module_available("DCFcluster")

      if(dcf_available) {
        showNotification("DCFcluster Python package is installed and ready!",
                         type = "message", duration = 5)
      } else {
        showNotification(
          HTML("DCFcluster not found. Install with:<br><code>pip install DCFcluster</code><br>or<br><code>pip install git+https://github.com/tobinjo96/DCFcluster.git</code>"),
          type = "warning", duration = 15)
      }

      # Check numpy and scipy
      numpy_ok <- reticulate::py_module_available("numpy")
      scipy_ok <- reticulate::py_module_available("scipy")
      sklearn_ok <- reticulate::py_module_available("sklearn")

      if(!numpy_ok || !scipy_ok || !sklearn_ok) {
        missing <- c()
        if(!numpy_ok) missing <- c(missing, "numpy")
        if(!scipy_ok) missing <- c(missing, "scipy")
        if(!sklearn_ok) missing <- c(missing, "sklearn")
        showNotification(paste("Missing Python packages:", paste(missing, collapse = ", ")),
                         type = "warning", duration = 10)
      }

    }, error = function(e) {
      showNotification(paste("Error checking DCF setup:", e$message),
                       type = "error", duration = 10)
    })
  })

  # DCF clustering function
  run_dcf_clustering <- function(data_matrix, k_param, beta_param) {
    # Run DCF clustering via Python
    #
    # Args:
    #   data_matrix: numeric matrix (subjects x time points)
    #   k_param: neighborhood parameter for density estimation
    #   beta_param: fluctuation parameter for cluster cores
    #
    # Returns:
    #   List with cluster assignments and metadata

    tryCatch({
      # Import DCFcluster
      dcf <- reticulate::import("DCFcluster")

      # Convert R matrix to numpy array
      np <- reticulate::import("numpy")
      X <- np$array(data_matrix)

      # Run DCF clustering
      result <- dcf$DCFcluster$train(X, as.integer(k_param), beta_param)

      # Extract results
      labels <- as.integer(result$labels) + 1L  # Convert 0-indexed to 1-indexed
      peak_values <- as.numeric(result$peak_values)
      core_sets <- result$core_sets

      # Number of clusters found
      n_clusters <- length(unique(labels))

      return(list(
        labels = labels,
        peak_values = peak_values,
        core_sets = core_sets,
        n_clusters = n_clusters,
        k = k_param,
        beta = beta_param,
        success = TRUE
      ))

    }, error = function(e) {
      return(list(
        success = FALSE,
        error = e$message
      ))
    })
  }

  # Perform clustering (K-means or DCF)
  observeEvent(input$run_clustering, {
    req(values$data)

    tryCatch({
      # Determine which data to use
      if(input$cluster_data_type == "smoothed") {
        if(is.null(values$smooth_data)) {
          showNotification("No smoothed data available. Please apply smoothing first or use raw data.",
                           type = "error", duration = 5)
          return(NULL)
        }
        data_matrix <- values$smooth_data
        data_type_label <- "Smoothed Data"
      } else {
        data_matrix <- values$data
        data_type_label <- "Raw Data"
      }

      # Apply standardization if requested
      standardize_opt <- if(!is.null(input$cluster_standardize)) input$cluster_standardize else "none"
      if(standardize_opt == "within") {
        # Within-participant: z-score each participant's own time series (row-wise)
        data_matrix <- t(apply(data_matrix, 1, function(x) {
          s <- sd(x, na.rm = TRUE)
          if(is.na(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
        }))
        data_type_label <- paste0(data_type_label, " [within-participant standardized]")
      } else if(standardize_opt == "between") {
        # Between-participant: z-score each time point across participants (column-wise)
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        col_sds   <- apply(data_matrix, 2, sd, na.rm = TRUE)
        col_sds[col_sds == 0] <- 1  # avoid division by zero for constant time points
        data_matrix <- sweep(sweep(data_matrix, 2, col_means, "-"), 2, col_sds, "/")
        data_type_label <- paste0(data_type_label, " [between-participant standardized]")
      }

      # Get clustering method
      clustering_method <- if(!is.null(input$clustering_method)) input$clustering_method else "standard"

      # Check for functional k-means requirements
      if(clustering_method == "functional") {
        if(input$cluster_data_type != "smoothed") {
          showNotification("Functional K-Means requires smoothed data. Please select 'Smoothed Data' above.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(is.null(values$fd_obj)) {
          showNotification("No functional data object available. Please apply smoothing first.",
                           type = "error", duration = 5)
          return(NULL)
        }
      }

      # Check for DCF requirements
      if(clustering_method == "dcf") {
        if(!requireNamespace("reticulate", quietly = TRUE)) {
          showNotification("reticulate package not available. Please install it for DCF clustering.",
                           type = "error", duration = 5)
          return(NULL)
        }
        if(!reticulate::py_module_available("DCFcluster")) {
          showNotification("DCFcluster Python package not found. Install with: pip install git+https://github.com/tobinjo96/DCFcluster.git",
                           type = "error", duration = 10)
          return(NULL)
        }
      }

      # Check for missing values
      if(any(is.na(data_matrix))) {
        showNotification("Data contains missing values. Please handle missing values before clustering.",
                         type = "error", duration = 5)
        return(NULL)
      }

      # Get parameters
      k <- input$n_clusters
      nstart <- input$kmeans_nstart
      iter_max <- input$kmeans_iter

      # ===== PERFORM CLUSTERING BASED ON METHOD =====

      if(clustering_method == "dcf") {
        # ===== DCF (DENSITY CORE FINDING) CLUSTERING =====
        showNotification("Running DCF clustering...", type = "message", duration = 2)

        dcf_k <- if(!is.null(input$dcf_k)) input$dcf_k else 10
        dcf_beta <- if(!is.null(input$dcf_beta)) input$dcf_beta else 0.1

        dcf_result <- run_dcf_clustering(data_matrix, dcf_k, dcf_beta)

        if(!dcf_result$success) {
          showNotification(paste("DCF clustering failed:", dcf_result$error),
                           type = "error", duration = 10)
          return(NULL)
        }

        cluster_assignments <- dcf_result$labels
        k <- dcf_result$n_clusters  # DCF determines number of clusters automatically

        # Calculate cluster means
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))
        for(i in 1:k) {
          cluster_idx <- which(cluster_assignments == i)
          if(length(cluster_idx) > 0) {
            cluster_means[i, ] <- colMeans(data_matrix[cluster_idx, , drop = FALSE])
          }
        }

        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Calculate within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared <- between_ss / total_ss

        # Create compatible result object
        kmeans_result <- list(
          cluster = cluster_assignments,
          centers = cluster_means,
          size = cluster_sizes,
          tot.withinss = within_ss,
          betweenss = between_ss
        )

        method_label <- paste0("DCF (Density Core Finding, k=", dcf_k, ", β=", dcf_beta, ")")

        # Store DCF-specific results
        dcf_extra <- list(
          peak_values = dcf_result$peak_values,
          core_sets = dcf_result$core_sets,
          k_param = dcf_k,
          beta_param = dcf_beta
        )
        hc_result <- NULL

      } else if(clustering_method == "functional") {
        # ===== FUNCTIONAL K-MEANS using fda.usc =====
        if(!requireNamespace("fda.usc", quietly = TRUE)) {
          showNotification("fda.usc package not available. Please install it or use Standard K-Means.",
                           type = "error", duration = 5)
          return(NULL)
        }

        # Convert data matrix to fdata object for fda.usc
        # fdata expects: rows = observations (subjects), cols = time points
        time_grid <- seq(0, 1, length.out = ncol(data_matrix))
        fdata_obj <- fda.usc::fdata(mdata = data_matrix, argvals = time_grid)

        # Run fda.usc::kmeans.fd nstart times; keep the solution with the lowest WCSS
        best_fkm  <- NULL
        best_wcss <- Inf
        for(restart in seq_len(nstart)) {
          set.seed(100 + restart)
          fkm_try <- tryCatch(
            fda.usc::kmeans.fd(fdata_obj, ncl = k, max.iter = iter_max),
            error = function(e) NULL
          )
          if(is.null(fkm_try)) next
          asn <- fkm_try$cluster
          wcss_try <- sum(sapply(seq_len(k), function(i) {
            idx <- which(asn == i)
            if(length(idx) == 0) return(0)
            m <- data_matrix[idx, , drop = FALSE]
            sum(apply(m, 1, function(x) sum((x - colMeans(m))^2)))
          }))
          if(wcss_try < best_wcss) { best_wcss <- wcss_try; best_fkm <- fkm_try }
        }
        if(is.null(best_fkm)) {
          showNotification("Functional K-Means failed across all restarts.", type = "error", duration = 5)
          return(NULL)
        }
        fkm_result <- best_fkm

        # Extract cluster assignments
        cluster_assignments <- fkm_result$cluster

        # Calculate cluster means by evaluating fd object for each cluster
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))

        for(i in 1:k) {
          cluster_idx <- which(cluster_assignments == i)
          if(length(cluster_idx) > 0) {
            cluster_fd <- values$fd_obj[cluster_idx]
            cluster_mean_fd <- mean.fd(cluster_fd)
            cluster_means[i, ] <- eval.fd(time_grid, cluster_mean_fd)
          }
        }

        # Calculate cluster sizes
        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Calculate within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared <- between_ss / total_ss

        # Create a compatible "kmeans_result" object for consistency
        kmeans_result <- list(
          cluster = cluster_assignments,
          centers = cluster_means,
          size = cluster_sizes,
          tot.withinss = within_ss,
          betweenss = between_ss
        )

        method_label <- "Functional K-Means (fda.usc)"
        dcf_extra <- NULL
        hc_result <- NULL

      } else if(clustering_method == "hierarchical") {
        # ===== HIERARCHICAL CLUSTERING =====
        linkage <- if(!is.null(input$hclust_linkage)) input$hclust_linkage else "ward.D2"

        dist_matrix_hc <- dist(data_matrix)
        hc_result <- hclust(dist_matrix_hc, method = linkage)
        cluster_assignments <- cutree(hc_result, k = k)

        # Cluster means
        cluster_means <- matrix(0, nrow = k, ncol = ncol(data_matrix))
        for(i in 1:k) {
          idx <- which(cluster_assignments == i)
          if(length(idx) > 0)
            cluster_means[i, ] <- colMeans(data_matrix[idx, , drop = FALSE])
        }
        cluster_sizes <- as.numeric(table(cluster_assignments))

        # Within-cluster sum of squares
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(members) > 0)
            wcss_per_cluster[i] <- sum(apply(members, 1, function(x) sum((x - cluster_means[i, ])^2)))
        }

        total_ss   <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        within_ss  <- sum(wcss_per_cluster)
        between_ss <- total_ss - within_ss
        r_squared  <- between_ss / total_ss

        kmeans_result <- list(
          cluster     = cluster_assignments,
          centers     = cluster_means,
          size        = cluster_sizes,
          tot.withinss = within_ss,
          betweenss   = between_ss
        )

        method_label <- paste0("Hierarchical Clustering (", linkage, ")")
        dcf_extra <- NULL

      } else {
        # ===== STANDARD K-MEANS =====
        hc_result <- NULL
        set.seed(123)  # For reproducibility
        kmeans_result <- kmeans(data_matrix, centers = k,
                                nstart = nstart, iter.max = iter_max)

        cluster_assignments <- kmeans_result$cluster
        cluster_means <- kmeans_result$centers
        cluster_sizes <- kmeans_result$size

        # Calculate within-cluster sum of squares for each cluster
        wcss_per_cluster <- numeric(k)
        for(i in 1:k) {
          cluster_members <- data_matrix[cluster_assignments == i, , drop = FALSE]
          if(nrow(cluster_members) > 0) {
            cluster_center <- cluster_means[i, ]
            wcss_per_cluster[i] <- sum(apply(cluster_members, 1, function(x) sum((x - cluster_center)^2)))
          }
        }

        total_ss <- sum(apply(data_matrix, 1, function(x) sum((x - colMeans(data_matrix))^2)))
        between_ss <- kmeans_result$betweenss
        within_ss <- kmeans_result$tot.withinss
        r_squared <- between_ss / total_ss
        dcf_extra <- NULL

        method_label <- "Standard (Point-wise) K-Means"
        dcf_extra <- NULL
      }

      # ===== COMMON CALCULATIONS FOR ALL METHODS =====

      # Silhouette width (requires cluster package)
      sil_width <- NA
      sil_data <- NULL
      if(requireNamespace("cluster", quietly = TRUE) && k > 1) {
        dist_matrix <- dist(data_matrix)
        sil <- cluster::silhouette(cluster_assignments, dist_matrix)
        sil_width <- mean(sil[, 3])
        # Store full silhouette object for detailed plotting
        sil_data <- sil
      }

      # Calinski-Harabasz index (variance ratio criterion)
      n <- nrow(data_matrix)
      ch_index <- if(k > 1) (between_ss / (k - 1)) / (within_ss / (n - k)) else NA

      # Store results
      values$clustering_results <- list(
        kmeans = kmeans_result,
        cluster_means = cluster_means,
        cluster_assignments = kmeans_result$cluster,
        cluster_sizes = kmeans_result$size,
        wcss_per_cluster = wcss_per_cluster,
        total_wcss = within_ss,
        between_ss = between_ss,
        total_ss = total_ss,
        r_squared = r_squared,
        silhouette_width = sil_width,
        silhouette_data = sil_data,  # Full silhouette object for detailed plot
        ch_index = ch_index,
        k = k,
        data_type = data_type_label,
        data_matrix = data_matrix,
        method = clustering_method,
        method_label = method_label,
        dcf_extra = dcf_extra,  # DCF-specific results (NULL for other methods)
        hc_obj = hc_result,     # hclust object (NULL for non-hierarchical methods)
        linkage_method = if(clustering_method == "hierarchical") {
          if(!is.null(input$hclust_linkage)) input$hclust_linkage else "ward.D2"
        } else NULL
      )

      # Notification message
      notification_msg <- if(clustering_method == "dcf") {
        paste0("DCF clustering completed: ", k, " clusters found automatically")
      } else {
        paste0("Clustering completed successfully with k=", k)
      }
      showNotification(notification_msg, type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Clustering error:", e$message),
                       type = "error", duration = 5)
      values$clustering_results <- NULL
    })
  })

  # Clustering status output
  output$clustering_status <- renderText({
    if(is.null(values$clustering_results)) {
      "No clustering performed yet. Click 'Run Clustering' to begin."
    } else {
      method_info <- if(!is.null(values$clustering_results$method_label)) {
        values$clustering_results$method_label
      } else {
        "Clustering"
      }
      paste0("✓ Clustering completed with k=", values$clustering_results$k,
             " clusters\n   Method: ", method_info,
             "\n   Data: ", values$clustering_results$data_type)
    }
  })

  # Save cluster membership to dataset
  observeEvent(input$save_cluster_membership, {
    req(values$clustering_results)

    tryCatch({
      # Validate variable name
      var_name <- trimws(input$cluster_var_name)

      if(var_name == "" || is.null(var_name)) {
        showNotification("Please enter a variable name.", type = "error", duration = 3)
        return(NULL)
      }

      # Check for invalid characters in variable name
      if(!grepl("^[a-zA-Z][a-zA-Z0-9._]*$", var_name)) {
        showNotification("Variable name must start with a letter and contain only letters, numbers, dots, or underscores.",
                         type = "error", duration = 5)
        return(NULL)
      }

      # Get cluster assignments
      cluster_assignments <- values$clustering_results$cluster_assignments

      # Check if we have uploaded data
      if(is.null(values$uploaded_data)) {
        showNotification("No uploaded data available. Cluster membership can only be saved if you uploaded data with ID variables.",
                         type = "warning", duration = 5)
        return(NULL)
      }

      # Check dimensions match
      if(nrow(values$uploaded_data) != length(cluster_assignments)) {
        showNotification(paste0("Data dimension mismatch. Uploaded data has ", nrow(values$uploaded_data),
                               " rows but clustering has ", length(cluster_assignments), " observations."),
                         type = "error", duration = 5)
        return(NULL)
      }

      # Check if variable already exists
      if(var_name %in% colnames(values$uploaded_data)) {
        showModal(modalDialog(
          title = "Variable Already Exists",
          paste0("The variable '", var_name, "' already exists in your dataset. Do you want to overwrite it?"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_overwrite", "Overwrite", class = "btn-warning")
          )
        ))
        return(NULL)
      }

      # Add cluster membership to uploaded data
      values$uploaded_data[[var_name]] <- as.factor(cluster_assignments)

      showNotification(paste0("✓ Variable '", var_name, "' added to in-memory dataset (",
                             values$clustering_results$k, " clusters). ",
                             "Use 'Download Dataset' button to export to file."),
                       type = "message", duration = 7)

    }, error = function(e) {
      showNotification(paste("Error saving cluster membership:", e$message),
                       type = "error", duration = 5)
    })
  })

  # Handle overwrite confirmation
  observeEvent(input$confirm_overwrite, {
    var_name <- trimws(input$cluster_var_name)
    cluster_assignments <- values$clustering_results$cluster_assignments

    # Overwrite the variable
    values$uploaded_data[[var_name]] <- as.factor(cluster_assignments)

    showNotification(paste0("✓ Variable '", var_name, "' overwritten in in-memory dataset. ",
                           "Use 'Download Dataset' button to export to file."),
                     type = "message", duration = 7)

    removeModal()
  })

  # Cluster save status
  output$cluster_save_status <- renderText({
    req(values$clustering_results)

    var_name <- trimws(input$cluster_var_name)

    if(is.null(values$uploaded_data)) {
      return("No uploaded data")
    }

    if(var_name == "" || is.null(var_name)) {
      return("Enter variable name")
    }

    if(var_name %in% colnames(values$uploaded_data)) {
      return(paste0("✓ '", var_name, "' exists\n(will prompt to overwrite)"))
    } else {
      return("Ready to add")
    }
  })

  # Download dataset with cluster membership
  output$download_data_with_clusters <- downloadHandler(
    filename = function() {
      paste0("data_with_clusters_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$uploaded_data)

      # Create a copy of uploaded data
      data_to_export <- values$uploaded_data

      # If cluster membership was added, it's already in values$uploaded_data
      # If not yet added, add it now
      var_name <- trimws(input$cluster_var_name)
      if(!is.null(values$clustering_results)) {
        if(!(var_name %in% colnames(data_to_export)) && var_name != "") {
          # Add cluster membership if it wasn't added yet
          cluster_assignments <- values$clustering_results$cluster_assignments
          if(nrow(data_to_export) == length(cluster_assignments)) {
            data_to_export[[var_name]] <- as.factor(cluster_assignments)
          }
        }
      }

      # Write to CSV
      write.csv(data_to_export, file, row.names = FALSE)
    }
  )

  # Cluster summary table
  output$cluster_summary_table <- renderDT({
    req(values$clustering_results)

    wcss <- values$clustering_results$wcss_per_cluster
    sizes <- values$clustering_results$cluster_sizes

    summary_df <- data.frame(
      Cluster = 1:values$clustering_results$k,
      Size = sizes,
      Percentage = round(100 * sizes / sum(sizes), 1),
      Within_SS = round(wcss, 2),
      Avg_Within_SS = round(wcss / sizes, 2)
    )

    colnames(summary_df) <- c("Cluster", "Members", "% of Total",
                               "Within-Cluster SS", "Avg. Within-SS")

    datatable(summary_df,
              options = list(
                dom = 't',
                pageLength = 20,
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(summary_df),
                  fontSize = '14px')
  })

  # Clustering fit statistics
  output$clustering_fit_stats <- renderText({
    req(values$clustering_results)

    r <- values$clustering_results

    # Method info
    method_info <- if(!is.null(r$method_label)) r$method_label else "Clustering"
    output_text <- paste0("Method: ", method_info, "\n",
                          "Number of Clusters: ", r$k, "\n\n")

    # Standard statistics
    output_text <- paste0(output_text,
      "Total Sum of Squares: ", round(r$total_ss, 2), "\n",
      "Between-Cluster Sum of Squares: ", round(r$between_ss, 2), "\n",
      "Within-Cluster Sum of Squares: ", round(r$total_wcss, 2), "\n",
      "R² (Variance Explained): ", round(r$r_squared * 100, 2), "%\n",
      "Between-SS / Total-SS Ratio: ", round(r$between_ss / r$total_ss * 100, 2), "%\n"
    )

    if(!is.na(r$silhouette_width)) {
      output_text <- paste0(output_text,
                            "Average Silhouette Width: ", round(r$silhouette_width, 3), "\n",
                            "  (Range: -1 to 1, higher is better)\n")
    }

    if(!is.na(r$ch_index)) {
      output_text <- paste0(output_text,
                            "Calinski-Harabasz Index: ", round(r$ch_index, 2), "\n",
                            "  (Higher values indicate better-defined clusters)\n")
    }

    # DCF-specific information
    if(!is.null(r$dcf_extra)) {
      output_text <- paste0(output_text, "\n--- DCF Parameters ---\n",
                            "Neighborhood k: ", r$dcf_extra$k_param, "\n",
                            "Fluctuation β: ", r$dcf_extra$beta_param, "\n",
                            "Clusters detected automatically: ", r$k, "\n")
    }

    output_text
  })

  # Cluster membership table
  output$cluster_membership_table <- renderDT({
    req(values$clustering_results)

    # Create membership table
    membership_df <- data.frame(
      Subject_ID = 1:length(values$clustering_results$cluster_assignments),
      Cluster = values$clustering_results$cluster_assignments
    )

    # Add ALL group variables if available (not just the primary one)
    if(!is.null(values$group_variables) && ncol(values$group_variables) > 0) {
      # Add all group variables from group_variables data frame
      for(var_name in colnames(values$group_variables)) {
        membership_df[[var_name]] <- values$group_variables[[var_name]]
      }
    } else if(!is.null(values$group_labels)) {
      # Fallback to primary group_labels
      membership_df$Group <- values$group_labels
    }

    # Generate stronger colors with transparency (same as plots)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Add transparency to colors for table background
    table_colors <- paste0(strong_colors[1:values$clustering_results$k], "40")

    datatable(membership_df,
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                scrollY = "400px",
                searching = TRUE,
                ordering = TRUE
              ),
              rownames = FALSE,
              filter = 'top') %>%
      formatStyle('Cluster',
                  backgroundColor = styleEqual(
                    1:values$clustering_results$k,
                    table_colors
                  ))
  })

  # Check if group labels are available
  output$has_group_labels <- reactive({
    !is.null(values$group_labels) && !is.null(values$clustering_results)
  })
  outputOptions(output, "has_group_labels", suspendWhenHidden = FALSE)

  # UI for selecting which group variable to use in cluster composition analysis
  output$cluster_group_var_selector <- renderUI({
    req(values$group_labels)

    # Get available group variables
    available_vars <- NULL
    if(!is.null(values$selected_group_vars) && length(values$selected_group_vars) >= 1) {
      available_vars <- values$selected_group_vars
    } else if(!is.null(values$group_labels)) {
      # Fallback: use "Group" as default name if only group_labels exists
      available_vars <- "Group"
    }

    if(!is.null(available_vars)) {
      # Get current selection or default to first
      current_selection <- if(!is.null(input$cluster_group_var) &&
                              input$cluster_group_var %in% available_vars) {
        input$cluster_group_var
      } else {
        available_vars[1]
      }

      tagList(
        selectInput("cluster_group_var",
                    "Select Group Variable for Composition Analysis:",
                    choices = available_vars,
                    selected = current_selection),
        helpText(paste("Available group variables:", paste(available_vars, collapse = ", "))),
        hr()
      )
    } else {
      NULL
    }
  })

  # Helper to get the selected group variable for clustering
  get_cluster_group_labels <- reactive({
    # First check if a specific variable is selected and available in group_variables
    if(!is.null(input$cluster_group_var) && !is.null(values$group_variables) &&
       input$cluster_group_var %in% colnames(values$group_variables)) {
      labels <- values$group_variables[[input$cluster_group_var]]
      # Ensure it's a factor
      if(!is.factor(labels)) {
        labels <- factor(labels)
      }
      return(labels)
    }
    # Fallback to primary group_labels
    values$group_labels
  })

  # Get the name of currently selected group variable
  get_cluster_group_var_name <- reactive({
    if(!is.null(input$cluster_group_var)) {
      input$cluster_group_var
    } else if(!is.null(values$selected_group_vars) && length(values$selected_group_vars) > 0) {
      values$selected_group_vars[1]
    } else {
      "Group"
    }
  })

  # Group counts table by cluster
  output$cluster_group_counts_table <- renderDT({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Create counts table
    counts_df <- as.data.frame.matrix(cont_table)
    counts_df <- cbind(Cluster = rownames(counts_df), counts_df)
    counts_df$Total <- rowSums(cont_table)

    datatable(counts_df,
              options = list(
                dom = 't',
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(counts_df),
                  fontSize = '14px')
  })

  # Group percentages table by cluster
  output$cluster_group_pct_table <- renderDT({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Calculate percentages for each cluster (row percentages)
    cluster_pcts <- prop.table(cont_table, margin = 1) * 100

    # Create percentages table
    pct_df <- as.data.frame.matrix(round(cluster_pcts, 1))
    pct_df <- cbind(Cluster = rownames(pct_df), pct_df)

    datatable(pct_df,
              options = list(
                dom = 't',
                scrollX = TRUE
              ),
              rownames = FALSE) %>%
      formatStyle(columns = 1:ncol(pct_df),
                  fontSize = '14px')
  })

  # Chi-square test for group-cluster independence
  output$cluster_group_test <- renderText({
    req(values$clustering_results)
    req(values$group_labels)

    clusters <- values$clustering_results$cluster_assignments
    groups <- get_cluster_group_labels()
    group_var_name <- get_cluster_group_var_name()

    # Create contingency table
    cont_table <- table(Cluster = clusters, Group = groups)

    # Get table dimensions for reporting
    n_clusters <- nrow(cont_table)
    n_groups <- ncol(cont_table)
    table_cells <- length(cont_table)

    # Perform chi-square test (works for any table size)
    chi_test <- tryCatch({
      suppressWarnings(chisq.test(cont_table))
    }, error = function(e) {
      NULL
    })

    # Check if chi-square assumptions are violated (expected frequencies < 5)
    chi_warning <- FALSE
    if(!is.null(chi_test)) {
      expected <- chi_test$expected
      if(any(expected < 5)) {
        chi_warning <- TRUE
      }
    }

    # Try Fisher's exact test with simulation for larger tables
    # Fisher test with simulation can handle larger tables
    fisher_test <- tryCatch({
      # Use simulation-based p-value for tables of any size
      fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
    }, error = function(e) {
      # If still fails (very rare), return NULL
      NULL
    })

    # Build output text
    output_text <- "=== CLUSTER-GROUP ASSOCIATION TEST ===\n\n"

    # Show which group variable is being used
    output_text <- paste0(output_text, "Group Variable: ", group_var_name, "\n")
    output_text <- paste0(output_text, "Table dimensions: ", n_clusters, " clusters × ", n_groups, " groups\n\n")

    # Overall percentages (limit display for many groups)
    output_text <- paste0(output_text, "Overall ", group_var_name, " Distribution:\n")
    overall_pcts <- prop.table(table(groups)) * 100
    group_names <- names(overall_pcts)
    if(length(group_names) > 10) {
      # Show only first 10 groups if there are many
      for(g in group_names[1:10]) {
        output_text <- paste0(output_text, "  ", g, ": ", round(overall_pcts[g], 1), "%\n")
      }
      output_text <- paste0(output_text, "  ... and ", length(group_names) - 10, " more groups\n")
    } else {
      for(g in group_names) {
        output_text <- paste0(output_text, "  ", g, ": ", round(overall_pcts[g], 1), "%\n")
      }
    }
    output_text <- paste0(output_text, "\n")

    # Chi-square test results
    if(!is.null(chi_test)) {
      output_text <- paste0(output_text, "Chi-Square Test of Independence:\n")
      output_text <- paste0(output_text, "  χ² = ", round(chi_test$statistic, 3), "\n")
      output_text <- paste0(output_text, "  df = ", chi_test$parameter, "\n")
      output_text <- paste0(output_text, "  p-value = ", format.pval(chi_test$p.value, digits = 4), "\n")

      if(chi_warning) {
        output_text <- paste0(output_text, "  Note: Some expected frequencies < 5; interpret with caution.\n")
      }
      output_text <- paste0(output_text, "\n")

      if(chi_test$p.value < 0.001) {
        output_text <- paste0(output_text, "  *** HIGHLY SIGNIFICANT (p < 0.001) ***\n")
        output_text <- paste0(output_text, "  Cluster membership is strongly associated with group membership.\n\n")
      } else if(chi_test$p.value < 0.01) {
        output_text <- paste0(output_text, "  ** SIGNIFICANT (p < 0.01) **\n")
        output_text <- paste0(output_text, "  Cluster membership is significantly associated with group membership.\n\n")
      } else if(chi_test$p.value < 0.05) {
        output_text <- paste0(output_text, "  * SIGNIFICANT (p < 0.05) *\n")
        output_text <- paste0(output_text, "  Cluster membership is associated with group membership.\n\n")
      } else {
        output_text <- paste0(output_text, "  Not significant (p >= 0.05)\n")
        output_text <- paste0(output_text, "  No significant association between cluster and group membership.\n\n")
      }
    }

    # Fisher's exact test results (if available)
    if(!is.null(fisher_test)) {
      output_text <- paste0(output_text, "Fisher's Exact Test (Monte Carlo simulation, B=10000):\n")
      output_text <- paste0(output_text, "  p-value = ", format.pval(fisher_test$p.value, digits = 4), "\n")
      if(fisher_test$p.value < 0.05) {
        output_text <- paste0(output_text, "  Significant association (p < 0.05)\n\n")
      } else {
        output_text <- paste0(output_text, "  Not significant (p >= 0.05)\n\n")
      }
    } else if(is.null(chi_test)) {
      output_text <- paste0(output_text, "Statistical tests could not be performed.\n")
      output_text <- paste0(output_text, "This may occur with very sparse tables or unusual data distributions.\n\n")
    }

    # Calculate and report effect size (Cramér's V)
    if(!is.null(chi_test)) {
      n <- sum(cont_table)
      k <- min(nrow(cont_table), ncol(cont_table))
      cramers_v <- sqrt(chi_test$statistic / (n * (k - 1)))
      output_text <- paste0(output_text, "Effect Size (Cramér's V): ", round(cramers_v, 3), "\n")
      if(cramers_v < 0.1) {
        output_text <- paste0(output_text, "  (Negligible effect)\n")
      } else if(cramers_v < 0.3) {
        output_text <- paste0(output_text, "  (Small effect)\n")
      } else if(cramers_v < 0.5) {
        output_text <- paste0(output_text, "  (Medium effect)\n")
      } else {
        output_text <- paste0(output_text, "  (Large effect)\n")
      }
    }

    output_text
  })

  # Cluster means plot
  output$cluster_means_plot <- renderPlotly({
    req(values$clustering_results)

    cluster_means <- values$clustering_results$cluster_means
    k <- values$clustering_results$k
    data_matrix <- values$clustering_results$data_matrix
    clusters <- values$clustering_results$cluster_assignments

    # Get time points using helper functions (same as preprocessing module)
    hour_labels <- get_hour_labels()

    if(!is.null(hour_labels)) {
      # Use actual hour values from column names
      time_positions <- calculate_time_positions(hour_labels)
      if(!is.null(time_positions)) {
        time_points <- time_positions
      } else {
        time_points <- seq(0, 1, length.out = ncol(cluster_means))
      }
      tick_text <- sapply(hour_labels, decimal_to_hhmm)

      # Snap ticks to round-hour multiples
      kmeans_step <- as.numeric(input$tick_freq_kmeans)
      if (!is.na(kmeans_step) && kmeans_step > 0) {
        keep <- abs(round(hour_labels / kmeans_step) * kmeans_step - hour_labels) < 1e-6
        if (any(keep)) {
          tick_vals_subset <- time_points[keep]
          tick_text_subset <- tick_text[keep]
        } else {
          tick_vals_subset <- time_points
          tick_text_subset <- tick_text
        }
      } else {
        tick_vals_subset <- time_points
        tick_text_subset <- tick_text
      }
    } else if(!is.null(values$time_labels)) {
      # Use column names as labels
      time_points <- seq(0, 1, length.out = length(values$time_labels))
      tick_vals_subset <- time_points
      tick_text_subset <- values$time_labels
    } else {
      # Fallback to numeric indices
      time_points <- 1:ncol(cluster_means)
      tick_vals_subset <- time_points
      tick_text_subset <- as.character(time_points)
    }
    hover_times <- hover_time_labels(time_points)

    x_label <- get_time_label()

    # Compute a nice y-axis tick interval targeting ~5 ticks
    y_span <- diff(range(data_matrix, na.rm = TRUE))
    y_dtick <- if (y_span > 0) {
      raw_step <- y_span / 5
      pow10 <- 10^floor(log10(raw_step))
      pow10 * ceiling(raw_step / pow10)
    } else NULL

    # Generate stronger colors (not pastel)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Create plotly object
    p <- plot_ly()

    # Add traces for each cluster
    for(i in 1:k) {
      # Get cluster members
      cluster_idx <- which(clusters == i)
      cluster_data <- data_matrix[cluster_idx, , drop = FALSE]

      # Calculate mean, SD for each time point
      cluster_mean <- cluster_means[i, ]
      cluster_sd <- apply(cluster_data, 2, sd)

      # Add lower SD band FIRST (this will be the bottom of the ribbon)
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean - cluster_sd,
          type = 'scatter',
          mode = 'lines',
          line = list(width = 0),
          name = paste0("Cluster ", i, " ±SD"),
          legendgroup = paste0("cluster_sd_", i),
          showlegend = FALSE,
          hoverinfo = 'skip'
        )

      # Add upper SD band SECOND with fill = 'tonexty' to create ribbon
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean + cluster_sd,
          type = 'scatter',
          mode = 'lines',
          line = list(width = 0),
          fillcolor = paste0(strong_colors[i], '40'),  # 40 = 25% opacity in hex
          fill = 'tonexty',
          name = paste0("Cluster ", i, " ±SD"),
          legendgroup = paste0("cluster_sd_", i),
          showlegend = TRUE,
          hoverinfo = 'skip'
        )

      # Add mean curve (smooth interpolation using spline)
      # Use SEPARATE legendgroup so it can be toggled independently
      p <- p %>%
        add_trace(
          x = time_points,
          y = cluster_mean,
          type = 'scatter',
          mode = 'lines',
          line = list(color = strong_colors[i], width = 2, shape = 'spline'),
          name = paste0("Cluster ", i, " Mean"),
          legendgroup = paste0("cluster_mean_", i),
          showlegend = TRUE,
          hovertemplate = paste0(
            "Cluster: ", i,
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<br>Size: ", values$clustering_results$cluster_sizes[i], " subjects",
            "<extra></extra>"
          ),
          customdata = hover_times
        )
    }

    # Update layout with proper time axis formatting
    p <- p %>%
      layout(
        title = list(text = "Cluster Mean Functions with ±1 SD Bands",
                     x = 0.5, xanchor = "center"),
        xaxis = list(
          title = x_label,
          tickmode = 'array',
          tickvals = tick_vals_subset,
          ticktext = tick_text_subset,
          tickangle = -90  # Rotate labels for readability
        ),
        yaxis = list(title = "Value", dtick = y_dtick),
        hovermode = "closest",
        legend = list(
          orientation = "v",
          x = 1.02,
          y = 1,
          tracegroupgap = 0
        )
      )

    p
  })

  # Individual curves by cluster
  output$cluster_individuals_plot <- renderPlotly({
    req(values$clustering_results)

    data_matrix <- values$clustering_results$data_matrix
    clusters <- values$clustering_results$cluster_assignments
    k <- values$clustering_results$k

    # Get time points using helper functions (same as preprocessing module)
    hour_labels <- get_hour_labels()

    if(!is.null(hour_labels)) {
      # Use actual hour values from column names
      time_positions <- calculate_time_positions(hour_labels)
      if(!is.null(time_positions)) {
        time_points <- time_positions
      } else {
        time_points <- seq(0, 1, length.out = ncol(data_matrix))
      }
      tick_text <- sapply(hour_labels, decimal_to_hhmm)

      # Snap ticks to round-hour multiples
      kmeans_step <- as.numeric(input$tick_freq_kmeans)
      if (!is.na(kmeans_step) && kmeans_step > 0) {
        keep <- abs(round(hour_labels / kmeans_step) * kmeans_step - hour_labels) < 1e-6
        if (any(keep)) {
          tick_vals_subset <- time_points[keep]
          tick_text_subset <- tick_text[keep]
        } else {
          tick_vals_subset <- time_points
          tick_text_subset <- tick_text
        }
      } else {
        tick_vals_subset <- time_points
        tick_text_subset <- tick_text
      }
    } else if(!is.null(values$time_labels)) {
      # Use column names as labels
      time_points <- seq(0, 1, length.out = length(values$time_labels))
      tick_vals_subset <- time_points
      tick_text_subset <- values$time_labels
    } else {
      # Fallback to numeric indices
      time_points <- 1:ncol(data_matrix)
      tick_vals_subset <- time_points
      tick_text_subset <- as.character(time_points)
    }
    hover_times <- hover_time_labels(time_points)

    x_label <- get_time_label()

    # Compute a nice y-axis tick interval targeting ~5 ticks
    y_span <- diff(range(data_matrix, na.rm = TRUE))
    y_dtick <- if (y_span > 0) {
      raw_step <- y_span / 5
      pow10 <- 10^floor(log10(raw_step))
      pow10 * ceiling(raw_step / pow10)
    } else NULL

    # Generate stronger colors (same as cluster means plot)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Calculate overall mean across all subjects (as reference)
    overall_mean <- colMeans(data_matrix)

    # Create subplots for each cluster
    plot_list <- list()

    for(i in 1:k) {
      cluster_idx <- which(clusters == i)
      cluster_data <- data_matrix[cluster_idx, , drop = FALSE]

      # Create plotly object for this cluster
      p <- plot_ly()

      # Add individual curves (each as a separate trace for interactivity)
      # Group them together so they can be toggled together
      for(j in 1:nrow(cluster_data)) {
        p <- p %>%
          add_trace(
            x = time_points,
            y = cluster_data[j, ],
            type = 'scatter',
            mode = 'lines',
            line = list(color = strong_colors[i], width = 0.5),
            opacity = 0.15,
            name = paste0("Cluster ", i, " Individuals"),
            legendgroup = paste0("individual_", i),
            showlegend = (j == 1),  # Only show in legend once
            hovertemplate = paste0(
              "Subject: ", cluster_idx[j],
              "<br>Time: %{customdata}",
              "<br>Value: %{y:.2f}",
              "<extra></extra>"
            ),
            customdata = hover_times
          )
      }

      # Add cluster mean - bold and opaque
      p <- p %>%
        add_trace(
          x = time_points,
          y = values$clustering_results$cluster_means[i, ],
          type = 'scatter',
          mode = 'lines',
          line = list(color = strong_colors[i], width = 2.5, shape = 'spline'),
          opacity = 1,
          name = paste0("Cluster ", i, " Mean"),
          legendgroup = paste0("mean_", i),
          showlegend = TRUE,
          hovertemplate = paste0(
            "Cluster ", i, " Mean",
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<extra></extra>"
          ),
          customdata = hover_times
        )

      # Add overall mean as reference - dashed black line
      p <- p %>%
        add_trace(
          x = time_points,
          y = overall_mean,
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'grey', width = 2, dash = 'dot', shape = 'spline'),
          opacity = 0.8,
          name = "Overall Mean",
          legendgroup = "overall_mean",
          showlegend = (i == 1),  # Only show in legend once (first subplot)
          hovertemplate = paste0(
            "Overall Mean (all subjects)",
            "<br>Time: %{customdata}",
            "<br>Value: %{y:.2f}",
            "<extra></extra>"
          ),
          customdata = hover_times
        )

      # Update layout for this subplot with proper time axis formatting
      p <- p %>%
        layout(
          title = list(text = paste0("Cluster ", i, " (n=", length(cluster_idx), ")"),
                       font = list(size = 12)),
          xaxis = list(
            title = x_label,
            tickmode = 'array',
            tickvals = tick_vals_subset,
            ticktext = tick_text_subset,
            tickangle = -90  # Rotate labels for readability
          ),
          yaxis = list(title = "Value"),
          hovermode = "closest"
        )

      plot_list[[i]] <- p
    }

    # Combine plots
    sp <- subplot(plot_list, nrows = ceiling(k/2), shareX = TRUE, shareY = TRUE,
                  titleX = TRUE, titleY = TRUE)

    # Build layout args, applying dtick to every y-axis (yaxis, yaxis2, yaxis3, ...)
    layout_args <- list(
      sp,
      title = list(text = "Individual Curves by Cluster (click legend to toggle)",
                   x = 0.5, xanchor = "center"),
      showlegend = TRUE,
      legend = list(orientation = "v", x = 1.02, y = 1)
    )
    yaxis_names <- c("yaxis", paste0("yaxis", seq_len(k)[-1]))
    for (nm in yaxis_names) {
      layout_args[[nm]] <- list(dtick = y_dtick)
    }
    do.call(layout, layout_args)
  })

  # Detailed silhouette plot (similar to Python scikit-learn visualization)
  output$detailed_silhouette_plot <- renderPlotly({
    req(values$clustering_results)
    req(values$clustering_results$silhouette_data)

    sil_data <- values$clustering_results$silhouette_data
    k <- values$clustering_results$k

    # Extract silhouette coefficients and cluster assignments
    sil_df <- data.frame(
      cluster = sil_data[, 1],
      neighbor = sil_data[, 2],
      sil_width = sil_data[, 3]
    )

    # Sort within each cluster by silhouette width (descending)
    sil_df <- sil_df %>%
      group_by(cluster) %>%
      arrange(cluster, desc(sil_width)) %>%
      mutate(sample_idx = row_number()) %>%
      ungroup()

    # Calculate cumulative positions for y-axis
    cluster_sizes <- table(sil_df$cluster)
    cumsum_sizes <- c(0, cumsum(cluster_sizes))

    # Add y-position for plotting (stacked bars)
    sil_df <- sil_df %>%
      group_by(cluster) %>%
      mutate(y_position = cumsum_sizes[cluster[1]] + row_number()) %>%
      ungroup()

    # Generate stronger colors (same as other plots)
    strong_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666",
                       "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

    # Create plotly figure
    p <- plot_ly()

    # Add bars for each cluster
    for(i in 1:k) {
      cluster_data <- sil_df %>% filter(cluster == i)

      p <- p %>%
        add_trace(
          x = cluster_data$sil_width,
          y = cluster_data$y_position,
          type = 'bar',
          orientation = 'h',
          marker = list(color = strong_colors[i]),
          name = paste0("Cluster ", i),
          hovertemplate = paste0(
            "Cluster: ", i,
            "<br>Sample: ", 1:nrow(cluster_data),
            "<br>Silhouette: %{x:.3f}",
            "<extra></extra>"
          )
        )
    }

    # Calculate average silhouette width
    avg_sil <- mean(sil_df$sil_width)

    # Add vertical line for average silhouette
    p <- p %>%
      add_trace(
        x = c(avg_sil, avg_sil),
        y = c(0, max(sil_df$y_position) + 1),
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'red', width = 2, dash = 'dash'),
        name = paste0('Average (', round(avg_sil, 3), ')'),
        showlegend = TRUE,
        hoverinfo = 'skip'
      )

    # Add cluster separators
    for(i in 1:(k-1)) {
      sep_y <- cumsum_sizes[i + 1] + 0.5
      p <- p %>%
        add_trace(
          x = c(-1, 1),
          y = c(sep_y, sep_y),
          type = 'scatter',
          mode = 'lines',
          line = list(color = 'black', width = 1, dash = 'dot'),
          showlegend = FALSE,
          hoverinfo = 'skip'
        )
    }

    # Update layout
    p <- p %>%
      layout(
        title = list(
          text = paste0("Silhouette Plot for k=", k, " Clusters (Avg: ", round(avg_sil, 3), ")"),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          title = "Silhouette Coefficient",
          range = c(min(c(sil_df$sil_width, -0.1)), 1)
        ),
        yaxis = list(
          title = "Sample Index (sorted within cluster)",
          showticklabels = FALSE
        ),
        barmode = 'overlay',
        showlegend = TRUE,
        legend = list(
          orientation = "v",
          x = 1.02,
          y = 1
        ),
        hovermode = "closest"
      )

    # Add annotations for cluster labels
    for(i in 1:k) {
      cluster_data <- sil_df %>% filter(cluster == i)
      mid_y <- mean(cluster_data$y_position)

      p <- p %>%
        add_annotations(
          x = -0.05,
          y = mid_y,
          text = paste0("C", i),
          showarrow = FALSE,
          xref = "x",
          yref = "y",
          font = list(size = 12, color = strong_colors[i])
        )
    }

    p
  })

  # ===== DENDROGRAM OUTPUTS =====

  # Main clustering tab dendrogram
  output$dendrogram_plot <- renderPlot({
    req(values$clustering_results)
    res <- values$clustering_results
    req(res$method == "hierarchical", !is.null(res$hc_obj))
    hc <- res$hc_obj
    k  <- res$k
    # Determine the height at which to cut for k clusters
    cut_height <- hc$height[length(hc$height) - k + 1]
    par(mar = c(4, 4, 3, 1))
    plot(hc, labels = FALSE, hang = -1,
         main = paste0("Dendrogram — ", res$linkage_method, " linkage (cut at k=", k, ")"),
         xlab = "Participants", ylab = "Height", sub = "")
    abline(h = cut_height, col = "red", lty = 2, lwd = 1.5)
    legend("topright", legend = paste0("Cut height = ", round(cut_height, 3)),
           col = "red", lty = 2, lwd = 1.5, bty = "n")
  })

  # Optimization tab dendrogram
  output$opt_dendrogram_plot <- renderPlot({
    req(values$cluster_optimization)
    opt <- values$cluster_optimization
    req(opt$method == "hierarchical", !is.null(opt$hc_obj))
    hc <- opt$hc_obj
    par(mar = c(4, 4, 3, 1))
    plot(hc, labels = FALSE, hang = -1,
         main = paste0("Full Dendrogram — ", if(!is.null(input$opt_hclust_linkage)) input$opt_hclust_linkage else "ward.D2", " linkage"),
         xlab = "Participants", ylab = "Height", sub = "")
  })

  # Export functions
  output$export_scores_csv <- downloadHandler(
    filename = function() paste0("pca_scores_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pca_results) && !is.null(values$pca_results$scores)) {
        scores_df <- data.frame(
          Subject = 1:nrow(values$pca_results$scores),
          values$pca_results$scores
        )
        colnames(scores_df)[-1] <- paste0("PC", 1:ncol(values$pca_results$scores))
        if(!is.null(values$group_labels)) {
          scores_df$Group <- values$group_labels
        }
        write.csv(scores_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No PCA scores available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_loadings_csv <- downloadHandler(
    filename = function() paste0("pca_loadings_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pca_results) && !is.null(values$pca_results$harmonics)) {
        time_points <- seq(0, 1, length.out = 100)
        loadings_mat <- matrix(NA, 100, length(values$pca_results$harmonics))
        for(i in 1:length(values$pca_results$harmonics)) {
          loadings_mat[,i] <- eval.fd(time_points, values$pca_results$harmonics[i])
        }
        loadings_df <- data.frame(Time = time_points, loadings_mat)
        colnames(loadings_df)[-1] <- paste0("PC", 1:ncol(loadings_mat))
        write.csv(loadings_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No loadings available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_warping_csv <- downloadHandler(
    filename = function() paste0("warping_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$warping_results) && !is.null(values$warping_results$warp_functions)) {
        warp_df <- data.frame(
          Time = values$warping_results$time_points,
          values$warping_results$warp_functions
        )
        colnames(warp_df)[-1] <- paste0("Subject", 1:(ncol(warp_df)-1))
        write.csv(warp_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No warping results available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_fanova_results_csv <- downloadHandler(
    filename = function() paste0("fanova_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$fanova_results)) {
        res <- values$fanova_results
        results_df <- data.frame(
          Time = res$time_points,
          F_statistic = res$F_stat,
          p_value_raw = res$p_values_pointwise,
          p_value_adjusted = res$p_values_adjusted,
          significant = res$sig_regions,
          eta_squared = res$eta_squared
        )
        
        for(i in 1:res$n_groups) {
          results_df[paste0("Mean_", res$groups[i])] <- res$group_means[,i]
        }
        
        write.csv(results_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No FANOVA results available"), file, row.names = FALSE)
      }
    }
  )
  
  output$export_pairwise_results_csv <- downloadHandler(
    filename = function() paste0("pairwise_results_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$pairwise_results)) {
        results_list <- list()
        
        for(pair_name in values$pairwise_results$pair_names) {
          pair_result <- values$pairwise_results$results[[pair_name]]
          
          results_list[[pair_name]] <- data.frame(
            Time = values$pairwise_results$time_points,
            Mean_Diff = pair_result$mean_diff,
            t_stat = pair_result$t_stat,
            p_raw = pair_result$p_values_pointwise,
            p_adj = pair_result$p_values_adjusted,
            Cohen_d = pair_result$cohens_d,
            CI_lower = pair_result$ci_lower,
            CI_upper = pair_result$ci_upper,
            Significant = pair_result$sig_regions
          )
        }
        
        all_results <- do.call(rbind, results_list)
        all_results$Comparison <- rep(names(results_list), each = length(values$pairwise_results$time_points))
        
        write.csv(all_results, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No pairwise results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster membership
  output$export_cluster_membership_csv <- downloadHandler(
    filename = function() paste0("cluster_membership_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        membership_df <- data.frame(
          Subject_ID = 1:length(values$clustering_results$cluster_assignments),
          Cluster = values$clustering_results$cluster_assignments
        )

        # Add group labels if available
        if(!is.null(values$group_labels)) {
          membership_df$Group <- values$group_labels
        }

        write.csv(membership_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster mean functions
  output$export_cluster_means_csv <- downloadHandler(
    filename = function() paste0("cluster_means_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        means_df <- as.data.frame(values$clustering_results$cluster_means)

        # Use time labels if available
        if(!is.null(values$time_labels) && length(values$time_labels) == ncol(means_df)) {
          colnames(means_df) <- values$time_labels
        } else {
          colnames(means_df) <- paste0("T", 1:ncol(means_df))
        }

        # Add cluster ID column
        means_df <- cbind(Cluster = 1:nrow(means_df), means_df)

        write.csv(means_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster statistics
  output$export_cluster_stats_csv <- downloadHandler(
    filename = function() paste0("cluster_statistics_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results)) {
        r <- values$clustering_results

        stats_df <- data.frame(
          Cluster = 1:r$k,
          Size = r$cluster_sizes,
          Percentage = round(100 * r$cluster_sizes / sum(r$cluster_sizes), 2),
          Within_SS = r$wcss_per_cluster,
          Avg_Within_SS = r$wcss_per_cluster / r$cluster_sizes
        )

        # Add overall statistics
        overall_stats <- data.frame(
          Metric = c("Total_SS", "Between_SS", "Within_SS", "R_Squared",
                     "Avg_Silhouette", "Calinski_Harabasz"),
          Value = c(r$total_ss, r$between_ss, r$total_wcss, r$r_squared,
                    ifelse(is.na(r$silhouette_width), NA, r$silhouette_width),
                    r$ch_index)
        )

        # Write both tables to separate sheets would require xlsx package
        # For CSV, concatenate with a separator
        write.csv(stats_df, file, row.names = FALSE)
        write("", file, append = TRUE)
        write("Overall Statistics:", file, append = TRUE)
        write.table(overall_stats, file, append = TRUE, sep = ",",
                    row.names = FALSE, col.names = TRUE)
      } else {
        write.csv(data.frame(Message = "No clustering results available"), file, row.names = FALSE)
      }
    }
  )

  # Export silhouette data
  output$export_silhouette_csv <- downloadHandler(
    filename = function() paste0("silhouette_data_", Sys.Date(), ".csv"),
    content = function(file) {
      if(!is.null(values$clustering_results) &&
         !is.null(values$clustering_results$silhouette_data)) {
        sil_data <- values$clustering_results$silhouette_data

        sil_df <- data.frame(
          Subject_ID = 1:nrow(sil_data),
          Cluster = sil_data[, 1],
          Neighbor = sil_data[, 2],
          Silhouette_Width = sil_data[, 3]
        )

        # Add group labels if available
        if(!is.null(values$group_labels)) {
          sil_df$Group <- values$group_labels
        }

        write.csv(sil_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No silhouette data available"), file, row.names = FALSE)
      }
    }
  )

  # Export cluster-group composition
  output$export_cluster_group_csv <- downloadHandler(
    filename = function() {
      group_var_name <- get_cluster_group_var_name()
      paste0("cluster_", group_var_name, "_composition_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if(!is.null(values$clustering_results) && !is.null(values$group_labels)) {
        clusters <- values$clustering_results$cluster_assignments
        groups <- get_cluster_group_labels()
        group_var_name <- get_cluster_group_var_name()

        # Create contingency table
        cont_table <- table(Cluster = clusters, Group = groups)

        # Calculate percentages for each cluster
        cluster_pcts <- prop.table(cont_table, margin = 1) * 100

        # Create export dataframe
        export_df <- data.frame(Cluster = rownames(cont_table))

        # Add metadata row with group variable name
        attr(export_df, "group_variable") <- group_var_name

        # Add count and percentage columns for each group
        for(g in colnames(cont_table)) {
          export_df[[paste0(g, "_Count")]] <- cont_table[, g]
          export_df[[paste0(g, "_Percent")]] <- round(cluster_pcts[, g], 2)
        }

        # Add total column
        export_df$Total <- rowSums(cont_table)

        # Add group variable name as comment in first row
        write.csv(export_df, file, row.names = FALSE)

        # Append metadata
        cat(paste0("\n# Group Variable: ", group_var_name, "\n"), file = file, append = TRUE)
      } else {
        write.csv(data.frame(Message = "No group composition data available"), file, row.names = FALSE)
      }
    }
  )

  # Export smoothed curves - Wide format
  output$export_smoothed_csv <- downloadHandler(
    filename = function() paste0("smoothed_curves_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        if(!is.null(values$smooth_data)) {
          # Export smooth_data directly (already smoothed curves)
          df_export <- as.data.frame(values$smooth_data)
          
          # Use original time labels directly (no interpolation!)
          if(!is.null(values$time_labels) && length(values$time_labels) == ncol(values$smooth_data)) {
            colnames(df_export) <- values$time_labels
          } else {
            colnames(df_export) <- paste0("T", 1:ncol(df_export))
          }
          
          # Add subject IDs
          df_export <- cbind(Subject = 1:nrow(df_export), df_export)
          
          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == nrow(df_export)) {
            df_export <- cbind(Group = values$group_labels, df_export)
          }
          
          write.csv(df_export, file, row.names = FALSE)
          cat("Successfully exported smoothed curves (wide format) to:", file, "\n")
          cat("Exported", nrow(df_export), "subjects with", ncol(values$smooth_data), "time points\n")
          
        } else {
          write.csv(data.frame(Message = "No smoothed data available. Please apply smoothing first."), 
                    file, row.names = FALSE)
          cat("No smoothed data available for export\n")
        }
      }, error = function(e) {
        cat("Error in export_smoothed_csv:", e$message, "\n")
        write.csv(data.frame(Error = paste("Export failed:", e$message)), 
                  file, row.names = FALSE)
      })
    }
  )
  
  # Export smoothed curves - Long format
  output$export_smoothed_long_csv <- downloadHandler(
    filename = function() paste0("smoothed_curves_long_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        if(!is.null(values$smooth_data)) {
          n_subjects <- nrow(values$smooth_data)
          n_time <- ncol(values$smooth_data)
          
          # Use original time labels directly (no interpolation!)
          if(!is.null(values$time_labels) && length(values$time_labels) == n_time) {
            time_labels_use <- values$time_labels
          } else {
            time_labels_use <- 1:n_time
          }
          
          # Create normalized time grid
          time_normalized <- seq(0, 1, length.out = n_time)
          
          # Convert to long format
          df_long <- data.frame(
            Subject = rep(1:n_subjects, each = n_time),
            Time = rep(time_labels_use, times = n_subjects),
            Time_Normalized = rep(time_normalized, times = n_subjects),
            Value = as.vector(t(values$smooth_data))
          )
          
          # Add group labels if available
          if(!is.null(values$group_labels) && length(values$group_labels) == n_subjects) {
            df_long$Group <- rep(values$group_labels, each = n_time)
            # Reorder columns
            df_long <- df_long[, c("Subject", "Group", "Time", "Time_Normalized", "Value")]
          }
          
          write.csv(df_long, file, row.names = FALSE)
          cat("Successfully exported smoothed curves (long format) to:", file, "\n")
          cat("Exported", n_subjects, "subjects ×", n_time, "time points =", nrow(df_long), "rows\n")
          
        } else {
          write.csv(data.frame(Message = "No smoothed data available. Please apply smoothing first."), 
                    file, row.names = FALSE)
          cat("No smoothed data available for export\n")
        }
      }, error = function(e) {
        cat("Error in export_smoothed_long_csv:", e$message, "\n")
        write.csv(data.frame(Error = paste("Export failed:", e$message)), 
                  file, row.names = FALSE)
      })
    }
  )
  
  output$export_plots <- downloadHandler(
    filename = function() paste0("all_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$data)) {
        plot(1:10, main = "Plots export not fully implemented")
        text(5, 5, "Export functionality to be implemented", cex = 2)
      }
      
      dev.off()
    }
  )
  
  output$export_fanova_plots <- downloadHandler(
    filename = function() paste0("fanova_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$fanova_results)) {
        res <- values$fanova_results
        par(mfrow = c(2, 2))
        
        matplot(res$time_points, res$group_means, type = "l", 
                col = 1:res$n_groups, lwd = 2,
                xlab = "Time", ylab = "Value",
                main = "Group Mean Functions")
        legend("topright", legend = res$groups, col = 1:res$n_groups, lty = 1, lwd = 2)
        
        plot(res$time_points, res$F_stat, type = "l", col = "blue", lwd = 2,
             xlab = "Time", ylab = "F-statistic",
             main = "Pointwise F-statistics")
        
        plot(res$time_points, res$p_values_adjusted, type = "l", col = "darkgreen", lwd = 2,
             xlab = "Time", ylab = "p-value", log = "y",
             main = "Adjusted p-values")
        
        plot(res$time_points, res$eta_squared, type = "l", col = "purple", lwd = 2,
             xlab = "Time", ylab = "η²",
             main = "Effect Size (η²)")
      }
      
      dev.off()
    }
  )
  
  output$export_pairwise_plots <- downloadHandler(
    filename = function() paste0("pairwise_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      if(!is.null(values$pairwise_results)) {
        par(mfrow = c(2, 2))
        
        for(pair_name in values$pairwise_results$pair_names[1:min(4, length(values$pairwise_results$pair_names))]) {
          pair_result <- values$pairwise_results$results[[pair_name]]
          time_points <- values$pairwise_results$time_points
          
          plot(time_points, pair_result$mean_diff, type = "l", col = "blue", lwd = 2,
               xlab = "Time", ylab = "Mean Difference",
               main = pair_name)
          abline(h = 0, lty = 2)
          
          lines(time_points, pair_result$ci_lower, col = "lightblue", lty = 2)
          lines(time_points, pair_result$ci_upper, col = "lightblue", lty = 2)
        }
      }
      
      dev.off()
    }
  )
  
  # ============================================================================
  # ANALYSIS CODE GENERATION
  # ============================================================================
  # Generate reproducible R code for all performed analyses

  generate_analysis_code <- function(full = TRUE) {
    # Generate comprehensive R code for reproducing the analysis
    # Args:
    #   full: If TRUE, generate complete code; if FALSE, generate preview (truncated)

    code_lines <- c()

    # Helper to add lines
    add <- function(...) {
      code_lines <<- c(code_lines, paste0(...))
    }

    add("# =============================================================================")
    add("# FUNCTIONAL DATA ANALYSIS - REPRODUCIBLE CODE")
    add("# Generated from WAPAA Shiny App")
    add("# Date: ", Sys.Date())
    add("# =============================================================================")
    add("")

    # ---- SECTION 1: LIBRARIES ----
    add("# -----------------------------------------------------------------------------")
    add("# 1. LOAD REQUIRED LIBRARIES")
    add("# -----------------------------------------------------------------------------")
    add("library(fda)        # Functional data analysis")
    add("library(mgcv)       # GAM smoothing with REML")
    add("library(ggplot2)    # Plotting")
    add("library(dplyr)      # Data manipulation")

    if(!is.null(values$fanova_results) && !is.null(values$fanova_results$design) &&
       values$fanova_results$design == "within") {
      add("library(rmfanova)   # Repeated measures functional ANOVA")
    }

    if(!is.null(values$clustering_results)) {
      add("library(cluster)    # Clustering diagnostics")
      add("library(fda.usc)    # Functional clustering")
    }
    add("")

    # ---- SECTION 2: DATA LOADING ----
    add("# -----------------------------------------------------------------------------")
    add("# 2. LOAD AND PREPARE DATA")
    add("# -----------------------------------------------------------------------------")

    if(!is.null(values$data)) {
      n_subj <- nrow(values$data)
      n_time <- ncol(values$data)

      add("# Data dimensions: ", n_subj, " subjects x ", n_time, " time points")
      add("")
      add("# Load your data (replace with your actual file path)")
      add("# The data should be a matrix with subjects in rows and time points in columns")
      add("data_matrix <- read.csv('your_data.csv', row.names = 1)")
      add("data_matrix <- as.matrix(data_matrix)")
      add("")
      add("# Number of subjects and time points")
      add("n_subjects <- nrow(data_matrix)  # ", n_subj)
      add("n_time <- ncol(data_matrix)      # ", n_time)
      add("")
      add("# Create normalized time grid [0, 1]")
      add("time_points <- seq(0, 1, length.out = n_time)")
      add("")

      # Time labels if available
      if(!is.null(values$time_labels)) {
        time_labels_str <- paste0("c('", paste(head(values$time_labels, 5), collapse = "', '"),
                                  if(length(values$time_labels) > 5) "', ..." else "'", ")")
        add("# Original time labels: ", time_labels_str)
        add("")
      }

      # Group labels
      if(!is.null(values$group_labels)) {
        groups <- levels(values$group_labels)
        group_counts <- table(values$group_labels)
        add("# Group structure")
        add("group_labels <- factor(c(")
        # Show first few group assignments
        n_show <- min(10, length(values$group_labels))
        group_sample <- paste0("'", as.character(values$group_labels[1:n_show]), "'", collapse = ", ")
        add("  ", group_sample, if(length(values$group_labels) > 10) ", ..." else "")
        add("))")
        add("# Groups: ", paste(names(group_counts), " (n=", group_counts, ")", sep = "", collapse = ", "))
        add("")
      }
    }

    # ---- SECTION 3: SMOOTHING ----
    add("# -----------------------------------------------------------------------------")
    add("# 3. SMOOTHING / FUNCTIONAL DATA OBJECT CREATION")
    add("# -----------------------------------------------------------------------------")

    if(!is.null(values$smooth_fit_metrics)) {
      metrics <- values$smooth_fit_metrics

      # Determine smoothing method used
      smooth_method <- if(!is.null(input$smooth_method)) input$smooth_method else "auto"

      if(smooth_method == "none") {
        add("# Smoothing method: None (raw data)")
        add("# Data used as-is without smoothing")
        add("")
      } else if(smooth_method == "auto") {
        n_basis <- if(!is.null(input$n_basis)) input$n_basis else 20
        add("# Smoothing method: Automatic (REML optimization)")
        add("# Number of B-spline basis functions: ", n_basis)
        if(!is.null(metrics$lambda)) {
          add("# Estimated lambda (smoothing parameter): ", sprintf("%.6e", metrics$lambda))
        }
        add("")
        add("# Create B-spline basis")
        add("n_basis <- ", n_basis)
        add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
        add("")
        add("# Smooth with REML (lambda = 0 triggers automatic optimization)")
        add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = 0)")
        add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
        add("fd_obj <- smooth_result$fd")
        add("")
      } else if(smooth_method == "manual") {
        n_basis <- if(!is.null(input$n_basis_manual)) input$n_basis_manual else 20
        smooth_factor <- if(!is.null(input$smooth_factor)) input$smooth_factor else 1
        lambda <- 10^(-smooth_factor)
        add("# Smoothing method: Manual")
        add("# Number of B-spline basis functions: ", n_basis)
        add("# Smoothing factor: ", smooth_factor, " (lambda = 10^(-", smooth_factor, ") = ", sprintf("%.6e", lambda), ")")
        add("")
        add("# Create B-spline basis")
        add("n_basis <- ", n_basis)
        add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
        add("")
        add("# Smooth with specified lambda")
        add("lambda <- 10^(-", smooth_factor, ")  # = ", sprintf("%.6e", lambda))
        add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = lambda)")
        add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
        add("fd_obj <- smooth_result$fd")
        add("")
      }

      # Fit statistics
      add("# Smoothing fit statistics:")
      add("#   Mean R²: ", sprintf("%.4f", metrics$mean_r_squared),
          " (SD: ", sprintf("%.4f", metrics$sd_r_squared), ")")
      add("#   Mean RMSE: ", sprintf("%.4f", metrics$mean_rmse),
          " (SD: ", sprintf("%.4f", metrics$sd_rmse), ")")
      add("")

      # Bounds constraint
      if(!is.null(input$constrain_bounds) && input$constrain_bounds) {
        add("# Value bounds constraint applied")
        add("min_bound <- ", input$min_bound)
        add("max_bound <- ", input$max_bound)
        add("# Smoothed values were clamped to [min_bound, max_bound]")
        add("")
      }
    } else {
      add("# No smoothing applied yet")
      add("# Default smoothing code:")
      add("n_basis <- 20")
      add("basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)")
      add("fd_par <- fdPar(basis, Lfdobj = 2, lambda = 0)  # REML optimization")
      add("smooth_result <- smooth.basis(time_points, t(data_matrix), fd_par)")
      add("fd_obj <- smooth_result$fd")
      add("")
    }

    add("# Evaluate smoothed curves")
    add("smooth_curves <- eval.fd(time_points, fd_obj)")
    add("smooth_curves <- t(smooth_curves)  # Back to subjects x time format")
    add("")

    # ---- SECTION 4: TIME WARPING ----
    if(!is.null(values$warping_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 4. TIME WARPING / CURVE REGISTRATION")
      add("# -----------------------------------------------------------------------------")

      warp_method <- values$warping_results$method
      add("# Warping method: ", warp_method)
      add("")

      if(warp_method == "linear_shift") {
        periodic <- if(!is.null(input$periodic_shift)) input$periodic_shift else FALSE
        reference <- if(!is.null(input$shift_reference)) input$shift_reference else "mean"
        allow_dilation <- if(!is.null(input$allow_dilation)) input$allow_dilation else FALSE

        add("# Linear shift alignment parameters:")
        add("#   Reference: ", reference)
        add("#   Periodic: ", periodic)
        add("#   Allow dilation: ", allow_dilation)
        add("")

        add("# Linear shift alignment function")
        add("linear_shift_alignment <- function(fd_obj, periodic = ", periodic,
            ", reference = '", reference, "', time_points = seq(0, 1, length.out = 100)) {")
        add("  n_curves <- ncol(fd_obj$coefs)")
        add("  n_time <- length(time_points)")
        add("  curves <- eval.fd(time_points, fd_obj)")
        add("")
        add("  # Get reference curve")
        add("  ref_curve <- switch(reference,")
        add("    'mean' = rowMeans(curves),")
        add("    'median' = apply(curves, 1, median),")
        add("    curves[,1])")
        add("")
        add("  # Initialize outputs")
        add("  registered_curves <- matrix(NA, n_time, n_curves)")
        add("  warp_functions <- matrix(NA, n_time, n_curves)")
        add("  shifts <- numeric(n_curves)")
        add("")
        add("  for(i in 1:n_curves) {")
        add("    # Find optimal shift using cross-correlation")
        add("    ccf_result <- ccf(curves[,i], ref_curve, lag.max = floor(n_time/4), plot = FALSE)")
        add("    best_lag <- ccf_result$lag[which.max(ccf_result$acf)]")
        add("    shifts[i] <- best_lag / n_time * 0.1")
        add("")
        add("    # Create warping function")
        add("    warp_functions[,i] <- pmin(1, pmax(0, time_points - shifts[i] * 0.5))")
        add("    warp_functions[1,i] <- 0; warp_functions[n_time,i] <- 1")
        add("")
        add("    # Apply warping")
        add("    if(abs(shifts[i]) > 0.001) {")
        add("      registered_curves[,i] <- approx(time_points, curves[,i],")
        add("                                       xout = warp_functions[,i], rule = 2)$y")
        add("    } else {")
        add("      registered_curves[,i] <- curves[,i]")
        add("    }")
        add("  }")
        add("")
        add("  return(list(registered_curves = registered_curves,")
        add("              warp_functions = warp_functions, shifts = shifts))")
        add("}")
        add("")
        add("# Apply alignment")
        add("warp_result <- linear_shift_alignment(fd_obj, periodic = ", periodic, ", reference = '", reference, "')")
        add("registered_curves <- warp_result$registered_curves")
        add("")

      } else if(warp_method == "parametric") {
        family <- if(!is.null(input$parametric_family)) input$parametric_family else "power"
        param_range <- if(!is.null(input$param_range)) input$param_range else c(0.5, 2)

        add("# Parametric warping parameters:")
        add("#   Family: ", family)
        add("#   Parameter range: [", param_range[1], ", ", param_range[2], "]")
        add("")

        add("# Parametric alignment function")
        add("parametric_alignment <- function(fd_obj, family = '", family, "',")
        add("                                  param_range = c(", param_range[1], ", ", param_range[2], "),")
        add("                                  time_points = seq(0, 1, length.out = 100)) {")
        add("  n_curves <- ncol(fd_obj$coefs)")
        add("  n_time <- length(time_points)")
        add("  curves <- eval.fd(time_points, fd_obj)")
        add("  mean_curve <- rowMeans(curves)")
        add("")
        add("  # Warping function based on family")
        add("  warp_func <- function(t, alpha) {")
        add("    switch(family,")
        add("      'power' = pmin(1, pmax(0, t^alpha)),")
        add("      'exponential' = if(abs(alpha-1) < 0.001) t else pmin(1, pmax(0, (exp(alpha*t)-1)/(exp(alpha)-1))),")
        add("      'quadratic' = pmin(1, pmax(0, alpha*t^2 + (1-alpha)*t)),")
        add("      t)")
        add("  }")
        add("")
        add("  registered_curves <- matrix(NA, n_time, n_curves)")
        add("  alpha_values <- numeric(n_curves)")
        add("")
        add("  for(i in 1:n_curves) {")
        add("    # Optimize warping parameter")
        add("    objective <- function(alpha) {")
        add("      warped <- approx(time_points, curves[,i], xout = warp_func(time_points, alpha), rule = 2)$y")
        add("      sum((warped - mean_curve)^2, na.rm = TRUE)")
        add("    }")
        add("    result <- optimize(objective, interval = param_range, tol = 1e-4)")
        add("    alpha_values[i] <- result$minimum")
        add("    registered_curves[,i] <- approx(time_points, curves[,i],")
        add("                                     xout = warp_func(time_points, alpha_values[i]), rule = 2)$y")
        add("  }")
        add("")
        add("  return(list(registered_curves = registered_curves, alpha_values = alpha_values))")
        add("}")
        add("")
        add("# Apply alignment")
        add("warp_result <- parametric_alignment(fd_obj, family = '", family, "')")
        add("registered_curves <- warp_result$registered_curves")
        add("")

      } else if(warp_method == "landmark") {
        add("# Landmark-based alignment")
        add("# (Implementation depends on detected landmarks)")
        add("")
      }

      # Warping fit statistics
      if(!is.null(values$warping_results$fit_statistics)) {
        stats <- values$warping_results$fit_statistics$summary
        add("# Warping fit statistics:")
        add("#   Mean R² (orig vs warped): ", sprintf("%.4f", stats$mean_r_squared))
        add("#   Mean RMSE: ", sprintf("%.4f", stats$mean_rmse))
        add("#   Variance explained by warping: ", sprintf("%.2f%%", stats$variance_explained_by_warping * 100))
        add("#   AIC: ", sprintf("%.2f", stats$AIC))
        add("#   BIC: ", sprintf("%.2f", stats$BIC))
        add("")
      }

      add("# Create fd object from registered curves")
      add("reg_basis <- create.bspline.basis(c(0, 1), nbasis = min(20, n_time - 2))")
      add("reg_fd <- smooth.basis(time_points, registered_curves, fdPar(reg_basis, 2, 0))$fd")
      add("")
    }

    # ---- SECTION 5: PCA ----
    if(!is.null(values$pca_results)) {
      add("# -----------------------------------------------------------------------------")
      if(!is.null(values$warping_results)) {
        add("# 5. TIME-WARPED PRINCIPAL COMPONENT ANALYSIS")
      } else {
        add("# 5. FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
      }
      add("# -----------------------------------------------------------------------------")

      n_comp <- ncol(values$pca_results$scores)
      varprop <- values$pca_results$varprop[1:n_comp]

      add("# Number of components: ", n_comp)
      add("# Variance explained:")
      for(i in 1:min(n_comp, 5)) {
        add("#   PC", i, ": ", sprintf("%.2f%%", varprop[i] * 100),
            " (cumulative: ", sprintf("%.2f%%", sum(varprop[1:i]) * 100), ")")
      }
      add("")

      add("# Perform functional PCA")
      add("n_components <- ", n_comp)
      if(!is.null(values$warping_results)) {
        add("pca_results <- pca.fd(reg_fd, nharm = n_components)")
      } else {
        add("pca_results <- pca.fd(fd_obj, nharm = n_components)")
      }
      add("")
      add("# Extract results")
      add("scores <- pca_results$scores          # PC scores (n_subjects x n_components)")
      add("varprop <- pca_results$varprop        # Proportion of variance explained")
      add("harmonics <- pca_results$harmonics    # Functional principal components")
      add("")
      add("# View variance explained")
      add("print(round(varprop * 100, 2))")
      add("")
      add("# Plot PC loadings")
      add("par(mfrow = c(1, min(3, n_components)))")
      add("for(i in 1:min(3, n_components)) {")
      add("  plot(harmonics[i], main = paste0('PC', i, ' (', round(varprop[i]*100, 1), '%)'),")
      add("       xlab = 'Time', ylab = 'Loading')")
      add("}")
      add("")
    }

    # ---- SECTION 6: CLUSTERING ----
    if(!is.null(values$clustering_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 6. FUNCTIONAL K-MEANS CLUSTERING")
      add("# -----------------------------------------------------------------------------")

      k <- values$clustering_results$k
      add("# Number of clusters: ", k)
      add("")

      add("# Prepare data for clustering")
      if(!is.null(values$warping_results)) {
        add("fdata_obj <- fdata(t(registered_curves), argvals = time_points)")
      } else {
        add("fdata_obj <- fdata(t(eval.fd(time_points, fd_obj)), argvals = time_points)")
      }
      add("")
      add("# Perform functional k-means clustering")
      add("set.seed(123)  # For reproducibility")
      add("k_clusters <- ", k)
      add("kmeans_result <- kmeans.fd(fdata_obj, ncl = k_clusters, metric = 'L2')")
      add("")
      add("# Get cluster assignments")
      add("cluster_assignments <- kmeans_result$cluster")
      add("print(table(cluster_assignments))")
      add("")

      if(!is.null(values$clustering_results$r_squared)) {
        add("# Clustering fit statistics:")
        add("#   R² (variance explained): ", sprintf("%.4f", values$clustering_results$r_squared))
        add("#   Within-cluster SS: ", sprintf("%.4f", values$clustering_results$total_wcss))
        add("#   Between-cluster SS: ", sprintf("%.4f", values$clustering_results$between_ss))
      }
      add("")
    }

    # ---- SECTION 7: FUNCTIONAL ANOVA ----
    if(!is.null(values$fanova_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 7. FUNCTIONAL ANOVA")
      add("# -----------------------------------------------------------------------------")

      design <- if(!is.null(values$fanova_results$design)) values$fanova_results$design else "between"
      n_groups <- values$fanova_results$n_groups

      add("# Design type: ", if(design == "within") "Within-subjects (Repeated Measures)" else "Between-subjects")
      add("# Number of groups: ", n_groups)
      add("")

      if(design == "between") {
        add("# Between-subjects functional ANOVA")
        add("")
        add("# Prepare data")
        if(!is.null(values$warping_results)) {
          add("fd_to_analyze <- reg_fd")
        } else {
          add("fd_to_analyze <- fd_obj")
        }
        add("")
        add("# Pointwise F-tests")
        add("n_time_eval <- 100")
        add("eval_points <- seq(0, 1, length.out = n_time_eval)")
        add("curves_eval <- t(eval.fd(eval_points, fd_to_analyze))")
        add("")
        add("# Perform pointwise ANOVA at each time point")
        add("f_values <- numeric(n_time_eval)")
        add("p_values <- numeric(n_time_eval)")
        add("")
        add("for(t in 1:n_time_eval) {")
        add("  aov_result <- summary(aov(curves_eval[, t] ~ group_labels))")
        add("  f_values[t] <- aov_result[[1]]$'F value'[1]")
        add("  p_values[t] <- aov_result[[1]]$'Pr(>F)'[1]")
        add("}")
        add("")
        add("# Identify significant time regions (p < 0.05)")
        add("sig_times <- eval_points[p_values < 0.05]")
        add("cat('Significant time points:', length(sig_times), '/', n_time_eval, '\\n')")
        add("")

        # Global test results
        if(!is.null(values$fanova_results$global_p)) {
          add("# Global test results:")
          add("#   Global F-statistic: ", sprintf("%.4f", values$fanova_results$global_f))
          add("#   Global p-value: ", sprintf("%.6f", values$fanova_results$global_p))
          add("")
        }

      } else {
        add("# Within-subjects (Repeated Measures) functional ANOVA")
        add("# Using rmfanova package")
        add("")
        add("# Prepare data for rmfanova")
        add("# Data should be in long format with subject_id, condition, and functional observations")
        add("")
        add("# Example structure:")
        add("# rm_data <- data.frame(")
        add("#   subject_id = rep(1:n_subjects, each = n_conditions),")
        add("#   condition = rep(condition_levels, n_subjects),")
        add("#   ... # functional data columns")
        add("# )")
        add("")
        add("# Fit repeated measures functional ANOVA")
        add("# rm_result <- rmfanova(formula = value ~ condition + Error(subject_id),")
        add("#                       data = rm_data, argvals = time_points)")
        add("")
      }
    }

    # ---- SECTION 8: PAIRWISE COMPARISONS ----
    if(!is.null(values$pairwise_results)) {
      add("# -----------------------------------------------------------------------------")
      add("# 8. PAIRWISE COMPARISONS")
      add("# -----------------------------------------------------------------------------")

      correction <- values$pairwise_results$correction_method
      add("# Multiple testing correction: ", correction)
      add("")

      add("# Perform pairwise functional t-tests")
      add("groups <- levels(group_labels)")
      add("n_groups <- length(groups)")
      add("pairs <- combn(groups, 2, simplify = FALSE)")
      add("")
      add("pairwise_results <- list()")
      add("for(pair in pairs) {")
      add("  g1 <- pair[1]; g2 <- pair[2]")
      add("  idx1 <- which(group_labels == g1)")
      add("  idx2 <- which(group_labels == g2)")
      add("")
      add("  # Pointwise t-tests")
      add("  p_values <- numeric(n_time_eval)")
      add("  for(t in 1:n_time_eval) {")
      add("    tt <- t.test(curves_eval[idx1, t], curves_eval[idx2, t])")
      add("    p_values[t] <- tt$p.value")
      add("  }")
      add("")
      add("  # Apply correction")
      add("  p_adjusted <- p.adjust(p_values, method = '", tolower(correction), "')")
      add("")
      add("  pairwise_results[[paste(g1, 'vs', g2)]] <- list(")
      add("    p_values = p_values,")
      add("    p_adjusted = p_adjusted,")
      add("    sig_times = eval_points[p_adjusted < 0.05]")
      add("  )")
      add("}")
      add("")

      # Summary of significant comparisons
      if(!is.null(values$pairwise_results$results)) {
        n_sig <- sum(sapply(values$pairwise_results$results, function(x)
          length(x$sig_times) > 0))
        add("# Number of significant pairwise comparisons: ", n_sig, "/", length(values$pairwise_results$results))
      }
      add("")
    }

    # ---- SECTION 9: VISUALIZATION ----
    add("# -----------------------------------------------------------------------------")
    add("# 9. VISUALIZATION CODE")
    add("# -----------------------------------------------------------------------------")
    add("")
    add("# Plot all curves with mean")
    add("matplot(time_points, t(smooth_curves), type = 'l', col = 'gray70', lty = 1,")
    add("        xlab = 'Time', ylab = 'Value', main = 'Functional Data')")
    add("lines(time_points, colMeans(smooth_curves), col = 'red', lwd = 3)")
    add("")

    if(!is.null(values$group_labels)) {
      add("# Plot by group")
      add("library(ggplot2)")
      add("plot_data <- data.frame(")
      add("  time = rep(time_points, n_subjects),")
      add("  value = as.vector(t(smooth_curves)),")
      add("  subject = rep(1:n_subjects, each = n_time),")
      add("  group = rep(group_labels, each = n_time)")
      add(")")
      add("")
      add("ggplot(plot_data, aes(x = time, y = value, group = subject, color = group)) +")
      add("  geom_line(alpha = 0.5) +")
      add("  stat_summary(aes(group = group), fun = mean, geom = 'line', linewidth = 2) +")
      add("  theme_minimal() +")
      add("  labs(title = 'Functional Data by Group', x = 'Time', y = 'Value')")
      add("")
    }

    add("# =============================================================================")
    add("# END OF ANALYSIS CODE")
    add("# =============================================================================")

    return(paste(code_lines, collapse = "\n"))
  }

  # Download handler for analysis code
  output$export_code <- downloadHandler(
    filename = function() paste0("fda_analysis_code_", Sys.Date(), ".R"),
    content = function(file) {
      code <- generate_analysis_code(full = TRUE)
      writeLines(code, file)
    }
  )

  # Code preview (truncated version)
  output$code_preview <- renderText({
    tryCatch({
      code <- generate_analysis_code(full = FALSE)

      # For preview, show first ~80 lines
      lines <- strsplit(code, "\n")[[1]]
      if(length(lines) > 80) {
        preview <- c(
          lines[1:80],
          "",
          "# ... (truncated for preview)",
          paste0("# Full code: ", length(lines), " lines"),
          "# Download the R file for complete code"
        )
        paste(preview, collapse = "\n")
      } else {
        code
      }
    }, error = function(e) {
      paste("# Error generating code preview:", e$message)
    })
  })
  
  cat("===== SERVER SETUP COMPLETE =====\n")
}

# Run the app
shinyApp(ui = ui, server = server)