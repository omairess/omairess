# ==========================================================================
# ui/72_harmonic.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 424-586  (Harmonic (cosinor) regression)
# ==========================================================================
ui_tab_harmonic <- tabItem(
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
      )
