# ==========================================================================
# ui/40_settings.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 351-460  (fPCA / time-warped PCA settings)
# ==========================================================================
ui_tab_settings <- tabItem(
        tabName = "settings",
        fluidRow(
          box(
            title = "PCA Type",
            collapsible = TRUE, collapsed = FALSE,
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
              collapsible = TRUE, collapsed = FALSE,
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
                  # AUDIT (P0.8): allow_dilation, dilation_range and symmetric_warp were read
                  # from the UI and passed as arguments that no function body ever
                  # referenced. Disabled rather than deleted so the intent is on
                  # record; a dilation term needs an estimator, not a checkbox.
                  checkboxInput("allow_dilation", "Allow slight scaling (dilation) - not implemented", FALSE),
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
                  checkboxInput("symmetric_warp", "Force symmetric warping - not implemented", FALSE)
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
            collapsible = TRUE, collapsed = FALSE,
            status = "success",
            solidHeader = TRUE,
            width = 12,
            numericInput("n_components", "Number of components to extract:",
                         value = 5, min = 1, max = 10),
            helpText(HTML("Was three. Functional PCA components are nested \u2014
                           extracting five leaves PC1\u2013PC3 numerically
                           identical and simply makes PC4 and PC5 available, so
                           the display controls downstream are not silently
                           capped at three.")),
            actionButton("run_analysis", "Run Analysis", class = "btn-success")
          )
        )
      )
