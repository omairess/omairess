# ==========================================================================
# ui/70_fosr.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 275-343  (Function-on-Scalar regression)
# ==========================================================================
ui_tab_fosr <- tabItem(
        tabName = "fosr",
        fluidRow(
          box(
            title = "FoSR Settings (Response: Curve)", status = "primary", solidHeader = TRUE, width = 4,
            collapsible = TRUE, collapsed = FALSE,
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
            collapsible = TRUE, collapsed = FALSE,
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
      )
