# ==========================================================================
# ui/71_sofr.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 346-421  (Scalar-on-Function regression)
# ==========================================================================
ui_tab_sofr <- tabItem(
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
      )
