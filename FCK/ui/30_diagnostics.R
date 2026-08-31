# ==========================================================================
# ui/30_diagnostics.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 232-348  (Smoothing Diagnostics tab)
# ==========================================================================
ui_tab_smooth_diag <- tabItem(
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
      )
