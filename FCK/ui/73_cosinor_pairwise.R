# ==========================================================================
# ui/73_cosinor_pairwise.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 589-654  (cosinor pairwise group tests; ids prefixed hp_)
# ==========================================================================
ui_tab_cosinor_pairwise <- tabItem(
        tabName = "harm_pairwise",
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
            selectInput("hp_param", "Parameter to Compare:",
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
            selectInput("hp_correction", "Multiple Comparison Correction:",
                       choices = c("None" = "none",
                                   "Bonferroni" = "bonferroni",
                                   "Holm" = "holm",
                                   "Hochberg" = "hochberg",
                                   "Hommel" = "hommel",
                                   "Benjamini-Hochberg (FDR)" = "BH",
                                   "Benjamini-Yekutieli (FDR)" = "BY"),
                       selected = "holm"),
            hr(),
            checkboxInput("hp_show_effect_size", "Show Effect Sizes (Cohen's d)", value = TRUE),
            checkboxInput("hp_show_ci", "Show 95% Confidence Intervals", value = TRUE),
            hr(),
            actionButton("hp_run", "Run Pairwise Comparisons", class = "btn-warning", icon = icon("play"))
          ),
          box(
            title = "Pairwise Test Results", status = "success", solidHeader = TRUE, width = 8,
            verbatimTextOutput("hp_results")
          )
        ),
        fluidRow(
          box(
            title = "Visualization", status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("hp_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Comparison Matrix", status = "info", solidHeader = TRUE, width = 6,
            uiOutput("hp_matrix_help"),
            tableOutput("hp_matrix")
          ),
          box(
            title = "Download", status = "warning", solidHeader = TRUE, width = 6,
            downloadButton("hp_export_results", "Download Pairwise Results (CSV)", class = "btn-primary"),
            br(), br(),
            downloadButton("hp_export_plot", "Download Plot (PNG)", class = "btn-success")
          )
        )
      )
