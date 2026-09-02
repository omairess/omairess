# ==========================================================================
# ui/51_posthoc.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 697-789  (fANOVA post-hoc tests)
# ==========================================================================
ui_tab_posthoc <- tabItem(
        tabName = "pairwise",
        fluidRow(
          box(
            title = "Pairwise Comparison Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.fanova_completed",
              # MERGED APP: choose what these tests compare (default: whatever
              # the omnibus fANOVA used). See FCK/server/52_posthoc_source.R.
              uiOutput("pairwise_source_ui"),
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
      )
