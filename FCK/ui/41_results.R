# ==========================================================================
# ui/41_results.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 463-576  (Functional PCA results)
# ==========================================================================
ui_tab_results <- tabItem(
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
      )
