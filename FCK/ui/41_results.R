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
            collapsible = TRUE, collapsed = FALSE,
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
              collapsible = TRUE, collapsed = FALSE,
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
              # AUDIT (P8.5): these headings were left over from the version
              # whose statistics they described. The server was corrected at
              # P5.2/P5.3 -- there is no AIC/BIC any more, and the dispersion
              # numbers are explicitly NOT a variance decomposition -- so the UI
              # was announcing "EFDA methodology", "Variance Decomposition" and
              # "Model Selection Criteria" directly above output that says none
              # of those things are on offer. A label that contradicts its own
              # panel is worse than no label.
              h4("Registration Diagnostics"),
              p(HTML(paste(
                "Descriptive measures of how much registration transformed the",
                "curves, and how between-curve dispersion changed. <b>These are",
                "not likelihood-based model-selection criteria</b> and there is",
                "no AIC or BIC here: no probability model for the observed",
                "curves under a candidate registration exists in this module.",
                "Choose a registration method from what you know about the data."))),
              wellPanel(
                style = "background-color: #f8f9fa;",
                fluidRow(
                  column(6,
                    h5("Transformation summary (averaged over subjects)"),
                    verbatimTextOutput("warping_fit_summary")
                  ),
                  column(6,
                    h5("Pre/post registration dispersion"),
                    verbatimTextOutput("warping_variance_decomposition")
                  )
                ),
                hr(),
                fluidRow(
                  column(12,
                    h5("Registration summary"),
                    verbatimTextOutput("warping_model_criteria")
                  )
                )
              ),
              hr(),
              h4("Per-Subject Transformation Diagnostics"),
              p(HTML(paste(
                "How far registration moved each subject's curve. <b>R&sup2;, RMSE",
                "and MAE here compare each curve with its own registered version",
                "&mdash; they measure the MAGNITUDE OF THE TRANSFORMATION, not",
                "how well the sample was aligned.</b> Doing nothing gives a",
                "perfect R&sup2; and zero RMSE, which is not a good registration;",
                "read the dispersion panel above for whether the curves actually",
                "came into alignment."))),
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
            collapsible = TRUE, collapsed = FALSE,
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
            fluidRow(
              column(6,
                     sliderInput("effect_size", "Effect Size Multiplier:",
                                 min = 0.5, max = 3, value = 1, step = 0.1)),
              column(6,
                     sliderInput("effect_n_comp", "Components to show:",
                                 min = 1, max = 10, value = 3, step = 1))
              # The ceiling is retuned to nharm each time the PCA runs; see
              # the observer in server/40_fpca.R.
            ),
            helpText(HTML("Was fixed at three. Beyond the fifth component hue
                           alone stops separating them reliably \u2014 solid/dash is
                           already spoken for by the sign of the deviation \u2014 so
                           read a long list with the legend, or show fewer.")),
            plotlyOutput("scores_plot", height = "400px"),
            hr(),
            DTOutput("scores_table")
          )
        ),

        # ====================================================================
        # NEW (2026-09-03): group comparisons on the component scores.
        #
        # The tab produced an n x k score matrix and stopped. This asks whether
        # the groups differ on each component, with both multiplicity families
        # corrected -- pairwise within a component, and omnibus across
        # components. Server: server/42_fpca_anova.R; arithmetic:
        # server/09_helpers_pcanova.R.
        # ====================================================================
        fluidRow(
          box(
            title = "Group Comparison on Component Scores (ANOVA per PC)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(4,
                     uiOutput("pca_anova_controls")
              ),
              column(8,
                     plotlyOutput("pca_anova_plot", height = "420px"),
                     br(),
                     downloadButton("download_pca_anova",
                                    "Download the full table (CSV)", class = "btn-sm")
              )
            ),
            hr(),
            verbatimTextOutput("pca_anova_results")
          )
        )
      )
