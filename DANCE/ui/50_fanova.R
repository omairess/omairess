# ==========================================================================
# ui/50_fanova.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 579-694  (Functional ANOVA)
# ==========================================================================
ui_tab_fanova <- tabItem(
        tabName = "fanova",
        fluidRow(
          box(
            title = "Functional ANOVA Settings",
            collapsible = TRUE, collapsed = FALSE,
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
              helpText(HTML("For within-subjects designs, specify subject ID and the
                             repeated measures factor (e.g. visit, condition, time).<br>
                             <b>Method:</b> pointwise repeated-measures ANOVA at each
                             time point, with p-values from within-subject permutation
                             of the condition labels. One relabelling is drawn per
                             subject per replicate and applied across the whole
                             trajectory, so the permuted curves stay temporally
                             coherent. This does not use the <code>rmfanova</code>
                             package.")),
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
            conditionalPanel(
              condition = "input.fanova_design == 'within'",
              checkboxInput("rm_global_test",
                            "Also run the global test (rmfanova package)", FALSE),
              helpText(HTML("The pointwise test above says <i>where</i> the conditions
                             differ; it computes no overall p-value. This adds one, via
                             <code>rmfanova</code> (Kurylo &amp; Smaga 2023).<br>
                             It needs a <b>complete balanced design</b> &mdash; every
                             subject in every condition exactly once &mdash; and reports
                             how many subjects it had to drop.<br>
                             Only three of its fifteen outputs are shown: in a simulation
                             of 400 nulls, two were anti-conservative (17.5% and 14.2%
                             rejection at a nominal 5%) and two had no power at all.
                             The three shown were closest to nominal with full power."))
            ),
            numericInput("n_permutations", "Number of permutations for testing:",
                         value = 5000, min = 200, max = 20000),
            helpText(HTML("Was 200, which resolves a p-value only to 0.005 and
                           leaves nothing for an FDR correction across time
                           points to work with. 5,000 is the minimum worth
                           reporting; raise it for a final analysis.")),
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
            collapsible = TRUE, collapsed = FALSE,
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
            collapsible = TRUE, collapsed = FALSE,
            status = "info",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(4,
                checkboxInput("fanova_show_sd_bands", "Show \u00b11 SD bands", value = TRUE),
                helpText(HTML("Bootstrap intervals shown on the post-hoc curves are
                               <b>pointwise</b> 95% percentile intervals from 2,000
                               replicates. They are not simultaneous functional bands:
                               a region where the interval excludes zero is not a
                               family-wise-controlled claim about the whole curve.
                               Use the FDR-adjusted pointwise p-values for that."))
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
            collapsible = TRUE, collapsed = FALSE,
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
            collapsible = TRUE, collapsed = FALSE,
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("fanova_effect_size_plot", height = "400px"),
            hr(),
            verbatimTextOutput("fanova_effect_summary")
          )
        )
      )
