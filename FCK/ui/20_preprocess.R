# ==============================================================================
# ui/20_preprocess.R — SHARED preprocessing / smoothing (tab "preprocess")
#
# Merged by hand from the two source apps' smoothing tabs:
#   WaPaa1_3.R lines 145-229  (separate n_basis for auto vs manual, the
#                              GCV-vs-n-basis sweep, the interactive curve
#                              plot with click-to-select, fit statistics)
#   CIRCAREG.R lines 111-154  (the cyclic/Fourier basis option and the compact
#                              fit-metrics panel)
#
# Both apps ran the SAME smoothing algorithm here (WaPaa's code says so in its
# own comments: "Applying smoothing using CIRCAREG method").  The union of the
# controls is offered once, and its output — values$smooth_data and
# values$fd_obj — is what every analysis tab in the app reads.
# ==============================================================================

ui_tab_preprocess <- tabItem(
  tabName = "preprocess",
  fluidRow(
    box(
      title = "Smoothing Options",
      collapsible = TRUE, collapsed = FALSE,
      status = "warning",
      solidHeader = TRUE,
      width = 6,
      radioButtons("smooth_method", "Smoothing Method:",
                   choices = list("Raw data (no smoothing)" = "none",
                                  "Automatic smoothing (GCV)" = "auto",
                                  "Manual smoothing" = "manual"),
                   selected = "auto"),

      # From CIRCAREG: periodic data gets a Fourier basis instead of B-splines.
      checkboxInput("is_cyclic", "Is data cyclic? (e.g. 24h clock)", FALSE),

      # Both source apps smoothed against the column INDEX (1, 2, 3, ...),
      # which is only right when the columns are evenly spaced in real time.
      # Off by default so existing results are reproduced exactly.
      checkboxInput("use_real_time",
                    "Space time points by their real clock times", FALSE),
      conditionalPanel(
        condition = "input.use_real_time == true",
        helpText(HTML(
          "Uses the hours parsed from the column names (e.g. <code>Base9h</code>,
           <code>KSS_9u_dag1</code>, <code>20h15</code>, <code>08:00</code>),
           unwrapped across midnight, as the smoothing argument instead of
           1, 2, 3, &hellip;<br>
           <b>Tick this if your measurements are unevenly spaced</b> — hourly by
           day and 2-hourly at night, say — otherwise a 2-hour gap is smoothed as
           though it were as short as a 1-hour one.<br>
           The roughness penalty is then per <i>hour</i> rather than per column,
           so the smoothing factor will need re-tuning; the Data Import tab
           reports whether clock times were parsed at all."))
      ),
      conditionalPanel(
        condition = "input.is_cyclic == true",
        helpText(HTML(
          "A Fourier basis is used instead of B-splines, with <b>min(n_time, 13)</b>
           basis functions — the rule CIRCAREG used. The number of B-spline bases
           below is ignored while this is ticked."))
      ),

      checkboxInput("constrain_bounds", "Constrain smoothed values to specific range", FALSE),
      conditionalPanel(
        condition = "input.constrain_bounds == true",
        numericInput("min_bound", "Minimum value:", value = 0, step = 1),
        numericInput("max_bound", "Maximum value:", value = 100, step = 1),
        helpText("Smoothed values will be clamped to [min, max] range.")
      ),
      conditionalPanel(
        condition = "input.smooth_method == 'manual'",
        sliderInput("smooth_factor", "Smoothing Factor:",
                    min = 0.1, max = 10, value = 1, step = 0.1),
        numericInput("n_basis_manual", "Number of B-spline basis functions:",
                     value = 20, min = 4, max = 100),
        helpText(HTML("<b>Smoothing Factor explained:</b><br>
                       Lambda = 10^(-Smoothing Factor)<br>
                       - Lower factor (0.1-2) → Higher lambda → Smoother curves<br>
                       - Higher factor (6-10) → Lower lambda → More flexible curves<br>
                       <br><b>Number of B-splines:</b><br>
                       - More bases → More potential flexibility<br>
                       - Typical: 10-20 for most datasets<br>
                       - Max recommended: n_time_points - 2<br>
                       <br><b>💡 Tip:</b> Use 'Smoothing Diagnostics' tab to find optimal values!")),
        conditionalPanel(
          condition = "output.diagnostics_available",
          hr(),
          actionButton("use_diagnostic_lambda", "📊 Use Diagnostic Results",
                       class = "btn-info btn-sm"),
          helpText("Sets the manual smoothing factor from the diagnostics tab. Note that tab fits an mgcv model, whose smoothing parameter is not on the same scale as fda's lambda; treat the transferred value as a starting point, not a selection.")
        )
      ),
      conditionalPanel(
        condition = "input.smooth_method == 'auto'",
        numericInput("n_basis", "Number of B-spline basis functions:",
                     value = 20, min = 4, max = 100),
        helpText(HTML("Lambda is chosen by minimising the mean GCV score across
                       subjects, using the same smoother that then does the
                       smoothing. This control previously said REML and set
                       lambda to zero, which applied no penalty at all &mdash; on a
                       basis the size of your time grid that is exact
                       interpolation, not a smooth."))
      ),
      conditionalPanel(
        condition = "input.smooth_method != 'none' && input.is_cyclic == false",
        hr(),
        actionButton("run_nbasis_analysis", "Analyse Optimal n-Basis",
                     class = "btn-sm btn-info", icon = icon("chart-line")),
        helpText("Tests a range of n-basis values and plots the mean GCV score. Lower GCV = better fit–complexity trade-off."),
        plotlyOutput("nbasis_gcv_plot", height = "280px")
      ),
      hr(),
      h5("Missing values"),
      helpText(HTML(
        "Smoothing fills every gap: a subject measured at 8 of 20 times comes out
         with 20 values. Filling <i>between</i> two of that subject's measurements
         is what smoothing is for. Filling <i>beyond</i> their first or last one is
         not — a spline carried past its data follows its end polynomial, so a
         short recording in a long protocol can come back with arbitrary values.")),
      checkboxInput("allow_extrapolation",
                    "Let the fitted curve run beyond each subject's observed range",
                    FALSE),
      conditionalPanel(
        condition = "input.allow_extrapolation == false",
        helpText("Values past a subject's first/last measurement are held flat at",
                 "that measurement. They are still invented — the map below marks them.")
      ),
      conditionalPanel(
        condition = "input.allow_extrapolation == true",
        helpText(HTML("<b>The curve is extrapolated.</b> Check the map below for how
                       far past the data that reaches before using these curves."))
      ),
      hr(),
      actionButton("apply_smooth", "Apply Smoothing", class = "btn-warning"),
      helpText("This is the only smoothing step in the app: every analysis tab",
               "reads the curves it produces."),
      hr(),
      # CIRCAREG's compact metrics panel...
      uiOutput("smooth_fit_display"),
      # ...and WaPaa's full statistics printout.
      h4("Smoothing Fit Statistics:"),
      verbatimTextOutput("smoothing_fit_summary")
    ),

    box(
      title = "Data Visualization",
      collapsible = TRUE, collapsed = FALSE,
      status = "info",
      solidHeader = TRUE,
      width = 6,
      selectInput("data_plot_n", "Curves to draw:",
                  choices = c("First 10" = 10, "First 25" = 25, "First 50" = 50,
                              "First 100" = 100, "First 250" = 250,
                              "All" = 0),
                  selected = 50, width = "200px"),
      helpText("Drawing every curve of a large set is slow and turns the plot",
               "into a solid block; the first 50 is usually enough to judge the",
               "smoothing. The title says how many of how many are shown."),
      selectInput("tick_freq_preprocess", "X-axis tick frequency:",
                  choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                 "Every hour" = 1, "Every 2 hours" = 2,
                                 "Every 4 hours" = 4, "All" = 0),
                  selected = 1, width = "200px"),
      plotlyOutput("data_plot", height = "400px"),
      hr(),
      wellPanel(
        style = "background-color: #f0f8ff; padding: 10px;",
        h5("Selected Curve Info"),
        p("Click on a curve in the plot to select it.", style = "font-size: 11px; color: #666;"),
        verbatimTextOutput("selected_curve_info"),
        actionButton("clear_curve_selection", "Clear Selection", class = "btn-sm btn-default")
      )
    )
  ),

  fluidRow(
    box(
      title = "Missing data & filled points",
      status = "danger", solidHeader = TRUE, width = 12, collapsible = TRUE,
      helpText("Which values in the smoothed curves are measurements, and which",
               "the smoother supplied. Apply smoothing first."),
      verbatimTextOutput("missing_headline"),
      uiOutput("missing_legend"),
      hr(),
      fluidRow(
        column(8,
          selectInput("missing_sort", "Order rows by:",
                      choices = c("As in the file" = "file",
                                  "Fewest measured points first" = "observed",
                                  "Most filled-in points first" = "filled"),
                      selected = "file", width = "260px"),
          plotlyOutput("missing_map", height = "420px")
        ),
        column(4,
          uiOutput("missing_subject_ui"),
          plotlyOutput("missing_curve", height = "360px")
        )
      ),
      hr(),
      DTOutput("missing_table"),
      br(),
      downloadButton("export_fill_status_csv",
                     "Download the observed/interpolated/extrapolated map (CSV)",
                     class = "btn-primary")
    )
  )
)
