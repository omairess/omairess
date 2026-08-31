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
      status = "warning",
      solidHeader = TRUE,
      width = 6,
      radioButtons("smooth_method", "Smoothing Method:",
                   choices = list("Raw data (no smoothing)" = "none",
                                  "Automatic smoothing (REML)" = "auto",
                                  "Manual smoothing" = "manual"),
                   selected = "auto"),

      # From CIRCAREG: periodic data gets a Fourier basis instead of B-splines.
      checkboxInput("is_cyclic", "Is data cyclic? (e.g. 24h clock)", FALSE),
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
          helpText("Click to automatically set smoothing factor based on REML/CV analysis")
        )
      ),
      conditionalPanel(
        condition = "input.smooth_method == 'auto'",
        numericInput("n_basis", "Number of B-spline basis functions:",
                     value = 20, min = 4, max = 100),
        helpText("Lambda = 0 uses REML optimization automatically")
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
      status = "info",
      solidHeader = TRUE,
      width = 6,
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
  )
)
