# ==========================================================================
# ui/72_harmonic.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   CIRCAREG.R lines 424-586  (Harmonic (cosinor) regression)
#
# CHANGELOG - 2026-09-03 cosinor audit (this file is NO LONGER verbatim)
# Adds the DV name/units/bounds inputs (1.7), the raw-vs-smoothed data source
# (2.1), the time-origin selector (1.4.3 / 2.2), the nested-model comparison and
# fixed-tau controls (2.4 / 2.2), and renames the MESOR bound inputs to what
# they actually bound: the intercept (1.4).
# ==========================================================================
ui_tab_harmonic <- tabItem(
        tabName = "harmonic",
        fluidRow(
          box(
            title = "Harmonic Regression Settings", status = "success", solidHeader = TRUE, width = 4,
            collapsible = TRUE, collapsed = FALSE,
            h4("Data Selection"),
            uiOutput("harmonic_var_select_ui"),

            # ================================================================
            # AUDIT 1.7 / 2.1 / 1.4.3: what is being modelled, and on what
            # ================================================================
            hr(),
            h4("The dependent variable"),
            helpText(HTML("The report never named the DV, so nobody reading it could
                           check whether the fitted curve was even admissible.")),
            textInput("harmonic_dv_name", "Name:", value = "", placeholder = "e.g. KSS sleepiness"),
            textInput("harmonic_dv_units", "Units:", value = "", placeholder = "e.g. KSS points"),
            fluidRow(
              column(6, numericInput("harmonic_dv_min", "Lowest possible value:", value = NA, step = 1)),
              column(6, numericInput("harmonic_dv_max", "Highest possible value:", value = NA, step = 1))
            ),
            helpText(HTML("Fitted values outside this range are flagged. With
                           <i>M + A\u2081 + A\u2082</i> the trough of a reported pooled fit can
                           sit below zero before the trend is added, which is
                           structurally impossible for a non-negative scale.")),

            hr(),
            h4("Data source"),
            radioButtons("harmonic_data_source", NULL,
                         choices = c("Raw (cosinor handles gaps natively)" = "raw",
                                     "Smoothed (FDA-interpolated)" = "smoothed"),
                         # P6.7: raw is the default. Cosinor is a regression on
                         # the observations and handles gaps natively; smoothing
                         # first buys nothing here and costs the inflation the
                         # note below describes. Defaulting to the option that
                         # inflates R2 and anticonservatively biases the
                         # zero-amplitude test was the wrong way round.
                         selected = "raw"),
            helpText(HTML("Fitting on smoothed data removes independent noise and
                           induces residual autocorrelation: R\u00b2 is inflated, LOOCV is
                           optimistic because a held-out point is partly rebuilt from
                           its neighbours, and the zero-amplitude F test is
                           anticonservative. <b>Run both</b> and compare \u2014 the gap is
                           the inflation.")),

            hr(),
            h4("Model Specification"),
            radioButtons("harmonic_time_origin", "Time origin (t = 0 at):",
                         choices = c("First observation" = "first_observation",
                                     "Midnight" = "midnight"),
                         # P6.7: the first observation is the default. With a
                         # saturating trend, anchoring S(t) at the first
                         # observation while the harmonics sit at midnight gives
                         # an intercept that is the value at neither origin;
                         # one origin makes it interpretable and conditions the
                         # A_sat/tau pair better. Midnight remains available for
                         # continuity with earlier runs.
                         selected = "first_observation"),
            helpText(HTML("With a saturating trend the fitter anchors <i>S(t)</i> at the
                           first observation while the harmonics stay anchored at
                           midnight \u2014 two origins, one constant, so the intercept is
                           the value at neither. Re-anchoring both makes the intercept
                           interpretable and improves the conditioning of the
                           <i>A_sat</i>/<i>\u03c4</i> pair.")),
            numericInput("harmonic_period", "Fundamental Period (τ):", value = 24, min = 1, max = 168, step = 1),
            helpText("Period in the same units as your time variable (e.g., 24 for circadian)."),
            sliderInput("n_harmonics", "Number of Harmonics:", min = 1, max = 8, value = 1, step = 1),
            helpText("1 = fundamental (τ), 2 = adds τ/2, 3 = adds τ/3, etc."),
            uiOutput("harmonic_warning_ui"),
            selectInput("harmonic_trend_type", "Homeostatic Trend Model:",
                        choices = c("None (circadian only)" = "none",
                                    "Linear" = "linear",
                                    "Logarithmic" = "log",
                                    "Saturating exponential" = "exp_sat"),
                        selected = "none"),
            conditionalPanel(
              condition = "input.harmonic_trend_type != 'none'",
              helpText(HTML("<b>Two-Process Model:</b> Separates homeostatic sleep pressure (trend) from circadian modulation (sinusoidal)."))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'linear'",
              helpText(HTML("Y(t) = M + β·t + harmonics<br><small>Constant rate of change</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'log'",
              helpText(HTML("Y(t) = M + β·log(t+1) + harmonics<br><small>Rapid initial change, slowing over time</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'exp_sat'",
              helpText(HTML("Y(t) = M + A·(1-e<sup>-t/τ</sup>) + harmonics<br>
                            <small>Asymptotic saturation (classic Process S)</small><br>
                            <small><em>⚠️ Uses nonlinear least squares. If fit is poor, try Log trend instead.</em></small>"))
            ),
            hr(),
            h4("Group Analysis"),
            uiOutput("harmonic_group_var_ui"),
            helpText("Optional: Select a group variable to compare rhythms between groups."),
            hr(),
            h4("Which fits enter the summaries"),
            checkboxInput("harmonic_include_boundary",
                          "Include fits that hit a parameter bound", TRUE),
            helpText(HTML("A fit pinned to a constraint converged to the <b>edge of
                           the feasible region</b>, not to an interior optimum: the
                           value is where the optimiser was stopped and its standard
                           error is meaningless there. Including them keeps the whole
                           sample and lets the bound show through the mean; excluding
                           them gives a cleaner summary of a smaller, possibly biased
                           set. Either way the report tables <b>which</b> bound each
                           fit hit and <b>which fits hit more than one</b> \u2014 two
                           pinned parameters usually means a ridge.<br><br>
                           Non-converged fits are always excluded: there is no
                           solution to average.")),

            hr(),
            h4("Diagnostics"),
            checkboxInput("harmonic_model_selection",
                          "Compare nested models (\u0394AICc table)", FALSE),
            helpText(HTML("Fits trend \u2208 {none, linear, saturating} \u00d7 harmonics \u2208 {1,2,3}
                           and reports \u0394AICc with Akaike weights. Absolute AIC/AICc/BIC
                           with no competing model are constant offsets of one another
                           and carry no information. <b>Slow:</b> 9 fits per subject.")),
            numericInput("harmonic_tau_fixed",
                         "\u03c4 held at (h), for the free-vs-fixed \u0394AIC:", value = 18,
                         min = 1, max = 72, step = 0.5),
            helpText(HTML("Daan, Beersma & Borb\u00e9ly (1984) give \u03c4_rise \u2248 18 h under
                           extended wakefulness. If free \u03c4 does not beat fixed \u03c4 on AIC,
                           \u03c4 is not identified by these data and a between-group
                           comparison of it is a comparison of where the optimiser
                           stopped on a ridge.")),

            hr(),
            h4("Bootstrap Options"),
            checkboxInput("harmonic_bootstrap", "Compute Bootstrap CIs", FALSE),
            conditionalPanel(
              condition = "input.harmonic_bootstrap == true",
              numericInput("harmonic_n_boot", "Bootstrap Iterations (B):", value = 500, min = 100, max = 2000)
            ),
            hr(),
            h4("Parameter Constraints"),
            checkboxInput("harmonic_use_bounds", "Enable Parameter Bounding", FALSE),
            conditionalPanel(
              condition = "input.harmonic_use_bounds == true",
              helpText("Constrain parameters to stay within plausible ranges based on your data scale."),
              uiOutput("harmonic_bounds_hints"),
              hr(),
              h5("Common Parameters"),
              numericInput("harmonic_mesor_min", "Intercept (\u03b2\u2080) Min:", value = NA, step = 0.1),
              numericInput("harmonic_mesor_max", "Intercept (\u03b2\u2080) Max:", value = NA, step = 0.1),
              numericInput("harmonic_amplitude_min", "Amplitude Min:", value = 0, step = 0.1),
              numericInput("harmonic_amplitude_max", "Amplitude Max:", value = NA, step = 0.1),
              conditionalPanel(
                condition = "input.harmonic_trend_type == 'exp_sat'",
                hr(),
                h5("Exponential Saturation Parameters"),
                numericInput("harmonic_A_sat_min", "A_sat (Asymptote) Min:", value = NA, step = 0.1),
                numericInput("harmonic_A_sat_max", "A_sat (Asymptote) Max:", value = NA, step = 0.1),
                numericInput("harmonic_tau_min", "τ (Time Constant) Min:", value = 0.5, step = 0.1),
                numericInput("harmonic_tau_max", "τ (Time Constant) Max:", value = NA, step = 1),
                helpText(HTML("<small>τ controls how fast saturation is reached. Smaller τ = faster saturation.</small>"))
              ),
              helpText(HTML("<small>Leave blank (NA) for no bound. Defaults based on data range shown in hints above.</small>"))
            ),
            hr(),
            actionButton("run_harmonic", "Run Harmonic Regression", class = "btn-success", icon = icon("play"))
          ),
          # The summary grew a great deal during the audit -- fit outcomes, bound
          # tables, commonality, conditioning, per-group detail. All of it earns
          # its place when you are checking a model, and none of it does when you
          # are not, so the box collapses and its text scrolls inside a fixed
          # height instead of pushing the rest of the tab off the screen.
          box(
            title = "Model Summary", status = "info", solidHeader = TRUE, width = 8,
            collapsible = TRUE, collapsed = FALSE,
            fluidRow(
              column(7,
                     sliderInput("harmonic_summary_height", "Panel height (px):",
                                 min = 200, max = 2000, value = 600, step = 100,
                                 width = "100%")),
              column(5,
                     br(),
                     downloadButton("download_harmonic_summary",
                                    "Save the summary as text", class = "btn-sm"))
            ),
            helpText(HTML("<small>Collapse the box with the <b>\u2212</b> in its header.
                           The summary is long by design \u2014 it is the audit trail \u2014
                           so it scrolls here rather than pushing the tab down.</small>")),
            uiOutput("harmonic_summary_box")
          )
        ),
        fluidRow(
          tabBox(
            title = "Harmonic Regression Results", id = "harmonic_tabs", width = 12,
            tabPanel("1. Fitted Curves", icon = icon("chart-line"),
                     fluidRow(
                       column(8, plotlyOutput("harmonic_fit_plot", height = "500px")),
                       column(4, 
                              h4("Display Options"),
                              uiOutput("harmonic_subject_selector"),
                              checkboxInput("harmonic_show_ci", "Show Confidence Bands", TRUE),
                              # P6.7: off by default. On a few hundred subjects
                              # the points bury the fitted curve they are there
                              # to support.
                              checkboxInput("harmonic_show_data", "Show Raw Data Points", FALSE),
                              checkboxInput("harmonic_show_components", "Show Harmonic Components", FALSE)
                       )
                     )
            ),
            tabPanel("2. Polar Plot (Acrophase)", icon = icon("compass"),
                     fluidRow(
                       # Same height as the Polar Density tab next door, so the
                       # two dials are directly comparable when you flip between
                       # them rather than one appearing to have a tighter spread
                       # purely because it is drawn smaller.
                       column(8, plotlyOutput("harmonic_polar_plot", height = "540px")),
                       column(4,
                              h4("Polar Plot Settings"),
                              uiOutput("harmonic_selector_polar"),
                              helpText("Acrophase displayed in polar coordinates. Radius = Amplitude, Angle = Acrophase."),
                              checkboxInput("polar_show_mean", "Show Population Mean Vector", TRUE),
                              checkboxInput("polar_show_ellipse", "Show Confidence Ellipse", TRUE)
                       )
                     )
            ),
            # MERGED APP: the same circle as a DENSITY, oriented as a clock
            # face -- noon at the top, midnight at the bottom, so the upper
            # half is daytime and the lower half night.  See
            # FCK/server/07_helpers_circular.R and 74_polar_density.R.
            tabPanel("2b. Polar Density", icon = icon("circle-notch"),
                     fluidRow(
                       column(8,
                              plotlyOutput("harmonic_density_plot", height = "540px"),
                              verbatimTextOutput("harmonic_density_note")),
                       column(4,
                              h4("Density Settings"),
                              helpText(HTML(
                                "A filled shape around a clock face. The radius
                                 carries either a <b>von Mises kernel density of
                                 the acrophases</b> -- the circular analogue of a
                                 KDE, so mass near midnight wraps instead of
                                 splitting between the ends of a histogram -- or
                                 the <b>signal itself averaged over the
                                 clock</b>.<br><br>
                                 <b>Noon is at the top, midnight at the bottom:</b>
                                 the upper half of the circle is daytime
                                 (06:00-18:00), the lower half night. Hours run
                                 clockwise.")),
                              uiOutput("density_controls_ui")
                       )
                     )
            ),
            tabPanel("3. Parameter Distribution", icon = icon("chart-bar"),
                     fluidRow(
                       column(12, 
                              uiOutput("harmonic_selector_dist"),
                              hr()
                       )
                     ),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_amplitude_hist", height = "400px")),
                       column(6, plotlyOutput("harmonic_acrophase_hist", height = "400px"))
                     ),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_mesor_plot", height = "300px")),
                       column(6, plotlyOutput("harmonic_trend_hist", height = "300px"))
                     )
            ),
            tabPanel("4. Individual Results", icon = icon("table"),
                     fluidRow(
                       column(12,
                              h4("Individual Subject Parameters"),
                              helpText("R-squared and other parameters for each subject. Use the export button to download."),
                              DT::dataTableOutput("harmonic_individual_table"),
                              br(),
                              downloadButton("export_harmonic_individual", "Export Individual Parameters (CSV)", class = "btn-success")
                       )
                     )
            ),
            tabPanel("5. Residual Diagnostics", icon = icon("microscope"),
                     fluidRow(
                       column(6, plotlyOutput("harmonic_residual_plot", height = "400px")),
                       column(6, plotlyOutput("harmonic_qq_plot", height = "400px"))
                     ),
                     fluidRow(
                       column(12, verbatimTextOutput("harmonic_gof_stats"))
                     )
            ),
            tabPanel("6. Group Comparison", icon = icon("users"),
                     fluidRow(
                       column(12, 
                              uiOutput("harmonic_selector_group"),
                              hr()
                       )
                     ),
                     fluidRow(
                       column(12, plotlyOutput("harmonic_group_comparison_plot", height = "500px"))
                     ),
                     hr(),
                     verbatimTextOutput("harmonic_group_test_results")
            )
          )
        )
      )
