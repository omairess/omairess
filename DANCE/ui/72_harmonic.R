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
            helpText(HTML("Fitted values outside this range are flagged: a curve whose
                           trough dips below zero, for instance, is impossible on a
                           non-negative scale.")),

            hr(),
            # ================================================================
            # THE APPROACH. One model over every observation, or one fit per
            # participant summarised afterwards. The same questions are asked
            # either way and every tab answers them from whichever was chosen,
            # so nothing on screen mixes the two estimators without saying so.
            # ================================================================
            h4("Approach"),
            radioButtons("harmonic_approach", NULL,
                         choices = c("Mixed-effects cosinor (one model, all participants)" = "mixed",
                                     "Two-stage cosinor (one fit per participant, then compared)" = "two_stage"),
                         selected = "mixed"),
            helpText(HTML("<b>Mixed-effects</b>: participants share a group-level rhythm and
                           each carries their own deviation from it, estimated jointly;
                           group comparisons come from the model's own contrasts and respect
                           a repeated-measures design. Individual results are the model's
                           <i>shrunken</i> predictions for each participant.<br>
                           <b>Two-stage</b>: an independent cosinor per participant, then the
                           point estimates are compared. Simple and classical; it treats
                           every estimate as exact and cannot express a within-participant
                           factor.")),

            hr(),
            h4("Model Specification"),
            # The time origin is not a choice any more: t = 0 is the first
            # observation, always. Both the trend and the harmonics are anchored
            # there, so the intercept is the model value at the start of the
            # recording; axes and every acrophase still read in clock time.
            # The period is P, not tau: tau is the saturating trend's time
            # constant, and one symbol for two quantities made "tau held at
            # 18 h" readable as a period.
            numericInput("harmonic_period", "Fundamental period P:", value = 24, min = 1, max = 168, step = 1),
            helpText("Period in the same units as your time variable (e.g., 24 for circadian)."),
            sliderInput("n_harmonics", "Number of Harmonics:", min = 1, max = 8, value = 1, step = 1),
            helpText("1 = fundamental (P), 2 = adds P/2, 3 = adds P/3, etc."),
            uiOutput("harmonic_warning_ui"),
            selectInput("harmonic_trend_type", "Non-periodic Trend Model:",
                        choices = c("None (circadian only)" = "none",
                                    "Linear" = "linear",
                                    "Logarithmic" = "log",
                                    "Saturating exponential" = "exp_sat"),
                        selected = "none"),
            conditionalPanel(
              condition = "input.harmonic_trend_type != 'none'",
              # The control used to be called "Homeostatic Trend Model" and this
              # note called the result a "Two-Process Model". A cosinor plus a
              # trend separates a PERIODIC component from a NON-PERIODIC one.
              # That may be informative about circadian modulation and
              # homeostatic accumulation, but it is not the mechanistic Borbely
              # two-process model, and naming it so asserts a mechanism the fit
              # does not contain.
              helpText(HTML("Separates a <b>periodic</b> component (the harmonics) from a
                             <b>non-periodic</b> change over time (the trend). The trend may
                             be informative about homeostatic processes, but this is a
                             descriptive decomposition, not a mechanistic
                             two-process model."))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'linear'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·t + harmonics<br><small>Constant rate of
                             change: β is in response units per hour</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'log'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·log(t+1) + harmonics<br><small>Rapid initial
                             change, slowing over time</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'exp_sat'",
              helpText(HTML("Y(t) = β<sub>0</sub> + A·(1-e<sup>-t/τ</sup>) + harmonics<br>
                            <small>Asymptotic saturation \u2014 the shape a homeostatic process is
                            often assumed to take, fitted descriptively. A is the rise above
                            β<sub>0</sub> the trend approaches, in response units; τ is its
                            time constant in hours.</small>")),
              # THE TIME CONSTANT: one value every fit holds, or a range every
              # estimate stays inside. The same value and range serve both
              # approaches; what "estimate" means differs and is said below.
              # The mixed model is linear in its coefficients only for a given
              # tau, so its trajectory model shares ONE tau across cells; a tau
              # per cell is a nonlinear model of its own, offered below.
              conditionalPanel(
                condition = "input.harmonic_approach == 'mixed'",
                radioButtons("harmonic_traj_tau_mode", "Time constant \u03c4:",
                             choices = c("Hold \u03c4 at a value" = "fixed",
                                         "Estimate one \u03c4 for all cells, within a range" = "profile"),
                             selected = "fixed")),
              conditionalPanel(
                condition = "input.harmonic_approach == 'two_stage'",
                radioButtons("harmonic_ts_tau_mode", "Time constant \u03c4:",
                             choices = c("Estimate \u03c4 for each participant, within a range" = "estimate",
                                         "Hold \u03c4 at a value for every participant" = "fixed"),
                             selected = "estimate")),
              conditionalPanel(
                condition = paste("(input.harmonic_approach == 'mixed' && input.harmonic_traj_tau_mode == 'fixed') ||",
                                  "(input.harmonic_approach == 'two_stage' && input.harmonic_ts_tau_mode == 'fixed')"),
                numericInput("harmonic_traj_tau", "\u03c4 held at (h):", value = 18,
                             min = 0.5, max = 168, step = 0.5)),
              conditionalPanel(
                condition = paste("(input.harmonic_approach == 'mixed' && (input.harmonic_traj_tau_mode == 'profile' || input.harmonic_traj_sat_cells)) ||",
                                  "(input.harmonic_approach == 'two_stage' && input.harmonic_ts_tau_mode == 'estimate')"),
                fluidRow(
                  column(6, numericInput("harmonic_tau_lo", "Search \u03c4 from (h):", value = 1,
                                         min = 0.1, step = 0.5)),
                  column(6, numericInput("harmonic_tau_hi", "to (h):", value = 72,
                                         min = 0.5, step = 1)))),
              helpText(HTML("<small>Daan, Beersma &amp; Borb\u00e9ly (1984) give
                             τ<sub>rise</sub> \u2248 18 h under extended wakefulness. An
                             estimate is never allowed outside the range; one that ends AT
                             either end is a bound, not an estimate, and is flagged.</small>")),
              conditionalPanel(
                condition = "input.harmonic_approach == 'mixed'",
                helpText(HTML("<small>Every cell shares one τ and the model estimates the rise A
                               per cell. Estimating τ profiles the likelihood at 24 values
                               across the range, which takes roughly 15 to 20 times as long as
                               one fit.</small>")),
                checkboxInput("harmonic_traj_sat_cells",
                              "Also estimate A and \u03c4 in every cell (nonlinear mixed model)",
                              FALSE),
                helpText(HTML("<small>Fits A, τ and the initial slope A/τ separately in each
                               design cell (nlme), τ within the same range, with intervals and
                               tests, and tests whether one τ is enough. Shown on the
                               comparison tab.</small>"))),
              conditionalPanel(
                condition = "input.harmonic_approach == 'two_stage'",
                helpText(HTML("<small><em>Each participant's A and τ are fitted by nonlinear least
                               squares; with τ held, A is a least-squares fit. If fits are poor,
                               try the log trend.</em></small>")))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type != 'none'",
              # The control used to be called "Homeostatic Trend Model" and this
              # note called the result a "Two-Process Model". A cosinor plus a
              # trend separates a PERIODIC component from a NON-PERIODIC one.
              # That may be informative about circadian modulation and
              # homeostatic accumulation, but it is not the mechanistic Borbely
              # two-process model, and naming it so asserts a mechanism the fit
              # does not contain.
              helpText(HTML("Separates a <b>periodic</b> component (the harmonics) from a
                             <b>non-periodic</b> change over time (the trend). The trend may
                             be informative about homeostatic processes, but this is a
                             descriptive decomposition, not a mechanistic
                             two-process model."))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'linear'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·t + harmonics<br><small>Constant rate of
                             change: β is in response units per hour</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'log'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·log(t+1) + harmonics<br><small>Rapid initial
                             change, slowing over time</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type == 'exp_sat'",
              helpText(HTML("Y(t) = β<sub>0</sub> + A·(1-e<sup>-t/τ</sup>) + harmonics<br>
                            <small>Asymptotic saturation \u2014 the shape a homeostatic process is
                            often assumed to take, fitted descriptively. A is the rise above
                            β<sub>0</sub> the trend approaches, in response units; τ is its
                            time constant in hours.</small>")),
              conditionalPanel(
                condition = "input.harmonic_approach == 'two_stage'",
                helpText(HTML("<small><em>Two-stage estimates A and τ for every participant by
                               nonlinear least squares. If fits are poor, try the log
                               trend.</em></small>"))
              ),
              # The mixed model is linear in its coefficients only for a given
              # tau, so the trajectory model shares ONE tau across cells: held at
              # a value, or estimated by profiling the likelihood. Estimating tau
              # per cell is a nonlinear model of its own, offered below it.
              # Visible here because tau IS part of the model -- it used to be
              # read from a hidden two-stage field.
              conditionalPanel(
                condition = "input.harmonic_approach == 'mixed'",
                radioButtons("harmonic_traj_tau_mode", "Time constant \u03c4:",
                             choices = c("Hold \u03c4 at a value" = "fixed",
                                         "Estimate one \u03c4 for all cells (profile likelihood)" = "profile"),
                             selected = "fixed"),
                conditionalPanel(
                  condition = "input.harmonic_traj_tau_mode == 'fixed'",
                  numericInput("harmonic_traj_tau", "\u03c4 held at (h):", value = 18,
                               min = 0.5, max = 168, step = 0.5)),
                helpText(HTML("<small>Every cell shares one τ and the model estimates
                               the rise A per cell. Daan, Beersma &amp; Borb\u00e9ly
                               (1984) give τ<sub>rise</sub> \u2248 18 h under extended
                               wakefulness. Estimating τ refits the model at about 30
                               values of τ, which takes roughly 15 to 20 times as long
                               as one fit.</small>")),
                checkboxInput("harmonic_traj_sat_cells",
                              "Also estimate A and \u03c4 in every cell (nonlinear mixed model)",
                              FALSE),
                helpText(HTML("<small>Fits A, τ and the initial slope A/τ separately
                               in each design cell (nlme), with intervals and tests,
                               and tests whether one τ is enough. Shown on the
                               comparison tab.</small>"))
              )
            ),
            conditionalPanel(
              condition = "input.harmonic_trend_type != 'none'",
              helpText(HTML("<small>β<sub>0</sub> is the <b>intercept</b>: the midline at the
                             first observation, where the trend is zero. The <b>MESOR</b> is the
                             midline β<sub>0</sub> + trend averaged over the recording. Without
                             a trend the two are the same number.</small>"))
            ),
            # ================================================================
            # STUDY DESIGN (P21 phase 4)
            # ================================================================
            # This replaces a single optional "Group Variable". That control
            # could express exactly one design -- one between-participant factor
            # -- and a repeated-measures study had nowhere to say so, which is
            # how the same participant's conditions ended up treated as
            # independent. The design is now named explicitly, because it
            # decides the random-effects structure and therefore the standard
            # errors of every comparison.
            hr(),
            h4("Study Design"),
            radioButtons("harmonic_design", NULL,
                         choices = c("Between-subject (independent groups)" = "between",
                                     "Within-subject (repeated conditions)" = "within",
                                     "Mixed (between x within)" = "mixed"),
                         selected = "between"),
            conditionalPanel(
              condition = "input.harmonic_design == 'between'",
              uiOutput("harmonic_between_var_ui"),
              helpText("Each participant contributes one curve, in one group. Any number of levels.")
            ),
            conditionalPanel(
              condition = "input.harmonic_design == 'within'",
              uiOutput("harmonic_within_var_ui"),
              helpText("Every participant is measured in every condition. Any number of levels.")
            ),
            conditionalPanel(
              condition = "input.harmonic_design == 'mixed'",
              uiOutput("harmonic_between_var_ui2"),
              uiOutput("harmonic_within_var_ui2"),
              helpText("Any size: 2 x 2, 2 x 3, 3 x 4 \u2014 nothing is hard-coded.")
            ),
            uiOutput("harmonic_design_readout"),
            conditionalPanel(
              condition = "input.harmonic_approach == 'two_stage'",
              hr(),
              h4("Uncertainty"),
              checkboxInput("harmonic_bootstrap", "Compute Bootstrap CIs", FALSE),
              conditionalPanel(
                condition = "input.harmonic_bootstrap == true",
                numericInput("harmonic_n_boot", "Bootstrap Iterations (B):", value = 500, min = 100, max = 2000)
              )
            ),
            hr(),
            # ================================================================
            # ADVANCED (P21 phase 4, brief SS13 and SS17)
            # ================================================================
            # Model comparison, bounding and the boundary-fits rule all stay --
            # they are genuinely useful and nothing about them changed. They
            # move BEHIND a collapsed panel because they were sitting in the
            # main column ahead of the Run button, which made model selection
            # look like a step you had to complete before you could fit
            # anything. It is not, and a first-time user should reach Run
            # without reading nine paragraphs.
            tags$details(
              tags$summary(tags$b("Advanced"),
                           style = "cursor:pointer; padding:6px 0; color:#555"),
              div(style = "border-left:2px solid #e3e3e3; padding-left:10px; margin-top:6px",
              conditionalPanel(
              condition = "input.harmonic_approach == 'two_stage'",
              hr(),
              h5("Which fits enter the summaries"),
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
                             solution to average."))
              ),

              conditionalPanel(
              condition = "input.harmonic_approach == 'mixed'",
              hr(),
              # ============================================================
                # DEGREES OF FREEDOM (P21 amendment 5)
                # ============================================================
                # There is NO sample-size theorem here. Kenward-Roger is the
                # better small-sample approximation at every n, and the only
                # thing that changes with n is how long it takes: it refits the
                # reduced model and computes an adjusted covariance, so one KR
                # test costs about as much as the original fit. Benchmarked on
                # 1052 participants with two harmonics: the fit 21.2 s, each KR
                # test 18.6 s, each Satterthwaite test under 0.1 s.
                #
                # So this is a COMPUTATIONAL choice, offered as one, with the
                # measured cost stated. It is not a statistical threshold.
                h5("Degrees of freedom"),
              radioButtons("harmonic_df_method", NULL,
                           choices = c("Kenward-Roger (preferred)" = "kr",
                                       "Satterthwaite (much faster on large samples)" = "satterthwaite"),
                           selected = "kr"),
              helpText(HTML("Kenward\u2013Roger is the better small-sample approximation at
                             <b>every</b> sample size \u2014 there is no n above which
                             Satterthwaite becomes correct instead. What changes with n is the
                             <b>cost</b>: each KR test refits the reduced model, so on ~1000
                             participants with two harmonics it takes about 19 s per block
                             against under 0.1 s for Satterthwaite. Choose Satterthwaite to
                             explore quickly; switch back before you report. The method that
                             ran is always named in the output."))
              ),

              conditionalPanel(
              condition = "input.harmonic_approach == 'two_stage'",
              hr(),
              h5("Model comparison and diagnostics"),
              # AUDIT (P15.2). These are TWO INDEPENDENT diagnostics that sat under
              # one heading with nothing separating them, so the tau field read as
              # if it fed the nested table. It does not: the nested set estimates
              # tau FREELY in its saturating-exponential cells. They are now
              # labelled as separate checks.
              tags$p(tags$b("1. Which specification?"), style = "margin-bottom:4px"),
              checkboxInput("harmonic_model_selection",
                            "Compare nested models (\u0394AICc table)", FALSE),
              # P15.2: rendered from the current setting rather than written out by
              # hand. The old text said "{none, linear, saturating} x {1,2,3} ...
              # 9 fits per subject", which stopped being true when the set began
              # following the harmonic count and gained the logarithmic trend.
              uiOutput("harmonic_model_selection_help"),

              tags$p(tags$b("2. Is \u03c4 identified?"), style = "margin-top:12px; margin-bottom:4px"),
              numericInput("harmonic_tau_fixed",
                           "\u03c4 held at (h), for the free-vs-fixed \u0394AIC:", value = 18,
                           min = 1, max = 72, step = 0.5),
              helpText(HTML("A <b>separate</b> check, and not part of the \u0394AICc table
                             above: that table estimates \u03c4 freely in every
                             saturating-exponential cell. Here the same model is refitted
                             with \u03c4 <i>held</i> at the value above and the two are
                             compared on AIC. Daan, Beersma &amp; Borb\u00e9ly (1984) give
                             \u03c4_rise \u2248 18 h under extended wakefulness. If free
                             \u03c4 does not beat fixed \u03c4, \u03c4 is not identified by
                             these data and a between-group comparison of it is a
                             comparison of where the optimiser stopped on a ridge.
                             <b>Clear the field to skip this check</b> \u2014 the readout
                             then says it was not run, rather than omitting it silently.")),

              hr(),
              h5("Parameter constraints"),
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
                  # tau's range moved to Model Specification, where it applies
                  # whether or not these bounds are switched on
                  helpText(HTML("<small>τ's range is set with the trend under Model
                                 Specification.</small>"))
                ),
                helpText(HTML("<small>Leave blank (NA) for no bound. Defaults based on data range shown in hints above.</small>"))
              )
              ),   # end of the two-stage advanced block
              )
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
                              conditionalPanel(
                                condition = "input.harmonic_approach == 'two_stage'",
                                checkboxInput("harmonic_show_components", "Show Harmonic Components", FALSE),
                                helpText(HTML("<small>Two-stage bands: the group's mean-coefficient
                                               curve &plusmn; z&middot;SE(t), where SE(t) is the standard
                                               error across participants of their own fitted curves.
                                               Pointwise, not simultaneous, and the SAME band
                                               <b>6. Group / Condition Comparison</b> draws.</small>"))),
                              conditionalPanel(
                                condition = "input.harmonic_approach == 'mixed'",
                                helpText(HTML("<small>The cell curves here are the SAME curves as
                                               <b>6. Group / Condition Comparison</b> draws:
                                               one model, one set of trajectories. Choose
                                               <i>All participants</i> above to overlay each
                                               participant's shrunken prediction.</small>")))
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
            # DANCE/server/07_helpers_circular.R and 74_polar_density.R.
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
            # ==========================================================
            # 6. GROUP / CONDITION COMPARISON  (P21 phase 4)
            # ==========================================================
            # Primary result first, decomposition below it, then the pictures,
            # then post-hoc. The order is the order it should be read in: the
            # full-trajectory test is the answer to "do these differ", and
            # everything under it explains why.
            tabPanel("6. Group / Condition Comparison", icon = icon("users"),
                     uiOutput("harmonic_traj_header"),

                     h4("Primary \u2014 full trajectory difference"),
                     helpText(HTML("Does the fitted trajectory differ? Intercept, non-periodic trend
                                    and every harmonic, tested <b>jointly</b>. A group that is
                                    simply higher all day counts as different even when its
                                    amplitude and acrophase are identical.")),
                     uiOutput("harmonic_traj_primary"),

                     hr(),
                     h4("Secondary \u2014 trajectory components"),
                     helpText("Why they differ. Only rows that apply to the fitted model are shown."),
                     uiOutput("harmonic_traj_components"),

                     hr(),
                     h4("Estimated trajectories"),
                     fluidRow(
                       column(4, selectInput("harmonic_traj_component", "Component:",
                                             choices = c("Full fitted trajectory" = "full",
                                                         "Harmonics only (centred on zero)" = "harmonics",
                                                         "Intercept + harmonics (trend removed)" = "baseline_harm",
                                                         "Non-periodic change from origin (first observation)" = "trend"),
                                             selected = "full")),
                       column(4, selectInput("harmonic_traj_band", "Confidence band:",
                                             choices = c("Pointwise" = "pointwise",
                                                         "Simultaneous (Scheffe)" = "simultaneous"),
                                             selected = "pointwise")),
                       column(4, checkboxInput("harmonic_traj_raw", "Overlay raw data", FALSE))
                     ),
                     plotlyOutput("harmonic_traj_curves", height = "440px"),
                     uiOutput("harmonic_traj_band_note"),

                     # The level of each trajectory as numbers: the intercept,
                     # the MESOR and (with a trend) the non-periodic change over
                     # the recording -- per cell, or as pairwise differences with
                     # the same tests as the pairwise table below.
                     h5(tags$b("Intercept, MESOR and change"), style = "margin-top:14px"),
                     radioButtons("harmonic_traj_level_view", NULL,
                                  choices = c("Estimates per cell" = "cells",
                                              "Pairwise differences" = "pairs"),
                                  selected = "cells", inline = TRUE),
                     plotlyOutput("harmonic_traj_levels_plot", height = "340px"),
                     uiOutput("harmonic_traj_levels_note"),

                     # The saturating trend's own parameters per cell: the rise A,
                     # the time constant tau and the initial slope A/tau -- from
                     # the per-cell nonlinear model (mixed) or from the
                     # per-participant fits (two-stage).
                     conditionalPanel(
                       condition = "input.harmonic_trend_type == 'exp_sat'",
                       h5(tags$b("Saturating trend: A, \u03c4 and A/\u03c4"), style = "margin-top:14px"),
                       uiOutput("harmonic_traj_sat_header"),
                       radioButtons("harmonic_traj_sat_view", NULL,
                                    choices = c("Estimates per cell" = "cells",
                                                "Pairwise differences" = "pairs",
                                                "Trend curves" = "curves",
                                                "Profile of \u03c4" = "profile"),
                                    selected = "cells", inline = TRUE),
                       plotlyOutput("harmonic_traj_sat_plot", height = "340px"),
                       uiOutput("harmonic_traj_sat_tests")),

                     hr(),
                     h4("Difference between two cells"),
                     fluidRow(
                       column(5, uiOutput("harmonic_traj_diff_a_ui")),
                       column(5, uiOutput("harmonic_traj_diff_b_ui")),
                       column(2, br(), checkboxInput("harmonic_traj_diff_sim", "Simultaneous", TRUE))
                     ),
                     plotlyOutput("harmonic_traj_diff_plot", height = "340px"),
                     uiOutput("harmonic_traj_diff_note"),

                     hr(),
                     h4("Pairwise comparisons"),
                     fluidRow(
                       column(4, uiOutput("harmonic_traj_effect_ui")),
                       # The choices are replaced from the server once a model is
                       # fitted: the omnibus blocks it can answer depend on the
                       # model, and offering a trend comparison to a model with
                       # no trend is the same mistake the component selector had.
                       column(4, selectInput("harmonic_traj_what", "Compare:",
                                             choices = c("Intercept" = "intercept",
                                                         "Amplitude — H1" = "amplitude1",
                                                         "Acrophase — H1" = "phase1"),
                                             selected = "intercept")),
                       column(4, selectInput("harmonic_traj_adjust", "Multiplicity:",
                                             choices = c("Holm" = "holm",
                                                         "Bonferroni" = "bonferroni",
                                                         "None (family is then wrong)" = "none"),
                                             selected = "holm"))
                     ),
                     uiOutput("harmonic_traj_pairwise"),

            )
          )
        )
      )
