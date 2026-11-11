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
# fixed-tau controls (2.4 / 2.2), and renames the time-averaged midline level bound inputs to what
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
            # THE SCALE OF THE ANALYSIS: a choice the user makes, stated
            # everywhere it applies, never applied behind their back.
            radioButtons("harmonic_scale", "Scale of the analysis:",
                         choices = c("Original values" = "original",
                                     "Natural log, ln(y + c)" = "log"),
                         selected = "original", inline = TRUE),
            conditionalPanel(
              condition = "input.harmonic_scale == 'log'",
              numericInput("harmonic_log_offset", "c, added before the log (0 unless there are zeros):",
                           value = 0, min = 0, step = 0.1),
              helpText(HTML("For positive, right-skewed outcomes such as cortisol or melatonin.
                             The model, its tests, bounds and tables are on the log scale; on the
                             original scale the fitted curve can never fall below &minus;c,
                             the time-averaged midline level becomes a geometric mean (when c = 0) and an amplitude
                             A the peak-to-trough ratio e<sup>2A</sup> of its harmonic."))),
            h5(tags$b("Unsampled intervals")),
            numericInput("harmonic_gap_h", "A gap between sampling times longer than (h):",
                         value = NA, min = 0, step = 0.5),
            helpText(HTML("Blank: twice the median interval between sampling times; 0: never mark.
                           A gap is two consecutive sampling times (columns with at least one
                           observation) further apart than this. The curve is drawn dashed over
                           it, it is shaded, and a time reported there (an acrophase, the curve's
                           minimum) is flagged as extrapolated. The summary names the two
                           columns either side.")),

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
            # THE TIME ORIGIN, t = 0: the first observation by default, or a
            # clock time before it (e.g. wake time). Both the trend and the
            # harmonics are anchored there, so the intercept is the model value
            # at the origin; axes and every acrophase still read in clock time.
            # The period is P, not tau: tau is the saturating trend's time
            # constant, and one symbol for two quantities made "tau held at
            # 18 h" readable as a period.
            numericInput("harmonic_period", "Fundamental period P:", value = 24, min = 1, max = 168, step = 1),
            helpText("Period in the same units as your time variable (e.g., 24 for circadian)."),
            radioButtons("harmonic_origin", "Time origin (t = 0):",
                         choices = c("The first observation" = "first",
                                     "A clock time before it (e.g. wake time)" = "clock"),
                         selected = "first"),
            conditionalPanel(
              condition = "input.harmonic_origin == 'clock'",
              textInput("harmonic_origin_clock", "Origin at (hh:mm):", value = "", placeholder = "e.g. 08:00"),
              helpText(HTML("<small>One origin for everyone: the last time the clock shows this
                             at or before the first observation. The intercept and the trend
                             are anchored there (the intercept is then an <b>extrapolation</b>
                             to a time nobody was sampled at), and a saturating trend's
                             1&minus;e<sup>&minus;t/&tau;</sup> counts from it. The time-averaged
                             levels, the change and the rates are still taken over the sampled
                             window. With a linear or saturating trend the fitted curve, &tau; and the
                             acrophases do not change: only the intercept and the trend's
                             coefficients are re-expressed (A at the origin is
                             A<sub>first</sub>&middot;e<sup>a/&tau;</sup>, a hours earlier). A logarithmic
                             trend, log(1 + t), changes shape with its origin, so its fit changes too.</small>")),
              checkboxInput("harmonic_draw_from_origin",
                            "Draw the fitted curves from the origin (dotted and shaded before the first observation)",
                            TRUE)),
            # THE MODEL STRUCTURE: three analyses, chosen here, carried into
            # the fit, the tests, the plots, the comparison, the report and the
            # export. A rhythm-only model has no trend; a trend-only model has
            # no harmonics (and so no amplitude, acrophase or rhythm test).
            radioButtons("harmonic_structure", "Model structure:",
                         choices = c("Rhythmic component only" = "rhythm",
                                     "Non-periodic trend only" = "trend",
                                     "Rhythmic component and non-periodic trend" = "combined"),
                         selected = "rhythm"),
            conditionalPanel(
              condition = "input.harmonic_structure != 'trend'",
              sliderInput("n_harmonics", "Number of Harmonics:", min = 1, max = 8, value = 1, step = 1),
              helpText("1 = fundamental (P), 2 = adds P/2, 3 = adds P/3, etc."),
              uiOutput("harmonic_warning_ui")),
            conditionalPanel(
              condition = "input.harmonic_structure != 'rhythm'",
              selectInput("harmonic_trend_type", "Non-periodic Trend Model:",
                          choices = c("Linear" = "linear",
                                      "Logarithmic" = "log",
                                      "Saturating exponential" = "exp_sat"),
                          selected = "linear")),
            conditionalPanel(
              condition = "input.harmonic_structure != 'rhythm'",
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
              condition = "input.harmonic_structure != 'rhythm' && input.harmonic_trend_type == 'linear'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·t + harmonics<br><small>Constant rate of
                             change: β is in response units per hour</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_structure != 'rhythm' && input.harmonic_trend_type == 'log'",
              helpText(HTML("Y(t) = β<sub>0</sub> + β·log(t+1) + harmonics<br><small>Rapid initial
                             change, slowing over time</small>"))
            ),
            conditionalPanel(
              condition = "input.harmonic_structure != 'rhythm' && input.harmonic_trend_type == 'exp_sat'",
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
              # The estimate IS the profile likelihood, restricted to the range:
              # one method, named in its label.
              conditionalPanel(
                condition = "input.harmonic_approach == 'mixed'",
                radioButtons("harmonic_traj_tau_mode", "Time constant \u03c4:",
                             choices = c("Hold \u03c4 at a value" = "fixed",
                                         "Estimate one \u03c4 for all cells (profile likelihood), within a range" = "profile"),
                             selected = "fixed")),
              conditionalPanel(
                condition = "input.harmonic_approach == 'two_stage'",
                radioButtons("harmonic_ts_tau_mode", "Time constant \u03c4:",
                             choices = c("Estimate \u03c4 for each participant (nonlinear least squares), within a range" = "estimate",
                                         "Hold \u03c4 at a value for every participant" = "fixed"),
                             selected = "estimate")),
              # Under the mixed approach the value is shown in BOTH modes: when
              # tau is estimated it is what the model holds tau at if the
              # profile is flat, and a value that decides a fit is not hidden.
              conditionalPanel(
                condition = paste("input.harmonic_approach == 'mixed' ||",
                                  "(input.harmonic_approach == 'two_stage' && input.harmonic_ts_tau_mode == 'fixed')"),
                numericInput("harmonic_traj_tau", "\u03c4 held at (h):", value = 18,
                             min = 0.5, max = 168, step = 0.5)),
              conditionalPanel(
                condition = "input.harmonic_approach == 'mixed' && input.harmonic_traj_tau_mode == 'profile'",
                helpText(HTML("<small>Used only if the profile is flat: \u03c4 is then not
                               identified by the data, and the model holds it at this value,
                               which must lie inside the range.</small>"))),
              conditionalPanel(
                condition = paste("(input.harmonic_approach == 'mixed' && (input.harmonic_traj_tau_mode == 'profile' || input.harmonic_traj_sat_cells)) ||",
                                  "(input.harmonic_approach == 'two_stage' && input.harmonic_ts_tau_mode == 'estimate')"),
                fluidRow(
                  column(6, numericInput("harmonic_tau_lo", "Search \u03c4 from (h):", value = 1,
                                         min = 0.1, step = 0.5)),
                  column(6, numericInput("harmonic_tau_hi", "to (h):", value = 72,
                                         min = 0.5, step = 1)))),
              # what is wrong with the setting, as it is typed, not at Run
              uiOutput("harmonic_tau_check"),
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
                               comparison tab.</small>")),
                conditionalPanel(
                  condition = "input.harmonic_traj_sat_cells",
                  checkboxInput("harmonic_traj_sat_hier",
                                "Advanced: let \u03c4 vary between participants too (hierarchical \u03c4)",
                                FALSE),
                  helpText(HTML("<small>Adds a participant-level random effect on τ (on its
                                 carrying scale) to the per-cell model: the population τ per cell,
                                 and how much participants scatter around it. It needs many
                                 participants and a recording that reaches well into the rise;
                                 when it does not converge the next, simpler structure is used
                                 and the summary says so.</small>")))),
              conditionalPanel(
                condition = "input.harmonic_approach == 'two_stage'",
                helpText(HTML("<small><em>Each participant's A and τ are fitted by nonlinear least
                               squares; with τ held, A is a least-squares fit. If fits are poor,
                               try the log trend.</em></small>")))
            ),
            # THE RESIDUALS (mixed approach). Independent is the classical mixed
            # cosinor. Densely sampled curves usually leave residuals that are
            # correlated with their neighbours in time; the continuous-time AR(1)
            # models that, with uneven gaps and missing observations as they are.
            conditionalPanel(
              condition = "input.harmonic_approach == 'mixed'",
              radioButtons("harmonic_traj_resid", "Residuals within a curve:",
                           choices = c("Independent" = "independent",
                                       "Autocorrelated: continuous-time AR(1), estimated" = "car1"),
                           selected = "independent"),
              helpText(HTML("<small><b>Independent</b> is the usual model in sleep and
                             circadian cosinor work, with the participant-level random effects
                             carrying between-person differences. Use the <b>AR(1)</b> for intensive
                             longitudinal sampling (readings close enough in time to share
                             short-term fluctuations), or as a sensitivity analysis: it is the
                             continuous-time form, so uneven or gapped sampling and missing
                             observations are handled as they are, and on an even grid it is the
                             ordinary AR(1).</small>")),
              conditionalPanel(
                condition = "input.harmonic_traj_resid == 'car1'",
                helpText(HTML("<small>Two residuals of one curve Δt hours apart correlate
                               φ<sup>Δt</sup>. φ is estimated by profile likelihood
                               (REML; by maximum likelihood together with τ when τ is
                               estimated) and the model is refitted at φ&#770;; the tests treat
                               φ as known. Estimating φ refits the model about 25 times (and
                               repeats the τ profile when τ is estimated too).</small>")))
            ),
            # TWO-STAGE: the same control, greyed out. Each participant's cosinor
            # has independent residuals; the group comparison uses only the
            # spread of the participants' estimates, which autocorrelation does
            # not bias (Bingham et al. 1982; Cornelissen 2014). Plain disabled
            # radios with no name: not a Shiny input, so nothing can be set here.
            conditionalPanel(
              condition = "input.harmonic_approach == 'two_stage'",
              tags$div(class = "form-group dance-resid-disabled", style = "opacity:0.5;cursor:not-allowed",
                tags$label(class = "control-label", "Residuals within a curve:"),
                tags$div(
                  tags$div(class = "radio", tags$label(
                    tags$input(type = "radio", disabled = NA, checked = NA), tags$span("Independent"))),
                  tags$div(class = "radio", tags$label(
                    tags$input(type = "radio", disabled = NA),
                    tags$span("Autocorrelated: continuous-time AR(1), estimated"))))),
              helpText(HTML("<small>Not available for the two-stage cosinor: each participant is
                             fitted with independent residuals. Autocorrelation within a participant
                             does not bias the group comparison, which uses only the spread of the
                             participants' estimates (Bingham et al., 1982; Cornelissen, 2014); it
                             makes each participant's own p-value and intervals too optimistic. To
                             model it, use the mixed-effects approach.</small>"))),
            conditionalPanel(
              condition = "input.harmonic_structure != 'rhythm'",
              helpText(HTML("<small>β<sub>0</sub> is the <b>intercept</b>: the midline at the
                             first observation, where the trend is zero. The <b>Time-averaged midline level</b> is the
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

              hr(),
              h5("Model comparison and diagnostics"),
              # AUDIT (P15.2). These are TWO INDEPENDENT diagnostics that sat under
              # one heading with nothing separating them, so the tau field read as
              # if it fed the nested table. It does not: the nested set estimates
              # tau FREELY in its saturating-exponential cells. They are now
              # labelled as separate checks.
              # The comparison runs under BOTH approaches: per participant (AICc
              # averaged) for two-stage, by maximum likelihood for mixed-effects.
              tags$p(tags$b("1. Which specification?"), style = "margin-bottom:4px"),
              checkboxInput("harmonic_model_selection",
                            "Compare models across every structure (\u0394AICc table)", FALSE),
              # P15.2: rendered from the current setting rather than written out by
              # hand. The old text said "{none, linear, saturating} x {1,2,3} ...
              # 9 fits per subject", which stopped being true when the set began
              # following the harmonic count and gained the logarithmic trend.
              uiOutput("harmonic_model_selection_help"),
              conditionalPanel(
              condition = "input.harmonic_approach == 'two_stage'",

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
                  condition = "input.harmonic_structure != 'rhythm' && input.harmonic_trend_type == 'exp_sat'",
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
                              checkboxInput("harmonic_show_ci", "Show the 95% confidence bands", TRUE),
                              conditionalPanel(
                                condition = "input.harmonic_scale == 'log'",
                                checkboxInput("harmonic_show_original",
                                              "Show on the original scale (back-transformed, exp(fit) \u2212 c)", TRUE)),
                              # P6.7: off by default. On a few hundred subjects
                              # the points bury the fitted curve they are there
                              # to support.
                              checkboxInput("harmonic_show_data", "Show raw data", FALSE),
                              conditionalPanel(
                                condition = "input.harmonic_show_data == true",
                                radioButtons("harmonic_raw_style", NULL,
                                             choices = c("Points" = "points",
                                                         "Lines (one per curve; broken at missing values)" = "lines"),
                                             selected = "points")),
                              checkboxInput("harmonic_show_components",
                                            "Show the components (each harmonic and the trend)", FALSE),
                              conditionalPanel(
                                condition = "input.harmonic_approach == 'two_stage'",
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
                              uiOutput("harmonic_individual_eq"),
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
                     ),
                     # ONE TRAJECTORY (no between or within factor): its
                     # estimates with intervals, its tests against zero and its
                     # variance components; with factors, a pointer to tab 6
                     uiOutput("harmonic_single_traj")
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

                     # THE STREAMLINED VIEW, the default: the primary test, the
                     # secondary comparisons by question, the rhythm descriptors
                     # and the quality checks, each value with its equation one
                     # click away. Everything else is under Detailed statistics.
                     uiOutput("harmonic_streamlined"),
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
                                                         "Simultaneous (Scheffe)" = "simultaneous",
                                                         "None (curves only)" = "none"),
                                             selected = "pointwise")),
                       column(4, checkboxInput("harmonic_traj_raw", "Overlay raw data", FALSE),
                              conditionalPanel(
                                condition = "input.harmonic_scale == 'log'",
                                checkboxInput("harmonic_traj_original",
                                              "Full trajectory on the original scale (back-transformed)", TRUE)),
                              conditionalPanel(
                                condition = "input.harmonic_traj_raw == true",
                                radioButtons("harmonic_traj_raw_style", NULL,
                                             choices = c("Points" = "points",
                                                         "Lines (one per curve; broken at missing values)" = "lines"),
                                             selected = "points", inline = TRUE)))
                     ),
                     # THE THREE GRAPHS are this plot's three views: the full
                     # fitted trajectory (overall level), the harmonics (the
                     # rhythm) and the non-periodic change from the first
                     # observation (accumulation). The annotations are read off
                     # the same outcome table as the "By question" panel.
                     checkboxGroupInput("harmonic_traj_annot", "Annotate (in the view each belongs to):",
                                        choices = c("Time-averaged fitted response (full view)" = "mean",
                                                    "Total change and average rate (trend view)" = "change",
                                                    "Initial trend rate, as a tangent (trend view)" = "tangent",
                                                    "Asymptote A and the 63.2% point at tau (trend view)" = "asymptote",
                                                    "Acrophase and amplitude (rhythm views)" = "acro"),
                                        selected = c("mean", "change", "acro"), inline = TRUE),
                     plotlyOutput("harmonic_traj_curves", height = "440px"),
                     uiOutput("harmonic_traj_annot_note"),
                     uiOutput("harmonic_traj_band_note"),

                     tags$details(class = "dance-detailed", style = "margin-top:16px",
                       tags$summary(style = "cursor:pointer;font-size:16px;font-weight:600;color:#1f5f7a",
                                    "Detailed statistics"),
                       tags$div(style = "margin-top:10px",
                     h4("Primary \u2014 full trajectory difference"),
                     helpText(HTML("Does the fitted trajectory differ? Intercept, non-periodic trend
                                    and every harmonic, tested <b>jointly</b>. A group that is
                                    simply higher all day counts as different even when its
                                    amplitude and acrophase are identical.")),
                     uiOutput("harmonic_traj_primary"),

                     hr(),
                     # THE THREE QUESTIONS, each with its tests and its group
                     # estimates: overall level, rhythmicity, accumulation. The
                     # same tables go to the CSV and to the APA report.
                     h4("By question \u2014 overall level, rhythmicity, accumulation"),
                     helpText("Why they differ, one scientific question at a time. Primary outcomes
                               are marked; only what the fitted model contains is shown."),
                     uiOutput("harmonic_outcomes"),
                     downloadButton("download_harmonic_outcomes",
                                    "Download outcomes and tests by question (CSV)", class = "btn-sm"),
                     tags$details(style = "margin-top:10px",
                       tags$summary(style = "cursor:pointer;color:#555",
                                    "Technical: every component block test, and where the fitted curve peaks"),
                       uiOutput("harmonic_traj_components")),

                     hr(),

                     # The level of each trajectory as numbers: the intercept,
                     # the time-averaged midline level and (with a trend) the non-periodic change over
                     # the recording -- per cell, or as pairwise differences with
                     # the same tests as the pairwise table below.
                     h5(tags$b("Intercept, midline level and change"), style = "margin-top:14px"),
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
                       condition = "input.harmonic_structure != 'rhythm' && input.harmonic_trend_type == 'exp_sat'",
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
                     downloadButton("download_harmonic_pairwise",
                                    "Download every pairwise table (CSV)",
                                    class = "btn-sm"),
                     )),

            )
          )
        )
      )
