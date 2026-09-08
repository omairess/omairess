# ==============================================================================
# ui/55_mixed.R — mixed (between x within) designs
# ==============================================================================
# The two existing fANOVA designs are one-factor. This tab is for a design with
# BOTH: one factor varying between subjects and one varying within them, and an
# interaction between the two. See server/06_helpers_mixed.R for why neither
# existing kernel can stand in for it.
# ==============================================================================

ui_tab_mixed <- tabItem(
  tabName = "mixed",
  fluidRow(
    box(
      title = "Mixed design (between x within)", status = "primary",
      solidHeader = TRUE, width = 4,

      helpText(HTML(
        "For a design where one factor varies <b>between</b> subjects and another
         <b>within</b> them &mdash; e.g. group as between, condition as within.
         Needs a participant identifier so the repeated measurement can be
         modelled rather than treated as extra independent curves.")),

      uiOutput("mixed_subject_ui"),
      uiOutput("mixed_between_ui"),
      uiOutput("mixed_within_ui"),
      hr(),

      radioButtons("mixed_analysis", "Analysis:",
                   choices = list(
                     "Functional: smooth per cell + subject random smooth" = "fanova",
                     "Cosinor: mixed harmonic regression" = "cosinor"),
                   selected = "fanova"),

      conditionalPanel(
        condition = "input.mixed_analysis == 'fanova'",
        sliderInput("mixed_k_time", "Basis size for the per-cell smooth (k):",
                    min = 4, max = 20, value = 12, step = 1),
        sliderInput("mixed_k_subject", "Basis size for each subject's curve (k):",
                    min = 3, max = 12, value = 6, step = 1),
        helpText(HTML(
          "The subject term is a nuisance: give it the same flexibility as the
           cell smooths and it will absorb the effect you are testing. Keep it
           coarser."))
      ),

      conditionalPanel(
        condition = "input.mixed_analysis == 'cosinor'",
        numericInput("mixed_period", "Period (time units):", value = 24,
                     min = 0.1, step = 0.5),
        sliderInput("mixed_harmonics", "Harmonics:", min = 1, max = 3, value = 1, step = 1),
        helpText(HTML(
          "Fits <code>y ~ (cos + sin) * between * within + (1 + cos + sin | subject)</code>:
           a MESOR, amplitude and acrophase per cell, and a rhythm per subject.
           This is a genuine mixed model, not the two-stage route on the
           Cosinor tabs, which cannot respect pairing."))
      ),

      helpText(HTML(
        "<b>These models are fitted to the raw observations</b>, not to the
         smoothed curves from the Preprocessing tab. That is deliberate: the
         model estimates the temporal structure itself, so feeding it
         pre-smoothed data would smooth twice and distort the residual model.
         Re-smoothing therefore does not invalidate a mixed fit, and the app
         does not clear one when you re-smooth.")),

      checkboxInput("mixed_real_time",
                    "Use real elapsed clock time as the time axis", TRUE),
      helpText(HTML(
        "Leave this on unless your columns are evenly spaced. A smooth of the
         column index on unevenly spaced data is a smooth of the wrong variable.")),

      hr(),
      actionButton("run_mixed", "Run mixed analysis",
                   icon = icon("play"), class = "btn-primary"),
      br(), br(),
      uiOutput("mixed_design_summary")
    ),

    box(
      title = "Results", status = "info", solidHeader = TRUE, width = 8,
      verbatimTextOutput("mixed_results"),
      hr(),
      plotlyOutput("mixed_plot", height = "420px")
    )
  )
)
