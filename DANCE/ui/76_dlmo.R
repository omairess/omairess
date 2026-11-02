# ==========================================================================
# ui/76_dlmo.R -- dim-light melatonin onset (DLMO)
#
# Works on whatever was selected at Data Import: each row a recording, each
# column a sample. Nothing here assumes a particular protocol, number of
# samples or number of nights.
# ==========================================================================
ui_tab_dlmo <- tabItem(
  tabName = "dlmo",
  fluidRow(
    box(
      title = "DLMO Settings", status = "success", solidHeader = TRUE, width = 4,
      collapsible = TRUE,
      helpText(HTML("Uses the curves and time points selected at <b>Data Import</b>: each row a
                     recording, each column a melatonin sample.")),
      h4("Sampling times"),
      radioButtons("dlmo_time_src", NULL,
                   choices = c("Clock times parsed at import" = "_shared_",
                               "Enter the times" = "_manual_"),
                   selected = "_shared_"),
      conditionalPanel("input.dlmo_time_src == '_manual_'",
        textInput("dlmo_manual_times", "Times (hours or HH:MM, comma-separated):", "",
                  placeholder = "e.g. 18, 18.5, 19, ..., 23:30, 00:00"),
        helpText("A time earlier than the one before it is read as after midnight.")),
      textInput("dlmo_units", "Concentration units:", "pg/mL"),

      hr(),
      h4("Profiles, nights and groups"),
      uiOutput("dlmo_var_ui"),
      numericInput("dlmo_gap_h", "A gap between scheduled samples longer than (h) starts a new profile:",
                   value = 6, min = 1, step = 1),
      helpText(HTML("Two nights in one row (e.g. baseline and sleep deprivation) are split here
                     and numbered Night 1, Night 2. With one row per night, choose the column
                     that names the night.")),

      hr(),
      h4("Methods"),
      checkboxGroupInput("dlmo_methods", NULL,
                         choices = c("Fixed threshold" = "fixed",
                                     "Variable threshold (baseline mean + k SD)" = "variable",
                                     "Hockey-stick (dlmoR package)" = "hockey"),
                         selected = "fixed"),
      uiOutput("dlmo_hockey_status"),
      numericInput("dlmo_threshold", "Threshold (fixed threshold and hockey-stick):",
                   value = 3, min = 0, step = 0.5),
      helpText(HTML("3&ndash;4 pg/mL is usual for saliva, 10 pg/mL for plasma.")),
      conditionalPanel("input.dlmo_methods.indexOf('variable') > -1",
        fluidRow(
          column(6, numericInput("dlmo_n_base", "Baseline samples:", value = 3, min = 2, step = 1)),
          column(6, numericInput("dlmo_k_sd", "k (SD above mean):", value = 2, min = 0, step = 0.5))),
        helpText("The first samples of each profile, which should precede the rise.")),
      conditionalPanel("input.dlmo_methods.indexOf('hockey') > -1",
        numericInput("dlmo_interval_limit", "Hockey-stick: repeated crossings within (h) are one rise:",
                     value = 2, min = 0.25, step = 0.25),
        checkboxInput("dlmo_hockey_fine", "Fine grid search after the coarse one (dlmoR's default)", TRUE),
        helpText(HTML("The hockey-stick takes about 1&ndash;2 minutes per profile with the fine search."))),

      h5(tags$b("Crossing rules (threshold methods)")),
      numericInput("dlmo_confirm", "Samples after the crossing that must also be above the threshold:",
                   value = 1, min = 0, max = 5, step = 1),
      numericInput("dlmo_max_gap", "Interpolate only between samples at most this far apart (h):",
                   value = 2, min = 0.1, step = 0.5),

      h5(tags$b("Limit of quantification")),
      numericInput("dlmo_lloq", "LLOQ (leave empty if none):", value = NA, min = 0, step = 0.1),
      conditionalPanel("input.dlmo_lloq !== null && input.dlmo_lloq !== ''",
        radioButtons("dlmo_lloq_rule", "Values below the LLOQ are set to:",
                     choices = c("LLOQ / 2" = "half", "the LLOQ" = "lloq", "0" = "zero"),
                     selected = "half", inline = TRUE)),

      hr(),
      h4("Habitual bedtime (optional)"),
      radioButtons("dlmo_bed_src", NULL,
                   choices = c("Not known" = "none", "The same for everyone" = "value",
                               "From a column" = "column"),
                   selected = "none"),
      conditionalPanel("input.dlmo_bed_src == 'value'",
        textInput("dlmo_bed_value", "Bedtime (HH:MM):", "23:00")),
      conditionalPanel("input.dlmo_bed_src == 'column'", uiOutput("dlmo_bed_var_ui")),
      helpText(HTML("The phase angle (bedtime &minus; DLMO) is computed only when a bedtime is given.")),

      hr(),
      actionButton("run_dlmo", "Estimate DLMO", icon = icon("play"), class = "btn-success")
    ),
    tabBox(
      title = "DLMO Results", id = "dlmo_tabs", width = 8,
      tabPanel("1. Profiles", icon = icon("chart-line"),
               uiOutput("dlmo_profile_select_ui"),
               plotlyOutput("dlmo_profile_plot", height = "480px"),
               uiOutput("dlmo_profile_note")),
      tabPanel("2. Onsets", icon = icon("table"),
               verbatimTextOutput("dlmo_summary"),
               downloadButton("dlmo_download", "Download onsets (CSV)"),
               br(), br(),
               tableOutput("dlmo_table")),
      tabPanel("3. Comparisons", icon = icon("users"),
               uiOutput("dlmo_compare_ui"))
    )
  )
)
