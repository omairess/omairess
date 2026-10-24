# ============================================================================
# app.R — CoDA sleep-diary app: multilevel compositional data analysis
# (multilevelcoda; Le, Stanford, Dumuid & Wiley, 2025) for sleep diaries.
#
# Tabs: Data -> Diary & composition -> Model -> Substitution -> Export.
# One step recorder feeds both the pipeline log and the reproducible script.
# Run with:  shiny::runApp("multilevelcoda-app")
# ============================================================================

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

ensure_packages()   # installs missing packages only \u2014 never auto-updates

HELPERS_PATH <- normalizePath(file.path("R", "02_diary_helpers.R"))

about_panel <- shiny::fluidPage(
  shiny::h3("What the app does"),
  shiny::tags$ol(
    shiny::tags$li("Reads one or many diary files (Excel, CSV, JSON, SPSS, RDS) in long, items-in-rows or wide layout and stacks them into one row per participant-night."),
    shiny::tags$li("Derives time in bed (TIB), total sleep time (TST), sleep-onset latency (SOL), wake after sleep onset (WASO), early-morning awakening (EMA) and sleep efficiency per night using Consensus Sleep Diary definitions."),
    shiny::tags$li("Builds a composition, e.g. the 24-h day as Sleep / Awake in bed / Out of bed, and computes between- and within-person isometric log-ratio (ilr) coordinates with multilevelcoda::complr()."),
    shiny::tags$li("Fits Bayesian multilevel models with multilevelcoda::brmcoda() (brms / Stan), with the composition as predictor of e.g. QS, QR or FJ, or as outcome."),
    shiny::tags$li("Estimates what reallocating minutes between parts would do to the outcome with multilevelcoda::substitution()."),
    shiny::tags$li("Exports a stand-alone R script, the pipeline, the data and all results.")),
  shiny::h3("Points to keep in mind"),
  shiny::tags$ul(
    shiny::tags$li("Parts must be strictly positive: log-ratios of 0 are undefined. Zeros (often WASO or SOL) are replaced multiplicatively or the night is dropped."),
    shiny::tags$li("Between-person effects compare people's usual compositions; within-person effects compare a person's nights with their own average. With few participants, between-person estimates are imprecise."),
    shiny::tags$li("A single patient's diary gives a single-level model of their nights. Nights are serially dependent, so treat those estimates as a description of that patient."),
    shiny::tags$li("Substitution estimates are model-based contrasts, not causal effects of changing behaviour, unless the design supports that."),
    shiny::tags$li("The model uses brms default priors; they are printed in the exported script (prior_summary).")),
  shiny::h3("References"),
  shiny::tags$ul(
    shiny::tags$li("Le, F., Stanford, T. E., Dumuid, D., & Wiley, J. F. (2025). Bayesian multilevel compositional data analysis: Introduction, evaluation, and application. Psychological Methods. doi:10.1037/met0000750"),
    shiny::tags$li("Le, F., Dumuid, D., Stanford, T. E., & Wiley, J. F. (2025). Bayesian multilevel compositional data analysis with the R package multilevelcoda. Multivariate Behavioral Research. doi:10.1080/00273171.2025.2565598"),
    shiny::tags$li("Carney, C. E., Buysse, D. J., Ancoli-Israel, S., Edinger, J. D., Krystal, A. D., Lichstein, K. L., & Morin, C. M. (2012). The Consensus Sleep Diary: Standardizing prospective sleep self-monitoring. Sleep, 35(2), 287\u2013302. doi:10.5665/sleep.1642"),
    shiny::tags$li("Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society: Series B, 44(2), 139\u2013177."),
    shiny::tags$li("Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcel\u00f3-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. Mathematical Geology, 35(3), 279\u2013300."),
    shiny::tags$li("Mart\u00edn-Fern\u00e1ndez, J. A., Barcel\u00f3-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. Mathematical Geology, 35(3), 253\u2013278."),
    shiny::tags$li("Dumuid, D., Pedi\u0161i\u0107, \u017d., Stanford, T. E., et al. (2019). The compositional isotemporal substitution model: A method for estimating changes in a health outcome for reallocation of time between sleep, physical activity and sedentary behaviour. Statistical Methods in Medical Research, 28(3), 846\u2013857."),
    shiny::tags$li("B\u00fcrkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan. Journal of Statistical Software, 80(1), 1\u201328. doi:10.18637/jss.v080.i01"))
)

ui <- shiny::navbarPage(
  title = "CoDA sleep diaries",
  windowTitle = "CoDA sleep diaries \u2014 multilevelcoda",
  shiny::tabPanel("1 \u00b7 Data", dataModuleUI("data")),
  shiny::tabPanel("2 \u00b7 Diary & composition", diaryModuleUI("diary")),
  shiny::tabPanel("3 \u00b7 Model", modelModuleUI("model")),
  shiny::tabPanel("4 \u00b7 Substitution", substitutionModuleUI("sub")),
  shiny::tabPanel("5 \u00b7 Pipeline & export", exportsModuleUI("export")),
  shiny::tabPanel("About & references", about_panel)
)

server <- function(input, output, session) {
  rec      <- new_recorder()
  data_bus <- dataModuleServer("data", rec)
  comp     <- diaryModuleServer("diary", data_bus, rec)
  model    <- modelModuleServer("model", comp, rec)
  sub      <- substitutionModuleServer("sub", model, rec)
  exportsModuleServer("export", rec, data_bus, comp, model, sub, HELPERS_PATH)
}

shiny::shinyApp(ui, server)
