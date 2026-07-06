# ============================================================================
# app.R — merged network-psychometrics app (BootSON + DAGger + PsychoNetrix)
#
# Assembly skeleton (stage 2). One shared recorder, one data contract, three
# namespaced analysis tabs. Governing conventions: shiny-house-style (ten
# rules) + network-psychometrics-correctness.
#
# DECISION: the three source apps become three shiny::moduleServer tabs —
#   NS() namespacing structurally kills the input$run / input$nodeBorderColor
#   / input$layout ID collisions found in stage 1, and lets each tab keep its
#   own layout vocabulary (spring/circle/groups/ega/pca vs spring/circle/
#   cascade) without semantic clashes.
# DECISION: one app-wide seed input, stored in the recorder and re-emitted
#   before every stochastic step in both the app and the exported script.
# ============================================================================

for (f in sort(list.files("R", full.names = TRUE))) source(f)

ensure_packages()   # install missing only — never auto-update (rule 0)

ui <- shiny::navbarPage(
  "Network Psychometrics Workbench",
  shiny::tabPanel("Data",       dataModuleUI("data"),
                  shiny::numericInput("seed", "Master seed (all bootstraps)",
                                      value = 12345, min = 1)),
  shiny::tabPanel("PMRFs (bootnet)",       bootnetTabUI("bootnet")),
  shiny::tabPanel("DAGs (bnlearn)",        daggerTabUI("dag")),
  shiny::tabPanel("LNMs (psychonetrics)",  psynetTabUI("psynet")),
  shiny::tabPanel("Compare networks (NCT)", nctUI("nct")),
  shiny::tabPanel(
    "Pipeline & Export",
    shiny::h4("Analysis pipeline (house rule 8)"),
    shiny::verbatimTextOutput("pipeline_log"),
    shiny::downloadButton("dl_script", "Download reproducible R script"),
    shiny::downloadButton("dl_bundle", "Download data + results bundle (.rds)"),
    shiny::helpText("Per-figure PNG/PDF/SVG export is on each analysis tab",
                    "(next to its colour/size controls).")
    # TODO(stage3-exports): CSV exports of edge lists / centrality / arc
    # tables (beyond the .rds bundle above) — not yet implemented.
  )
)

server <- function(input, output, session) {

  rec  <- new_recorder(seed = 12345L)
  vers <- record_setup_step(rec)

  shiny::observeEvent(input$seed, {
    rec_set_seed(rec, input$seed)
  }, ignoreInit = TRUE)

  data_bus <- dataModuleServer("data", rec)

  bootnet_out <- bootnetTabServer("bootnet", data_bus, rec)
  dag_out     <- daggerTabServer("dag", data_bus, rec)
  psynet_out  <- psynetTabServer("psynet", data_bus, rec)
  # NCT is now its OWN top-level tab. It still consumes the bootnet tab's GGM
  # settings reactive — one source of truth for estimator/corMethod/gamma, so
  # the comparison uses the same settings as the networks on the GGM tab.
  nctServer("nct", data_bus, rec, gg_settings = bootnet_out$settings)

  # --- Rule 8: live pipeline log -------------------------------------------
  output$pipeline_log <- shiny::renderText({
    rec$rev()                      # re-render on every recorder change
    build_pipeline_log(rec)
  })

  # --- Rule 6: standalone script, always current ---------------------------
  output$dl_script <- shiny::downloadHandler(
    filename = function() sprintf("analysis_%s.R", format(Sys.Date())),
    content = function(file) {
      writeLines(build_script(rec, APP_PACKAGES, vers), file)
    }
  )

  # --- Rule 7: everything-bundle -------------------------------------------
  output$dl_bundle <- shiny::downloadHandler(
    filename = function() sprintf("bundle_%s.rds", format(Sys.Date())),
    content = function(file) {
      saveRDS(list(
        raw_data     = tryCatch(data_bus$raw(),  error = function(e) NULL),
        wide_data    = tryCatch(data_bus$wide(), error = function(e) NULL),
        meta         = tryCatch(data_bus$meta(), error = function(e) NULL),
        bootnet_net  = bootnet_out$net(),
        dag_avg      = dag_out$avg(),
        psynet_model = psynet_out$model(),
        pipeline     = build_pipeline_log(rec),
        script       = build_script(rec, APP_PACKAGES, vers),
        pkg_versions = vers,
        seed         = rec_seed(rec),
        timestamp    = Sys.time()
      ), file)
    }
  )
}

shiny::shinyApp(ui, server)
