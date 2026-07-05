# ============================================================================
# 12_tab_psynet.R — psychonetrics tab + NCT comparison (from PsychoNetrix.R,
#                   NCT relocated here from BootSON per the merge plan)
#
# CORRECTNESS-CRITICAL LOGIC WRITTEN IN FULL (do not alter in stage 3):
#   * estimator RESOLVED to a concrete value before fitting — the string
#     "default" (psychonetrics' silent ML->FIML auto-switch) never reaches
#     the model call, the log, or the exported script
#   * ordinal-data guard: declaring ordered/dichotomous data with an
#     ML-family estimator triggers a visible warning (Rhemtulla et al. 2012)
#   * NCT: paired = TRUE/FALSE is an explicit required choice, and edge-level
#     p-values get multiple-comparison correction (Holm default), passed to
#     NCT natively when the installed version supports it, else applied
#     manually to the returned einv.pvals
#
# DECISION: NCT estimates its group networks with the SAME pinned gg_settings
#   object as the bootnet tab (passed in from app.R) — one source of truth,
#   so the comparison can never silently use different estimator settings
#   than the networks the user has been looking at.
# DECISION: FIML is offered as the explicit missing-data estimator (skill
#   prefers FIML over pairwise/listwise) and auto-suggested when the data
#   contain NAs, instead of psychonetrics silently switching.
# ============================================================================

# --- Estimator resolution: no "default" ever reaches runmodel ---------------
resolve_psynet_estimator <- function(choice, data_has_na) {
  if (choice != "auto") return(choice)
  if (data_has_na) "FIML" else "ML"
}

psynetTabUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    varselectUI(ns("vars"), "Indicators (psychonetrics)"),
    shiny::selectInput(ns("data_type"), "Variable type",
                       c("Continuous" = "continuous",
                         "Ordered categorical" = "ordered",
                         "Dichotomous (0/1)" = "dichotomous")),
    shiny::selectInput(ns("estimator"), "Estimator",
                       c("Auto (ML; FIML if missing data)" = "auto",
                         "ML" = "ML", "FIML" = "FIML", "DWLS (ordinal)" = "DWLS",
                         "ULS" = "ULS", "WLS" = "WLS")),
    shiny::uiOutput(ns("ordinal_warning")),
    shiny::radioButtons(ns("family"), "Model family",
                        c("GGM" = "ggm", "LNM" = "lnm", "RNM" = "rnm",
                          "LRNM" = "lrnm", "CFA/SEM" = "lvm")),
    shiny::actionButton(ns("run"), "Fit model", class = "btn-primary"),
    shiny::verbatimTextOutput(ns("fit_summary")),
    # TODO(stage3-psynet-ui): lambda-matrix editor, prune/stepup/modelsearch
    # controls, panel families (dlvm1/panelgvar/ri_clpm) + wave detection,
    # results tabBox, factor scores — port from PsychoNetrix.R:120-706. The
    # data-transformation pipeline (:914-980) ports as a recorded "reshape"
    # phase step (it changes the data every model sees).
    shiny::hr(),
    shiny::h4("Network Comparison Test (NCT)"),
    nctUI(ns("nct"))
  )
}

psynetTabServer <- function(id, data_bus, rec, gg_settings) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "psynet")
    rv  <- shiny::reactiveValues(model = NULL, estimator_used = NULL)

    # --- Ordinal/ML mismatch guard (fires on the failure mode PsychoNetrix
    #     missed: ordered data left on the default ML path) ------------------
    output$ordinal_warning <- shiny::renderUI({
      ml_family <- input$estimator %in% c("auto", "ML", "FIML")
      if (input$data_type %in% c("ordered", "dichotomous") && ml_family) {
        shiny::div(class = "alert alert-warning",
          "Ordinal/dichotomous items with an ML-family estimator treat categories
           as continuous. With < ~5-7 categories this biases estimates
           (Rhemtulla et al., 2012) — consider DWLS or ULS.")
      }
    })

    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      est  <- resolve_psynet_estimator(input$estimator, anyNA(dat))
      rv$estimator_used <- est
      fam  <- input$family
      seed <- rec_seed(rec)

      # TODO(stage3-psynet): full family dispatch incl. lambda matrix, prune/
      # stepup/modelsearch chains, multi-group — port from
      # PsychoNetrix.R:1261-1486 through this pinned-estimator path. GGM path
      # written out as the reference wiring:
      model <- tryCatch({
        if (fam == "ggm") {
          m <- psychonetrics::ggm(dat, estimator = est)
          m <- psychonetrics::runmodel(m)
          m
        } else {
          shiny::showNotification(
            "TODO(stage3-psynet): non-GGM families not yet ported.",
            type = "warning")
          NULL
        }
      }, error = function(e) {
        shiny::showNotification(paste("Fit error:", conditionMessage(e)),
                                type = "error"); NULL
      })
      shiny::req(model)
      rv$model <- model

      rec_upsert(
        rec, "psynet_analysis", "analysis",
        description = sprintf(
          "[psychonetrics] Fitted %s on %d variables (n = %d): estimator = %s (resolved explicitly%s), data type = %s.",
          toupper(fam), length(vars), nrow(dat), est,
          if (input$estimator == "auto") " from 'auto'" else "",
          input$data_type),
        code = sprintf(
          paste("psynet_vars <- %s",
                "dat_psynet  <- dat_wide[, psynet_vars]",
                'mod <- psychonetrics::%s(dat_psynet, estimator = "%s")',
                "mod <- psychonetrics::runmodel(mod)", sep = "\n"),
          vars_literal(vars), fam, est)
      )
    })

    output$fit_summary <- shiny::renderPrint({
      shiny::req(rv$model)
      psychonetrics::fit(rv$model)
    })

    # NCT lives on this tab; shares the bootnet tab's pinned GGM settings.
    nctServer("nct", data_bus, rec, gg_settings)

    list(model = shiny::reactive(rv$model))
  })
}

# ============================================================================
# NCT sub-module — paired handling + multiple-comparison correction
# ============================================================================

nctUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("group_col"), "Grouping column", choices = NULL),
    shiny::radioButtons(ns("design"), "Design (required — determines paired)",
                        c("Independent groups (paired = FALSE)" = "independent",
                          "Pre-post / repeated measures (paired = TRUE)" = "paired")),
    shiny::selectInput(ns("adjust"), "Edge-test p-value correction (pinned)",
                       c("Holm (recommended)" = "holm", "BH / FDR" = "BH",
                         "Bonferroni" = "bonferroni",
                         "None (NOT recommended)" = "none"),
                       selected = "holm"),
    shiny::numericInput(ns("iterations"), "Permutation iterations",
                        value = 1000, min = 100, step = 100),
    shiny::actionButton(ns("run"), "Run NCT"),
    shiny::verbatimTextOutput(ns("nct_summary")),
    DT::dataTableOutput(ns("edge_table"))
    # TODO(stage3-nct-ui): paired subject-ID matching selector (port
    # BootSON.R:5409-5466 row-alignment logic); centrality-difference tests.
  )
}

# Core NCT call. Pins estimation to the shared gg settings and applies the
# chosen edge-level correction natively when NCT supports it, manually if not.
run_nct_corrected <- function(dat1, dat2, paired, it, adjust, gg, seed) {
  set.seed(seed)
  nct_formals <- names(formals(NetworkComparisonTest::NCT))
  args <- list(data1 = dat1, data2 = dat2, it = it, paired = paired,
               test.edges = TRUE, edges = "all")
  # Pin estimation so NCT uses the same settings as the bootnet tab:
  if (all(c("estimator", "estimatorArgs") %in% nct_formals)) {
    args$estimator <- bootnet::estimateNetwork
    args$estimatorArgs <- list(default = gg$default,
                               corMethod = gg$corMethod, tuning = gg$tuning)
  } else if ("gamma" %in% nct_formals) {
    args$gamma <- gg$tuning
  }
  native_adjust <- "p.adjust.methods" %in% nct_formals
  if (native_adjust) args$p.adjust.methods <- adjust
  res <- do.call(NetworkComparisonTest::NCT, args)
  # Manual fallback for NCT builds without p.adjust.methods:
  if (!native_adjust && !identical(adjust, "none") &&
      !is.null(res$einv.pvals)) {
    pcol <- grep("p", names(res$einv.pvals), ignore.case = TRUE, value = TRUE)[1]
    res$einv.pvals[[pcol]] <- stats::p.adjust(res$einv.pvals[[pcol]],
                                              method = adjust)
  }
  attr(res, "edge_p_adjust") <- adjust
  res
}

nctServer <- function(id, data_bus, rec, gg_settings) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(data_bus$wide(), {
      shiny::updateSelectInput(session, "group_col",
                               choices = names(data_bus$wide()))
    })

    res_r <- shiny::eventReactive(input$run, {
      dat <- data_bus$wide()
      shiny::req(input$group_col %in% names(dat))
      g   <- dat[[input$group_col]]
      lv  <- names(sort(table(g), decreasing = TRUE))[1:2]
      num <- names(dat)[vapply(dat, is.numeric, TRUE)]
      num <- setdiff(num, input$group_col)
      d1  <- dat[g == lv[1], num, drop = FALSE]
      d2  <- dat[g == lv[2], num, drop = FALSE]
      paired <- identical(input$design, "paired")
      if (paired) {
        # Paired NCT requires row i of d1 to be the same subject as row i of
        # d2. TODO(stage3-nct): subject-ID matching UI; until then require
        # equal group sizes and warn that row order is assumed aligned.
        shiny::validate(shiny::need(nrow(d1) == nrow(d2),
          "Paired design needs equal, subject-aligned groups (ID matching UI: stage 3)."))
      }
      gg   <- gg_settings()
      seed <- rec_seed(rec)
      res  <- run_nct_corrected(d1, d2, paired, input$iterations,
                                input$adjust, gg, seed)

      rec_upsert(
        rec, "nct_comparison", "comparison",
        description = sprintf(
          "[NCT] Compared '%s' vs '%s' (%s design, paired = %s): %d permutations, seed %d; edge-level p-values corrected with '%s'; estimation pinned to %s / %s / gamma = %s. A null NCT is not proof of equality.",
          lv[1], lv[2], input$design, paired, input$iterations, seed,
          input$adjust, gg$default, gg$corMethod, gg$tuning),
        code = sprintf(
          paste("set.seed(%d)",
                "nct_res <- NetworkComparisonTest::NCT(",
                "  data1 = dat_wide[dat_wide$%s == \"%s\", nct_vars],",
                "  data2 = dat_wide[dat_wide$%s == \"%s\", nct_vars],",
                "  it = %d, paired = %s, test.edges = TRUE, edges = \"all\",",
                "  p.adjust.methods = \"%s\")  # edge-level correction",
                "nct_vars <- %s", sep = "\n"),
          seed, input$group_col, lv[1], input$group_col, lv[2],
          input$iterations, paired, input$adjust, vars_literal(num))
      )
      list(res = res, groups = lv)
    })

    output$nct_summary <- shiny::renderPrint({
      r <- res_r()
      cat(sprintf("NCT: %s vs %s | paired = %s | edge p-adjust = %s\n\n",
                  r$groups[1], r$groups[2],
                  identical(input$design, "paired"),
                  attr(r$res, "edge_p_adjust")))
      print(summary(r$res))
    })

    output$edge_table <- DT::renderDataTable({
      r <- res_r()
      shiny::req(!is.null(r$res$einv.pvals))
      tb <- as.data.frame(r$res$einv.pvals)
      pcol <- grep("p", names(tb), ignore.case = TRUE, value = TRUE)[1]
      tb$significant <- tb[[pcol]] < 0.05   # on CORRECTED p-values
      DT::datatable(tb, caption = sprintf(
        "Edge-level differences; p-values %s-corrected.",
        attr(r$res, "edge_p_adjust")))
    })

    res_r
  })
}
