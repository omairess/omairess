# ============================================================================
# 10_tab_bootnet.R — bootnet GGM tab (from BootSON.R)
#
# CORRECTNESS-CRITICAL LOGIC WRITTEN IN FULL (do not alter in stage 3):
#   * estimator / corMethod / EBIC-gamma pinned, corMethod default cor_auto
#   * nonparametric + case-dropping bootstraps run UNCONDITIONALLY with
#     estimation (BootSON's opt-in checkbox is gone)
#   * CS-coefficient gate: CS < 0.25 BLOCKS the ordered centrality plot;
#     values remain visible only unordered (alphabetical) under a red banner
#   * expected influence force-included and listed first whenever the
#     estimated network contains any negative edge
#
# DECISION: bootstraps are not opt-in — "Estimate & validate" runs
#   estimate -> nonparametric boot -> case-dropping boot as one recorded
#   pipeline, because the correctness skill makes them mandatory for any
#   edge/centrality interpretation.
# DECISION: gate is PER STATISTIC (corStability returns one CS per measure);
#   a stable strength does not license interpreting an unstable betweenness.
# DECISION: centrality (strength, EI) is computed transparently from the
#   weights matrix (colSums), not via a helper whose defaults may drift.
# DECISION: corMethod default is "cor_auto" (polychoric auto-detect), fixing
#   BootSON's Pearson-by-default bug for Likert-type items.
# ============================================================================

# Per-statistic interpretation gate (Epskamp, Borsboom & Fried, 2018).
cs_gate <- function(cs) {
  if (!is.finite(cs)) "unknown"
  else if (cs < 0.25) "blocked"
  else if (cs < 0.50) "caution"
  else "ok"
}

bootnetTabUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    varselectUI(ns("vars"), "Node set (bootnet GGM)"),
    shiny::selectInput(ns("bootnet_default"), "Estimator (pinned)",
                       c("EBICglasso", "ggmModSelect", "pcor", "cor")),
    shiny::selectInput(ns("bootnet_cormethod"), "Correlation method (pinned)",
                       c("cor_auto (polychoric auto-detect)" = "cor_auto",
                         "spearman" = "spearman",
                         "cor (Pearson — continuous data only)" = "cor"),
                       selected = "cor_auto"),
    shiny::numericInput(ns("bootnet_tuning"), "EBIC gamma (pinned)",
                        value = 0.5, min = 0, max = 1, step = 0.05),
    shiny::numericInput(ns("nboots"), "Bootstrap samples (both types)",
                        value = 1000, min = 250, step = 250),
    shiny::actionButton(ns("run"), "Estimate & validate",
                        class = "btn-primary"),
    shiny::helpText("Runs estimation + nonparametric bootstrap (edge CIs) +",
                    "case-dropping bootstrap (CS-coefficients) as one step.",
                    "Centrality stays locked until all three complete."),
    shiny::hr(),
    shiny::uiOutput(ns("stability_banner")),
    shiny::fluidRow(
      shiny::column(8,
        shiny::plotOutput(ns("network_plot")),
        shiny::plotOutput(ns("edge_ci_plot")),
        shiny::uiOutput(ns("centrality_ui"))
      ),
      shiny::column(4, appearanceUI(ns("look"), signed = TRUE))
    )
    # TODO(stage3-bootnet-ui): network sub-tabs; port Bridge-Symptoms + NSON
    # sub-tabs from BootSON.R:1633-2224 unchanged (Bridge already has its own
    # case-drop + corStability — keep it).
  )
}

bootnetTabServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "bootnet")

    # Pinned GGM settings: SINGLE SOURCE, exported to the NCT module via the
    # return value so the comparison can never drift from this tab.
    gg_settings <- shiny::reactive({
      list(
        default   = input$bootnet_default   %||% "EBICglasso",
        corMethod = input$bootnet_cormethod %||% "cor_auto",
        tuning    = input$bootnet_tuning    %||% 0.5
      )
    })
    rv  <- shiny::reactiveValues(net = NULL, boot_np = NULL, boot_cd = NULL,
                                 cs = NULL, has_neg = FALSE, plot_fn = NULL)

    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      s    <- gg_settings()
      nb   <- max(250L, as.integer(input$nboots))
      if (nb < 500) shiny::showNotification(
        "nBoots < 500: stability estimates will be rough; 1000 recommended.",
        type = "warning")
      seed <- rec_seed(rec)

      shiny::withProgress(message = "Estimating + validating network", {

        # -- 1. Estimation, every setting pinned explicitly ------------------
        shiny::incProgress(0.1, detail = "estimateNetwork")
        net <- bootnet::estimateNetwork(
          dat,
          default   = s$default,
          corMethod = s$corMethod,
          tuning    = s$tuning
        )
        wmat <- qgraph::getWmat(net)
        rv$has_neg <- any(wmat[upper.tri(wmat)] < 0)

        rec_upsert(
          rec, "bootnet_analysis", "analysis",
          description = sprintf(
            "[bootnet] Estimated GGM on %d variables (n = %d): default = %s, corMethod = %s, EBIC gamma = %s.%s",
            length(vars), nrow(dat), s$default, s$corMethod, s$tuning,
            if (rv$has_neg) " Network contains negative edges -> expected influence reported first."
            else ""),
          code = sprintf(
            paste("bootnet_vars <- %s",
                  "dat_bootnet  <- dat_wide[, bootnet_vars]",
                  'net <- bootnet::estimateNetwork(dat_bootnet,',
                  '  default   = "%s",',
                  '  corMethod = "%s",',
                  '  tuning    = %s)', sep = "\n"),
            vars_literal(vars), s$default, s$corMethod, s$tuning)
        )

        # -- 2. Nonparametric bootstrap: edge-weight accuracy (mandatory) ----
        shiny::incProgress(0.3, detail = "nonparametric bootstrap (edge CIs)")
        set.seed(seed)
        boot_np <- bootnet::bootnet(net, nBoots = nb, type = "nonparametric")

        # -- 3. Case-dropping bootstrap: CS-coefficients (mandatory) ---------
        # EI included in the case bootstrap whenever negative edges exist.
        cd_stats <- c("strength",
                      if (rv$has_neg) "expectedInfluence",
                      "closeness", "betweenness")
        shiny::incProgress(0.3, detail = "case-dropping bootstrap (CS)")
        set.seed(seed)
        boot_cd <- bootnet::bootnet(net, nBoots = nb, type = "case",
                                    statistics = cd_stats)
        cs <- bootnet::corStability(boot_cd)

        rec_upsert(
          rec, "bootnet_stability", "stability",
          description = sprintf(
            "[bootnet] Ran %d nonparametric (edge CIs) and %d case-dropping bootstraps (seed %d); CS-coefficients: %s. CS < 0.25 blocks centrality-order interpretation; >= 0.5 preferred.",
            nb, nb, seed,
            paste(sprintf("%s = %.2f", names(cs), cs), collapse = ", ")),
          code = sprintf(
            paste("set.seed(%d)",
                  'boot_np <- bootnet::bootnet(net, nBoots = %d, type = "nonparametric")',
                  "set.seed(%d)",
                  'boot_cd <- bootnet::bootnet(net, nBoots = %d, type = "case",',
                  "  statistics = %s)",
                  "cs <- bootnet::corStability(boot_cd)",
                  "# Overlapping edge CIs in plot(boot_np) mean edge ORDER is not interpretable.",
                  sep = "\n"),
            seed, nb, seed, nb, vars_literal(cd_stats))
        )

        rv$net <- net; rv$boot_np <- boot_np; rv$boot_cd <- boot_cd; rv$cs <- cs
      })
    })

    # --- Transparent centrality from the weights matrix ---------------------
    centrality_tbl <- shiny::reactive({
      shiny::req(rv$net)
      w  <- qgraph::getWmat(rv$net)
      tb <- data.frame(
        node               = colnames(w),
        strength           = colSums(abs(w)),
        expected_influence = colSums(w),
        row.names = NULL
      )
      # DECISION: EI column listed first when negative edges exist — strength
      # can mislead there because opposite-sign edges inflate it.
      if (rv$has_neg) tb[, c("node", "expected_influence", "strength")] else tb
    })

    # --- The CS gate ---------------------------------------------------------
    # Gate the ORDERED display on the CS of the leading (interpreted) measure.
    lead_cs <- shiny::reactive({
      shiny::req(rv$cs)
      key <- if (rv$has_neg) "expectedInfluence" else "strength"
      if (key %in% names(rv$cs)) unname(rv$cs[[key]]) else NA_real_
    })

    output$stability_banner <- shiny::renderUI({
      shiny::req(rv$cs)
      g <- cs_gate(lead_cs())
      msg <- sprintf(
        "CS-coefficient (%s) = %.2f. %s",
        if (rv$has_neg) "expected influence" else "strength", lead_cs(),
        switch(g,
          blocked = "BELOW 0.25 — centrality order is NOT interpretable at this sample size; ordered centrality output is disabled.",
          caution = "Between 0.25 and 0.5 — interpret centrality order with caution.",
          ok      = "At or above 0.5 — centrality order interpretable.",
          unknown = "Could not be computed."))
      shiny::div(class = switch(g, blocked = "alert alert-danger",
                                caution = "alert alert-warning",
                                "alert alert-success"), msg)
    })

    output$centrality_ui <- shiny::renderUI({
      ns <- session$ns
      shiny::req(rv$net)
      shiny::validate(shiny::need(
        !is.null(rv$cs),
        "Centrality locked: the mandatory case-dropping bootstrap has not completed."))
      if (cs_gate(lead_cs()) == "blocked") {
        # BLOCK ordered interpretation: values shown alphabetically only,
        # no ordered plot, explicit refusal text — a gate, not a footnote.
        shiny::tagList(
          shiny::div(class = "alert alert-danger",
            sprintf("CS = %.2f < 0.25: the ordered centrality plot is disabled (Epskamp, Borsboom & Fried, 2018). Values below are in ALPHABETICAL order and their ranking must not be interpreted.",
                    lead_cs())),
          shiny::tableOutput(ns("cent_table_unordered"))
        )
      } else {
        shiny::tagList(
          if (cs_gate(lead_cs()) == "caution") shiny::div(
            class = "alert alert-warning",
            "CS between 0.25 and 0.5 — interpret ordering cautiously."),
          shiny::plotOutput(ns("cent_plot")),
          shiny::tableOutput(ns("cent_table_ordered"))
        )
      }
    })

    output$cent_table_unordered <- shiny::renderTable({
      tb <- centrality_tbl()
      tb[order(tb$node), ]                      # alphabetical — never by value
    })
    output$cent_table_ordered <- shiny::renderTable({
      tb <- centrality_tbl()
      tb[order(-tb[[2]]), ]                     # by leading measure, gate passed
    })
    output$cent_plot <- shiny::renderPlot({
      shiny::req(cs_gate(lead_cs()) != "blocked")
      tb  <- centrality_tbl()
      tb  <- tb[order(tb[[2]]), ]
      graphics::dotchart(tb[[2]], labels = tb$node,
                         xlab = names(tb)[2], pch = 19)
    })

    output$edge_ci_plot <- shiny::renderPlot({
      shiny::req(rv$boot_np)
      # bootnet's CI plot; overlapping CIs = edge order not interpretable.
      plot(rv$boot_np, labels = FALSE, order = "sample")
    })

    look <- appearanceServer("look", plot_closure = shiny::reactive(rv$plot_fn))

    output$network_plot <- shiny::renderPlot({
      shiny::req(rv$net)
      s <- look()
      fn <- function() qgraph::qgraph(
        qgraph::getWmat(rv$net), layout = "spring",
        posCol = s$pos_edge, negCol = s$neg_edge,
        color = s$node_fill, border.color = s$node_border,
        vsize = s$vsize, esize = s$esize, label.cex = s$label_cex)
      rv$plot_fn <- fn

      rec_upsert(
        rec, "bootnet_plot", "plot",
        description = "[bootnet] Plotted the estimated network (qgraph, spring layout) with the chosen colours/sizes.",
        code = sprintf(
          paste('qgraph::qgraph(qgraph::getWmat(net), layout = "spring",',
                '  posCol = "%s", negCol = "%s",',
                '  color = "%s", border.color = "%s",',
                "  vsize = %s, esize = %s, label.cex = %s)", sep = "\n"),
          s$pos_edge, s$neg_edge, s$node_fill, s$node_border,
          s$vsize, s$esize, s$label_cex)
      )
      fn()
    })

    # TODO(stage3-bootnet): port from BootSON.R — split-group estimation
    # (:2763-2833), bootstrap edge-significance thresholding (:2903-2945,
    # keep p.adjust), Bridge-Symptoms tab (:6028-6089, keep its existing
    # corStability), NSON sub-app (:6357-6687), estimator branches for
    # huge/IsingFit/mgm (:3516-3741) — each new estimator branch MUST pin its
    # settings and register its own recorder fragment like the EBICglasso
    # path above. bootnet::differenceTest UI for edge/centrality pairs.

    list(net = shiny::reactive(rv$net),
         plot_fn = shiny::reactive(rv$plot_fn),
         settings = gg_settings)
  })
}
