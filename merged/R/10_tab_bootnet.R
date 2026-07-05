# ============================================================================
# 10_tab_bootnet.R — bootnet network tab (from BootSON.R)
#
# CORRECTNESS-CRITICAL LOGIC (do not alter casually):
#   * every estimator's settings pinned explicitly, never package defaults;
#     corMethod defaults to cor_auto for the correlation-based estimators
#   * bootstraps run ON DEMAND (user decision 2026-07-05, overriding the
#     original run-with-estimation design) — BUT the correctness gate stands:
#     centrality and edge-CI outputs stay LOCKED until the nonparametric +
#     case-dropping bootstraps have been run for the CURRENT network, and
#     re-estimating drops any stale stability results/fragments
#   * CS-coefficient gate: CS < 0.25 BLOCKS the ordered centrality display;
#     values remain visible only alphabetically under a red banner
#   * expected influence force-included and listed first whenever the
#     estimated network contains any negative edge
#
# DECISION: all estimators (EBICglasso / ggmModSelect / pcor / cor / huge /
#   IsingFit / mgm / relimp) go through bootnet::estimateNetwork so the same
#   bootstrap machinery, CS gate, and recorder wiring covers every branch.
# DECISION: relimp (relative-importance) networks are DIRECTED and
#   non-negative; the plot flips to directed for that estimator.
# DECISION: EGA is used here for COMMUNITIES/COLOURING/LAYOUT only (a plot
#   step in the recorder), not as the estimator — the reported network is
#   always the pinned estimateNetwork result.
# ============================================================================

# Per-statistic interpretation gate (Epskamp, Borsboom & Fried, 2018).
cs_gate <- function(cs) {
  if (!is.finite(cs)) "unknown"
  else if (cs < 0.25) "blocked"
  else if (cs < 0.50) "caution"
  else "ok"
}

# Pinned arguments per estimator. Returns list(args, desc); args go straight
# into estimateNetwork AND (deparsed) into the exported fragment, so the app
# and the script cannot diverge.
bootnet_pinned <- function(def, input) {
  switch(def,
    EBICglasso = list(
      args = list(corMethod = input$bootnet_cormethod,
                  tuning    = input$bootnet_tuning),
      desc = sprintf("corMethod = %s, EBIC gamma = %s",
                     input$bootnet_cormethod, input$bootnet_tuning)),
    ggmModSelect = list(
      args = list(corMethod = input$bootnet_cormethod,
                  tuning    = input$gms_tuning,
                  stepwise  = isTRUE(input$gms_stepwise)),
      desc = sprintf("corMethod = %s, gamma = %s, stepwise = %s",
                     input$bootnet_cormethod, input$gms_tuning,
                     isTRUE(input$gms_stepwise))),
    pcor = list(
      args = list(corMethod = input$bootnet_cormethod,
                  threshold = input$pcor_threshold,
                  alpha     = input$pcor_alpha),
      desc = sprintf("corMethod = %s, threshold = %s, alpha = %s",
                     input$bootnet_cormethod, input$pcor_threshold,
                     input$pcor_alpha)),
    cor = list(
      args = list(corMethod = input$bootnet_cormethod,
                  threshold = input$pcor_threshold,
                  alpha     = input$pcor_alpha),
      desc = sprintf("corMethod = %s, threshold = %s, alpha = %s",
                     input$bootnet_cormethod, input$pcor_threshold,
                     input$pcor_alpha)),
    huge = list(
      args = list(tuning    = input$huge_tuning,
                  criterion = input$huge_criterion,
                  npn       = isTRUE(input$huge_npn)),
      desc = sprintf("EBIC gamma = %s, criterion = %s, nonparanormal = %s",
                     input$huge_tuning, input$huge_criterion,
                     isTRUE(input$huge_npn))),
    IsingFit = list(
      args = list(tuning = input$ising_tuning, rule = input$ising_rule),
      desc = sprintf("gamma = %s, rule = %s (binary data required)",
                     input$ising_tuning, input$ising_rule)),
    mgm = list(
      args = list(tuning    = input$mgm_tuning,
                  criterion = input$mgm_criterion,
                  rule      = input$mgm_rule),
      desc = sprintf("EBIC gamma = %s, criterion = %s, rule = %s (mixed types auto-detected)",
                     input$mgm_tuning, input$mgm_criterion, input$mgm_rule)),
    relimp = list(
      args = list(normalized = isTRUE(input$relimp_normalized)),
      desc = sprintf("relative-importance network (directed, non-negative), normalized = %s",
                     isTRUE(input$relimp_normalized)))
  )
}

# Deparse an args list into "  name = value" lines for the exported script.
args_code <- function(args) {
  paste(vapply(names(args), function(nm)
    sprintf("  %s = %s", nm, deparse(args[[nm]])), ""), collapse = ",\n")
}

# Node-size-by-column-mean helper (shared with the DAG tab; the exported
# fragment emits this exact function so app and script match).
scale_vsize_by_mean <- function(means, vmin, vmax) {
  m <- means
  if (all(is.na(m))) return(rep((vmin + vmax) / 2, length(m)))
  m[is.na(m)] <- mean(m, na.rm = TRUE)
  r <- range(m)
  n <- if (diff(r) > 0) (m - r[1]) / diff(r) else rep(0.5, length(m))
  vmin + n * (vmax - vmin)
}

SCALE_VSIZE_FRAGMENT <- paste(
  "scale_vsize <- function(m, vmin, vmax) {",
  "  if (all(is.na(m))) return(rep((vmin + vmax) / 2, length(m)))",
  "  m[is.na(m)] <- mean(m, na.rm = TRUE); r <- range(m)",
  "  n <- if (diff(r) > 0) (m - r[1]) / diff(r) else rep(0.5, length(m))",
  "  vmin + n * (vmax - vmin)",
  "}", sep = "\n")

bootnetTabUI <- function(id) {
  ns <- shiny::NS(id)
  cond <- function(vals) paste(sprintf("input['%s'] == '%s'",
                                       ns("bootnet_default"), vals),
                               collapse = " || ")
  shiny::tagList(
    varselectUI(ns("vars"), "Node set (bootnet network)"),
    shiny::selectInput(ns("bootnet_default"), "Estimator (pinned)",
                       c("EBICglasso (regularized GGM)"        = "EBICglasso",
                         "ggmModSelect (non-regularized GGM)"  = "ggmModSelect",
                         "pcor (partial correlations)"          = "pcor",
                         "cor (marginal correlations)"          = "cor",
                         "huge (glasso via huge)"               = "huge",
                         "IsingFit (binary data)"               = "IsingFit",
                         "mgm (mixed graphical model)"          = "mgm",
                         "relimp (relative importance network)" = "relimp")),

    # -- per-estimator pinned settings ---------------------------------------
    shiny::conditionalPanel(cond(c("EBICglasso", "ggmModSelect", "pcor", "cor")),
      shiny::selectInput(ns("bootnet_cormethod"), "Correlation method (pinned)",
                         c("cor_auto (polychoric auto-detect)" = "cor_auto",
                           "spearman" = "spearman",
                           "cor (Pearson — continuous data only)" = "cor"),
                         selected = "cor_auto")),
    shiny::conditionalPanel(cond("EBICglasso"),
      shiny::numericInput(ns("bootnet_tuning"), "EBIC gamma (pinned)",
                          value = 0.5, min = 0, max = 1, step = 0.05)),
    shiny::conditionalPanel(cond("ggmModSelect"),
      shiny::numericInput(ns("gms_tuning"), "gamma (pinned; 0 = BIC selection)",
                          value = 0, min = 0, max = 1, step = 0.05),
      shiny::checkboxInput(ns("gms_stepwise"), "Stepwise search", value = TRUE)),
    shiny::conditionalPanel(cond(c("pcor", "cor")),
      shiny::selectInput(ns("pcor_threshold"), "Thresholding",
                         c("none", "sig", "holm", "bonferroni", "BH"),
                         selected = "sig"),
      shiny::numericInput(ns("pcor_alpha"), "alpha",
                          value = 0.05, min = 0.001, max = 0.2, step = 0.005)),
    shiny::conditionalPanel(cond("huge"),
      shiny::numericInput(ns("huge_tuning"), "EBIC gamma (pinned)",
                          value = 0.5, min = 0, max = 1, step = 0.05),
      shiny::selectInput(ns("huge_criterion"), "Selection criterion",
                         c("ebic", "ric", "stars"), selected = "ebic"),
      shiny::checkboxInput(ns("huge_npn"),
                           "Nonparanormal transform (huge.npn)", value = TRUE)),
    shiny::conditionalPanel(cond("IsingFit"),
      shiny::numericInput(ns("ising_tuning"), "gamma (pinned)",
                          value = 0.25, min = 0, max = 1, step = 0.05),
      shiny::selectInput(ns("ising_rule"), "Rule", c("OR", "AND")),
      shiny::helpText("Requires binary (two-valued) variables. Dichotomize",
                      "on the psychonetrics tab's transform step if needed.")),
    shiny::conditionalPanel(cond("mgm"),
      shiny::numericInput(ns("mgm_tuning"), "EBIC gamma (pinned)",
                          value = 0.25, min = 0, max = 1, step = 0.05),
      shiny::selectInput(ns("mgm_criterion"), "Criterion", c("EBIC", "CV")),
      shiny::selectInput(ns("mgm_rule"), "Rule", c("AND", "OR"))),
    shiny::conditionalPanel(cond("relimp"),
      shiny::checkboxInput(ns("relimp_normalized"),
                           "Normalize relative importance", value = TRUE)),

    shiny::actionButton(ns("run"), "Estimate network", class = "btn-primary"),
    shiny::hr(),

    # -- on-demand validation (correctness gate still enforced) --------------
    shiny::numericInput(ns("nboots"), "Bootstrap samples (both types)",
                        value = 1000, min = 250, step = 250),
    shiny::actionButton(ns("run_boot"),
                        "Run bootstrap validation (required before interpreting centrality / edge order)",
                        class = "btn-warning"),
    shiny::helpText("On demand: runs the nonparametric bootstrap (edge CIs) +",
                    "case-dropping bootstrap (CS-coefficients) for the current",
                    "network. Centrality and edge-CI panels stay locked until",
                    "it has been run; re-estimating clears old results."),
    shiny::hr(),
    shiny::uiOutput(ns("stability_banner")),
    shiny::fluidRow(
      shiny::column(8,
        shiny::plotOutput(ns("network_plot")),
        shiny::plotOutput(ns("edge_ci_plot")),
        shiny::uiOutput(ns("centrality_ui"))
      ),
      shiny::column(4,
        shiny::selectInput(ns("layout_type"), "Layout",
                           c("Spring" = "spring", "Circle" = "circle",
                             "EGA communities" = "ega",
                             "PCA (first two components)" = "pca")),
        appearanceUI(ns("look"), signed = TRUE)
      )
    )
    # TODO(port): Bridge-Symptoms + NSON sub-tabs from BootSON.R:1633-2224
    # unchanged (Bridge already has its own case-drop + corStability — keep
    # it); split-group estimation; bootstrap edge-significance thresholding;
    # bootnet::differenceTest UI for edge/centrality pairs.
  )
}

bootnetTabServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "bootnet")

    # Pinned GGM settings: SINGLE SOURCE, exported to the NCT module via the
    # return value. NCT always compares EBICglasso-style GGMs, so it reads
    # these EBICglasso inputs regardless of the estimator chosen here.
    gg_settings <- shiny::reactive({
      list(
        default   = "EBICglasso",
        corMethod = input$bootnet_cormethod %||% "cor_auto",
        tuning    = input$bootnet_tuning    %||% 0.5
      )
    })

    rv <- shiny::reactiveValues(net = NULL, def = NULL, dat = NULL,
                                boot_np = NULL, boot_cd = NULL, cs = NULL,
                                has_neg = FALSE, plot_fn = NULL)

    # ---- 1. Estimation (instant, no bootstraps) -----------------------------
    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      def  <- input$bootnet_default

      if (def == "IsingFit") {
        n_uniq <- vapply(dat, function(x) length(unique(stats::na.omit(x))), 0L)
        if (any(n_uniq > 2)) {
          shiny::showNotification(
            sprintf("IsingFit needs binary variables; not binary: %s",
                    paste(vars[n_uniq > 2], collapse = ", ")),
            type = "error", duration = 10)
          return()
        }
      }

      pin <- bootnet_pinned(def, input)

      net <- tryCatch(
        do.call(bootnet::estimateNetwork,
                c(list(data = dat, default = def), pin$args)),
        error = function(e) {
          shiny::showNotification(paste("Estimation error:", conditionMessage(e)),
                                  type = "error", duration = 10)
          NULL
        })
      shiny::req(net)

      wmat <- qgraph::getWmat(net)
      rv$has_neg <- any(wmat[upper.tri(wmat)] < 0) || any(wmat > 0 & t(wmat) < 0)
      rv$net <- net; rv$def <- def; rv$dat <- dat

      # A new/changed network voids any previous stability results — the CS
      # gate must never be satisfied by bootstraps of a DIFFERENT network.
      rv$boot_np <- NULL; rv$boot_cd <- NULL; rv$cs <- NULL
      rec_drop(rec, "bootnet_stability")

      rec_upsert(
        rec, "bootnet_analysis", "analysis",
        description = sprintf(
          "[bootnet] Estimated %s network on %d variables (n = %d): %s.%s Bootstrap validation not yet run for this network.",
          def, length(vars), nrow(dat), pin$desc,
          if (rv$has_neg) " Contains negative edges -> expected influence reported first."
          else ""),
        code = sprintf(
          paste("bootnet_vars <- %s",
                "dat_bootnet  <- dat_wide[, bootnet_vars]",
                "net <- bootnet::estimateNetwork(dat_bootnet,",
                '  default = "%s",',
                "%s)", sep = "\n"),
          vars_literal(vars), def, args_code(pin$args))
      )
    })

    # ---- 2. On-demand bootstrap validation ----------------------------------
    shiny::observeEvent(input$run_boot, {
      shiny::req(rv$net)
      nb   <- max(250L, as.integer(input$nboots))
      if (nb < 500) shiny::showNotification(
        "nBoots < 500: stability estimates will be rough; 1000 recommended.",
        type = "warning")
      seed <- rec_seed(rec)

      shiny::withProgress(message = "Bootstrap validation", {
        shiny::incProgress(0.1, detail = "nonparametric bootstrap (edge CIs)")
        set.seed(seed)
        boot_np <- bootnet::bootnet(rv$net, nBoots = nb, type = "nonparametric")

        cd_stats <- c("strength",
                      if (rv$has_neg) "expectedInfluence",
                      "closeness", "betweenness")
        shiny::incProgress(0.5, detail = "case-dropping bootstrap (CS)")
        set.seed(seed)
        boot_cd <- bootnet::bootnet(rv$net, nBoots = nb, type = "case",
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

        rv$boot_np <- boot_np; rv$boot_cd <- boot_cd; rv$cs <- cs
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
      if (rv$has_neg) tb[, c("node", "expected_influence", "strength")] else tb
    })

    # --- The CS gate ---------------------------------------------------------
    lead_cs <- shiny::reactive({
      shiny::req(rv$cs)
      key <- if (rv$has_neg) "expectedInfluence" else "strength"
      if (key %in% names(rv$cs)) unname(rv$cs[[key]]) else NA_real_
    })

    output$stability_banner <- shiny::renderUI({
      shiny::req(rv$net)
      if (is.null(rv$cs)) {
        return(shiny::div(class = "alert alert-info",
          "Bootstrap validation has not been run for this network — centrality
           and edge-CI panels are locked until you run it (button above)."))
      }
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
        "Centrality locked: run the bootstrap validation first — centrality must not be interpreted without the case-dropping bootstrap (CS-coefficient)."))
      if (cs_gate(lead_cs()) == "blocked") {
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
      shiny::req(rv$net)
      shiny::validate(shiny::need(!is.null(rv$boot_np),
        "Edge-weight CIs locked: run the bootstrap validation to see them."))
      # bootnet's CI plot; overlapping CIs = edge order not interpretable.
      plot(rv$boot_np, labels = FALSE, order = "sample")
    })

    # ---- Plot: layouts (spring/circle/EGA/PCA) + node-size-by-mean ----------
    look <- appearanceServer("look", plot_closure = shiny::reactive(rv$plot_fn))

    output$network_plot <- shiny::renderPlot({
      shiny::req(rv$net)
      s    <- look()
      wmat <- qgraph::getWmat(rv$net)
      dat  <- rv$dat
      directed <- identical(rv$def, "relimp")
      lt   <- input$layout_type %||% "spring"

      # Layout / grouping
      groups_list <- NULL
      node_cols   <- s$node_fill
      layout_arg  <- "spring"
      layout_code <- 'layout_arg <- "spring"'
      if (lt == "circle") {
        layout_arg  <- "circle"
        layout_code <- 'layout_arg <- "circle"'
      } else if (lt == "pca") {
        pc <- tryCatch(stats::prcomp(scale(stats::na.omit(dat))),
                       error = function(e) NULL)
        if (!is.null(pc)) {
          layout_arg  <- pc$rotation[, 1:2]
          layout_code <- paste(
            "pca <- stats::prcomp(scale(na.omit(dat_bootnet)))",
            "layout_arg <- pca$rotation[, 1:2]  # nodes at their PC1/PC2 loadings",
            sep = "\n")
        } else shiny::showNotification("PCA layout failed; using spring.",
                                       type = "warning")
      } else if (lt == "ega") {
        ega <- tryCatch(
          EGAnet::EGA(dat, model = "glasso", algorithm = "walktrap",
                      plot.EGA = FALSE),
          error = function(e) NULL)
        if (!is.null(ega) && !is.null(ega$wc)) {
          groups_list <- split(colnames(wmat), ega$wc)
          node_cols   <- house_group_colors(length(groups_list))
          layout_code <- paste(
            "# EGA used for communities/colouring only; the reported network",
            "# is still the pinned estimateNetwork result above.",
            'ega <- EGAnet::EGA(dat_bootnet, model = "glasso",',
            '                   algorithm = "walktrap", plot.EGA = FALSE)',
            "groups_list <- split(colnames(dat_bootnet), ega$wc)",
            'layout_arg <- "spring"', sep = "\n")
        } else shiny::showNotification("EGA failed; using spring layout.",
                                       type = "warning")
      }

      # Node size: fixed slider or scaled by column means (user option)
      if (isTRUE(s$scale_nodes)) {
        means <- vapply(dat, function(x)
          if (is.numeric(x)) mean(x, na.rm = TRUE) else NA_real_, numeric(1))
        vsize_arg  <- scale_vsize_by_mean(means, s$vsize_min, s$vsize_max)
        vsize_code <- paste(
          SCALE_VSIZE_FRAGMENT,
          "node_means <- vapply(dat_bootnet, function(x)",
          "  if (is.numeric(x)) mean(x, na.rm = TRUE) else NA_real_, numeric(1))",
          sprintf("vsize_arg <- scale_vsize(node_means, %s, %s)",
                  s$vsize_min, s$vsize_max), sep = "\n")
      } else {
        vsize_arg  <- s$vsize
        vsize_code <- sprintf("vsize_arg <- %s", s$vsize)
      }

      fn <- function() qgraph::qgraph(
        wmat, directed = directed, layout = layout_arg,
        groups = groups_list,
        posCol = s$pos_edge, negCol = s$neg_edge,
        color = node_cols, border.color = s$node_border,
        vsize = vsize_arg, esize = s$esize, label.cex = s$label_cex)
      rv$plot_fn <- fn

      rec_upsert(
        rec, "bootnet_plot", "plot",
        description = sprintf(
          "[bootnet] Plotted the %s network (qgraph, %s layout%s%s).",
          rv$def, lt,
          if (lt == "ega") ", coloured by EGA walktrap communities" else "",
          if (isTRUE(s$scale_nodes)) ", node size scaled by column means" else ""),
        code = paste(
          layout_code, vsize_code,
          sprintf(
            paste("qgraph::qgraph(qgraph::getWmat(net), directed = %s,",
                  "  layout = layout_arg,%s",
                  '  posCol = "%s", negCol = "%s",',
                  '  color = %s, border.color = "%s",',
                  "  vsize = vsize_arg, esize = %s, label.cex = %s)",
                  sep = "\n"),
            directed,
            if (lt == "ega") "\n  groups = groups_list," else "",
            s$pos_edge, s$neg_edge,
            # bake literal colours — the script must run without app helpers
            if (lt == "ega") vars_literal(node_cols)
            else sprintf('"%s"', s$node_fill),
            s$node_border, s$esize, s$label_cex),
          sep = "\n")
      )
      fn()
    })

    list(net = shiny::reactive(rv$net),
         plot_fn = shiny::reactive(rv$plot_fn),
         settings = gg_settings)
  })
}
