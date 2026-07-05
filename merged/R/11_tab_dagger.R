# ============================================================================
# 11_tab_dagger.R — bnlearn DAG tab (from Dagger_zero.R)
#
# CORRECTNESS-CRITICAL LOGIC WRITTEN IN FULL (do not alter in stage 3):
#   * boot.strength (R >= 500) + averaged.network is the PRIMARY result;
#     the single algorithm run is demoted to "exploratory single fit"
#   * CPDAG / Markov-equivalence caveat computed and surfaced in the UI
#     whenever a DAG is displayed: arcs whose direction is not identified
#     are listed by name, and cextend()'ed graphs are labelled as one
#     arbitrary member of the equivalence class
#
# DECISION: the equivalence-class analysis runs on the averaged network
#   BEFORE any cextend() call — cextend manufactures directions, so a caveat
#   computed after it would understate what the data cannot identify.
# DECISION: score and algorithm are always passed explicitly (Dagger already
#   complied); the recorder fragment bakes algorithm, score, R, threshold,
#   and seed so the exported script is settings-complete.
# DECISION: unqualified causal UI copy from Dagger_zero ("gold standard for
#   causal discovery" etc.) is NOT carried over; stage 3 must port help text
#   through the cpdag caveat wording below.
# ============================================================================

# Equivalence-class audit for an averaged network. Returns the CPDAG, the
# list of direction-unidentified arcs, and boot direction confidence.
dag_equivalence_info <- function(avg_net, boot_str) {
  cpd <- bnlearn::cpdag(avg_net)
  und <- bnlearn::undirected.arcs(cpd)
  # undirected arcs appear in both directions in bnlearn; deduplicate
  if (nrow(und) > 0) {
    key <- apply(und, 1, function(r) paste(sort(r), collapse = "~"))
    und <- und[!duplicated(key), , drop = FALSE]
  }
  # direction confidence from the bootstrap (0.5 = coin flip)
  arcs <- bnlearn::arcs(avg_net)
  dir_conf <- if (nrow(arcs) > 0) {
    m <- merge(as.data.frame(arcs), boot_str,
               by = c("from", "to"), all.x = TRUE)
    m[order(m$direction), c("from", "to", "strength", "direction")]
  } else NULL
  list(cpdag = cpd, undirected = und, dir_conf = dir_conf,
       n_arcs = nrow(arcs), n_unidentified = nrow(und))
}

daggerTabUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    varselectUI(ns("vars"), "Node set (Bayesian network / DAG)"),
    shiny::selectInput(ns("algorithm"), "Structure-learning algorithm (pinned)",
                       c("Hill-Climbing" = "hc", "Tabu Search" = "tabu",
                         "MMHC" = "mmhc", "PC Stable" = "pc.stable",
                         "Grow-Shrink" = "gs", "IAMB" = "iamb",
                         "RSMAX2" = "rsmax2")),
    shiny::selectInput(ns("score"), "Network score (pinned)",
                       c("BIC (Gaussian)" = "bic-g", "BGe" = "bge",
                         "BIC (discrete)" = "bic", "BDe" = "bde")),
    # TODO(stage3-dagger-ui): restore Dagger's full 26-score list + per-score
    # help text (Dagger_zero.R:51-125), gated by detected data type.
    shiny::numericInput(ns("boot_r"), "Bootstrap replicates (boot.strength)",
                        value = 500, min = 500, step = 100),
    shiny::sliderInput(ns("threshold"), "Arc-strength inclusion threshold",
                       min = 0.05, max = 0.95, value = 0.50, step = 0.05),
    shiny::actionButton(ns("run"), "Learn & validate structure",
                        class = "btn-primary"),
    shiny::helpText("The reported network is the bootstrap-AVERAGED structure;",
                    "a single algorithm run is shown only as an exploratory fit."),
    shiny::hr(),
    shiny::uiOutput(ns("cpdag_caveat")),
    shiny::fluidRow(
      shiny::column(8,
        shiny::plotOutput(ns("dag_plot")),
        DT::dataTableOutput(ns("arc_table"))
      ),
      # signed = FALSE: bnlearn arcs have no sign; the "positive edge" picker
      # is relabelled "Arc colour" and the negative-edge picker is hidden.
      shiny::column(4, appearanceUI(ns("look"), signed = FALSE))
    )
    # TODO(stage3-dagger-ui): blacklist/whitelist constraint builder
    # (Dagger_zero.R:2312-2379, 4140-4189), split analysis, DAG-diff plot,
    # folded temporal graph tabs — port as-is (they already consume the
    # underscore wide naming this app's data contract guarantees).
  )
}

daggerTabServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "dag")
    rv  <- shiny::reactiveValues(single = NULL, boot_str = NULL, avg = NULL,
                                 eq = NULL, plot_fn = NULL)

    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      # TODO(stage3-dagger): port prepare_mixed_data() type handling
      # (Dagger_zero.R:228-285) — continuous as numeric, discrete as factor,
      # zero-variance checks — before this point.
      dat  <- stats::na.omit(dat)
      algo <- input$algorithm
      scr  <- input$score
      R    <- max(500L, as.integer(input$boot_r))   # skill floor: R >= 500
      thr  <- input$threshold
      seed <- rec_seed(rec)

      shiny::withProgress(message = "Structure learning + bootstrap validation", {

        # -- 1. Exploratory single fit (NOT the reported model) --------------
        shiny::incProgress(0.1, detail = "single algorithm run")
        algo_fn <- get(algo, envir = asNamespace("bnlearn"))
        single  <- if (algo %in% c("hc", "tabu"))
                     algo_fn(dat, score = scr) else algo_fn(dat)

        # -- 2. Bootstrap model averaging: the reported model -----------------
        shiny::incProgress(0.4, detail = sprintf("boot.strength (R = %d)", R))
        set.seed(seed)
        boot_str <- bnlearn::boot.strength(
          dat, R = R, algorithm = algo,
          algorithm.args = if (algo %in% c("hc", "tabu"))
                             list(score = scr) else list())
        avg <- bnlearn::averaged.network(boot_str, threshold = thr)

        # -- 3. Markov-equivalence audit BEFORE any cextend -------------------
        eq <- dag_equivalence_info(avg, boot_str)

        rec_upsert(
          rec, "dag_analysis", "analysis",
          description = sprintf(
            "[DAG] Learned structure on %d variables (n = %d): algorithm = %s, score = %s; validated by boot.strength (R = %d, seed %d) + averaged.network (threshold = %.2f). Reported model is the bootstrap average, not a single run.",
            length(vars), nrow(dat), algo, scr, R, seed, thr),
          code = sprintf(
            paste("dag_vars <- %s",
                  "dat_dag  <- na.omit(dat_wide[, dag_vars])",
                  "set.seed(%d)",
                  'boot_str <- bnlearn::boot.strength(dat_dag, R = %d,',
                  '  algorithm = "%s", algorithm.args = %s)',
                  "avg_net  <- bnlearn::averaged.network(boot_str, threshold = %.2f)",
                  "# Single exploratory fit (NOT the reported model):",
                  '# single <- bnlearn::%s(dat_dag%s)', sep = "\n"),
            vars_literal(vars), seed, R, algo,
            if (algo %in% c("hc", "tabu"))
              sprintf('list(score = "%s")', scr) else "list()",
            thr, algo,
            if (algo %in% c("hc", "tabu")) sprintf(', score = "%s"', scr) else "")
        )

        rec_upsert(
          rec, "dag_stability", "stability",
          description = sprintf(
            "[DAG] Equivalence-class audit: %d of %d arcs have direction NOT identified by the data (undirected in the CPDAG); shown directions for those arcs are algorithmic convention, not evidence. Cross-sectional DAGs are identified only up to a Markov equivalence class.",
            eq$n_unidentified, eq$n_arcs),
          code = paste(
            "cpd <- bnlearn::cpdag(avg_net)",
            "undirected_in_cpdag <- bnlearn::undirected.arcs(cpd)",
            "# Arcs listed here are reversible within the Markov equivalence",
            "# class: the data cannot distinguish their direction. Do not make",
            "# causal-direction claims about them.", sep = "\n")
        )

        rv$single <- single; rv$boot_str <- boot_str
        rv$avg <- avg; rv$eq <- eq
      })
    })

    # --- The CPDAG caveat, surfaced with every displayed DAG -----------------
    output$cpdag_caveat <- shiny::renderUI({
      shiny::req(rv$eq)
      eq <- rv$eq
      und_txt <- if (eq$n_unidentified > 0)
        paste(apply(eq$undirected, 1, paste, collapse = " — "),
              collapse = "; ")
      else "none"
      shiny::div(
        class = if (eq$n_unidentified > 0) "alert alert-warning"
                else "alert alert-info",
        shiny::strong("Markov-equivalence caveat: "),
        sprintf(
          "this cross-sectional DAG is identified only up to its equivalence class (CPDAG). %d of %d arcs have statistically unidentified direction: %s. Arrows on those arcs are an algorithmic convention, not causal evidence. Arc-direction confidence from the bootstrap is listed in the table (direction near 0.5 = undecidable).",
          eq$n_unidentified, eq$n_arcs, und_txt)
      )
    })

    look <- appearanceServer("look", plot_closure = shiny::reactive(rv$plot_fn))

    output$dag_plot <- shiny::renderPlot({
      shiny::req(rv$avg)
      s    <- look()
      amat <- bnlearn::amat(rv$avg)

      # Edge width reflects bootstrap arc strength (falls back to a flat
      # width if a bootstrap strength isn't available for some reason).
      width_mat <- matrix(s$esize * 0.3, nrow(amat), ncol(amat),
                          dimnames = dimnames(amat))
      dc <- rv$eq$dir_conf
      if (!is.null(dc) && nrow(dc) > 0) {
        for (i in seq_len(nrow(dc)))
          width_mat[dc$from[i], dc$to[i]] <- s$esize * dc$strength[i]
      }

      # Dashed styling for arcs whose direction is NOT identified by the
      # data (undirected in the CPDAG) — visually surfaces the caveat above.
      lty_mat <- matrix(1, nrow(amat), ncol(amat), dimnames = dimnames(amat))
      und <- rv$eq$undirected
      if (!is.null(und) && nrow(und) > 0) {
        for (i in seq_len(nrow(und))) {
          a <- und[i, "from"]; b <- und[i, "to"]
          if (amat[a, b] == 1) lty_mat[a, b] <- 2
          if (amat[b, a] == 1) lty_mat[b, a] <- 2
        }
      }

      fn <- function() qgraph::qgraph(
        amat, directed = TRUE, layout = "spring",
        edge.color = s$pos_edge,      # arcs are unsigned: pos picker only
        color = s$node_fill, border.color = s$node_border,
        vsize = s$vsize, esize = s$esize, label.cex = s$label_cex,
        edge.width = width_mat, lty = lty_mat)
      rv$plot_fn <- fn

      rec_upsert(
        rec, "dag_plot", "plot",
        description = "[DAG] Plotted the averaged network (qgraph, spring layout); edge width reflects bootstrap arc strength, dashed arcs have unidentified direction per the CPDAG.",
        code = sprintf(
          paste("amat <- bnlearn::amat(avg_net)",
                "# width_mat: bootstrap arc strength scaled by esize; dashed",
                "# (lty = 2) arcs are undirected in cpdag(avg_net) (see dag_stability step).",
                'qgraph::qgraph(amat, directed = TRUE, layout = "spring",',
                '  edge.color = "%s", color = "%s", border.color = "%s",',
                "  vsize = %s, esize = %s, label.cex = %s,",
                "  edge.width = width_mat, lty = lty_mat)", sep = "\n"),
          s$pos_edge, s$node_fill, s$node_border, s$vsize, s$esize, s$label_cex)
      )
      fn()
    })

    output$arc_table <- DT::renderDataTable({
      shiny::req(rv$eq$dir_conf)
      DT::datatable(rv$eq$dir_conf,
                    caption = "Arc strength & direction confidence (bootstrap). direction ~ 0.5 means the data cannot decide the arrow.",
                    options = list(pageLength = 15)) |>
        DT::formatRound(c("strength", "direction"), 3)
    })

    # TODO(stage3-dagger): port split-group analysis (Dagger_zero.R:4194-4272)
    # + common layout (:288-406), DAG-diff plot (:409-542), folded temporal
    # graph (:1369-1710), score/algorithm help text (scrubbed of unqualified
    # causal language), cextend option — if exposed, its output MUST carry the
    # label "one arbitrary member of the equivalence class" and rv$eq must
    # keep being computed pre-cextend.

    list(avg = shiny::reactive(rv$avg), plot_fn = shiny::reactive(rv$plot_fn))
  })
}
