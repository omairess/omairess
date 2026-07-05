# ============================================================================
# 11_tab_dagger.R — bnlearn DAG tab (from Dagger_zero.R)
#
# CORRECTNESS-CRITICAL LOGIC (do not alter casually):
#   * boot.strength (R >= 500) + averaged.network is the PRIMARY result;
#     the single algorithm run is demoted to "exploratory single fit"
#   * CPDAG / Markov-equivalence caveat computed and surfaced in the UI
#     whenever a DAG is displayed: arcs whose direction is not identified
#     are listed by name; the audit runs on the averaged network BEFORE any
#     cextend() (cextend manufactures directions, so a caveat computed after
#     it would understate what the data cannot identify)
#
# DECISION: blacklist/whitelist constraints are threaded into
#   algorithm.args, so they apply to EVERY bootstrap replicate (not just a
#   single fit) — the validated model respects the constraints. Whitelisted
#   arcs that would create a cycle make bnlearn error; caught and surfaced.
# DECISION: the "advanced" arc metrics are bnlearn's own boot.strength
#   outputs: strength = P(arc present), direction = P(this orientation |
#   present). Edge width / labels can show strength, direction, or their
#   product; the arc table always lists both so nothing is hidden.
# ============================================================================

# Equivalence-class audit for an averaged network. Returns the CPDAG, the
# list of direction-unidentified arcs, and boot direction confidence.
dag_equivalence_info <- function(avg_net, boot_str) {
  cpd <- bnlearn::cpdag(avg_net)
  und <- bnlearn::undirected.arcs(cpd)
  if (nrow(und) > 0) {
    key <- apply(und, 1, function(r) paste(sort(r), collapse = "~"))
    und <- und[!duplicated(key), , drop = FALSE]
  }
  arcs <- bnlearn::arcs(avg_net)
  dir_conf <- if (nrow(arcs) > 0) {
    m <- merge(as.data.frame(arcs), boot_str,
               by = c("from", "to"), all.x = TRUE)
    m$combined <- m$strength * m$direction
    m[order(-m$strength), c("from", "to", "strength", "direction", "combined")]
  } else NULL
  list(cpdag = cpd, undirected = und, dir_conf = dir_conf,
       n_arcs = nrow(arcs), n_unidentified = nrow(und))
}

# Build a from/to constraint data.frame from a cross product of node sets,
# dropping self-arcs. Returns NULL if empty (bnlearn treats NULL as "none").
constraint_df <- function(from_vars, to_vars) {
  if (!length(from_vars) || !length(to_vars)) return(NULL)
  g <- expand.grid(from = from_vars, to = to_vars, stringsAsFactors = FALSE)
  g <- g[g$from != g$to, , drop = FALSE]
  if (!nrow(g)) return(NULL)
  unique(g)
}

# Deparse a constraint data.frame into a literal for the exported script.
constraint_code <- function(df, name) {
  if (is.null(df) || !nrow(df))
    return(sprintf("%s <- NULL  # no constraints", name))
  sprintf("%s <- data.frame(from = %s, to = %s, stringsAsFactors = FALSE)",
          name, vars_literal(df$from), vars_literal(df$to))
}

daggerTabUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    varselectUI(ns("vars"), "Node set (Bayesian network / DAG)"),
    shiny::selectInput(ns("algorithm"), "Structure-learning algorithm",
                       c("Hill-Climbing" = "hc", "Tabu Search" = "tabu",
                         "MMHC" = "mmhc", "PC Stable" = "pc.stable",
                         "Grow-Shrink" = "gs", "IAMB" = "iamb",
                         "RSMAX2" = "rsmax2")),
    shiny::conditionalPanel(
      sprintf("['hc','tabu','mmhc','rsmax2'].includes(input['%s'])", ns("algorithm")),
      shiny::selectInput(ns("score"), "Network score",
                         c("BIC (Gaussian)" = "bic-g", "BGe" = "bge",
                           "BIC (discrete)" = "bic", "BDe" = "bde",
                           "AIC (Gaussian)" = "aic-g", "AIC (discrete)" = "aic"))
    ),
    shiny::numericInput(ns("boot_r"), "Bootstrap replicates (boot.strength)",
                        value = 500, min = 500, step = 100),
    shiny::sliderInput(ns("threshold"), "Arc-strength inclusion threshold",
                       min = 0.05, max = 0.95, value = 0.50, step = 0.05),

    # -- Constraint builder (blacklist / whitelist) --------------------------
    shiny::hr(),
    shiny::h5("Constraints (optional)"),
    shiny::helpText("Blacklist = arcs forbidden. Whitelist = arcs forced.",
                    "Pick FROM and TO node sets, then add the cross-product."),
    shiny::fluidRow(
      shiny::column(6, shinyWidgets::pickerInput(ns("con_from"), "FROM",
        choices = NULL, multiple = TRUE,
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, liveSearch = TRUE))),
      shiny::column(6, shinyWidgets::pickerInput(ns("con_to"), "TO",
        choices = NULL, multiple = TRUE,
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, liveSearch = TRUE)))
    ),
    shiny::actionButton(ns("add_bl"), "Add to blacklist", class = "btn-danger btn-sm"),
    shiny::actionButton(ns("add_wl"), "Add to whitelist", class = "btn-success btn-sm"),
    shiny::actionButton(ns("clear_con"), "Clear all", class = "btn-default btn-sm"),
    shiny::fluidRow(
      shiny::column(6, shiny::strong("Blacklist"), DT::dataTableOutput(ns("bl_table"))),
      shiny::column(6, shiny::strong("Whitelist"), DT::dataTableOutput(ns("wl_table")))
    ),

    shiny::hr(),
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
      shiny::column(4,
        shiny::selectInput(ns("layout_type"), "Layout",
                           c("Spring" = "spring", "Circle" = "circle",
                             "Hierarchical (Sugiyama)" = "tree")),
        shiny::selectInput(ns("edge_metric"), "Arc width / label reflects",
                           c("Strength  P(arc present)"          = "strength",
                             "Direction  P(orientation | present)" = "direction",
                             "Strength x Direction (combined)"     = "combined")),
        shiny::checkboxInput(ns("show_edge_labels"), "Show arc values", FALSE),
        # signed = FALSE: bnlearn arcs have no sign; the "positive edge" picker
        # is relabelled "Arc colour" and the negative-edge picker is hidden.
        appearanceUI(ns("look"), signed = FALSE)
      )
    )
    # TODO(port): split analysis + DAG-diff plot (Dagger_zero.R:409-542,
    # 4194-4272), folded temporal graph (:1369-1710), full 26-score list.
  )
}

daggerTabServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "dag")
    rv  <- shiny::reactiveValues(single = NULL, boot_str = NULL, avg = NULL,
                                 eq = NULL, plot_fn = NULL,
                                 bl = NULL, wl = NULL)

    # keep constraint pickers in sync with the selected node set
    shiny::observeEvent(sel$vars(), {
      shinyWidgets::updatePickerInput(session, "con_from", choices = sel$vars())
      shinyWidgets::updatePickerInput(session, "con_to",   choices = sel$vars())
    })

    # -- constraint accumulation --------------------------------------------
    arc_key <- function(df) if (is.null(df)) character(0)
                            else paste(df$from, df$to, sep = "->")
    shiny::observeEvent(input$add_bl, {
      new <- constraint_df(input$con_from, input$con_to)
      if (is.null(new)) return()
      rv$bl <- unique(rbind(rv$bl, new))
      # a blacklisted arc cannot also be whitelisted
      if (!is.null(rv$wl))
        rv$wl <- rv$wl[!arc_key(rv$wl) %in% arc_key(rv$bl), , drop = FALSE]
    })
    shiny::observeEvent(input$add_wl, {
      new <- constraint_df(input$con_from, input$con_to)
      if (is.null(new)) return()
      rv$wl <- unique(rbind(rv$wl, new))
      if (!is.null(rv$bl))
        rv$bl <- rv$bl[!arc_key(rv$bl) %in% arc_key(rv$wl), , drop = FALSE]
    })
    shiny::observeEvent(input$clear_con, { rv$bl <- NULL; rv$wl <- NULL })

    output$bl_table <- DT::renderDataTable(
      DT::datatable(rv$bl %||% data.frame(from = character(), to = character()),
                    options = list(dom = "tp", pageLength = 5), rownames = FALSE))
    output$wl_table <- DT::renderDataTable(
      DT::datatable(rv$wl %||% data.frame(from = character(), to = character()),
                    options = list(dom = "tp", pageLength = 5), rownames = FALSE))

    # -- estimation + validation --------------------------------------------
    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      dat  <- stats::na.omit(dat)
      algo <- input$algorithm
      uses_score <- algo %in% c("hc", "tabu", "mmhc", "rsmax2")
      scr  <- if (uses_score) input$score else NA_character_
      R    <- max(500L, as.integer(input$boot_r))   # skill floor: R >= 500
      thr  <- input$threshold
      seed <- rec_seed(rec)
      bl   <- rv$bl; wl <- rv$wl

      # algorithm.args: score (if score-based) + constraints, threaded into
      # EVERY bootstrap replicate so the validated model respects them.
      aargs <- list()
      if (uses_score) aargs$score <- scr
      if (!is.null(bl)) aargs$blacklist <- bl
      if (!is.null(wl)) aargs$whitelist <- wl

      res <- tryCatch(
        shiny::withProgress(message = "Structure learning + bootstrap validation", {
          shiny::incProgress(0.1, detail = "single algorithm run")
          algo_fn <- get(algo, envir = asNamespace("bnlearn"))
          single  <- do.call(algo_fn, c(list(x = dat), aargs))

          shiny::incProgress(0.4, detail = sprintf("boot.strength (R = %d)", R))
          set.seed(seed)
          boot_str <- bnlearn::boot.strength(dat, R = R, algorithm = algo,
                                             algorithm.args = aargs)
          avg <- bnlearn::averaged.network(boot_str, threshold = thr)
          eq  <- dag_equivalence_info(avg, boot_str)
          list(single = single, boot_str = boot_str, avg = avg, eq = eq)
        }),
        error = function(e) {
          shiny::showNotification(
            paste("Structure learning failed:", conditionMessage(e),
                  "(a whitelist that forces a cycle is a common cause)"),
            type = "error", duration = 12)
          NULL
        })
      shiny::req(res)

      con_desc <- sprintf("%d blacklisted, %d whitelisted arcs",
                          if (is.null(bl)) 0L else nrow(bl),
                          if (is.null(wl)) 0L else nrow(wl))

      rec_upsert(
        rec, "dag_analysis", "analysis",
        description = sprintf(
          "[DAG] Learned structure on %d variables (n = %d): algorithm = %s%s; %s; validated by boot.strength (R = %d, seed %d) + averaged.network (threshold = %.2f). Reported model is the bootstrap average, not a single run.",
          length(vars), nrow(dat), algo,
          if (uses_score) sprintf(", score = %s", scr) else "",
          con_desc, R, seed, thr),
        code = paste(
          sprintf("dag_vars <- %s", vars_literal(vars)),
          "dat_dag  <- na.omit(dat_wide[, dag_vars])",
          constraint_code(bl, "dag_blacklist"),
          constraint_code(wl, "dag_whitelist"),
          sprintf("dag_aargs <- list(%s)",
            paste(c(if (uses_score) sprintf('score = "%s"', scr),
                    "blacklist = dag_blacklist", "whitelist = dag_whitelist"),
                  collapse = ", ")),
          sprintf("set.seed(%d)", seed),
          sprintf('boot_str <- bnlearn::boot.strength(dat_dag, R = %d, algorithm = "%s",',
                  R, algo),
          "                                   algorithm.args = dag_aargs)",
          sprintf("avg_net  <- bnlearn::averaged.network(boot_str, threshold = %.2f)", thr),
          "# Single exploratory fit (NOT the reported model):",
          sprintf('# single <- do.call(bnlearn::%s, c(list(x = dat_dag), dag_aargs))', algo),
          sep = "\n")
      )

      rec_upsert(
        rec, "dag_stability", "stability",
        description = sprintf(
          "[DAG] Equivalence-class audit: %d of %d arcs have direction NOT identified by the data (undirected in the CPDAG); shown directions for those arcs are algorithmic convention, not evidence. Cross-sectional DAGs are identified only up to a Markov equivalence class.",
          res$eq$n_unidentified, res$eq$n_arcs),
        code = paste(
          "cpd <- bnlearn::cpdag(avg_net)",
          "undirected_in_cpdag <- bnlearn::undirected.arcs(cpd)",
          "# Arcs listed here are reversible within the Markov equivalence",
          "# class: the data cannot distinguish their direction. Do not make",
          "# causal-direction claims about them.", sep = "\n")
      )

      rv$single <- res$single; rv$boot_str <- res$boot_str
      rv$avg <- res$avg; rv$eq <- res$eq
    })

    # --- The CPDAG caveat, surfaced with every displayed DAG -----------------
    output$cpdag_caveat <- shiny::renderUI({
      shiny::req(rv$eq)
      eq <- rv$eq
      und_txt <- if (eq$n_unidentified > 0)
        paste(apply(eq$undirected, 1, paste, collapse = " — "), collapse = "; ")
      else "none"
      shiny::div(
        class = if (eq$n_unidentified > 0) "alert alert-warning" else "alert alert-info",
        shiny::strong("Markov-equivalence caveat: "),
        sprintf(
          "this cross-sectional DAG is identified only up to its equivalence class (CPDAG). %d of %d arcs have statistically unidentified direction: %s. Arrows on those arcs are an algorithmic convention, not causal evidence. Arc-direction confidence from the bootstrap is listed in the table (direction near 0.5 = undecidable).",
          eq$n_unidentified, eq$n_arcs, und_txt)
      )
    })

    look <- appearanceServer("look", plot_closure = shiny::reactive(rv$plot_fn))

    output$dag_plot <- shiny::renderPlot({
      shiny::req(rv$avg)
      s      <- look()
      amat   <- bnlearn::amat(rv$avg)
      metric <- input$edge_metric %||% "strength"
      dc     <- rv$eq$dir_conf

      # metric value per retained arc (strength / direction / combined)
      metric_mat <- matrix(0, nrow(amat), ncol(amat), dimnames = dimnames(amat))
      if (!is.null(dc) && nrow(dc) > 0) {
        for (i in seq_len(nrow(dc)))
          metric_mat[dc$from[i], dc$to[i]] <- dc[[metric]][i]
      }
      # width scales with the chosen metric; floor so thin arcs stay visible
      width_mat <- pmax(0.5, metric_mat * s$esize)

      # dashed = direction unidentified in the CPDAG (reinforces the caveat)
      lty_mat <- matrix(1, nrow(amat), ncol(amat), dimnames = dimnames(amat))
      und <- rv$eq$undirected
      if (!is.null(und) && nrow(und) > 0) {
        for (i in seq_len(nrow(und))) {
          a <- und[i, "from"]; b <- und[i, "to"]
          if (amat[a, b] == 1) lty_mat[a, b] <- 2
          if (amat[b, a] == 1) lty_mat[b, a] <- 2
        }
      }

      elabels <- if (isTRUE(input$show_edge_labels)) round(metric_mat, 2) else FALSE

      # node size: fixed slider or scaled by column means (shared option)
      if (isTRUE(s$scale_nodes)) {
        means <- vapply(data_bus$wide()[, colnames(amat), drop = FALSE],
                        function(x) if (is.numeric(x)) mean(x, na.rm = TRUE)
                                    else NA_real_, numeric(1))
        vsize_arg <- scale_vsize_by_mean(means, s$vsize_min, s$vsize_max)
      } else vsize_arg <- s$vsize

      # layout
      layout_arg <- switch(input$layout_type %||% "spring",
        circle = "circle",
        tree   = tryCatch({
          g  <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
          igraph::layout_with_sugiyama(g)$layout
        }, error = function(e) "spring"),
        "spring")

      fn <- function() qgraph::qgraph(
        amat, directed = TRUE, layout = layout_arg,
        edge.color = s$pos_edge,      # arcs are unsigned: pos picker only
        color = s$node_fill, border.color = s$node_border,
        vsize = vsize_arg, esize = s$esize, label.cex = s$label_cex,
        edge.width = width_mat, lty = lty_mat, edge.labels = elabels,
        edge.label.cex = 0.9)
      rv$plot_fn <- fn

      rec_upsert(
        rec, "dag_plot", "plot",
        description = sprintf(
          "[DAG] Plotted the averaged network (qgraph, %s layout); arc width/labels show %s; dashed arcs have unidentified direction per the CPDAG.",
          input$layout_type %||% "spring", metric),
        code = paste(
          "amat <- bnlearn::amat(avg_net)",
          "# arc-metric matrix from boot.strength (see arc table): strength /",
          "# direction / combined; dashed (lty=2) arcs are undirected in cpdag(avg_net).",
          sprintf('# metric shown: %s', metric),
          "arc_info <- merge(as.data.frame(bnlearn::arcs(avg_net)), boot_str,",
          '                  by = c("from","to"), all.x = TRUE)',
          "arc_info$combined <- arc_info$strength * arc_info$direction",
          sprintf('qgraph::qgraph(amat, directed = TRUE, layout = "%s",',
                  if (identical(input$layout_type, "tree")) "spring" else (input$layout_type %||% "spring")),
          sprintf('  edge.color = "%s", color = "%s", border.color = "%s",',
                  s$pos_edge, s$node_fill, s$node_border),
          sprintf("  vsize = %s, esize = %s, label.cex = %s)  # + edge.width by %s",
                  if (isTRUE(s$scale_nodes)) "vsize_arg" else as.character(s$vsize),
                  s$esize, s$label_cex, metric),
          sep = "\n")
      )
      fn()
    })

    output$arc_table <- DT::renderDataTable({
      shiny::req(rv$eq$dir_conf)
      DT::datatable(rv$eq$dir_conf,
                    caption = "Arc strength, direction confidence, and their product (bootstrap). strength = P(arc present); direction = P(this orientation | present); combined = strength x direction. direction ~ 0.5 means the data cannot decide the arrow.",
                    options = list(pageLength = 15)) |>
        DT::formatRound(c("strength", "direction", "combined"), 3)
    })

    list(avg = shiny::reactive(rv$avg), plot_fn = shiny::reactive(rv$plot_fn))
  })
}
