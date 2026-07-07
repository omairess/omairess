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

# --- Cascade layout (ported from Dagger_zero.R:1073-1114) -------------------
# Layered top-down placement: y = Katz centrality band (sources high, sinks
# low), x spread evenly within each band. Deterministic and overlap-free.
dagger_cascade_coords <- function(amat) {
  n <- nrow(amat)
  g <- igraph::graph_from_adjacency_matrix((amat > 0) * 1, mode = "directed",
                                           diag = FALSE)
  katz <- tryCatch(igraph::alpha_centrality(g, alpha = 0.1),
                   error = function(e) igraph::degree(g, mode = "in"))
  yn <- if (length(unique(katz)) > 1)
          1 - (katz - min(katz)) / (max(katz) - min(katz))
        else rep(0.5, n)
  y_lev  <- round(yn * 8) / 8            # snap to 8 bands
  coords <- matrix(0, n, 2, dimnames = list(rownames(amat), NULL))
  for (lv in unique(y_lev)) {
    idx <- which(abs(y_lev - lv) < 0.07)
    coords[idx, 1] <- if (length(idx) > 1)
                        seq(-1, 1, length.out = length(idx)) else 0
    coords[idx, 2] <- lv                 # same band = exactly same y
  }
  coords
}

# --- Cascade edge curving (ported from Dagger_zero.R:924-968) ----------------
# Bends any edge that would pass through an intermediate node, so the layered
# view stays readable. Returns a per-edge curvature matrix for qgraph.
compute_cascade_curves <- function(adj, coords) {
  n <- nrow(adj)
  curve_mat <- matrix(0, n, n)
  coord_range <- max(diff(range(coords[, 1])), diff(range(coords[, 2])), 0.5)
  node_radius <- coord_range * 0.10
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (adj[i, j] <= 0) next
    x1 <- coords[i, 1]; y1 <- coords[i, 2]
    x2 <- coords[j, 1]; y2 <- coords[j, 2]
    dx <- x2 - x1;      dy <- y2 - y1
    len_sq <- dx^2 + dy^2
    if (len_sq < 1e-10) next
    y_lo <- min(y1, y2); y_hi <- max(y1, y2)
    for (k in seq_len(n)) {
      if (k == i || k == j) next
      xk <- coords[k, 1]; yk <- coords[k, 2]
      if (yk < y_lo - node_radius || yk > y_hi + node_radius) next
      t <- max(0, min(1, ((xk - x1) * dx + (yk - y1) * dy) / len_sq))
      dist <- sqrt((xk - x1 - t * dx)^2 + (yk - y1 - t * dy)^2)
      if (dist < node_radius) {
        cross <- dx * (yk - y1) - dy * (xk - x1)
        curve_mat[i, j] <- if (cross > 0) 0.4 else -0.4
        break
      }
    }
  }
  curve_mat
}

# --- Score/data-type reconciliation ------------------------------------------
# bnlearn is strict: discrete scores (bic/aic/bde) need FACTOR columns,
# Gaussian scores (bic-g/aic-g/bge) need numeric. Binary 0/1 columns read
# from a file arrive numeric, so "BIC (discrete)" used to error while
# "BIC (Gaussian)" silently treated binary data as continuous. This coerces
# the data to what the chosen score needs (auto-detect for constraint-based
# algorithms: <= 5 unique values -> factor).
DAG_DISCRETE_SCORES <- c("bic", "aic", "bde")
DAG_GAUSSIAN_SCORES <- c("bic-g", "aic-g", "bge")

dag_prepare_types <- function(dat, score = NA_character_) {
  to_factor  <- function(d) { d[] <- lapply(d, function(x) as.factor(x)); d }
  to_numeric <- function(d) {
    d[] <- lapply(d, function(x)
      if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x))))
    d
  }
  if (!is.na(score) && score %in% DAG_DISCRETE_SCORES)
    return(list(dat = to_factor(dat), coerced = "factor",
                code = "dat_dag[] <- lapply(dat_dag, as.factor)  # discrete score needs factors"))
  if (!is.na(score) && score %in% DAG_GAUSSIAN_SCORES)
    return(list(dat = to_numeric(dat), coerced = "numeric",
                code = "dat_dag[] <- lapply(dat_dag, function(x) as.numeric(as.character(x)))  # Gaussian score needs numeric"))
  # constraint-based: auto-detect (few unique values -> categorical tests)
  uq <- vapply(dat, function(x) length(unique(stats::na.omit(x))), 0L)
  if (all(uq <= 5))
    list(dat = to_factor(dat), coerced = "factor (auto: <= 5 unique values)",
         code = "dat_dag[] <- lapply(dat_dag, as.factor)  # auto: few unique values -> categorical tests")
  else
    list(dat = to_numeric(dat), coerced = "numeric (auto)",
         code = "dat_dag[] <- lapply(dat_dag, function(x) as.numeric(as.character(x)))")
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
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'hc'", ns("algorithm")),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("restarts"),
          "Random restarts", value = 0, min = 0, step = 1)),
        shiny::column(6, shiny::numericInput(ns("perturb"),
          "Perturbations per restart", value = 1, min = 1, step = 1))
      ),
      shiny::helpText("Restarts re-run hill-climbing from perturbed graphs to",
                      "escape local optima (bnlearn's restart/perturb).")
    ),
    shiny::numericInput(ns("boot_r"), "Bootstrap replicates (boot.strength)",
                        value = 500, min = 500, step = 100),
    shiny::sliderInput(ns("threshold"), "Arc-strength inclusion threshold",
                       min = 0, max = 0.95, value = 0.50, step = 0.05),
    shiny::sliderInput(ns("dir_threshold"), "Direction inclusion threshold",
                       min = 0, max = 0.95, value = 0.50, step = 0.05),
    shiny::helpText("Displayed arcs need bootstrap strength >= the first",
                    "threshold AND direction confidence >= the second.",
                    "Set BOTH to 0 to see every arc that ever appeared",
                    "in the bootstrap (as in DAGger)."),

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
    shiny::helpText("Add as many FROM/TO batches as you like — each click",
                    "APPENDS to the list. Click rows in a table to select",
                    "them, then remove just those."),
    shiny::fluidRow(
      shiny::column(6, shiny::strong("Blacklist"),
        DT::dataTableOutput(ns("bl_table")),
        shiny::actionButton(ns("rm_bl"), "Remove selected", class = "btn-default btn-xs"),
        shiny::actionButton(ns("clear_bl"), "Clear blacklist", class = "btn-default btn-xs")),
      shiny::column(6, shiny::strong("Whitelist"),
        DT::dataTableOutput(ns("wl_table")),
        shiny::actionButton(ns("rm_wl"), "Remove selected", class = "btn-default btn-xs"),
        shiny::actionButton(ns("clear_wl"), "Clear whitelist", class = "btn-default btn-xs"))
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
        shiny::uiOutput(ns("dag_plot_ui")),        # height follows slider
        DT::dataTableOutput(ns("arc_table"))
      ),
      shiny::column(4,
        shiny::selectInput(ns("layout_type"), "Layout",
                           c("Spring" = "spring", "Circle" = "circle",
                             "Cascade (layered, curved edges)" = "cascade",
                             "Layered / dot-like (with controls)" = "dot",
                             "Stress majorization (neato-like)" = "kk",
                             "Force-directed (fdp-like)" = "fr",
                             "Manual (edit positions)" = "manual")),
        shiny::helpText("dot/neato/fdp are igraph equivalents (Sugiyama /",
                        "Kamada-Kawai / Fruchterman-Reingold) — Rgraphviz is",
                        "Bioconductor-only, which the app's CRAN installer",
                        "cannot fetch. Arrows and edge widths/labels are",
                        "independent of the layout and stay visible."),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'dot'", ns("layout_type")),
          shiny::sliderInput(ns("ranksep"), "Rank separation (vertical stretch)",
                             min = 0.2, max = 4, value = 1, step = 0.1),
          shiny::sliderInput(ns("nodesep"), "Node separation (horizontal stretch)",
                             min = 0.2, max = 4, value = 1, step = 0.1),
          shiny::sliderInput(ns("organic"),
                             "Relax rank alignment (0 = strict ranks, 1 = organic)",
                             min = 0, max = 1, value = 0.3, step = 0.05),
          shiny::checkboxInput(ns("layout_weights"),
                               "Use arc strength as layout weights", value = TRUE),
          shiny::helpText("Relaxing blends the layered coordinates with a",
                          "stress layout, so nodes stop sitting exactly",
                          "below one another while ranks stay recognisable.")
        ),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'manual'", ns("layout_type")),
          shiny::helpText("Double-click a cell to edit x/y (positions are",
                          "pinned across replots). Starts from the last",
                          "computed layout."),
          DT::dataTableOutput(ns("pos_table")),
          shiny::actionButton(ns("reset_pos"), "Reset to last computed layout",
                              class = "btn-default btn-xs")
        ),
        shiny::selectInput(ns("edge_metric"), "Arc width / label reflects",
                           c("Strength  P(arc present)"          = "strength",
                             "Direction  P(orientation | present)" = "direction",
                             "Strength x Direction (combined)"     = "combined")),
        # signed = FALSE: bnlearn arcs have no sign; the "positive edge" picker
        # is relabelled "Arc colour" and the negative-edge picker is hidden.
        # default_esize = 2, slider capped at 4 (arc widths multiply by
        # bootstrap strength, so large esize values are useless here).
        appearanceUI(ns("look"), signed = FALSE, default_esize = 2,
                     esize_max = 4)
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
    shiny::observeEvent(input$clear_bl, { rv$bl <- NULL })
    shiny::observeEvent(input$clear_wl, { rv$wl <- NULL })
    shiny::observeEvent(input$rm_bl, {
      sel <- input$bl_table_rows_selected
      if (length(sel) && !is.null(rv$bl)) {
        rv$bl <- rv$bl[-sel, , drop = FALSE]
        if (!nrow(rv$bl)) rv$bl <- NULL
      }
    })
    shiny::observeEvent(input$rm_wl, {
      sel <- input$wl_table_rows_selected
      if (length(sel) && !is.null(rv$wl)) {
        rv$wl <- rv$wl[-sel, , drop = FALSE]
        if (!nrow(rv$wl)) rv$wl <- NULL
      }
    })

    output$bl_table <- DT::renderDataTable(
      DT::datatable(rv$bl %||% data.frame(from = character(), to = character()),
                    options = list(dom = "tp", pageLength = 5), rownames = FALSE,
                    selection = "multiple"))
    output$wl_table <- DT::renderDataTable(
      DT::datatable(rv$wl %||% data.frame(from = character(), to = character()),
                    options = list(dom = "tp", pageLength = 5), rownames = FALSE,
                    selection = "multiple"))

    # -- estimation + validation --------------------------------------------
    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      # Per-module transform: after variable selection, before learning.
      trans <- sel$transform()
      dat  <- apply_house_transform(dat, trans)
      dat  <- stats::na.omit(dat)
      algo <- input$algorithm
      uses_score <- algo %in% c("hc", "tabu", "mmhc", "rsmax2")
      scr  <- if (uses_score) input$score else NA_character_
      if (!is.null(sel$group_var())) shiny::showNotification(
        "Grouping on the DAG tab is not yet supported (split analysis pending); the grouping variable is ignored here.",
        type = "warning", duration = 8)
      # Coerce column types to what the chosen score/algorithm needs
      # (fixes: binary 0/1 data + discrete BIC previously errored).
      prep <- dag_prepare_types(dat, scr)
      dat  <- prep$dat
      shiny::showNotification(sprintf("DAG variables treated as %s.", prep$coerced),
                              type = "message", duration = 5)
      R    <- max(500L, as.integer(input$boot_r))   # skill floor: R >= 500
      thr  <- input$threshold
      seed <- rec_seed(rec)
      bl   <- rv$bl; wl <- rv$wl

      # algorithm.args: score (if score-based) + constraints + hc restarts,
      # threaded into EVERY bootstrap replicate so the validated model
      # respects them.
      aargs <- list()
      if (uses_score) aargs$score <- scr
      if (algo == "hc") {
        aargs$restart <- max(0L, as.integer(input$restarts %||% 0))
        aargs$perturb <- max(1L, as.integer(input$perturb %||% 1))
      }
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
          "[DAG] Learned structure on %d variables (n = %d): algorithm = %s%s%s; %s; transform: %s; validated by boot.strength (R = %d, seed %d) + averaged.network (threshold = %.2f). Reported model is the bootstrap average, not a single run.",
          length(vars), nrow(dat), algo,
          if (uses_score) sprintf(", score = %s", scr) else "",
          if (algo == "hc") sprintf(", restarts = %d (perturb = %d)",
                                    aargs$restart, aargs$perturb) else "",
          con_desc, names(TRANSFORM_LABELS)[TRANSFORM_LABELS == trans],
          R, seed, thr),
        code = paste(
          sprintf("dag_vars <- %s", vars_literal(vars)),
          "dat_dag  <- dat_wide[, dag_vars]",
          transform_code_fragment(trans, "dat_dag"),
          "dat_dag  <- na.omit(dat_dag)",
          prep$code,
          constraint_code(bl, "dag_blacklist"),
          constraint_code(wl, "dag_whitelist"),
          sprintf("dag_aargs <- list(%s)",
            paste(c(if (uses_score) sprintf('score = "%s"', scr),
                    if (algo == "hc") sprintf("restart = %d, perturb = %d",
                                              aargs$restart, aargs$perturb),
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

    # Resizable plot window: height follows the appearance slider.
    output$dag_plot_ui <- shiny::renderUI({
      shiny::plotOutput(session$ns("dag_plot"),
                        height = sprintf("%dpx", look()$plot_height))
    })

    # --- Manual node positions (pin mode) ------------------------------------
    layout_store <- shiny::reactiveVal(NULL)   # last computed coords (named)
    manual_pos   <- shiny::reactiveVal(NULL)   # user-pinned coords

    output$pos_table <- DT::renderDataTable({
      co <- manual_pos() %||% layout_store()
      shiny::validate(shiny::need(!is.null(co),
        "Plot once with any other layout, then switch to Manual."))
      DT::datatable(data.frame(node = rownames(co),
                               x = round(co[, 1], 2), y = round(co[, 2], 2)),
                    editable = list(target = "cell", disable = list(columns = 0)),
                    options = list(pageLength = 25, dom = "tp"),
                    rownames = FALSE, selection = "none")
    })
    shiny::observeEvent(input$pos_table_cell_edit, {
      ed <- input$pos_table_cell_edit
      co <- manual_pos() %||% layout_store()
      shiny::req(co)
      val <- suppressWarnings(as.numeric(ed$value))
      # DT columns are 0-based with rownames = FALSE: 0 = node (locked),
      # 1 = x, 2 = y — which conveniently match the matrix columns 1/2.
      if (!is.na(val) && ed$col %in% c(1, 2)) {
        co[ed$row, ed$col] <- val
        manual_pos(co)
      }
    })
    shiny::observeEvent(input$reset_pos, manual_pos(NULL))

    # Rescale columns to [-1, 1] so different engines can be blended.
    .norm_coords <- function(co) {
      apply(co, 2, function(v) {
        r <- range(v)
        if (diff(r) < 1e-9) rep(0, length(v))
        else -1 + 2 * (v - r[1]) / diff(r)
      })
    }

    # DISPLAY arc set (request #4): every bootstrap arc passing BOTH the
    # strength and the direction thresholds, applied live at plot time.
    # Both at 0 = every arc that ever appeared (DAGger behaviour). The
    # REPORTED model (recorder / CPDAG audit) stays averaged.network at the
    # strength threshold chosen when "Learn & validate" was clicked.
    display_arcs <- shiny::reactive({
      shiny::req(rv$boot_str)
      bs <- as.data.frame(rv$boot_str)
      bs <- bs[bs$strength > 0 &
               bs$strength  >= input$threshold &
               bs$direction >= input$dir_threshold, , drop = FALSE]
      bs$combined <- bs$strength * bs$direction
      bs[order(-bs$strength), ]
    })

    output$dag_plot <- shiny::renderPlot({
      shiny::req(rv$avg)
      s      <- look()
      da     <- display_arcs()
      vars   <- bnlearn::nodes(rv$avg)
      metric <- input$edge_metric %||% "strength"

      # adjacency + metric matrices from the display arc set
      amat <- matrix(0, length(vars), length(vars),
                     dimnames = list(vars, vars))
      metric_mat <- amat
      if (nrow(da) > 0) for (i in seq_len(nrow(da))) {
        amat[da$from[i], da$to[i]]       <- 1
        metric_mat[da$from[i], da$to[i]] <- da[[metric]][i]
      }
      # width scales with the chosen metric; floor so thin arcs stay visible
      width_mat <- pmax(0.5, metric_mat * s$esize)

      # dashed = direction unidentified in the CPDAG (reinforces the caveat)
      lty_mat <- matrix(1, nrow(amat), ncol(amat), dimnames = dimnames(amat))
      und <- rv$eq$undirected
      if (!is.null(und) && nrow(und) > 0) {
        for (i in seq_len(nrow(und))) {
          a <- und[i, "from"]; b <- und[i, "to"]
          if (a %in% vars && b %in% vars) {
            if (amat[a, b] == 1) lty_mat[a, b] <- 2
            if (amat[b, a] == 1) lty_mat[b, a] <- 2
          }
        }
      }

      # edge labels: show the chosen arc metric when the appearance toggle is on
      elabels <- if (isTRUE(s$show_edge_labels)) round(metric_mat, 2) else FALSE

      # node size: fixed slider or scaled by column means (shared option)
      if (isTRUE(s$scale_nodes)) {
        means <- vapply(data_bus$wide()[, vars, drop = FALSE],
                        function(x) if (is.numeric(x)) mean(x, na.rm = TRUE)
                                    else NA_real_, numeric(1))
        vsize_arg <- scale_vsize_by_mean(means, s$vsize_min, s$vsize_max)
      } else vsize_arg <- s$vsize

      # layout — engines:
      #   cascade: DAGger's layered coords + per-edge curves
      #   dot:     Sugiyama with ranksep/nodesep controls and an "organic"
      #            blend toward a stress layout, breaking exact rank
      #            alignment while keeping the hierarchy readable
      #   kk/fr:   stress majorization (neato-like) / force-directed (fdp-like)
      #   manual:  user-pinned coordinates from the editable table
      lt <- input$layout_type %||% "spring"
      curve_arg <- 0.2; curve_all <- FALSE
      g_lay <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
      lay_w <- if (isTRUE(input$layout_weights) && nrow(da) > 0) {
        # heavier arcs pull nodes closer (graphviz `weight` analogue)
        e <- igraph::as_edgelist(g_lay)
        vapply(seq_len(nrow(e)), function(i) {
          hit <- which(da$from == e[i, 1] & da$to == e[i, 2])
          if (length(hit)) max(0.1, da$strength[hit[1]]) else 0.1
        }, numeric(1))
      } else NULL
      layout_arg <- switch(lt,
        circle  = "circle",
        cascade = {
          coords <- dagger_cascade_coords(amat)
          curve_mat <- compute_cascade_curves(amat, coords)
          curve_arg <- curve_mat; curve_all <- TRUE
          coords
        },
        dot = tryCatch({
          co <- igraph::layout_with_sugiyama(g_lay, weights = lay_w)$layout
          co <- .norm_coords(co)
          b  <- input$organic %||% 0
          if (b > 0) {                        # relax strict rank alignment
            # KK weights are desired DISTANCES -> invert so strong arcs attract
            co_kk <- .norm_coords(igraph::layout_with_kk(
              g_lay, weights = if (!is.null(lay_w)) 1 / lay_w))
            co <- (1 - b) * co + b * co_kk
          }
          # ranksep/nodesep: uniform scaling would be cancelled by qgraph's
          # own rescale, so they work as a STRETCH RATIO (y vs x), preserved
          # by plotting with aspect = TRUE (see qgraph call below).
          co[, 1] <- co[, 1] * (input$nodesep %||% 1)
          co[, 2] <- co[, 2] * (input$ranksep %||% 1)
          rownames(co) <- colnames(amat)
          co
        }, error = function(e) "spring"),
        kk = tryCatch({
          co <- .norm_coords(igraph::layout_with_kk(
            g_lay, weights = if (!is.null(lay_w)) 1 / lay_w))
          rownames(co) <- colnames(amat); co
        }, error = function(e) "spring"),
        fr = tryCatch({
          set.seed(rec_seed(rec))             # deterministic force layout
          co <- .norm_coords(igraph::layout_with_fr(g_lay, weights = lay_w))
          rownames(co) <- colnames(amat); co
        }, error = function(e) "spring"),
        manual = {
          co <- manual_pos() %||% layout_store()
          if (is.null(co)) "spring" else co[colnames(amat), , drop = FALSE]
        },
        "spring")

      # remember computed coordinates so Manual mode can start from them
      if (is.matrix(layout_arg) && lt != "manual") {
        co_store <- layout_arg
        rownames(co_store) <- colnames(amat)
        layout_store(co_store)
      }

      fn <- function() qgraph::qgraph(
        amat, directed = TRUE, layout = layout_arg,
        aspect = identical(lt, "dot"),  # preserve the ranksep/nodesep ratio
        labels = colnames(amat), label.scale = FALSE,   # equal-size node labels
        label.font = if (isTRUE(s$label_bold)) 2 else 1,
        edge.color = s$pos_edge,      # arcs are unsigned: pos picker only
        color = s$node_fill, border.color = s$node_border,
        vsize = vsize_arg, esize = s$esize, asize = s$asize,
        label.cex = s$label_cex,
        minimum = s$min_edge,
        edge.width = width_mat, lty = lty_mat, edge.labels = elabels,
        edge.label.cex = s$edge_label_cex,
        curve = curve_arg, curveAll = curve_all)
      rv$plot_fn <- fn

      rec_upsert(
        rec, "dag_plot", "plot",
        description = sprintf(
          "[DAG] Plotted arcs with bootstrap strength >= %.2f AND direction >= %.2f (%d arcs, %s layout); width/labels show %s; dashed arcs have unidentified direction per the CPDAG. The reported model remains the averaged network from the analysis step.",
          input$threshold, input$dir_threshold, nrow(da), lt, metric),
        code = paste(
          "# DISPLAY arc set: strength & direction thresholds applied to boot_str",
          "bs <- as.data.frame(boot_str)",
          sprintf("bs <- bs[bs$strength > 0 & bs$strength >= %.2f & bs$direction >= %.2f, ]",
                  input$threshold, input$dir_threshold),
          "bs$combined <- bs$strength * bs$direction",
          "vars <- unique(c(bs$from, bs$to, colnames(dat_dag)))",
          "amat <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))",
          "for (i in seq_len(nrow(bs))) amat[bs$from[i], bs$to[i]] <- 1",
          sprintf('# metric shown on arcs: %s', metric),
          if (lt == "manual" && is.matrix(layout_arg)) sprintf(
            "layout_coords <- matrix(c(%s), ncol = 2)  # pinned positions",
            paste(round(as.vector(layout_arg), 3), collapse = ", ")),
          sprintf('qgraph::qgraph(amat, directed = TRUE, layout = %s,',
                  if (lt == "manual" && is.matrix(layout_arg)) "layout_coords"
                  else if (lt %in% c("spring", "circle")) sprintf('"%s"', lt)
                  else '"spring"  # engine layout approximated in export'),
          "  labels = colnames(amat), label.scale = FALSE,",
          sprintf('  edge.color = "%s", color = "%s", border.color = "%s",',
                  s$pos_edge, s$node_fill, s$node_border),
          sprintf("  minimum = %s, edge.label.cex = %s,", s$min_edge, s$edge_label_cex),
          sprintf("  vsize = %s, esize = %s, label.cex = %s)  # + edge.width by %s",
                  if (isTRUE(s$scale_nodes)) "vsize_arg" else as.character(s$vsize),
                  s$esize, s$label_cex, metric),
          sep = "\n")
      )
      fn()
    })

    output$arc_table <- DT::renderDataTable({
      da <- display_arcs()
      DT::datatable(da[, c("from", "to", "strength", "direction", "combined")],
                    caption = "Arcs passing the current strength AND direction thresholds. strength = P(arc present); direction = P(this orientation | present); combined = their product. direction ~ 0.5 means the data cannot decide the arrow.",
                    options = list(pageLength = 15)) |>
        DT::formatRound(c("strength", "direction", "combined"), 3)
    })

    list(avg = shiny::reactive(rv$avg), plot_fn = shiny::reactive(rv$plot_fn))
  })
}
