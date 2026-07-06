# ============================================================================
# 05_appearance_module.R — shared colours / sizes / figure export (rules 4,5,9)
#
# One appearance contract for all tabs. Returns a reactive list:
#   node_fill, node_border, pos_edge, neg_edge   (rule 4 colour pickers)
#   vsize, esize, label_cex                      (rule 5 sliders)
#   export_w, export_h, export_res               (rule 5 export sizing)
#
# DECISION: Dagger's bnlearn arcs are UNSIGNED, so its tab maps pos_edge onto
#   arc colour and ignores neg_edge (with a UI note) — we do not fake a sign
#   concept bnlearn doesn't have.
# DECISION: figure export re-plots the stored zero-argument plot closure into
#   a fresh device (never grabs the screen device), giving BootSON plot
#   export it never had and adding SVG everywhere (svglite).
# ============================================================================

# signed = FALSE for tabs whose edges have no sign (bnlearn arcs): the
# positive-edge picker is relabelled "Arc colour" and the negative-edge
# picker is hidden rather than faking a sign concept that doesn't exist.
# default_esize / default_edge_labels let a tab set its own starting values
# (DAG wants thinner edges; the RI network wants edge labels on).
appearanceUI <- function(id, signed = TRUE, default_esize = 8,
                         default_edge_labels = FALSE) {
  ns  <- shiny::NS(id)
  pal <- house_pastel()
  shiny::tagList(
    shiny::h5("Colours"),
    colourpicker::colourInput(ns("node_fill"), "Node colour",
                              value = pal$node_fill, palette = "square"),
    colourpicker::colourInput(ns("node_border"), "Node border colour",
                              value = pal$node_border, palette = "square"),
    colourpicker::colourInput(ns("pos_edge"),
                              if (signed) "Positive edges" else "Arc colour",
                              value = pal$pos_edge, palette = "square"),
    if (signed) colourpicker::colourInput(ns("neg_edge"), "Negative edges",
                                          value = pal$neg_edge, palette = "square"),
    shiny::h5("Sizes"),
    shiny::checkboxInput(ns("scale_nodes"),
                         "Scale node size by column mean", value = FALSE),
    shiny::conditionalPanel(
      sprintf("!input['%s']", ns("scale_nodes")),
      shiny::sliderInput(ns("vsize"), "Node size", min = 2, max = 20,
                         value = 8, step = 0.5)
    ),
    shiny::conditionalPanel(
      sprintf("input['%s']", ns("scale_nodes")),
      # node size ranges linearly between these when scaled by column mean
      shiny::sliderInput(ns("vsize_range"), "Node size range (min-max)",
                         min = 2, max = 20, value = c(4, 14), step = 0.5),
      shiny::helpText("Smallest node = lowest column mean, largest = highest.")
    ),
    shiny::sliderInput(ns("esize"), "Edge scaling", min = 1, max = 20,
                       value = default_esize, step = 0.5),
    shiny::sliderInput(ns("label_cex"), "Node label size", min = 0.3, max = 3,
                       value = 1, step = 0.1),
    shiny::helpText("Node labels are drawn at one fixed size (not scaled per",
                    "node), so every label is equally legible."),
    shiny::h5("Edges shown"),
    # minimum absolute weight for an edge to be DRAWN (qgraph `minimum`)
    shiny::numericInput(ns("min_edge"),
                        "Minimum |value| to show an edge", value = 0,
                        min = 0, step = 0.01),
    shiny::checkboxInput(ns("show_edge_labels"), "Show edge values as labels",
                         value = default_edge_labels),
    shiny::conditionalPanel(
      sprintf("input['%s']", ns("show_edge_labels")),
      shiny::sliderInput(ns("edge_label_cex"), "Edge label size",
                         min = 0.3, max = 2, value = 0.8, step = 0.1)
    ),
    shiny::sliderInput(ns("plot_height"), "Plot window height (px)",
                       min = 300, max = 1400, value = 520, step = 20),
    shiny::h5("Export"),
    shiny::fluidRow(
      shiny::column(4, shiny::numericInput(ns("export_w"), "Width (in)", 10, min = 2)),
      shiny::column(4, shiny::numericInput(ns("export_h"), "Height (in)", 8, min = 2)),
      shiny::column(4, shiny::numericInput(ns("export_res"), "PNG dpi", 300, min = 72))
    ),
    shiny::downloadButton(ns("dl_png"), "PNG"),
    shiny::downloadButton(ns("dl_pdf"), "PDF"),
    shiny::downloadButton(ns("dl_svg"), "SVG")
  )
}

appearanceServer <- function(id, plot_closure) {
  # plot_closure: reactive returning a zero-arg function that redraws the
  # tab's current figure (each tab stores one after plotting).
  shiny::moduleServer(id, function(input, output, session) {
    pal <- house_pastel()

    settings <- shiny::reactive({
      vr <- input$vsize_range %||% c(4, 14)
      list(
        node_fill        = input$node_fill   %||% pal$node_fill,
        node_border      = input$node_border %||% pal$node_border,
        pos_edge         = input$pos_edge    %||% pal$pos_edge,
        neg_edge         = input$neg_edge    %||% pal$neg_edge,
        scale_nodes      = isTRUE(input$scale_nodes),
        vsize            = input$vsize       %||% 8,
        vsize_min        = vr[1],
        vsize_max        = vr[2],
        esize            = input$esize       %||% 8,
        label_cex        = input$label_cex   %||% 1,
        min_edge         = input$min_edge    %||% 0,
        show_edge_labels = isTRUE(input$show_edge_labels),
        edge_label_cex   = input$edge_label_cex %||% 0.8,
        plot_height      = input$plot_height %||% 520,
        export_w         = input$export_w    %||% 10,
        export_h         = input$export_h    %||% 8,
        export_res       = input$export_res  %||% 300
      )
    })

    # Shared export factory — fresh device per download (rule 5).
    make_plot_download <- function(fmt) {
      shiny::downloadHandler(
        filename = function() sprintf("network_%s.%s",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), fmt),
        content = function(file) {
          s  <- settings()
          fn <- plot_closure()
          shiny::req(is.function(fn))
          switch(fmt,
            png = grDevices::png(file, width = s$export_w, height = s$export_h,
                                 units = "in", res = s$export_res),
            pdf = grDevices::pdf(file, width = s$export_w, height = s$export_h),
            svg = svglite::svglite(file, width = s$export_w, height = s$export_h)
          )
          on.exit(grDevices::dev.off(), add = TRUE)
          fn()
        }
      )
    }
    output$dl_png <- make_plot_download("png")
    output$dl_pdf <- make_plot_download("pdf")
    output$dl_svg <- make_plot_download("svg")

    settings
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Shared qgraph arguments built from an appearance `settings` list, so every
# network in the app obeys the same rules:
#   * label.scale = FALSE       -> node labels are ALL the same size (rule #2)
#   * labels = colnames(W)      -> real names, never numbered-with-legend (#4)
#   * edge.labels + edge.label.cex from the toggle/slider           (#3)
#   * minimum = min |value| to draw an edge                          (#9)
# `directed` flips arrows on; `node_col`/`groups` allow community colouring.
house_qgraph_args <- function(W, s, directed = FALSE,
                              node_col = NULL, groups = NULL,
                              vsize = NULL) {
  labs <- colnames(W); if (is.null(labs)) labs <- rownames(W)
  list(
    directed        = directed,
    labels          = labs,             # explicit names (no numbered legend)
    label.scale     = FALSE,            # all node labels equal size
    label.cex       = s$label_cex,
    vsize           = vsize %||% s$vsize,
    esize           = s$esize,
    minimum         = s$min_edge,       # hide edges below this |value|
    edge.labels     = isTRUE(s$show_edge_labels),
    edge.label.cex  = s$edge_label_cex,
    posCol          = s$pos_edge,
    negCol          = s$neg_edge,
    color           = node_col %||% s$node_fill,
    border.color    = s$node_border,
    groups          = groups
  )
}

# The exported-script twin: renders the same arguments as R source text.
# `wobj` is the name of the weight-matrix variable in the script (e.g. "W").
house_qgraph_args_code <- function(s, directed, wobj = "W",
                                   node_col_expr = NULL, groups_expr = NULL,
                                   vsize_expr = NULL) {
  lines <- c(
    sprintf("  directed = %s,", directed),
    sprintf("  labels = colnames(%s), label.scale = FALSE, label.cex = %s,",
            wobj, s$label_cex),
    sprintf("  vsize = %s, esize = %s, minimum = %s,",
            vsize_expr %||% as.character(s$vsize), s$esize, s$min_edge),
    sprintf("  edge.labels = %s, edge.label.cex = %s,",
            isTRUE(s$show_edge_labels), s$edge_label_cex),
    sprintf('  posCol = "%s", negCol = "%s",', s$pos_edge, s$neg_edge),
    sprintf("  color = %s, border.color = \"%s\"",
            node_col_expr %||% sprintf('"%s"', s$node_fill), s$node_border),
    if (!is.null(groups_expr)) sprintf("  , groups = %s", groups_expr))
  paste(lines, collapse = "\n")
}
