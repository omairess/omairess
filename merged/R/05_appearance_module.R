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
appearanceUI <- function(id, signed = TRUE) {
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
    shiny::sliderInput(ns("vsize"), "Node size", min = 2, max = 20,
                       value = 8, step = 0.5),
    shiny::sliderInput(ns("esize"), "Edge scaling", min = 1, max = 20,
                       value = 8, step = 0.5),
    shiny::sliderInput(ns("label_cex"), "Label size", min = 0.3, max = 3,
                       value = 1, step = 0.1),
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
      list(
        node_fill   = input$node_fill   %||% pal$node_fill,
        node_border = input$node_border %||% pal$node_border,
        pos_edge    = input$pos_edge    %||% pal$pos_edge,
        neg_edge    = input$neg_edge    %||% pal$neg_edge,
        vsize       = input$vsize       %||% 8,
        esize       = input$esize       %||% 8,
        label_cex   = input$label_cex   %||% 1,
        export_w    = input$export_w    %||% 10,
        export_h    = input$export_h    %||% 8,
        export_res  = input$export_res  %||% 300
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
