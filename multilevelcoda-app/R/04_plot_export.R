# ============================================================================
# 04_plot_export.R — colour pickers, size controls and PNG/PDF/SVG export for
# ggplot figures (house rules 4, 5, 9). Pastel defaults; every figure is
# re-drawn into a fresh device at the requested size for export.
# ============================================================================

partColoursUI <- function(ns, id = "colours") shiny::uiOutput(ns(id))

# Colour pickers, one per compositional part, with pastel defaults.
render_part_colours <- function(output, session, parts_r, id = "colours") {
  output[[id]] <- shiny::renderUI({
    parts <- parts_r()
    shiny::req(length(parts))
    pal <- pastel_palette(length(parts))
    shiny::tagList(
      shiny::tags$label("Part colours"),
      shiny::fluidRow(lapply(seq_along(parts), function(i)
        shiny::column(4, colourpicker::colourInput(
          session$ns(paste0(id, "_", i)), parts[i], pal[i])))))
  })
}

read_part_colours <- function(input, parts, id = "colours") {
  pal <- pastel_palette(length(parts))
  cols <- vapply(seq_along(parts), function(i)
    input[[paste0(id, "_", i)]] %||% pal[i], character(1))
  stats::setNames(cols, parts)
}

plotExportUI <- function(ns, prefix, w = 2400, h = 1600) {
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns(paste0(prefix, "_w")), "Width (px)", w, min = 300)),
      shiny::column(3, shiny::numericInput(ns(paste0(prefix, "_h")), "Height (px)", h, min = 300)),
      shiny::column(2, shiny::numericInput(ns(paste0(prefix, "_res")), "dpi", 300, min = 72)),
      shiny::column(4, shiny::radioButtons(ns(paste0(prefix, "_fmt")), "Format",
                                           c("png", "pdf", "svg"), inline = TRUE))),
    shiny::downloadButton(ns(paste0(prefix, "_dl")), "Download figure"),
    shiny::helpText("Copying a plot to the clipboard is not reliable from Shiny;",
                    "download it instead (or right-click the on-screen image).")
  )
}

plotExportServer <- function(input, output, prefix, plot_r, basename) {
  output[[paste0(prefix, "_dl")]] <- shiny::downloadHandler(
    filename = function() paste0(basename, ".", input[[paste0(prefix, "_fmt")]]),
    content = function(file) {
      w <- input[[paste0(prefix, "_w")]]; h <- input[[paste0(prefix, "_h")]]
      res <- input[[paste0(prefix, "_res")]]
      switch(input[[paste0(prefix, "_fmt")]],
        png = grDevices::png(file, width = w, height = h, res = res),
        pdf = grDevices::pdf(file, width = w / res, height = h / res),
        svg = svglite::svglite(file, width = w / res, height = h / res))
      print(plot_r())
      grDevices::dev.off()
    })
}

colours_literal <- function(cols) {
  sprintf("c(%s)", paste(sprintf('%s = "%s"', names(cols), cols), collapse = ", "))
}
