# ==============================================================================
# server/22_missing_views.R — showing which points are measurements
#
# Three views of the same matrix (values$fill_status), because they answer
# different questions:
#   the map        which subjects and which times are affected, at a glance
#   the table      exactly how bad each subject is, sortable so the worst rise
#   the inspector  for one subject: the curve, its real points, and its filled
#                  ones, so you can see what the fit is doing across a gap
#
# Colour never carries this alone: the map has a legend and per-cell hover text
# naming the state, and the table repeats the counts as numbers.
# ==============================================================================

# --- the map -----------------------------------------------------------------
output$missing_map <- renderPlotly({
  req(values$fill_status)
  st   <- values$fill_status
  labs <- values$time_labels
  if (is.null(labs) || length(labs) != ncol(st)) labs <- paste0("T", seq_len(ncol(st)))
  rows <- fck_row_labels(values)

  ord <- switch(input$missing_sort %||% "file",
    file     = seq_len(nrow(st)),
    observed = order(rowSums(st == FCK_FILL_OBSERVED)),
    filled   = order(-rowSums(st != FCK_FILL_OBSERVED)))
  st <- st[ord, , drop = FALSE]; rows <- rows[ord]

  hover <- matrix(sprintf("%s<br>%s<br><b>%s</b>",
                          rep(rows, times = ncol(st)),
                          rep(labs, each  = nrow(st)),
                          FCK_FILL_LABELS[as.vector(st) + 1L]),
                  nrow = nrow(st))

  # a discrete three-step scale: plotly wants a continuous one, so each state
  # gets a flat band rather than a gradient between them
  scale <- list(
    list(0.000, unname(FCK_FILL_COLORS["Observed"])),
    list(0.333, unname(FCK_FILL_COLORS["Observed"])),
    list(0.333, unname(FCK_FILL_COLORS["Interpolated"])),
    list(0.667, unname(FCK_FILL_COLORS["Interpolated"])),
    list(0.667, unname(FCK_FILL_COLORS["Extrapolated"])),
    list(1.000, unname(FCK_FILL_COLORS["Extrapolated"])))

  plot_ly(z = st, x = labs, y = rows, type = "heatmap",
          colorscale = scale, zmin = 0, zmax = 2, showscale = FALSE,
          xgap = 1, ygap = 1,                       # a surface gap between cells
          hoverinfo = "text", text = hover) %>%
    layout(
      title  = list(text = "Every cell: measured, or filled in by the smoother?",
                    font = list(size = 14)),
      xaxis  = list(title = "", tickangle = -90, tickfont = list(size = 9)),
      yaxis  = list(title = "", tickfont = list(size = 8), autorange = "reversed"),
      margin = list(l = 110, b = 70, t = 40))
})

output$missing_legend <- renderUI({
  swatch <- function(state, note) {
    tags$span(style = "display:inline-block; margin-right:22px; white-space:nowrap;",
      tags$span(style = sprintf(
        "display:inline-block; width:13px; height:13px; background:%s; border:1px solid #999; vertical-align:middle; margin-right:6px;",
        FCK_FILL_COLORS[[state]])),
      tags$b(state), tags$span(style = "color:#52514e;", paste0(" — ", note)))
  }
  tagList(
    swatch("Observed",     "a real measurement"),
    swatch("Interpolated", "filled between two observed points of that subject"),
    swatch("Extrapolated", "filled beyond that subject's first/last observation")
  )
})

output$missing_headline <- renderText({
  req(values$fill_status)
  st <- values$fill_status
  hps <- fck_hours_per_step(values)
  extrap_rows <- sum(rowSums(st == FCK_FILL_EXTRAPOLATED) > 0)
  half_rows   <- sum(rowMeans(st != FCK_FILL_OBSERVED) > 0.5)

  paste0(
    fck_fill_headline(st), "\n",
    sprintf("%d of %d rows have values carried beyond their observed range",
            extrap_rows, nrow(st)),
    if (is.finite(hps)) sprintf(" (up to %.0f h past the last measurement)",
                                max(rowSums(st == FCK_FILL_EXTRAPOLATED)) * hps) else "",
    ".\n",
    if (half_rows > 0)
      sprintf("%d row%s more than half filled in — check %s in the table below before trusting %s.\n",
              half_rows, if (half_rows == 1) " is" else "s are",
              if (half_rows == 1) "it" else "them",
              if (half_rows == 1) "it" else "them")
    else "No row is more than half filled in.\n")
})

# --- the table ---------------------------------------------------------------
output$missing_table <- renderDT({
  req(values$fill_status)
  st  <- values$fill_status
  hps <- fck_hours_per_step(values)
  tab <- fck_fill_per_subject(st, hps)
  tab$Row <- fck_row_labels(values)
  names(tab)[1] <- "Subject"

  labs <- values$time_labels
  if (!is.null(labs) && length(labs) == ncol(st)) {
    tab$`First observed` <- labs[tab$`First observed`]
    tab$`Last observed`  <- labs[tab$`Last observed`]
  }

  datatable(tab, rownames = FALSE,
            options = list(pageLength = 10, scrollX = TRUE,
                           order = list(list(3, "desc"))),   # worst extrapolation first
            caption = "Sorted by extrapolated points. 'Longest gap' is the widest run of consecutive filled points inside the observed range.") %>%
    formatStyle("Extrapolated",
                background = styleInterval(0, c("", FCK_FILL_COLORS[["Extrapolated"]])),
                color      = styleInterval(0, c("", "white")))
})

# --- the inspector -----------------------------------------------------------
output$missing_subject_ui <- renderUI({
  req(values$data)
  rows <- fck_row_labels(values)
  st <- values$fill_status
  # default to the worst offender: that is the one worth looking at first
  sel <- if (!is.null(st)) which.max(rowSums(st != FCK_FILL_OBSERVED)) else 1L
  selectInput("missing_subject", "Inspect one subject:",
              choices = stats::setNames(seq_along(rows), rows), selected = sel)
})

output$missing_curve <- renderPlotly({
  req(values$data, values$fill_status, input$missing_subject)
  i  <- as.integer(input$missing_subject)
  st <- values$fill_status[i, ]
  raw <- values$data[i, ]
  fit <- if (!is.null(values$smooth_data)) values$smooth_data[i, ] else rep(NA_real_, length(raw))

  x <- fck_cumulative_hours(values$time_labels)
  xlab <- "Hours from the first time point"
  if (is.null(x) || length(x) != length(raw)) { x <- seq_along(raw); xlab <- "Time point" }
  labs <- values$time_labels
  if (is.null(labs) || length(labs) != length(raw)) labs <- paste0("T", seq_along(raw))

  p <- plot_ly()
  if (any(is.finite(fit))) {
    p <- add_trace(p, x = x, y = fit, type = "scatter", mode = "lines",
                   name = "Smoothed curve",
                   line = list(color = "#52514e", width = 2),
                   hoverinfo = "text",
                   text = sprintf("%s<br>fitted %.2f", labs, fit))
  }
  for (state in c("Observed", "Interpolated", "Extrapolated")) {
    k <- which(st == (match(state, FCK_FILL_LABELS) - 1L))
    if (!length(k)) next
    y <- if (state == "Observed") raw[k] else fit[k]
    p <- add_trace(p, x = x[k], y = y, type = "scatter", mode = "markers",
                   name = state,
                   marker = list(
                     size   = if (state == "Observed") 10 else 9,
                     color  = if (state == "Observed") FCK_FILL_COLORS[["Observed"]] else "rgba(0,0,0,0)",
                     symbol = if (state == "Observed") "circle" else
                              if (state == "Interpolated") "circle-open" else "x",
                     line   = list(width = 2, color = unname(FCK_FILL_COLORS[[state]]))),
                   hoverinfo = "text",
                   text = sprintf("%s<br>%s<br>%.2f", labs[k], state, y))
  }
  n_obs <- sum(st == FCK_FILL_OBSERVED)
  plotly::layout(p,
    title = list(text = sprintf("%s — %d of %d points measured",
                                fck_row_labels(values)[i], n_obs, length(st)),
                 font = list(size = 14)),
    xaxis = list(title = xlab),
    yaxis = list(title = "Value"),
    legend = list(orientation = "h", y = -0.18),
    margin = list(b = 70))
})

# --- export ------------------------------------------------------------------
output$export_fill_status_csv <- downloadHandler(
  filename = function() paste0("fill_status_", Sys.Date(), ".csv"),
  content = function(file) {
    req(values$fill_status)
    st <- values$fill_status
    df <- as.data.frame(matrix(FCK_FILL_LABELS[as.vector(st) + 1L], nrow = nrow(st)))
    labs <- values$time_labels
    names(df) <- if (!is.null(labs) && length(labs) == ncol(st)) labs else paste0("T", seq_len(ncol(st)))
    df <- cbind(Subject = fck_row_labels(values), df)
    write.csv(df, file, row.names = FALSE)
  }
)
