# ==============================================================================
# server/76_dlmo.R -- dim-light melatonin onset (DLMO)
#
# The estimators are in server/08l_helpers_dlmo.R; this file reads the inputs,
# runs them and draws the results. A run is stored in values$dlmo_results, and
# the report (server/93_apa_report.R) and the script (server/90_export.R) read
# it from there.
# ==============================================================================

dlmo_covariates <- reactive(names(values$covariates) %||% character(0))

output$dlmo_var_ui <- renderUI({
  cv <- dlmo_covariates()
  tagList(
    selectInput("dlmo_night_var", "Within participants: night / condition column (one row per night or condition):",
                choices = c("None: one row per participant" = "_none_", cv),
                selected = isolate(input$dlmo_night_var) %||% "_none_"),
    selectInput("dlmo_group_var", "Between participants: group column:",
                choices = c("None" = "_none_", cv),
                selected = isolate(input$dlmo_group_var) %||% "_none_"))
})

output$dlmo_bed_var_ui <- renderUI({
  tagList(
    selectInput("dlmo_bed_var", "Bedtime column:", choices = dlmo_covariates(),
                selected = isolate(input$dlmo_bed_var)),
    helpText("HH:MM, decimal hours (23.5), or an Excel time (a day fraction below 1)."))
})

output$dlmo_hockey_status <- renderUI({
  if (requireNamespace("dlmoR", quietly = TRUE))
    dance_status(sprintf("Hockey-stick: dlmoR %s found.", as.character(utils::packageVersion("dlmoR"))))
  else dance_status(HTML(paste0(
    "Hockey-stick: the <b>dlmoR</b> package is not installed, so this method is unavailable ",
    "(it is not reimplemented here). dlmoR is on GitHub only: ",
    "<code>remotes::install_github(\"tscnlab/dlmoR\")</code>, then restart the app.")),
    style = "color:#8a5a12")
})

# Hours from "18", "18.5" or "18:30"; a time earlier than the one before is
# after midnight (dance_unwrap_clock).
dlmo_times <- function() {
  n <- ncol(values$data)
  if (identical(input$dlmo_time_src, "_manual_")) {
    s <- trimws(strsplit(input$dlmo_manual_times %||% "", ",", fixed = TRUE)[[1]])
    s <- s[nzchar(s)]
    h <- vapply(s, function(x) {
      if (grepl(":", x, fixed = TRUE)) {
        p <- suppressWarnings(as.numeric(strsplit(x, ":", fixed = TRUE)[[1]]))
        if (length(p) < 2L || any(!is.finite(p))) NA_real_ else sum(p * c(1, 1 / 60, 1 / 3600)[seq_along(p)])
      } else suppressWarnings(as.numeric(x))
    }, numeric(1), USE.NAMES = FALSE)
    if (length(h) != n || any(!is.finite(h)))
      stop(sprintf("%d readable time(s) entered for %d sample columns.", sum(is.finite(h)), n))
  } else {
    h <- values$time_clock
    if (is.null(h) || length(h) != n)
      stop("No clock times could be read from the column names: choose 'Enter the times'.")
  }
  t <- dance_unwrap_clock(as.numeric(h))$times
  # two columns at the same time cannot both be samples of one profile
  if (any(diff(t) <= 0)) {
    k <- which(diff(t) <= 0)[1]
    stop(sprintf(paste("Columns %d and %d read as the same time (%s): the hour in each column name is",
                       "taken from its first number. Enter the times instead."),
                 k, k + 1L, dance_clock_label(t[k], 24, show_day = FALSE)))
  }
  t
}

dlmo_bedtime <- function() {
  n <- nrow(values$data)
  switch(input$dlmo_bed_src %||% "none",
    none = NULL,
    value = {
      b <- dance_dlmo_clock_hours(input$dlmo_bed_value %||% "")
      if (!is.finite(b)) stop("The bedtime could not be read: enter it as HH:MM.")
      rep(b, n)
    },
    column = {
      v <- values$covariates[[input$dlmo_bed_var %||% ""]]
      if (is.null(v)) stop("Choose the bedtime column.")
      b <- dance_dlmo_clock_hours(v)
      if (!any(is.finite(b))) stop("No bedtime in that column could be read as a time of day.")
      b
    })
}

observeEvent(input$run_dlmo, {
  req(values$data)
  res <- tryCatch({
    times <- dlmo_times()
    meth <- intersect(c("fixed", "variable", "hockey"), input$dlmo_methods %||% character(0))
    if (!length(meth)) stop("Choose at least one method.")
    if ("hockey" %in% meth && !requireNamespace("dlmoR", quietly = TRUE))
      stop("The hockey-stick method needs the dlmoR package: remotes::install_github(\"tscnlab/dlmoR\").")
    thr <- input$dlmo_threshold
    if (!isTRUE(is.finite(thr)) || thr <= 0) stop("Enter a threshold above 0.")
    nv <- input$dlmo_night_var %||% "_none_"; gv <- input$dlmo_group_var %||% "_none_"
    ids <- values$subject_ids %||% rownames(values$data) %||% as.character(seq_len(nrow(values$data)))
    gap_h <- if (isTRUE(is.finite(input$dlmo_gap_h))) input$dlmo_gap_h else 6
    prof <- dance_dlmo_profiles(values$data, times, ids,
                                night = if (nv != "_none_") values$covariates[[nv]],
                                group = if (gv != "_none_") values$covariates[[gv]], gap_h = gap_h)
    lloq <- input$dlmo_lloq
    opt <- list(threshold = thr, confirm = as.integer(input$dlmo_confirm %||% 1L),
                max_gap = input$dlmo_max_gap %||% 2, n_base = as.integer(input$dlmo_n_base %||% 3L),
                k_sd = input$dlmo_k_sd %||% 2, lloq = if (isTRUE(is.finite(lloq))) lloq else NA_real_,
                lloq_rule = input$dlmo_lloq_rule %||% "half",
                interval_limit = input$dlmo_interval_limit %||% 2,
                fine = isTRUE(input$dlmo_hockey_fine %||% TRUE))
    bed <- dlmo_bedtime()
    tab <- withProgress(message = "Estimating DLMO", value = 0, {
      dance_dlmo_run(prof, meth, opt, bed, progress = function(i, n)
        setProgress(i / n, detail = sprintf("profile %d of %d", i, n)))
    })
    anchor <- attr(tab, "anchor")
    prof <- lapply(prof, function(p) { p$t <- p$t - dance_dlmo_day_shift(min(p$t), anchor); p })
    list(ok = TRUE, tab = tab, profiles = prof, methods = meth, opt = opt, anchor = anchor,
         times = times, n_rows = nrow(values$data), gap_h = gap_h,
         time_src = input$dlmo_time_src %||% "_shared_", units = input$dlmo_units %||% "pg/mL",
         night_var = if (nv != "_none_") nv, group_var = if (gv != "_none_") gv,
         has_night_var = nv != "_none_", has_night = length(unique(tab$night)) > 1L,
         has_group = gv != "_none_", has_bed = !is.null(bed),
         bed_src = input$dlmo_bed_src %||% "none", bed_value = input$dlmo_bed_value,
         bed_var = input$dlmo_bed_var, bedtime = bed,
         dlmoR_version = if ("hockey" %in% meth) as.character(utils::packageVersion("dlmoR")))
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  values$dlmo_results <- res
  if (!isTRUE(res$ok)) showNotification(res$message, type = "error", duration = 10)
  else showNotification(sprintf("DLMO estimated in %d profile(s).", length(res$profiles)), type = "message")
})

dlmo_result <- reactive({
  r <- values$dlmo_results
  validate(need(!is.null(r), "Choose the settings and press 'Estimate DLMO'."),
           need(isTRUE(r$ok), r$message %||% "The last run failed."))
  r
})

# ---- 1. PROFILES ------------------------------------------------------------------
output$dlmo_profile_select_ui <- renderUI({
  r <- values$dlmo_results
  if (!isTRUE(r$ok)) return(NULL)
  ids <- unique(vapply(r$profiles, `[[`, "", "id"))
  gl <- if (isTRUE(r$has_group)) sort(unique(stats::na.omit(r$tab$group))) else character(0)
  ch <- c("All participants" = "_all_",
          stats::setNames(paste0("_g_", gl), sprintf("Group: %s", gl)), stats::setNames(ids, ids))
  fluidRow(
    column(6, selectInput("dlmo_profile_id", "Show:", choices = ch,
                          selected = isolate(input$dlmo_profile_id) %||% ids[1])),
    column(6, radioButtons("dlmo_summary_stat", "Onset summary (all participants / a group):",
                           choices = c("Mean" = "mean", "Median" = "median", "None" = "none"),
                           selected = isolate(input$dlmo_summary_stat) %||% "mean", inline = TRUE)))
})

DLMO_METHOD_DASH <- c(fixed = "solid", variable = "dash", hockey = "dot")

# What the plot shows: one participant, everyone, or one group; the cells it is
# coloured by; and the onset summary of those cells (server/08l_helpers_dlmo.R).
dlmo_plot_view <- reactive({
  r <- dlmo_result()
  sel <- input$dlmo_profile_id %||% r$profiles[[1]]$id
  grp_sel <- if (startsWith(sel, "_g_")) sub("^_g_", "", sel) else NULL
  many <- identical(sel, "_all_") || !is.null(grp_sel)
  pr <- if (identical(sel, "_all_")) r$profiles
        else if (!is.null(grp_sel)) Filter(function(p) identical(p$group, grp_sel), r$profiles)
        else Filter(function(p) identical(p$id, sel), r$profiles)
  # with everyone and a group variable, a cell is group x night; otherwise the night
  by_group <- isTRUE(r$has_group) && identical(sel, "_all_")
  cell_of <- function(p) if (by_group) paste(p$group %||% "(no group)", p$night, sep = " / ") else p$night
  stat <- input$dlmo_summary_stat %||% "mean"
  summ <- NULL
  if (many && stat != "none") {
    tb <- r$tab
    if (!is.null(grp_sel)) tb <- tb[!is.na(tb$group) & tb$group == grp_sel, , drop = FALSE]
    summ <- dance_dlmo_onset_summary(tb, stat, by_group = by_group)
  }
  list(r = r, sel = sel, grp = grp_sel, many = many, profiles = pr, by_group = by_group,
       cell_of = cell_of, stat = stat, summary = summ)
})

output$dlmo_profile_plot <- renderPlotly({
  v <- dlmo_plot_view(); r <- v$r; pr <- v$profiles; all <- v$many
  validate(need(length(pr) > 0, "No profile for this selection."))
  cells <- unique(vapply(if (v$by_group) r$profiles else pr, v$cell_of, ""))
  cols <- dance_group_colors(cells)
  names(cols) <- cells
  u <- r$units %||% "pg/mL"
  p <- plotly::plot_ly()
  shown <- character(0)
  ymax <- max(unlist(lapply(pr, function(q) q$y)), na.rm = TRUE)
  for (q in pr) {
    ce <- v$cell_of(q); cl <- cols[[ce]]
    ok <- is.finite(q$y)
    p <- plotly::add_trace(p, x = q$t[ok], y = q$y[ok], type = "scatter",
                           mode = if (all) "lines" else "lines+markers",
                           line = list(color = cl, width = if (all) 1 else 2), marker = list(color = cl),
                           opacity = if (all && !is.null(v$summary)) 0.5 else 1,
                           name = ce, legendgroup = ce, showlegend = !(ce %in% shown),
                           text = sprintf("%s, %s<br>%s: %.2f %s", q$id, ce,
                                          dance_clock_label(q$t[ok], 24, show_day = FALSE), q$y[ok], u),
                           hoverinfo = "text")
    shown <- c(shown, ce)
    rows <- r$tab[r$tab$id == q$id & r$tab$night == q$night, , drop = FALSE]
    for (j in seq_len(nrow(rows))) {
      rw <- rows[j, ]
      if (!all && rw$method == "variable" && is.finite(rw$threshold))
        p <- plotly::add_segments(p, x = min(q$t), xend = max(q$t), y = rw$threshold, yend = rw$threshold,
                                  line = list(color = cl, dash = "dash", width = 1), showlegend = FALSE,
                                  hoverinfo = "text", text = sprintf("variable threshold %.2f %s", rw$threshold, u))
      if (rw$status != "determined" || !is.finite(rw$dlmo_h)) next
      lab <- sprintf("%s, %s<br>DLMO (%s): %s", q$id, q$night, DANCE_DLMO_METHOD_LABEL[[rw$method]], rw$dlmo_clock)
      p <- if (all)
        plotly::add_markers(p, x = rw$dlmo_h, y = rw$threshold %||% r$opt$threshold, marker = list(color = cl, size = 7,
                            symbol = switch(rw$method, fixed = "circle", variable = "diamond", "triangle-up")),
                            showlegend = FALSE, hoverinfo = "text", text = lab)
      else plotly::add_segments(p, x = rw$dlmo_h, xend = rw$dlmo_h, y = 0, yend = max(q$y, na.rm = TRUE),
                                line = list(color = cl, dash = DLMO_METHOD_DASH[[rw$method]], width = 2),
                                showlegend = FALSE, hoverinfo = "text", text = lab)
    }
    if (!all && r$has_bed && is.finite(rows$bedtime_h[1])) {
      bt <- rows$bedtime_h[1] + 24 * round((mean(range(q$t)) - rows$bedtime_h[1]) / 24)
      p <- plotly::add_segments(p, x = bt, xend = bt, y = 0, yend = max(q$y, na.rm = TRUE),
                                line = list(color = "#555", dash = "dashdot", width = 1), showlegend = FALSE,
                                hoverinfo = "text", text = sprintf("habitual bedtime %s",
                                                                   dance_clock_label(bt, 24, show_day = FALSE)))
    }
  }
  # THE ONSET SUMMARY of each cell, one line per method in the method's dash
  sm <- v$summary
  if (!is.null(sm)) for (i in which(is.finite(sm$value))) {
    s <- sm[i, ]; cl <- cols[[s$cell]] %||% "#333"
    lab <- sprintf("%s onset, %s (%s): %s<br>from %d of %d profiles", tools::toTitleCase(s$stat), s$cell,
                   DANCE_DLMO_METHOD_LABEL[[s$method]], s$clock, s$n, s$n_profiles)
    p <- plotly::add_segments(p, x = s$value, xend = s$value, y = 0, yend = ymax,
                              line = list(color = cl, dash = DLMO_METHOD_DASH[[s$method]], width = 3),
                              name = sprintf("%s onset", tools::toTitleCase(s$stat)), legendgroup = "summary",
                              showlegend = FALSE, hoverinfo = "text", text = lab)
  }
  xr <- range(unlist(lapply(pr, `[[`, "t")), na.rm = TRUE)
  tk <- dance_clock_ticks(xr)
  shapes <- if ("fixed" %in% r$methods || "hockey" %in% r$methods)
    list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = r$opt$threshold, y1 = r$opt$threshold,
              line = list(color = "#b22222", dash = "dot", width = 1)))
  plotly::layout(p, shapes = shapes,
                 xaxis = c(list(title = "Clock time"), if (!is.null(tk)) list(tickvals = tk$vals, ticktext = tk$text)),
                 yaxis = list(title = sprintf("Melatonin (%s)", u), rangemode = "tozero"),
                 legend = list(title = list(text = if (v$by_group) "Group / night" else "Night")))
})

output$dlmo_profile_note <- renderUI({
  r <- values$dlmo_results
  if (!isTRUE(r$ok)) return(NULL)
  v <- tryCatch(dlmo_plot_view(), error = function(e) NULL)
  sm <- v$summary
  tagList(
    helpText(HTML(paste0(
      "Red dotted line: the threshold (", format(r$opt$threshold), " ", r$units, "). ",
      "Vertical lines: the onset by each method &mdash; solid fixed threshold, dashed variable threshold ",
      "(its per-profile threshold drawn as a thin dashed line), dotted hockey-stick.",
      if (r$has_bed) " Dash-dot: habitual bedtime." else "",
      " Each profile is drawn on the same evening clock. With all participants or a group, the markers ",
      "are the individual onsets and the thick lines each cell's ",
      if (!is.null(sm)) sm$stat[1] else "mean or median",
      " onset, of the determined onsets only (the same numbers as the Comparisons tab)."))),
    if (!is.null(sm) && nrow(sm))
      tags$table(style = "font-size:12px;margin-top:4px",
        tags$thead(tags$tr(lapply(c("Cell", "Method", tools::toTitleCase(sm$stat[1]), "Determined"),
                                  function(h) tags$th(style = "padding:2px 14px 2px 0", h)))),
        tags$tbody(lapply(seq_len(nrow(sm)), function(i) tags$tr(
          tags$td(style = "padding:2px 14px 2px 0", sm$cell[i]),
          tags$td(style = "padding:2px 14px 2px 0", DANCE_DLMO_METHOD_LABEL[[sm$method[i]]]),
          tags$td(style = "padding:2px 14px 2px 0;font-family:monospace", if (is.finite(sm$value[i])) sm$clock[i] else "--"),
          tags$td(style = "padding:2px 14px 2px 0", sprintf("%d of %d", sm$n[i], sm$n_profiles[i])))))))
}
)

# ---- 2. ONSETS --------------------------------------------------------------------
dlmo_table_frame <- function(r) {
  t <- r$tab
  data.frame(Participant = t$id, Night = t$night, Group = t$group,
             Method = unname(DANCE_DLMO_METHOD_LABEL[t$method]),
             Threshold = round(t$threshold, 3), DLMO = ifelse(t$status == "determined", t$dlmo_clock, ""),
             DLMO_h = round(t$dlmo_h, 3), Status = unname(DANCE_DLMO_STATUS_LABEL[t$status]),
             Samples = t$n_samples, Below_LLOQ = t$n_below_lloq,
             Bedtime = if (r$has_bed) dance_clock_label(t$bedtime_h, 24, show_day = FALSE) else NA,
             Phase_angle_h = if (r$has_bed) round(t$phase_angle_h, 3) else NA,
             Note = t$note, stringsAsFactors = FALSE, check.names = FALSE)
}

output$dlmo_summary <- renderText({
  r <- dlmo_result()
  st <- lapply(r$methods, function(m) {
    s <- table(factor(r$tab$status[r$tab$method == m], levels = names(DANCE_DLMO_STATUS_LABEL)))
    s <- s[s > 0]
    sprintf("  %s: %s", DANCE_DLMO_METHOD_LABEL[[m]],
            paste(sprintf("%d %s", as.integer(s), DANCE_DLMO_STATUS_LABEL[names(s)]), collapse = "; "))
  })
  paste(c("WHAT THIS RUN DID", paste0("  ", dance_dlmo_settings_text(r)), "", "ONSETS", unlist(st)),
        collapse = "\n")
})

output$dlmo_table <- renderTable({
  r <- dlmo_result()
  d <- dlmo_table_frame(r)
  keep <- vapply(d, function(x) !all(is.na(x) | x == ""), logical(1))
  d[, keep, drop = FALSE]
}, striped = TRUE, na = "")

output$dlmo_download <- downloadHandler(
  filename = function() sprintf("dlmo_onsets_%s.csv", Sys.Date()),
  content = function(file) {
    r <- values$dlmo_results
    req(isTRUE(r$ok))
    utils::write.csv(dlmo_table_frame(r), file, row.names = FALSE, na = "")
  })

# ---- 3. COMPARISONS ---------------------------------------------------------------
dlmo_desc_frame <- function(ds, clock = TRUE) {
  f <- function(x) if (clock) dance_clock_label(x, 24, show_day = FALSE) else sprintf("%.2f", x)
  data.frame(Cell = ds$cell, Profiles = ds$n_profiles, Determined = ds$n_determined,
             `Before window` = ds$n_before, `No onset` = ds$n_after, Other = ds$n_other,
             Mean = f(ds$mean), `SD (h)` = sprintf("%.2f", ds$sd), Median = f(ds$median),
             Range = paste(f(ds$min), "to", f(ds$max)), check.names = FALSE, stringsAsFactors = FALSE)
}
dlmo_test_frame <- function(cp) {
  tt <- cp$tests
  data.frame(Effect = tt$effect, Test = tt$test,
             Difference = ifelse(is.finite(tt$estimate), dance_dlmo_fmt_diff(tt$estimate), ""),
             `95% CI (h)` = ifelse(is.finite(tt$lo), sprintf("[%.2f, %.2f]", tt$lo, tt$hi), ""),
             Statistic = vapply(seq_len(nrow(tt)), function(i) dance_dlmo_test_text(tt[i, ]), ""),
             n = tt$n, check.names = FALSE, stringsAsFactors = FALSE) -> out
  # the Holm-adjusted p of the pairwise contrasts, when there are any
  if (any(is.finite(tt$p_adj)))
    out[["p (Holm)"]] <- ifelse(is.finite(tt$p_adj), sub("^0", "", sprintf("%.3f", tt$p_adj)), "")
  out
}

output$dlmo_compare_ui <- renderUI({
  r <- dlmo_result()
  tagList(
    # the equation of each method and of the phase angle (server/08n_helpers_equations.R)
    tags$div(style = "margin-bottom:8px;color:#555", "Equations:",
             lapply(r$methods, function(m) tags$span(style = "margin-right:12px",
               dance_eq_link(paste0("dlmo_", m), DANCE_DLMO_METHOD_LABEL[[m]]))),
             if (r$has_bed) dance_eq_link("phase_angle", "Phase angle")),
    if (length(r$methods) > 1L)
      selectInput("dlmo_compare_method", "Method:", choices = setNames(r$methods, DANCE_DLMO_METHOD_LABEL[r$methods]),
                  selected = isolate(input$dlmo_compare_method)),
    h4("DLMO"), uiOutput("dlmo_design_note"), tableOutput("dlmo_desc_table"), uiOutput("dlmo_test_note"),
    tableOutput("dlmo_test_table"),
    if (r$has_bed) tagList(h4("Phase angle (bedtime - DLMO, h)"), tableOutput("dlmo_pa_desc_table"),
                           uiOutput("dlmo_pa_test_note"), tableOutput("dlmo_pa_test_table")),
    if (length(r$methods) > 1L) tagList(h4("Agreement between methods"), tableOutput("dlmo_agree_table")),
    helpText(HTML(paste(
      "Only determined onsets enter the means and tests. If being undetermined depends on the condition",
      "(an onset pushed past the window, for instance), the comparison is biased; the counts per cell",
      "are reported for that reason. Nights / conditions are compared within participant (the",
      "participant identifier must repeat across a participant's rows): a paired t test for two",
      "conditions, otherwise a linear mixed model with a random intercept per participant (lmerTest,",
      "Satterthwaite df), with every pair of conditions as a contrast, Holm-adjusted within each group.",
      "Groups alone: Welch's t test, or Welch's one-way ANOVA."))))
})

dlmo_cmp <- reactive({
  r <- dlmo_result()
  m <- input$dlmo_compare_method %||% r$methods[1]
  if (!m %in% r$methods) m <- r$methods[1]
  dance_dlmo_comparisons(r, m)
})
output$dlmo_design_note <- renderUI(dance_status(dlmo_cmp()$design))
output$dlmo_desc_table <- renderTable(dlmo_desc_frame(dlmo_cmp()$desc), na = "")
output$dlmo_test_note <- renderUI({
  cp <- dlmo_cmp()$test
  dance_status(if (isTRUE(cp$ok)) paste0("Test: ", cp$method,
                                     if (!is.null(cp$n_pairs)) sprintf(" (%d participants with an onset on both nights; %d without)",
                                                                       cp$n_pairs, cp$n_unpaired) else "")
           else cp$message)
})
output$dlmo_test_table <- renderTable({
  cp <- dlmo_cmp()$test
  if (isTRUE(cp$ok)) dlmo_test_frame(cp)
}, na = "")
output$dlmo_pa_desc_table <- renderTable({
  pd <- dlmo_cmp()$pa_desc
  if (!is.null(pd)) dlmo_desc_frame(pd, clock = FALSE)
}, na = "")
output$dlmo_pa_test_note <- renderUI({
  cp <- dlmo_cmp()$pa_test
  if (!is.null(cp)) dance_status(if (isTRUE(cp$ok)) paste("Test:", cp$method) else cp$message)
})
output$dlmo_pa_test_table <- renderTable({
  cp <- dlmo_cmp()$pa_test
  if (isTRUE(cp$ok)) dlmo_test_frame(cp)
}, na = "")
output$dlmo_agree_table <- renderTable({
  ag <- dance_dlmo_agreement(dlmo_result()$tab)
  if (!is.null(ag)) data.frame(Comparison = ag$comparison, `Profiles with both` = ag$n_both,
                               `Mean difference` = dance_dlmo_fmt_diff(ag$mean_diff_h),
                               `SD (h)` = sprintf("%.2f", ag$sd_diff_h), `Largest |difference| (h)` = sprintf("%.2f", ag$max_abs_diff_h),
                               check.names = FALSE)
}, na = "")
