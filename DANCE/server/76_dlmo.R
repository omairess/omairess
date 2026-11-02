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
    selectInput("dlmo_night_var", "Night / condition column (one row per night):",
                choices = c("None: one row per participant" = "_none_", cv),
                selected = isolate(input$dlmo_night_var) %||% "_none_"),
    selectInput("dlmo_group_var", "Group (between participants):",
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
    helpText(sprintf("Hockey-stick: dlmoR %s found.", as.character(utils::packageVersion("dlmoR"))))
  else helpText(HTML(paste0(
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
  selectInput("dlmo_profile_id", "Participant:", choices = c("All participants" = "_all_", ids),
              selected = isolate(input$dlmo_profile_id) %||% ids[1])
})

DLMO_METHOD_DASH <- c(fixed = "solid", variable = "dash", hockey = "dot")

output$dlmo_profile_plot <- renderPlotly({
  r <- dlmo_result()
  sel <- input$dlmo_profile_id %||% r$profiles[[1]]$id
  all <- identical(sel, "_all_")
  pr <- if (all) r$profiles else Filter(function(p) identical(p$id, sel), r$profiles)
  validate(need(length(pr) > 0, "No profile for this participant."))
  nights <- unique(vapply(r$profiles, `[[`, "", "night"))
  cols <- dance_group_colors(nights)
  names(cols) <- nights
  u <- r$units %||% "pg/mL"
  p <- plotly::plot_ly()
  shown <- character(0)
  for (q in pr) {
    cl <- cols[[q$night]]
    ok <- is.finite(q$y)
    p <- plotly::add_trace(p, x = q$t[ok], y = q$y[ok], type = "scatter",
                           mode = if (all) "lines" else "lines+markers",
                           line = list(color = cl, width = if (all) 1 else 2), marker = list(color = cl),
                           name = q$night, legendgroup = q$night, showlegend = !(q$night %in% shown),
                           text = sprintf("%s, %s<br>%s: %.2f %s", q$id, q$night,
                                          dance_clock_label(q$t[ok], 24, show_day = FALSE), q$y[ok], u),
                           hoverinfo = "text")
    shown <- c(shown, q$night)
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
  xr <- range(unlist(lapply(pr, `[[`, "t")), na.rm = TRUE)
  tk <- dance_clock_ticks(xr)
  shapes <- if ("fixed" %in% r$methods || "hockey" %in% r$methods)
    list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = r$opt$threshold, y1 = r$opt$threshold,
              line = list(color = "#b22222", dash = "dot", width = 1)))
  plotly::layout(p, shapes = shapes,
                 xaxis = c(list(title = "Clock time"), if (!is.null(tk)) list(tickvals = tk$vals, ticktext = tk$text)),
                 yaxis = list(title = sprintf("Melatonin (%s)", u), rangemode = "tozero"),
                 legend = list(title = list(text = "Night")))
})

output$dlmo_profile_note <- renderUI({
  r <- values$dlmo_results
  if (!isTRUE(r$ok)) return(NULL)
  helpText(HTML(paste0(
    "Red dotted line: the threshold (", format(r$opt$threshold), " ", r$units, "). ",
    "Vertical lines: the onset by each method &mdash; solid fixed threshold, dashed variable threshold ",
    "(its per-profile threshold drawn as a thin dashed line), dotted hockey-stick.",
    if (r$has_bed) " Dash-dot: habitual bedtime." else "",
    " Each profile is drawn on the same evening clock. With all participants, the markers are the onsets.")))
})

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
             n = tt$n, check.names = FALSE, stringsAsFactors = FALSE)
}

output$dlmo_compare_ui <- renderUI({
  r <- dlmo_result()
  tagList(
    if (length(r$methods) > 1L)
      selectInput("dlmo_compare_method", "Method:", choices = setNames(r$methods, DANCE_DLMO_METHOD_LABEL[r$methods]),
                  selected = isolate(input$dlmo_compare_method)),
    h4("DLMO"), tableOutput("dlmo_desc_table"), uiOutput("dlmo_test_note"), tableOutput("dlmo_test_table"),
    if (r$has_bed) tagList(h4("Phase angle (bedtime - DLMO, h)"), tableOutput("dlmo_pa_desc_table"),
                           uiOutput("dlmo_pa_test_note"), tableOutput("dlmo_pa_test_table")),
    if (length(r$methods) > 1L) tagList(h4("Agreement between methods"), tableOutput("dlmo_agree_table")),
    helpText(HTML(paste(
      "Only determined onsets enter the means and tests. If being undetermined depends on the condition",
      "(an onset pushed past the window, for instance), the comparison is biased; the counts per cell",
      "are reported for that reason. Nights are compared within participant: a paired t test for two",
      "nights, otherwise a linear mixed model with a random intercept per participant (lmerTest,",
      "Satterthwaite df). Groups alone: Welch's t test, or Welch's one-way ANOVA."))))
})

dlmo_cmp <- reactive({
  r <- dlmo_result()
  m <- input$dlmo_compare_method %||% r$methods[1]
  if (!m %in% r$methods) m <- r$methods[1]
  dance_dlmo_comparisons(r, m)
})
output$dlmo_desc_table <- renderTable(dlmo_desc_frame(dlmo_cmp()$desc), na = "")
output$dlmo_test_note <- renderUI({
  cp <- dlmo_cmp()$test
  helpText(if (isTRUE(cp$ok)) paste0("Test: ", cp$method,
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
  if (!is.null(cp)) helpText(if (isTRUE(cp$ok)) paste("Test:", cp$method) else cp$message)
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
