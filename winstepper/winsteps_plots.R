# =============================================================================
# winsteps_plots.R  --  WINSTEPS-style figures, base graphics.
#
# Every plot takes a `style` list so that colours (house rule 4), sizes
# (house rule 5) and the pastel default palette (house rule 9) are controlled
# from the UI, and every plot can be re-drawn into any graphics device for
# export at arbitrary size / resolution.
# =============================================================================

#' Default (pastel) style. Overridden by the UI colour pickers / sliders.
ws_style <- function(...) {
  s <- list(
    col_person   = "#B3CDE3",   # Pastel1 blue
    col_item     = "#FBB4AE",   # Pastel1 red
    col_border   = "#666666",
    col_line     = "#4D4D4D",
    col_ref      = "#CCCCCC",
    col_pos      = "#CCEBC5",
    col_neg      = "#FDDAEC",
    palette      = "Pastel1",
    cex_label    = 0.85,
    cex_point    = 1.2,
    lwd          = 2,
    bins         = 30,
    show_grid    = TRUE,
    show_msT     = TRUE       # Winsteps M / S / T markers
  )
  m <- list(...)
  s[names(m)] <- m
  s
}

.pal <- function(n, which = "Pastel1") {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    mx <- RColorBrewer::brewer.pal.info[which, "maxcolors"]
    base <- RColorBrewer::brewer.pal(mx, which)
    if (n <= mx) return(base[seq_len(n)])
    return(grDevices::colorRampPalette(base)(n))
  }
  grDevices::hcl.colors(n, "Pastel 1")
}

.grid <- function(style, h = TRUE) {
  if (isTRUE(style$show_grid))
    graphics::grid(nx = if (h) NA else NULL, ny = if (h) NULL else NA,
                   col = style$col_ref, lty = 3)
}

# ---------------------------------------------------------------------------
# Wright map (person-item map; WINSTEPS Tables 1 / 12 / 16)
# ---------------------------------------------------------------------------

#' @param what "items" plots item difficulties, "thresholds" plots
#'   Rasch-Thurstone thresholds (item.category), "expected" plots the measures
#'   at which the expected score on the item is at each half-score point.
plot_wright <- function(fit, style = ws_style(), what = c("items", "thresholds"),
                        item_labels = NULL, max_labels = 60, main = "Wright map (person-item map)") {
  what <- match.arg(what)
  pm <- fit$theta[fit$keep_p & !is.na(fit$theta)]
  if (what == "items") {
    im <- fit$delta[fit$keep_i & !is.na(fit$delta)]
    lab <- (item_labels %||% fit$item_id)[fit$keep_i & !is.na(fit$delta)]
  } else {
    td <- threshold_data(fit)
    im <- td$measure; lab <- td$label
  }
  rng <- range(c(pm, im), na.rm = TRUE)
  rng <- rng + c(-1, 1) * max(diff(rng) / 20, 0.2)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, nrow = 1), widths = c(1, 1.25))
  graphics::par(mar = c(4, 4.2, 3.5, 0.5), oma = c(0, 0, 2, 0))

  ## --- left: person distribution (same logit scale as the right panel) -----
  brks <- seq(rng[1], rng[2], length.out = max(style$bins, 5) + 1)
  h <- graphics::hist(pm, breaks = brks, plot = FALSE)
  xmax <- max(h$counts) * 1.06
  graphics::plot(NA, xlim = c(xmax, 0), ylim = rng, axes = FALSE,
                 xlab = "", ylab = "", main = "PERSONS",
                 cex.main = style$cex_label + 0.15, xaxs = "i")
  if (isTRUE(style$show_grid))
    graphics::abline(h = pretty(rng, 8), col = style$col_ref, lty = 3)
  graphics::rect(0, h$breaks[-length(h$breaks)], h$counts, h$breaks[-1],
                 col = style$col_person, border = style$col_border)
  graphics::axis(2, at = pretty(rng, 8), las = 1, cex.axis = style$cex_label)
  graphics::axis(1, at = pretty(c(0, xmax), 4), cex.axis = style$cex_label * 0.9)
  graphics::mtext("logits", side = 2, line = 2.8, cex = style$cex_label)
  graphics::mtext("persons (count)", side = 1, line = 2.3, cex = style$cex_label * 0.9)
  if (isTRUE(style$show_msT)) {
    mm <- mean(pm); ss <- stats::sd(pm)
    graphics::abline(h = mm, col = style$col_line, lwd = style$lwd)
    graphics::abline(h = c(mm - ss, mm + ss), col = style$col_line, lty = 2)
    graphics::abline(h = c(mm - 2 * ss, mm + 2 * ss), col = style$col_line, lty = 3)
    graphics::mtext(c("T", "S", "M", "S", "T"), side = 2, line = -0.9, las = 1,
                    at = c(mm - 2 * ss, mm - ss, mm, mm + ss, mm + 2 * ss),
                    cex = style$cex_label * 0.8)
  }

  ## --- right: item ladder ---------------------------------------------------
  graphics::par(mar = c(4, 0.5, 3.5, 2.6))
  graphics::plot(NA, xlim = c(0, 1), ylim = rng, axes = FALSE, xlab = "", ylab = "",
                 main = if (what == "items") "ITEMS" else "ITEM THRESHOLDS",
                 cex.main = style$cex_label + 0.15)
  graphics::axis(4, at = pretty(rng, 8), las = 1, cex.axis = style$cex_label)
  if (isTRUE(style$show_grid))
    graphics::abline(h = pretty(rng, 8), col = style$col_ref, lty = 3)
  ord <- order(im)
  im <- im[ord]; lab <- lab[ord]
  if (length(im) > max_labels) {
    keep <- round(seq(1, length(im), length.out = max_labels))
    im <- im[keep]; lab <- lab[keep]
  }
  # stack labels that would overlap
  spacing <- diff(rng) / 45
  xoff <- numeric(length(im)); last <- -Inf; col_i <- 0
  for (j in seq_along(im)) {
    if (im[j] - last < spacing) col_i <- col_i + 1 else col_i <- 0
    xoff[j] <- col_i; last <- im[j]
  }
  xpos <- 0.06 + xoff * 0.12
  graphics::points(xpos, im, pch = 22, bg = style$col_item, col = style$col_border,
                   cex = style$cex_point)
  graphics::text(xpos + 0.03, im, lab, adj = 0, cex = style$cex_label * 0.9)
  graphics::mtext("logits", side = 4, line = 2.4, cex = style$cex_label)
  graphics::mtext(main, side = 3, line = 0.2, outer = TRUE, cex = style$cex_label + 0.35, font = 2)
  if (isTRUE(style$show_msT)) {
    mm <- mean(im); ss <- stats::sd(im)
    graphics::abline(h = mm, col = style$col_line, lwd = style$lwd, lty = 1)
    graphics::abline(h = c(mm - ss, mm + ss), col = style$col_line, lty = 2)
    graphics::abline(h = c(mm - 2 * ss, mm + 2 * ss), col = style$col_line, lty = 3)
    graphics::mtext(c("T", "S", "M", "S", "T"), side = 4, line = 0.2,
                    at = c(mm - 2 * ss, mm - ss, mm, mm + ss, mm + 2 * ss),
                    cex = style$cex_label * 0.8, las = 1)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Rating-scale functioning graphs (WINSTEPS Graphs menu)
# ---------------------------------------------------------------------------

#' Category probability curves for one item group (optionally at one item's
#' difficulty). This is WINSTEPS "Category probability curves".
plot_category_probs <- function(fit, group, style = ws_style(), item = NULL,
                                from = -6, to = 6, show_thresholds = TRUE) {
  de <- if (is.null(item)) 0 else fit$delta[if (is.character(item)) match(item, fit$item_id) else item]
  if (is.na(de)) de <- 0
  cd <- curve_data(fit, group, delta = de, from = from + de, to = to + de)
  m <- ncol(cd$prob) - 1
  cols <- .pal(m + 1, style$palette)
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(NA, xlim = range(cd$measure), ylim = c(0, 1),
                 xlab = "Person measure - item difficulty (logits)",
                 ylab = "Category probability",
                 main = sprintf("Category probability curves - group '%s'%s", group,
                                if (is.null(item)) "" else paste0(" (item ", item, ")")),
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  for (k in 0:m)
    graphics::lines(cd$measure, cd$prob[, k + 1], col = cols[k + 1], lwd = style$lwd + 0.6)
  if (show_thresholds && m >= 2) {
    tau <- fit$tau[[group]]
    graphics::abline(v = de + tau, col = style$col_ref, lty = 2)
    graphics::axis(3, at = de + tau, labels = paste0("t", seq_along(tau)),
                   cex.axis = style$cex_label * 0.85, tick = FALSE, line = -0.8)
  }
  graphics::legend("top", legend = paste("Category", 0:m), col = cols, lwd = style$lwd + 0.6,
                   horiz = TRUE, bty = "n", cex = style$cex_label * 0.9)
  invisible(NULL)
}

#' Cumulative category probabilities P(X >= k), with the Rasch-Thurstone
#' 50% thresholds marked (WINSTEPS "Cumulative probabilities").
plot_cumulative <- function(fit, group, style = ws_style(), item = NULL,
                            from = -6, to = 6) {
  de <- if (is.null(item)) 0 else fit$delta[if (is.character(item)) match(item, fit$item_id) else item]
  if (is.na(de)) de <- 0
  cd <- curve_data(fit, group, delta = de, from = from + de, to = to + de)
  m <- ncol(cd$cum)
  cols <- .pal(max(m, 3), style$palette)[seq_len(m)]
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(NA, xlim = range(cd$measure), ylim = c(0, 1),
                 xlab = "Measure (logits)", ylab = "P(X >= k)",
                 main = sprintf("Cumulative category probabilities - group '%s'", group),
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  graphics::abline(h = 0.5, col = style$col_line, lty = 2, lwd = style$lwd)
  for (k in seq_len(m))
    graphics::lines(cd$measure, cd$cum[, k], col = cols[k], lwd = style$lwd + 0.6)
  th <- vapply(seq_len(m), function(k) .thurstone(k, fit$tau[[group]]), numeric(1)) + de
  graphics::points(th, rep(0.5, m), pch = 21, bg = style$col_item,
                   col = style$col_border, cex = style$cex_point)
  graphics::legend("topright", legend = paste0("P(X>=", seq_len(m), ")"), col = cols,
                   lwd = style$lwd + 0.6, bty = "n", cex = style$cex_label * 0.9)
  invisible(NULL)
}

#' Model expected score curve with empirical ICC overlay (WINSTEPS
#' "Expected score ICC" + "Empirical ICC").
plot_expected_score <- function(fit, item, style = ws_style(), nbins = 10,
                                from = -6, to = 6, ci = TRUE) {
  i <- if (is.character(item)) match(item, fit$item_id) else item
  g <- fit$groups[i]; de <- fit$delta[i]
  cd <- curve_data(fit, g, delta = de, from = de + from, to = de + to)
  m  <- fit$max_cat[[g]]
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(NA, xlim = range(cd$measure), ylim = c(0, m),
                 xlab = "Person measure (logits)", ylab = "Score on item",
                 main = sprintf("Expected score curve - item %s (measure %.2f)", fit$item_id[i], de),
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  if (ci) {
    se <- sqrt(cd$information)
    graphics::polygon(c(cd$measure, rev(cd$measure)),
                      pmin(pmax(c(cd$expected + se, rev(cd$expected - se)), 0), m),
                      col = grDevices::adjustcolor(style$col_item, alpha.f = .3), border = NA)
  }
  graphics::lines(cd$measure, cd$expected, col = style$col_line, lwd = style$lwd + 1)
  graphics::abline(v = de, col = style$col_ref, lty = 2)
  ei <- empirical_icc(fit, i, nbins)
  if (!is.null(ei)) {
    graphics::points(ei$measure, ei$observed, pch = 21, bg = style$col_person,
                     col = style$col_border, cex = style$cex_point * sqrt(ei$n / mean(ei$n)))
    graphics::lines(ei$measure, ei$observed, col = style$col_person, lwd = style$lwd, lty = 2)
  }
  graphics::legend("topleft", c("Model expected score", "Empirical (binned)", "+/- 1 model SD"),
                   col = c(style$col_line, style$col_person, style$col_item),
                   lwd = c(style$lwd + 1, style$lwd, 6), lty = c(1, 2, 1),
                   bty = "n", cex = style$cex_label * 0.85)
  invisible(NULL)
}

#' Item information function(s).
plot_item_info <- function(fit, items = NULL, style = ws_style(), from = -6, to = 6) {
  idx <- if (is.null(items)) which(!is.na(fit$delta)) else
    if (is.character(items)) match(items, fit$item_id) else items
  x <- seq(from, to, length.out = 401)
  cols <- .pal(max(length(idx), 3), style$palette)
  info <- lapply(idx, function(i) {
    P <- .probs_group(x, fit$delta[i], fit$tau[[fit$groups[i]]])
    .moments(P)$W[, 1]
  })
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(NA, xlim = range(x), ylim = c(0, max(unlist(info)) * 1.05),
                 xlab = "Person measure (logits)", ylab = "Item information",
                 main = "Item information functions",
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  for (j in seq_along(idx))
    graphics::lines(x, info[[j]], col = cols[(j - 1) %% length(cols) + 1], lwd = style$lwd)
  if (length(idx) <= 12)
    graphics::legend("topright", fit$item_id[idx], col = cols[seq_along(idx)],
                     lwd = style$lwd, bty = "n", cex = style$cex_label * 0.8)
  invisible(NULL)
}

#' Test characteristic curve and test information (WINSTEPS Table 20 graphs).
plot_test_curves <- function(fit, style = ws_style(), what = c("both", "tcc", "info"),
                             from = -6, to = 6) {
  what <- match.arg(what)
  tc <- test_curves(fit, from, to)
  if (what == "both") { op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op)) }
  if (what %in% c("both", "tcc")) {
    graphics::par(mar = c(4.5, 4.5, 3.5, 1))
    graphics::plot(tc$measure, tc$expected_score, type = "l", col = style$col_line,
                   lwd = style$lwd + 1, xlab = "Person measure (logits)",
                   ylab = "Expected raw score", main = "Test characteristic curve",
                   cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
    .grid(style, h = FALSE); .grid(style, h = TRUE)
  }
  if (what %in% c("both", "info")) {
    graphics::par(mar = c(4.5, 4.5, 3.5, 4.5))
    graphics::plot(tc$measure, tc$information, type = "l", col = style$col_line,
                   lwd = style$lwd + 1, xlab = "Person measure (logits)",
                   ylab = "Test information", main = "Test information & standard error",
                   cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
    .grid(style, h = FALSE); .grid(style, h = TRUE)
    graphics::par(new = TRUE)
    graphics::plot(tc$measure, tc$SE, type = "l", col = style$col_item, lwd = style$lwd,
                   lty = 2, axes = FALSE, xlab = "", ylab = "",
                   ylim = c(0, stats::quantile(tc$SE, .95)))
    graphics::axis(4, cex.axis = style$cex_label); graphics::mtext("SE (logits)", side = 4, line = 2.6,
                                                                  cex = style$cex_label)
    graphics::legend("topright", c("Information", "SE"), col = c(style$col_line, style$col_item),
                     lty = c(1, 2), lwd = style$lwd, bty = "n", cex = style$cex_label * 0.85)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Fit diagnostics
# ---------------------------------------------------------------------------

#' Bond & Fox "pathway" bubble chart: measure vs fit, bubble size = SE.
plot_pathway <- function(fit, style = ws_style(), margin = c("items", "persons"),
                         stat = c("infit_zstd", "outfit_zstd", "infit_mnsq", "outfit_mnsq")) {
  margin <- match.arg(margin); stat <- match.arg(stat)
  tab <- if (margin == "items") item_table(fit, FALSE) else person_table(fit)
  lab <- if (margin == "items") tab$Item else tab$Person
  keep <- if (margin == "items") fit$keep_i else fit$keep_p
  y <- switch(stat,
              infit_zstd = tab$Infit_ZSTD, outfit_zstd = tab$Outfit_ZSTD,
              infit_mnsq = tab$Infit_MNSQ, outfit_mnsq = tab$Outfit_MNSQ)
  x <- tab$Measure; se <- tab$Model_SE
  ok <- keep & is.finite(x) & is.finite(y)
  cexv <- style$cex_point * (0.6 + 2.2 * (se[ok] - min(se[ok], na.rm = TRUE)) /
                               max(diff(range(se[ok], na.rm = TRUE)), 1e-9))
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(x[ok], y[ok], pch = 21,
                 bg = if (margin == "items") style$col_item else style$col_person,
                 col = style$col_border, cex = cexv,
                 xlab = "Measure (logits)",
                 ylab = switch(stat, infit_zstd = "Infit ZSTD", outfit_zstd = "Outfit ZSTD",
                               infit_mnsq = "Infit MNSQ", outfit_mnsq = "Outfit MNSQ"),
                 main = sprintf("Pathway (bubble) chart - %s; bubble size = S.E.", margin),
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  ref <- if (grepl("zstd", stat)) c(-2, 2) else c(0.5, 1.5)
  graphics::abline(h = ref, col = style$col_line, lty = 2, lwd = style$lwd)
  graphics::abline(h = if (grepl("zstd", stat)) 0 else 1, col = style$col_ref)
  if (sum(ok) <= 80)
    graphics::text(x[ok], y[ok], lab[ok], pos = 4, cex = style$cex_label * 0.75, offset = 0.4)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Dimensionality
# ---------------------------------------------------------------------------

#' Scree plot of the residual contrasts, with the random-level reference band
#' (Smith & Miao 1994: 1.4; Raiche 2005: up to ~2.0).
plot_scree <- function(pca, style = ws_style(), simulated = NULL) {
  ev <- pca$eigenvalues[seq_len(min(10, length(pca$eigenvalues)))]
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(seq_along(ev), ev, type = "b", pch = 21, bg = style$col_item,
                 col = style$col_border, lwd = style$lwd, cex = style$cex_point,
                 xlab = "Contrast (PCA component of standardized residuals)",
                 ylab = "Eigenvalue (item units)",
                 ylim = c(0, max(ev, 2.2, simulated$P95, na.rm = TRUE) * 1.1),
                 main = "Scree plot of residual contrasts",
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = TRUE)
  graphics::rect(0, 1.4, length(ev) + 1, 2.0,
                 col = grDevices::adjustcolor(style$col_ref, alpha.f = .45), border = NA)
  graphics::abline(h = 2, col = style$col_line, lty = 2, lwd = style$lwd)
  if (!is.null(simulated))
    graphics::lines(simulated$Contrast, simulated$P95, col = style$col_person,
                    lwd = style$lwd, lty = 3, type = "b", pch = 4)
  graphics::legend("topright",
                   c("Observed contrast", "Random level (1.4-2.0)",
                     if (!is.null(simulated)) "Simulated 95th pct"),
                   col = c(style$col_item, style$col_ref,
                           if (!is.null(simulated)) style$col_person),
                   pch = c(21, 15, if (!is.null(simulated)) 4), bty = "n",
                   cex = style$cex_label * 0.85)
  invisible(NULL)
}

#' Contrast loading plot (WINSTEPS Table 23.2): loading vs item measure.
plot_pca_contrast <- function(pca, contrast = 1, style = ws_style()) {
  ld <- pca$loadings
  y <- ld[[paste0("Contrast", contrast)]]
  cl <- ld$Cluster
  cols <- .pal(3, style$palette)
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(ld$Measure, y, pch = 21, bg = cols[cl], col = style$col_border,
                 cex = style$cex_point + 0.3,
                 xlab = "Item measure (logits)",
                 ylab = sprintf("Loading on contrast %d", contrast),
                 main = sprintf("Contrast %d loading plot (eigenvalue %.2f)",
                                contrast, pca$eigenvalues[contrast]),
                 ylim = range(c(y, -0.6, 0.6)),
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = FALSE); .grid(style, h = TRUE)
  graphics::abline(h = 0, col = style$col_line, lwd = style$lwd)
  graphics::abline(h = c(-0.4, 0.4), col = style$col_ref, lty = 2)
  graphics::text(ld$Measure, y, ld$Item, pos = 4, cex = style$cex_label * 0.75, offset = 0.4)
  graphics::legend("bottomright", paste("Cluster", 1:3), pt.bg = cols, pch = 21,
                   bty = "n", cex = style$cex_label * 0.85)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# DIF
# ---------------------------------------------------------------------------

#' DIF measures by person class (WINSTEPS Table 30 plot).
plot_dif <- function(dif, style = ws_style(), errbars = TRUE) {
  lv <- unique(c(dif$Class_A, dif$Class_B))
  loc <- attr(dif, "local_measures")
  se  <- attr(dif, "local_se")
  items <- rownames(loc)
  cols <- .pal(max(length(lv), 3), style$palette)[seq_along(lv)]
  graphics::par(mar = c(7, 4.5, 3.5, 1))
  yl <- range(c(loc - 2 * se, loc + 2 * se), na.rm = TRUE)
  graphics::plot(NA, xlim = c(0.5, nrow(loc) + 0.5), ylim = yl, xaxt = "n",
                 xlab = "", ylab = "Local item difficulty (logits)",
                 main = "DIF measures by person class (Rasch-Welch)",
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  .grid(style, h = TRUE)
  graphics::axis(1, at = seq_len(nrow(loc)), labels = items, las = 2,
                 cex.axis = style$cex_label * 0.85)
  for (j in seq_along(lv)) {
    x <- seq_len(nrow(loc)) + (j - (length(lv) + 1) / 2) * 0.12
    if (errbars)
      graphics::segments(x, loc[, lv[j]] - 1.96 * se[, lv[j]],
                         x, loc[, lv[j]] + 1.96 * se[, lv[j]],
                         col = style$col_border)
    graphics::points(x, loc[, lv[j]], pch = 21, bg = cols[j], col = style$col_border,
                     cex = style$cex_point)
    graphics::lines(x, loc[, lv[j]], col = cols[j], lwd = style$lwd, lty = 3)
  }
  graphics::legend("topleft", lv, pt.bg = cols, pch = 21, bty = "n",
                   cex = style$cex_label * 0.85, horiz = TRUE)
  invisible(NULL)
}

#' DIF contrast (size) plot with the ETS B / C reference lines.
plot_dif_contrast <- function(dif, style = ws_style()) {
  d <- dif[order(dif$Entry), ]
  cols <- c(A = style$col_ref, B = style$col_person, C = style$col_item)
  graphics::par(mar = c(7, 4.5, 3.5, 1))
  bp <- graphics::barplot(d$DIF_Contrast, col = cols[d$ETS_Class], border = style$col_border,
                          ylab = "DIF contrast (logits)", names.arg = d$Item, las = 2,
                          cex.names = style$cex_label * 0.85,
                          main = sprintf("DIF contrast: %s minus %s (ETS classification)",
                                         d$Class_A[1], d$Class_B[1]),
                          cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2,
                          ylim = range(c(d$DIF_Contrast, -0.8, 0.8), na.rm = TRUE))
  graphics::abline(h = c(-0.64, -0.43, 0, 0.43, 0.64),
                   col = c(style$col_line, style$col_ref, "black", style$col_ref, style$col_line),
                   lty = c(2, 3, 1, 3, 2), lwd = style$lwd)
  graphics::legend("topleft", c("A negligible", "B slight-moderate", "C moderate-large"),
                   fill = cols, bty = "n", cex = style$cex_label * 0.85, horiz = TRUE)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Category structure summary figure (Table 3.2 companion)
# ---------------------------------------------------------------------------

#' Observed vs expected average measure per category, and threshold ordering.
plot_category_diagnostics <- function(ct, group, style = ws_style()) {
  d <- ct[ct$Group == group, ]
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1)); on.exit(graphics::par(op))
  graphics::plot(d$Category, d$Obsvd_Avrge, type = "b", pch = 21, bg = style$col_person,
                 col = style$col_border, lwd = style$lwd, cex = style$cex_point,
                 ylim = range(c(d$Obsvd_Avrge, d$Sample_Expect), na.rm = TRUE),
                 xlab = "Category", ylab = "Average measure (Bn - Di)",
                 main = "Observed vs expected average measure",
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label)
  .grid(style, h = TRUE)
  graphics::lines(d$Category, d$Sample_Expect, type = "b", pch = 22,
                  bg = style$col_item, col = style$col_border, lty = 2, lwd = style$lwd)
  graphics::legend("topleft", c("Observed", "Expected"), pt.bg = c(style$col_person, style$col_item),
                   pch = c(21, 22), bty = "n", cex = style$cex_label * 0.85)

  tau <- d$Andrich_Threshold[-1]
  graphics::barplot(tau, names.arg = paste0("t", seq_along(tau)),
                    col = ifelse(c(TRUE, diff(tau) > 0), style$col_pos, style$col_neg),
                    border = style$col_border, ylab = "Andrich threshold (logits)",
                    main = "Andrich thresholds (disordered = highlighted)",
                    cex.lab = style$cex_label + 0.1, cex.main = style$cex_label)
  graphics::abline(h = 0, col = style$col_line)
  invisible(NULL)
}
