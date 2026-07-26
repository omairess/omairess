# =============================================================================
# winstepper_extras.R  --  WINSTEPPER-only additions on top of the audited
# rasch_engine.R / winsteps_plots.R. Kept separate so the reused engine and
# plot files stay byte-identical to R-Winsteps.
#
# Adds:
#   * keyform_data() / plot_keyform()  -- WINSTEPS Table 2.2 general keyform
#     (expected-score), with Rasch-Thurstone (2.3) and modal (2.1) variants.
#   * dgf_analysis() / plot_dgf()      -- WINSTEPS Table 33 Differential Group
#     Functioning: item-class x person-class interaction, sitting alongside DIF.
#
# Source AFTER rasch_engine.R and winsteps_plots.R (it uses their internals:
# category_table(), .probs_group()-style maths, .ets_class(), ws_style(),
# .pal(), .grid(), %||%).
# =============================================================================

# ---------------------------------------------------------------------------
# Table 2.2 : General keyform
# ---------------------------------------------------------------------------

#' Long data frame of the category "markers" that make up a keyform.
#'
#' @param kind "expected" (Table 2.2, measure at which the expected score on the
#'   item equals each category), "thurstone" (Table 2.3, Rasch-Thurstone 50%
#'   cumulative thresholds) or "modal" (Table 2.1, Andrich thresholds = points
#'   where adjacent categories are equally probable).
#' @return data.frame(Entry, Item, Item_Measure, Category, Measure, Label)
keyform_data <- function(fit, kind = c("expected", "thurstone", "modal"),
                         include_extreme = FALSE) {
  kind <- match.arg(kind)
  ct <- category_table(fit)
  keep <- if (include_extreme) rep(TRUE, length(fit$delta)) else fit$keep_i
  rows <- list()
  for (i in which(keep & !is.na(fit$delta))) {
    g <- fit$groups[i]
    d <- ct[ct$Group == g, , drop = FALSE]
    de <- fit$delta[i]
    meas <- switch(kind,
                   expected  = de + d$Category_Measure,
                   thurstone = de + d$Thurstone_50pct,
                   modal     = de + d$Andrich_Threshold)
    ok <- is.finite(meas)
    if (!any(ok)) next
    rows[[length(rows) + 1]] <- data.frame(
      Entry = i, Item = fit$item_id[i], Item_Measure = de,
      Category = d$Category[ok], Measure = meas[ok],
      Label = as.character(d$Category[ok]), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  attr(out, "kind") <- kind
  out
}

#' WINSTEPS Table 2.2 general keyform. Items are rows ordered by measure; the
#' category numbers are placed along the shared logit axis where they "key".
plot_keyform <- function(fit, style = ws_style(), kind = c("expected", "thurstone", "modal"),
                         max_items = 80, main = NULL) {
  kind <- match.arg(kind)
  kf <- keyform_data(fit, kind)
  if (!nrow(kf)) { graphics::plot.new(); graphics::text(0.5, 0.5, "No non-extreme items to plot."); return(invisible(NULL)) }
  im <- tapply(kf$Item_Measure, kf$Item, function(v) v[1])
  items <- names(im)[order(im)]
  if (length(items) > max_items)
    items <- items[round(seq(1, length(items), length.out = max_items))]
  kf <- kf[kf$Item %in% items, , drop = FALSE]
  yof <- match(kf$Item, items)
  xr <- range(kf$Measure, na.rm = TRUE); xr <- xr + c(-1, 1) * max(diff(xr) / 15, 0.3)

  op <- graphics::par(mar = c(4.5, 8, 3.5, 1)); on.exit(graphics::par(op))
  ttl <- main %||% sprintf("General keyform - %s",
                           c(expected = "expected score (Table 2.2)",
                             thurstone = "Rasch-Thurstone 50% (Table 2.3)",
                             modal = "most probable / modal (Table 2.1)")[kind])
  graphics::plot(NA, xlim = xr, ylim = c(0.5, length(items) + 0.5), yaxt = "n",
                 xlab = "Measure (logits)", ylab = "", main = ttl,
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  if (isTRUE(style$show_grid))
    graphics::abline(v = pretty(xr, 10), col = style$col_ref, lty = 3)
  graphics::axis(2, at = seq_along(items), labels = items, las = 1,
                 cex.axis = style$cex_label * 0.8)
  for (j in seq_along(items)) {
    seg <- kf[yof == j, , drop = FALSE]
    graphics::segments(min(seg$Measure), j, max(seg$Measure), j,
                       col = style$col_ref, lwd = style$lwd)
    graphics::points(seg$Measure, rep(j, nrow(seg)), pch = 21, bg = style$col_item,
                     col = style$col_border, cex = style$cex_point * 0.85)
    graphics::text(seg$Measure, rep(j, nrow(seg)), seg$Label, pos = 3,
                   cex = style$cex_label * 0.75, offset = 0.3)
  }
  graphics::abline(v = mean(fit$delta, na.rm = TRUE), col = style$col_line, lty = 2, lwd = style$lwd)
  graphics::mtext("item difficulty mean", side = 3, line = -1,
                  at = mean(fit$delta, na.rm = TRUE), cex = style$cex_label * 0.7)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Table 33 : Differential Group Functioning (DGF)
# ---------------------------------------------------------------------------

#' Elementwise expected score and model variance for paired (theta, delta).
#' Mirrors the engine's .probs_group() maths but for vectors of paired values.
.escore_pairs <- function(t, del, tau) {
  m <- length(tau)
  ct <- c(0, cumsum(tau))
  td <- t - del
  psi <- vapply(0:m, function(k) k * td - ct[k + 1], numeric(length(td)))
  psi <- matrix(psi, ncol = m + 1)
  mx <- apply(psi, 1, max)
  ep <- exp(psi - mx)
  s <- rowSums(ep)
  P <- ep / s
  E <- as.vector(P %*% (0:m))
  V <- as.vector(P %*% ((0:m)^2)) - E^2
  list(E = E, W = V)
}

#' Uniform difficulty shift for one item-class x person-class cell, with person
#' measures, baseline item difficulties and thresholds all anchored.
.dgf_cell <- function(x, t, del, grp, tau_list, tol = 1e-8, maxit = 200) {
  d <- 0
  for (it in seq_len(maxit)) {
    E <- W <- numeric(length(x))
    for (g in unique(grp)) {
      sel <- grp == g
      es <- .escore_pairs(t[sel], del[sel] + d, tau_list[[g]])
      E[sel] <- es$E; W[sel] <- es$W
    }
    step <- (sum(E) - sum(x)) / max(sum(W), 1e-8)   # Newton: d <- d - (sumx-sumE)/sumW
    step <- max(min(step, 1), -1)
    d <- max(min(d + step, 20), -20)                # guard degenerate (near-extreme) cells
    if (abs(step) < tol) break
  }
  E <- W <- numeric(length(x))
  for (g in unique(grp)) {
    sel <- grp == g
    es <- .escore_pairs(t[sel], del[sel] + d, tau_list[[g]])
    E[sel] <- es$E; W[sel] <- es$W
  }
  list(shift = d, se = 1 / sqrt(max(sum(W), 1e-8)), n = length(x),
       obsexp = mean(x - E))
}

#' WINSTEPS Table 33 DGF: the interaction of an item classification with a
#' person classification. For each (item-class, person-class) cell a single
#' uniform difficulty shift is estimated on top of the baseline calibration;
#' the DGF contrast is the difference of shifts across person classes.
#'
#' @param person_class factor length nrow(fit$X): the person grouping (DIF=)
#' @param item_group   vector length ncol(fit$X): the item classification
#'   (DIF@= / item groups). Defaults to the model's rating-scale groups.
#' @param min_n minimum observations per cell.
dgf_analysis <- function(fit, person_class, item_group = fit$groups, min_n = 5) {
  person_class <- as.factor(person_class)
  if (length(person_class) != nrow(fit$X))
    stop("person_class must have one entry per person.")
  ig <- as.character(item_group)
  if (length(ig) != ncol(fit$X))
    stop("item_group must have one entry per item.")
  plv <- levels(droplevels(person_class[fit$keep_p]))
  if (length(plv) < 2)
    stop("DGF needs at least two person classes with non-extreme persons.")
  iglv <- sort(unique(ig[fit$keep_i]))

  shift <- se <- nobs <- oe <-
    matrix(NA_real_, length(iglv), length(plv), dimnames = list(iglv, plv))
  for (cc in iglv) {
    icols <- which(ig == cc & fit$keep_i)
    if (!length(icols)) next
    for (l in plv) {
      prows <- which(fit$keep_p & !is.na(person_class) & person_class == l)
      if (!length(prows)) next
      sub <- fit$mask[prows, icols, drop = FALSE]
      idx <- which(sub, arr.ind = TRUE)
      nobs[cc, l] <- nrow(idx)
      if (nrow(idx) < min_n) next
      pr <- prows[idx[, 1]]; ic <- icols[idx[, 2]]
      cell <- .dgf_cell(x = fit$X[cbind(pr, ic)], t = fit$theta[pr],
                        del = fit$delta[ic], grp = fit$groups[ic],
                        tau_list = fit$tau)
      shift[cc, l] <- cell$shift; se[cc, l] <- cell$se; oe[cc, l] <- cell$obsexp
    }
  }

  pairs <- utils::combn(plv, 2, simplify = FALSE)
  rows <- list()
  for (pp in pairs) {
    A <- pp[1]; B <- pp[2]
    contrast <- shift[, A] - shift[, B]
    jse <- sqrt(se[, A]^2 + se[, B]^2)
    tval <- contrast / jse
    df <- pmax(nobs[, A] + nobs[, B] - 2, 1)
    pval <- 2 * stats::pt(-abs(tval), df = df)
    rows[[length(rows) + 1]] <- data.frame(
      Item_Class = iglv,
      Class_A = A, ObsExp_A = oe[, A], Shift_A = shift[, A], SE_A = se[, A], N_A = nobs[, A],
      Class_B = B, ObsExp_B = oe[, B], Shift_B = shift[, B], SE_B = se[, B], N_B = nobs[, B],
      DGF_Contrast = contrast, Joint_SE = jse, t = tval, df = df, p = pval,
      ETS_Class = .ets_class(contrast, pval),
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  attr(res, "shift") <- shift
  attr(res, "se") <- se
  rownames(res) <- NULL
  res
}

#' DGF plot: item-class difficulty shift by person class (Table 33 companion).
plot_dgf <- function(dgf, style = ws_style(), errbars = TRUE) {
  sh <- attr(dgf, "shift"); se <- attr(dgf, "se")
  cls <- colnames(sh); ic <- rownames(sh)
  cols <- .pal(max(length(cls), 3), style$palette)[seq_along(cls)]
  op <- graphics::par(mar = c(7, 4.5, 3.5, 1)); on.exit(graphics::par(op))
  yl <- range(c(sh - 2 * se, sh + 2 * se), na.rm = TRUE)
  if (!all(is.finite(yl))) yl <- c(-1, 1)
  graphics::plot(NA, xlim = c(0.5, nrow(sh) + 0.5), ylim = yl, xaxt = "n",
                 xlab = "", ylab = "Item-class difficulty shift (logits)",
                 main = "DGF: item class x person class (Table 33)",
                 cex.lab = style$cex_label + 0.1, cex.main = style$cex_label + 0.2)
  if (isTRUE(style$show_grid))
    graphics::abline(h = pretty(yl, 8), col = style$col_ref, lty = 3)
  graphics::axis(1, at = seq_len(nrow(sh)), labels = ic, las = 2,
                 cex.axis = style$cex_label * 0.85)
  for (j in seq_along(cls)) {
    x <- seq_len(nrow(sh)) + (j - (length(cls) + 1) / 2) * 0.14
    if (errbars)
      graphics::segments(x, sh[, j] - 1.96 * se[, j], x, sh[, j] + 1.96 * se[, j],
                         col = style$col_border)
    graphics::points(x, sh[, j], pch = 21, bg = cols[j], col = style$col_border,
                     cex = style$cex_point)
    graphics::lines(x, sh[, j], col = cols[j], lwd = style$lwd, lty = 3)
  }
  graphics::abline(h = 0, col = style$col_line, lwd = style$lwd)
  graphics::legend("topleft", cls, pt.bg = cols, pch = 21, bty = "n",
                   cex = style$cex_label * 0.85, horiz = TRUE)
  invisible(NULL)
}
