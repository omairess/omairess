# =============================================================================
# winstepper_extras.R  --  WINSTEPPER-only additions on top of the audited
# rasch_engine.R / winsteps_plots.R. Kept separate so the reused engine and
# plot files stay byte-identical to R-Winsteps.
#
# Adds:
#   * threshold_advance_min() / category_diagnostics()  -- category-count-aware
#     minimum threshold advance. DELIBERATELY OVERRIDES the engine's
#     category_diagnostics() (see note below).
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
# Minimum threshold advance (replaces the flat "1.4 logits" rule)
# ---------------------------------------------------------------------------

#' Minimum required advance between adjacent Andrich thresholds.
#'
#' The often-quoted "thresholds must advance by 1.4 logits" is NOT a general
#' constant: 1.4 is the three-category case only. The criterion comes from the
#' binomial boundary condition, under which the threshold locations for m
#' thresholds (m + 1 categories) are
#'
#'     tau_k = ln( k / (m + 1 - k) ),        k = 1 .. m
#'
#' so the minimum required advance between adjacent thresholds is
#'
#'     Delta_k = ln( (k + 1)(m + 1 - k) / (k (m - k)) ),   k = 1 .. m - 1
#'
#' For 3 categories this gives ln(2) - ln(1/2) = 2 ln 2 = 1.386 ~ 1.4; for more
#' categories the required advance is smaller (and Linacre's rule of thumb in
#' RMT 2006, 20:1, p. 1052 is "thresholds must advance by one logit").
#'
#' @param m number of thresholds (= number of categories - 1)
#' @return numeric vector of length m - 1 (empty when m < 2)
#' @references Linacre, J. M. (2006). Rasch Measurement Transactions, 20(1), 1052.
threshold_advance_min <- function(m) {
  if (!is.finite(m) || m < 2) return(numeric(0))
  k <- seq_len(m - 1)
  log((k + 1) * (m + 1 - k) / (k * (m - k)))
}

#' Linacre's (1999, 2002, 2006) rating-scale quality guidelines.
#'
#' NOTE: this deliberately overrides `category_diagnostics()` from
#' rasch_engine.R, which tested every threshold advance against a flat 1.4
#' logits. That is only correct for three categories. Everything else in this
#' function is unchanged from the engine version. The engine file is left
#' byte-identical to R-Winsteps; the fix lives here.
category_diagnostics <- function(ct) {
  out <- list()
  for (g in unique(ct$Group)) {
    d <- ct[ct$Group == g, ]
    m <- nrow(d) - 1                      # number of thresholds
    tau <- d$Andrich_Threshold[-1]
    chk <- function(name, ok, detail) data.frame(Group = g, Guideline = name,
                                                 Status = ifelse(ok, "OK", "REVIEW"),
                                                 Detail = detail, stringsAsFactors = FALSE)
    out[[length(out) + 1]] <- chk(
      "At least 10 observations per category",
      all(d$Count >= 10),
      paste0("min count = ", min(d$Count)))
    out[[length(out) + 1]] <- chk(
      "Regular / uniform category distribution",
      all(d$Percent > 0),
      paste0("percentages: ", paste(sprintf("%.1f", d$Percent), collapse = ", ")))
    oa <- d$Obsvd_Avrge
    out[[length(out) + 1]] <- chk(
      "Observed average measures advance with category",
      all(diff(oa[!is.na(oa)]) > 0),
      paste0("observed averages: ", paste(sprintf("%.2f", oa), collapse = ", ")))
    out[[length(out) + 1]] <- chk(
      "Category outfit mean-squares < 2.0",
      all(d$Outfit_MNSQ < 2, na.rm = TRUE),
      paste0("max outfit = ", sprintf("%.2f", max(d$Outfit_MNSQ, na.rm = TRUE))))
    if (m >= 2) {
      out[[length(out) + 1]] <- chk(
        "Andrich thresholds advance (ordered)",
        all(diff(tau) > 0),
        paste0("thresholds: ", paste(sprintf("%.2f", tau), collapse = ", ")))
      dif <- diff(tau)
      req <- threshold_advance_min(m)
      out[[length(out) + 1]] <- chk(
        sprintf("Threshold advances >= binomial minimum (%d categories: %s logits)",
                m + 1, paste(sprintf("%.2f", req), collapse = ", ")),
        all(dif >= req, na.rm = TRUE),
        paste0("advances: ", paste(sprintf("%.2f", dif), collapse = ", "),
               "  |  required: ", paste(sprintf("%.2f", req), collapse = ", "),
               "  |  shortfall: ",
               paste(sprintf("%+.2f", dif - req), collapse = ", ")))
      out[[length(out) + 1]] <- chk(
        "Threshold advances <= 5.0 logits (no gaps in the variable)",
        all(dif <= 5),
        paste0("advances: ", paste(sprintf("%.2f", dif), collapse = ", ")))
    }
    out[[length(out) + 1]] <- chk(
      "Coherence M->C >= 40% (Linacre 2002)",
      all(d$Coherence_M_to_C >= 40, na.rm = TRUE),
      paste0("M->C: ", paste(sprintf("%.0f%%", d$Coherence_M_to_C), collapse = ", ")))
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Table 3.1 : NA-robust summary blocks
# ---------------------------------------------------------------------------

#' Summary block for Table 3.1, tolerant of non-estimable measures.
#'
#' DELIBERATELY OVERRIDES `.summary_block()` from rasch_engine.R. The engine
#' version does `sep <- if (rmse > 0) ...`, which throws
#' "missing value where TRUE/FALSE needed" as soon as any element of the block
#' is NA. That happens routinely on real data: an extreme person whose observed
#' responses all fall on extreme items never receives a measure (in
#' rasch_jmle() the extreme-person loop runs before the extreme-item loop, so
#' those item difficulties are still NA and the person is skipped), and the
#' "(all)" rows of Table 3.1 include exactly those cases. Long rating scales
#' (e.g. 0-10) hit this often because floor/ceiling items are common.
#'
#' This version summarises the estimable cases and reports how many were
#' dropped, instead of failing. `summary_table()` from the engine picks it up
#' automatically (both live in the global environment; extras is sourced last).
.summary_block <- function(meas, se_model, se_real, infit, outfit, score, count, label) {
  n_all <- length(meas)
  ok <- is.finite(meas) & is.finite(se_model) & is.finite(se_real)
  ok[is.na(ok)] <- FALSE
  meas <- meas[ok]; se_model <- se_model[ok]; se_real <- se_real[ok]
  infit <- infit[ok]; outfit <- outfit[ok]; score <- score[ok]; count <- count[ok]
  n <- length(meas)
  if (n < 2) return(NULL)

  sd_pop <- function(v) stats::sd(v) * sqrt((length(v) - 1) / length(v))
  sdm <- sd_pop(meas)
  rmse_m <- sqrt(mean(se_model^2)); rmse_r <- sqrt(mean(se_real^2))
  tsd_m <- sqrt(max(sdm^2 - rmse_m^2, 0)); tsd_r <- sqrt(max(sdm^2 - rmse_r^2, 0))
  # isTRUE() so an NA that slips through yields NA rather than an error.
  sep_m <- if (isTRUE(rmse_m > 0)) tsd_m / rmse_m else NA_real_
  sep_r <- if (isTRUE(rmse_r > 0)) tsd_r / rmse_r else NA_real_
  rel_m <- if (isTRUE(is.finite(sep_m))) sep_m^2 / (1 + sep_m^2) else NA_real_
  rel_r <- if (isTRUE(is.finite(sep_r))) sep_r^2 / (1 + sep_r^2) else NA_real_
  data.frame(
    Statistic = label, N = n, N_Not_Estimable = n_all - n,
    Mean_Score = mean(score), Mean_Count = mean(count),
    Mean_Measure = mean(meas), SD_Measure = sdm,
    Mean_Model_SE = mean(se_model), Mean_Real_SE = mean(se_real),
    Mean_Infit_MNSQ = mean(infit, na.rm = TRUE),
    Mean_Outfit_MNSQ = mean(outfit, na.rm = TRUE),
    Model_RMSE = rmse_m, Real_RMSE = rmse_r,
    True_SD_Model = tsd_m, True_SD_Real = tsd_r,
    Separation_Model = sep_m, Separation_Real = sep_r,
    Reliability_Model = rel_m, Reliability_Real = rel_r,
    Strata_Model = (4 * sep_m + 1) / 3, Strata_Real = (4 * sep_r + 1) / 3,
    stringsAsFactors = FALSE
  )
}

#' Which persons / items could not be measured, and why. Use this to explain a
#' non-zero `N_Not_Estimable` in Table 3.1.
measure_health <- function(fit) {
  bad_p <- which(!is.finite(fit$theta))
  bad_i <- which(!is.finite(fit$delta))
  n_resp_p <- rowSums(fit$mask)
  n_resp_i <- colSums(fit$mask)
  data.frame(
    Quantity = c("Persons", "  extreme (min/max score)", "  measure not estimable",
                 "Items", "  extreme (min/max score)", "  measure not estimable"),
    N = c(nrow(fit$X), length(fit$extreme_persons), length(bad_p),
          ncol(fit$X), length(fit$extreme_items), length(bad_i)),
    Detail = c("", "", if (length(bad_p))
      paste0(paste(utils::head(fit$person_id[bad_p], 10), collapse = ", "),
             if (length(bad_p) > 10) ", ..." else "",
             sprintf("  (median responses = %g)", stats::median(n_resp_p[bad_p]))) else "",
      "", "", if (length(bad_i))
        paste0(paste(utils::head(fit$item_id[bad_i], 10), collapse = ", "),
               if (length(bad_i) > 10) ", ..." else "",
               sprintf("  (median responses = %g)", stats::median(n_resp_i[bad_i]))) else ""),
    stringsAsFactors = FALSE)
}

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
#'
#' @param show_persons if TRUE, a person-measure histogram is drawn underneath,
#'   on exactly the same logit axis (so the person distribution can be read
#'   against the category keys).
plot_keyform <- function(fit, style = ws_style(), kind = c("expected", "thurstone", "modal"),
                         max_items = 80, show_persons = TRUE, main = NULL) {
  kind <- match.arg(kind)
  kf <- keyform_data(fit, kind)
  if (!nrow(kf)) { graphics::plot.new(); graphics::text(0.5, 0.5, "No non-extreme items to plot."); return(invisible(NULL)) }
  im <- tapply(kf$Item_Measure, kf$Item, function(v) v[1])
  items <- names(im)[order(im)]
  if (length(items) > max_items)
    items <- items[round(seq(1, length(items), length.out = max_items))]
  kf <- kf[kf$Item %in% items, , drop = FALSE]
  yof <- match(kf$Item, items)

  pm <- fit$theta[fit$keep_p & !is.na(fit$theta)]
  # One shared logit axis across both panels (cf. the Wright map convention).
  xr <- range(c(kf$Measure, if (show_persons) pm), na.rm = TRUE)
  xr <- xr + c(-1, 1) * max(diff(xr) / 15, 0.3)

  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::layout(1); graphics::par(op) })
  if (show_persons && length(pm) > 1)
    graphics::layout(matrix(1:2, nrow = 2), heights = c(3.2, 1))
  else show_persons <- FALSE

  ttl <- main %||% sprintf("General keyform - %s",
                           c(expected = "expected score (Table 2.2)",
                             thurstone = "Rasch-Thurstone 50% (Table 2.3)",
                             modal = "most probable / modal (Table 2.1)")[kind])
  graphics::par(mar = c(if (show_persons) 0.6 else 4.5, 8, 3.5, 1))
  graphics::plot(NA, xlim = xr, ylim = c(0.5, length(items) + 0.5), yaxt = "n",
                 xaxt = if (show_persons) "n" else "s",
                 xlab = if (show_persons) "" else "Measure (logits)", ylab = "", main = ttl,
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
  dm <- mean(fit$delta, na.rm = TRUE)
  graphics::abline(v = dm, col = style$col_line, lty = 2, lwd = style$lwd)
  graphics::mtext("item difficulty mean", side = 3, line = -1, at = dm,
                  cex = style$cex_label * 0.7)

  ## --- bottom panel: person distribution on the same axis ------------------
  if (show_persons) {
    graphics::par(mar = c(4.5, 8, 0.6, 1))
    brks <- seq(xr[1], xr[2], length.out = max(style$bins, 5) + 1)
    h <- graphics::hist(pm[pm >= xr[1] & pm <= xr[2]], breaks = brks, plot = FALSE)
    graphics::plot(NA, xlim = xr, ylim = c(0, max(h$counts, 1) * 1.08),
                   xlab = "Measure (logits)", ylab = "Persons", yaxs = "i",
                   cex.lab = style$cex_label + 0.1, las = 1,
                   cex.axis = style$cex_label * 0.9)
    if (isTRUE(style$show_grid))
      graphics::abline(v = pretty(xr, 10), col = style$col_ref, lty = 3)
    graphics::rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
                   col = style$col_person, border = style$col_border)
    if (isTRUE(style$show_msT)) {
      mm <- mean(pm); ss <- stats::sd(pm)
      graphics::abline(v = mm, col = style$col_line, lwd = style$lwd)
      graphics::abline(v = c(mm - ss, mm + ss), col = style$col_line, lty = 2)
      graphics::mtext(c("S", "M", "S"), side = 1, line = -1.1,
                      at = c(mm - ss, mm, mm + ss), cex = style$cex_label * 0.75)
    }
    graphics::mtext(sprintf("n = %d non-extreme persons", length(pm)), side = 3,
                    line = -1.2, adj = 0.99, cex = style$cex_label * 0.7)
  }
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
