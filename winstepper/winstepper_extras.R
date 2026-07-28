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
# Suggested rescore: propose a NEWSCORE= collapse for a misbehaving scale
# ---------------------------------------------------------------------------

#' Format a numeric vector the way the CODES= / NEWSCORE= boxes expect.
fmt_codes <- function(v) paste(v, collapse = ",")

#' Suggest a category collapse for one item group.
#'
#' Purely a suggestion: it returns the CODES= and NEWSCORE= strings and a
#' plain-language reason for each merge. Nothing is applied and nothing is
#' re-estimated -- the user confirms or edits it.
#'
#' Two passes over the group's rows of `category_table()`, each stopping before
#' the scale would drop below `min_cat` categories:
#'
#'  1. Sparse categories: while any block has fewer than `min_count`
#'     observations, merge the sparsest into whichever neighbour has the smaller
#'     count (at either end, merge inward). Counts are exactly additive under
#'     merging, so this pass needs no re-estimation and is exact.
#'  2. Narrow / disordered thresholds: a block whose entry and exit thresholds
#'     advance by less than `threshold_advance_min()` allows -- including the
#'     disordered case, where the advance is negative -- is merged into its
#'     smaller neighbour.
#'
#' Pass 2 necessarily reasons from the CURRENT calibration: thresholds are
#' re-estimated once categories change, so treat the result as a starting point
#' and re-read the diagnostics after applying it.
#'
#' @param ct output of `category_table(fit)`
#' @param group item group to work on (defaults to the first)
#' @param min_count minimum observations per category (Linacre 1999: 10)
#' @param min_cat never propose fewer than this many categories
#' @return list(codes, newscore, codes_txt, newscore_txt, reasons, n_cat, changed)
suggest_collapse <- function(ct, group = NULL, min_count = 10, min_cat = 3) {
  if (is.null(ct) || !nrow(ct)) return(NULL)
  if (is.null(group) || !group %in% ct$Group) group <- ct$Group[1]
  d <- ct[ct$Group == group, , drop = FALSE]
  d <- d[order(d$Category), , drop = FALSE]
  cats <- d$Category; counts <- d$Count; tau <- d$Andrich_Threshold
  n <- nrow(d)
  blocks <- as.list(seq_len(n))          # each block = indices of merged categories
  reasons <- character(0)
  lab <- function(b) paste(cats[b], collapse = "+")
  bcount <- function(b) sum(counts[b], na.rm = TRUE)

  ## --- pass 1: sparse categories (exact; counts add under merging) ---------
  repeat {
    if (length(blocks) <= min_cat) break
    bc <- vapply(blocks, bcount, numeric(1))
    if (all(bc >= min_count)) break
    j <- which.min(bc)
    k <- if (j == 1) 2L else if (j == length(blocks)) length(blocks) - 1L else
      if (bc[j - 1] <= bc[j + 1]) j - 1L else j + 1L
    reasons <- c(reasons, sprintf(
      "Category %s has %d observation%s (fewer than %d); merged with category %s.",
      lab(blocks[[j]]), bc[j], if (bc[j] == 1) "" else "s", min_count, lab(blocks[[k]])))
    lo <- min(j, k); hi <- max(j, k)
    blocks[[lo]] <- sort(c(blocks[[lo]], blocks[[hi]]))
    blocks[[hi]] <- NULL
  }

  ## --- pass 2: narrow or disordered thresholds (from the current solution) --
  repeat {
    if (length(blocks) <= max(min_cat, 3)) break
    # entry threshold of each block after the first = tau of its first category
    bnd <- vapply(blocks[-1], function(b) tau[b[1]], numeric(1))
    if (!all(is.finite(bnd))) break                 # cannot judge; stop quietly
    adv <- diff(bnd)                                # width of blocks 2..(k-1)
    req <- threshold_advance_min(length(bnd))       # category-count dependent
    if (!length(adv) || !length(req)) break
    short <- adv - req
    bad <- which(short < 0)
    if (!length(bad)) break
    i <- bad[which.min(short[bad])]                 # worst offender
    j <- i + 1L                                     # the too-narrow block
    bc <- vapply(blocks, bcount, numeric(1))
    k <- if (j == 1) 2L else if (j == length(blocks)) length(blocks) - 1L else
      if (bc[j - 1] <= bc[j + 1]) j - 1L else j + 1L
    reasons <- c(reasons, sprintf(
      "Category %s spans only %.2f logits where %.2f is required for %d categories%s; merged with category %s.",
      lab(blocks[[j]]), adv[i], req[i], length(bnd) + 1L,
      if (adv[i] <= 0) " (thresholds disordered)" else "", lab(blocks[[k]])))
    lo <- min(j, k); hi <- max(j, k)
    blocks[[lo]] <- sort(c(blocks[[lo]], blocks[[hi]]))
    blocks[[hi]] <- NULL
  }

  newscore <- integer(n)
  for (i in seq_along(blocks)) newscore[blocks[[i]]] <- i - 1L
  list(codes = cats, newscore = newscore,
       codes_txt = fmt_codes(cats), newscore_txt = fmt_codes(newscore),
       reasons = reasons, n_cat = length(blocks),
       changed = length(blocks) < n, group = group)
}

# ---------------------------------------------------------------------------
# Table 3.2 : NA-robust category structure
# ---------------------------------------------------------------------------

#' Category structure table, tolerant of NA measures and unbounded zones.
#'
#' DELIBERATELY OVERRIDES `category_table()` from rasch_engine.R, which has two
#' `if (<NA>)` faults that abort with "missing value where TRUE/FALSE needed":
#'
#'  * line 683, `if (sum(w) > 0)` with `w <- P[[k+1]] * selg`. In R
#'    `NA * FALSE` is `NA`, not 0, so multiplying by the logical mask does NOT
#'    zero out the NA probabilities belonging to persons with no estimable
#'    measure -- `sum(w)` becomes NA. Same for `sum(bd * w)`.
#'  * line 705, `if (sum(inzone) > 0)`. The zone bounds come from
#'    `.score_to_measure()`, which returns NA when the expected score cannot be
#'    bracketed in [-25, 25]. `TRUE & NA` is NA, so `inzone` is poisoned. Long
#'    rating scales (0-10) hit this because the thresholds spread far enough
#'    that the half-score points fall outside the bracket, especially when some
#'    categories are sparsely observed.
#'
#' Here the mask is applied by assignment rather than multiplication, and a
#' category whose zone is not computable reports NA coherence instead of
#' failing. The root-search bracket is also widened from +/-25 to +/-60 logits
#' (see .sm_wide): the expected score is strictly monotone in the measure, so
#' the root is unique and widening can only turn an NA into the value that was
#' always there -- it cannot change a number the engine already produced.
#' Everything else matches the engine.
#' Wider-bracket wrappers. The expected score / cumulative probability are
#' strictly monotone in the measure, so the root is unique; widening the search
#' interval only recovers roots that lie outside the engine's +/-25 default,
#' which happens on long rating scales with widely spread thresholds.
.sm_wide <- function(target, tau) .score_to_measure(target, tau, lo = -60, hi = 60)
.th_wide <- function(k, tau)      .thurstone(k, tau, lo = -60, hi = 60)

#' Rasch-Thurstone item thresholds, using the wider bracket (overrides the
#' engine's threshold_data(), which is what the Wright map's "thresholds" view
#' plots).
threshold_data <- function(fit) {
  out <- list()
  for (i in seq_along(fit$delta)) {
    g <- fit$groups[i]; tau <- fit$tau[[g]]; m <- fit$max_cat[[g]]
    if (is.na(fit$delta[i])) next
    th <- vapply(1:m, function(k) .th_wide(k, tau), numeric(1))
    out[[length(out) + 1]] <- data.frame(
      item = fit$item_id[i], threshold = 1:m,
      measure = fit$delta[i] + th,
      label = paste0(fit$item_id[i], ".", 1:m), stringsAsFactors = FALSE)
  }
  if (!length(out)) return(data.frame(item = character(0), threshold = integer(0),
                                      measure = numeric(0), label = character(0)))
  do.call(rbind, out)
}

category_table <- function(fit, group = NULL, cat_extreme = 0.25) {
  gl <- if (is.null(group)) names(fit$tau) else group
  res <- list()
  for (g in gl) {
    cols <- which(fit$groups == g)
    m <- fit$max_cat[[g]]
    tau <- fit$tau[[g]]
    sel <- fit$mask
    sel[!fit$keep_p, ] <- FALSE
    sel[, !fit$keep_i] <- FALSE
    selg <- sel[, cols, drop = FALSE]
    xg   <- fit$X[, cols, drop = FALSE]
    bd   <- outer(fit$theta, fit$delta[cols], "-")   # Bn - Di
    Eg   <- fit$E[, cols, drop = FALSE]
    Wg   <- fit$W[, cols, drop = FALSE]
    P    <- .probs_group(fit$theta, fit$delta[cols], tau)
    # bd is NA exactly where theta/delta are non-estimable; those cells are
    # always outside selg, so zeroing them cannot bias a weighted mean.
    bdz <- bd; bdz[is.na(bdz)] <- 0

    cnt <- oa <- ea <- outms <- infms <- numeric(m + 1)
    for (k in 0:m) {
      inc <- selg & !is.na(xg) & xg == k
      cnt[k + 1] <- sum(inc)
      oa[k + 1]  <- if (cnt[k + 1] > 0) mean(bd[inc]) else NA_real_
      w <- P[[k + 1]]
      w[is.na(w)] <- 0
      w[!selg] <- 0                       # mask by assignment, not by NA * FALSE
      sw <- sum(w)
      ea[k + 1] <- if (isTRUE(sw > 0)) sum(bdz * w) / sw else NA_real_
      if (cnt[k + 1] > 0) {
        y <- (xg - Eg)[inc]; wv <- Wg[inc]
        outms[k + 1] <- mean(y^2 / wv)
        infms[k + 1] <- sum(y^2) / sum(wv)
      } else { outms[k + 1] <- NA_real_; infms[k + 1] <- NA_real_ }
    }
    andrich <- c(NA_real_, tau)
    se_and  <- c(NA_real_, fit$se_tau[[g]])
    catmeas <- vapply(0:m, function(k) {
      tgt <- if (k == 0) cat_extreme else if (k == m) m - cat_extreme else k
      .sm_wide(tgt, tau)
    }, numeric(1))
    thur <- c(NA_real_, vapply(1:m, function(k) .th_wide(k, tau), numeric(1)))
    zlo <- vapply(0:m, function(k) if (k == 0) -Inf else .sm_wide(k - 0.5, tau), numeric(1))
    zhi <- vapply(0:m, function(k) if (k == m)  Inf else .sm_wide(k + 0.5, tau), numeric(1))

    mc <- cm <- rep(NA_real_, m + 1)
    for (k in 0:m) {
      lo <- zlo[k + 1]; hi <- zhi[k + 1]
      if (is.na(lo) || is.na(hi)) next    # zone not computable -> coherence NA
      inzone <- selg & !is.na(xg) & bd >= lo & bd < hi
      inzone[is.na(inzone)] <- FALSE
      incat  <- selg & !is.na(xg) & xg == k
      both <- sum(inzone & incat)
      mc[k + 1] <- if (sum(inzone) > 0) 100 * both / sum(inzone) else NA_real_
      cm[k + 1] <- if (sum(incat)  > 0) 100 * both / sum(incat)  else NA_real_
    }
    tot <- sum(cnt)
    res[[g]] <- data.frame(
      Group = g, Category = 0:m, Score = 0:m, Count = cnt,
      Percent = 100 * cnt / max(tot, 1),
      Obsvd_Avrge = oa, Sample_Expect = ea,
      Infit_MNSQ = infms, Outfit_MNSQ = outms,
      Andrich_Threshold = andrich, Threshold_SE = se_and,
      Category_Measure = catmeas,
      Thurstone_50pct = thur,
      Zone_Lo = zlo, Zone_Hi = zhi,
      Coherence_M_to_C = mc, Coherence_C_to_M = cm,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, res)
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
  if (!nrow(kf)) {
    graphics::plot.new()
    graphics::text(0.5, 0.58, "No keyform coordinates could be computed.", font = 2)
    graphics::text(0.5, 0.42, paste(
      "Every item's category measures are undefined for this keyform type.",
      "Try the 'most probable / modal' type, which uses the Andrich thresholds",
      "directly, and check the category structure on the Rating scale tab.",
      sep = "\n"), cex = 0.9)
    return(invisible(NULL))
  }
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
