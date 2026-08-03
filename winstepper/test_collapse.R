# =============================================================================
# test_collapse.R  --  checks for suggest_collapse(), threshold_advance_min()
#                      and the person-item barchart (pi_threshold_data/plot_pi_map)
#
# Run from inside the winstepper folder:   Rscript test_collapse.R
# No Shiny, no network, no API. Pure R.
# =============================================================================

source("rasch_engine.R")
source("winsteps_plots.R")
source("winstepper_extras.R")

ok <- 0L; bad <- 0L
check <- function(label, cond, detail = "") {
  if (isTRUE(cond)) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { bad <<- bad + 1L; cat(sprintf("  FAIL  %s %s\n", label, detail)) }
}

# Build a minimal category_table()-shaped data.frame.
mk <- function(counts, tau) {
  m <- length(counts) - 1
  data.frame(Group = "R1", Category = 0:m, Score = 0:m, Count = counts,
             Andrich_Threshold = c(NA_real_, tau), stringsAsFactors = FALSE)
}

cat("\n1. threshold_advance_min() vs Linacre (RMT 2006, 20:1, p.1052)\n")
ref <- list("3" = 1.39, "4" = c(1.10, 1.10), "5" = c(0.98, 0.81, 0.98),
            "6" = c(0.92, 0.69, 0.69, 0.92),
            "7" = c(0.88, 0.63, 0.58, 0.63, 0.88))
for (nc in names(ref)) {
  got <- round(threshold_advance_min(as.integer(nc) - 1), 2)
  check(sprintf("%s categories -> %s", nc, paste(got, collapse = ", ")),
        isTRUE(all.equal(got, ref[[nc]], tolerance = 1e-8)),
        sprintf("(expected %s)", paste(ref[[nc]], collapse = ", ")))
}
check("dichotomy yields no advance requirement", length(threshold_advance_min(1)) == 0)

cat("\n2. Sparse categories are merged (11-category scale)\n")
counts <- c(50, 3, 2, 40, 60, 55, 4, 30, 20, 2, 1)
s <- suggest_collapse(mk(counts, seq(-3, 3, length.out = 10)), min_count = 10, min_cat = 3)
cat(sprintf("     CODES:    %s\n     NEWSCORE: %s\n     %d -> %d categories\n",
            s$codes_txt, s$newscore_txt, length(s$codes), s$n_cat))
check("newscore has one entry per code", length(s$newscore) == length(s$codes))
check("newscore starts at 0", s$newscore[1] == 0)
check("newscore is monotone non-decreasing", all(diff(s$newscore) >= 0))
check("newscore steps by at most 1", all(diff(s$newscore) <= 1))
check("uses every score 0..n_cat-1 with no gaps",
      identical(sort(unique(s$newscore)), as.integer(seq_len(s$n_cat) - 1)))
check("collapsed categories all reach the minimum count",
      all(tapply(counts, s$newscore, sum) >= 10),
      sprintf("(got %s)", paste(tapply(counts, s$newscore, sum), collapse = ", ")))
check("never drops below min_cat", s$n_cat >= 3)

cat("\n3. A healthy scale is left alone\n")
h <- suggest_collapse(mk(c(80, 70, 75, 90), c(-1.5, 0, 1.5)),
                      min_count = 10, min_sep = 0.5, min_cat = 3)
check("changed flag is FALSE", !isTRUE(h$changed))
check("newscore equals codes", identical(as.integer(h$newscore), as.integer(h$codes)))

cat("\n4. REGRESSION: a compressed 10-category scale must NOT collapse to a\n")
cat("   two-poles-and-one-huge-middle solution (the 0 | 1-8 | 9 failure).\n")
# Real marginal distribution from NSGGMgameCONT.xlsx, 17 items x 1025 persons,
# 1-10 shifted to 0-9. No category is sparse; the problem is compression.
real <- c(3802, 1886, 1959, 1413, 1565, 1489, 1675, 1601, 988, 1047)
r <- suggest_collapse(mk(real, seq(-0.55, 0.55, length.out = 9)),
                      min_count = 10, min_sep = 1.0, min_cat = 3)
cat(sprintf("     NEWSCORE: %s   (%d -> %d categories)\n",
            r$newscore_txt, length(r$codes), r$n_cat))
cat(sprintf("     bands:    %s\n", paste(r$bands$Categories, collapse = " | ")))
check("keeps more than 3 categories", r$n_cat > 3, sprintf("(got %d)", r$n_cat))
widest <- max(r$bands$Pct)
check("no band swallows more than half the responses", widest <= 50,
      sprintf("(widest band = %.1f%%)", widest))
check("the floor category stays its own level",
      identical(r$bands$Categories[1], "0"),
      sprintf("(first band = %s)", r$bands$Categories[1]))
check("retained thresholds meet the requested separation",
      length(r$sep) == 0 || min(r$sep) >= 1.0,
      sprintf("(min separation = %.2f)", if (length(r$sep)) min(r$sep) else Inf))
check("an alternative is offered", !is.null(r$alt_newscore_txt))

cat("\n5. Raising the separation requirement yields a coarser scale\n")
coarse <- suggest_collapse(mk(real, seq(-0.55, 0.55, length.out = 9)),
                           min_count = 10, min_sep = 1.4, min_cat = 3)
cat(sprintf("     min_sep 1.4 -> %s  (%d categories)\n",
            coarse$newscore_txt, coarse$n_cat))
check("coarser than the 1.0-logit solution", coarse$n_cat <= r$n_cat)
check("still above the floor", coarse$n_cat >= 3)

cat("\n6. Degenerate: no collapse can reach min_count (group too small)\n")
f <- suggest_collapse(mk(c(1, 1, 1, 1, 1, 1), seq(-2, 2, length.out = 5)),
                      min_count = 10, min_cat = 3)
check("still returns a suggestion rather than nothing", !is.null(f))
check("does not fall below min_cat = 3", !is.null(f) && f$n_cat >= 3,
      if (is.null(f)) "(returned NULL)" else sprintf("(got %d)", f$n_cat))
check("explains that the count floor is unreachable",
      !is.null(f) && any(grepl("in total", f$notes)),
      if (is.null(f)) "" else sprintf("(notes: %s)", paste(f$notes, collapse = " | ")))

cat("\n7. Person-item barchart threshold data (all four WINSTEPS types)\n")
set.seed(1)
dat  <- demo_data(n = 250, k = 8, m = 3)
prep <- rasch_prep(as.matrix(dat[, grep("^Item", names(dat))]), recode = "shift")
fit  <- rasch_jmle(prep, maxit = 150, conv = 1e-3, rconv = 1e-2)
nI   <- sum(fit$keep_i)
for (w in c("andrich", "thurstone", "halfpoint", "fullpoint")) {
  td <- pi_threshold_data(fit, w)
  check(sprintf("%-10s returns m thresholds for every non-extreme item", w),
        nrow(td) == nI * 3 && all(is.finite(td$Measure)),
        sprintf("(rows = %d, expected %d; %d non-finite)",
                nrow(td), nI * 3, sum(!is.finite(td$Measure))))
}
td <- pi_threshold_data(fit, "thurstone")
check("thresholds ascend within each item",
      all(tapply(td$Measure, td$Item, function(v) all(diff(v) > 0))))
check("plot_pi_map() runs and produces a plot", {
  pf <- tempfile(fileext = ".png")
  grDevices::png(pf, width = 1400, height = 900, res = 110)
  r2 <- try(plot_pi_map(fit, ws_style(), person_class = dat$Sex,
                        what = "andrich", sort_by = "measure"), silent = TRUE)
  grDevices::dev.off()
  !inherits(r2, "try-error") && file.exists(pf) && file.info(pf)$size > 5000
})

cat("\n8. Labels always fit the space reserved for them\n")
# The figures size their margins from .str_in(), which estimates text width from
# par("cin") because strwidth() cannot be called before par(mar=) is set. If that
# estimate ever under-shoots the real rendered width, labels get clipped - which
# is exactly the bug the helpers were added to fix. So: measure both.
labs <- c("a", "Item01", "GamingObsession_prob", "WWWWWWWWWWWW", "@@@@@@@@@@",
          "mmmmmmmmmmmm", "OOOOOOOOOO", strrep("i", 30),
          "Mixed Case With Spaces 123", "Item_04 (reverse-scored)")
lo <- Inf; hi <- 0
pf <- tempfile(fileext = ".png")
grDevices::png(pf, width = 1200, height = 800, res = 110)
graphics::plot.new()
for (cx in c(0.5, 0.75, 1, 1.5, 2.5)) for (l in labs) {
  r <- .str_in(l, cx) / graphics::strwidth(l, units = "inches", cex = cx)
  lo <- min(lo, r); hi <- max(hi, r)
}
lh <- graphics::par("csi"); fw <- graphics::par("fin")[1L]
fl_ok <- .fit_labels(labs, 1)
enough <- (fl_ok$lines - 1.1) * lh >= max(graphics::strwidth(fl_ok$labels, "inches", cex = 1))
long <- strrep("W", 200)
fl_tr <- .fit_labels(long, 1, side = 2L, max_frac = 0.30)
capped <- (fl_tr$lines - 1.1) * lh <= 0.30 * fw + 1e-9
trunc_fits <- graphics::strwidth(fl_tr$labels, "inches", cex = 1) <= 0.30 * fw
grDevices::dev.off()

check("estimated width never under-shoots the rendered width", lo >= 1,
      sprintf("(worst estimate/actual ratio = %.3f)", lo))
check("and is not wastefully generous", hi < 2.2,
      sprintf("(largest estimate/actual ratio = %.3f)", hi))
check(".fit_labels() reserves enough margin for the longest label", enough)
check("a pathological label is truncated with an ellipsis",
      nchar(fl_tr$labels) < 200L && grepl("…$", fl_tr$labels) &&
        isTRUE(fl_tr$truncated))
check("the truncated label really fits the cap", isTRUE(trunc_fits))
check("the margin is capped at max_frac of the figure", isTRUE(capped),
      sprintf("(%.3f in vs cap %.3f in)", (fl_tr$lines - 1.1) * lh, 0.30 * fw))
check("margin grows with the label size", {
  a <- .fit_labels("GamingObsession_prob", 0.6)$lines
  b <- .fit_labels("GamingObsession_prob", 1.8)$lines
  b > a * 1.5
})

cat("\n9. REGRESSION: long item names are not clipped in the keyform\n")
# The reported failure: plot_keyform() hard-coded an 8-line left margin, so
# 'GamingObsession_prob' rendered as 'gObsession_prob'.
long_names <- c("GamingObsession_prob", "GamingWithdrawal_prob",
                "GamingToleranceIncrease_prob", "Mood", "Conflict_prob",
                "SalienceOfGamingActivity_prob", "RelapseAfterAbstinence_prob",
                "ProblemsAtWorkOrSchool_prob")
fit2 <- fit
colnames(fit2$X) <- fit2$item_id <- long_names
check("plot_keyform() runs with long names at a large label size", {
  kf <- tempfile(fileext = ".png")
  grDevices::png(kf, width = 1400, height = 900, res = 110)
  r3 <- try(plot_keyform(fit2, ws_style(cex_label = 2.0), kind = "expected"),
            silent = TRUE)
  grDevices::dev.off()
  !inherits(r3, "try-error") && file.exists(kf) && file.info(kf)$size > 5000
})
check("the left margin scales with the names and the label size", {
  kf <- tempfile(fileext = ".png")
  grDevices::png(kf, width = 1400, height = 900, res = 110)
  graphics::plot.new()
  small <- .fit_labels(long_names, 0.85 * 0.8, side = 2L)$lines
  big   <- .fit_labels(long_names, 2.00 * 0.8, side = 2L)$lines
  grDevices::dev.off()
  small > 8 && big > small          # 8 = the old hard-coded value that clipped
})

cat("\n10. Every figure survives an extreme label size\n")
for (nm in c("plot_wright", "plot_pathway", "plot_dif_contrast", "plot_pi_map")) {
  check(sprintf("%s() at cex_label = 2.5", nm), {
    ff <- tempfile(fileext = ".png")
    grDevices::png(ff, width = 1400, height = 900, res = 110)
    st <- ws_style(cex_label = 2.5)
    rr <- try(switch(nm,
      plot_wright       = plot_wright(fit2, st),
      plot_pathway      = plot_pathway(fit2, st),
      plot_dif_contrast = plot_dif_contrast(dif_analysis(fit2, dat$Sex), st),
      plot_pi_map       = plot_pi_map(fit2, st, what = "andrich")), silent = TRUE)
    grDevices::dev.off()
    !inherits(rr, "try-error") && file.exists(ff) && file.info(ff)$size > 5000
  })
}

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
