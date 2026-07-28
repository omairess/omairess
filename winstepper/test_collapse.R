# =============================================================================
# test_collapse.R  --  checks for suggest_collapse() and threshold_advance_min()
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

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
