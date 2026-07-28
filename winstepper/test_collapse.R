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
for (r in s$reasons) cat("       - ", r, "\n", sep = "")
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
check("a reason is given for every merge", length(s$reasons) == length(s$codes) - s$n_cat)

cat("\n3. A healthy scale is left alone\n")
h <- suggest_collapse(mk(c(80, 70, 75, 90), c(-1.5, 0, 1.5)), min_count = 10, min_cat = 3)
check("changed flag is FALSE", !isTRUE(h$changed))
check("newscore equals codes", identical(as.integer(h$newscore), as.integer(h$codes)))
check("no reasons", length(h$reasons) == 0)

cat("\n4. Disordered thresholds trigger a merge\n")
d <- suggest_collapse(mk(c(60, 55, 50, 65, 70), c(-1.0, -1.2, 0.5, 1.4)),
                      min_count = 10, min_cat = 3)
check("a merge was proposed", isTRUE(d$changed))
check("reason mentions the threshold problem",
      any(grepl("logits|disordered", d$reasons)),
      sprintf("(reasons: %s)", paste(d$reasons, collapse = " | ")))

cat("\n5. min_cat floor is respected even when everything is sparse\n")
f <- suggest_collapse(mk(c(1, 1, 1, 1, 1, 1), seq(-2, 2, length.out = 5)),
                      min_count = 10, min_cat = 3)
check("stops at min_cat = 3", f$n_cat == 3, sprintf("(got %d)", f$n_cat))
f2 <- suggest_collapse(mk(c(1, 1, 1, 1, 1, 1), seq(-2, 2, length.out = 5)),
                       min_count = 10, min_cat = 2)
check("honours min_cat = 2", f2$n_cat == 2, sprintf("(got %d)", f2$n_cat))

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
