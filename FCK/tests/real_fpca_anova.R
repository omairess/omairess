# ==============================================================================
# tests/real_fpca_anova.R — the new per-component ANOVA, on the real data
#
# Runs the smoothing the app runs, then fPCA, then the component ANOVA, using
# the shipped helper module. This is the feature demonstrated end to end rather
# than unit-tested in isolation.
#
#   Rscript tests/real_fpca_anova.R <xlsx>
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressWarnings(suppressMessages({library(readxl); library(fda)}))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))
source(file.path(app_dir, "server/09_helpers_pcanova.R"))

xlsx <- commandArgs(trailingOnly = TRUE)[1]
d <- suppressWarnings(suppressMessages(read_excel(xlsx, sheet = "slaperigheid", .name_repair = "unique")))
a <- suppressWarnings(suppressMessages(read_excel(xlsx, sheet = "Algemene variabelen", .name_repair = "unique")))
tc <- c("8H","9H","10H","11H","12H","14H","16H","18H","20H","21H","22H","23H","0H","2H","4H","6H")
Y <- as.matrix(d[, tc]); storage.mode(Y) <- "double"
tlin <- c(8,9,10,11,12,14,16,18,20,21,22,23,24,26,28,30)
grp  <- as.character(a$AGEcategory)          # row-position join, as the app does
ids  <- as.character(d$ID)

# drop rows with too few observations, as the import does
keep <- rowSums(!is.na(Y)) >= 8
Y <- Y[keep, , drop = FALSE]; grp <- grp[keep]; ids <- ids[keep]
cat(sprintf("Curves: %d, time points: %d\n", nrow(Y), ncol(Y)))

# ---- smoothing, on the REAL time axis (the app's opt-in real-time option) ----
rng   <- range(tlin)
basis <- create.bspline.basis(rangeval = rng, breaks = tlin, norder = 4)
fdP   <- fdPar(basis, 2, 1e-1)
Ys <- matrix(NA_real_, nrow(Y), length(tlin))
for (i in seq_len(nrow(Y))) {
  ok <- !is.na(Y[i, ])
  if (sum(ok) < 4) next
  s <- tryCatch(smooth.basis(tlin[ok], Y[i, ok], fdP), error = function(e) NULL)
  if (!is.null(s)) Ys[i, ] <- as.vector(eval.fd(tlin, s$fd))
}
good <- !apply(is.na(Ys), 1, any)
Ys <- Ys[good, , drop = FALSE]; grp <- grp[good]; ids <- ids[good]
cat(sprintf("Smoothed: %d curves\n", nrow(Ys)))

fdobj <- smooth.basis(tlin, t(Ys), fdP)$fd
pca <- pca.fd(fdobj, nharm = 5)
cat("Variance explained:", paste(sprintf("%.1f%%", 100 * pca$varprop), collapse = ", "), "\n\n")

res <- fck_pc_anova_all(
  scores = pca$scores, g = factor(grp),
  eigenvalues = pca$values, varprop = pca$varprop,
  posthoc_method = "games-howell",
  across_pc_correction = "holm",
  posthoc_gate = 0.05, subject_ids = ids)

cat(strrep("=", 76), "\nCOMPONENT ANOVA — AGEcategory\n", strrep("=", 76), "\n", sep = "")
ic <- res$independence
if (!isTRUE(ic$ok))
  cat(sprintf("\n! %d score rows from %d participants (%d repeats): these tests treat\n  correlated rows as independent and are ANTICONSERVATIVE.\n",
              ic$n_rows, ic$n_subjects, ic$n_repeated))

cat(sprintf("\n%-5s %8s %11s %12s %11s %10s  %s\n",
            "PC", "var%", "test", "p", "p adj (Holm)", "omega2", "separation"))
cat(strrep("-", 82), "\n")
for (j in seq_len(res$k)) {
  o <- res$omnibus[[j]]; if (is.null(o)) next
  usedw <- is.finite(o$welch_p) && isTRUE(all.equal(res$p_omnibus[j], o$welch_p))
  cat(sprintf("%-5s %7s%% %11s %12s %11s %10s  %s\n", paste0("PC", j),
              fmt1(100 * res$varprop[j]), if (usedw) "Welch F" else "Fisher F",
              format.pval(res$p_omnibus[j], digits = 3, eps = 1e-16),
              format.pval(res$p_omnibus_adj[j], digits = 3, eps = 1e-16),
              fmt3(o$omega2), res$separation$stability[j]))
}
for (j in seq_len(res$k)) {
  ph <- res$posthoc[[j]]; if (is.null(ph)) next
  o <- res$omnibus[[j]]
  cat(sprintf("\nPC%d post-hoc (Games-Howell, family-wise by construction):\n", j))
  for (lv in o$levels)
    cat(sprintf("   %-12s mean %10s (SD %8s, n = %d)\n", lv,
                fmt2(o$means[lv]), fmt2(o$sds[lv]), o$ns[lv]))
  cat(sprintf("   %-12s %-12s %10s %22s %10s %9s\n",
              "A", "B", "diff", "95% CI", "p", "Hedges g"))
  for (i in seq_len(nrow(ph))) {
    r <- ph[i, ]
    cat(sprintf("   %-12s %-12s %10s %22s %10s %9s\n", r$a, r$b, fmt2(r$diff),
                sprintf("[%s, %s]", fmt2(r$ci_lo), fmt2(r$ci_hi)),
                format.pval(r$p_adj, digits = 3, eps = 1e-16), fmt2(r$hedges_g)))
  }
}
cat("\nComponents with no post-hoc table were closed by the across-component\ncorrection, not tested and dismissed.\n")
