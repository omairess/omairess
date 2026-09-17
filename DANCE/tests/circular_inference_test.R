# ==============================================================================
# tests/circular_inference_test.R — the circular inference helpers, pinned
# ==============================================================================
# P21 phase 0 moved these out of the server body of 72_harmonic.R without
# changing a line of arithmetic (finding A8). "Without changing a line" is a
# claim, and this file is what makes it checkable: every value below was
# MEASURED from the pre-extraction code -- sliced out of git HEAD and evaluated
# beside the new file -- and all eight comparisons came back identical at
# tolerance = 0. The numbers are then hard-coded here so the pin survives
# without git, and a later edit that moves them fails rather than passes.
#
# (The first draft of this file carried values I had reasoned out rather than
# run. Six of the eight were wrong. They are measured now.)
#
# Run with:  Rscript tests/circular_inference_test.R      (from the DANCE dir)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/08b_helpers_popcosinor.R"), local = e)
source(file.path(app_dir, "server/08c_helpers_circstat.R"), local = e)

near <- function(a, b, tol = 1e-9) is.finite(a) && is.finite(b) && abs(a - b) < tol

cat("-- circular mean, resultant length, SD, SE ---------------------------\n")
# a cluster straddling the wrap: the case a linear mean gets wrong
ang <- c(0.05, 6.20, 0.10, 6.25, 0.00)
m <- e$dance_circular_mean(ang)
chk(near(m %% (2*pi), 0.0067224788, 1e-9),
    sprintf("circular mean of a wrap-straddling cluster is %.6f rad", m %% (2*pi)),
    sprintf("circular mean moved: %.6f", m %% (2*pi)))
chk(near(mean(ang), 2.52, 1e-9),
    "a LINEAR mean of the same angles is 2.52 rad -- which is why this helper exists",
    "the linear-mean contrast no longer holds")
chk(near(e$dance_mean_resultant_length(ang), 0.9979717381, 1e-9),
    sprintf("r-bar = %.7f", e$dance_mean_resultant_length(ang)),
    "r-bar changed")
chk(near(e$dance_circular_sd(ang), 0.0637231760, 1e-9),
    sprintf("circular SD = %.7f rad", e$dance_circular_sd(ang)),
    "circular SD changed")
chk(near(e$dance_circular_se(ang), 0.4481225053, 1e-9),
    sprintf("circular SE = %.7f", e$dance_circular_se(ang)),
    "circular SE changed")

cat("\n-- the degenerate inputs the reactives actually hit -------------------\n")
chk(is.nan(e$dance_circular_mean(numeric(0))),
    "circular mean of an empty vector is NaN, not an error",
    "the empty-vector behaviour changed")
chk(is.na(e$dance_circular_sd(c(NA_real_, NA_real_))),
    "circular SD of all-NA is NA",
    "all-NA no longer returns NA")
chk(identical(e$dance_circular_sd(rep(1.2, 5)), 0),
    "identical angles give a circular SD of exactly 0",
    "the zero-dispersion case changed")

cat("\n-- Watson-Williams ---------------------------------------------------\n")
set.seed(11)
g1 <- rnorm(20, 0.4, 0.30); g2 <- rnorm(20, 1.1, 0.30); g3 <- rnorm(20, 0.5, 0.30)
ww2 <- e$dance_watson_williams_test(list(g1, g2))
chk(near(ww2$F, 90.8681279306, 1e-8) && ww2$df1 == 1 && ww2$df2 == 38,
    sprintf("two groups: F(%d, %d) = %.6f, p = %.3g", ww2$df1, ww2$df2, ww2$F, ww2$p),
    sprintf("two-group Watson-Williams changed: F = %.6f on (%s, %s)",
            ww2$F, ww2$df1, ww2$df2))
ww3 <- e$dance_watson_williams_test(list(g1, g2, g3))
chk(near(ww3$F, 42.8806964228, 1e-8) && ww3$df1 == 2 && ww3$df2 == 57,
    sprintf("three groups: F(%d, %d) = %.6f, p = %.3g", ww3$df1, ww3$df2, ww3$F, ww3$p),
    sprintf("three-group Watson-Williams changed: F = %.6f on (%s, %s)",
            ww3$F, ww3$df1, ww3$df2))
chk(!is.null(e$dance_watson_williams_test(list(g1))$message),
    "a single group is refused",
    "a single group was accepted")

cat("\n-- Bingham parameter test / Hotelling T-squared -----------------------\n")
set.seed(11)
mats <- lapply(1:3, function(i) cbind(rnorm(20, i * 0.4), rnorm(20, -i * 0.3)))
ht <- e$dance_hotelling_t2(lapply(mats, function(m) m[, 1]),
                           lapply(mats, function(m) m[, 2]))
chk(near(ht$F, 3.644338, 1e-5) && ht$df1 == 4 && ht$df2 == 112,
    sprintf("three groups of 20: F(%d, %d) = %.6f, Wilks = %.6f",
            ht$df1, ht$df2, ht$F, ht$lambda),
    sprintf("the multigroup Bingham test changed: F = %.6f on (%s, %s)",
            ht$F, ht$df1, ht$df2))
# and it must still reduce to the exact two-sample T-squared
set.seed(5); n1 <- 14; n2 <- 11
m1 <- cbind(rnorm(n1, 0.9), rnorm(n1, 0.2)); m2 <- cbind(rnorm(n2), rnorm(n2))
h2 <- e$dance_hotelling_t2(list(m1[, 1], m2[, 1]), list(m1[, 2], m2[, 2]))
N <- n1 + n2
Sp <- ((n1 - 1) * cov(m1) + (n2 - 1) * cov(m2)) / (N - 2)
d  <- colMeans(m1) - colMeans(m2)
T2 <- (n1 * n2) / N * as.numeric(t(d) %*% solve(Sp) %*% d)
chk(near(h2$F, (N - 3) / ((N - 2) * 2) * T2, 1e-9) && h2$df1 == 2 && h2$df2 == N - 3,
    "two groups still reduce exactly to the two-sample Hotelling T-squared",
    "the two-group reduction broke")

cat("\n-- the extraction itself ---------------------------------------------\n")
h72 <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
             collapse = "\n")
chk(!any(vapply(c("circular_mean", "circular_sd", "circular_se",
                  "mean_resultant_length", "watson_williams_test", "hotelling_t2"),
                function(fn) grepl(sprintf("\n  %s <- function", fn), h72, fixed = TRUE),
                logical(1))),
    "no circular helper is defined inside the observer file any more",
    "a circular helper is still defined inside server/72_harmonic.R")

cat("\n")
if (bad) { cat(sprintf("Circular inference tests FAILED  (%d passed, %d failed)\n", ok_n, bad)); quit(status = 1) }
cat(sprintf("Circular inference tests PASSED  (%d passed, 0 failed)\n", ok_n))
