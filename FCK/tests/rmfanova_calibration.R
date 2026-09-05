# ==============================================================================
# tests/rmfanova_calibration.R -- which of rmfanova's fifteen outputs to trust
#
# rmfanova reports three statistics under five resampling schemes. They are not
# interchangeable. This measures Type-I error and power for all fifteen in one
# plausible setting, which is what FCK_RMFANOVA_DEFAULT is chosen from.
#
# One scenario is not a verdict on the package -- the authors characterise these
# tests across many more -- but it is enough to justify showing three of fifteen
# rather than all of them.
#
#   Rscript tests/rmfanova_calibration.R      (~4 minutes)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
if (!requireNamespace("rmfanova", quietly = TRUE))
  stop("rmfanova is not installed; this calibration cannot run.", call. = FALSE)
suppressMessages(library(rmfanova))

n <- 15; p <- 20; l <- 3; t <- seq(0, 1, length.out = p)
gen <- function(shift_fun) {
  subj <- matrix(rnorm(n, 0, 1.5), n, p)
  lapply(seq_len(l), function(j)
    subj + matrix(rep(shift_fun(j), each = n), n, p) +
           matrix(rnorm(n * p, 0, 0.6), n, p))
}
run <- function(shift_fun, R) {
  set.seed(99)
  P <- NULL
  for (r in seq_len(R)) {
    x <- gen(shift_fun)
    v <- tryCatch(as.numeric(rmfanova(x, n_perm = 200, n_boot = 200)$p_values),
                  error = function(e) stop("rmfanova failed: ", conditionMessage(e),
                                           call. = FALSE))
    P <- rbind(P, v)
  }
  colnames(P) <- c("Cn_P1","Cn_P2","Cn_B1","Cn_B2","Cn_B3",
                   "Dn_P1","Dn_P2","Dn_B1","Dn_B2","Dn_B3",
                   "En_P1","En_P2","En_B1","En_B2","En_B3")
  P
}
R <- 400
cat("=== CALIBRATION: true null (all conditions identical),", R, "replicates ===\n")
cat("nominal 5%; a usable test lands near 5%\n\n")
P0 <- run(function(j) 0 * t, R)
bump <- function(j) if (j == 3) 1.0 * exp(-((t - 0.5)^2) / 0.02) else 0 * t
cat("computing power...\n")
P1 <- run(bump, R)

res <- data.frame(
  method = colnames(P0),
  type1_5pct = round(100 * colMeans(P0 < 0.05, na.rm = TRUE), 1),
  power_5pct = round(100 * colMeans(P1 < 0.05, na.rm = TRUE), 1))
res$verdict <- ifelse(res$type1_5pct > 9, "INFLATED",
               ifelse(res$type1_5pct < 1.5, "over-conservative",
               ifelse(res$power_5pct < 30, "no power", "usable")))
print(res, row.names = FALSE)

# ---- assertions -------------------------------------------------------------
app_dir <- if (dir.exists("server")) "." else "FCK"
source(file.path(app_dir, "server/09b_helpers_rmfanova.R"))
def <- res[res$method %in% FCK_RMFANOVA_DEFAULT, ]
bad <- res[res$method %in% FCK_RMFANOVA_INFLATED, ]
stopifnot(
  "every default method must be near nominal" = all(def$type1_5pct < 9),
  "every default method must have power"      = all(def$power_5pct > 90),
  "the flagged methods really are inflated"   = all(bad$type1_5pct > 9)
)
cat("\nCalibration assertions passed: the curated default is the calibrated subset.\n")
