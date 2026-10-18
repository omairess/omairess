# ==============================================================================
# tests/traj_car1_calibration.R -- does the mixed cosinor's test keep its size
# with autocorrelated residuals, modelled or not?
#
# A validation study, not a routine test (it takes tens of minutes). Data:
# two groups of 12 participants, 16 uneven times from 08:00 to 15:00 the next
# day, one harmonic, a random intercept and random rhythm per participant, and
# residuals that are either independent or a continuous-time AR(1) with
# phi = 0.6 per hour (two hours apart: 0.36). Under the null the groups share
# their rhythm; under the alternative group B's acrophase is 0.75 h later.
#
# For each data set the circadian block test of Group (Kenward-Roger) is run on
#   - the model with independent residuals (the default), and
#   - the model with continuous-time AR(1) residuals (phi estimated, REML).
# Reported: rejection rate at .05 (size under the null, power under the
# alternative), and for the AR(1) model the mean phi-hat and how often its 95%
# profile interval covers the true phi.
#
# Run with:   Rscript tests/traj_car1_calibration.R [reps]   (from DANCE)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else "DANCE"
for (f in c("08d_helpers_traj.R", "08e_helpers_trajfit.R", "08f_helpers_trajinf.R",
            "08g_helpers_trajband.R", "08i_helpers_trajsat.R"))
  eval(parse(file.path(app_dir, "server", f), encoding = "UTF-8"), envir = globalenv())
suppressMessages(suppressPackageStartupMessages({ library(lme4); library(parallel) }))

reps <- as.integer(commandArgs(TRUE)[1] %||% 150)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 15)
lin <- clock; for (i in 2:length(lin)) if (lin[i] < lin[i - 1]) lin[i:length(lin)] <- lin[i:length(lin)] + 24
tt <- lin - lin[1]
grp <- rep(c("A", "B"), each = 12)

simulate <- function(seed, phi, shift_h) {
  set.seed(seed)
  L <- if (phi > 0) t(chol(phi^abs(outer(tt, tt, "-")))) else diag(length(tt))
  Y <- t(vapply(seq_along(grp), function(i) {
    acro <- 16 + if (grp[i] == "B") shift_h else 0
    20 + rnorm(1, 0, 1.5) + (4 + rnorm(1, 0, 0.8)) * cos(2 * pi * (lin - acro) / 24) +
      rnorm(1, 0, 0.8) * sin(2 * pi * lin / 24) + as.numeric(L %*% rnorm(length(tt), 0, 1.2))
  }, numeric(length(tt))))
  Y[cbind(sample(nrow(Y), 4), sample(ncol(Y), 4))] <- NA        # a few missing points
  dance_traj_long(Y, tt, sprintf("S%02d", seq_along(grp)), list(Group = grp))
}

one <- function(seed, phi, shift_h) {
  d <- simulate(seed, phi, shift_h)
  out <- c(p_ind = NA_real_, p_car1 = NA_real_, phi_hat = NA_real_, cover = NA_real_,
           rung_ind = NA_real_, rung_car1 = NA_real_)
  fi <- tryCatch(dance_traj_run(d, 24, 1, "none", roles = c(Group = "between"), t0 = 0),
                 error = function(e) NULL)
  if (isTRUE(fi$ok)) {
    ti <- dance_traj_marginal_test(fi, "Group", "circadian", df_method = "kr")
    if (isTRUE(ti$ok)) out["p_ind"] <- ti$p
    out["rung_ind"] <- fi$re_rung
  }
  fa <- tryCatch(dance_traj_run(d, 24, 1, "none", roles = c(Group = "between"), t0 = 0,
                                residuals = "car1"), error = function(e) NULL)
  if (isTRUE(fa$ok)) {
    ta <- dance_traj_marginal_test(fa, "Group", "circadian", df_method = "kr")
    if (isTRUE(ta$ok)) out["p_car1"] <- ta$p
    out["phi_hat"] <- fa$residual_ar$phi
    out["cover"] <- as.numeric(fa$residual_ar$ci[1] <= phi && fa$residual_ar$ci[2] >= phi)
    out["rung_car1"] <- fa$re_rung
  }
  out
}

scenarios <- list(
  list(name = "null, independent residuals", phi = 0,   shift = 0),
  list(name = "null, AR(1) residuals",       phi = 0.6, shift = 0),
  list(name = "0.75 h phase shift, AR(1)",   phi = 0.6, shift = 0.75))
cores <- max(1L, min(3L, parallel::detectCores() - 1L))
t_start <- Sys.time()
for (sc in scenarios) {
  res <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r)
    one(1000 + r, sc$phi, sc$shift), mc.cores = cores))
  rate <- function(p) { p <- p[is.finite(p)]; sprintf("%.3f (n = %d)", mean(p < 0.05), length(p)) }
  cat(sprintf("\n%s, %d replicates\n", sc$name, reps))
  cat(sprintf("  rejection at .05, independent-residual model: %s\n", rate(res[, "p_ind"])))
  cat(sprintf("  rejection at .05, AR(1)-residual model:       %s\n", rate(res[, "p_car1"])))
  cat(sprintf("  AR(1): mean phi-hat %.3f (true %.1f); interval covers the truth in %.3f\n",
              mean(res[, "phi_hat"], na.rm = TRUE), sc$phi, mean(res[, "cover"], na.rm = TRUE)))
  cat(sprintf("  random structure at rung 1: independent %.2f, AR(1) %.2f\n",
              mean(res[, "rung_ind"] == 1, na.rm = TRUE), mean(res[, "rung_car1"] == 1, na.rm = TRUE)))
}
cat(sprintf("\n%.1f min\n", as.numeric(Sys.time() - t_start, units = "mins")))
