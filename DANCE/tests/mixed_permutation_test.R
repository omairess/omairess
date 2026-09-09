# ==============================================================================
# tests/mixed_permutation_test.R — is the mixed permutation scheme EXACT?
#
# The one-way kernels' authority comes from exchangeability. The mixed kernel
# claims the same for all three effects, including the interaction, which is
# usually said to have no exact scheme. A claim like that is worth nothing
# argued and everything measured, so this file measures it: simulate under each
# null and check the rejection rate against the nominal level.
#
# The interaction null is the demanding one. It is simulated with STRONG main
# effects in both factors present, because a scheme that is only valid when
# nothing else is going on is not useful.
#
# Run with:   Rscript tests/mixed_permutation_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/06_helpers_mixed.R"), local = e)
source(file.path(app_dir, "server/07_helpers_mixed_perm.R"), local = e)

mk <- function(n_g, n_t, n_cond, within_eff, between_eff, inter_eff, seed) {
  set.seed(seed)
  g <- rep(c("g1", "g2"), each = n_g); ns <- length(g)
  tp <- seq(0, 1, length.out = n_t)
  subj <- rnorm(ns, 0, 3)
  wave <- sin(2 * pi * tp)
  rows <- list(); k <- 1L
  for (i in seq_len(ns)) for (ci in seq_len(n_cond)) {
    cc <- (ci - 1) / max(1, n_cond - 1)
    y <- subj[i] +
      within_eff  * cc * wave +                                  # within main
      between_eff * (g[i] == "g2") * cos(2 * pi * tp) +          # between main
      inter_eff   * cc * (g[i] == "g2") * wave +                 # interaction
      rnorm(n_t, 0, 2)
    rows[[k]] <- data.frame(subject = sprintf("S%02d", i), within = paste0("c", ci),
                            between = g[i], t = tp, y = y, stringsAsFactors = FALSE)
    k <- k + 1L
  }
  d <- do.call(rbind, rows)
  d$subject <- factor(d$subject); d$within <- factor(d$within); d$between <- factor(d$between)
  d$cell <- droplevels(interaction(d$within, d$between, sep = " x ", drop = TRUE))
  d
}

calibrate <- function(effect, n_cond, nsim, B, ...) {
  p <- vapply(seq_len(nsim), function(s) {
    d <- mk(n_g = 10, n_t = 10, n_cond = n_cond, seed = 1000 + s, ...)
    r <- e$dance_mixed_permutation(d, n_permutations = B, effects = effect)
    r[[effect]]$global_p
  }, numeric(1))
  list(rate = mean(p <= 0.05), ks = suppressWarnings(stats::ks.test(p, "punif")$p.value))
}

cat("-- calibration under the null (nominal .05) ----------------------------\n")
NS <- 400; B <- 299

# INTERACTION null, with BOTH main effects strongly present
r <- calibrate("interaction", n_cond = 2, nsim = NS, B = B,
               within_eff = 8, between_eff = 6, inter_eff = 0)
chk(abs(r$rate - .05) < .025 && r$ks > .01,
    sprintf("interaction, 2 conditions, strong main effects: rate %.3f, KS p %.3f",
            r$rate, r$ks),
    sprintf("interaction test is MIScalibrated: rate %.3f (KS p %.3f)", r$rate, r$ks))

# and with THREE within levels, which the 2-level difference argument does not cover
r3 <- calibrate("interaction", n_cond = 3, nsim = NS, B = B,
                within_eff = 8, between_eff = 6, inter_eff = 0)
chk(abs(r3$rate - .05) < .025 && r3$ks > .01,
    sprintf("interaction, 3 conditions: rate %.3f, KS p %.3f", r3$rate, r3$ks),
    sprintf("the 3-condition interaction test is MIScalibrated: rate %.3f (KS p %.3f)",
            r3$rate, r3$ks))

# WITHIN null, with a between effect present
rw <- calibrate("within", n_cond = 2, nsim = NS, B = B,
                within_eff = 0, between_eff = 6, inter_eff = 0)
chk(abs(rw$rate - .05) < .025 && rw$ks > .01,
    sprintf("within main effect: rate %.3f, KS p %.3f", rw$rate, rw$ks),
    sprintf("the within test is MIScalibrated: rate %.3f (KS p %.3f)", rw$rate, rw$ks))

# BETWEEN null, with a within effect present
rb <- calibrate("between", n_cond = 2, nsim = NS, B = B,
                within_eff = 8, between_eff = 0, inter_eff = 0)
chk(abs(rb$rate - .05) < .025 && rb$ks > .01,
    sprintf("between main effect: rate %.3f, KS p %.3f", rb$rate, rb$ks),
    sprintf("the between test is MIScalibrated: rate %.3f (KS p %.3f)", rb$rate, rb$ks))

cat("\n-- power: a real effect must be found --------------------------------\n")
d <- mk(10, 10, 2, within_eff = 8, between_eff = 6, inter_eff = 7, seed = 99)
r <- e$dance_mixed_permutation(d, n_permutations = 499)
for (nm in c("within", "between", "interaction"))
  chk(r[[nm]]$global_p < .05,
      sprintf("%-12s planted effect detected (global p = %.4f)", nm, r[[nm]]$global_p),
      sprintf("%-12s planted effect MISSED (global p = %.4f)", nm, r[[nm]]$global_p))

cat("\n-- the promises the output makes -------------------------------------\n")
chk(all(r$interaction$p_values >= r$p_floor - 1e-12),
    sprintf("no p falls below the Monte Carlo floor 1/(B+1) = %.4f", r$p_floor),
    "a p-value below the permutation floor was reported")
chk(all(r$interaction$p_adjusted >= r$interaction$p_values - 1e-12),
    "adjusted p values are never smaller than raw ones",
    "an adjusted p is smaller than its raw p")
chk(length(r$interaction$statistic) == length(unique(d$t)),
    "one pointwise statistic per evaluation point",
    "the pointwise statistic has the wrong length")

# ---- missing observations: the scheme stays exact, a missing CELL does not --
# Requiring no missing value at all would refuse ordinary data (the dataset this
# was built for is 1.2% missing). What exactness needs is that relabelling be a
# symmetry, and it is per time point provided a participant contributes there
# only when every condition is observed. A missing CELL is different: it is a
# missing half of the contrast, and no relabelling can invent it.
set.seed(4)
d_gap <- d; drop <- sample(nrow(d_gap), 40)
d_gap <- d_gap[-drop, ]
arr <- e$dance_mixed_array(d_gap)
chk(isTRUE(arr$complete) && any(!arr$usable),
    sprintf("scattered missing values leave every cell present (%d unusable time cells)",
            sum(!arr$usable)),
    "the gappy fixture did not produce the case under test")
r_gap <- e$dance_mixed_permutation(d_gap, n_permutations = 199)
chk(isTRUE(r_gap$ok) && isTRUE(r_gap$any_missing),
    "a dataset with missing observations is analysed and flagged as such",
    "missing observations were refused, which is stricter than exactness requires")
chk(isTRUE(r_gap$ok) && all(is.finite(r_gap$interaction$statistic)),
    "every pointwise statistic is finite despite the gaps",
    "a pointwise statistic came back non-finite on gappy data")

cal_gap <- local({
  p <- vapply(seq_len(200), function(s) {
    dd <- mk(n_g = 10, n_t = 10, n_cond = 2, seed = 7000 + s,
             within_eff = 8, between_eff = 6, inter_eff = 0)
    set.seed(8000 + s); dd <- dd[-sample(nrow(dd), floor(0.03 * nrow(dd))), ]
    e$dance_mixed_permutation(dd, n_permutations = 299, effects = "interaction")$interaction$global_p
  }, numeric(1))
  mean(p <= .05)
})
chk(abs(cal_gap - .05) < .035,
    sprintf("with 3%% of observations missing the interaction still rejects at %.3f", cal_gap),
    sprintf("missing data breaks calibration: rate %.3f", cal_gap))

# an incomplete design must be REFUSED, not silently approximated
d2 <- d[!(d$subject == levels(d$subject)[1] & d$within == levels(d$within)[2]), ]
r2 <- e$dance_mixed_permutation(d2, n_permutations = 49)
chk(!isTRUE(r2$ok) && grepl("every level", r2$message),
    "an incomplete design is refused, with the reason and the alternative",
    "an incomplete design was accepted, so the scheme is not exact for it")

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (bad == 0) "Mixed permutation tests PASSED" else "Mixed permutation tests FAILED",
            ok_n, bad))
quit(status = if (bad == 0) 0 else 1)
