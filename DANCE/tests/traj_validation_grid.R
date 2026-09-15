# ==============================================================================
# tests/traj_validation_grid.R — THE PRE-PHASE-4 VALIDATION GATE
# ==============================================================================
# THIS IS NOT A UNIT TEST. It is the simulation study that has to clear before
# the trajectory framework may become the default analysis, and it is deliberately
# separate from the fast suites so that neither can be mistaken for the other.
#
# WHY THE PILOT IS NOT ENOUGH. The phase-2 fix took the circadian block's type-I
# error from 0.275 to 0.020-0.050 across four null designs. Those were 60-100
# simulations per cell, and the exact binomial interval around 4/100 runs
# [0.011, 0.099]. That rules out 0.275. It does NOT distinguish a test sitting
# at .05 from one sitting at .09, and it says nothing at all about two
# harmonics, a trend term, unbalanced cells, or any block other than circadian.
# Reporting "calibrated" on that evidence would be the same kind of overclaim
# the 0.275 caveat existed to prevent.
#
# WHAT CLEARS THE GATE. For each cell below: at least N_MAIN simulations, the
# rejection rate, its exact (Clopper-Pearson) interval, and -- where an effect
# is planted -- power with its own interval. A cell PASSES when the exact
# interval for the null rate contains .05 and its upper limit is below .08.
# Cells that fail stay labelled validated = FALSE for that configuration; they
# do not block the cells that pass, and they are not argued away.
#
#   Rscript tests/traj_validation_grid.R                 # full gate, N = 1000
#   DANCE_GRID_N=50 Rscript tests/traj_validation_grid.R # smoke run of the grid
#   DANCE_GRID_CELLS=k1_nt_mixed,k2_nt_mixed Rscript ... # named cells only
#
# RUNTIME. At N = 1000 the full grid is thousands of mixed-model fits and is
# measured in hours, not minutes. It is meant to be run deliberately, once,
# before promotion -- and again whenever layer B or C changes.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

N_MAIN <- suppressWarnings(as.integer(Sys.getenv("DANCE_GRID_N")))
if (is.na(N_MAIN) || N_MAIN < 1) N_MAIN <- 1000L
ONLY <- Sys.getenv("DANCE_GRID_CELLS")
ONLY <- if (nzchar(ONLY)) trimws(strsplit(ONLY, ",")[[1]]) else NULL
PASS_UPPER <- 0.08

# ------------------------------------------------------------------------------
# ONE GENERATOR, covering every axis the grid varies
# ------------------------------------------------------------------------------
# Random variation is placed at BOTH the participant and the curve level,
# because that is what repeated-measures data looks like and it is the harder
# case: it is the configuration that measured 0.283 before the curve-level rung
# existed.
# An irregular observation protocol: 19 times over just past one period, no two
# gaps alike. Fixed across replicates because that is how a study's schedule
# behaves -- awkward, but the same awkward for everyone.
TIMES_IRREG <- c(0, 1.4, 2.1, 4.3, 5.0, 7.2, 8.9, 9.5, 11.1, 12.8,
                 13.3, 15.6, 16.2, 17.9, 19.1, 20.4, 21.0, 22.3, 23.1)

gen <- function(seed, n_per_group = 8, groups = c("a", "b"), conds = c("p", "q"),
                nt = 12, amp = 4, amp2 = 0, slope = 0, noise = 1,
                shift = 0, shift_cell = NULL, level = 0, level_group = NULL,
                drop_frac = 0, times = NULL) {
  set.seed(seed)
  tp <- times %||% seq(0, 22, length.out = nt)
  subj <- sprintf("S%03d", seq_len(n_per_group * length(groups)))
  sg <- rep(groups, each = n_per_group)
  rows <- list(); k <- 1L
  meta <- list(subject = character(0), Group = character(0), Condition = character(0))
  for (i in seq_along(subj)) {
    sb <- rnorm(1, 0, 2); sa <- rnorm(1, 0, .4); sp <- rnorm(1, 0, .2)
    for (cc in conds) {
      base <- 10 + sb + rnorm(1, 0, 2)
      a <- amp + sa + rnorm(1, 0, .4)
      ph <- 2 + sp + rnorm(1, 0, .2)
      if (!is.null(shift_cell) && paste(sg[i], cc) == shift_cell) ph <- ph + shift
      lv <- if (!is.null(level_group) && sg[i] == level_group) level else 0
      y <- base + lv + a * cos(2 * pi * tp / 24 - ph) + slope * tp +
           amp2 * cos(2 * pi * 2 * tp / 24 - 1) + rnorm(length(tp), 0, noise)
      rows[[k]] <- y
      meta$subject <- c(meta$subject, subj[i]); meta$Group <- c(meta$Group, sg[i])
      meta$Condition <- c(meta$Condition, cc); k <- k + 1L
    }
  }
  Y <- do.call(rbind, rows)
  # MODERATE IMBALANCE, as missing observations rather than missing curves:
  # that is how real data arrives, and it is what makes cells unbalanced.
  if (drop_frac > 0) {
    n_drop <- round(drop_frac * length(Y))
    Y[sample.int(length(Y), n_drop)] <- NA_real_
  }
  list(Y = Y, t = tp, meta = meta)
}
build <- function(g, within = TRUE) {
  fl <- list(Group = g$meta$Group)
  if (within) fl$Condition <- g$meta$Condition
  e$dance_traj_long(g$Y, g$t, g$meta$subject, fl)
}

# ------------------------------------------------------------------------------
# THE GRID
# ------------------------------------------------------------------------------
# Each cell names: the generator arguments, the model (harmonics, trend), the
# block to test, and whether an effect was planted (so the number is power
# rather than size).
cell <- function(id, desc, args = list(), K = 1, trend = "none",
                 block = "circadian", within = TRUE, planted = FALSE, n = NULL)
  list(id = id, desc = desc, args = args, K = K, trend = trend, block = block,
       within = within, planted = planted, n = n)

GRID <- list(
  # --- harmonics -------------------------------------------------------------
  cell("k1_nt_mixed",   "K=1, no trend, 2x2 mixed, 16 participants"),
  cell("k2_nt_mixed",   "K=2, no trend, 2x2 mixed, 16 participants",
       list(amp2 = 2), K = 2),
  # --- trend -----------------------------------------------------------------
  cell("k1_lin_mixed",  "K=1, linear trend, 2x2 mixed, 16 participants",
       list(slope = 0.08), trend = "linear"),
  cell("k2_lin_mixed",  "K=2, linear trend, 2x2 mixed, 16 participants",
       list(amp2 = 2, slope = 0.08), K = 2, trend = "linear"),
  # --- design shape ----------------------------------------------------------
  cell("k1_nt_between", "K=1, no trend, between-participant only, 16 participants",
       list(conds = "p"), within = FALSE),
  cell("k1_nt_within",  "K=1, no trend, within-participant only (3 conditions)",
       list(groups = "a", conds = c("p", "q", "r"), n_per_group = 16)),
  # --- sample size -----------------------------------------------------------
  cell("k1_nt_n8",      "K=1, no trend, 2x2 mixed, 8 participants",
       list(n_per_group = 4)),
  cell("k1_nt_n32",     "K=1, no trend, 2x2 mixed, 32 participants",
       list(n_per_group = 16)),
  cell("k1_nt_n64",     "K=1, no trend, 2x2 mixed, 64 participants",
       list(n_per_group = 32), n = max(200L, N_MAIN %/% 4L)),
  # --- balance ---------------------------------------------------------------
  cell("k1_nt_unbal",   "K=1, no trend, 2x2 mixed, 10% of observations missing",
       list(drop_frac = 0.10)),
  cell("k1_nt_unbal25", "K=1, no trend, 2x2 mixed, 25% of observations missing",
       list(drop_frac = 0.25)),
  # --- the other two omnibus blocks ------------------------------------------
  cell("k1_lin_shape",  "SHAPE block, K=1, linear trend, 2x2 mixed",
       list(slope = 0.08), trend = "linear", block = "shape"),
  cell("k1_lin_full",   "FULL block, K=1, linear trend, 2x2 mixed",
       list(slope = 0.08), trend = "linear", block = "full"),
  cell("k1_nt_level",   "LEVEL block, K=1, no trend, 2x2 mixed", block = "level"),
  # --- the configuration a real study actually has ---------------------------
  # Asked for directly: the app warned that a user's fit was OUTSIDE this grid
  # on five counts at once -- full block, K = 2, a linear trend, unbalanced
  # cells and irregular times -- and the honest answer to "then simulate THAT"
  # is these cells. Every axis is varied one at a time above; real data varies
  # them together, and approximations that each survive alone can still fail in
  # combination, so the combination has to be its own cell.
  #
  # Deliberately at 32 participants and not at the ~1300 of the study that
  # prompted it. Kenward-Roger and Satterthwaite are SMALL-sample corrections:
  # the approximation is under most strain when the variance components are
  # poorly determined, and it only improves as participants are added. A cell
  # that holds at n = 32 therefore covers the same structure at n = 1305, and
  # costs minutes instead of days. If it fails at 32 the honest next step is to
  # climb the sample size until it holds, not to assume the study is fine.
  cell("real_k2_lin_full",  "FULL block, K=2, linear trend, 4 between-groups, irregular + unbalanced",
       list(groups = c("a", "b", "c", "d"), conds = "p", n_per_group = 8,
            amp2 = 2, slope = 0.08, drop_frac = 0.10, times = TIMES_IRREG),
       K = 2, trend = "linear", block = "full", within = FALSE,
       n = max(200L, N_MAIN %/% 4L)),
  cell("real_k2_lin_circ",  "CIRCADIAN block, same configuration as real_k2_lin_full",
       list(groups = c("a", "b", "c", "d"), conds = "p", n_per_group = 8,
            amp2 = 2, slope = 0.08, drop_frac = 0.10, times = TIMES_IRREG),
       K = 2, trend = "linear", block = "circadian", within = FALSE,
       n = max(200L, N_MAIN %/% 4L)),
  # --- power, so size is never read on its own -------------------------------
  cell("pow_phase_05",  "POWER: 0.5 rad phase shift in one cell",
       list(shift = 0.5, shift_cell = "a q"), planted = TRUE,
       n = max(200L, N_MAIN %/% 4L)),
  cell("pow_level_3",   "POWER: level shift of 3 in one group, FULL block",
       list(level = 3, level_group = "a"), block = "full", planted = TRUE,
       n = max(200L, N_MAIN %/% 4L))
)
if (!is.null(ONLY)) GRID <- Filter(function(c) c$id %in% ONLY, GRID)

# ------------------------------------------------------------------------------
run_cell <- function(cl) {
  n <- cl$n %||% N_MAIN
  p <- rep(NA_real_, n); rung <- sing <- rep(NA_real_, n)
  t0 <- Sys.time()
  for (i in seq_len(n)) {
    g <- do.call(gen, c(list(seed = 90000 + i), cl$args))
    d <- build(g, within = cl$within)
    sp <- e$dance_traj_spec(d, 24, cl$K, cl$trend)
    if (!isTRUE(sp$ok)) next
    ff <- try(e$dance_traj_fit(sp), silent = TRUE)
    if (inherits(ff, "try-error") || !isTRUE(ff$ok)) next
    rung[i] <- ff$re_rung; sing[i] <- isTRUE(ff$singular)
    o <- try(e$dance_traj_omnibus(ff, cl$block), silent = TRUE)
    if (!inherits(o, "try-error") && isTRUE(o$ok)) p[i] <- o$p %||% NA_real_
  }
  pp <- p[is.finite(p)]
  k <- sum(pp <= .05); m <- length(pp)
  ci <- if (m > 0) stats::binom.test(k, m)$conf.int else c(NA_real_, NA_real_)
  rate <- if (m > 0) k / m else NA_real_
  verdict <- if (cl$planted) "power"
             else if (m == 0) "NO FITS"
             else if (ci[1] <= .05 && ci[2] <= PASS_UPPER) "PASS"
             else if (ci[1] > .05) "FAIL (too liberal)"
             else if (ci[2] > PASS_UPPER) "INCONCLUSIVE (interval too wide)"
             else "FAIL"
  data.frame(cell = cl$id, description = cl$desc, block = cl$block,
             n_sim = m, n_failed = n - m, rate = rate,
             ci_lo = ci[1], ci_hi = ci[2], verdict = verdict,
             top_rung_pct = 100 * mean(rung == 1, na.rm = TRUE),
             singular_pct = 100 * mean(sing, na.rm = TRUE),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
             stringsAsFactors = FALSE)
}

cat(sprintf("THE PRE-PHASE-4 VALIDATION GATE\n%s\nN = %d per cell, %d cell(s). A null cell PASSES when its exact\nbinomial interval contains .05 and its upper limit is <= %.2f.\n\n",
            strrep("=", 76), N_MAIN, length(GRID), PASS_UPPER))
res <- do.call(rbind, suppressWarnings(suppressMessages(lapply(GRID, function(cl) {
  r <- run_cell(cl)
  cat(sprintf("%-16s %-52s %s\n  rate %.4f  exact 95%% CI [%.4f, %.4f]  n = %4d (%d failed to fit)  rung1 %3.0f%%  singular %3.0f%%  %.1f min\n\n",
              r$cell, substr(r$description, 1, 52), r$verdict,
              r$rate, r$ci_lo, r$ci_hi, r$n_sim, r$n_failed,
              r$top_rung_pct, r$singular_pct, r$minutes))
  r
}))))

out <- file.path(if (dir.exists("tests")) "tests" else ".", "traj_validation_grid_result.csv")
utils::write.csv(res, out, row.names = FALSE)
nulls <- res[res$verdict != "power", , drop = FALSE]
cat(sprintf("%s\n%d of %d null cells PASS; %d fail; %d inconclusive.\nWritten to %s\n",
            strrep("=", 76), sum(nulls$verdict == "PASS"), nrow(nulls),
            sum(grepl("^FAIL", nulls$verdict)),
            sum(grepl("INCONCLUSIVE", nulls$verdict)), out))
cat(paste(
  "\nA cell that fails or is inconclusive does NOT block the cells that pass.",
  "\nIt means dance_traj_calibration() keeps returning validated = FALSE for that",
  "\nconfiguration, with the reason, until the cause is found and fixed --",
  "\nwhich is what the 0.275 episode established as the way to handle this.\n"))
