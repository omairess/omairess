# ==============================================================================
# tests/testthat/test-traj-participant-cell.R
#
# dance_traj_participant_curves() added every participant's conditional modes to
# mean(co$a) -- the average over ALL design cells -- instead of to their own
# cell's fixed-effect rhythm. With two groups of different amplitude, every
# individual in the high group read low and every individual in the low group
# read high, all pulled toward a grand mean nobody belongs to.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk <- function(nper = 14, seed = 12) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  grp <- rep(c("lo", "hi"), each = nper)
  amp <- ifelse(grp == "hi", 12, 5)                  # planted: 12 vs 5
  Y <- t(sapply(seq_along(grp), function(i)
    30 + rnorm(1, 0, 3) + (amp[i] + rnorm(1, 0, 0.8)) * cos(2*pi*(tp - 3)/24) +
      rnorm(length(tp), 0, 1.5)))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 1, "none")
  list(fit = e$dance_traj_fit(sp), grp = grp)
}

test_that("each participant deviates from their OWN cell's rhythm", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk(); skip_if_not(isTRUE(m$fit$ok))
  pc <- e$dance_traj_participant_curves(m$fit, 1)
  skip_if_not(isTRUE(pc$ok), pc$message)
  tab <- pc$table
  expect_true(all(c("subject", "curve", "cell", "amplitude") %in% names(tab)))
  hi <- tab$amplitude[tab$cell == "hi"]; lo <- tab$amplitude[tab$cell == "lo"]
  expect_gt(length(hi), 0); expect_gt(length(lo), 0)
  # centred on their own group, not on the grand mean (8.5)
  expect_gt(mean(hi), 10.5); expect_lt(mean(lo), 6.5)
  co <- e$dance_traj_cell_coefs(m$fit, 1)
  A_cell <- sqrt(co$a^2 + co$b^2); names(A_cell) <- co$cells
  # the mean conditional mode is ~0 by construction, so the mean of a cell's
  # participants sits on that cell's own fixed-effect amplitude
  expect_lt(abs(mean(hi) - A_cell[["hi"]]), 1.0)
  expect_lt(abs(mean(lo) - A_cell[["lo"]]), 1.0)
  # and the old construction would NOT satisfy this: a grand-mean baseline
  # puts both groups' means at the same place
  expect_gt(abs(mean(hi) - mean(lo)), 4)
})

test_that("with no design factor there is one cell and nothing changes", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  set.seed(2)
  tp <- seq(0, 22, length.out = 12)
  Y <- t(sapply(1:16, function(i) 30 + rnorm(1,0,3) + 8*cos(2*pi*(tp-3)/24) + rnorm(12,0,1.5)))
  sp <- e$dance_traj_spec(e$dance_traj_long(Y, tp, sprintf("S%03d", 1:16), list()), 24, 1, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  pc <- e$dance_traj_participant_curves(ff, 1)
  skip_if_not(isTRUE(pc$ok), pc$message)
  expect_equal(nrow(pc$table), 16L)
  expect_equal(length(unique(pc$table$cell)), 1L)
  expect_lt(abs(mean(pc$table$amplitude) - 8), 1.5)
})
