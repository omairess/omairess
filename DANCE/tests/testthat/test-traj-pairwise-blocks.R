# ==============================================================================
# tests/testthat/test-traj-pairwise-blocks.R
#
# Asked for directly: the secondary table reports an omnibus test for each
# component block, and the pairwise panel could only compare level, amplitude
# and acrophase. So a "Rhythmic block p = .029" had no follow-up telling you
# WHICH pair of cells carried it.
#
# The pairwise block contrast is cell i minus cell j over exactly the components
# the omnibus tests, run through the same L-matrix machinery -- not a second
# implementation. These pin that relationship, because two notions of what a
# block IS is how a significant omnibus and an empty pairwise panel end up in
# the same figure.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk <- function(K = 2, trend = "linear", nper = 10, seed = 9) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  grp <- rep(c("A", "B", "C"), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    20 + rnorm(1, 0, 3) + (8 + 3 * (g == "C")) * cos(2*pi*(tp - 2 - 2*(g == "B"))/24) +
      2 * cos(2*pi*2*tp/24) + 0.1 * (g == "C") * tp + rnorm(length(tp), 0, 2)
  }))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, K, trend)
  ff <- e$dance_traj_fit(sp); ff$clock_origin <- 8; ff
}

test_that("the pairwise vocabulary follows the fit, with no duplicate blocks", {
  skip_if_not_installed("lme4")
  f2 <- mk(); skip_if_not(isTRUE(f2$ok))
  ws <- e$dance_traj_pair_whats(f2)
  expect_true(all(c("level", "amplitude", "phase", "full", "shape", "circadian",
                    "trend") %in% ws))
  f1 <- mk(K = 1, trend = "none"); skip_if_not(isTRUE(f1$ok))
  w1 <- e$dance_traj_pair_whats(f1)
  expect_false("trend" %in% w1)
  # with no trend, shape and circadian are the SAME columns: one entry, not two
  expect_equal(sum(c("shape", "circadian") %in% w1), 1L)
  expect_true("full" %in% w1)
  expect_true(all(ws %in% names(e$DANCE_TRAJ_PAIR_LABEL)))
})

test_that("a pair contrast uses the same components as the omnibus it follows", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  grid <- e$dance_traj_cell_grid(ff$spec)
  for (b in c("full", "shape", "circadian", "trend", "level")) {
    comps <- e$DANCE_TRAJ_BLOCK_COMPONENTS(ff$spec, b)
    Lp <- e$dance_traj_pair_L(ff, grid$.cell[1], grid$.cell[2], b)
    expect_equal(nrow(Lp), length(comps), info = b)
    # each row is cell i minus cell j on ONE component, so it is the difference
    # of the two cells' functionals and nothing else
    fx <- e$dance_traj_cell_functionals(ff)
    get1 <- function(ci, cp) if (identical(cp, ".level@t0")) fx$cells[[ci]]$level0
                             else fx$cells[[ci]]$coef[[cp]]
    for (k in seq_along(comps)) {
      want <- get1(1, comps[k]) - get1(2, comps[k])
      want[abs(want) < 1e-12] <- 0
      expect_equal(unname(Lp[k, ]), unname(want[colnames(Lp)]), info = paste(b, comps[k]))
    }
  }
})

test_that("the numerator df is one cell-contrast worth of the omnibus", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  grid <- e$dance_traj_cell_grid(ff$spec)
  nlev <- nrow(grid)
  for (b in c("full", "shape", "circadian", "trend", "level")) {
    om <- e$dance_traj_marginal_test(ff, character(0), b, df_method = "satterthwaite")
    pw <- e$dance_traj_pair_block_test(ff, grid$.cell[1], grid$.cell[2], b,
                                       df_method = "satterthwaite")
    expect_true(isTRUE(om$ok) && isTRUE(pw$ok), info = b)
    # the omnibus spans (cells - 1) contrasts per component; a pair spans one
    expect_equal(om$df1, pw$df1 * (nlev - 1L), info = b)
  }
})

test_that("a one-component block comes back as a scalar with an interval", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  grid <- e$dance_traj_cell_grid(ff$spec)
  for (b in c("trend", "level")) {
    r <- e$dance_traj_pair_block_test(ff, grid$.cell[1], grid$.cell[2], b,
                                      df_method = "satterthwaite")
    expect_equal(r$df1, 1, info = b)
    expect_true(is.finite(r$estimate) && is.finite(r$se), info = b)
    expect_true(is.finite(r$lo) && is.finite(r$hi), info = b)
    expect_lt(r$lo, r$estimate); expect_gt(r$hi, r$estimate)
    # F on 1 df is t squared, and the interval excludes zero exactly when p < .05
    expect_equal(sqrt(r$statistic), abs(r$estimate / r$se), tolerance = 1e-8, info = b)
    expect_equal(r$p < 0.05, (r$lo > 0 || r$hi < 0), info = b)
  }
  # a multi-component block has no single number and says so
  rj <- e$dance_traj_pair_block_test(ff, grid$.cell[1], grid$.cell[2], "circadian",
                                     df_method = "satterthwaite")
  expect_gt(rj$df1, 1)
  expect_true(is.na(rj$lo) && is.na(rj$hi))
})

test_that("the LEVEL block pairwise agrees with the existing level contrast", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  # The scalar "Level (at t = 0)" entry evaluates the design rows at t0 and
  # differences them; the block path differences the cells' level functionals.
  # Two constructions of one quantity, so they must agree exactly -- in size AND
  # in sign, or the two panels would disagree about which cell is higher.
  ct <- e$dance_traj_contrasts(ff, "level", adjust = "none")
  expect_true(isTRUE(ct$ok))
  for (k in seq_len(nrow(ct$table))) {
    z <- ct$table[k, ]
    r <- e$dance_traj_pair_block_test(ff, z$cell1, z$cell2, "level",
                                      df_method = "satterthwaite")
    expect_equal(r$estimate, z$estimate, tolerance = 1e-8)
    expect_equal(r$se, z$se, tolerance = 1e-8)
    expect_equal(sign(r$estimate), sign(z$estimate))
  }
})

test_that("the contrast table renders for every block the fit offers", {
  skip_if_not_installed("lme4")
  ff <- mk(); skip_if_not(isTRUE(ff$ok))
  for (w in setdiff(e$dance_traj_pair_whats(ff), c("amplitude", "phase"))) {
    ct <- e$dance_traj_contrasts(ff, w, adjust = "holm", df_method = "satterthwaite")
    expect_true(isTRUE(ct$ok), info = w)
    expect_equal(nrow(ct$table), 3L, info = w)       # 3 cells -> 3 pairs
    expect_true(all(is.finite(ct$table$p_raw)), info = w)
    expect_true(all(ct$table$p_adj >= ct$table$p_raw - 1e-12), info = w)
    expect_true(nzchar(ct$note), info = w)
    if (w %in% e$DANCE_TRAJ_BLOCKS) {
      expect_true(all(is.finite(ct$table$statistic)), info = w)
      expect_identical(ct$joint, any(ct$table$df1 > 1), info = w)
    }
  }
})

test_that("a quantity this fit cannot compare is refused by name", {
  skip_if_not_installed("lme4")
  f1 <- mk(K = 1, trend = "none"); skip_if_not(isTRUE(f1$ok))
  bad <- e$dance_traj_contrasts(f1, "nonsense")
  expect_false(isTRUE(bad$ok))
  expect_match(bad$message, "not a comparable quantity")
  expect_match(bad$message, "Available:")
})

test_that("simple effects answer the same questions as the whole design", {
  skip_if_not_installed("lme4")
  set.seed(3)
  tp <- seq(0, 22, length.out = 12); nper <- 8
  sub <- sprintf("S%03d", seq_len(nper * 2))
  sg <- rep(c("g1", "g2"), each = nper)
  rows <- list(); meta <- list(subject = character(), Group = character(), Cond = character())
  k <- 1L
  for (i in seq_along(sub)) for (cc in c("p", "q")) {
    rows[[k]] <- 20 + rnorm(1, 0, 3) + 8 * cos(2*pi*(tp - 2 - (sg[i] == "g2"))/24) +
      rnorm(length(tp), 0, 2)
    meta$subject <- c(meta$subject, sub[i]); meta$Group <- c(meta$Group, sg[i])
    meta$Cond <- c(meta$Cond, cc); k <- k + 1L
  }
  sp <- e$dance_traj_spec(e$dance_traj_long(do.call(rbind, rows), tp, meta$subject,
                            list(Group = meta$Group, Cond = meta$Cond)), 24, 1, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  se <- e$dance_traj_simple_effects(ff, "Group", at = list(Cond = "p"),
                                    what = "circadian", df_method = "satterthwaite")
  expect_true(isTRUE(se$ok))
  expect_equal(nrow(se$table), 1L)
  expect_gt(se$table$df1[1], 1)
  expect_true(isTRUE(se$joint))
})
