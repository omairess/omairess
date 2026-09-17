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
  # amplitude and phase are now per harmonic -- see the harmonic block below
  expect_true(all(c("level", "amplitude1", "phase1", "full", "shape", "circadian",
                    "trend") %in% ws))
  f1 <- mk(K = 1, trend = "none"); skip_if_not(isTRUE(f1$ok))
  w1 <- e$dance_traj_pair_whats(f1)
  expect_false("trend" %in% w1)
  # with no trend, shape and circadian are the SAME columns: one entry, not two
  expect_equal(sum(c("shape", "circadian") %in% w1), 1L)
  expect_true("full" %in% w1)
  # every entry has a label, block or harmonic alike, and none is empty
  labs <- vapply(ws, e$dance_traj_pair_label, character(1),
                 period = f2$spec$period, n_harmonics = f2$spec$n_harmonics)
  expect_true(all(nzchar(labs)))
  expect_false(any(labs %in% ws))          # a label, not the raw key echoed back
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

# ==============================================================================
# WHICH HARMONIC. "Amplitude" and "Acrophase" are properties of ONE harmonic, and
# the panel took harmonic = 1 from a default argument the UI never set -- so a
# two-harmonic model silently compared the 24 h component under a bare label, and
# the 12 h component could not be compared at all. The harmonic is now part of
# the name and of every label.
# ==============================================================================

mk_h2 <- function(nper = 10, seed = 9) {
  # H1 amplitude differs A vs B; H2 amplitude differs A vs C. If the panel can
  # only ever show H1 these two planted effects are indistinguishable.
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  grp <- rep(c("A", "B", "C"), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    20 + rnorm(1, 0, 3) + (8 + 4 * (g == "B")) * cos(2*pi*(tp - 2)/24) +
      (2 + 4 * (g == "C")) * cos(2*pi*2*(tp - 5)/24) + rnorm(length(tp), 0, 2)
  }))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 2, "none")
  ff <- e$dance_traj_fit(sp); ff$clock_origin <- 8; ff
}

test_that("every fitted harmonic gets its own amplitude and acrophase entry", {
  skip_if_not_installed("lme4")
  ff <- mk_h2(); skip_if_not(isTRUE(ff$ok))
  ws <- e$dance_traj_pair_whats(ff)
  expect_true(all(c("amplitude1", "amplitude2", "phase1", "phase2") %in% ws))
  expect_false("amplitude" %in% ws)          # nothing offered is bare
  expect_false("phase" %in% ws)
  # and each label names the harmonic AND its period, which is what makes the
  # number interpretable
  expect_match(e$dance_traj_pair_label("amplitude1", 24, 2), "^Amplitude .* H1 \\(24 h\\)$")
  expect_match(e$dance_traj_pair_label("phase2", 24, 2), "^Acrophase .* H2 \\(12 h\\)$")
  expect_match(e$dance_traj_pair_label("phase1", 12, 3), "H1 \\(12 h\\)")
  expect_match(e$dance_traj_pair_label("amplitude3", 24, 3), "H3 \\(8 h\\)")
  expect_equal(e$dance_traj_pair_label("level", 24, 2), "Level (at t = 0)")
})

test_that("the harmonic in the name is the harmonic compared", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  ff <- mk_h2(); skip_if_not(isTRUE(ff$ok))
  a1 <- e$dance_traj_contrasts(ff, "amplitude1", adjust = "none")
  a2 <- e$dance_traj_contrasts(ff, "amplitude2", adjust = "none")
  expect_true(isTRUE(a1$ok) && isTRUE(a2$ok))
  expect_false(isTRUE(all.equal(a1$table$estimate, a2$table$estimate)))
  g <- function(ct, c1, c2) ct$table$estimate[ct$table$cell1 == c1 & ct$table$cell2 == c2]
  # H1 separates A from B and not A from C; H2 does the opposite. Planted at 4.
  expect_gt(g(a1, "A", "B"), 3); expect_lt(abs(g(a1, "A", "C")), 1.5)
  expect_gt(g(a2, "A", "C"), 3); expect_lt(abs(g(a2, "A", "B")), 1.5)
  # the unit line names the harmonic and its period, so a printed table cannot
  # be read as being about the wrong rhythm
  expect_match(a1$unit, "H1"); expect_match(a1$unit, "24")
  expect_match(a2$unit, "H2"); expect_match(a2$unit, "12")
  p2 <- e$dance_traj_contrasts(ff, "phase2", adjust = "none")
  expect_match(p2$unit, "H2 acrophase")
  expect_match(p2$label, "Acrophase .* H2")
})

test_that("a bare name still means H1, so existing callers are unaffected", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  ff <- mk_h2(); skip_if_not(isTRUE(ff$ok))
  for (nm in c("amplitude", "phase")) {
    bare <- e$dance_traj_contrasts(ff, nm, adjust = "none")
    one  <- e$dance_traj_contrasts(ff, paste0(nm, "1"), adjust = "none")
    expect_true(isTRUE(bare$ok), info = nm)
    expect_equal(bare$table$estimate, one$table$estimate, info = nm)
    expect_equal(bare$table$p_raw, one$table$p_raw, info = nm)
  }
  expect_equal(e$dance_traj_pair_parse("amplitude")$harmonic, 1L)
  expect_equal(e$dance_traj_pair_parse("phase3")$harmonic, 3L)
  expect_true(is.na(e$dance_traj_pair_parse("circadian")$harmonic))
})

test_that("a harmonic the model does not fit is refused by number, not by name", {
  skip_if_not_installed("lme4")
  ff <- mk_h2(); skip_if_not(isTRUE(ff$ok))
  bad <- e$dance_traj_contrasts(ff, "amplitude5")
  expect_false(isTRUE(bad$ok))
  expect_match(bad$message, "fits 2 harmonic\\(s\\)")
  expect_match(bad$message, "no H5")
  # an unknown name still lists what IS available
  other <- e$dance_traj_contrasts(ff, "nonsense")
  expect_match(other$message, "not a comparable quantity")
  expect_match(other$message, "amplitude2")
})

test_that("simple effects take the harmonic too", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  set.seed(4)
  tp <- seq(0, 22, length.out = 12); nper <- 8
  sub <- sprintf("S%03d", seq_len(nper * 2)); sg <- rep(c("g1", "g2"), each = nper)
  rows <- list(); meta <- list(subject = character(), Group = character(), Cond = character())
  k <- 1L
  for (i in seq_along(sub)) for (cc in c("p", "q")) {
    rows[[k]] <- 20 + rnorm(1, 0, 3) + 8 * cos(2*pi*(tp - 2)/24) +
      (2 + 3 * (sg[i] == "g2")) * cos(2*pi*2*tp/24) + rnorm(length(tp), 0, 2)
    meta$subject <- c(meta$subject, sub[i]); meta$Group <- c(meta$Group, sg[i])
    meta$Cond <- c(meta$Cond, cc); k <- k + 1L
  }
  sp <- e$dance_traj_spec(e$dance_traj_long(do.call(rbind, rows), tp, meta$subject,
                            list(Group = meta$Group, Cond = meta$Cond)), 24, 2, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  se1 <- e$dance_traj_simple_effects(ff, "Group", at = list(Cond = "p"), what = "amplitude1")
  se2 <- e$dance_traj_simple_effects(ff, "Group", at = list(Cond = "p"), what = "amplitude2")
  expect_true(isTRUE(se1$ok) && isTRUE(se2$ok))
  expect_false(isTRUE(all.equal(se1$table$estimate, se2$table$estimate)))
  expect_match(se2$unit, "H2")
  # the H2 difference is the planted one
  expect_gt(abs(se2$table$estimate[1]), abs(se1$table$estimate[1]))
})
