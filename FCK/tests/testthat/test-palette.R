# ==============================================================================
# tests/testthat/test-palette.R — one group palette, keyed by name
#
# Seven palettes were in use across the app, so the same age band was a
# different colour in every figure. Worse, every one of them was indexed by
# POSITION: filter a group out, or sort the levels differently in one tab, and
# the survivors are repainted -- "the blue group" then means different things in
# two figures on the same screen. That is the "colour follows the entity, never
# its rank" rule, and it was broken everywhere.
#
# The palette itself was also re-stepped. The old one FAILED the dataviz
# validator in --pairs all mode at four groups (amber vs orange, normal-vision
# dE 13.7 against a floor of 15), which is exactly the reporter's case: four age
# bands, compared all-pairs in a scatter and an overlay.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/02b_helpers_palette.R"))

AGE <- c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")

# ------------------------------------------- colour follows the entity
test_that("a group keeps its colour when another is filtered out", {
  full <- fck_group_colors(AGE)
  sub  <- fck_group_colors(c("YOUTH", "ELDERLY"), all_levels = AGE)
  expect_equal(unname(sub[["YOUTH"]]),   unname(full[["YOUTH"]]))
  expect_equal(unname(sub[["ELDERLY"]]), unname(full[["ELDERLY"]]))
})

test_that("a group keeps its colour when the levels arrive in another order", {
  a <- fck_group_colors(AGE)
  b <- fck_group_colors(rev(AGE))
  for (g in AGE) expect_equal(unname(a[[g]]), unname(b[[g]]), info = g)
})

test_that("the assignment is by name, so it survives a differently-sorted tab", {
  a <- fck_group_colors(c("ADULT", "YOUTH"))
  b <- fck_group_colors(c("YOUTH", "ADULT"))
  expect_equal(unname(a[["YOUTH"]]), unname(b[["YOUTH"]]))
  expect_equal(unname(a[["ADULT"]]), unname(b[["ADULT"]]))
  # and the two groups are NOT the same colour as each other
  expect_false(identical(unname(a[["YOUTH"]]), unname(a[["ADULT"]])))
})

test_that("every group gets a distinct colour up to the validated depth", {
  for (n in 2:5) {
    cols <- fck_group_colors(LETTERS[seq_len(n)])
    expect_equal(length(unique(unname(cols))), n, info = paste("n =", n))
  }
})

test_that("NA and empty levels are dropped, not coloured", {
  cols <- fck_group_colors(c("YOUTH", NA, "", "ADULT"))
  expect_equal(sort(names(cols)), c("ADULT", "YOUTH"))
})

test_that("an unknown level reads gray rather than borrowing a colour", {
  cols <- fck_group_colors(c("YOUTH", "MYSTERY"), all_levels = AGE)
  expect_equal(unname(cols[["MYSTERY"]]), "#6f6f6f")
  expect_false(unname(cols[["YOUTH"]]) == "#6f6f6f")
})

# ------------------------------------------------- the validated palette
test_that("the validated prefix is exactly the five all-pairs-separable slots", {
  expect_equal(length(FCK_GROUP_COLORS), 5)
  expect_true(all(grepl("^#[0-9a-fA-F]{6}$", FCK_GROUP_COLORS)))
  # the amber that failed the all-pairs normal-vision check is gone
  expect_false("#eda100" %in% FCK_GROUP_COLORS)
})

test_that("dash is required past the validated depth, and only past it", {
  expect_false(fck_needs_dash(5))
  expect_true(fck_needs_dash(6))
})

test_that("the ramp extends past five rather than erroring", {
  expect_equal(length(fck_group_ramp(8)), 8)
  expect_equal(fck_group_ramp(3), FCK_GROUP_COLORS[1:3])
  expect_equal(length(fck_group_ramp(12)), 12)   # recycles, with dash carrying identity
})

test_that("dashes are keyed the same way as the colours", {
  d1 <- fck_group_dashes(AGE)
  d2 <- fck_group_dashes(rev(AGE))
  for (g in AGE) expect_equal(unname(d1[[g]]), unname(d2[[g]]), info = g)
})

# ------------------------------------------------------------- derivatives
test_that("fills and rgba derive from the colour, so they cannot drift", {
  hex <- FCK_GROUP_COLORS[1]
  expect_equal(fck_group_fill(hex, 0.25), paste0(hex, "40"))
  expect_equal(fck_group_fill(hex, 1), paste0(hex, "FF"))
  expect_equal(fck_group_fill(hex, 0), paste0(hex, "00"))
  expect_true(grepl("^rgba\\(", fck_group_rgba(hex, 0.4)))
  # the rgba channels really are this colour's
  m <- grDevices::col2rgb(hex)
  expect_true(grepl(sprintf("rgba\\(%d,%d,%d", m[1], m[2], m[3]),
                    fck_group_rgba(hex, 0.4)))
})

test_that("alpha is clamped rather than producing invalid hex", {
  hex <- FCK_GROUP_COLORS[1]
  expect_equal(fck_group_fill(hex, 5), paste0(hex, "FF"))
  expect_equal(fck_group_fill(hex, -1), paste0(hex, "00"))
})

# --------------------------------------------- no stray palettes survive
test_that("rainbow() is gone from the pairwise plots", {
  src <- paste(readLines(file.path(app_dir, "server/73_cosinor_pairwise.R"), warn = FALSE),
               collapse = "\n")
  expect_false(grepl("rainbow(n_groups", src, fixed = TRUE))
})

test_that("the harmonic tab no longer hard-codes its own hues", {
  src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
               collapse = "\n")
  for (dead in c("'#B22222', '#4682B4'", "'firebrick', 'steelblue'", "'#FF6B6B', '#87CEEB'"))
    expect_false(grepl(dead, src, fixed = TRUE), info = dead)
  expect_true(grepl("fck_group_colors(", src, fixed = TRUE))
})

test_that("the polar density draws from the shared palette", {
  src <- paste(readLines(file.path(app_dir, "server/74_polar_density.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl("FCK_DENSITY_COLORS <- c(FCK_GROUP_COLORS", src, fixed = TRUE))
  # and keys on name, not position
  expect_false(grepl("FCK_DENSITY_COLORS[((i - 1)", src, fixed = TRUE))
})

# ------------------------------------------------------ warping parameters
source(file.path(app_dir, "server/08_helpers_cosinor.R"))
source(file.path(app_dir, "server/09_helpers_pcanova.R"))

test_that("warping parameters are extracted for a shift-based method", {
  tp <- seq(0, 1, length.out = 20)
  sh <- c(-0.05, 0, 0.05, 0.1)
  wf <- vapply(sh, function(s) tp - s, numeric(length(tp)))
  wp <- fck_warp_params(list(shifts = sh, warp_functions = wf, time_points = tp))
  expect_false(is.null(wp))
  expect_true(all(c("shift", "warp_amplitude") %in% names(wp)))
  expect_equal(wp$shift, sh)
  # RMS distance from identity is |shift| for a pure shift
  expect_equal(wp$warp_amplitude, abs(sh), tolerance = 1e-9)
})

test_that("a warp that moved nothing yields nothing to test", {
  tp <- seq(0, 1, length.out = 20)
  wf <- matrix(rep(tp, 4), ncol = 4)
  expect_null(fck_warp_params(list(shifts = rep(0, 4), warp_functions = wf, time_points = tp)))
  expect_null(fck_warp_params(NULL))
})

test_that("warp amplitude works for a nonlinear warp with no single shift", {
  tp <- seq(0, 1, length.out = 50)
  wf <- vapply(c(0.8, 1, 1.3), function(pw) tp^pw, numeric(length(tp)))
  wp <- fck_warp_params(list(warp_functions = wf, time_points = tp))
  expect_false(is.null(wp))
  expect_equal(names(wp), "warp_amplitude")     # no shift column exists
  expect_equal(wp$warp_amplitude[2], 0, tolerance = 1e-12)   # the identity warp
  expect_true(all(wp$warp_amplitude[c(1, 3)] > 0))
})

test_that("the warping parameters run through the same ANOVA machinery", {
  set.seed(4)
  g <- factor(rep(c("A", "B"), each = 40))
  wp <- data.frame(shift = c(rnorm(40, 0, 0.02), rnorm(40, 0.06, 0.02)))
  r <- fck_pc_anova_all(as.matrix(wp), g, across_pc_correction = "holm")
  expect_equal(r$k, 1)
  expect_lt(r$p_omnibus_adj[1], 1e-6)           # a planted phase difference
  expect_false(is.null(r$posthoc[[1]]))
})
