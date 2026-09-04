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
# The palette itself is now the reporter's "Alternating Light/Dark Pairs" set:
# four hue families (blue, pink, green, orange), each with a light and a medium
# member. The MEDIUMS take slots 1-4 and the LIGHTS slots 5-8 -- taken in the
# order supplied, slots 1 and 2 would be the light and the medium blue, two
# series dE 14.9 apart that read as one group split in two. See the header of
# server/02b_helpers_palette.R for the full validator report.
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
  for (n in 2:8) {
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
test_that("the prefix is the four medium members, one per hue family", {
  expect_equal(length(FCK_GROUP_COLORS), 4)
  expect_true(all(grepl("^#[0-9a-fA-F]{6}$", FCK_GROUP_COLORS)))
  expect_equal(FCK_GROUP_COLORS, c("#00b0be", "#f45f74", "#98c127", "#ffb255"))
  # not the order they were listed in: that would put the two blues adjacent
  expect_false(identical(fck_group_ramp(2), c("#8fd7d7", "#00b0be")))
})

test_that("slots 5-8 are the light partners of slots 1-4", {
  expect_equal(FCK_GROUP_COLORS_EXTRA, c("#8fd7d7", "#ff8ca1", "#bdd373", "#ffcd8e"))
  # every hex the reporter supplied is somewhere in the first eight slots
  supplied <- c("#8fd7d7", "#00b0be", "#ff8ca1", "#f45f74",
                "#bdd373", "#98c127", "#ffcd8e", "#ffb255")
  expect_setequal(fck_group_ramp(8), supplied)
})

test_that("the cvd_safe variant re-steps only the green and the orange", {
  # switching the one constant must not disturb the blue or the pink
  expect_equal(.FCK_MEDIUMS$cvd_safe[1:2], .FCK_MEDIUMS$supplied[1:2])
  expect_false(identical(.FCK_MEDIUMS$cvd_safe[3], .FCK_MEDIUMS$supplied[3]))
  expect_false(identical(.FCK_MEDIUMS$cvd_safe[4], .FCK_MEDIUMS$supplied[4]))
  expect_equal(.FCK_MEDIUMS$cvd_safe, c("#00b0be", "#f45f74", "#9ec055", "#d38400"))
  # an unrecognised variant name falls back rather than erroring at load
  expect_true(FCK_PALETTE_VARIANT %in% names(.FCK_MEDIUMS))
})

test_that("dash is required past the validated depth, and only past it", {
  expect_false(fck_needs_dash(4))
  expect_true(fck_needs_dash(5))
  # eight distinct patterns, so slot k and its light partner at k+4 differ
  expect_equal(length(unique(FCK_GROUP_DASH)), 6)
  d <- fck_group_dashes(LETTERS[1:8])
  for (k in 1:4) expect_false(identical(unname(d[[k]]), unname(d[[k + 4]])), info = k)
})

test_that("the ramp extends past the palette rather than erroring", {
  expect_equal(length(fck_group_ramp(8)), 8)
  expect_equal(fck_group_ramp(3), FCK_GROUP_COLORS[1:3])
  expect_equal(length(fck_group_ramp(12)), 12)   # recycles, with dash carrying identity
})

test_that("components draw from the same vocabulary as groups", {
  expect_equal(fck_component_colors(4), fck_group_ramp(4))
  expect_equal(length(fck_component_colors(8)), 8)
})

test_that("emphasis and reference roles sit OUTSIDE the categorical palette", {
  # a red population-mean line would read as a fifth group next to the pink
  full <- fck_group_ramp(8)
  expect_false(FCK_EMPHASIS %in% full)
  expect_false(FCK_NEUTRAL %in% full)
  expect_equal(FCK_SERIES1, unname(FCK_GROUP_COLORS[1]))
  expect_equal(FCK_SERIES2, unname(FCK_GROUP_COLORS[2]))
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

# ============================================================== figure layout
# Source-level guards for four defects the reporter hit on screen. They are
# greps because the alternative is a headless browser: what is asserted is that
# the layout property that caused each overlap or reset is still set the way the
# fix set it.

test_that("the acrophase dial does not stack its radial axis on 12 o'clock", {
  src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
               collapse = "\n")
  blk <- sub(".*output\\$harmonic_polar_plot", "", src)
  blk <- substr(blk, 1, 12000)
  # angle = 90 put the radial title and its ticks through the topmost angular
  # label and up into the subtitle; 45 deg is midway between two 30 deg ticks
  expect_false(grepl("tickangle = 0, angle = 90", blk, fixed = TRUE))
  expect_true(grepl("angle = 45", blk, fixed = TRUE))
  # and the legend is under the dial, not eating into it from the right
  expect_true(grepl('legend = list(orientation = "h"', blk, fixed = TRUE))
  # the subtitle is a title line, not a paper-anchored annotation: as an
  # annotation it had no idea where the dial ended and the topmost angular
  # label -- drawn OUTSIDE the polar domain -- was written through it
  expect_true(grepl("<br><sub>%s</sub>", blk, fixed = TRUE))
  expect_false(grepl('x = 0.5, y = 0.925, xref = "paper"', blk, fixed = TRUE))
  expect_false(grepl('x = 0.5, y = 0.945, xref = "paper"', blk, fixed = TRUE))
})

test_that("both polar figures survive a re-render with the view intact", {
  for (f in c("server/72_harmonic.R", "server/74_polar_density.R")) {
    src <- paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
    expect_true(grepl("uirevision", src, fixed = TRUE), info = f)
  }
  # the density legend gets its OWN revision, keyed on which rings exist, so
  # hiding a group persists across a knob but not across a regrouping
  dsrc <- paste(readLines(file.path(app_dir, "server/74_polar_density.R"), warn = FALSE),
                collapse = "\n")
  expect_true(grepl('uirevision = "fck-density"', dsrc, fixed = TRUE))
  expect_true(grepl("uirevision = paste(nm, collapse", dsrc, fixed = TRUE))
})

test_that("the component-scores slider cannot promise more than the PCA kept", {
  src <- paste(readLines(file.path(app_dir, "server/40_fpca.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl('updateSliderInput(session, "effect_n_comp"', src, fixed = TRUE))
  expect_true(grepl("max = n_avail", src, fixed = TRUE))
  # and the figure names the remedy when the two still disagree
  expect_true(grepl("Number of components to extract", src, fixed = TRUE))
  # extraction default raised so the common case is not capped at three
  ui <- paste(readLines(file.path(app_dir, "ui/40_settings.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('"n_components", "Number of components to extract:",\n                         value = 5',
                    ui))
})

test_that("no unvalidated CSS colour names are left carrying series identity", {
  for (f in c("server/72_harmonic.R", "server/40_fpca.R")) {
    src <- paste(readLines(file.path(app_dir, f), warn = FALSE), collapse = "\n")
    for (dead in c("color = 'red'", 'color = "red"', "color = 'steelblue'",
                   "color = 'firebrick'", "color = 'blue'"))
      expect_false(grepl(dead, src, fixed = TRUE), info = paste(f, dead))
  }
  h <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl('c("green", "orange", "purple", "brown"', h, fixed = TRUE))
  expect_true(grepl("fck_component_colors(", h, fixed = TRUE))
})
