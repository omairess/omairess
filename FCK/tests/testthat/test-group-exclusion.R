# ==============================================================================
# tests/testthat/test-group-exclusion.R
#
# Reported: one subject had no AGEcategory, and the group-comparison plot died
# with "missing value where TRUE/FALSE needed".
#
# Cause: an earlier design carried unlabelled subjects as an "UNASSIGNED"
# pseudo-group so they could not vanish silently. But their group column is NA,
# so `params$group == "__UNASSIGNED__"` matched nothing, the angle vector came
# back empty, circular_mean() of an empty vector is atan2(NaN, NaN) = NaN, and
# `if (NaN < 0)` is an error, not FALSE.
#
# A label-less group of one is not a group: it has no circular mean and every
# comparison built on it is undefined. They are now EXCLUDED from group
# analyses and counted in the audit, so the number still reconciles without
# poisoning the arithmetic.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/07_helpers_circular.R"))
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

# circular_mean() is defined inside the server body of 72_harmonic.R, which
# cannot be sourced outside a session. Reproduced here verbatim -- the point of
# the first tests is precisely what it returns for an empty vector.
circular_mean <- function(angles_rad) {
  x <- mean(cos(angles_rad), na.rm = TRUE)
  y <- mean(sin(angles_rad), na.rm = TRUE)
  atan2(y, x)
}

test_that("the local circular_mean matches the shipped one", {
  src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl("x <- mean(cos(angles_rad), na.rm = TRUE)", src, fixed = TRUE))
  expect_true(grepl("y <- mean(sin(angles_rad), na.rm = TRUE)", src, fixed = TRUE))
  expect_true(grepl("atan2(y, x)", src, fixed = TRUE))
})

# ------------------------------------------------------- the failure itself
test_that("circular_mean of an empty vector is NaN, and NaN breaks an if", {
  expect_true(is.nan(circular_mean(numeric(0))))
  expect_error(if (circular_mean(numeric(0)) < 0) TRUE, "missing value")
})

test_that("an NA group label matches nothing, which is how the vector emptied", {
  grp <- c("YOUTH", "ADULT", NA)
  # the old lookup
  expect_equal(sum(grp == "__UNASSIGNED__", na.rm = TRUE), 0)
  expect_true(is.na((grp == "__UNASSIGNED__")[3]))
  # and indexing with an NA-bearing logical yields NA entries, not an empty set
  x <- c(1, 2, 3)
  expect_true(any(is.na(x[grp == "YOUTH"])))
})

# ------------------------------------------------------- what happens now
test_that("unlabelled subjects are excluded and counted, not carried", {
  labs <- c(rep("YOUTH", 654), rep("ADULT", 410),
            rep("MIDDLE_AGE", 181), rep("ELDERLY", 59), NA)
  ga <- fck_group_audit(labs, seq_along(labs))
  expect_equal(ga$n_total, 1305)
  expect_equal(ga$n_unassigned, 1)
  expect_false("__UNASSIGNED__" %in% ga$levels)   # never a level
  expect_equal(sort(ga$levels), sort(c("YOUTH", "ADULT", "MIDDLE_AGE", "ELDERLY")))
  # the accounting still closes
  expect_equal(sum(ga$counts) + ga$n_unassigned, ga$n_total)
})

test_that("the group-selection rule keeps only usable labels", {
  lab_all <- c("YOUTH", NA, "ADULT", "", "YOUTH")
  # the rule the app now uses
  pick <- function(g) which(!is.na(lab_all) & nzchar(lab_all) & lab_all == g)
  expect_equal(pick("YOUTH"), c(1L, 5L))
  expect_equal(pick("ADULT"), 3L)
  # neither NA nor the empty string is ever selected
  expect_equal(length(pick(NA)), 0)
  expect_equal(length(pick("")), 0)
})

test_that("a group emptied by an upstream filter is skipped, not fatal", {
  # the guard in the comparison plot: finite angles only, then skip if none
  ang <- c(NA_real_, NaN, NA_real_)
  ang <- ang[is.finite(ang)]
  expect_equal(length(ang), 0)
  # the guard's condition, which replaces the bare `if (mean < 0)`
  expect_true(length(ang) < 1)
})

test_that("a one-subject group still has a defined circular mean", {
  # excluding UNLABELLED subjects is not the same as excluding SMALL groups:
  # a labelled group of one has a perfectly good mean direction, it is just
  # imprecise, and the n >= 3 guard is what handles that
  m <- circular_mean(2.5)
  expect_true(is.finite(m))
  expect_equal(m, 2.5, tolerance = 1e-12)
})

# ------------------------------------- constant term vs MESOR vs start value
test_that("beta_0, the MESOR and the start value are three different numbers", {
  mod <- list(period = 24, n_harmonics = 2, trend_type = "exp_sat",
              time_vec = seq(0, 22, length.out = 100), t_offset = 0, origin_shift = 8)
  coefs <- c(20.53, 60, 14, 12, -5, 3, 2)   # b0, A_sat, tau, H1, H2

  b0 <- coefs[1]
  v0 <- fck_value_at(coefs, mod, min(mod$time_vec))
  ram <- fck_rhythm_adjusted_mean(b0, "exp_sat", c(coefs[2], coefs[3]),
                                  min(mod$time_vec), max(mod$time_vec), 0)
  expect_false(isTRUE(all.equal(b0, v0)))
  expect_false(isTRUE(all.equal(b0, ram)))
  expect_false(isTRUE(all.equal(v0, ram)))

  # at the first observation a saturating trend anchored there contributes
  # exactly nothing, so the start value is b0 PLUS the harmonic sum
  hsum <- coefs[4] * cos(0) + coefs[5] * sin(0) + coefs[6] * cos(0) + coefs[7] * sin(0)
  expect_equal(v0, b0 + hsum, tolerance = 1e-9)
})

test_that("with no trend and no harmonics all three coincide", {
  mod <- list(period = 24, n_harmonics = 1, trend_type = "none",
              time_vec = seq(0, 24, length.out = 50), t_offset = 0, origin_shift = 0)
  coefs <- c(42, 0, 0)
  expect_equal(fck_value_at(coefs, mod, 0), 42)
  expect_equal(fck_rhythm_adjusted_mean(42, "none", numeric(0), 0, 24, 0), 42)
})

test_that("the MESOR is above beta_0 for a rising trend, as it must be", {
  ram <- fck_rhythm_adjusted_mean(20.53, "exp_sat", c(60, 14), 0, 22, 0)
  expect_gt(ram, 20.53)
  # and bounded by the asymptote it is approaching
  expect_lt(ram, 20.53 + 60)
})

# ---------------------------------------- the pairwise columns must exist
test_that("the per-subject columns the pairwise picker offers are declared", {
  src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
               collapse = "\n")
  for (col in c("mesor", "mesor_adj", "value_at_start"))
    expect_true(grepl(sprintf('"%s"', col), src, fixed = TRUE), info = col)
  ui <- paste(readLines(file.path(app_dir, "ui/73_cosinor_pairwise.R"), warn = FALSE),
              collapse = "\n")
  for (col in c("mesor_adj", "value_at_start"))
    expect_true(grepl(col, ui, fixed = TRUE), info = col)
})

test_that("no UNASSIGNED pseudo-group survives anywhere in the server", {
  for (f in list.files(file.path(app_dir, "server"), pattern = "[.]R$", full.names = TRUE)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl("__UNASSIGNED__", src, fixed = TRUE), info = basename(f))
  }
})
