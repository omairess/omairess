# ==============================================================================
# tests/testthat/test-bounds.R — which bounds a fit sits on, and who hits two
#
# The convergence gate used to exclude every fit pinned to a parameter bound and
# report one number ("938 on a bound"). That is the wrong shape of answer twice
# over: it silently shrank the sample from 1305 to 367, and it did not say WHICH
# bound, so a badly chosen constraint was indistinguishable from a sample full
# of odd subjects.
#
# On the real Circaflex data the per-bound table says it in one line: 851 fits
# (65.2%) at tau's UPPER bound and 87 (6.7%) at its lower, A_sat never, and no
# fit on two at once. tau's ceiling there is t_max * 5 = 110 h -- a numerical
# guard, not a physiological claim -- so two thirds of the sample want a tau the
# constraint will not let them have, which is the identifiability problem
# stated as a fact rather than an inference.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

lo <- c(mesor = -Inf, A_sat = -Inf, tau = 0.5,  b_cos1 = -Inf, b_sin1 = -Inf)
up <- c(mesor =  Inf, A_sat =  Inf, tau = 110,  b_cos1 =  Inf, b_sin1 =  Inf)

test_that("a parameter exactly on a bound is detected, by name and side", {
  expect_equal(fck_bounds_hit(c(mesor = 20, A_sat = 30, tau = 110, b_cos1 = 1, b_sin1 = 2), lo, up),
               "tau (upper)")
  expect_equal(fck_bounds_hit(c(mesor = 20, A_sat = 30, tau = 0.5, b_cos1 = 1, b_sin1 = 2), lo, up),
               "tau (lower)")
})

test_that("an interior fit reports nothing", {
  expect_equal(fck_bounds_hit(c(mesor = 20, A_sat = 30, tau = 14, b_cos1 = 1, b_sin1 = 2), lo, up),
               character(0))
})

test_that("infinite bounds are never counted as hit", {
  # A_sat is unbounded here; a huge value is not a bound hit
  expect_false(any(grepl("A_sat",
    fck_bounds_hit(c(mesor = 20, A_sat = 1e9, tau = 14, b_cos1 = 1, b_sin1 = 2), lo, up))))
  expect_equal(fck_bounds_hit(c(a = 1e300), c(a = -Inf), c(a = Inf)), character(0))
})

test_that("both bounds on one parameter cannot both fire", {
  # a degenerate box where lower == upper: the value is at both, and both are
  # reported, because that IS what happened
  hit <- fck_bounds_hit(c(tau = 5), c(tau = 5), c(tau = 5))
  expect_equal(length(hit), 2)
})

test_that("several parameters can be pinned at once", {
  hit <- fck_bounds_hit(c(mesor = 0, A_sat = 500, tau = 110),
                        c(mesor = 0, A_sat = -Inf, tau = 0.5),
                        c(mesor = Inf, A_sat = 500, tau = 110))
  expect_equal(sort(hit), sort(c("mesor (lower)", "A_sat (upper)", "tau (upper)")))
})

test_that("the tolerance is relative, so a big and a small bound bite alike", {
  # tau = 110 * (1 + 1e-9) is on the bound; 110 * (1 + 1e-3) is not
  expect_equal(fck_bounds_hit(c(tau = 110 * (1 + 1e-9)), c(tau = 0.5), c(tau = 110)),
               "tau (upper)")
  expect_equal(fck_bounds_hit(c(tau = 110 * (1 + 1e-3)), c(tau = 0.5), c(tau = 110)),
               character(0))
  # the same sensitivity at the small end
  expect_equal(fck_bounds_hit(c(tau = 0.5 * (1 + 1e-9)), c(tau = 0.5), c(tau = 110)),
               "tau (lower)")
})

test_that("an unnamed or empty coefficient vector is handled", {
  expect_equal(fck_bounds_hit(c(1, 2, 3), lo, up), character(0))
  expect_equal(fck_bounds_hit(numeric(0), lo, up), character(0))
})

test_that("NA and non-finite coefficients never report a hit", {
  expect_equal(fck_bounds_hit(c(tau = NA_real_), c(tau = 0.5), c(tau = 110)), character(0))
  expect_equal(fck_bounds_hit(c(tau = Inf), c(tau = 0.5), c(tau = 110)), character(0))
})

# ------------------------------------------------------------------ summaries
mk <- function() list(
  character(0),
  "tau (upper)",
  "tau (upper)",
  c("tau (upper)", "A_sat (upper)"),
  "tau (lower)",
  c("tau (lower)", "mesor (lower)", "A_sat (upper)"))

test_that("the per-bound table counts every bound, not every fit", {
  s <- fck_bounds_summary(mk())
  expect_equal(s$n, 6)
  expect_equal(s$n_any, 5)                     # one fit is interior
  pb <- s$per_bound
  expect_equal(pb$n[pb$bound == "tau (upper)"], 3)
  expect_equal(pb$n[pb$bound == "A_sat (upper)"], 2)
  expect_equal(pb$n[pb$bound == "tau (lower)"], 2)
  # a fit hitting two bounds contributes to two rows, so these exceed n_any
  expect_gt(sum(pb$n), s$n_any)
  # sorted with the worst offender first
  expect_equal(pb$bound[1], "tau (upper)")
})

test_that("the per-count table partitions the fits exactly", {
  s <- fck_bounds_summary(mk())
  expect_equal(sum(s$per_count$n_subjects), s$n)          # every fit counted once
  expect_equal(sum(s$per_count$pct), 100, tolerance = 1e-9)
  # the fixture holds: 1 interior, 3 with one bound, 1 with two, 1 with three
  expect_equal(s$per_count$n_subjects[s$per_count$n_bounds == 0], 1)
  expect_equal(s$per_count$n_subjects[s$per_count$n_bounds == 1], 3)
  expect_equal(s$per_count$n_subjects[s$per_count$n_bounds == 2], 1)
  expect_equal(s$per_count$n_subjects[s$per_count$n_bounds == 3], 1)
})

test_that("fits on more than one bound are listed with which bounds", {
  s <- fck_bounds_summary(mk(), subject_ids = letters[1:6])
  expect_equal(s$n_multi, 2)
  expect_equal(nrow(s$multi), 2)
  expect_equal(s$multi$subject, c("d", "f"))
  expect_equal(s$multi$n_bounds, c(2L, 3L))
  expect_true(grepl("A_sat", s$multi$bounds[1]))
  expect_true(all(s$multi$n_bounds >= 2))
})

test_that("no fit on a bound gives an empty per-bound table, not an error", {
  s <- fck_bounds_summary(list(character(0), character(0)))
  expect_equal(s$n_any, 0)
  expect_null(s$per_bound)
  expect_null(s$multi)
  expect_equal(s$per_count$n_subjects[s$per_count$n_bounds == 0], 2)
})

test_that("subject ids fall back to the row number when none are supplied", {
  s <- fck_bounds_summary(mk())
  expect_equal(s$multi$subject, c("4", "6"))
  expect_equal(s$multi$row, c(4L, 6L))
})

test_that("the real Circaflex pattern reproduces from its counts", {
  # 851 at tau's ceiling, 87 at its floor, 367 interior, none on two
  bl <- c(replicate(851, "tau (upper)", simplify = FALSE),
          replicate(87, "tau (lower)", simplify = FALSE),
          replicate(367, character(0), simplify = FALSE))
  s <- fck_bounds_summary(bl)
  expect_equal(s$n, 1305)
  expect_equal(s$n_any, 938)
  expect_equal(s$n_multi, 0)
  expect_null(s$multi)
  expect_equal(s$per_bound$n, c(851, 87))
  expect_equal(s$per_bound$pct[1], 100 * 851 / 1305, tolerance = 1e-9)
  # tau is pinned for 71.9% of the sample and A_sat for none of it: the
  # constraint is on tau, and it is the identifiability problem, not the data
  expect_false(any(grepl("A_sat", s$per_bound$bound)))
})

# ------------------------------------------------------- what the gate does
test_that("including bound-hit fits keeps the whole converged sample", {
  conv  <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
  bound <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
  # the shipped rule: keep_rows <- conv & (include | !bound)
  keep_incl <- conv & (TRUE  | !bound)
  keep_excl <- conv & (FALSE | !bound)
  expect_equal(sum(keep_incl), 4)      # every converged fit
  expect_equal(sum(keep_excl), 2)      # only the interior ones
  # a non-converged fit is excluded either way: there is no solution to average
  expect_false(keep_incl[5])
  expect_false(keep_excl[5])
})
